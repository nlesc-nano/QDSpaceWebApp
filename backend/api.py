# backend/api.py
from pathlib import Path
import os, glob, shlex, shutil, tempfile, subprocess, sys, uuid
from typing import List, Dict, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field

# ---------- Paths ----------
PROPS_ROOT = Path(
    os.environ.get(
        "PROPS_ROOT",
        Path(__file__).resolve().parent.parent / "docs" / "library"
    )
).resolve()

PROPS_SUBDIR = os.environ.get("PROPS_SUBDIR", "properties")

# Temp output directory for large HTML plots (served via /api/plot_file)
TMP_PLOTS_DIR = Path(__file__).resolve().parent / "_tmp_plots"
TMP_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# ---------- App ----------
app = FastAPI(title="miniCAT backend")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
# gzip big HTML responses
app.add_middleware(GZipMiddleware, minimum_size=1024)

# ---------- Models ----------
class Job(BaseModel):
    ligands: List[str] = Field(..., min_items=1)
    dummy: str
    dist: str  # "r1:r2:...:mode" where mode ∈ {random, segmented}

class MiniCATRequest(BaseModel):
    xyztext: str
    out_prefix: str = "final_passivated_dot"
    jobs: List[Job] = Field(..., min_items=1)

class LegacyAttachRequest(BaseModel):
    xyztext: str
    smiles: str
    split: bool = True  # split→random, not split→segmented

class PlotRequest(BaseModel):
    folder: str                 # e.g. "II-VI/CdTe/HLE17/12ang"
    fuzzy: str                  # e.g. "fuzzy_data.npz"
    pdos: str                   # e.g. "pdos_data.csv"
    coop: str                   # e.g. "coop_data.csv"
    out: Optional[str] = None   # ignored; backend uses temp file
    ef: Optional[float] = None
    title: Optional[str] = None
    normalize_coop: bool = True

# ---------- Health ----------
@app.get("/")
def root():
    return {"message": "miniCAT backend is alive. POST /attach"}

# ---------- miniCAT attach ----------
@app.post("/attach")
def attach(payload: Dict):
    # Parse new schema first
    req_jobs: List[Job]
    out_prefix = "final_passivated_dot"
    try:
        new = MiniCATRequest(**payload)
        xyztext = new.xyztext
        out_prefix = new.out_prefix
        req_jobs = new.jobs
    except Exception:
        # Fallback legacy: 1 ligand, single ratio
        try:
            old = LegacyAttachRequest(**payload)
        except Exception:
            raise HTTPException(status_code=400, detail="Invalid request body")
        xyztext = old.xyztext
        mode = "random" if old.split else "segmented"
        dist = f"1.0:{mode}"
        req_jobs = [Job(ligands=[old.smiles], dummy="Cl", dist=dist)]

    tmp = tempfile.mkdtemp(prefix="minicat_")
    try:
        core = os.path.join(tmp, "initial_dot.xyz")
        with open(core, "w") as f:
            f.write(xyztext)

        cmd = ["miniCAT", "--qd", "initial_dot.xyz", "--out_prefix", out_prefix]
        for j in req_jobs:
            cmd += ["--job-ligands", *j.ligands, "--job-dummy", j.dummy, "--job-dist", j.dist]
        cmd_str = " ".join(shlex.quote(c) for c in cmd)

        try:
            proc = subprocess.run(cmd, cwd=tmp, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            raise HTTPException(status_code=500, detail="miniCAT executable not found on PATH")
        except subprocess.CalledProcessError as e:
            raise HTTPException(status_code=500, detail=e.stderr or e.stdout or "miniCAT failed")

        outs = sorted(glob.glob(os.path.join(tmp, f"{out_prefix}*.xyz")))
        if not outs:
            raise HTTPException(status_code=500, detail="miniCAT produced no .xyz files")

        results = []
        for p in outs:
            with open(p, "r") as f:
                results.append({"filename": os.path.basename(p), "xyz": f.read()})

        return {
            "message": f"miniCAT OK ({len(results)} file(s))",
            "results": results,
            "cmd": cmd_str,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        }
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

# ---------- Properties plotting ----------
@app.post("/plot")
def plot_interactive(req: PlotRequest):
    # Resolve target folder safely under PROPS_ROOT
    base = (PROPS_ROOT / (req.folder or "")).resolve()
    if not str(base).startswith(str(PROPS_ROOT)):
        raise HTTPException(status_code=400, detail="Invalid folder path")

    workdir = base / PROPS_SUBDIR
    if not workdir.is_dir():
        raise HTTPException(status_code=400, detail=f"'properties' folder not found: {workdir}")

    # Absolute input paths (script will read these)
    fuzzy_path = (workdir / req.fuzzy).resolve()
    pdos_path  = (workdir / req.pdos ).resolve()
    coop_path  = (workdir / req.coop ).resolve()
    for p in (fuzzy_path, pdos_path, coop_path):
        if not p.is_file():
            raise HTTPException(status_code=400, detail=f"Missing input file: {p}")

    # Always use the backend plotter (single source of truth)
    script = Path(__file__).resolve().parent / "plot_interactive.py"
    if not script.is_file():
        raise HTTPException(status_code=500, detail="backend/plot_interactive.py not found")

    # Unique temp filename to avoid caching & collisions
    fid = f"{uuid.uuid4().hex}.html"
    out_path = (TMP_PLOTS_DIR / fid).resolve()

    # Build command; run in TMP_PLOTS_DIR (keeps datasets read-only)
    cmd = [
        sys.executable, str(script),
        "--fuzzy", str(fuzzy_path),
        "--pdos",  str(pdos_path),
        "--coop",  str(coop_path),
        "--out",   str(out_path),
    ]
    if req.normalize_coop:
        cmd.append("--normalize-coop")
    if req.ef is not None:
        cmd += ["--ef", str(req.ef)]
    if req.title:
        cmd += ["--title", req.title]

    try:
        run = subprocess.run(
            cmd,
            cwd=str(TMP_PLOTS_DIR),
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        err = (e.stderr or e.stdout or "").strip()
        raise HTTPException(status_code=500, detail=f"plot_interactive failed: {err}")

    if not out_path.is_file():
        raise HTTPException(status_code=500, detail=f"Output HTML not found: {out_path}")

    # Return a link to embed in an <iframe>
    href = f"/api/plot_file?fid={fid}"
    return {
        "message": "OK",
        "href": href,
        "cmd": " ".join(cmd),
        "stdout": run.stdout or "",
        "stderr": run.stderr or "",
    }

@app.get("/plot_file", response_class=HTMLResponse)
def plot_file(fid: str):
    # Safe lookup inside temp dir
    fp = (TMP_PLOTS_DIR / fid).resolve()
    if not str(fp).startswith(str(TMP_PLOTS_DIR)) or not fp.is_file():
        raise HTTPException(status_code=404, detail="Plot HTML not found")
    html = fp.read_text(encoding="utf-8")
    # Avoid browser caching
    return HTMLResponse(content=html, headers={"Cache-Control": "no-store"})

