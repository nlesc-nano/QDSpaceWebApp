# backend/api.py
from pathlib import Path
import os, sys, tempfile, subprocess, shutil, re, traceback, glob, shlex, uuid, gzip, io
from typing import List, Dict, Optional
import boto3
from botocore.exceptions import ClientError

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from pydantic import BaseModel, Field

# -------------------------------------------------------------------
# S3 + Config
# -------------------------------------------------------------------
PROPS_S3_BUCKET = os.environ.get("PROPS_S3_BUCKET", "qdspace-frontend-v1")
PROPS_S3_PREFIX = os.environ.get("PROPS_S3_PREFIX", "library")
AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")
_s3 = boto3.client("s3", region_name=AWS_REGION)

TMP_PLOTS_DIR = Path("/tmp/plots")
TMP_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
def _safe_folder(folder: str) -> str:
    # Only allow A–Z, a–z, 0–9, slash, underscore, dash, and dot
    return re.sub(r"[^A-Za-z0-9/_\.\-]", "_", (folder or "")).strip("/")


def _s3_exists(bucket: str, key: str) -> bool:
    try:
        _s3.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError as e:
        code = (e.response or {}).get("Error", {}).get("Code", "")
        if code in ("404", "NotFound", "NoSuchKey"):
            return False
        raise

def _s3_get_bytes(bucket: str, key: str) -> bytes:
    obj = _s3.get_object(Bucket=bucket, Key=key)
    return obj["Body"].read()

def _presign_get(bucket: str, key: str, seconds: int = 3600) -> str:
    return _s3.generate_presigned_url(
        "get_object", Params={"Bucket": bucket, "Key": key}, ExpiresIn=seconds
    )

# -------------------------------------------------------------------
# FastAPI app
# -------------------------------------------------------------------
app = FastAPI(title="QDSpace Backend")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
app.add_middleware(GZipMiddleware, minimum_size=1024)

# -------------------------------------------------------------------
# Models
# -------------------------------------------------------------------
class PlotRequest(BaseModel):
    folder: str
    fuzzy: str
    pdos: str
    coop: str
    out: Optional[str] = None
    ef: Optional[float] = None
    title: Optional[str] = None
    normalize_coop: bool = True

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
# -------------------------------------------------------------------
# Root health
# -------------------------------------------------------------------
@app.get("/")
def root():
    return {"message": "QDSpace backend is alive"}

# -------------------------------------------------------------------
# Attach endpoint
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Plot endpoint
# -------------------------------------------------------------------
@app.post("/plot")
def plot_interactive(req: PlotRequest):
    """
    Compute-once / reuse-forever plot:
      • Looks for deterministic files in …/library/<folder>/properties/: plot.html + plot.png
      • If HTML exists -> return it (fast path).
      • If missing -> download inputs, render both PNG and HTML once, upload to that folder.
    """
    folder = req.folder or ""
    safe_folder = _safe_folder(folder)
    base_props_key = f"{PROPS_S3_PREFIX.rstrip('/')}/{safe_folder}/properties"

    fuzzy_key = f"{base_props_key}/{req.fuzzy or 'fuzzy_data.npz'}"
    pdos_key  = f"{base_props_key}/{req.pdos  or 'pdos_data.csv'}"
    coop_key  = f"{base_props_key}/{req.coop  or 'coop_data.csv'}"

    out_html_key = f"{base_props_key}/plot.html"
    out_png_key  = f"{base_props_key}/plot.png"

    cdn_base = os.environ.get("FRONTEND_CDN", "").rstrip("/")
    cdn_html = f"{cdn_base}/{out_html_key}" if cdn_base else None

    # --- 1. Fast path: if HTML already exists, just return it
    if _s3_exists(PROPS_S3_BUCKET, out_html_key):
        href = cdn_html or _presign_get(PROPS_S3_BUCKET, out_html_key, 3600)
        return {"message": "OK (cache)", "href": href, "key": out_html_key}

    # --- 2. Validate that required inputs exist
    have_fuzzy = _s3_exists(PROPS_S3_BUCKET, fuzzy_key)
    have_pdos  = _s3_exists(PROPS_S3_BUCKET,  pdos_key)
    have_coop  = _s3_exists(PROPS_S3_BUCKET,  coop_key)

    # --- 2. Validate that required inputs exist before build
    missing_inputs = []
    if not have_fuzzy: missing_inputs.append(fuzzy_key)
    if not have_pdos:  missing_inputs.append(pdos_key)
    if not have_coop:  missing_inputs.append(coop_key)

    if missing_inputs:
        raise HTTPException(
            status_code=400,
            detail={
                "code": "PLOT_INPUTS_MISSING",
                "missing": missing_inputs,
                "folder": base_props_key,
            },
        )

    # --- 3. Slow path: build and upload
    try:
        TMP_PLOTS_DIR.mkdir(parents=True, exist_ok=True)
        fuzzy_local = TMP_PLOTS_DIR / "fuzzy_data.npz"
        pdos_local  = TMP_PLOTS_DIR / "pdos_data.csv"
        coop_local  = TMP_PLOTS_DIR / "coop_data.csv"
        fuzzy_local.write_bytes(_s3_get_bytes(PROPS_S3_BUCKET, fuzzy_key))
        pdos_local.write_bytes(_s3_get_bytes(PROPS_S3_BUCKET,  pdos_key))
        coop_local.write_bytes(_s3_get_bytes(PROPS_S3_BUCKET,  coop_key))

        local_html = TMP_PLOTS_DIR / "plot.html"
        local_png  = TMP_PLOTS_DIR / "plot.png"
        script = Path(__file__).resolve().parent / "plot_interactive.py"
        if not script.is_file():
            raise HTTPException(status_code=500, detail="plot_interactive.py not found")

        cmd = [
            sys.executable, str(script),
            "--fuzzy", str(fuzzy_local),
            "--pdos",  str(pdos_local),
            "--coop",  str(coop_local),
            "--out",   str(local_html),
            "--png",   str(local_png),
        ]
        if req.normalize_coop: cmd.append("--normalize-coop")
        if req.ef is not None: cmd += ["--ef", str(req.ef)]
        if req.title:          cmd += ["--title", req.title]

        subprocess.run(cmd, cwd=str(TMP_PLOTS_DIR), check=True, capture_output=True, text=True)

        if not local_html.is_file():
            raise HTTPException(status_code=500, detail="Plot build failed: HTML not produced")
  
        # Upload PNG first (so UI can show preview)
        if local_png.is_file():
            with open(local_png, "rb") as f:
                _s3.put_object(
                    Bucket=PROPS_S3_BUCKET,
                    Key=out_png_key,
                    Body=f,
                    ContentType="image/png",
                    CacheControl="public, max-age=31536000, immutable",
                )

        # Upload gzipped HTML next
        with open(local_html, "rb") as f:
            raw = f.read()
        gzbuf = io.BytesIO()
        with gzip.GzipFile(fileobj=gzbuf, mode="wb") as gz:
            gz.write(raw)
        gzbuf.seek(0)
        _s3.put_object(
            Bucket=PROPS_S3_BUCKET,
            Key=out_html_key,
            Body=gzbuf,
            ContentType="text/html; charset=utf-8",
            ContentEncoding="gzip",
            CacheControl="public, max-age=31536000, immutable",
        )

        href = cdn_html or _presign_get(PROPS_S3_BUCKET, out_html_key, 3600)
        return {"message": "OK", "href": href, "key": out_html_key}

    except subprocess.CalledProcessError as e:
        err = (e.stderr or e.stdout or "").strip()
        raise HTTPException(status_code=500, detail={"code": "PLOT_BUILD_FAILED", "error": err})
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail={"code": "PLOT_BUILD_FAILED", "error": str(e)})
