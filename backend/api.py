# backend/api.py
from pathlib import Path
import os, sys, tempfile, subprocess, shutil, re, traceback, glob, shlex, uuid, gzip, io
from typing import List, Dict, Optional
import boto3
from botocore.exceptions import ClientError

import json
import queue
import threading
import asyncio
import contextlib

from fastapi.responses import StreamingResponse
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from pydantic import BaseModel, Field

import io
import argparse
from ase.io import read as ase_read, write as ase_write

# Import directly from your installed miniCAT library!
from minicat.main import (
    _sanitize_simple_acid_smiles, smiles_to_3d_mol, 
    build_per_site_variant, execute_passivation_job
)
from minicat.functional_groups_class import get_fg_registry, detect_fg_matches_neutral

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
class LineQueueWriter(io.TextIOBase):
    def __init__(self, q): self.q, self.buf = q, ""
    def write(self, s):
        self.buf += s
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1)
            self.q.put(line)
        return len(s)
    def flush(self):
        if self.buf: self.q.put(self.buf); self.buf = ""

# -------------------------------------------------------------------
# Attach endpoint (Upgraded for Real-Time Streaming)
# -------------------------------------------------------------------
class LineQueueWriter(io.TextIOBase):
    def __init__(self, q): self.q, self.buf = q, ""
    def write(self, s):
        self.buf += s
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1)
            self.q.put(line)
        return len(s)
    def flush(self):
        if self.buf: self.q.put(self.buf); self.buf = ""

@app.post("/attach")
async def attach(payload: Dict):
    # 1. Parse request schema
    try:
        new = MiniCATRequest(**payload)
        xyztext, out_prefix, req_jobs = new.xyztext, new.out_prefix, new.jobs
    except Exception:
        try:
            old = LegacyAttachRequest(**payload)
            xyztext, mode = old.xyztext, "random" if old.split else "segmented"
            req_jobs = [Job(ligands=[old.smiles], dummy="Cl", dist=f"1.0:{mode}")]
            out_prefix = "final_passivated_dot"
        except Exception:
            raise HTTPException(status_code=400, detail="Invalid request body")

    async def _stream():
        # Let the frontend know we started successfully
        yield json.dumps({"event": "status", "line": "accepted"}) + "\n"
        
        log_q = queue.Queue()
        qw = LineQueueWriter(log_q)
        res_box = {}

        def _thread_run():
            try:
                # Redirect all miniCAT print() statements to our queue
                with contextlib.redirect_stdout(qw), contextlib.redirect_stderr(qw):
                    print("[cap] Starting in-memory passivation engine...")
                    args = argparse.Namespace(
                        seed=1337, ff="uff", refinement_passes=3, neigh=8,
                        coarse_step_deg=20.0, sterics_mode="vdw", sterics_margin=0.25,
                        warn_tol=1.4, adaptive_offset_steps=4, adaptive_offset_step=0.15,
                        anchor_mode="dummy", offset_out=0.0, bond_len_override=None,
                        neighbor_repulsion=0.5
                    )
                    
                    reg = get_fg_registry()
                    prepared_jobs = []
                    
                    # Pre-compute ligands
                    for req_job in req_jobs:
                        dist_parts = req_job.dist.split(':')
                        strategy, ratios = dist_parts[-1], [float(r) for r in dist_parts[:-1]]
                        
                        precomputed_ligands = []
                        role = None
                        for smiles in req_job.ligands:
                            print(f"[cap] Building 3D geometry for {smiles}...")
                            mol_neutral = smiles_to_3d_mol(_sanitize_simple_acid_smiles(smiles), args.seed, args.ff)
                            fg_matches = detect_fg_matches_neutral(mol_neutral, 'anion', reg) or detect_fg_matches_neutral(mol_neutral, 'cation', reg)
                            
                            if not fg_matches:
                                print(f"[error] No recognizable group for {smiles}")
                                raise ValueError(f"No recognizable group for {smiles}")
                            
                            ionic_mol, (fg_active, match_active) = build_per_site_variant(fg_matches, mol_neutral, 0, args.ff)
                            precomputed_ligands.append({"mol": ionic_mol, "fg": fg_active, "match": match_active})
                            if role is None:
                                role = fg_matches[0][0].role

                        prepared_jobs.append({
                            "smiles": req_job.ligands, "dummy": req_job.dummy, "role": role,
                            "ratios": ratios, "strategy": strategy, "precomputed_ligands": precomputed_ligands
                        })

                    # Run Passivation
                    current_qd = ase_read(io.StringIO(xyztext), format="xyz")
                    for job in prepared_jobs:
                        print(f"[cap] Passivating {job['dummy']} sites based on '{job['strategy']}' distribution...")
                        current_qd = execute_passivation_job(current_qd, job, args)
                        
                    out_xyz = io.StringIO()
                    ase_write(out_xyz, current_qd, format="xyz")
                    res_box['final_str'] = out_xyz.getvalue()
                    print("[cap] Passivation complete!")
                    
            except Exception as e:
                import traceback
                print(f"[fatal] {traceback.format_exc()}", file=qw)
                res_box['error'] = str(e)
            finally:
                qw.flush()
                log_q.put(None)

        # Spin up the thread
        t = threading.Thread(target=_thread_run)
        t.start()

        # Stream lines to the frontend as soon as they appear
        while True:
            try:
                line = log_q.get_nowait()
                if line is None: break
                if line.strip(): yield json.dumps({"event": "log", "line": line.strip()}) + "\n"
            except queue.Empty:
                await asyncio.sleep(0.01)

        t.join()

        # Send the final geometry back in the same structure the Library expects
        if 'error' in res_box:
            yield json.dumps({"event": "result", "status": "failed", "error": res_box['error']}) + "\n"
        else:
            yield json.dumps({
                "event": "result",
                "status": "success",
                "message": "miniCAT (In-Memory) OK",
                "results": [{"filename": f"{out_prefix}_final.xyz", "xyz": res_box['final_str']}],
                "cmd": "Executed via direct memory import",
                "stdout": "Success",
                "stderr": ""
            }) + "\n"

    return StreamingResponse(
        _stream(),
        media_type="text/plain; charset=utf-8",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"}
    )

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
