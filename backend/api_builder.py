from __future__ import annotations

import io
import json, logging, os, re, shutil, subprocess, sys, tempfile

import yaml
import traceback

import pty, tty, fcntl, select
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import threading
import queue
import asyncio
import contextlib
from builder.main import main as nc_builder_main

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from pydantic import BaseModel, Field

import argparse
from ase.io import read as ase_read, write as ase_write
from minicat.main import (
    _sanitize_simple_acid_smiles, smiles_to_3d_mol, 
    build_per_site_variant, execute_passivation_job
)
from minicat.functional_groups_class import get_fg_registry, detect_fg_matches_neutral

import numpy as np

# Core inorganic elements for each material (ignores organic ligands like MA/FA)
MATERIAL_ELEMENTS = {
    "CSPBCL3": ["Cs", "Pb", "Cl"], "CSPBBR3": ["Cs", "Pb", "Br"], "CSPBI3":  ["Cs", "Pb", "I"],
    "MAPBI3":  ["Pb", "I"], "FAPBI3":  ["Pb", "I"], 
    "ZNS":     ["Zn", "S"], "ZNSE":    ["Zn", "Se"], "ZNTE":    ["Zn", "Te"],
    "CDS":     ["Cd", "S"], "CDSE":    ["Cd", "Se"], "CDTE":    ["Cd", "Te"],
    "HGS":     ["Hg", "S"], "HGSE":    ["Hg", "Se"], "HGTE":    ["Hg", "Te"],
    "ALP":     ["Al", "P"], "ALAS":    ["Al", "As"], "ALSB":    ["Al", "Sb"],
    "GAP":     ["Ga", "P"], "GAAS":    ["Ga", "As"], "GASB":    ["Ga", "Sb"],
    "INP":     ["In", "P"], "INAS":    ["In", "As"], "INSB":    ["In", "Sb"],
    "PBS":     ["Pb", "S"], "PBSE":    ["Pb", "Se"]
}

def get_cluster_size_metrics(coords_ang, atom_symbols=None, material_name=None):
    """Calculates size using Axis-Aligned Bounding Box (ideal for lattice-cut QDs)."""
    coords = np.asarray(coords_ang, dtype=float)

    if atom_symbols is not None and material_name is not None:
        m_name = material_name.upper()
        if 'MATERIAL_ELEMENTS' in globals() and m_name in MATERIAL_ELEMENTS:
            core_elements = [el.lower() for el in MATERIAL_ELEMENTS[m_name]]
            core_coords = [coords[i] for i, sym in enumerate(atom_symbols) if sym.lower() in core_elements]
            if len(core_coords) > 0:
                coords = np.array(core_coords)

    if len(coords) < 2:
        return {'R_eff_hull': 1.0, 'diameter_hull': 2.0}

    # Blazing fast: Calculate peak-to-peak distance directly along X, Y, and Z axes
    spans = np.ptp(coords, axis=0)

    # Add 2.5 Å for the physical outer electron cloud (van der Waals radii)
    spans += 2.5

    # The effective diameter is the average of length, width, and height
    avg_diameter = np.mean(spans)
    R_eff = avg_diameter / 2.0

    return {'R_eff_hull': float(R_eff), 'diameter_hull': float(avg_diameter)}

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ---------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------

class Facet(BaseModel):
    hkl: str
    gamma: float


class ShellLayer(BaseModel):
    material_cif: str
    aspect: List[float] = [1.0, 1.0, 1.0]
    facets: List[Facet] = []


class LigJob(BaseModel):
    smiles: str
    ratio: float = 1.0
    dummy: str = ""

class BuildOptions(BaseModel):
    radius_A: float
    core_cif_filename: str
    aspect: List[float] = [1.0, 1.0, 1.0]
    facets: List[Facet] = []
    shells: List[ShellLayer] = []

    # legacy (kept for compatibility elsewhere)
    ligand: Optional[str] = None
    surf_tol: Optional[float] = None

    # SMILES passivation
    cap_distribution: Optional[str] = "uniform"  # 'uniform' | 'segmented' | 'random'
    cap_anionic_jobs: Optional[List[LigJob]] = None  # replace Cl
    cap_cationic_jobs: Optional[List[LigJob]] = None  # replace Rb
    skip_core_build: bool = False
    xyz_unpassivated: Optional[str] = None

# ---------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------

app = FastAPI(
    title="QD_Builder API (Auto-Passivation)",
    version="5.3.6",
    description="Builds nanocrystals with automatic ligand selection for charge neutrality.",
)

allowed_origins = [
    os.getenv("FRONTEND_CDN", "https://dmq59n0f96mxz.cloudfront.net"),
    "http://localhost:8000",
    "http://127.0.0.1:8000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------------------------
# Static
# ---------------------------------------------------------------------


@app.get("/", response_class=HTMLResponse)
async def read_index():
    index_path = Path(__file__).resolve().parent.parent / "docs" / "builder" / "index.html" 
    if not index_path.exists():
        return HTMLResponse(
            content="<h1>index.html not found</h1><p>Make sure the index.html file is in the same directory as api.py.</p>",
            status_code=404,
        )
    return HTMLResponse(content=index_path.read_text())


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

LIG_PLACEHOLDERS = {"Cl", "Rb"}
LOG_KEEP = 300
_LOG_RE = re.compile(r"(Build(ing)?|Wulff|facet|shell|core|Atoms?|Total charge|Export|Passivation|Done|Warn|Error)$", re.I)
_SKIP_RE = re.compile(r"(downgrade|CN=|unique,\s*depth=|q=\d|neighbor)", re.I)


def _filter_log(text: str) -> str:
    out: List[str] = []
    last = None
    for raw in text.splitlines():
        line = raw.strip()
        if _SKIP_RE.search(line):
            continue
        if not _LOG_RE.search(line):
            continue
        if len(line) > 140:
            line = line[:140] + " …"
        if line == last:
            continue
        out.append(line)
        last = line
    return "\n".join(out[-LOG_KEEP:])

def _xyz_count_symbols(xyz_text: str, targets=("Cl","Rb")) -> dict[str,int]:
    counts = {t: 0 for t in targets}
    lines = xyz_text.splitlines()
    if len(lines) < 3:
        return counts
    try:
        n = int(lines[0].strip())
    except Exception:
        n = len(lines) - 2
    start = 2
    end = min(len(lines), start + max(0, n))
    for i in range(start, end):
        parts = (lines[i] or "").split()
        if parts:
            sym = parts[0]
            if sym in counts:
                counts[sym] += 1
    return counts

def _allocate_counts_by_ratio(total: int, jobs) -> dict[str, int]:
    """Split `total` across jobs by their ratios (largest-remainder)."""
    items = []
    for j in jobs or []:
        if isinstance(j, dict):
            smiles = (j.get("smiles") or "").strip()
            ratio  = float(j.get("ratio") or 0)
        else:
            smiles = (getattr(j, "smiles", "") or "").strip()
            ratio  = float(getattr(j, "ratio", 0) or 0)
        if smiles:
            items.append([smiles, max(0.0, ratio)])
    if not items or total <= 0:
        return {}
    sw = sum(w for _, w in items)
    if sw <= 0:
        base, rem = divmod(total, len(items))
        out = {s: base for s,_ in items}
        for i in range(rem): out[items[i%len(items)][0]] += 1
        return out
    # proportional + largest remainder
    prov, acc = [], 0
    for s,w in items:
        exact = total * (w / sw)
        f = int(exact // 1)
        acc += f
        prov.append([s, f, exact - f])
    rem = total - acc
    prov.sort(key=lambda x: x[2], reverse=True)
    for i in range(rem): prov[i][1] += 1
    return {s:c for s,c,_ in prov}


def _count_grouped_xyz(xyz_text: str, core_elems: set[str], shell_elems: set[str]) -> dict:
    counts = {"core": Counter(), "shell": Counter(), "ligand": Counter()}
    lines = xyz_text.splitlines()
    try:
        n = int(lines[0].strip())
    except Exception:
        n = 0

    for i in range(2, 2 + n):
        if i >= len(lines):
            break
        el = (lines[i].split() or ["?"])[0]
        if el in LIG_PLACEHOLDERS:
            grp = "ligand"
        elif el in core_elems:
            grp = "core"
        elif el in shell_elems:
            grp = "shell"
        else:
            grp = "ligand"
        counts[grp][el] += 1

    total = sum(sum(c.values()) for c in counts.values())
    fmt = lambda c: {"total": sum(c.values()), "by_element": dict(c)}  # noqa: E731
    return {
        "total_atoms": total,
        "core": fmt(counts["core"]),
        "shell": fmt(counts["shell"]),
        "ligand": fmt(counts["ligand"]),
    }


def _charge_of_xyz(xyz_text: str, charges: Dict[str, float]) -> float:
    lines = xyz_text.splitlines()
    try:
        n = int(lines[0].strip())
    except Exception:
        return 0.0

    total = 0.0
    for i in range(2, 2 + n):
        if i >= len(lines):
            break
        el = (lines[i].split() or ["?"])[0]
        total += charges.get(el, 0.0)
    return total


def run_cmd(cmd: List[str], cwd: Path) -> Tuple[str, str]:
    if cmd[0] == "nc-builder":
        exe = shutil.which("nc-builder")
        if not exe:
            raise RuntimeError("nc-builder not found on PATH for this server process.")
        cmd[0] = exe
    p = subprocess.run(cmd, cwd=str(cwd), check=True, text=True, capture_output=True)
    return p.stdout, p.stderr


def parse_cif_oxidation_numbers(text: str) -> Dict[str, float]:
    lines = text.splitlines()
    i = 0
    charges: Dict[str, float] = {}

    def read_loop(start: int):
        hdrs, rows, j = [], [], start + 1
        while j < len(lines) and lines[j].lstrip().startswith("_"):
            hdrs.append(lines[j].strip())
            j += 1
        while j < len(lines):
            t = lines[j].strip()
            if not t or t.startswith(("loop_", "data_", "_")):
                break
            rows.append(t)
            j += 1
        return hdrs, rows, j

    _KEY_ALIASES = {
        "_atom_type.symbol": "_atom_type_symbol",
        "_atom_type.oxidation_number": "_atom_type_oxidation_number",
        "_atom_type.charge": "_atom_type_charge",
        "_atom_site.type_symbol": "_atom_site_type_symbol",
        "_atom_site.oxidation_number": "_atom_site_oxidation_number",
        "_atom_site.charge": "_atom_site_charge",
        "_atom_site.label": "_atom_site_label",
    }

    def _norm_key(k: str) -> str:
        return _KEY_ALIASES.get(k.strip().lower(), k.strip().lower())

    def _elem_from_label(s: str) -> str:
        m = re.match(r"([A-Z][a-z]?)", s)
        return m.group(1) if m else s

    def _parse_charge_str(v: str) -> float:
        v = str(v).strip().replace(",", ".").replace(" ", "")
        if v.endswith("+"):
            v = v[:-1]
        elif v.endswith("-"):
            v = "-" + v[:-1]
        if v.startswith("+"):
            v = v[1:]
        return float(v) if v else 0.0

    while i < len(lines):
        if lines[i].strip().lower().startswith("loop_"):
            hdrs, rows, j = read_loop(i)
            headers = [_norm_key(h) for h in hdrs]
            sym_candidates = ["_atom_type_symbol", "_atom_site_type_symbol", "_atom_site_label"]
            ox_candidates = [
                "_atom_type_oxidation_number",
                "_atom_site_oxidation_number",
                "_atom_type_charge",
                "_atom_site_charge",
            ]
            si = next((headers.index(k) for k in sym_candidates if k in headers), None)
            oi = next((headers.index(k) for k in ox_candidates if k in headers), None)
            if si is not None and oi is not None:
                for row in rows:
                    toks = re.findall(r"(?:'[^']*'|\"[^\"]*\"|\S+)", row)
                    if len(toks) <= max(si, oi):
                        continue
                    raw_sym, raw_ox = toks[si].strip("'\""), toks[oi].strip("'\"")
                    if raw_sym in ".?" or raw_ox in ".?":
                        continue
                    sym = _elem_from_label(raw_sym)
                    try:
                        charges[sym] = _parse_charge_str(raw_ox)
                    except (ValueError, IndexError):
                        continue
            i = j
        else:
            i += 1

    return charges


def execute_builder_command(
    cmd: List[str], out_dir: str, primary_output_filename: str
) -> Tuple[bool, str, Dict[str, Path | None]]:
    out_dir_path = Path(out_dir)
    final_path = out_dir_path / primary_output_filename
    output_stem = final_path.stem
    cut_path = out_dir_path / f"{output_stem}_cut.xyz"

    logging.info(f"Running command in {out_dir_path}: {' '.join(cmd)}")

    try:
        p = subprocess.run(
            cmd,
            cwd=str(out_dir_path),
            capture_output=True,
            env=os.environ.copy(),
            check=False,
            text=True,
        )
        stdout = p.stdout or ""
        stderr = p.stderr or ""
        log_output = stdout + ("\n" + stderr if stderr else "")

        ok = p.returncode == 0 and final_path.exists() and final_path.stat().st_size > 0
        if not ok:
            logging.error(
                f"nc-builder execution check failed (return code: {p.returncode}).\nLogs:\n{log_output}"
            )

        found_files = {
            "final": final_path if final_path.exists() else None,
            "cut": cut_path if cut_path.exists() else None,
        }
        return ok, log_output, found_files

    except FileNotFoundError:
        msg = "Error: 'nc-builder' command not found. Is your conda environment activated?"
        logging.error(msg)
        return False, msg, {}

    except Exception as e:
        msg = f"An unexpected error occurred: {e}"
        logging.error(msg, exc_info=True)
        return False, msg, {}

def run_quiet(cmd, cwd=None, env=None):
    # Fastest: no stdout/err shipping over HTTP; write a log file for post-mortem.
    log_path = os.path.join(cwd or ".", "build.log")
    with open(log_path, "wb") as lf:
        p = subprocess.run(
            cmd,
            cwd=cwd,
            env=env,
            stdout=lf,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return p.returncode, log_path

def popen_stream_tty(cmd, cwd=None, env=None):
    """
    Spawn under a pseudo-terminal to force line-buffering in the child.
    Yields decoded lines ASAP. Slightly more CPU than pipes, but best interactivity.
    """
    master_fd, slave_fd = pty.openpty()
    # Optional: put the tty in raw-ish mode
    try:
        tty.setraw(master_fd)
    except Exception:
        pass

    # Make non-blocking
    fl = fcntl.fcntl(master_fd, fcntl.F_GETFL)
    fcntl.fcntl(master_fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)

    # Important: direct child's stdout/err to the slave tty
    p = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        stdin=slave_fd,
        stdout=slave_fd,
        stderr=slave_fd,
        bufsize=0,
        close_fds=True,
    )
    os.close(slave_fd)

    decoder = None
    try:
        while True:
            r, _, _ = select.select([master_fd], [], [], 0.25)
            if master_fd in r:
                try:
                    chunk = os.read(master_fd, 65536)
                except BlockingIOError:
                    chunk = b""
                if not chunk:
                    # could be just EAGAIN; fall through to poll
                    pass
                else:
                    if decoder is None:
                        decoder = (lambda b: b.decode("utf-8", "replace"))
                    text = decoder(chunk)
                    for line in text.splitlines():
                        yield line
            if p.poll() is not None:
                # Drain any residual
                try:
                    while True:
                        chunk = os.read(master_fd, 65536)
                        if not chunk:
                            break
                        text = (decoder or (lambda b: b.decode("utf-8","replace")))(chunk)
                        for line in text.splitlines():
                            yield line
                except Exception:
                    pass
                break
    finally:
        try: os.close(master_fd)
        except Exception: pass

def popen_stream(cmd: List[str], cwd: Path, env: Optional[Dict[str, str]] = None):
    """
    Yield text lines (stdout+stderr merged) from a long-running command.
    """
    if cmd[0] == "nc-builder":
        nc_builder_path = shutil.which("nc-builder")
        if not nc_builder_path:
            raise RuntimeError("nc-builder not found on PATH for this server process.")
        python_exe = sys.executable
        cmd = [python_exe, "-u", nc_builder_path] + cmd[1:]

    base_env = os.environ.copy()
    if env:
        base_env.update(env)

    # Ensure unbuffered UTF-8 streams for child
    base_env["PYTHONUNBUFFERED"] = "1"
    base_env["PYTHONIOENCODING"] = "utf-8"
    base_env["NCBUILDER_FLUSH"] = "1"
    base_env["QD_BUILDER_UNBUFFERED"] = "1"
    # ---------------------
    p = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=base_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        text=True,                # <- text mode (aka universal_newlines=True)
        encoding="utf-8",
        errors="replace",         # avoid crashes on odd bytes
    )
    try:
        for line in p.stdout:
            # line already decoded text with trailing newline if any
            yield line
    finally:
        if p.stdout:
            p.stdout.close()
        p.wait()

def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


def format_facets_to_dict(facets: List[Facet]) -> Dict[str, float]:
    return {f.hkl: f.gamma for f in facets}


def _dist_string(jobs: List[LigJob], mode: str) -> str:
    xs = [f"{float(max(0.0, j.ratio)):.6g}" for j in jobs if j.smiles.strip()]
    if not xs:
        xs = ["1"]
    mode = (mode or "uniform").lower()
    if mode not in {"uniform", "segmented", "random"}:
        mode = "uniform"
    return ":".join(xs + [mode])

def run_capper_cli(script_path: Path, xyz_text: str, jobs: List[LigJob],
                   dummy: str, dist_mode: str, workdir: Path) -> Tuple[str, str, str]:
    """Return (capped_xyz, cap_log, out_name)."""
    qd_path = workdir / "qd_for_capping.xyz"
    qd_path.write_text(xyz_text, encoding="utf-8")

    out_prefix = workdir / f"capped_{dummy.lower()}"
    smiles = [j.smiles.strip() for j in jobs if j.smiles.strip()]
    if not smiles:
        return xyz_text, "[cap] No ligands provided; skipping.", ""

    # Use system-wide miniCAT CLI (on PATH) instead of the bundled script
    exe = shutil.which("miniCAT")
    if not exe:
        raise RuntimeError("miniCAT executable not found on PATH. Activate the env where miniCAT is installed.")

    # We assume miniCAT exposes a 'passivate' subcommand compatible with attach_from_smiles
    cmd = [
        exe, 
        "--qd", str(qd_path),
        "--out_prefix", str(out_prefix),
        "--job-ligands", *smiles,
        "--job-dummy", dummy,
        "--job-dist", _dist_string(jobs, dist_mode),
    ]
    p = subprocess.run(cmd, cwd=str(workdir), text=True, capture_output=True)
    log = (p.stdout or "") + ("\n" + p.stderr if p.stderr else "")

    candidates = [f"{out_prefix}_final.xyz", f"{out_prefix}.xyz", f"{out_prefix}_capped.xyz"]
    xyz_out, out_name = None, ""
    for cand in candidates:
        cp = Path(cand)
        if cp.exists() and cp.stat().st_size > 0:
            xyz_out = cp.read_text(encoding="utf-8")
            out_name = cp.name
            break

    if p.returncode != 0 or not xyz_out:
        raise RuntimeError(f"Capper failed (rc={p.returncode}).\n{log}")
    return xyz_out, log, out_name

def run_in_memory_capper(xyz_text: str, jobs: List[LigJob], default_dummy: str, dist_mode: str, out_stream=sys.stdout) -> str:
    """Runs miniCAT natively in memory and streams its logs to out_stream."""
    if not jobs:
        print("[cap] No ligands provided.", file=out_stream)
        return xyz_text
    
    valid_jobs = [j for j in jobs if j.smiles.strip()]
    if not valid_jobs:
        print("[cap] No valid SMILES.", file=out_stream)
        return xyz_text

    args = argparse.Namespace(
        seed=1337, ff="uff", refinement_passes=3, neigh=8,
        coarse_step_deg=20.0, sterics_mode="vdw", sterics_margin=0.25,
        warn_tol=1.4, adaptive_offset_steps=4, adaptive_offset_step=0.15,
        anchor_mode="dummy", offset_out=0.0, bond_len_override=None,
        neighbor_repulsion=0.5
    )
    reg = get_fg_registry()

    smiles_list = [j.smiles.strip() for j in valid_jobs]
    ratios = [j.ratio for j in valid_jobs]
    dummy = valid_jobs[0].dummy.strip() if valid_jobs[0].dummy else default_dummy
    
    # Redirect prints directly to our stream queue
    with contextlib.redirect_stdout(out_stream), contextlib.redirect_stderr(out_stream):
        print(f"[cap] Passivating {dummy} sites natively with {smiles_list} (mode={dist_mode})")
        precomputed_ligands = []
        role = None
        for sm in smiles_list:
            mol_neutral = smiles_to_3d_mol(_sanitize_simple_acid_smiles(sm), args.seed, args.ff)
            fg_matches = detect_fg_matches_neutral(mol_neutral, 'anion', reg) or detect_fg_matches_neutral(mol_neutral, 'cation', reg)
            if not fg_matches:
                print(f"No recognizable functional group for {sm}")
                continue
            ionic_mol, (fg_active, match_active) = build_per_site_variant(fg_matches, mol_neutral, 0, args.ff)
            precomputed_ligands.append({"mol": ionic_mol, "fg": fg_active, "match": match_active})
            if role is None:
                role = fg_matches[0][0].role

        job_dict = {
            "smiles": smiles_list, "dummy": dummy, "role": role,
            "ratios": ratios, "strategy": dist_mode, "precomputed_ligands": precomputed_ligands
        }
        
        current_qd = ase_read(io.StringIO(xyz_text), format="xyz")
        current_qd = execute_passivation_job(current_qd, job_dict, args)
        
        out_xyz = io.StringIO()
        ase_write(out_xyz, current_qd, format="xyz")
        
    return out_xyz.getvalue()


# ---------------------------------------------------------------------
# Route: build
# ---------------------------------------------------------------------
@app.post("/api/build", response_class=JSONResponse)
async def build_nanocrystal(files: List[UploadFile] = File(...), options: str = Form(...)):
    try:
        opts = BuildOptions.parse_raw(options)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Invalid options JSON: {e}")

    tmpdir = tempfile.mkdtemp(prefix="qdb_")
    logging.info(f"Temporary directory created at: {tmpdir}")

    try:
        tmp_path = Path(tmpdir)

        # Save uploads
        file_map: Dict[str, Path] = {}
        for uploaded_file in files:
            s_filename = safe_filename(uploaded_file.filename)
            file_path = tmp_path / s_filename
            file_path.write_bytes(await uploaded_file.read())
            file_map[s_filename] = file_path

        # Core
        core_path = file_map.get(safe_filename(opts.core_cif_filename))
        if not core_path:
            raise HTTPException(404, f"Core file '{opts.core_cif_filename}' not found.")

        charges = parse_cif_oxidation_numbers(core_path.read_text("utf-8", "ignore"))
        core_elements = set(charges.keys())
        shell_elements: set[str] = set()

        import yaml

        # ---------------- First pass (info) ----------------
        logging.info("Running first pass to gather info...")
        first_pass_charges = charges.copy()
        if "Cl" not in first_pass_charges:
            first_pass_charges["Cl"] = -1.0

        temp_yaml_dict = {
            "shape": {"aspect": opts.aspect},
            "facets": [f.dict() for f in opts.facets],
            "charges": first_pass_charges,
            "passivation": {"ligand": "Cl", "surf_tol": 2.0},
        }
        temp_yaml_path = tmp_path / "temp_config.yml"
        temp_yaml_path.write_text(yaml.safe_dump(temp_yaml_dict, sort_keys=False), encoding="utf-8")

        first_pass_output_file = "final.xyz"
        first_pass_cmd = [
            "nc-builder",
            str(core_path.resolve()),
            str(temp_yaml_path.resolve()),
            "-r",
            f"{opts.radius_A:.4f}",
            "-o",
            str((tmp_path / first_pass_output_file).resolve()),
            "--center",
            "--write-all",
            "--verbose",
        ]
        p1_ok, p1_log, _ = execute_builder_command(first_pass_cmd, tmpdir, first_pass_output_file)
        if not p1_ok:
            return JSONResponse(status_code=500, content={"status": "failed", "log": p1_log})

        # ---------------- Second pass (final build) ----------------
        final_charges = charges.copy()
        if "Cl" not in final_charges:
            final_charges["Cl"] = -1.0
        if "Rb" not in final_charges:
            final_charges["Rb"] = 1.0

        passivation_block = {"ligand": "Cl", "surf_tol": 2.0}

        if not opts.shells:
            yaml_dict = {
                "shape": {"aspect": opts.aspect},
                "facets": [f.dict() for f in opts.facets],
                "charges": final_charges,
                "passivation": passivation_block,
            }
            cmd = ["nc-builder", str(core_path.resolve())]
        else:
            materials = [
                {
                    "name": "core",
                    "cif": str(core_path.resolve()),
                    "facets": format_facets_to_dict(opts.facets),
                    "shape": {"aspect": opts.aspect},
                }
            ]
            outermost_shell_path = core_path

            for i, shell in enumerate(opts.shells):
                shell_cif_path = file_map.get(safe_filename(shell.material_cif))
                if not shell_cif_path:
                    raise HTTPException(404, f"Shell material '{shell.material_cif}' not found.")

                materials.append(
                    {
                        "name": f"shell{i+1}",
                        "cif": str(shell_cif_path.resolve()),
                        "facets": format_facets_to_dict(shell.facets) or "inherit",
                        "shape": {"aspect": shell.aspect},
                    }
                )

                shell_charges = parse_cif_oxidation_numbers(shell_cif_path.read_text("utf-8", "ignore"))
                final_charges.update(shell_charges)
                shell_elements.update(shell_charges.keys())
                outermost_shell_path = shell_cif_path

            yaml_dict = {
                "materials": materials,
                "charges": final_charges,
                "symmetry": {"proper_rotations_only": True},
                "facet_options": {"pair_opposites": True},
                "passivation": passivation_block,
            }
            cmd = ["nc-builder", str(outermost_shell_path.resolve())]

        final_yaml_path = tmp_path / "config.yml"
        final_yaml_path.write_text(yaml.safe_dump(yaml_dict, sort_keys=False), encoding="utf-8")
        logging.info(f"Generated FINAL YAML config:\n{final_yaml_path.read_text()}")

        final_output_file = "final.xyz"
        cmd.extend(
            [
                str(final_yaml_path.resolve()),
                "-r",
                f"{opts.radius_A:.4f}",
                "-o",
                str((tmp_path / final_output_file).resolve()),
                "--write-all",
                "--center",
                "--verbose",
            ]
        )
        if opts.shells:
            cmd.extend(["--core-lattice-fit", "--core-strain-width", "2.5", "--core-center", "com"])

        ok, log, out_files = execute_builder_command(cmd, tmpdir, final_output_file)

        xyz_passivated_content = out_files.get("final").read_text() if out_files.get("final") else None
        xyz_unpassivated_content = out_files.get("cut").read_text() if out_files.get("cut") else None

        full_log = p1_log + "\n--- Second Pass ---\n" + log
        raw_log = full_log

        if not ok:
            return JSONResponse(status_code=500, content={"status": "failed", "log": raw_log})

        elements, total_charge = "unknown", 0
        grouped_counts: Dict[str, Any] = {}



        if xyz_passivated_content:
            try:
                from ase.io import read as ase_read

                atoms = ase_read(io.StringIO(xyz_passivated_content), format="xyz")
                symbols = atoms.get_chemical_symbols()
                elements = ",".join(sorted(set(symbols)))
                total_charge = sum(final_charges.get(s, 0.0) for s in symbols)
                grouped_counts = _count_grouped_xyz(xyz_passivated_content, core_elements, shell_elements)
            except Exception as e:
                logging.error(f"Could not parse passivated XYZ or calculate charge: {e}")

        # ---------- Organic passivation from SMILES (optional) ----------
        CAP_SCRIPT = None 
 
        cap_logs: List[str] = []
        current_xyz = xyz_passivated_content or xyz_unpassivated_content or ""
        download_name = None 
        try:
            if opts.cap_anionic_jobs:
                current_xyz, log1, name1 = run_capper_cli(
                    CAP_SCRIPT, current_xyz, opts.cap_anionic_jobs or [],
                    dummy="Cl", dist_mode=(opts.cap_distribution or "uniform"),
                    workdir=tmp_path
                )
                cap_logs.append(log1)
                download_name = name1 or download_name
        
            if opts.cap_cationic_jobs:
                current_xyz, log2, name2 = run_capper_cli(
                    CAP_SCRIPT, current_xyz, opts.cap_cationic_jobs or [],
                    dummy="Rb", dist_mode=(opts.cap_distribution or "uniform"),
                    workdir=tmp_path
                )
                cap_logs.append(log2)
                download_name = name2 or download_name
        
            if cap_logs:
                xyz_passivated_content = current_xyz
                # recompute metadata on capped structure
                from ase.io import read as ase_read

                atoms = ase_read(io.StringIO(current_xyz), format="xyz")
                symbols = atoms.get_chemical_symbols()
                elements = ",".join(sorted(set(symbols)))
                total_charge = sum(final_charges.get(s, 0.0) for s in symbols)
                grouped_counts = _count_grouped_xyz(current_xyz, core_elements, shell_elements)
                raw_log += "\n--- Organic capping ---\n" + "\n".join(cap_logs)

        except Exception as e:
            logging.exception("Capping step failed")
            raw_log += f"\n[cap][error] {e}"
        # -----------------------------------------------------------------

        # ---- compact summary ----
        init_charge = 0.0
        if xyz_unpassivated_content:
            init_charge = _charge_of_xyz(xyz_unpassivated_content, final_charges)

        cl_n = grouped_counts.get("ligand", {}).get("by_element", {}).get("Cl", 0)
        rb_n = grouped_counts.get("ligand", {}).get("by_element", {}).get("Rb", 0)
        core_comp = grouped_counts.get("core", {}).get("by_element", {})
        shell_comp = grouped_counts.get("shell", {}).get("by_element", {})

        summary_lines = [
            "Summary:",
            f"- Atoms total: {grouped_counts.get('total_atoms', 'NA')}",
            f"- Core:  {grouped_counts.get('core', {}).get('total', 0)}  | "
            + " ".join(f"{k}:{v}" for k, v in sorted(core_comp.items())),
            f"- Shell: {grouped_counts.get('shell', {}).get('total', 0)} | "
            + " ".join(f"{k}:{v}" for k, v in sorted(shell_comp.items())),
            f"- Ligand placeholders: Cl={cl_n}, Rb={rb_n}",
            f"- Charge: initial≈{round(init_charge)}, final≈{round(total_charge)}",
        ]
        raw_log = (raw_log + "\n" + "\n".join(summary_lines)).strip()

        return JSONResponse(
            content={
                "status": "success",
                "log": raw_log,
                "elements": elements,
                "xyz_passivated": xyz_passivated_content,
                "xyz_unpassivated": xyz_unpassivated_content,
                "total_charge": round(total_charge),
                "grouped_counts": grouped_counts,
                "download_name": download_name or "final.xyz"
            }
        )

    finally:
        logging.info(f"Cleaning up temporary directory: {tmpdir}")
        # shutil.rmtree(tmpdir)

class CapOnlyOptions(BaseModel):
    cap_distribution: Optional[str] = "uniform"
    cap_anionic_jobs: Optional[List[LigJob]] = None
    cap_cationic_jobs: Optional[List[LigJob]] = None

@app.post("/api/passivate")
async def passivate_endpoint(
    xyz_file: UploadFile = File(...),
    options: str = Form(...),
):
    """
    JSON in, JSON out. Response shape:
      {
        "status": "success" | "failed",
        "cmd": "python -u miniCAT --input ... --distribution ... --anion ... --ratio ...",
        "xyz_passivated": "<XYZ text>",                # on success
        "download_name": "capped_final.xyz",           # on success
        "log": "...",                                  # optional text log
        "error": "message", "traceback": "..."         # on failure
      }
    """

    # --- helper: normalize incoming jobs (dicts/objects) to attribute-style objects ---
    def _normalize_jobs(jobs):
        """
        Accepts a list of dicts or objects with (smiles, ratio) and
        returns a list of SimpleNamespace(smiles=..., ratio=float).
        Drops entries with empty smiles.
        """
        norm = []
        for j in jobs or []:
            if isinstance(j, dict):
                s = (j.get("smiles") or "").strip()
                r = float(j.get("ratio") or 0)
            else:
                s = (getattr(j, "smiles", "") or "").strip()
                r = float(getattr(j, "ratio", 0) or 0)
            if s:
                norm.append(SimpleNamespace(smiles=s, ratio=r))
        return norm

    try:
        cfg = json.loads(options or "{}")
    except Exception as e:
        return JSONResponse({"status": "failed", "error": f"Invalid options JSON: {e}"}, status_code=400)

    tmpdir = tempfile.mkdtemp(prefix="qdb_pass_")
    tmp = Path(tmpdir)
    pseudo_cmd = None  # filled before running so we can echo it even on failure

    try:
        # --- write incoming XYZ to disk ---
        xyz_path = tmp / (xyz_file.filename or "input.xyz")
        xyz_bytes = await xyz_file.read()
        xyz_path.write_bytes(xyz_bytes)

        # --- read + normalize options ---
        raw_anionic  = cfg.get("cap_anionic_jobs")  or []
        raw_cationic = cfg.get("cap_cationic_jobs") or []
        dist         = (cfg.get("cap_distribution") or "uniform").lower()
        if dist not in {"uniform", "segmented", "random"}:
            dist = "uniform"

        anionic  = _normalize_jobs(raw_anionic)
        cationic = _normalize_jobs(raw_cationic)

        # --- build a human-readable command we will echo back ---
        CAP_SCRIPT = None  

        exe_shown = shutil.which("miniCAT") or "miniCAT"
        
        # show two sequential miniCAT calls (anionic then cationic), matching how run_capper_cli is invoked
        def _dist_str(jobs):  # reuse same allocation logic for the display string
            return _dist_string(jobs, dist)
        
        an_cmd = f'{exe_shown} passivate --qd "{xyz_path}" --out_prefix "{(tmp / "capped_cl")}" --job-ligands ' + \
                 " ".join(j.smiles for j in anionic) + f' --job-dummy Cl --job-dist "{_dist_str(anionic)}"' if anionic else ""
        ca_cmd = f'{exe_shown} passivate --qd "{xyz_path}" --out_prefix "{(tmp / "capped_rb")}" --job-ligands ' + \
                 " ".join(j.smiles for j in cationic) + f' --job-dummy Rb --job-dist "{_dist_str(cationic)}"' if cationic else ""
        
        pseudo_cmd = " && ".join([c for c in (an_cmd, ca_cmd) if c])
        if not pseudo_cmd:
            pseudo_cmd = f"{exe_shown} passivate (no ligands provided)"

        # --- run passivation using your helper (no external shell needed) ---
        current_xyz = xyz_path.read_text()
        log_chunks = []

        def run_cap(jobs, ion_label):
            nonlocal current_xyz
            if not jobs:
                return None
            # run_capper_cli must accept (script_path, xyz_text, jobs, ion_symbol, dist, workdir)
            xyz_out, cap_log, out_name = run_capper_cli(
                CAP_SCRIPT, current_xyz, jobs, ion_label, dist, tmp
            )
            current_xyz = xyz_out
            if cap_log:
                log_chunks.append(cap_log)
            return out_name

        name1 = run_cap(anionic,  "Cl")
        name2 = run_cap(cationic, "Rb")
        dl_name = name2 or name1 or "capped_final.xyz"

        # --- corrected total charge (ligands: anion = -1, cation = +1) ---
        before_counts = _xyz_count_symbols(xyz_path.read_text(), targets=("Cl","Rb"))
        after_counts  = _xyz_count_symbols(current_xyz,       targets=("Cl","Rb"))

        n_anionic  = max(0, before_counts.get("Cl", 0) - after_counts.get("Cl", 0))
        n_cationic = max(0, before_counts.get("Rb", 0) - after_counts.get("Rb", 0))

        elem_charge = {"Cl": -1.0, "Rb": +1.0}
        raw_total = 0.0
        for line in current_xyz.splitlines()[2:]:
            parts = line.split()
            if parts:
                raw_total += elem_charge.get(parts[0], 0.0)

        total_charge = int(round(raw_total + (-1)*n_anionic + (+1)*n_cationic))

        ligand_detail = {
            "anionic":  _allocate_counts_by_ratio(n_anionic,  anionic),
            "cationic": _allocate_counts_by_ratio(n_cationic, cationic),
            "total": int(n_anionic + n_cationic),
        }
        
        return JSONResponse({
            "status": "success",
            "cmd": pseudo_cmd,
            "xyz_passivated": current_xyz,
            "download_name": dl_name,
            "log": "\n".join(log_chunks)[:200000],
            "total_charge": total_charge,   # ← include it
            "ligand_detail": ligand_detail, 
        })

    except Exception as e:
        tb = traceback.format_exc(limit=8)
        return JSONResponse({
            "status": "failed",
            "error": str(e),
            "traceback": tb,
            "cmd": pseudo_cmd or "N/A",
        }, status_code=500)
    finally:
        # Keep tmpdir for debugging? Comment next line to keep artifacts.
        # shutil.rmtree(tmpdir, ignore_errors=True)
        pass


# #####################################################################
# ###                  CORRECTED STREAMING ROUTE                    ###
# #####################################################################
#
# The placeholder function has been removed and the @app.post decorator
# is now on the real implementation below.
@app.post("/api/build_stream")

async def build_nanocrystal_stream(
    files: List[UploadFile] = File(...),
    options: str = Form(...),
    mode: str = Form("quiet"), # "quiet" (fast, default), "live-tty", or "live-pipe"
    positive_q_mode: str = Form("remove")
):
    """
    NDJSON stream:
      {"event":"status","line":"accepted"}           # flushed immediately
      {"event":"diag","files":[...],"shells":[...]}  # what server sees
      {"event":"cmd","line":"nc-builder ..."}        # the single run
      {"event":"log","line":"..."}                   # only in live modes
      {"event":"result", ...payload...}
    """

    # ---- safe facet formatter (fallback if helper not present) ----
    def format_facets_to_dict_safe(facets):
        """
        Accepts a list of objects/dicts like {"hkl": "100", "gamma": 1.0}
        and returns [{"hkl":"100","gamma":1.0}, ...] with normalized types.
        """
        out = []
        for f in facets or []:
            if hasattr(f, "hkl"):
                hkl = getattr(f, "hkl")
                gamma = float(getattr(f, "gamma"))
            else:
                hkl = f.get("hkl")
                gamma = float(f.get("gamma", 0.0))
            out.append({"hkl": str(hkl), "gamma": gamma})
        return out

    try:
        opts = BuildOptions.parse_raw(options)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Invalid options JSON: {e}")

    tmpdir = tempfile.mkdtemp(prefix="qdb_stream_")
    tmp_path = Path(tmpdir)

    # ---- READ & BUFFER UPLOADS BEFORE streaming (prevents closed-body errors) ----
    buffered_uploads: List[Tuple[str, bytes]] = []
    for uf in files:
        fname = safe_filename(uf.filename)
        data = await uf.read()
        buffered_uploads.append((fname, data))

    async def _run():
        # Always flush headers immediately so frontend won't see "Failed to fetch"
        yield json.dumps({"event": "status", "line": "accepted"}) + "\n"

        try:
            # ---- WRITE buffered uploads to disk inside the generator ----
            file_map: Dict[str, Path] = {}
            for fname, data in buffered_uploads:
                fp = tmp_path / fname
                fp.write_bytes(data)
                file_map[fname] = fp

            # Small diagnostic to help catch filename mismatches for shells
            try:
                diag = {
                    "files": sorted(list(file_map.keys())),
                    "core_expected": getattr(opts, "core_cif_filename", None),
                    "shells_expected": [getattr(sh, "material_cif", None) for sh in (opts.shells or [])],
                }
                yield json.dumps({"event": "diag", **diag}) + "\n"
            except Exception:
                pass

            core_basename = safe_filename(opts.core_cif_filename)
            core_path = file_map.get(core_basename)
            if not core_path:
                yield json.dumps({"event":"log","line": f"[error] Core file '{opts.core_cif_filename}' not found among uploaded: {', '.join(file_map)}"}) + "\n"
                yield json.dumps({"event":"result","status":"failed"}) + "\n"
                return

            # ---- Charges and bookkeeping ----
            charges = parse_cif_oxidation_numbers(core_path.read_text("utf-8", "ignore"))
            core_elements = set(charges.keys())
            shell_elements: set[str] = set()

            # ---- Build ONE YAML (core-only OR core@shell) ----
            final_charges = charges.copy()
            final_charges.setdefault("Cl", -1.0)
            final_charges.setdefault("Rb",  1.0)

            if not opts.shells:
                # Core-only schema
                facets_input = [f.dict() if hasattr(f, "dict") else f for f in (opts.facets or [])]
                if not facets_input:
                    # Default low-index cubic facets if user didn’t provide any
                    facets_input = [
                        {"hkl": "100", "gamma": 1.0},
                        {"hkl": "110", "gamma": 1.0},
                        {"hkl": "111", "gamma": 1.0},
                    ]

                yaml_dict = {
                    "shape":   {"aspect": opts.aspect},
                    "facets":  format_facets_to_dict_safe(facets_input),
                    "charges": final_charges,
                    "passivation": {"ligand": "Cl", "surf_tol": 2.0},
                }
                root_path = core_path

            else:
                # Core@shell schema
                materials = [{
                    "name": "core",
                    "cif": str(core_path.resolve()),
                    "facets": format_facets_to_dict_safe([f.dict() if hasattr(f, "dict") else f for f in opts.facets]),
                    "shape": {"aspect": opts.aspect},
                }]
                outer = core_path
                for i, sh in enumerate(opts.shells):
                    shell_name = getattr(sh, "material_cif", "")
                    sp = file_map.get(safe_filename(shell_name))
                    if not sp:
                        yield json.dumps({"event":"log","line": f"[error] Shell '{shell_name}' not found among uploaded: {', '.join(file_map)}"}) + "\n"
                        yield json.dumps({"event":"result","status":"failed"}) + "\n"
                        return
                    materials.append({
                        "name":  f"shell{i+1}",
                        "cif":   str(sp.resolve()),
                        "facets": format_facets_to_dict_safe([f.dict() if hasattr(f, "dict") else f for f in (sh.facets or [])]) or "inherit",
                        "shape": {"aspect": getattr(sh, "aspect", None) or [1.0, 1.0, 1.0]},
                    })
                    sc = parse_cif_oxidation_numbers(sp.read_text("utf-8", "ignore"))
                    final_charges.update(sc)
                    shell_elements.update(sc.keys())
                    outer = sp

                yaml_dict = {
                    "materials": materials,
                    "charges": final_charges,
                    "symmetry": {"proper_rotations_only": True},
                    "facet_options": {"pair_opposites": True},
                    "passivation": {"ligand": "Cl", "surf_tol": 2.0},
                }
                root_path = outer

            final_yaml = tmp_path / "config.yml"
            final_yaml.write_text(yaml.safe_dump(yaml_dict, sort_keys=False), encoding="utf-8")

            # ---- Single nc-builder command ----
            final_out = tmp_path / "final.xyz"
            posq = (positive_q_mode or "remove").strip().lower()
            if posq not in ("remove", "add"):
                posq = "remove"

            cmd = [
                "nc-builder",
                str(root_path.resolve()),
                str(final_yaml.resolve()),
                "-r", f"{opts.radius_A:.4f}",
                "-o", str(final_out.resolve()),
                "--write-all", "--center", "--verbose",
                "--positive-q-mode", posq,
            ]
            if opts.shells:
                cmd += ["--core-lattice-fit", "--core-strain-width", "2.5", "--core-center", "com"]

            yield json.dumps({"event": "cmd", "line": " ".join(cmd)}) + "\n"

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

            # ==========================================================
            # SMART BUILD EXECUTION (Skip if core unchanged)
            # ==========================================================
            if getattr(opts, "skip_core_build", False) and getattr(opts, "xyz_unpassivated", None):
                yield json.dumps({"event": "log", "line": "[system] Core parameters unchanged. Bypassing nc-builder and re-using existing core..."}) + "\n"
                xyz_un = opts.xyz_unpassivated
                xyz_pass = xyz_un
                cmd = ["(skipped nc-builder)"]
            else:
                # Run the heavy nc-builder
                final_out = tmp_path / "final.xyz"
                posq = (positive_q_mode or "remove").strip().lower()
                if posq not in ("remove", "add"): posq = "remove"

                cmd = [
                    "nc-builder", str(root_path.resolve()), str(final_yaml.resolve()),
                    "-r", f"{opts.radius_A:.4f}", "-o", str(final_out.resolve()),
                    "--write-all", "--center", "--verbose", "--positive-q-mode", posq,
                ]
                if opts.shells:
                    cmd += ["--core-lattice-fit", "--core-strain-width", "2.5", "--core-center", "com"]

                yield json.dumps({"event": "cmd", "line": " ".join(cmd)}) + "\n"

                log_queue = queue.Queue()

                def _run_builder_sync():
                    qw = LineQueueWriter(log_queue)
                    argv = cmd[1:] 
                    old_cwd = os.getcwd()  # Save current dir
                    with contextlib.redirect_stdout(qw), contextlib.redirect_stderr(qw):
                        try: 
                            os.chdir(tmpdir)  # FORCE nc-builder outputs into the tmp folder
                            nc_builder_main(argv)
                        except SystemExit as e:
                            if e.code != 0 and e.code is not None: print(f"[error] nc-builder exited with code {e.code}")
                        except Exception as e: print(f"[fatal] {traceback.format_exc()}")
                        finally: 
                            os.chdir(old_cwd) # Restore dir safely
                            qw.flush()
                            log_queue.put(None)

                builder_thread = threading.Thread(target=_run_builder_sync)
                builder_thread.start()

                while True:
                    try:
                        line = log_queue.get_nowait()
                        if line is None: break
                        if line.strip(): yield json.dumps({"event": "log", "line": line}) + "\n"
                    except queue.Empty:
                        await asyncio.sleep(0.01)
                
                builder_thread.join()

                if not final_out.exists() or final_out.stat().st_size == 0:
                    yield json.dumps({"event": "result", "status": "failed", "log": "nc-builder did not produce final.xyz"}) + "\n"
                    return

                xyz_pass = final_out.read_text()
                xyz_cut = (tmp_path / f"{final_out.stem}_cut.xyz")
                xyz_un = xyz_cut.read_text() if xyz_cut.exists() else xyz_pass

            elements, total_charge, grouped_counts = "unknown", 0, {}
            try:
                from ase.io import read as ase_read
                atoms = ase_read(io.StringIO(xyz_pass), format="xyz")
                symbols = atoms.get_chemical_symbols()
                elements = ",".join(sorted(set(symbols)))
                total_charge = sum(final_charges.get(s, 0.0) for s in symbols)
                grouped_counts = _count_grouped_xyz(xyz_pass, core_elements, shell_elements)
            except Exception as e:
                logging.error(f"XYZ parse fail: {e}")

            # ---- Optional SMILES capping (Native Memory Integration) ----
            # ---- Optional SMILES capping (Native Memory Integration) ----
            current_xyz = xyz_pass or xyz_un or ""
            download_name = None

            def _dist():
                m = (opts.cap_distribution or "uniform").lower()
                return m if m in {"uniform", "segmented", "random"} else "uniform"

            # Async helper to thread the capper and stream logs to frontend
            async def stream_capper(jobs, dummy_label):
                nonlocal current_xyz
                log_q = queue.Queue()
                qw = LineQueueWriter(log_q)
                res_box = []

                def _thread_cap():
                    try:
                        res = run_in_memory_capper(current_xyz, jobs, dummy_label, _dist(), out_stream=qw)
                        res_box.append(res)
                    except Exception as e:
                        print(f"[cap][error] {e}", file=qw)
                    finally:
                        qw.flush()
                        log_q.put(None)

                t = threading.Thread(target=_thread_cap)
                t.start()
                
                while True:
                    try:
                        line = log_q.get_nowait()
                        if line is None: break
                        if line.strip(): yield json.dumps({"event": "log", "line": line.strip()}) + "\n"
                    except queue.Empty:
                        await asyncio.sleep(0.01)
                t.join()
                
                if res_box:
                    current_xyz = res_box[0]

            if opts.cap_anionic_jobs:
                async for chunk in stream_capper(opts.cap_anionic_jobs, "Cl"):
                    yield chunk

            if opts.cap_cationic_jobs:
                async for chunk in stream_capper(opts.cap_cationic_jobs, "Rb"):
                    yield chunk

            # Recompute quick metadata after capping (best-effort)
            try:
                if current_xyz:
                    from ase.io import read as ase_read
                    atoms = ase_read(io.StringIO(current_xyz), format="xyz")
                    symbols = atoms.get_chemical_symbols()
                    elements = ",".join(sorted(set(symbols)))
            
                    # Raw charge from element map (Cl = -1, Rb = +1, organics 0…)
                    raw_total = sum(final_charges.get(s, 0.0) for s in symbols)
            
                    # Infer how many placeholders were replaced by ligands:
                    # compare counts before passivation (xyz_pass) vs after (current_xyz).
                    before_counts = _xyz_count_symbols(xyz_pass, targets=("Cl","Rb"))
                    after_counts  = _xyz_count_symbols(current_xyz, targets=("Cl","Rb"))
                    n_anionic  = max(0, before_counts.get("Cl", 0) - after_counts.get("Cl", 0))  # ligands placed on anionic sites
                    n_cationic = max(0, before_counts.get("Rb", 0) - after_counts.get("Rb", 0))  # ligands placed on cationic sites
            
                    # Each anionic ligand is −1, each cationic ligand is +1
                    ligand_correction = (-1) * n_anionic + (+1) * n_cationic
                    total_charge = int(round(raw_total + ligand_correction))

                    # Per-SMILES counts estimated from ratios
                    ligand_detail = {
                        "anionic":  _allocate_counts_by_ratio(n_anionic,  opts.cap_anionic_jobs),
                        "cationic": _allocate_counts_by_ratio(n_cationic, opts.cap_cationic_jobs),
                        "total": int(n_anionic + n_cationic),
                    }
                                
                    grouped_counts = _count_grouped_xyz(current_xyz, core_elements, shell_elements)
            except Exception as e:
                logging.error(f"Post-cap metadata recompute failed: {e}")
            
            size_metrics = None
            if current_xyz:
                try:
                    from ase.io import read as ase_read
                    # Read the XYZ string into an ASE atoms object
                    tmp_atoms = ase_read(io.StringIO(current_xyz), format="xyz")
                    c_coords = tmp_atoms.get_positions()
                    c_symbols = tmp_atoms.get_chemical_symbols()
                    
                    # Calculate metrics
                    size_metrics = get_cluster_size_metrics(c_coords, c_symbols)
                except Exception as e:
                    logging.error(f"Failed to calculate size metrics: {e}")           

            payload = {
                "status": "success",
                "elements": elements,
                "xyz_passivated": current_xyz,
                "xyz_unpassivated": xyz_un,
                "xyz_dummy": xyz_pass,
                "total_charge": round(total_charge),
                "grouped_counts": grouped_counts,
                "download_name": download_name or "final.xyz",
                "last_command": " ".join(cmd),
                "ligand_detail": ligand_detail,
                "size_metrics": size_metrics, 
            }
            yield json.dumps({"event": "result", **payload}) + "\n"

        except Exception as e:
            # Stream a readable error instead of letting the connection drop
            tb = traceback.format_exc(limit=5)
            yield json.dumps({"event": "log", "line": f"[fatal] {e}\n{tb}"}) + "\n"
            yield json.dumps({"event": "result", "status": "failed"}) + "\n"
        finally:
            logging.info(f"Cleaning up temporary directory: {tmpdir}")
            # shutil.rmtree(tmpdir, ignore_errors=True)
    
    return StreamingResponse(
        _run(),
        # OLD: media_type="application/x-ndjson",
        media_type="text/plain; charset=utf-8",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"}
    )

# ---------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------

if __name__ == "__main__":
    import uvicorn

    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=False)

