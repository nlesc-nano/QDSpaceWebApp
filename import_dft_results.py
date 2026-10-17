#!/usr/bin/env python3
"""
Publish one dot's DFT chain output (CP2K PBE relaxation + QDEX) into its
library record folder.

    <chain_dir>/relaxed.xyz                  PBE-relaxed geometry
    <chain_dir>/qdex/webapp.json             QDEX summary, band-edge excitons, spectra
    <chain_dir>/qdex/spatial_MO_*.cube       HOMO-1, HOMO, LUMO, LUMO+1
    <chain_dir>/qdex/*.html                  QDEX dashboards
    <chain_dir>/qdex_trap/                   optional: QDEX with the trap filter on (transitions with surface
                                             trap states excluded); its exciton pages and summary are
                                             published as qdex/*_traps_filtered.{html,json}

are copied to <record_dir>/dft/ (same layout; "+" in names becomes "_plus") and <record_dir>/dft/manifest.json
lists every file with its size and whether it is S3-only.  The cubes and the
dashboards are large: they are git-ignored and served from S3 only
(upload_large_files.py); git keeps relaxed.xyz, webapp.json and the manifest.
ingest_library.py reads the manifest, never the S3-only files, which a CI
checkout does not have.

Usage
  python import_dft_results.py <chain_dir> <record_dir>
"""
from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

PROVENANCE = {"code": "CP2K", "functional": "PBE", "basis": "DZVP-MOLOPT"}
ORBITALS = ("HOMO-1", "HOMO", "LUMO", "LUMO+1")
DASHBOARDS = (
    "exciton_analysis_sf.html", "exciton_analysis_soc.html",
    "fuzzy_dashboard_sf.html", "fuzzy_dashboard_soc.html",
    "fuzzy_dashboard_sf_qp.html", "fuzzy_dashboard_soc_qp.html",
)


def element_counts(xyz: Path) -> Counter:
    """Element counts of the first frame of an xyz file."""
    lines = xyz.read_text().splitlines()
    n = int(lines[0].split()[0])
    return Counter(line.split()[0] for line in lines[2:2 + n])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("chain_dir", type=Path, help="DFT chain output (relaxed.xyz, qdex/)")
    ap.add_argument("record_dir", type=Path, help="library record folder (record.json, start.xyz)")
    ap.add_argument("--force", action="store_true",
                    help="publish even if the relaxed composition differs from start.xyz")
    args = ap.parse_args(argv)
    chain, record = args.chain_dir, args.record_dir
    qdex = chain / "qdex"

    if not (record / "record.json").is_file():
        ap.error(f"{record} has no record.json")
    for need in (chain / "relaxed.xyz", qdex / "webapp.json"):
        if not need.is_file():
            ap.error(f"missing {need}")
    # Relaxation never changes the composition: a mismatch means the wrong record.
    if (record / "start.xyz").is_file() and not args.force:
        got, want = element_counts(chain / "relaxed.xyz"), element_counts(record / "start.xyz")
        if got != want:
            ap.error(f"composition of {chain / 'relaxed.xyz'} ({dict(got)}) differs from "
                     f"{record / 'start.xyz'} ({dict(want)}); --force to publish anyway")

    # (source, path inside dft/, S3-only, kind, orbital label)
    files = [(chain / "relaxed.xyz", "relaxed.xyz", False, "geometry", None),
             (qdex / "webapp.json", "qdex/webapp.json", False, "summary", None)]
    for label in ORBITALS:
        src = qdex / f"spatial_MO_{label}.cube"
        if src.is_file():
            # S3 reads "+" in a URL path as a space: LUMO+1 -> spatial_MO_LUMO_plus1.cube
            files.append((src, f"qdex/{src.name.replace('+', '_plus')}", True, "orbital", label))
    for name in DASHBOARDS:
        if (qdex / name).is_file():
            files.append((qdex / name, f"qdex/{name}", True, "dashboard", None))
    trap = chain / "qdex_trap"               # same orbitals, BSE without trap-state transitions
    if (trap / "webapp.json").is_file():
        files.append((trap / "webapp.json", "qdex/webapp_traps_filtered.json", False, "summary", None))
        for name in ("exciton_analysis_sf.html", "exciton_analysis_soc.html"):
            if (trap / name).is_file():
                files.append((trap / name, f"qdex/{name[:-5]}_traps_filtered.html", True, "dashboard", None))

    out = record / "dft"
    entries = []
    for src, rel, s3_only, kind, label in files:
        dest = out / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
        entry = {"path": rel, "size": dest.stat().st_size, "s3_only": s3_only, "kind": kind}
        if label:
            entry["orbital"] = label
        entries.append(entry)

    summary = json.loads((qdex / "webapp.json").read_text()).get("summary", {})
    manifest = {
        "schema_version": 1,
        **PROVENANCE,
        "qdex_commit": summary.get("qdex_commit"),
        "imported": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "files": entries,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    n_s3 = sum(e["s3_only"] for e in entries)
    print(f"{len(entries)} files -> {out} ({n_s3} S3-only, "
          f"{sum(e['size'] for e in entries if e['s3_only']) / 1e6:.1f} MB); manifest -> {out / 'manifest.json'}")
    missing = [lab for lab in ORBITALS if not any(e.get("orbital") == lab for e in entries)]
    if missing:
        print(f"warning: no cube for {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
