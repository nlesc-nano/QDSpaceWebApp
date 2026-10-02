#!/usr/bin/env python3
"""
Build the unified library index (one record per structure).

Sources
  * the existing library tree (qd-frontend/public/<family>/<material>/...):
    files are grouped into structures (start / geo_opt / md stages plus a
    properties/ folder). A `meta.yaml` next to the stage folders describes
    what cannot be read from coordinates (source, surface, DFT settings).
    Without one, it is prefilled from the legacy filename parsing
    (make_metadata.parse_metadata) and, with --write-meta, written out.
  * builder series produced by QD_builder
    (python -m builder.scripts.generate_library ...), one record.json each.

Everything computable from coordinates (formula, centre species, size,
fingerprint) comes from builder.library_record, the same code the generator
uses. DFT structures whose start geometry has the fingerprint of a builder
structure are merged with it (the record gains the recipe and the DFT stages).

Usage
  python ingest_library.py --builder <generator out>/CdSe --out library_index.json \
      --report library_report.md [--write-meta]
"""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import yaml

from builder.library_record import (
    SCHEMA_VERSION,
    describe_structure,
    make_record,
    radial_signature,
    read_xyz_first_frame,
    signature_distance,
)
from make_metadata import MATERIAL_ELEMENTS, find_xyz_files, parse_metadata

PUBLIC = Path("qd-frontend/public")
STAGE_DIRS = {"start": "start", "geo_opt": "geo_opt", "geo_op": "geo_opt", "md": "md"}
# Twin = same material, formula, centre and surface; the mean radial-profile
# deviation (Å) only guards against a different cut.  DFT start files are
# pre-relaxed copies of builder cuts (0.2-0.4 Å on small dots).
TWIN_TOL = 0.5
SAME_START_TOL = 0.02
FORMAL_CHARGES = {
    "Cd": 2, "Zn": 2, "Pb": 2, "Hg": 2, "In": 3, "Ga": 3, "Al": 3, "Cs": 1, "Rb": 1,
    "S": -2, "Se": -2, "Te": -2, "O": -2, "P": -3, "As": -3, "Sb": -3,
    "F": -1, "Cl": -1, "Br": -1, "I": -1,
}
FAMILY_PHASE = {
    "II-VI": "zinc-blende", "III-V": "zinc-blende", "IV-VI": "rock-salt",
    "ABX3": "cubic perovskite", "II-VI@II-VI": "zinc-blende",
}
# Prefill for meta.yaml only (the meta.yaml is the source of truth afterwards):
# manually built, reconstructed DFT structures.
LEGACY_SURFACE_OVERRIDES = {
    "II-VI/CdSe/HLE17/34ang": "reconstructed",
    "II-VI/CdSe/HLE17/40ang": "reconstructed",
}


def native_elements(material: str) -> List[str]:
    """Native (inorganic core) elements, in formula order; e.g. CdSe_ZnS -> Cd, Se, Zn, S."""
    out: List[str] = []
    for p in material.replace("@", "_").split("_"):
        elems = MATERIAL_ELEMENTS.get(p.upper()) or [
            e for e in re.findall(r"[A-Z][a-z]?", p) if e in FORMAL_CHARGES
        ]
        for e in elems:
            if e not in out:
                out.append(e)
    return out


def group_legacy_files(files: List[str]) -> Dict[str, Dict[str, List[str]]]:
    """structure dir -> {stage: [relpaths]}; files outside stage folders stand alone."""
    groups: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    for rel in files:
        parts = rel.split("/")
        if len(parts) >= 2 and parts[-2] in STAGE_DIRS:
            groups["/".join(parts[:-2])][STAGE_DIRS[parts[-2]]].append(rel)
        else:
            low = parts[-1].lower()
            stage = "md" if "pos" in low else ("geo_opt" if "opt" in low else "start")
            groups[rel][stage].append(rel)
    return groups


def legacy_meta(group: str, stages: Dict[str, List[str]]) -> dict:
    first = (stages.get("geo_opt") or stages.get("start") or stages.get("md"))[0]
    legacy = parse_metadata(first)
    family = legacy["system_type"]
    return {
        "source": "dft",
        "surface": LEGACY_SURFACE_OVERRIDES.get(group, "clean"),
        "family": family,
        "material": legacy["material"],
        "phase": FAMILY_PHASE.get(family, "unknown"),
        "dft": {"code": legacy.get("code") or None, "functional": legacy.get("functional") or None,
                "basis": legacy.get("basis") or None},
        "notes": "",
    }


def legacy_record(group: str, stages: Dict[str, List[str]], meta: dict):
    ref = (stages.get("start") or stages.get("geo_opt") or stages.get("md"))[0]
    symbols, pts = read_xyz_first_frame(str(PUBLIC / ref))
    sig = radial_signature(symbols, pts)
    desc = describe_structure(symbols, pts, native_order=native_elements(meta["material"]),
                              charges=FORMAL_CHARGES)
    dft = meta.get("dft", {})
    stage_list = []
    for stage in ("start", "geo_opt", "md"):
        for rel in sorted(stages.get(stage, [])):
            entry = {"stage": stage, "file": rel}
            if stage != "start":
                entry.update({k: v for k, v in dft.items() if v})
            stage_list.append(entry)
    prop_dir = PUBLIC / group / "properties" if (PUBLIC / group).is_dir() else None
    properties = {}
    if prop_dir and prop_dir.is_dir():
        for f in sorted(prop_dir.iterdir()):
            if f.suffix in (".html", ".gz", ".png"):
                properties[f.name] = str(f.relative_to(PUBLIC))
    rec = make_record(
        material=meta["material"], family=meta["family"], phase=meta["phase"],
        surface=meta["surface"], source=meta["source"], description=desc, stages=stage_list,
        origin={"legacy_paths": [group], "notes": meta.get("notes", "")},
        extra={"properties": properties},
    )
    return rec, sig, bool(stages.get("start"))


def load_builder_records(dirs: List[str], prefix: str):
    """
    Builder records with their start-geometry radial signatures.  Series
    stored inside the site root get paths relative to it; others are
    addressed as <prefix>/<material>/<id>.
    """
    recs = []
    public = PUBLIC.resolve()
    for d in dirs:
        for rec_path in sorted(Path(d).glob("*/record.json")):
            rec = json.loads(rec_path.read_text())
            symbols, pts = read_xyz_first_frame(str(rec_path.parent / "start.xyz"))
            folder = rec_path.parent.resolve()
            if public in folder.parents:
                base = str(folder.relative_to(public))
            else:
                base = f"{prefix.rstrip('/')}/{rec['material']}/{rec['id']}"
            for st in rec["stages"]:
                st["file"] = f"{base}/{st['file']}"
            rec.setdefault("properties", {})
            _attach_qdprops(rec, rec_path.parent, base)
            recs.append((rec, radial_signature(symbols, pts)))
    return recs


def _attach_qdprops(rec: dict, folder: Path, base: str) -> None:
    """Fold QD_builder `qdprops` results (<id>/props/) into a builder record."""
    props = folder / "props"
    summary_path = props / "properties.json"
    if not summary_path.is_file():
        return
    summary = json.loads(summary_path.read_text())
    head = (summary.get("provenance", {}).get("relax") or {}).get("head", "")
    if (props / "relaxed.xyz").is_file():
        rec["stages"].append({"stage": "geo_opt", "file": f"{base}/props/relaxed.xyz",
                              "functional": f"MACE-MH-1 ({head})" if head else "MACE-MH-1", "code": "MACE"})
    for name in ("ground_state.html", "ground_state.png", "synthesis.html"):
        if (props / name).is_file():
            rec["properties"][name] = f"{base}/props/{name}"
    rec["computed"] = {"schema_version": summary.get("schema_version"), "summary": summary.get("summary", {}),
                       "provenance": summary.get("provenance", {})}


def _best_match(rec: dict, sig, pool, tol: float):
    """Closest (record, distance) in pool with the same material/phase/formula/centre/surface."""
    best = None
    for other, other_sig in pool:
        if (other["material"], other.get("phase"), other["formula"], other["centre"], other["surface"]) != (
            rec["material"], rec.get("phase"), rec["formula"], rec["centre"], rec["surface"]
        ):
            continue
        d = signature_distance(sig, other_sig)
        if d <= tol and (best is None or d < best[1]):
            best = (other, d)
    return best


def _merge_into(target: dict, rec: dict, *, keep_start: bool) -> None:
    for st in rec["stages"]:
        if st["stage"] == "start" and not keep_start:
            continue
        if st not in target["stages"]:
            target["stages"].append(st)
    target["properties"].update(rec.get("properties", {}))
    paths = target.setdefault("origin", {}).setdefault("legacy_paths", [])
    for p in rec["origin"].get("legacy_paths", []):
        if p not in paths:
            paths.append(p)


def flags(rec: dict) -> dict:
    stages = {s["stage"] for s in rec["stages"]}
    return {"optimized": "geo_opt" in stages, "md": "md" in stages,
            "properties": bool(rec.get("properties"))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--builder", nargs="*", default=None,
                    help="generator output dirs (…/<material>); default: every */builder dir in the site root")
    ap.add_argument("--builder-prefix", default="library_builder",
                    help="path prefix under the site root where builder series are served")
    ap.add_argument("--out", default=str(PUBLIC / "library_index.json"))
    ap.add_argument("--report", default="library_report.md")
    ap.add_argument("--write-meta", action="store_true", help="write prefilled meta.yaml files")
    args = ap.parse_args(argv)

    # builder/ holds the default recipe; builder_<tag>/ holds variant recipes
    # (e.g. builder_100).  Sorted, so the default series comes first.
    builder_dirs = args.builder if args.builder is not None else [
        str(p) for p in sorted(PUBLIC.glob("*/*/builder*")) if p.is_dir()
    ]
    builder = load_builder_records(builder_dirs, args.builder_prefix)
    # A variant recipe can cut the same structure as the default one: keep the
    # first and note the other recipe on it.
    by_fp: Dict[str, dict] = {}
    unique_builder = []
    builder_dupes: List[str] = []
    for rec, sig in builder:
        first = by_fp.get(rec["fingerprint"])
        if first is not None:
            first["origin"].setdefault("also_from_recipes", []).append(rec["origin"].get("recipe"))
            builder_dupes.append(f"{rec['origin'].get('config')}: {rec['id']} -> {first['id']}")
            continue
        by_fp[rec["fingerprint"]] = rec
        unique_builder.append((rec, sig))
    builder = unique_builder

    legacy: List[tuple] = []   # (record, signature, has_start) after merging identical starts
    report = {"builder_dupes": builder_dupes, "legacy_groups": 0, "merged": [], "same_start": [], "meta_written": 0,
              "collisions": [], "errors": []}
    groups = group_legacy_files([
        f for f in find_xyz_files(str(PUBLIC))
        if not (PUBLIC / f).with_name("record.json").is_file()   # builder series: record-driven
        and "/props/" not in f"/{f}"                              # qdprops outputs of a record
    ])
    for group, stages in sorted(groups.items()):
        report["legacy_groups"] += 1
        meta_path = PUBLIC / group / "meta.yaml"
        if meta_path.is_file():
            meta = yaml.safe_load(meta_path.read_text())
        else:
            meta = legacy_meta(group, stages)
            if args.write_meta and (PUBLIC / group).is_dir():
                meta_path.write_text(yaml.safe_dump(meta, sort_keys=False))
                report["meta_written"] += 1
        try:
            rec, sig, has_start = legacy_record(group, stages, meta)
        except Exception as exc:
            report["errors"].append(f"{group}: {type(exc).__name__}: {exc}")
            continue
        # Same start geometry computed with another functional -> one structure.
        same = _best_match(rec, sig, [(r, s_) for r, s_, h in legacy if h], SAME_START_TOL) if has_start else None
        if same is not None:
            _merge_into(same[0], rec, keep_start=False)
            report["same_start"].append(f"{group} -> {same[0]['origin']['legacy_paths'][0]}")
            continue
        legacy.append((rec, sig, has_start))

    records: List[dict] = []
    for rec, sig, has_start in legacy:
        twin = _best_match(rec, sig, builder, TWIN_TOL) if has_start else None
        if twin is not None:
            target, dist = twin
            for st in rec["stages"]:
                if st["stage"] == "start":
                    st["dft_start"] = True
            _merge_into(target, rec, keep_start=True)
            target["source"] = "builder+dft"
            target["origin"]["dft_twin_deviation_A"] = round(dist, 3)
            report["merged"].append(
                f"{', '.join(rec['origin']['legacy_paths'])} -> {target['id']} (Δr̄ = {dist:.3f} Å)")
            continue
        records.append(rec)
    # Builder records first: they keep the canonical ids on formula clashes.
    records = [r for r, _ in builder] + records

    seen: Dict[str, int] = defaultdict(int)
    for rec in records:
        seen[rec["id"]] += 1
        if seen[rec["id"]] > 1:
            new_id = f"{rec['id']}-v{seen[rec['id']]}"
            where = ", ".join(rec["origin"].get("legacy_paths", [])) or "builder"
            report["collisions"].append(f"{rec['id']} -> {new_id} ({where})")
            rec["id"] = new_id
        rec["flags"] = flags(rec)

    index = {"schema_version": SCHEMA_VERSION,
             "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "structures": records}
    Path(args.out).write_text(json.dumps(index, indent=1))
    _write_report(Path(args.report), records, report)
    print(f"{len(records)} structures -> {args.out}; report -> {args.report}")
    return 0


def _write_report(path: Path, records: List[dict], rep: dict) -> None:
    by_mat = defaultdict(list)
    for r in records:
        by_mat[(r["family"], r["material"])].append(r)
    lines = ["# Library index report", "",
             f"- builder variant-recipe duplicates (merged): {len(rep['builder_dupes'])}",
             f"- legacy structure groups: {rep['legacy_groups']}",
             f"- same start, other functional (merged): {len(rep['same_start'])}",
             f"- merged with builder twins: {len(rep['merged'])}",
             f"- meta.yaml written: {rep['meta_written']}",
             f"- id collisions (suffixed): {len(rep['collisions'])}",
             f"- errors: {len(rep['errors'])}", "",
             "| family | material | structures | builder | dft | builder+dft | reconstructed | centres |",
             "|---|---|---|---|---|---|---|---|"]
    for (fam, mat), rs in sorted(by_mat.items()):
        src = defaultdict(int)
        for r in rs:
            src[r["source"]] += 1
        centres = ", ".join(f"{c}:{n}" for c, n in sorted(
            {c: sum(1 for r in rs if r["centre"] == c) for c in {r["centre"] for r in rs}}.items()))
        lines.append(f"| {fam} | {mat} | {len(rs)} | {src['builder']} | {src['dft']} | {src['builder+dft']} | "
                     f"{sum(1 for r in rs if r['surface'] == 'reconstructed')} | {centres} |")
    for title, key in (("Builder variant-recipe duplicates", "builder_dupes"),
                       ("Merged with builder twins", "merged"),
                       ("Same start geometry, other functional", "same_start"),
                       ("Id collisions", "collisions"),
                       ("Errors", "errors")):
        if rep[key]:
            lines += ["", f"## {title}", ""] + [f"- {x}" for x in rep[key]]
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
