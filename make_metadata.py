import os
import json
import re
import itertools
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


def extract_coords_and_symbols(xyz_path, material_name=None):
    symbols = []
    coords = []
    
    # 1. Look up the allowed core elements for this specific material
    allowed_elements = None
    if material_name:
        m_name = material_name.upper()
        if m_name in MATERIAL_ELEMENTS:
            allowed_elements = {el.lower() for el in MATERIAL_ELEMENTS[m_name]}

    try:
        with open(xyz_path, 'r') as f:
            first = f.readline()
            if not first: return symbols, coords
            
            n_atoms = int(first.strip())
            f.readline() # skip comment
            
            for _ in range(n_atoms):
                line = f.readline()
                if not line: break
                parts = line.split()
                
                if len(parts) >= 4:
                    sym = parts[0]
                    
                    # 2. BLAZING FAST FILTER: Only parse if it's in the allowed elements (or if no material matched)
                    if allowed_elements is None or sym.lower() in allowed_elements:
                        symbols.append(sym)
                        # We only pay the cost of float() conversion for core atoms!
                        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                        
    except Exception as e:
        print(f"Error parsing {xyz_path}: {e}")
        
    return symbols, coords


def parse_metadata(relpath):
    """
    Given a relative path like "II-VI/ZnSe/HLE17/28ang/geo_opt/StartXYZ.xyz",
    extract:
      - system_type   → "II-VI"
      - material      → "ZnSe"
      - filename      → "StartXYZ.xyz"
      - size (in nm)  → parsed from “28ang” → 2.8
      - functional    → "HLE17" if present
      - basis         → parsed from filename if “DZVP”/“TZVP” or default to "DZVP" when functional=="HLE17"
      - run_type      → "Geometry Optimization" if any folder is "geo_opt",
                         "Molecular Dynamics" if any folder is "md",
                         "Start" if any folder is "start" or if neither geo_opt nor md appear
      - code          → "ORCA" if “orca” in filename; else default "CP2K"
    """
    parts = relpath.split('/')
    if parts and parts[0] == 'library':
        parts = parts[1:]
    filename = os.path.basename(relpath)
    metadata = {
        "system_type": parts[0] if len(parts) > 0 else "",
        "material": parts[1] if len(parts) > 1 else "",
        "filename": filename
    }

    # ─── Size in nm ───────────────────────────────────────────────────────
    nm_match  = re.search(r'(\d+(\.\d+)?)\s*nm', filename, re.IGNORECASE)
    ang_match = re.search(r'(\d+(\.\d+)?)\s*ang', filename, re.IGNORECASE)
    if nm_match:
        metadata["size"] = float(nm_match.group(1))
    elif ang_match:
        metadata["size"] = round(float(ang_match.group(1)) / 10.0, 3)
    else:
        metadata["size"] = None

    # ─── Functional ───────────────────────────────────────────────────────
    func_match = re.search(r'(HLE17|PBE|B3LYP|HSE06)', filename, re.IGNORECASE)
    metadata["functional"] = func_match.group(1).upper() if func_match else ""

    # ─── Basis Set ────────────────────────────────────────────────────────
    basis_match = re.search(r'(DZVP|TZVP)', filename, re.IGNORECASE)
    if basis_match:
        metadata["basis"] = basis_match.group(1).upper()
    elif metadata["functional"] == "HLE17":
        metadata["basis"] = "DZVP"
    else:
        metadata["basis"] = ""

    # ─── Run Type ─────────────────────────────────────────────────────────
    run_type = ""
    for part in parts:
        low = part.lower()
        if low == "geo_opt":
            run_type = "Geometry Optimization"
            break
        elif low == "md":
            run_type = "Molecular Dynamics"
            break
        elif low == "start":
            run_type = "Start"
            break
    # If neither geo_opt, md, nor start was found, default to "Start"
    if not run_type:
        run_type = "Start"
    metadata["run_type"] = run_type

    # ─── DFT Code ──────────────────────────────────────────────────────────
    if re.search(r'orca', filename, re.IGNORECASE):
        metadata["code"] = "ORCA"
    else:
        metadata["code"] = "CP2K"

    return metadata

def count_atoms(xyz_path):
    """
    Count atoms only from the first frame of an XYZ file.
    That way, for an MD “pos” file with multiple frames, we only count the first frame.
    """
    counts = {}
    try:
        with open(xyz_path, 'r') as f:
            first = f.readline()
            if not first:
                return counts
            try:
                n_atoms = int(first.strip())
            except ValueError:
                # If header is malformed, fall back to counting all lines after line 2
                lines = [first] + f.readlines()
                for line in lines[2:]:
                    parts = line.strip().split()
                    if parts:
                        el = parts[0]
                        counts[el] = counts.get(el, 0) + 1
                return counts

            # Skip comment line
            f.readline()

            # Read exactly n_atoms lines for the first frame
            for _ in range(n_atoms):
                line = f.readline()
                if not line:
                    break
                parts = line.strip().split()
                if parts:
                    el = parts[0]
                    counts[el] = counts.get(el, 0) + 1
    except Exception:
        pass

    return counts

def compute_all_ratios(counts):
    ratios = {}
    elements = [el for el in counts if counts[el] > 0]
    for el1, el2 in itertools.combinations(elements, 2):
        n1, n2 = counts.get(el1, 0), counts.get(el2, 0)
        if n2:
            ratios[f"{el1}/{el2}"] = round(n1 / n2, 3)
        if n1:
            ratios[f"{el2}/{el1}"] = round(n2 / n1, 3)
    return ratios

def find_xyz_files(root):
    xyz_paths = []
    for dirpath, dirnames, filenames in os.walk(root):
        parts = dirpath.split(os.sep)
        in_md_folder = any(p.lower() == "md" for p in parts)
        for f in filenames:
            if not f.lower().endswith(".xyz"):
                continue
            rel = os.path.relpath(os.path.join(dirpath, f), root).replace("\\", "/")
            if in_md_folder:
                if "pos" in f.lower():
                    xyz_paths.append(rel)
            else:
                xyz_paths.append(rel)
    xyz_paths.sort()
    return xyz_paths

def main():
    # Update target to the Svelte public folder
    target_dir = "qd-frontend/public"
    metadata_out = os.path.join(target_dir, "metadata.json")

    # 1. Load existing metadata to avoid recomputing
    existing_meta = {}
    if os.path.exists(metadata_out):
        try:
            with open(metadata_out, 'r') as f:
                existing_meta = json.load(f)
            print(f"Found existing metadata with {len(existing_meta)} entries. Will skip recomputing known files.")
        except Exception as e:
            print(f"Warning: Could not read existing metadata.json: {e}")

    xyz_files = find_xyz_files(target_dir)
    meta = {}
    
    computed_count = 0
    skipped_count = 0

    print("Processing structures...")

    for relpath in xyz_files:
        entry = parse_metadata(relpath)
        full_path = os.path.join(target_dir, relpath)
        
        # Uncomment if you have your count_atoms functions in the file
        # atom_counts = count_atoms(full_path)
        # entry["stoichiometry"] = atom_counts
        # entry["ratios"] = compute_all_ratios(atom_counts)
        
        # Check if we already have the size metrics for this exact file (from caching)
        if relpath in existing_meta and "size_metrics" in existing_meta[relpath]:
            entry["size_metrics"] = existing_meta[relpath]["size_metrics"]
            skipped_count += 1
        else:
            # Pass the material to the extraction function to filter early!
            material = entry.get("material")
            symbols, coords = extract_coords_and_symbols(full_path, material)
            
            if coords:
                # The coordinates are ALREADY filtered down to just the core now
                metrics = get_cluster_size_metrics(coords, symbols, material)
                entry["size_metrics"] = metrics
                computed_count += 1
                
                if computed_count % 10 == 0:
                    print(f"  ...computed {computed_count} new size metrics...")

        meta[relpath] = entry

    # 4. Final summary and save
    with open(metadata_out, "w") as out:
        json.dump(meta, out, indent=2)
        
    print("\n--- Metadata Generation Complete ---")
    print(f"Structures skipped (reused metrics): {skipped_count}")
    print(f"Structures computed (new metrics):   {computed_count}")
    print(f"Generated {metadata_out} with {len(meta)} entries.")

if __name__ == "__main__":
    main()

