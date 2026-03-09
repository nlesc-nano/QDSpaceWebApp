import os
import json

# Point to the new Svelte public folder
target_dir = "qd-frontend/public"
out_file = os.path.join(target_dir, "file_list.json")

xyz_files = []
for root, dirs, files in os.walk(target_dir):
    parts = root.split(os.sep)
    in_md_folder = any(part.lower() == "md" for part in parts)

    for f in files:
        if not f.lower().endswith(".xyz"):
            continue
        relpath = os.path.relpath(os.path.join(root, f), target_dir).replace("\\", "/")
        
        if in_md_folder:
            if "pos" in f.lower():
                xyz_files.append(relpath)
        else:
            xyz_files.append(relpath)

xyz_files.sort()

# Save as clean JSON
with open(out_file, "w") as out:
    json.dump(xyz_files, out, indent=2)

print(f"Generated {out_file} with {len(xyz_files)} .xyz files.")


