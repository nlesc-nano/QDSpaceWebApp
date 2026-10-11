#!/usr/bin/env python3
"""
Upload the S3-only DFT files (orbital cubes, QDEX dashboards) of the library.

Every dft/manifest.json under qd-frontend/public lists its files; those marked
"s3_only" are git-ignored and live only in the site bucket, under the same key
as their path in public/ (e.g. II-VI/CdSe/builder/<id>/dft/qdex/spatial_MO_HOMO.cube).
Each is gzipped and stored with Content-Encoding: gzip, so browsers get plain
text from fetch() and iframes.  The deploy sync (deploy-frontend.yml) excludes
these keys, so its --delete leaves them alone.

Dry run by default (lists what would be uploaded, does not contact S3);
--upload performs it and skips objects whose ETag already matches.  gzip runs
with mtime=0, so the same file always gives the same bytes and ETag.

Usage
  python upload_large_files.py [--upload] [--bucket qdspace-frontend-v1] [paths ...]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

PUBLIC = Path("qd-frontend/public")
BUCKET = "qdspace-frontend-v1"
CONTENT_TYPES = {".html": "text/html; charset=utf-8", ".cube": "text/plain; charset=utf-8"}
CACHE_CONTROL = "public, max-age=86400"


def s3_only_files(roots):
    """(local path, bucket key) of every S3-only file in the manifests under roots."""
    out = []
    for root in roots:
        for manifest in sorted(Path(root).glob("**/dft/manifest.json")):
            for entry in json.loads(manifest.read_text()).get("files", []):
                if entry.get("s3_only"):
                    path = manifest.parent / entry["path"]
                    out.append((path, path.resolve().relative_to(PUBLIC.resolve()).as_posix(), entry))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="*", default=[str(PUBLIC)],
                    help="folders searched for dft/manifest.json (default: the whole site root)")
    ap.add_argument("--upload", action="store_true", help="upload (default: dry run)")
    ap.add_argument("--bucket", default=BUCKET)
    args = ap.parse_args(argv)

    files = s3_only_files(args.paths)
    s3 = None
    if args.upload:
        import boto3
        from botocore.exceptions import ClientError
        s3 = boto3.client("s3")
    n_up = n_skip = n_missing = 0
    total = 0
    for path, key, entry in files:
        if not path.is_file():
            # Only the machine that ran import_dft_results.py has the file.
            print(f"missing  {path} (listed in its manifest)")
            n_missing += 1
            continue
        raw = path.read_bytes()
        if len(raw) != entry.get("size", len(raw)):
            print(f"warning  {path}: {len(raw)} bytes, manifest says {entry['size']}")
        body = gzip.compress(raw, compresslevel=9, mtime=0)
        etag = hashlib.md5(body).hexdigest()
        ctype = CONTENT_TYPES.get(path.suffix, "application/octet-stream")
        if s3 is None:
            print(f"would upload  s3://{args.bucket}/{key}  {len(raw) / 1e6:.2f} MB -> {len(body) / 1e6:.2f} MB gz  {ctype}")
            n_up += 1
            total += len(body)
            continue
        try:
            head = s3.head_object(Bucket=args.bucket, Key=key)
            if head.get("ETag", "").strip('"') == etag and head.get("ContentEncoding") == "gzip":
                n_skip += 1
                continue
        except ClientError as exc:
            if exc.response.get("Error", {}).get("Code") not in ("404", "NoSuchKey", "NotFound"):
                raise
        # Single-part put (the largest dashboard gzips to a few MB): ETag = MD5 of the body.
        s3.put_object(Bucket=args.bucket, Key=key, Body=body, ContentType=ctype, ContentEncoding="gzip",
                      CacheControl=CACHE_CONTROL)
        print(f"uploaded  s3://{args.bucket}/{key}  {len(body) / 1e6:.2f} MB gz")
        n_up += 1
        total += len(body)
    verb = "uploaded" if s3 is not None else "would upload (dry run; --upload to perform)"
    print(f"{n_up} files {verb}, {total / 1e6:.1f} MB gz; {n_skip} unchanged; {n_missing} missing locally")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
