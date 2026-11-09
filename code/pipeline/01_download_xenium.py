#!/usr/bin/env python3
"""Download Xenium matrices/cells for the Nat Immunol 2026 deposit.

Uses the Harvard Dataverse API; treats `-1` file variants as separate sections.
"""
import hashlib
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

ROOT = "/data/workplace/wyx/perio_tls_20261107"
OUT = os.path.join(ROOT, "data/raw/xenium")
LOG = os.path.join(ROOT, "logs/downloads/xenium_download.log")
DOI = "doi:10.7910/DVN/F2ICVP"
API = "https://dataverse.harvard.edu/api"

os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.dirname(LOG), exist_ok=True)


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def fetch_manifest():
    url = f"{API}/datasets/:persistentId/?persistentId={DOI}"
    r = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=120)
    r.raise_for_status()
    data = r.json()["data"]["latestVersion"]["files"]
    rows = []
    for item in data:
        df = item["dataFile"]
        rows.append({"id": df["id"], "name": df["filename"], "size": df.get("filesize", 0)})
    return rows


def parse_sections(rows):
    sections = {}
    pat = re.compile(r"^(.+?)_matrix\.mtx(-1)?\.gz$")
    for row in rows:
        m = pat.match(row["name"])
        if not m:
            continue
        key, variant = m.group(1), (m.group(2) or "")
        sid = key + variant
        sections[sid] = {"key": key, "variant": variant, "files": {}}
    # attach exact files
    for sid, meta in sections.items():
        key, variant = meta["key"], meta["variant"]
        for typ in ["features.tsv", "barcodes.tsv", "matrix.mtx", "cells.csv", "cell_boundaries.csv"]:
            # actual filename: key_typ + variant + .gz, e.g. P1_cells.csv-1.gz
            want = f"{key}_{typ}{variant}.gz"
            matches = [r for r in rows if r["name"] == want]
            if not matches:
                # sometimes variant is absent on some auxiliary files; try base
                want2 = f"{key}_{typ}.gz"
                matches = [r for r in rows if r["name"] == want2]
            if matches:
                meta["files"][typ] = matches[0]
    return sections


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024 * 16), b""):
            h.update(chunk)
    return h.hexdigest()


def download(sid, typ, row):
    d = os.path.join(OUT, sid)
    os.makedirs(d, exist_ok=True)
    # normalized local name
    local_name = typ + ".gz"
    out = os.path.join(d, local_name)
    if os.path.exists(out) and os.path.getsize(out) == row["size"]:
        return sid, typ, out, row["size"], "exists"
    url = f"{API}/access/datafile/{row['id']}"
    for attempt in range(1, 6):
        try:
            with requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, stream=True, timeout=300) as r:
                r.raise_for_status()
                tmp = out + ".part"
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            fh.write(chunk)
                os.replace(tmp, out)
            got = os.path.getsize(out)
            if got != row["size"]:
                raise IOError(f"size mismatch {got} != {row['size']}")
            return sid, typ, out, got, "downloaded"
        except Exception as e:
            log(f"RETRY {sid} {typ} attempt {attempt}: {e}")
            time.sleep(3 * attempt)
    raise RuntimeError(f"failed {sid} {typ}")


def main():
    rows = fetch_manifest()
    log(f"manifest rows={len(rows)}")
    sections = parse_sections(rows)
    log(f"sections={len(sections)}: {sorted(sections)}")
    jobs = []
    for sid, meta in sorted(sections.items()):
        for typ, row in meta["files"].items():
            jobs.append((sid, typ, row))
        missing = set(["features.tsv", "barcodes.tsv", "matrix.mtx", "cells.csv", "cell_boundaries.csv"]) - set(meta["files"])
        if missing:
            log(f"WARN {sid} missing {sorted(missing)}")
    log(f"download jobs={len(jobs)}")
    manifest = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(download, *job): job for job in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            sid, typ, out, size, status = fut.result()
            manifest.append({"section": sid, "type": typ, "path": out, "size": size, "status": status})
            if i % 20 == 0 or i == len(jobs):
                log(f"progress {i}/{len(jobs)}")
    # hashes
    log("hashing")
    for rec in manifest:
        rec["sha256"] = sha256(rec["path"])
    with open(os.path.join(ROOT, "data/manifests/xenium_files.tsv"), "w") as fh:
        fh.write("section\ttype\tpath\tsize\tstatus\tsha256\n")
        for rec in manifest:
            fh.write(f"{rec['section']}\t{rec['type']}\t{rec['path']}\t{rec['size']}\t{rec['status']}\t{rec['sha256']}\n")
    log(f"done files={len(manifest)}")


if __name__ == "__main__":
    main()
