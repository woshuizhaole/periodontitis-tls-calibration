#!/usr/bin/env python3
import os, hashlib, time
import pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"
OUT=os.path.join(ROOT,"data/manifests/DATA_MANIFEST.tsv"); os.makedirs(os.path.dirname(OUT),exist_ok=True)
rows=[]
def add(name, src, method, path, size, status="present", notes=""):
    rows.append(dict(dataset_name=name, source_url=src, download_method=method, downloaded_at=time.strftime("%Y-%m-%dT%H:%M:%S"),
                     raw_path=path, file_size=size, md5_or_sha256="", expected_files="", actual_files="", status=status, notes=notes))
# external symlinks
for f in sorted(os.listdir(os.path.join(ROOT,"data/external"))):
    p=os.path.join(ROOT,"data/external",f)
    if os.path.islink(p):
        tgt=os.path.realpath(p); add(f, tgt, "symlink", p, os.path.getsize(tgt) if os.path.exists(tgt) else 0,
                                     "present" if os.path.exists(tgt) else "broken", f"target={tgt}")
# geo tars
for root,dirs,files in os.walk(os.path.join(ROOT,"data/raw/geo")):
    for f in files:
        if f.endswith(".tar"):
            p=os.path.join(root,f); add(f, "", "wget", p, os.path.getsize(p))
# xenium files
for root,dirs,files in os.walk(os.path.join(ROOT,"data/raw/xenium")):
    for f in files:
        p=os.path.join(root,f); add("xenium_"+os.path.basename(root), "", "dataverse_api", p, os.path.getsize(p))
# ibex
for f in os.listdir(os.path.join(ROOT,"data/raw/ibex")):
    p=os.path.join(ROOT,"data/raw/ibex",f)
    if os.path.isfile(p): add("ibex_"+f, "", "dataverse_api", p, os.path.getsize(p))
# tonsil
for f in os.listdir(os.path.join(ROOT,"data/raw/tonsil")):
    p=os.path.join(ROOT,"data/raw/tonsil",f)
    if os.path.isfile(p): add("tonsil_"+f, "", "wget", p, os.path.getsize(p))
pd_rows=pd.DataFrame(rows)
pd_rows.to_csv(OUT,sep="\t",index=False)
print(OUT,len(pd_rows))
