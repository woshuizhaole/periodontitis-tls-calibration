#!/usr/bin/env python3
import os, hashlib, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"
OUT=os.path.join(ROOT,"metadata/checksums.tsv")
rows=[]
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as fh:
        for chunk in iter(lambda: fh.read(1024*1024*16), b""): h.update(chunk)
    return h.hexdigest()
# xenium hashes already computed
xp=os.path.join(ROOT,"data/manifests/xenium_files.tsv")
if os.path.exists(xp):
    x=pd.read_csv(xp,sep="\t")
    for _,r in x.iterrows(): rows.append(dict(dataset="xenium",path=r.path,size=r["size"],sha256=r.sha256))
# geo tars
for root,dirs,files in os.walk(os.path.join(ROOT,"data/raw/geo")):
    for f in files:
        if f.endswith(".tar"):
            p=os.path.join(root,f); rows.append(dict(dataset="geo",path=p,size=os.path.getsize(p),sha256=sha(p)))
# ibex
for f in os.listdir(os.path.join(ROOT,"data/raw/ibex")):
    p=os.path.join(ROOT,"data/raw/ibex",f)
    if os.path.isfile(p) and f.endswith(".ims"):
        rows.append(dict(dataset="ibex",path=p,size=os.path.getsize(p),sha256=sha(p)))
# tonsil
for f in os.listdir(os.path.join(ROOT,"data/raw/tonsil")):
    p=os.path.join(ROOT,"data/raw/tonsil",f)
    if os.path.isfile(p):
        rows.append(dict(dataset="tonsil",path=p,size=os.path.getsize(p),sha256=sha(p)))
pd.DataFrame(rows).to_csv(OUT,sep="\t",index=False)
print(OUT,len(rows))
