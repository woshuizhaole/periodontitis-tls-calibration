#!/usr/bin/env python3
import os, hashlib, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"
def sha_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as fh:
        for chunk in iter(lambda: fh.read(1024*1024*16), b""): h.update(chunk)
    return h.hexdigest()
def sha_tree(d):
    h=hashlib.sha256(); items=[]
    for root,dirs,files in os.walk(d):
        for f in files:
            p=os.path.join(root,f); rel=os.path.relpath(p,d)
            items.append((rel,os.path.getsize(p),sha_file(p)))
    for rel,size,s in sorted(items):
        h.update(f"{rel}\0{size}\0{s}\n".encode())
    return h.hexdigest(), len(items), sum(x[1] for x in items)
d=pd.read_csv(os.path.join(ROOT,"metadata/datasets.tsv"),sep="\t"); d["checksum"]=d["checksum"].astype(object)
for g,sub in d.groupby("dataset_id"):
    p=os.path.join(ROOT,str(sub.file_path.iloc[0]))
    if os.path.isdir(p):
        h,n,sz=sha_tree(p); val=f"tree_sha256:{h};files:{n};bytes:{sz}"
    elif os.path.isfile(p):
        val=f"sha256:{sha_file(p)}"
    else:
        val="MISSING"
    d.loc[d.dataset_id==g,"checksum"]=val
    d.loc[d.dataset_id==g,"file_size"]=sz if os.path.isdir(p) else (os.path.getsize(p) if os.path.isfile(p) else 0)
d.to_csv(os.path.join(ROOT,"metadata/datasets.tsv"),sep="\t",index=False)
# fill manifest md5 column where possible
m=pd.read_csv(os.path.join(ROOT,"data/manifests/DATA_MANIFEST.tsv"),sep="\t"); m["md5_or_sha256"]=m["md5_or_sha256"].astype(object)
cmap={}
for _,r in d.drop_duplicates("dataset_id").iterrows():
    cmap[r.file_path]=r.checksum
for i,r in m.iterrows():
    rel=os.path.relpath(r.raw_path,ROOT) if str(r.raw_path).startswith(ROOT) else str(r.raw_path)
    if rel in cmap: m.at[i,"md5_or_sha256"]=cmap[rel]
m.to_csv(os.path.join(ROOT,"data/manifests/DATA_MANIFEST.tsv"),sep="\t",index=False)
print(d[["dataset_id","file_size","checksum"]].drop_duplicates().to_string(index=False))
