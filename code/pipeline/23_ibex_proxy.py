#!/usr/bin/env python3
"""Stage 2d: IBEX image-level organization proxy (CD20/Ki67/CD138/CD31 area fractions)."""
import os, glob, re, time, warnings
import numpy as np, pandas as pd, h5py
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage2/ibex"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage2_ibex.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
def decode(x):
    if isinstance(x,np.ndarray):
        try: return b"".join(x.tolist()).decode(errors="ignore")
        except Exception: return str(x)
    return str(x)
markers={"CD20":"CD20","Ki67":"Ki67","CD138":"CD138","CD31":"CD31","HLA-DR":"HLA-DR","CD3":"CD3","CD68":"CD68","PanCK":"PanCK","Vimentin":"Vimentin"}
rows=[]
for p in sorted(glob.glob(os.path.join(ROOT,"data/raw/ibex/*.ims"))):
    name=os.path.basename(p); stem=os.path.splitext(name)[0]; m=re.search(r"([HP]\d+)",stem); donor=m.group(1) if m else stem; group=donor[0]
    try:
        with h5py.File(p,"r") as f:
            names=[]
            for k in sorted([x for x in f["DataSetInfo"].keys() if x.startswith("Channel")],key=lambda s:int(s.split()[-1])):
                names.append(decode(f["DataSetInfo"][k].attrs.get("Name","")))
            vals={}
            for mk in markers:
                idx=[i for i,n in enumerate(names) if re.search(rf"\b{re.escape(mk)}\b",n,re.I)]
                if not idx: continue
                arr=f[f"DataSet/ResolutionLevel 0/TimePoint 0/Channel {idx[0]}/Data"][0,::8,::8].astype(np.float32)
                # normalize by 99th percentile, threshold at 75th percentile of positive signal
                q99=np.percentile(arr,99)
                if q99<=0: vals[mk]=0.0; continue
                arr=arr/q99
                vals[mk]=float((arr>0.25).mean())
            vals.update(file=name,donor=donor,group=group); rows.append(vals)
            log(f"{name} {donor} {vals}")
    except Exception as e: log(f"ERROR {name}: {e}")
df=pd.DataFrame(rows); df.to_csv(os.path.join(OUT,"image_proxy.tsv"),sep="\t",index=False)
# composite organization proxy
cols=[c for c in ["CD20","Ki67","CD138","CD31"] if c in df]
if cols:
    z=df[cols].apply(lambda x:(x-x.mean())/(x.std()+1e-9))
    df["organization_proxy"]=z.mean(1); df.to_csv(os.path.join(OUT,"image_proxy.tsv"),sep="\t",index=False)
print(df.head().to_string(index=False)); print("rows",len(df)); log("DONE")
