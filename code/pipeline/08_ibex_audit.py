#!/usr/bin/env python3
import os, glob, re, h5py, numpy as np, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"metadata/ibex_audit.tsv")
def decode(x):
    if isinstance(x,np.ndarray):
        try: return b"".join(x.tolist()).decode(errors="ignore")
        except Exception: return str(x)
    return str(x)
rows=[]
for p in sorted(glob.glob(os.path.join(ROOT,"data/raw/ibex/*.ims"))):
    name=os.path.basename(p); stem=os.path.splitext(name)[0]; m=re.search(r"([HP]\d+)",stem); donor=m.group(1) if m else stem; group=donor[0]
    try:
        with h5py.File(p,"r") as f:
            names=[]
            for k in sorted([x for x in f["DataSetInfo"].keys() if x.startswith("Channel")], key=lambda s:int(s.split()[-1])):
                names.append(decode(f["DataSetInfo"][k].attrs.get("Name","")))
            d=f["DataSet/ResolutionLevel 0/TimePoint 0/Channel 0/Data"]
            shape=tuple(d.shape); z=shape[0]; xdim=shape[2]
            ext=None
            for key in ["ExtMax0","ExtMax1","ExtMax2"]:
                if key in f["DataSetInfo/Image"].attrs:
                    ext=float(decode(f["DataSetInfo/Image"].attrs[key]))
                    break
            pixel=(ext/xdim) if ext and xdim else np.nan
            rows.append(dict(file=name, donor=donor, group=group, size=os.path.getsize(p), shape="x".join(map(str,shape)),
                             z_planes=z, pixel_size_um=pixel, n_channels=len(names), channels=";".join(names)))
    except Exception as e:
        rows.append(dict(file=name, donor=donor, group=group, size=os.path.getsize(p), shape="", z_planes=0, pixel_size_um=np.nan, n_channels=0, channels="", error=str(e)))
pd.DataFrame(rows).to_csv(OUT,sep="\t",index=False)
print(pd.DataFrame(rows)[["file","donor","group","z_planes","pixel_size_um","n_channels"]].to_string(index=False))
