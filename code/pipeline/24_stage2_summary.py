#!/usr/bin/env python3
import os, glob, pandas as pd, numpy as np
from scipy.stats import spearmanr
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage2")
rows=[]
for f in sorted(glob.glob(os.path.join(OUT,"visium","*_scores.tsv"))):
    d=pd.read_csv(f,sep="\t",index_col=0); sid=os.path.basename(f).replace("_scores.tsv","")
    for target in ["Bcell","Tcell","Plasma"]:
        r,p=spearmanr(d.OTLS_core,d[target])
        rows.append(dict(dataset="visium",sample=sid,target=target,rho=float(r),p=float(p),n=len(d)))
# IBEX vs Xenium donor OTLS
ib=pd.read_csv(os.path.join(OUT,"ibex","image_proxy.tsv"),sep="\t")
ib=ib.groupby("donor").agg(organization_proxy=("organization_proxy","mean"),CD20=("CD20","mean"),Ki67=("Ki67","mean"),CD138=("CD138","mean"),CD31=("CD31","mean")).reset_index()
st=pd.read_csv(os.path.join(OUT,"xenium","section_summary.tsv"),sep="\t")
st=st[st.main_quant].dropna(subset=["OTLS"])
xd=st.groupby("donor").apply(lambda g:np.average(g.OTLS,weights=g.roi_um2)).rename("OTLS_xenium").reset_index()
m=ib.merge(xd,on="donor")
for target in ["organization_proxy","CD20","Ki67","CD138","CD31"]:
    if len(m)>=5:
        r,p=spearmanr(m[target],m.OTLS_xenium); rows.append(dict(dataset="ibex_xenium",sample="donor",target=target,rho=float(r),p=float(p),n=len(m)))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"stage2_summary.tsv"),sep="\t",index=False)
print(pd.DataFrame(rows).to_string(index=False))
print("IBEX-Xenium merged donors",len(m)); print(m.to_string(index=False))
