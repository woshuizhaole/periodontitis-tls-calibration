#!/usr/bin/env python3
"""P2.3: GSE16134 OTLS effect within B/Plasma/T composition strata."""
import os, pandas as pd, numpy as np
from scipy.stats import wilcoxon
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/narrative_B/bulk"); os.makedirs(OUT,exist_ok=True)
sc=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE16134_scores.tsv"),sep="\t",index_col=0)
meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t").set_index("gsm")
d=sc.join(meta[["patient_id","site_status"]],how="inner")
for c in ["OTLS_core","Bcell","Plasma","Tcell"]: d[c]=(d[c]-d[c].mean())/(d[c].std()+1e-12)
d["BPT"]=d[["Bcell","Plasma","Tcell"]].mean(1)
# patient-level BPT (mean across sites) and OTLS difference
piv=d.pivot_table(index="patient_id",columns="site_status",values=["OTLS_core","BPT"],aggfunc="mean").dropna()
piv["BPT_mean"]=piv[("BPT","affected")]+piv[("BPT","unaffected")]; piv["BPT_mean"]/=2
piv["OTLS_diff"]=piv[("OTLS_core","affected")]-piv[("OTLS_core","unaffected")]
piv["BPT_diff"]=piv[("BPT","affected")]-piv[("BPT","unaffected")]
piv["stratum"]=pd.qcut(piv.BPT_mean,4,labels=["Q1_low","Q2","Q3","Q4_high"])
rows=[]
for s,g in piv.groupby("stratum",observed=True):
    rows.append(dict(stratum=str(s),n=len(g),OTLS_diff_mean=float(g.OTLS_diff.mean()),OTLS_diff_median=float(g.OTLS_diff.median()),BPT_diff_mean=float(g.BPT_diff.mean()),p=float(wilcoxon(g.OTLS_diff).pvalue)))
allrow=dict(stratum="all",n=len(piv),OTLS_diff_mean=float(piv.OTLS_diff.mean()),OTLS_diff_median=float(piv.OTLS_diff.median()),BPT_diff_mean=float(piv.BPT_diff.mean()),p=float(wilcoxon(piv.OTLS_diff).pvalue))
res=pd.DataFrame([allrow]+rows); res.to_csv(os.path.join(OUT,"BPT_stratified_OTLS.tsv"),sep="\t",index=False)
print(res.to_string(index=False)); print(piv.groupby("stratum",observed=True).BPT_diff.mean().to_string())
