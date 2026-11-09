import pandas as pd, numpy as np, os
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/narrative_B/spatial")
d=pd.read_csv(os.path.join(ROOT,"diagnostics/gse297924/spot_scores.tsv"),sep="\t",index_col=0)
for c in ["Bcell","Plasma","Tcell"]: d[c]=(d[c]-d[c].mean())/(d[c].std()+1e-12)
d["BPT"]=d[["Bcell","Plasma","Tcell"]].mean(1)
rows=[]
for s,g in d.groupby("sample"):
    rows.append(dict(dataset="GSE297924",sample=s,n=len(g),BPT_mean=float(g.BPT.mean()),OTLS_core_mean=float(g.OTLS_core.mean()),rho_BPT_OTLS=float(g.BPT.corr(g.OTLS_core,method="spearman")),BPT_top5_OTLS=float(g.loc[g.BPT>=g.BPT.quantile(.95),"OTLS_core"].mean()),rest_OTLS=float(g.loc[g.BPT<g.BPT.quantile(.95),"OTLS_core"].mean())))
r=pd.DataFrame(rows); r.to_csv(os.path.join(OUT,"GSE297924_BPT_summary.tsv"),sep="\t",index=False); print(r.to_string(index=False))
