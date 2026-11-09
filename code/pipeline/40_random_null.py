#!/usr/bin/env python3
"""Stage 4: random gene-set null for the GSE16134 paired effect."""
import os, time, warnings
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage4"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage4_random_null.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
expr=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_expr_gene.tsv"),sep="\t",index_col=0)
meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t").set_index("gsm")
sc=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE16134_scores.tsv"),sep="\t",index_col=0)
genes=expr.index.astype(str).values; n=len(genes); rng=np.random.default_rng(42)
sets=[rng.choice(n,size=30,replace=False) for _ in range(100)]
# precompute ranks per sample once
X=expr.values.astype(np.float32)
rows=[]
for j in range(100):
    gs=set(genes[sets[j]])
    hitmask=np.array([g in gs for g in genes])
    vals=[]
    for c in range(X.shape[1]):
        order=np.argsort(-X[:,c]); hit=hitmask[order]
        if hit.sum()<2: vals.append(np.nan); continue
        ranks=np.arange(1,n+1,dtype=float); w=np.where(hit,ranks**0.25,0.0)
        es=np.cumsum(w/w.sum()-((~hit).astype(float)/(~hit).sum()))
        vals.append(float(es[np.argmax(np.abs(es))]))
    s=pd.Series(vals,index=expr.columns)
    d=sc.join(meta[["patient_id","site_status"]],how="inner")
    d["null"]=s.reindex(d.index).values
    piv=d.pivot_table(index="patient_id",columns="site_status",values="null",aggfunc="mean").dropna()
    if len(piv)>=5:
        diff=piv["affected"]-piv["unaffected"]
        rows.append(dict(set=j,n_pairs=len(piv),mean_diff=float(diff.mean()),cohen_dz=float(diff.mean()/(diff.std()+1e-12))))
    if (j+1)%10==0: log(f"null {j+1}/100")
null=pd.DataFrame(rows); null.to_csv(os.path.join(OUT,"random_null_gse16134.tsv"),sep="\t",index=False)
obs=pd.read_csv(os.path.join(ROOT,"results/stage3/stats/gse16134_paired.tsv"),sep="\t").set_index("signature").loc["OTLS_core"]
p=(1+(null.mean_diff.abs()>=abs(obs.mean_diff)).sum())/(len(null)+1)
pd.DataFrame([dict(observed=float(obs.mean_diff),null_mean=float(null.mean_diff.mean()),null_sd=float(null.mean_diff.std()),
                   empirical_p_two_sided=float(p),n_null=len(null))]).to_csv(os.path.join(OUT,"random_null_summary.tsv"),sep="\t",index=False)
log(f"observed={obs.mean_diff:.3f} null mean={null.mean_diff.mean():.3f} sd={null.mean_diff.std():.3f} p={p:.3g}")
log("DONE")
