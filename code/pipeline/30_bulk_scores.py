#!/usr/bin/env python3
"""Stage 3a: ssGSEA scores for all bulk cohorts and signatures."""
import os, gzip, time, warnings, yaml
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage3/scores"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage3_bulk_scores.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
def sig(n): return yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
SIGS=["OTLS_core","OTLS_ext","OTLS_panel","OTLS_core_independent","Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]
def ssgsea(x, gene_sets, alpha=0.25):
    x=np.asarray(x,dtype=float); n=len(x); order=np.argsort(-x); genes=np.asarray(gene_sets["_genes"])
    out={}
    for name,gs in gene_sets.items():
        if name=="_genes": continue
        hit=np.isin(genes[order],list(gs))
        if hit.sum()<2: out[name]=np.nan; continue
        ranks=np.arange(1,n+1,dtype=float)
        w=np.where(hit,ranks**alpha,0.0); wh=w/w.sum()
        miss=(~hit).astype(float); wm=miss/miss.sum()
        es=np.cumsum(wh-wm); out[name]=float(es[np.argmax(np.abs(es))])
    return out
def score_matrix(df, cohort, meta=None):
    # df: genes x samples
    genes=df.index.astype(str).values
    sets={"_genes":genes}
    for s in SIGS:
        gs=[g for g in sig(s) if g in set(genes)]
        sets[s]=gs; log(f"{cohort} {s} genes={len(gs)}")
    rows=[]
    for c in df.columns:
        d=ssgsea(df[c].values,sets); d["sample"]=c; rows.append(d)
    sc=pd.DataFrame(rows).set_index("sample")
    sc.to_csv(os.path.join(OUT,f"{cohort}_scores.tsv"),sep="\t")
    return sc
# GSE16134/10334
for cohort,f in [("GSE16134","data/external/GSE16134_expr_gene.tsv"),("GSE10334","data/external/GSE10334_expr_gene.tsv")]:
    log(f"load {cohort}")
    df=pd.read_csv(os.path.join(ROOT,f),sep="\t",index_col=0)
    log(f"{cohort} {df.shape}")
    score_matrix(df,cohort)
# GSE173078 logCPM
log("load GSE173078")
df=pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_logCPM.tsv"),sep="\t",index_col=0)
log(f"GSE173078 {df.shape}"); score_matrix(df,"GSE173078")
# GSE223328 series matrix + GPL21185 probe->gene
log("load GSE223328")
mp={}
with gzip.open(os.path.join(ROOT,"data/external/GPL21185.txt"),"rt",errors="ignore") as fh:
    tab=False
    for line in fh:
        if line.startswith("!platform_table_begin"): tab=True; hdr=next(fh).rstrip("\n").split("\t"); continue
        if not tab: continue
        if line.startswith("!platform_table_end"): break
        parts=line.rstrip("\n").split("\t")
        if len(parts)<=5: continue
        g=parts[5].strip()
        if g and g!="": mp[parts[0]]=g
log(f"probe map={len(mp)}")
lines=[]; tab=False
with gzip.open(os.path.join(ROOT,"data/raw/geo/GSE223328/GSE223328_series_matrix.txt.gz"),"rt",errors="ignore") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"): tab=True; lines.append(next(fh)); continue
        if not tab: continue
        if line.startswith("!series_matrix_table_end"): break
        lines.append(line)
import io
df=pd.read_csv(io.StringIO("".join(lines)),sep="\t",index_col=0)
df.index=[mp.get(i,i) for i in df.index]
df=df.groupby(level=0).mean()
log(f"GSE223328 {df.shape}")
score_matrix(df,"GSE223328")
# metadata merged for later
meta=pd.read_csv(os.path.join(ROOT,"metadata/datasets.tsv"),sep="\t")
m=meta[meta.dataset_id=="GSE223328"][["sample_id","disease_state","notes"]].rename(columns={"sample_id":"sample","disease_state":"condition"})
m.to_csv(os.path.join(OUT,"GSE223328_meta.tsv"),sep="\t",index=False)
log("DONE")
