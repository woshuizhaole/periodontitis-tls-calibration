#!/usr/bin/env python3
"""Stage 1 replication/holdout using per-sample GSE171213 counts (fast) and GSE152042 h5ad."""
import os, glob, re, time, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, yaml
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage1"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage1_replication.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
def sig(n): return yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
def mean_z_score(adata, genes, name):
    genes=[g for g in genes if g in adata.var_names]
    if len(genes)<2: return None
    vals=[]
    for g in genes:
        x=np.asarray(adata[:,g].X.todense() if sp.issparse(adata[:,g].X) else adata[:,g].X).ravel().astype(float)
        sd=x.std()
        if sd>0: vals.append((x-x.mean())/sd)
    adata.obs[name]=np.mean(vals,axis=0); return len(genes)
# ---------- GSE171213 per sample ----------
log("GSE171213 per-sample load")
parts=[]
for f in sorted(glob.glob(os.path.join(ROOT,"data/external/GSE171213_raw","GSM*_*.counts.tsv.gz"))):
    if "AllSample" in f: continue
    sample=os.path.basename(f).split("_",1)[1].replace(".counts.tsv.gz","")
    df=pd.read_csv(f,sep="\t",index_col=0)
    X=sp.csr_matrix(df.values.T.astype(np.float32))
    obs=pd.DataFrame({"sample":sample,"donor":sample,"condition":("healthy" if sample.startswith("HC") else ("treated" if sample.startswith("PDT") else "periodontitis"))},index=df.columns)
    a=ad.AnnData(X=X,obs=obs,var=pd.DataFrame(index=df.index.astype(str))); a.var_names_make_unique()
    sc.pp.filter_cells(a,min_genes=200); sc.pp.filter_genes(a,min_cells=3)
    a.var["mt"]=a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True)
    a=a[(a.obs.pct_counts_mt<20)&(a.obs.n_genes_by_counts<6000)].copy()
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    parts.append(a); log(f"{sample} {a.shape}")
a=ad.concat(parts,join="outer",merge="same"); a.obs_names_make_unique(); log(f"concat {a.shape}")
for n in ["OTLS_core","OTLS_ext","OTLS_panel","Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]:
    k=mean_z_score(a,sig(n),n); log(f"score {n} genes={k}")
a.obs.to_csv(os.path.join(OUT,"GSE171213_cell_scores.tsv"),sep="\t")
summ=a.obs.groupby(["sample","condition"]).agg(n=("OTLS_core","size"),OTLS_core=("OTLS_core","mean"),OTLS_panel=("OTLS_panel","mean"),Bcell=("Bcell","mean"),Plasma=("Plasma","mean"),Tcell=("Tcell","mean")).reset_index()
summ.to_csv(os.path.join(OUT,"GSE171213_sample_summary.tsv"),sep="\t",index=False); log(summ.to_string(index=False))
# ---------- GSE152042 ----------
log("GSE152042 load")
a2=ad.read_h5ad(os.path.join(ROOT,"data/external/GSE152042_labeled.h5ad"))
for n in ["OTLS_core","OTLS_ext","OTLS_panel","Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]:
    k=mean_z_score(a2,sig(n),n); log(f"score {n} genes={k}")
a2.obs.to_csv(os.path.join(OUT,"GSE152042_cell_scores.tsv"),sep="\t")
summ2=a2.obs.groupby(["sample","condition"]).agg(n=("OTLS_core","size"),OTLS_core=("OTLS_core","mean"),OTLS_panel=("OTLS_panel","mean"),Bcell=("Bcell","mean"),Plasma=("Plasma","mean"),Tcell=("Tcell","mean")).reset_index()
summ2.to_csv(os.path.join(OUT,"GSE152042_sample_summary.tsv"),sep="\t",index=False); log(summ2.to_string(index=False))
log("DONE")
