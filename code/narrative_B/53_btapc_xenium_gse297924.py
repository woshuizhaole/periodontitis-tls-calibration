#!/usr/bin/env python3
"""B score in Xenium and GSE297924 full Visium."""
import os, glob, gzip, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, anndata as ad, scanpy as sc
from scipy.stats import spearmanr, mannwhitneyu
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/narrative_B/spatial"); os.makedirs(OUT,exist_ok=True)
APC=["HLA-DRB5","HLA-DMB","HLA-DQA2","HLA-DQB2","CD14","CD68","CD163","CD1C","FCER1A","LAMP3","CCR7","CD83"]
# Xenium
nb=pd.read_parquet(os.path.join(ROOT,"diagnostics/xenium/cache/nb.parquet"))
nb["B_T_APC_min"]=np.minimum.reduce([nb.nb_B_50,nb.nb_T_50,nb.nb_APC_50])
nb["B_T_APC_mean"]=(nb.nb_B_50+nb.nb_T_50+nb.nb_APC_50)/3
cells=pd.concat([pd.read_csv(f,sep="\t") for f in glob.glob(os.path.join(ROOT,"results/stage2/xenium/cells/*.tsv.gz"))],ignore_index=True)
d=cells.merge(nb[["section","cell_id","B_T_APC_min","B_T_APC_mean"]],on=["section","cell_id"],how="left")
dm=pd.read_csv(os.path.join(ROOT,"metadata/xenium_donor_map.tsv"),sep="\t"); d=d.merge(dm[["section","donor","group","main_quant"]],on="section",how="left")
sec=d.groupby(["section","group","donor","main_quant"]).agg(B_T_APC_min=("B_T_APC_min","mean"),B_T_APC_mean=("B_T_APC_mean","mean"),OTLS_panel=("OTLS_panel","mean")).reset_index()
sec.to_csv(os.path.join(OUT,"xenium_BTAPC_section.tsv"),sep="\t",index=False)
rows=[]
for s in ["B_T_APC_min","B_T_APC_mean"]:
    x=sec[sec.main_quant]; h=x[x.group=="H"][s]; p=x[x.group=="P"][s]
    rows.append(dict(score=s,rho_OTLS=float(spearmanr(x[s],x.OTLS_panel).correlation),p_rho=float(spearmanr(x[s],x.OTLS_panel).pvalue),H_mean=float(h.mean()),P_mean=float(p.mean()),p_HvP=float(mannwhitneyu(h,p).pvalue),n=len(x)))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"xenium_BTAPC_HvP.tsv"),sep="\t",index=False)
# GSE297924
rows=[]
# GSE297924 BPT is in 58_gse297924_bpt.py (uses existing spot_scores)
from sklearn.metrics import roc_auc_score
import yaml
sig=lambda n: yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
def mean_z(a,genes,name):
    gs=[g for g in genes if g in a.var_names]
    if len(gs)<2: return
    vals=[]
    for g in gs:
        x=np.asarray(a[:,g].X.todense() if sp.issparse(a[:,g].X) else a[:,g].X).ravel().astype(float); sd=x.std()
        if sd>0: vals.append((x-x.mean())/sd)
    a.obs[name]=np.mean(vals,axis=0)
for dpath in sorted(glob.glob(os.path.join(ROOT,"diagnostics/data/tmp/gse297924/*/"))):
    sid=os.path.basename(dpath.rstrip("/")); mtx=os.path.join(dpath,"matrix.mtx.gz")
    if not os.path.exists(mtx): continue
    with gzip.open(mtx,"rb") as fh: X=scipy.io.mmread(fh).tocsr()
    f=pd.read_csv(os.path.join(dpath,"features.tsv.gz"),sep="\t",header=None); b=pd.read_csv(os.path.join(dpath,"barcodes.tsv.gz"),header=None)[0].astype(str).values
    a=ad.AnnData(X=X.T.tocsr().astype(np.float32),obs=pd.DataFrame(index=b),var=pd.DataFrame(index=f[1].astype(str).values)); a.var_names_make_unique()
    sc.pp.filter_cells(a,min_genes=200); sc.pp.filter_genes(a,min_cells=3); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    for n,gs in {"Bcell":["MS4A1","CD79A","CD19","BANK1"],"Tcell":["CD3D","CD3E","CD2","TRAC"],"APC":APC,"Plasma":sig("Plasma")}.items(): mean_z(a,gs,n)
    a.obs["BPTA_min"]=np.minimum.reduce([a.obs.Bcell,a.obs.Tcell,a.obs.APC,a.obs.Plasma]); a.obs["BPTA_mean"]=(a.obs.Bcell+a.obs.Tcell+a.obs.APC+a.obs.Plasma)/4
    a.obs.to_csv(os.path.join(OUT,f"GSE297924_{sid}_BPTA.tsv"),sep="\t")
    rows.append(dict(sample=sid,n=a.n_obs,BPTA_min=float(a.obs.BPTA_min.mean()),BPTA_mean=float(a.obs.BPTA_mean.mean()),OTLS_core=float(a.obs.OTLS_core.mean()) if "OTLS_core" in a.obs else np.nan))
print(pd.DataFrame(rows).to_string(index=False)); print(pd.read_csv(os.path.join(OUT,"xenium_BTAPC_HvP.tsv"),sep="\t").to_string(index=False))
