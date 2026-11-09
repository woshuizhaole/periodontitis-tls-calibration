#!/usr/bin/env python3
"""Narrative B: B/T/APC co-localization score in tonsil and GSE206621 Visium."""
import os, glob, gzip, warnings, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, scanpy as sc, anndata as ad, yaml
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/narrative_B/spatial"); os.makedirs(OUT,exist_ok=True)
APC=["HLA-DRB5","HLA-DMB","HLA-DQA2","HLA-DQB2","CD14","CD68","CD163","CD1C","FCER1A","LAMP3","CCR7","CD83"]
B=["MS4A1","CD79A","CD19","BANK1"]; T=["CD3D","CD3E","CD2","TRAC"]
def sig(n): return yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
def mean_z(a,genes,name):
    genes=[g for g in genes if g in a.var_names]
    if len(genes)<2: a.obs[name]=np.nan; return 0
    vals=[]
    for g in genes:
        x=np.asarray(a[:,g].X.todense() if sp.issparse(a[:,g].X) else a[:,g].X).ravel().astype(float)
        sd=x.std()
        if sd>0: vals.append((x-x.mean())/sd)
    a.obs[name]=np.mean(vals,axis=0); return len(genes)
# tonsil
a=sc.read_10x_h5(os.path.join(ROOT,"data/raw/tonsil/tonsil_expr.h5")); a.var_names_make_unique()
sc.pp.filter_cells(a,min_genes=200); a.var["mt"]=a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True); a=a[a.obs.pct_counts_mt<20].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
for n,gs in {"Bcell":B,"Tcell":T,"APC":APC,"OTLS_core":sig("OTLS_core"),"OTLS_panel":sig("OTLS_panel")}.items(): mean_z(a,gs,n)
a.obs["BTAPC_min"]=np.minimum.reduce([a.obs.Bcell,a.obs.Tcell,a.obs.APC]); a.obs["BTAPC_mean"]=(a.obs.Bcell+a.obs.Tcell+a.obs.APC)/3
a.obs["follicle_proxy"]=a.obs.Bcell>=a.obs.Bcell.quantile(0.90); a.obs.to_csv(os.path.join(OUT,"tonsil_scores.tsv"),sep="\t")
rows=[]
for s in ["OTLS_core","OTLS_panel","Bcell","Tcell","APC","BTAPC_min","BTAPC_mean"]:
    rows.append(dict(dataset="tonsil",score=s,AUC=float(roc_auc_score(a.obs.follicle_proxy,a.obs[s])),n=len(a),TLS_mean=float(a.obs.loc[a.obs.follicle_proxy,s].mean()),rest=float(a.obs.loc[~a.obs.follicle_proxy,s].mean())))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"tonsil_BTAPC.tsv"),sep="\t",index=False)
# GSE206621
meta=pd.read_csv(os.path.join(ROOT,"metadata/datasets.tsv"),sep="\t"); meta=meta[meta.dataset_id=="GSE206621"][["subject_id","disease_state"]].rename(columns={"subject_id":"sample_id"})
rows=[]
for f in sorted(glob.glob(os.path.join(ROOT,"results/stage2/visium","*_scores.tsv"))):
    sid=os.path.basename(f).replace("_scores.tsv",""); d=pd.read_csv(f,sep="\t",index_col=0)
    # APC not in score table: load raw matrix and compute
    raw=os.path.join(ROOT,"data/external/GSE206621_raw",sid); mtx=glob.glob(os.path.join(raw,"*matrix.mtx*")) or glob.glob(os.path.join(raw,"matrix.mtx*"))
    feat=glob.glob(os.path.join(raw,"*features.tsv*")) or glob.glob(os.path.join(raw,"features.tsv*"))
    if not mtx or not feat: continue
    with (gzip.open if mtx[0].endswith(".gz") else open)(mtx[0],"rb") as fh: X=scipy.io.mmread(fh).tocsr()
    f=pd.read_csv(feat[0],sep="\t",header=None); b=pd.read_csv(glob.glob(os.path.join(raw,"*barcodes.tsv*"))[0] or glob.glob(os.path.join(raw,"barcodes.tsv*"))[0],header=None)[0].astype(str).values
    aa=ad.AnnData(X=X.T.tocsr().astype(np.float32),obs=pd.DataFrame(index=b),var=pd.DataFrame(index=f[1].astype(str).values)); aa.var_names_make_unique()
    sc.pp.filter_cells(aa,min_genes=200); sc.pp.normalize_total(aa,target_sum=1e4); sc.pp.log1p(aa)
    mean_z(aa,APC,"APC"); d=d.join(aa.obs[["APC"]],how="left")
    d["BTAPC_min"]=np.minimum.reduce([d.Bcell,d.Tcell,d.APC]); d["BTAPC_mean"]=(d.Bcell+d.Tcell+d.APC)/3
    d.to_csv(os.path.join(OUT,f"GSE206621_{sid}_BTAPC.tsv"),sep="\t")
    lab=meta[meta.sample_id==sid].disease_state.iloc[0]
    rows.append(dict(sample=sid,label=lab,n=len(d),OTLS_core=float(d.OTLS_core.mean()),BTAPC_min=float(d.BTAPC_min.mean()),BTAPC_mean=float(d.BTAPC_mean.mean()),Bcell=float(d.Bcell.mean()),Tcell=float(d.Tcell.mean()),APC=float(d.APC.mean())))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"GSE206621_summary.tsv"),sep="\t",index=False)
print(pd.DataFrame(rows).to_string(index=False)); print(pd.read_csv(os.path.join(OUT,"tonsil_BTAPC.tsv"),sep="\t").to_string(index=False))
