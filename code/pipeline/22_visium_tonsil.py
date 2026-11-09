#!/usr/bin/env python3
"""Stage 2c: Visium attenuation cohort + tonsil positive control."""
import os, glob, gzip, time, warnings
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, scanpy as sc, anndata as ad, yaml
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage2"); os.makedirs(OUT+"/visium",exist_ok=True); os.makedirs(OUT+"/tonsil",exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage2_visium_tonsil.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
def sig(n): return yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
def mean_z(adata,genes,name):
    genes=[g for g in genes if g in adata.var_names]
    if len(genes)<2: return
    vals=[]
    for g in genes:
        x=np.asarray(adata[:,g].X.todense() if sp.issparse(adata[:,g].X) else adata[:,g].X).ravel().astype(float)
        sd=x.std()
        if sd>0: vals.append((x-x.mean())/sd)
    adata.obs[name]=np.mean(vals,axis=0)
# ---- Visium ----
vis_rows=[]
for d in sorted(glob.glob(os.path.join(ROOT,"data/external/GSE206621_raw","*"))):
    if not os.path.isdir(d): continue
    sid=os.path.basename(d)
    mtx=glob.glob(os.path.join(d,"*matrix.mtx*")); feats=glob.glob(os.path.join(d,"*features.tsv*")); bars=glob.glob(os.path.join(d,"*barcodes.tsv*")); pos=glob.glob(os.path.join(d,"*tissue_positions_list.csv*")) or glob.glob(os.path.join(d,"spatial","*tissue_positions_list.csv*"))
    if not (mtx and feats and bars and pos): log(f"skip {sid}"); continue
    with (gzip.open if mtx[0].endswith(".gz") else open)(mtx[0],"rb") as fh: X=scipy.io.mmread(fh).tocsr()
    f=pd.read_csv(feats[0],sep="\t",header=None); b=pd.read_csv(bars[0],header=None)[0].astype(str).values
    p=pd.read_csv(pos[0],header=None,names=["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"])
    a=ad.AnnData(X=X.T.tocsr().astype(np.float32),obs=pd.DataFrame(index=b),var=pd.DataFrame(index=f[1].astype(str).values))
    a.var_names_make_unique(); p=p.set_index("barcode").reindex(a.obs_names)
    a.obs["in_tissue"]=p["in_tissue"].values; a=a[a.obs.in_tissue==1].copy()
    a.obsm["spatial"]=p.loc[a.obs_names,["pxl_col","pxl_row"]].values
    a.var["mt"]=a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True)
    sc.pp.filter_cells(a,min_genes=200); a=a[a.obs.pct_counts_mt<20].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    for n in ["OTLS_core","OTLS_panel","OTLS_core_independent","Bcell","Tcell","Plasma","12CK","TLS39"]:
        mean_z(a,sig(n),n)
    a.obs["sample"]=sid; a.obs.to_csv(os.path.join(OUT,"visium",sid+"_scores.tsv"),sep="\t")
    vis_rows.append(dict(sample=sid,n_spots=a.n_obs,OTLS_core=float(a.obs.OTLS_core.mean()),OTLS_panel=float(a.obs.OTLS_panel.mean()),
                         Bcell=float(a.obs.Bcell.mean()),Tcell=float(a.obs.Tcell.mean()),Plasma=float(a.obs.Plasma.mean())))
    log(f"visium {sid} {a.shape}")
pd.DataFrame(vis_rows).to_csv(os.path.join(OUT,"visium","sample_summary.tsv"),sep="\t",index=False)
# ---- tonsil ----
log("tonsil load")
a=sc.read_10x_h5(os.path.join(ROOT,"data/raw/tonsil/tonsil_expr.h5"))
a.var_names_make_unique()
pos=pd.read_csv(os.path.join(ROOT,"data/raw/tonsil/spatial/tissue_positions.csv"))
pos.columns=[c.lower() for c in pos.columns]
bcol="barcode" if "barcode" in pos.columns else pos.columns[0]
pos=pos.set_index(bcol)
a.obs=a.obs.join(pos)
if "in_tissue" in a.obs: a=a[a.obs.in_tissue.astype(str).isin(["1","True","TRUE"])].copy()
a.obsm["spatial"]=a.obs[["pxl_col_in_fullres","pxl_row_in_fullres"]].values if "pxl_col_in_fullres" in a.obs else a.obsm.get("spatial")
a.var["mt"]=a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True)
sc.pp.filter_cells(a,min_genes=200); a=a[a.obs.pct_counts_mt<20].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
for n in ["OTLS_core","OTLS_panel","OTLS_core_independent","Bcell","Tcell","Plasma","12CK","TLS39"]:
    mean_z(a,sig(n),n)
# follicle proxy: top 10% Bcell
thr=a.obs.Bcell.quantile(0.90); a.obs["follicle_proxy"]=a.obs.Bcell>=thr
a.obs.to_csv(os.path.join(OUT,"tonsil","tonsil_scores.tsv"),sep="\t")
from sklearn.metrics import roc_auc_score
auc=roc_auc_score(a.obs.follicle_proxy,a.obs.OTLS_core)
auc_p=roc_auc_score(a.obs.follicle_proxy,a.obs.OTLS_panel)
pd.DataFrame([dict(n_spots=a.n_obs,follicle_frac=float(a.obs.follicle_proxy.mean()),AUC_OTLS_core=float(auc),AUC_OTLS_panel=float(auc_p),
                   OTLS_core_follicle=float(a.obs[a.obs.follicle_proxy].OTLS_core.mean()),OTLS_core_rest=float(a.obs[~a.obs.follicle_proxy].OTLS_core.mean()))]).to_csv(os.path.join(OUT,"tonsil","positive_control.tsv"),sep="\t",index=False)
log(f"tonsil {a.shape} auc={auc:.3f}/{auc_p:.3f}")
log("DONE")
