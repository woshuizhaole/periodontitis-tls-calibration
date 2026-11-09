#!/usr/bin/env python3
"""Diagnostic: independent full-transcriptome Visium GSE297924 (4 sections, filtered spots)."""
import os, glob, gzip, time, warnings, yaml
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, scanpy as sc, anndata as ad
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; IN=os.path.join(ROOT,"diagnostics/data/tmp/gse297924"); OUT=os.path.join(ROOT,"diagnostics/gse297924"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(OUT,"visium.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
def sig(n): return yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
SIGS=["OTLS_core","OTLS_ext","OTLS_panel","OTLS_core_independent","Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]
EXTRA={"OTLS_stromal":["CXCL13","CCL19","CCL21","LTB","TNFSF13B","CXCL12","PTGDS","PLVAP","MMRN1","SELE","ICAM1","VCAM1","VWF"],
"OTLS_lymphoid":["BCL6","TNFRSF13B","CD27","CXCR5","PDCD1","ICOS","TOX","CD200","TIGIT","LAMP3","CCR7","CD83","CD40","IRF4","CCL20","RBP5","BCL2A1"]}
def mean_z(a,gs,name):
    gs=[g for g in gs if g in a.var_names]
    if len(gs)<2: return 0
    vals=[]
    for g in gs:
        x=np.asarray(a[:,g].X.todense() if sp.issparse(a[:,g].X) else a[:,g].X).ravel().astype(float)
        sd=x.std()
        if sd>0: vals.append((x-x.mean())/sd)
    a.obs[name]=np.mean(vals,axis=0); return len(gs)
rows=[]; allobs=[]
for d in sorted(glob.glob(os.path.join(IN,"*/"))):
    name=os.path.basename(d.rstrip("/")); mtx=os.path.join(d,"matrix.mtx.gz")
    if not os.path.exists(mtx): continue
    with gzip.open(mtx,"rb") as fh: X=scipy.io.mmread(fh).tocsr()
    f=pd.read_csv(os.path.join(d,"features.tsv.gz"),sep="\t",header=None); b=pd.read_csv(os.path.join(d,"barcodes.tsv.gz"),header=None)[0].astype(str).values
    a=ad.AnnData(X=X.T.tocsr().astype(np.float32),obs=pd.DataFrame(index=b),var=pd.DataFrame(index=f[1].astype(str).values))
    a.var_names_make_unique(); sc.pp.filter_cells(a,min_genes=200); sc.pp.filter_genes(a,min_cells=3)
    a.var["mt"]=a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True)
    a=a[a.obs.pct_counts_mt<20].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    for s in SIGS: mean_z(a,sig(s),s)
    for s,gs in EXTRA.items(): mean_z(a,gs,s)
    # B/T/APC modules for TLS-like spot definition (full transcriptome)
    for mod,gs in {"Bmod":["MS4A1","CD79A","CD19","BANK1"],"Tmod":["CD3D","CD3E","CD2","TRAC"],"APCmod":["HLA-DRB5","HLA-DMB","HLA-DQA2","CD14","CD68","CD163","CD1C","FCER1A","LAMP3","CCR7","CD83"]}.items(): mean_z(a,gs,mod)
    z=a.obs[["Bmod","Tmod","APCmod"]].apply(lambda x:(x-x.mean())/(x.std()+1e-9))
    a.obs["BTAPC"]=z.mean(1); thr=a.obs.BTAPC.quantile(0.95); a.obs["TLS_like_spot"]=a.obs.BTAPC>=thr
    a.obs["sample"]=name; allobs.append(a.obs.copy())
    rows.append(dict(sample=name,n_spots=a.n_obs,n_genes=a.n_vars,OTLS_core=float(a.obs.OTLS_core.mean()),OTLS_panel=float(a.obs.OTLS_panel.mean()),
                     Bcell=float(a.obs.Bcell.mean()),Plasma=float(a.obs.Plasma.mean()),Tcell=float(a.obs.Tcell.mean()),
                     TLS_OTLS=float(a.obs.loc[a.obs.TLS_like_spot,"OTLS_core"].mean()),nonTLS_OTLS=float(a.obs.loc[~a.obs.TLS_like_spot,"OTLS_core"].mean())))
    log(f"{name} {a.shape} TLS_like={a.obs.TLS_like_spot.sum()}")
pd.DataFrame(rows).to_csv(os.path.join(OUT,"sample_summary.tsv"),sep="\t",index=False)
pd.concat(allobs).to_csv(os.path.join(OUT,"spot_scores.tsv"),sep="\t")
# correlations / AUC
from sklearn.metrics import roc_auc_score
obs=pd.concat(allobs)
out=[]
for s in ["OTLS_core","OTLS_panel","OTLS_core_independent","OTLS_stromal","OTLS_lymphoid","Bcell","Plasma","Tcell"]:
    if s in obs:
        r=obs[s].corr(obs.BTAPC,method="spearman")
        try: auc=roc_auc_score(obs.TLS_like_spot.astype(int),obs[s])
        except Exception: auc=np.nan
        out.append(dict(score=s,spearman_BTAPC=float(r),AUC_TLSlike=float(auc),TLS_mean=float(obs.loc[obs.TLS_like_spot,s].mean()),nonTLS_mean=float(obs.loc[~obs.TLS_like_spot,s].mean())))
pd.DataFrame(out).to_csv(os.path.join(OUT,"score_vs_BTAPC.tsv"),sep="\t",index=False)
print(pd.DataFrame(rows).to_string(index=False)); print(pd.DataFrame(out).to_string(index=False)); log("DONE")
