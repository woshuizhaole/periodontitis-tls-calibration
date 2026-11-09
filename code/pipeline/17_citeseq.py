#!/usr/bin/env python3
"""Stage 1 CITE-seq validation: GSE296447 GEX matrices, 4 periodontitis patients."""
import os, glob, re, gzip, time, warnings
import numpy as np, pandas as pd, scipy.io, scanpy as sc, anndata as ad, yaml
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage1"); RAW=os.path.join(ROOT,"data/processed/GSE296447")
os.makedirs(RAW,exist_ok=True); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage1_citeseq.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
tar=os.path.join(ROOT,"data/raw/geo/GSE296447/GSE296447_RAW.tar")
if not glob.glob(os.path.join(RAW,"*_mRNA_matrix.mtx.gz")):
    log("extracting tar")
    os.system(f"tar -xf {tar} -C {RAW}")
def sig(n): return yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
rows=[]; allobs=[]
for mtx in sorted(glob.glob(os.path.join(RAW,"*_mRNA_matrix.mtx.gz"))):
    pref=mtx.replace("_mRNA_matrix.mtx.gz","")
    gsm=os.path.basename(pref).split("_")[0]
    m=re.search(r"_(P\d+)_raw_mRNA",os.path.basename(mtx)); pid=m.group(1) if m else gsm
    features=pref+"_mRNA_features.tsv.gz"; barcodes=pref+"_mRNA_barcodes.tsv.gz"
    if not (os.path.exists(features) and os.path.exists(barcodes)): continue
    log(f"read {pid} {gsm}")
    with gzip.open(mtx,"rb") as fh: X=scipy.io.mmread(fh).tocsr().T
    f=pd.read_csv(features,sep="\t",header=None); b=pd.read_csv(barcodes,header=None)[0].astype(str).values
    a=ad.AnnData(X=X.astype(np.float32), obs=pd.DataFrame(index=b), var=pd.DataFrame(index=f[1].astype(str).values))
    a.var_names_make_unique(); sc.pp.filter_cells(a,min_genes=200); sc.pp.filter_genes(a,min_cells=3)
    a.var["mt"]=a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True)
    a=a[(a.obs.pct_counts_mt<20)&(a.obs.n_genes_by_counts<6000)].copy()
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    for n in ["OTLS_core","OTLS_ext","OTLS_panel","Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]:
        gs=[g for g in sig(n) if g in a.var_names]
        if len(gs)>=2: sc.tl.score_genes(a,gs,score_name=n,ctrl_size=50,n_bins=25)
    a.obs["patient"]=pid; a.obs["gsm"]=gsm; allobs.append(a.obs.copy())
    rows.append(dict(patient=pid,gsm=gsm,n_cells=a.n_obs,n_genes=a.n_vars,
                     OTLS_core=float(a.obs.OTLS_core.mean()),OTLS_panel=float(a.obs.OTLS_panel.mean()),
                     Bcell=float(a.obs.Bcell.mean()),Plasma=float(a.obs.Plasma.mean()),Tcell=float(a.obs.Tcell.mean())))
    log(f"done {pid} {a.shape}")
pd.DataFrame(rows).to_csv(os.path.join(OUT,"GSE296447_patient_summary.tsv"),sep="\t",index=False)
if allobs: pd.concat(allobs).to_csv(os.path.join(OUT,"GSE296447_cell_scores.tsv"),sep="\t")
log("DONE")
