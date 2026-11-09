#!/usr/bin/env python3
"""Stage 2a: annotate Xenium cells and score OTLS panel / B/T/APC modules per section."""
import os, glob, gzip, time, warnings, re
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, scanpy as sc, anndata as ad, yaml
from concurrent.futures import ProcessPoolExecutor, as_completed
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; XEN=os.path.join(ROOT,"data/raw/xenium"); OUT=os.path.join(ROOT,"results/stage2/xenium/cells")
os.makedirs(OUT,exist_ok=True); os.makedirs(os.path.join(ROOT,"logs/analysis"),exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage2_xenium_cells.log")
MODULES={
 "B":["MS4A1","CD79A","CD19","BANK1"],
 "Plasma":["MZB1","SDC1","PRDM1","TNFRSF17","IGHA1","DERL3"],
 "T":["CD3D","CD3E","CD2","TRAC","CD4","CD8A"],
 "APC":["HLA-DRB5","HLA-DMB","HLA-DQA2","HLA-DQB2","CD14","CD68","CD163","CD1C","FCER1A","LAMP3","CCR7","CD83"],
 "Endo":["PECAM1","VWF","ACKR1","MMRN1"],
 "Fib":["PDGFRA","PDGFRB","COL5A2","COL8A1","POSTN","THY1","DPT","ASPN","VCAN","COL16A1"],
 "Epi":["KRT5","KRT19","KRT7","EPCAM","COL17A1","CLCA2","ODAM"],
 "Neut":["CSF3R","FPR1","MPO","S100A12","CAMP"],
}
BTAPC=["MS4A1","CD79A","CD19","CD3D","CD3E","HLA-DRB5","HLA-DMB","CD14","CD68"]
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
def sig(n): return yaml.safe_load(open(os.path.join(ROOT,"signatures",n+".yaml")))["genes"]
def read_mtx(p):
    with open(p,"rb") as fh: magic=fh.read(2)
    if magic==b"\x1f\x8b":
        with gzip.open(p,"rb") as fh: return scipy.io.mmread(fh).tocsr()
    with open(p,"rb") as fh: return scipy.io.mmread(fh).tocsr()
def process(sid):
    d=os.path.join(XEN,sid)
    feats=pd.read_csv(os.path.join(d,"features.tsv.gz"),sep="\t",header=None)
    mask=feats[2].astype(str).eq("Gene Expression").values
    genes=feats.loc[mask,1].astype(str).values
    bars=pd.read_csv(os.path.join(d,"barcodes.tsv.gz"),header=None)[0].astype(str).values
    X=read_mtx(os.path.join(d,"matrix.mtx.gz"))[mask,:]
    cells=pd.read_csv(os.path.join(d,"cells.csv.gz"),index_col=0).loc[bars]
    total=np.asarray(X.sum(0)).ravel().astype(float); keep=total>=10
    X=X[:,keep]; bars=bars[keep]; cells=cells.loc[bars]; total=total[keep]
    Xn=X.multiply(1.0/np.maximum(total,1)*100).tocsr(); Xn.data=np.log1p(Xn.data)
    a=ad.AnnData(X=Xn.T.tocsr(), obs=pd.DataFrame(index=bars), var=pd.DataFrame(index=genes))
    for name,gs in MODULES.items():
        gs=[g for g in gs if g in a.var_names]
        if len(gs)>=2: sc.tl.score_genes(a,gs,score_name=name,ctrl_size=20,n_bins=10)
    for name,gs in [("OTLS_panel",sig("OTLS_panel")),("OTLS_core_independent",sig("OTLS_core_independent")),("BTAPC",BTAPC),("EpiNeg",MODULES["Epi"])]:
        gs=[g for g in gs if g in a.var_names]
        if len(gs)>=2: sc.tl.score_genes(a,gs,score_name=name,ctrl_size=20,n_bins=10)
    score_cols=[c for c in ["B","Plasma","T","APC","Endo","Fib","Epi","Neut"] if c in a.obs]
    Z=a.obs[score_cols].apply(lambda x:(x-x.mean())/(x.std()+1e-9))
    a.obs["cell_type"]=Z.idxmax(axis=1)
    a.obs.loc[Z.max(axis=1)<0.5,"cell_type"]="Other"
    out=pd.DataFrame({"section":sid,"cell_id":a.obs_names,"x":cells["x_centroid"].values,"y":cells["y_centroid"].values,
                      "total_counts":total,"cell_type":a.obs["cell_type"].values})
    for c in ["B","Plasma","T","APC","Endo","Fib","Epi","Neut","OTLS_panel","OTLS_core_independent","BTAPC","EpiNeg"]:
        if c in a.obs: out[c]=a.obs[c].values
    out.to_csv(os.path.join(OUT,sid+".tsv.gz"),sep="\t",index=False,compression="gzip")
    return sid,len(out),out.cell_type.value_counts().to_dict()
def main():
    sids=sorted(os.path.basename(d) for d in glob.glob(os.path.join(XEN,"*")) if os.path.isdir(d) and os.path.exists(os.path.join(d,"matrix.mtx.gz")))
    log(f"sections={len(sids)}")
    rows=[]
    with ProcessPoolExecutor(max_workers=10) as ex:
        futs={ex.submit(process,s):s for s in sids}
        for i,f in enumerate(as_completed(futs),1):
            try:
                sid,n,ct=f.result(); rows.append(dict(section=sid,n_cells=n)); log(f"{i}/{len(sids)} {sid} n={n} {ct}")
            except Exception as e: log(f"ERROR {futs[f]}: {e}")
    pd.DataFrame(rows).to_csv(os.path.join(OUT,"..","cell_counts.tsv"),sep="\t",index=False)
    log("DONE")
if __name__=="__main__": main()
