#!/usr/bin/env python3
"""Stage 0 Xenium audit: cell counts, negative controls, gene coverage, quick B/T/APC aggregate scan."""
import os, re, json, math, time
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import gzip
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT="/data/workplace/wyx/perio_tls_20261107"
XEN=os.path.join(ROOT,"data/raw/xenium")
OUT=os.path.join(ROOT,"metadata/xenium_section_audit.tsv")
LOG=os.path.join(ROOT,"logs/preprocessing/xenium_audit.log")

CANDIDATES = """MS4A1 CD79A CD79B CD19 BANK1 CD74 HLA-DRA MZB1 JCHAIN SDC1 PRDM1 XBP1 TNFRSF17 IGHG1 IGHA1
CXCR5 PDCD1 ICOS BCL6 TOX TCF7 IL21 CD40LG CD200 CXCL13 CR2 FDCSP LTB TNFSF13B CCL19 CCL21
PLVAP MMRN1 ACKR1 VWF SELE ICAM1 VCAM1 CD34 LAMP3 CCR7 CD83 CD40 IRF4 CLEC9A
COL1A1 PDGFRA FAP THY1 CXCL12 IL33 AICDA MEF2B MME BCL2A1 CD27 TNFRSF13B
HLA-DRB5 HLA-DMB HLA-DQA2 HLA-DQB2 CD14 CD68 CD163 CD1C FCER1A CD3D CD3E CD2 TRAC CD4 CD8A
KRT5 KRT19 KRT7 EPCAM COL17A1 CLCA2 ODAM CSF3R FPR1 MPO S100A12 CAMP""".split()

B_MARK=["MS4A1","CD79A","CD19","BANK1"]
T_MARK=["CD3D","CD3E","CD2","TRAC"]
APC_MARK=["HLA-DRB5","HLA-DMB","HLA-DQA2","HLA-DQB2","CD14","CD68","CD163","CD1C","FCER1A"]
EPI_MARK=["KRT5","KRT19","KRT7","EPCAM","COL17A1","CLCA2","ODAM"]
PIXEL_UM=0.2125
RADIUS_PX=50/PIXEL_UM


def log(msg):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"; print(line,flush=True)
    open(LOG,"a").write(line+"\n")


def process_section(sid):
    d=os.path.join(XEN,sid)
    try:
        feats=pd.read_csv(os.path.join(d,"features.tsv.gz"),sep="\t",header=None)
        mask=feats[2].astype(str).eq("Gene Expression").values
        genes=feats.loc[mask,1].astype(str).values
        bars=pd.read_csv(os.path.join(d,"barcodes.tsv.gz"),header=None)[0].astype(str).values
        mtx=os.path.join(d,"matrix.mtx.gz")
        with open(mtx,"rb") as _fh: magic=_fh.read(2)
        if magic==b"\x1f\x8b":
            with gzip.open(mtx,"rb") as _fh: X=scipy.io.mmread(_fh).tocsr()
        else:
            with open(mtx,"rb") as _fh: X=scipy.io.mmread(_fh).tocsr()
        # matrix is genes x cells; subset GE
        X=X[mask,:]
        cells=pd.read_csv(os.path.join(d,"cells.csv.gz"),index_col=0)
        cells=cells.loc[bars]
        total=np.asarray(X.sum(0)).ravel().astype(float)
        n_deposited=int(X.shape[1])
        neg_all=(cells.get("control_probe_counts",0).values + cells.get("control_codeword_counts",0).values)
        neg_section=float(neg_all.sum()/max(cells["total_counts"].values.sum(),1))
        keep=total>=10
        X=X[:,keep]; bars=bars[keep]; cells=cells.loc[bars]; total=total[keep]
        # negative control fraction from cells.csv
        neg=(cells.get("control_probe_counts",0).values + cells.get("control_codeword_counts",0).values)
        neg_frac=np.divide(neg, np.maximum(cells["total_counts"].values,1))
        # normalize
        Xn=X.multiply(1.0/np.maximum(total,1)*100).tocsr()
        Xn.data=np.log1p(Xn.data)
        def score(markers):
            idx=[i for i,g in enumerate(genes) if g in markers]
            if not idx: return np.zeros(Xn.shape[1])
            sub=Xn[idx,:]
            return np.asarray(sub.mean(0)).ravel()
        b=score(B_MARK); t=score(T_MARK); a=score(APC_MARK); e=score(EPI_MARK)
        # local fractions via KDTree
        xy=cells[["x_centroid","y_centroid"]].values
        tree=cKDTree(xy)
        neigh=tree.query_ball_point(xy, RADIUS_PX)
        b_hi=b>np.quantile(b,0.90); t_hi=t>np.quantile(t,0.90); a_hi=a>np.quantile(a,0.90)
        nbscore=np.zeros(len(xy))
        for i,idx in enumerate(neigh):
            if not idx: continue
            nbscore[i]=(b_hi[idx].mean()+t_hi[idx].mean()+a_hi[idx].mean())/3
        # high neighborhoods, cluster
        thr=np.quantile(nbscore,0.95)
        high=np.where(nbscore>=thr)[0]
        regions=0; region_sizes=[]
        null_regions=[]
        def count_regions(mask):
            idx=np.where(mask)[0]
            if len(idx)<20: return 0
            lab=DBSCAN(eps=RADIUS_PX,min_samples=20).fit_predict(xy[idx])
            return sum(1 for k in set(lab) if k!=-1 and int((lab==k).sum())>=20)
        regions=count_regions(nbscore>=thr)
        if regions:
            region_sizes=[int(nbscore[mask].size) for mask in [nbscore>=thr]]
        rng=np.random.default_rng(42)
        for _ in range(3):
            b2=rng.permutation(b_hi); t2=rng.permutation(t_hi); a2=rng.permutation(a_hi)
            nb2=np.zeros(len(xy))
            for i,idx in enumerate(neigh):
                if idx: nb2[i]=(b2[idx].mean()+t2[idx].mean()+a2[idx].mean())/3
            null_regions.append(count_regions(nb2>=np.quantile(nb2,0.95)))
        # gene coverage
        cov=sum(1 for g in CANDIDATES if g in set(genes))
        return dict(section=sid, group="P" if sid.startswith("P") else "H",
                    n_cells=int(len(bars)), n_cells_deposited=int(n_deposited), median_total=float(np.median(cells["total_counts"].values)),
                    neg_frac_median=float(np.nanmedian(neg_frac)), neg_frac_section=neg_section, neg_frac_max_cell=float(np.nanmax(neg_frac)), neg_frac_gt5=int((neg_frac>0.05).sum()),
                    n_ge_genes=int(len(genes)), candidate_genes_present=int(cov),
                    btapc_rich_regions=int(regions), btapc_rich_cells=int(sum(region_sizes)),
                    null_btapc_rich_regions_mean=float(np.mean(null_regions)) if null_regions else 0.0,
                    analysis_note="relative per-section thresholds; B/T/APC-rich only, not TLS evidence",
                    radius_px=float(RADIUS_PX))
    except Exception as e:
        return dict(section=sid, group=("P" if sid.startswith("P") else "H"), error=str(e))


def main():
    sections=sorted([d for d in os.listdir(XEN) if os.path.isdir(os.path.join(XEN,d))])
    log(f"sections={len(sections)}")
    rows=[]
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs={ex.submit(process_section,s):s for s in sections}
        for i,fut in enumerate(as_completed(futs),1):
            r=fut.result(); rows.append(r); log(f"{i}/{len(sections)} {r.get('section')} cells={r.get('n_cells')} err={r.get('error','')}")
    df=pd.DataFrame(rows)
    df.to_csv(OUT,sep="\t",index=False)
    log(f"wrote {OUT} rows={len(df)}")

if __name__=="__main__":
    main()
