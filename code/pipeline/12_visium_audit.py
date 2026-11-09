#!/usr/bin/env python3
import os, gzip, glob, pandas as pd, numpy as np, scipy.io
ROOT="/data/workplace/wyx/perio_tls_20261107"; RAW=os.path.join(ROOT,"data/external/GSE206621_raw"); OUT=os.path.join(ROOT,"metadata/visium_sample_audit.tsv")
rows=[]
for d in sorted(glob.glob(os.path.join(RAW,"*"))):
    if not os.path.isdir(d): continue
    sid=os.path.basename(d)
    pos=glob.glob(os.path.join(d,"*tissue_positions_list.csv*"))
    if not pos:
        # maybe spatial subdir
        pos=glob.glob(os.path.join(d,"spatial","*tissue_positions_list.csv*"))
    n_int=0; pos_rows=0
    if pos:
        p=pos[0]
        op=gzip.open if p.endswith(".gz") else open
        with op(p,"rt") as fh:
            for line in fh:
                parts=line.rstrip().split(",")
                pos_rows+=1
                if len(parts)>1 and parts[1].strip() in ("1","True","TRUE"): n_int+=1
    feats=glob.glob(os.path.join(d,"*features.tsv*"))
    bars=glob.glob(os.path.join(d,"*barcodes.tsv*"))
    mats=glob.glob(os.path.join(d,"*matrix.mtx*"))
    n_barcodes=0; genes=set(); total_med=np.nan; ms4a1=cd3d=np.nan
    if bars:
        n_barcodes=sum(1 for _ in (gzip.open if bars[0].endswith(".gz") else open)(bars[0],"rt"))
    if feats:
        f=pd.read_csv(feats[0],sep="\t",header=None)
        genes=set(f[1].astype(str))
    if mats and bars and feats:
        op=gzip.open if mats[0].endswith(".gz") else open
        with op(mats[0],"rb") as fh:
            X=scipy.io.mmread(fh).tocsr()
        tot=np.asarray(X.sum(0)).ravel()
        total_med=float(np.median(tot)) if len(tot) else np.nan
        for g in ["MS4A1","CD3D","CD79A","CD3E"]:
            if g in genes:
                i=list(f[1].astype(str)).index(g)
                frac=float((X.getrow(i).toarray().ravel()>0).mean()) if X.shape[0]>i else np.nan
                if g=="MS4A1": ms4a1=frac
                if g=="CD3D": cd3d=frac
    rows.append(dict(sample=sid, n_positions_rows=pos_rows, n_in_tissue=n_int, n_barcodes=n_barcodes,
                     genes=len(genes), median_counts=total_med, MS4A1_detect_frac=ms4a1, CD3D_detect_frac=cd3d,
                     low_quality=(n_int<500)))
df=pd.DataFrame(rows); df.to_csv(OUT,sep="\t",index=False); print(df.to_string(index=False))
