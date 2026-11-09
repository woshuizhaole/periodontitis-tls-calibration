#!/usr/bin/env python3
"""P1.1: coarse pixel-level B/T/APC context in public 40-plex PhenoCycler (4 PD/4 health)."""
import os, glob, time, warnings
import numpy as np, pandas as pd, tifffile
from scipy import ndimage
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; IN=os.path.join(ROOT,"data/external/phenocycler"); OUT=os.path.join(ROOT,"results/narrative_B/spatial/phenocycler"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(OUT,"phenocycler.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
CH=["DAPI","Keratin 14","CD31","CD4","AQP3","HLA-A","PD-L1","CD34","CD20","CD40","Podoplanin","SMA","CD68","CD45RO","Vimentin","CD45","IFNG","IDO1","CD8","FOXP3","Mac2/Galectin-3","HLA-DR","CD14","CD3e","Granzyme B","Collagen IV","Ki67","PD-1","Tryptase","Caveolin","ICOS","CD141","Bcl-2","CD38","ST2/IL1RL1","MPO","EpCAM","CD39","Pan-Cytokeratin"]
B=["CD20"]; T=["CD3e","CD4","CD8"]; APC=["HLA-DR","CD14","CD68","CD141"]; HEV=["CD31","Podoplanin","Collagen IV"]; PRO=["Ki67"]
def block8(a):
    h,w=a.shape; h-=h%8; w-=w%8; return a[:h,:w].reshape(h//8,8,w//8,8).mean((1,3)).astype(np.float32)
def pos(vals):
    vals=vals[np.isfinite(vals)]
    if len(vals)==0: return None
    q=np.quantile(vals,[0.5,0.99]); return float(q[1]*0.5) if q[1]>q[0] else np.nan
rows=[]
for p in sorted(glob.glob(os.path.join(IN,"*.ome.tif"))):
    name=os.path.basename(p).replace(".ome.tif","").replace(" ","_"); grp="PD" if name.startswith("Grade") else "H"
    log(f"read {name}")
    arrs=[]
    with tifffile.TiffFile(p) as tf:
        n=len(tf.pages)
        for i in range(min(n,len(CH))):
            arrs.append(block8(tf.pages[i].asarray()))
    A=np.stack(arrs,0);  # C,H,W
    dapi=A[CH.index("DAPI")]; tissue=dapi>np.quantile(dapi,0.35)
    # threshold each marker at 50% of its 99th percentile in tissue
    masks={}
    for m in set(B+T+APC+HEV+PRO):
        v=A[CH.index(m)]; thr=pos(v[tissue])
        masks[m]=(v>thr)&tissue if thr and np.isfinite(thr) else np.zeros_like(tissue,bool)
    def frac(mask): return float(mask.sum()/max(tissue.sum(),1))
    def combo(ms): return np.mean([masks[m] for m in ms],0)
    bt=combo(B); tt=combo(T); at=combo(APC); ht=combo(HEV)
    # smooth and co-localization: all three smoothed masks positive
    bs=ndimage.gaussian_filter(bt.astype(float),2); ts=ndimage.gaussian_filter(tt.astype(float),2); as_=ndimage.gaussian_filter(at.astype(float),2)
    coloc=(bs>0.15)&(ts>0.15)&(as_>0.15)&tissue
    # B aggregates: connected components of B mask after closing
    bm=ndimage.binary_closing(masks["CD20"],iterations=2); lab,nlab=ndimage.label(bm)
    sizes=np.bincount(lab.ravel())[1:]; big=sizes[sizes>=50]
    rows.append(dict(sample=name,group=grp,n_channels=n,shape=f"{A.shape[1]}x{A.shape[2]}",tissue_frac=float(tissue.mean()),
                     B_area=frac(masks["CD20"]),T_area=float(tt.mean()),APC_area=float(at.mean()),HEV_area=float(ht.mean()),Ki67_area=frac(masks["Ki67"]),
                     BTAPC_coloc=float(coloc.mean()),T_B_ratio=float(tt.mean()/(bt.mean()+1e-9)),
                     B_aggregates=int(len(big)),B_aggregate_area=float(big.sum()/max(tissue.sum(),1)),
                     B_max_agg=float(big.max()/max(tissue.sum(),1)) if len(big) else 0.0))
    log(f"done {name} tissue={tissue.mean():.3f} B={rows[-1]['B_area']:.4f} T={tt.mean():.4f} APC={at.mean():.4f} coloc={coloc.mean():.5f}")
df=pd.DataFrame(rows); df.to_csv(os.path.join(OUT,"phenocycler_metrics.tsv"),sep="\t",index=False)
from scipy.stats import mannwhitneyu
res=[]
for c in ["B_area","T_area","APC_area","HEV_area","Ki67_area","BTAPC_coloc","T_B_ratio","B_aggregates","B_aggregate_area","B_max_agg"]:
    h=df[df.group=="H"][c]; p=df[df.group=="PD"][c]
    res.append(dict(metric=c,H_mean=float(h.mean()),PD_mean=float(p.mean()),H_median=float(h.median()),PD_median=float(p.median()),p=float(mannwhitneyu(h,p).pvalue),nH=len(h),nPD=len(p)))
pd.DataFrame(res).to_csv(os.path.join(OUT,"phenocycler_group_tests.tsv"),sep="\t",index=False)
print(df.to_string(index=False)); print(pd.DataFrame(res).to_string(index=False)); log("DONE")
