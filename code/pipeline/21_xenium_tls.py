#!/usr/bin/env python3
"""Stage 2b: define B/T/APC-rich TLS-like regions, density-match controls, test OTLS_panel."""
import os, glob, time, warnings
import numpy as np, pandas as pd
from scipy.spatial import cKDTree, ConvexHull
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from scipy.stats import wilcoxon, mannwhitneyu
import statsmodels.api as sm
import statsmodels.formula.api as smf
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; IN=os.path.join(ROOT,"results/stage2/xenium/cells"); OUT=os.path.join(ROOT,"results/stage2/xenium")
os.makedirs(OUT,exist_ok=True); LOG=os.path.join(ROOT,"logs/analysis/stage2_xenium_tls.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
PIX=0.2125; RAD=50/PIX
# load cells
files=sorted(glob.glob(os.path.join(IN,"*.tsv.gz")))
cells=pd.concat([pd.read_csv(f,sep="\t") for f in files],ignore_index=True)
log(f"cells={cells.shape}")
dm=pd.read_csv(os.path.join(ROOT,"metadata/xenium_donor_map.tsv"),sep="\t")
cells=cells.merge(dm[["section","donor","group","main_quant"]],on="section",how="left")
# neighborhood fractions
frac=[]
for sid,g in cells.groupby("section"):
    xy=g[["x","y"]].values; tree=cKDTree(xy); idx=tree.query_ball_point(xy,RAD,workers=8)
    d={}
    for ct in ["B","T","APC","Plasma","Epi","Endo","Fib"]:
        flag=(g.cell_type.values==ct).astype(float)
        d["nb_"+ct]=np.array([flag[ii].mean() if len(ii) else 0 for ii in idx])
    f=pd.DataFrame(d,index=g.index); f["section"]=sid; frac.append(f)
frac=pd.concat(frac); cells=cells.join(frac.drop(columns="section"))
for c in ["nb_B","nb_T","nb_APC","nb_Epi"]:
    cells["z_"+c]=(cells[c]-cells[c].mean())/(cells[c].std()+1e-9)
cells["tls_nb_score"]=cells.z_nb_B+cells.z_nb_T+cells.z_nb_APC-0.5*cells.z_nb_Epi
# stringent, co-presence definition: top 5% neighborhood score AND B/T/APC each >=3% of neighbors
thr=float(np.quantile(cells.tls_nb_score,0.95))
cells["tls_like"]=(cells.tls_nb_score>=thr)&(cells.nb_B>=0.03)&(cells.nb_T>=0.03)&(cells.nb_APC>=0.03)
log(f"tls_nb_thr={thr:.3f} tls_like cells={cells.tls_like.sum()} / {len(cells)}")
# regions per section
regions=[]; controls=[]; region_cells={}
cells["region_id"]=""
rng=np.random.default_rng(42)
for sid,g in cells.groupby("section"):
    hi=g[g.tls_like]
    if len(hi)<20: continue
    lab=DBSCAN(eps=RAD,min_samples=20).fit_predict(hi[["x","y"]].values)
    for k in sorted(set(lab)):
        if k==-1: continue
        idx=hi.index[lab==k]
        if len(idx)<20: continue
        r=g.loc[idx]; xy=r[["x","y"]].values
        region_cells[(sid,k)]=idx
        cells.loc[idx,"region_id"]=f"{sid}_{k}"
        try:
            area=ConvexHull(xy).volume*PIX*PIX
        except Exception:
            area=len(r)*np.pi*(50**2)
        regions.append(dict(section=sid,region=k,n_cells=len(r),area_um2=float(area),
                            B=float(r.nb_B.mean()),T=float(r.nb_T.mean()),APC=float(r.nb_APC.mean()),Plasma=float(r.nb_Plasma.mean()),
                            OTLS_panel=float(r.OTLS_panel.mean()),BTAPC=float(r.BTAPC.mean()),
                            OTLS_ind=float(r.OTLS_core_independent.mean())))
        # controls: random non-TLS centers, n_cells nearest
        pool=g[~g.tls_like]
        if len(pool)<len(r): continue
        centers=pool.sample(n=min(5,len(pool)),random_state=int(rng.integers(1e9)))
        tree=cKDTree(g[["x","y"]].values)
        for _,c in centers.iterrows():
            d,ii=tree.query([c.x,c.y],k=min(len(r),len(g)))
            rr=g.iloc[ii]
            controls.append(dict(section=sid,region=k,n_cells=len(rr),area_um2=float(area),
                                 B=float(rr.nb_B.mean()),T=float(rr.nb_T.mean()),APC=float(rr.nb_APC.mean()),Plasma=float(rr.nb_Plasma.mean()),
                                 OTLS_panel=float(rr.OTLS_panel.mean()),BTAPC=float(rr.BTAPC.mean()),
                                 OTLS_ind=float(rr.OTLS_core_independent.mean()),control_id=f"{sid}_{k}_{c.name}",overlap_tls=float(rr.tls_like.mean())))
reg=pd.DataFrame(regions); ctl=pd.DataFrame(controls)
reg.to_csv(os.path.join(OUT,"tls_regions.tsv"),sep="\t",index=False)
ctl.to_csv(os.path.join(OUT,"control_regions.tsv"),sep="\t",index=False)
log(f"regions={len(reg)} controls={len(ctl)}")
# density matching TLS -> control
pairs=[]
if len(reg) and len(ctl):
    X=reg[["B","T","APC"]].values; Y=ctl[["B","T","APC"]].values
    nn=NearestNeighbors(n_neighbors=1).fit(Y); dist,idx=nn.kneighbors(X)
    for i,r in reg.iterrows():
        j=idx[i,0]; c=ctl.iloc[j]
        pairs.append(dict(section=r.section,region=r.region,control=c.control_id,density_dist=float(dist[i,0]),
                          OTLS_tls=r.OTLS_panel,OTLS_ctl=c.OTLS_panel,BTAPC_tls=r.BTAPC,BTAPC_ctl=c.BTAPC,
                          OTLSind_tls=r.OTLS_ind,OTLSind_ctl=c.OTLS_ind))
pairs=pd.DataFrame(pairs); pairs.to_csv(os.path.join(OUT,"matched_pairs.tsv"),sep="\t",index=False)
# tests
res=[]
if len(pairs):
    for name,a,b in [("OTLS_panel","OTLS_tls","OTLS_ctl"),("BTAPC","BTAPC_tls","BTAPC_ctl"),("OTLS_independent","OTLSind_tls","OTLSind_ctl")]:
        d=pairs[a]-pairs[b]; p=wilcoxon(pairs[a],pairs[b]).pvalue if len(pairs)>=5 else np.nan
        res.append(dict(test=f"matched_{name}",n=len(pairs),mean_diff=float(d.mean()),median_diff=float(d.median()),
                        cohen_d=float(d.mean()/(d.std()+1e-12)),p=float(p)))
# distance decay: each region uses its own cells (not all TLS-like cells in the section)
dec=[]; slopes_by_region={}
for (sid,rid),idx in region_cells.items():
    sub=cells[cells.section==sid]; regcells=cells.loc[idx]
    if len(regcells)<5: continue
    tree=cKDTree(regcells[["x","y"]].values); d,_=tree.query(sub[["x","y"]].values)
    m=(d>0)&(d<=200/PIX)
    if m.sum()<20: continue
    dd=d[m]*PIX; y=sub.OTLS_panel.values[m]; yb=sub.BTAPC.values[m]
    bins=[0,25,50,100,150,200]
    for lo,hi in zip(bins[:-1],bins[1:]):
        mm=(dd>=lo)&(dd<hi)
        if mm.sum()>=5:
            dec.append(dict(section=sid,region=rid,bin=f"{lo}-{hi}",n=int(mm.sum()),
                            OTLS=float(np.mean(y[mm])),BTAPC=float(np.mean(yb[mm]))))
dec=pd.DataFrame(dec); dec.to_csv(os.path.join(OUT,"distance_decay.tsv"),sep="\t",index=False)
if len(dec):
    for name,col in [("OTLS","OTLS"),("BTAPC","BTAPC")]:
        slopes=[]
        for (sid,rid),g in dec.groupby(["section","region"]):
            mid=np.array([np.mean(list(map(float,b.split("-")))) for b in g.bin])
            sl=float(np.polyfit(mid,g[col],1)[0]); slopes.append(sl); slopes_by_region.setdefault((sid,rid),{})[name]=sl
        res.append(dict(test=f"decay_slope_{name}",n=len(slopes),mean_diff=float(np.mean(slopes)),median_diff=float(np.median(slopes)),cohen_d=np.nan,p=np.nan))
    both=pd.DataFrame([dict(section=k[0],region=k[1],**v) for k,v in slopes_by_region.items() if "OTLS" in v and "BTAPC" in v])
    if len(both)>=5:
        p=wilcoxon(both.OTLS,both.BTAPC).pvalue
        res.append(dict(test="decay_slope_paired",n=len(both),mean_diff=float((both.OTLS-both.BTAPC).mean()),median_diff=float((both.OTLS-both.BTAPC).median()),cohen_d=float((both.OTLS-both.BTAPC).mean()/((both.OTLS-both.BTAPC).std()+1e-12)),p=float(p)))
# H vs P region counts / area / OTLS
section_tab=[]
for sid,g in cells.groupby("section"):
    rr=reg[reg.section==sid]
    # ROI area from hull of all cells
    try: roi=ConvexHull(g[["x","y"]].values).volume*PIX*PIX
    except Exception: roi=np.nan
    section_tab.append(dict(section=sid,group=g.group.iloc[0],donor=g.donor.iloc[0],main_quant=bool(g.main_quant.iloc[0]),
                            n_regions=len(rr),roi_um2=roi,area_frac=float(rr.area_um2.sum()/roi) if roi and len(rr) else 0,
                            OTLS=float(np.average(rr.OTLS_panel,weights=rr.area_um2)) if len(rr) else np.nan,
                            BTAPC=float(np.average(rr.BTAPC,weights=rr.area_um2)) if len(rr) else np.nan))
st=pd.DataFrame(section_tab); st.to_csv(os.path.join(OUT,"section_summary.tsv"),sep="\t",index=False)
# group tests on main_quant
for subset,name in [(st[st.main_quant],"main_quant"),(st,"all")]:
    if subset.group.nunique()==2 and len(subset)>=10:
        sub=subset.dropna(subset=["roi_um2"]).copy(); sub["log_area"]=np.log(sub.roi_um2/1e6)
        try:
            m=smf.glm("n_regions ~ group",data=sub,family=sm.families.Poisson(),offset=sub.log_area).fit()
            res.append(dict(test=f"region_count_{name}",n=len(sub),mean_diff=float(np.exp(m.params.get("group[T.P]",np.nan))),median_diff=np.nan,cohen_d=np.nan,p=float(m.pvalues.get("group[T.P]",np.nan))))
        except Exception as e: log(f"count model {name} failed {e}")
        try:
            x=sub.dropna(subset=["OTLS"])
            d=x[x.group=="P"].OTLS.mean()-x[x.group=="H"].OTLS.mean()
            p=mannwhitneyu(x[x.group=="P"].OTLS,x[x.group=="H"].OTLS).pvalue
            res.append(dict(test=f"OTLS_HvP_{name}",n=len(x),mean_diff=float(d),median_diff=np.nan,cohen_d=float(d/(x.OTLS.std()+1e-12)),p=float(p)))
        except Exception: pass
cells[["section","cell_id","x","y","tls_like","region_id"]].to_csv(os.path.join(OUT,"cells_tls.tsv.gz"),sep="\t",index=False,compression="gzip")
pd.DataFrame(res).to_csv(os.path.join(OUT,"tests.tsv"),sep="\t",index=False)
log("DONE")
