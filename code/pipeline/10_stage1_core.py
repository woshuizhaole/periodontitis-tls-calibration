#!/usr/bin/env python3
"""Stage 1: annotate GSE164241, define TLS-like states, derive OTLS core/ext/panel signatures."""
import os, re, json, time, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import yaml
import scipy.sparse as sp

warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"
OUT=os.path.join(ROOT,"results/stage1"); os.makedirs(OUT,exist_ok=True)
os.makedirs(os.path.join(OUT,"de_results"),exist_ok=True)
os.makedirs(os.path.join(ROOT,"signatures"),exist_ok=True)
os.makedirs(os.path.join(ROOT,"logs/analysis"),exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage1_core.log")

def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")

# ---------- load ----------
h5=os.path.join(ROOT,"data/external/GSE164241_integrated.h5ad")
log(f"loading {h5}")
a=ad.read_h5ad(h5)
log(f"loaded {a.shape}")
# use log-normalized X; counts retained in layers
cl=a.obs["leiden_r1.0"].astype(str)
cluster_map={
 "0":"Fibroblast","1":"Fibroblast","2":"Fibroblast","5":"Fibroblast","10":"Fibroblast","18":"Fibroblast",
 "6":"Endothelial","7":"Endothelial","4":"Epithelial","14":"B","11":"Plasma","22":"Plasma",
 "15":"T","16":"T","17":"NK","13":"Myeloid","20":"Mast","23":"Neutrophil",
 "3":"Mixed","8":"Mixed","9":"Mixed","12":"Mixed","19":"Mixed","21":"Mixed"}
a.obs["cell_type"]=cl.map(cluster_map).fillna("Mixed").astype(str)
log("cell types: "+str(a.obs["cell_type"].value_counts().to_dict()))

# ---------- state scores ----------
score_sets={
 "B_TLS":["BCL6","AICDA","MEF2B","MME","CD38","TNFRSF13B","CD27"],
 "T_TLS":["CXCR5","PDCD1","ICOS","IL21","BCL6","TOX","CD200"],
 "Fib_TLS":["CXCL13","CCL19","CCL21","LTB","TNFSF13B"],
 "Endo_TLS":["PLVAP","MMRN1","ACKR1","SELE","ICAM1"],
 "Myeloid_TLS":["LAMP3","CCR7","CD83","CD40"],
 "B_lineage":["MS4A1","CD79A","CD19","BANK1"],
 "Plasma_lineage":["MZB1","JCHAIN","SDC1","PRDM1","XBP1"],
 "T_lineage":["CD3D","CD3E","CD2","TRAC"],
}
for k,gs in score_sets.items():
    gs=[g for g in gs if g in a.var_names]
    sc.tl.score_genes(a,gs,score_name=k,ctrl_size=50,n_bins=25)
log("state scores done")

# high/low within lineage
def add_state(lineage, score_col, hi_q=0.75, lo_q=0.50, min_cells=50):
    m=(a.obs["cell_type"]==lineage).values
    s=a.obs.loc[m,score_col].values
    if len(s)<min_cells: return None
    hi=np.quantile(s,hi_q); lo=np.quantile(s,lo_q)
    state=np.full(a.n_obs,"other",dtype=object)
    state[m & (a.obs[score_col].values>=hi)]="high"
    state[m & (a.obs[score_col].values<=lo)]="low"
    a.obs[f"state_{lineage}"]=state
    log(f"{lineage} state high={int((state=='high').sum())} low={int((state=='low').sum())} q={hi_q}/{lo_q}")
    return {"lineage":lineage,"score":score_col,"n_high":int((state=='high').sum()),"n_low":int((state=='low').sum()),"hi_q":hi_q,"lo_q":lo_q}

summ=[]
summ.append(add_state("B","B_TLS",0.75,0.50,50))
summ.append(add_state("T","T_TLS",0.75,0.50,100))
summ.append(add_state("Fibroblast","Fib_TLS",0.75,0.50,200))
summ.append(add_state("Endothelial","Endo_TLS",0.75,0.50,100))
summ.append(add_state("Myeloid","Myeloid_TLS",0.75,0.50,100))
pd.DataFrame([x for x in summ if x]).to_csv(os.path.join(OUT,"tls_state_summary.tsv"),sep="\t",index=False)

# ---------- DE ----------
exclude_re=re.compile(r"^(IG[HKL][VDJ]|IGKC|IGLC|RPL|RPS|MT-|XIST|RPS4Y1|MKI67|TOP2A|KRT|S100A8|S100A9|FCGR3B|CXCR2|HBB|HBA|MALAT1|NEAT1)")
pure_lineage=set("MS4A1 CD79A CD79B CD19 BANK1 CD3D CD3E CD2 TRAC CD4 CD8A MZB1 JCHAIN SDC1 PRDM1 XBP1 COL1A1 COL1A2 PDGFRA DCN LUM PECAM1 VWF KRT5 KRT14 KRT76 KRT13 CSTB".split())
de_rows=[]
state_summary={}
for lineage in ["B","T","Fibroblast","Endothelial","Myeloid"]:
    col=f"state_{lineage}"
    if col not in a.obs: continue
    sub=a[(a.obs[col].isin(["high","low"]))].copy()
    sub.obs[col]=sub.obs[col].astype("category")
    if sub.obs[col].value_counts().min()<30:
        log(f"skip DE {lineage}: too few"); continue
    log(f"DE {lineage} n={sub.n_obs}")
    sc.tl.rank_genes_groups(sub,col,groups=["high"],reference="low",method="wilcoxon",n_genes=1000,pts=True)
    df=sc.get.rank_genes_groups_df(sub,group="high")
    df=df.rename(columns={"names":"gene","logfoldchanges":"logFC","pvals_adj":"padj"})
    df["pct_high"]=df.get("pct_nz_group", np.nan)
    df["pct_low"]=df.get("pct_nz_reference", np.nan)
    df["lineage"]=lineage
    df.to_csv(os.path.join(OUT,"de_results",f"{lineage}_high_vs_low.tsv"),sep="\t",index=False)
    sel=df[(df.padj<0.05)&(df.logFC>0.25)&(df.pct_high>0.10)&(~df.gene.str.match(exclude_re))]
    log(f"{lineage} selected={len(sel)}")
    for _,r in sel.iterrows():
        de_rows.append(dict(gene=r.gene,lineage=lineage,logFC=float(r.logFC),padj=float(r.padj),pct_high=float(r.get("pct_high",np.nan)),pct_low=float(r.get("pct_low",np.nan))))
de=pd.DataFrame(de_rows)
de.to_csv(os.path.join(OUT,"de_all_selected.tsv"),sep="\t",index=False)
log(f"DE table rows={len(de)}")

# cross-lineage candidates
agg=de.groupby("gene").agg(n_lineages=("lineage","nunique"),mean_logFC=("logFC","mean"),max_logFC=("logFC","max"),
                           min_padj=("padj","min"),max_pct=("pct_high","max"),lineages=("lineage",lambda x:",".join(sorted(set(x))))).reset_index()
agg=agg.sort_values(["n_lineages","mean_logFC"],ascending=False)
agg.to_csv(os.path.join(OUT,"candidate_genes_cross_lineage.tsv"),sep="\t",index=False)
log(f"candidates={len(agg)} multi-lineage={int((agg.n_lineages>=2).sum())}")

# ---------- platform gene sets ----------
g570=set(pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_expr_gene.tsv"),sep="\t",usecols=[0]).iloc[:,0].astype(str))
g173=set(pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_logCPM.tsv"),sep="\t",usecols=[0]).iloc[:,0].astype(str))
xen_feat=pd.read_csv(os.path.join(ROOT,"data/raw/xenium/H1/features.tsv.gz"),sep="\t",header=None)
xen=set(xen_feat.loc[xen_feat[2].astype(str).eq("Gene Expression"),1].astype(str))
# candidate pool: multi-lineage first, then single-lineage
pool=agg[agg.n_lineages>=2].copy()
if len(pool)<30:
    extra=agg[agg.n_lineages==1].head(80-len(pool)); pool=pd.concat([pool,extra])
pool=pool[pool.gene.isin(g570)&pool.gene.isin(g173)]
pool=pool[~pool.gene.isin(pure_lineage)]
# core: 40 genes, prefer multi-lineage
core=pool.head(40).gene.tolist()
# ensure >=30 if possible
if len(core)<30:
    extra=agg[(agg.n_lineages>=1)&(~agg.gene.isin(core))&agg.gene.isin(g570)&agg.gene.isin(g173)&(~agg.gene.isin(pure_lineage))].head(30-len(core))
    core+=extra.gene.tolist()
ext=pool.head(120).gene.tolist()
if len(ext)<80:
    extra=agg[(agg.n_lineages>=1)&(~agg.gene.isin(ext))&agg.gene.isin(g570)&agg.gene.isin(g173)&(~agg.gene.isin(pure_lineage))].head(80-len(ext))
    ext+=extra.gene.tolist()
panel=[g for g in core if g in xen]
if len(panel)<8:
    panel += [g for g in ext if g in xen and g not in panel][: max(0,8-len(panel))]
log(f"core={len(core)} ext={len(ext)} panel={len(panel)}")

# save YAML
def save_yaml(name, genes, source):
    d={"name":name,"version":"1.0","genome":"GRCh38","source":source,"genes":genes}
    with open(os.path.join(ROOT,"signatures",name+".yaml"),"w") as fh: yaml.safe_dump(d,fh,sort_keys=False)
save_yaml("OTLS_core",core,"GSE164241 DE high-vs-low TLS states; cross-lineage; filtered for GPL570/RNA-seq")
save_yaml("OTLS_ext",ext,"GSE164241 DE extended")
save_yaml("OTLS_panel",panel,"OTLS_core intersect Xenium 451 panel (supplemented from ext if <8)")
controls={
 "Bcell":["MS4A1","CD79A","CD79B","CD19","BANK1"],
 "Plasma":["MZB1","JCHAIN","SDC1","PRDM1","XBP1","TNFRSF17","IGHG1","IGHA1"],
 "Tcell":["CD3D","CD3E","CD2","TRAC","CD52","IL7R"],
 "12CK":["CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13"],
 "Cabrita9":["PTGDS","RBP5","EIF1AY","CETP","SKAP1","LAT","CCR6","CD1D","CD79B"],
 "TLS39":["CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1","CD4","CCR5","CXCR3","CSF2","IGSF6","IL2RA","CD38","CD40","CD5","MS4A1","SDC1","GFI1","IL1R1","IL1R2","IL10","CCL20","IRF4","TRAF6","STAT5A","TNFRSF17"],
 "Inflammation":["IFNG","STAT1","CXCL9","CXCL10","IL6","CXCL8","NFKB1","RELA"],
}
for k,v in controls.items(): save_yaml(k,v,"control signature")
with open(os.path.join(ROOT,"signatures","signature_evidence.tsv"),"w") as fh:
    fh.write("gene\tsignature\tn_lineages\tmean_logFC\tlineages\n")
    for g in core:
        r=agg[agg.gene==g].iloc[0]; fh.write(f"{g}\tcore\t{r.n_lineages}\t{r.mean_logFC:.3f}\t{r.lineages}\n")
    for g in ext:
        if g in core: continue
        r=agg[agg.gene==g].iloc[0]; fh.write(f"{g}\text\t{r.n_lineages}\t{r.mean_logFC:.3f}\t{r.lineages}\n")

# ---------- score discovery cells ----------
for k,gs in [("OTLS_core",core),("OTLS_ext",ext),("OTLS_panel",panel)]+[(k,v) for k,v in controls.items()]:
    gs=[g for g in gs if g in a.var_names]
    if len(gs)>=2:
        sc.tl.score_genes(a,gs,score_name=k,ctrl_size=50,n_bins=25)
# targeted LR support
lr=[("CXCL13","CXCR5"),("CCL19","CCR7"),("CCL21","CCR7"),("LTB","LTBR"),("TNFSF13B","TNFRSF13B"),("TNFSF13B","TNFRSF13C"),("ICAM1","ITGAL"),("SELE","SELPLG")]
lr_rows=[]
for lig,rec in lr:
    if lig in a.var_names and rec in a.var_names:
        for lineage in ["B","T","Fibroblast","Endothelial"]:
            col=f"state_{lineage}"
            if col not in a.obs: continue
            for st in ["high","low"]:
                m=(a.obs[col]==st).values
                if m.sum()<10: continue
                lr_rows.append(dict(ligand=lig,receptor=rec,lineage=lineage,state=st,n=int(m.sum()),
                                    ligand_mean=float(np.asarray(a[m,lig].X.todense() if sp.issparse(a[m,lig].X) else a[m,lig].X).ravel().mean()),
                                    receptor_mean=float(np.asarray(a[m,rec].X.todense() if sp.issparse(a[m,rec].X) else a[m,rec].X).ravel().mean())))
pd.DataFrame(lr_rows).to_csv(os.path.join(OUT,"ligand_receptor_support.tsv"),sep="\t",index=False)
# save obs table
obs=a.obs.copy(); obs.to_csv(os.path.join(OUT,"celltype_metadata.tsv"),sep="\t")
# save signatures summary
pd.DataFrame({"signature":["OTLS_core","OTLS_ext","OTLS_panel"],"n_genes":[len(core),len(ext),len(panel)],
              "genes":[";".join(core),";".join(ext),";".join(panel)]}).to_csv(os.path.join(OUT,"signature_summary.tsv"),sep="\t",index=False)
# symlink integrated h5ad for output
link=os.path.join(OUT,"integrated.h5ad")
if not os.path.exists(link): os.symlink(h5,link)
log("DONE")
