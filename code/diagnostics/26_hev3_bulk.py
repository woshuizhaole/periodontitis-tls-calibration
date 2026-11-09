#!/usr/bin/env python3
"""Diagnostic: HEV3 (CCL19/MMRN1/SELE) bulk association, independent of full panel."""
import os, gzip, io, numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy.stats import wilcoxon, spearmanr
from statsmodels.miscmodels.ordinal_model import OrderedModel
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"diagnostics/bulk/hev3"); os.makedirs(OUT,exist_ok=True)
GENES=["CCL19","MMRN1","SELE"]
def ssgsea(x, genes):
    x=np.asarray(x,float); n=len(x); order=np.argsort(-x); gs=np.asarray(genes)
    hit=np.isin(gs[order],GENES); ranks=np.arange(1,n+1,dtype=float)
    w=np.where(hit,ranks**0.25,0); es=np.cumsum(w/w.sum()-((~hit).astype(float)/(~hit).sum()))
    return float(es[np.argmax(np.abs(es))])
def add(df, cohort):
    s=pd.Series({c:ssgsea(df[c].values,df.index) for c in df.columns},name="HEV3")
    s=(s-s.mean())/(s.std()+1e-12); s.to_csv(os.path.join(OUT,f"{cohort}_hev3.tsv"),sep="\t")
    return s
res=[]
# GSE16134
sc=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE16134_scores.tsv"),sep="\t",index_col=0)
e=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_expr_gene.tsv"),sep="\t",index_col=0)
h=add(e,"GSE16134"); meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t").set_index("gsm")
d=sc.join(meta[["patient_id","site_status"]],how="inner"); d["HEV3"]=h.reindex(d.index).values
piv=d.pivot_table(index="patient_id",columns="site_status",values="HEV3",aggfunc="mean").dropna()
diff=piv.affected-piv.unaffected
res.append(dict(cohort="GSE16134",analysis="paired",n=len(piv),effect=float(diff.mean()),cohen_dz=float(diff.mean()/(diff.std()+1e-12)),p=float(wilcoxon(piv.affected,piv.unaffected).pvalue)))
for adj in [[],["Bcell"],["Bcell","Plasma","Tcell"],["Bcell","Plasma","Tcell","stromal_hev"]]:
    if "stromal_hev" in adj and "stromal_hev" not in d: continue
    f="HEV3 ~ C(site_status)"+("" if not adj else "+"+"+".join(adj))
    m=smf.mixedlm(f,d,groups=d.patient_id).fit()
    b=m.params.get("C(site_status)[T.unaffected]",np.nan); p=m.pvalues.get("C(site_status)[T.unaffected]",np.nan)
    res.append(dict(cohort="GSE16134",analysis="LMM_"+("none" if not adj else "+".join(adj)),n=len(d),effect=float(b),cohen_dz=np.nan,p=float(p)))
# GSE10334 technical
e2=pd.read_csv(os.path.join(ROOT,"data/external/GSE10334_expr_gene.tsv"),sep="\t",index_col=0); h2=add(e2,"GSE10334")
# GSE173078 ordinal
sc3=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE173078_scores.tsv"),sep="\t",index_col=0)
e3=pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_logCPM.tsv"),sep="\t",index_col=0); h3=add(e3,"GSE173078")
m3=pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_meta.tsv"),sep="\t")
d3=sc3.join(m3.set_index("sample_id")[["site_status"]],how="inner"); d3["HEV3"]=h3.reindex(d3.index).values
d3=d3.dropna(subset=["HEV3","site_status"]); d3["sev"]=pd.Categorical(d3.site_status,categories=["healthy","gingivitis","periodontitis"],ordered=True).codes.astype(int)
X=d3[["HEV3"]].values; m=OrderedModel(d3.sev.values,X,distr="logit").fit(method="bfgs",disp=False,maxiter=500)
res.append(dict(cohort="GSE173078",analysis="ordinal",n=len(d3),effect=float(np.exp(m.params[0])),cohen_dz=float(spearmanr(d3.HEV3,d3.sev).correlation),p=float(m.pvalues[0])))
# GSE223328
s4=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE223328_scores.tsv"),sep="\t",index_col=0)
# parse series matrix
lines=[]; tab=False
with gzip.open(os.path.join(ROOT,"data/raw/geo/GSE223328/GSE223328_series_matrix.txt.gz"),"rt",errors="ignore") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"): tab=True; lines.append(next(fh)); continue
        if not tab: continue
        if line.startswith("!series_matrix_table_end"): break
        lines.append(line)
x4=pd.read_csv(io.StringIO("".join(lines)),sep="\t",index_col=0)
mp={}
with gzip.open(os.path.join(ROOT,"data/external/GPL21185.txt"),"rt",errors="ignore") as fh:
    tab=False
    for line in fh:
        if line.startswith("!platform_table_begin"): tab=True; next(fh); continue
        if not tab: continue
        if line.startswith("!platform_table_end"): break
        p=line.rstrip("\n").split("\t")
        if len(p)>5 and p[5]: mp[p[0]]=p[5]
x4.index=[mp.get(i,i) for i in x4.index]; x4=x4.groupby(level=0).mean(); h4=add(x4,"GSE223328")
meta4=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE223328_meta.tsv"),sep="\t")
d4=s4.join(meta4.set_index("sample"),how="inner"); d4["HEV3"]=h4.reindex(d4.index).values; d4["status"]=(d4.condition=="periodontitis").astype(int); d4["ebv"]=d4.notes.str.contains("positive").astype(int)
X=sm.add_constant(d4[["HEV3","ebv"]]); m=sm.Logit(d4.status,X).fit(disp=0)
res.append(dict(cohort="GSE223328",analysis="logit_ebv",n=len(d4),effect=float(np.exp(m.params["HEV3"])),cohen_dz=np.nan,p=float(m.pvalues["HEV3"])))
pd.DataFrame(res).to_csv(os.path.join(OUT,"hev3_association.tsv"),sep="\t",index=False)
print(pd.DataFrame(res).to_string(index=False))
