#!/usr/bin/env python3
"""Narrative B: B/T/APC co-localization score in bulk cohorts."""
import os, gzip, io, time, warnings, yaml
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy.stats import wilcoxon, spearmanr
from statsmodels.miscmodels.ordinal_model import OrderedModel
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/narrative_B/bulk"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(OUT,"bulk.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
APC=["HLA-DRB5","HLA-DMB","HLA-DQA2","HLA-DQB2","CD14","CD68","CD163","CD1C","FCER1A","LAMP3","CCR7","CD83"]
def ssgsea(x, genes, gene_set):
    x=np.asarray(x,float); n=len(x); order=np.argsort(-x); gs=np.asarray(genes)
    hit=np.isin(gs[order],gene_set); ranks=np.arange(1,n+1,dtype=float)
    w=np.where(hit,ranks**0.25,0); es=np.cumsum(w/w.sum()-((~hit).astype(float)/(~hit).sum()))
    return float(es[np.argmax(np.abs(es))])
def add_apc(df, cohort):
    if not any(g in df.index for g in APC): return None
    s=pd.Series({c:ssgsea(df[c].values,df.index,APC) for c in df.columns},name="APC")
    s=(s-s.mean())/(s.std()+1e-12); s.to_csv(os.path.join(OUT,f"{cohort}_APC.tsv"),sep="\t")
    log(f"{cohort} APC genes present={sum(g in df.index for g in APC)}/12")
    return s
def zscore(s): return (s-s.mean())/(s.std()+1e-12)
# GSE16134
sc=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE16134_scores.tsv"),sep="\t",index_col=0)
e=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_expr_gene.tsv"),sep="\t",index_col=0); apc=add_apc(e,"GSE16134")
meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t").set_index("gsm")
d=sc.join(meta[["patient_id","site_status"]],how="inner"); d["APC"]=apc.reindex(d.index).values; d["APC"]=zscore(d.APC)
for c in ["OTLS_core","OTLS_core_independent","TLS39","Bcell","Plasma","Tcell"]: d[c]=zscore(d[c])
d["BTAPC_min"]=np.minimum.reduce([d.Bcell,d.Tcell,d.APC]); d["BTAPC_mean"]=(d.Bcell+d.Tcell+d.APC)/3
d["BPT"]=d[["Bcell","Plasma","Tcell"]].mean(1)
rows=[]
for s in ["OTLS_core","OTLS_core_independent","TLS39","Bcell","Tcell","APC","BTAPC_min","BTAPC_mean","BPT"]:
    piv=d.pivot_table(index="patient_id",columns="site_status",values=s,aggfunc="mean").dropna()
    diff=piv.affected-piv.unaffected
    rows.append(dict(cohort="GSE16134",analysis="paired",score=s,n=len(piv),effect=float(diff.mean()),cohen_dz=float(diff.mean()/(diff.std()+1e-12)),p=float(wilcoxon(piv.affected,piv.unaffected).pvalue)))
    try:
        m=smf.mixedlm(f"{s} ~ C(site_status)",d,groups=d.patient_id).fit()
        b=m.params.get("C(site_status)[T.unaffected]",np.nan); p=m.pvalues.get("C(site_status)[T.unaffected]",np.nan)
        rows.append(dict(cohort="GSE16134",analysis="LMM_status",score=s,n=len(d),effect=float(b),cohen_dz=np.nan,p=float(p)))
    except Exception as e: log(f"LMM {s} {e}")
# AUC patient-grouped for status (sample-level, but group by patient)
y=(d.site_status=="affected").astype(int).values; groups=d.patient_id.astype("category").cat.codes.values
for s in ["OTLS_core","BPT","BTAPC_min","BTAPC_mean","TLS39"]:
    X=d[[s]].values.astype(float)
    try:
        pred=cross_val_predict(LogisticRegression(max_iter=500),X,y,groups=groups,cv=GroupKFold(n_splits=5),method="predict_proba")[:,1]
        rows.append(dict(cohort="GSE16134",analysis="AUC_grouped5",score=s,n=len(d),effect=float(roc_auc_score(y,pred)),cohen_dz=np.nan,p=np.nan))
    except Exception as e: log(f"AUC {s} {e}")
# GSE10334 technical
e2=pd.read_csv(os.path.join(ROOT,"data/external/GSE10334_expr_gene.tsv"),sep="\t",index_col=0); apc2=add_apc(e2,"GSE10334")
s2=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE10334_scores.tsv"),sep="\t",index_col=0); m2=pd.read_csv(os.path.join(ROOT,"data/external/GSE10334_meta.tsv"),sep="\t")
# GSE173078
e3=pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_logCPM.tsv"),sep="\t",index_col=0); apc3=add_apc(e3,"GSE173078")
s3=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE173078_scores.tsv"),sep="\t",index_col=0); m3=pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_meta.tsv"),sep="\t")
d3=s3.join(m3.set_index("sample_id")[["site_status"]],how="inner");
for c in ["OTLS_core","TLS39","Bcell","Plasma","Tcell"]: d3[c]=zscore(d3[c])
d3["APC"]=zscore(apc3.reindex(d3.index)); d3["BTAPC_min"]=np.minimum.reduce([zscore(d3.Bcell),zscore(d3.Tcell),d3.APC]); d3["BTAPC_mean"]=(zscore(d3.Bcell)+zscore(d3.Tcell)+d3.APC)/3; d3["BPT"]=d3[["Bcell","Plasma","Tcell"]].mean(1)
d3=d3.dropna(subset=["site_status"]); d3["sev"]=pd.Categorical(d3.site_status,categories=["healthy","gingivitis","periodontitis"],ordered=True).codes.astype(int)
for s in ["OTLS_core","BPT","BTAPC_min","BTAPC_mean","APC"]:
    try:
        m=OrderedModel(d3.sev.values,d3[[s]].values,distr="logit").fit(method="bfgs",disp=False,maxiter=500)
        rows.append(dict(cohort="GSE173078",analysis="ordinal",score=s,n=len(d3),effect=float(np.exp(m.params[0])),cohen_dz=float(spearmanr(d3[s],d3.sev).correlation),p=float(m.pvalues[0])))
    except Exception as e: log(f"ordinal {s} {e}")
# GSE223328
lines=[]; tab=False
with gzip.open(os.path.join(ROOT,"data/raw/geo/GSE223328/GSE223328_series_matrix.txt.gz"),"rt",errors="ignore") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"): tab=True; lines.append(next(fh)); continue
        if not tab: continue
        if line.startswith("!series_matrix_table_end"): break
        lines.append(line)
x4=pd.read_csv(io.StringIO("".join(lines)),sep="\t",index_col=0); mp={}
with gzip.open(os.path.join(ROOT,"data/external/GPL21185.txt"),"rt",errors="ignore") as fh:
    tab=False
    for line in fh:
        if line.startswith("!platform_table_begin"): tab=True; next(fh); continue
        if not tab: continue
        if line.startswith("!platform_table_end"): break
        p=line.rstrip("\n").split("\t")
        if len(p)>5 and p[5]: mp[p[0]]=p[5]
x4.index=[mp.get(i,i) for i in x4.index]; x4=x4.groupby(level=0).mean(); apc4=add_apc(x4,"GSE223328")
s4=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE223328_scores.tsv"),sep="\t",index_col=0); m4=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE223328_meta.tsv"),sep="\t")
d4=s4.join(m4.set_index("sample"),how="inner");
for c in ["OTLS_core","TLS39","Bcell","Plasma","Tcell"]: d4[c]=zscore(d4[c])
d4["APC"]=zscore(apc4.reindex(d4.index)); d4["BTAPC_min"]=np.minimum.reduce([zscore(d4.Bcell),zscore(d4.Tcell),d4.APC]); d4["BTAPC_mean"]=(zscore(d4.Bcell)+zscore(d4.Tcell)+d4.APC)/3; d4["BPT"]=d4[["Bcell","Plasma","Tcell"]].mean(1)
d4["status"]=(d4.condition=="periodontitis").astype(int); d4["ebv"]=d4.notes.str.contains("positive").astype(int)
for s in ["OTLS_core","BPT","BTAPC_min","BTAPC_mean","APC"]:
    try:
        X=sm.add_constant(d4[[s,"ebv"]]); m=sm.Logit(d4.status,X).fit(disp=0)
        rows.append(dict(cohort="GSE223328",analysis="logit_ebv",score=s,n=len(d4),effect=float(np.exp(m.params[s])),cohen_dz=np.nan,p=float(m.pvalues[s])))
    except Exception as e: log(f"logit {s} {e}")
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"btapc_association.tsv"),sep="\t",index=False)
print(res.to_string(index=False)); log("DONE")
