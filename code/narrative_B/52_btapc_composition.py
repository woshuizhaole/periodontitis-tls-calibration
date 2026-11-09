#!/usr/bin/env python3
"""Compare B/Plasma/T/APC composition scores in bulk; refine narrative B."""
import os, pandas as pd, numpy as np
from scipy.stats import wilcoxon
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/narrative_B/bulk"); os.makedirs(OUT,exist_ok=True)
def z(s): return (s-s.mean())/(s.std()+1e-12)
rows=[]
# GSE16134
sc=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE16134_scores.tsv"),sep="\t",index_col=0)
apc=pd.read_csv(os.path.join(OUT,"GSE16134_APC.tsv"),sep="\t",index_col=0).iloc[:,0]
meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t").set_index("gsm")
d=sc.join(meta[["patient_id","site_status"]],how="inner"); d["APC"]=z(apc.reindex(d.index)); d["BPTA_mean"]=d[["Bcell","Plasma","Tcell","APC"]].apply(z).mean(1); d["BPTA_min"]=np.minimum.reduce([z(d[c]) for c in ["Bcell","Plasma","Tcell","APC"]])
for s in ["BPTA_mean","BPTA_min"]:
    piv=d.pivot_table(index="patient_id",columns="site_status",values=s,aggfunc="mean").dropna(); diff=piv.affected-piv.unaffected
    rows.append(dict(cohort="GSE16134",score=s,analysis="paired",n=len(piv),effect=float(diff.mean()),cohen_dz=float(diff.mean()/(diff.std()+1e-12)),p=float(wilcoxon(piv.affected,piv.unaffected).pvalue)))
y=(d.site_status=="affected").astype(int).values; groups=d.patient_id.astype("category").cat.codes.values
for s in ["BPTA_mean","BPTA_min"]:
    pred=cross_val_predict(LogisticRegression(max_iter=500),d[[s]].values,y,groups=groups,cv=GroupKFold(n_splits=5),method="predict_proba")[:,1]
    rows.append(dict(cohort="GSE16134",score=s,analysis="AUC_grouped5",n=len(d),effect=float(roc_auc_score(y,pred)),cohen_dz=np.nan,p=np.nan))
# GSE173078
sc3=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE173078_scores.tsv"),sep="\t",index_col=0); a3=pd.read_csv(os.path.join(OUT,"GSE173078_APC.tsv"),sep="\t",index_col=0).iloc[:,0]
m3=pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_meta.tsv"),sep="\t"); d3=sc3.join(m3.set_index("sample_id")[["site_status"]],how="inner"); d3["APC"]=z(a3.reindex(d3.index)); d3["BPTA_mean"]=d3[["Bcell","Plasma","Tcell","APC"]].apply(z).mean(1); d3["BPTA_min"]=np.minimum.reduce([z(d3[c]) for c in ["Bcell","Plasma","Tcell","APC"]]); d3=d3.dropna(subset=["site_status"])
from statsmodels.miscmodels.ordinal_model import OrderedModel
d3["sev"]=pd.Categorical(d3.site_status,categories=["healthy","gingivitis","periodontitis"],ordered=True).codes.astype(int)
for s in ["BPTA_mean","BPTA_min"]:
    m=OrderedModel(d3.sev.values,d3[[s]].values,distr="logit").fit(method="bfgs",disp=False,maxiter=500)
    rows.append(dict(cohort="GSE173078",score=s,analysis="ordinal",n=len(d3),effect=float(np.exp(m.params[0])),cohen_dz=np.nan,p=float(m.pvalues[0])))
# GSE223328
s4=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE223328_scores.tsv"),sep="\t",index_col=0); a4=pd.read_csv(os.path.join(OUT,"GSE223328_APC.tsv"),sep="\t",index_col=0).iloc[:,0]
m4=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE223328_meta.tsv"),sep="\t"); d4=s4.join(m4.set_index("sample"),how="inner"); d4["APC"]=z(a4.reindex(d4.index)); d4["BPTA_mean"]=d4[["Bcell","Plasma","Tcell","APC"]].apply(z).mean(1); d4["BPTA_min"]=np.minimum.reduce([z(d4[c]) for c in ["Bcell","Plasma","Tcell","APC"]]); d4["status"]=(d4.condition=="periodontitis").astype(int); d4["ebv"]=d4.notes.str.contains("positive").astype(int)
import statsmodels.api as sm
for s in ["BPTA_mean","BPTA_min"]:
    X=sm.add_constant(d4[[s,"ebv"]]); m=sm.Logit(d4.status,X).fit(disp=0)
    rows.append(dict(cohort="GSE223328",score=s,analysis="logit_ebv",n=len(d4),effect=float(np.exp(m.params[s])),cohen_dz=np.nan,p=float(m.pvalues[s])))
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"BPTA_composition.tsv"),sep="\t",index=False); print(res.to_string(index=False))
