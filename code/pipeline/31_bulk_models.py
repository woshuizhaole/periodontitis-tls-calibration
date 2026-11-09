#!/usr/bin/env python3
"""Stage 3b: bulk association models, technical replication, signature benchmark, random null."""
import os, time, warnings, yaml
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, ttest_rel, spearmanr, pearsonr
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.miscmodels.ordinal_model import OrderedModel
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; S=os.path.join(ROOT,"results/stage3/scores"); OUT=os.path.join(ROOT,"results/stage3/stats"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(ROOT,"logs/analysis/stage3_models.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
SIGS=["OTLS_core","OTLS_ext","OTLS_panel","OTLS_core_independent","Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]
def z(df,cols):
    for c in cols:
        if c in df: df[c]=(df[c]-df[c].mean())/(df[c].std()+1e-12)
    return df
# ---------- GSE16134 paired ----------
log("GSE16134")
sc=pd.read_csv(os.path.join(S,"GSE16134_scores.tsv"),sep="\t",index_col=0)
meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t")
df=sc.join(meta.set_index("gsm")[["patient_id","site_status","sample_id"]],how="inner")
df=z(df,SIGS)
rows=[]
for s in SIGS:
    piv=df.pivot_table(index="patient_id",columns="site_status",values=s,aggfunc="mean")
    piv=piv.dropna(subset=["affected","unaffected"])
    if len(piv)<5: continue
    d=piv["affected"]-piv["unaffected"]
    p=wilcoxon(piv["affected"],piv["unaffected"]).pvalue
    pt=ttest_rel(piv["affected"],piv["unaffected"]).pvalue
    rows.append(dict(signature=s,n_pairs=len(piv),mean_diff=float(d.mean()),median_diff=float(d.median()),
                     cohen_dz=float(d.mean()/(d.std()+1e-12)),p_wilcoxon=float(p),p_paired_t=float(pt),
                     affected_mean=float(piv["affected"].mean()),unaffected_mean=float(piv["unaffected"].mean())))
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"gse16134_paired.tsv"),sep="\t",index=False); log(res.to_string(index=False))
# LMM for OTLS_core and controls
for s in ["OTLS_core","OTLS_panel","OTLS_core_independent","BTAPC" if "BTAPC" in df else "Tcell"]:
    if s not in df: continue
    try:
        m=smf.mixedlm(f"{s} ~ C(site_status)",df,groups=df["patient_id"]).fit()
        b=m.params.get("C(site_status)[T.unaffected]",np.nan); p=m.pvalues.get("C(site_status)[T.unaffected]",np.nan)
        log(f"LMM {s} beta(unaffected vs affected)={b:.4f} p={p:.4g}")
    except Exception as e: log(f"LMM {s} failed {e}")
# adjusted model: status + B/Plasma/T
for adj in [[],["Bcell"],["Bcell","Plasma","Tcell"]]:
    if "OTLS_core" not in df: break
    f="OTLS_core ~ C(site_status)"+("" if not adj else "+"+"+".join(adj))
    try:
        m=smf.mixedlm(f,df,groups=df["patient_id"]).fit()
        b=m.params.get("C(site_status)[T.unaffected]",np.nan); p=m.pvalues.get("C(site_status)[T.unaffected]",np.nan)
        log(f"adjusted {adj} beta={b:.4f} p={p:.4g}")
    except Exception as e: log(f"adjusted {adj} failed {e}")
# ---------- GSE173078 ordinal ----------
log("GSE173078")
sc=pd.read_csv(os.path.join(S,"GSE173078_scores.tsv"),sep="\t",index_col=0)
meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_meta.tsv"),sep="\t")
df=sc.join(meta.set_index("sample_id")[["site_status"]],how="inner")
df=z(df,SIGS)
rows=[]
for s in SIGS:
    if s not in df: continue
    d=df[[s,"site_status"]].dropna(); d["sev"]=pd.Categorical(d.site_status,categories=["healthy","gingivitis","periodontitis"],ordered=True)
    d=d.dropna(subset=["sev"])
    if d.sev.nunique()<3: continue
    y=d.sev.cat.codes.values.astype(int); X=d[[s]].values
    try:
        m=OrderedModel(y,X,distr="logit").fit(method="bfgs",disp=False,maxiter=500)
        b=m.params[0]; se=m.bse[0]; p=m.pvalues[0]
        rows.append(dict(signature=s,n=len(d),OR_per_SD=float(np.exp(b)),CI_low=float(np.exp(b-1.96*se)),CI_high=float(np.exp(b+1.96*se)),p=float(p),
                         spearman=float(spearmanr(d[s],y).correlation),spearman_p=float(spearmanr(d[s],y).pvalue)))
    except Exception as e: log(f"ordinal {s} failed {e}")
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"gse173078_ordinal.tsv"),sep="\t",index=False); log(res.to_string(index=False))
# ---------- GSE223328 logistic ----------
log("GSE223328")
sc=pd.read_csv(os.path.join(S,"GSE223328_scores.tsv"),sep="\t",index_col=0)
meta=pd.read_csv(os.path.join(S,"GSE223328_meta.tsv"),sep="\t")
df=sc.join(meta.set_index("sample"),how="inner"); df=z(df,SIGS)
df["status"]=(df.condition=="periodontitis").astype(int); df["ebv"]=df.notes.str.contains("positive").astype(int)
rows=[]
for s in SIGS:
    if s not in df or df.status.nunique()<2: continue
    X=sm.add_constant(df[[s,"ebv"]]); 
    try:
        m=sm.Logit(df.status,X).fit(disp=0)
        b=m.params[s]; se=m.bse[s]; p=m.pvalues[s]
        # 4-fold stratified CV AUC (LOOCV unstable with n=16 and prior shifts)
        cv=StratifiedKFold(n_splits=4,shuffle=True,random_state=42)
        pred=cross_val_predict(LogisticRegression(max_iter=500),df[[s]],df.status,cv=cv,method="predict_proba")[:,1]
        auc=roc_auc_score(df.status,pred)
        rows.append(dict(signature=s,n=len(df),OR_per_SD=float(np.exp(b)),CI_low=float(np.exp(b-1.96*se)),CI_high=float(np.exp(b+1.96*se)),p=float(p),CV4_AUC=float(auc)))
    except Exception as e: log(f"logit {s} failed {e}")
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"gse223328_logistic.tsv"),sep="\t",index=False); log(res.to_string(index=False))
# ---------- GSE10334 technical replication ----------
log("GSE10334")
s1=pd.read_csv(os.path.join(S,"GSE16134_scores.tsv"),sep="\t",index_col=0)
s2=pd.read_csv(os.path.join(S,"GSE10334_scores.tsv"),sep="\t",index_col=0)
m1=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t")
m2=pd.read_csv(os.path.join(ROOT,"data/external/GSE10334_meta.tsv"),sep="\t")
for m in [m1,m2]:
    m["patient"]=m.patient_id.astype(str); m["sample_no"]=m.sample_id.astype(str).str.extract(r"sample_(\d+)")[0]; m["key"]=m.patient+"|"+m.sample_no+"|"+m.site_status
rows=[]
for s in SIGS:
    a=s1[[s]].rename(columns={s:"gse16134"}); b=s2[[s]].rename(columns={s:"gse10334"})
    a=a.join(m1.set_index("gsm")["key"]); b=b.join(m2.set_index("gsm")["key"])
    j=a.dropna().merge(b.dropna(),on="key",suffixes=("_1","_2"))
    if len(j)>=10:
        rows.append(dict(signature=s,n=len(j),pearson=float(pearsonr(j.gse16134,j.gse10334)[0]),spearman=float(spearmanr(j.gse16134,j.gse10334).correlation)))
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"gse10334_replication.tsv"),sep="\t",index=False); log(res.to_string(index=False))
log("DONE")
