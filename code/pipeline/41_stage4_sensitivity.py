#!/usr/bin/env python3
"""Stage 4: sensitivity analyses and multiplicity summary."""
import os, pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage4"); os.makedirs(OUT,exist_ok=True)
# GSE16134 all samples (unpaired) vs 69 complete pairs
sc=pd.read_csv(os.path.join(ROOT,"results/stage3/scores/GSE16134_scores.tsv"),sep="\t",index_col=0)
meta=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t")
df=sc.join(meta.set_index("gsm")[["patient_id","site_status"]],how="inner")
rows=[]
for s in ["OTLS_core","OTLS_panel","OTLS_core_independent"]:
    a=df[df.site_status=="affected"][s].dropna(); u=df[df.site_status=="unaffected"][s].dropna()
    p=mannwhitneyu(a,u).pvalue; d=(a.mean()-u.mean())/(df[s].std()+1e-12)
    rows.append(dict(analysis="GSE16134_all_unpaired",signature=s,n_affected=len(a),n_unaffected=len(u),effect=float(d),p=float(p)))
sens=pd.DataFrame(rows); sens.to_csv(os.path.join(OUT,"gse16134_all_unpaired.tsv"),sep="\t",index=False)
# BH on GSE16134 paired p
p=pd.read_csv(os.path.join(ROOT,"results/stage3/stats/gse16134_paired.tsv"),sep="\t")
p["p_BH"]=multipletests(p.p_wilcoxon,method="fdr_bh")[1]; p.to_csv(os.path.join(OUT,"gse16134_paired_BH.tsv"),sep="\t",index=False)
# Xenium all vs main from tests
t=pd.read_csv(os.path.join(ROOT,"results/stage2/xenium/tests.tsv"),sep="\t")
t.to_csv(os.path.join(OUT,"xenium_tests_reference.tsv"),sep="\t",index=False)
# Stage 3 model summary
summary=[]
for f,key in [("gse16134_paired.tsv","signature"),("gse173078_ordinal.tsv","signature"),("gse223328_logistic.tsv","signature")]:
    d=pd.read_csv(os.path.join(ROOT,"results/stage3/stats",f),sep="\t")
    summary.append(d.assign(table=f))
pd.concat(summary,ignore_index=True).to_csv(os.path.join(OUT,"stage3_model_summary.tsv"),sep="\t",index=False)
print(sens.to_string(index=False)); print(p[["signature","p_wilcoxon","p_BH"]].head().to_string(index=False))
