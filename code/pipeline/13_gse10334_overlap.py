#!/usr/bin/env python3
import os, re, numpy as np, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"metadata/gse10334_gse16134_overlap.tsv")
m1=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_meta.tsv"),sep="\t")
m2=pd.read_csv(os.path.join(ROOT,"data/external/GSE10334_meta.tsv"),sep="\t")
def parse(m):
    d=m.copy()
    d["patient"]=d["patient_id"].astype(str)
    d["sample_no"]=d["sample_id"].astype(str).str.extract(r"sample_(\d+)")[0]
    d["site"]=d["site_status"].astype(str)
    return d
a=parse(m1); b=parse(m2)
a["key"]=a.patient+"|"+a.sample_no+"|"+a.site; b["key"]=b.patient+"|"+b.sample_no+"|"+b.site
common=sorted(set(a.key)&set(b.key))
print("GSE16134",a.shape,"GSE10334",b.shape,"common keys",len(common))
# expression
e1=pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_expr_gene.tsv"),sep="\t",index_col=0)
e2=pd.read_csv(os.path.join(ROOT,"data/external/GSE10334_expr_gene.tsv"),sep="\t",index_col=0)
probes=e1.index.intersection(e2.index); print("common genes",len(probes))
rows=[]
for k in common:
    ra=a[a.key==k].iloc[0]; rb=b[b.key==k].iloc[0]
    x=e1.loc[probes,ra.gsm].values.astype(float); y=e2.loc[probes,rb.gsm].values.astype(float)
    ok=np.isfinite(x)&np.isfinite(y)
    rows.append(dict(patient=ra.patient,sample_no=ra.sample_no,site=ra.site,
                     gsm_16134=ra.gsm,gsm_10334=rb.gsm,n_genes=int(ok.sum()),
                     pearson=float(np.corrcoef(x[ok],y[ok])[0,1])))
df=pd.DataFrame(rows); df.to_csv(OUT,sep="\t",index=False)
print(df.pearson.describe().to_string())
print("patients shared", len(set(a.patient)&set(b.patient)), "GSE10334 patients",b.patient.nunique())
print("unmatched GSE10334 keys",len(set(b.key)-set(a.key)))
print(df.sort_values("pearson").head(8).to_string(index=False))
