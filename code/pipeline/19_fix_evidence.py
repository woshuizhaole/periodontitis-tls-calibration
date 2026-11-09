#!/usr/bin/env python3
import pandas as pd, yaml, os
ROOT="/data/workplace/wyx/perio_tls_20261107"
seed={"Fibroblast":["CXCL13","CCL19","CCL21","LTB","TNFSF13B","CXCL12","FDCSP","CR2","PTGDS","VCAN","FBN1"],
"Endothelial":["PLVAP","MMRN1","ACKR1","SELE","ICAM1","VCAM1","VWF"],
"B":["BCL6","TNFRSF13B","CD27","CD38","MME","BCL2A1"],
"T":["CXCR5","PDCD1","ICOS","IL21","TOX","CD200","TIGIT","CTLA4"],
"Myeloid":["LAMP3","CCR7","CD83","CD40","IRF4","CLEC9A"],"Other":["CCL20","IL6ST"]}
expected={g:lin for lin,gs in seed.items() for g in gs}
e=pd.read_csv(os.path.join(ROOT,"results/stage1/signature_evidence_v2.tsv"),sep="\t")
for c in ["expected_lineage","support_lineage"]: e[c]=e[c].astype(object)
for i,r in e.iterrows():
    if r.signature=="core":
        e.at[i,"expected_lineage"]=expected.get(r.gene,"NA")
        if r.support_lineage in expected.values() and r.support_lineage!=expected.get(r.gene):
            # if expected lineage also has significant DE, use expected lineage's statistics
            lin=expected.get(r.gene)
            p=f"{ROOT}/results/stage1/de_results/{lin}_high_vs_low.tsv"
            if os.path.exists(p):
                d=pd.read_csv(p,sep="\t"); x=d[d.gene==r.gene]
                if len(x) and float(x.padj.iloc[0])<0.05 and float(x.logFC.iloc[0])>0.5:
                    e.at[i,"support_lineage"]=lin; e.at[i,"logFC"]=float(x.logFC.iloc[0]); e.at[i,"padj"]=float(x.padj.iloc[0]); e.at[i,"pct_high"]=float(x.get("pct_high",pd.Series([float("nan")])).iloc[0])
e.to_csv(os.path.join(ROOT,"results/stage1/signature_evidence_v2.tsv"),sep="\t",index=False)
print(e[e.gene.isin(["LAMP3","CD200","BCL6","CCL20","RBP5","BCL2A1"])].to_string(index=False))
