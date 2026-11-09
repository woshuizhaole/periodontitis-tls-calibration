#!/usr/bin/env python3
import os, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"
META=os.path.join(ROOT,"metadata"); os.makedirs(META,exist_ok=True)
CAND="""MS4A1 CD79A CD79B CD19 BANK1 CD74 HLA-DRA MZB1 JCHAIN SDC1 PRDM1 XBP1 TNFRSF17 IGHG1 IGHA1
CXCR5 PDCD1 ICOS BCL6 TOX TCF7 IL21 CD40LG CD200 CXCL13 CR2 FDCSP LTB TNFSF13B CCL19 CCL21
PLVAP MMRN1 ACKR1 VWF SELE ICAM1 VCAM1 CD34 LAMP3 CCR7 CD83 CD40 IRF4 CLEC9A
COL1A1 PDGFRA FAP THY1 CXCL12 IL33 AICDA MEF2B MME BCL2A1 CD27 TNFRSF13B
HLA-DRB5 HLA-DMB HLA-DQA2 HLA-DQB2 CD14 CD68 CD163 CD1C FCER1A CD3D CD3E CD2 TRAC CD4 CD8A
KRT5 KRT19 KRT7 EPCAM COL17A1 CLCA2 ODAM CSF3R FPR1 MPO S100A12 CAMP
CCL2 CCL3 CCL4 CCL5 CCL8 CCL18 CXCL9 CXCL10 CXCL11 PTGDS RBP5 EIF1AY CETP SKAP1 LAT CCR6 CD1D
FBLN7 SGPP2 SH2D1A TIGIT IGSF6 IL2RA CD38 CD5 GFI1 IL1R1 IL1R2 IL10 CCL20 TRAF6 STAT5A IFNG STAT1 IL6 CXCL8 NFKB1 RELA
TCL1A IGHD IGHM MEF2B""".split()
CAND=sorted(set(CAND))
# gene universes
g16134=set(pd.read_csv(os.path.join(ROOT,"data/external/GSE16134_expr_gene.tsv"),sep="\t",usecols=[0]).iloc[:,0].astype(str))
g173078=set(pd.read_csv(os.path.join(ROOT,"data/external/GSE173078_logCPM.tsv"),sep="\t",usecols=[0]).iloc[:,0].astype(str))
# xenium features from first complete section
xen=None
for d in sorted(os.listdir(os.path.join(ROOT,"data/raw/xenium"))):
    p=os.path.join(ROOT,"data/raw/xenium",d,"features.tsv.gz")
    if os.path.exists(p):
        f=pd.read_csv(p,sep="\t",header=None)
        xen=set(f.loc[f[2].astype(str).eq("Gene Expression"),1].astype(str)); break
# visium features
vis=set()
for f in os.listdir(os.path.join(ROOT,"data/external/GSE206621_raw")):
    if f.endswith("_features.tsv.gz"):
        v=pd.read_csv(os.path.join(ROOT,"data/external/GSE206621_raw",f),sep="\t",header=None)
        vis=set(v[1].astype(str)); break
rows=[]
for g in CAND:
    rows.append(dict(gene=g, GPL570_gene=g in g16134, GSE173078_detected=g in g173078,
                     Xenium_panel=(xen is not None and g in xen), Visium_detected=g in vis))
df=pd.DataFrame(rows)
df["n_platforms"]=df[["GPL570_gene","GSE173078_detected","Xenium_panel","Visium_detected"]].sum(1)
df.to_csv(os.path.join(META,"gene_coverage.tsv"),sep="\t",index=False)
print(df.sum(numeric_only=True))
print(df.head(10).to_string(index=False))
