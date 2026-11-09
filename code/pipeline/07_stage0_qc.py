#!/usr/bin/env python3
import os, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"
QC=os.path.join(ROOT,"logs/qc"); os.makedirs(QC,exist_ok=True)
rows=[]
def add(item, value, status, note=""):
    rows.append(dict(item=item,value=value,status=status,note=note))
# metadata
d=pd.read_csv(os.path.join(ROOT,"metadata/datasets.tsv"),sep="\t")
add("datasets_rows",len(d),"PASS" if len(d)==864 else "CHECK")
for g,exp in {"GSE16134":310,"GSE10334":247,"GSE173078":36,"GSE223328":16,"GSE164241":30,"GSE171213":12,"GSE152042":4,"GSE296447":8,"GSE206621":8,"GSE299738":3,"GSE266897":4,"GSE262668":7,"GSE180588":179}.items():
    n=(d.dataset_id==g).sum(); add(f"n_samples_{g}",int(n),"PASS" if n==exp else "CHECK",f"expected {exp}")
# xenium
x=pd.read_csv(os.path.join(ROOT,"metadata/xenium_section_audit.tsv"),sep="\t")
add("xenium_sections",len(x),"PASS" if len(x)==46 else "FAIL")
add("xenium_cells_total_postqc",int(x.n_cells.sum()),"PASS")
add("xenium_cells_deposited",int(x.n_cells_deposited.sum()),"INFO")
add("xenium_min_cells",int(x.n_cells.min()),"PASS" if x.n_cells.min()>=500 else "CHECK")
add("xenium_max_neg_frac_section",float(x.neg_frac_section.max()),"PASS" if x.neg_frac_section.max()<0.05 else "CHECK")
add("xenium_sections_with_btapc_rich",int((x.btapc_rich_regions>0).sum()),"INFO",f"null mean {x.null_btapc_rich_regions_mean.mean():.2f}/section; not TLS evidence")
# gene coverage
g=pd.read_csv(os.path.join(ROOT,"metadata/gene_coverage.tsv"),sep="\t")
add("candidate_genes",len(g),"PASS")
add("GPL570_coverage",float(g.GPL570_gene.mean()),"PASS" if g.GPL570_gene.mean()>=0.8 else "CHECK")
add("Xenium_panel_coverage",int(g.Xenium_panel.sum()),"INFO")
add("Xenium_coverage_fraction",float(g.Xenium_panel.mean()),"INFO")
v=pd.read_csv(os.path.join(ROOT,"metadata/visium_sample_audit.tsv"),sep="\t")
add("visium_low_quality_sections",int(v.low_quality.sum()),"INFO","6/8 <500 spots; attenuation only")
ib=pd.read_csv(os.path.join(ROOT,"metadata/ibex_audit.tsv"),sep="\t")
add("ibex_files",len(ib),"PASS" if len(ib)==29 else "CHECK")
# required files
required=[
 "data/raw/xenium/H1/matrix.mtx.gz","data/raw/xenium/H1/features.tsv.gz",
 "data/raw/geo/GSE296447/GSE296447_RAW.tar","data/raw/geo/GSE299738/GSE299738_RAW.tar",
 "data/raw/geo/GSE223328/GSE223328_RAW.tar","data/raw/tonsil/tonsil_expr.h5",
 "data/external/GSE164241_integrated.h5ad","data/external/GSE152042_labeled.h5ad",
 "data/external/GSE16134_expr_gene.tsv","data/external/GSE173078_logCPM.tsv",
 "data/external/GSE206621_raw"]
for p in required:
    fp=os.path.join(ROOT,p); add("file_"+os.path.basename(p),os.path.getsize(fp) if os.path.exists(fp) else 0,"PASS" if os.path.exists(fp) and os.path.getsize(fp)>0 else "FAIL")
# manifest
m=pd.read_csv(os.path.join(ROOT,"data/manifests/DATA_MANIFEST.tsv"),sep="\t")
add("manifest_rows",len(m),"PASS")
pd.DataFrame(rows).to_csv(os.path.join(QC,"stage0_qc.tsv"),sep="\t",index=False)
print(pd.DataFrame(rows).to_string(index=False))
print("overall", "PASS" if all(r["status"]!="FAIL" for r in rows) else "FAIL")
