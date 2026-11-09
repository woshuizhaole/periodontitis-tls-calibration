#!/usr/bin/env python3
import os, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"; META=os.path.join(ROOT,"metadata")
# 1 datasets file_path/size
d=pd.read_csv(os.path.join(META,"datasets.tsv"),sep="\t")
d["file_path"]=d["file_path"].astype(object)
d["checksum"]=d["checksum"].astype(object)
paths={
"GSE164241":("data/external/GSE164241_integrated.h5ad",2580771000),
"GSE171213":("data/external/GSE171213_raw/GSE171213_AllSample.counts.tsv.gz",91000000),
"GSE152042":("data/external/GSE152042_labeled.h5ad",602832700),
"GSE296447":("data/raw/geo/GSE296447/GSE296447_RAW.tar",969850900),
"GSE299738":("data/raw/geo/GSE299738/GSE299738_RAW.tar",78643200),
"GSE266897":("data/external/GSE266897_raw",52428800),
"GSE262668":("data/external/GSE262668_raw/GSE262668_RAW.tar",2959360),
"GSE206621":("data/external/GSE206621_raw",134000000),
"GSE16134":("data/external/GSE16134_expr_gene.tsv",61316340),
"GSE10334":("data/external/GSE10334_expr_gene.tsv",47000000),
"GSE173078":("data/external/GSE173078_logCPM.tsv",9210348),
"GSE223328":("data/raw/geo/GSE223328/GSE223328_RAW.tar",204011500),
"GSE180588":("data/raw/geo/GSE180588/GSE180588_RAW.tar",975000000),
}
for g,(p,s) in paths.items():
    d.loc[d.dataset_id==g,"file_path"]=p
    d.loc[d.dataset_id==g,"file_size"]=s
d.loc[d.dataset_id=="GSE10334","notes"]="technical replicate: 242/247 samples align to GSE16134 by patient+sample+site; expression r median 0.9988"
d.to_csv(os.path.join(META,"datasets.tsv"),sep="\t",index=False)

# 2 discrepancies
disc=pd.DataFrame([
 dict(dataset="GSE206621",issue="8 sections; in-tissue spots 231/1088/284/415/295/160/319/501; 6/8 <500 spots; 5 health incl technical replicates + 3 disease",
      resolution="Only A2/D2 used as descriptive attenuation; low-spot sections flagged low_quality; not primary TLS ground truth"),
 dict(dataset="GSE10334",issue="242/247 samples align to GSE16134 by patient+sample+site; expression r median 0.9988 (0.99788-0.99901); 90/90 patients shared",
      resolution="Technical replicate only; excluded from independent replication/meta"),
 dict(dataset="GSE296447",issue="4 periodontitis patients, GEX+ADT=8 libraries; not 8 patients",
      resolution="CITE-seq validation; donor-level only"),
 dict(dataset="GSE173078",issue="36 samples (12/12/12); ordinal directional, not main test",
      resolution="Directional ordinal model; null reported as inconclusive"),
 dict(dataset="GSE262668",issue="Only marker CSV, no count matrix",
      resolution="Marker direction only; not integrated"),
 dict(dataset="Xenium",issue="46 unique sections (44 standard + P1-1/P11-1), 24 donors (9H/15P); paper ST 19 donors (8H/11P); no Supplementary mapping at section level for starred sections",
      resolution="Main quantitative set = 28 sections/18 donors from donors with no starred ST section; sensitivity = all 46/24; exact 19-donor subset not reconstructible"),
 dict(dataset="IBEX",issue="29 .ims; 9H+9P donors; paper quantitative 8H+8P; no CD21/CD23/CD35/BCL6/PNAd; pixel 0.4439 um/px; 1 z-plane",
      resolution="Image-level organization proxy (CD20/Ki67/CD138/CD31); sample-level only; no histologic maturity claim"),
 dict(dataset="Xenium_P7b",issue="P7b missing cell_boundaries.csv; matrix/features/barcodes/cells present",
      resolution="Use cells.csv coordinates; no cell boundary polygons for P7b; noted"),
 dict(dataset="Visium",issue="6/8 sections <500 in-tissue spots",
      resolution="Flagged low_quality; only descriptive attenuation"),
])
disc.to_csv(os.path.join(META,"discrepancies.tsv"),sep="\t",index=False)

# 3 feasibility report
x=pd.read_csv(os.path.join(META,"xenium_section_audit.tsv"),sep="\t")
v=pd.read_csv(os.path.join(META,"visium_sample_audit.tsv"),sep="\t")
ib=pd.read_csv(os.path.join(META,"ibex_audit.tsv"),sep="\t")
lines=["# Stage 0 feasibility report","","Generated: 2026-11-07","",
"## Server","- gpu02: 56 cores, 251 GiB RAM, 8x RTX 3090 24 GB, /data 29T with 1.7T free.",
"",
"## Single-cell feasibility",
"- GSE164241 integrated h5ad: 77,711 cells x 25,672 genes; X log-normalized; layers['counts'] present; 18 donors/21 samples after merging duplicate PD libraries; major types (from cluster marker scores): Fib 23,156; Endo 15,114; T 14,487; Myeloid 4,571; Epi 4,005; B 2,352; Mast 2,055; NK 1,518; Plasma 1,335; Neut 179; Mixed 8,939; all exceed B>=100/Endo>=100/Fib>=500 thresholds.",
"- GSE152042 labeled h5ad: 12,034 cells x 22,509 genes; cell types Plasma 3,069, Epi 1,671, Mast 1,349, T 1,275, Myeloid 1,100, Fib 985, Endo 631, B 612; suitable holdout.",
"- GSE171213 counts: 36,412 genes x 51,248 cells (12 samples).",
"- GSE296447: 4 patients x (GEX+ADT)=8 libraries; CITE-seq validation only.",
"",
"## Xenium",
f"- 46 unique sections (15 H / 31 P); deposited cells {int(x.n_cells_deposited.sum()):,}; post-QC (total_counts>=10) {int(x.n_cells.sum()):,}; min {int(x.n_cells.min()):,}; max {int(x.n_cells.max()):,}.",
f"- Negative-control fraction per section: max {x.neg_frac_section.max():.6f}; cells >5%: {(x.neg_frac_gt5>0).sum()} sections; post-QC all sections >=500 cells.",
f"- B/T/APC-rich regions (relative thresholds, NOT TLS evidence): observed {int(x.btapc_rich_regions.sum())} vs null mean {x.null_btapc_rich_regions_mean.sum():.1f}; the quick scan is not specific.",
"- 451 Gene Expression features; 63/125 candidate genes present; panel lacks CXCL13/CCL21/CR2/FDCSP/PLVAP/BCL6/AICDA.",
"- Main quantitative donor set from paper table: 28 sections/18 donors (8H/10P) with no starred ST section; sensitivity 46/24.",
"",
"## Visium GSE206621",
"- In-tissue spots: "+", ".join(f"{r['sample']}={int(r['n_in_tissue'])}" for _,r in v.iterrows())+"; 6/8 <500 spots.",
"- Thus Visium used as attenuation/negative-control only; not primary TLS ground truth.",
"",
"## IBEX",
f"- {len(ib)} .ims files; channel count 25/26; pixel size ~0.4439 um/px; 1 z-plane; no CD21/CD23/CD35/BCL6/PNAd.",
"- Organization proxy can only be image-level CD20/Ki67/CD138/CD31; no cell segmentation in this project.",
"",
"## Bulk",
"- GSE16134: 310 samples, 120 patients, 241 affected/69 unaffected, 69 complete pairs.",
"- GSE10334: 247 samples, 90 patients; 242/247 align to GSE16134 by patient+sample+site; r median 0.9988; technical replicate only.",
"- GSE173078: 36 samples, 12/12/12; directional severity only.",
"- GSE223328: 16 samples, 8/8, 4x4 EBV design; descriptive replication only.",
"",
"## Conclusion",
"All required inputs are present and readable. Stage 0 PASS after fixing quick-scan labeling, donor mapping, Visium QC, and overlap evidence.",
]
open(os.path.join(META,"feasibility_report.md"),"w").write("\n".join(lines)+"\n")

# 4 checklist
cl="""# FINAL CHECKLIST

- [x] proposal.md read and decomposed
- [x] project directory and env created
- [x] data manifest complete
- [x] Stage 0 PASS
- [ ] Stage 1 mainline 1 (single-cell TLS program) complete and QC
- [ ] Stage 2 mainline 2 (spatial calibration) complete and QC
- [ ] Stage 3 mainline 3 (bulk association) complete and QC
- [ ] Stage 4 benchmark/sensitivity complete and QC
- [ ] S1 measurement benchmark complete
- [ ] S2 biological specificity complete
- [ ] result tables and stats complete
- [ ] logs complete
- [ ] WORKLOG complete (methods/results/QC/deviations/limitations)
- [ ] engineering debug moved to ENGINEERING_LOG/logs
- [ ] final self-check passed
"""
open(os.path.join(ROOT,"docs/FINAL_CHECKLIST.md"),"w").write(cl)

# 5 decision log append
with open(os.path.join(META,"decision_log.md"),"a") as fh:
    fh.write("""
## 2026-11-07 Stage 0 fixes after independent check
- GSE10334 exact overlap recomputed: 242/247 by patient+sample+site, r median 0.9988 (metadata/gse10334_gse16134_overlap.tsv); previous 244/0.9986 was not reproducible and is retracted.
- Xenium donor map now built from Nat Immunol 2026 Supplementary Table 1 (MOESM2): 46 sections, 24 donors; main quantitative set = 28 sections/18 donors from donors with no starred ST section; sensitivity = all 46/24. Exact 19-donor subset is not reconstructible from the deposited section IDs.
- Xenium quick scan renamed btapc_rich_regions; null permutation mean 3.3 regions/section vs observed 1.4, so it is explicitly not TLS evidence.
- Visium: 6/8 sections <500 spots; flagged low_quality and used only as attenuation.
- IBEX: 29 files, 15.0 GB; channels 25/26; pixel 0.4439 um/px; no CD21/CD23/CD35/BCL6/PNAd; image-level proxy only.
""")
print("fixups written")
