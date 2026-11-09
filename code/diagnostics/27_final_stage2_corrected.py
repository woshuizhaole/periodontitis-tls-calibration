#!/usr/bin/env python3
"""Final corrected Stage 2: cell-level composition-matched null; no increment beyond B/T/APC."""
import os, pandas as pd, numpy as np
ROOT="/data/workplace/wyx/perio_tls_20261107"; X=os.path.join(ROOT,"results/stage2/xenium"); D=os.path.join(ROOT,"diagnostics/xenium/results"); OUT=os.path.join(X,"corrected"); os.makedirs(OUT,exist_ok=True)
# audit original matched control overlap
mp=pd.read_csv(os.path.join(X,"matched_pairs.tsv"),sep="\t"); cr=pd.read_csv(os.path.join(X,"control_regions.tsv"),sep="\t")
j=mp.merge(cr[["control_id","overlap_tls"]],left_on="control",right_on="control_id",how="left")
pd.DataFrame([dict(n_pairs=len(j),n_overlap_gt5=int((j.overlap_tls>0.05).sum()),n_overlap_gt20=int((j.overlap_tls>0.20).sum()),
                   median_overlap=float(j.overlap_tls.median()),mean_overlap=float(j.overlap_tls.mean()),max_overlap=float(j.overlap_tls.max()),
                   note="original density-matched controls contaminated; p=0.766 retracted")]).to_csv(os.path.join(OUT,"control_overlap_audit.tsv"),sep="\t",index=False)
# final cell-level composition-matched tests
c=pd.read_csv(os.path.join(D,"matched_cell_region_tests.tsv"),sep="\t")
rows=[]
for _,r in c.iterrows():
    rows.append(dict(test=f"matched_composition_{r.score}",n=int(r.n),mean_diff=float(r["mean"]),cohen_d=float(r.cohen_d),
                     p=float(r.p_wilcoxon),q=float(r.q_bh),analysis="cell_level_1to1_composition_matched"))
# region-level clean controls (composition-confounded, reported separately)
try:
    rt=pd.read_csv(os.path.join(D,"score_tests.tsv"),sep="\t")
    for _,r in rt[rt.score.isin(["OTLS_panel","lymphoid7","stromal7","hev3","BTAPC","B","T","APC"])].iterrows():
        rows.append(dict(test=f"region_clean_{r.score}",n=int(r.n),mean_diff=float(r.mean_diff),cohen_d=float(r.cohen_d),
                         p=float(r.p_wilcoxon),q=float(r.q_bh),analysis="clean_region_controls_but_composition_diff"))
except Exception: pass
# unchanged region-level disease/decay
main=pd.read_csv(os.path.join(X,"tests_v1_contaminated.tsv"),sep="\t")
for _,r in main.iterrows():
    if r.test.startswith("region_count") or r.test.startswith("OTLS_HvP") or r.test.startswith("decay"):
        rows.append(dict(test=r.test,n=int(r.n),mean_diff=float(r.mean_diff),cohen_d=float(r.cohen_d),p=float(r.p),q=np.nan,analysis="unchanged_region_level"))
cor=pd.DataFrame(rows); cor.to_csv(os.path.join(OUT,"corrected_tests.tsv"),sep="\t",index=False); cor.to_csv(os.path.join(X,"tests.tsv"),sep="\t",index=False)
print(pd.read_csv(os.path.join(OUT,"control_overlap_audit.tsv"),sep="\t").to_string(index=False))
print(cor[cor.test.str.startswith("matched_composition")].to_string(index=False))
