#!/usr/bin/env python3
"""Correct Stage 2 matched-control analysis: audit contaminated controls, import strict-control tests."""
import os, pandas as pd, numpy as np
ROOT="/data/workplace/wyx/perio_tls_20261107"; X=os.path.join(ROOT,"results/stage2/xenium"); OUT=os.path.join(X,"corrected"); os.makedirs(OUT,exist_ok=True)
main=pd.read_csv(os.path.join(X,"tests.tsv"),sep="\t")
diag=pd.read_csv(os.path.join(ROOT,"diagnostics/xenium/results/score_tests.tsv"),sep="\t")
mp=pd.read_csv(os.path.join(X,"matched_pairs.tsv"),sep="\t")
cr=pd.read_csv(os.path.join(X,"control_regions.tsv"),sep="\t")
# audit overlap of the controls actually used by the original matched analysis
j=mp.merge(cr[["control_id","overlap_tls"]],left_on="control",right_on="control_id",how="left")
audit=pd.DataFrame([dict(n_pairs=len(j),n_overlap_gt5=int((j.overlap_tls>0.05).sum()),
                         median_overlap=float(j.overlap_tls.median()),mean_overlap=float(j.overlap_tls.mean()),
                         max_overlap=float(j.overlap_tls.max()),
                         note="original matched controls were contaminated by TLS-like cells")])
audit.to_csv(os.path.join(OUT,"control_overlap_audit.tsv"),sep="\t",index=False)
# corrected matched tests from strict-control diagnostic
want=["OTLS_panel","lymphoid7","stromal7","BTAPC","B","T","APC","OTLS_core_independent"]
rows=[]
for _,r in diag[diag.score.isin(want)].iterrows():
    rows.append(dict(test=f"matched_{r.score}",n=int(r.n),mean_diff=float(r.mean_diff),cohen_d=float(r.cohen_d),
                     p=float(r.p_wilcoxon),q=float(r.q_bh),analysis="strict_controls_overlap_lt5pct"))
# keep unaffected main tests
for _,r in main.iterrows():
    if r.test.startswith("region_count") or r.test.startswith("OTLS_HvP") or r.test.startswith("decay"):
        rows.append(dict(test=r.test,n=int(r.n),mean_diff=float(r.mean_diff),cohen_d=float(r.cohen_d),p=float(r.p),q=np.nan,analysis="unchanged_region_level"))
cor=pd.DataFrame(rows); cor.to_csv(os.path.join(OUT,"corrected_tests.tsv"),sep="\t",index=False)
# summary of decisive contrast
sel=cor[cor.test.isin(["matched_OTLS_panel","matched_lymphoid7","matched_stromal7","matched_BTAPC"])]
sel.to_csv(os.path.join(OUT,"decisive_contrasts.tsv"),sep="\t",index=False)
print(audit.to_string(index=False)); print(sel.to_string(index=False))
