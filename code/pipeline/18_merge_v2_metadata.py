#!/usr/bin/env python3
"""Merge v2 OTLS scores into cell metadata (provenance for the four OTLS columns)."""
import os, pandas as pd
ROOT="/data/workplace/wyx/perio_tls_20261107"
m=pd.read_csv(os.path.join(ROOT,"results/stage1/celltype_metadata.tsv"),sep="\t",index_col=0)
s=pd.read_csv(os.path.join(ROOT,"results/stage1/OTLS_v2_cell_scores.tsv"),sep="\t",index_col=0)
m=m.drop(columns=[c for c in ["OTLS_core","OTLS_ext","OTLS_panel","OTLS_core_independent"] if c in m.columns])
m=m.join(s)
m.to_csv(os.path.join(ROOT,"results/stage1/celltype_metadata.tsv"),sep="\t")
print(m.shape, [c for c in m.columns if c.startswith("OTLS")])
