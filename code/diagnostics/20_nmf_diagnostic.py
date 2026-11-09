#!/usr/bin/env python3
"""Diagnostic: NMF co-expression modules in GSE164241, test TLS association vs lineage/inflammation."""
import os, time, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
from sklearn.decomposition import NMF
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"diagnostics/nmf"); os.makedirs(OUT,exist_ok=True)
LOG=os.path.join(OUT,"nmf.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
log("load")
a=sc.read_h5ad(os.path.join(ROOT,"data/external/GSE164241_integrated.h5ad"))
meta=pd.read_csv(os.path.join(ROOT,"results/stage1/celltype_metadata.tsv"),sep="\t",index_col=0)
a=a[a.obs_names.isin(meta.index)].copy(); a.obs=meta.loc[a.obs_names]
log(f"{a.shape}")
# HVGs on log data
sc.pp.highly_variable_genes(a,n_top_genes=3000,batch_key="donor")
hv=a.var_names[a.var.highly_variable].tolist()
# pseudobulk cell_type x donor
rows=[]; groups=[]
X=a[:,hv].X
if sp.issparse(X): X=X.tocsr()
for (ct,donor),idx in a.obs.groupby(["cell_type","donor"]).groups.items():
    ii=a.obs.index.get_indexer(idx)
    if len(ii)<20: continue
    rows.append(np.asarray(X[ii].mean(0)).ravel()); groups.append((ct,donor))
P=np.vstack(rows).astype(np.float32); log(f"pseudobulk {P.shape}")
for k in [15,20,25]:
    nmf=NMF(n_components=k,init="nndsvda",random_state=42,max_iter=600,beta_loss="frobenius")
    W=nmf.fit_transform(P); H=nmf.components_
    np.save(os.path.join(OUT,f"W_k{k}.npy"),W); np.save(os.path.join(OUT,f"H_k{k}.npy"),H)
    pd.DataFrame(H,columns=hv).to_csv(os.path.join(OUT,f"H_k{k}.tsv"),sep="\t")
    log(f"k={k} reconstruction={nmf.reconstruction_err_:.3f}")
# use k=20 top genes
H=pd.read_csv(os.path.join(OUT,"H_k20.tsv"),sep="\t",index_col=0)
top={i:H.loc[i].sort_values(ascending=False).head(30).index.tolist() for i in H.index}
# module scores per cell (mean z of top genes)
genes=set(g for gs in top.values() for g in gs); genes=[g for g in genes if g in a.var_names]
M=a[:,genes].X
if sp.issparse(M): M=M.toarray()
M=np.asarray(M,dtype=np.float32); Mz=(M-M.mean(0))/(M.std(0)+1e-9); gidx={g:i for i,g in enumerate(genes)}
scores={}
for i,gs in top.items():
    idx=[gidx[g] for g in gs if g in gidx]
    scores[f"module_{i}"]=Mz[:,idx].mean(1)
S=pd.DataFrame(scores,index=a.obs_names); S.to_csv(os.path.join(OUT,"module_scores.tsv"),sep="\t")
# correlations with TLS/lineage/inflammation
meta2=meta.loc[a.obs_names]
targets=["B_TLS","T_TLS","Fib_TLS","Endo_TLS","Myeloid_TLS","OTLS_core","OTLS_panel","Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]
corr=pd.DataFrame({m:[S[m].corr(meta2[t],method="spearman") for t in targets] for m in S.columns},index=targets).T
corr["max_TLS"]=corr[["B_TLS","T_TLS","Fib_TLS","Endo_TLS","Myeloid_TLS"]].max(1)
corr["max_control"]=corr[["Bcell","Plasma","Tcell","12CK","Cabrita9","TLS39","Inflammation"]].abs().max(1)
corr["specificity"]=corr.max_TLS-corr.max_control
corr.sort_values("specificity",ascending=False).to_csv(os.path.join(OUT,"module_correlations.tsv"),sep="\t")
for i in H.index:
    gs=top[i]; log(f"module_{i}: "+",".join(gs[:15]))
log("DONE")
