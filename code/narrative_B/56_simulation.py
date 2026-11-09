#!/usr/bin/env python3
"""P0.2 simple simulation: composition vs independent TLS program."""
import os, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/narrative_B/simulation"); os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(42); n=3000
# composition (B,T,APC) positive, normalized
comp=rng.dirichlet([1.2,1.2,1.2],size=n); b,t,a=comp.T
# genes: 30 OTLS; half composition-driven, half TLS-program-driven
G=30; X=np.zeros((n,G))
for j in range(G):
    if j<15:
        X[:,j]=0.6*b+0.6*t+0.6*a+rng.normal(0,0.5,n)
    else:
        X[:,j]=0.2*(b+t+a)+rng.normal(0,0.5,n)
# independent TLS program (stromal/HEV) exists in half the samples, orthogonal to composition
program=rng.binomial(1,0.5,size=n)
for j in range(15,G):
    X[:,j]+=0.8*program
# observed scores
OTLS=X.mean(1); BPTA=comp.mean(1)
# composition-only vs program scenarios: evaluate against program label
rows=[]
for label,score in [("OTLS",OTLS),("BPTA",BPTA)]:
    # residualize on composition
    resid=score-LinearRegression().fit(comp,score).predict(comp)
    rows.append(dict(score=label,AUC_program=float(roc_auc_score(program,score)),AUC_resid=float(roc_auc_score(program,resid)),rho_comp=float(np.corrcoef(score,BPTA)[0,1])))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"simulation_scores.tsv"),sep="\t",index=False)
print(pd.DataFrame(rows).to_string(index=False))
