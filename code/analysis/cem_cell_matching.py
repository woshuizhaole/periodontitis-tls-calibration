import pandas as pd, numpy as np, glob, os
files=sorted(glob.glob("results/stage2/xenium/cells/*.tsv.gz"))
frames=[]
for f in files:
    sec=os.path.basename(f).replace(".tsv.gz","")
    d=pd.read_csv(f,sep="\t",usecols=["cell_id","B","Plasma","T","APC","OTLS_panel","BTAPC"])
    d["section"]=sec
    frames.append(d)
d=pd.concat(frames,ignore_index=True)
for c in ["B","Plasma","T","APC","OTLS_panel","BTAPC"]: d[c]=pd.to_numeric(d[c],errors="coerce")
t=pd.read_csv("results/stage2/xenium/cells_tls.tsv.gz",sep="\t",usecols=["section","cell_id","tls_like"])
d=d.merge(t,on=["section","cell_id"],how="left")
d["tls_like"]=d.tls_like.fillna(False)
for c in ["B","T","APC"]:
    d[c+"_bin"]=d.groupby("section")[c].transform(lambda x: pd.qcut(x,3,labels=False,duplicates="drop"))
d=d.dropna(subset=["B_bin","T_bin","APC_bin"])
rows=[]; rng=np.random.default_rng(11)
for (sec,b,t_,a),g in d.groupby(["section","B_bin","T_bin","APC_bin"]):
    tr=g[g.tls_like]; co=g[~g.tls_like]
    if len(tr)==0 or len(co)==0: continue
    n=min(len(tr),len(co))
    tr=tr.sample(n=n,random_state=int(rng.integers(1e9))); co=co.sample(n=n,random_state=int(rng.integers(1e9)))
    rows.append(pd.DataFrame({"section":sec,"B_bin":b,"T_bin":t_,"APC_bin":a,
        "d_OTLS_panel":tr["OTLS_panel"].to_numpy(float)-co["OTLS_panel"].to_numpy(float),
        "d_BTAPC":tr["BTAPC"].to_numpy(float)-co["BTAPC"].to_numpy(float),
        "d_B":tr["B"].to_numpy(float)-co["B"].to_numpy(float),
        "d_T":tr["T"].to_numpy(float)-co["T"].to_numpy(float),
        "d_APC":tr["APC"].to_numpy(float)-co["APC"].to_numpy(float)}))
r=pd.concat(rows,ignore_index=True)
print("cells",len(r),"strata",len(rows),"sections",r.section.nunique())
for c in ["d_OTLS_panel","d_BTAPC","d_B","d_T","d_APC"]:
    x=r[c].to_numpy(float); bs=[np.mean(rng.choice(x,len(x),replace=True)) for _ in range(1000)]
    print(c,len(x),round(float(x.mean()),4),round(float(x.mean()/(x.std(ddof=1)+1e-12)),4),[round(float(v),4) for v in np.percentile(bs,[2.5,97.5])])
r.to_csv("diagnostics/xenium/results/cem_cell_matched.tsv",sep="\t",index=False)
