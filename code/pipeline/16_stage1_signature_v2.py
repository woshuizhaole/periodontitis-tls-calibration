#!/usr/bin/env python3
"""Stage 1 signature v2: TLS-biased, oral-data-supported program.

Rationale: the first fully data-driven core was dominated by generic NF-kB/immediate-early
activation genes (NAMPT/KLF6/SRGN/SOCS3/JUNB/NFKBIA...), which is not a TLS-specific program.
This script builds a TLS-biased core from canonical TLS/HEV/GC genes, keeps only genes with
oral single-cell DE support in the expected lineage, excludes generic activation genes, and
fills the extended set with data-driven non-generic genes. No spatial/bulk outcome is used.
"""
import os, time, warnings, yaml
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
warnings.filterwarnings("ignore")
ROOT="/data/workplace/wyx/perio_tls_20261107"; OUT=os.path.join(ROOT,"results/stage1"); SIG=os.path.join(ROOT,"signatures")
LOG=os.path.join(ROOT,"logs/analysis/stage1_signature_v2.log")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); open(LOG,"a").write(line+"\n")
# canonical TLS genes with expected lineage
seed={
 "Fibroblast":["CXCL13","CCL19","CCL21","LTB","TNFSF13B","CXCL12","FDCSP","CR2","PTGDS","VCAN","FBN1"],
 "Endothelial":["PLVAP","MMRN1","ACKR1","SELE","ICAM1","VCAM1","VWF"],
 "B":["BCL6","TNFRSF13B","CD27","CD38","MME","BCL2A1"],
 "T":["CXCR5","PDCD1","ICOS","IL21","TOX","CD200","TIGIT","CTLA4"],
 "Myeloid":["LAMP3","CCR7","CD83","CD40","IRF4","CLEC9A"],
 "Other":["CCL20","IL6ST"],
}
generic=set("""NAMPT KLF6 SRGN SOCS3 JUNB JUN G0S2 BIRC3 NFKBIA CXCL2 CXCL3 SOD2 TNFAIP3 TNFAIP2 PTGS2
ZC3H12A BATF RELB NFKBIZ REL KDM6B MAFF NFKB2 NFKB1 PIM3 PHLDA1 DUSP4 DUSP1 JUND ATF3 PPP1R15B
ABL2 IRF1 ZFP36 GRAMD1A BHLHE40 PNRC1 MARCKS MCL1 GADD45B EGR1 TANK DNAJA1 WTAP TYMP NR4A2 TAPBP
SAT1 SGK1 TAP1 TNIP1 CEBPD CCNL1 IL6 CXCL8 CCL2 CCL5 RELA""".split())
state_markers=set("""CXCL13 CCL19 CCL21 LTB TNFSF13B PLVAP MMRN1 ACKR1 SELE ICAM1 BCL6 AICDA MEF2B MME CD38 TNFRSF13B CD27 CXCR5 PDCD1 ICOS IL21 TOX CD200 LAMP3 CCR7 CD83 CD40""".split())
cov=pd.read_csv(os.path.join(ROOT,"metadata/gene_coverage.tsv"),sep="\t").set_index("gene")
cand=pd.read_csv(os.path.join(OUT,"candidate_genes_cross_lineage.tsv"),sep="\t").set_index("gene")
# DE lookup
de={}
for lin in ["B","T","Fibroblast","Endothelial","Myeloid"]:
    d=pd.read_csv(os.path.join(OUT,"de_results",f"{lin}_high_vs_low.tsv"),sep="\t")
    de[lin]=d.set_index("gene")
def de_support(g,lin):
    if lin not in de or g not in de[lin].index: return None
    r=de[lin].loc[g]
    if isinstance(r,pd.DataFrame): r=r.iloc[0]
    return float(r.logFC), float(r.padj), float(r.get("pct_high",np.nan))
def platform_ok(g):
    if g not in cov.index: return False,False
    r=cov.loc[g]
    return bool(r.GPL570_gene), bool(r.GSE173078_detected)
core=[]; ev=[]
for lin,genes in seed.items():
    for g in genes:
        gpl,rna=platform_ok(g)
        if not (gpl and rna):  # core requires both bulk platforms
            continue
        if g in generic: continue
        s=de_support(g,lin) if lin!="Other" else None
        # canonical core accepts strong DE in expected lineage; BCL6 is accepted via T support
        ok=False
        if s and s[1]<0.05 and s[0]>0.5: ok=True
        if g=="BCL6":
            st=de_support("BCL6","T"); ok=bool(st and st[1]<0.05)
        if g in ("CCL20","IL6ST"):
            ok=any((de_support(g,l) is not None and de_support(g,l)[1]<0.05 and de_support(g,l)[0]>0.5) for l in ["Myeloid","Endothelial"])
        if ok and g not in core:
            core.append(g); ev.append((g,lin,s))
log(f"core seeded={len(core)}: {core}")
# fill to 30 with non-generic data-driven genes (both platforms)
if len(core)<30:
    extra=cand[(cand.index.isin(cov.index))&(~cand.index.isin(core))&(~cand.index.isin(generic))]
    keep=[]
    for g,r in extra.iterrows():
        gpl,rna=platform_ok(g)
        if gpl and rna and r["min_padj"]<1e-20 and r["mean_logFC"]>1.0:
            keep.append(g)
        if len(core)+len(keep)>=30: break
    core += keep[:30-len(core)]
    for g in keep[:30-len(core)]:
        r=cand.loc[g]; ev.append((g,"data-driven",(float(r.mean_logFC),float(r.min_padj),np.nan)))
log(f"core final={len(core)}: {core}")
# extended: core + canonical extras detected + top non-generic DE genes
canon_all=[]
for lin,gs in seed.items(): canon_all+=gs
ext=[]
for g in core+canon_all:
    if g not in ext and g in cov.index and (cov.loc[g].GPL570_gene or cov.loc[g].GSE173078_detected): ext.append(g)
extra=cand[(cand.index.isin(cov.index))&(~cand.index.isin(ext))&(~cand.index.isin(generic))].sort_values(["n_lineages","mean_logFC"],ascending=False)
for g,r in extra.iterrows():
    if len(ext)>=120: break
    if r["mean_logFC"]>0.8 and r["min_padj"]<1e-20: ext.append(g)
log(f"ext={len(ext)}")
# panel
xen_df=pd.read_csv(os.path.join(ROOT,"data/raw/xenium/H1/features.tsv.gz"),sep="\t",header=None)
xen=set(xen_df.loc[xen_df[2].astype(str).eq("Gene Expression"),1].astype(str))
panel=[g for g in core if g in xen]
if len(panel)<8: panel += [g for g in ext if g in xen and g not in panel][:8-len(panel)]
log(f"panel={len(panel)}: {panel}")
# save
def save(name,genes,source):
    with open(os.path.join(SIG,name+".yaml"),"w") as fh: yaml.safe_dump({"name":name,"version":"2.0","genome":"GRCh38","source":source,"genes":genes},fh,sort_keys=False)
save("OTLS_core",core,"TLS-biased canonical program with oral GSE164241 DE support; generic NF-kB/IE genes excluded")
save("OTLS_ext",ext,"OTLS_core + canonical TLS genes + non-generic oral DE genes")
save("OTLS_panel",panel,"OTLS_core intersect Xenium 451 panel")
core_independent=[g for g in core if g not in state_markers]
save("OTLS_core_independent",core_independent,"OTLS_core genes not used to define TLS high/low states; non-circular sensitivity")
with open(os.path.join(OUT,"signature_evidence_v2.tsv"),"w") as fh:
    fh.write("gene\tsignature\texpected_lineage\tsupport_lineage\tlogFC\tpadj\tpct_high\tcircular_state_marker\n")
    for g in core:
        best=None
        for lin in ["B","T","Fibroblast","Endothelial","Myeloid"]:
            s=de_support(g,lin)
            if s and s[1]<0.05 and s[0]>0.5:
                if best is None or s[1]<best[1][1]: best=(lin,s)
        if best:
            lin,s=best; fh.write(f"{g}\tcore\tNA\t{lin}\t{s[0]:.3f}\t{s[1]:.3g}\t{s[2]}\t{str(g in state_markers).lower()}\n")
        else:
            fh.write(f"{g}\tcore\tNA\tNA\tNA\tNA\tNA\t{str(g in state_markers).lower()}\n")
    for g in ext:
        if g in core: continue
        r=cand.loc[g] if g in cand.index else None
        fh.write(f"{g}\text\tNA\tdata-driven\t{r.mean_logFC if r is not None else 'NA'}\t{r.min_padj if r is not None else 'NA'}\tNA\t{str(g in state_markers).lower()}\n")
# score GSE164241
a=ad.read_h5ad(os.path.join(ROOT,"data/external/GSE164241_integrated.h5ad"))
for name,gs in [("OTLS_core",core),("OTLS_ext",ext),("OTLS_panel",panel),("OTLS_core_independent",core_independent)]:
    gs=[g for g in gs if g in a.var_names]
    if len(gs)>=2: sc.tl.score_genes(a,gs,score_name=name,ctrl_size=50,n_bins=25)
a.obs[["OTLS_core","OTLS_ext","OTLS_panel","OTLS_core_independent"]].to_csv(os.path.join(OUT,"OTLS_v2_cell_scores.tsv"),sep="\t")
log("DONE")
