#!/usr/bin/env python3
"""Build Xenium section->donor map from Nat Immunol 2026 Supplementary Table 1 (MOESM2)."""
import pandas as pd, re
ROOT="/data/workplace/wyx/perio_tls_20261107"
x=pd.read_csv(f"{ROOT}/metadata/xenium_section_audit.tsv",sep="\t")
def donor(s):
    s=s.replace("-1",""); s=re.sub(r"[a-z]$","",s); return s
x["donor"]=x.section.map(donor); x["group"]=x.section.str[0]
# ST column counts: number of x and number of * (parsed from MOESM2 Human_subject_metadata)
tab={
"H1":(1,0),"H2":(2,0),"H3":(2,2),"H4":(3,0),"H5":(1,0),"H6":(1,0),"H7":(2,0),"H8":(1,0),"H9":(2,0),
"P1":(2,0),"P2":(4,2),"P3":(1,0),"P4":(2,0),"P5":(3,6),"P6":(1,0),"P7":(2,0),"P8":(1,0),"P9":(2,0),
"P10":(1,0),"P11":(1,0),"P12":(2,2),"P13":(3,3),"P14":(2,0),"P15":(1,3)}
x["paper_st_sections"]=x.donor.map(lambda d:tab.get(d,(0,0))[0])
x["paper_starred_sections"]=x.donor.map(lambda d:tab.get(d,(0,0))[1])
x["main_quant"]=(x.paper_starred_sections==0)
# deposit has 5 P9 sections; paper table lists 2; extra 3 are not in the paper quantitative set
x.loc[(x.donor=="P9") & (~x.section.isin(["P9a","P9b"])),"main_quant"]=False
def note(r):
    if r.paper_starred_sections>0:
        return "paper table: donor has starred/excluded sections; exact section-level star status not resolvable from deposit"
    if r.donor=="P9" and r.section not in ["P9a","P9b"]:
        return "paper table lists 2 P9 sections; this deposit section is extra and excluded from main"
    return "paper table: donor has no starred ST section"
x["note"]=x.apply(note,axis=1)
x[["section","donor","group","n_cells","paper_st_sections","paper_starred_sections","main_quant","note"]].to_csv(f"{ROOT}/metadata/xenium_donor_map.tsv",sep="\t",index=False)
print(x.groupby(["group","main_quant"]).agg(sections=("section","size"),donors=("donor","nunique")).to_string())
print("main sections",int(x.main_quant.sum()),"donors",x[x.main_quant].donor.nunique())
print("paper ST total",int(x.drop_duplicates("donor").paper_st_sections.sum()),"star total",int(x.drop_duplicates("donor").paper_starred_sections.sum()))
