#!/usr/bin/env python3
"""Stage 0: download GEO MINiML metadata, build sample/donor tables and feasibility report."""
import gzip
import os
import re
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
from datetime import datetime

import pandas as pd

ROOT = "/data/workplace/wyx/perio_tls_20261107"
META = os.path.join(ROOT, "metadata")
RAW = os.path.join(ROOT, "data/raw/geo_metadata")
os.makedirs(META, exist_ok=True)
os.makedirs(RAW, exist_ok=True)

DATASETS = {
    "GSE164241": dict(path="GSE164nnn", role="single-cell discovery", platform="10x scRNA-seq"),
    "GSE171213": dict(path="GSE171nnn", role="single-cell discovery", platform="10x scRNA-seq"),
    "GSE152042": dict(path="GSE152nnn", role="single-cell holdout", platform="10x scRNA-seq"),
    "GSE296447": dict(path="GSE296nnn", role="CITE-seq validation", platform="10x CITE-seq"),
    "GSE299738": dict(path="GSE299nnn", role="single-cell single-arm", platform="10x scRNA-seq"),
    "GSE266897": dict(path="GSE266nnn", role="single-cell single-arm", platform="10x scRNA-seq"),
    "GSE262668": dict(path="GSE262nnn", role="marker-only", platform="scRNA marker CSV"),
    "GSE206621": dict(path="GSE206nnn", role="Visium attenuation", platform="10x Visium"),
    "GSE16134": dict(path="GSE16nnn", role="bulk main", platform="Affymetrix HG-U133 Plus 2.0"),
    "GSE10334": dict(path="GSE10nnn", role="bulk technical replicate", platform="Affymetrix HG-U133 Plus 2.0"),
    "GSE173078": dict(path="GSE173nnn", role="bulk ordinal directional", platform="RNA-seq"),
    "GSE223328": dict(path="GSE223nnn", role="bulk descriptive replication", platform="Agilent GPL21185"),
    "GSE180588": dict(path="GSE180nnn", role="optional NHP dynamics", platform="microarray"),
}


def fetch_miniml(gse, path):
    xml = os.path.join(RAW, f"{gse}_family.xml")
    if os.path.exists(xml) and os.path.getsize(xml) > 1000:
        return xml
    tgz = os.path.join(RAW, f"{gse}_family.xml.tgz")
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{path}/{gse}/miniml/{gse}_family.xml.tgz"
    if not os.path.exists(tgz) or os.path.getsize(tgz) < 1000:
        subprocess.run(["curl", "-L", "--retry", "3", "-m", "180", "-o", tgz, url], check=True)
    with tarfile.open(tgz, "r:gz") as tf:
        tf.extractall(RAW)
    return xml


def strip_ns(tag):
    return tag.split("}")[-1]


def text(el):
    return "".join(el.itertext()).strip() if el is not None else ""


def parse_family(gse, xml):
    tree = ET.parse(xml)
    root = tree.getroot()
    samples = []
    for smp in root.iter():
        if strip_ns(smp.tag) != "Sample":
            continue
        title = ""
        chars = {}
        gsm = ""
        platform = ""
        for child in smp.iter():
            tag = strip_ns(child.tag)
            if tag == "Title" and not title:
                title = text(child)
            elif tag == "Accession" and not gsm:
                gsm = text(child)
            elif tag == "Characteristics":
                chars[child.attrib.get("tag", "")] = text(child)
            elif tag == "Platform-Ref" and not platform:
                platform = child.attrib.get("ref", "")
        samples.append(dict(gsm=gsm, title=title, platform=platform, chars=chars))
    return samples


def derive(gse, s):
    chars = {k.lower(): v for k, v in s["chars"].items()}
    title = s["title"]
    low = title.lower()
    d = dict(
        dataset_id=gse,
        accession=gse,
        platform=DATASETS[gse]["platform"],
        sample_id=s["gsm"] or title,
        subject_id="",
        tissue=chars.get("tissue", chars.get("tissue type", "")),
        disease_state=chars.get("disease state", chars.get("phenotype", chars.get("subject status", ""))),
        severity_label="",
        site_status="",
        age=chars.get("age", ""),
        sex=chars.get("sex", chars.get("gender", "")),
        batch=chars.get("batch", ""),
        file_path="",
        file_size="",
        checksum="",
        notes=title,
    )
    if gse in ("GSE16134", "GSE10334"):
        m = re.search(r"patient\s*(\d+)", title, re.I)
        d["subject_id"] = f"P{m.group(1)}" if m else ""
        d["tissue"] = "gingiva"
        d["site_status"] = "unaffected" if "unaffected" in low else ("affected" if "affected" in low else "")
        d["disease_state"] = "healthy_site" if d["site_status"] == "unaffected" else "diseased_site"
    elif gse == "GSE164241":
        base = title
        base = re.sub(r"(b|c)$", "", base)
        d["subject_id"] = base
        d["tissue"] = chars.get("tissue", "").lower()
        d["disease_state"] = chars.get("disease state", "").lower()
    elif gse == "GSE171213":
        d["subject_id"] = title.split("_")[0]
        grp = chars.get("group", "")
        d["disease_state"] = grp
        d["severity_label"] = "severe_chronic" if "Severe" in grp else ("treated" if "Treatment" in grp else "healthy")
    elif gse == "GSE152042":
        d["subject_id"] = title
        d["disease_state"] = chars.get("disease state", "")
        d["tissue"] = "gingiva"
    elif gse == "GSE296447":
        m = re.search(r"(P\d+)", title)
        d["subject_id"] = m.group(1) if m else title
        d["tissue"] = "gingiva"
        d["disease_state"] = "periodontitis"
        d["notes"] = title
    elif gse == "GSE206621":
        m = re.search(r"\[([A-D][12])\]", title)
        d["subject_id"] = m.group(1) if m else title
        d["tissue"] = "oral_mucosa"
        d["disease_state"] = chars.get("subject status", "")
    elif gse == "GSE173078":
        d["subject_id"] = title.split("_")[0]
        d["tissue"] = "gingiva"
        d["disease_state"] = chars.get("phenotype", "")
        d["severity_label"] = chars.get("phenotype", "")
    elif gse == "GSE223328":
        d["subject_id"] = title.split()[0]
        d["tissue"] = chars.get("tissue", "")
        d["disease_state"] = "periodontitis" if "Periodontitis" in chars.get("tissue", "") else "healthy"
        d["notes"] = f"EBV={chars.get('ebv status','')}"
    elif gse == "GSE180588":
        m = re.search(r"([A-Za-z]+ NHP \d+)", title)
        d["subject_id"] = m.group(1) if m else title
        d["tissue"] = "gingiva"
        d["disease_state"] = chars.get("time point", "")
    elif gse in ("GSE299738", "GSE266897", "GSE262668"):
        d["subject_id"] = title
        d["tissue"] = chars.get("tissue", "gingiva")
        d["disease_state"] = chars.get("subject status", chars.get("group", ""))
    return d


def main():
    all_rows = []
    for gse, cfg in DATASETS.items():
        try:
            xml = fetch_miniml(gse, cfg["path"])
            samples = parse_family(gse, xml)
            rows = [derive(gse, s) for s in samples]
            all_rows.extend(rows)
            print(f"{gse}: {len(rows)} samples")
        except Exception as e:
            print(f"ERROR {gse}: {e}", file=sys.stderr)
    df = pd.DataFrame(all_rows)
    df["role"] = df["dataset_id"].map({k: v["role"] for k, v in DATASETS.items()})
    df.to_csv(os.path.join(META, "datasets.tsv"), sep="\t", index=False)

    # donor map
    dm = df[["dataset_id", "sample_id", "subject_id", "disease_state", "severity_label", "site_status", "tissue"]].copy()
    dm["notes"] = ""
    dm.to_csv(os.path.join(META, "donor_map.tsv"), sep="\t", index=False)

    # discrepancies
    disc = [
        dict(dataset="GSE206621", issue="GEO has 5 health (incl. technical replicates) + 3 disease sections; paper reports 4 disease samples",
             resolution="Use GEO-deposited 3 disease sections; effective biological n ~3 health + 3 disease; attenuation only"),
        dict(dataset="GSE10334", issue="242/247 samples overlap GSE16134 by patient+sample+site; expression r median 0.9988",
             resolution="Treat as technical replicate; not independent; excluded from meta-analysis"),
        dict(dataset="GSE296447", issue="4 periodontitis patients, GEX+ADT = 8 libraries; not 8 patients",
             resolution="Use as CITE-seq validation; donor-level analysis"),
        dict(dataset="GSE173078", issue="Only 36 samples (12/12/12); ordinal severity directional, not main test",
             resolution="Directional ordinal model; zero result reported as inconclusive"),
        dict(dataset="GSE262668", issue="Only marker CSV, no count matrix",
             resolution="Marker direction check only; not integrated"),
        dict(dataset="Xenium", issue="Deposit has 46 unique sections (44 standard + P1-1 and P11-1), 24 donors (9H/15P); paper ST analysis 19 donors (8H/11P)",
             resolution="Main analysis uses 28 sections/18 donors if mapping available; 46 sections/24 donors for sensitivity"),
        dict(dataset="IBEX", issue="No CD21/CD23/CD35/BCL6; 29 images, deposited 9H+9P, analyzed 8H+8P",
             resolution="Use HEV/Ki67/CD138/CD31 organization proxy; no histologic maturity claim"),
    ]
    pd.DataFrame(disc).to_csv(os.path.join(META, "discrepancies.tsv"), sep="\t", index=False)

    # feasibility report
    lines = ["# Stage 0 feasibility report", "", f"Generated: {datetime.now().isoformat(timespec='seconds')}", ""]
    lines.append("## Sample counts by dataset")
    for gse, sub in df.groupby("dataset_id"):
        lines.append(f"- {gse}: {len(sub)} samples; disease states: {sub['disease_state'].value_counts().to_dict()}")
    lines.append("")
    lines.append("## Data roles")
    for gse, cfg in DATASETS.items():
        lines.append(f"- {gse}: {cfg['role']} ({cfg['platform']})")
    lines.append("")
    lines.append("## Decisions")
    lines.append("- Xenium section/donor audit is in metadata/xenium_section_audit.tsv; donor mapping from Supplementary Table 1 where available.")
    lines.append("- GSE10334 is a technical replicate of GSE16134; not independent.")
    lines.append("- GSE173078 is directional only; primary disease-status test is GSE16134 69 complete pairs.")
    lines.append("- IBEX organization proxy uses Ki67/CD138/CD31; no FDC/GC maturity claim.")
    lines.append("- All data are open; no controlled-access data.")
    with open(os.path.join(META, "feasibility_report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # analysis plan
    plan = """# Analysis plan v1

## Primary (confirmatory)
1. Xenium: OTLS_panel enrichment in B/T/APC density-matched TLS vs non-TLS regions (section/donor random effects; section-level bootstrap/permutation).
2. Xenium: distance decay of OTLS vs BTAPC (torus/block permutation).
3. Xenium: incremental value over BTAPC by section-level cross-validated AUC/effect size.
4. Xenium: organization proxy (HEV/Ki67/CD138/CD31) including organized vs unorganized.
5. Xenium: health vs disease (region count per mm2, area fraction, OTLS level).
6. Bulk: GSE16134 69 complete affected/unaffected pairs (patient-level paired test; LMM as supporting).

## Secondary / directional
- GSE173078 ordinal severity (n=36; effect size/CI; inconclusive if null).
- GSE223328 descriptive binary replication (Firth; LOOCV AUC; OTLS x EBV exploratory).
- Design-stratified forest; no cross-design pooled estimate.
- 12CK, Cabrita9, TLS39 comparisons with pre-specified overlap sensitivity.

## Exploratory
- OLP/NHP; bacterial burden; EBV interaction; TLS subtypes.

## Statistical policy
- BH within each pre-specified family.
- Report effect sizes and 95% CI; no arbitrary success threshold.
- Analysis plan can be revised with reasons in decision_log.md.
"""
    with open(os.path.join(META, "analysis_plan_v1.md"), "w") as fh:
        fh.write(plan)

    # decision log
    dlog = """# Decision log

## 2026-11-07 Stage 0
- Reused processed GSE164241/GSE152042 h5ad and bulk GSE16134/GSE10334/GSE173078 matrices from prior in-house projects after checking dimensions and sample labels; raw GEO files remain available for audit.
- GSE10334 treated as technical replicate (242/247 overlap, r median 0.9988).
- GSE206621 used as attenuation/negative-control Visium cohort, not primary TLS ground truth.
- Xenium deposit has 46 unique sections; primary quantitative set 28 sections/18 donors if donor map available, otherwise section-level with donor sensitivity.
- IBEX lacks FDC/GC markers; organization proxy only.
"""
    with open(os.path.join(META, "decision_log.md"), "w") as fh:
        fh.write(dlog)

    print("wrote metadata tables")


if __name__ == "__main__":
    main()
