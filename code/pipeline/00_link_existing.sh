#!/usr/bin/env bash
set -euo pipefail
P=/data/workplace/wyx/perio_tls_20261107
mkdir -p "$P"/data/external "$P"/data/raw/geo "$P"/data/raw/xenium "$P"/data/raw/ibex "$P"/data/raw/tonsil "$P"/data/raw/phenocycler
link() { src="$1"; dst="$2"; if [ -e "$src" ]; then ln -sfn "$src" "$dst"; echo "LINK $dst -> $src"; else echo "MISSING $src"; fi; }
link /data/workplace/wyx/periodontitis_field_defect/data_raw/GSE164241 "$P/data/external/GSE164241_raw"
link /data/workplace/wyx/oral_deconv_reference_20260707/data/interim/gse164241_gingiva_integrated.h5ad "$P/data/external/GSE164241_integrated.h5ad"
link /data/workplace/wyx/oral_deconv_reference_20260707/data/interim/gse152042_labeled.h5ad "$P/data/external/GSE152042_labeled.h5ad"
link /data/workplace/wyx/periodontitis_network_inversion_20260504/data/raw/GSE152042 "$P/data/external/GSE152042_raw"
link /data/workplace/wyx/exposome-to-endotype/data_raw/sc_st/GSE171213 "$P/data/external/GSE171213_raw"
link /data/workplace/wyx/periodontitis_network_inversion_20260504/data/raw/GSE266897 "$P/data/external/GSE266897_raw"
link /data/workplace/wyx/oral_deconv_reference_20260707/data/raw/GSE206621 "$P/data/external/GSE206621_raw"
link /data/workplace/wyx/oral_deconv_reference_20260707/data/raw/GSE262668 "$P/data/external/GSE262668_raw"
link /data/workplace/wyx/periodontitis_field_defect/data_processed/GSE16134_expr_gene.tsv "$P/data/external/GSE16134_expr_gene.tsv"
link /data/workplace/wyx/periodontitis_field_defect/metadata/GSE16134_meta.tsv "$P/data/external/GSE16134_meta.tsv"
link /data/workplace/wyx/periodontitis_field_defect/data_processed/GSE10334_expr_gene.tsv "$P/data/external/GSE10334_expr_gene.tsv"
link /data/workplace/wyx/periodontitis_field_defect/metadata/GSE10334_meta.tsv "$P/data/external/GSE10334_meta.tsv"
link /data/workplace/wyx/periodontitis_field_defect/data_processed/GSE173078_logCPM.tsv "$P/data/external/GSE173078_logCPM.tsv"
link /data/workplace/wyx/periodontitis_field_defect/metadata/GSE173078_meta.tsv "$P/data/external/GSE173078_meta.tsv"
link /data/workplace/wyx/periodontitis_shared_cell_state_liability_20260507/data/raw/cellxgene/periodontitis_atlas.h5ad "$P/data/external/periodontitis_atlas.h5ad"
