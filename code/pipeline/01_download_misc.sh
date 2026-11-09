#!/usr/bin/env bash
set -euo pipefail
P=/data/workplace/wyx/perio_tls_20261107
cd "$P"
mkdir -p data/raw/tonsil data/raw/phenocycler logs/downloads
# 10x tonsil Visium FFPE (direct)
if [ ! -s data/raw/tonsil/tonsil_expr.h5 ]; then
  wget -c -O data/raw/tonsil/tonsil_expr.h5 "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil_AddOns/CytAssist_FFPE_Protein_Expression_Human_Tonsil_AddOns_filtered_feature_bc_matrix.h5"
fi
if [ ! -s data/raw/tonsil/tonsil_spatial.tar.gz ]; then
  wget -c -O data/raw/tonsil/tonsil_spatial.tar.gz "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil_AddOns/CytAssist_FFPE_Protein_Expression_Human_Tonsil_AddOns_spatial.tar.gz"
fi
mkdir -p data/raw/tonsil/spatial
tar -xzf data/raw/tonsil/tonsil_spatial.tar.gz -C data/raw/tonsil/ || true
