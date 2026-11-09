#!/usr/bin/env bash
set -euo pipefail
P=/data/workplace/wyx/perio_tls_20261107
cd "$P"
mkdir -p data/raw/geo/{GSE296447,GSE299738,GSE223328,GSE180588} logs/downloads
dl() { url="$1"; out="$2"; if [ -s "$out" ]; then echo "EXISTS $out"; return; fi; echo "DOWNLOAD $url"; wget -c --tries=5 --timeout=60 -O "$out" "$url"; }
dl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296447/suppl/GSE296447_RAW.tar data/raw/geo/GSE296447/GSE296447_RAW.tar
dl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299738/suppl/GSE299738_RAW.tar data/raw/geo/GSE299738/GSE299738_RAW.tar
dl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223328/suppl/GSE223328_RAW.tar data/raw/geo/GSE223328/GSE223328_RAW.tar
# optional NHP
dl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180588/suppl/GSE180588_RAW.tar data/raw/geo/GSE180588/GSE180588_RAW.tar || true
