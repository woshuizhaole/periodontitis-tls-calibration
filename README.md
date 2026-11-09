# periodontitis-tls-calibration

Analysis code for the study **TLS expression scores measure lymphoid infiltration in periodontitis**.

## Data

The analysis uses public bulk transcriptomes GSE16134, GSE10334, GSE173078 and GSE223328, a Xenium spatial dataset, full-transcriptome Visium data (GSE297924) and 40-plex protein images. Raw data are not redistributed here; accession-level metadata and processed score tables are in the manuscript and its supplementary tables.

## Code layout

- `code/pipeline/` — preprocessing, score construction, bulk models, Xenium/Visium processing and sensitivity analyses.
- `code/narrative_B/` — composition, bulk and spatial analyses used in the manuscript.
- `code/diagnostics/` — diagnostic and corrected-design analyses.

Scripts assume the project root layout described in the manuscript (or the path given by `PROJECT_ROOT`). 

## Environment

Python 3.11; principal packages: numpy, pandas, scipy, scikit-learn, statsmodels, anndata, scanpy, squidpy, tifffile. Install the versions available at the time of analysis from the manuscript's Methods section.

## Citation

Cui H, Liang H, Wang X, Jiao X, Wu L. TLS expression scores measure lymphoid infiltration in periodontitis. Manuscript submitted to *Clinical Oral Investigations*.
