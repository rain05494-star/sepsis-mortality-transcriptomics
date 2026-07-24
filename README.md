# Sepsis mortality whole-blood transcriptomic signature

This repository contains the analysis code accompanying the manuscript:

**Two-cohort discovery and external assessment of a mortality-associated whole-blood transcriptomic signature in sepsis.**

The study identifies an operationally defined whole-blood transcriptomic signature associated with short-term sepsis mortality using two discovery cohorts and evaluates its directional reproducibility, biological context, overlap with prior cross-cohort evidence, and descriptive MARS endotype-related variation in an external cohort.
## Repository contents
```text
sepsis-mortality-transcriptomics/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── code/
│   └── sepsis_mortality_transcriptomics_analysis.R
├── data/
│   └── README.md
├── metadata/
│   ├── GEO_accessions.txt
│   └── Giannini_E3_gene_set.csv
├── environment/
│   ├── sessionInfo.txt
│   └── package_versions.csv
└── docs/
    └── analysis_workflow.md
```
Only the analysis code, metadata needed to reproduce the analyses, and environment information are included in this repository. Public GEO expression data are not redistributed here and should be downloaded directly from GEO.
## Study overview

The analysis uses public whole-blood transcriptomic datasets from the Gene Expression Omnibus (GEO):

- **GSE272769**: discovery cohort; critically ill patients with sepsis; 30-day mortality endpoint.
- **GSE95233**: discovery cohort; day-1 samples only; 28-day mortality endpoint.
- **GSE65682**: external assessment cohort; 28-day mortality endpoint and published MARS endotype annotations.
The code is intended to reproduce the main tables, supplementary tables, and figure source outputs used in the manuscript. Minor numerical differences may occur if online databases are updated after the stated access dates.
## Data availability

The datasets analyzed in this study are publicly available from GEO:

| Dataset | GEO accession | Role in this study |
|---|---|---|
| GSE272769 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272769 | Discovery cohort |
| GSE95233 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95233 | Discovery cohort; day-1 samples |
| GSE65682 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682 | External assessment cohort |

No new human participant data were generated for this study.
## Required input files

Download the GEO Series Matrix files and the corresponding platform annotation files from GEO.

Place downloaded files under:

```text
data/
├── GSE272769_series_matrix.txt.gz
├── GSE95233_series_matrix.txt.gz
├── GSE65682_series_matrix.txt.gz
├── GPLxxxx_annotation_file_for_GSE272769.txt.gz
├── GPLxxxx_annotation_file_for_GSE95233.txt.gz
└── GPL13667_annotation_file_for_GSE65682.txt.gz
```
# The exact GEO filenames may differ depending on the download route. If filenames differ from those expected in the script, update the file path variables in the configuration section of:
```text
code/sepsis_mortality_transcriptomics_analysis.R
```
## External metadata

The overlap-exclusion sensitivity analysis uses the published Giannini E3 mortality-associated gene set.
```text
metadata/Giannini_E3_gene_set.csv
```
Recommended columns are:

```text
gene_symbol
giannini_messi_log2fc
giannini_mars_log2fc
```
If the script uses a different column naming convention, use the format specified in the header comments of the analysis script.
This comparison is a sensitivity analysis of gene-set dependence and should not be interpreted as independent external replication, because the present study and the Giannini study share public datasets.
## Software requirements

The analysis was performed in R.

Recommended R version:

```text
R >= 4.5.0
The following R packages are required:
```text
GEOquery
limma
metafor
clusterProfiler
org.Hs.eg.db
fgsea
msigdbr
xCell
enrichR
pROC
httr
jsonlite
rentrez
ggplot2
ggrepel
dplyr
tidyr
reshape2
pheatmap
patchwork
VennDiagram
```
The exact software environment used for the manuscript is recorded in:

```text
environment/sessionInfo.txt
environment/package_versions.csv
```
Before running the full analysis, install missing packages using Bioconductor or CRAN as appropriate.
Example:
```text
install.packages(c(
  "ggplot2", "ggrepel", "dplyr", "tidyr", "reshape2",
  "pheatmap", "patchwork", "VennDiagram", "metafor",
  "pROC", "httr", "jsonlite", "rentrez"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db",
  "fgsea", "msigdbr"
))
```
The xCell package may require installation according to its current distribution route.
## Running the analysis

Clone the repository:

```bash
git clone https://github.com/rain05494-star/sepsis-mortality-transcriptomics.git
cd sepsis-mortality-transcriptomics
```
Place all required GEO input files in:
```text
data/
Then run:
Rscript code/sepsis_mortality_transcriptomics_analysis.R
Alternatively, from an R session opened at the repository root:
source("code/sepsis_mortality_transcriptomics_analysis.R")
The script assumes that the working directory is the repository root. If running from another directory, provide or modify the project root path in the configuration section of the script.
## Expected outputs

The script generates output folders under:

```text
results_meta/
├── 01_preprocessing/
├── 02_discovery_meta/
├── 03_functional_annotation/
├── 04_GSE65682/
├── 05_triple_cohort/
├── 06_drug/
├── 07_tables/
└── 08_Giannini_sensitivity/
```
Depending on the final script version, folder names may differ slightly. The main output files include:

- Discovery cohort differential expression results.
- Fisher combined P-value and FDR tables.
- 32-gene signature table.
- External assessment tables for GSE65682.
- Outcome-label permutation results.
- PC1 and random-gene reference outputs.
- Hallmark GSEA tables.
- xCell enrichment-score tables.
- MARS-stratified gene-level and module-level summaries.
- Giannini E3 overlap and residual-gene sensitivity-analysis tables.
- Exploratory Enrichr, DGIdb, and PubMed annotation outputs.
## Randomness and reproducibility

Random seeds are set in the analysis script for:

- Outcome-label permutation analyses.
- Random-gene PC1 reference analyses.
- Any other simulation-based or resampling-based procedures.

The main deterministic analyses should reproduce exactly when the same input files, software versions, and local settings are used.

Analyses using external online resources may show minor changes over time because the underlying databases are updated.
## External resources and access dates

The following public resources were used:

| Resource | Purpose | Access or cutoff date |
|---|---|---|
| GEO | Public expression datasets and platform annotations | July 19, 2026 |
| MSigDB Hallmark via `msigdbr` | Hallmark GSEA gene sets | July 19, 2026 |
| Enrichr | Transcription factor target-set enrichment and exploratory perturbational annotation | [TODO: date] |
| DGIdb v5 | Curated drug-gene interaction records | July 19, 2026 |
| PubMed E-utilities | Bibliometric counts for sepsis-related compound mentions | Publication-date cutoff: July 19, 2026 |

Because Enrichr, DGIdb, PubMed, MSigDB, and gene annotation resources are periodically updated, rerunning the analysis at a later date may produce minor differences in database-derived annotations.
