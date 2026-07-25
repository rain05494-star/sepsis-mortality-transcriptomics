# Data directory

This directory is intended for locally downloaded public GEO input files.

The raw and processed GEO files are not redistributed in this repository. Users should download the required files directly from the Gene Expression Omnibus (GEO).

## Required GEO datasets

The analysis uses the following public GEO datasets:

- GSE272769
- GSE95233
- GSE65682

## Download links

- GSE272769: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272769
- GSE95233: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95233
- GSE65682: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682

## Expected local placement

After downloading the GEO Series Matrix files and the corresponding platform annotation files, place them in this directory.

A typical local structure is:

```text
data/
├── GSE272769_series_matrix.txt.gz
├── GSE95233_series_matrix.txt.gz
├── GSE65682_series_matrix.txt.gz
├── GPLxxxx_annotation_file_for_GSE272769.txt.gz
├── GPLxxxx_annotation_file_for_GSE95233.txt.gz
└── GPL13667_annotation_file_for_GSE65682.txt.gz
