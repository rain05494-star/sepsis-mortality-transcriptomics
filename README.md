# Sepsis mortality transcriptomics and single-cell reference mapping

Analysis code and data-source documentation for the manuscript:

**Cross-cohort and single-cell analyses characterize proliferative and granule-associated components of mortality-associated transcription in sepsis**

## Study overview

This study integrates gene-level evidence from two whole-blood discovery
cohorts, assesses the resulting signature in a third mortality cohort,
and maps the signature to a single-cell reference.

The analysis distinguishes relative program scores from
captured-transcript contributions across reference cell states.

The main analytical components are:

- Cross-cohort identification of a 32-gene mortality-associated signature.
- External assessment of gene-level effect directions and effect estimates.
- Functional annotation of investigator-defined proliferative
  (P-program) and granule-associated (G-program) gene sets.
- Comparison with the published Giannini E3 gene set.
- Donor-aware single-cell reference mapping of gene expression,
  program scores, and captured-transcript contributions.
- Sensitivity analyses of gene-set composition and reference-state mapping.
- Supplementary acute–convalescent comparisons.
- Exploratory mortality contrasts across MARS endotypes.

## Repository organization

Scripts are organized by analysis type:

| Location | Contents |
| --- | --- |
| [code/bulk](code/bulk/) | Bulk transcriptomic analysis scripts |
| [code/SCRNA](code/SCRNA/) | Single-cell reference analysis scripts |
| [data/GEO_accessions.txt](data/GEO_accessions.txt) | Dataset accessions and study roles |
| [data/README.md](data/README.md) | Data sources and input preparation |
| [environment](environment/) | Software environment records |
| [CITATION.cff](CITATION.cff) | Software citation metadata |
| [LICENSE](LICENSE) | Repository license |

Script filenames, their roles, and execution order should be documented
in the README files for the corresponding code directories.

## Datasets

| Dataset | Role in this study | Outcome or reference information |
| --- | --- | --- |
| [GSE272769](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272769) | Discovery cohort | 30-day mortality |
| [GSE95233](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95233) | Discovery cohort; admission day-1 samples | 28-day mortality |
| [GSE65682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682) | External assessment and exploratory endotype analyses | 28-day mortality and published MARS annotations |
| [GSE216009](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216009) | Single-cell reference mapping | Source-defined cell-state and sample annotations |

Primary single-cell mapping uses 151,837 acute-infection cells from
26 donors. Supplementary analyses also use convalescent,
healthy-control, and surgery-control samples.

Mortality outcomes were unavailable in the single-cell reference used
for this study.

No new primary datasets were generated.

## Input preparation

### Bulk transcriptomic analyses

Required inputs include:

- Processed expression matrices.
- Platform annotations for probe-to-gene mapping.
- Sample metadata and mortality outcomes.
- Published MARS endotype annotations for the relevant GSE65682 samples.
- The published Giannini E3 reference gene set.

The GSE95233 discovery analysis uses admission day-1 samples only.

### Single-cell reference analyses

The analyses require the processed GSE216009 reference object,
including the RNA counts and source-defined cell-state and sample metadata.

The local object filename used in the analysis scripts is:

`GSE216009_rhapsody_wholeblood_sobj.rds.gz`

This is the local filename expected by the scripts, not a verified
GEO download filename.

The supplied single-cell scripts use metadata fields including
`fine_annot`, `sample_id`, and `diagnosis`.

A raw count matrix alone does not contain all annotations required
by the mapping workflow.

### Giannini E3 reference set

The E3 reference set was obtained from:

`TableE3_Transcohort_DE_Genes.xlsx`

This supplementary file accompanies the Giannini study cited in the
manuscript. Use the original supplementary file and any preprocessing
steps required by the analysis scripts.

The E3 derivation datasets overlap with cohorts used in this study;
the E3 comparison therefore has a shared-cohort context.

Additional input guidance is provided in [data/README.md](data/README.md).

## Software environment

The analyses were performed in R.

Environment records are stored under [environment](environment/).
Use the record corresponding to the relevant analysis stage and
execution environment.

Bulk and single-cell analyses may require different package sets.
Relevant packages include:

- Bulk analyses: GEOquery, limma, metafor, clusterProfiler,
  org.Hs.eg.db, fgsea, msigdbr, xCell, enrichR, and pROC.
- Single-cell analyses: SeuratObject, Matrix, and UCell.
- Visualization and data handling: ggplot2 and additional packages
  imported by the individual scripts.

This list summarizes major dependencies; the scripts and
stage-specific environment records determine the complete requirements.

Do not substitute current package versions for the versions used
to generate the manuscript results without checking compatibility.

## Running the analyses

Clone the repository:

```bash
git clone https://github.com/rain05494-star/sepsis-mortality-transcriptomics.git
cd sepsis-mortality-transcriptomics
```

Before running scripts:

1. Obtain the required input files.
2. Identify the scripts corresponding to the intended analysis stage.
3. Review the software requirements.
4. Update project, input, and output paths in the configuration sections.
5. Confirm that upstream outputs required by downstream scripts exist.

Some scripts contain server-specific paths. Moving a script into
`code/bulk/` or `code/SCRNA/` does not automatically update its
internal file paths.

Placing input files in the repository's `data/` directory also does
not automatically configure every script.

### Analysis stages

The bulk workflow covers:

1. Input processing and gene mapping.
2. Cohort-specific mortality contrasts and cross-cohort integration.
3. Functional annotation and enrichment.
4. External assessment and E3 comparisons.
5. Exploratory MARS endotype analyses.
6. Figure generation and table assembly.

The single-cell workflow covers:

1. Reference-object inspection, gene mapping, and detectability.
2. Donor-state gene-expression and transcript-contribution summaries.
3. Program scores and program-level contribution profiles.
4. E3-defined subset comparisons and leave-one-gene-out analyses.
5. Technical sensitivity and acute–convalescent analyses.
6. Figure generation and table assembly.

These stages describe the analytical workflow. Follow the verified
script-level instructions in each code directory rather than treating
every script in the directory as a sequential pipeline.

## Outputs

Depending on the script, outputs include:

- Gene-level differential-expression and signature tables.
- Functional enrichment results.
- External-assessment and MARS endotype summaries.
- Gene- and program-level single-cell reference mappings.
- Donor-aware program-score and transcript-contribution summaries.
- Sensitivity-analysis and quality-control tables.
- Figure files and supplementary-table source files.

Output paths and filenames are defined by the individual scripts.

Keep outputs from distinct runs separate and ensure that plotting
and table-assembly scripts use the intended upstream results.

## Reproducibility

Reproduction requires matching:

- Input data and sample selections.
- Gene identifiers and platform annotations.
- Reference-object annotations.
- Gene-set definitions and analysis thresholds.
- Software versions.
- Random seeds for resampling procedures.
- Annotation-library versions or fixed term–gene mappings.

Where an analysis uses cached annotation resources or fixed mappings,
reuse those resources when reproducing the reported results.
A fresh online query may return different annotations.

Retain run logs and provenance files to connect each output to its
input files, script version, and execution environment.

## Versioning and citation

The `main` branch contains ongoing repository updates.

For manuscript reproduction, use the tagged release cited in the
submitted manuscript. Changes committed to `main` are not
automatically included in an earlier tagged release.

Available releases:
https://github.com/rain05494-star/sepsis-mortality-transcriptomics/releases

Software citation metadata is provided in [CITATION.cff](CITATION.cff).
Use the citation metadata associated with the release being used.

## License

Repository code is distributed under the [MIT License](LICENSE).

Third-party datasets and reference resources remain subject to
the terms of their original providers.
