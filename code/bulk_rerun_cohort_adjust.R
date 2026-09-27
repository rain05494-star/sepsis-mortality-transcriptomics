graphics.off()
rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(pROC)
  library(fgsea)
  library(msigdbr)
})

# ============================================================
# 0. Project paths
# ============================================================

command_args <- commandArgs(
  trailingOnly = TRUE
)

if (
  length(command_args) >= 1L &&
  nzchar(command_args[1])
) {
  PROJECT_DIR <- normalizePath(
    command_args[1],
    winslash = "/",
    mustWork = TRUE
  )
} else {
  PROJECT_DIR <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}

DATA_DIR <- file.path(
  PROJECT_DIR,
  "data"
)

ROOT_DIR <- file.path(
  PROJECT_DIR,
  "results_meta"
)

S2_DIR <- file.path(
  ROOT_DIR,
  "02_meta_analysis"
)

# 本轮所有新结果放入独立目录，不覆盖旧结果
REVISION_DIR <- file.path(
  ROOT_DIR,
  "08_priority1_revision"
)

S4_DIR <- file.path(
  REVISION_DIR,
  "04_GSE65682"
)

S5_DIR <- file.path(
  REVISION_DIR,
  "05_GSEA"
)

MARS_INTERACTION_DIR <- file.path(
  REVISION_DIR,
  "09_MARS_module_interaction"
)

output_directories <- c(
  REVISION_DIR,
  S4_DIR,
  S5_DIR,
  MARS_INTERACTION_DIR
)

for (
  output_directory in
  output_directories
) {
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 1. Required local inputs
# ============================================================

checkpoint_s2_file <- file.path(
  S2_DIR,
  "checkpoint_s2.rds"
)

gse65682_file <- file.path(
  DATA_DIR,
  "GSE65682_series_matrix.txt.gz"
)

gpl13667_file <- file.path(
  DATA_DIR,
  "GPL13667.soft.gz"
)

required_files <- c(
  checkpoint_s2_file,
  gse65682_file,
  gpl13667_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(sprintf(
    "Missing required files:\n%s",
    paste(
      paste0(
        "  - ",
        missing_files
      ),
      collapse = "\n"
    )
  ))
}

cat(
  "\nProject directory:\n",
  PROJECT_DIR,
  "\n\n"
)

cat(
  "Revision output directory:\n",
  REVISION_DIR,
  "\n\n"
)

cat(
  "Required input files:\n"
)

print(
  data.frame(
    file = required_files,
    exists = file.exists(
      required_files
    ),
    size_bytes = file.info(
      required_files
    )$size,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

# ============================================================
# 2. Load discovery checkpoint
# ============================================================

checkpoint_s2 <- readRDS(
  checkpoint_s2_file
)

required_checkpoint_objects <- c(
  "meta_df",
  "fisher_out",
  "gmat_272769",
  "gmat_95233",
  "gmat_272769_linear",
  "gmat_95233_linear",
  "group_272769",
  "group_95233",
  "common_genes",
  "module_map",
  "res_272769",
  "res_95233"
)

missing_checkpoint_objects <- setdiff(
  required_checkpoint_objects,
  names(checkpoint_s2)
)

if (
  length(
    missing_checkpoint_objects
  ) > 0L
) {
  stop(sprintf(
    "checkpoint_s2 is missing objects: %s",
    paste(
      missing_checkpoint_objects,
      collapse = ", "
    )
  ))
}

list2env(
  checkpoint_s2,
  envir = .GlobalEnv
)

rm(checkpoint_s2)

# ============================================================
# 3. Restore fixed analysis constants
# ============================================================

LFC_CUT <- 0.5
FDR_CUT <- 0.05

logFC_COL_95233 <- "logFC_95233 (D01)"
pval_COL_95233 <- "pval_95233 (D01)"
fdr_COL_95233 <- "fdr_95233 (D01)"
t_COL_95233 <- "t_95233 (D01)"
ave_COL_95233 <- "AveExpr_95233 (D01)"

MODULE_NEUTRO <- names(
  module_map
)[
  module_map ==
    "Neutrophil degranulation"
]

MODULE_CELL_CYCLE <- names(
  module_map
)[
  module_map ==
    "Cell cycle / proliferation"
]

MODULE_DOWN_DEFENSE <- names(
  module_map
)[
  module_map ==
    "Inflammatory / Down"
]

MODULE_OTHER <- names(
  module_map
)[
  module_map ==
    "Other"
]

MODULE_LEVELS <- c(
  "Neutrophil degranulation",
  "Cell cycle / proliferation",
  "Inflammatory / Down",
  "Other",
  "Unclassified"
)

# ============================================================
# 4. Checkpoint integrity audit
# ============================================================

required_meta_columns <- c(
  "gene",
  "logFC_272769",
  "t_272769",
  "pval_272769",
  logFC_COL_95233,
  t_COL_95233,
  pval_COL_95233
)

missing_meta_columns <- setdiff(
  required_meta_columns,
  colnames(meta_df)
)

if (
  length(
    missing_meta_columns
  ) > 0L
) {
  stop(sprintf(
    "meta_df is missing columns: %s",
    paste(
      missing_meta_columns,
      collapse = ", "
    )
  ))
}

stopifnot(
  nrow(fisher_out) == 32L,
  length(unique(fisher_out$gene)) == 32L,
  sum(
    fisher_out$direction == "Up"
  ) == 29L,
  sum(
    fisher_out$direction == "Down"
  ) == 3L,
  length(MODULE_NEUTRO) == 10L,
  length(MODULE_CELL_CYCLE) == 11L,
  length(MODULE_DOWN_DEFENSE) == 3L,
  length(MODULE_OTHER) == 8L,
  ncol(gmat_272769) ==
    length(group_272769),
  ncol(gmat_95233) ==
    length(group_95233)
)

cat(
  "\n========== CHECKPOINT S2 AUDIT ==========\n"
)

cat(sprintf(
  "Signature genes: %d\n",
  nrow(fisher_out)
))

cat(sprintf(
  "Upregulated: %d\n",
  sum(
    fisher_out$direction == "Up"
  )
))

cat(sprintf(
  "Downregulated: %d\n",
  sum(
    fisher_out$direction == "Down"
  )
))

cat(sprintf(
  "M2 genes: %d\n",
  length(MODULE_NEUTRO)
))

cat(sprintf(
  "M1 genes: %d\n",
  length(MODULE_CELL_CYCLE)
))

cat(sprintf(
  "Inflammatory/downregulated genes: %d\n",
  length(MODULE_DOWN_DEFENSE)
))

cat(sprintf(
  "Other genes: %d\n",
  length(MODULE_OTHER)
))

cat(sprintf(
  "GSE272769: %d genes × %d samples\n",
  nrow(gmat_272769),
  ncol(gmat_272769)
))

cat(sprintf(
  "GSE95233: %d genes × %d samples\n",
  nrow(gmat_95233),
  ncol(gmat_95233)
))

cat(sprintf(
  "GSE95233 t-statistic column: %s\n",
  t_COL_95233
))

cat(
  "\nCheckpoint S2 loaded successfully.\n"
)

cat(
  "No online query has been executed.\n"
)

cat(
  "=========================================\n"
)
# ============================================================
# Step 2: GSE65682 metadata and cohort audit
# ============================================================

cat(
  "\n========== STEP 2: GSE65682 METADATA AUDIT ==========\n"
)

# ------------------------------------------------------------
# 2.1 Read the local Series Matrix
# ------------------------------------------------------------

gset <- GEOquery::getGEO(
  filename = gse65682_file,
  getGPL = FALSE
)

pdata <- Biobase::pData(
  gset
)

common_samples_65682 <- intersect(
  colnames(
    Biobase::exprs(gset)
  ),
  rownames(pdata)
)

if (
  length(common_samples_65682) == 0L
) {
  stop(
    "No common sample IDs were found between expression data and metadata."
  )
}

gset <- gset[
  ,
  common_samples_65682
]

pdata <- pdata[
  common_samples_65682,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(
      Biobase::exprs(gset)
    ),
    rownames(pdata)
  ),
  ncol(
    Biobase::exprs(gset)
  ) == nrow(pdata)
)

cat(sprintf(
  "GSE65682 deposited records: %d\n",
  nrow(pdata)
))

cat(sprintf(
  "Expression matrix: %d probes × %d records\n",
  nrow(
    Biobase::exprs(gset)
  ),
  ncol(
    Biobase::exprs(gset)
  )
))

# ------------------------------------------------------------
# 2.2 Verify required metadata schema
# ------------------------------------------------------------

mortality_col <- "mortality_event_28days:ch1"
mars_col <- "endotype_class:ch1"
cohort_col <- "endotype_cohort:ch1"

required_metadata_columns <- c(
  mortality_col,
  mars_col,
  cohort_col
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(pdata)
)

if (
  length(
    missing_metadata_columns
  ) > 0L
) {
  cohort_candidates <- grep(
    "cohort",
    colnames(pdata),
    value = TRUE,
    ignore.case = TRUE
  )
  
  mortality_candidates <- grep(
    "mortality",
    colnames(pdata),
    value = TRUE,
    ignore.case = TRUE
  )
  
  endotype_candidates <- grep(
    "endotype",
    colnames(pdata),
    value = TRUE,
    ignore.case = TRUE
  )
  
  stop(sprintf(
    paste0(
      "Missing required metadata columns: %s\n",
      "Cohort-like columns: %s\n",
      "Mortality-like columns: %s\n",
      "Endotype-like columns: %s"
    ),
    paste(
      missing_metadata_columns,
      collapse = ", "
    ),
    paste(
      cohort_candidates,
      collapse = ", "
    ),
    paste(
      mortality_candidates,
      collapse = ", "
    ),
    paste(
      endotype_candidates,
      collapse = ", "
    )
  ))
}

cat(
  "Required metadata fields verified:\n"
)

print(
  required_metadata_columns
)

# ------------------------------------------------------------
# 2.3 Normalize mortality, MARS and cohort annotations
# ------------------------------------------------------------

mortality_raw <- trimws(
  as.character(
    pdata[[mortality_col]]
  )
)

mars_raw <- trimws(
  as.character(
    pdata[[mars_col]]
  )
)

cohort_raw <- trimws(
  as.character(
    pdata[[cohort_col]]
  )
)

valid_mortality <- mortality_raw %in%
  c(
    "0",
    "1"
  )

valid_mars <- mars_raw %in%
  c(
    "Mars1",
    "Mars2",
    "Mars3",
    "Mars4"
  )

cohort_lower <- tolower(
  cohort_raw
)

cohort_missing <- is.na(cohort_raw) |
  cohort_lower %in% c(
    "",
    "na",
    "n/a",
    "unknown"
  )

cohort_std <- rep(
  NA_character_,
  length(cohort_raw)
)

cohort_std[
  !cohort_missing &
    cohort_lower == "discovery"
] <- "Discovery"

cohort_std[
  !cohort_missing &
    cohort_lower == "validation"
] <- "Validation"

unrecognized_cohort <- !cohort_missing &
  is.na(cohort_std)

if (
  any(
    unrecognized_cohort
  )
) {
  stop(sprintf(
    "Unrecognized endotype-cohort values: %s",
    paste(
      unique(
        cohort_raw[
          unrecognized_cohort
        ]
      ),
      collapse = ", "
    )
  ))
}

valid_cohort <- !is.na(
  cohort_std
)

# Analysis eligibility remains defined by mortality and MARS.
include_65682 <- valid_mortality &
  valid_mars

if (
  any(
    include_65682 &
    !valid_cohort
  )
) {
  missing_ids <- rownames(pdata)[
    include_65682 &
      !valid_cohort
  ]
  
  stop(sprintf(
    paste0(
      "%d mortality/MARS-eligible samples lacked ",
      "a valid cohort assignment: %s"
    ),
    length(missing_ids),
    paste(
      head(
        missing_ids,
        20L
      ),
      collapse = ", "
    )
  ))
}

# ------------------------------------------------------------
# 2.4 Build eligibility audit
# ------------------------------------------------------------

selection_summary_65682 <- data.frame(
  Item = c(
    "Deposited GSE65682 sample records",
    "Samples with valid 28-day mortality",
    "Samples with valid MARS classification",
    "Samples with valid endotype-cohort classification",
    "Samples with mortality and MARS",
    "Mortality available but MARS unavailable",
    "MARS available but mortality unavailable",
    "Eligible samples missing cohort",
    "Final outcome-analysis population"
  ),
  N = c(
    nrow(pdata),
    sum(valid_mortality),
    sum(valid_mars),
    sum(valid_cohort),
    sum(
      valid_mortality &
        valid_mars
    ),
    sum(
      valid_mortality &
        !valid_mars
    ),
    sum(
      !valid_mortality &
        valid_mars
    ),
    sum(
      include_65682 &
        !valid_cohort
    ),
    sum(include_65682)
  ),
  stringsAsFactors = FALSE
)

sample_audit_65682 <- data.frame(
  sample_id = rownames(pdata),
  mortality_event_28days =
    mortality_raw,
  endotype_class =
    mars_raw,
  endotype_cohort_raw =
    cohort_raw,
  endotype_cohort =
    cohort_std,
  valid_mortality =
    valid_mortality,
  valid_MARS =
    valid_mars,
  valid_endotype_cohort =
    valid_cohort,
  included =
    include_65682,
  stringsAsFactors = FALSE
)

print(
  selection_summary_65682,
  row.names = FALSE
)

# ------------------------------------------------------------
# 2.5 Construct the 479-sample analysis population
# ------------------------------------------------------------

keep_idx <- which(
  include_65682
)

gset_sepsis <- gset[
  ,
  keep_idx
]

pdata_sepsis <- pdata[
  keep_idx,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(
      Biobase::exprs(
        gset_sepsis
      )
    ),
    rownames(
      pdata_sepsis
    )
  )
)

group_val <- factor(
  ifelse(
    mortality_raw[
      keep_idx
    ] == "1",
    "NonSurvivor",
    "Survivor"
  ),
  levels = c(
    "Survivor",
    "NonSurvivor"
  )
)

endo_val <- mars_raw[
  keep_idx
]

cohort_val <- factor(
  cohort_std[
    keep_idx
  ],
  levels = c(
    "Discovery",
    "Validation"
  )
)

sample_ids_val <- rownames(
  pdata_sepsis
)

names(group_val) <- sample_ids_val
names(endo_val) <- sample_ids_val
names(cohort_val) <- sample_ids_val

stopifnot(
  length(group_val) == 479L,
  length(endo_val) == 479L,
  length(cohort_val) == 479L,
  !anyNA(group_val),
  !anyNA(endo_val),
  !anyNA(cohort_val),
  identical(
    names(group_val),
    sample_ids_val
  ),
  identical(
    names(endo_val),
    sample_ids_val
  ),
  identical(
    names(cohort_val),
    sample_ids_val
  )
)

# ------------------------------------------------------------
# 2.6 Fixed count checks
# ------------------------------------------------------------

outcome_counts <- table(
  group_val
)

cohort_counts <- table(
  cohort_val
)

mars_counts <- table(
  factor(
    endo_val,
    levels = c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    )
  )
)

stopifnot(
  as.integer(
    outcome_counts[
      "Survivor"
    ]
  ) == 365L,
  as.integer(
    outcome_counts[
      "NonSurvivor"
    ]
  ) == 114L,
  as.integer(
    cohort_counts[
      "Discovery"
    ]
  ) == 263L,
  as.integer(
    cohort_counts[
      "Validation"
    ]
  ) == 216L
)

mars_cohort_counts <- table(
  MARS = factor(
    endo_val,
    levels = c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    )
  ),
  Cohort = cohort_val
)

expected_mars_cohort_counts <- matrix(
  c(
    72L, 60L,
    97L, 79L,
    60L, 58L,
    34L, 19L
  ),
  nrow = 4L,
  byrow = TRUE,
  dimnames = list(
    MARS = c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    ),
    Cohort = c(
      "Discovery",
      "Validation"
    )
  )
)

print(
  mars_cohort_counts
)

print(
  expected_mars_cohort_counts
)

stopifnot(
  identical(
    dimnames(mars_cohort_counts),
    dimnames(expected_mars_cohort_counts)
  ),
  identical(
    as.integer(mars_cohort_counts),
    as.integer(expected_mars_cohort_counts)
  )
)

mars_cohort_outcome_counts <- table(
  MARS = factor(
    endo_val,
    levels = c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    )
  ),
  Cohort = cohort_val,
  Outcome = group_val
)

if (
  any(
    mars_cohort_outcome_counts ==
    0L
  )
) {
  stop(paste0(
    "At least one MARS × cohort × outcome cell was empty:\n",
    paste(
      capture.output(
        print(
          mars_cohort_outcome_counts
        )
      ),
      collapse = "\n"
    )
  ))
}

# Known sparse cell, retained as an auditable fact.
stopifnot(
  as.integer(
    mars_cohort_outcome_counts[
      "Mars4",
      "Validation",
      "NonSurvivor"
    ]
  ) == 1L
)

# ------------------------------------------------------------
# 2.7 Export metadata audit
# ------------------------------------------------------------

write.csv(
  selection_summary_65682,
  file.path(
    S4_DIR,
    "GSE65682_selection_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_audit_65682,
  file.path(
    S4_DIR,
    "GSE65682_sample_inclusion_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    mars_cohort_counts
  ),
  file.path(
    S4_DIR,
    "GSE65682_MARS_by_cohort_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    mars_cohort_outcome_counts
  ),
  file.path(
    S4_DIR,
    "GSE65682_MARS_by_cohort_by_outcome_counts.csv"
  ),
  row.names = FALSE
)

metadata_checkpoint <- list(
  pdata_sepsis = pdata_sepsis,
  sample_ids = sample_ids_val,
  group_val = group_val,
  endo_val = endo_val,
  cohort_val = cohort_val,
  selection_summary =
    selection_summary_65682,
  sample_audit =
    sample_audit_65682,
  mars_cohort_counts =
    mars_cohort_counts,
  mars_cohort_outcome_counts =
    mars_cohort_outcome_counts
)

saveRDS(
  metadata_checkpoint,
  file.path(
    S4_DIR,
    "checkpoint_step2_metadata.rds"
  )
)

cat(
  "\nOutcome counts:\n"
)

print(
  outcome_counts
)

cat(
  "\nCohort counts:\n"
)

print(
  cohort_counts
)

cat(
  "\nMARS counts:\n"
)

print(
  mars_counts
)

cat(
  "\nMARS × cohort counts:\n"
)

print(
  mars_cohort_counts
)

cat(
  "\nMARS × cohort × outcome counts:\n"
)

print(
  mars_cohort_outcome_counts
)

cat(
  "\nSTEP 2 completed successfully.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 3: GPL13667 annotation and gene-level expression matrix
# ============================================================

cat(
  "\n========== STEP 3: GSE65682 GENE-LEVEL MATRIX ==========\n"
)

# ------------------------------------------------------------
# 3.1 Required-object audit
# ------------------------------------------------------------

required_step3_objects <- c(
  "gpl13667_file",
  "gset_sepsis",
  "pdata_sepsis",
  "sample_ids_val",
  "group_val",
  "endo_val",
  "cohort_val",
  "fisher_out",
  "S4_DIR"
)

missing_step3_objects <- required_step3_objects[
  !vapply(
    required_step3_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (
  length(missing_step3_objects) > 0L
) {
  stop(sprintf(
    "Missing objects required for Step 3: %s",
    paste(
      missing_step3_objects,
      collapse = ", "
    )
  ))
}

if (
  !file.exists(gpl13667_file)
) {
  stop(sprintf(
    "Local GPL13667 annotation file was not found: %s",
    gpl13667_file
  ))
}

stopifnot(
  identical(
    colnames(
      Biobase::exprs(gset_sepsis)
    ),
    sample_ids_val
  ),
  identical(
    rownames(pdata_sepsis),
    sample_ids_val
  ),
  identical(
    names(group_val),
    sample_ids_val
  ),
  identical(
    names(endo_val),
    sample_ids_val
  ),
  identical(
    names(cohort_val),
    sample_ids_val
  ),
  length(sample_ids_val) == 479L
)

# ------------------------------------------------------------
# 3.2 Read local GPL13667 annotation
# ------------------------------------------------------------

gpl13667 <- GEOquery::getGEO(
  filename = gpl13667_file,
  GSEMatrix = FALSE
)

tab13667 <- GEOquery::Table(
  gpl13667
)

gene_symbol_col <- "Gene Symbol"

if (
  !gene_symbol_col %in%
  colnames(tab13667)
) {
  gene_symbol_candidates <- grep(
    "gene.*symbol|symbol",
    colnames(tab13667),
    value = TRUE,
    ignore.case = TRUE
  )
  
  stop(sprintf(
    paste0(
      "Expected GPL13667 column '%s' was not found.\n",
      "Symbol-like columns: %s"
    ),
    gene_symbol_col,
    paste(
      gene_symbol_candidates,
      collapse = ", "
    )
  ))
}

if (
  !"ID" %in%
  colnames(tab13667)
) {
  stop(
    "GPL13667 annotation does not contain the required ID column."
  )
}

cat(sprintf(
  "GPL13667 annotation: %d rows × %d columns\n",
  nrow(tab13667),
  ncol(tab13667)
))

cat(sprintf(
  "Gene-symbol field: %s\n",
  gene_symbol_col
))

# ------------------------------------------------------------
# 3.3 Parse probe-to-gene annotation
# ------------------------------------------------------------

invalid_gene_step3 <- function(x) {
  x <- trimws(
    as.character(x)
  )
  
  is.na(x) |
    !nzchar(x) |
    x %in% c(
      "---",
      "NA",
      "N/A"
    ) |
    nchar(x) > 30L
}

parse_gene_symbol_step3 <- function(x) {
  vapply(
    trimws(
      as.character(x)
    ),
    function(z) {
      if (
        is.na(z) ||
        !nzchar(z)
      ) {
        return(
          NA_character_
        )
      }
      
      parts <- unlist(
        strsplit(
          z,
          "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*",
          perl = TRUE
        )
      )
      
      parts <- trimws(parts)
      parts <- parts[
        !invalid_gene_step3(parts)
      ]
      
      if (
        length(parts) == 0L
      ) {
        NA_character_
      } else {
        parts[1L]
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
}

map13667 <- data.frame(
  probe = trimws(
    as.character(
      tab13667$ID
    )
  ),
  gene = parse_gene_symbol_step3(
    tab13667[[gene_symbol_col]]
  ),
  stringsAsFactors = FALSE
)

map13667 <- map13667[
  !is.na(map13667$probe) &
    nzchar(map13667$probe) &
    !invalid_gene_step3(map13667$gene),
  ,
  drop = FALSE
]

map13667 <- map13667[
  !duplicated(
    map13667[
      ,
      c(
        "probe",
        "gene"
      )
    ]
  ),
  ,
  drop = FALSE
]

cat(sprintf(
  "Valid probe-to-gene annotation rows: %d\n",
  nrow(map13667)
))

# ------------------------------------------------------------
# 3.4 Check expression scale and sample order
# ------------------------------------------------------------

exprs_raw_65682 <- Biobase::exprs(
  gset_sepsis
)

stopifnot(
  identical(
    colnames(exprs_raw_65682),
    sample_ids_val
  )
)

finite_expression_values <- exprs_raw_65682[
  is.finite(exprs_raw_65682)
]

if (
  length(finite_expression_values) == 0L
) {
  stop(
    "GSE65682 expression matrix contains no finite values."
  )
}

range_65682 <- range(
  finite_expression_values
)

is_log2_65682 <-
  range_65682[2L] < 30 &&
  range_65682[1L] > -2

cat(sprintf(
  "GSE65682 expression range: %.4f to %.4f\n",
  range_65682[1L],
  range_65682[2L]
))

cat(sprintf(
  "Log2-like expression scale: %s\n",
  is_log2_65682
))

if (
  !is_log2_65682
) {
  stop(
    paste0(
      "GSE65682 does not satisfy the expected log2-like ",
      "processed-expression range. Do not transform automatically."
    )
  )
}

# ------------------------------------------------------------
# 3.5 Collapse probes by maximum across-sample variance
# ------------------------------------------------------------

collapse_probes_step3 <- function(
    expression_matrix,
    probe_map
) {
  expression_matrix <- as.matrix(
    expression_matrix
  )
  
  common_probes <- intersect(
    rownames(expression_matrix),
    probe_map$probe
  )
  
  if (
    length(common_probes) == 0L
  ) {
    stop(
      "No GPL13667 probes matched the GSE65682 expression matrix."
    )
  }
  
  expression_subset <- expression_matrix[
    common_probes,
    ,
    drop = FALSE
  ]
  
  matched_map <- probe_map[
    match(
      common_probes,
      probe_map$probe
    ),
    ,
    drop = FALSE
  ]
  
  gene_groups <- split(
    seq_len(
      nrow(expression_subset)
    ),
    matched_map$gene
  )
  
  selected_indices <- vapply(
    gene_groups,
    function(index) {
      probe_variance <- apply(
        expression_subset[
          index,
          ,
          drop = FALSE
        ],
        1,
        stats::var,
        na.rm = TRUE
      )
      
      probe_variance[
        !is.finite(probe_variance)
      ] <- -Inf
      
      if (
        all(probe_variance == -Inf)
      ) {
        index[1L]
      } else {
        index[
          which.max(probe_variance)
        ]
      }
    },
    integer(1)
  )
  
  gene_matrix <- expression_subset[
    selected_indices,
    ,
    drop = FALSE
  ]
  
  rownames(gene_matrix) <- names(
    gene_groups
  )
  
  selected_variance <- vapply(
    selected_indices,
    function(index) {
      current_variance <- stats::var(
        expression_subset[
          index,
          ,
          drop = TRUE
        ],
        na.rm = TRUE
      )
      
      if (
        is.finite(current_variance)
      ) {
        current_variance
      } else {
        NA_real_
      }
    },
    numeric(1)
  )
  
  selection_audit <- data.frame(
    gene = names(gene_groups),
    selected_probe =
      rownames(expression_subset)[
        selected_indices
      ],
    candidate_probe_n =
      lengths(gene_groups),
    selected_probe_variance =
      selected_variance,
    stringsAsFactors = FALSE
  )
  
  list(
    matrix = gene_matrix,
    selection_audit = selection_audit,
    matched_probe_n =
      length(common_probes)
  )
}

collapse_result_65682 <- collapse_probes_step3(
  expression_matrix =
    exprs_raw_65682,
  probe_map =
    map13667
)

gmat_val <- collapse_result_65682$matrix

probe_selection_audit_65682 <-
  collapse_result_65682$selection_audit

stopifnot(
  ncol(gmat_val) == 479L,
  identical(
    colnames(gmat_val),
    sample_ids_val
  ),
  !anyDuplicated(
    rownames(gmat_val)
  )
)

cat(sprintf(
  "Expression probes matched to annotation: %d\n",
  collapse_result_65682$matched_probe_n
))

cat(sprintf(
  "Gene-level expression matrix: %d genes × %d samples\n",
  nrow(gmat_val),
  ncol(gmat_val)
))

# ------------------------------------------------------------
# 3.6 Genome-wide non-finite-value audit
# ------------------------------------------------------------

nonfinite_total_65682 <- sum(
  !is.finite(gmat_val)
)

nonfinite_gene_n_65682 <- sum(
  rowSums(
    !is.finite(gmat_val)
  ) > 0L
)

nonfinite_sample_n_65682 <- sum(
  colSums(
    !is.finite(gmat_val)
  ) > 0L
)

if (
  nonfinite_total_65682 > 0L
) {
  warning(sprintf(
    paste0(
      "The complete gene-level matrix contains %d non-finite ",
      "values across %d genes and %d samples. ",
      "Strict finite-value checks will be applied to ",
      "the evaluable signature genes."
    ),
    nonfinite_total_65682,
    nonfinite_gene_n_65682,
    nonfinite_sample_n_65682
  ))
}

# ------------------------------------------------------------
# 3.7 Signature-gene coverage audit
# ------------------------------------------------------------

signature_genes <- unique(
  as.character(
    fisher_out$gene
  )
)

if (
  length(signature_genes) != 32L
) {
  stop(sprintf(
    "Expected 32 signature genes but found %d.",
    length(signature_genes)
  ))
}

signature_adjusted <- signature_genes[
  signature_genes %in%
    rownames(gmat_val)
]

missing_signature_genes <- setdiff(
  signature_genes,
  rownames(gmat_val)
)

expected_missing_signature_genes <- c(
  "HIST1H3B",
  "HIST1H2BM"
)

if (
  length(signature_adjusted) != 30L ||
  !setequal(
    missing_signature_genes,
    expected_missing_signature_genes
  )
) {
  stop(sprintf(
    paste0(
      "Unexpected signature coverage in GSE65682.\n",
      "Detected: %d/32\n",
      "Missing: %s"
    ),
    length(signature_adjusted),
    paste(
      missing_signature_genes,
      collapse = ", "
    )
  ))
}

signature_matrix_65682 <- gmat_val[
  signature_adjusted,
  ,
  drop = FALSE
]

bad_signature_cells <- which(
  !is.finite(signature_matrix_65682),
  arr.ind = TRUE
)

if (
  nrow(bad_signature_cells) > 0L
) {
  bad_examples <- apply(
    head(
      bad_signature_cells,
      20L
    ),
    1,
    function(index) {
      paste0(
        rownames(signature_matrix_65682)[
          index[1L]
        ],
        "@",
        colnames(signature_matrix_65682)[
          index[2L]
        ]
      )
    }
  )
  
  stop(sprintf(
    paste0(
      "Non-finite values were detected in evaluable ",
      "signature genes: %s"
    ),
    paste(
      bad_examples,
      collapse = ", "
    )
  ))
}

signature_sd_65682 <- apply(
  signature_matrix_65682,
  1,
  stats::sd
)

if (
  any(
    !is.finite(signature_sd_65682) |
    signature_sd_65682 <= 0
  )
) {
  stop(sprintf(
    "Zero-variance or invalid signature genes: %s",
    paste(
      names(signature_sd_65682)[
        !is.finite(signature_sd_65682) |
          signature_sd_65682 <= 0
      ],
      collapse = ", "
    )
  ))
}

direction_col_step3 <- intersect(
  c(
    "direction",
    "Direction"
  ),
  colnames(fisher_out)
)

discovery_direction_step3 <- if (
  length(direction_col_step3) == 1L
) {
  as.character(
    fisher_out[
      match(
        signature_genes,
        fisher_out$gene
      ),
      direction_col_step3
    ]
  )
} else {
  rep(
    NA_character_,
    length(signature_genes)
  )
}

signature_detection_audit_65682 <- data.frame(
  gene = signature_genes,
  discovery_direction =
    discovery_direction_step3,
  detected_in_GSE65682 =
    signature_genes %in%
    rownames(gmat_val),
  selected_probe =
    probe_selection_audit_65682$
    selected_probe[
      match(
        signature_genes,
        probe_selection_audit_65682$gene
      )
    ],
  candidate_probe_n =
    probe_selection_audit_65682$
    candidate_probe_n[
      match(
        signature_genes,
        probe_selection_audit_65682$gene
      )
    ],
  stringsAsFactors = FALSE
)

cat(sprintf(
  "Signature genes detected: %d/32\n",
  length(signature_adjusted)
))

cat(sprintf(
  "Signature genes not detected: %s\n",
  paste(
    missing_signature_genes,
    collapse = ", "
  )
))

# ------------------------------------------------------------
# 3.8 Export preprocessing audit and checkpoint
# ------------------------------------------------------------

expression_qc_65682 <- data.frame(
  Item = c(
    "GSE65682 analysis samples",
    "GPL13667 annotation rows",
    "Valid probe-to-gene annotation rows",
    "Matched expression probes",
    "Gene-level expression features",
    "Detected signature genes",
    "Missing signature genes",
    "Genome-wide non-finite values",
    "Genes containing non-finite values",
    "Samples containing non-finite values",
    "Expression minimum",
    "Expression maximum",
    "Log2-like expression scale"
  ),
  Value = c(
    ncol(gmat_val),
    nrow(tab13667),
    nrow(map13667),
    collapse_result_65682$matched_probe_n,
    nrow(gmat_val),
    length(signature_adjusted),
    length(missing_signature_genes),
    nonfinite_total_65682,
    nonfinite_gene_n_65682,
    nonfinite_sample_n_65682,
    range_65682[1L],
    range_65682[2L],
    is_log2_65682
  ),
  stringsAsFactors = FALSE
)

write.csv(
  expression_qc_65682,
  file.path(
    S4_DIR,
    "GSE65682_expression_preprocessing_QC.csv"
  ),
  row.names = FALSE
)

write.csv(
  signature_detection_audit_65682,
  file.path(
    S4_DIR,
    "GSE65682_signature_detection_audit.csv"
  ),
  row.names = FALSE
)

expression_checkpoint <- list(
  gmat_val = gmat_val,
  signature_genes = signature_genes,
  signature_adjusted =
    signature_adjusted,
  missing_signature_genes =
    missing_signature_genes,
  probe_selection_audit =
    probe_selection_audit_65682,
  expression_qc =
    expression_qc_65682
)

saveRDS(
  expression_checkpoint,
  file.path(
    S4_DIR,
    "checkpoint_step3_expression.rds"
  )
)

cat(
  "\nExpression preprocessing QC:\n"
)

print(
  expression_qc_65682,
  row.names = FALSE
)

cat(
  "\nSignature detection audit:\n"
)

print(
  signature_detection_audit_65682,
  row.names = FALSE
)

cat(
  "\nSTEP 3 completed successfully.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 4: Cohort-adjusted and cohort-specific DE analysis
# ============================================================

cat(
  "\n========== STEP 4: GSE65682 DIFFERENTIAL EXPRESSION ==========\n"
)

# ------------------------------------------------------------
# 4.1 Required-object and design audit
# ------------------------------------------------------------

required_step4_objects <- c(
  "gmat_val",
  "group_val",
  "cohort_val",
  "sample_ids_val",
  "signature_adjusted",
  "missing_signature_genes",
  "fisher_out",
  "LFC_CUT",
  "FDR_CUT",
  "S4_DIR"
)

missing_step4_objects <- required_step4_objects[
  !vapply(
    required_step4_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (
  length(missing_step4_objects) > 0L
) {
  stop(sprintf(
    "Missing objects required for Step 4: %s",
    paste(
      missing_step4_objects,
      collapse = ", "
    )
  ))
}

stopifnot(
  ncol(gmat_val) == 479L,
  identical(
    colnames(gmat_val),
    sample_ids_val
  ),
  identical(
    names(group_val),
    sample_ids_val
  ),
  identical(
    names(cohort_val),
    sample_ids_val
  ),
  identical(
    levels(group_val),
    c(
      "Survivor",
      "NonSurvivor"
    )
  ),
  identical(
    levels(cohort_val),
    c(
      "Discovery",
      "Validation"
    )
  ),
  length(signature_adjusted) == 30L,
  setequal(
    missing_signature_genes,
    c(
      "HIST1H2BM",
      "HIST1H3B"
    )
  )
)

cohort_outcome_counts_step4 <- table(
  Cohort = cohort_val,
  Outcome = group_val
)

cat(
  "\nCohort × outcome counts:\n"
)

print(
  cohort_outcome_counts_step4
)

if (
  any(
    cohort_outcome_counts_step4 == 0L
  )
) {
  stop(
    "At least one cohort × outcome cell is empty."
  )
}

# ------------------------------------------------------------
# 4.2 Unadjusted limma function
# ------------------------------------------------------------

run_limma_unadjusted_step4 <- function(
    expression_matrix,
    outcome,
    analysis_name
) {
  outcome <- droplevels(
    factor(
      as.character(outcome),
      levels = c(
        "Survivor",
        "NonSurvivor"
      )
    )
  )
  
  stopifnot(
    ncol(expression_matrix) ==
      length(outcome),
    !anyNA(outcome),
    nlevels(outcome) == 2L
  )
  
  outcome_counts <- table(
    outcome
  )
  
  if (
    any(
      outcome_counts < 2L
    )
  ) {
    stop(sprintf(
      "%s has insufficient outcome counts: %s",
      analysis_name,
      paste(
        names(outcome_counts),
        outcome_counts,
        sep = "=",
        collapse = ", "
      )
    ))
  }
  
  design <- stats::model.matrix(
    ~ 0 + outcome
  )
  
  expected_columns <- c(
    "outcomeSurvivor",
    "outcomeNonSurvivor"
  )
  
  if (
    !all(
      expected_columns %in%
      colnames(design)
    )
  ) {
    stop(sprintf(
      "%s design columns were unexpected: %s",
      analysis_name,
      paste(
        colnames(design),
        collapse = ", "
      )
    ))
  }
  
  if (
    qr(design)$rank !=
    ncol(design)
  ) {
    stop(sprintf(
      "%s design matrix was not full rank.",
      analysis_name
    ))
  }
  
  contrast <- limma::makeContrasts(
    outcomeNonSurvivor -
      outcomeSurvivor,
    levels = design
  )
  
  fit <- limma::lmFit(
    expression_matrix,
    design
  )
  
  fit <- limma::contrasts.fit(
    fit,
    contrast
  )
  
  fit <- limma::eBayes(
    fit,
    trend = TRUE
  )
  
  result <- limma::topTable(
    fit,
    coef = 1L,
    number = Inf,
    sort.by = "none"
  )
  
  result$gene <- rownames(
    result
  )
  
  nonfinite_result <- !is.finite(
    result$logFC
  ) |
    !is.finite(
      result$t
    )
  
  if (
    any(nonfinite_result)
  ) {
    warning(sprintf(
      "%s produced non-finite logFC/t values for %d genes.",
      analysis_name,
      sum(nonfinite_result)
    ))
  }
  
  cat(sprintf(
    "%s: %d genes | Survivor=%d | NonSurvivor=%d\n",
    analysis_name,
    nrow(result),
    outcome_counts["Survivor"],
    outcome_counts["NonSurvivor"]
  ))
  
  list(
    results = result,
    design = design,
    contrast = contrast,
    outcome_counts = outcome_counts
  )
}

# ------------------------------------------------------------
# 4.3 Cohort-adjusted limma function
# ------------------------------------------------------------

run_limma_adjusted_step4 <- function(
    expression_matrix,
    outcome,
    cohort,
    analysis_name
) {
  outcome <- factor(
    as.character(outcome),
    levels = c(
      "Survivor",
      "NonSurvivor"
    )
  )
  
  cohort <- factor(
    as.character(cohort),
    levels = c(
      "Discovery",
      "Validation"
    )
  )
  
  stopifnot(
    ncol(expression_matrix) ==
      length(outcome),
    ncol(expression_matrix) ==
      length(cohort),
    !anyNA(outcome),
    !anyNA(cohort),
    nlevels(outcome) == 2L,
    nlevels(cohort) == 2L
  )
  
  cross_counts <- table(
    Cohort = cohort,
    Outcome = outcome
  )
  
  if (
    any(
      cross_counts == 0L
    )
  ) {
    stop(sprintf(
      paste0(
        "%s contains an empty cohort × outcome cell:\n%s"
      ),
      analysis_name,
      paste(
        capture.output(
          print(cross_counts)
        ),
        collapse = "\n"
      )
    ))
  }
  
  design <- stats::model.matrix(
    ~ 0 + outcome + cohort
  )
  
  expected_columns <- c(
    "outcomeSurvivor",
    "outcomeNonSurvivor",
    "cohortValidation"
  )
  
  if (
    !all(
      expected_columns %in%
      colnames(design)
    )
  ) {
    stop(sprintf(
      "%s design columns were unexpected: %s",
      analysis_name,
      paste(
        colnames(design),
        collapse = ", "
      )
    ))
  }
  
  if (
    qr(design)$rank !=
    ncol(design)
  ) {
    stop(sprintf(
      paste0(
        "%s design matrix was not full rank: ",
        "rank=%d, columns=%d."
      ),
      analysis_name,
      qr(design)$rank,
      ncol(design)
    ))
  }
  
  contrast <- limma::makeContrasts(
    outcomeNonSurvivor -
      outcomeSurvivor,
    levels = design
  )
  
  fit <- limma::lmFit(
    expression_matrix,
    design
  )
  
  fit <- limma::contrasts.fit(
    fit,
    contrast
  )
  
  fit <- limma::eBayes(
    fit,
    trend = TRUE
  )
  
  result <- limma::topTable(
    fit,
    coef = 1L,
    number = Inf,
    sort.by = "none"
  )
  
  result$gene <- rownames(
    result
  )
  
  nonfinite_result <- !is.finite(
    result$logFC
  ) |
    !is.finite(
      result$t
    )
  
  if (
    any(nonfinite_result)
  ) {
    warning(sprintf(
      "%s produced non-finite logFC/t values for %d genes.",
      analysis_name,
      sum(nonfinite_result)
    ))
  }
  
  cat(sprintf(
    paste0(
      "%s: %d genes | design rank=%d/%d | ",
      "Survivor=%d | NonSurvivor=%d\n"
    ),
    analysis_name,
    nrow(result),
    qr(design)$rank,
    ncol(design),
    sum(outcome == "Survivor"),
    sum(outcome == "NonSurvivor")
  ))
  
  list(
    results = result,
    design = design,
    contrast = contrast,
    cross_counts = cross_counts
  )
}

# ------------------------------------------------------------
# 4.4 Fit combined analyses
# ------------------------------------------------------------

fit_val_unadjusted <- run_limma_unadjusted_step4(
  expression_matrix = gmat_val,
  outcome = group_val,
  analysis_name =
    "GSE65682 combined unadjusted"
)

fit_val_adjusted <- run_limma_adjusted_step4(
  expression_matrix = gmat_val,
  outcome = group_val,
  cohort = cohort_val,
  analysis_name =
    "GSE65682 combined cohort-adjusted"
)

res_val_unadjusted <-
  fit_val_unadjusted$results

res_val_adjusted <-
  fit_val_adjusted$results

# Formal external-assessment result used downstream.
res_val <- res_val_adjusted

cat(
  "\nAdjusted design columns:\n"
)

print(
  colnames(
    fit_val_adjusted$design
  )
)

# ------------------------------------------------------------
# 4.5 Fit the two deposited subcohorts separately
# ------------------------------------------------------------

discovery_index_step4 <- which(
  cohort_val == "Discovery"
)

validation_index_step4 <- which(
  cohort_val == "Validation"
)

stopifnot(
  length(discovery_index_step4) == 263L,
  length(validation_index_step4) == 216L,
  identical(
    colnames(
      gmat_val[
        ,
        discovery_index_step4,
        drop = FALSE
      ]
    ),
    sample_ids_val[
      discovery_index_step4
    ]
  ),
  identical(
    colnames(
      gmat_val[
        ,
        validation_index_step4,
        drop = FALSE
      ]
    ),
    sample_ids_val[
      validation_index_step4
    ]
  )
)

fit_val_discovery <- run_limma_unadjusted_step4(
  expression_matrix = gmat_val[
    ,
    discovery_index_step4,
    drop = FALSE
  ],
  outcome = group_val[
    discovery_index_step4
  ],
  analysis_name =
    "GSE65682 deposited discovery cohort"
)

fit_val_validation <- run_limma_unadjusted_step4(
  expression_matrix = gmat_val[
    ,
    validation_index_step4,
    drop = FALSE
  ],
  outcome = group_val[
    validation_index_step4
  ],
  analysis_name =
    "GSE65682 deposited validation cohort"
)

res_val_discovery <-
  fit_val_discovery$results

res_val_validation <-
  fit_val_validation$results

# ------------------------------------------------------------
# 4.6 Strict checks restricted to the 30 evaluable genes
# ------------------------------------------------------------

result_sets_step4 <- list(
  combined_unadjusted =
    res_val_unadjusted,
  combined_adjusted =
    res_val_adjusted,
  discovery_cohort =
    res_val_discovery,
  validation_cohort =
    res_val_validation
)

for (
  result_name in
  names(result_sets_step4)
) {
  current_result <-
    result_sets_step4[[
      result_name
    ]]
  
  matched_index <- match(
    signature_adjusted,
    current_result$gene
  )
  
  if (
    anyNA(matched_index)
  ) {
    stop(sprintf(
      "%s was missing evaluable signature genes: %s",
      result_name,
      paste(
        signature_adjusted[
          is.na(matched_index)
        ],
        collapse = ", "
      )
    ))
  }
  
  current_signature_result <-
    current_result[
      matched_index,
      ,
      drop = FALSE
    ]
  
  required_numeric_columns <- c(
    "logFC",
    "t",
    "P.Value",
    "adj.P.Val"
  )
  
  finite_matrix <- as.matrix(
    current_signature_result[
      ,
      required_numeric_columns,
      drop = FALSE
    ]
  )
  
  if (
    any(
      !is.finite(finite_matrix)
    )
  ) {
    stop(sprintf(
      "%s contained non-finite statistics among the 30 genes.",
      result_name
    ))
  }
}

# ------------------------------------------------------------
# 4.7 Full-genome adjusted DE counts
# ------------------------------------------------------------

n_val_up <- sum(
  res_val_adjusted$logFC >
    LFC_CUT &
    res_val_adjusted$adj.P.Val <
    FDR_CUT,
  na.rm = TRUE
)

n_val_down <- sum(
  res_val_adjusted$logFC <
    -LFC_CUT &
    res_val_adjusted$adj.P.Val <
    FDR_CUT,
  na.rm = TRUE
)

cat(sprintf(
  paste0(
    "\nCohort-adjusted genome-wide DE counts: ",
    "Up=%d | Down=%d\n"
  ),
  n_val_up,
  n_val_down
))

# Deliberately no assertion against historical 36/9 counts.

# ------------------------------------------------------------
# 4.8 Build the 30-gene concordance table
# ------------------------------------------------------------

effect_direction_step4 <- function(log_fc) {
  ifelse(
    log_fc > 0,
    "Up",
    ifelse(
      log_fc < 0,
      "Down",
      "Zero"
    )
  )
}

extract_signature_statistics_step4 <- function(
    result_table,
    prefix
) {
  matched_index <- match(
    signature_adjusted,
    result_table$gene
  )
  
  current <- result_table[
    matched_index,
    ,
    drop = FALSE
  ]
  
  output <- data.frame(
    gene = signature_adjusted,
    stringsAsFactors = FALSE
  )
  
  output[[
    paste0(
      prefix,
      "_logFC"
    )
  ]] <- current$logFC
  
  output[[
    paste0(
      prefix,
      "_t"
    )
  ]] <- current$t
  
  output[[
    paste0(
      prefix,
      "_P"
    )
  ]] <- current$P.Value
  
  output[[
    paste0(
      prefix,
      "_FDR"
    )
  ]] <- current$adj.P.Val
  
  output[[
    paste0(
      prefix,
      "_direction"
    )
  ]] <- effect_direction_step4(
    current$logFC
  )
  
  output
}

combined_unadjusted_stats <-
  extract_signature_statistics_step4(
    res_val_unadjusted,
    "combined_unadjusted"
  )

combined_adjusted_stats <-
  extract_signature_statistics_step4(
    res_val_adjusted,
    "combined_adjusted"
  )

discovery_cohort_stats <-
  extract_signature_statistics_step4(
    res_val_discovery,
    "discovery_cohort"
  )

validation_cohort_stats <-
  extract_signature_statistics_step4(
    res_val_validation,
    "validation_cohort"
  )

discovery_signature_direction <-
  as.character(
    fisher_out$direction[
      match(
        signature_adjusted,
        fisher_out$gene
      )
    ]
  )

if (
  anyNA(
    discovery_signature_direction
  ) ||
  !all(
    discovery_signature_direction %in%
    c(
      "Up",
      "Down"
    )
  )
) {
  stop(
    "Discovery signature directions were missing or invalid."
  )
}

per_gene_three_way_concordance <- data.frame(
  gene = signature_adjusted,
  discovery_signature_direction =
    discovery_signature_direction,
  stringsAsFactors = FALSE
)

statistics_tables_step4 <- list(
  combined_unadjusted_stats,
  combined_adjusted_stats,
  discovery_cohort_stats,
  validation_cohort_stats
)

for (
  current_statistics in
  statistics_tables_step4
) {
  stopifnot(
    identical(
      current_statistics$gene,
      per_gene_three_way_concordance$gene
    )
  )
  
  columns_to_add <- setdiff(
    colnames(current_statistics),
    "gene"
  )
  
  per_gene_three_way_concordance[
    ,
    columns_to_add
  ] <- current_statistics[
    ,
    columns_to_add,
    drop = FALSE
  ]
}

per_gene_three_way_concordance$
  combined_unadjusted_agrees_with_signature <-
  per_gene_three_way_concordance$
  combined_unadjusted_direction ==
  per_gene_three_way_concordance$
  discovery_signature_direction

per_gene_three_way_concordance$
  combined_adjusted_agrees_with_signature <-
  per_gene_three_way_concordance$
  combined_adjusted_direction ==
  per_gene_three_way_concordance$
  discovery_signature_direction

per_gene_three_way_concordance$
  discovery_cohort_agrees_with_signature <-
  per_gene_three_way_concordance$
  discovery_cohort_direction ==
  per_gene_three_way_concordance$
  discovery_signature_direction

per_gene_three_way_concordance$
  validation_cohort_agrees_with_signature <-
  per_gene_three_way_concordance$
  validation_cohort_direction ==
  per_gene_three_way_concordance$
  discovery_signature_direction

per_gene_three_way_concordance$
  discovery_validation_same_direction <-
  per_gene_three_way_concordance$
  discovery_cohort_direction ==
  per_gene_three_way_concordance$
  validation_cohort_direction

per_gene_three_way_concordance$
  both_subcohorts_agree_with_signature <-
  per_gene_three_way_concordance$
  discovery_cohort_agrees_with_signature &
  per_gene_three_way_concordance$
  validation_cohort_agrees_with_signature

per_gene_three_way_concordance$
  all_three_GSE65682_estimates_agree_with_signature <-
  per_gene_three_way_concordance$
  combined_adjusted_agrees_with_signature &
  per_gene_three_way_concordance$
  discovery_cohort_agrees_with_signature &
  per_gene_three_way_concordance$
  validation_cohort_agrees_with_signature

agreement_columns_step4 <- c(
  "combined_unadjusted_agrees_with_signature",
  "combined_adjusted_agrees_with_signature",
  "discovery_cohort_agrees_with_signature",
  "validation_cohort_agrees_with_signature",
  "discovery_validation_same_direction",
  "both_subcohorts_agree_with_signature",
  "all_three_GSE65682_estimates_agree_with_signature"
)

three_way_concordance_summary <- data.frame(
  Metric = agreement_columns_step4,
  N = vapply(
    agreement_columns_step4,
    function(column_name) {
      as.integer(
        sum(
          per_gene_three_way_concordance[[column_name]],
          na.rm = TRUE
        )
      )
    },
    integer(1)
  ),
  Denominator = nrow(
    per_gene_three_way_concordance
  ),
  stringsAsFactors = FALSE
)

three_way_concordance_summary$Percent <-
  100 *
  three_way_concordance_summary$N /
  three_way_concordance_summary$Denominator

cat(
  "\nThirty-gene concordance summary:\n"
)

print(
  three_way_concordance_summary,
  row.names = FALSE
)

cat(
  "\nPer-gene directional results:\n"
)

print(
  per_gene_three_way_concordance[
    ,
    c(
      "gene",
      "discovery_signature_direction",
      "combined_adjusted_logFC",
      "combined_adjusted_direction",
      "discovery_cohort_logFC",
      "discovery_cohort_direction",
      "validation_cohort_logFC",
      "validation_cohort_direction",
      "all_three_GSE65682_estimates_agree_with_signature"
    )
  ],
  row.names = FALSE
)

# No 30/30 directional-concordance assertion is imposed.

# ------------------------------------------------------------
# 4.9 Export and checkpoint
# ------------------------------------------------------------

write.csv(
  res_val_unadjusted,
  file.path(
    S4_DIR,
    "GSE65682_DE_combined_unadjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  res_val_adjusted,
  file.path(
    S4_DIR,
    "GSE65682_DE_combined_cohort_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  res_val_discovery,
  file.path(
    S4_DIR,
    "GSE65682_DE_discovery_cohort.csv"
  ),
  row.names = FALSE
)

write.csv(
  res_val_validation,
  file.path(
    S4_DIR,
    "GSE65682_DE_validation_cohort.csv"
  ),
  row.names = FALSE
)

write.csv(
  per_gene_three_way_concordance,
  file.path(
    S4_DIR,
    "GSE65682_per_gene_three_way_concordance.csv"
  ),
  row.names = FALSE
)

write.csv(
  three_way_concordance_summary,
  file.path(
    S4_DIR,
    "GSE65682_three_way_concordance_summary.csv"
  ),
  row.names = FALSE
)

de_checkpoint <- list(
  res_val_unadjusted =
    res_val_unadjusted,
  res_val_adjusted =
    res_val_adjusted,
  res_val_discovery =
    res_val_discovery,
  res_val_validation =
    res_val_validation,
  per_gene_three_way_concordance =
    per_gene_three_way_concordance,
  three_way_concordance_summary =
    three_way_concordance_summary,
  n_val_up = n_val_up,
  n_val_down = n_val_down,
  adjusted_design =
    fit_val_adjusted$design,
  adjusted_contrast =
    fit_val_adjusted$contrast
)

saveRDS(
  de_checkpoint,
  file.path(
    S4_DIR,
    "checkpoint_step4_de.rds"
  )
)

cat(
  "\nSTEP 4 completed successfully.\n"
)

cat(
  "Formal downstream GSE65682 result: res_val_adjusted.\n"
)

cat(
  "No historical 36/9 or 30/30 assertion was imposed.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 5: Within-cohort outcome-label permutation test
# ============================================================

cat(
  "\n========== STEP 5: WITHIN-COHORT PERMUTATION ==========\n"
)

# ------------------------------------------------------------
# 5.1 Required-object audit
# ------------------------------------------------------------

required_step5_objects <- c(
  "gmat_val",
  "group_val",
  "cohort_val",
  "signature_adjusted",
  "fisher_out",
  "res_val_adjusted",
  "S4_DIR"
)

missing_step5_objects <- required_step5_objects[
  !vapply(
    required_step5_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (
  length(missing_step5_objects) > 0L
) {
  stop(sprintf(
    "Missing objects required for Step 5: %s",
    paste(
      missing_step5_objects,
      collapse = ", "
    )
  ))
}

genes_tested <- signature_adjusted

disc_direction_tested <- as.character(
  fisher_out$direction[
    match(
      genes_tested,
      fisher_out$gene
    )
  ]
)

disc_sign_tested <- ifelse(
  disc_direction_tested == "Up",
  1L,
  ifelse(
    disc_direction_tested == "Down",
    -1L,
    NA_integer_
  )
)

observed_result_index <- match(
  genes_tested,
  res_val_adjusted$gene
)

if (
  anyNA(observed_result_index) ||
  anyNA(disc_sign_tested)
) {
  stop(
    "The 30 evaluable genes could not be matched to the formal results."
  )
}

observed_logFC <- res_val_adjusted$logFC[
  observed_result_index
]

if (
  any(
    !is.finite(observed_logFC)
  )
) {
  stop(
    "Observed cohort-adjusted logFC contained non-finite values."
  )
}

observed_sign <- sign(
  observed_logFC
)

n_comparable <- length(
  genes_tested
)

actual_concordant <- sum(
  observed_sign ==
    disc_sign_tested
)

actual_concord_rate <-
  actual_concordant /
  n_comparable

up_index <- disc_sign_tested == 1L
down_index <- disc_sign_tested == -1L

n_up_disc <- sum(up_index)
n_down_disc <- sum(down_index)

actual_concord_up <- sum(
  observed_sign[up_index] == 1L
)

actual_concord_down <- sum(
  observed_sign[down_index] == -1L
)

cat(sprintf(
  paste0(
    "Observed cohort-adjusted concordance: ",
    "Up=%d/%d | Down=%d/%d | Overall=%d/%d (%.1f%%)\n"
  ),
  actual_concord_up,
  n_up_disc,
  actual_concord_down,
  n_down_disc,
  actual_concordant,
  n_comparable,
  100 * actual_concord_rate
))

# No 30/30 assertion is imposed.

# ------------------------------------------------------------
# 5.2 Prepare expression and permutation strata
# ------------------------------------------------------------

expr_tested <- gmat_val[
  genes_tested,
  ,
  drop = FALSE
]

group_chr <- as.character(
  group_val
)

cohort_chr <- as.character(
  cohort_val
)

stopifnot(
  nrow(expr_tested) == 30L,
  ncol(expr_tested) == 479L,
  identical(
    rownames(expr_tested),
    genes_tested
  ),
  identical(
    colnames(expr_tested),
    names(group_val)
  ),
  identical(
    names(group_val),
    names(cohort_val)
  ),
  all(
    is.finite(expr_tested)
  ),
  all(
    group_chr %in%
      c(
        "Survivor",
        "NonSurvivor"
      )
  ),
  all(
    cohort_chr %in%
      c(
        "Discovery",
        "Validation"
      )
  )
)

permutation_stratum_counts <- table(
  Cohort = cohort_chr,
  Outcome = group_chr
)

cat(
  "\nOutcome counts retained within each permutation stratum:\n"
)

print(
  permutation_stratum_counts
)

if (
  any(
    permutation_stratum_counts == 0L
  )
) {
  stop(
    "At least one cohort lacks one of the two outcome levels."
  )
}

# ------------------------------------------------------------
# 5.3 Construct the contrast template
# ------------------------------------------------------------

po_template <- factor(
  group_chr,
  levels = c(
    "Survivor",
    "NonSurvivor"
  )
)

pc_template <- factor(
  cohort_chr,
  levels = c(
    "Discovery",
    "Validation"
  )
)

design_template <- stats::model.matrix(
  ~ 0 + po_template + pc_template
)

expected_permutation_columns <- c(
  "po_templateSurvivor",
  "po_templateNonSurvivor",
  "pc_templateValidation"
)

if (
  !identical(
    colnames(design_template),
    expected_permutation_columns
  )
) {
  stop(sprintf(
    "Unexpected permutation-template columns: %s",
    paste(
      colnames(design_template),
      collapse = ", "
    )
  ))
}

if (
  qr(design_template)$rank !=
  ncol(design_template)
) {
  stop(
    "The permutation-template design matrix is not full rank."
  )
}

contrast_template <- limma::makeContrasts(
  po_templateNonSurvivor -
    po_templateSurvivor,
  levels = design_template
)

# ------------------------------------------------------------
# 5.4 Run 10,000 within-cohort permutations
# ------------------------------------------------------------

set.seed(2024L)

n_perm <- 10000L

perm_concordant_count <- integer(
  n_perm
)

permutation_start_time <- as.numeric(
  proc.time()["elapsed"]
)

for (
  permutation_index in
  seq_len(n_perm)
) {
  perm_group <- group_chr
  
  for (
    current_cohort in
    c(
      "Discovery",
      "Validation"
    )
  ) {
    current_index <- which(
      cohort_chr ==
        current_cohort
    )
    
    perm_group[current_index] <- sample(
      group_chr[current_index],
      replace = FALSE
    )
  }
  
  po <- factor(
    perm_group,
    levels = c(
      "Survivor",
      "NonSurvivor"
    )
  )
  
  pc <- factor(
    cohort_chr,
    levels = c(
      "Discovery",
      "Validation"
    )
  )
  
  design_perm <- stats::model.matrix(
    ~ 0 + po + pc
  )
  
  expected_current_columns <- c(
    "poSurvivor",
    "poNonSurvivor",
    "pcValidation"
  )
  
  if (
    !identical(
      colnames(design_perm),
      expected_current_columns
    )
  ) {
    stop(sprintf(
      "Unexpected design columns at permutation %d: %s",
      permutation_index,
      paste(
        colnames(design_perm),
        collapse = ", "
      )
    ))
  }
  
  contrast_perm <- limma::makeContrasts(
    poNonSurvivor -
      poSurvivor,
    levels = design_perm
  )
  
  fit_perm <- limma::lmFit(
    expr_tested,
    design_perm
  )
  
  fit_perm <- limma::contrasts.fit(
    fit_perm,
    contrast_perm
  )
  
  # eBayes is intentionally omitted:
  # it does not alter the estimated contrast coefficients.
  perm_logFC <- as.numeric(
    fit_perm$coefficients[
      ,
      1L
    ]
  )
  
  if (
    any(
      !is.finite(perm_logFC)
    )
  ) {
    stop(sprintf(
      "Non-finite permutation coefficient at iteration %d.",
      permutation_index
    ))
  }
  
  perm_concordant_count[
    permutation_index
  ] <- sum(
    sign(perm_logFC) ==
      disc_sign_tested
  )
  
  if (
    permutation_index %% 1000L ==
    0L
  ) {
    cat(sprintf(
      "Completed %d/%d permutations\n",
      permutation_index,
      n_perm
    ))
  }
}

permutation_elapsed_seconds <-
  as.numeric(
    proc.time()["elapsed"]
  ) -
  permutation_start_time

# ------------------------------------------------------------
# 5.5 Empirical P value and null summary
# ------------------------------------------------------------

perm_concord_rate <-
  perm_concordant_count /
  n_comparable

n_exceed <- sum(
  perm_concordant_count >=
    actual_concordant
)

permutation_empirical_P <-
  (n_exceed + 1) /
  (n_perm + 1)

null_quantiles <- stats::quantile(
  perm_concordant_count,
  probs = c(
    0.025,
    0.25,
    0.50,
    0.75,
    0.975
  ),
  names = FALSE,
  type = 7
)

permutation_summary <- data.frame(
  permutation_scheme =
    "Outcome labels permuted within deposited cohort",
  formal_model =
    "expression ~ mortality + endotype_cohort",
  contrast =
    "NonSurvivor - Survivor",
  seed = 2024L,
  evaluable_genes =
    n_comparable,
  discovery_up_genes =
    n_up_disc,
  discovery_down_genes =
    n_down_disc,
  observed_concordant_up =
    actual_concord_up,
  observed_concordant_down =
    actual_concord_down,
  observed_concordant_total =
    actual_concordant,
  observed_concordance_rate =
    actual_concord_rate,
  permutations =
    n_perm,
  permutations_at_least_as_extreme =
    n_exceed,
  empirical_P =
    permutation_empirical_P,
  null_mean_concordant =
    mean(perm_concordant_count),
  null_SD_concordant =
    stats::sd(perm_concordant_count),
  null_min_concordant =
    min(perm_concordant_count),
  null_q025_concordant =
    null_quantiles[1L],
  null_q25_concordant =
    null_quantiles[2L],
  null_median_concordant =
    null_quantiles[3L],
  null_q75_concordant =
    null_quantiles[4L],
  null_q975_concordant =
    null_quantiles[5L],
  null_max_concordant =
    max(perm_concordant_count),
  elapsed_seconds =
    permutation_elapsed_seconds,
  stringsAsFactors = FALSE
)

permutation_source <- data.frame(
  permutation =
    seq_len(n_perm),
  concordant_genes =
    perm_concordant_count,
  concordance_rate =
    perm_concord_rate,
  stringsAsFactors = FALSE
)

cat(
  "\nWithin-cohort permutation summary:\n"
)

print(
  permutation_summary,
  row.names = FALSE
)

cat(
  "\nNull distribution of concordant-gene counts:\n"
)

print(
  table(
    perm_concordant_count
  )
)

# ------------------------------------------------------------
# 5.6 Export and checkpoint
# ------------------------------------------------------------

write.csv(
  permutation_summary,
  file.path(
    S4_DIR,
    "GSE65682_within_cohort_permutation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  permutation_source,
  file.path(
    S4_DIR,
    "GSE65682_within_cohort_permutation_source_data.csv"
  ),
  row.names = FALSE
)

permutation_checkpoint <- list(
  genes_tested = genes_tested,
  discovery_direction =
    disc_direction_tested,
  observed_logFC =
    observed_logFC,
  observed_concordant =
    actual_concordant,
  permutation_counts =
    perm_concordant_count,
  permutation_summary =
    permutation_summary,
  seed = 2024L
)

saveRDS(
  permutation_checkpoint,
  file.path(
    S4_DIR,
    "checkpoint_step5_permutation.rds"
  )
)

cat(
  "\nSTEP 5 completed successfully.\n"
)

cat(
  "Permutation scheme: outcomes reshuffled within cohort.\n"
)

cat(
  "Every permutation effect passed through the mortality contrast.\n"
)

cat(
  "No eBayes call was used inside the permutation loop.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 6: Adjusted external effect-size assessment
# ============================================================

cat(
  "\n========== STEP 6: ADJUSTED EXTERNAL EFFECT ASSESSMENT ==========\n"
)

# ------------------------------------------------------------
# 6.1 Required-object audit
# ------------------------------------------------------------

required_step6_objects <- c(
  "signature_adjusted",
  "fisher_out",
  "res_val_unadjusted",
  "res_val_adjusted",
  "res_val_discovery",
  "res_val_validation",
  "S4_DIR"
)

missing_step6_objects <- required_step6_objects[
  !vapply(
    required_step6_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (
  length(missing_step6_objects) > 0L
) {
  stop(sprintf(
    "Missing objects required for Step 6: %s",
    paste(
      missing_step6_objects,
      collapse = ", "
    )
  ))
}

required_discovery_columns <- c(
  "gene",
  "direction",
  "re_estimate"
)

missing_discovery_columns <- setdiff(
  required_discovery_columns,
  colnames(fisher_out)
)

if (
  length(missing_discovery_columns) > 0L
) {
  stop(sprintf(
    "fisher_out is missing columns: %s",
    paste(
      missing_discovery_columns,
      collapse = ", "
    )
  ))
}

# ------------------------------------------------------------
# 6.2 Match the 30 evaluable genes
# ------------------------------------------------------------

gene_index_discovery <- match(
  signature_adjusted,
  fisher_out$gene
)

gene_index_unadjusted <- match(
  signature_adjusted,
  res_val_unadjusted$gene
)

gene_index_adjusted <- match(
  signature_adjusted,
  res_val_adjusted$gene
)

gene_index_deposited_discovery <- match(
  signature_adjusted,
  res_val_discovery$gene
)

gene_index_deposited_validation <- match(
  signature_adjusted,
  res_val_validation$gene
)

all_match_indices <- list(
  gene_index_discovery,
  gene_index_unadjusted,
  gene_index_adjusted,
  gene_index_deposited_discovery,
  gene_index_deposited_validation
)

if (
  any(
    vapply(
      all_match_indices,
      anyNA,
      logical(1)
    )
  )
) {
  stop(
    "At least one evaluable gene was missing from an effect-result table."
  )
}

module_step6 <- rep(
  "Unclassified",
  length(signature_adjusted)
)

if (
  exists("module_map") &&
  !is.null(
    names(module_map)
  )
) {
  matched_module <- as.character(
    module_map[
      signature_adjusted
    ]
  )
  
  module_step6[
    !is.na(matched_module)
  ] <- matched_module[
    !is.na(matched_module)
  ]
}

effect_source_step6 <- data.frame(
  gene = signature_adjusted,
  module = module_step6,
  discovery_direction = as.character(
    fisher_out$direction[
      gene_index_discovery
    ]
  ),
  discovery_random_effects_logFC =
    fisher_out$re_estimate[
      gene_index_discovery
    ],
  combined_unadjusted_logFC =
    res_val_unadjusted$logFC[
      gene_index_unadjusted
    ],
  combined_unadjusted_t =
    res_val_unadjusted$t[
      gene_index_unadjusted
    ],
  combined_unadjusted_P =
    res_val_unadjusted$P.Value[
      gene_index_unadjusted
    ],
  combined_unadjusted_FDR =
    res_val_unadjusted$adj.P.Val[
      gene_index_unadjusted
    ],
  combined_adjusted_logFC =
    res_val_adjusted$logFC[
      gene_index_adjusted
    ],
  combined_adjusted_t =
    res_val_adjusted$t[
      gene_index_adjusted
    ],
  combined_adjusted_P =
    res_val_adjusted$P.Value[
      gene_index_adjusted
    ],
  combined_adjusted_FDR =
    res_val_adjusted$adj.P.Val[
      gene_index_adjusted
    ],
  deposited_discovery_logFC =
    res_val_discovery$logFC[
      gene_index_deposited_discovery
    ],
  deposited_discovery_t =
    res_val_discovery$t[
      gene_index_deposited_discovery
    ],
  deposited_discovery_P =
    res_val_discovery$P.Value[
      gene_index_deposited_discovery
    ],
  deposited_discovery_FDR =
    res_val_discovery$adj.P.Val[
      gene_index_deposited_discovery
    ],
  deposited_validation_logFC =
    res_val_validation$logFC[
      gene_index_deposited_validation
    ],
  deposited_validation_t =
    res_val_validation$t[
      gene_index_deposited_validation
    ],
  deposited_validation_P =
    res_val_validation$P.Value[
      gene_index_deposited_validation
    ],
  deposited_validation_FDR =
    res_val_validation$adj.P.Val[
      gene_index_deposited_validation
    ],
  stringsAsFactors = FALSE
)

numeric_effect_columns_step6 <- setdiff(
  colnames(effect_source_step6),
  c(
    "gene",
    "module",
    "discovery_direction"
  )
)

if (
  any(
    !is.finite(
      as.matrix(
        effect_source_step6[
          ,
          numeric_effect_columns_step6,
          drop = FALSE
        ]
      )
    )
  )
) {
  stop(
    "Non-finite statistics were detected in the 30-gene effect table."
  )
}

# ------------------------------------------------------------
# 6.3 Direction and magnitude annotations
# ------------------------------------------------------------

direction_from_effect_step6 <- function(x) {
  ifelse(
    x > 0,
    "Up",
    ifelse(
      x < 0,
      "Down",
      "Zero"
    )
  )
}

effect_source_step6$
  combined_unadjusted_direction <-
  direction_from_effect_step6(
    effect_source_step6$
      combined_unadjusted_logFC
  )

effect_source_step6$
  combined_adjusted_direction <-
  direction_from_effect_step6(
    effect_source_step6$
      combined_adjusted_logFC
  )

effect_source_step6$
  deposited_discovery_direction <-
  direction_from_effect_step6(
    effect_source_step6$
      deposited_discovery_logFC
  )

effect_source_step6$
  deposited_validation_direction <-
  direction_from_effect_step6(
    effect_source_step6$
      deposited_validation_logFC
  )

effect_source_step6$
  combined_adjusted_directionally_concordant <-
  effect_source_step6$
  combined_adjusted_direction ==
  effect_source_step6$
  discovery_direction

effect_source_step6$
  combined_adjusted_above_0_5 <-
  abs(
    effect_source_step6$
      combined_adjusted_logFC
  ) > 0.5

effect_source_step6$
  adjustment_delta_logFC <-
  effect_source_step6$
  combined_adjusted_logFC -
  effect_source_step6$
  combined_unadjusted_logFC

# ------------------------------------------------------------
# 6.4 Cross-gene metric function
# ------------------------------------------------------------

summarize_cross_gene_effect_step6 <- function(
    analysis_name,
    external_logFC,
    external_P,
    external_FDR
) {
  discovery_logFC <-
    effect_source_step6$
    discovery_random_effects_logFC
  
  complete_index <-
    is.finite(discovery_logFC) &
    is.finite(external_logFC)
  
  x <- discovery_logFC[
    complete_index
  ]
  
  y <- external_logFC[
    complete_index
  ]
  
  current_direction <- direction_from_effect_step6(
    y
  )
  
  reference_direction <-
    effect_source_step6$
    discovery_direction[
      complete_index
    ]
  
  pearson_test <- stats::cor.test(
    x,
    y,
    method = "pearson"
  )
  
  spearman_test <- suppressWarnings(
    stats::cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  ols_fit <- stats::lm(
    y ~ x
  )
  
  ols_coefficients <- stats::coef(
    ols_fit
  )
  
  ols_confidence_interval <- stats::confint(
    ols_fit,
    level = 0.95
  )
  
  data.frame(
    analysis = analysis_name,
    n_genes = length(x),
    directionally_concordant_n = sum(
      current_direction ==
        reference_direction
    ),
    external_abs_logFC_gt_0_5_n = sum(
      abs(y) > 0.5
    ),
    external_nominal_P_lt_0_05_n = sum(
      external_P[
        complete_index
      ] < 0.05
    ),
    external_FDR_lt_0_05_n = sum(
      external_FDR[
        complete_index
      ] < 0.05
    ),
    Pearson_r = unname(
      pearson_test$estimate
    ),
    Pearson_P = pearson_test$p.value,
    Spearman_rho = unname(
      spearman_test$estimate
    ),
    Spearman_P = spearman_test$p.value,
    OLS_intercept = unname(
      ols_coefficients[1L]
    ),
    OLS_slope_external_on_discovery = unname(
      ols_coefficients[2L]
    ),
    OLS_slope_CI95_lower =
      ols_confidence_interval[
        2L,
        1L
      ],
    OLS_slope_CI95_upper =
      ols_confidence_interval[
        2L,
        2L
      ],
    stringsAsFactors = FALSE
  )
}

cross_gene_effect_metrics_step6 <- do.call(
  rbind,
  list(
    summarize_cross_gene_effect_step6(
      analysis_name =
        "Combined unadjusted",
      external_logFC =
        effect_source_step6$
        combined_unadjusted_logFC,
      external_P =
        effect_source_step6$
        combined_unadjusted_P,
      external_FDR =
        effect_source_step6$
        combined_unadjusted_FDR
    ),
    summarize_cross_gene_effect_step6(
      analysis_name =
        "Combined cohort-adjusted formal",
      external_logFC =
        effect_source_step6$
        combined_adjusted_logFC,
      external_P =
        effect_source_step6$
        combined_adjusted_P,
      external_FDR =
        effect_source_step6$
        combined_adjusted_FDR
    ),
    summarize_cross_gene_effect_step6(
      analysis_name =
        "Deposited discovery cohort",
      external_logFC =
        effect_source_step6$
        deposited_discovery_logFC,
      external_P =
        effect_source_step6$
        deposited_discovery_P,
      external_FDR =
        effect_source_step6$
        deposited_discovery_FDR
    ),
    summarize_cross_gene_effect_step6(
      analysis_name =
        "Deposited validation cohort",
      external_logFC =
        effect_source_step6$
        deposited_validation_logFC,
      external_P =
        effect_source_step6$
        deposited_validation_P,
      external_FDR =
        effect_source_step6$
        deposited_validation_FDR
    )
  )
)

# ------------------------------------------------------------
# 6.5 Formal adjusted-result sensitivity analyses
# ------------------------------------------------------------

formal_x <- effect_source_step6$
  discovery_random_effects_logFC

formal_y <- effect_source_step6$
  combined_adjusted_logFC

without_cx3cr1 <- effect_source_step6$gene !=
  "CX3CR1"

upregulated_only <- effect_source_step6$
  discovery_direction == "Up"

pearson_without_cx3cr1 <- stats::cor(
  formal_x[
    without_cx3cr1
  ],
  formal_y[
    without_cx3cr1
  ],
  method = "pearson"
)

pearson_upregulated_only <- stats::cor(
  formal_x[
    upregulated_only
  ],
  formal_y[
    upregulated_only
  ],
  method = "pearson"
)

leave_one_gene_out_pearson <- vapply(
  seq_along(formal_x),
  function(index_to_remove) {
    stats::cor(
      formal_x[
        -index_to_remove
      ],
      formal_y[
        -index_to_remove
      ],
      method = "pearson"
    )
  },
  numeric(1)
)

formal_sensitivity_step6 <- data.frame(
  analysis = c(
    "Pearson excluding CX3CR1",
    "Pearson upregulated genes only",
    "Leave-one-gene-out Pearson minimum",
    "Leave-one-gene-out Pearson maximum",
    "Mean absolute logFC change after cohort adjustment",
    "Maximum absolute logFC change after cohort adjustment",
    "Genes changing direction after cohort adjustment"
  ),
  n_genes = c(
    sum(without_cx3cr1),
    sum(upregulated_only),
    length(formal_x) - 1L,
    length(formal_x) - 1L,
    length(formal_x),
    length(formal_x),
    length(formal_x)
  ),
  estimate = c(
    pearson_without_cx3cr1,
    pearson_upregulated_only,
    min(
      leave_one_gene_out_pearson
    ),
    max(
      leave_one_gene_out_pearson
    ),
    mean(
      abs(
        effect_source_step6$
          adjustment_delta_logFC
      )
    ),
    max(
      abs(
        effect_source_step6$
          adjustment_delta_logFC
      )
    ),
    sum(
      effect_source_step6$
        combined_adjusted_direction !=
        effect_source_step6$
        combined_unadjusted_direction
    )
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 6.6 Print and export
# ------------------------------------------------------------

cat(
  "\nCross-gene effect metrics:\n"
)

print(
  cross_gene_effect_metrics_step6,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nFormal adjusted-result sensitivity analyses:\n"
)

print(
  formal_sensitivity_step6,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nFormal 30-gene adjusted effect table:\n"
)

print(
  effect_source_step6[
    ,
    c(
      "gene",
      "module",
      "discovery_random_effects_logFC",
      "combined_unadjusted_logFC",
      "combined_adjusted_logFC",
      "adjustment_delta_logFC",
      "combined_adjusted_P",
      "combined_adjusted_FDR",
      "combined_adjusted_directionally_concordant",
      "combined_adjusted_above_0_5"
    )
  ],
  row.names = FALSE,
  digits = 6
)

write.csv(
  effect_source_step6,
  file.path(
    S4_DIR,
    "GSE65682_adjusted_external_effect_source_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  cross_gene_effect_metrics_step6,
  file.path(
    S4_DIR,
    "GSE65682_cross_gene_effect_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  formal_sensitivity_step6,
  file.path(
    S4_DIR,
    "GSE65682_adjusted_effect_sensitivity.csv"
  ),
  row.names = FALSE
)

effect_assessment_checkpoint <- list(
  effect_source =
    effect_source_step6,
  cross_gene_metrics =
    cross_gene_effect_metrics_step6,
  formal_sensitivity =
    formal_sensitivity_step6
)

saveRDS(
  effect_assessment_checkpoint,
  file.path(
    S4_DIR,
    "checkpoint_step6_effect_assessment.rds"
  )
)

cat(
  "\nSTEP 6 completed successfully.\n"
)

cat(
  "Formal external effects were taken from res_val_adjusted.\n"
)

cat(
  "OLS direction: external effect regressed on discovery effect.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 7: Cohort-adjusted MARS-stratified differential expression
# ============================================================

cat(
  "\n========== STEP 7: MARS-STRATIFIED DIFFERENTIAL EXPRESSION ==========\n"
)

# ------------------------------------------------------------
# 7.1 Required-object and cell-count audit
# ------------------------------------------------------------

required_step7_objects <- c(
  "gmat_val",
  "group_val",
  "endo_val",
  "cohort_val",
  "signature_adjusted",
  "fisher_out",
  "S4_DIR"
)

missing_step7_objects <- required_step7_objects[
  !vapply(
    required_step7_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (
  length(missing_step7_objects) > 0L
) {
  stop(sprintf(
    "Missing objects required for Step 7: %s",
    paste(
      missing_step7_objects,
      collapse = ", "
    )
  ))
}

mars_levels_step7 <- c(
  "Mars1",
  "Mars2",
  "Mars3",
  "Mars4"
)

mars_factor_step7 <- factor(
  as.character(endo_val),
  levels = mars_levels_step7
)

outcome_factor_step7 <- factor(
  as.character(group_val),
  levels = c(
    "Survivor",
    "NonSurvivor"
  )
)

cohort_factor_step7 <- factor(
  as.character(cohort_val),
  levels = c(
    "Discovery",
    "Validation"
  )
)

stopifnot(
  ncol(gmat_val) == 479L,
  length(mars_factor_step7) == 479L,
  length(outcome_factor_step7) == 479L,
  length(cohort_factor_step7) == 479L,
  !anyNA(mars_factor_step7),
  !anyNA(outcome_factor_step7),
  !anyNA(cohort_factor_step7),
  identical(
    colnames(gmat_val),
    names(group_val)
  ),
  identical(
    names(group_val),
    names(endo_val)
  ),
  identical(
    names(group_val),
    names(cohort_val)
  )
)

mars_cohort_outcome_step7 <- table(
  MARS = mars_factor_step7,
  Cohort = cohort_factor_step7,
  Outcome = outcome_factor_step7
)

cat(
  "\nMARS × cohort × outcome audit:\n"
)

print(
  mars_cohort_outcome_step7
)

if (
  any(
    mars_cohort_outcome_step7 == 0L
  )
) {
  stop(
    paste0(
      "At least one MARS × cohort × outcome cell is empty:\n",
      paste(
        capture.output(
          print(
            mars_cohort_outcome_step7
          )
        ),
        collapse = "\n"
      )
    )
  )
}

# ------------------------------------------------------------
# 7.2 Fit all four MARS strata exactly once
# ------------------------------------------------------------

fit_mars_stratified_de_step7 <- function(
    expression_matrix,
    outcome,
    mars,
    cohort,
    genes
) {
  mars_levels <- c(
    "Mars1",
    "Mars2",
    "Mars3",
    "Mars4"
  )
  
  result_list <- setNames(
    vector(
      "list",
      length(mars_levels)
    ),
    mars_levels
  )
  
  signature_result_list <- setNames(
    vector(
      "list",
      length(mars_levels)
    ),
    mars_levels
  )
  
  logFC_matrix <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(mars_levels),
    dimnames = list(
      genes,
      mars_levels
    )
  )
  
  t_matrix <- logFC_matrix
  P_matrix <- logFC_matrix
  FDR_matrix <- logFC_matrix
  
  design_audit_rows <- vector(
    "list",
    length(mars_levels)
  )
  
  for (
    mars_index in
    seq_along(mars_levels)
  ) {
    current_mars <- mars_levels[
      mars_index
    ]
    
    current_sample_index <- which(
      mars == current_mars
    )
    
    ge <- factor(
      as.character(
        outcome[
          current_sample_index
        ]
      ),
      levels = c(
        "Survivor",
        "NonSurvivor"
      )
    )
    
    ce <- factor(
      as.character(
        cohort[
          current_sample_index
        ]
      ),
      levels = c(
        "Discovery",
        "Validation"
      )
    )
    
    current_counts <- table(
      Cohort = ce,
      Outcome = ge
    )
    
    if (
      anyNA(ge) ||
      anyNA(ce) ||
      any(
        current_counts == 0L
      )
    ) {
      stop(sprintf(
        paste0(
          "%s failed the cohort/outcome coverage audit:\n%s"
        ),
        current_mars,
        paste(
          capture.output(
            print(current_counts)
          ),
          collapse = "\n"
        )
      ))
    }
    
    current_expression <- expression_matrix[
      ,
      current_sample_index,
      drop = FALSE
    ]
    
    stopifnot(
      ncol(current_expression) ==
        length(ge),
      ncol(current_expression) ==
        length(ce)
    )
    
    design_e <- stats::model.matrix(
      ~ 0 + ge + ce
    )
    
    expected_design_columns <- c(
      "geSurvivor",
      "geNonSurvivor",
      "ceValidation"
    )
    
    if (
      !identical(
        colnames(design_e),
        expected_design_columns
      )
    ) {
      stop(sprintf(
        "%s produced unexpected design columns: %s",
        current_mars,
        paste(
          colnames(design_e),
          collapse = ", "
        )
      ))
    }
    
    design_rank <- qr(
      design_e
    )$rank
    
    if (
      design_rank !=
      ncol(design_e)
    ) {
      stop(sprintf(
        paste0(
          "%s design was not full rank: ",
          "rank=%d, columns=%d.\n%s"
        ),
        current_mars,
        design_rank,
        ncol(design_e),
        paste(
          capture.output(
            print(current_counts)
          ),
          collapse = "\n"
        )
      ))
    }
    
    contrast_e <- limma::makeContrasts(
      geNonSurvivor -
        geSurvivor,
      levels = design_e
    )
    
    fit_e <- limma::lmFit(
      current_expression,
      design_e
    )
    
    fit_e <- limma::contrasts.fit(
      fit_e,
      contrast_e
    )
    
    fit_e <- limma::eBayes(
      fit_e,
      trend = TRUE
    )
    
    result_e <- limma::topTable(
      fit_e,
      coef = 1L,
      number = Inf,
      sort.by = "none"
    )
    
    result_e$gene <- rownames(
      result_e
    )
    
    genome_nonfinite <- !is.finite(
      result_e$logFC
    ) |
      !is.finite(
        result_e$t
      )
    
    if (
      any(genome_nonfinite)
    ) {
      warning(sprintf(
        "%s produced non-finite logFC/t for %d genome-wide genes.",
        current_mars,
        sum(genome_nonfinite)
      ))
    }
    
    gene_index <- match(
      genes,
      result_e$gene
    )
    
    if (
      anyNA(gene_index)
    ) {
      stop(sprintf(
        "%s was missing signature genes: %s",
        current_mars,
        paste(
          genes[
            is.na(gene_index)
          ],
          collapse = ", "
        )
      ))
    }
    
    signature_result_e <- result_e[
      gene_index,
      ,
      drop = FALSE
    ]
    
    signature_numeric <- as.matrix(
      signature_result_e[
        ,
        c(
          "logFC",
          "t",
          "P.Value",
          "adj.P.Val"
        ),
        drop = FALSE
      ]
    )
    
    if (
      any(
        !is.finite(signature_numeric)
      )
    ) {
      stop(sprintf(
        "%s produced non-finite signature-gene statistics.",
        current_mars
      ))
    }
    
    logFC_matrix[
      ,
      current_mars
    ] <- signature_result_e$logFC
    
    t_matrix[
      ,
      current_mars
    ] <- signature_result_e$t
    
    P_matrix[
      ,
      current_mars
    ] <- signature_result_e$P.Value
    
    FDR_matrix[
      ,
      current_mars
    ] <- signature_result_e$adj.P.Val
    
    result_list[[current_mars]] <-
      result_e
    
    signature_result_list[[current_mars]] <-
      signature_result_e
    
    design_audit_rows[[mars_index]] <- data.frame(
      MARS = current_mars,
      samples =
        length(current_sample_index),
      survivors =
        sum(ge == "Survivor"),
      non_survivors =
        sum(ge == "NonSurvivor"),
      discovery_survivors =
        as.integer(
          current_counts[
            "Discovery",
            "Survivor"
          ]
        ),
      discovery_non_survivors =
        as.integer(
          current_counts[
            "Discovery",
            "NonSurvivor"
          ]
        ),
      validation_survivors =
        as.integer(
          current_counts[
            "Validation",
            "Survivor"
          ]
        ),
      validation_non_survivors =
        as.integer(
          current_counts[
            "Validation",
            "NonSurvivor"
          ]
        ),
      design_rank =
        design_rank,
      design_columns =
        ncol(design_e),
      stringsAsFactors = FALSE
    )
    
    cat(sprintf(
      paste0(
        "%s completed: n=%d | S=%d | NS=%d | ",
        "design rank=%d/%d\n"
      ),
      current_mars,
      length(current_sample_index),
      sum(ge == "Survivor"),
      sum(ge == "NonSurvivor"),
      design_rank,
      ncol(design_e)
    ))
  }
  
  design_audit <- do.call(
    rbind,
    design_audit_rows
  )
  
  stopifnot(
    identical(
      rownames(logFC_matrix),
      genes
    ),
    identical(
      colnames(logFC_matrix),
      mars_levels
    ),
    identical(
      dimnames(t_matrix),
      dimnames(logFC_matrix)
    ),
    identical(
      dimnames(P_matrix),
      dimnames(logFC_matrix)
    ),
    identical(
      dimnames(FDR_matrix),
      dimnames(logFC_matrix)
    ),
    all(
      is.finite(logFC_matrix)
    ),
    all(
      is.finite(t_matrix)
    ),
    all(
      is.finite(P_matrix)
    ),
    all(
      is.finite(FDR_matrix)
    )
  )
  
  list(
    logFC = logFC_matrix,
    t = t_matrix,
    P = P_matrix,
    FDR = FDR_matrix,
    full_results = result_list,
    signature_results =
      signature_result_list,
    design_audit = design_audit
  )
}

mars_de_step7 <- fit_mars_stratified_de_step7(
  expression_matrix = gmat_val,
  outcome = outcome_factor_step7,
  mars = mars_factor_step7,
  cohort = cohort_factor_step7,
  genes = signature_adjusted
)

# Canonical downstream objects: do not recompute later.
per_endo_fc <- mars_de_step7$logFC
per_endo_t <- mars_de_step7$t
per_endo_P <- mars_de_step7$P
per_endo_FDR <- mars_de_step7$FDR
per_endo_results <- mars_de_step7$full_results
per_endo_signature_results <-
  mars_de_step7$signature_results

endo_levels <- mars_levels_step7

mars3_logfc <- per_endo_fc[
  ,
  "Mars3"
]

# ------------------------------------------------------------
# 7.3 Directional concordance summary
# ------------------------------------------------------------

discovery_direction_step7 <- as.character(
  fisher_out$direction[
    match(
      signature_adjusted,
      fisher_out$gene
    )
  ]
)

discovery_sign_step7 <- ifelse(
  discovery_direction_step7 == "Up",
  1L,
  ifelse(
    discovery_direction_step7 == "Down",
    -1L,
    NA_integer_
  )
)

if (
  anyNA(
    discovery_sign_step7
  )
) {
  stop(
    "Invalid discovery directions in Step 7."
  )
}

mars_outcome_counts_step7 <- table(
  MARS = mars_factor_step7,
  Outcome = outcome_factor_step7
)

endo_results <- do.call(
  rbind,
  lapply(
    mars_levels_step7,
    function(current_mars) {
      current_agreement <-
        sign(
          per_endo_fc[
            ,
            current_mars
          ]
        ) ==
        discovery_sign_step7
      
      data.frame(
        Endotype = current_mars,
        S = as.integer(
          mars_outcome_counts_step7[
            current_mars,
            "Survivor"
          ]
        ),
        NS = as.integer(
          mars_outcome_counts_step7[
            current_mars,
            "NonSurvivor"
          ]
        ),
        agree = sum(
          current_agreement
        ),
        total = length(
          current_agreement
        ),
        pct = round(
          mean(
            current_agreement
          ) * 100,
          1
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

# ------------------------------------------------------------
# 7.4 Long-form signature results
# ------------------------------------------------------------

module_step7 <- rep(
  "Unclassified",
  length(signature_adjusted)
)

if (
  exists("module_map") &&
  !is.null(
    names(module_map)
  )
) {
  matched_module_step7 <- as.character(
    module_map[
      signature_adjusted
    ]
  )
  
  module_step7[
    !is.na(matched_module_step7)
  ] <- matched_module_step7[
    !is.na(matched_module_step7)
  ]
}

per_mars_signature_long <- do.call(
  rbind,
  lapply(
    mars_levels_step7,
    function(current_mars) {
      current_table <-
        per_endo_signature_results[[current_mars]]
      
      current_direction <- ifelse(
        current_table$logFC > 0,
        "Up",
        ifelse(
          current_table$logFC < 0,
          "Down",
          "Zero"
        )
      )
      
      data.frame(
        MARS = current_mars,
        gene = signature_adjusted,
        module = module_step7,
        discovery_direction =
          discovery_direction_step7,
        logFC = current_table$logFC,
        AveExpr = current_table$AveExpr,
        t = current_table$t,
        P = current_table$P.Value,
        FDR = current_table$adj.P.Val,
        B = current_table$B,
        MARS_direction =
          current_direction,
        directionally_concordant =
          current_direction ==
          discovery_direction_step7,
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  per_mars_signature_long
) <- NULL

# ------------------------------------------------------------
# 7.5 Print results
# ------------------------------------------------------------

cat(
  "\nPer-MARS design audit:\n"
)

print(
  mars_de_step7$design_audit,
  row.names = FALSE
)

cat(
  "\nPer-MARS directional concordance:\n"
)

print(
  endo_results,
  row.names = FALSE
)

cat(
  "\nPer-MARS 30-gene logFC matrix:\n"
)

print(
  per_endo_fc,
  digits = 5
)

# ------------------------------------------------------------
# 7.6 Export
# ------------------------------------------------------------

write.csv(
  mars_de_step7$design_audit,
  file.path(
    S4_DIR,
    "GSE65682_per_MARS_design_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  endo_results,
  file.path(
    S4_DIR,
    "GSE65682_per_MARS_direction_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = rownames(per_endo_fc),
    per_endo_fc,
    check.names = FALSE
  ),
  file.path(
    S4_DIR,
    "GSE65682_per_MARS_logFC_matrix.csv"
  ),
  row.names = FALSE
)

write.csv(
  per_mars_signature_long,
  file.path(
    S4_DIR,
    "GSE65682_per_MARS_signature_statistics.csv"
  ),
  row.names = FALSE
)

for (
  current_mars in
  mars_levels_step7
) {
  write.csv(
    per_endo_results[[current_mars]],
    file.path(
      S4_DIR,
      paste0(
        "GSE65682_",
        current_mars,
        "_cohort_adjusted_DE_all_genes.csv"
      )
    ),
    row.names = FALSE
  )
}

mars_de_checkpoint <- list(
  per_endo_fc =
    per_endo_fc,
  per_endo_t =
    per_endo_t,
  per_endo_P =
    per_endo_P,
  per_endo_FDR =
    per_endo_FDR,
  per_endo_results =
    per_endo_results,
  per_endo_signature_results =
    per_endo_signature_results,
  endo_results =
    endo_results,
  per_mars_signature_long =
    per_mars_signature_long,
  design_audit =
    mars_de_step7$design_audit,
  mars_cohort_outcome_counts =
    mars_cohort_outcome_step7
)

saveRDS(
  mars_de_checkpoint,
  file.path(
    S4_DIR,
    "checkpoint_step7_per_MARS_DE.rds"
  )
)

cat(
  "\nSTEP 7 completed successfully.\n"
)

cat(
  "All four MARS strata were retained.\n"
)

cat(
  "Each MARS model adjusted for deposited cohort.\n"
)

cat(
  "per_endo_fc is the canonical downstream logFC matrix.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 8: Cohort-adjusted MARS-by-module interaction
# ============================================================

cat(
  "\n========== STEP 8: COHORT-ADJUSTED MODULE × MARS INTERACTION ==========\n"
)

# ------------------------------------------------------------
# 8.1 Required-object audit
# ------------------------------------------------------------

required_step8_objects <- c(
  "gmat_val",
  "group_val",
  "endo_val",
  "cohort_val",
  "MODULE_CELL_CYCLE",
  "MODULE_NEUTRO",
  "REVISION_DIR"
)

missing_step8_objects <- required_step8_objects[
  !vapply(
    required_step8_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (
  length(missing_step8_objects) > 0L
) {
  stop(sprintf(
    "Missing objects required for Step 8: %s",
    paste(
      missing_step8_objects,
      collapse = ", "
    )
  ))
}

MARS_INTERACTION_DIR <- file.path(
  REVISION_DIR,
  "09_MARS_module_interaction"
)

dir.create(
  MARS_INTERACTION_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

mars_levels_step8 <- c(
  "Mars1",
  "Mars2",
  "Mars3",
  "Mars4"
)

module_sets_step8 <- list(
  "Cell cycle / proliferation" =
    MODULE_CELL_CYCLE,
  "Neutrophil degranulation" =
    MODULE_NEUTRO
)

stopifnot(
  length(MODULE_CELL_CYCLE) == 11L,
  length(MODULE_NEUTRO) == 10L,
  ncol(gmat_val) == 479L,
  identical(
    colnames(gmat_val),
    names(group_val)
  ),
  identical(
    names(group_val),
    names(endo_val)
  ),
  identical(
    names(group_val),
    names(cohort_val)
  )
)

# ------------------------------------------------------------
# 8.2 Module-score function
# ------------------------------------------------------------

calculate_module_score_step8 <- function(
    expression_matrix,
    module_genes
) {
  detected_genes <- intersect(
    module_genes,
    rownames(expression_matrix)
  )
  
  missing_genes <- setdiff(
    module_genes,
    detected_genes
  )
  
  if (
    length(missing_genes) > 0L
  ) {
    stop(sprintf(
      "Module genes missing from GSE65682: %s",
      paste(
        missing_genes,
        collapse = ", "
      )
    ))
  }
  
  module_matrix <- expression_matrix[
    detected_genes,
    ,
    drop = FALSE
  ]
  
  gene_sd <- apply(
    module_matrix,
    1,
    stats::sd
  )
  
  invalid_gene <- !is.finite(
    gene_sd
  ) |
    gene_sd <= 0
  
  if (
    any(invalid_gene)
  ) {
    stop(sprintf(
      "Invalid module-gene variance: %s",
      paste(
        names(gene_sd)[
          invalid_gene
        ],
        collapse = ", "
      )
    ))
  }
  
  gene_z <- t(
    scale(
      t(module_matrix),
      center = TRUE,
      scale = TRUE
    )
  )
  
  if (
    any(
      !is.finite(gene_z)
    )
  ) {
    stop(
      "Module gene-level z-score matrix contains non-finite values."
    )
  }
  
  module_score <- colMeans(
    gene_z
  )
  
  if (
    any(
      !is.finite(module_score)
    )
  ) {
    stop(
      "Module score contains non-finite values."
    )
  }
  
  list(
    score = module_score,
    detected_genes = detected_genes,
    missing_genes = missing_genes
  )
}

# ------------------------------------------------------------
# 8.3 Fit one cohort-adjusted interaction model
# ------------------------------------------------------------

fit_module_interaction_step8 <- function(
    module_name,
    module_genes
) {
  score_object <- calculate_module_score_step8(
    expression_matrix = gmat_val,
    module_genes = module_genes
  )
  
  analysis_data <- data.frame(
    sample_id = colnames(gmat_val),
    outcome = as.integer(
      as.character(group_val) ==
        "NonSurvivor"
    ),
    outcome_label = factor(
      as.character(group_val),
      levels = c(
        "Survivor",
        "NonSurvivor"
      )
    ),
    MARS = factor(
      as.character(endo_val),
      levels = mars_levels_step8
    ),
    cohort = factor(
      as.character(cohort_val),
      levels = c(
        "Discovery",
        "Validation"
      )
    ),
    module_score =
      score_object$score,
    stringsAsFactors = FALSE
  )
  
  analysis_data$score_z <- as.numeric(
    scale(
      analysis_data$module_score
    )
  )
  
  required_analysis_columns <- c(
    "outcome",
    "MARS",
    "cohort",
    "score_z"
  )
  
  if (
    any(
      !complete.cases(
        analysis_data[
          ,
          required_analysis_columns,
          drop = FALSE
        ]
      )
    )
  ) {
    stop(sprintf(
      "%s contains incomplete model data.",
      module_name
    ))
  }
  
  stopifnot(
    nrow(analysis_data) == 479L,
    nlevels(analysis_data$MARS) == 4L,
    nlevels(analysis_data$cohort) == 2L,
    all(
      is.finite(
        analysis_data$score_z
      )
    )
  )
  
  analysis_counts <- table(
    MARS = analysis_data$MARS,
    Cohort = analysis_data$cohort,
    Outcome = analysis_data$outcome_label
  )
  
  if (
    any(
      analysis_counts == 0L
    )
  ) {
    stop(sprintf(
      paste0(
        "%s has an empty MARS × cohort × outcome cell:\n%s"
      ),
      module_name,
      paste(
        capture.output(
          print(analysis_counts)
        ),
        collapse = "\n"
      )
    ))
  }
  
  main_effect_model <- stats::glm(
    outcome ~
      score_z +
      MARS +
      cohort,
    family = stats::binomial(),
    data = analysis_data
  )
  
  interaction_model <- stats::glm(
    outcome ~
      score_z * MARS +
      cohort,
    family = stats::binomial(),
    data = analysis_data
  )
  
  if (
    !main_effect_model$converged ||
    !interaction_model$converged
  ) {
    stop(sprintf(
      "%s logistic model failed to converge.",
      module_name
    ))
  }
  
  interaction_lrt <- stats::anova(
    main_effect_model,
    interaction_model,
    test = "LRT"
  )
  
  if (
    nrow(interaction_lrt) != 2L ||
    !is.finite(
      interaction_lrt$Deviance[2L]
    ) ||
    !is.finite(
      interaction_lrt$`Pr(>Chi)`[2L]
    )
  ) {
    stop(sprintf(
      "%s produced an invalid interaction LRT.",
      module_name
    ))
  }
  
  global_result <- data.frame(
    module = module_name,
    submitted_genes =
      length(module_genes),
    detected_genes =
      length(
        score_object$detected_genes
      ),
    analyzed_samples =
      nrow(analysis_data),
    non_survivors =
      sum(
        analysis_data$outcome == 1L
      ),
    survivors =
      sum(
        analysis_data$outcome == 0L
      ),
    interaction_df =
      interaction_lrt$Df[2L],
    interaction_chisq =
      interaction_lrt$Deviance[2L],
    interaction_P =
      interaction_lrt$`Pr(>Chi)`[2L],
    main_model_converged =
      main_effect_model$converged,
    interaction_model_converged =
      interaction_model$converged,
    AIC_main_effect =
      stats::AIC(
        main_effect_model
      ),
    AIC_interaction =
      stats::AIC(
        interaction_model
      ),
    model_formula_main =
      "outcome ~ score_z + MARS + cohort",
    model_formula_interaction =
      "outcome ~ score_z * MARS + cohort",
    stringsAsFactors = FALSE
  )
  
  # --------------------------------------------------------
  # Model-derived score slope within each MARS
  # --------------------------------------------------------
  
  beta <- stats::coef(
    interaction_model
  )
  
  variance_matrix <- stats::vcov(
    interaction_model
  )
  
  if (
    any(
      !is.finite(beta)
    ) ||
    any(
      !is.finite(variance_matrix)
    )
  ) {
    stop(sprintf(
      "%s contains non-finite model coefficients.",
      module_name
    ))
  }
  
  model_terms <- stats::delete.response(
    stats::terms(
      interaction_model
    )
  )
  
  endotype_slopes <- do.call(
    rbind,
    lapply(
      mars_levels_step8,
      function(current_mars) {
        new_zero <- data.frame(
          score_z = 0,
          MARS = factor(
            current_mars,
            levels =
              mars_levels_step8
          ),
          cohort = factor(
            "Discovery",
            levels = c(
              "Discovery",
              "Validation"
            )
          )
        )
        
        new_one <- data.frame(
          score_z = 1,
          MARS = factor(
            current_mars,
            levels =
              mars_levels_step8
          ),
          cohort = factor(
            "Discovery",
            levels = c(
              "Discovery",
              "Validation"
            )
          )
        )
        
        matrix_zero <- stats::model.matrix(
          model_terms,
          data = new_zero
        )
        
        matrix_one <- stats::model.matrix(
          model_terms,
          data = new_one
        )
        
        if (
          !all(
            names(beta) %in%
            colnames(matrix_zero)
          ) ||
          !all(
            names(beta) %in%
            colnames(matrix_one)
          )
        ) {
          stop(sprintf(
            "%s/%s model-matrix columns did not match coefficients.",
            module_name,
            current_mars
          ))
        }
        
        contrast_vector <-
          matrix_one[
            1L,
            names(beta)
          ] -
          matrix_zero[
            1L,
            names(beta)
          ]
        
        log_or <- sum(
          contrast_vector *
            beta
        )
        
        log_or_variance <- as.numeric(
          t(contrast_vector) %*%
            variance_matrix %*%
            contrast_vector
        )
        
        if (
          !is.finite(log_or_variance) ||
          log_or_variance <= 0
        ) {
          stop(sprintf(
            "%s/%s produced invalid slope variance.",
            module_name,
            current_mars
          ))
        }
        
        log_or_se <- sqrt(
          log_or_variance
        )
        
        z_value <- log_or /
          log_or_se
        
        data.frame(
          module = module_name,
          endotype =
            current_mars,
          log_OR_per_1SD =
            log_or,
          SE =
            log_or_se,
          OR_per_1SD =
            exp(log_or),
          CI95_lower =
            exp(
              log_or -
                1.96 *
                log_or_se
            ),
          CI95_upper =
            exp(
              log_or +
                1.96 *
                log_or_se
            ),
          Wald_Z =
            z_value,
          Wald_P =
            2 *
            stats::pnorm(
              -abs(z_value)
            ),
          cohort_reference_for_contrast =
            "Discovery",
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  # --------------------------------------------------------
  # Descriptive score summaries
  # --------------------------------------------------------
  
  score_summary <- do.call(
    rbind,
    lapply(
      mars_levels_step8,
      function(current_mars) {
        do.call(
          rbind,
          lapply(
            c(
              "Survivor",
              "NonSurvivor"
            ),
            function(current_outcome) {
              current_values <-
                analysis_data$score_z[
                  analysis_data$MARS ==
                    current_mars &
                    analysis_data$outcome_label ==
                    current_outcome
                ]
              
              data.frame(
                module =
                  module_name,
                endotype =
                  current_mars,
                outcome =
                  current_outcome,
                n =
                  length(
                    current_values
                  ),
                mean_score_z =
                  mean(
                    current_values
                  ),
                SD_score_z =
                  stats::sd(
                    current_values
                  ),
                median_score_z =
                  stats::median(
                    current_values
                  ),
                Q1_score_z =
                  unname(
                    stats::quantile(
                      current_values,
                      0.25
                    )
                  ),
                Q3_score_z =
                  unname(
                    stats::quantile(
                      current_values,
                      0.75
                    )
                  ),
                stringsAsFactors = FALSE
              )
            }
          )
        )
      }
    )
  )
  
  gene_table <- data.frame(
    module = module_name,
    gene = module_genes,
    detected_in_GSE65682 =
      module_genes %in%
      score_object$detected_genes,
    stringsAsFactors = FALSE
  )
  
  sample_scores <- data.frame(
    module = module_name,
    analysis_data,
    stringsAsFactors = FALSE
  )
  
  list(
    global = global_result,
    slopes = endotype_slopes,
    score_summary = score_summary,
    gene_table = gene_table,
    sample_scores = sample_scores,
    main_effect_model =
      main_effect_model,
    interaction_model =
      interaction_model
  )
}

# ------------------------------------------------------------
# 8.4 Run both formal module tests
# ------------------------------------------------------------

module_interaction_results_step8 <- lapply(
  names(module_sets_step8),
  function(current_module) {
    fit_module_interaction_step8(
      module_name = current_module,
      module_genes =
        module_sets_step8[[current_module]]
    )
  }
)

names(module_interaction_results_step8) <-
  names(module_sets_step8)

module_interaction_global_step8 <- do.call(
  rbind,
  lapply(
    module_interaction_results_step8,
    function(x) {
      x$global
    }
  )
)

module_interaction_global_step8$
  interaction_FDR <- stats::p.adjust(
    module_interaction_global_step8$
      interaction_P,
    method = "BH"
  )

module_interaction_slopes_step8 <- do.call(
  rbind,
  lapply(
    module_interaction_results_step8,
    function(x) {
      x$slopes
    }
  )
)

module_interaction_slopes_step8$
  Wald_FDR <- stats::p.adjust(
    module_interaction_slopes_step8$
      Wald_P,
    method = "BH"
  )

module_score_summary_step8 <- do.call(
  rbind,
  lapply(
    module_interaction_results_step8,
    function(x) {
      x$score_summary
    }
  )
)

module_gene_table_step8 <- do.call(
  rbind,
  lapply(
    module_interaction_results_step8,
    function(x) {
      x$gene_table
    }
  )
)

module_sample_scores_step8 <- do.call(
  rbind,
  lapply(
    module_interaction_results_step8,
    function(x) {
      x$sample_scores
    }
  )
)

rownames(module_interaction_global_step8) <- NULL
rownames(module_interaction_slopes_step8) <- NULL
rownames(module_score_summary_step8) <- NULL
rownames(module_gene_table_step8) <- NULL
rownames(module_sample_scores_step8) <- NULL

# ------------------------------------------------------------
# 8.5 Print and export
# ------------------------------------------------------------

cat(
  "\nGlobal cohort-adjusted module × MARS tests:\n"
)

print(
  module_interaction_global_step8,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nCohort-adjusted model-derived slopes:\n"
)

print(
  module_interaction_slopes_step8,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nDescriptive module-score summaries:\n"
)

print(
  module_score_summary_step8,
  row.names = FALSE,
  digits = 6
)

write.csv(
  module_interaction_global_step8,
  file.path(
    MARS_INTERACTION_DIR,
    "Table_MARS_module_interaction_global_cohort_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_interaction_slopes_step8,
  file.path(
    MARS_INTERACTION_DIR,
    "Table_MARS_module_interaction_endotype_slopes_cohort_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_score_summary_step8,
  file.path(
    MARS_INTERACTION_DIR,
    "Table_MARS_module_score_summary_cohort_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_gene_table_step8,
  file.path(
    MARS_INTERACTION_DIR,
    "Table_MARS_module_gene_membership.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_sample_scores_step8,
  file.path(
    MARS_INTERACTION_DIR,
    "MARS_module_sample_scores_cohort_adjusted.csv"
  ),
  row.names = FALSE
)

saveRDS(
  module_interaction_results_step8,
  file.path(
    MARS_INTERACTION_DIR,
    "MARS_module_interaction_models_cohort_adjusted.rds"
  )
)

cat(
  "\nSTEP 8 completed successfully.\n"
)

cat(
  "Both interaction models included cohort as a main effect.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 9: Cohort-adjusted PC1 outcome association
#         and descriptive AUC assessment
# ============================================================

cat(
  "\n========== STEP 9: PC1 OUTCOME ASSESSMENT ==========\n"
)

# ------------------------------------------------------------
# 9.1 Required-object and package audit
# ------------------------------------------------------------

required_step9_objects <- c(
  "gmat_val",
  "group_val",
  "endo_val",
  "cohort_val",
  "signature_adjusted",
  "REVISION_DIR"
)

missing_step9_objects <- required_step9_objects[
  !vapply(
    required_step9_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_step9_objects) > 0L) {
  stop(sprintf(
    "Missing objects required for Step 9: %s",
    paste(
      missing_step9_objects,
      collapse = ", "
    )
  ))
}

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop(
    "Package 'pROC' is required but is not installed."
  )
}

PC1_DIR <- file.path(
  REVISION_DIR,
  "10_PC1_outcome_assessment"
)

dir.create(
  PC1_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  ncol(gmat_val) == 479L,
  identical(
    colnames(gmat_val),
    names(group_val)
  ),
  identical(
    names(group_val),
    names(endo_val)
  ),
  identical(
    names(group_val),
    names(cohort_val)
  )
)

signature_genes_step9 <- as.character(
  signature_adjusted
)

stopifnot(
  length(signature_genes_step9) == 30L,
  length(unique(signature_genes_step9)) == 30L,
  all(signature_genes_step9 %in% rownames(gmat_val))
)

# ------------------------------------------------------------
# 9.2 Calculate the signature PC1 once
# ------------------------------------------------------------

signature_expression_step9 <- gmat_val[
  signature_genes_step9,
  ,
  drop = FALSE
]

if (any(!is.finite(signature_expression_step9))) {
  stop(
    "Signature expression matrix contains non-finite values."
  )
}

signature_gene_sd_step9 <- apply(
  signature_expression_step9,
  1L,
  stats::sd
)

if (
  any(!is.finite(signature_gene_sd_step9)) ||
  any(signature_gene_sd_step9 <= 0)
) {
  stop(
    "At least one signature gene has invalid variance."
  )
}

pca_32 <- stats::prcomp(
  t(signature_expression_step9),
  center = TRUE,
  scale. = TRUE
)

stopifnot(
  identical(
    rownames(pca_32$x),
    colnames(gmat_val)
  )
)

sig_score_raw <- as.numeric(
  pca_32$x[, 1L]
)

names(sig_score_raw) <- colnames(gmat_val)

group_character_step9 <- as.character(
  group_val
)

mean_raw_non_survivor <- mean(
  sig_score_raw[
    group_character_step9 == "NonSurvivor"
  ]
)

mean_raw_survivor <- mean(
  sig_score_raw[
    group_character_step9 == "Survivor"
  ]
)

if (
  !is.finite(mean_raw_non_survivor) ||
  !is.finite(mean_raw_survivor) ||
  mean_raw_non_survivor == mean_raw_survivor
) {
  stop(
    "PC1 direction could not be determined."
  )
}

# Orient once: higher PC1 corresponds to higher mean in non-survivors.
pc1_orientation_multiplier <- if (
  mean_raw_non_survivor >
  mean_raw_survivor
) {
  1
} else {
  -1
}

sig_score <- sig_score_raw *
  pc1_orientation_multiplier

names(sig_score) <- colnames(gmat_val)

PC1_z <- as.numeric(
  scale(sig_score)
)

names(PC1_z) <- colnames(gmat_val)

if (any(!is.finite(PC1_z))) {
  stop(
    "Standardized PC1 scores contain non-finite values."
  )
}

stopifnot(
  mean(
    sig_score[
      group_character_step9 == "NonSurvivor"
    ]
  ) >
    mean(
      sig_score[
        group_character_step9 == "Survivor"
      ]
    )
)

pc1_variance_explained_pct <- 100 *
  pca_32$sdev[1L]^2 /
  sum(pca_32$sdev^2)

# ------------------------------------------------------------
# 9.3 Construct the single analysis dataset
# ------------------------------------------------------------

pc1_data <- data.frame(
  Sample = colnames(gmat_val),
  outcome_binary = as.integer(
    group_character_step9 ==
      "NonSurvivor"
  ),
  Outcome = factor(
    group_character_step9,
    levels = c(
      "Survivor",
      "NonSurvivor"
    )
  ),
  Endotype = factor(
    as.character(endo_val),
    levels = c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    )
  ),
  Cohort = factor(
    as.character(cohort_val),
    levels = c(
      "Discovery",
      "Validation"
    )
  ),
  Signature_Score = sig_score,
  PC1_z = PC1_z,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(pc1_data) == 479L,
  length(unique(pc1_data$Sample)) == 479L,
  sum(pc1_data$outcome_binary) == 114L,
  all(complete.cases(pc1_data))
)

# ------------------------------------------------------------
# 9.4 ROC helper functions
# ------------------------------------------------------------

build_roc_step9 <- function(
    outcome_vector,
    predictor_vector
) {
  if (
    length(outcome_vector) !=
    length(predictor_vector)
  ) {
    stop(
      "ROC outcome and predictor lengths differ."
    )
  }
  
  if (
    length(unique(as.character(
      outcome_vector
    ))) != 2L
  ) {
    stop(
      "ROC analysis requires both outcome classes."
    )
  }
  
  if (any(!is.finite(predictor_vector))) {
    stop(
      "ROC predictor contains non-finite values."
    )
  }
  
  pROC::roc(
    response = outcome_vector,
    predictor = predictor_vector,
    levels = c(
      "Survivor",
      "NonSurvivor"
    ),
    direction = "<",
    quiet = TRUE
  )
}

build_auc_row_step9 <- function(
    analysis_name,
    roc_object,
    outcome_vector
) {
  current_ci <- as.numeric(
    pROC::ci.auc(
      roc_object,
      method = "delong"
    )
  )
  
  data.frame(
    analysis = analysis_name,
    n = length(outcome_vector),
    survivors = sum(
      as.character(outcome_vector) ==
        "Survivor"
    ),
    non_survivors = sum(
      as.character(outcome_vector) ==
        "NonSurvivor"
    ),
    AUC = as.numeric(
      pROC::auc(roc_object)
    ),
    AUC_CI95_lower = current_ci[1L],
    AUC_CI95_upper = current_ci[3L],
    direction_definition =
      "Higher score or fitted probability indicates higher mortality",
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 9.5 Combined and cohort-specific descriptive PC1 AUC
# ------------------------------------------------------------

roc_obj <- build_roc_step9(
  outcome_vector = pc1_data$Outcome,
  predictor_vector = pc1_data$PC1_z
)

auc_val <- as.numeric(
  pROC::auc(roc_obj)
)

auc_ci <- as.numeric(
  pROC::ci.auc(
    roc_obj,
    method = "delong"
  )
)

descriptive_auc_rows_step9 <- list(
  Combined = build_auc_row_step9(
    analysis_name = "Combined",
    roc_object = roc_obj,
    outcome_vector = pc1_data$Outcome
  )
)

cohort_roc_objects_step9 <- list()

for (current_cohort in levels(pc1_data$Cohort)) {
  
  current_index <- pc1_data$Cohort == current_cohort
  
  current_roc <- build_roc_step9(
    outcome_vector = pc1_data$Outcome[current_index],
    predictor_vector = pc1_data$PC1_z[current_index]
  )
  
  cohort_roc_objects_step9[[current_cohort]] <- current_roc
  
  descriptive_auc_rows_step9[[current_cohort]] <- build_auc_row_step9(
    analysis_name = current_cohort,
    roc_object = current_roc,
    outcome_vector = pc1_data$Outcome[current_index]
  )
}

pc1_descriptive_auc_step9 <- do.call(
  rbind,
  descriptive_auc_rows_step9
)

rownames(pc1_descriptive_auc_step9) <- NULL

# ------------------------------------------------------------
# 9.6 Cohort-adjusted PC1 association
# ------------------------------------------------------------

model_cohort_only_step9 <- stats::glm(
  outcome_binary ~ Cohort,
  family = stats::binomial(),
  data = pc1_data
)

model_cohort_pc1_step9 <- stats::glm(
  outcome_binary ~ Cohort + PC1_z,
  family = stats::binomial(),
  data = pc1_data
)

if (
  !model_cohort_only_step9$converged ||
  !model_cohort_pc1_step9$converged
) {
  stop(
    "At least one PC1 logistic model failed to converge."
  )
}

pc1_model_lrt_step9 <- stats::anova(
  model_cohort_only_step9,
  model_cohort_pc1_step9,
  test = "LRT"
)

if (
  nrow(pc1_model_lrt_step9) != 2L ||
  !is.finite(
    pc1_model_lrt_step9$
    `Pr(>Chi)`[2L]
  )
) {
  stop(
    "Invalid likelihood-ratio test for PC1."
  )
}

pc1_coefficient_step9 <- summary(
  model_cohort_pc1_step9
)$coefficients["PC1_z", ]

pc1_beta_step9 <- unname(
  pc1_coefficient_step9["Estimate"]
)

pc1_se_step9 <- unname(
  pc1_coefficient_step9["Std. Error"]
)

pc1_wald_z_step9 <- unname(
  pc1_coefficient_step9["z value"]
)

pc1_wald_p_step9 <- unname(
  pc1_coefficient_step9["Pr(>|z|)"]
)

pc1_association_step9 <- data.frame(
  predictor = "PC1_z",
  scaling = "Per 1-SD increase",
  analyzed_samples = nrow(pc1_data),
  survivors = sum(
    pc1_data$Outcome == "Survivor"
  ),
  non_survivors = sum(
    pc1_data$Outcome == "NonSurvivor"
  ),
  beta = pc1_beta_step9,
  SE = pc1_se_step9,
  OR_per_1SD = exp(pc1_beta_step9),
  CI95_lower = exp(
    pc1_beta_step9 -
      1.96 * pc1_se_step9
  ),
  CI95_upper = exp(
    pc1_beta_step9 +
      1.96 * pc1_se_step9
  ),
  Wald_Z = pc1_wald_z_step9,
  Wald_P = pc1_wald_p_step9,
  LRT_df =
    pc1_model_lrt_step9$Df[2L],
  LRT_chisq =
    pc1_model_lrt_step9$Deviance[2L],
  LRT_P =
    pc1_model_lrt_step9$
    `Pr(>Chi)`[2L],
  AIC_cohort_only = stats::AIC(
    model_cohort_only_step9
  ),
  AIC_cohort_plus_PC1 = stats::AIC(
    model_cohort_pc1_step9
  ),
  model_formula_null =
    "outcome ~ cohort",
  model_formula_PC1 =
    "outcome ~ cohort + PC1_z",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 9.7 Descriptive fitted-probability AUC
# ------------------------------------------------------------

pc1_data$Predicted_risk_cohort_only <-
  as.numeric(
    stats::predict(
      model_cohort_only_step9,
      type = "response"
    )
  )

pc1_data$Predicted_risk_cohort_plus_PC1 <-
  as.numeric(
    stats::predict(
      model_cohort_pc1_step9,
      type = "response"
    )
  )

roc_model_cohort_only_step9 <-
  build_roc_step9(
    outcome_vector = pc1_data$Outcome,
    predictor_vector =
      pc1_data$Predicted_risk_cohort_only
  )

roc_model_cohort_pc1_step9 <-
  build_roc_step9(
    outcome_vector = pc1_data$Outcome,
    predictor_vector =
      pc1_data$
      Predicted_risk_cohort_plus_PC1
  )

pc1_model_auc_step9 <- rbind(
  build_auc_row_step9(
    analysis_name =
      "Cohort-only fitted probability",
    roc_object =
      roc_model_cohort_only_step9,
    outcome_vector = pc1_data$Outcome
  ),
  build_auc_row_step9(
    analysis_name =
      "Cohort-plus-PC1 fitted probability",
    roc_object =
      roc_model_cohort_pc1_step9,
    outcome_vector = pc1_data$Outcome
  )
)

rownames(pc1_model_auc_step9) <- NULL

# ------------------------------------------------------------
# 9.8 Reuse the same combined ROC object for Youden metrics
# ------------------------------------------------------------

youden_coordinates_step9 <- pROC::coords(
  roc_obj,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity",
    "youden"
  ),
  transpose = FALSE
)

if (is.null(dim(youden_coordinates_step9))) {
  youden_table_step9 <- as.data.frame(
    as.list(youden_coordinates_step9),
    check.names = FALSE
  )
} else {
  youden_table_step9 <- as.data.frame(
    youden_coordinates_step9,
    check.names = FALSE
  )
}

required_youden_columns_step9 <- c(
  "threshold",
  "sensitivity",
  "specificity"
)

if (
  !all(
    required_youden_columns_step9 %in%
    colnames(youden_table_step9)
  )
) {
  stop(
    "Youden output did not contain the expected columns."
  )
}

youden_table_step9$analysis <-
  "Combined internally derived PC1 score"

youden_table_step9$
  optimal_threshold_index <-
  seq_len(nrow(youden_table_step9))

youden_table_step9 <- youden_table_step9[
  ,
  c(
    "analysis",
    "optimal_threshold_index",
    setdiff(
      colnames(youden_table_step9),
      c(
        "analysis",
        "optimal_threshold_index"
      )
    )
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 9.9 PCA summary and oriented loadings
# ------------------------------------------------------------

pc1_pca_summary_step9 <- data.frame(
  submitted_signature_genes = 32L,
  detected_signature_genes =
    length(signature_genes_step9),
  analyzed_samples = nrow(pc1_data),
  PC1_variance_explained_pct =
    pc1_variance_explained_pct,
  raw_PC1_mean_survivor =
    mean_raw_survivor,
  raw_PC1_mean_non_survivor =
    mean_raw_non_survivor,
  orientation_multiplier =
    pc1_orientation_multiplier,
  oriented_PC1_mean_survivor =
    mean(
      pc1_data$Signature_Score[
        pc1_data$Outcome ==
          "Survivor"
      ]
    ),
  oriented_PC1_mean_non_survivor =
    mean(
      pc1_data$Signature_Score[
        pc1_data$Outcome ==
          "NonSurvivor"
      ]
    ),
  orientation_rule =
    "PC1 sign oriented once so that the non-survivor mean was higher",
  stringsAsFactors = FALSE
)

pc1_loadings_step9 <- data.frame(
  gene = rownames(pca_32$rotation),
  loading_raw =
    pca_32$rotation[, 1L],
  loading_oriented =
    pca_32$rotation[, 1L] *
    pc1_orientation_multiplier,
  stringsAsFactors = FALSE
)

pc1_loadings_step9$absolute_loading_rank <-
  rank(
    -abs(
      pc1_loadings_step9$
        loading_oriented
    ),
    ties.method = "first"
  )

pc1_loadings_step9 <- pc1_loadings_step9[
  order(
    pc1_loadings_step9$
      absolute_loading_rank
  ),
  ,
  drop = FALSE
]

rownames(pc1_loadings_step9) <- NULL

# Legacy-compatible combined AUC summary.
auc_metrics_step9 <- data.frame(
  Metric = c(
    "PC1_variance_explained_pct",
    "AUC",
    "AUC_95CI_low",
    "AUC_95CI_high",
    "Youden_threshold_first",
    "Sensitivity_first",
    "Specificity_first",
    "N_optimal_thresholds",
    "N_total",
    "N_Survivor",
    "N_NonSurvivor"
  ),
  Value = c(
    pc1_variance_explained_pct,
    auc_val,
    auc_ci[1L],
    auc_ci[3L],
    youden_table_step9$threshold[1L],
    youden_table_step9$sensitivity[1L],
    youden_table_step9$specificity[1L],
    nrow(youden_table_step9),
    nrow(pc1_data),
    sum(pc1_data$Outcome == "Survivor"),
    sum(
      pc1_data$Outcome ==
        "NonSurvivor"
    )
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 9.10 Print and export
# ------------------------------------------------------------

cat(
  "\nPC1 calculation summary:\n"
)

print(
  pc1_pca_summary_step9,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nDescriptive PC1 AUCs:\n"
)

print(
  pc1_descriptive_auc_step9,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nCohort-adjusted PC1 association:\n"
)

print(
  pc1_association_step9,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nDescriptive fitted-probability AUCs:\n"
)

print(
  pc1_model_auc_step9,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nYouden coordinate(s), using the same combined ROC object:\n"
)

print(
  youden_table_step9,
  row.names = FALSE,
  digits = 6
)

sample_pc1_export_step9 <- pc1_data[
  ,
  c(
    "Sample",
    "Signature_Score",
    "PC1_z",
    "Outcome",
    "Endotype",
    "Cohort",
    "Predicted_risk_cohort_only",
    "Predicted_risk_cohort_plus_PC1"
  ),
  drop = FALSE
]

write.csv(
  sample_pc1_export_step9,
  file.path(
    PC1_DIR,
    "Table_S7_signature_score_AUC.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc1_pca_summary_step9,
  file.path(
    PC1_DIR,
    "Table_PC1_calculation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc1_loadings_step9,
  file.path(
    PC1_DIR,
    "Table_PC1_oriented_loadings.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc1_descriptive_auc_step9,
  file.path(
    PC1_DIR,
    "Table_PC1_descriptive_AUC_by_cohort.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc1_association_step9,
  file.path(
    PC1_DIR,
    "Table_PC1_cohort_adjusted_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc1_model_auc_step9,
  file.path(
    PC1_DIR,
    "Table_PC1_model_prediction_AUC.csv"
  ),
  row.names = FALSE
)

write.csv(
  youden_table_step9,
  file.path(
    PC1_DIR,
    "Table_PC1_Youden_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  auc_metrics_step9,
  file.path(
    PC1_DIR,
    "Table_AUC_metrics.csv"
  ),
  row.names = FALSE
)

pc1_checkpoint_step9 <- list(
  pca = pca_32,
  orientation_multiplier =
    pc1_orientation_multiplier,
  pc1_data = pc1_data,
  pca_summary =
    pc1_pca_summary_step9,
  loadings =
    pc1_loadings_step9,
  descriptive_auc =
    pc1_descriptive_auc_step9,
  association =
    pc1_association_step9,
  model_auc =
    pc1_model_auc_step9,
  youden =
    youden_table_step9,
  cohort_only_model =
    model_cohort_only_step9,
  cohort_plus_pc1_model =
    model_cohort_pc1_step9,
  combined_pc1_roc = roc_obj,
  cohort_specific_rocs =
    cohort_roc_objects_step9,
  cohort_only_probability_roc =
    roc_model_cohort_only_step9,
  cohort_plus_pc1_probability_roc =
    roc_model_cohort_pc1_step9
)

saveRDS(
  pc1_checkpoint_step9,
  file.path(
    PC1_DIR,
    "checkpoint_step9_PC1_outcome_assessment.rds"
  )
)

cat(
  "\nSTEP 9 completed successfully.\n"
)

cat(
  "PC1 was calculated and oriented once.\n"
)

cat(
  "No figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)



cat(
  "\n========== STEP 10 PRECHECK: LOCAL GSEA RESOURCES ==========\n"
)

project_dir_step10 <- if (
  exists("PROJECT_DIR")
) {
  PROJECT_DIR
} else {
  "/home/sunshine/predicate/test03"
}

object_audit_step10 <- data.frame(
  object = c(
    "meta_df",
    "res_val_adjusted",
    "t_COL_95233",
    "hallmark_list",
    "msig_h"
  ),
  exists = vapply(
    c(
      "meta_df",
      "res_val_adjusted",
      "t_COL_95233",
      "hallmark_list",
      "msig_h"
    ),
    exists,
    logical(1),
    inherits = TRUE
  ),
  stringsAsFactors = FALSE
)

print(
  object_audit_step10,
  row.names = FALSE
)

cat("\nPackage availability:\n")

package_audit_step10 <- data.frame(
  package = c(
    "fgsea",
    "msigdbr"
  ),
  installed = c(
    requireNamespace(
      "fgsea",
      quietly = TRUE
    ),
    requireNamespace(
      "msigdbr",
      quietly = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

print(
  package_audit_step10,
  row.names = FALSE
)

if (exists("meta_df")) {
  cat("\nmeta_df dimensions:\n")
  print(dim(meta_df))
  
  cat("\nCandidate moderated-t columns in meta_df:\n")
  print(
    grep(
      "^t_",
      colnames(meta_df),
      value = TRUE
    )
  )
}

if (exists("res_val_adjusted")) {
  cat("\nres_val_adjusted dimensions:\n")
  print(dim(res_val_adjusted))
  
  cat("\nRequired GSE65682 columns present:\n")
  print(
    c(
      gene = "gene" %in%
        colnames(res_val_adjusted),
      t = "t" %in%
        colnames(res_val_adjusted)
    )
  )
}

if (exists("t_COL_95233")) {
  cat("\nt_COL_95233:\n")
  print(t_COL_95233)
}

if (exists("hallmark_list")) {
  cat("\nhallmark_list pathway count:\n")
  print(length(hallmark_list))
  
  cat("\nFirst Hallmark pathway names:\n")
  print(head(names(hallmark_list)))
}

local_gsea_files_step10 <- list.files(
  project_dir_step10,
  pattern = "hallmark|msigdb|\\.gmt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\nCandidate local Hallmark/MSigDB/GMT files:\n")

if (length(local_gsea_files_step10) == 0L) {
  cat("None found in project directory.\n")
} else {
  print(
    head(
      local_gsea_files_step10,
      50L
    )
  )
}

cat(
  "\nNo online query was executed.\n"
)

cat(
  "====================================================\n"
)
cat(
  "\n========== STEP 10B: MSIGDB LOCAL CACHE AUDIT ==========\n"
)

msigdbr_version_step10 <- as.character(
  utils::packageVersion("msigdbr")
)

msigdbr_package_dir_step10 <- system.file(
  package = "msigdbr"
)

cat("\nInstalled msigdbr version:\n")
print(msigdbr_version_step10)

cat("\nInstalled msigdbr package directory:\n")
print(msigdbr_package_dir_step10)

# ------------------------------------------------------------
# Files bundled inside the installed package
# ------------------------------------------------------------

msigdbr_package_files_step10 <- list.files(
  msigdbr_package_dir_step10,
  recursive = TRUE,
  full.names = TRUE
)

msigdbr_package_candidates_step10 <-
  msigdbr_package_files_step10[
    grepl(
      paste0(
        "msig|hallmark|",
        "\\.gmt(\\.gz)?$|",
        "\\.rds$|\\.rda$|\\.rdata$"
      ),
      msigdbr_package_files_step10,
      ignore.case = TRUE
    )
  ]

cat("\nCandidate data files bundled in msigdbr:\n")

if (
  length(msigdbr_package_candidates_step10) ==
  0L
) {
  cat("None detected.\n")
} else {
  print(
    head(
      msigdbr_package_candidates_step10,
      100L
    )
  )
}

# ------------------------------------------------------------
# Standard msigdbr cache locations
# ------------------------------------------------------------

msigdbr_cache_dirs_step10 <- unique(
  c(
    tryCatch(
      tools::R_user_dir(
        "msigdbr",
        which = "cache"
      ),
      error = function(e) {
        character(0)
      }
    ),
    path.expand(
      "~/.cache/R/msigdbr"
    ),
    path.expand(
      "~/.cache/msigdbr"
    ),
    path.expand(
      "~/.local/share/R/msigdbr"
    )
  )
)

msigdbr_cache_audit_step10 <- data.frame(
  directory =
    msigdbr_cache_dirs_step10,
  exists =
    dir.exists(
      msigdbr_cache_dirs_step10
    ),
  stringsAsFactors = FALSE
)

cat("\nCandidate cache directories:\n")

print(
  msigdbr_cache_audit_step10,
  row.names = FALSE
)

existing_cache_dirs_step10 <-
  msigdbr_cache_dirs_step10[
    dir.exists(
      msigdbr_cache_dirs_step10
    )
  ]

if (
  length(existing_cache_dirs_step10) >
  0L
) {
  msigdbr_cache_files_step10 <- unique(
    unlist(
      lapply(
        existing_cache_dirs_step10,
        function(current_directory) {
          list.files(
            current_directory,
            recursive = TRUE,
            full.names = TRUE
          )
        }
      ),
      use.names = FALSE
    )
  )
} else {
  msigdbr_cache_files_step10 <-
    character(0)
}

cat("\nFiles found in msigdbr cache directories:\n")

if (
  length(msigdbr_cache_files_step10) ==
  0L
) {
  cat("None detected.\n")
} else {
  cache_file_information_step10 <- file.info(
    msigdbr_cache_files_step10
  )
  
  cache_file_table_step10 <- data.frame(
    file =
      msigdbr_cache_files_step10,
    size_bytes =
      cache_file_information_step10$size,
    modified =
      cache_file_information_step10$mtime,
    stringsAsFactors = FALSE
  )
  
  print(
    cache_file_table_step10,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# Inspect relevant installed namespace functions
# without calling them
# ------------------------------------------------------------

msigdbr_namespace_names_step10 <- ls(
  envir = asNamespace("msigdbr"),
  all.names = TRUE
)

cat("\nPotential data-loading functions in msigdbr namespace:\n")

print(
  grep(
    "download|cache|load|msigdb",
    msigdbr_namespace_names_step10,
    value = TRUE,
    ignore.case = TRUE
  )
)

cat(
  "\nNo msigdbr data function was called.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 10C: Load the local MSigDB cache and construct
#           moderated-t ranked statistics
# ============================================================

cat(
  "\n========== STEP 10C: LOCAL HALLMARK AND RANKED STATISTICS ==========\n"
)

required_step10_objects <- c(
  "meta_df",
  "res_val_adjusted",
  "t_COL_95233",
  "REVISION_DIR"
)

missing_step10_objects <- required_step10_objects[
  !vapply(
    required_step10_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_step10_objects) > 0L) {
  stop(sprintf(
    "Missing objects required for Step 10C: %s",
    paste(
      missing_step10_objects,
      collapse = ", "
    )
  ))
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("Package 'fgsea' is unavailable.")
}

GSEA_DIR <- file.path(
  REVISION_DIR,
  "11_Hallmark_GSEA_moderated_t"
)

dir.create(
  GSEA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

MSIGDB_CACHE_FILE <- paste0(
  "/home/sunshine/.cache/R/msigdbr/",
  "msigdb.2025.1.Hs.rds"
)

if (!file.exists(MSIGDB_CACHE_FILE)) {
  stop(sprintf(
    "Local MSigDB cache file was not found: %s",
    MSIGDB_CACHE_FILE
  ))
}

# ------------------------------------------------------------
# 10C.1 Read the existing local cache directly
# ------------------------------------------------------------

msigdb_hs_step10 <- readRDS(
  MSIGDB_CACHE_FILE
)

if (!is.data.frame(msigdb_hs_step10)) {
  stop(
    "The cached MSigDB object was not a data frame."
  )
}

required_msigdb_columns_step10 <- c(
  "db_version",
  "db_target_species",
  "gs_collection",
  "gs_name",
  "db_gene_symbol"
)

missing_msigdb_columns_step10 <- setdiff(
  required_msigdb_columns_step10,
  colnames(msigdb_hs_step10)
)

if (
  length(missing_msigdb_columns_step10) >
  0L
) {
  stop(sprintf(
    "Cached MSigDB file lacks required columns: %s",
    paste(
      missing_msigdb_columns_step10,
      collapse = ", "
    )
  ))
}

msigdb_version_step10 <- unique(
  as.character(
    msigdb_hs_step10$db_version
  )
)

msigdb_species_step10 <- unique(
  as.character(
    msigdb_hs_step10$db_target_species
  )
)

if (
  length(msigdb_version_step10) != 1L ||
  length(msigdb_species_step10) != 1L
) {
  stop(
    "MSigDB version or target species was not unique."
  )
}

cat("\nMSigDB database version:\n")
print(msigdb_version_step10)

cat("\nMSigDB target species:\n")
print(msigdb_species_step10)

# ------------------------------------------------------------
# 10C.2 Construct the Hallmark membership list
# ------------------------------------------------------------

msig_h_step10 <- msigdb_hs_step10[
  as.character(
    msigdb_hs_step10$gs_collection
  ) == "H",
  ,
  drop = FALSE
]

hallmark_membership_step10 <- unique(
  data.frame(
    pathway = trimws(
      as.character(
        msig_h_step10$gs_name
      )
    ),
    gene = trimws(
      as.character(
        msig_h_step10$db_gene_symbol
      )
    ),
    stringsAsFactors = FALSE
  )
)

valid_membership_step10 <-
  !is.na(hallmark_membership_step10$pathway) &
  !is.na(hallmark_membership_step10$gene) &
  nzchar(hallmark_membership_step10$pathway) &
  nzchar(hallmark_membership_step10$gene)

hallmark_membership_step10 <-
  hallmark_membership_step10[
    valid_membership_step10,
    ,
    drop = FALSE
  ]

hallmark_membership_step10$db_version <-
  msigdb_version_step10

hallmark_list <- split(
  hallmark_membership_step10$gene,
  hallmark_membership_step10$pathway
)

hallmark_list <- lapply(
  hallmark_list,
  function(current_genes) {
    sort(
      unique(
        as.character(current_genes)
      )
    )
  }
)

hallmark_list <- hallmark_list[
  sort(names(hallmark_list))
]

hallmark_pathway_sizes_step10 <- data.frame(
  pathway = names(hallmark_list),
  database_gene_count = vapply(
    hallmark_list,
    length,
    integer(1)
  ),
  stringsAsFactors = FALSE
)

if (length(hallmark_list) != 50L) {
  stop(sprintf(
    "Expected 50 Hallmark pathways but obtained %d.",
    length(hallmark_list)
  ))
}

if (
  any(
    hallmark_pathway_sizes_step10$
    database_gene_count <= 0L
  )
) {
  stop(
    "At least one Hallmark pathway contained no genes."
  )
}

cat("\nHallmark pathway count:\n")
print(length(hallmark_list))

cat("\nHallmark database gene-count range:\n")
print(
  range(
    hallmark_pathway_sizes_step10$
      database_gene_count
  )
)

# ------------------------------------------------------------
# 10C.3 Function for moderated-t ranked statistics
# ------------------------------------------------------------

build_ranked_stat_step10 <- function(
    data_frame,
    gene_column,
    statistic_column,
    dataset_name
) {
  if (
    !gene_column %in% colnames(data_frame) ||
    !statistic_column %in%
    colnames(data_frame)
  ) {
    stop(sprintf(
      "%s lacks gene or statistic column.",
      dataset_name
    ))
  }
  
  gene_vector <- trimws(
    as.character(
      data_frame[[gene_column]]
    )
  )
  
  statistic_vector <- suppressWarnings(
    as.numeric(
      data_frame[[statistic_column]]
    )
  )
  
  valid_record <- !is.na(gene_vector) &
    nzchar(gene_vector) &
    is.finite(statistic_vector)
  
  ranking_frame <- data.frame(
    gene = gene_vector[valid_record],
    statistic =
      statistic_vector[valid_record],
    original_order =
      which(valid_record),
    stringsAsFactors = FALSE
  )
  
  duplicate_rows_removed <- sum(
    duplicated(ranking_frame$gene)
  )
  
  if (duplicate_rows_removed > 0L) {
    ranking_frame <- ranking_frame[
      order(
        ranking_frame$gene,
        -abs(ranking_frame$statistic),
        ranking_frame$original_order
      ),
      ,
      drop = FALSE
    ]
    
    ranking_frame <- ranking_frame[
      !duplicated(ranking_frame$gene),
      ,
      drop = FALSE
    ]
  }
  
  ranked_statistic <- ranking_frame$statistic
  
  names(ranked_statistic) <-
    ranking_frame$gene
  
  ranked_statistic <- sort(
    ranked_statistic,
    decreasing = TRUE
  )
  
  if (
    length(ranked_statistic) < 1000L ||
    any(!is.finite(ranked_statistic)) ||
    anyDuplicated(names(ranked_statistic))
  ) {
    stop(sprintf(
      "%s produced an invalid ranked statistic.",
      dataset_name
    ))
  }
  
  audit <- data.frame(
    dataset = dataset_name,
    gene_column = gene_column,
    statistic_column =
      statistic_column,
    ranking_metric =
      "limma moderated t",
    input_rows = nrow(data_frame),
    valid_rows_before_deduplication =
      sum(valid_record),
    invalid_rows_removed =
      sum(!valid_record),
    duplicate_rows_removed =
      duplicate_rows_removed,
    final_ranked_genes =
      length(ranked_statistic),
    positive_statistics =
      sum(ranked_statistic > 0),
    negative_statistics =
      sum(ranked_statistic < 0),
    zero_statistics =
      sum(ranked_statistic == 0),
    minimum_statistic =
      min(ranked_statistic),
    maximum_statistic =
      max(ranked_statistic),
    tied_statistic_entries =
      sum(duplicated(ranked_statistic)),
    stringsAsFactors = FALSE
  )
  
  list(
    stats = ranked_statistic,
    audit = audit
  )
}

# ------------------------------------------------------------
# 10C.4 Build all three rankings
# ------------------------------------------------------------

stopifnot(
  identical(
    t_COL_95233,
    "t_95233 (D01)"
  ),
  "gene" %in% colnames(meta_df),
  "t_272769" %in% colnames(meta_df),
  t_COL_95233 %in% colnames(meta_df),
  "gene" %in%
    colnames(res_val_adjusted),
  "t" %in%
    colnames(res_val_adjusted)
)

rank_object_272769_step10 <-
  build_ranked_stat_step10(
    data_frame = meta_df,
    gene_column = "gene",
    statistic_column = "t_272769",
    dataset_name = "GSE272769"
  )

rank_object_95233_step10 <-
  build_ranked_stat_step10(
    data_frame = meta_df,
    gene_column = "gene",
    statistic_column =
      t_COL_95233,
    dataset_name = "GSE95233_D01"
  )

rank_object_65682_step10 <-
  build_ranked_stat_step10(
    data_frame =
      res_val_adjusted,
    gene_column = "gene",
    statistic_column = "t",
    dataset_name =
      "GSE65682_cohort_adjusted"
  )

rank_272769 <-
  rank_object_272769_step10$stats

rank_95233 <-
  rank_object_95233_step10$stats

rank_65682 <-
  rank_object_65682_step10$stats

gsea_rank_audit_step10 <- rbind(
  rank_object_272769_step10$audit,
  rank_object_95233_step10$audit,
  rank_object_65682_step10$audit
)

rownames(gsea_rank_audit_step10) <- NULL

cat("\nModerated-t ranking audit:\n")

print(
  gsea_rank_audit_step10,
  row.names = FALSE,
  digits = 6
)

# ------------------------------------------------------------
# 10C.5 Cache provenance and export
# ------------------------------------------------------------

cache_file_information_step10 <- file.info(
  MSIGDB_CACHE_FILE
)

msigdb_cache_manifest_step10 <- data.frame(
  msigdbr_package_version =
    as.character(
      utils::packageVersion("msigdbr")
    ),
  MSigDB_database_version =
    msigdb_version_step10,
  database_target_species =
    msigdb_species_step10,
  cache_file =
    MSIGDB_CACHE_FILE,
  cache_size_bytes =
    cache_file_information_step10$size,
  cache_modified =
    as.character(
      cache_file_information_step10$mtime
    ),
  cache_MD5 =
    unname(
      tools::md5sum(
        MSIGDB_CACHE_FILE
      )
    ),
  Hallmark_pathways =
    length(hallmark_list),
  Hallmark_membership_rows =
    nrow(
      hallmark_membership_step10
    ),
  stringsAsFactors = FALSE
)

write.csv(
  msigdb_cache_manifest_step10,
  file.path(
    GSEA_DIR,
    "MSigDB_local_cache_manifest.csv"
  ),
  row.names = FALSE
)

write.csv(
  hallmark_membership_step10,
  file.path(
    GSEA_DIR,
    paste0(
      "Hallmark_gene_set_membership_",
      "MSigDB_2025.1.Hs.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  hallmark_pathway_sizes_step10,
  file.path(
    GSEA_DIR,
    "Hallmark_pathway_database_sizes.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_rank_audit_step10,
  file.path(
    GSEA_DIR,
    "GSEA_moderated_t_rank_audit.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    hallmark_list =
      hallmark_list,
    hallmark_membership =
      hallmark_membership_step10,
    pathway_sizes =
      hallmark_pathway_sizes_step10,
    rank_272769 =
      rank_272769,
    rank_95233 =
      rank_95233,
    rank_65682 =
      rank_65682,
    rank_audit =
      gsea_rank_audit_step10,
    cache_manifest =
      msigdb_cache_manifest_step10
  ),
  file.path(
    GSEA_DIR,
    "checkpoint_step10C_GSEA_inputs.rds"
  )
)

cat(
  "\nSTEP 10C completed successfully.\n"
)

cat(
  "Hallmark memberships were read directly from the existing local RDS cache.\n"
)

cat(
  "No GSEA was run yet.\n"
)

cat(
  "No figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 10D: Triple-cohort Hallmark GSEA using
#           limma moderated-t ranked statistics
# ============================================================

cat(
  "\n========== STEP 10D: MODERATED-T HALLMARK GSEA ==========\n"
)

required_step10d_objects <- c(
  "hallmark_list",
  "rank_272769",
  "rank_95233",
  "rank_65682",
  "gsea_rank_audit_step10",
  "msigdb_version_step10",
  "msigdb_cache_manifest_step10",
  "GSEA_DIR"
)

missing_step10d_objects <- required_step10d_objects[
  !vapply(
    required_step10d_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_step10d_objects) > 0L) {
  stop(sprintf(
    "Missing objects required for Step 10D: %s",
    paste(
      missing_step10d_objects,
      collapse = ", "
    )
  ))
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("Package 'fgsea' is unavailable.")
}

stopifnot(
  length(hallmark_list) == 50L,
  length(rank_272769) == 18230L,
  length(rank_95233) == 18230L,
  length(rank_65682) == 11518L,
  !anyDuplicated(names(rank_272769)),
  !anyDuplicated(names(rank_95233)),
  !anyDuplicated(names(rank_65682)),
  all(is.finite(rank_272769)),
  all(is.finite(rank_95233)),
  all(is.finite(rank_65682))
)

# ------------------------------------------------------------
# 10D.1 Run one reproducible fgseaMultilevel analysis
# ------------------------------------------------------------

run_hallmark_fgsea_step10 <- function(
    ranked_statistic,
    dataset_name,
    random_seed
) {
  set.seed(random_seed)
  
  fgsea_result <- fgsea::fgseaMultilevel(
    pathways = hallmark_list,
    stats = ranked_statistic,
    minSize = 10L,
    maxSize = 500L,
    eps = 0,
    nPermSimple = 10000L,
    nproc = 1L
  )
  
  fgsea_result <- as.data.frame(
    fgsea_result
  )
  
  required_fgsea_columns <- c(
    "pathway",
    "pval",
    "padj",
    "ES",
    "NES",
    "size",
    "leadingEdge"
  )
  
  missing_fgsea_columns <- setdiff(
    required_fgsea_columns,
    colnames(fgsea_result)
  )
  
  if (length(missing_fgsea_columns) > 0L) {
    stop(sprintf(
      "%s fgsea output lacks columns: %s",
      dataset_name,
      paste(
        missing_fgsea_columns,
        collapse = ", "
      )
    ))
  }
  
  if (nrow(fgsea_result) != 50L) {
    stop(sprintf(
      "%s returned %d rather than 50 Hallmark pathways.",
      dataset_name,
      nrow(fgsea_result)
    ))
  }
  
  invalid_result <- !is.finite(
    fgsea_result$NES
  ) |
    !is.finite(
      fgsea_result$pval
    ) |
    !is.finite(
      fgsea_result$padj
    )
  
  if (any(invalid_result)) {
    stop(sprintf(
      "%s produced non-finite NES/P/FDR for %d pathways.",
      dataset_name,
      sum(invalid_result)
    ))
  }
  
  fgsea_result$leadingEdge <- vapply(
    fgsea_result$leadingEdge,
    function(current_genes) {
      paste(
        as.character(current_genes),
        collapse = ";"
      )
    },
    character(1)
  )
  
  fgsea_result$dataset <- dataset_name
  fgsea_result$rank_metric <-
    "limma moderated t"
  fgsea_result$MSigDB_version <-
    msigdb_version_step10
  fgsea_result$random_seed <-
    random_seed
  
  fgsea_result <- fgsea_result[
    order(
      fgsea_result$padj,
      fgsea_result$pval,
      -abs(fgsea_result$NES)
    ),
    ,
    drop = FALSE
  ]
  
  rownames(fgsea_result) <- NULL
  
  cat(sprintf(
    paste0(
      "%s: %d pathways | ",
      "nominal P<0.05=%d | ",
      "FDR<0.05=%d | ",
      "positive NES=%d | ",
      "negative NES=%d\n"
    ),
    dataset_name,
    nrow(fgsea_result),
    sum(fgsea_result$pval < 0.05),
    sum(fgsea_result$padj < 0.05),
    sum(fgsea_result$NES > 0),
    sum(fgsea_result$NES < 0)
  ))
  
  fgsea_result
}

# ------------------------------------------------------------
# 10D.2 Run all three cohorts
# ------------------------------------------------------------

fg_h_272769 <- run_hallmark_fgsea_step10(
  ranked_statistic = rank_272769,
  dataset_name = "GSE272769",
  random_seed = 4201L
)

fg_h_95233 <- run_hallmark_fgsea_step10(
  ranked_statistic = rank_95233,
  dataset_name = "GSE95233",
  random_seed = 4202L
)

fg_h_65682 <- run_hallmark_fgsea_step10(
  ranked_statistic = rank_65682,
  dataset_name = "GSE65682",
  random_seed = 4203L
)

fg_list <- list(
  GSE272769 = fg_h_272769,
  GSE95233 = fg_h_95233,
  GSE65682 = fg_h_65682
)

# ------------------------------------------------------------
# 10D.3 Per-cohort summary
# ------------------------------------------------------------

gsea_cohort_summary_step10 <- do.call(
  rbind,
  lapply(
    names(fg_list),
    function(current_dataset) {
      current_result <-
        fg_list[[current_dataset]]
      
      data.frame(
        dataset = current_dataset,
        pathways = nrow(current_result),
        nominal_P_lt_0_05 =
          sum(
            current_result$pval <
              0.05
          ),
        FDR_lt_0_05 =
          sum(
            current_result$padj <
              0.05
          ),
        positive_NES =
          sum(
            current_result$NES >
              0
          ),
        negative_NES =
          sum(
            current_result$NES <
              0
          ),
        median_absolute_NES =
          stats::median(
            abs(
              current_result$NES
            )
          ),
        maximum_absolute_NES =
          max(
            abs(
              current_result$NES
            )
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(gsea_cohort_summary_step10) <- NULL

# ------------------------------------------------------------
# 10D.4 Merge the three cohort results
# ------------------------------------------------------------

extract_gsea_summary_step10 <- function(
    fgsea_result,
    dataset_name
) {
  result_summary <- fgsea_result[
    ,
    c(
      "pathway",
      "ES",
      "NES",
      "pval",
      "padj",
      "size"
    ),
    drop = FALSE
  ]
  
  colnames(result_summary)[
    -1L
  ] <- paste0(
    colnames(result_summary)[-1L],
    "_",
    dataset_name
  )
  
  result_summary
}

hall_3 <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "pathway",
      all = TRUE
    )
  },
  list(
    extract_gsea_summary_step10(
      fg_h_272769,
      "GSE272769"
    ),
    extract_gsea_summary_step10(
      fg_h_95233,
      "GSE95233"
    ),
    extract_gsea_summary_step10(
      fg_h_65682,
      "GSE65682"
    )
  )
)

hall_3 <- as.data.frame(
  hall_3
)

expected_nes_columns_step10 <- c(
  "NES_GSE272769",
  "NES_GSE95233",
  "NES_GSE65682"
)

expected_pval_columns_step10 <- c(
  "pval_GSE272769",
  "pval_GSE95233",
  "pval_GSE65682"
)

expected_padj_columns_step10 <- c(
  "padj_GSE272769",
  "padj_GSE95233",
  "padj_GSE65682"
)

stopifnot(
  nrow(hall_3) == 50L,
  all(
    expected_nes_columns_step10 %in%
      colnames(hall_3)
  ),
  all(
    expected_pval_columns_step10 %in%
      colnames(hall_3)
  ),
  all(
    expected_padj_columns_step10 %in%
      colnames(hall_3)
  )
)

nes_matrix_step10 <- as.matrix(
  hall_3[
    ,
    expected_nes_columns_step10,
    drop = FALSE
  ]
)

pval_matrix_step10 <- as.matrix(
  hall_3[
    ,
    expected_pval_columns_step10,
    drop = FALSE
  ]
)

padj_matrix_step10 <- as.matrix(
  hall_3[
    ,
    expected_padj_columns_step10,
    drop = FALSE
  ]
)

if (
  any(!is.finite(nes_matrix_step10)) ||
  any(!is.finite(pval_matrix_step10)) ||
  any(!is.finite(padj_matrix_step10))
) {
  stop(
    "Triple-cohort GSEA summary contains non-finite values."
  )
}

hall_3$mean_abs_NES <- rowMeans(
  abs(nes_matrix_step10)
)

hall_3$all_same_dir <- apply(
  nes_matrix_step10,
  1L,
  function(current_nes) {
    all(current_nes != 0) &&
      length(
        unique(
          sign(current_nes)
        )
      ) == 1L
  }
)

hall_3$n_same_pairs <- apply(
  nes_matrix_step10,
  1L,
  function(current_nes) {
    direction_vector <- sign(
      current_nes
    )
    
    agreement_matrix <- outer(
      direction_vector,
      direction_vector,
      "=="
    )
    
    sum(
      agreement_matrix[
        upper.tri(
          agreement_matrix
        )
      ]
    )
  }
)

hall_3$total_pairs <- 3L

hall_3$concordant_direction <- ifelse(
  hall_3$all_same_dir &
    hall_3$NES_GSE272769 > 0,
  "Positive",
  ifelse(
    hall_3$all_same_dir &
      hall_3$NES_GSE272769 < 0,
    "Negative",
    "Mixed"
  )
)

hall_3$n_cohorts_nominal_P_lt_0_05 <-
  rowSums(
    pval_matrix_step10 < 0.05
  )

hall_3$n_cohorts_FDR_lt_0_05 <-
  rowSums(
    padj_matrix_step10 < 0.05
  )

hall_3$all_three_nominal_P_lt_0_05 <-
  hall_3$
  n_cohorts_nominal_P_lt_0_05 ==
  3L

hall_3$all_three_FDR_lt_0_05 <-
  hall_3$
  n_cohorts_FDR_lt_0_05 ==
  3L

hall_3$display_rank_mean_abs_NES <- rank(
  -hall_3$mean_abs_NES,
  ties.method = "first"
)

hall_3 <- hall_3[
  order(hall_3$pathway),
  ,
  drop = FALSE
]

rownames(hall_3) <- NULL

# ------------------------------------------------------------
# 10D.5 Pairwise NES concordance
# ------------------------------------------------------------

gsea_pairs_step10 <- list(
  GSE272769_vs_GSE95233 =
    c(
      "NES_GSE272769",
      "NES_GSE95233"
    ),
  GSE272769_vs_GSE65682 =
    c(
      "NES_GSE272769",
      "NES_GSE65682"
    ),
  GSE95233_vs_GSE65682 =
    c(
      "NES_GSE95233",
      "NES_GSE65682"
    )
)

gsea_pairwise_concordance_step10 <- do.call(
  rbind,
  lapply(
    names(gsea_pairs_step10),
    function(pair_name) {
      current_columns <-
        gsea_pairs_step10[[pair_name]]
      
      x <- as.numeric(
        hall_3[[current_columns[1L]]]
      )
      
      y <- as.numeric(
        hall_3[[current_columns[2L]]]
      )
      
      data.frame(
        comparison = pair_name,
        pathways = length(x),
        same_direction =
          sum(
            sign(x) ==
              sign(y)
          ),
        same_direction_fraction =
          mean(
            sign(x) ==
              sign(y)
          ),
        Pearson_r =
          stats::cor(
            x,
            y,
            method = "pearson"
          ),
        Spearman_rho =
          stats::cor(
            x,
            y,
            method = "spearman"
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  gsea_pairwise_concordance_step10
) <- NULL

# ------------------------------------------------------------
# 10D.6 Run manifest
# ------------------------------------------------------------

gsea_run_manifest_step10 <- data.frame(
  dataset = c(
    "GSE272769",
    "GSE95233",
    "GSE65682"
  ),
  rank_metric =
    "limma moderated t",
  cohort_adjustment = c(
    "Not applicable",
    "Not applicable",
    "Included"
  ),
  ranked_genes = c(
    length(rank_272769),
    length(rank_95233),
    length(rank_65682)
  ),
  Hallmark_pathways = 50L,
  minSize = 10L,
  maxSize = 500L,
  eps = 0,
  nPermSimple = 10000L,
  nproc = 1L,
  random_seed = c(
    4201L,
    4202L,
    4203L
  ),
  fgsea_version =
    as.character(
      utils::packageVersion("fgsea")
    ),
  MSigDB_version =
    msigdb_version_step10,
  MSigDB_cache_MD5 =
    msigdb_cache_manifest_step10$
    cache_MD5,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10D.7 Print results
# ------------------------------------------------------------

cat("\nPer-cohort Hallmark GSEA summary:\n")

print(
  gsea_cohort_summary_step10,
  row.names = FALSE,
  digits = 6
)

cat(sprintf(
  paste0(
    "\nTriple-cohort Hallmark GSEA: ",
    "%d pathways | ",
    "all-three directionally concordant=%d (%.1f%%) | ",
    "all-three nominal P<0.05=%d | ",
    "all-three FDR<0.05=%d\n"
  ),
  nrow(hall_3),
  sum(hall_3$all_same_dir),
  100 * mean(hall_3$all_same_dir),
  sum(
    hall_3$
      all_three_nominal_P_lt_0_05
  ),
  sum(
    hall_3$
      all_three_FDR_lt_0_05
  )
))

cat("\nPairwise NES concordance:\n")

print(
  gsea_pairwise_concordance_step10,
  row.names = FALSE,
  digits = 6
)

hallmark_top15_step10 <- hall_3[
  order(
    -hall_3$mean_abs_NES
  ),
  c(
    "pathway",
    expected_nes_columns_step10,
    "mean_abs_NES",
    "all_same_dir",
    "concordant_direction",
    "n_cohorts_nominal_P_lt_0_05",
    "n_cohorts_FDR_lt_0_05"
  ),
  drop = FALSE
]

hallmark_top15_step10 <-
  head(
    hallmark_top15_step10,
    15L
  )

cat("\nTop 15 pathways by mean absolute NES:\n")

print(
  hallmark_top15_step10,
  row.names = FALSE,
  digits = 5
)

# ------------------------------------------------------------
# 10D.8 Export
# ------------------------------------------------------------

write.csv(
  fg_h_272769,
  file.path(
    GSEA_DIR,
    "GSEA_Hallmark_GSE272769_moderated_t.csv"
  ),
  row.names = FALSE
)

write.csv(
  fg_h_95233,
  file.path(
    GSEA_DIR,
    "GSEA_Hallmark_GSE95233_D01_moderated_t.csv"
  ),
  row.names = FALSE
)

write.csv(
  fg_h_65682,
  file.path(
    GSEA_DIR,
    paste0(
      "GSEA_Hallmark_GSE65682_",
      "cohort_adjusted_moderated_t.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  hall_3,
  file.path(
    GSEA_DIR,
    paste0(
      "GSEA_Hallmark_triple_cohort_",
      "moderated_t.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  gsea_cohort_summary_step10,
  file.path(
    GSEA_DIR,
    "GSEA_Hallmark_cohort_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_pairwise_concordance_step10,
  file.path(
    GSEA_DIR,
    "GSEA_Hallmark_pairwise_NES_concordance.csv"
  ),
  row.names = FALSE
)

write.csv(
  hallmark_top15_step10,
  file.path(
    GSEA_DIR,
    "GSEA_Hallmark_top15_display_source.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_run_manifest_step10,
  file.path(
    GSEA_DIR,
    "GSEA_Hallmark_run_manifest.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    fg_h_272769 =
      fg_h_272769,
    fg_h_95233 =
      fg_h_95233,
    fg_h_65682 =
      fg_h_65682,
    fg_list =
      fg_list,
    hall_3 =
      hall_3,
    cohort_summary =
      gsea_cohort_summary_step10,
    pairwise_concordance =
      gsea_pairwise_concordance_step10,
    top15 =
      hallmark_top15_step10,
    run_manifest =
      gsea_run_manifest_step10
  ),
  file.path(
    GSEA_DIR,
    "checkpoint_step10D_Hallmark_GSEA_results.rds"
  )
)

cat(
  "\nSTEP 10D completed successfully.\n"
)

cat(
  "All three rankings used limma moderated t statistics.\n"
)

cat(
  "GSE65682 used the cohort-adjusted moderated t statistic.\n"
)

cat(
  "No figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 10E: Compare the previous GSEA summary with the
#           moderated-t revision
# ============================================================

cat(
  "\n========== STEP 10E: OLD-VERSUS-REVISED GSEA COMPARISON ==========\n"
)

required_step10e_objects <- c(
  "hall_3",
  "expected_nes_columns_step10",
  "GSEA_DIR"
)

missing_step10e_objects <- required_step10e_objects[
  !vapply(
    required_step10e_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_step10e_objects) > 0L) {
  stop(sprintf(
    "Missing objects required for Step 10E: %s",
    paste(
      missing_step10e_objects,
      collapse = ", "
    )
  ))
}

project_dir_step10e <- if (
  exists("PROJECT_DIR")
) {
  PROJECT_DIR
} else {
  "/home/sunshine/predicate/test03"
}

old_gsea_file_step10e <- file.path(
  project_dir_step10e,
  "results_meta",
  "05_triple_cohort",
  "GSEA_Hallmark_triple_cohort.csv"
)

if (!file.exists(old_gsea_file_step10e)) {
  stop(sprintf(
    "Previous GSEA summary was not found: %s",
    old_gsea_file_step10e
  ))
}

old_hall_3_step10e <- read.csv(
  old_gsea_file_step10e,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_comparison_columns_step10e <- c(
  "pathway",
  expected_nes_columns_step10,
  "mean_abs_NES",
  "all_same_dir"
)

missing_old_columns_step10e <- setdiff(
  required_comparison_columns_step10e,
  colnames(old_hall_3_step10e)
)

missing_new_columns_step10e <- setdiff(
  required_comparison_columns_step10e,
  colnames(hall_3)
)

if (
  length(missing_old_columns_step10e) >
  0L ||
  length(missing_new_columns_step10e) >
  0L
) {
  stop(sprintf(
    paste0(
      "Missing comparison columns. ",
      "Old: %s; New: %s"
    ),
    paste(
      missing_old_columns_step10e,
      collapse = ", "
    ),
    paste(
      missing_new_columns_step10e,
      collapse = ", "
    )
  ))
}

stopifnot(
  nrow(old_hall_3_step10e) == 50L,
  nrow(hall_3) == 50L,
  !anyDuplicated(
    old_hall_3_step10e$pathway
  ),
  !anyDuplicated(
    hall_3$pathway
  )
)

old_comparison_step10e <-
  old_hall_3_step10e[
    ,
    required_comparison_columns_step10e,
    drop = FALSE
  ]

new_comparison_step10e <-
  hall_3[
    ,
    required_comparison_columns_step10e,
    drop = FALSE
  ]

colnames(old_comparison_step10e) <- c(
  "pathway",
  paste0(
    expected_nes_columns_step10,
    "_old"
  ),
  "mean_abs_NES_old",
  "all_same_dir_old"
)

colnames(new_comparison_step10e) <- c(
  "pathway",
  paste0(
    expected_nes_columns_step10,
    "_new"
  ),
  "mean_abs_NES_new",
  "all_same_dir_new"
)

gsea_old_new_comparison_step10e <- merge(
  old_comparison_step10e,
  new_comparison_step10e,
  by = "pathway",
  all = TRUE
)

if (
  nrow(gsea_old_new_comparison_step10e) !=
  50L
) {
  stop(
    "Old/new pathway reconciliation did not return 50 pathways."
  )
}

dataset_names_step10e <- c(
  "GSE272769",
  "GSE95233",
  "GSE65682"
)

for (
  current_dataset in
  dataset_names_step10e
) {
  old_column <- paste0(
    "NES_",
    current_dataset,
    "_old"
  )
  
  new_column <- paste0(
    "NES_",
    current_dataset,
    "_new"
  )
  
  delta_column <- paste0(
    "NES_delta_",
    current_dataset
  )
  
  direction_column <- paste0(
    "direction_changed_",
    current_dataset
  )
  
  gsea_old_new_comparison_step10e[[delta_column]] <-
    gsea_old_new_comparison_step10e[[new_column]] -
    gsea_old_new_comparison_step10e[[old_column]]
  
  gsea_old_new_comparison_step10e[[direction_column]] <-
    sign(
      gsea_old_new_comparison_step10e[[new_column]]
    ) !=
    sign(
      gsea_old_new_comparison_step10e[[old_column]]
    )
}

gsea_old_new_summary_step10e <- do.call(
  rbind,
  lapply(
    dataset_names_step10e,
    function(current_dataset) {
      old_column <- paste0(
        "NES_",
        current_dataset,
        "_old"
      )
      
      new_column <- paste0(
        "NES_",
        current_dataset,
        "_new"
      )
      
      old_nes <-
        gsea_old_new_comparison_step10e[[old_column]]
      
      new_nes <-
        gsea_old_new_comparison_step10e[[new_column]]
      
      data.frame(
        dataset = current_dataset,
        pathways = length(old_nes),
        Pearson_r =
          stats::cor(
            old_nes,
            new_nes,
            method = "pearson"
          ),
        Spearman_rho =
          stats::cor(
            old_nes,
            new_nes,
            method = "spearman"
          ),
        mean_absolute_NES_change =
          mean(
            abs(
              new_nes -
                old_nes
            )
          ),
        maximum_absolute_NES_change =
          max(
            abs(
              new_nes -
                old_nes
            )
          ),
        direction_changed =
          sum(
            sign(old_nes) !=
              sign(new_nes)
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(gsea_old_new_summary_step10e) <- NULL

logical_from_csv_step10e <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  toupper(
    trimws(
      as.character(x)
    )
  ) == "TRUE"
}

old_same_direction_step10e <-
  logical_from_csv_step10e(
    gsea_old_new_comparison_step10e$
      all_same_dir_old
  )

new_same_direction_step10e <-
  logical_from_csv_step10e(
    gsea_old_new_comparison_step10e$
      all_same_dir_new
  )

old_top15_step10e <- head(
  old_hall_3_step10e$pathway[
    order(
      -old_hall_3_step10e$
        mean_abs_NES
    )
  ],
  15L
)

new_top15_step10e <- head(
  hall_3$pathway[
    order(
      -hall_3$mean_abs_NES
    )
  ],
  15L
)

gsea_old_new_comparison_step10e$
  displayed_in_old_top15 <-
  gsea_old_new_comparison_step10e$
  pathway %in%
  old_top15_step10e

gsea_old_new_comparison_step10e$
  displayed_in_new_top15 <-
  gsea_old_new_comparison_step10e$
  pathway %in%
  new_top15_step10e

gsea_global_revision_summary_step10e <-
  data.frame(
    old_all_three_directionally_concordant =
      sum(
        old_same_direction_step10e,
        na.rm = TRUE
      ),
    new_all_three_directionally_concordant =
      sum(
        new_same_direction_step10e,
        na.rm = TRUE
      ),
    concordant_count_change =
      sum(
        new_same_direction_step10e,
        na.rm = TRUE
      ) -
      sum(
        old_same_direction_step10e,
        na.rm = TRUE
      ),
    old_top15_pathways = 15L,
    new_top15_pathways = 15L,
    top15_overlap =
      length(
        intersect(
          old_top15_step10e,
          new_top15_step10e
        )
      ),
    pathways_entering_top15 =
      length(
        setdiff(
          new_top15_step10e,
          old_top15_step10e
        )
      ),
    pathways_leaving_top15 =
      length(
        setdiff(
          old_top15_step10e,
          new_top15_step10e
        )
      ),
    stringsAsFactors = FALSE
  )

cat("\nOld-versus-revised NES comparison:\n")

print(
  gsea_old_new_summary_step10e,
  row.names = FALSE,
  digits = 6
)

cat("\nGlobal display comparison:\n")

print(
  gsea_global_revision_summary_step10e,
  row.names = FALSE
)

cat("\nPathways entering the revised Top 15:\n")

print(
  setdiff(
    new_top15_step10e,
    old_top15_step10e
  )
)

cat("\nPathways leaving the revised Top 15:\n")

print(
  setdiff(
    old_top15_step10e,
    new_top15_step10e
  )
)

write.csv(
  gsea_old_new_comparison_step10e,
  file.path(
    GSEA_DIR,
    "GSEA_old_vs_moderated_t_pathway_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_old_new_summary_step10e,
  file.path(
    GSEA_DIR,
    "GSEA_old_vs_moderated_t_cohort_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_global_revision_summary_step10e,
  file.path(
    GSEA_DIR,
    "GSEA_old_vs_moderated_t_global_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\nSTEP 10E completed successfully.\n"
)

cat(
  "No GSEA was rerun in this comparison step.\n"
)

cat(
  "No figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)
# ============================================================
# Step 11A: Revision-script completeness and syntax audit
# ============================================================

cat(
  "\n========== STEP 11A: REVISION SCRIPT AUDIT ==========\n"
)

SCRIPT_FILE_STEP11 <- paste0(
  "/home/sunshine/predicate/test03/",
  "05bulk_rerun_priority1_bulk.R"
)

if (!file.exists(SCRIPT_FILE_STEP11)) {
  stop(sprintf(
    "Revision script was not found: %s",
    SCRIPT_FILE_STEP11
  ))
}

script_lines_step11 <- readLines(
  SCRIPT_FILE_STEP11,
  warn = FALSE,
  encoding = "UTF-8"
)

required_script_markers_step11 <- c(
  "Step 8: Cohort-adjusted MARS-by-module interaction",
  "Step 9: Cohort-adjusted PC1 outcome association",
  "Step 10C: Load the local MSigDB cache",
  "Step 10D: Triple-cohort Hallmark GSEA",
  "Step 10E: Compare the previous GSEA summary"
)

script_marker_audit_step11 <- data.frame(
  marker = required_script_markers_step11,
  present = vapply(
    required_script_markers_step11,
    function(current_marker) {
      any(
        grepl(
          current_marker,
          script_lines_step11,
          fixed = TRUE
        )
      )
    },
    logical(1)
  ),
  stringsAsFactors = FALSE
)

cat("\nRequired section markers:\n")

print(
  script_marker_audit_step11,
  row.names = FALSE
)

missing_script_markers_step11 <-
  script_marker_audit_step11$marker[
    !script_marker_audit_step11$present
  ]

if (
  length(missing_script_markers_step11) >
  0L
) {
  stop(sprintf(
    paste0(
      "The following completed analysis blocks ",
      "have not been saved in the script: %s"
    ),
    paste(
      missing_script_markers_step11,
      collapse = "; "
    )
  ))
}

# Detect the earlier broken [[ syntax pattern.
broken_double_bracket_lines_step11 <- grep(
  paste0(
    "cohort_roc_objects_step9\\s*\\[$|",
    "descriptive_auc_rows_step9\\s*\\[$"
  ),
  trimws(script_lines_step11),
  value = TRUE
)

cat("\nBroken double-bracket patterns found:\n")

print(
  broken_double_bracket_lines_step11
)

if (
  length(broken_double_bracket_lines_step11) >
  0L
) {
  stop(
    "The saved script still contains the broken Step 9 bracket syntax."
  )
}

parse_error_step11 <- tryCatch(
  {
    parse(
      file = SCRIPT_FILE_STEP11,
      encoding = "UTF-8"
    )
    
    NA_character_
  },
  error = function(e) {
    conditionMessage(e)
  }
)

cat("\nFull-script syntax result:\n")

if (is.na(parse_error_step11)) {
  cat("PASS: the complete script parsed successfully.\n")
} else {
  cat(
    paste0(
      "FAIL: ",
      parse_error_step11,
      "\n"
    )
  )
  
  stop(
    "The complete revision script contains a syntax error."
  )
}

write.csv(
  script_marker_audit_step11,
  file.path(
    REVISION_DIR,
    "revision_script_section_audit.csv"
  ),
  row.names = FALSE
)

cat(
  "\nSTEP 11A completed successfully.\n"
)

cat(
  "All required revision blocks are saved in the script.\n"
)

cat(
  "The full revision script is syntactically valid.\n"
)

cat(
  "No analysis was rerun.\n"
)

cat(
  "No figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)







# ============================================================
# Illustrator-compatible PDF device
# ============================================================

ai_pdf_device_step12 <- function(
    filename,
    width,
    height,
    bg = "white",
    ...
) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    family = "Helvetica",
    useDingbats = FALSE,
    bg = bg,
    ...
  )
}
# ============================================================
# Step 12A: Revised Figure 4B
# External assessment of effect directions and magnitudes
# ============================================================

cat(
  "\n========== STEP 12A: REVISED FIGURE 4B ==========\n"
)

suppressPackageStartupMessages(
  library(ggplot2)
)

suppressPackageStartupMessages(
  library(ggrepel)
)

# ------------------------------------------------------------
# 12A.1 Required-object audit
# ------------------------------------------------------------

required_objects_step12a <- c(
  "REVISION_DIR",
  "fisher_out",
  "res_val_adjusted"
)

missing_objects_step12a <- required_objects_step12a[
  !vapply(
    required_objects_step12a,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects_step12a) > 0L) {
  stop(
    sprintf(
      "Missing required objects: %s",
      paste(
        missing_objects_step12a,
        collapse = ", "
      )
    )
  )
}

required_discovery_columns_step12a <- c(
  "gene",
  "re_estimate",
  "direction"
)

required_external_columns_step12a <- c(
  "gene",
  "logFC",
  "t",
  "P.Value",
  "adj.P.Val"
)

stopifnot(
  all(
    required_discovery_columns_step12a %in%
      colnames(fisher_out)
  ),
  all(
    required_external_columns_step12a %in%
      colnames(res_val_adjusted)
  )
)

# ------------------------------------------------------------
# 12A.2 Output directory
# ------------------------------------------------------------

figure_root_step12a <- file.path(
  REVISION_DIR,
  "12_revised_figures"
)

figure_dir_step12a <- file.path(
  figure_root_step12a,
  "Figure_4B_external_assessment"
)

dir.create(
  figure_dir_step12a,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 12A.3 Construct plotting data explicitly from the
# adjusted GSE65682 analysis
# ------------------------------------------------------------

signature_genes_step12a <- unique(
  as.character(fisher_out$gene)
)

evaluable_genes_step12a <-
  signature_genes_step12a[
    signature_genes_step12a %in%
      res_val_adjusted$gene
  ]

missing_genes_step12a <- setdiff(
  signature_genes_step12a,
  evaluable_genes_step12a
)

discovery_index_step12a <- match(
  evaluable_genes_step12a,
  fisher_out$gene
)

external_index_step12a <- match(
  evaluable_genes_step12a,
  res_val_adjusted$gene
)

figure4b_data_step12a <- data.frame(
  gene = evaluable_genes_step12a,
  
  discovery_logFC =
    fisher_out$re_estimate[
      discovery_index_step12a
    ],
  
  external_adjusted_logFC =
    res_val_adjusted$logFC[
      external_index_step12a
    ],
  
  external_moderated_t =
    res_val_adjusted$t[
      external_index_step12a
    ],
  
  external_P =
    res_val_adjusted$P.Value[
      external_index_step12a
    ],
  
  external_FDR =
    res_val_adjusted$adj.P.Val[
      external_index_step12a
    ],
  
  discovery_direction =
    fisher_out$direction[
      discovery_index_step12a
    ],
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 12A.4 Reproducibility assertions
# ------------------------------------------------------------

stopifnot(
  length(signature_genes_step12a) == 32L,
  nrow(figure4b_data_step12a) == 30L,
  setequal(
    missing_genes_step12a,
    c(
      "HIST1H2BM",
      "HIST1H3B"
    )
  ),
  !anyNA(
    figure4b_data_step12a$discovery_logFC
  ),
  !anyNA(
    figure4b_data_step12a$
      external_adjusted_logFC
  ),
  all(
    is.finite(
      figure4b_data_step12a$
        discovery_logFC
    )
  ),
  all(
    is.finite(
      figure4b_data_step12a$
        external_adjusted_logFC
    )
  )
)

# ------------------------------------------------------------
# 12A.5 Direction and external magnitude classifications
# ------------------------------------------------------------

figure4b_data_step12a$
  direction_retained <-
  sign(
    figure4b_data_step12a$
      discovery_logFC
  ) ==
  sign(
    figure4b_data_step12a$
      external_adjusted_logFC
  )

external_logFC_cut_step12a <- 0.5

figure4b_data_step12a$
  external_above_threshold <-
  abs(
    figure4b_data_step12a$
      external_adjusted_logFC
  ) >
  external_logFC_cut_step12a

n_signature_step12a <-
  length(signature_genes_step12a)

n_evaluable_step12a <-
  nrow(figure4b_data_step12a)

n_retained_step12a <- sum(
  figure4b_data_step12a$
    direction_retained
)

n_above_threshold_step12a <- sum(
  figure4b_data_step12a$
    external_above_threshold
)

n_nominal_step12a <- sum(
  figure4b_data_step12a$external_P < 0.05,
  na.rm = TRUE
)

n_fdr_step12a <- sum(
  figure4b_data_step12a$external_FDR < 0.05,
  na.rm = TRUE
)

stopifnot(
  n_retained_step12a == 30L,
  n_above_threshold_step12a == 10L,
  n_nominal_step12a == 21L,
  n_fdr_step12a == 14L
)

high_class_step12a <- sprintf(
  "|External log2FC| > 0.5 (n = %d)",
  n_above_threshold_step12a
)

lower_class_step12a <- sprintf(
  "|External log2FC| <= 0.5 (n = %d)",
  n_evaluable_step12a -
    n_above_threshold_step12a
)

figure4b_data_step12a$
  external_effect_class <- factor(
    ifelse(
      figure4b_data_step12a$
        external_above_threshold,
      high_class_step12a,
      lower_class_step12a
    ),
    levels = c(
      high_class_step12a,
      lower_class_step12a
    )
  )

# Label the 10 genes exceeding the display threshold.
figure4b_data_step12a$label <- ifelse(
  figure4b_data_step12a$
    external_above_threshold,
  figure4b_data_step12a$gene,
  ""
)

# ------------------------------------------------------------
# 12A.6 Descriptive cross-gene statistics
# ------------------------------------------------------------

pearson_r_step12a <- cor(
  figure4b_data_step12a$
    discovery_logFC,
  figure4b_data_step12a$
    external_adjusted_logFC,
  method = "pearson"
)

spearman_rho_step12a <- cor(
  figure4b_data_step12a$
    discovery_logFC,
  figure4b_data_step12a$
    external_adjusted_logFC,
  method = "spearman"
)

ols_fit_step12a <- lm(
  external_adjusted_logFC ~
    discovery_logFC,
  data = figure4b_data_step12a
)

ols_intercept_step12a <- unname(
  coef(ols_fit_step12a)[1]
)

ols_slope_step12a <- unname(
  coef(ols_fit_step12a)[2]
)

ols_ci_step12a <- confint(
  ols_fit_step12a,
  level = 0.95
)

ols_slope_ci_low_step12a <-
  unname(
    ols_ci_step12a[
      "discovery_logFC",
      1
    ]
  )

ols_slope_ci_high_step12a <-
  unname(
    ols_ci_step12a[
      "discovery_logFC",
      2
    ]
  )

# Expected values from the completed adjusted analysis.
stopifnot(
  abs(
    pearson_r_step12a - 0.813635
  ) < 0.0001,
  abs(
    ols_slope_step12a - 0.696881
  ) < 0.0001
)

# ------------------------------------------------------------
# 12A.7 Export source data and metrics
# ------------------------------------------------------------

figure4b_metrics_step12a <- data.frame(
  Metric = c(
    "Signature genes",
    "Evaluable genes",
    "Direction retained",
    "External absolute log2FC above 0.5",
    "External nominal P below 0.05",
    "External FDR below 0.05",
    "Descriptive cross-gene Pearson r",
    "Descriptive cross-gene Spearman rho",
    "Descriptive OLS intercept",
    "Descriptive OLS slope",
    "OLS slope lower 95% CI",
    "OLS slope upper 95% CI"
  ),
  
  Value = c(
    n_signature_step12a,
    n_evaluable_step12a,
    n_retained_step12a,
    n_above_threshold_step12a,
    n_nominal_step12a,
    n_fdr_step12a,
    pearson_r_step12a,
    spearman_rho_step12a,
    ols_intercept_step12a,
    ols_slope_step12a,
    ols_slope_ci_low_step12a,
    ols_slope_ci_high_step12a
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  figure4b_data_step12a,
  file.path(
    figure_dir_step12a,
    "Figure_4B_source_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  figure4b_metrics_step12a,
  file.path(
    figure_dir_step12a,
    "Figure_4B_metrics.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12A.8 Common axis range
# ------------------------------------------------------------

common_limits_step12a <- range(
  c(
    figure4b_data_step12a$
      discovery_logFC,
    figure4b_data_step12a$
      external_adjusted_logFC,
    -external_logFC_cut_step12a,
    external_logFC_cut_step12a,
    0
  ),
  finite = TRUE
)

limit_padding_step12a <- max(
  diff(common_limits_step12a) * 0.08,
  0.05
)

common_limits_step12a <- c(
  common_limits_step12a[1] -
    limit_padding_step12a,
  common_limits_step12a[2] +
    limit_padding_step12a
)

effect_colors_step12a <- setNames(
  c(
    "#D94841",
    "#3977B7"
  ),
  c(
    high_class_step12a,
    lower_class_step12a
  )
)

# ------------------------------------------------------------
# 12A.9 Draw final revised Figure 4B
# ------------------------------------------------------------

# Shorter legend labels.
figure4b_data_step12a$
  external_effect_class <- factor(
    ifelse(
      figure4b_data_step12a$
        external_above_threshold,
      "Above threshold",
      "At or below threshold"
    ),
    levels = c(
      "Above threshold",
      "At or below threshold"
    )
  )

effect_colors_step12a <- c(
  "Above threshold" = "#D94841",
  "At or below threshold" = "#3977B7"
)

# Fixed, reproducible label layout for the 10 highlighted genes.
label_layout_step12a <- data.frame(
  gene = c(
    "CA2",
    "RHAG",
    "CTSG",
    "CD24",
    "ELANE",
    "MS4A3",
    "CEACAM6",
    "CEACAM8",
    "DEFA4",
    "CX3CR1"
  ),
  
  label_dx = c(
    -0.09,
    -0.09,
    -0.07,
    0.06,
    -0.05,
    -0.07,
    0.12,
    0.11,
    0.08,
    -0.07
  ),
  
  label_dy = c(
    0.09,
    -0.10,
    0.11,
    -0.07,
    -0.14,
    0.13,
    0.10,
    -0.15,
    0.12,
    -0.10
  ),
  
  label_hjust = c(
    1,
    1,
    1,
    0,
    1,
    1,
    0,
    0,
    0,
    1
  ),
  
  stringsAsFactors = FALSE
)

label_data_step12a <- merge(
  figure4b_data_step12a,
  label_layout_step12a,
  by = "gene",
  all = FALSE,
  sort = FALSE
)

label_data_step12a <- label_data_step12a[
  match(
    label_layout_step12a$gene,
    label_data_step12a$gene
  ),
  ,
  drop = FALSE
]

stopifnot(
  nrow(label_data_step12a) == 10L,
  identical(
    label_data_step12a$gene,
    label_layout_step12a$gene
  ),
  all(
    label_data_step12a$
      external_above_threshold
  )
)

label_data_step12a$label_x <-
  label_data_step12a$discovery_logFC +
  label_data_step12a$label_dx

label_data_step12a$label_y <-
  label_data_step12a$
  external_adjusted_logFC +
  label_data_step12a$label_dy

p_figure4b_step12a <- ggplot(
  figure4b_data_step12a,
  aes(
    x = discovery_logFC,
    y = external_adjusted_logFC
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    color = "grey58"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey58"
  ) +
  geom_hline(
    yintercept = c(
      -external_logFC_cut_step12a,
      external_logFC_cut_step12a
    ),
    linetype = "dotted",
    linewidth = 0.4,
    color = "grey72"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.55,
    color = "grey58"
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey15",
    fill = "grey55",
    linewidth = 0.7,
    linetype = "dashed",
    alpha = 0.18
  ) +
  geom_segment(
    data = label_data_step12a,
    aes(
      x = discovery_logFC,
      y = external_adjusted_logFC,
      xend = label_x,
      yend = label_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.28,
    color = "grey48"
  ) +
  geom_point(
    aes(
      fill = external_effect_class
    ),
    shape = 21,
    size = 4.0,
    stroke = 0.35,
    color = "white",
    alpha = 0.94
  ) +
  geom_text(
    data = label_data_step12a,
    aes(
      x = label_x,
      y = label_y,
      label = gene,
      hjust = label_hjust
    ),
    inherit.aes = FALSE,
    size = 3.05,
    vjust = 0.5,
    color = "grey12",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = effect_colors_step12a,
    breaks = c(
      "Above threshold",
      "At or below threshold"
    ),
    labels = c(
      sprintf(
        "> 0.5 (n = %d)",
        n_above_threshold_step12a
      ),
      sprintf(
        "\u2264 0.5 (n = %d)",
        n_evaluable_step12a -
          n_above_threshold_step12a
      )
    ),
    drop = FALSE
  ) +
  coord_equal(
    xlim = common_limits_step12a,
    ylim = common_limits_step12a,
    clip = "off"
  ) +
  labs(
    title = paste(
      "External assessment of effect",
      "directions and magnitudes"
    ),
    
    subtitle = sprintf(
      paste0(
        "%d/%d evaluable genes retained direction\n",
        "Descriptive cross-gene Pearson r = %.3f; ",
        "external-on-discovery slope = %.3f"
      ),
      n_retained_step12a,
      n_evaluable_step12a,
      pearson_r_step12a,
      ols_slope_step12a
    ),
    
    x = paste0(
      "Discovery random-effects log\u2082FC\n",
      "(non-survivors vs survivors)"
    ),
    
    y = paste0(
      "GSE65682 cohort-adjusted log\u2082FC\n",
      "(non-survivors vs survivors)"
    ),
    
    fill = paste0(
      "Cohort-adjusted external |log\u2082FC|"
    )
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid = element_blank(),
    
    plot.title = element_text(
      size = 13,
      face = "bold",
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      size = 9.5,
      lineheight = 1.18,
      margin = margin(b = 8)
    ),
    
    axis.title = element_text(
      size = 10.2
    ),
    
    axis.text = element_text(
      size = 9.2,
      color = "grey15"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 8.8,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8.7
    ),
    
    legend.key.width = grid::unit(
      0.55,
      "cm"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 20,
      b = 8,
      l = 8
    )
  )

print(
  p_figure4b_step12a
)

# ------------------------------------------------------------
# 12A.10 Save final standalone panel files
# ------------------------------------------------------------

pdf_device_step12a <- ai_pdf_device_step12

ggsave(
  filename = file.path(
    figure_dir_step12a,
    paste0(
      "Figure_4B_external_effect_",
      "assessment_cohort_adjusted.pdf"
    )
  ),
  plot = p_figure4b_step12a,
  device = pdf_device_step12a,
  width = 7.8,
  height = 6.5,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir_step12a,
    paste0(
      "Figure_4B_external_effect_",
      "assessment_cohort_adjusted.png"
    )
  ),
  plot = p_figure4b_step12a,
  width = 7.8,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir_step12a,
    paste0(
      "Figure_4B_external_effect_",
      "assessment_cohort_adjusted.tiff"
    )
  ),
  plot = p_figure4b_step12a,
  device = "tiff",
  width = 7.8,
  height = 6.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

saveRDS(
  p_figure4b_step12a,
  file.path(
    figure_dir_step12a,
    "Figure_4B_plot_object.rds"
  )
)



# ============================================================
# Step 12B: Triple-cohort Hallmark GSEA heatmap
# Original pheatmap design using revised moderated-t results
# ============================================================

cat(
  "\n========== STEP 12B: HALLMARK GSEA HEATMAP ==========\n"
)

suppressPackageStartupMessages(
  library(pheatmap)
)

# ------------------------------------------------------------
# 12B.1 Required-object and schema audit
# ------------------------------------------------------------

required_objects_step12b <- c(
  "REVISION_DIR",
  "hall_3"
)

missing_objects_step12b <- required_objects_step12b[
  !vapply(
    required_objects_step12b,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects_step12b) > 0L) {
  stop(
    sprintf(
      "Missing required objects: %s",
      paste(
        missing_objects_step12b,
        collapse = ", "
      )
    )
  )
}

required_columns_step12b <- c(
  "pathway",
  "NES_GSE272769",
  "NES_GSE95233",
  "NES_GSE65682",
  "padj_GSE272769",
  "padj_GSE95233",
  "padj_GSE65682",
  "mean_abs_NES",
  "all_same_dir",
  "all_three_nominal_P_lt_0_05",
  "all_three_FDR_lt_0_05",
  "n_cohorts_nominal_P_lt_0_05",
  "n_cohorts_FDR_lt_0_05"
)

stopifnot(
  is.data.frame(hall_3),
  nrow(hall_3) == 50L,
  all(
    required_columns_step12b %in%
      colnames(hall_3)
  )
)

NES_columns_step12b <- c(
  "NES_GSE272769",
  "NES_GSE95233",
  "NES_GSE65682"
)

FDR_columns_step12b <- c(
  "padj_GSE272769",
  "padj_GSE95233",
  "padj_GSE65682"
)

stopifnot(
  all(
    vapply(
      hall_3[
        ,
        NES_columns_step12b,
        drop = FALSE
      ],
      is.numeric,
      logical(1)
    )
  ),
  all(
    vapply(
      hall_3[
        ,
        FDR_columns_step12b,
        drop = FALSE
      ],
      is.numeric,
      logical(1)
    )
  ),
  all(
    is.finite(
      as.matrix(
        hall_3[
          ,
          NES_columns_step12b,
          drop = FALSE
        ]
      )
    )
  ),
  all(
    is.finite(
      as.matrix(
        hall_3[
          ,
          FDR_columns_step12b,
          drop = FALSE
        ]
      )
    )
  )
)

# ------------------------------------------------------------
# 12B.2 Output directory
# ------------------------------------------------------------

figure_dir_step12b <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_Hallmark_triple_cohort_heatmap"
)

dir.create(
  figure_dir_step12b,
  recursive = TRUE,
  showWarnings = FALSE
)

# Remove files produced by the abandoned ggplot heatmap design.
obsolete_files_step12b <- file.path(
  figure_dir_step12b,
  c(
    paste0(
      "Figure_Hallmark_triple_cohort_",
      "NES_heatmap_moderated_t.pdf"
    ),
    paste0(
      "Figure_Hallmark_triple_cohort_",
      "NES_heatmap_moderated_t.png"
    ),
    paste0(
      "Figure_Hallmark_triple_cohort_",
      "NES_heatmap_moderated_t.tiff"
    ),
    "Figure_Hallmark_heatmap_cell_source.csv"
  )
)

obsolete_files_step12b <-
  obsolete_files_step12b[
    file.exists(obsolete_files_step12b)
  ]

if (length(obsolete_files_step12b) > 0L) {
  unlink(
    obsolete_files_step12b,
    force = TRUE
  )
}

# ------------------------------------------------------------
# 12B.3 Select the revised Top 15 pathways
# ------------------------------------------------------------

hallmark_order_step12b <- order(
  -hall_3$mean_abs_NES,
  hall_3$pathway
)

hallmark_display_wide_step12b <- head(
  hall_3[
    hallmark_order_step12b,
    ,
    drop = FALSE
  ],
  15L
)

hallmark_display_wide_step12b$
  display_rank <- seq_len(
    nrow(hallmark_display_wide_step12b)
  )

expected_top15_step12b <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_HEME_METABOLISM",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_SPERMATOGENESIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)

stopifnot(
  nrow(hallmark_display_wide_step12b) == 15L,
  identical(
    as.character(
      hallmark_display_wide_step12b$pathway
    ),
    expected_top15_step12b
  ),
  sum(hall_3$all_same_dir) == 23L,
  sum(
    hall_3$
      all_three_nominal_P_lt_0_05
  ) == 8L,
  sum(
    hall_3$
      all_three_FDR_lt_0_05
  ) == 7L,
  sum(
    hallmark_display_wide_step12b$
      all_same_dir
  ) == 12L,
  sum(
    hallmark_display_wide_step12b$
      all_three_FDR_lt_0_05
  ) == 7L
)

# ------------------------------------------------------------
# 12B.4 Construct the display matrix
# ------------------------------------------------------------

hallmark_matrix_step12b <- as.matrix(
  hallmark_display_wide_step12b[
    ,
    NES_columns_step12b,
    drop = FALSE
  ]
)

storage.mode(
  hallmark_matrix_step12b
) <- "numeric"

rownames(
  hallmark_matrix_step12b
) <- gsub(
  "^HALLMARK_",
  "",
  hallmark_display_wide_step12b$pathway
)

colnames(
  hallmark_matrix_step12b
) <- c(
  "GSE272769",
  "GSE95233 D01",
  "GSE65682\n(cohort-adjusted)"
)

# No row reversal is performed.
# The highest-ranked pathway therefore appears at the top.

stopifnot(
  nrow(hallmark_matrix_step12b) == 15L,
  ncol(hallmark_matrix_step12b) == 3L,
  all(
    is.finite(
      hallmark_matrix_step12b
    )
  ),
  identical(
    rownames(hallmark_matrix_step12b)[1],
    "E2F_TARGETS"
  ),
  identical(
    rownames(hallmark_matrix_step12b)[15],
    "ESTROGEN_RESPONSE_LATE"
  )
)

# ------------------------------------------------------------
# 12B.5 Original symmetric blue-white-red colour design
# ------------------------------------------------------------

maximum_absolute_NES_step12b <- max(
  abs(hallmark_matrix_step12b),
  na.rm = TRUE
)

heatmap_colours_step12b <- colorRampPalette(
  c(
    "#4E79A7",
    "white",
    "#E15759"
  )
)(100)

heatmap_breaks_step12b <- seq(
  -maximum_absolute_NES_step12b,
  maximum_absolute_NES_step12b,
  length.out = 101
)

# ------------------------------------------------------------
# 12B.6 Create the pheatmap object
# ------------------------------------------------------------

p_hallmark_heatmap_step12b <- pheatmap::pheatmap(
  hallmark_matrix_step12b,
  
  color = heatmap_colours_step12b,
  breaks = heatmap_breaks_step12b,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  display_numbers = TRUE,
  number_format = "%.2f",
  
  fontsize = 10,
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 10,
  
  angle_col = 0,
  border_color = NA,
  
  main = "GSEA Hallmark - Triple-cohort",
  
  legend = TRUE,
  treeheight_row = 0,
  treeheight_col = 0,
  
  silent = TRUE
)

stopifnot(
  !is.null(
    p_hallmark_heatmap_step12b$gtable
  )
)

# Display in the RStudio Plot pane.
grid::grid.newpage()

grid::grid.draw(
  p_hallmark_heatmap_step12b$gtable
)

# ------------------------------------------------------------
# 12B.7 Export figure source data and audit metrics
# ------------------------------------------------------------

hallmark_matrix_export_step12b <- data.frame(
  pathway = rownames(
    hallmark_matrix_step12b
  ),
  
  hallmark_matrix_step12b,
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)

hallmark_heatmap_metrics_step12b <- data.frame(
  Metric = c(
    "Hallmark pathways evaluated",
    "All-three directionally concordant",
    "All-three nominal P below 0.05",
    "All-three FDR below 0.05",
    "Pathways displayed",
    paste(
      "Displayed pathways",
      "directionally concordant"
    ),
    paste(
      "Displayed pathways with",
      "all-three FDR below 0.05"
    )
  ),
  
  Value = c(
    nrow(hall_3),
    
    sum(
      hall_3$all_same_dir
    ),
    
    sum(
      hall_3$
        all_three_nominal_P_lt_0_05
    ),
    
    sum(
      hall_3$
        all_three_FDR_lt_0_05
    ),
    
    nrow(
      hallmark_display_wide_step12b
    ),
    
    sum(
      hallmark_display_wide_step12b$
        all_same_dir
    ),
    
    sum(
      hallmark_display_wide_step12b$
        all_three_FDR_lt_0_05
    )
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  hallmark_display_wide_step12b,
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_top15_selection_source.csv"
  ),
  row.names = FALSE
)

write.csv(
  hallmark_matrix_export_step12b,
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_heatmap_matrix.csv"
  ),
  row.names = FALSE
)

write.csv(
  hallmark_heatmap_metrics_step12b,
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_heatmap_metrics.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12B.8 Save PDF, PNG and TIFF files
# ------------------------------------------------------------

figure_width_step12b <- 8

figure_height_step12b <- min(
  8,
  nrow(hallmark_matrix_step12b) *
    0.35 + 2
)

figure_base_step12b <- file.path(
  figure_dir_step12b,
  "Figure_GSEA_triple_heatmap_moderated_t"
)

# PDF
if (capabilities("cairo")) {
  
  ai_pdf_device_step12(
    filename = paste0(
      figure_base_step12b,
      ".pdf"
    ),
    width = figure_width_step12b,
    height = figure_height_step12b
  )
  
} else {
  
  grDevices::pdf(
    file = paste0(
      figure_base_step12b,
      ".pdf"
    ),
    width = figure_width_step12b,
    height = figure_height_step12b
  )
}

grid::grid.newpage()

grid::grid.draw(
  p_hallmark_heatmap_step12b$gtable
)

grDevices::dev.off()

# PNG
grDevices::png(
  filename = paste0(
    figure_base_step12b,
    ".png"
  ),
  width = figure_width_step12b,
  height = figure_height_step12b,
  units = "in",
  res = 300,
  bg = "white"
)

grid::grid.newpage()

grid::grid.draw(
  p_hallmark_heatmap_step12b$gtable
)

grDevices::dev.off()

# TIFF
grDevices::tiff(
  filename = paste0(
    figure_base_step12b,
    ".tiff"
  ),
  width = figure_width_step12b,
  height = figure_height_step12b,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)

grid::grid.newpage()

grid::grid.draw(
  p_hallmark_heatmap_step12b$gtable
)

grDevices::dev.off()

saveRDS(
  p_hallmark_heatmap_step12b,
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_heatmap_plot_object.rds"
  )
)

# ------------------------------------------------------------
# 12B.9 Output-file audit
# ------------------------------------------------------------

expected_figure_files_step12b <- c(
  paste0(
    figure_base_step12b,
    ".pdf"
  ),
  paste0(
    figure_base_step12b,
    ".png"
  ),
  paste0(
    figure_base_step12b,
    ".tiff"
  ),
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_top15_selection_source.csv"
  ),
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_heatmap_matrix.csv"
  ),
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_heatmap_metrics.csv"
  ),
  file.path(
    figure_dir_step12b,
    "Figure_Hallmark_heatmap_plot_object.rds"
  )
)

output_file_audit_step12b <- data.frame(
  file = expected_figure_files_step12b,
  
  exists = file.exists(
    expected_figure_files_step12b
  ),
  
  size_bytes = ifelse(
    file.exists(
      expected_figure_files_step12b
    ),
    file.info(
      expected_figure_files_step12b
    )$size,
    NA_real_
  ),
  
  stringsAsFactors = FALSE
)

print(
  output_file_audit_step12b,
  row.names = FALSE
)

stopifnot(
  all(
    output_file_audit_step12b$exists
  ),
  all(
    output_file_audit_step12b$
      size_bytes > 0
  )
)

# ------------------------------------------------------------
# 12B.10 Console audit
# ------------------------------------------------------------

cat("\nHallmark heatmap metrics:\n")

print(
  hallmark_heatmap_metrics_step12b,
  row.names = FALSE
)

cat(
  "\nDisplayed pathways and revised NES values:\n"
)

print(
  hallmark_matrix_export_step12b,
  row.names = FALSE,
  digits = 5
)

cat(
  sprintf(
    paste0(
      "\nMaximum absolute NES: %.5f\n",
      "Output directory:\n%s\n"
    ),
    maximum_absolute_NES_step12b,
    figure_dir_step12b
  )
)

cat(
  "\nSTEP 12B completed successfully.\n"
)

cat(
  "The original pheatmap visual structure was retained.\n"
)

cat(
  "The highest-ranked pathway is displayed at the top.\n"
)

cat(
  "GSE65682 is labelled as cohort-adjusted.\n"
)

cat(
  "The heatmap uses the revised moderated-t GSEA results.\n"
)

cat(
  "No GSEA analysis was rerun.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "No multipanel figure was assembled.\n"
)

cat(
  "====================================================\n"
)


# ============================================================
# Step 12C: Pairwise Hallmark NES concordance plots
# Original plot design using revised moderated-t results
# ============================================================

cat(
  "\n========== STEP 12C: PAIRWISE HALLMARK NES PLOTS ==========\n"
)

suppressPackageStartupMessages(
  library(ggplot2)
)

# ------------------------------------------------------------
# 12C.1 Required-object and schema audit
# ------------------------------------------------------------

stopifnot(
  exists("REVISION_DIR"),
  exists("hall_3"),
  is.data.frame(hall_3),
  nrow(hall_3) == 50L
)

required_columns_step12c <- c(
  "pathway",
  "NES_GSE272769",
  "NES_GSE95233",
  "NES_GSE65682"
)

stopifnot(
  all(
    required_columns_step12c %in%
      colnames(hall_3)
  ),
  all(
    is.finite(
      as.matrix(
        hall_3[
          ,
          required_columns_step12c[-1L],
          drop = FALSE
        ]
      )
    )
  )
)

# ------------------------------------------------------------
# 12C.2 Output directory
# ------------------------------------------------------------

figure_dir_step12c <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_Hallmark_pairwise_NES_concordance"
)

dir.create(
  figure_dir_step12c,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 12C.3 Prespecify the three pairwise comparisons
# ------------------------------------------------------------

pair_specification_step12c <- data.frame(
  comparison = c(
    "GSE272769 vs GSE95233 D01",
    "GSE272769 vs GSE65682",
    "GSE95233 D01 vs GSE65682"
  ),
  
  x_column = c(
    "NES_GSE272769",
    "NES_GSE272769",
    "NES_GSE95233"
  ),
  
  y_column = c(
    "NES_GSE95233",
    "NES_GSE65682",
    "NES_GSE65682"
  ),
  
  x_label = c(
    "GSE272769 NES",
    "GSE272769 NES",
    "GSE95233 D01 NES"
  ),
  
  y_label = c(
    "GSE95233 D01 NES",
    "GSE65682 cohort-adjusted NES",
    "GSE65682 cohort-adjusted NES"
  ),
  
  file_stem = c(
    "GSE272769_vs_GSE95233_D01",
    "GSE272769_vs_GSE65682",
    "GSE95233_D01_vs_GSE65682"
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 12C.4 Common symmetric axis range
# ------------------------------------------------------------

all_NES_values_step12c <- unlist(
  hall_3[
    ,
    required_columns_step12c[-1L],
    drop = FALSE
  ],
  use.names = FALSE
)

maximum_absolute_NES_step12c <- max(
  abs(all_NES_values_step12c),
  na.rm = TRUE
)

common_axis_limit_step12c <- (
  ceiling(
    maximum_absolute_NES_step12c * 10
  ) / 10
) + 0.1

axis_breaks_step12c <- c(
  -4,
  -2,
  0,
  2,
  4
)

axis_breaks_step12c <- axis_breaks_step12c[
  axis_breaks_step12c >=
    -common_axis_limit_step12c &
    axis_breaks_step12c <=
    common_axis_limit_step12c
]

direction_colours_step12c <- c(
  "Concordant" = "#4E79A7",
  "Discordant" = "#B07AA1"
)

# ------------------------------------------------------------
# 12C.5 Construct data, calculate metrics and draw panels
# ------------------------------------------------------------

pair_data_list_step12c <- vector(
  mode = "list",
  length = nrow(
    pair_specification_step12c
  )
)

pair_plot_list_step12c <- vector(
  mode = "list",
  length = nrow(
    pair_specification_step12c
  )
)

pair_metric_list_step12c <- vector(
  mode = "list",
  length = nrow(
    pair_specification_step12c
  )
)

for (
  i in seq_len(
    nrow(pair_specification_step12c)
  )
) {
  
  current_x_column_step12c <-
    pair_specification_step12c$
    x_column[i]
  
  current_y_column_step12c <-
    pair_specification_step12c$
    y_column[i]
  
  current_x_step12c <- as.numeric(
    getElement(
      hall_3,
      current_x_column_step12c
    )
  )
  
  current_y_step12c <- as.numeric(
    getElement(
      hall_3,
      current_y_column_step12c
    )
  )
  
  current_valid_step12c <-
    is.finite(current_x_step12c) &
    is.finite(current_y_step12c)
  
  current_pair_data_step12c <- data.frame(
    comparison =
      pair_specification_step12c$
      comparison[i],
    
    pathway =
      hall_3$pathway[
        current_valid_step12c
      ],
    
    x_cohort =
      current_x_column_step12c,
    
    y_cohort =
      current_y_column_step12c,
    
    x_NES =
      current_x_step12c[
        current_valid_step12c
      ],
    
    y_NES =
      current_y_step12c[
        current_valid_step12c
      ],
    
    stringsAsFactors = FALSE
  )
  
  stopifnot(
    nrow(
      current_pair_data_step12c
    ) == 50L
  )
  
  current_pair_data_step12c$
    direction_status <- ifelse(
      sign(
        current_pair_data_step12c$x_NES
      ) ==
        sign(
          current_pair_data_step12c$y_NES
        ),
      "Concordant",
      "Discordant"
    )
  
  current_pair_data_step12c$
    direction_status <- factor(
      current_pair_data_step12c$
        direction_status,
      levels = c(
        "Concordant",
        "Discordant"
      )
    )
  
  current_total_step12c <- nrow(
    current_pair_data_step12c
  )
  
  current_concordant_step12c <- sum(
    current_pair_data_step12c$
      direction_status == "Concordant"
  )
  
  current_fraction_step12c <-
    current_concordant_step12c /
    current_total_step12c
  
  current_pearson_step12c <- cor(
    current_pair_data_step12c$x_NES,
    current_pair_data_step12c$y_NES,
    method = "pearson"
  )
  
  current_spearman_step12c <- cor(
    current_pair_data_step12c$x_NES,
    current_pair_data_step12c$y_NES,
    method = "spearman"
  )
  
  pair_metric_list_step12c[i] <- list(
    data.frame(
      comparison =
        pair_specification_step12c$
        comparison[i],
      
      pathways =
        current_total_step12c,
      
      same_direction =
        current_concordant_step12c,
      
      same_direction_fraction =
        current_fraction_step12c,
      
      Pearson_r =
        current_pearson_step12c,
      
      Spearman_rho =
        current_spearman_step12c,
      
      stringsAsFactors = FALSE
    )
  )
  
  pair_data_list_step12c[i] <- list(
    current_pair_data_step12c
  )
  
  current_annotation_step12c <- sprintf(
    paste0(
      "Same direction: %d/%d (%.0f%%)\n",
      "Cross-pathway Pearson r = %.3f"
    ),
    current_concordant_step12c,
    current_total_step12c,
    100 * current_fraction_step12c,
    current_pearson_step12c
  )
  
  current_plot_step12c <- ggplot(
    current_pair_data_step12c,
    aes(
      x = x_NES,
      y = y_NES,
      color = direction_status
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      color = "grey70"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      color = "grey70"
    ) +
    geom_point(
      size = 2.6,
      alpha = 0.84
    ) +
    annotate(
      geom = "label",
      x = -Inf,
      y = Inf,
      label =
        current_annotation_step12c,
      hjust = -0.05,
      vjust = 1.08,
      size = 3.0,
      color = "grey20",
      fill = "white",
      linewidth = 0.25,
      label.padding = grid::unit(
        0.18,
        "lines"
      )
    ) +
    scale_color_manual(
      values =
        direction_colours_step12c,
      breaks = c(
        "Concordant",
        "Discordant"
      ),
      drop = FALSE,
      name = NULL
    ) +
    scale_x_continuous(
      breaks = axis_breaks_step12c
    ) +
    scale_y_continuous(
      breaks = axis_breaks_step12c
    ) +
    coord_equal(
      xlim = c(
        -common_axis_limit_step12c,
        common_axis_limit_step12c
      ),
      ylim = c(
        -common_axis_limit_step12c,
        common_axis_limit_step12c
      ),
      clip = "off"
    ) +
    labs(
      title =
        pair_specification_step12c$
        comparison[i],
      
      x =
        pair_specification_step12c$
        x_label[i],
      
      y =
        pair_specification_step12c$
        y_label[i]
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid = element_blank(),
      
      plot.title = element_text(
        size = 11.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 6
        )
      ),
      
      axis.title = element_text(
        size = 10
      ),
      
      axis.text = element_text(
        size = 9,
        color = "grey15"
      ),
      
      legend.position = "bottom",
      
      legend.text = element_text(
        size = 8.8
      ),
      
      legend.key.width = grid::unit(
        0.6,
        "cm"
      ),
      
      plot.margin = margin(
        t = 8,
        r = 10,
        b = 8,
        l = 8
      )
    )
  
  pair_plot_list_step12c[i] <- list(
    current_plot_step12c
  )
}

names(
  pair_data_list_step12c
) <- pair_specification_step12c$file_stem

names(
  pair_plot_list_step12c
) <- pair_specification_step12c$file_stem

# ------------------------------------------------------------
# 12C.6 Combine source data and metrics
# ------------------------------------------------------------

pairwise_source_data_step12c <- do.call(
  rbind,
  pair_data_list_step12c
)

rownames(
  pairwise_source_data_step12c
) <- NULL

pairwise_metrics_step12c <- do.call(
  rbind,
  pair_metric_list_step12c
)

rownames(
  pairwise_metrics_step12c
) <- NULL

# Expected revised results from Step 10D.
stopifnot(
  identical(
    as.integer(
      pairwise_metrics_step12c$
        same_direction
    ),
    c(
      34L,
      36L,
      26L
    )
  ),
  max(
    abs(
      pairwise_metrics_step12c$
        same_direction_fraction -
        c(
          0.68,
          0.72,
          0.52
        )
    )
  ) < 0.000001,
  max(
    abs(
      pairwise_metrics_step12c$
        Pearson_r -
        c(
          0.543472,
          0.739185,
          0.282979
        )
    )
  ) < 0.00001,
  max(
    abs(
      pairwise_metrics_step12c$
        Spearman_rho -
        c(
          0.532197,
          0.814742,
          0.300264
        )
    )
  ) < 0.00001
)

# ------------------------------------------------------------
# 12C.7 Export source data and metrics
# ------------------------------------------------------------

write.csv(
  pairwise_source_data_step12c,
  file.path(
    figure_dir_step12c,
    "Figure_GSEA_pairwise_source_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_metrics_step12c,
  file.path(
    figure_dir_step12c,
    "Figure_GSEA_pairwise_metrics.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12C.8 Save three standalone panels
# ------------------------------------------------------------

pdf_device_step12c <- ai_pdf_device_step12

for (
  i in seq_len(
    nrow(pair_specification_step12c)
  )
) {
  
  current_file_stem_step12c <-
    pair_specification_step12c$
    file_stem[i]
  
  current_plot_step12c <- getElement(
    pair_plot_list_step12c,
    current_file_stem_step12c
  )
  
  current_file_base_step12c <- file.path(
    figure_dir_step12c,
    paste0(
      "Figure_GSEA_pairwise_",
      current_file_stem_step12c
    )
  )
  
  ggsave(
    filename = paste0(
      current_file_base_step12c,
      ".pdf"
    ),
    plot = current_plot_step12c,
    device = pdf_device_step12c,
    width = 5.5,
    height = 5.5,
    units = "in",
    bg = "white"
  )
  
  ggsave(
    filename = paste0(
      current_file_base_step12c,
      ".png"
    ),
    plot = current_plot_step12c,
    width = 5.5,
    height = 5.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(
      current_file_base_step12c,
      ".tiff"
    ),
    plot = current_plot_step12c,
    device = "tiff",
    width = 5.5,
    height = 5.5,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

saveRDS(
  list(
    plots =
      pair_plot_list_step12c,
    
    source_data =
      pairwise_source_data_step12c,
    
    metrics =
      pairwise_metrics_step12c,
    
    pair_specification =
      pair_specification_step12c
  ),
  file.path(
    figure_dir_step12c,
    "Figure_GSEA_pairwise_plot_objects.rds"
  )
)

# ------------------------------------------------------------
# 12C.9 Display the first panel in RStudio
# ------------------------------------------------------------

print(
  getElement(
    pair_plot_list_step12c,
    "GSE272769_vs_GSE95233_D01"
  )
)

# ------------------------------------------------------------
# 12C.10 Output-file audit
# ------------------------------------------------------------

expected_image_files_step12c <- unlist(
  lapply(
    pair_specification_step12c$
      file_stem,
    function(current_stem) {
      
      current_base <- file.path(
        figure_dir_step12c,
        paste0(
          "Figure_GSEA_pairwise_",
          current_stem
        )
      )
      
      paste0(
        current_base,
        c(
          ".pdf",
          ".png",
          ".tiff"
        )
      )
    }
  ),
  use.names = FALSE
)

expected_output_files_step12c <- c(
  expected_image_files_step12c,
  
  file.path(
    figure_dir_step12c,
    "Figure_GSEA_pairwise_source_data.csv"
  ),
  
  file.path(
    figure_dir_step12c,
    "Figure_GSEA_pairwise_metrics.csv"
  ),
  
  file.path(
    figure_dir_step12c,
    "Figure_GSEA_pairwise_plot_objects.rds"
  )
)

output_file_audit_step12c <- data.frame(
  file = expected_output_files_step12c,
  
  exists = file.exists(
    expected_output_files_step12c
  ),
  
  size_bytes = ifelse(
    file.exists(
      expected_output_files_step12c
    ),
    file.info(
      expected_output_files_step12c
    )$size,
    NA_real_
  ),
  
  stringsAsFactors = FALSE
)

print(
  output_file_audit_step12c,
  row.names = FALSE
)

stopifnot(
  all(
    output_file_audit_step12c$exists
  ),
  all(
    output_file_audit_step12c$
      size_bytes > 0
  )
)

# ------------------------------------------------------------
# 12C.11 Console audit
# ------------------------------------------------------------

cat("\nPairwise Hallmark NES metrics:\n")

print(
  pairwise_metrics_step12c,
  row.names = FALSE,
  digits = 6
)

cat(
  sprintf(
    paste0(
      "\nCommon symmetric axis limit: +/- %.2f\n",
      "Output directory:\n%s\n"
    ),
    common_axis_limit_step12c,
    figure_dir_step12c
  )
)

cat(
  "\nSTEP 12C completed successfully.\n"
)

cat(
  "Three standalone panels were generated.\n"
)

cat(
  "No multipanel figure was assembled.\n"
)

cat(
  "All three rankings used limma moderated t statistics.\n"
)

cat(
  "No GSEA analysis was rerun.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)



# ============================================================
# Step 12D: Cohort-adjusted MARS-stratified logFC heatmap
# Original pheatmap design using canonical Step 7 results
# ============================================================

cat(
  "\n========== STEP 12D: MARS-STRATIFIED LOGFC HEATMAP ==========\n"
)

suppressPackageStartupMessages(
  library(pheatmap)
)

# ------------------------------------------------------------
# 12D.1 Required-object and schema audit
# ------------------------------------------------------------

required_objects_step12d <- c(
  "REVISION_DIR",
  "per_endo_fc",
  "per_endo_t",
  "per_endo_P",
  "per_endo_FDR",
  "endo_results",
  "module_map"
)

missing_objects_step12d <- required_objects_step12d[
  !vapply(
    required_objects_step12d,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects_step12d) > 0L) {
  stop(
    sprintf(
      "Missing required objects: %s",
      paste(
        missing_objects_step12d,
        collapse = ", "
      )
    )
  )
}

mars_levels_step12d <- c(
  "Mars1",
  "Mars2",
  "Mars3",
  "Mars4"
)

stopifnot(
  is.matrix(per_endo_fc),
  nrow(per_endo_fc) == 30L,
  ncol(per_endo_fc) == 4L,
  identical(
    colnames(per_endo_fc),
    mars_levels_step12d
  ),
  identical(
    dimnames(per_endo_t),
    dimnames(per_endo_fc)
  ),
  identical(
    dimnames(per_endo_P),
    dimnames(per_endo_fc)
  ),
  identical(
    dimnames(per_endo_FDR),
    dimnames(per_endo_fc)
  ),
  all(
    is.finite(per_endo_fc)
  ),
  all(
    is.finite(per_endo_t)
  ),
  all(
    is.finite(per_endo_P)
  ),
  all(
    is.finite(per_endo_FDR)
  )
)

# Verify the canonical Step 7 directional summary.
endo_results_step12d <- endo_results[
  match(
    mars_levels_step12d,
    endo_results$Endotype
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    as.character(
      endo_results_step12d$Endotype
    ),
    mars_levels_step12d
  ),
  identical(
    as.integer(
      endo_results_step12d$agree
    ),
    c(
      28L,
      30L,
      13L,
      28L
    )
  ),
  identical(
    as.integer(
      endo_results_step12d$total
    ),
    rep(
      30L,
      4L
    )
  )
)

# ------------------------------------------------------------
# 12D.2 Output directory
# ------------------------------------------------------------

figure_dir_step12d <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_MARS_stratified_logFC_heatmap"
)

dir.create(
  figure_dir_step12d,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 12D.3 Assign genes to the fixed functional modules
# ------------------------------------------------------------

genes_step12d <- rownames(
  per_endo_fc
)

stopifnot(
  length(genes_step12d) == 30L,
  !anyDuplicated(genes_step12d)
)

module_step12d <- rep(
  "Unclassified",
  length(genes_step12d)
)

names(
  module_step12d
) <- genes_step12d

matched_module_step12d <- as.character(
  module_map[
    genes_step12d
  ]
)

module_step12d[
  !is.na(matched_module_step12d)
] <- matched_module_step12d[
  !is.na(matched_module_step12d)
]

default_module_levels_step12d <- c(
  "Neutrophil degranulation",
  "Cell cycle / proliferation",
  "Inflammatory / Down",
  "Other",
  "Unclassified"
)

if (
  exists("MODULE_LEVELS") &&
  length(MODULE_LEVELS) > 0L
) {
  
  module_levels_step12d <- unique(
    c(
      as.character(MODULE_LEVELS),
      "Unclassified"
    )
  )
  
} else {
  
  module_levels_step12d <-
    default_module_levels_step12d
}

missing_module_levels_step12d <- setdiff(
  unique(module_step12d),
  module_levels_step12d
)

module_levels_step12d <- c(
  module_levels_step12d,
  missing_module_levels_step12d
)

module_factor_step12d <- factor(
  module_step12d,
  levels = module_levels_step12d
)

stopifnot(
  !anyNA(module_factor_step12d)
)

# ------------------------------------------------------------
# 12D.4 Sort genes by module and gene symbol
# ------------------------------------------------------------

gene_order_step12d <- order(
  module_factor_step12d,
  genes_step12d
)

heatmap_matrix_step12d <- per_endo_fc[
  gene_order_step12d,
  mars_levels_step12d,
  drop = FALSE
]

ordered_modules_step12d <-
  module_factor_step12d[
    gene_order_step12d
  ]

row_annotation_step12d <- data.frame(
  Module = ordered_modules_step12d,
  row.names = rownames(
    heatmap_matrix_step12d
  )
)

row_annotation_step12d$Module <- droplevels(
  row_annotation_step12d$Module
)

stopifnot(
  identical(
    rownames(row_annotation_step12d),
    rownames(heatmap_matrix_step12d)
  )
)

# ------------------------------------------------------------
# 12D.5 Module colours
# ------------------------------------------------------------

module_colours_step12d <- c(
  "Neutrophil degranulation" =
    "#E15759",
  
  "Cell cycle / proliferation" =
    "#4E79A7",
  
  "Inflammatory / Down" =
    "#59A14F",
  
  "Other" =
    "#B07AA1",
  
  "Unclassified" =
    "grey70"
)

modules_without_colour_step12d <- setdiff(
  levels(
    row_annotation_step12d$Module
  ),
  names(module_colours_step12d)
)

if (
  length(
    modules_without_colour_step12d
  ) > 0L
) {
  
  additional_module_colours_step12d <- rep(
    "grey70",
    length(
      modules_without_colour_step12d
    )
  )
  
  names(
    additional_module_colours_step12d
  ) <- modules_without_colour_step12d
  
  module_colours_step12d <- c(
    module_colours_step12d,
    additional_module_colours_step12d
  )
}

annotation_colours_step12d <- list(
  Module =
    module_colours_step12d[
      levels(
        row_annotation_step12d$Module
      )
    ]
)

# ------------------------------------------------------------
# 12D.6 Symmetric logFC colour scale
# ------------------------------------------------------------

maximum_absolute_logFC_step12d <- max(
  abs(heatmap_matrix_step12d),
  na.rm = TRUE
)
colour_scale_limit_step12d <- ceiling(
  maximum_absolute_logFC_step12d * 2
) / 2

stopifnot(
  colour_scale_limit_step12d == 1.5
)

heatmap_colours_step12d <- colorRampPalette(
  c(
    "#4E79A7",
    "white",
    "#E15759"
  )
)(100)

heatmap_breaks_step12d <- seq(
  -colour_scale_limit_step12d,
  colour_scale_limit_step12d,
  length.out = 101
)

legend_breaks_step12d <- seq(
  -colour_scale_limit_step12d,
  colour_scale_limit_step12d,
  by = 0.5
)

legend_labels_step12d <- sprintf(
  "%.1f",
  legend_breaks_step12d
)

# Display the cohort-adjusted logFC values without
# direction-reversal stars.
heatmap_numbers_step12d <- matrix(
  sprintf(
    "%.2f",
    heatmap_matrix_step12d
  ),
  nrow = nrow(
    heatmap_matrix_step12d
  ),
  ncol = ncol(
    heatmap_matrix_step12d
  ),
  dimnames = dimnames(
    heatmap_matrix_step12d
  )
)

# ------------------------------------------------------------
# 12D.7 Create the pheatmap object
# ------------------------------------------------------------

p_mars_heatmap_step12d <- pheatmap::pheatmap(
  heatmap_matrix_step12d,
  
  color = heatmap_colours_step12d,
  breaks = heatmap_breaks_step12d,
  
  legend_breaks = legend_breaks_step12d,
  legend_labels = legend_labels_step12d,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  annotation_row =
    row_annotation_step12d,
  
  annotation_colors =
    annotation_colours_step12d,
  
  annotation_names_row = FALSE,
  
  display_numbers =
    heatmap_numbers_step12d,
  
  number_color = "grey20",
  
  fontsize = 10,
  fontsize_number = 6.8,
  fontsize_row = 8.5,
  fontsize_col = 10,
  
  cellheight = 13,
  cellwidth = 42,
  
  angle_col = 0,
  border_color = NA,
  
  main = paste0(
    "Cohort-adjusted log2FC across MARS classes\n",
    "(non-survivors vs survivors)"
  ),
  
  na_col = "grey90",
  
  treeheight_row = 0,
  treeheight_col = 0,
  
  silent = TRUE
)

stopifnot(
  !is.null(
    p_mars_heatmap_step12d$gtable
  )
)

# Display in the RStudio Plot pane.
grid::grid.newpage()

grid::grid.draw(
  p_mars_heatmap_step12d$gtable
)

# ------------------------------------------------------------
# 12D.8 Construct figure source data
# ------------------------------------------------------------

figure_source_wide_step12d <- data.frame(
  gene = rownames(
    heatmap_matrix_step12d
  ),
  
  module = as.character(
    row_annotation_step12d$Module
  ),
  
  Mars1 =
    heatmap_matrix_step12d[
      ,
      "Mars1"
    ],
  
  Mars2 =
    heatmap_matrix_step12d[
      ,
      "Mars2"
    ],
  
  Mars3 =
    heatmap_matrix_step12d[
      ,
      "Mars3"
    ],
  
  Mars4 =
    heatmap_matrix_step12d[
      ,
      "Mars4"
    ],
  
  stringsAsFactors = FALSE
)

figure_statistics_long_step12d <- do.call(
  rbind,
  lapply(
    mars_levels_step12d,
    function(current_mars) {
      
      current_gene_order <- match(
        rownames(
          heatmap_matrix_step12d
        ),
        rownames(per_endo_fc)
      )
      
      data.frame(
        gene = rownames(
          heatmap_matrix_step12d
        ),
        
        module = as.character(
          row_annotation_step12d$Module
        ),
        
        MARS = current_mars,
        
        logFC =
          per_endo_fc[
            current_gene_order,
            current_mars
          ],
        
        moderated_t =
          per_endo_t[
            current_gene_order,
            current_mars
          ],
        
        P =
          per_endo_P[
            current_gene_order,
            current_mars
          ],
        
        FDR =
          per_endo_FDR[
            current_gene_order,
            current_mars
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  figure_statistics_long_step12d
) <- NULL

stopifnot(
  nrow(
    figure_statistics_long_step12d
  ) == 120L
)

write.csv(
  figure_source_wide_step12d,
  file.path(
    figure_dir_step12d,
    "Figure_MARS_logFC_heatmap_source.csv"
  ),
  row.names = FALSE
)

write.csv(
  figure_statistics_long_step12d,
  file.path(
    figure_dir_step12d,
    "Figure_MARS_gene_statistics_long.csv"
  ),
  row.names = FALSE
)

write.csv(
  endo_results_step12d,
  file.path(
    figure_dir_step12d,
    "Figure_MARS_direction_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12D.9 Save PDF, PNG and TIFF
# ------------------------------------------------------------

figure_width_step12d <- 7.4
figure_height_step12d <- 7.8

figure_base_step12d <- file.path(
  figure_dir_step12d,
  paste0(
    "Figure_MARS_stratified_",
    "cohort_adjusted_logFC_heatmap"
  )
)

# PDF
if (capabilities("cairo")) {
  
  ai_pdf_device_step12(
    filename = paste0(
      figure_base_step12d,
      ".pdf"
    ),
    width = figure_width_step12d,
    height = figure_height_step12d
  )
  
} else {
  
  grDevices::pdf(
    file = paste0(
      figure_base_step12d,
      ".pdf"
    ),
    width = figure_width_step12d,
    height = figure_height_step12d
  )
}

grid::grid.newpage()

grid::grid.draw(
  p_mars_heatmap_step12d$gtable
)

grDevices::dev.off()

# PNG
grDevices::png(
  filename = paste0(
    figure_base_step12d,
    ".png"
  ),
  width = figure_width_step12d,
  height = figure_height_step12d,
  units = "in",
  res = 300,
  bg = "white"
)

grid::grid.newpage()

grid::grid.draw(
  p_mars_heatmap_step12d$gtable
)

grDevices::dev.off()

# TIFF
grDevices::tiff(
  filename = paste0(
    figure_base_step12d,
    ".tiff"
  ),
  width = figure_width_step12d,
  height = figure_height_step12d,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)

grid::grid.newpage()

grid::grid.draw(
  p_mars_heatmap_step12d$gtable
)

grDevices::dev.off()

saveRDS(
  p_mars_heatmap_step12d,
  file.path(
    figure_dir_step12d,
    "Figure_MARS_heatmap_plot_object.rds"
  )
)

# ------------------------------------------------------------
# 12D.10 Output-file audit
# ------------------------------------------------------------

expected_output_files_step12d <- c(
  paste0(
    figure_base_step12d,
    ".pdf"
  ),
  
  paste0(
    figure_base_step12d,
    ".png"
  ),
  
  paste0(
    figure_base_step12d,
    ".tiff"
  ),
  
  file.path(
    figure_dir_step12d,
    "Figure_MARS_logFC_heatmap_source.csv"
  ),
  
  file.path(
    figure_dir_step12d,
    "Figure_MARS_gene_statistics_long.csv"
  ),
  
  file.path(
    figure_dir_step12d,
    "Figure_MARS_direction_summary.csv"
  ),
  
  file.path(
    figure_dir_step12d,
    "Figure_MARS_heatmap_plot_object.rds"
  )
)

output_file_audit_step12d <- data.frame(
  file = expected_output_files_step12d,
  
  exists = file.exists(
    expected_output_files_step12d
  ),
  
  size_bytes = ifelse(
    file.exists(
      expected_output_files_step12d
    ),
    file.info(
      expected_output_files_step12d
    )$size,
    NA_real_
  ),
  
  stringsAsFactors = FALSE
)

print(
  output_file_audit_step12d,
  row.names = FALSE
)

stopifnot(
  all(
    output_file_audit_step12d$exists
  ),
  all(
    output_file_audit_step12d$
      size_bytes > 0
  )
)

# ------------------------------------------------------------
# 12D.11 Console audit
# ------------------------------------------------------------

cat("\nMARS directional summary:\n")

print(
  endo_results_step12d,
  row.names = FALSE
)

cat(
  sprintf(
    paste0(
      "\nMaximum absolute cohort-adjusted logFC: %.5f\n",
      "Output directory:\n%s\n"
    ),
    maximum_absolute_logFC_step12d,
    figure_dir_step12d
  )
)

cat(
  "\nSTEP 12D completed successfully.\n"
)

cat(
  "The canonical Step 7 per_endo_fc matrix was used.\n"
)

cat(
  "No MARS-stratified differential expression was rerun.\n"
)

cat(
  "No direction-reversal stars were added.\n"
)

cat(
  "No multipanel figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)


# ============================================================
# Step 12E: Cohort-adjusted module-by-MARS interaction plots
# ============================================================

cat(
  "\n========== STEP 12E: MODULE × MARS INTERACTION PLOTS ==========\n"
)

suppressPackageStartupMessages(
  library(ggplot2)
)

# ------------------------------------------------------------
# 12E.1 Required-object and schema audit
# ------------------------------------------------------------

required_objects_step12e <- c(
  "REVISION_DIR",
  "module_interaction_global_step8",
  "module_interaction_slopes_step8"
)

missing_objects_step12e <- required_objects_step12e[
  !vapply(
    required_objects_step12e,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects_step12e) > 0L) {
  stop(
    sprintf(
      "Missing required objects: %s",
      paste(
        missing_objects_step12e,
        collapse = ", "
      )
    )
  )
}

required_global_columns_step12e <- c(
  "module",
  "submitted_genes",
  "detected_genes",
  "analyzed_samples",
  "non_survivors",
  "survivors",
  "interaction_df",
  "interaction_chisq",
  "interaction_P",
  "interaction_FDR"
)

required_slope_columns_step12e <- c(
  "module",
  "endotype",
  "log_OR_per_1SD",
  "SE",
  "OR_per_1SD",
  "CI95_lower",
  "CI95_upper",
  "Wald_Z",
  "Wald_P",
  "Wald_FDR"
)

stopifnot(
  is.data.frame(
    module_interaction_global_step8
  ),
  is.data.frame(
    module_interaction_slopes_step8
  ),
  all(
    required_global_columns_step12e %in%
      colnames(
        module_interaction_global_step8
      )
  ),
  all(
    required_slope_columns_step12e %in%
      colnames(
        module_interaction_slopes_step8
      )
  )
)

# ------------------------------------------------------------
# 12E.2 Fixed module and MARS definitions
# ------------------------------------------------------------

module_levels_step12e <- c(
  "Cell cycle / proliferation",
  "Neutrophil degranulation"
)

mars_levels_step12e <- c(
  "Mars1",
  "Mars2",
  "Mars3",
  "Mars4"
)

module_titles_step12e <- c(
  "Cell cycle / proliferation" =
    "Cell-cycle/proliferation module (M1)",
  
  "Neutrophil degranulation" =
    "Neutrophil-degranulation module (M2)"
)

module_colours_step12e <- c(
  "Cell cycle / proliferation" =
    "#4E79A7",
  
  "Neutrophil degranulation" =
    "#E15759"
)

module_file_stems_step12e <- c(
  "Cell cycle / proliferation" =
    "M1_cell_cycle_proliferation",
  
  "Neutrophil degranulation" =
    "M2_neutrophil_degranulation"
)

stopifnot(
  setequal(
    module_interaction_global_step8$module,
    module_levels_step12e
  ),
  setequal(
    module_interaction_slopes_step8$module,
    module_levels_step12e
  ),
  setequal(
    module_interaction_slopes_step8$endotype,
    mars_levels_step12e
  ),
  nrow(
    module_interaction_global_step8
  ) == 2L,
  nrow(
    module_interaction_slopes_step8
  ) == 8L
)

# ------------------------------------------------------------
# 12E.3 Numerical-result audit
# ------------------------------------------------------------

global_audit_step12e <-
  module_interaction_global_step8[
    match(
      module_levels_step12e,
      module_interaction_global_step8$module
    ),
    ,
    drop = FALSE
  ]

stopifnot(
  max(
    abs(
      global_audit_step12e$
        interaction_P -
        c(
          0.307803,
          0.396725
        )
    )
  ) < 0.00001,
  
  max(
    abs(
      global_audit_step12e$
        interaction_FDR -
        c(
          0.396725,
          0.396725
        )
    )
  ) < 0.00001,
  
  identical(
    as.integer(
      global_audit_step12e$
        interaction_df
    ),
    c(
      3L,
      3L
    )
  ),
  
  all(
    module_interaction_slopes_step8$
      OR_per_1SD > 0
  ),
  
  all(
    module_interaction_slopes_step8$
      CI95_lower > 0
  ),
  
  all(
    module_interaction_slopes_step8$
      CI95_upper >
      module_interaction_slopes_step8$
      CI95_lower
  ),
  
  sum(
    module_interaction_slopes_step8$
      Wald_FDR < 0.05
  ) == 1L
)

mars3_M1_step12e <-
  module_interaction_slopes_step8[
    module_interaction_slopes_step8$
      module ==
      "Cell cycle / proliferation" &
      module_interaction_slopes_step8$
      endotype ==
      "Mars3",
    ,
    drop = FALSE
  ]

mars3_M2_step12e <-
  module_interaction_slopes_step8[
    module_interaction_slopes_step8$
      module ==
      "Neutrophil degranulation" &
      module_interaction_slopes_step8$
      endotype ==
      "Mars3",
    ,
    drop = FALSE
  ]

stopifnot(
  nrow(mars3_M1_step12e) == 1L,
  nrow(mars3_M2_step12e) == 1L,
  
  abs(
    mars3_M1_step12e$OR_per_1SD -
      0.885813
  ) < 0.00001,
  
  abs(
    mars3_M2_step12e$OR_per_1SD -
      1.135621
  ) < 0.00001
)

# ------------------------------------------------------------
# 12E.4 Output directory
# ------------------------------------------------------------

figure_dir_step12e <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_MARS_module_interaction_forest"
)

dir.create(
  figure_dir_step12e,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 12E.5 Construct two standalone forest plots
# ------------------------------------------------------------

module_plot_list_step12e <- vector(
  mode = "list",
  length = length(
    module_levels_step12e
  )
)

module_source_list_step12e <- vector(
  mode = "list",
  length = length(
    module_levels_step12e
  )
)

for (
  i in seq_along(
    module_levels_step12e
  )
) {
  
  current_module_step12e <-
    module_levels_step12e[i]
  
  current_global_step12e <-
    module_interaction_global_step8[
      module_interaction_global_step8$
        module ==
        current_module_step12e,
      ,
      drop = FALSE
    ]
  
  current_slopes_step12e <-
    module_interaction_slopes_step8[
      module_interaction_slopes_step8$
        module ==
        current_module_step12e,
      ,
      drop = FALSE
    ]
  
  current_slopes_step12e <-
    current_slopes_step12e[
      match(
        mars_levels_step12e,
        current_slopes_step12e$endotype
      ),
      ,
      drop = FALSE
    ]
  
  stopifnot(
    nrow(current_global_step12e) == 1L,
    nrow(current_slopes_step12e) == 4L,
    identical(
      as.character(
        current_slopes_step12e$endotype
      ),
      mars_levels_step12e
    )
  )
  
  # Reverse factor levels so Mars1 appears at the top.
  current_slopes_step12e$
    endotype_display <- factor(
      current_slopes_step12e$endotype,
      levels = rev(
        mars_levels_step12e
      )
    )
  
  current_slopes_step12e$
    OR_CI_label <- sprintf(
      "%.2f (%.2f-%.2f)",
      current_slopes_step12e$
        OR_per_1SD,
      current_slopes_step12e$
        CI95_lower,
      current_slopes_step12e$
        CI95_upper
    )
  
  current_slopes_step12e$
    label_x <- pmin(
      current_slopes_step12e$
        CI95_upper * 1.09,
      5.85
    )
  
  current_colour_step12e <- unname(
    module_colours_step12e[
      current_module_step12e
    ]
  )
  
  current_title_step12e <- unname(
    module_titles_step12e[
      current_module_step12e
    ]
  )
  
  current_subtitle_step12e <- sprintf(
    paste0(
      "Global module-score \u00D7 MARS: ",
      "\u03C7\u00B2(%d) = %.2f; ",
      "P = %.3f; FDR = %.3f"
    ),
    as.integer(
      current_global_step12e$
        interaction_df
    ),
    current_global_step12e$
      interaction_chisq,
    current_global_step12e$
      interaction_P,
    current_global_step12e$
      interaction_FDR
  )
  
  current_plot_step12e <- ggplot(
    current_slopes_step12e,
    aes(
      y = endotype_display
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.55,
      color = "grey48"
    ) +
    geom_segment(
      aes(
        x = CI95_lower,
        xend = CI95_upper,
        yend = endotype_display
      ),
      linewidth = 0.85,
      color = current_colour_step12e
    ) +
    geom_point(
      aes(
        x = OR_per_1SD
      ),
      shape = 21,
      size = 4.0,
      stroke = 0.4,
      fill = current_colour_step12e,
      color = "white"
    ) +
    geom_text(
      aes(
        x = label_x,
        label = OR_CI_label
      ),
      hjust = 0,
      size = 3.15,
      color = "grey15"
    ) +
    scale_x_log10(
      limits = c(
        0.35,
        7.5
      ),
      breaks = c(
        0.5,
        1,
        2,
        4
      ),
      labels = c(
        "0.5",
        "1",
        "2",
        "4"
      )
    ) +
    scale_y_discrete(
      drop = FALSE
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    labs(
      title = current_title_step12e,
      
      subtitle =
        current_subtitle_step12e,
      
      x = paste0(
        "Odds ratio for 28-day mortality\n",
        "per 1-SD increase in module score"
      ),
      
      y = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.major.y =
        element_blank(),
      
      panel.grid.minor =
        element_blank(),
      
      plot.title = element_text(
        size = 12,
        face = "bold",
        margin = margin(
          b = 5
        )
      ),
      
      plot.subtitle = element_text(
        size = 9.4,
        lineheight = 1.08,
        margin = margin(
          b = 9
        )
      ),
      
      axis.title.x = element_text(
        size = 10
      ),
      
      axis.text.x = element_text(
        size = 9,
        color = "grey15"
      ),
      
      axis.text.y = element_text(
        size = 10,
        color = "grey10"
      ),
      
      plot.margin = margin(
        t = 8,
        r = 30,
        b = 8,
        l = 8
      )
    )
  
  module_plot_list_step12e[i] <- list(
    current_plot_step12e
  )
  
  current_slopes_step12e$
    global_interaction_chisq <-
    current_global_step12e$
    interaction_chisq
  
  current_slopes_step12e$
    global_interaction_df <-
    current_global_step12e$
    interaction_df
  
  current_slopes_step12e$
    global_interaction_P <-
    current_global_step12e$
    interaction_P
  
  current_slopes_step12e$
    global_interaction_FDR <-
    current_global_step12e$
    interaction_FDR
  
  module_source_list_step12e[i] <- list(
    current_slopes_step12e
  )
}

names(
  module_plot_list_step12e
) <- module_levels_step12e

names(
  module_source_list_step12e
) <- module_levels_step12e

# ------------------------------------------------------------
# 12E.6 Export figure source data
# ------------------------------------------------------------

module_forest_source_step12e <- do.call(
  rbind,
  module_source_list_step12e
)

rownames(
  module_forest_source_step12e
) <- NULL

write.csv(
  module_forest_source_step12e,
  file.path(
    figure_dir_step12e,
    "Figure_MARS_module_forest_source.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_interaction_global_step8,
  file.path(
    figure_dir_step12e,
    "Figure_MARS_module_global_tests.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12E.7 Save two standalone panels
# ------------------------------------------------------------

pdf_device_step12e <- ai_pdf_device_step12

for (
  i in seq_along(
    module_levels_step12e
  )
) {
  
  current_module_step12e <-
    module_levels_step12e[i]
  
  current_plot_step12e <- getElement(
    module_plot_list_step12e,
    current_module_step12e
  )
  
  current_file_stem_step12e <- unname(
    module_file_stems_step12e[
      current_module_step12e
    ]
  )
  
  current_file_base_step12e <- file.path(
    figure_dir_step12e,
    paste0(
      "Figure_MARS_interaction_",
      current_file_stem_step12e
    )
  )
  
  ggsave(
    filename = paste0(
      current_file_base_step12e,
      ".pdf"
    ),
    plot = current_plot_step12e,
    device = pdf_device_step12e,
    width = 7.0,
    height = 4.5,
    units = "in",
    bg = "white"
  )
  
  ggsave(
    filename = paste0(
      current_file_base_step12e,
      ".png"
    ),
    plot = current_plot_step12e,
    width = 7.0,
    height = 4.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(
      current_file_base_step12e,
      ".tiff"
    ),
    plot = current_plot_step12e,
    device = "tiff",
    width = 7.0,
    height = 4.5,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

saveRDS(
  list(
    plots =
      module_plot_list_step12e,
    
    source_data =
      module_forest_source_step12e,
    
    global_tests =
      module_interaction_global_step8
  ),
  file.path(
    figure_dir_step12e,
    "Figure_MARS_module_forest_plot_objects.rds"
  )
)

# ------------------------------------------------------------
# 12E.8 Display the M1 panel in RStudio
# ------------------------------------------------------------

print(
  getElement(
    module_plot_list_step12e,
    "Cell cycle / proliferation"
  )
)

# ------------------------------------------------------------
# 12E.9 Output-file audit
# ------------------------------------------------------------

expected_image_files_step12e <- unlist(
  lapply(
    module_levels_step12e,
    function(current_module) {
      
      current_stem <- unname(
        module_file_stems_step12e[
          current_module
        ]
      )
      
      current_base <- file.path(
        figure_dir_step12e,
        paste0(
          "Figure_MARS_interaction_",
          current_stem
        )
      )
      
      paste0(
        current_base,
        c(
          ".pdf",
          ".png",
          ".tiff"
        )
      )
    }
  ),
  use.names = FALSE
)

expected_output_files_step12e <- c(
  expected_image_files_step12e,
  
  file.path(
    figure_dir_step12e,
    "Figure_MARS_module_forest_source.csv"
  ),
  
  file.path(
    figure_dir_step12e,
    "Figure_MARS_module_global_tests.csv"
  ),
  
  file.path(
    figure_dir_step12e,
    "Figure_MARS_module_forest_plot_objects.rds"
  )
)

output_file_audit_step12e <- data.frame(
  file = expected_output_files_step12e,
  
  exists = file.exists(
    expected_output_files_step12e
  ),
  
  size_bytes = ifelse(
    file.exists(
      expected_output_files_step12e
    ),
    file.info(
      expected_output_files_step12e
    )$size,
    NA_real_
  ),
  
  stringsAsFactors = FALSE
)

print(
  output_file_audit_step12e,
  row.names = FALSE
)

stopifnot(
  all(
    output_file_audit_step12e$exists
  ),
  all(
    output_file_audit_step12e$
      size_bytes > 0
  )
)

# ------------------------------------------------------------
# 12E.10 Console audit
# ------------------------------------------------------------

cat("\nGlobal module-score × MARS tests:\n")

print(
  module_interaction_global_step8[
    ,
    c(
      "module",
      "interaction_df",
      "interaction_chisq",
      "interaction_P",
      "interaction_FDR"
    )
  ],
  row.names = FALSE,
  digits = 6
)

cat(
  "\nMARS-specific OR estimates:\n"
)

print(
  module_interaction_slopes_step8[
    ,
    c(
      "module",
      "endotype",
      "OR_per_1SD",
      "CI95_lower",
      "CI95_upper",
      "Wald_P",
      "Wald_FDR"
    )
  ],
  row.names = FALSE,
  digits = 6
)

cat(
  sprintf(
    "\nOutput directory:\n%s\n",
    figure_dir_step12e
  )
)

cat(
  "\nSTEP 12E completed successfully.\n"
)

cat(
  "Two standalone module panels were generated.\n"
)

cat(
  "No multipanel figure was assembled.\n"
)

cat(
  "No interaction model was rerun.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)


# ============================================================
# Step 12F: Revised Figure 3C
#           GSE65682-derived PC1 and 28-day mortality
#           Standalone panel; no multi-panel assembly
# ============================================================

cat(
  "\n========== STEP 12F: PC1 OUTCOME PANEL ==========\n"
)

# ------------------------------------------------------------
# 12F.1 Package and object audit
# ------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required but is not installed.")
}

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Package 'pROC' is required but is not installed.")
}

if (!exists("REVISION_DIR", inherits = TRUE)) {
  stop("Object 'REVISION_DIR' was not found.")
}

required_objects_step12f <- c(
  "pc1_data",
  "roc_obj",
  "pc1_descriptive_auc_step9",
  "pc1_association_step9",
  "pc1_pca_summary_step9"
)

required_objects_present_step12f <- vapply(
  required_objects_step12f,
  exists,
  logical(1),
  inherits = TRUE
)

# If the R session was restarted, restore Step 9 objects
# from the saved checkpoint without rerunning the analysis.
if (!all(required_objects_present_step12f)) {
  
  checkpoint_file_step12f <- file.path(
    REVISION_DIR,
    "10_PC1_outcome_assessment",
    "checkpoint_step9_PC1_outcome_assessment.rds"
  )
  
  if (!file.exists(checkpoint_file_step12f)) {
    stop(
      paste0(
        "Required Step 9 objects were missing and the checkpoint ",
        "file was not found:\n",
        checkpoint_file_step12f
      )
    )
  }
  
  checkpoint_step12f <- readRDS(
    checkpoint_file_step12f
  )
  
  pc1_data <- getElement(
    checkpoint_step12f,
    "pc1_data"
  )
  
  roc_obj <- getElement(
    checkpoint_step12f,
    "combined_pc1_roc"
  )
  
  pc1_descriptive_auc_step9 <- getElement(
    checkpoint_step12f,
    "descriptive_auc"
  )
  
  pc1_association_step9 <- getElement(
    checkpoint_step12f,
    "association"
  )
  
  pc1_pca_summary_step9 <- getElement(
    checkpoint_step12f,
    "pca_summary"
  )
  
  cat(
    "Step 9 objects were restored from the saved checkpoint.\n"
  )
}

stopifnot(
  is.data.frame(pc1_data),
  is.data.frame(pc1_descriptive_auc_step9),
  is.data.frame(pc1_association_step9),
  is.data.frame(pc1_pca_summary_step9),
  nrow(pc1_data) == 479L,
  length(unique(pc1_data$Sample)) == 479L,
  sum(pc1_data$Outcome == "Survivor") == 365L,
  sum(pc1_data$Outcome == "NonSurvivor") == 114L
)

# ------------------------------------------------------------
# 12F.2 Extract audited statistics
# ------------------------------------------------------------

combined_auc_row_step12f <-
  pc1_descriptive_auc_step9[
    as.character(
      getElement(
        pc1_descriptive_auc_step9,
        "analysis"
      )
    ) == "Combined",
    ,
    drop = FALSE
  ]

pc1_association_row_step12f <-
  pc1_association_step9[
    as.character(
      getElement(
        pc1_association_step9,
        "predictor"
      )
    ) == "PC1_z",
    ,
    drop = FALSE
  ]

stopifnot(
  nrow(combined_auc_row_step12f) == 1L,
  nrow(pc1_association_row_step12f) == 1L,
  nrow(pc1_pca_summary_step9) == 1L
)

auc_step12f <- as.numeric(
  getElement(
    combined_auc_row_step12f,
    "AUC"
  )
)

auc_lower_step12f <- as.numeric(
  getElement(
    combined_auc_row_step12f,
    "AUC_CI95_lower"
  )
)

auc_upper_step12f <- as.numeric(
  getElement(
    combined_auc_row_step12f,
    "AUC_CI95_upper"
  )
)

or_step12f <- as.numeric(
  getElement(
    pc1_association_row_step12f,
    "OR_per_1SD"
  )
)

or_lower_step12f <- as.numeric(
  getElement(
    pc1_association_row_step12f,
    "CI95_lower"
  )
)

or_upper_step12f <- as.numeric(
  getElement(
    pc1_association_row_step12f,
    "CI95_upper"
  )
)

wald_p_step12f <- as.numeric(
  getElement(
    pc1_association_row_step12f,
    "Wald_P"
  )
)

pc1_variance_step12f <- as.numeric(
  getElement(
    pc1_pca_summary_step9,
    "PC1_variance_explained_pct"
  )
)

# Audit against the completed Step 9 analysis.
stopifnot(
  abs(auc_step12f - 0.619106) < 1e-5,
  abs(auc_lower_step12f - 0.562267) < 1e-5,
  abs(auc_upper_step12f - 0.675945) < 1e-5,
  abs(or_step12f - 1.45242) < 1e-5,
  abs(or_lower_step12f - 1.18200) < 1e-4,
  abs(or_upper_step12f - 1.78470) < 1e-4,
  abs(wald_p_step12f - 0.000384167) < 1e-8,
  abs(pc1_variance_step12f - 51.3128) < 1e-4,
  abs(
    as.numeric(pROC::auc(roc_obj)) -
      auc_step12f
  ) < 1e-8
)

# ------------------------------------------------------------
# 12F.3 Construct ROC source data
# ------------------------------------------------------------

roc_source_step12f <- data.frame(
  specificity = as.numeric(
    getElement(
      roc_obj,
      "specificities"
    )
  ),
  sensitivity = as.numeric(
    getElement(
      roc_obj,
      "sensitivities"
    )
  ),
  threshold = as.numeric(
    getElement(
      roc_obj,
      "thresholds"
    )
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(roc_source_step12f) > 2L,
  all(
    is.finite(
      roc_source_step12f$specificity
    )
  ),
  all(
    is.finite(
      roc_source_step12f$sensitivity
    )
  ),
  all(
    roc_source_step12f$specificity >= 0 &
      roc_source_step12f$specificity <= 1
  ),
  all(
    roc_source_step12f$sensitivity >= 0 &
      roc_source_step12f$sensitivity <= 1
  )
)

roc_source_step12f <- roc_source_step12f[
  order(
    roc_source_step12f$specificity,
    -roc_source_step12f$sensitivity
  ),
  ,
  drop = FALSE
]

rownames(roc_source_step12f) <- NULL

# ------------------------------------------------------------
# 12F.4 Labels
# ------------------------------------------------------------

auc_label_step12f <- sprintf(
  "Descriptive AUC = %.3f (95%% CI %.3f-%.3f)",
  auc_step12f,
  auc_lower_step12f,
  auc_upper_step12f
)

or_label_step12f <- sprintf(
  paste0(
    "Cohort-adjusted OR per 1-SD PC1 = %.2f ",
    "(95%% CI %.2f-%.2f)"
  ),
  or_step12f,
  or_lower_step12f,
  or_upper_step12f
)

p_label_step12f <- sprintf(
  "P = %.6f",
  wald_p_step12f
)

statistic_annotation_step12f <- paste(
  auc_label_step12f,
  or_label_step12f,
  p_label_step12f,
  sep = "\n"
)

subtitle_step12f <- sprintf(
  paste0(
    "30 evaluable signature genes | ",
    "n = %d | PC1 variance explained = %.1f%%"
  ),
  nrow(pc1_data),
  pc1_variance_step12f
)

# ------------------------------------------------------------
# 12F.5 Draw the standalone ROC panel
# ------------------------------------------------------------

pc1_roc_plot_step12f <- ggplot2::ggplot(
  roc_source_step12f,
  ggplot2::aes(
    x = specificity,
    y = sensitivity
  )
) +
  ggplot2::geom_abline(
    intercept = 1,
    slope = -1,
    colour = "#B3B3B3",
    linewidth = 0.65
  ) +
  ggplot2::geom_step(
    colour = "#D95F5F",
    linewidth = 1.05,
    direction = "vh",
    lineend = "round"
  ) +
  ggplot2::annotate(
    geom = "text",
    x = 0.035,
    y = 0.045,
    label = statistic_annotation_step12f,
    hjust = 1,
    vjust = 0,
    size = 3.25,
    lineheight = 1.08,
    colour = "#202020"
  ) +
  ggplot2::scale_x_reverse(
    limits = c(1, 0),
    breaks = seq(
      1,
      0,
      by = -0.2
    ),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    expand = c(0, 0)
  ) +
  ggplot2::coord_equal(
    ratio = 1,
    clip = "off"
  ) +
  ggplot2::labs(
    title =
      "GSE65682-derived signature PC1 and 28-day mortality",
    subtitle = subtitle_step12f,
    x = "Specificity",
    y = "Sensitivity"
  ) +
  ggplot2::theme_classic(
    base_size = 12,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 14.5,
      face = "bold",
      colour = "#111111",
      margin = ggplot2::margin(
        b = 5
      )
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5,
      colour = "#333333",
      margin = ggplot2::margin(
        b = 10
      )
    ),
    axis.title = ggplot2::element_text(
      size = 11.5,
      colour = "#111111"
    ),
    axis.text = ggplot2::element_text(
      size = 10,
      colour = "#222222"
    ),
    axis.line = ggplot2::element_line(
      linewidth = 0.55,
      colour = "#333333"
    ),
    axis.ticks = ggplot2::element_line(
      linewidth = 0.45,
      colour = "#333333"
    ),
    plot.margin = ggplot2::margin(
      t = 10,
      r = 14,
      b = 10,
      l = 10
    )
  )

print(pc1_roc_plot_step12f)

# ------------------------------------------------------------
# 12F.6 Export standalone panel and source data
# ------------------------------------------------------------

FIGURE_3C_DIR_STEP12F <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_3C_PC1_outcome_association"
)

dir.create(
  FIGURE_3C_DIR_STEP12F,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_base_step12f <- file.path(
  FIGURE_3C_DIR_STEP12F,
  "Figure_3C_GSE65682_PC1_outcome_association"
)

# PDF: use the Illustrator-compatible PDF device
ggplot2::ggsave(
  filename = paste0(
    figure_base_step12j,
    ".pdf"
  ),
  plot = mars3_plot_step12j,
  device = ai_pdf_device_step12,
  width = 10.2,
  height = 6.2,
  units = "in",
  bg = "white"
)

# PNG: use the PNG device
ggplot2::ggsave(
  filename = paste0(
    figure_base_step12j,
    ".png"
  ),
  plot = mars3_plot_step12j,
  device = "png",
  width = 10.2,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

# TIFF: use the TIFF device
ggplot2::ggsave(
  filename = paste0(
    figure_base_step12j,
    ".tiff"
  ),
  plot = mars3_plot_step12j,
  device = "tiff",
  width = 10.2,
  height = 6.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(
    filename = paste0(
      figure_base_step12f,
      ".svg"
    ),
    plot = pc1_roc_plot_step12f,
    device = svglite::svglite,
    width = 5.7,
    height = 5.5,
    units = "in",
    bg = "white"
  )
}

write.csv(
  roc_source_step12f,
  file.path(
    FIGURE_3C_DIR_STEP12F,
    "Figure_3C_source_data_ROC_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc1_descriptive_auc_step9,
  file.path(
    FIGURE_3C_DIR_STEP12F,
    "Figure_3C_descriptive_AUC_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  pc1_association_step9,
  file.path(
    FIGURE_3C_DIR_STEP12F,
    "Figure_3C_cohort_adjusted_PC1_association.csv"
  ),
  row.names = FALSE
)

figure_metadata_step12f <- data.frame(
  item = c(
    "Analyzed samples",
    "Survivors",
    "Non-survivors",
    "Evaluable signature genes",
    "PC1 variance explained (%)",
    "Descriptive combined AUC",
    "AUC 95% CI lower",
    "AUC 95% CI upper",
    "Cohort-adjusted OR per 1-SD PC1",
    "OR 95% CI lower",
    "OR 95% CI upper",
    "Wald P"
  ),
  value = c(
    nrow(pc1_data),
    sum(pc1_data$Outcome == "Survivor"),
    sum(pc1_data$Outcome == "NonSurvivor"),
    30L,
    pc1_variance_step12f,
    auc_step12f,
    auc_lower_step12f,
    auc_upper_step12f,
    or_step12f,
    or_lower_step12f,
    or_upper_step12f,
    wald_p_step12f
  ),
  stringsAsFactors = FALSE
)

write.csv(
  figure_metadata_step12f,
  file.path(
    FIGURE_3C_DIR_STEP12F,
    "Figure_3C_plot_metadata.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    plot = pc1_roc_plot_step12f,
    roc_source = roc_source_step12f,
    descriptive_auc =
      pc1_descriptive_auc_step9,
    cohort_adjusted_association =
      pc1_association_step9,
    figure_metadata =
      figure_metadata_step12f
  ),
  file.path(
    FIGURE_3C_DIR_STEP12F,
    "Figure_3C_plot_objects.rds"
  )
)

output_files_step12f <- list.files(
  FIGURE_3C_DIR_STEP12F,
  full.names = TRUE
)

stopifnot(
  file.exists(
    paste0(
      figure_base_step12f,
      ".pdf"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12f,
      ".png"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12f,
      ".tiff"
    )
  )
)

cat(
  "\nFigure 3C statistics:\n"
)

print(
  figure_metadata_step12f,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nExported files:\n"
)

print(output_files_step12f)

cat(
  "\nSTEP 12F completed successfully.\n"
)

cat(
  "The existing Step 9 PC1 and ROC objects were reused.\n"
)

cat(
  "No predictive model was refitted.\n"
)

cat(
  "No multi-panel figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)


# ============================================================
# Step 12G: Revised Figure 2C right panel
#           ChEA 2022 results for three downregulated genes
#           Standalone panel; no multi-panel assembly
# ============================================================

cat(
  "\n========== STEP 12G: DOWNREGULATED TF FDR PANEL ==========\n"
)

# ------------------------------------------------------------
# 12G.1 Package and project audit
# ------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required but is not installed.")
}

PROJECT_DIR_STEP12G <- if (
  exists("PROJECT_DIR", inherits = TRUE)
) {
  PROJECT_DIR
} else {
  "/home/sunshine/predicate/test03"
}

REVISION_DIR_STEP12G <- if (
  exists("REVISION_DIR", inherits = TRUE)
) {
  REVISION_DIR
} else {
  file.path(
    PROJECT_DIR_STEP12G,
    "results_meta",
    "08_priority1_revision"
  )
}

stopifnot(
  dir.exists(PROJECT_DIR_STEP12G),
  dir.exists(REVISION_DIR_STEP12G)
)

# ------------------------------------------------------------
# 12G.2 Locate the existing local ChEA result
# ------------------------------------------------------------

expected_tf_file_step12g <- file.path(
  PROJECT_DIR_STEP12G,
  "results_meta",
  "05_triple_cohort",
  "TF_ChEA_2022_Down.csv"
)

if (file.exists(expected_tf_file_step12g)) {
  
  tf_file_step12g <- expected_tf_file_step12g
  
} else {
  
  candidate_tf_files_step12g <- list.files(
    PROJECT_DIR_STEP12G,
    pattern = "^TF_ChEA_2022_Down\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )
  
  if (length(candidate_tf_files_step12g) != 1L) {
    
    cat(
      "\nCandidate local ChEA downregulated files:\n"
    )
    
    print(candidate_tf_files_step12g)
    
    stop(
      paste0(
        "Expected exactly one local ",
        "TF_ChEA_2022_Down.csv file, but found ",
        length(candidate_tf_files_step12g),
        "."
      )
    )
  }
  
  tf_file_step12g <- candidate_tf_files_step12g[1L]
}

cat(
  "\nUsing local TF result file:\n",
  tf_file_step12g,
  "\n",
  sep = ""
)

# No Enrichr request is made here.
tf_down_raw_step12g <- read.csv(
  tf_file_step12g,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 12G.3 Schema and numeric audit
# ------------------------------------------------------------

required_tf_columns_step12g <- c(
  "Term",
  "Overlap",
  "P.value",
  "Adjusted.P.value",
  "Odds.Ratio",
  "Combined.Score",
  "Genes"
)

missing_tf_columns_step12g <- setdiff(
  required_tf_columns_step12g,
  colnames(tf_down_raw_step12g)
)

if (length(missing_tf_columns_step12g) > 0L) {
  stop(
    sprintf(
      "Missing ChEA columns: %s",
      paste(
        missing_tf_columns_step12g,
        collapse = ", "
      )
    )
  )
}

tf_down_raw_step12g$P.value <- suppressWarnings(
  as.numeric(tf_down_raw_step12g$P.value)
)

tf_down_raw_step12g$Adjusted.P.value <- suppressWarnings(
  as.numeric(tf_down_raw_step12g$Adjusted.P.value)
)

tf_down_raw_step12g$Odds.Ratio <- suppressWarnings(
  as.numeric(tf_down_raw_step12g$Odds.Ratio)
)

tf_down_raw_step12g$Combined.Score <- suppressWarnings(
  as.numeric(tf_down_raw_step12g$Combined.Score)
)

tf_down_raw_step12g$Overlap_n <- suppressWarnings(
  as.integer(
    sub(
      "/.*$",
      "",
      tf_down_raw_step12g$Overlap
    )
  )
)

tf_down_raw_step12g$TF_symbol <- sub(
  "\\s+.*$",
  "",
  trimws(tf_down_raw_step12g$Term)
)

tf_down_clean_step12g <- tf_down_raw_step12g[
  !is.na(tf_down_raw_step12g$Adjusted.P.value) &
    is.finite(tf_down_raw_step12g$Adjusted.P.value) &
    tf_down_raw_step12g$Adjusted.P.value > 0 &
    !is.na(tf_down_raw_step12g$Combined.Score) &
    is.finite(tf_down_raw_step12g$Combined.Score) &
    !is.na(tf_down_raw_step12g$TF_symbol) &
    nzchar(tf_down_raw_step12g$TF_symbol),
  ,
  drop = FALSE
]

rownames(tf_down_clean_step12g) <- NULL

tested_term_n_step12g <- nrow(
  tf_down_clean_step12g
)

significant_term_n_step12g <- sum(
  tf_down_clean_step12g$Adjusted.P.value < 0.05
)

minimum_fdr_step12g <- min(
  tf_down_clean_step12g$Adjusted.P.value
)

cat(
  "\nChEA 2022 downregulated-query audit:\n"
)

print(
  data.frame(
    tested_target_sets =
      tested_term_n_step12g,
    FDR_below_0_05 =
      significant_term_n_step12g,
    minimum_FDR =
      minimum_fdr_step12g
  ),
  row.names = FALSE,
  digits = 8
)

# Audit against the current submitted Supplementary Table S5.
stopifnot(
  tested_term_n_step12g == 124L,
  significant_term_n_step12g == 0L,
  abs(
    minimum_fdr_step12g -
      0.212063252346456
  ) < 1e-12
)

# ------------------------------------------------------------
# 12G.4 Preserve the original Combined Score bubble design
# ------------------------------------------------------------

# Preserve the original selection rule:
# one record per TF symbol, retaining the largest Combined Score.
tf_split_step12g <- split(
  tf_down_clean_step12g,
  tf_down_clean_step12g$TF_symbol
)

tf_representative_step12g <- do.call(
  rbind,
  lapply(
    tf_split_step12g,
    function(current_tf_table) {
      
      selected_index <- which.max(
        current_tf_table$Combined.Score
      )
      
      current_tf_table[
        selected_index,
        ,
        drop = FALSE
      ]
    }
  )
)

rownames(tf_representative_step12g) <- NULL

tf_representative_step12g <-
  tf_representative_step12g[
    order(
      -tf_representative_step12g$Combined.Score,
      tf_representative_step12g$TF_symbol
    ),
    ,
    drop = FALSE
  ]

tf_plot_data_step12g <- head(
  tf_representative_step12g,
  15L
)

tf_plot_data_step12g$log10_combined_score <- log10(
  tf_plot_data_step12g$Combined.Score + 1
)

tf_plot_data_step12g$FDR_status <-
  "FDR >= 0.05"

# Largest Combined Score appears at the top.
tf_plot_data_step12g$TF_symbol <- factor(
  tf_plot_data_step12g$TF_symbol,
  levels = rev(
    tf_plot_data_step12g$TF_symbol
  )
)

stopifnot(
  nrow(tf_plot_data_step12g) == 15L,
  all(
    tf_plot_data_step12g$Adjusted.P.value >=
      0.05
  ),
  as.character(
    tf_plot_data_step12g$TF_symbol[1L]
  ) == "FLI1"
)

# ------------------------------------------------------------
# 12G.5 Draw the revised original-style bubble panel
# ------------------------------------------------------------

maximum_x_step12g <- max(
  tf_plot_data_step12g$log10_combined_score
)

tf_down_plot_step12g <- ggplot2::ggplot(
  tf_plot_data_step12g,
  ggplot2::aes(
    x = log10_combined_score,
    y = TF_symbol,
    size = Overlap_n
  )
) +
  ggplot2::geom_point(
    shape = 21,
    fill = "#9DBBD4",
    colour = "#4E79A7",
    stroke = 0.45,
    alpha = 0.95
  ) +
  ggplot2::annotate(
    geom = "label",
    x = maximum_x_step12g * 0.98,
    y = 1.4,
    label = paste0(
      "FDR-significant target sets: 0 / ",
      tested_term_n_step12g,
      "\nMinimum FDR = ",
      sprintf(
        "%.3f",
        minimum_fdr_step12g
      )
    ),
    hjust = 1,
    vjust = 0,
    size = 3.15,
    lineheight = 1.05,
    colour = "#3F3F3F",
    fill = "white",
    linewidth = 0.35,
    label.padding = grid::unit(
      0.18,
      "lines"
    )
  ) +
  ggplot2::scale_size_continuous(
    name = "n Genes",
    range = c(
      3.0,
      7.0
    ),
    breaks = sort(
      unique(
        tf_plot_data_step12g$Overlap_n
      )
    )
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +
  ggplot2::labs(
    title =
      "ChEA 2022: Downregulated genes (3 genes)",
    subtitle = sprintf(
      paste0(
        "No target set met FDR < 0.05 ",
        "(minimum FDR = %.3f)"
      ),
      minimum_fdr_step12g
    ),
    x = expression(
      log[10]("Combined Score + 1")
    ),
    y = NULL
  ) +
  ggplot2::theme_bw(
    base_size = 11,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 13,
      face = "bold",
      colour = "#111111",
      margin = ggplot2::margin(
        b = 5
      )
    ),
    plot.subtitle = ggplot2::element_text(
      size = 9.8,
      colour = "#444444",
      margin = ggplot2::margin(
        b = 8
      )
    ),
    axis.title.x = ggplot2::element_text(
      size = 10.5,
      margin = ggplot2::margin(
        t = 6
      )
    ),
    axis.text.x = ggplot2::element_text(
      size = 9.5,
      colour = "#333333"
    ),
    axis.text.y = ggplot2::element_text(
      size = 9.7,
      colour = "#222222"
    ),
    panel.grid.major = ggplot2::element_line(
      colour = "#E6E6E6",
      linewidth = 0.35
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "#777777",
      linewidth = 0.45
    ),
    legend.title = ggplot2::element_text(
      size = 9.5
    ),
    legend.text = ggplot2::element_text(
      size = 9
    ),
    plot.margin = ggplot2::margin(
      t = 10,
      r = 12,
      b = 10,
      l = 10
    )
  )

print(tf_down_plot_step12g)
# ------------------------------------------------------------
# 12G.6 Export standalone panel and source data
# ------------------------------------------------------------

FIGURE_2C_DIR_STEP12G <- file.path(
  REVISION_DIR_STEP12G,
  "12_revised_figures",
  "Figure_2C_downregulated_TF_FDR"
)

dir.create(
  FIGURE_2C_DIR_STEP12G,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_base_step12g <- file.path(
  FIGURE_2C_DIR_STEP12G,
  "Figure_2C_ChEA2022_downregulated_FDR"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12g,
    ".pdf"
  ),
  plot = tf_down_plot_step12g,
  device = ai_pdf_device_step12,
  width = 6.2,
  height = 5.6,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12g,
    ".png"
  ),
  plot = tf_down_plot_step12g,
  width = 6.2,
  height = 5.6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12g,
    ".tiff"
  ),
  plot = tf_down_plot_step12g,
  device = "tiff",
  width = 6.2,
  height = 5.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(
    filename = paste0(
      figure_base_step12g,
      ".svg"
    ),
    plot = tf_down_plot_step12g,
    device = svglite::svglite,
    width = 6.2,
    height = 5.6,
    units = "in",
    bg = "white"
  )
}

write.csv(
  tf_plot_data_step12g,
  file.path(
    FIGURE_2C_DIR_STEP12G,
    "Figure_2C_displayed_TF_target_sets.csv"
  ),
  row.names = FALSE
)

write.csv(
  tf_down_clean_step12g,
  file.path(
    FIGURE_2C_DIR_STEP12G,
    "Figure_2C_complete_ChEA2022_downregulated_results.csv"
  ),
  row.names = FALSE
)

tf_summary_step12g <- data.frame(
  query = "Three downregulated signature genes",
  database = "ChEA 2022",
  tested_target_sets = tested_term_n_step12g,
  FDR_below_0_05 = significant_term_n_step12g,
  minimum_FDR = minimum_fdr_step12g,
  display_rule = paste0(
    "Fifteen lowest-FDR TF symbols; ",
    "one representative target set per TF symbol"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  tf_summary_step12g,
  file.path(
    FIGURE_2C_DIR_STEP12G,
    "Figure_2C_TF_FDR_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    plot = tf_down_plot_step12g,
    displayed_data = tf_plot_data_step12g,
    complete_results = tf_down_clean_step12g,
    summary = tf_summary_step12g
  ),
  file.path(
    FIGURE_2C_DIR_STEP12G,
    "Figure_2C_TF_FDR_plot_objects.rds"
  )
)

output_files_step12g <- list.files(
  FIGURE_2C_DIR_STEP12G,
  full.names = TRUE
)

stopifnot(
  file.exists(
    paste0(
      figure_base_step12g,
      ".pdf"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12g,
      ".png"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12g,
      ".tiff"
    )
  )
)

cat(
  "\nExported files:\n"
)

print(output_files_step12g)

cat(
  "\nSTEP 12G completed successfully.\n"
)

cat(
  "The existing local ChEA 2022 result was reused.\n"
)

cat(
  "No TF enrichment analysis was rerun.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "No multi-panel figure was assembled.\n"
)

cat(
  "====================================================\n"
)




# ============================================================
# Step 12H: Revised Figure 3A
#           Cohort-adjusted GSE65682 volcano plot
#           Standalone panel; no multi-panel assembly
# ============================================================

cat(
  "\n========== STEP 12H: COHORT-ADJUSTED VOLCANO ==========\n"
)

# ------------------------------------------------------------
# 12H.1 Package and object audit
# ------------------------------------------------------------

required_packages_step12h <- c(
  "ggplot2",
  "ggrepel"
)

missing_packages_step12h <- required_packages_step12h[
  !vapply(
    required_packages_step12h,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_step12h) > 0L) {
  stop(
    sprintf(
      "Missing required package(s): %s",
      paste(
        missing_packages_step12h,
        collapse = ", "
      )
    )
  )
}

if (!exists("REVISION_DIR", inherits = TRUE)) {
  stop("Object 'REVISION_DIR' was not found.")
}

# Restore Step 4 result only if the current R session no longer
# contains the cohort-adjusted DE table.
if (!exists("res_val_adjusted", inherits = TRUE)) {
  
  checkpoint_candidates_step12h <- list.files(
    REVISION_DIR,
    pattern = "^checkpoint_step4_de\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(checkpoint_candidates_step12h) != 1L) {
    
    cat(
      "\nStep 4 checkpoint candidates:\n"
    )
    
    print(checkpoint_candidates_step12h)
    
    stop(
      paste0(
        "Expected exactly one checkpoint_step4_de.rds file, ",
        "but found ",
        length(checkpoint_candidates_step12h),
        "."
      )
    )
  }
  
  checkpoint_step12h <- readRDS(
    checkpoint_candidates_step12h[1L]
  )
  
  res_val_adjusted <- getElement(
    checkpoint_step12h,
    "res_val_adjusted"
  )
  
  cat(
    "The cohort-adjusted Step 4 result was restored from checkpoint.\n"
  )
}

stopifnot(
  is.data.frame(res_val_adjusted),
  nrow(res_val_adjusted) == 11518L
)

required_columns_step12h <- c(
  "gene",
  "logFC",
  "P.Value",
  "adj.P.Val"
)

missing_columns_step12h <- setdiff(
  required_columns_step12h,
  colnames(res_val_adjusted)
)

if (length(missing_columns_step12h) > 0L) {
  stop(
    sprintf(
      "Missing DE columns: %s",
      paste(
        missing_columns_step12h,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# 12H.2 Construct the plotting table
# ------------------------------------------------------------

LFC_CUT_STEP12H <- 0.5
FDR_CUT_STEP12H <- 0.05

de_plot_data_step12h <- data.frame(
  gene = as.character(
    res_val_adjusted$gene
  ),
  logFC = as.numeric(
    res_val_adjusted$logFC
  ),
  P_value = as.numeric(
    res_val_adjusted$P.Value
  ),
  FDR = as.numeric(
    res_val_adjusted$adj.P.Val
  ),
  stringsAsFactors = FALSE
)

de_plot_data_step12h <- de_plot_data_step12h[
  !is.na(de_plot_data_step12h$gene) &
    nzchar(de_plot_data_step12h$gene) &
    is.finite(de_plot_data_step12h$logFC) &
    is.finite(de_plot_data_step12h$P_value) &
    is.finite(de_plot_data_step12h$FDR),
  ,
  drop = FALSE
]

de_plot_data_step12h$neg_log10_FDR <- -log10(
  pmax(
    de_plot_data_step12h$FDR,
    .Machine$double.xmin
  )
)

de_plot_data_step12h$Status <- "Not significant"

de_plot_data_step12h$Status[
  de_plot_data_step12h$logFC >
    LFC_CUT_STEP12H &
    de_plot_data_step12h$FDR <
    FDR_CUT_STEP12H
] <- "Upregulated"

de_plot_data_step12h$Status[
  de_plot_data_step12h$logFC <
    -LFC_CUT_STEP12H &
    de_plot_data_step12h$FDR <
    FDR_CUT_STEP12H
] <- "Downregulated"

de_plot_data_step12h$Status <- factor(
  de_plot_data_step12h$Status,
  levels = c(
    "Upregulated",
    "Downregulated",
    "Not significant"
  )
)

n_up_step12h <- sum(
  de_plot_data_step12h$Status ==
    "Upregulated"
)

n_down_step12h <- sum(
  de_plot_data_step12h$Status ==
    "Downregulated"
)

n_not_significant_step12h <- sum(
  de_plot_data_step12h$Status ==
    "Not significant"
)

cat(
  "\nCohort-adjusted volcano counts:\n"
)

print(
  data.frame(
    category = c(
      "Upregulated",
      "Downregulated",
      "Not significant"
    ),
    n = c(
      n_up_step12h,
      n_down_step12h,
      n_not_significant_step12h
    )
  ),
  row.names = FALSE
)

stopifnot(
  nrow(de_plot_data_step12h) == 11518L,
  n_up_step12h == 36L,
  n_down_step12h == 9L,
  n_not_significant_step12h ==
    11518L - 36L - 9L
)

# ------------------------------------------------------------
# 12H.3 Select labels
# ------------------------------------------------------------

upregulated_data_step12h <- de_plot_data_step12h[
  de_plot_data_step12h$Status ==
    "Upregulated",
  ,
  drop = FALSE
]

downregulated_data_step12h <- de_plot_data_step12h[
  de_plot_data_step12h$Status ==
    "Downregulated",
  ,
  drop = FALSE
]

upregulated_data_step12h <- upregulated_data_step12h[
  order(
    upregulated_data_step12h$FDR,
    -abs(
      upregulated_data_step12h$logFC
    ),
    upregulated_data_step12h$gene
  ),
  ,
  drop = FALSE
]

downregulated_data_step12h <- downregulated_data_step12h[
  order(
    downregulated_data_step12h$FDR,
    -abs(
      downregulated_data_step12h$logFC
    ),
    downregulated_data_step12h$gene
  ),
  ,
  drop = FALSE
]

# Label the 15 strongest upregulated results and all nine
# downregulated results.
label_data_step12h <- rbind(
  head(
    upregulated_data_step12h,
    15L
  ),
  head(
    downregulated_data_step12h,
    9L
  )
)

label_data_step12h <- label_data_step12h[
  !duplicated(
    label_data_step12h$gene
  ),
  ,
  drop = FALSE
]

de_plot_data_step12h$Labelled <- (
  de_plot_data_step12h$gene %in%
    label_data_step12h$gene
)

stopifnot(
  nrow(label_data_step12h) == 24L,
  sum(label_data_step12h$Status == "Upregulated") ==
    15L,
  sum(label_data_step12h$Status == "Downregulated") ==
    9L
)

# ------------------------------------------------------------
# 12H.4 Draw the cohort-adjusted volcano plot
# ------------------------------------------------------------

status_colours_step12h <- c(
  "Upregulated" = "#E15759",
  "Downregulated" = "#4E79A7",
  "Not significant" = "#D2D2D2"
)

status_labels_step12h <- c(
  "Upregulated" = sprintf(
    "Upregulated (n = %d)",
    n_up_step12h
  ),
  "Downregulated" = sprintf(
    "Downregulated (n = %d)",
    n_down_step12h
  ),
  "Not significant" = "Not significant"
)

volcano_plot_step12h <- ggplot2::ggplot(
  de_plot_data_step12h,
  ggplot2::aes(
    x = logFC,
    y = neg_log10_FDR,
    colour = Status
  )
) +
  ggplot2::geom_hline(
    yintercept = -log10(
      FDR_CUT_STEP12H
    ),
    colour = "#8C8C8C",
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  ggplot2::geom_vline(
    xintercept = c(
      -LFC_CUT_STEP12H,
      LFC_CUT_STEP12H
    ),
    colour = "#8C8C8C",
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 1.25,
    alpha = 0.78
  ) +
  ggrepel::geom_text_repel(
    data = label_data_step12h,
    ggplot2::aes(
      label = gene,
      colour = Status
    ),
    size = 3.0,
    seed = 20260831,
    box.padding = 0.32,
    point.padding = 0.16,
    min.segment.length = 0,
    segment.colour = "#A0A0A0",
    segment.linewidth = 0.35,
    max.overlaps = Inf,
    force = 0.8,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = status_colours_step12h,
    breaks = c(
      "Upregulated",
      "Downregulated",
      "Not significant"
    ),
    labels = status_labels_step12h,
    drop = FALSE
  ) +
  ggplot2::labs(
    title =
      "Cohort-adjusted differential expression in GSE65682",
    subtitle = paste0(
      "Non-survivors versus survivors | ",
      "FDR < 0.05 and |log2FC| > 0.5"
    ),
    x = expression(
      "Cohort-adjusted log"[2] *
        " fold change"
    ),
    y = expression(
      -log[10]("Benjamini-Hochberg FDR")
    ),
    colour = NULL
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 3.2,
        alpha = 1
      )
    )
  ) +
  ggplot2::theme_bw(
    base_size = 11,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 14,
      face = "bold",
      colour = "#111111",
      margin = ggplot2::margin(
        b = 5
      )
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.2,
      colour = "#333333",
      margin = ggplot2::margin(
        b = 9
      )
    ),
    axis.title = ggplot2::element_text(
      size = 10.8,
      colour = "#111111"
    ),
    axis.text = ggplot2::element_text(
      size = 9.3,
      colour = "#333333"
    ),
    panel.grid.major = ggplot2::element_line(
      colour = "#E7E7E7",
      linewidth = 0.35
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "#555555",
      linewidth = 0.55
    ),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(
      size = 9.2
    ),
    legend.key.width = grid::unit(
      0.8,
      "lines"
    ),
    legend.spacing.x = grid::unit(
      0.25,
      "cm"
    ),
    plot.margin = ggplot2::margin(
      t = 10,
      r = 14,
      b = 8,
      l = 10
    )
  ) +
  ggplot2::coord_cartesian(
    clip = "off"
  )

print(volcano_plot_step12h)

# ------------------------------------------------------------
# 12H.5 Export standalone panel and source data
# ------------------------------------------------------------

FIGURE_3A_DIR_STEP12H <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_3A_GSE65682_cohort_adjusted_volcano"
)

dir.create(
  FIGURE_3A_DIR_STEP12H,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_base_step12h <- file.path(
  FIGURE_3A_DIR_STEP12H,
  "Figure_3A_GSE65682_cohort_adjusted_volcano"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12h,
    ".pdf"
  ),
  plot = volcano_plot_step12h,
  device = ai_pdf_device_step12,
  width = 7.2,
  height = 5.8,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12h,
    ".png"
  ),
  plot = volcano_plot_step12h,
  width = 7.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12h,
    ".tiff"
  ),
  plot = volcano_plot_step12h,
  device = "tiff",
  width = 7.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(
    filename = paste0(
      figure_base_step12h,
      ".svg"
    ),
    plot = volcano_plot_step12h,
    device = svglite::svglite,
    width = 7.2,
    height = 5.8,
    units = "in",
    bg = "white"
  )
}

write.csv(
  de_plot_data_step12h,
  file.path(
    FIGURE_3A_DIR_STEP12H,
    "Figure_3A_complete_cohort_adjusted_DE_source.csv"
  ),
  row.names = FALSE
)

write.csv(
  label_data_step12h,
  file.path(
    FIGURE_3A_DIR_STEP12H,
    "Figure_3A_labelled_genes.csv"
  ),
  row.names = FALSE
)

volcano_summary_step12h <- data.frame(
  analysis =
    "GSE65682 cohort-adjusted differential expression",
  contrast =
    "Non-survivors versus survivors",
  analyzed_genes =
    nrow(de_plot_data_step12h),
  log2FC_threshold =
    LFC_CUT_STEP12H,
  FDR_threshold =
    FDR_CUT_STEP12H,
  upregulated =
    n_up_step12h,
  downregulated =
    n_down_step12h,
  not_significant =
    n_not_significant_step12h,
  model =
    "expression ~ outcome + deposited cohort",
  stringsAsFactors = FALSE
)

write.csv(
  volcano_summary_step12h,
  file.path(
    FIGURE_3A_DIR_STEP12H,
    "Figure_3A_volcano_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    plot = volcano_plot_step12h,
    complete_source =
      de_plot_data_step12h,
    labelled_genes =
      label_data_step12h,
    summary =
      volcano_summary_step12h
  ),
  file.path(
    FIGURE_3A_DIR_STEP12H,
    "Figure_3A_volcano_plot_objects.rds"
  )
)

output_files_step12h <- list.files(
  FIGURE_3A_DIR_STEP12H,
  full.names = TRUE
)

stopifnot(
  file.exists(
    paste0(
      figure_base_step12h,
      ".pdf"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12h,
      ".png"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12h,
      ".tiff"
    )
  )
)

cat(
  "\nExported files:\n"
)

print(output_files_step12h)

cat(
  "\nSTEP 12H completed successfully.\n"
)

cat(
  "The existing cohort-adjusted Step 4 DE result was reused.\n"
)

cat(
  "No differential-expression model was refitted.\n"
)

cat(
  "No multi-panel figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)


# ============================================================
# Step 12I: Revised Figure 5A
#           Direction concordance across MARS classes
#           Standalone panel; no multi-panel assembly
# ============================================================

cat(
  "\n========== STEP 12I: MARS DIRECTION CONCORDANCE ==========\n"
)

# ------------------------------------------------------------
# 12I.1 Package and object audit
# ------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required but is not installed.")
}

if (!exists("REVISION_DIR", inherits = TRUE)) {
  stop("Object 'REVISION_DIR' was not found.")
}

required_objects_step12i <- c(
  "endo_results",
  "per_endo_fc"
)

required_objects_present_step12i <- vapply(
  required_objects_step12i,
  exists,
  logical(1),
  inherits = TRUE
)

# Restore Step 7 results if the R session has been restarted.
if (!all(required_objects_present_step12i)) {
  
  checkpoint_candidates_step12i <- list.files(
    REVISION_DIR,
    pattern = "^checkpoint_step7_per_MARS_DE\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(checkpoint_candidates_step12i) != 1L) {
    
    cat(
      "\nStep 7 checkpoint candidates:\n"
    )
    
    print(checkpoint_candidates_step12i)
    
    stop(
      paste0(
        "Expected exactly one ",
        "checkpoint_step7_per_MARS_DE.rds file, but found ",
        length(checkpoint_candidates_step12i),
        "."
      )
    )
  }
  
  checkpoint_step12i <- readRDS(
    checkpoint_candidates_step12i[1L]
  )
  
  endo_results <- getElement(
    checkpoint_step12i,
    "endo_results"
  )
  
  per_endo_fc <- getElement(
    checkpoint_step12i,
    "per_endo_fc"
  )
  
  cat(
    "Step 7 objects were restored from checkpoint.\n"
  )
}

stopifnot(
  is.data.frame(endo_results),
  is.matrix(per_endo_fc),
  nrow(per_endo_fc) == 30L,
  identical(
    colnames(per_endo_fc),
    c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    )
  )
)

required_summary_columns_step12i <- c(
  "Endotype",
  "S",
  "NS",
  "agree",
  "total",
  "pct"
)

missing_summary_columns_step12i <- setdiff(
  required_summary_columns_step12i,
  colnames(endo_results)
)

if (length(missing_summary_columns_step12i) > 0L) {
  stop(
    sprintf(
      "Missing MARS summary columns: %s",
      paste(
        missing_summary_columns_step12i,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# 12I.2 Prepare and audit the plotting data
# ------------------------------------------------------------

mars_concordance_step12i <- endo_results[
  match(
    c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    ),
    endo_results$Endotype
  ),
  ,
  drop = FALSE
]

rownames(mars_concordance_step12i) <- NULL

stopifnot(
  identical(
    as.character(
      mars_concordance_step12i$Endotype
    ),
    c(
      "Mars1",
      "Mars2",
      "Mars3",
      "Mars4"
    )
  ),
  identical(
    as.integer(
      mars_concordance_step12i$agree
    ),
    c(
      28L,
      30L,
      13L,
      28L
    )
  ),
  identical(
    as.integer(
      mars_concordance_step12i$total
    ),
    rep(
      30L,
      4L
    )
  ),
  identical(
    as.integer(
      mars_concordance_step12i$S
    ),
    c(
      87L,
      138L,
      97L,
      43L
    )
  ),
  identical(
    as.integer(
      mars_concordance_step12i$NS
    ),
    c(
      45L,
      38L,
      21L,
      10L
    )
  ),
  all(
    abs(
      as.numeric(
        mars_concordance_step12i$pct
      ) -
        c(
          93.3,
          100.0,
          43.3,
          93.3
        )
    ) < 0.05
  )
)

mars_concordance_step12i$Endotype <- factor(
  mars_concordance_step12i$Endotype,
  levels = c(
    "Mars1",
    "Mars2",
    "Mars3",
    "Mars4"
  )
)

mars_concordance_step12i$label <- ifelse(
  abs(
    mars_concordance_step12i$pct -
      round(
        mars_concordance_step12i$pct
      )
  ) < 0.05,
  sprintf(
    "%d/%d\n(%.0f%%)",
    mars_concordance_step12i$agree,
    mars_concordance_step12i$total,
    mars_concordance_step12i$pct
  ),
  sprintf(
    "%d/%d\n(%.1f%%)",
    mars_concordance_step12i$agree,
    mars_concordance_step12i$total,
    mars_concordance_step12i$pct
  )
)

mars_concordance_step12i$sample_label <- sprintf(
  "S = %d; NS = %d",
  mars_concordance_step12i$S,
  mars_concordance_step12i$NS
)

cat(
  "\nFigure 5A source data:\n"
)

print(
  mars_concordance_step12i,
  row.names = FALSE
)

# ------------------------------------------------------------
# 12I.3 Draw Figure 5A
# ------------------------------------------------------------

mars_colours_step12i <- c(
  "Mars1" = "#66C2A5",
  "Mars2" = "#FC8D62",
  "Mars3" = "#8DA0CB",
  "Mars4" = "#E78AC3"
)

mars_concordance_plot_step12i <- ggplot2::ggplot(
  mars_concordance_step12i,
  ggplot2::aes(
    x = Endotype,
    y = pct,
    fill = Endotype
  )
) +
  ggplot2::geom_hline(
    yintercept = 50,
    colour = "#8C8C8C",
    linewidth = 0.55,
    linetype = "dotted"
  ) +
  ggplot2::geom_col(
    width = 0.60,
    alpha = 0.88
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = label
    ),
    vjust = -0.35,
    size = 3.65,
    lineheight = 0.95,
    colour = "#222222"
  ) +
  ggplot2::scale_fill_manual(
    values = mars_colours_step12i,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      110
    ),
    breaks = c(
      0,
      30,
      60,
      90
    ),
    expand = c(
      0,
      0
    )
  ) +
  ggplot2::labs(
    title =
      "Discovery-signature direction concordance by MARS class",
    subtitle =
      "Thirty evaluable genes; cohort-adjusted within-class estimates",
    x = NULL,
    y = "Direction concordance (%)"
  ) +
  ggplot2::theme_bw(
    base_size = 11,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 14,
      face = "bold",
      colour = "#111111",
      margin = ggplot2::margin(
        b = 5
      )
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.2,
      colour = "#333333",
      margin = ggplot2::margin(
        b = 9
      )
    ),
    axis.title.y = ggplot2::element_text(
      size = 10.8,
      colour = "#111111",
      margin = ggplot2::margin(
        r = 7
      )
    ),
    axis.text.x = ggplot2::element_text(
      size = 10.5,
      colour = "#333333",
      margin = ggplot2::margin(
        t = 5
      )
    ),
    axis.text.y = ggplot2::element_text(
      size = 9.5,
      colour = "#333333"
    ),
    panel.grid.major.x =
      ggplot2::element_blank(),
    panel.grid.minor =
      ggplot2::element_blank(),
    panel.grid.major.y =
      ggplot2::element_line(
        colour = "#E6E6E6",
        linewidth = 0.35
      ),
    panel.border = ggplot2::element_rect(
      colour = "#555555",
      linewidth = 0.55
    ),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      t = 10,
      r = 12,
      b = 10,
      l = 10
    )
  )

print(mars_concordance_plot_step12i)

# ------------------------------------------------------------
# 12I.4 Export standalone panel and source data
# ------------------------------------------------------------

FIGURE_5A_DIR_STEP12I <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_5A_MARS_direction_concordance"
)

dir.create(
  FIGURE_5A_DIR_STEP12I,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_base_step12i <- file.path(
  FIGURE_5A_DIR_STEP12I,
  "Figure_5A_MARS_direction_concordance"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12i,
    ".pdf"
  ),
  plot = mars_concordance_plot_step12i,
  device = ai_pdf_device_step12,
  width = 6.8,
  height = 5.0,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12i,
    ".png"
  ),
  plot = mars_concordance_plot_step12i,
  width = 6.8,
  height = 5.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    figure_base_step12i,
    ".tiff"
  ),
  plot = mars_concordance_plot_step12i,
  device = "tiff",
  width = 6.8,
  height = 5.0,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(
    filename = paste0(
      figure_base_step12i,
      ".svg"
    ),
    plot = mars_concordance_plot_step12i,
    device = svglite::svglite,
    width = 6.8,
    height = 5.0,
    units = "in",
    bg = "white"
  )
}

write.csv(
  mars_concordance_step12i,
  file.path(
    FIGURE_5A_DIR_STEP12I,
    "Figure_5A_source_data.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    plot = mars_concordance_plot_step12i,
    source_data =
      mars_concordance_step12i,
    colours =
      mars_colours_step12i
  ),
  file.path(
    FIGURE_5A_DIR_STEP12I,
    "Figure_5A_plot_objects.rds"
  )
)

output_files_step12i <- list.files(
  FIGURE_5A_DIR_STEP12I,
  full.names = TRUE
)

stopifnot(
  file.exists(
    paste0(
      figure_base_step12i,
      ".pdf"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12i,
      ".png"
    )
  ),
  file.exists(
    paste0(
      figure_base_step12i,
      ".tiff"
    )
  )
)

cat(
  "\nExported files:\n"
)

print(output_files_step12i)

cat(
  "\nSTEP 12I completed successfully.\n"
)

cat(
  "The existing cohort-adjusted Step 7 results were reused.\n"
)

cat(
  "No MARS-stratified DE model was refitted.\n"
)

cat(
  "No inferential significance symbols were added.\n"
)

cat(
  "No multi-panel figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)



# ============================================================
# Step 12J: Revised Figure 5B
#           Mars3 cohort-adjusted gene-level log2FC
#           Standalone panel; no multi-panel assembly
# ============================================================

cat(
  "\n========== STEP 12J: MARS3 GENE-LEVEL EFFECTS ==========\n"
)

# ------------------------------------------------------------
# 12J.1 Package and object audit
# ------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required but is not installed.")
}

if (!exists("REVISION_DIR", inherits = TRUE)) {
  stop("Object 'REVISION_DIR' was not found.")
}

if (
  !exists(
    "per_mars_signature_long",
    inherits = TRUE
  )
) {
  
  checkpoint_candidates_step12j <- list.files(
    REVISION_DIR,
    pattern = "^checkpoint_step7_per_MARS_DE\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(checkpoint_candidates_step12j) != 1L) {
    
    cat(
      "\nStep 7 checkpoint candidates:\n"
    )
    
    print(checkpoint_candidates_step12j)
    
    stop(
      paste0(
        "Expected exactly one ",
        "checkpoint_step7_per_MARS_DE.rds file, but found ",
        length(checkpoint_candidates_step12j),
        "."
      )
    )
  }
  
  checkpoint_step12j <- readRDS(
    checkpoint_candidates_step12j[1L]
  )
  
  per_mars_signature_long <- getElement(
    checkpoint_step12j,
    "per_mars_signature_long"
  )
  
  cat(
    "Step 7 Mars-stratified results were restored from checkpoint.\n"
  )
}

stopifnot(
  is.data.frame(per_mars_signature_long)
)

required_columns_step12j <- c(
  "MARS",
  "gene",
  "discovery_direction",
  "logFC",
  "directionally_concordant"
)

missing_columns_step12j <- setdiff(
  required_columns_step12j,
  colnames(per_mars_signature_long)
)

if (length(missing_columns_step12j) > 0L) {
  stop(
    sprintf(
      "Missing Mars3 source columns: %s",
      paste(
        missing_columns_step12j,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# 12J.2 Extract and classify the Mars3 estimates
# ------------------------------------------------------------

mars3_plot_data_step12j <-
  per_mars_signature_long[
    as.character(
      per_mars_signature_long$MARS
    ) == "Mars3",
    required_columns_step12j,
    drop = FALSE
  ]

rownames(mars3_plot_data_step12j) <- NULL

mars3_plot_data_step12j$gene <- as.character(
  mars3_plot_data_step12j$gene
)

mars3_plot_data_step12j$discovery_direction <-
  as.character(
    mars3_plot_data_step12j$discovery_direction
  )

mars3_plot_data_step12j$logFC <- as.numeric(
  mars3_plot_data_step12j$logFC
)

mars3_plot_data_step12j$directionally_concordant <-
  as.logical(
    mars3_plot_data_step12j$directionally_concordant
  )

stopifnot(
  nrow(mars3_plot_data_step12j) == 30L,
  length(
    unique(
      mars3_plot_data_step12j$gene
    )
  ) == 30L,
  all(
    mars3_plot_data_step12j$discovery_direction %in%
      c(
        "Up",
        "Down"
      )
  ),
  all(
    is.finite(
      mars3_plot_data_step12j$logFC
    )
  ),
  !anyNA(
    mars3_plot_data_step12j$
      directionally_concordant
  )
)

MARS3_EFFECT_THRESHOLD_STEP12J <- 0.2

mars3_plot_data_step12j$effect_group <- ifelse(
  abs(
    mars3_plot_data_step12j$logFC
  ) <= MARS3_EFFECT_THRESHOLD_STEP12J,
  "Near-zero",
  ifelse(
    mars3_plot_data_step12j$
      directionally_concordant,
    "Direction-concordant",
    "Direction-reversed"
  )
)

mars3_plot_data_step12j$effect_group <- factor(
  mars3_plot_data_step12j$effect_group,
  levels = c(
    "Near-zero",
    "Direction-concordant",
    "Direction-reversed"
  )
)

near_zero_n_step12j <- sum(
  mars3_plot_data_step12j$effect_group ==
    "Near-zero"
)

concordant_n_step12j <- sum(
  mars3_plot_data_step12j$effect_group ==
    "Direction-concordant"
)

reversed_n_step12j <- sum(
  mars3_plot_data_step12j$effect_group ==
    "Direction-reversed"
)

reversed_genes_step12j <- mars3_plot_data_step12j$gene[
  mars3_plot_data_step12j$effect_group ==
    "Direction-reversed"
]

stopifnot(
  identical(
    c(
      near_zero_n_step12j,
      concordant_n_step12j,
      reversed_n_step12j
    ),
    c(
      21L,
      8L,
      1L
    )
  ),
  identical(
    reversed_genes_step12j,
    "CA2"
  ),
  abs(
    mars3_plot_data_step12j$logFC[
      mars3_plot_data_step12j$gene ==
        "CX3CR1"
    ] -
      -0.7952066
  ) < 1e-6,
  abs(
    mars3_plot_data_step12j$logFC[
      mars3_plot_data_step12j$gene ==
        "DEFA4"
    ] -
      0.6296314
  ) < 1e-6
)

# ------------------------------------------------------------
# 12J.3 Define left-to-right gene ordering
# ------------------------------------------------------------

mars3_plot_data_step12j <-
  mars3_plot_data_step12j[
    order(
      abs(
        mars3_plot_data_step12j$logFC
      ),
      mars3_plot_data_step12j$logFC,
      mars3_plot_data_step12j$gene
    ),
    ,
    drop = FALSE
  ]

# The former top-to-bottom ordering is now shown
# from left to right in the landscape waterfall plot.
gene_levels_step12j <- as.character(
  mars3_plot_data_step12j$gene
)

mars3_plot_data_step12j$gene <- factor(
  as.character(
    mars3_plot_data_step12j$gene
  ),
  levels = gene_levels_step12j
)

mars3_plot_data_step12j$display_order_left_to_right <-
  seq_len(
    nrow(mars3_plot_data_step12j)
  )

# Four display categories preserve the two grey shades
# used for near-zero estimates.
mars3_plot_data_step12j$display_group <- ifelse(
  abs(
    mars3_plot_data_step12j$logFC
  ) <= MARS3_EFFECT_THRESHOLD_STEP12J,
  ifelse(
    mars3_plot_data_step12j$
      directionally_concordant,
    "Near-zero, concordant",
    "Near-zero, discordant"
  ),
  ifelse(
    mars3_plot_data_step12j$
      directionally_concordant,
    "Direction-concordant",
    "Direction-reversed"
  )
)

mars3_plot_data_step12j$display_group <- factor(
  mars3_plot_data_step12j$display_group,
  levels = c(
    "Near-zero, discordant",
    "Near-zero, concordant",
    "Direction-concordant",
    "Direction-reversed"
  )
)

cat(
  "\nMars3 effect-group summary:\n"
)

mars3_effect_summary_step12j <- data.frame(
  category = c(
    "Near-zero",
    "Direction-concordant",
    "Direction-reversed"
  ),
  definition = c(
    "|log2FC| <= 0.2",
    paste0(
      "|log2FC| > 0.2 and ",
      "discovery-direction concordant"
    ),
    paste0(
      "|log2FC| > 0.2 and ",
      "discovery-direction reversed"
    )
  ),
  n = c(
    near_zero_n_step12j,
    concordant_n_step12j,
    reversed_n_step12j
  ),
  stringsAsFactors = FALSE
)

print(
  mars3_effect_summary_step12j,
  row.names = FALSE
)

cat(
  "\nDirection-reversed gene(s):\n"
)

print(
  reversed_genes_step12j
)

# ------------------------------------------------------------
# 12J.4 Draw the landscape Figure 5B
#      x = gene; y = Mars3 log2FC
# ------------------------------------------------------------

mars3_colours_step12j <- c(
  "Near-zero, discordant" = "#8C8C8C",
  "Near-zero, concordant" = "#D0D0D0",
  "Direction-concordant" = "#4E79A7",
  "Direction-reversed" = "#E15759"
)

mars3_axis_limit_step12j <- max(
  0.85,
  ceiling(
    max(
      abs(
        mars3_plot_data_step12j$logFC
      )
    ) * 10
  ) / 10 + 0.05
)

mars3_subtitle_step12j <- sprintf(
  paste0(
    "Near-zero: %d (|log2FC| <= 0.2)  |  ",
    "|log2FC| > 0.2: concordant = %d; reversed = %d"
  ),
  near_zero_n_step12j,
  concordant_n_step12j,
  reversed_n_step12j
)

mars3_plot_step12j <- ggplot2::ggplot(
  mars3_plot_data_step12j,
  ggplot2::aes(
    x = gene,
    y = logFC,
    fill = display_group
  )
) +
  ggplot2::geom_col(
    width = 0.78,
    colour = NA
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    colour = "black"
  ) +
  ggplot2::geom_hline(
    yintercept = c(
      -MARS3_EFFECT_THRESHOLD_STEP12J,
      MARS3_EFFECT_THRESHOLD_STEP12J
    ),
    linewidth = 0.4,
    linetype = "dashed",
    colour = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = mars3_colours_step12j,
    drop = FALSE
  ) +
  ggplot2::scale_x_discrete(
    drop = FALSE,
    expand = ggplot2::expansion(
      add = c(
        0.5,
        0.5
      )
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      -mars3_axis_limit_step12j,
      mars3_axis_limit_step12j
    ),
    breaks = seq(
      -0.8,
      0.8,
      by = 0.2
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  ggplot2::labs(
    title = expression(
      "Mars3 gene-level " * log[2] * "FC"
    ),
    subtitle = mars3_subtitle_step12j,
    x = NULL,
    y = expression(
      log[2] * "FC in Mars3"
    )
  ) +
  ggplot2::theme_minimal(
    base_size = 10,
    base_family = "sans"
  ) +
  ggplot2::theme(
    legend.position = "none",
    
    axis.line.x = ggplot2::element_line(
      colour = "black",
      linewidth = 0.55
    ),
    axis.line.y = ggplot2::element_line(
      colour = "black",
      linewidth = 0.55
    ),
    axis.ticks.x = ggplot2::element_line(
      colour = "black",
      linewidth = 0.4
    ),
    axis.ticks.y = ggplot2::element_line(
      colour = "black",
      linewidth = 0.4
    ),
    axis.ticks.length = grid::unit(
      2.5,
      "pt"
    ),
    
    panel.grid.major.x =
      ggplot2::element_blank(),
    panel.grid.minor =
      ggplot2::element_blank(),
    panel.grid.major.y =
      ggplot2::element_line(
        colour = "#E4E4E4",
        linewidth = 0.35
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 55,
        hjust = 1,
        vjust = 1,
        size = 8,
        face = "italic",
        colour = "#333333"
      ),
    axis.text.y =
      ggplot2::element_text(
        size = 8.5,
        colour = "#333333"
      ),
    axis.title.y =
      ggplot2::element_text(
        size = 10,
        margin = ggplot2::margin(
          r = 7
        )
      ),
    plot.title =
      ggplot2::element_text(
        size = 15,
        face = "bold",
        hjust = 0
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 10,
        hjust = 0,
        margin = ggplot2::margin(
          b = 7
        )
      ),
    plot.margin =
      ggplot2::margin(
        t = 10,
        r = 12,
        b = 10,
        l = 10
      )
  )
print(
  mars3_plot_step12j
)

# ------------------------------------------------------------
# 12J.5 Export standalone panel and source data
# ------------------------------------------------------------

if (
  !exists(
    "ai_pdf_device_step12",
    mode = "function",
    inherits = TRUE
  )
) {
  stop(
    "Function 'ai_pdf_device_step12' was not found."
  )
}

FIGURE_5B_DIR_STEP12J <- file.path(
  REVISION_DIR,
  "12_revised_figures",
  "Figure_5B_Mars3_gene_level_effects"
)

dir.create(
  FIGURE_5B_DIR_STEP12J,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_base_step12j <- file.path(
  FIGURE_5B_DIR_STEP12J,
  "Figure_5B_Mars3_cohort_adjusted_logFC"
)

# PDF
ggplot2::ggsave(
  filename = paste0(
    figure_base_step12j,
    ".pdf"
  ),
  plot = mars3_plot_step12j,
  device = ai_pdf_device_step12,
  width = 10.2,
  height = 6.2,
  units = "in",
  bg = "white"
)

# PNG
ggplot2::ggsave(
  filename = paste0(
    figure_base_step12j,
    ".png"
  ),
  plot = mars3_plot_step12j,
  device = "png",
  width = 10.2,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

# TIFF
ggplot2::ggsave(
  filename = paste0(
    figure_base_step12j,
    ".tiff"
  ),
  plot = mars3_plot_step12j,
  device = "tiff",
  width = 10.2,
  height = 6.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# SVG
if (
  requireNamespace(
    "svglite",
    quietly = TRUE
  )
) {
  ggplot2::ggsave(
    filename = paste0(
      figure_base_step12j,
      ".svg"
    ),
    plot = mars3_plot_step12j,
    device = svglite::svglite,
    width = 10.2,
    height = 6.2,
    units = "in",
    bg = "white"
  )
}

# Source data
write.csv(
  mars3_plot_data_step12j,
  file.path(
    FIGURE_5B_DIR_STEP12J,
    "Figure_5B_source_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  mars3_effect_summary_step12j,
  file.path(
    FIGURE_5B_DIR_STEP12J,
    "Figure_5B_effect_group_summary.csv"
  ),
  row.names = FALSE
)

# Plot object
saveRDS(
  list(
    plot =
      mars3_plot_step12j,
    source_data =
      mars3_plot_data_step12j,
    effect_summary =
      mars3_effect_summary_step12j,
    colours =
      mars3_colours_step12j,
    threshold =
      MARS3_EFFECT_THRESHOLD_STEP12J
  ),
  file.path(
    FIGURE_5B_DIR_STEP12J,
    "Figure_5B_plot_objects.rds"
  )
)

output_files_step12j <- list.files(
  FIGURE_5B_DIR_STEP12J,
  full.names = TRUE
)

required_output_files_step12j <- c(
  paste0(
    figure_base_step12j,
    ".pdf"
  ),
  paste0(
    figure_base_step12j,
    ".png"
  ),
  paste0(
    figure_base_step12j,
    ".tiff"
  )
)

stopifnot(
  all(
    file.exists(
      required_output_files_step12j
    )
  ),
  all(
    file.info(
      required_output_files_step12j
    )$size > 0
  )
)

cat(
  "\nExported files:\n"
)

print(
  output_files_step12j
)

cat(
  "\nSTEP 12J completed successfully.\n"
)

cat(
  "The existing cohort-adjusted Step 7 Mars3 estimates were reused.\n"
)

cat(
  "No MARS-stratified model was refitted.\n"
)

cat(
  "No inferential significance symbols were added.\n"
)

cat(
  "No multi-panel figure was assembled.\n"
)

cat(
  "No online query was executed.\n"
)

cat(
  "====================================================\n"
)



# ============================================================
# Step 12K: Revised Figure 1B
# Wide vertical bar chart with all 32 genes
# ============================================================

cat(
  # ------------------------------------------------------------
  # 12K.3 Fixed display order matching the original wide panel
  # ------------------------------------------------------------
)
  expected_gene_order_step12k <- c(
    "BUB1",
    "SECTM1",
    "ASPM",
    "MKI67",
    "ANLN",
    "CEACAM6",
    "IL1B",
    "CA2",
    "CX3CR1",
    "RHAG",
    "CEACAM8",
    "ELANE",
    "MS4A3",
    "CTSG",
    "KIF11",
    "TPX2",
    "TYMS",
    "OLFM4",
    "DEFA4",
    "TCN1",
    "RNASE3",
    "HIST1H2BM",
    "CCNA2",
    "MPO",
    "HIST1H3B",
    "SERPINB10",
    "CD24",
    "CHI3L1",
    "RRM2",
    "TOP2A",
    "CKS2",
    "ABCA13"
  )
  # Remove factor attributes left from an earlier plotting run.
  i2_plot_data_step12k$gene <- trimws(
    as.character(
      i2_plot_data_step12k$gene
    )
  )
  missing_genes_step12k <- setdiff(
    expected_gene_order_step12k,
    i2_plot_data_step12k$gene
  )
  
  extra_genes_step12k <- setdiff(
    i2_plot_data_step12k$gene,
    expected_gene_order_step12k
  )
  
  if (
    length(missing_genes_step12k) > 0L ||
    length(extra_genes_step12k) > 0L
  ) {
    stop(
      paste0(
        "Figure 1B gene-set mismatch.\n",
        "Missing genes: ",
        paste(
          missing_genes_step12k,
          collapse = ", "
        ),
        "\nExtra genes: ",
        paste(
          extra_genes_step12k,
          collapse = ", "
        )
      )
    )
  }
  
  i2_plot_data_step12k <- i2_plot_data_step12k[
    match(
      expected_gene_order_step12k,
      i2_plot_data_step12k$gene
    ),
    ,
    drop = FALSE
  ]
  
  rownames(i2_plot_data_step12k) <- NULL
  
  stopifnot(
    identical(
      as.character(
        i2_plot_data_step12k$gene
      ),
      expected_gene_order_step12k
    ),
    sum(i2_plot_data_step12k$I2 == 0) == 28L,
    sum(i2_plot_data_step12k$I2 > 0) == 4L,
    identical(
      tail(
        as.character(
          i2_plot_data_step12k$gene
        ),
        4L
      ),
      c(
        "RRM2",
        "TOP2A",
        "CKS2",
        "ABCA13"
      )
    )
  )
  
  i2_plot_data_step12k$gene <- factor(
    i2_plot_data_step12k$gene,
    levels = expected_gene_order_step12k
  )
  stopifnot(
    identical(
      levels(
        i2_plot_data_step12k$gene
      ),
      expected_gene_order_step12k
    ),
    identical(
      as.character(
        i2_plot_data_step12k$gene
      ),
      expected_gene_order_step12k
    )
  )
  # ------------------------------------------------------------
  # 12K.4 I2 display categories
  # ------------------------------------------------------------
  
  i2_plot_data_step12k$I2_category <- cut(
    i2_plot_data_step12k$I2,
    breaks = c(
      -Inf,
      25,
      50,
      75,
      Inf
    ),
    right = FALSE,
    labels = c(
      "I\u00b2 < 25%",
      "25 \u2264 I\u00b2 < 50%",
      "50 \u2264 I\u00b2 < 75%",
      "I\u00b2 \u2265 75%"
    )
  )
  
  i2_plot_data_step12k$I2_category <- factor(
    i2_plot_data_step12k$I2_category,
    levels = c(
      "I\u00b2 < 25%",
      "25 \u2264 I\u00b2 < 50%",
      "50 \u2264 I\u00b2 < 75%",
      "I\u00b2 \u2265 75%"
    )
  )
  
  i2_palette_step12k <- c(
    "I\u00b2 < 25%" = "#4E79A7",
    "25 \u2264 I\u00b2 < 50%" = "#D6A43A",
    "50 \u2264 I\u00b2 < 75%" = "#E15759",
    "I\u00b2 \u2265 75%" = "#8E6C8A"
  )
  
  # ------------------------------------------------------------
  # 12K.5 Wide vertical bar chart matching the reference design
  # ------------------------------------------------------------
  
  i2_plot_step12k <- ggplot2::ggplot(
    i2_plot_data_step12k,
    ggplot2::aes(
      x = gene,
      y = I2,
      fill = I2_category
    )
  ) +
    ggplot2::geom_col(
      width = 0.90,
      colour = NA
    ) +
    
    # Heterogeneity thresholds
    ggplot2::geom_hline(
      yintercept = c(
        25,
        50,
        75
      ),
      colour = "#4F4F4F",
      linewidth = 0.42,
      linetype = "dashed"
    ) +
    
    # Visible zero reference line
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "#4E79A7",
      linewidth = 0.42,
      linetype = "dashed"
    ) +
    
    ggplot2::scale_fill_manual(
      values = i2_palette_step12k,
      name = NULL,
      drop = TRUE
    ) +
    
    ggplot2::scale_x_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(
        add = c(
          0.12,
          0.12
        )
      )
    ) +
    
    ggplot2::scale_y_continuous(
      limits = c(
        -3,
        78
      ),
      breaks = c(
        0,
        20,
        40,
        60
      ),
      minor_breaks = seq(
        0,
        70,
        by = 10
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    ggplot2::labs(
      title = "32-gene signature: I\u00b2 heterogeneity",
      subtitle = paste0(
        "k = 2 cohorts | ",
        "31/32 genes had I\u00b2 < 25%"
      ),
      x = NULL,
      y = expression(I^2 ~ "(%)")
    ) +
    
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          colour = NA
        )
      )
    ) +
    
    ggplot2::theme_bw(
      base_size = 7
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 9.8,
        face = "bold",
        hjust = 0,
        colour = "#151515",
        margin = ggplot2::margin(
          b = 2
        )
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 7.3,
        hjust = 0,
        colour = "#303030",
        margin = ggplot2::margin(
          b = 6
        )
      ),
      
      plot.title.position = "plot", 
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = "#303030",
        linewidth = 0.45
      ),
      
      panel.grid.major = ggplot2::element_line(
        colour = "#E2E2E2",
        linewidth = 0.32
      ),
      
      panel.grid.minor = ggplot2::element_line(
        colour = "#F0F0F0",
        linewidth = 0.25
      ),
      
      axis.line = ggplot2::element_blank(),
      
      axis.ticks = ggplot2::element_line(
        colour = "#303030",
        linewidth = 0.35
      ),
      
      axis.ticks.length = grid::unit(
        1.3,
        "mm"
      ),
      
      # Required rotated x-axis labels
      axis.text.x = ggplot2::element_text(
        angle = 50,
        hjust = 1,
        vjust = 1,
        size = 6.2,
        colour = "#303030",
        margin = ggplot2::margin(
          t = 2
        )
      ),
      
      axis.text.y = ggplot2::element_text(
        size = 6.6,
        colour = "#303030"
      ),
      
      axis.title.y = ggplot2::element_text(
        size = 7.4,
        margin = ggplot2::margin(
          r = 5
        )
      ),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      
      legend.text = ggplot2::element_text(
        size = 6.7,
        colour = "#303030"
      ),
      
      legend.key.width = grid::unit(
        4.5,
        "mm"
      ),
      
      legend.key.height = grid::unit(
        3.8,
        "mm"
      ),
      
      legend.margin = ggplot2::margin(
        t = 4,
        b = 0
      ),
      
      plot.margin = ggplot2::margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  
  print(i2_plot_step12k)
  # ------------------------------------------------------------
  # 12K.6 Final figure export
  # ------------------------------------------------------------
  
  stopifnot(
    exists("REVISION_DIR"),
    exists("i2_plot_step12k"),
    inherits(
      i2_plot_step12k,
      "ggplot"
    ),
    exists("i2_plot_data_step12k"),
    nrow(i2_plot_data_step12k) == 32L,
    exists("expected_gene_order_step12k"),
    identical(
      as.character(
        i2_plot_data_step12k$gene
      ),
      expected_gene_order_step12k
    ),
    exists("ai_pdf_device_step12"),
    is.function(ai_pdf_device_step12)
  )
  
  # ------------------------------------------------------------
  # 12K.6.1 Output directory
  # ------------------------------------------------------------
  
  FIGURE_1B_DIR_STEP12K <- file.path(
    REVISION_DIR,
    "12_revised_figures",
    "Figure_1B_I2_heterogeneity_vertical"
  )
  
  dir.create(
    FIGURE_1B_DIR_STEP12K,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  stopifnot(
    dir.exists(
      FIGURE_1B_DIR_STEP12K
    )
  )
  
  figure_1b_base_step12k <- file.path(
    FIGURE_1B_DIR_STEP12K,
    "Figure_1B_I2_heterogeneity_vertical"
  )
  
  # Final physical size for a wide main-figure panel
  FIG1B_WIDTH_STEP12K <- 7.2
  FIG1B_HEIGHT_STEP12K <- 3.45
  FIG1B_DPI_STEP12K <- 600
  
  # ------------------------------------------------------------
  # 12K.6.2 PDF
  # ------------------------------------------------------------
  
  figure_1b_pdf_step12k <- paste0(
    figure_1b_base_step12k,
    ".pdf"
  )
  
  ggplot2::ggsave(
    filename = figure_1b_pdf_step12k,
    plot = i2_plot_step12k,
    device = ai_pdf_device_step12,
    width = FIG1B_WIDTH_STEP12K,
    height = FIG1B_HEIGHT_STEP12K,
    units = "in",
    limitsize = FALSE
  )
  
  # Do not add compression = "lzw" to PDF.
  
  # ------------------------------------------------------------
  # 12K.6.3 PNG
  # ------------------------------------------------------------
  
  figure_1b_png_step12k <- paste0(
    figure_1b_base_step12k,
    ".png"
  )
  
  ggplot2::ggsave(
    filename = figure_1b_png_step12k,
    plot = i2_plot_step12k,
    device = "png",
    width = FIG1B_WIDTH_STEP12K,
    height = FIG1B_HEIGHT_STEP12K,
    units = "in",
    dpi = FIG1B_DPI_STEP12K,
    bg = "white",
    limitsize = FALSE
  )
  
  # ------------------------------------------------------------
  # 12K.6.4 TIFF
  # ------------------------------------------------------------
  
  figure_1b_tiff_step12k <- paste0(
    figure_1b_base_step12k,
    ".tiff"
  )
  
  ggplot2::ggsave(
    filename = figure_1b_tiff_step12k,
    plot = i2_plot_step12k,
    device = "tiff",
    width = FIG1B_WIDTH_STEP12K,
    height = FIG1B_HEIGHT_STEP12K,
    units = "in",
    dpi = FIG1B_DPI_STEP12K,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
  
  # ------------------------------------------------------------
  # 12K.6.5 SVG
  # ------------------------------------------------------------
  
  if (!requireNamespace(
    "svglite",
    quietly = TRUE
  )) {
    stop(
      "Package 'svglite' is required for SVG export."
    )
  }
  
  figure_1b_svg_step12k <- paste0(
    figure_1b_base_step12k,
    ".svg"
  )
  
  ggplot2::ggsave(
    filename = figure_1b_svg_step12k,
    plot = i2_plot_step12k,
    device = svglite::svglite,
    width = FIG1B_WIDTH_STEP12K,
    height = FIG1B_HEIGHT_STEP12K,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  
  # ------------------------------------------------------------
  # 12K.6.6 Complete 32-gene source data
  # ------------------------------------------------------------
  
  i2_source_data_export_step12k <- data.frame(
    display_order = seq_len(
      nrow(i2_plot_data_step12k)
    ),
    gene = as.character(
      i2_plot_data_step12k$gene
    ),
    I2 = i2_plot_data_step12k$I2,
    I2_category = as.character(
      i2_plot_data_step12k$I2_category
    ),
    source_order = i2_plot_data_step12k$source_order,
    stringsAsFactors = FALSE
  )
  
  stopifnot(
    nrow(i2_source_data_export_step12k) == 32L,
    !anyDuplicated(
      i2_source_data_export_step12k$gene
    ),
    identical(
      i2_source_data_export_step12k$gene,
      expected_gene_order_step12k
    )
  )
  
  figure_1b_csv_step12k <- file.path(
    FIGURE_1B_DIR_STEP12K,
    "Figure_1B_I2_source_data.csv"
  )
  
  write.csv(
    i2_source_data_export_step12k,
    figure_1b_csv_step12k,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  # ------------------------------------------------------------
  # 12K.6.7 Output audit
  # ------------------------------------------------------------
  
  expected_output_files_step12k <- c(
    figure_1b_pdf_step12k,
    figure_1b_png_step12k,
    figure_1b_tiff_step12k,
    figure_1b_svg_step12k,
    figure_1b_csv_step12k
  )
  
  stopifnot(
    all(
      file.exists(
        expected_output_files_step12k
      )
    )
  )
  
  output_file_info_step12k <- file.info(
    expected_output_files_step12k
  )
  
  output_audit_step12k <- data.frame(
    file = expected_output_files_step12k,
    size_bytes = output_file_info_step12k$size,
    modified = output_file_info_step12k$mtime,
    stringsAsFactors = FALSE
  )
  
  print(
    output_audit_step12k,
    row.names = FALSE
  )
  
  cat(
    "\nFigure 1B saved successfully.\n",
    "Output directory:\n",
    normalizePath(
      FIGURE_1B_DIR_STEP12K,
      mustWork = TRUE
    ),
    "\n\nFiles exported:\n",
    paste(
      basename(
        expected_output_files_step12k
      ),
      collapse = "\n"
    ),
    "\n====================================================\n"
  )
  
  
  # ============================================================
  # Step 12L: E3 stratification and module composition
  # ============================================================
  
  cat(
    "\n========== STEP 12L: E3 LAYER MODULE COMPOSITION ==========\n"
  )
  
  # ------------------------------------------------------------
  # 12L.1 Package audit
  # ------------------------------------------------------------
  
  required_packages_step12l <- c(
    "ggplot2",
    "dplyr",
    "tidyr",
    "svglite",
    "ragg"
  )
  
  missing_packages_step12l <- required_packages_step12l[
    !vapply(
      required_packages_step12l,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(missing_packages_step12l) > 0L) {
    stop(
      paste0(
        "Missing required packages: ",
        paste(missing_packages_step12l, collapse = ", ")
      )
    )
  }
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # ------------------------------------------------------------
  # 12L.2 Output directory
  # ------------------------------------------------------------
  
  project_dir_step12l <- if (exists("PROJECT_DIR")) {
    PROJECT_DIR
  } else {
    "/home/sunshine/predicate/test03"
  }
  
  revision_dir_step12l <- if (exists("REVISION_DIR")) {
    REVISION_DIR
  } else {
    file.path(
      project_dir_step12l,
      "results_meta",
      "08_priority1_revision"
    )
  }
  
  figure_dir_step12l <- file.path(
    revision_dir_step12l,
    "12_revised_figures",
    "Figure_E3_layer_module_composition"
  )
  
  dir.create(
    figure_dir_step12l,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ------------------------------------------------------------
  # 12L.3 Frozen source values from Supplementary Table S7
  # Sheet: Layer_module_composition
  # ------------------------------------------------------------
  
  module_order_step12l <- c(
    "M1 cell cycle/proliferation",
    "M2 neutrophil degranulation",
    "Inflammatory/downregulated",
    "Other"
  )
  
  layer_wide_step12l <- data.frame(
    Layer = c(
      "Full signature",
      "E3-overlap",
      "E3-residual"
    ),
    `M1 cell cycle/proliferation` = c(11L, 5L, 6L),
    `M2 neutrophil degranulation` = c(10L, 6L, 4L),
    `Inflammatory/downregulated` = c(3L, 2L, 1L),
    Other = c(8L, 2L, 6L),
    Total = c(32L, 15L, 17L),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # ------------------------------------------------------------
  # 12L.4 Numerical assertions
  # ------------------------------------------------------------
  
  stopifnot(
    all(
      rowSums(
        layer_wide_step12l[, module_order_step12l]
      ) == layer_wide_step12l$Total
    ),
    layer_wide_step12l$Total[
      layer_wide_step12l$Layer == "Full signature"
    ] == 32L,
    layer_wide_step12l$Total[
      layer_wide_step12l$Layer == "E3-overlap"
    ] == 15L,
    layer_wide_step12l$Total[
      layer_wide_step12l$Layer == "E3-residual"
    ] == 17L
  )
  
  full_module_counts_step12l <-
    as.integer(
      layer_wide_step12l[
        layer_wide_step12l$Layer == "Full signature",
        module_order_step12l
      ]
    )
  
  partitioned_module_counts_step12l <-
    as.integer(
      layer_wide_step12l[
        layer_wide_step12l$Layer == "E3-overlap",
        module_order_step12l
      ]
    ) +
    as.integer(
      layer_wide_step12l[
        layer_wide_step12l$Layer == "E3-residual",
        module_order_step12l
      ]
    )
  
  stopifnot(
    identical(
      full_module_counts_step12l,
      partitioned_module_counts_step12l
    ),
    round(100 * 15 / 32, 1) == 46.9,
    round(100 * 17 / 32, 1) == 53.1
  )
  
  # ------------------------------------------------------------
  # 12L.5 Long-format plotting data
  # ------------------------------------------------------------
  
  layer_long_step12l <- layer_wide_step12l |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(module_order_step12l),
      names_to = "Module",
      values_to = "Genes"
    )
  
  layer_levels_step12l <- c(
    "E3-residual",
    "E3-overlap",
    "Full signature"
  )
  
  layer_long_step12l$Layer <- factor(
    layer_long_step12l$Layer,
    levels = layer_levels_step12l
  )
  
  layer_long_step12l$Module <- factor(
    layer_long_step12l$Module,
    levels = module_order_step12l
  )
  
  # ------------------------------------------------------------
  # 12L.6 Figure 1C module-colour contract
  # ------------------------------------------------------------
  
  module_palette_step12l <- c(
    "M1 cell cycle/proliferation" = "#B07AA1",
    "M2 neutrophil degranulation" = "#E15759",
    "Inflammatory/downregulated" = "#4E79A7",
    "Other" = "#B6992D"
  )
  
  segment_text_palette_step12l <- c(
    "M1 cell cycle/proliferation" = "white",
    "M2 neutrophil degranulation" = "white",
    "Inflammatory/downregulated" = "white",
    "Other" = "#1F1F1F"
  )
  
  # ------------------------------------------------------------
  # 12L.7 Right-side annotations
  # ------------------------------------------------------------
  
  right_annotations_step12l <- data.frame(
    Layer = factor(
      c("E3-overlap", "E3-residual"),
      levels = layer_levels_step12l
    ),
    x = c(15.7, 17.7),
    label = c(
      "15/32 genes overlapping E3\n(46.9%)",
      "17/32 genes retained after\noverlap exclusion (53.1%)"
    ),
    stringsAsFactors = FALSE
  )
  
  direction_callout_step12l <- data.frame(
    Layer = factor(
      "E3-overlap",
      levels = layer_levels_step12l
    ),
    x = 21.5,
    label = paste0(
      "15/15 overlap genes retained the discovery-defined direction\n",
      "in both MESSI and MARS"
    ),
    stringsAsFactors = FALSE
  )
  
  layer_axis_labels_step12l <- c(
    "Full signature" = "Full signature\n(n = 32)",
    "E3-overlap" = "E3-overlap\n(n = 15)",
    "E3-residual" = "E3-residual\n(n = 17)"
  )
  
  # ------------------------------------------------------------
  # 12L.8 Plot
  # ------------------------------------------------------------
  
  e3_module_plot_step12l <- ggplot(
    layer_long_step12l,
    aes(
      x = Genes,
      y = Layer,
      fill = Module
    )
  ) +
    geom_col(
      width = 0.56,
      colour = "white",
      linewidth = 0.40,
      position = position_stack(reverse = TRUE)
    ) +
    geom_text(
      aes(
        label = Genes,
        colour = Module
      ),
      position = position_stack(
        vjust = 0.5,
        reverse = TRUE
      ),
      size = 3.35,
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_text(
      data = right_annotations_step12l,
      aes(
        x = x,
        y = Layer,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = 3.05,
      lineheight = 0.95,
      colour = "#252525"
    ) +
    geom_label(
      data = direction_callout_step12l,
      aes(
        x = x,
        y = Layer,
        label = label
      ),
      inherit.aes = FALSE,
      nudge_y = 0.48,
      hjust = 0.5,
      vjust = 0.5,
      size = 2.85,
      lineheight = 0.95,
      fontface = "bold",
      colour = "#3A3A3A",
      fill = "#F5F5F5",
      label.size = 0.25,
      label.padding = grid::unit(0.15, "lines"),
      label.r = grid::unit(0.08, "lines")
    ) +
    scale_fill_manual(
      values = module_palette_step12l,
      breaks = module_order_step12l,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = segment_text_palette_step12l,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = c(0, 8, 16, 24, 32),
      minor_breaks = NULL,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(
      labels = layer_axis_labels_step12l,
      expand = expansion(add = c(0.38, 0.50))
    ) +
    coord_cartesian(
      xlim = c(0, 32),
      clip = "off"
    ) +
    labs(
      title = "E3 stratification and module composition of the 32-gene signature",
      subtitle = paste0(
        "Both E3-overlap and E3-residual layers retained ",
        "cell-cycle/proliferation and neutrophil-degranulation components"
      ),
      x = "Number of genes",
      y = NULL,
      fill = "Functional module"
    ) +
    theme_classic(
      base_size = 10.2,
      base_family = "sans"
    ) +
    theme(
      plot.title = element_text(
        size = 12.4,
        face = "bold",
        colour = "black",
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9.4,
        colour = "#3F3F3F",
        margin = margin(b = 16)
      ),
      axis.title.x = element_text(
        size = 9.6,
        margin = margin(t = 8)
      ),
      axis.text.x = element_text(
        size = 8.8,
        colour = "#303030"
      ),
      axis.text.y = element_text(
        size = 9.1,
        colour = "#202020",
        lineheight = 0.93,
        margin = margin(r = 8)
      ),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(
        linewidth = 0.40,
        colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.45,
        colour = "black"
      ),
      panel.grid.major.x = element_line(
        colour = "#E4E4E4",
        linewidth = 0.35
      ),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(
        size = 8.8,
        face = "bold"
      ),
      legend.text = element_text(
        size = 8.3
      ),
      legend.key.height = grid::unit(4.2, "mm"),
      legend.key.width = grid::unit(6.2, "mm"),
      legend.margin = margin(t = 5),
      plot.margin = margin(
        t = 12,
        r = 82,
        b = 10,
        l = 12
      )
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          colour = NA
        )
      )
    )
  
  print(e3_module_plot_step12l)
  
  # ------------------------------------------------------------
  # 12L.9 Export source data
  # ------------------------------------------------------------
  
  layer_long_export_step12l <- layer_long_step12l |>
    mutate(
      Source_table = "Supplementary Table S7",
      Source_sheet = "Layer_module_composition"
    )
  
  overlap_summary_export_step12l <- data.frame(
    Metric = c(
      "Current signature genes",
      "Giannini E3 genes",
      "Overlapping genes",
      "Overlap percentage",
      "Directionally concordant overlap genes",
      "E3-residual genes",
      "E3-residual percentage"
    ),
    Value = c(
      "32",
      "199",
      "15",
      "46.9%",
      "15/15",
      "17",
      "53.1%"
    ),
    Source_sheet = "Overlap_summary",
    stringsAsFactors = FALSE
  )
  
  write.csv(
    layer_long_export_step12l,
    file.path(
      figure_dir_step12l,
      "Figure_E3_layer_module_composition_source_data.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    overlap_summary_export_step12l,
    file.path(
      figure_dir_step12l,
      "Figure_E3_overlap_summary_source_data.csv"
    ),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # 12L.10 Export figure
  # ------------------------------------------------------------
  
  figure_base_step12l <- file.path(
    figure_dir_step12l,
    "Figure_E3_layer_module_composition"
  )
  
  figure_width_step12l <- 7.20
  figure_height_step12l <- 4.35
  figure_dpi_step12l <- 600
  
  # PDF: reuse the Adobe-Illustrator-compatible device established earlier.
  if (exists("ai_pdf_device_step12")) {
    
    ai_pdf_device_step12(
      paste0(figure_base_step12l, ".pdf"),
      width = figure_width_step12l,
      height = figure_height_step12l
    )
    
  } else {
    
    grDevices::pdf(
      file = paste0(figure_base_step12l, ".pdf"),
      width = figure_width_step12l,
      height = figure_height_step12l,
      family = "Helvetica",
      useDingbats = FALSE,
      onefile = FALSE
    )
  }
  
  print(e3_module_plot_step12l)
  grDevices::dev.off()
  
  # SVG
  svglite::svglite(
    file = paste0(figure_base_step12l, ".svg"),
    width = figure_width_step12l,
    height = figure_height_step12l,
    bg = "white"
  )
  
  print(e3_module_plot_step12l)
  grDevices::dev.off()
  
  # PNG
  ragg::agg_png(
    filename = paste0(figure_base_step12l, ".png"),
    width = figure_width_step12l,
    height = figure_height_step12l,
    units = "in",
    res = figure_dpi_step12l,
    background = "white"
  )
  
  print(e3_module_plot_step12l)
  grDevices::dev.off()
  
  # TIFF
  ragg::agg_tiff(
    filename = paste0(figure_base_step12l, ".tiff"),
    width = figure_width_step12l,
    height = figure_height_step12l,
    units = "in",
    res = figure_dpi_step12l,
    background = "white",
    compression = "lzw"
  )
  
  print(e3_module_plot_step12l)
  grDevices::dev.off()
  
  # ------------------------------------------------------------
  # 12L.11 Output audit
  # ------------------------------------------------------------
  
  expected_figure_files_step12l <- paste0(
    figure_base_step12l,
    c(".pdf", ".svg", ".png", ".tiff")
  )
  
  stopifnot(
    all(file.exists(expected_figure_files_step12l)),
    file.exists(
      file.path(
        figure_dir_step12l,
        "Figure_E3_layer_module_composition_source_data.csv"
      )
    ),
    file.exists(
      file.path(
        figure_dir_step12l,
        "Figure_E3_overlap_summary_source_data.csv"
      )
    )
  )
  
  cat("\nSTEP 12L completed successfully.\n")
  cat(
    "Figure directory:\n",
    normalizePath(
      figure_dir_step12l,
      mustWork = TRUE
    ),
    "\n"
  )
  cat("No upstream analysis was rerun.\n")
  cat("No online query was executed.\n")
  cat("====================================================\n")
  
  
  # ============================================================
  # Step 12M: E3-residual external directional assessment
  # Step 12N: Functional enrichment after E3-overlap exclusion
  # ============================================================
  
  cat(
    "\n========== STEPS 12M-12N: E3-RESIDUAL PANELS ==========\n"
  )
  
  # ------------------------------------------------------------
  # 12M.1 Package audit
  # ------------------------------------------------------------
  
  required_packages_step12mn <- c(
    "ggplot2",
    "ggrepel",
    "svglite",
    "ragg"
  )
  
  missing_packages_step12mn <- required_packages_step12mn[
    !vapply(
      required_packages_step12mn,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(missing_packages_step12mn) > 0L) {
    stop(
      paste0(
        "Missing required packages: ",
        paste(missing_packages_step12mn, collapse = ", ")
      )
    )
  }
  
  library(ggplot2)
  
  # ------------------------------------------------------------
  # 12M.2 Output directories
  # ------------------------------------------------------------
  
  project_dir_step12mn <- if (exists("PROJECT_DIR")) {
    PROJECT_DIR
  } else {
    "/home/sunshine/predicate/test03"
  }
  
  revision_dir_step12mn <- if (exists("REVISION_DIR")) {
    REVISION_DIR
  } else {
    file.path(
      project_dir_step12mn,
      "results_meta",
      "08_priority1_revision"
    )
  }
  
  scatter_dir_step12m <- file.path(
    revision_dir_step12mn,
    "12_revised_figures",
    "Figure_E3_residual_external_scatter"
  )
  
  dumbbell_dir_step12n <- file.path(
    revision_dir_step12mn,
    "12_revised_figures",
    "Figure_E3_residual_functional_dumbbell"
  )
  
  dir.create(
    scatter_dir_step12m,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    dumbbell_dir_step12n,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ------------------------------------------------------------
  # 12M.3 Common style and export function
  # ------------------------------------------------------------
  
  base_family_step12mn <- "Helvetica"
  
  module_palette_step12mn <- c(
    "Cell cycle / proliferation" = "#B07AA1",
    "Neutrophil degranulation" = "#E15759",
    "Inflammatory / Down" = "#4E79A7",
    "Other" = "#B6992D"
  )
  
  theme_panel_step12mn <- function() {
    theme_classic(
      base_size = 8.5,
      base_family = base_family_step12mn
    ) +
      theme(
        axis.line = element_line(
          colour = "black",
          linewidth = 0.4
        ),
        axis.ticks = element_line(
          colour = "black",
          linewidth = 0.35
        ),
        axis.title = element_text(
          colour = "black",
          size = 8.8
        ),
        axis.text = element_text(
          colour = "black",
          size = 7.8
        ),
        plot.title = element_text(
          size = 10.0,
          face = "bold",
          hjust = 0,
          margin = margin(b = 6)
        ),
        plot.caption = element_text(
          size = 6.6,
          colour = "#444444",
          hjust = 0,
          margin = margin(t = 5)
        ),
        legend.title = element_text(
          size = 7.5,
          face = "bold"
        ),
        legend.text = element_text(
          size = 7.0
        ),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = 1),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = margin(8, 10, 6, 8)
      )
  }
  
  render_device_step12mn <- function(open_device, plot_object) {
    
    open_device()
    opened_device <- grDevices::dev.cur()
    
    on.exit(
      {
        active_devices <- grDevices::dev.list()
        
        if (
          !is.null(active_devices) &&
          opened_device %in% active_devices
        ) {
          grDevices::dev.off(which = opened_device)
        }
      },
      add = TRUE
    )
    
    print(plot_object)
    grDevices::dev.off(which = opened_device)
    
    invisible(NULL)
  }
  
  save_panel_step12mn <- function(
    plot_object,
    file_base,
    width_mm,
    height_mm,
    dpi = 600
  ) {
    
    width_in <- width_mm / 25.4
    height_in <- height_mm / 25.4
    
    # Base PDF device: avoids Cairo/Nimbus/DejaVu substitution.
    render_device_step12mn(
      function() {
        grDevices::pdf(
          file = paste0(file_base, ".pdf"),
          width = width_in,
          height = height_in,
          family = "Helvetica",
          encoding = "WinAnsi.enc",
          useDingbats = FALSE,
          onefile = TRUE
        )
      },
      plot_object
    )
    
    render_device_step12mn(
      function() {
        svglite::svglite(
          file = paste0(file_base, ".svg"),
          width = width_in,
          height = height_in
        )
      },
      plot_object
    )
    
    render_device_step12mn(
      function() {
        ragg::agg_png(
          filename = paste0(file_base, ".png"),
          width = width_in,
          height = height_in,
          units = "in",
          res = dpi,
          background = "white"
        )
      },
      plot_object
    )
    
    render_device_step12mn(
      function() {
        ragg::agg_tiff(
          filename = paste0(file_base, ".tiff"),
          width = width_in,
          height = height_in,
          units = "in",
          res = dpi,
          background = "white"
        )
      },
      plot_object
    )
  }
  
  # ============================================================
  # Step 12M: E3-residual external directional assessment
  # ============================================================
  
  # ------------------------------------------------------------
  # 12M.4 Audited values from Supplementary Table S8
  # ------------------------------------------------------------
  
  residual_gene_step12m <- data.frame(
    signature_gene = c(
      "ASPM",
      "ANLN",
      "RRM2",
      "KIF11",
      "BUB1",
      "MKI67",
      "IL1B",
      "OLFM4",
      "CEACAM8",
      "RNASE3",
      "TCN1",
      "ABCA13",
      "CD24",
      "SERPINB10",
      "CHI3L1"
    ),
    Module = c(
      rep("Cell cycle / proliferation", 6),
      "Inflammatory / Down",
      rep("Neutrophil degranulation", 4),
      rep("Other", 4)
    ),
    discovery_RE_log2FC = c(
      0.662665023361024,
      0.639377774792075,
      0.617122326036678,
      0.596860627307398,
      0.586609950257587,
      0.579807082298657,
      -0.578504119362258,
      1.230905037299770,
      1.064737185498790,
      0.751024746701887,
      0.744217507778250,
      0.990035344508998,
      0.737311983032069,
      0.684541634074426,
      0.651121136864208
    ),
    log2FC_GSE65682_cohort_adjusted = c(
      0.106004084416722,
      0.0979709224695586,
      0.313095898967182,
      0.240873373813039,
      0.0619909945744288,
      0.0850025042484521,
      -0.305319555882648,
      0.340372660621589,
      0.728620928026702,
      0.429354271819953,
      0.286390998922093,
      0.186384559574593,
      0.567227724910582,
      0.177946414600648,
      0.0351022342528173
    ),
    P_GSE65682_cohort_adjusted = c(
      0.0578470769081707,
      0.110086767955495,
      0.00900411833532524,
      0.0170529492189881,
      0.406263244176974,
      0.182917325872608,
      0.018581239005381,
      0.240026840309680,
      0.00195870280526869,
      0.00553079928726837,
      0.105811871338132,
      0.0525820706996692,
      0.00335602495424866,
      0.0565804014944423,
      0.860390037493294
    ),
    FDR_GSE65682_cohort_adjusted = c(
      0.188215432719862,
      0.269446635829721,
      0.0724733997108848,
      0.099799825616097,
      0.588050564711171,
      0.364312407844052,
      0.104883603732538,
      0.427754220008544,
      0.0351955365227532,
      0.0576723737044071,
      0.264139820995363,
      0.178823382982689,
      0.0471398724671171,
      0.186038556783610,
      0.921343664173276
    ),
    stringsAsFactors = FALSE
  )
  
  residual_gene_step12m$direction_concordant <-
    sign(residual_gene_step12m$discovery_RE_log2FC) ==
    sign(
      residual_gene_step12m$
        log2FC_GSE65682_cohort_adjusted
    )
  
  residual_gene_step12m$external_fdr_group <- ifelse(
    residual_gene_step12m$
      FDR_GSE65682_cohort_adjusted < 0.05,
    "FDR < 0.05",
    "FDR >= 0.05"
  )
  
  residual_gene_step12m$external_fdr_group <- factor(
    residual_gene_step12m$external_fdr_group,
    levels = c(
      "FDR >= 0.05",
      "FDR < 0.05"
    )
  )
  
  residual_gene_step12m$Module <- factor(
    residual_gene_step12m$Module,
    levels = names(module_palette_step12mn)
  )
  
  # ------------------------------------------------------------
  # 12M.5 Numerical assertions
  # ------------------------------------------------------------
  
  pearson_r_step12m <- cor(
    residual_gene_step12m$discovery_RE_log2FC,
    residual_gene_step12m$
      log2FC_GSE65682_cohort_adjusted,
    method = "pearson"
  )
  
  stopifnot(
    nrow(residual_gene_step12m) == 15L,
    all(residual_gene_step12m$direction_concordant),
    sum(
      residual_gene_step12m$
        FDR_GSE65682_cohort_adjusted < 0.05
    ) == 2L,
    sum(
      residual_gene_step12m$
        P_GSE65682_cohort_adjusted < 0.05
    ) == 6L,
    sum(
      abs(
        residual_gene_step12m$
          log2FC_GSE65682_cohort_adjusted
      ) > 0.5
    ) == 2L,
    abs(pearson_r_step12m - 0.732) < 0.001
  )
  
  cat("\nStep 12M Pearson correlation:\n")
  print(pearson_r_step12m)
  
  # ------------------------------------------------------------
  # 12M.6 Scatter plot: manually controlled labels
  # ------------------------------------------------------------
  
  label_positions_step12m <- data.frame(
    signature_gene = c(
      "ASPM", "ANLN", "RRM2", "KIF11", "BUB1",
      "MKI67", "IL1B", "OLFM4", "CEACAM8",
      "RNASE3", "TCN1", "ABCA13", "CD24",
      "SERPINB10", "CHI3L1"
    ),
    label_x = c(
      0.80, 0.53, 0.53, 0.49, 0.53,
      0.49, -0.61, 1.16, 1.03,
      0.82, 0.84, 1.03, 0.66,
      0.78, 0.73
    ),
    label_y = c(
      -0.04, 0.14, 0.42, 0.32, -0.20,
      -0.05, -0.41, 0.45, 0.83,
      0.61, 0.35, 0.10, 0.67,
      0.19, -0.16
    ),
    label_hjust = c(
      0, 1, 1, 1, 1,
      1, 1, 0, 0.5,
      0, 0, 0, 1,
      0, 0
    ),
    stringsAsFactors = FALSE
  )
  label_match_step12m <- match(
    residual_gene_step12m$signature_gene,
    label_positions_step12m$signature_gene
  )
  
  stopifnot(
    !anyNA(label_match_step12m),
    length(unique(label_positions_step12m$signature_gene)) == 15L
  )
  
  residual_gene_step12m$label_x <-
    label_positions_step12m$label_x[label_match_step12m]
  
  residual_gene_step12m$label_y <-
    label_positions_step12m$label_y[label_match_step12m]
  
  residual_gene_step12m$label_hjust <-
    label_positions_step12m$label_hjust[label_match_step12m]
  
  scatter_info_step12m <- paste0(
    "15/17 evaluable; 15/15 directionally concordant\n",
    "Pearson r = 0.732 (P = 0.00194)\n",
    "External FDR < 0.05: 2/15"
  )
  
  scatter_plot_step12m <- ggplot(
    residual_gene_step12m,
    aes(
      x = discovery_RE_log2FC,
      y = log2FC_GSE65682_cohort_adjusted
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      colour = "#A0A0A0",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = 0,
      colour = "#555555",
      linewidth = 0.32,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = 0,
      colour = "#555555",
      linewidth = 0.32,
      linetype = "dashed"
    ) +
    geom_segment(
      aes(
        xend = label_x,
        yend = label_y
      ),
      colour = "#969696",
      linewidth = 0.20,
      lineend = "round"
    ) +
    geom_point(
      aes(
        fill = Module,
        shape = external_fdr_group
      ),
      colour = "black",
      stroke = 0.45,
      size = 3.0
    ) +
    geom_text(
      data = residual_gene_step12m,
      aes(
        x = label_x,
        y = label_y,
        label = signature_gene,
        hjust = label_hjust
      ),
      inherit.aes = FALSE,
      size = 2.45,
      family = base_family_step12mn,
      colour = "#202020"
    ) +
    annotate(
      geom = "label",
      x = -0.71,
      y = 1.08,
      label = scatter_info_step12m,
      hjust = 0,
      vjust = 1,
      size = 2.25,
      family = base_family_step12mn,
      colour = "black",
      fill = "white",
      linewidth = 0.25,
      label.padding = grid::unit(0.18, "lines")
    ) +
    scale_fill_manual(
      name = "Module",
      values = module_palette_step12mn,
      drop = FALSE
    ) +
    scale_shape_manual(
      name = "External FDR",
      values = c(
        "FDR >= 0.05" = 21,
        "FDR < 0.05" = 23
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq(-0.5, 1.25, by = 0.25),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      breaks = seq(-0.25, 1.00, by = 0.25),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    coord_equal(
      xlim = c(-0.75, 1.42),
      ylim = c(-0.48, 1.12),
      clip = "off"
    ) +
    labs(
      title = paste0(
        "E3-residual genes retained discovery-defined\n",
        "directions in GSE65682"
      ),
      x = expression(
        paste("Discovery random-effects ", log[2], "FC")
      ),
      y = expression(
        paste("GSE65682 cohort-adjusted ", log[2], "FC")
      )
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        title.position = "left",
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          shape = 21,
          colour = "black",
          size = 2.6
        )
      ),
      shape = guide_legend(
        order = 2,
        title.position = "left",
        nrow = 1,
        override.aes = list(
          fill = "white",
          colour = "black",
          size = 2.6
        )
      )
    ) +
    theme_panel_step12mn() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 9.6,
        face = "bold",
        lineheight = 1.03
      ),
      legend.title = element_text(
        size = 7.0,
        face = "bold"
      ),
      legend.text = element_text(size = 6.7),
      legend.key.width = grid::unit(0.40, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.margin = margin(7, 9, 5, 7)
    )
  
  write.csv(
    residual_gene_step12m,
    paste0(
      scatter_base_step12m,
      "_source_data.csv"
    ),
    row.names = FALSE
  )
  
  residual_summary_step12m <- data.frame(
    Metric = c(
      "E3-residual signature genes",
      "Externally evaluable residual genes",
      "Directionally concordant residual genes",
      "Pearson correlation",
      "Pearson correlation P value",
      "Spearman correlation",
      "Spearman correlation P value",
      "External nominal P below 0.05",
      "External FDR below 0.05",
      "External absolute log2FC above 0.5"
    ),
    Value = c(
      "17",
      "15",
      "15/15",
      "0.732",
      "0.00194",
      "0.729",
      "0.00207",
      "6/15",
      "2/15",
      "2/15"
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    residual_summary_step12m,
    paste0(
      scatter_base_step12m,
      "_summary.csv"
    ),
    row.names = FALSE
  )
  
  save_panel_step12mn(
    plot_object = scatter_plot_step12m,
    file_base = scatter_base_step12m,
    width_mm = 135,
    height_mm = 96,
    dpi = 600
  )
  
  # ============================================================
  # Step 12N: Functional enrichment after overlap exclusion
  # ============================================================
  
  # ------------------------------------------------------------
  # 12N.1 Audited values from Supplementary Tables S4 and S9
  # ------------------------------------------------------------
  
  enrichment_step12n <- data.frame(
    ID = c(
      "GO:0000280",
      "GO:0048285",
      "GO:0140014",
      "hsa05322",
      "hsa04613"
    ),
    Term = c(
      "Nuclear division",
      "Organelle fission",
      "Mitotic nuclear division",
      "Systemic lupus erythematosus\u2020",
      "Neutrophil extracellular trap formation"
    ),
    Database = c(
      "GO Biological Process",
      "GO Biological Process",
      "GO Biological Process",
      "KEGG",
      "KEGG"
    ),
    Full_FDR = c(
      3.79831060924414e-4,
      4.02505520531409e-4,
      1.11914052981892e-2,
      1.85303624768699e-4,
      1.85303624768699e-4
    ),
    Residual_FDR = c(
      1.44062969065282e-2,
      1.44062969065282e-2,
      1.89381410431842e-2,
      4.70865775199584e-2,
      6.68379365333511e-2
    ),
    Full_gene_count = c(
      8L,
      8L,
      5L,
      4L,
      5L
    ),
    Residual_gene_count = c(
      5L,
      5L,
      4L,
      2L,
      2L
    ),
    stringsAsFactors = FALSE
  )
  
  enrichment_step12n$Full_minus_log10_FDR <-
    -log10(enrichment_step12n$Full_FDR)
  
  enrichment_step12n$Residual_minus_log10_FDR <-
    -log10(enrichment_step12n$Residual_FDR)
  
  enrichment_step12n$Residual_status <- ifelse(
    enrichment_step12n$Residual_FDR < 0.05,
    "Retained",
    "Not retained"
  )
  
  enrichment_step12n$gene_count_change <- paste0(
    enrichment_step12n$Full_gene_count,
    " -> ",
    enrichment_step12n$Residual_gene_count
  )
  
  term_order_step12n <- enrichment_step12n$Term
  
  enrichment_step12n$Term <- factor(
    enrichment_step12n$Term,
    levels = rev(term_order_step12n)
  )
  
  # ------------------------------------------------------------
  # 12N.2 Numerical assertions
  # ------------------------------------------------------------
  
  stopifnot(
    nrow(enrichment_step12n) == 5L,
    sum(enrichment_step12n$Residual_FDR < 0.05) == 4L,
    enrichment_step12n$Residual_FDR[
      enrichment_step12n$ID == "hsa04613"
    ] > 0.05,
    enrichment_step12n$Residual_FDR[
      enrichment_step12n$ID == "hsa05322"
    ] < 0.05,
    identical(
      enrichment_step12n$Full_gene_count,
      c(8L, 8L, 5L, 4L, 5L)
    ),
    identical(
      enrichment_step12n$Residual_gene_count,
      c(5L, 5L, 4L, 2L, 2L)
    )
  )
  
  fdr_threshold_step12n <- -log10(0.05)
  # ------------------------------------------------------------
  # 12N.3 Dumbbell plot: full-width layout
  # ------------------------------------------------------------
  
  term_display_step12n <- c(
    "GO:0000280" = "Nuclear division",
    "GO:0048285" = "Organelle fission",
    "GO:0140014" = "Mitotic nuclear division",
    "hsa05322" = "Systemic lupus erythematosus\u2020",
    "hsa04613" = "Neutrophil extracellular\ntrap formation"
  )
  
  enrichment_step12n$Term_display <- unname(
    term_display_step12n[enrichment_step12n$ID]
  )
  
  stopifnot(
    !anyNA(enrichment_step12n$Term_display)
  )
  
  enrichment_step12n$Term_display <- factor(
    enrichment_step12n$Term_display,
    levels = rev(unname(term_display_step12n))
  )
  
  enrichment_plot_step12n <- ggplot(
    enrichment_step12n,
    aes(y = Term_display)
  ) +
    geom_vline(
      xintercept = fdr_threshold_step12n,
      colour = "#B2473E",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_segment(
      aes(
        x = Full_minus_log10_FDR,
        xend = Residual_minus_log10_FDR,
        yend = Term_display
      ),
      colour = "#C7C7C7",
      linewidth = 0.85
    ) +
    geom_point(
      aes(
        x = Full_minus_log10_FDR,
        colour = "Full NS-higher signature"
      ),
      shape = 16,
      size = 3.0
    ) +
    geom_point(
      aes(
        x = Residual_minus_log10_FDR,
        colour = "E3-residual NS-higher subset"
      ),
      shape = 16,
      size = 3.0
    ) +
    geom_text(
      aes(
        x = 4.17,
        label = gene_count_change
      ),
      hjust = 0,
      size = 2.55,
      family = base_family_step12mn,
      colour = "#303030"
    ) +
    annotate(
      geom = "text",
      x = fdr_threshold_step12n,
      y = 5.62,
      label = "FDR = 0.05",
      hjust = 0.5,
      vjust = 0,
      size = 2.40,
      family = base_family_step12mn,
      colour = "#B2473E"
    ) +
    annotate(
      geom = "text",
      x = 4.17,
      y = 5.62,
      label = "Overlap genes\nfull -> residual",
      hjust = 0,
      vjust = 0,
      size = 2.30,
      lineheight = 0.92,
      family = base_family_step12mn,
      colour = "#303030"
    ) +
    scale_colour_manual(
      name = NULL,
      values = c(
        "Full NS-higher signature" = "#A7A7A7",
        "E3-residual NS-higher subset" = "#365F8D"
      ),
      breaks = c(
        "Full NS-higher signature",
        "E3-residual NS-higher subset"
      )
    ) +
    scale_x_continuous(
      breaks = 0:4,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(
      expand = expansion(add = c(0.45, 1.10))
    ) +
    coord_cartesian(
      xlim = c(0, 5.05),
      clip = "off"
    ) +
    labs(
      title = paste0(
        "Overlap exclusion retained cell-cycle enrichment\n",
        "while attenuating NET enrichment"
      ),
      x = expression(-log[10](FDR)),
      y = NULL,
      caption = "\u2020 KEGG database annotation term."
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2.8)
      )
    ) +
    theme_panel_step12mn() +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(
        size = 9.5,
        face = "bold",
        lineheight = 1.03
      ),
      panel.grid.major.x = element_line(
        colour = "#E4E4E4",
        linewidth = 0.30
      ),
      panel.grid.minor = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(
        size = 7.5,
        lineheight = 0.92,
        hjust = 1
      ),
      legend.position = "bottom",
      legend.text = element_text(size = 7.0),
      legend.key.width = grid::unit(0.42, "cm"),
      legend.spacing.x = grid::unit(0.15, "cm"),
      plot.caption = element_text(
        size = 6.5,
        colour = "#4A4A4A",
        hjust = 0
      ),
      plot.margin = margin(7, 9, 5, 7)
    )
  
  dumbbell_base_step12n <- file.path(
    dumbbell_dir_step12n,
    "Figure_E3_residual_functional_dumbbell"
  )
  
  write.csv(
    enrichment_step12n,
    paste0(
      dumbbell_base_step12n,
      "_source_data.csv"
    ),
    row.names = FALSE
  )
  
  save_panel_step12mn(
    plot_object = enrichment_plot_step12n,
    file_base = dumbbell_base_step12n,
    width_mm = 170,
    height_mm = 92,
    dpi = 600
  )
  
  # ------------------------------------------------------------
  # 12M-12N.4 Final output audit
  # ------------------------------------------------------------
  
  expected_outputs_step12mn <- c(
    paste0(scatter_base_step12m, ".pdf"),
    paste0(scatter_base_step12m, ".svg"),
    paste0(scatter_base_step12m, ".png"),
    paste0(scatter_base_step12m, ".tiff"),
    paste0(
      scatter_base_step12m,
      "_source_data.csv"
    ),
    paste0(
      scatter_base_step12m,
      "_summary.csv"
    ),
    paste0(dumbbell_base_step12n, ".pdf"),
    paste0(dumbbell_base_step12n, ".svg"),
    paste0(dumbbell_base_step12n, ".png"),
    paste0(dumbbell_base_step12n, ".tiff"),
    paste0(
      dumbbell_base_step12n,
      "_source_data.csv"
    )
  )
  
  output_audit_step12mn <- data.frame(
    file = expected_outputs_step12mn,
    exists = file.exists(expected_outputs_step12mn),
    size_bytes = file.info(
      expected_outputs_step12mn
    )$size,
    stringsAsFactors = FALSE
  )
  
  print(
    output_audit_step12mn,
    row.names = FALSE
  )
  
  stopifnot(
    all(output_audit_step12mn$exists),
    all(output_audit_step12mn$size_bytes > 0)
  )
  
  cat("\nStep 12M scatter directory:\n")
  cat(
    normalizePath(
      scatter_dir_step12m,
      mustWork = TRUE
    ),
    "\n"
  )
  
  cat("\nStep 12N dumbbell directory:\n")
  cat(
    normalizePath(
      dumbbell_dir_step12n,
      mustWork = TRUE
    ),
    "\n"
  )
  
  cat(
    "\nSTEPS 12M-12N completed successfully.\n"
  )
  
  cat(
    "Two independent panels were generated; no main figure was assembled.\n"
  )
  
  cat(
    "No online query was executed.\n"
  )
  
  cat(
    "====================================================\n"
  )