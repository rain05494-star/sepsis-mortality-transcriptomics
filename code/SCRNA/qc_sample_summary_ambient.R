
# ============================================================================
# SCRIPT -- Single-cell QC summary + ambient PROXY audit (Methods backing)
#   Produces the sample-layer QC summary and a DESCRIPTIVE ambient-RNA proxy
#   audit comparing granule / residual-signature marker detection between
#   neutrophil-axis, other granulocyte-like, and STRICT non-granulocyte cells.
#   It DOES NOT perform decontamination (no SoupX/CellBender/etc.).
#
#   CHANGELOG
#     v1.1: strict non-granulocyte compartments; residual-neutrophil panel
#           added (CEACAM8/OLFM4/RNASE3/TCN1).
#     v1.2: donor-level ambient correlation split by marker panel; pooled
#           proxy retained only for backward compatibility.
#     v1.3: donor-panel ambient table reports donor-level denominator cell
#           counts; duplicate ambient markers across panels are fail-closed.
#     v1.4: donor-state eligibility (MIN_CELLS_PER_DONOR_STATE) applied;
#           bootstrap CI + three-way verdict + BH-FDR + depth-adjusted
#           sensitivity; full-object invariants; residual panel renamed
#           residual_E3_granulocyte_panel; RUNNING/COMPLETE/FAILED status.
#     v1.5: main() wrapper with reliable on.exit status capture; fixed the
#           vapply integer/double mismatch; manifest check after sample_tab
#           and canonicalized; FORMAL_RUN gate; MIN_SHARED_STATES_PER_DONOR;
#           partial/UMI-fraction sensitivities gated on evaluability; BH-FDR
#           over formal panels only; sha256 with md5 fallback; script path
#           from commandArgs(); dead variables removed.
#     v1.6: removed undefined RESIDUAL_PANEL_SOURCE_NOTE (was the only
#           guaranteed runtime blocker); E3 registry verified by STRICT gene
#           set equality + duplicate check + non-empty values + single frozen
#           evidence tier + external EXPECTED_E3_REGISTRY_SHA256 (no in-file
#           self-reference); formal mode forbids md5 fallback and requires
#           real SHA-256 for RDS/script/manifest/registry; provenance records
#           frozen_manifest_hash + e3_registry_file/hash; RUN_PURPOSE is
#           mode-stamped (_TEST/_FORMAL) so outputs cannot be confused;
#           donor_cov records per-donor eligible_shared_state_names (state
#           sets may differ across donors); BH-FDR uses the pre-specified
#           four-panel family (n = 4), comment matches implementation.
#
#     v1.7: fixed two warning sources without changing any statistic:
#           (a) any(z) -> any(z > 0) in the three tapply calls of the
#               ambient_rows block (logical coercion warning);
#           (b) GetAssayData slot=counts -> SeuratObject>=5.0.0 version-gated
#               layer=counts branch (deprecation warning).
#   MAIN-pipeline context (pre-ACSL1 environment). ALL outputs under:
#     /home/sunshine/predicate/singlecell/outputs/GSE216009_32gene_analysis/
#       qc_sample_summary_ambient/<YYYYMMDD>_singlecell_qc_TEST[_N]
#       qc_sample_summary_ambient/<YYYYMMDD>_singlecell_qc_FORMAL[_N]
#   Nothing touches the ACSL1_verification/ folder.
#
#   Interpretive guardrails (v1.6):
#     - Bootstrap CIs are descriptive and NOT multiplicity-adjusted.
#     - BH-adjusted P-values test zero correlation only; they do NOT
#       adjudicate the rho >= 0.50 threshold verdict.
#     - Donor signal is state-standardized as the mean across THAT donor's
#       eligible shared strict non-gran states; eligible state sets may
#       differ across donors (names are recorded per donor).
#     - The script supports "no strong donor-level positive structure in this
#       proxy audit"; it CANNOT support "ambient RNA contamination has been
#       excluded".
#
#   Non-correction declaration: analyses use the supplied processed object;
#     NO additional ambient or doublet correction is applied. Doublet
#     sensitivity is handled separately as robustness (§3.9 / S18).
#
#   Cohort reconstruction is fail-closed (identical to Scripts 01-05):
#     48 samples / 39 donors / 26 acute / 9 paired / 151,837 primary / 24 states.
#
#   RUN: Rscript qc_sample_summary_ambient.R
# ============================================================================

options(stringsAsFactors = FALSE)

# ---- 0. CONFIG --------------------------------------------------------------
SCRIPT_VERSION <- "v1.7"
FORMAL_RUN     <- FALSE   # TRUE: hard-stop unless all frozen inputs configured

OUTPUT_ROOT <- "/home/sunshine/predicate/singlecell/outputs/GSE216009_32gene_analysis"
INPUT_RDS   <- file.path("/home/sunshine/predicate/singlecell",
                         "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
STAGE_NAME  <- "qc_sample_summary_ambient"
# v1.6: outputs are mode-stamped so test and formal runs cannot be confused.
RUN_PURPOSE <- if (FORMAL_RUN) "singlecell_qc_FORMAL" else "singlecell_qc_TEST"

COUNT_ASSAY      <- "RNA"
COUNT_LAYER      <- "counts"
FINE_STATE_FIELD <- "fine_annot"
SAMPLE_FIELD     <- "sample_id"
CONDITION_FIELD  <- "diagnosis"
ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")

MIN_CELLS_PER_DONOR_STATE <- 20L   # donor-state eligibility (Scripts 02/03)
MIN_SHARED_STATES_PER_DONOR <- 3L  # donor signal requires >= this many eligible shared states
MIN_EVALUABLE_DONORS      <- 20L   # minimum donors for a correlation verdict
STRONG_RHO_THRESHOLD      <- 0.50
BOOT_REPS                 <- 2000L
BOOT_SEED                 <- 20260813L
MIN_SHARED_STATE_DONOR_N  <- ceiling(26L / 2)   # shared-state standard = >=13/26 donors

EXPECTED_N_TOTAL_CELLS    <- 272993L
EXPECTED_N_SAMPLES        <- 48L
EXPECTED_N_DONORS         <- 39L
EXPECTED_N_ACUTE_DONORS   <- 26L
EXPECTED_N_HC_DONORS      <- 6L
EXPECTED_N_SURGERY_DONORS <- 7L
EXPECTED_N_CONV_SAMPLES   <- 9L
EXPECTED_N_PAIRED_DONORS  <- 9L
EXPECTED_N_PRIMARY_CELLS  <- 151837L
EXPECTED_N_FINE_STATES    <- 24L

# v1.4: PASTE the exact 24 fine-state names from the frozen Script-01 output
#   (sort(unique(fine_annot))). Required when FORMAL_RUN = TRUE.
FROZEN_ALL_STATES_24 <- NULL

# v1.4: path to the frozen 48-sample manifest CSV (columns sample_id, donor,
#   diagnosis). Required when FORMAL_RUN = TRUE.
FROZEN_MANIFEST_FILE <- NULL

# v1.5/v1.6: verifiable E3 provenance. Registry CSV requires columns
#   gene, panel, evidence_tier, source_file. The expected SHA-256 comes from
#   THIS independent config constant (NOT from a hash stored inside the
#   registry file, which would be self-referential). The amendment ID is YOUR
#   governance record value. All required in formal mode.
E3_RESIDUAL_REGISTRY_FILE <- NULL
EXPECTED_E3_REGISTRY_SHA256 <- NULL
E3_RESIDUAL_AMENDMENT_ID  <- NULL
E3_RESIDUAL_LOCK_DATE     <- "2026-08-12"

AXIS_STATES <- c(
  "Cycling_neutrophil_progenitors",
  "MPO+_immature_neutrophils_or_progenitors",
  "PADI4+_immature_neutrophils",
  "IL1R2+_immature_neutrophils",
  "S100A8-9_hi_neutrophils",
  "Mature_neutrophils",
  "Degranulating_neutrophils",
  "Apoptosing_neutrophils")

GRANULE_PANELS <- list(
  primary_azurophilic = c("MPO","ELANE","CTSG","DEFA4"),
  secondary_specific = c("LTF","LCN2"),
  tertiary_gelatinase = c("MMP8","MMP9"))

GRANULOCYTE_EXTRA_STATES <- c(
  "Eosinophils",
  "Mast_cells/eosiniophils"
)
GRANULOCYTE_LIKE_STATES <- unique(c(AXIS_STATES, GRANULOCYTE_EXTRA_STATES))

# residual E3 granulocyte panel = manuscript §3.9 residual ND four-gene subset
# (S7/S18, locked 2026-08-12). RNASE3 is eosinophil-associated, hence "E3
# granulocyte", not "neutrophil". If a gene is missing from the Rhapsody panel,
# the panel emits not_evaluable_incomplete_panel (no directional verdict).
AMBIENT_PANELS <- c(
  GRANULE_PANELS,
  list(
    residual_E3_granulocyte_panel = c("CEACAM8", "OLFM4", "RNASE3", "TCN1")
  )
)
REQUIRED_AMBIENT_PANELS <- c("primary_azurophilic", "secondary_specific",
                             "tertiary_gelatinase")

# ---- 1. PACKAGES + HELPERS --------------------------------------------------
needed <- c("SeuratObject","Matrix")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
HAS_GGPLOT2 <- requireNamespace("ggplot2", quietly = TRUE)
HAS_DIGEST  <- requireNamespace("digest", quietly = TRUE)

stop_msg <- function(...) stop(paste0(...), call. = FALSE)
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}
near_constant <- function(x, tol = 1e-12) {
  x <- x[is.finite(x)]
  length(x) < 3L || length(unique(round(x, 12))) < 2L || diff(range(x)) <= tol
}

write_csv_atomic <- function(x, path) {
  stopifnot(is.data.frame(x))
  tf <- tempfile(tmpdir = dirname(path), fileext = ".csv")
  on.exit(if (file.exists(tf)) unlink(tf), add = TRUE)
  con <- file(tf, open = "w", encoding = "UTF-8")
  con_closed <- FALSE
  on.exit(try(if (!con_closed && isOpen(con)) close(con), silent = TRUE), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  flush(con)
  close(con)
  con_closed <- TRUE
  ok <- file.rename(tf, path)
  if (!ok) stop_msg("Failed to atomically write: ", path)
  invisible(path)
}

make_run_dir <- function(root, stage, purpose) {
  base <- file.path(root, stage)
  stamp <- format(Sys.Date(), "%Y%m%d")
  d0 <- file.path(base, paste0(stamp, "_", purpose))
  if (!dir.exists(d0)) { dir.create(d0, recursive = TRUE); return(d0) }
  i <- 2L
  repeat {
    d <- file.path(base, paste0(stamp, "_", purpose, "_", i))
    if (!dir.exists(d)) { dir.create(d, recursive = TRUE); return(d) }
    i <- i + 1L
  }
}

# v1.5: resolve the running script's own path (works under
#   Rscript /full/path/script.R; md5 otherwise resolves to NA).
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  w <- grep("^--file=", args, value = TRUE)
  if (length(w)) return(tryCatch(normalizePath(sub("^--file=", "", w[1])),
                                 error = function(e) sub("^--file=", "", w[1])))
  NA_character_
}

file_hash <- function(path, algo = "sha256") {
  if (!file.exists(path)) return(NA_character_)
  if (algo == "sha256" && HAS_DIGEST) {
    tryCatch(digest::digest(path, algo = "sha256", file = TRUE),
             error = function(e) NA_character_)
  } else {
    tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  }
}

canonicalize_manifest <- function(x) {
  x <- data.frame(
    sample_id = as.character(x$sample_id),
    donor     = as.character(x$donor),
    diagnosis = as.character(x$diagnosis),
    stringsAsFactors = FALSE)
  x <- unique(x)
  x <- x[order(x$sample_id, x$donor, x$diagnosis), , drop = FALSE]
  rownames(x) <- NULL
  x
}

# ---- 2. MAIN -----------------------------------------------------------------
main <- function() {
  
  # ---- 2.1 OUTPUT ROOT + RUN STATUS ----------------------------------------
  OUT_ROOT  <- make_run_dir(OUTPUT_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  PROV_DIR  <- file.path(OUT_ROOT, "provenance")
  for (d in c(TABLE_DIR, FIG_DIR, PROV_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  STATUS_FILE <- file.path(OUT_ROOT, "STATUS.txt")
  writeLines("RUNNING", STATUS_FILE)
  run_status <- "RUNNING"
  on.exit(if (identical(run_status, "RUNNING")) {
    try(writeLines("FAILED", STATUS_FILE), silent = TRUE)
  }, add = TRUE)
  cat("[setup] QC output root (", RUN_PURPOSE, "):\n  ", OUT_ROOT, "\n", sep = "")
  
  # ---- 2.2 FORMAL-RUN GATE ---------------------------------------------------
  if (FORMAL_RUN) {
    if (is.null(FROZEN_ALL_STATES_24)) stop_msg("Formal run requires FROZEN_ALL_STATES_24.")
    if (is.null(FROZEN_MANIFEST_FILE)) stop_msg("Formal run requires FROZEN_MANIFEST_FILE.")
    if (is.null(E3_RESIDUAL_REGISTRY_FILE)) stop_msg("Formal run requires E3_RESIDUAL_REGISTRY_FILE.")
    if (is.null(EXPECTED_E3_REGISTRY_SHA256)) stop_msg("Formal run requires EXPECTED_E3_REGISTRY_SHA256.")
    if (is.null(E3_RESIDUAL_AMENDMENT_ID)) stop_msg("Formal run requires E3_RESIDUAL_AMENDMENT_ID.")
    if (!HAS_DIGEST) stop_msg("Formal run requires package 'digest' for SHA-256 hashing.")
  }
  
  # ---- 2.3 READ RDS + COUNT LAYER + OBJECT INVARIANTS ------------------------
  if (!file.exists(INPUT_RDS)) stop_msg("Input RDS not found: ", INPUT_RDS)
  input_rds_hash <- file_hash(INPUT_RDS)
  hash_algo_used <- if (HAS_DIGEST) "sha256" else "md5_fallback"
  obj <- readRDS(INPUT_RDS)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  if (ncol(obj) != EXPECTED_N_TOTAL_CELLS)
    stop_msg("Expected ", EXPECTED_N_TOTAL_CELLS, " total cells; observed ", ncol(obj))
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay = COUNT_ASSAY, layer = COUNT_LAYER, fast = FALSE)
  } else {
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, layer = COUNT_LAYER)
    } else {
      cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, slot = COUNT_LAYER)
    }
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  if (!identical(colnames(cnt), colnames(obj))) stop_msg("Count matrix colnames != object cells")
  
  if (!all(is.finite(cnt@x))) stop_msg("Count matrix contains non-finite values.")
  if (any(cnt@x < 0)) stop_msg("Count matrix contains negative values.")
  if (max(abs(cnt@x - round(cnt@x))) > 1e-9)
    stop_msg("Count matrix contains non-integer values.")
  
  lib_all <- as.numeric(Matrix::colSums(cnt))
  if (any(!is.finite(lib_all) | lib_all <= 0)) stop_msg("Non-positive cell library sizes found.")
  
  for (f in c(SAMPLE_FIELD, CONDITION_FIELD, FINE_STATE_FIELD))
    if (!(f %in% colnames(meta))) stop_msg("meta.data missing required column: ", f)
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  
  obs_states <- sort(unique(cell_state))
  if (length(obs_states) != EXPECTED_N_FINE_STATES)
    stop_msg("Expected ", EXPECTED_N_FINE_STATES, " fine states; observed ",
             length(obs_states), ". fine_annot may have changed.")
  axis_missing <- setdiff(AXIS_STATES, obs_states)
  if (length(axis_missing))
    stop_msg("Maturation-axis states not observed: ",
             paste(axis_missing, collapse = ", "))
  cat("[object] cells =", ncol(obj), "| samples =", length(unique(cell_sample)),
      "| fine states =", length(obs_states), "\n")
  
  # ---- 2.4 COMPARTMENT DEFINITION ---------------------------------------------
  gran_extra_missing <- setdiff(GRANULOCYTE_EXTRA_STATES, obs_states)
  if (length(gran_extra_missing))
    stop_msg("Granulocyte-like extra states absent from this object: ",
             paste(gran_extra_missing, collapse = ", "),
             ". They define the strict non-granulocyte boundary; do not proceed.")
  GRANULOCYTE_LIKE_STATES <- intersect(GRANULOCYTE_LIKE_STATES, obs_states)
  STRICT_NON_GRAN_STATES <- setdiff(obs_states, GRANULOCYTE_LIKE_STATES)
  if (!length(STRICT_NON_GRAN_STATES))
    stop_msg("No strict non-granulocyte states left after excluding granulocyte-like states.")
  cat("[compartment] axis =", length(AXIS_STATES),
      "| other granulocyte-like =", length(GRANULOCYTE_EXTRA_STATES),
      "| strict non-granulocyte =", length(STRICT_NON_GRAN_STATES), "\n")
  
  # ---- 2.5 FROZEN STATE-SET CHECK --------------------------------------------
  frozen_state_check <- "not_configured"
  if (!is.null(FROZEN_ALL_STATES_24)) {
    stopifnot(length(FROZEN_ALL_STATES_24) == EXPECTED_N_FINE_STATES)
    if (!setequal(obs_states, FROZEN_ALL_STATES_24))
      stop_msg("fine_annot state set does not match the frozen 24-state set.")
    frozen_state_check <- "verified"
  }
  
  # ---- 2.6 COHORT REBUILD (fail-closed; Scripts 02/03) ------------------------
  sample_tab <- unique(data.frame(sample_id = cell_sample, diagnosis = cell_diag,
                                  stringsAsFactors = FALSE))
  n_diag_per_sample <- tapply(sample_tab$diagnosis, sample_tab$sample_id,
                              function(z) length(unique(z)))
  if (any(n_diag_per_sample != 1L)) stop_msg("Samples with multiple diagnosis values.")
  sample_tab$donor  <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort <- NA_character_
  sample_tab$cohort[sample_tab$diagnosis %in% ACUTE_DIAG] <- "Acute_sepsis"
  sample_tab$cohort[sample_tab$diagnosis == "Conv"]        <- "Convalescent"
  sample_tab$cohort[sample_tab$diagnosis == "HV"]          <- "Healthy_control"
  sample_tab$cohort[sample_tab$diagnosis == "CS"]          <- "Surgery_control"
  if (anyNA(sample_tab$cohort))
    stop_msg("Unmapped diagnosis: ",
             paste(unique(sample_tab$diagnosis[is.na(sample_tab$cohort)]), collapse = ", "))
  stopifnot(identical(grepl("_CONV$", sample_tab$sample_id), sample_tab$diagnosis == "Conv"))
  
  n_samples <- length(unique(sample_tab$sample_id))
  n_donors  <- length(unique(sample_tab$donor))
  acute_donors   <- unique(sample_tab$donor[sample_tab$cohort == "Acute_sepsis"])
  hc_donors      <- unique(sample_tab$donor[sample_tab$cohort == "Healthy_control"])
  surgery_donors <- unique(sample_tab$donor[sample_tab$cohort == "Surgery_control"])
  conv_donors    <- unique(sample_tab$donor[sample_tab$cohort == "Convalescent"])
  paired_donors  <- intersect(acute_donors, conv_donors)
  stopifnot(n_samples == EXPECTED_N_SAMPLES, n_donors == EXPECTED_N_DONORS)
  stopifnot(length(acute_donors)   == EXPECTED_N_ACUTE_DONORS,
            length(hc_donors)      == EXPECTED_N_HC_DONORS,
            length(surgery_donors) == EXPECTED_N_SURGERY_DONORS,
            sum(sample_tab$cohort == "Convalescent") == EXPECTED_N_CONV_SAMPLES,
            length(paired_donors)  == EXPECTED_N_PAIRED_DONORS,
            all(conv_donors %in% acute_donors))
  
  # ---- 2.7 FROZEN MANIFEST CHECK (after sample_tab exists; v1.5) --------------
  frozen_manifest_check <- "not_configured"
  if (!is.null(FROZEN_MANIFEST_FILE)) {
    if (!file.exists(FROZEN_MANIFEST_FILE))
      stop_msg("FROZEN_MANIFEST_FILE not found: ", FROZEN_MANIFEST_FILE)
    fm <- read.csv(FROZEN_MANIFEST_FILE, stringsAsFactors = FALSE)
    req_cols <- c("sample_id", "donor", "diagnosis")
    if (!all(req_cols %in% colnames(fm)))
      stop_msg("Frozen manifest missing columns: ",
               paste(setdiff(req_cols, colnames(fm)), collapse = ", "))
    frozen_triples <- unique(data.frame(
      sample_id = fm$sample_id, donor = fm$donor, diagnosis = fm$diagnosis,
      stringsAsFactors = FALSE))
    obs_triples <- unique(data.frame(
      sample_id = sample_tab$sample_id, donor = sample_tab$donor,
      diagnosis = sample_tab$diagnosis, stringsAsFactors = FALSE))
    if (!identical(canonicalize_manifest(frozen_triples),
                   canonicalize_manifest(obs_triples)))
      stop_msg("Object sample_id/donor/diagnosis manifest does not match the frozen manifest.")
    frozen_manifest_check <- "verified"
  }
  cat("[frozen] state set:", frozen_state_check, "| manifest:", frozen_manifest_check, "\n")
  
  m <- match(cell_sample, sample_tab$sample_id)
  if (anyNA(m)) stop_msg("Sample->donor resolution failed.")
  donor_all  <- sample_tab$donor[m]
  cohort_all <- sample_tab$cohort[m]
  is_primary <- cohort_all == "Acute_sepsis"
  stopifnot(all(!is_primary[grepl("_CONV$", cell_sample)]))
  primary_idx <- which(is_primary)
  if (length(primary_idx) != EXPECTED_N_PRIMARY_CELLS)
    stop_msg("Expected ", EXPECTED_N_PRIMARY_CELLS, " primary cells; observed ",
             length(primary_idx))
  cat("[cohort] samples =", n_samples, "| donors =", n_donors,
      "| acute =", length(acute_donors), "| paired =", length(paired_donors),
      "| primary cells =", length(primary_idx), "\n")
  
  # ---- 2.8 E3 REGISTRY VERIFICATION (v1.6: strict) -----------------------------
  e3_registry_check <- "not_configured"
  e3_evidence_tier <- NA_character_
  e3_source_files  <- NA_character_
  if (!is.null(E3_RESIDUAL_REGISTRY_FILE)) {
    if (!file.exists(E3_RESIDUAL_REGISTRY_FILE))
      stop_msg("E3 registry file not found: ", E3_RESIDUAL_REGISTRY_FILE)
    reg <- read.csv(E3_RESIDUAL_REGISTRY_FILE, stringsAsFactors = FALSE)
    need_cols <- c("gene", "panel", "evidence_tier", "source_file")
    if (!all(need_cols %in% colnames(reg)))
      stop_msg("E3 registry missing columns: ",
               paste(setdiff(need_cols, colnames(reg)), collapse = ", "))
    reg_res <- reg[as.character(reg$panel) == "residual_E3_granulocyte_panel", , drop = FALSE]
    
    # 1) Strict set equality: no missing genes, no extra genes.
    expected_genes <- sort(AMBIENT_PANELS$residual_E3_granulocyte_panel)
    observed_genes <- sort(unique(trimws(as.character(reg_res$gene))))
    if (!identical(expected_genes, observed_genes))
      stop_msg("E3 registry gene set mismatch. Missing: ",
               paste(setdiff(expected_genes, observed_genes), collapse = ", "),
               "; extra: ", paste(setdiff(observed_genes, expected_genes), collapse = ", "))
    
    # 2) No duplicate genes / conflicting records.
    if (anyDuplicated(reg_res$gene))
      stop_msg("Duplicated genes in the residual E3 registry.")
    
    # 3) No missing/empty values in any required column.
    for (f in c("gene", "panel", "evidence_tier", "source_file")) {
      z <- trimws(as.character(reg_res[[f]]))
      if (anyNA(z) || any(!nzchar(z)))
        stop_msg("E3 registry contains missing/empty values in: ", f)
    }
    
    # 4) evidence_tier must be a single frozen value (not a semicolon-joined mix).
    if (length(unique(trimws(as.character(reg_res$evidence_tier)))) != 1L)
      stop_msg("E3 registry evidence_tier is not a single frozen value.")
    
    e3_registry_check <- "verified"
    e3_evidence_tier <- unique(trimws(as.character(reg_res$evidence_tier)))[1]
    e3_source_files  <- paste(unique(trimws(as.character(reg_res$source_file))),
                              collapse = ";")
    
    # 5) SHA-256 against the independent config constant (no in-file self-reference).
    if (!is.null(EXPECTED_E3_REGISTRY_SHA256)) {
      obs_sha <- file_hash(E3_RESIDUAL_REGISTRY_FILE, algo = "sha256")
      if (is.na(obs_sha) || !identical(tolower(obs_sha), tolower(EXPECTED_E3_REGISTRY_SHA256)))
        stop_msg("E3 registry SHA-256 mismatch.")
    }
    cat("[e3] registry:", e3_registry_check,
        "| amendment:", E3_RESIDUAL_AMENDMENT_ID, "\n")
  }
  
  # ---- 2.9 PER-CELL QC METRICS + EXACT CHECKS --------------------------------
  for (qc in c("nCount_RNA","nFeature_RNA"))
    if (!(qc %in% colnames(meta))) stop_msg("meta.data missing QC column: ", qc)
  umi_cell  <- as.numeric(meta[["nCount_RNA"]])
  gene_cell <- as.numeric(meta[["nFeature_RNA"]])
  if (any(!is.finite(umi_cell) | umi_cell < 0))
    stop_msg("nCount_RNA has non-finite or negative values.")
  if (any(!is.finite(gene_cell) | gene_cell < 0))
    stop_msg("nFeature_RNA has non-finite or negative values.")
  
  lib_chk <- abs(lib_all - umi_cell)
  if (max(lib_chk) != 0)
    stop_msg("nCount_RNA != colSums(counts) for some cells (max |diff| = ",
             max(lib_chk), "); metadata/counts layer mismatch.")
  gene_chk <- as.numeric(Matrix::colSums(cnt > 0))
  if (max(abs(gene_chk - gene_cell)) != 0)
    stop_msg("nFeature_RNA != colSums(counts > 0) for some cells (max |diff| = ",
             max(abs(gene_chk - gene_cell)), ").")
  cat("[qc] nCount_RNA == colSums(counts): OK | nFeature_RNA == colSums(counts>0): OK\n")
  
  # ---- 2.10 PART A -- SAMPLE-LAYER QC SUMMARY --------------------------------
  by_sample <- split(seq_len(ncol(obj)), cell_sample)
  sample_sum <- do.call(rbind, lapply(names(by_sample), function(s) {
    idx <- by_sample[[s]]
    data.frame(
      sample_id      = s,
      donor          = sample_tab$donor[match(s, sample_tab$sample_id)],
      cohort         = sample_tab$cohort[match(s, sample_tab$sample_id)],
      n_cells        = length(idx),
      median_UMI     = safe_median(umi_cell[idx]),
      median_genes   = safe_median(gene_cell[idx]),
      total_UMI      = sum(umi_cell[idx]),
      stringsAsFactors = FALSE)
  }))
  sample_sum <- sample_sum[order(sample_sum$cohort, -sample_sum$n_cells), ]
  write_csv_atomic(sample_sum, file.path(TABLE_DIR, "sample_level_qc_summary.csv"))
  
  donor_sum <- do.call(rbind, lapply(sort(unique(donor_all)), function(d) {
    idx <- which(donor_all == d)
    data.frame(
      donor         = d,
      cohort        = paste(unique(sample_tab$cohort[sample_tab$donor == d]), collapse = "&"),
      n_samples     = length(unique(cell_sample[idx])),
      n_cells       = length(idx),
      median_UMI    = safe_median(umi_cell[idx]),
      median_genes  = safe_median(gene_cell[idx]),
      stringsAsFactors = FALSE)
  }))
  write_csv_atomic(donor_sum, file.path(TABLE_DIR, "donor_level_qc_summary.csv"))
  
  primary_umi  <- umi_cell[primary_idx]
  primary_gene <- gene_cell[primary_idx]
  object_summary <- data.frame(
    metric = c("total_cells","total_samples","total_donors","fine_states",
               "acute_donors","paired_donors","primary_cells",
               "median_UMI_per_cell_all","median_genes_per_cell_all",
               "median_UMI_per_cell_primary","median_genes_per_cell_primary",
               "median_cells_per_sample","min_cells_per_sample","max_cells_per_sample"),
    value = c(ncol(obj), n_samples, n_donors, length(obs_states),
              length(acute_donors), length(paired_donors), length(primary_idx),
              safe_median(umi_cell), safe_median(gene_cell),
              safe_median(primary_umi), safe_median(primary_gene),
              safe_median(sample_sum$n_cells), min(sample_sum$n_cells), max(sample_sum$n_cells)),
    stringsAsFactors = FALSE)
  write_csv_atomic(object_summary, file.path(TABLE_DIR, "object_qc_summary.csv"))
  
  cat("\n[Part A] object-level QC\n")
  print(object_summary, row.names = FALSE)
  cat("\n[Part A] per-cohort sample cell counts (median, min-max)\n")
  cohort_n <- lapply(split(sample_sum$n_cells, sample_sum$cohort), function(z)
    c(median = safe_median(z), min = min(z), max = max(z)))
  print(do.call(rbind, cohort_n))
  
  # ---- 2.11 PART B -- AMBIENT PROXY AUDIT ------------------------------------
  panel_of <- setNames(
    rep(names(AMBIENT_PANELS), lengths(AMBIENT_PANELS)),
    unlist(AMBIENT_PANELS, use.names = FALSE))
  dup_panel_marker <- unique(names(panel_of)[duplicated(names(panel_of))])
  if (length(dup_panel_marker))
    stop_msg("Ambient markers assigned to multiple panels: ",
             paste(dup_panel_marker, collapse = ", "),
             ". Keep panels mutually exclusive or explicitly define how duplicated markers should be handled.")
  markers <- names(panel_of)
  feat <- rownames(cnt)
  gene_map <- data.frame(gene = markers, dataset_feature = markers,
                         panel = unname(panel_of[markers]),
                         present = markers %in% feat,
                         stringsAsFactors = FALSE)
  for (i in which(!gene_map$present)) {
    g <- gene_map$gene[i]
    ci <- feat[tolower(feat) == tolower(g)]
    if (length(ci) == 1L) { gene_map$dataset_feature[i] <- ci; gene_map$present[i] <- TRUE }
  }
  absent_g <- gene_map$gene[!gene_map$present]
  absent_gran <- absent_g[panel_of[absent_g] %in% REQUIRED_AMBIENT_PANELS]
  if (length(absent_gran))
    stop_msg("Required ambient markers absent from counts: ",
             paste(absent_gran, collapse = ", "))
  gene_map <- gene_map[gene_map$present, , drop = FALSE]
  markers <- gene_map$gene
  panel_complete_flag <- vapply(names(AMBIENT_PANELS), function(pn) {
    all(AMBIENT_PANELS[[pn]] %in% markers)
  }, logical(1))
  incomplete_panels <- names(AMBIENT_PANELS)[!panel_complete_flag]
  if (length(incomplete_panels))
    cat("[warn] incomplete panels (not_evaluable_incomplete_panel): ",
        paste(incomplete_panels, collapse = ", "), "\n", sep = "")
  gi <- match(gene_map$dataset_feature, feat)
  stopifnot(!anyNA(gi), !anyDuplicated(gi))
  Xm <- cnt[gi, primary_idx, drop = FALSE]
  rownames(Xm) <- markers
  
  state_p     <- cell_state[primary_idx]
  axis_idx_p  <- which(state_p %in% AXIS_STATES)
  gran_like_idx_p <- which(state_p %in% GRANULOCYTE_LIKE_STATES)
  non_gran_idx_p  <- which(state_p %in% STRICT_NON_GRAN_STATES)
  other_gl_idx_p  <- setdiff(gran_like_idx_p, axis_idx_p)   # Eos + Mast-eos only
  
  donor_p <- donor_all[primary_idx]
  ng_donor <- donor_p[non_gran_idx_p]
  ax_donor <- donor_p[axis_idx_p]
  gran_like_donor <- donor_p[gran_like_idx_p]
  other_gl_donor  <- donor_p[other_gl_idx_p]
  
  ambient_rows <- lapply(markers, function(g) {
    det <- as.numeric(Xm[g, ] > 0)
    ng_frac <- tapply(det[non_gran_idx_p], ng_donor, mean)
    ng_any  <- tapply(det[non_gran_idx_p], ng_donor, function(z) any(z > 0))
    ax_frac <- tapply(det[axis_idx_p], ax_donor, mean)
    ax_any  <- tapply(det[axis_idx_p], ax_donor, function(z) any(z > 0))
    og_frac <- tapply(det[other_gl_idx_p], other_gl_donor, mean)
    og_any  <- tapply(det[other_gl_idx_p], other_gl_donor, function(z) any(z > 0))
    data.frame(
      marker                            = g,
      panel                             = unname(panel_of[g]),
      axis_pooled_detection             = safe_mean(det[axis_idx_p]),
      axis_det_fraction_median          = safe_median(ax_frac),
      axis_donor_detection_rate         = safe_mean(ax_any),
      other_granulocyte_like_pooled_detection = safe_mean(det[other_gl_idx_p]),
      other_granulocyte_like_donor_detection_rate = safe_mean(og_any),
      strict_non_gran_pooled_detection  = safe_mean(det[non_gran_idx_p]),
      strict_non_gran_det_fraction_median = safe_median(ng_frac),
      strict_non_gran_donor_detection_rate = safe_mean(ng_any),
      stringsAsFactors = FALSE)
  })
  ambient_report <- do.call(rbind, ambient_rows)
  write_csv_atomic(ambient_report, file.path(TABLE_DIR, "ambient_granule_signal.csv"))
  
  ng_states <- STRICT_NON_GRAN_STATES
  presence_mat <- t(vapply(markers, function(g) {
    det <- as.numeric(Xm[g, ] > 0)
    vapply(ng_states, function(s) {
      sidx <- which(state_p == s)
      if (length(sidx)) safe_mean(det[sidx]) else NA_real_
    }, numeric(1))
  }, numeric(length(ng_states))))
  presence_long <- data.frame(
    marker            = rep(markers, each = length(ng_states)),
    panel             = rep(unname(panel_of[markers]), each = length(ng_states)),
    strict_non_gran_state = rep(ng_states, times = length(markers)),
    pooled_detection  = as.vector(t(presence_mat)),
    stringsAsFactors = FALSE)
  write_csv_atomic(presence_long, file.path(TABLE_DIR, "ambient_marker_presence_audit.csv"))
  
  # ---- 2.12 DONOR-STATE ELIGIBLE CELLS + DONOR COVARIATES --------------------
  ng_state_cells <- lapply(acute_donors, function(d) {
    cell_idx <- which(donor_p == d)
    lapply(STRICT_NON_GRAN_STATES, function(s) {
      idx <- cell_idx[state_p[cell_idx] == s]
      if (length(idx) >= MIN_CELLS_PER_DONOR_STATE) idx else integer(0)
    })
  })
  names(ng_state_cells) <- acute_donors
  
  state_eligible_donor_n <- vapply(STRICT_NON_GRAN_STATES, function(s) {
    j <- match(s, STRICT_NON_GRAN_STATES)
    sum(vapply(acute_donors, function(d)
      length(ng_state_cells[[d]][[j]]) >= MIN_CELLS_PER_DONOR_STATE, logical(1)))
  }, integer(1))
  shared_states <- STRICT_NON_GRAN_STATES[state_eligible_donor_n >= MIN_SHARED_STATE_DONOR_N]
  
  # v1.6: per-donor eligible shared-state NAMES, not just the count, so readers
  # can judge whether donors use comparable background state compositions.
  donor_cov <- data.frame(
    donor = acute_donors,
    primary_cell_n = vapply(acute_donors, function(d) sum(donor_p == d), integer(1)),
    strict_non_gran_cell_n = vapply(acute_donors, function(d)
      length(intersect(which(donor_p == d), non_gran_idx_p)), integer(1)),
    n_eligible_shared_states = vapply(acute_donors, function(d) {
      sum(vapply(shared_states, function(s) {
        j <- match(s, STRICT_NON_GRAN_STATES)
        length(ng_state_cells[[d]][[j]]) >= MIN_CELLS_PER_DONOR_STATE
      }, logical(1)))
    }, integer(1)),
    eligible_shared_state_names = vapply(acute_donors, function(d) {
      elig <- vapply(shared_states, function(s) {
        j <- match(s, STRICT_NON_GRAN_STATES)
        length(ng_state_cells[[d]][[j]]) >= MIN_CELLS_PER_DONOR_STATE
      }, logical(1))
      paste(shared_states[elig], collapse = ";")
    }, character(1)),
    neutrophil_fraction = vapply(acute_donors, function(d) {
      idxd <- which(donor_p == d)
      if (!length(idxd)) return(NA_real_)
      mean(state_p[idxd] %in% AXIS_STATES)
    }, numeric(1)),
    neutrophil_UMI_fraction = vapply(acute_donors, function(d) {
      idxd <- which(donor_p == d)
      if (!length(idxd)) return(NA_real_)
      sum(primary_umi[idxd][state_p[idxd] %in% AXIS_STATES]) / sum(primary_umi[idxd])
    }, numeric(1)),
    median_UMI = vapply(acute_donors, function(d) safe_median(primary_umi[donor_p == d]), numeric(1)),
    median_genes = vapply(acute_donors, function(d) safe_median(primary_gene[donor_p == d]), numeric(1)),
    stringsAsFactors = FALSE)
  dc_nd <- donor_cov[, setdiff(names(donor_cov), "donor"), drop = FALSE]
  
  # ---- 2.13 STATE-STANDARDIZED DONOR SIGNAL (v1.4/v1.5) ----------------------
  # Detection is computed per eligible (>= 20 cell) SHARED strict non-gran state
  # and averaged over the donor's eligible shared states. A donor with fewer
  # than MIN_SHARED_STATES_PER_DONOR eligible shared states yields NA (excluded
  # from correlation) so donors with tiny state coverage do not enter the test.
  # NOTE (v1.6): this is the mean across THAT donor's eligible shared states;
  # eligible state sets may differ across donors (recorded per donor).
  panel_state_signal_for_donor <- function(gset, d) {
    st_cells <- ng_state_cells[[d]]
    state_frac <- vapply(shared_states, function(s) {
      j <- match(s, STRICT_NON_GRAN_STATES)
      idx <- st_cells[[j]]
      if (length(idx) < MIN_CELLS_PER_DONOR_STATE) return(NA_real_)
      safe_mean(vapply(gset, function(g) mean(as.numeric(Xm[g, idx] > 0)), numeric(1)))
    }, numeric(1))
    elig <- is.finite(state_frac)
    if (sum(elig) < MIN_SHARED_STATES_PER_DONOR) return(NA_real_)
    mean(state_frac[elig])
  }
  
  donor_panel_signal <- function(panel_name, gset) {
    sig <- vapply(acute_donors, function(d) panel_state_signal_for_donor(gset, d), numeric(1))
    data.frame(
      donor = acute_donors,
      panel = panel_name,
      n_panel_genes_total = length(AMBIENT_PANELS[[panel_name]]),
      n_panel_genes_present = length(gset),
      missing_genes = paste(setdiff(AMBIENT_PANELS[[panel_name]], gset), collapse = ";"),
      panel_complete = panel_complete_flag[[panel_name]],
      dc_nd,
      strict_non_gran_detection_fraction = sig,
      stringsAsFactors = FALSE, row.names = NULL)
  }
  
  ambient_panel_presence_summary <- do.call(rbind, lapply(names(AMBIENT_PANELS), function(pn) {
    gset <- intersect(AMBIENT_PANELS[[pn]], markers)
    data.frame(
      panel = pn,
      n_total_genes = length(AMBIENT_PANELS[[pn]]),
      n_present_genes = length(gset),
      n_missing_genes = length(setdiff(AMBIENT_PANELS[[pn]], gset)),
      missing_genes = paste(setdiff(AMBIENT_PANELS[[pn]], gset), collapse = ";"),
      panel_complete = panel_complete_flag[[pn]],
      stringsAsFactors = FALSE)
  }))
  write_csv_atomic(ambient_panel_presence_summary,
                   file.path(TABLE_DIR, "ambient_panel_presence_summary.csv"))
  
  # ---- 2.14 CORRELATION ASSESSMENT (bootstrap CI + three-way verdict) --------
  ambient_cor_assess <- function(neut_fraction, signal, panel_complete = TRUE,
                                 min_donors = MIN_EVALUABLE_DONORS) {
    ok <- is.finite(neut_fraction) & is.finite(signal)
    base <- data.frame(n_donors_with_data = sum(ok),
                       rho = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                       p_value = NA_real_, padj = NA_real_,
                       verdict = NA_character_, stringsAsFactors = FALSE)
    if (!panel_complete) { base$verdict <- "not_evaluable_incomplete_panel"; return(base) }
    if (sum(ok) < min_donors) { base$verdict <- "not_evaluable_too_few_donors"; return(base) }
    if (near_constant(signal[ok])) {
      base$verdict <- "not_evaluable_zero_or_near_constant_non_gran_signal"; return(base)
    }
    if (near_constant(neut_fraction[ok])) {
      base$verdict <- "not_evaluable_near_constant_neutrophil_fraction"; return(base)
    }
    ct <- suppressWarnings(stats::cor.test(neut_fraction[ok], signal[ok],
                                           method = "spearman", exact = FALSE))
    rho <- unname(ct$estimate); pv <- ct$p.value
    x <- neut_fraction[ok]; y <- signal[ok]
    boot <- vapply(seq_len(BOOT_REPS), function(i) {
      b <- sample.int(length(x), replace = TRUE)
      if (length(unique(x[b])) < 2L || length(unique(y[b])) < 2L) return(NA_real_)
      stats::cor(rank(x[b]), rank(y[b]))
    }, numeric(1))
    ci <- stats::quantile(boot, c(0.025, 0.975), na.rm = TRUE)
    verdict <- if (ci[[2]] < STRONG_RHO_THRESHOLD) {
      "strong_positive_structure_not_supported"
    } else if (rho >= STRONG_RHO_THRESHOLD && ci[[1]] > 0) {
      "positive_structure_supported"
    } else "inconclusive"
    data.frame(n_donors_with_data = sum(ok), rho = rho, ci_low = ci[[1]],
               ci_high = ci[[2]], p_value = pv, padj = NA_real_,
               verdict = verdict, stringsAsFactors = FALSE)
  }
  
  spearman_pair <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10L) return(c(rho = NA_real_, p = NA_real_))
    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
  }
  depth_adjusted_spearman <- function(signal, neut_fraction, median_umi) {
    ok <- is.finite(signal) & is.finite(neut_fraction) & is.finite(median_umi)
    if (sum(ok) < 10L) return(c(rho = NA_real_, p = NA_real_))
    x <- rank(neut_fraction[ok]); y <- rank(signal[ok]); z <- rank(median_umi[ok])
    xr <- resid(stats::lm(x ~ z)); yr <- resid(stats::lm(y ~ z))
    ct <- suppressWarnings(stats::cor.test(yr, xr))
    c(rho = unname(ct$estimate), p = ct$p.value)
  }
  
  set.seed(BOOT_SEED)   # deterministic bootstrap across all panels
  ambient_panel_donor <- do.call(rbind, lapply(names(AMBIENT_PANELS), function(pn) {
    gset <- intersect(AMBIENT_PANELS[[pn]], markers)
    donor_panel_signal(pn, gset)
  }))
  write_csv_atomic(ambient_panel_donor,
                   file.path(TABLE_DIR, "ambient_vs_neutrophil_fraction_by_panel_donor.csv"))
  
  ambient_panel_cor <- do.call(rbind, lapply(
    split(ambient_panel_donor, ambient_panel_donor$panel), function(df) {
      cg <- ambient_cor_assess(df$neutrophil_fraction, df$strict_non_gran_detection_fraction,
                               panel_complete = unique(df$panel_complete))
      evaluable <- !is.na(cg$rho)   # gate sensitivities on the main verdict
      da <- if (evaluable) depth_adjusted_spearman(
        df$strict_non_gran_detection_fraction, df$neutrophil_fraction, df$median_UMI)
      else c(rho = NA_real_, p = NA_real_)
      uf <- if (evaluable) spearman_pair(
        df$neutrophil_UMI_fraction, df$strict_non_gran_detection_fraction)
      else c(rho = NA_real_, p = NA_real_)
      data.frame(
        panel = unique(df$panel),
        n_panel_genes_total = unique(df$n_panel_genes_total),
        n_panel_genes_present = unique(df$n_panel_genes_present),
        missing_genes = unique(df$missing_genes),
        panel_complete = unique(df$panel_complete),
        cg,
        rho_partial_umi_adjusted = unname(da["rho"]),
        p_partial_umi_adjusted = unname(da["p"]),
        rho_neut_umifrac_sensitivity = unname(uf["rho"]),
        p_neut_umifrac_sensitivity = unname(uf["p"]),
        stringsAsFactors = FALSE)
    }))
  
  # Pooled proxy over COMPLETE panels only; legacy, excluded from FDR family.
  complete_panels <- names(AMBIENT_PANELS)[panel_complete_flag]
  pooled_gset <- unlist(AMBIENT_PANELS[complete_panels], use.names = FALSE)
  sig_pooled <- vapply(acute_donors, function(d) panel_state_signal_for_donor(pooled_gset, d), numeric(1))
  ambient_pooled_donor <- data.frame(
    donor = acute_donors,
    panel = "all_ambient_markers_pooled_legacy",
    n_panel_genes_total = length(unique(unlist(AMBIENT_PANELS[complete_panels]))),
    n_panel_genes_present = length(pooled_gset),
    missing_genes = paste(setdiff(unique(unlist(AMBIENT_PANELS[complete_panels])), pooled_gset),
                          collapse = ";"),
    panel_complete = TRUE,
    dc_nd,
    strict_non_gran_detection_fraction = sig_pooled,
    stringsAsFactors = FALSE, row.names = NULL)
  cg_pooled <- ambient_cor_assess(ambient_pooled_donor$neutrophil_fraction,
                                  ambient_pooled_donor$strict_non_gran_detection_fraction,
                                  panel_complete = TRUE)
  evaluable_pooled <- !is.na(cg_pooled$rho)
  da_pooled <- if (evaluable_pooled) depth_adjusted_spearman(
    ambient_pooled_donor$strict_non_gran_detection_fraction,
    ambient_pooled_donor$neutrophil_fraction, ambient_pooled_donor$median_UMI)
  else c(rho = NA_real_, p = NA_real_)
  uf_pooled <- if (evaluable_pooled) spearman_pair(
    ambient_pooled_donor$neutrophil_UMI_fraction,
    ambient_pooled_donor$strict_non_gran_detection_fraction)
  else c(rho = NA_real_, p = NA_real_)
  ambient_pooled_cor <- data.frame(
    panel = "all_ambient_markers_pooled_legacy",
    n_panel_genes_total = length(unique(unlist(AMBIENT_PANELS[complete_panels]))),
    n_panel_genes_present = length(pooled_gset),
    missing_genes = paste(setdiff(unique(unlist(AMBIENT_PANELS[complete_panels])), pooled_gset),
                          collapse = ";"),
    panel_complete = TRUE,
    cg_pooled,
    rho_partial_umi_adjusted = unname(da_pooled["rho"]),
    p_partial_umi_adjusted = unname(da_pooled["p"]),
    rho_neut_umifrac_sensitivity = unname(uf_pooled["rho"]),
    p_neut_umifrac_sensitivity = unname(uf_pooled["p"]),
    stringsAsFactors = FALSE)
  
  ambient_panel_cor <- rbind(ambient_panel_cor, ambient_pooled_cor)
  # v1.6: BH-FDR over the PRE-SPECIFIED family of four formal panels (n = 4),
  #   matching the designated-comparison design. Not-evaluable panels carry
  #   p = NA and keep padj = NA; the remaining p-values are adjusted with the
  #   full pre-specified family size. The pooled legacy row is excluded.
  formal_panels <- c("primary_azurophilic", "secondary_specific",
                     "tertiary_gelatinase", "residual_E3_granulocyte_panel")
  fp <- ambient_panel_cor$panel %in% formal_panels & is.finite(ambient_panel_cor$p_value)
  ambient_panel_cor$padj[fp] <- stats::p.adjust(ambient_panel_cor$p_value[fp],
                                                method = "BH", n = length(formal_panels))
  write_csv_atomic(ambient_panel_cor,
                   file.path(TABLE_DIR, "ambient_vs_neutrophil_fraction_by_panel.csv"))
  
  # Preserve the v1.1 filename/columns for downstream compatibility (pooled).
  ambient_cor <- data.frame(
    donor = ambient_pooled_donor$donor,
    neutrophil_fraction = ambient_pooled_donor$neutrophil_fraction,
    non_gran_granule_detection_fraction = ambient_pooled_donor$strict_non_gran_detection_fraction,
    stringsAsFactors = FALSE)
  write_csv_atomic(ambient_cor, file.path(TABLE_DIR, "ambient_vs_neutrophil_fraction.csv"))
  
  # ---- 2.15 FIGURES (guarded; recorded in provenance) ------------------------
  figures_generated <- FALSE
  if (HAS_GGPLOT2) {
    suppressPackageStartupMessages(library(ggplot2))
    save_pdf_png <- function(device, path_base, width, height) {
      grDevices::pdf(file = paste0(path_base, ".pdf"), width = width, height = height)
      print(device); grDevices::dev.off()
      grDevices::png(file = paste0(path_base, ".png"), width = width, height = height,
                     units = "in", res = 300)
      print(device); grDevices::dev.off()
      invisible(TRUE)
    }
    sm <- sample_sum
    sm$sample_id <- factor(sm$sample_id, levels = sm$sample_id[order(sm$cohort, -sm$n_cells)])
    p1 <- ggplot(sm, aes(x = sample_id, y = n_cells, fill = cohort)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = c(Acute_sepsis = "#C8102E", Convalescent = "#F0A0A0",
                                   Healthy_control = "#4D9BE6", Surgery_control = "#B8B8B8")) +
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
            legend.position = "top") +
      labs(x = NULL, y = "Cells per sample",
           title = "GSE216009 sample-layer QC: retained cells per sample")
    save_pdf_png(p1, file.path(FIG_DIR, "fig1_cells_per_sample"), width = 8, height = 4)
    
    # Mutually exclusive compartments: axis | other granulocyte-like | strict non-gran
    ar <- ambient_report
    ar$marker <- factor(ar$marker, levels = markers)
    ar_plot <- data.frame(
      marker = rep(ar$marker, 3),
      panel  = rep(ar$panel, 3),
      compartment = rep(c("neutrophil axis",
                          "other granulocyte-like",
                          "strict non-granulocyte"), each = nrow(ar)),
      pooled_detection = c(ar$axis_pooled_detection,
                           ar$other_granulocyte_like_pooled_detection,
                           ar$strict_non_gran_pooled_detection),
      stringsAsFactors = FALSE)
    p2 <- ggplot(ar_plot, aes(x = marker, y = pooled_detection, color = compartment)) +
      geom_point(size = 2) +
      scale_color_manual(values = c("neutrophil axis" = "#C8102E",
                                    "other granulocyte-like" = "#B8B8B8",
                                    "strict non-granulocyte" = "#4D9BE6")) +
      theme_classic(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = NULL, y = "Pooled detection fraction (acute primary, 26 donors)",
           title = "Ambient proxy: granule/residual markers by compartment (mutually exclusive)")
    save_pdf_png(p2, file.path(FIG_DIR, "fig2_ambient_granule_signal"), width = 7, height = 4.5)
    figures_generated <- TRUE
  } else {
    cat("[warn] ggplot2 unavailable; figures skipped (recorded in provenance).\n")
  }
  
  # ---- 2.16 PROVENANCE -------------------------------------------------------
  script_path <- get_script_path()
  script_hash <- file_hash(script_path)
  frozen_manifest_hash <- if (!is.null(FROZEN_MANIFEST_FILE) && file.exists(FROZEN_MANIFEST_FILE))
    file_hash(FROZEN_MANIFEST_FILE) else NA_character_
  e3_registry_path <- if (!is.null(E3_RESIDUAL_REGISTRY_FILE) && file.exists(E3_RESIDUAL_REGISTRY_FILE))
    normalizePath(E3_RESIDUAL_REGISTRY_FILE) else NA_character_
  e3_registry_hash <- if (!is.null(E3_RESIDUAL_REGISTRY_FILE) && file.exists(E3_RESIDUAL_REGISTRY_FILE))
    file_hash(E3_RESIDUAL_REGISTRY_FILE) else NA_character_
  
  # v1.6: formal mode must have real SHA-256 for every required artifact.
  if (FORMAL_RUN) {
    if (is.na(input_rds_hash) || !nzchar(input_rds_hash)) stop_msg("Failed to compute input RDS SHA-256.")
    if (is.na(script_hash) || !nzchar(script_hash)) stop_msg("Failed to compute script SHA-256.")
    if (!is.null(FROZEN_MANIFEST_FILE) && (is.na(frozen_manifest_hash) || !nzchar(frozen_manifest_hash)))
      stop_msg("Failed to compute frozen manifest SHA-256.")
    if (!is.null(E3_RESIDUAL_REGISTRY_FILE) && (is.na(e3_registry_hash) || !nzchar(e3_registry_hash)))
      stop_msg("Failed to compute E3 registry SHA-256.")
  }
  
  prov_rows <- list()
  add_prov <- function(k, v) prov_rows[[length(prov_rows) + 1L]] <<- data.frame(
    key = k, value = as.character(v), stringsAsFactors = FALSE)
  add_prov("script", "qc_sample_summary_ambient.R")
  add_prov("version", SCRIPT_VERSION)
  add_prov("formal_run", FORMAL_RUN)
  add_prov("run_purpose", RUN_PURPOSE)
  add_prov("output_root", OUT_ROOT)
  add_prov("input_rds", INPUT_RDS)
  add_prov("input_rds_hash", input_rds_hash)
  add_prov("hash_algo_used", hash_algo_used)
  add_prov("script_path", script_path)
  add_prov("script_hash", script_hash)
  add_prov("run_status", "COMPLETE")
  add_prov("figures_generated", figures_generated)
  add_prov("frozen_state_set_check", frozen_state_check)
  add_prov("frozen_manifest_check", frozen_manifest_check)
  add_prov("frozen_manifest_file",
           if (is.null(FROZEN_MANIFEST_FILE)) NA_character_ else FROZEN_MANIFEST_FILE)
  add_prov("frozen_manifest_hash", frozen_manifest_hash)
  add_prov("e3_registry_check", e3_registry_check)
  add_prov("e3_registry_file", e3_registry_path)
  add_prov("e3_registry_hash", e3_registry_hash)
  add_prov("e3_registry_sha256_check",
           if (is.null(EXPECTED_E3_REGISTRY_SHA256)) "not_configured" else "verified")
  add_prov("e3_residual_amendment_id",
           if (is.null(E3_RESIDUAL_AMENDMENT_ID)) NA_character_ else E3_RESIDUAL_AMENDMENT_ID)
  add_prov("e3_residual_lock_date", E3_RESIDUAL_LOCK_DATE)
  add_prov("e3_evidence_tier", e3_evidence_tier)
  add_prov("e3_source_files", e3_source_files)
  add_prov("expected_total_cells", EXPECTED_N_TOTAL_CELLS)
  add_prov("total_cells_observed", ncol(obj))
  add_prov("n_samples", n_samples)
  add_prov("n_donors", n_donors)
  add_prov("n_acute_donors", length(acute_donors))
  add_prov("n_paired_donors", length(paired_donors))
  add_prov("n_primary_cells", length(primary_idx))
  add_prov("n_fine_states", length(obs_states))
  add_prov("ambient_panels", paste(names(AMBIENT_PANELS), collapse = ";"))
  add_prov("granulocyte_like_states_excluded_from_non_gran",
           paste(GRANULOCYTE_LIKE_STATES, collapse = ";"))
  add_prov("strict_non_gran_states", paste(STRICT_NON_GRAN_STATES, collapse = ";"))
  add_prov("shared_state_standard_n_donors_min", MIN_SHARED_STATE_DONOR_N)
  add_prov("shared_states", paste(shared_states, collapse = ";"))
  add_prov("min_cells_per_donor_state", MIN_CELLS_PER_DONOR_STATE)
  add_prov("min_shared_states_per_donor", MIN_SHARED_STATES_PER_DONOR)
  add_prov("min_evaluable_donors", MIN_EVALUABLE_DONORS)
  add_prov("strong_rho_threshold", STRONG_RHO_THRESHOLD)
  add_prov("bootstrap_reps", BOOT_REPS)
  add_prov("bootstrap_seed", BOOT_SEED)
  add_prov("residual_markers_not_evaluable",
           if (length(incomplete_panels))
             paste(setdiff(unlist(AMBIENT_PANELS[incomplete_panels]), markers), collapse = ";")
           else "none")
  for (pn in c("primary_azurophilic", "secondary_specific", "tertiary_gelatinase",
               "residual_E3_granulocyte_panel", "all_ambient_markers_pooled_legacy")) {
    r_ <- ambient_panel_cor[ambient_panel_cor$panel == pn, , drop = FALSE]
    if (!nrow(r_)) next
    r_ <- r_[1, , drop = FALSE]
    if (identical(pn, "all_ambient_markers_pooled_legacy")) {
      add_prov(paste0("verdict_", pn), paste0(
        r_$verdict, " (rho=", round(r_$rho, 4), ", CI ", round(r_$ci_low, 4), "-",
        round(r_$ci_high, 4), "); LEGACY backward-compat only, excluded from FDR"))
    } else {
      add_prov(paste0("verdict_", pn), paste0(
        r_$verdict, " (rho=", round(r_$rho, 4), ", CI ", round(r_$ci_low, 4), "-",
        round(r_$ci_high, 4), ", padj=", round(r_$padj, 4),
        ", rho_umi_adj=", round(r_$rho_partial_umi_adjusted, 4),
        ", rho_umifrac=", round(r_$rho_neut_umifrac_sensitivity, 4), ")"))
    }
  }
  add_prov("notes", paste0(
    "DESCRIPTIVE ambient proxy audit using the supplied processed object; no ",
    "additional ambient or doublet correction applied (non-correction ",
    "declaration); does NOT perform decontamination. Doublet sensitivity ",
    "handled as robustness (manuscript §3.9 / S18). Donor signal is the mean ",
    "across THAT donor's eligible shared strict non-gran states; eligible state ",
    "sets may differ across donors (per-donor state names in donor tables). ",
    "Bootstrap CIs are descriptive and NOT multiplicity-adjusted; BH-FDR uses the ",
    "pre-specified four-panel family (n = 4) and tests zero correlation only, not ",
    "the rho>=0.50 threshold verdict. Pooled proxy is legacy backward-compat only."))
  add_prov("R_version", R.version.string)
  prov <- do.call(rbind, prov_rows)
  write_csv_atomic(prov, file.path(PROV_DIR, "run_provenance.csv"))
  writeLines(capture.output(sessionInfo()), file.path(PROV_DIR, "sessionInfo.txt"))
  writeLines(paste("run_dir", OUT_ROOT, sep = "\t"), file.path(PROV_DIR, "run_dir.txt"))
  
  # ---- 2.17 CONSOLE SUMMARY ---------------------------------------------------
  cat("\n================ Single-cell QC summary + ambient PROXY audit ================\n")
  cat("  mode: ", if (FORMAL_RUN) "FORMAL" else "TEST", " | version ", SCRIPT_VERSION,
      "\n", sep = "")
  cat("-- Part A: sample-layer QC (whole atlas)\n")
  cat("  cells:", ncol(obj), "| samples:", n_samples, "| donors:", n_donors,
      "| states:", length(obs_states), "\n")
  cat("  median cells/sample:", object_summary$value[object_summary$metric == "median_cells_per_sample"],
      "(range", object_summary$value[object_summary$metric == "min_cells_per_sample"],
      "-", object_summary$value[object_summary$metric == "max_cells_per_sample"], ")\n")
  cat("  median UMI/cell:", object_summary$value[object_summary$metric == "median_UMI_per_cell_all"],
      "| median genes/cell:", object_summary$value[object_summary$metric == "median_genes_per_cell_all"], "\n")
  cat("  acute primary: ", length(primary_idx), " cells / ", length(acute_donors),
      " donors\n", sep = "")
  cat("\n-- Part B: ambient PROXY audit (acute primary, 26 donors)\n")
  cat("  shared-state standard: ", length(shared_states), "/", length(STRICT_NON_GRAN_STATES),
      " strict non-gran states; donor signal = mean across that donor's eligible ",
      "shared states (>= ", MIN_SHARED_STATES_PER_DONOR, "; names in tables)\n", sep = "")
  print(ambient_report, row.names = FALSE)
  cat("\n  Panel verdicts (rho | 95% CI | BH-padj | umi-adj partial rho | umifrac sens):\n")
  ambient_panel_cor <- ambient_panel_cor[order(ambient_panel_cor$panel), , drop = FALSE]
  for (i in seq_len(nrow(ambient_panel_cor))) {
    r_ <- ambient_panel_cor[i, ]
    tag <- if (grepl("legacy", r_$panel)) " [LEGACY]" else ""
    cat("    ", r_$panel, tag, ":\n", sep = "")
    cat("        n_donors=", r_$n_donors_with_data,
        " rho=", round(r_$rho, 4),
        " CI=", round(r_$ci_low, 4), "-", round(r_$ci_high, 4),
        " padj=", round(r_$padj, 4),
        " rho_umi_adj=", round(r_$rho_partial_umi_adjusted, 4),
        " rho_umifrac=", round(r_$rho_neut_umifrac_sensitivity, 4),
        " | ", r_$verdict, "\n", sep = "")
  }
  cat("\n  strict non-gran pooled detection by panel:\n")
  for (pn in names(AMBIENT_PANELS)) {
    v <- safe_mean(ambient_report$strict_non_gran_pooled_detection[
      ambient_report$panel == pn])
    cat("    [", pn, "] = ", round(v, 5), "\n", sep = "")
  }
  cat("\n  frozen checks: state set =", frozen_state_check,
      "| manifest =", frozen_manifest_check,
      "| e3 registry =", e3_registry_check,
      "| figures =", figures_generated, "\n")
  
  writeLines("COMPLETE", STATUS_FILE)
  run_status <- "COMPLETE"
  cat("\n[complete] Script QC/ambient ", SCRIPT_VERSION, " complete.\nRun dir: ",
      OUT_ROOT, "\n", sep = "")
  invisible(TRUE)
}

# ---- 3. RUN -----------------------------------------------------------------
main()