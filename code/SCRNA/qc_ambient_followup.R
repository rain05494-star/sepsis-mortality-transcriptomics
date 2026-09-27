# ============================================================================
# SCRIPT -- qc_ambient_followup.R  (v1.4)
#   Follow-up to qc_sample_summary_ambient.R (main QC; target v1.7).
#   Adjudicates what the donor-level positive ambient-proxy structure
#   (primary_azurophilic and residual-E3 panels) is:
#     Part 0  doublet-risk fail-closed validation + exclusion sets + manifest
#     Part 1  state-level localization (4 metrics)
#     Part 2  shared-state donor coverage
#     Part 3  leave-one-donor-out jackknife (20-donor standard)
#     Part 4  donor-level ambient rho, DUAL-BRANCH doublet sensitivity
#             (per_sample primary; global_acute_primary legacy)
#     Part 5  M1 count-UMI contribution fraction, DUAL-BRANCH doublet
#             sensitivity + PAIRED donor bootstrap (baseline replication
#             fail-closed; no arbitrary formal-equivalence threshold)
#
#   Frozen decisions (2026-08-13, v1.4 final):
#     - QC PREFLIGHT runs before output-dir creation and RDS loading
#       (fail-fast; no FAILED dir is left when QC_RUN_DIR/version/purpose/
#       formal_run/hash-algo/input-RDS-hash/script-path is wrong).
#     - QC_RUN_DIR must point at the directory of the ACTUAL v1.7 main-QC
#       TEST re-run. EXPECTED_QC_VERSION (v1.7), EXPECTED_QC_PURPOSE
#       (singlecell_qc_TEST), formal_run (FALSE) and hash_algo_used (sha256)
#       are enforced from the QC provenance. formal_run is validated as
#       exactly TRUE/FALSE before any comparison (a malformed value such as
#       "foo"/"not_configured" is a hard stop, never silently treated as FALSE).
#     - The audit must be launched with Rscript, not source(), because the
#       executing script path is SHA-256 hashed.
#     - per_sample is the PRIMARY doublet scope; global_acute_primary is the
#       legacy/manuscript-replication branch. Both always computed.
#     - 26 acute samples == 26 acute donors (1 acute sample per donor), so
#       per_sample scope == per_donor scope in this cohort (recorded).
#     - Exact-k removal: ceiling(n_group * fraction); ties broken by
#       order(-doublet_score, cell_barcode, method="radix"); NOT score>=cutoff.
#       Strict expected totals: per_sample_top5 7607, per_sample_top10 15200,
#       global_top5 7592, global_top10 15184 (fail-closed verified).
#     - M1 readout is count-UMI contribution fraction; the gene set is FROZEN
#       to exactly the 11 Script-02-verified genes (no manifest override).
#     - Part 5 FULL group MUST reproduce pre-audit Script 02/S17 baselines
#       progenitor=0.346059, broad_residual=0.348480 within 1e-6, or stop.
#     - Bootstrap is PAIRED: a single donor-resample matrix (sampling the
#       donors commonly evaluable in ALL groups) is shared by every group.
#       Point estimates are reported TWICE: diff_*_all_evaluable (that
#       group's own evaluable donors) and diff_*_common_donors (exactly the
#       donor set the paired CI is computed on). The paired CI always
#       corresponds to the common-donor estimate; if all groups are 26/26
#       the two coincide. Paired requires >= MIN_EVALUABLE_DONORS common
#       donors, else no CI is claimed.
#     - Bootstrap "P" is a bootstrap SUPPORT PROBABILITY, not a hypothesis P.
#     - doublet_score is a continuous doublet-RISK ranking score, not a
#       validated doublet label; predicted_doublet must be non-missing and
#       all-FALSE in acute-primary cells (else stop).
#     - Input hashing is fixed to SHA-256 (digest is required); the follow-up
#       only accepts a main-QC run that records hash_algo_used = sha256.
#     - Cross-check compares evaluable donor IDENTITY per panel against the
#       QC donor-level table (not just donor counts).
#     - FORMAL follow-up is BLOCKED until frozen state manifest / E3 registry /
#       expected SHA-256 / FORMAL main-QC provenance are configured.
#
#   CHANGELOG
#     v1.0  initial four parts + cross-check.
#     v1.1  20-donor minimum in spearman_rho; formal panel-completeness hard
#           stop; corrected state-contribution metrics; deletion-composition
#           reporting; hardened cross-check; scope param.
#     v1.2  DUAL-BRANCH doublet sensitivity; Step A fail-closed validation +
#           per-cell exclusion manifest; ceiling + radix tie-break frozen;
#           count-based M1 Part 5 with 03c bootstrap and fail-closed baseline
#           replication; RHO_TOL=1e-6.
#     v1.3  (first review) complete file structure; M1_GENES frozen+validated;
#           QC version gate; count-column order check; scope checks in all
#           modes; predicted_doublet non-missing/all-FALSE; Part 4 split into
#           panel-summary + donor-signals-long; exploratory stability band;
#           PAIRED bootstrap; exact-k totals verified; FORMAL blocked.
#     v1.4  (second + third review, RC -> frozen) QC PREFLIGHT fully moved
#           before output dir + RDS load (dir/STATUS/version/purpose/
#           formal_run/hash-algo/input-RDS-hash/script-path, fail-fast);
#           hash_algo_used=sha256 required; TEST identity (run_purpose +
#           formal_run) enforced; formal_run validated as exactly TRUE/FALSE
#           (malformed value = hard stop); Rscript required (script hashed
#           before any computation); M1 manifest override removed; paired
#           bootstrap restricted to donors commonly evaluable in all groups
#           with BOTH all-evaluable and common-donor point estimates reported
#           and the paired CI tied to the common-donor estimate; file_hash
#           fixed to SHA-256 via required digest; cross-check compares
#           evaluable donor IDENTITY per panel (one row per panel enforced);
#           jackknife all-NA protection; rank labels unified,
#           rank_order_all_evaluable / rank_order_common_donors recorded;
#           header RUN SEQUENCE syntax-check line single-line (no dangling
#           quote) so the file parses.
#
#   Non-correction declaration: read-only diagnostic (no cell removal from the
#   primary analysis, no counts modification, no re-clustering). Doublet
#   sensitivity is DESCRIPTIVE robustness, not a correction.
#
#   RUN SEQUENCE (required, in order):
#     (1) clean main QC script qc_sample_summary_ambient.R to v1.7 (provenance
#         must record run_purpose=singlecell_qc_TEST, formal_run=FALSE,
#         hash_algo_used=sha256), re-run its TEST, and note the console-printed
#         output directory;
#     (2) set QC_RUN_DIR below to that ACTUAL v1.7 directory (a pre-existing
#         .../20260813_singlecell_qc_TEST produced by v1.6 will correctly fail
#         the preflight);
#     (3) syntax check:
#         Rscript -e 'parse(file="~/predicate/singlecell/script/qc_ambient_followup.R"); cat("PARSE_OK\n")'
#     (4) run: Rscript ~/predicate/singlecell/script/qc_ambient_followup.R
#         (NOT source(); the script hash requires Rscript).
# ============================================================================

options(stringsAsFactors = FALSE)

# ---- 0. CONFIG --------------------------------------------------------------
SCRIPT_VERSION <- "v1.4"
FORMAL_RUN     <- FALSE

OUTPUT_ROOT  <- "/home/sunshine/predicate/singlecell/outputs/GSE216009_32gene_analysis"
INPUT_RDS    <- file.path("/home/sunshine/predicate/singlecell",
                          "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
STAGE_NAME   <- "qc_ambient_followup"
RUN_PURPOSE  <- if (FORMAL_RUN) "ambient_followup_FORMAL" else "ambient_followup_TEST"

# MUST point at the directory printed by the v1.7 main-QC TEST re-run (step 1).
# The existing .../20260813_singlecell_qc_TEST was produced by v1.6 and will
# (correctly) fail the preflight. Do not guess the suffix.
QC_RUN_DIR <- "/home/sunshine/predicate/singlecell/outputs/GSE216009_32gene_analysis/qc_sample_summary_ambient/20260813_singlecell_qc_TEST_5"
EXPECTED_QC_VERSION   <- "v1.7"
EXPECTED_QC_PURPOSE   <- "singlecell_qc_TEST"
EXPECTED_QC_FORMAL_RUN <- FALSE
EXPECTED_QC_HASH_ALGO <- "sha256"

DOUBLET_RISK_FIELD         <- "doublet_score"
DOUBLET_PRIMARY_SCOPE      <- "per_sample"
DOUBLET_SENSITIVITY_SCOPES <- c("per_sample", "global_acute_primary")
DOUBLET_TOP_FRACTIONS      <- c(0.05, 0.10)
EXACT_K_RULE               <- "ceiling"
TIE_BREAK_RULE             <- "doublet_score_desc_then_cell_barcode_asc"
M1_READOUT                 <- "count_UMI_contribution_fraction"
M1_UCELL_AVAILABLE         <- FALSE

# ---- 11 M1 cell-cycle genes: FROZEN (2026-08-13, user-verified vs Script 02) ----
# The frozen set is the ONLY entry point; there is no manifest override.
FROZEN_M1_GENES <- c("ANLN","ASPM","BUB1","CCNA2","CKS2","KIF11",
                     "MKI67","RRM2","TOP2A","TPX2","TYMS")
M1_GENES <- FROZEN_M1_GENES

# Pre-audit Script 02 / S17 baselines (user independent replication:
# progenitor = 0.346058858326, broad residual = 0.348480251181; the 6-dp
# constants below satisfy tolerance 1e-6 by ~1.4e-7 / ~2.5e-7).
M1_BASELINE <- c(progenitor_granulopoiesis = 0.346059, broad_residual = 0.348480)
M1_BASELINE_TOL <- 1e-6
M1_BOOT_REPS    <- 10000L
M1_BOOT_SEED    <- 20260804L
# EXPLORATORY only: a 5-percentage-point |observed diff change| flag. It is a
# descriptive column (exploratory_5pp_rule), NOT a formal equivalence verdict.
# It is evaluated on the PAIRED common-donor point estimate.
M1_EXPLORATORY_STABILITY_BAND <- 0.05

M1_COMPARTMENTS <- list(
  progenitor_granulopoiesis = c("HSPCs",
                                "Cycling_neutrophil_progenitors",
                                "MPO+_immature_neutrophils_or_progenitors"),
  immature_neutrophil = c("PADI4+_immature_neutrophils",
                          "IL1R2+_immature_neutrophils",
                          "S100A8-9_hi_neutrophils"),
  cycling_lymphoid = "Cycling_TNK")
# broad_residual = remaining states (computed from observed 24)

COUNT_ASSAY      <- "RNA"
COUNT_LAYER      <- "counts"
FINE_STATE_FIELD <- "fine_annot"
SAMPLE_FIELD     <- "sample_id"
CONDITION_FIELD  <- "diagnosis"
ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")

MIN_CELLS_PER_DONOR_STATE   <- 20L
MIN_SHARED_STATES_PER_DONOR <- 3L
MIN_EVALUABLE_DONORS        <- 20L
STRONG_RHO_THRESHOLD        <- 0.50
MIN_SHARED_STATE_DONOR_N    <- 13L
RHO_TOL                     <- 1e-6

EXPECTED_N_TOTAL_CELLS   <- 272993L
EXPECTED_N_FINE_STATES   <- 24L
EXPECTED_N_PRIMARY_CELLS <- 151837L
EXPECTED_N_ACUTE_DONORS  <- 26L

# Strict exact-k removal totals (deterministic; fail-closed verified at run)
EXPECTED_EXCLUSION_N <- c(
  per_sample_top5             = 7607L,
  per_sample_top10            = 15200L,
  global_acute_primary_top5   = 7592L,
  global_acute_primary_top10  = 15184L)

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
  secondary_specific  = c("LTF","LCN2"),
  tertiary_gelatinase = c("MMP8","MMP9"))

GRANULOCYTE_EXTRA_STATES <- c("Eosinophils", "Mast_cells/eosiniophils")
GRANULOCYTE_LIKE_STATES  <- unique(c(AXIS_STATES, GRANULOCYTE_EXTRA_STATES))

AMBIENT_PANELS <- c(
  GRANULE_PANELS,
  list(residual_E3_granulocyte_panel = c("CEACAM8","OLFM4","RNASE3","TCN1")))
REQUIRED_AMBIENT_PANELS <- c("primary_azurophilic","secondary_specific","tertiary_gelatinase")
FORMAL_PANELS <- c("primary_azurophilic","secondary_specific",
                   "tertiary_gelatinase","residual_E3_granulocyte_panel")

# ---- 1. PACKAGES + HELPERS --------------------------------------------------
needed <- c("SeuratObject","Matrix","digest")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({ library(SeuratObject); library(Matrix); library(digest) })

stop_msg <- function(...) stop(paste0(...), call. = FALSE)
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
  flush(con); close(con); con_closed <- TRUE
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
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  w <- grep("^--file=", args, value = TRUE)
  if (length(w)) return(tryCatch(normalizePath(sub("^--file=", "", w[1])),
                                 error = function(e) sub("^--file=", "", w[1])))
  NA_character_
}
# SHA-256 ONLY (fixed algorithm; digest is a required package)
file_hash <- function(path) {
  if (!file.exists(path)) stop_msg("Cannot hash missing file: ", path)
  h <- tryCatch(digest::digest(path, algo = "sha256", file = TRUE),
                error = function(e) NA_character_)
  if (is.na(h) || !nzchar(h)) stop_msg("Failed to compute SHA-256: ", path)
  h
}
spearman_rho <- function(x, y, min_donors = MIN_EVALUABLE_DONORS) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_donors) return(NA_real_)
  if (near_constant(x[ok]) || near_constant(y[ok])) return(NA_real_)
  unname(suppressWarnings(stats::cor.test(x[ok], y[ok],
                                          method = "spearman", exact = FALSE)$estimate))
}

# ---- 2. MAIN -----------------------------------------------------------------
main <- function() {
  
  # ---- 2.0 FROZEN SCOPE VALIDATION (all modes, not only FORMAL) ----------------
  if (!identical(DOUBLET_PRIMARY_SCOPE, "per_sample"))
    stop_msg("DOUBLET_PRIMARY_SCOPE must be 'per_sample'.")
  if (!identical(DOUBLET_SENSITIVITY_SCOPES, c("per_sample", "global_acute_primary")))
    stop_msg("Unexpected DOUBLET_SENSITIVITY_SCOPES: ",
             paste(DOUBLET_SENSITIVITY_SCOPES, collapse = ","))
  if (any(!is.finite(DOUBLET_TOP_FRACTIONS)) ||
      any(DOUBLET_TOP_FRACTIONS <= 0 | DOUBLET_TOP_FRACTIONS >= 1) ||
      anyDuplicated(DOUBLET_TOP_FRACTIONS))
    stop_msg("Invalid DOUBLET_TOP_FRACTIONS.")
  
  # ---- 2.1 FORMAL GATE ---------------------------------------------------------
  if (FORMAL_RUN) {
    stop_msg(paste0("FORMAL follow-up is BLOCKED until the frozen 24-state ",
                    "manifest, E3 registry/amendment ID, expected SHA-256 and ",
                    "a FORMAL main-QC provenance are configured. This version ",
                    "supports TEST only."))
  }
  
  # ---- 2.2 QC PREFLIGHT (fail-fast; BEFORE output dir and RDS loading) ----------
  # (a) QC_RUN_DIR + STATUS + provenance identity (TEST only, version, hash algo)
  if (!is.character(QC_RUN_DIR) || length(QC_RUN_DIR) != 1L || !nzchar(QC_RUN_DIR) ||
      grepl("^<.*>$", QC_RUN_DIR) || !dir.exists(QC_RUN_DIR))
    stop_msg("QC_RUN_DIR is not configured to an existing v1.7 QC directory: ", QC_RUN_DIR)
  status_path <- file.path(QC_RUN_DIR, "STATUS.txt")
  qc_prov_path <- file.path(QC_RUN_DIR, "provenance", "run_provenance.csv")
  if (!file.exists(status_path)) stop_msg("QC run STATUS.txt not found: ", status_path)
  if (!identical(trimws(paste(readLines(status_path, warn = FALSE), collapse = "")), "COMPLETE"))
    stop_msg("QC run is not COMPLETE: ", QC_RUN_DIR)
  if (!file.exists(qc_prov_path)) stop_msg("QC provenance not found: ", qc_prov_path)
  qp <- read.csv(qc_prov_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("key", "value") %in% names(qp)))
    stop_msg("QC provenance must contain key/value columns.")
  if (anyDuplicated(qp$key)) stop_msg("Duplicated keys in QC provenance.")
  get_qc_prov <- function(k) {
    z <- qp$value[qp$key == k]
    if (length(z) != 1L || is.na(z) || !nzchar(z))
      stop_msg("Required QC provenance key missing/invalid: ", k)
    z
  }
  qc_version    <- get_qc_prov("version")
  qc_purpose    <- get_qc_prov("run_purpose")
  qc_hash       <- get_qc_prov("input_rds_hash")
  qc_shared     <- get_qc_prov("shared_states")
  qc_hash_algo  <- get_qc_prov("hash_algo_used")
  qc_formal_run <- get_qc_prov("formal_run")
  if (!identical(qc_version, EXPECTED_QC_VERSION))
    stop_msg("QC version mismatch: expected ", EXPECTED_QC_VERSION,
             "; observed ", qc_version, ". Re-run main QC v1.7 and update QC_RUN_DIR.")
  if (!identical(qc_hash_algo, EXPECTED_QC_HASH_ALGO))
    stop_msg("QC hash_algo_used mismatch: expected ", EXPECTED_QC_HASH_ALGO,
             "; observed ", qc_hash_algo, ". Main QC v1.7 must record SHA-256.")
  if (!identical(qc_purpose, EXPECTED_QC_PURPOSE))
    stop_msg("QC run_purpose mismatch: expected ", EXPECTED_QC_PURPOSE,
             "; observed ", qc_purpose, ". Only a singlecell_qc_TEST run is accepted.")
  # formal_run: validate as exactly TRUE/FALSE BEFORE any comparison, so a
  # malformed value ("foo"/"not_configured") is a hard stop, never silently
  # coerced to FALSE.
  qc_formal_norm <- tolower(trimws(qc_formal_run))
  if (!(qc_formal_norm %in% c("true", "false")))
    stop_msg("QC formal_run must be exactly TRUE or FALSE; observed ", qc_formal_run)
  qc_formal_logical <- identical(qc_formal_norm, "true")
  if (!identical(qc_formal_logical, EXPECTED_QC_FORMAL_RUN))
    stop_msg("QC formal_run mismatch: expected ", EXPECTED_QC_FORMAL_RUN,
             "; observed ", qc_formal_run, ". Only a non-formal TEST run is accepted.")
  
  # (b) Rscript enforcement (executing script must be hashable, fail-fast)
  script_path <- get_script_path()
  if (length(script_path) != 1L || is.na(script_path) || !nzchar(script_path) ||
      !file.exists(script_path))
    stop_msg("This audit must be run with Rscript, not source(), ",
             "so the executing script can be SHA-256 hashed.")
  script_hash <- file_hash(script_path)
  
  # (c) input RDS hash BEFORE output-dir creation (no FAILED dir on mismatch)
  if (!file.exists(INPUT_RDS)) stop_msg("Input RDS not found: ", INPUT_RDS)
  input_rds_hash <- file_hash(INPUT_RDS)
  if (!identical(qc_hash, input_rds_hash))
    stop_msg("QC run input RDS hash differs from current object: expected ",
             qc_hash, ", observed ", input_rds_hash)
  
  # ---- 2.3 OUTPUT ROOT + STATUS -------------------------------------------------
  OUT_ROOT  <- make_run_dir(OUTPUT_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  PROV_DIR  <- file.path(OUT_ROOT, "provenance")
  for (d in c(TABLE_DIR, PROV_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  STATUS_FILE <- file.path(OUT_ROOT, "STATUS.txt")
  writeLines("RUNNING", STATUS_FILE)
  run_status <- "RUNNING"
  on.exit(if (identical(run_status, "RUNNING")) {
    try(writeLines("FAILED", STATUS_FILE), silent = TRUE)
  }, add = TRUE)
  cat("[setup] ", RUN_PURPOSE, " output root:\n  ", OUT_ROOT, "\n", sep = "")
  
  # ---- 2.4 READ OBJECT + INVARIANTS + COUNT ORDER (hash already verified) ---------
  obj <- readRDS(INPUT_RDS)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  if (ncol(obj) != EXPECTED_N_TOTAL_CELLS)
    stop_msg("Expected ", EXPECTED_N_TOTAL_CELLS, " total cells; observed ", ncol(obj))
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay = COUNT_ASSAY, layer = COUNT_LAYER, fast = FALSE)
  } else if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, layer = COUNT_LAYER)
  } else {
    cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, slot = COUNT_LAYER)
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  # count-matrix cell order MUST equal the Seurat object's cell order
  if (!identical(colnames(cnt), colnames(obj)))
    stop_msg("Count matrix columns are not identical to Seurat object cells.")
  for (f in c(SAMPLE_FIELD, CONDITION_FIELD, FINE_STATE_FIELD))
    if (!(f %in% colnames(meta))) stop_msg("meta.data missing: ", f)
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  obs_states <- sort(unique(cell_state))
  if (length(obs_states) != EXPECTED_N_FINE_STATES)
    stop_msg("Expected ", EXPECTED_N_FINE_STATES, " fine states; observed ", length(obs_states))
  
  # ---- 2.5 COMPARTMENTS -----------------------------------------------------------
  gran_extra_missing <- setdiff(GRANULOCYTE_EXTRA_STATES, obs_states)
  if (length(gran_extra_missing))
    stop_msg("Granulocyte-like extra states absent: ", paste(gran_extra_missing, collapse = ", "))
  GRANULOCYTE_LIKE_STATES <- intersect(GRANULOCYTE_LIKE_STATES, obs_states)
  STRICT_NON_GRAN_STATES <- setdiff(obs_states, GRANULOCYTE_LIKE_STATES)
  if (!length(STRICT_NON_GRAN_STATES)) stop_msg("No strict non-granulocyte states.")
  
  # ---- 2.6 COHORT REBUILD (fail-closed, identical to main) -------------------------
  sample_tab <- unique(data.frame(sample_id = cell_sample, diagnosis = cell_diag,
                                  stringsAsFactors = FALSE))
  n_diag_per_sample <- tapply(sample_tab$diagnosis, sample_tab$sample_id,
                              function(z) length(unique(z)))
  if (any(n_diag_per_sample != 1L)) stop_msg("Samples with multiple diagnosis values.")
  sample_tab$donor  <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort <- NA_character_
  sample_tab$cohort[sample_tab$diagnosis %in% ACUTE_DIAG] <- "Acute_sepsis"
  sample_tab$cohort[sample_tab$diagnosis == "Conv"] <- "Convalescent"
  sample_tab$cohort[sample_tab$diagnosis == "HV"]   <- "Healthy_control"
  sample_tab$cohort[sample_tab$diagnosis == "CS"]   <- "Surgery_control"
  if (anyNA(sample_tab$cohort))
    stop_msg("Unmapped diagnosis: ",
             paste(unique(sample_tab$diagnosis[is.na(sample_tab$cohort)]), collapse = ", "))
  acute_donors <- sort(unique(sample_tab$donor[sample_tab$cohort == "Acute_sepsis"]))
  if (length(acute_donors) != EXPECTED_N_ACUTE_DONORS)
    stop_msg("Expected ", EXPECTED_N_ACUTE_DONORS, " acute donors; observed ", length(acute_donors))
  m <- match(cell_sample, sample_tab$sample_id)
  if (anyNA(m)) stop_msg("Sample->donor resolution failed.")
  donor_all  <- sample_tab$donor[m]
  cohort_all <- sample_tab$cohort[m]
  is_primary <- cohort_all == "Acute_sepsis"
  primary_idx <- which(is_primary)
  if (length(primary_idx) != EXPECTED_N_PRIMARY_CELLS)
    stop_msg("Expected ", EXPECTED_N_PRIMARY_CELLS, " primary cells; observed ", length(primary_idx))
  n_primary <- length(primary_idx)
  
  state_p  <- cell_state[primary_idx]
  donor_p  <- donor_all[primary_idx]
  sample_of_primary <- cell_sample[primary_idx]
  cell_id_primary   <- rownames(meta)[primary_idx]
  # per_sample == per_donor holds iff exactly one acute sample per acute donor
  spd <- tapply(sample_of_primary, donor_p, function(z) length(unique(z)))
  if (length(spd) != EXPECTED_N_ACUTE_DONORS || any(spd != 1L))
    stop_msg("Expected 26 acute donors each with exactly one acute sample (per_sample != per_donor).")
  
  # ---- 2.7 MARKERS + Xm + PANEL COMPLETENESS --------------------------------------
  panel_of <- setNames(rep(names(AMBIENT_PANELS), lengths(AMBIENT_PANELS)),
                       unlist(AMBIENT_PANELS, use.names = FALSE))
  dup_panel_marker <- unique(names(panel_of)[duplicated(names(panel_of))])
  if (length(dup_panel_marker))
    stop_msg("Ambient markers in multiple panels: ", paste(dup_panel_marker, collapse = ", "))
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
    stop_msg("Required ambient markers absent: ", paste(absent_gran, collapse = ", "))
  gene_map <- gene_map[gene_map$present, , drop = FALSE]
  markers <- gene_map$gene
  panel_complete_flag <- vapply(names(AMBIENT_PANELS), function(pn)
    all(AMBIENT_PANELS[[pn]] %in% markers), logical(1))
  formal_incomplete <- FORMAL_PANELS[!panel_complete_flag[FORMAL_PANELS]]
  if (length(formal_incomplete))
    stop_msg("Formal follow-up panels incomplete: ",
             paste(formal_incomplete, collapse = ", "))
  gi <- match(gene_map$dataset_feature, feat)
  stopifnot(!anyNA(gi), !anyDuplicated(gi))
  Xm <- cnt[gi, primary_idx, drop = FALSE]
  rownames(Xm) <- markers
  
  non_gran_idx_p <- which(state_p %in% STRICT_NON_GRAN_STATES)
  is_axis <- state_p %in% AXIS_STATES
  is_ng   <- state_p %in% STRICT_NON_GRAN_STATES
  is_ogl  <- state_p %in% GRANULOCYTE_EXTRA_STATES
  
  # ---- 2.8 DONOR-STATE ELIGIBILITY + SHARED STATES + COVARIATES --------------------
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
  is_shared <- state_p %in% shared_states
  if (!identical(qc_shared, paste(shared_states, collapse = ";")))
    stop_msg("QC run shared states differ from current rebuild.")
  
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
    stringsAsFactors = FALSE)
  nf_full <- donor_cov$neutrophil_fraction
  
  compute_donor_signals <- function(keep_idx, gset) {
    vapply(acute_donors, function(d) {
      st_cells <- ng_state_cells[[d]]
      state_frac <- vapply(shared_states, function(s) {
        j <- match(s, STRICT_NON_GRAN_STATES)
        idx_full <- st_cells[[j]]
        if (!length(idx_full)) return(NA_real_)
        idx <- intersect(idx_full, keep_idx)
        if (length(idx) < MIN_CELLS_PER_DONOR_STATE) return(NA_real_)
        safe_mean(vapply(gset, function(g) mean(as.numeric(Xm[g, idx] > 0)), numeric(1)))
      }, numeric(1))
      elig <- is.finite(state_frac)
      if (sum(elig) < MIN_SHARED_STATES_PER_DONOR) return(NA_real_)
      mean(state_frac[elig])
    }, numeric(1))
  }
  
  # ---- 2.9 PART 0 -- DOUBLET-RISK VALIDATION + EXCLUSION SETS ---------------------
  ds_all <- as.numeric(as.character(meta[[DOUBLET_RISK_FIELD]]))
  if (length(ds_all) != nrow(meta) || anyNA(ds_all) || any(!is.finite(ds_all)))
    stop_msg("doublet_score is not fully numeric/finite.")
  ds_primary <- ds_all[primary_idx]
  if (near_constant(ds_primary)) stop_msg("doublet_score (near-)constant in primary.")
  
  if (!("predicted_doublet" %in% colnames(meta))) stop_msg("predicted_doublet missing.")
  pd_primary <- as.character(meta[["predicted_doublet"]][primary_idx])
  if (anyNA(pd_primary) || !all(pd_primary == "FALSE"))
    stop_msg("predicted_doublet must be non-missing and all-FALSE in acute-primary cells.")
  cat("[Part0] predicted_doublet in primary: all-FALSE (", sum(pd_primary %in% "FALSE"),
      " cells) -- dead/unusable, recorded only.\n", sep = "")
  
  sm_tab <- data.frame(
    sample = sort(unique(sample_of_primary)),
    score_median = vapply(sort(unique(sample_of_primary)), function(z)
      median(ds_primary[sample_of_primary == z]), numeric(1)),
    score_q95 = vapply(sort(unique(sample_of_primary)), function(z)
      unname(quantile(ds_primary[sample_of_primary == z], 0.95)), numeric(1)),
    stringsAsFactors = FALSE)
  write_csv_atomic(sm_tab, file.path(TABLE_DIR, "doublet_score_per_sample.csv"))
  
  batch_R2 <- NA_real_
  if ("batch" %in% colnames(meta)) {
    batch <- factor(meta[["batch"]][primary_idx])
    if (anyNA(batch)) stop_msg("batch contains NA in primary.")
    if (nlevels(batch) >= 2L)
      batch_R2 <- summary(lm(rank(ds_primary, ties.method = "average") ~ batch))$r.squared
  }
  
  build_exclusion <- function(scope, fraction) {
    groups <- if (scope == "per_sample")
      split(seq_len(n_primary), sample_of_primary) else
        list(all = seq_len(n_primary))
    removed <- integer(0); cutoffs <- numeric(); gnames <- character()
    for (g in names(groups)) {
      idx <- groups[[g]]
      k <- ceiling(length(idx) * fraction)
      if (k < 1L) next
      ord <- idx[order(-ds_primary[idx], cell_id_primary[idx], method = "radix")]
      sel <- ord[seq_len(k)]
      removed <- c(removed, sel)
      cutoffs <- c(cutoffs, min(ds_primary[sel]))
      gnames <- c(gnames, g)
    }
    names(cutoffs) <- gnames
    list(scope = scope, fraction = fraction,
         removed = removed, keep = setdiff(seq_len(n_primary), removed),
         n_removed = length(removed), cutoffs = cutoffs)
  }
  
  exclusion_sets <- list()
  for (sc in DOUBLET_SENSITIVITY_SCOPES)
    for (fr in DOUBLET_TOP_FRACTIONS)
      exclusion_sets[[paste0(sc, "_top", round(100 * fr))]] <- build_exclusion(sc, fr)
  
  # strict exact-k totals (fail-closed, not just printed)
  for (nm in names(EXPECTED_EXCLUSION_N)) {
    obs_n <- as.integer(exclusion_sets[[nm]]$n_removed)
    if (!identical(obs_n, as.integer(EXPECTED_EXCLUSION_N[[nm]])))
      stop_msg("Exclusion count mismatch for ", nm, ": expected ",
               EXPECTED_EXCLUSION_N[[nm]], ", observed ", obs_n)
  }
  cat("[Part0] exact-k exclusion totals verified: ",
      paste(sprintf("%s=%d", names(EXPECTED_EXCLUSION_N), EXPECTED_EXCLUSION_N), collapse=", "),
      "\n", sep = "")
  
  global_cut <- vapply(DOUBLET_TOP_FRACTIONS, function(fr) {
    k <- ceiling(n_primary * fr)
    ord <- seq_len(n_primary)[order(-ds_primary, cell_id_primary, method = "radix")]
    min(ds_primary[ord[seq_len(k)]])
  }, numeric(1))
  names(global_cut) <- paste0("top", round(100 * DOUBLET_TOP_FRACTIONS))
  
  # per-cell exclusion manifest (all primary cells x 4 scope/fraction combos)
  manifest_rows <- list()
  for (sc in DOUBLET_SENSITIVITY_SCOPES) {
    for (fr in DOUBLET_TOP_FRACTIONS) {
      es <- exclusion_sets[[paste0(sc, "_top", round(100 * fr))]]
      rem <- logical(n_primary); rem[es$removed] <- TRUE
      if (sc == "per_sample") {
        per_cut <- es$cutoffs[sample_of_primary]
        names(per_cut) <- NULL
        if (anyNA(per_cut)) stop_msg("per_sample cutoff missing for some sample.")
        man <- data.frame(
          cell_barcode = cell_id_primary, sample_id = sample_of_primary,
          donor = donor_p, doublet_score = ds_primary, scope = sc, fraction = fr,
          sample_cutoff = per_cut, global_cutoff = NA_real_,
          selected_for_exclusion = rem, tied_at_cutoff = ds_primary == per_cut,
          stringsAsFactors = FALSE)
      } else {
        gc <- global_cut[[paste0("top", round(100 * fr))]]
        man <- data.frame(
          cell_barcode = cell_id_primary, sample_id = sample_of_primary,
          donor = donor_p, doublet_score = ds_primary, scope = sc, fraction = fr,
          sample_cutoff = NA_real_, global_cutoff = gc,
          selected_for_exclusion = rem, tied_at_cutoff = ds_primary == gc,
          stringsAsFactors = FALSE)
      }
      manifest_rows[[length(manifest_rows) + 1L]] <- man
    }
  }
  manifest_df <- do.call(rbind, manifest_rows)
  write_csv_atomic(manifest_df, file.path(TABLE_DIR, "doublet_exclusion_cell_manifest.csv"))
  
  # deletion composition (compartment / donor / state) per exclusion set
  comp_rows <- list(); donor_rem_rows <- list(); state_rem_rows <- list()
  for (es in exclusion_sets) {
    rem <- es$removed
    comp_rows[[length(comp_rows) + 1L]] <- data.frame(
      scope = es$scope, fraction = es$fraction, n_removed_total = length(rem),
      n_removed_axis = sum(is_axis[rem]),
      frac_axis_removed = round(sum(is_axis[rem]) / sum(is_axis), 4),
      n_removed_strict_non_gran = sum(is_ng[rem]),
      frac_strict_non_gran_removed = round(sum(is_ng[rem]) / sum(is_ng), 4),
      n_removed_other_gl = sum(is_ogl[rem]),
      frac_other_gl_removed = round(sum(is_ogl[rem]) / sum(is_ogl), 4),
      n_removed_in_shared_states = sum(is_shared[rem]),
      frac_shared_state_cells_removed = round(sum(is_shared[rem]) / sum(is_shared), 4),
      n_donors_affected = length(unique(donor_p[rem])),
      stringsAsFactors = FALSE)
    donor_rem_rows[[length(donor_rem_rows) + 1L]] <- data.frame(
      scope = es$scope, fraction = es$fraction, donor = acute_donors,
      n_primary = vapply(acute_donors, function(d) sum(donor_p == d), integer(1)),
      n_removed = vapply(acute_donors, function(d) sum(donor_p[rem] == d), integer(1)),
      stringsAsFactors = FALSE)
    donor_rem_rows[[length(donor_rem_rows)]]$frac_removed <-
      round(donor_rem_rows[[length(donor_rem_rows)]]$n_removed /
              donor_rem_rows[[length(donor_rem_rows)]]$n_primary, 4)
    state_rem_rows[[length(state_rem_rows) + 1L]] <- data.frame(
      scope = es$scope, fraction = es$fraction, state = obs_states,
      n_cells = vapply(obs_states, function(s) sum(state_p == s), integer(1)),
      n_removed = vapply(obs_states, function(s) sum(state_p[rem] == s), integer(1)),
      stringsAsFactors = FALSE)
    state_rem_rows[[length(state_rem_rows)]]$frac_removed <-
      round(state_rem_rows[[length(state_rem_rows)]]$n_removed /
              state_rem_rows[[length(state_rem_rows)]]$n_cells, 4)
  }
  write_csv_atomic(do.call(rbind, comp_rows), file.path(TABLE_DIR, "doublet_deletion_composition.csv"))
  write_csv_atomic(do.call(rbind, donor_rem_rows), file.path(TABLE_DIR, "doublet_deletion_by_donor.csv"))
  write_csv_atomic(do.call(rbind, state_rem_rows), file.path(TABLE_DIR, "doublet_deletion_by_state.csv"))
  
  # ---- 2.10 PART 1 -- STATE-LEVEL LOCALIZATION (four metrics) ---------------------
  state_contrib <- do.call(rbind, lapply(markers, function(g) {
    det <- as.numeric(Xm[g, ] > 0)
    state_n <- vapply(STRICT_NON_GRAN_STATES, function(s) sum(state_p == s), integer(1))
    state_detected <- vapply(STRICT_NON_GRAN_STATES, function(s) {
      idx <- which(state_p == s)
      if (!length(idx)) return(0L)
      as.integer(sum(det[idx] > 0))
    }, integer(1))
    within_frac <- ifelse(state_n > 0, state_detected / state_n, NA_real_)
    detected_share <- if (sum(state_detected) > 0)
      state_detected / sum(state_detected) else rep(NA_real_, length(state_detected))
    data.frame(marker = g, panel = unname(panel_of[g]),
               strict_non_gran_state = STRICT_NON_GRAN_STATES,
               state_n_cells = state_n,
               state_n_detected = state_detected,
               within_state_detection_fraction = within_frac,
               detected_cell_share = detected_share,
               rank_within_marker_by_fraction = rank(-within_frac, ties.method = "min",
                                                     na.last = "keep"),
               stringsAsFactors = FALSE)
  }))
  write_csv_atomic(state_contrib, file.path(TABLE_DIR, "state_level_contribution.csv"))
  
  # ---- 2.11 PART 2 -- SHARED-STATE DONOR COVERAGE ---------------------------------
  coverage <- data.frame(
    donor = acute_donors,
    n_eligible_shared_states = donor_cov$n_eligible_shared_states,
    eligible_shared_state_names = donor_cov$eligible_shared_state_names,
    strict_non_gran_cell_n = donor_cov$strict_non_gran_cell_n,
    in_correlation = donor_cov$n_eligible_shared_states >= MIN_SHARED_STATES_PER_DONOR,
    stringsAsFactors = FALSE)
  coverage$excluded_reason <- ifelse(coverage$in_correlation, "",
                                     "fewer than MIN_SHARED_STATES_PER_DONOR eligible shared states")
  write_csv_atomic(coverage, file.path(TABLE_DIR, "shared_state_donor_coverage.csv"))
  
  # ---- 2.12 PART 3 -- LEAVE-ONE-DONOR-OUT JACKKNIFE (20-donor standard) -----------
  jack_rows <- list(); jack_summary <- list()
  for (pn in FORMAL_PANELS) {
    gset <- intersect(AMBIENT_PANELS[[pn]], markers)
    sig <- compute_donor_signals(seq_len(n_primary), gset)
    eval_mask <- is.finite(nf_full) & is.finite(sig)
    n_eval <- sum(eval_mask)
    if (n_eval < MIN_EVALUABLE_DONORS) {
      jack_summary[[pn]] <- data.frame(
        panel = pn, n_donors = n_eval, rho_full = NA_real_,
        rho_min_loo = NA_real_, rho_max_loo = NA_real_,
        n_loo_below_050 = NA_integer_, most_influential_donor = NA_character_,
        max_delta_rho = NA_real_, verdict = "not_evaluable_too_few_donors",
        stringsAsFactors = FALSE)
      next
    }
    rho_full <- spearman_rho(nf_full[eval_mask], sig[eval_mask])
    idx_e <- which(eval_mask)
    rlo <- numeric(length(idx_e)); delta <- numeric(length(idx_e))
    for (i in seq_along(idx_e)) {
      keep_ <- idx_e[-i]
      rlo[i] <- spearman_rho(nf_full[keep_], sig[keep_])
      delta[i] <- rho_full - rlo[i]
    }
    drop_below <- !is.na(rlo) & rlo < STRONG_RHO_THRESHOLD
    verdict <- if (is.na(rho_full)) "not_evaluable"
    else if (all(is.na(rlo))) "not_evaluable_leave_one_out"
    else if (rho_full >= STRONG_RHO_THRESHOLD) {
      if (any(drop_below)) "sensitive_to_single_donor" else "not_single_donor_driven"
    } else "structure_not_present_full_data"
    # all-NA protection for leave-one-out aggregates
    finite_loo <- is.finite(rlo)
    if (!any(finite_loo)) {
      rho_min_loo <- NA_real_; rho_max_loo <- NA_real_
      max_delta_rho <- NA_real_; influential_donor <- NA_character_
    } else {
      rho_min_loo <- min(rlo[finite_loo])
      rho_max_loo <- max(rlo[finite_loo])
      finite_delta <- is.finite(delta)
      if (any(finite_delta)) {
        jmax <- which(finite_delta)[which.max(delta[finite_delta])]
        max_delta_rho <- delta[jmax]
        influential_donor <- acute_donors[idx_e][jmax]
      } else {
        max_delta_rho <- NA_real_; influential_donor <- NA_character_
      }
    }
    jack_rows[[pn]] <- data.frame(
      panel = pn, donor = acute_donors[idx_e],
      rho_full = rho_full, rho_leave_one_out = rlo, delta_rho = delta,
      leave_one_out_below_050 = drop_below, stringsAsFactors = FALSE)
    jack_summary[[pn]] <- data.frame(
      panel = pn, n_donors = n_eval, rho_full = rho_full,
      rho_min_loo = rho_min_loo, rho_max_loo = rho_max_loo,
      n_loo_below_050 = sum(drop_below),
      most_influential_donor = influential_donor,
      max_delta_rho = max_delta_rho, verdict = verdict,
      stringsAsFactors = FALSE)
  }
  jack_df <- do.call(rbind, jack_rows)
  write_csv_atomic(jack_df, file.path(TABLE_DIR, "donor_jackknife.csv"))
  jack_sum_df <- do.call(rbind, jack_summary)
  write_csv_atomic(jack_sum_df, file.path(TABLE_DIR, "donor_jackknife_summary.csv"))
  
  # ---- 2.13 PART 4 -- DONOR-LEVEL AMBIENT RHO (DUAL-BRANCH) -------------------------
  keep_sets <- list(full = seq_len(n_primary))
  for (es in exclusion_sets) keep_sets[[paste0(es$scope, "_top", round(100 * es$fraction))]] <- es$keep
  
  donor_signals_long <- do.call(rbind, lapply(FORMAL_PANELS, function(pn) {
    gset <- intersect(AMBIENT_PANELS[[pn]], markers)
    do.call(rbind, lapply(names(keep_sets), function(ksn) {
      sig <- compute_donor_signals(keep_sets[[ksn]], gset)
      n_elig <- vapply(acute_donors, function(d) {
        st <- ng_state_cells[[d]]
        sum(vapply(shared_states, function(s) {
          j <- match(s, STRICT_NON_GRAN_STATES)
          length(intersect(st[[j]], keep_sets[[ksn]])) >= MIN_CELLS_PER_DONOR_STATE
        }, logical(1)))
      }, integer(1))
      scope <- if (ksn == "full") "full" else sub("_top[0-9]+$", "", ksn)
      fraction <- if (ksn == "full") NA_real_
      else as.numeric(sub("^.*_top([0-9]+)$", "\\1", ksn)) / 100
      data.frame(donor = acute_donors, panel = pn, scope = scope,
                 fraction = fraction,
                 neutrophil_fraction_full = nf_full,
                 ambient_signal = sig,
                 n_eligible_shared_states_after_exclusion = n_elig,
                 in_correlation = n_elig >= MIN_SHARED_STATES_PER_DONOR,
                 stringsAsFactors = FALSE)
    }))
  }))
  write_csv_atomic(donor_signals_long,
                   file.path(TABLE_DIR, "ambient_doublet_sensitivity_donor_signals_long.csv"))
  
  dt_rows <- lapply(FORMAL_PANELS, function(pn) {
    gset <- intersect(AMBIENT_PANELS[[pn]], markers)
    sig_full <- compute_donor_signals(seq_len(n_primary), gset)
    ok_full <- is.finite(nf_full) & is.finite(sig_full)
    rho_full <- spearman_rho(nf_full[ok_full], sig_full[ok_full])
    row <- list(panel = pn, n_donors_full = sum(ok_full), rho_full = rho_full)
    for (ksn in names(keep_sets)[-1]) {
      sig_k <- compute_donor_signals(keep_sets[[ksn]], gset)
      ok_k <- ok_full & is.finite(sig_k)
      row[[paste0("n_donors_", ksn)]] <- sum(ok_k)
      row[[paste0("rho_", ksn)]] <- spearman_rho(nf_full[ok_k], sig_k[ok_k])
    }
    r_ps <- row$rho_per_sample_top10
    r_gl <- row$rho_global_acute_primary_top10
    verdict <- if (is.na(rho_full) || rho_full < STRONG_RHO_THRESHOLD)
      "no_strong_positive_structure_at_full"
    else if (is.na(r_ps) || is.na(r_gl)) "not_evaluable_too_few_donors"
    else {
      s_ps <- r_ps >= STRONG_RHO_THRESHOLD
      s_gl <- r_gl >= STRONG_RHO_THRESHOLD
      if (s_ps && s_gl) "doublet_risk_robust"
      else if (s_ps && !s_gl) "global_scope_affected_primary_passes"
      else if (!s_ps && s_gl) "primary_fails_not_robust"
      else "doublet_risk_sensitive"
    }
    row$verdict_top10 <- verdict
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  dt_df <- do.call(rbind, dt_rows)
  write_csv_atomic(dt_df, file.path(TABLE_DIR, "ambient_doublet_sensitivity_panel_summary.csv"))
  
  # ---- 2.14 PART 5 -- M1 COUNT-BASED CONTRIBUTION (DUAL-BRANCH) + PAIRED BOOTSTRAP --
  # M1 genes: FROZEN set is the only entry (validated below; no manifest override)
  if (length(M1_GENES) != 11L || anyDuplicated(M1_GENES))
    stop_msg("M1_GENES must contain exactly 11 unique genes.")
  missing_m1 <- setdiff(M1_GENES, rownames(cnt))
  if (length(missing_m1))
    stop_msg("M1 genes absent from counts: ", paste(missing_m1, collapse = ", "))
  
  PROG_STATES <- M1_COMPARTMENTS$progenitor_granulopoiesis
  IMM_STATES  <- M1_COMPARTMENTS$immature_neutrophil
  CYC_STATES  <- M1_COMPARTMENTS$cycling_lymphoid
  named_states <- unique(c(PROG_STATES, IMM_STATES, CYC_STATES))
  miss_st <- setdiff(named_states, obs_states)
  if (length(miss_st))
    stop_msg("M1 compartment states absent from obs states: ",
             paste(miss_st, collapse = ", "))
  RESID_STATES <- setdiff(obs_states, named_states)
  if (length(RESID_STATES) != EXPECTED_N_FINE_STATES - length(named_states))
    stop_msg("broad_residual compartment size unexpected (expected ",
             EXPECTED_N_FINE_STATES - length(named_states), " states).")
  
  Xm1 <- cnt[match(M1_GENES, rownames(cnt)), primary_idx, drop = FALSE]
  m1_cell_umi <- as.numeric(Matrix::colSums(Xm1))
  
  compute_m1_group <- function(keep_idx) {
    do.call(rbind, lapply(acute_donors, function(d) {
      idxd <- intersect(which(donor_p == d), keep_idx)
      state_umi <- vapply(obs_states, function(s) {
        idxs <- idxd[state_p[idxd] == s]
        if (!length(idxs)) return(0)
        sum(m1_cell_umi[idxs])
      }, numeric(1))
      total <- sum(state_umi)
      if (total <= 0)
        return(data.frame(donor = d, prog = NA_real_, resid = NA_real_,
                          immature = NA_real_, cycling = NA_real_, total_umi = 0))
      data.frame(donor = d,
                 prog     = sum(state_umi[PROG_STATES]) / total,
                 resid    = sum(state_umi[RESID_STATES]) / total,
                 immature = sum(state_umi[IMM_STATES]) / total,
                 cycling  = sum(state_umi[CYC_STATES]) / total,
                 total_umi = total)
    }))
  }
  
  m1_full <- compute_m1_group(seq_len(n_primary))
  med_prog_full  <- median(m1_full$prog,  na.rm = TRUE)
  med_resid_full <- median(m1_full$resid, na.rm = TRUE)
  cat(sprintf("[Part5] M1 FULL baseline: progenitor=%.6f broad_residual=%.6f\n",
              med_prog_full, med_resid_full))
  if (abs(med_prog_full - M1_BASELINE[["progenitor_granulopoiesis"]]) > M1_BASELINE_TOL ||
      abs(med_resid_full - M1_BASELINE[["broad_residual"]]) > M1_BASELINE_TOL)
    stop_msg(sprintf(paste0("M1 baseline mismatch (tol=%.0e): got progenitor=%.6f (expect %.6f), ",
                            "broad_residual=%.6f (expect %.6f). M1 gene list / state names / ",
                            "compartments differ from pre-audit Script 02/S17."),
                     M1_BASELINE_TOL, med_prog_full, M1_BASELINE[["progenitor_granulopoiesis"]],
                     med_resid_full, M1_BASELINE[["broad_residual"]]))
  d_full <- med_prog_full - med_resid_full
  
  m1_groups <- list(full = m1_full)
  for (es in exclusion_sets)
    m1_groups[[paste0(es$scope, "_top", round(100 * es$fraction))]] <- compute_m1_group(es$keep)
  
  # ---- PAIRED bootstrap: sample ONLY donors commonly evaluable in ALL groups -------
  eval_in_all <- Reduce(`&`, lapply(m1_groups, function(g)
    is.finite(g$prog) & is.finite(g$resid) &
      is.finite(g$immature) & is.finite(g$cycling)))
  common_idx <- which(eval_in_all)
  n_eval_all <- length(common_idx)
  paired_ok <- n_eval_all >= MIN_EVALUABLE_DONORS
  if (paired_ok) {
    set.seed(M1_BOOT_SEED)
    boot_index <- replicate(M1_BOOT_REPS,
                            sample(common_idx, size = n_eval_all, replace = TRUE))
  } else {
    boot_index <- NULL
    cat("[Part5] Paired bootstrap not evaluable: ", n_eval_all,
        " common donors; requires at least ", MIN_EVALUABLE_DONORS,
        ". No CI/support will be claimed for any group.\n", sep = "")
  }
  
  # Common-donor point estimates (the estimand the paired CI belongs to);
  # all-evaluable estimates remain separate and always reported.
  d_full_common <- if (paired_ok) {
    median(m1_groups$full$prog[common_idx]) -
      median(m1_groups$full$resid[common_idx])
  } else {
    NA_real_
  }
  
  top_names <- c("progenitor_granulopoiesis","broad_residual",
                 "immature_neutrophil","cycling_lymphoid")
  rank_of_group <- function(gdat, idx = seq_len(nrow(gdat))) {
    meds <- c(progenitor_granulopoiesis = median(gdat$prog[idx], na.rm = TRUE),
              broad_residual = median(gdat$resid[idx], na.rm = TRUE),
              immature_neutrophil = median(gdat$immature[idx], na.rm = TRUE),
              cycling_lymphoid = median(gdat$cycling[idx], na.rm = TRUE))
    paste(names(sort(meds, decreasing = TRUE)), collapse = ">")
  }
  full_rank_all    <- rank_of_group(m1_groups$full)
  full_rank_common <- if (paired_ok) rank_of_group(m1_groups$full, common_idx) else NA_character_
  
  m1_summary <- list(); top_summary <- list(); full_deltas <- NULL
  for (gname in names(m1_groups)) {
    gdat <- m1_groups[[gname]]
    n_cells_donors <- sum(gdat$total_umi > 0)
    ok <- is.finite(gdat$prog) & is.finite(gdat$resid)
    n_eval <- sum(ok)
    uneval <- acute_donors[!ok]
    d_med_all <- median(gdat$prog, na.rm = TRUE) - median(gdat$resid, na.rm = TRUE)
    d_med_common <- if (paired_ok) {
      median(gdat$prog[common_idx]) - median(gdat$resid[common_idx])
    } else {
      NA_real_
    }
    boot_deltas <- NULL
    ci <- c(NA_real_, NA_real_); support_gt0 <- NA_real_; support_frac <- NA_real_
    if (paired_ok) {
      boot_deltas <- apply(boot_index, 2, function(idx)
        median(gdat$prog[idx]) - median(gdat$resid[idx]))
      ci <- unname(quantile(boot_deltas, c(0.025, 0.975)))
      support_gt0 <- mean(boot_deltas > 0)
      support_frac <- mean(boot_deltas > 0) + 0.5 * mean(boot_deltas == 0)
    }
    if (gname == "full" && paired_ok) full_deltas <- boot_deltas
    top_support <- setNames(rep(NA_real_, 4), top_names)
    if (paired_ok) {
      wins <- setNames(numeric(4), top_names)
      for (b in seq_len(M1_BOOT_REPS)) {
        idx <- boot_index[, b]
        meds <- c(median(gdat$prog[idx]), median(gdat$resid[idx]),
                  median(gdat$immature[idx]), median(gdat$cycling[idx]))
        top <- which(meds == max(meds))
        wins[top] <- wins[top] + 1 / length(top)
      }
      top_support <- wins / M1_BOOT_REPS
    }
    rank_all    <- rank_of_group(gdat)
    rank_common <- if (paired_ok) rank_of_group(gdat, common_idx) else NA_character_
    row <- data.frame(
      group = gname,
      n_donors_with_cells = n_cells_donors,
      n_donors_evaluable_M1 = n_eval,
      donors_unevaluable = paste(uneval, collapse = ";"),
      progenitor_median = median(gdat$prog,  na.rm = TRUE),
      broad_residual_median = median(gdat$resid, na.rm = TRUE),
      immature_median = median(gdat$immature, na.rm = TRUE),
      cycling_lymphoid_median = median(gdat$cycling, na.rm = TRUE),
      diff_prog_minus_resid_all_evaluable = d_med_all,
      diff_prog_minus_resid_common_donors = d_med_common,
      ci_common_lo = ci[1], ci_common_hi = ci[2],
      support_gt0 = support_gt0,
      support_fractional_ties = support_frac,
      rank_order_all_evaluable = rank_all,
      rank_order_common_donors = rank_common,
      stringsAsFactors = FALSE)
    row$d_delta_vs_full_all_evaluable <- if (gname == "full") NA_real_ else d_med_all - d_full
    row$d_delta_vs_full_common_donors <- if (gname == "full" || !paired_ok) NA_real_
    else d_med_common - d_full_common
    row$paired_diff_median  <- NA_real_
    row$paired_diff_ci_lo   <- NA_real_
    row$paired_diff_ci_hi   <- NA_real_
    row$paired_diff_contains_zero <- NA
    row$exploratory_5pp_rule <- NA
    row$rank_changed_vs_full <- if (gname == "full") FALSE
    else if (!paired_ok) NA
    else !identical(rank_common, full_rank_common)
    if (gname != "full" && paired_ok) {
      dd <- boot_deltas - full_deltas
      pci <- unname(quantile(dd, c(0.025, 0.975)))
      row$paired_diff_median <- median(dd)
      row$paired_diff_ci_lo <- pci[1]
      row$paired_diff_ci_hi <- pci[2]
      row$paired_diff_contains_zero <- (pci[1] <= 0) & (pci[2] >= 0)
      row$exploratory_5pp_rule <- abs(row$d_delta_vs_full_common_donors) <=
        M1_EXPLORATORY_STABILITY_BAND
    }
    m1_summary[[gname]] <- row
    top_summary[[gname]] <- data.frame(group = gname, t(top_support),
                                       check.names = FALSE, stringsAsFactors = FALSE)
  }
  m1_sum_df <- do.call(rbind, m1_summary)
  write_csv_atomic(m1_sum_df, file.path(TABLE_DIR, "m1_contribution_groups.csv"))
  write_csv_atomic(do.call(rbind, top_summary), file.path(TABLE_DIR, "m1_top_compartment_support.csv"))
  
  # per-donor x state cell counts (full + each exclusion set)
  full_counts <- do.call(rbind, lapply(acute_donors, function(d) {
    idxd <- which(donor_p == d)
    data.frame(donor = d, state = obs_states,
               n_full = vapply(obs_states, function(s) sum(state_p[idxd] == s), integer(1)),
               stringsAsFactors = FALSE)
  }))
  for (gname in names(m1_groups)[-1]) {
    keep_idx <- exclusion_sets[[gname]]$keep
    full_counts[[gname]] <- unlist(lapply(acute_donors, function(d) {
      idxd <- intersect(which(donor_p == d), keep_idx)
      vapply(obs_states, function(s) sum(state_p[idxd] == s), integer(1))
    }), use.names = FALSE)
  }
  write_csv_atomic(full_counts, file.path(TABLE_DIR, "m1_donor_state_cell_counts.csv"))
  
  # descriptive M1 summary (NOT a formal equivalence verdict)
  ps10 <- m1_sum_df[m1_sum_df$group == "per_sample_top10", , drop = FALSE]
  gl10 <- m1_sum_df[m1_sum_df$group == "global_acute_primary_top10", , drop = FALSE]
  m1_verdict <- "not_evaluable"
  if (paired_ok && nrow(ps10) == 1L && nrow(gl10) == 1L) {
    z_ps <- isTRUE(ps10$paired_diff_contains_zero)
    z_gl <- isTRUE(gl10$paired_diff_contains_zero)
    if (z_ps && z_gl) m1_verdict <- "descriptive_paired_ci_contains_zero"
    else if (z_ps && !z_gl) m1_verdict <- "descriptive_global_branch_paired_ci_excludes_zero"
    else if (!z_ps && z_gl) m1_verdict <- "descriptive_per_sample_branch_paired_ci_excludes_zero"
    else m1_verdict <- "descriptive_both_branches_paired_ci_exclude_zero"
  }
  
  # ---- 2.15 CROSS-CHECK vs QC RUN (preflight objects reused; RHO_TOL = 1e-6) --------
  # (QC_RUN_DIR, STATUS, version, purpose, formal_run, hash-algo, RDS hash,
  #  shared states already checked in preflight -- no provenance re-read here.)
  panel_csv <- file.path(QC_RUN_DIR, "tables", "ambient_vs_neutrophil_fraction_by_panel.csv")
  if (!file.exists(panel_csv)) stop_msg("QC run panel CSV not found: ", panel_csv)
  qc_tab <- read.csv(panel_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("panel" %in% names(qc_tab))) stop_msg("QC panel summary missing 'panel' column.")
  if (anyDuplicated(qc_tab$panel)) stop_msg("Duplicated panel rows in QC panel summary.")
  if (length(setdiff(FORMAL_PANELS, qc_tab$panel)))
    stop_msg("QC run missing formal panels: ",
             paste(setdiff(FORMAL_PANELS, qc_tab$panel), collapse = ", "))
  deltas <- numeric(length(FORMAL_PANELS)); names(deltas) <- FORMAL_PANELS
  for (pn in FORMAL_PANELS) {
    qr <- qc_tab[qc_tab$panel == pn, , drop = FALSE]
    dr <- dt_df[dt_df$panel == pn, , drop = FALSE]
    if (nrow(qr) != 1L || nrow(dr) != 1L)
      stop_msg("Expected exactly one row for panel: ", pn)
    if (is.na(qr$rho[1]) || is.na(dr$rho_full[1])) stop_msg("Panel rho NA in cross-check: ", pn)
    if ("n_donors_with_data" %in% names(qr) &&
        !is.na(qr$n_donors_with_data[1]) &&
        !identical(as.integer(qr$n_donors_with_data[1]), as.integer(dr$n_donors_full[1])))
      stop_msg("Donor-count mismatch for panel: ", pn)
    deltas[pn] <- abs(qr$rho[1] - dr$rho_full[1])
  }
  if (anyNA(deltas)) stop_msg("Cross-check requires all 4 formal panels.")
  if (max(deltas) > RHO_TOL)
    stop_msg(sprintf("Reproducibility MISMATCH vs QC run (tol=%.0e): max |drho| = %.7f",
                     RHO_TOL, max(deltas)))
  
  # evaluable donor IDENTITY comparison (not just donor counts)
  qc_donor_path <- file.path(QC_RUN_DIR, "tables",
                             "ambient_vs_neutrophil_fraction_by_panel_donor.csv")
  if (!file.exists(qc_donor_path))
    stop_msg("QC donor-level panel table not found: ", qc_donor_path)
  qc_donor <- read.csv(qc_donor_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_qc_donor_cols <- c("donor","panel","neutrophil_fraction",
                              "strict_non_gran_detection_fraction")
  miss_cols <- setdiff(required_qc_donor_cols, names(qc_donor))
  if (length(miss_cols))
    stop_msg("QC donor table missing columns: ", paste(miss_cols, collapse = ", "))
  current_full <- donor_signals_long[
    donor_signals_long$scope == "full" & is.na(donor_signals_long$fraction), , drop = FALSE]
  for (pn in FORMAL_PANELS) {
    q <- qc_donor[qc_donor$panel == pn, , drop = FALSE]
    d <- current_full[current_full$panel == pn, , drop = FALSE]
    q_ok <- is.finite(q$neutrophil_fraction) &
      is.finite(q$strict_non_gran_detection_fraction)
    d_ok <- is.finite(d$neutrophil_fraction_full) & is.finite(d$ambient_signal)
    q_ids <- sort(as.character(q$donor[q_ok]))
    d_ids <- sort(as.character(d$donor[d_ok]))
    if (!identical(q_ids, d_ids))
      stop_msg("Evaluable donor-ID mismatch for panel: ", pn)
  }
  cross_check <- "verified"
  cat("[cross-check] version=", qc_version, "; all 4 panels match QC run ",
      "(max |drho| = ", format(max(deltas), digits = 7), "); donor IDs verified.\n", sep = "")
  
  # ---- 2.16 PROVENANCE -------------------------------------------------------------
  # script_path and script_hash were computed in the preflight (Rscript required)
  prov_rows <- list()
  add_prov <- function(k, v) prov_rows[[length(prov_rows) + 1L]] <<- data.frame(
    key = k, value = as.character(v), stringsAsFactors = FALSE)
  add_prov("script", "qc_ambient_followup.R")
  add_prov("version", SCRIPT_VERSION)
  add_prov("formal_run", FORMAL_RUN)
  add_prov("run_purpose", RUN_PURPOSE)
  add_prov("output_root", OUT_ROOT)
  add_prov("qc_run_dir", QC_RUN_DIR)
  add_prov("qc_expected_version", EXPECTED_QC_VERSION)
  add_prov("qc_observed_version", qc_version)
  add_prov("qc_expected_purpose", EXPECTED_QC_PURPOSE)
  add_prov("qc_observed_purpose", qc_purpose)
  add_prov("qc_expected_formal_run", EXPECTED_QC_FORMAL_RUN)
  add_prov("qc_observed_formal_run", qc_formal_run)
  add_prov("qc_hash_algo_used", qc_hash_algo)
  add_prov("hash_algo_used", EXPECTED_QC_HASH_ALGO)
  add_prov("input_rds", INPUT_RDS)
  add_prov("input_rds_hash", input_rds_hash)
  add_prov("script_path", script_path)
  add_prov("script_hash", script_hash)
  add_prov("doublet_risk_field", DOUBLET_RISK_FIELD)
  add_prov("doublet_primary_scope", DOUBLET_PRIMARY_SCOPE)
  add_prov("doublet_sensitivity_scopes", paste(DOUBLET_SENSITIVITY_SCOPES, collapse = ";"))
  add_prov("exact_k_rule", EXACT_K_RULE)
  add_prov("tie_break_rule", TIE_BREAK_RULE)
  add_prov("expected_exclusion_n", paste(names(EXPECTED_EXCLUSION_N),
                                         EXPECTED_EXCLUSION_N, sep = "=", collapse = ";"))
  add_prov("per_sample_equals_per_donor", "TRUE (26 acute samples == 26 acute donors)")
  add_prov("predicted_doublet", "non-missing all-FALSE in primary; dead/unusable; not used")
  add_prov("batch_R2_score_ranks", if (is.na(batch_R2)) "NA" else round(batch_R2, 6))
  for (es in exclusion_sets) {
    lbl <- paste0(es$scope, "_top", round(100 * es$fraction))
    cut_r <- range(es$cutoffs)
    add_prov(paste0("exclusion_", lbl),
             sprintf("removed=%d cells; cutoff range %.6f-%.6f", es$n_removed,
                     cut_r[1], cut_r[2]))
  }
  add_prov("m1_readout", M1_READOUT)
  add_prov("m1_ucell_available", M1_UCELL_AVAILABLE)
  add_prov("frozen_m1_genes", paste(FROZEN_M1_GENES, collapse = ";"))
  add_prov("m1_genes", paste(M1_GENES, collapse = ";"))
  add_prov("m1_compartments", paste(vapply(M1_COMPARTMENTS, paste, character(1), collapse = ";"),
                                    collapse = " | "))
  add_prov("m1_broad_residual_states", paste(RESID_STATES, collapse = ";"))
  add_prov("m1_baseline_full",
           sprintf("progenitor=%.6f broad_residual=%.6f", med_prog_full, med_resid_full))
  add_prov("m1_baseline_tol", M1_BASELINE_TOL)
  add_prov("m1_boot_reps", M1_BOOT_REPS)
  add_prov("m1_boot_seed", M1_BOOT_SEED)
  add_prov("m1_paired_common_donors", n_eval_all)
  add_prov("m1_paired_common_donor_ids", paste(acute_donors[common_idx], collapse = ";"))
  add_prov("m1_paired_ok", paired_ok)
  add_prov("m1_exploratory_stability_band",
           sprintf("%s (exploratory_5pp_rule; NOT a formal equivalence threshold)",
                   M1_EXPLORATORY_STABILITY_BAND))
  add_prov("m1_verdict", m1_verdict)
  for (gname in names(m1_summary)) {
    r <- m1_sum_df[m1_sum_df$group == gname, , drop = FALSE]
    if (nrow(r))
      add_prov(paste0("m1_", gname),
               sprintf(paste0("prog=%.6f resid=%.6f diff_all=%.6f diff_common=%.6f ",
                              "ci_common=[%.6f,%.6f] support_gt0=%.4f support_frac=%.4f n_eval=%d | ",
                              "paired_diff_median=%.6f paired_ci=[%.6f,%.6f] zero_in=%s | ",
                              "explor_5pp=%s rank_changed=%s"),
                       r$progenitor_median[1], r$broad_residual_median[1],
                       r$diff_prog_minus_resid_all_evaluable[1],
                       r$diff_prog_minus_resid_common_donors[1],
                       r$ci_common_lo[1], r$ci_common_hi[1],
                       r$support_gt0[1], r$support_fractional_ties[1],
                       r$n_donors_evaluable_M1[1],
                       r$paired_diff_median[1], r$paired_diff_ci_lo[1],
                       r$paired_diff_ci_hi[1],
                       as.character(r$paired_diff_contains_zero[1]),
                       as.character(r$exploratory_5pp_rule[1]),
                       as.character(r$rank_changed_vs_full[1])))
  }
  add_prov("fixed_shared_states", paste(shared_states, collapse = ";"))
  add_prov("min_cells_per_donor_state", MIN_CELLS_PER_DONOR_STATE)
  add_prov("min_shared_states_per_donor", MIN_SHARED_STATES_PER_DONOR)
  add_prov("min_shared_state_donor_n", MIN_SHARED_STATE_DONOR_N)
  add_prov("min_evaluable_donors", MIN_EVALUABLE_DONORS)
  add_prov("strong_rho_threshold", STRONG_RHO_THRESHOLD)
  add_prov("rho_tol_cross_check", RHO_TOL)
  add_prov("shared_state_donor_n_evaluable", sum(coverage$in_correlation))
  for (pn in FORMAL_PANELS) {
    j <- jack_sum_df[jack_sum_df$panel == pn, , drop = FALSE]
    if (nrow(j))
      add_prov(paste0("jackknife_", pn),
               sprintf("%s (rho_full=%.4f, rho_loo %.4f-%.4f, n_loo<0.50=%d, influential=%s)",
                       j$verdict[1], j$rho_full[1], j$rho_min_loo[1], j$rho_max_loo[1],
                       j$n_loo_below_050[1], j$most_influential_donor[1]))
    d <- dt_df[dt_df$panel == pn, , drop = FALSE]
    if (nrow(d))
      add_prov(paste0("doublet_sens_", pn),
               sprintf(paste0("full=%.6f(n=%d) | ps5=%.6f ps10=%.6f | gl5=%.6f gl10=%.6f | ",
                              "verdict=%s"),
                       d$rho_full[1], d$n_donors_full[1],
                       d$rho_per_sample_top5[1], d$rho_per_sample_top10[1],
                       d$rho_global_acute_primary_top5[1], d$rho_global_acute_primary_top10[1],
                       d$verdict_top10[1]))
  }
  add_prov("cross_check_vs_qc_run", cross_check)
  add_prov("notes", paste0(
    "Dual-branch doublet sensitivity. per_sample is the PRIMARY scope ",
    "(sample-balanced; per-sample cutoffs and per-sample share of the global ",
    "top5 0.4-13.2% show cross-sample scale incomparability); ",
    "global_acute_primary is the legacy/manuscript-replication branch. ",
    "Exact-k removal with ceiling() and order(-score, cell_barcode, radix); ",
    "strict totals verified 7607/15200/7592/15184. M1 readout = count-UMI ",
    "contribution fraction (11 frozen genes), not UCell. Bootstrap is PAIRED ",
    "(single donor-resample matrix over donors commonly evaluable in all ",
    "groups); point estimates are reported both all-evaluable and common-donor, ",
    "and the paired CI always corresponds to the common-donor estimate; ",
    "'P' is a bootstrap SUPPORT PROBABILITY. paired_diff_ci_contains_zero and ",
    "exploratory_5pp_rule are descriptive statements, NOT a formal ",
    "equivalence test. doublet_score is a doublet-RISK ranking score, not a ",
    "validated doublet label. No cells were removed from any primary analysis; ",
    "these are descriptive robustness assessments (non-correction declaration)."))
  add_prov("R_version", R.version.string)
  prov <- do.call(rbind, prov_rows)
  write_csv_atomic(prov, file.path(PROV_DIR, "run_provenance.csv"))
  writeLines(capture.output(sessionInfo()), file.path(PROV_DIR, "sessionInfo.txt"))
  writeLines(paste("run_dir", OUT_ROOT, sep = "\t"), file.path(PROV_DIR, "run_dir.txt"))
  
  # ---- 2.17 CONSOLE SUMMARY ----------------------------------------------------------
  cat("\n========== Ambient follow-up ", SCRIPT_VERSION, " (", RUN_PURPOSE, ") ==========\n", sep = "")
  cat("\n-- QC preflight: ", QC_RUN_DIR, "\n", sep = "")
  cat("  version=", qc_version, " purpose=", qc_purpose,
      " formal_run=", qc_formal_run, " hash_algo=", qc_hash_algo, "\n", sep = "")
  cat("  script=", script_path, "\n", sep = "")
  cat("\n-- Part 0: doublet-risk (primary scope = ", DOUBLET_PRIMARY_SCOPE, ")\n", sep = "")
  cat("  predicted_doublet: non-missing, all-FALSE (recorded only)\n")
  cat("  batch R2 (score ranks): ", if (is.na(batch_R2)) "NA" else round(batch_R2, 4), "\n", sep = "")
  cat("  exclusions (exact-k verified):\n")
  for (es in exclusion_sets)
    cat(sprintf("    %s_top%d: removed=%d cells\n", es$scope, round(100 * es$fraction), es$n_removed))
  
  cat("\n-- Part 1: state-level localization (top states by within-state detection fraction)\n")
  for (g in c("MPO","RNASE3","ELANE","LTF","CEACAM8")) {
    if (!(g %in% markers)) { cat("  ", g, ": absent\n", sep = ""); next }
    sub <- state_contrib[state_contrib$marker == g, ]
    sub <- sub[order(sub$within_state_detection_fraction, decreasing = TRUE), ]
    top <- head(sub, 3)
    cat("  ", g, ": ", sep = "")
    cat(paste(sprintf("%s frac=%.4f n=%d share=%.1f%%",
                      top$strict_non_gran_state, top$within_state_detection_fraction,
                      top$state_n_detected, 100 * top$detected_cell_share), collapse = "; "), "\n")
  }
  
  cat("\n-- Part 2: shared-state donor coverage\n")
  cat("  shared states: ", length(shared_states), "/", length(STRICT_NON_GRAN_STATES),
      " | evaluable donors: ", sum(coverage$in_correlation), "/", length(acute_donors), "\n", sep = "")
  
  cat("\n-- Part 3: jackknife (leave-one-donor-out, 20-donor standard)\n")
  print(jack_sum_df, row.names = FALSE)
  
  cat("\n-- Part 4: donor-level ambient rho (dual-branch)\n")
  print(dt_df, row.names = FALSE)
  cat("  donor-level signals: tables/ambient_doublet_sensitivity_donor_signals_long.csv\n")
  
  cat("\n-- Part 5: M1 count-UMI contribution (dual-branch, PAIRED bootstrap)\n")
  cat("  paired common donors: ", n_eval_all, "/", length(acute_donors),
      " (paired_ok = ", paired_ok, ")\n", sep = "")
  if (paired_ok) cat("  common donor IDs: ", paste(acute_donors[common_idx], collapse = ", "), "\n", sep = "")
  print(m1_sum_df, row.names = FALSE)
  cat("  top-compartment support:\n")
  print(do.call(rbind, top_summary), row.names = FALSE)
  cat("  M1 descriptive summary: ", m1_verdict,
      " (diff_*_all_evaluable and diff_*_common_donors both reported; paired CI ",
      "corresponds to the common-donor estimate. paired_diff_ci_contains_zero / ",
      "exploratory_5pp_rule are DESCRIPTIVE, not a formal equivalence test)\n", sep = "")
  
  cat("\n-- cross-check vs QC run: ", cross_check, " (version ", qc_version, ")\n", sep = "")
  writeLines("COMPLETE", STATUS_FILE)
  run_status <- "COMPLETE"
  cat("\n[complete] Follow-up ", SCRIPT_VERSION, " complete.\nRun dir: ", OUT_ROOT, "\n", sep = "")
  invisible(TRUE)
}   # end main()

# ---- 3. RUN -----------------------------------------------------------------
main()