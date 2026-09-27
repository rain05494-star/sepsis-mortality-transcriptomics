# ============================================================================
# SCRIPT 01 v4 -- GSE216009 32-gene mapping + detectability audit
#
# Key fixes:
#   [1] Module assignment source of truth = Giannini E3 tables; tmod cross-checks
#   [2] 24 real fine_annot states mapped to display order/color/lineage
#   [3] pct_state uses detected cells: (X > 0) %*% t(S)
#   [4] donor/sample counting is sample-level, not cell-level
#   [5] donor_detection_rate denominator = number of donors
#   [6] count/metadata/cell-barcode alignment enforced
#   [7] gene mapping is fail-closed and auditable
#
# RUN:
#   Rscript 01_gene_mapping_detectability.R
# ============================================================================

options(stringsAsFactors = FALSE)

local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  PROJECT_ROOT <- "/home/sunshine/predicate/singlecell"
  RESULTS_META <- "/home/sunshine/predicate/test03/results_meta"
  INPUT_RDS    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  
  TABLE_32GENES <- file.path(
    RESULTS_META, "02_meta_analysis", "meta_significant_32_genes.csv"
  )
  TABLE_MODULES <- file.path(
    RESULTS_META, "09_MARS_module_interaction", "Table_MARS_module_gene_membership.csv"
  )
  TABLE_E3_OVERLAP <- file.path(
    RESULTS_META, "08_Giannini_sensitivity", "Table_Giannini_E3_overlap_direction.csv"
  )
  TABLE_E3_RESIDUAL <- file.path(
    RESULTS_META, "08_Giannini_sensitivity",
    "Table_Giannini_E3_nonoverlapping_signature_genes.csv"
  )
  
  # ---- meaningful output naming: {stage}/{YYYYMMDD}_{purpose} ----
  OUTPUTS_ROOT <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  STAGE_NAME   <- "01_gene_mapping_detectability"
  RUN_PURPOSE  <- "gene_mapping_detectability"
  
  make_run_dir <- function(root, stage, purpose) {
    stage_dir <- file.path(root, stage)
    dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
    base <- paste0(format(Sys.Date(), "%Y%m%d"), "_", purpose)
    run_dir <- file.path(stage_dir, base)
    if (dir.exists(run_dir)) {
      i <- 2L
      repeat {
        cand <- file.path(stage_dir, paste0(base, "_", i))
        if (!dir.exists(cand)) { run_dir <- cand; break }
        i <- i + 1L
      }
    }
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    run_dir
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  PROV_DIR  <- file.path(OUT_ROOT, "provenance")
  for (d in c(TABLE_DIR, FIG_DIR, PROV_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  PROV_DIR  <- file.path(OUT_ROOT, "provenance")
  
  for (d in c(TABLE_DIR, FIG_DIR, PROV_DIR)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  
  STABLE_MIN_DETECTED_CELLS <- 20L
  STABLE_MIN_DONOR_RATE     <- 0.50
  
  COUNT_ASSAY <- "RNA"
  COUNT_LAYER <- "counts"
  COUNT_TOTAL_METADATA_FIELD <- "nCount_RNA"
  
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  DONOR_FIELD      <- ""
  ALLOW_CONV_FALLBACK <- TRUE
  ENFORCE_CONV_PAIR_STRUCTURE <- TRUE
  
  EXPECTED_N_CELLS      <- 272993L
  EXPECTED_SAMPLE_COUNT <- 48L
  EXPECTED_DONOR_COUNT  <- 39L
  EXPECTED_PAIRED_DONOR <- 9L
  EXPECTED_FINE_STATE_COUNT <- 24L
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS + PACKAGES
  # ---------------------------------------------------------------------------
  needed <- c("SeuratObject", "Matrix", "ggplot2", "scales")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  suppressPackageStartupMessages({
    library(SeuratObject)
    library(Matrix)
    library(ggplot2)
    library(scales)
  })
  
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  
  require_cols <- function(x, cols, label) {
    missing <- setdiff(cols, names(x))
    if (length(missing)) {
      stop_msg("Missing required columns in ", label, ": ", paste(missing, collapse = ", "))
    }
    invisible(TRUE)
  }
  
  write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
    write.csv(x, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) {
      unlink(tmp)
      stop_msg("Failed to move temp CSV to final path: ", path)
    }
    invisible(path)
  }
  
  save_pdf_png <- function(p, base, dir = FIG_DIR, width = 7, height = 6) {
    ggsave(file.path(dir, paste0(base, ".pdf")), p, width = width, height = height, device = "pdf")
    ggsave(file.path(dir, paste0(base, ".png")), p, width = width, height = height, dpi = 300)
    invisible(TRUE)
  }
  
  required_files <- c(INPUT_RDS, TABLE_32GENES, TABLE_MODULES, TABLE_E3_OVERLAP, TABLE_E3_RESIDUAL)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) {
    stop_msg("Required input file(s) not found: ", paste(missing_files, collapse = "; "))
  }
  
  # ---------------------------------------------------------------------------
  # 2. FROZEN GENE DEFINITIONS
  # ---------------------------------------------------------------------------
  t32  <- read.csv(TABLE_32GENES, check.names = FALSE, stringsAsFactors = FALSE)
  tmod <- read.csv(TABLE_MODULES, check.names = FALSE, stringsAsFactors = FALSE)
  tE3o <- read.csv(TABLE_E3_OVERLAP, check.names = FALSE, stringsAsFactors = FALSE)
  tE3r <- read.csv(TABLE_E3_RESIDUAL, check.names = FALSE, stringsAsFactors = FALSE)
  
  module_short <- c(
    "Cell cycle / proliferation" = "Cell cycle",
    "Neutrophil degranulation"   = "Neutrophil",
    "Other"                      = "Other",
    "Inflammatory / Down"        = "Inflammatory"
  )
  module_order <- names(module_short)
  
  STATE_ORDER <- c(
    "HSPCs",
    "Cycling_neutrophil_progenitors",
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "IL1R2+_immature_neutrophils",
    "S100A8-9_hi_neutrophils",
    "Mature_neutrophils",
    "Degranulating_neutrophils",
    "Apoptosing_neutrophils",
    "Eosinophils",
    "Mast_cells/eosiniophils",
    "Classical_monocytes",
    "Non-classical_monocytes",
    "cDCs",
    "pDCs",
    "Naive_CD4_T_cells",
    "Memory_CD4_T_cells",
    "Naive_CD8_T_cells",
    "CD8_T_cells",
    "Cycling_TNK",
    "NK cells",
    "B_cells",
    "Plasmablasts",
    "Platelets"
  )
  
  STATE_LINEAGE <- c(
    "HSPCs" = "HSPC",
    "Cycling_neutrophil_progenitors" = "Neutrophil",
    "MPO+_immature_neutrophils_or_progenitors" = "Neutrophil",
    "PADI4+_immature_neutrophils" = "Neutrophil",
    "IL1R2+_immature_neutrophils" = "Neutrophil",
    "S100A8-9_hi_neutrophils" = "Neutrophil",
    "Mature_neutrophils" = "Neutrophil",
    "Degranulating_neutrophils" = "Neutrophil",
    "Apoptosing_neutrophils" = "Neutrophil",
    "Eosinophils" = "Eos/Mast",
    "Mast_cells/eosiniophils" = "Eos/Mast",
    "Classical_monocytes" = "Mono/DC",
    "Non-classical_monocytes" = "Mono/DC",
    "cDCs" = "Mono/DC",
    "pDCs" = "Mono/DC",
    "Naive_CD4_T_cells" = "Lymphoid",
    "Memory_CD4_T_cells" = "Lymphoid",
    "Naive_CD8_T_cells" = "Lymphoid",
    "CD8_T_cells" = "Lymphoid",
    "Cycling_TNK" = "Lymphoid",
    "NK cells" = "Lymphoid",
    "B_cells" = "Lymphoid",
    "Plasmablasts" = "Lymphoid",
    "Platelets" = "Platelets"
  )
  
  LINEAGE_LABEL <- c(
    "HSPC" = "HSPC",
    "Neutrophil" = "Neutrophil",
    "Eos/Mast" = "Eos/Mast",
    "Mono/DC" = "Mono/DC",
    "Lymphoid" = "Lymphoid",
    "Platelets" = "Platelets"
  )
  
  STATE_COLORS <- c(
    "HSPCs" = "#8968D4",
    "Cycling_neutrophil_progenitors" = "#5DBAFC",
    "MPO+_immature_neutrophils_or_progenitors" = "#26A9F1",
    "PADI4+_immature_neutrophils" = "#009AE4",
    "IL1R2+_immature_neutrophils" = "#008BD7",
    "S100A8-9_hi_neutrophils" = "#007DC9",
    "Mature_neutrophils" = "#006EBB",
    "Degranulating_neutrophils" = "#0063B0",
    "Apoptosing_neutrophils" = "#0058A0",
    "Eosinophils" = "#D25D96",
    "Mast_cells/eosiniophils" = "#A53C75",
    "Classical_monocytes" = "#3EB268",
    "Non-classical_monocytes" = "#009342",
    "cDCs" = "#00723B",
    "pDCs" = "#00582F",
    "Naive_CD4_T_cells" = "#E86518",
    "Memory_CD4_T_cells" = "#C33C00",
    "Naive_CD8_T_cells" = "#E9504D",
    "CD8_T_cells" = "#B9003D",
    "Cycling_TNK" = "#FF7D7C",
    "NK cells" = "#982000",
    "B_cells" = "#9B49AE",
    "Plasmablasts" = "#800D77",
    "Platelets" = "#CB9E00"
  )
  
  stopifnot(length(STATE_ORDER) == 24L)
  stopifnot(length(STATE_COLORS) == 24L)
  stopifnot(length(STATE_LINEAGE) == 24L)
  stopifnot(identical(STATE_ORDER, names(STATE_COLORS)))
  stopifnot(identical(STATE_ORDER, names(STATE_LINEAGE)))
  stopifnot(setequal(unique(unname(STATE_LINEAGE)), names(LINEAGE_LABEL)))
  
  require_cols(t32,  c("gene", "direction"), "meta_significant_32_genes.csv")
  require_cols(tmod, c("module", "gene"), "Table_MARS_module_gene_membership.csv")
  require_cols(tE3o, c("signature_gene", "canonical_gene", "module"),
               "Table_Giannini_E3_overlap_direction.csv")
  require_cols(tE3r, c("signature_gene", "canonical_gene", "module"),
               "Table_Giannini_E3_nonoverlapping_signature_genes.csv")
  
  sig <- data.frame(
    gene = t32$gene,
    direction = t32$direction,
    stringsAsFactors = FALSE
  )
  
  if (anyNA(sig$gene) || any(!nzchar(sig$gene))) {
    stop_msg("Empty or NA gene names in 32-gene table.")
  }
  if (anyDuplicated(sig$gene)) {
    stop_msg("Duplicate genes in 32-gene table: ",
             paste(unique(sig$gene[duplicated(sig$gene)]), collapse = ", "))
  }
  
  # ---- Module assignment: E3 tables are source of truth for all 32 genes ----
  e3_module <- rbind(
    tE3o[, c("signature_gene", "module"), drop = FALSE],
    tE3r[, c("signature_gene", "module"), drop = FALSE]
  )
  e3_module <- e3_module[e3_module$signature_gene %in% sig$gene, , drop = FALSE]
  
  e3_module_split <- split(e3_module$module, e3_module$signature_gene)
  e3_module_conflict <- names(e3_module_split)[vapply(
    e3_module_split,
    function(z) length(unique(z[!is.na(z)])) > 1L,
    logical(1)
  )]
  if (length(e3_module_conflict)) {
    stop_msg("Conflicting E3 module assignments: ",
             paste(e3_module_conflict, collapse = ", "))
  }
  
  e3_module_unique <- e3_module[!duplicated(e3_module$signature_gene), , drop = FALSE]
  sig$module <- e3_module_unique$module[match(sig$gene, e3_module_unique$signature_gene)]
  
  gap <- is.na(sig$module)
  if (any(gap)) {
    sig$module[gap] <- tmod$module[match(sig$gene[gap], tmod$gene)]
  }
  
  if (anyNA(sig$module)) {
    stop_msg("Some signature genes lack module assignment: ",
             paste(sig$gene[is.na(sig$module)], collapse = ", "))
  }
  
  # Cross-check where MARS table covers the gene.
  tmod_subset <- tmod[tmod$gene %in% sig$gene, c("gene", "module"), drop = FALSE]
  tmod_split <- split(tmod_subset$module, tmod_subset$gene)
  tmod_conflict <- names(tmod_split)[vapply(
    tmod_split,
    function(z) length(unique(z[!is.na(z)])) > 1L,
    logical(1)
  )]
  if (length(tmod_conflict)) {
    stop_msg("Conflicting MARS module assignments: ",
             paste(tmod_conflict, collapse = ", "))
  }
  
  tmod_unique <- tmod_subset[!duplicated(tmod_subset$gene), , drop = FALSE]
  tmod_cross <- tmod_unique$module[match(sig$gene, tmod_unique$gene)]
  tmod_discord <- !is.na(tmod_cross) & tmod_cross != sig$module
  if (any(tmod_discord)) {
    stop_msg("MARS module table disagrees with E3 module for: ",
             paste(sig$gene[tmod_discord], collapse = ", "))
  }
  
  missing_modules <- setdiff(unique(sig$module), module_order)
  if (length(missing_modules)) {
    stop_msg("Modules missing from module_short: ", paste(missing_modules, collapse = ", "))
  }
  
  e3_overlap <- unique(tE3o$signature_gene)
  e3_residual <- unique(tE3r$signature_gene)
  e3_both <- intersect(e3_overlap, e3_residual)
  if (length(e3_both)) {
    stop_msg("Genes appear in both E3 overlap and residual tables: ",
             paste(e3_both, collapse = ", "))
  }
  
  not_in_e3 <- setdiff(sig$gene, union(e3_overlap, e3_residual))
  if (length(not_in_e3)) {
    stop_msg("Some 32 genes are absent from both E3 tables: ",
             paste(not_in_e3, collapse = ", "))
  }
  
  canon_all <- rbind(
    tE3o[, c("signature_gene", "canonical_gene")],
    tE3r[, c("signature_gene", "canonical_gene")]
  )
  canon_all <- canon_all[canon_all$signature_gene %in% sig$gene, , drop = FALSE]
  
  canon_split <- split(canon_all$canonical_gene, canon_all$signature_gene)
  canon_conflict <- names(canon_split)[vapply(
    canon_split,
    function(z) length(unique(z[!is.na(z)])) > 1L,
    logical(1)
  )]
  if (length(canon_conflict)) {
    stop_msg("Conflicting canonical_gene assignments: ",
             paste(canon_conflict, collapse = ", "))
  }
  
  canon <- canon_all[!duplicated(canon_all$signature_gene), , drop = FALSE]
  sig$canonical_hgnc <- canon$canonical_gene[match(sig$gene, canon$signature_gene)]
  
  if (anyNA(sig$canonical_hgnc) || any(!nzchar(sig$canonical_hgnc))) {
    stop_msg("Missing canonical HGNC for: ",
             paste(sig$gene[is.na(sig$canonical_hgnc) | !nzchar(sig$canonical_hgnc)], collapse = ", "))
  }
  
  sig$E3_membership <- ifelse(
    sig$gene %in% e3_overlap,
    "E3_overlap_15",
    "E3_nonoverlap_residual_17"
  )
  
  stopifnot(nrow(sig) == 32L, length(unique(sig$gene)) == 32L)
  stopifnot(sum(sig$direction == "Up") == 29L, sum(sig$direction == "Down") == 3L)
  stopifnot(length(unique(sig$module)) == 4L)
  stopifnot(sum(sig$E3_membership == "E3_overlap_15") == 15L)
  stopifnot(sum(sig$E3_membership == "E3_nonoverlap_residual_17") == 17L)
  
  histone_map <- setNames(sig$canonical_hgnc, sig$gene)
  stopifnot(identical(unname(histone_map["HIST1H2BM"]), "H2BC14"))
  stopifnot(identical(unname(histone_map["HIST1H3B"]), "H3C2"))
  
  sig <- sig[order(match(sig$module, module_order), sig$gene), ]
  rownames(sig) <- NULL
  write_csv_atomic(sig, file.path(TABLE_DIR, "frozen_signature_input.csv"))
  
  # ---------------------------------------------------------------------------
  # 3. READ OBJECT + STRUCTURE AUDIT
  # ---------------------------------------------------------------------------
  obj <- readRDS(INPUT_RDS)
  meta <- obj@meta.data
  
  if (!is.data.frame(meta)) {
    stop_msg("obj@meta.data is not a data.frame.")
  }
  
  stopifnot(
    ncol(meta) > 0,
    !anyNA(colnames(meta)),
    !anyDuplicated(colnames(meta)),
    all(nzchar(colnames(meta)))
  )
  
  stopifnot(nrow(meta) == ncol(obj))
  stopifnot(identical(rownames(meta), colnames(obj)))
  
  write_csv_atomic(
    data.frame(
      n_cells = ncol(obj),
      n_metadata_fields = ncol(meta),
      assays = paste(Assays(obj), collapse = ";"),
      default_assay = DefaultAssay(obj),
      reductions = paste(Reductions(obj), collapse = ";"),
      stringsAsFactors = FALSE
    ),
    file.path(TABLE_DIR, "object_structure_audit.csv")
  )
  
  if (ncol(obj) != EXPECTED_N_CELLS) {
    stop_msg("Cell count mismatch: expected ", EXPECTED_N_CELLS, ", got ", ncol(obj),
             ". Investigate before proceeding.")
  }
  
  # ---------------------------------------------------------------------------
  # 4. RESOLVE METADATA FIELDS
  # ---------------------------------------------------------------------------
  field_audit <- data.frame(
    field = c(SAMPLE_FIELD, FINE_STATE_FIELD, COUNT_TOTAL_METADATA_FIELD, DONOR_FIELD),
    present = c(SAMPLE_FIELD, FINE_STATE_FIELD, COUNT_TOTAL_METADATA_FIELD, DONOR_FIELD) %in% colnames(meta),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(field_audit, file.path(TABLE_DIR, "configured_metadata_field_audit.csv"))
  
  stopifnot(SAMPLE_FIELD %in% colnames(meta))
  stopifnot(FINE_STATE_FIELD %in% colnames(meta))
  stopifnot(COUNT_TOTAL_METADATA_FIELD %in% colnames(meta))
  
  sample_id  <- as.character(meta[[SAMPLE_FIELD]])
  fine_annot <- as.character(meta[[FINE_STATE_FIELD]])
  ncount_meta <- as.numeric(meta[[COUNT_TOTAL_METADATA_FIELD]])
  
  stopifnot(!anyNA(sample_id), all(nzchar(sample_id)))
  stopifnot(!anyNA(fine_annot), all(nzchar(fine_annot)))
  stopifnot(!anyNA(ncount_meta), all(is.finite(ncount_meta)), all(ncount_meta >= 0))
  
  if (nzchar(DONOR_FIELD)) {
    stopifnot(DONOR_FIELD %in% colnames(meta))
    donor <- as.character(meta[[DONOR_FIELD]])
    donor_source <- paste0("field:", DONOR_FIELD)
  } else if (ALLOW_CONV_FALLBACK) {
    donor <- sub("_CONV$", "", sample_id)
    donor_source <- "fallback: sub('_CONV$','',sample_id)"
  } else {
    stop_msg("No donor source configured.")
  }
  
  stopifnot(!anyNA(donor), all(nzchar(donor)))
  
  sample_map <- unique(data.frame(
    sample_id = sample_id,
    donor = donor,
    stringsAsFactors = FALSE
  ))
  sample_map <- sample_map[order(sample_map$donor, sample_map$sample_id), ]
  rownames(sample_map) <- NULL
  
  if (anyDuplicated(sample_map$sample_id)) {
    dup_samples <- unique(sample_map$sample_id[duplicated(sample_map$sample_id)])
    stop_msg("Same sample_id maps to multiple donors: ", paste(dup_samples, collapse = ", "))
  }
  
  n_samples <- nrow(sample_map)
  n_donors <- length(unique(sample_map$donor))
  donor_sample_tab <- table(sample_map$donor)
  n_paired <- sum(donor_sample_tab == 2L)
  
  if (n_samples != EXPECTED_SAMPLE_COUNT ||
      n_donors != EXPECTED_DONOR_COUNT ||
      n_paired != EXPECTED_PAIRED_DONOR) {
    stop_msg("Donor resolution mismatch (", n_samples, " samples / ", n_donors,
             " donors / ", n_paired, " paired). Expected ",
             EXPECTED_SAMPLE_COUNT, "/", EXPECTED_DONOR_COUNT, "/",
             EXPECTED_PAIRED_DONOR, ".")
  }
  
  if (ENFORCE_CONV_PAIR_STRUCTURE && n_paired > 0) {
    two <- names(donor_sample_tab)[donor_sample_tab == 2L]
    ok2 <- vapply(two, function(d) {
      setequal(
        sample_map$sample_id[sample_map$donor == d],
        c(d, paste0(d, "_CONV"))
      )
    }, logical(1))
    
    if (!all(ok2)) {
      stop_msg("Paired donor sample IDs are not exactly {donor, donor_CONV}: ",
               paste(two[!ok2], collapse = ", "))
    }
  }
  
  donor_audit <- data.frame(
    donor = names(donor_sample_tab),
    n_samples = as.integer(donor_sample_tab),
    sample_ids = vapply(names(donor_sample_tab), function(d) {
      paste(sample_map$sample_id[sample_map$donor == d], collapse = ";")
    }, character(1)),
    source = donor_source,
    stringsAsFactors = FALSE
  )
  write_csv_atomic(donor_audit, file.path(TABLE_DIR, "sample_to_donor_mapping_audit.csv"))
  
  observed_states <- unique(fine_annot)
  missing_states <- setdiff(observed_states, names(STATE_COLORS))
  if (length(missing_states)) {
    stop_msg("Observed fine states not in color map: ", paste(missing_states, collapse = ", "))
  }
  
  extra_states <- setdiff(names(STATE_COLORS), observed_states)
  if (length(extra_states)) {
    stop_msg("Color map states not observed in data: ", paste(extra_states, collapse = ", "))
  }
  
  fine_states <- STATE_ORDER
  
  if (length(fine_states) != EXPECTED_FINE_STATE_COUNT) {
    warning(
      "Fine-state count ", length(fine_states),
      " != expected ", EXPECTED_FINE_STATE_COUNT,
      ". This is a dataset-identity warning."
    )
  }
  
  fine_state_counts <- table(factor(fine_annot, levels = fine_states))
  write_csv_atomic(
    data.frame(
      display_order = seq_along(fine_states),
      fine_state = fine_states,
      n_cells = as.integer(fine_state_counts),
      lineage = unname(STATE_LINEAGE[fine_states]),
      color = unname(STATE_COLORS[fine_states]),
      stringsAsFactors = FALSE
    ),
    file.path(TABLE_DIR, "fine_state_display_order_audit.csv")
  )
  
  # ---------------------------------------------------------------------------
  # 5. SELECT COUNT LAYER + VALIDATE
  # ---------------------------------------------------------------------------
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    lyr <- Layers(obj[[COUNT_ASSAY]])
    stopifnot(COUNT_LAYER %in% lyr)
    cnt <- LayerData(obj, assay = COUNT_ASSAY, layer = COUNT_LAYER, fast = FALSE)
    storage <- paste0("Assay5 layer:", COUNT_LAYER)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, slot = "counts")
    storage <- "legacy Assay counts"
  }
  
  stopifnot(is(cnt, "dgCMatrix"))
  stopifnot(nrow(cnt) > 0, ncol(cnt) > 0)
  stopifnot(!anyNA(rownames(cnt)), !anyDuplicated(rownames(cnt)), all(nzchar(rownames(cnt))))
  stopifnot(!anyNA(colnames(cnt)), !anyDuplicated(colnames(cnt)), all(nzchar(colnames(cnt))))
  stopifnot(identical(colnames(cnt), colnames(obj)))
  stopifnot(all(is.finite(cnt@x)), all(cnt@x >= 0), all(cnt@x == floor(cnt@x)))
  
  count_sums <- as.numeric(Matrix::colSums(cnt))
  max_diff <- max(abs(count_sums - ncount_meta))
  if (!is.finite(max_diff) || max_diff > 1e-8) {
    stop_msg("Count layer column sums deviate from ", COUNT_TOTAL_METADATA_FIELD,
             " (max abs diff ", format(max_diff, digits = 3), "). Wrong assay/layer?")
  }
  
  write_csv_atomic(
    data.frame(
      assay = COUNT_ASSAY,
      layer = COUNT_LAYER,
      storage = storage,
      n_features = nrow(cnt),
      n_cells = ncol(cnt),
      max_abs_diff_ncount = max_diff,
      stringsAsFactors = FALSE
    ),
    file.path(TABLE_DIR, "selected_count_layer_validation.csv")
  )
  
  # ---------------------------------------------------------------------------
  # 6. GENE MAPPING
  # ---------------------------------------------------------------------------
  feat <- rownames(cnt)
  
  mapping <- do.call(rbind, lapply(seq_len(nrow(sig)), function(i) {
    g <- sig$gene[i]
    cg <- sig$canonical_hgnc[i]
    query <- unique(c(g, cg))
    
    hits <- feat[feat %in% query]
    
    if (length(hits) == 0L) {
      ci <- feat[tolower(feat) %in% tolower(query)]
      mt <- if (length(ci) == 0L) "absent" else "case_insensitive_only"
      
      return(data.frame(
        signature_gene = g,
        module = sig$module[i],
        canonical_hgnc = cg,
        dataset_feature = if (mt == "absent") NA_character_ else paste(ci, collapse = ";"),
        mapping_type = mt,
        stringsAsFactors = FALSE
      ))
    }
    
    mt <- if (length(hits) == 1L) {
      if (identical(hits, g)) "exact" else "alias_canonical"
    } else {
      "ambiguous_multiple"
    }
    
    data.frame(
      signature_gene = g,
      module = sig$module[i],
      canonical_hgnc = cg,
      dataset_feature = paste(hits, collapse = ";"),
      mapping_type = mt,
      stringsAsFactors = FALSE
    )
  }))
  
  write_csv_atomic(mapping, file.path(TABLE_DIR, "gene_mapping_audit.csv"))
  write_csv_atomic(
    mapping[mapping$signature_gene %in% c("HIST1H2BM", "HIST1H3B"), , drop = FALSE],
    file.path(TABLE_DIR, "histone_alias_verification.csv")
  )
  
  mapped_rows_precheck <- mapping$mapping_type %in% c("exact", "alias_canonical")
  dup_feature <- names(which(table(mapping$dataset_feature[mapped_rows_precheck]) > 1))
  if (length(dup_feature)) {
    stop_msg("Same count feature serves multiple signature genes: ",
             paste(dup_feature, collapse = ", "))
  }
  
  if (any(mapping$mapping_type == "ambiguous_multiple")) {
    stop_msg("Ambiguous mapping, legacy and current names both present: ",
             paste(mapping$signature_gene[mapping$mapping_type == "ambiguous_multiple"], collapse = ", "))
  }
  
  if (any(mapping$mapping_type == "case_insensitive_only")) {
    stop_msg("Case-insensitive-only candidates unresolved; inspect gene_mapping_audit.csv: ",
             paste(mapping$signature_gene[mapping$mapping_type == "case_insensitive_only"], collapse = ", "))
  }
  
  mapped_rows <- mapping$mapping_type %in% c("exact", "alias_canonical")
  mapped_genes <- mapping$signature_gene[mapped_rows]
  mapped_feat  <- mapping$dataset_feature[mapped_rows]
  
  stopifnot(length(mapped_genes) == length(mapped_feat))
  stopifnot(!anyNA(mapped_feat))
  stopifnot(!any(grepl(";", mapped_feat, fixed = TRUE)))
  stopifnot(!anyDuplicated(mapped_feat))
  
  idx <- match(mapped_feat, feat)
  stopifnot(!anyNA(idx))
  
  X <- cnt[idx, , drop = FALSE]
  rownames(X) <- mapped_genes
  
  # ---------------------------------------------------------------------------
  # 7. DETECTION METRICS
  # ---------------------------------------------------------------------------
  cell_donor <- factor(donor, levels = unique(donor))
  cell_state <- factor(fine_annot, levels = fine_states)
  stopifnot(!anyNA(cell_donor), !anyNA(cell_state))
  
  D <- Matrix::sparseMatrix(
    i = as.integer(cell_donor),
    j = seq_along(cell_donor),
    x = 1,
    dims = c(nlevels(cell_donor), length(cell_donor)),
    dimnames = list(levels(cell_donor), colnames(obj))
  )
  
  S <- Matrix::sparseMatrix(
    i = as.integer(cell_state),
    j = seq_along(cell_state),
    x = 1,
    dims = c(nlevels(cell_state), length(cell_state)),
    dimnames = list(levels(cell_state), colnames(obj))
  )
  
  global_umi  <- Matrix::rowSums(X)
  n_cells_det <- Matrix::rowSums(X > 0)
  
  # Robust dense aggregation (k <= 32 -> cheap; dense %*% never drops dims)
  X_dense  <- as.matrix(X)                # k x n_cells
  D_dense  <- as.matrix(Matrix::t(D))     # n_cells x n_donors
  S_dense  <- as.matrix(Matrix::t(S))     # n_cells x n_states
  
  donor_umi <- X_dense %*% D_dense        # k x n_donors
  colnames(donor_umi) <- levels(cell_donor)
  rownames(donor_umi) <- mapped_genes
  stopifnot(nrow(donor_umi) == length(mapped_genes),
            ncol(donor_umi) == length(levels(cell_donor)))
  
  donor_det <- donor_umi > 0
  
  state_n_det <- (X_dense > 0) %*% S_dense   # k x n_states
  colnames(state_n_det) <- levels(cell_state)
  rownames(state_n_det) <- mapped_genes
  stopifnot(nrow(state_n_det) == length(mapped_genes),
            ncol(state_n_det) == length(levels(cell_state)))
  
  state_sizes <- as.numeric(table(cell_state)[fine_states])
  stopifnot(length(state_sizes) == length(fine_states), all(state_sizes > 0))
  
  pct_state <- sweep(state_n_det, 2, state_sizes, "/") * 100
  colnames(pct_state) <- fine_states
  rownames(pct_state) <- mapped_genes
  
  if (any(pct_state < -1e-8 | pct_state > 100 + 1e-8, na.rm = TRUE)) {
    stop_msg("Internal error: pct_state outside [0, 100].")
  }
  pct_state <- pmin(100, pmax(0, pct_state))
  pct_state <- matrix(as.numeric(pct_state),
                      nrow = length(mapped_genes),
                      dimnames = list(mapped_genes, fine_states))
  
  det <- data.frame(
    signature_gene = mapped_genes,
    global_total_umi = unname(global_umi),
    n_cells_detected = unname(n_cells_det),
    pct_cells_detected = 100 * unname(n_cells_det) / ncol(obj),
    n_donors_detected = rowSums(donor_det),
    donor_detection_rate = rowSums(donor_det) / nrow(D),
    stringsAsFactors = FALSE
  )
  
  absent_genes <- sig$gene[!sig$gene %in% mapped_genes]
  if (length(absent_genes)) {
    det <- rbind(
      det,
      data.frame(
        signature_gene = absent_genes,
        global_total_umi = NA_real_,
        n_cells_detected = NA_real_,
        pct_cells_detected = NA_real_,
        n_donors_detected = NA_real_,
        donor_detection_rate = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  det$module <- sig$module[match(det$signature_gene, sig$gene)]
  det$direction <- sig$direction[match(det$signature_gene, sig$gene)]
  det$E3 <- sig$E3_membership[match(det$signature_gene, sig$gene)]
  
  det <- det[order(match(det$module, module_order), det$signature_gene), ]
  rownames(det) <- NULL
  
  det$class <- ifelse(
    is.na(det$global_total_umi), "feature_missing",
    ifelse(
      det$global_total_umi == 0, "present_not_detected",
      ifelse(
        det$n_cells_detected >= STABLE_MIN_DETECTED_CELLS &
          det$donor_detection_rate >= STABLE_MIN_DONOR_RATE,
        "stable_detected", "low_detected"
      )
    )
  )
  
  det$class <- factor(
    det$class,
    levels = c("feature_missing", "present_not_detected", "low_detected", "stable_detected")
  )
  det$evaluable <- det$class == "stable_detected"
  
  write_csv_atomic(det, file.path(TABLE_DIR, "final_gene_evaluability_32.csv"))
  write_csv_atomic(
    det[det$class == "stable_detected", , drop = FALSE],
    file.path(TABLE_DIR, "primary_evaluable_genes_stable.csv")
  )
  
  k_summary <- data.frame(
    metric = c("mapped", "detected", "stable"),
    k = c(
      sum(!is.na(det$global_total_umi)),
      sum(!is.na(det$global_total_umi) & det$global_total_umi > 0),
      sum(det$class == "stable_detected")
    ),
    n_32 = rep(nrow(sig), 3L),
    stringsAsFactors = FALSE
  )
  k_summary$k_of_32 <- paste0(k_summary$k, " / ", k_summary$n_32)
  write_csv_atomic(k_summary, file.path(TABLE_DIR, "coverage_k_of_32_summary.csv"))
  
  donor_cols <- colnames(donor_umi)
  donor_df_wide <- as.data.frame(donor_umi, check.names = FALSE)
  donor_df_wide$signature_gene <- rownames(donor_umi)
  gd <- merge(
    det[, c("signature_gene","module","direction","E3","class",
            "global_total_umi","donor_detection_rate"), drop = FALSE],
    donor_df_wide, by = "signature_gene", all.x = TRUE, sort = FALSE)
  write_csv_atomic(gd, file.path(TABLE_DIR, "gene_by_donor_detectability.csv"))
  
  pct_df <- as.data.frame(pct_state, check.names = FALSE)
  pct_df$signature_gene <- rownames(pct_state)
  gs <- merge(
    det[, c("signature_gene","module","direction","E3","class"), drop = FALSE],
    pct_df, by = "signature_gene", all.x = TRUE, sort = FALSE)
  write_csv_atomic(gs, file.path(TABLE_DIR, "gene_by_fine_state_detectability.csv"))
  
  # ---------------------------------------------------------------------------
  # 8. VISUALIZATION
  # ---------------------------------------------------------------------------
  cell_colors <- STATE_COLORS
  lineage_of <- STATE_LINEAGE
  lineage_short <- LINEAGE_LABEL
  
  plot_fine_states <- STATE_ORDER
  
  det_colors <- c(
    feature_missing = "#9E9E9E",
    present_not_detected = "#CFCFCF",
    low_detected = "#2E6FD8",
    stable_detected = "#E03131"
  )
  
  class_labels <- c(
    feature_missing = "feature missing (NA)",
    present_not_detected = "present, 0 UMI",
    low_detected = "low detection",
    stable_detected = "stable detection"
  )
  
  raster_colors <- c(
    "detected" = "#E62E32",
    "0 UMI" = "#E0E0E0",
    "NA" = "#858585"
  )
  
  k_colors <- c(
    mapped = "#9E9E9E",
    detected = "#2E6FD8",
    stable = "#E03131"
  )
  
  class_shapes <- c(
    feature_missing = 4,
    present_not_detected = 17,
    low_detected = 15,
    stable_detected = 19
  )
  
  style_theme <- theme_classic(base_size = 11) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks = element_line(color = "grey30", linewidth = 0.3),
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, hjust = 0, color = "black", margin = margin(b = 6)),
      plot.margin = margin(8, 8, 8, 8),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key = element_blank(),
      legend.key.size = unit(3.5, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 9, color = "black", face = "bold")
    )
  
  # Fig 1
  d1 <- det[!is.na(det$donor_detection_rate), , drop = FALSE]
  d1$module_f <- factor(unname(module_short[as.character(d1$module)]),
                        levels = unname(module_short))
  d1$gene_f <- factor(d1$signature_gene, levels = rev(d1$signature_gene))
  
  p1 <- ggplot(d1, aes(x = gene_f, y = 100 * donor_detection_rate, fill = class)) +
    geom_col(width = 0.72) +
    facet_grid(module_f ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = det_colors, labels = class_labels, name = "Detection class") +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.02)),
      labels = function(x) paste0(x, "%")
    ) +
    coord_flip() +
    labs(x = NULL, y = "Donor detection rate", title = "Gene detection audit") +
    style_theme +
    theme(
      panel.spacing = unit(3, "mm"),
      strip.text.y = element_text(angle = 0, size = 7.5),
      plot.margin = margin(8, 80, 8, 8)
    )
  save_pdf_png(p1, "fig01_gene_detection_rates", width = 7, height = 8.5)
  
  # Fig 2
  det_mod <- factor(det$module, levels = module_order)
  mod_list <- split(det$signature_gene, det_mod, drop = TRUE)
  mods_rev <- rev(names(mod_list))
  mods_rev_short <- unname(module_short[mods_rev])
  
  gene_levels <- character(0)
  for (i in seq_along(mods_rev)) {
    gene_levels <- c(gene_levels, mod_list[[mods_rev[i]]])
    gene_levels <- c(gene_levels, paste0("spacer|", mods_rev_short[i]))
  }
  
  heat_df <- do.call(rbind, lapply(gene_levels, function(g) {
    data.frame(gene = g, fine_state = plot_fine_states, stringsAsFactors = FALSE)
  }))
  
  heat_df$pct <- vapply(seq_len(nrow(heat_df)), function(i) {
    g <- heat_df$gene[i]
    s <- heat_df$fine_state[i]
    if (grepl("^spacer", g)) {
      NA_real_
    } else if (g %in% rownames(pct_state)) {
      pct_state[g, s]
    } else {
      NA_real_
    }
  }, numeric(1))
  
  heat_df$fill_level <- as.character(cut(
    heat_df$pct,
    breaks = c(0, 20, 40, 60, 80, 100),
    include.lowest = TRUE
  ))
  heat_df$fill_level[grepl("^spacer", heat_df$gene)] <- "spacer"
  heat_df$fill_level[is.na(heat_df$pct) & !grepl("^spacer", heat_df$gene)] <- "feature missing"
  heat_df$rowtype <- "heat"
  
  strip_df <- data.frame(
    gene = "lineage",
    fine_state = plot_fine_states,
    fill_level = plot_fine_states,
    rowtype = "strip",
    stringsAsFactors = FALSE
  )
  
  state_lineage <- unname(lineage_of[plot_fine_states])
  stopifnot(!anyNA(state_lineage))
  
  lineage_order <- unique(unname(lineage_of))
  boundary_states <- character(0)
  if (length(lineage_order) > 1) {
    boundary_states <- vapply(lineage_order[-length(lineage_order)], function(lg) {
      st <- plot_fine_states[state_lineage == lg]
      if (length(st)) st[length(st)] else NA_character_
    }, character(1))
    boundary_states <- boundary_states[!is.na(boundary_states)]
  }
  
  comb <- rbind(
    heat_df[, c("gene", "fine_state", "fill_level", "rowtype")],
    strip_df[, c("gene", "fine_state", "fill_level", "rowtype")]
  )
  
  heat_fill_breaks <- c("[0,20]", "(20,40]", "(40,60]", "(60,80]", "(80,100]", "feature missing")
  heat_fill_values <- c(
    "[0,20]" = "#2E6FD8",
    "(20,40]" = "#8AB2DB",
    "(40,60]" = "#D5E3F6",
    "(60,80]" = "#F2A79B",
    "(80,100]" = "#E03131",
    "feature missing" = "#858585",
    "spacer" = "#FFFFFF"
  )
  fill_values <- c(heat_fill_values, cell_colors)
  
  comb$fill_level <- factor(
    comb$fill_level,
    levels = c(heat_fill_breaks, "spacer", plot_fine_states)
  )
  comb$gene <- factor(comb$gene, levels = c(gene_levels, "lineage"))
  comb$fine_state <- factor(comb$fine_state, levels = plot_fine_states)
  comb$rowtype <- factor(comb$rowtype, levels = c("heat", "strip"))
  comb$is_boundary <- comb$rowtype == "strip" &
    as.character(comb$fine_state) %in% boundary_states
  
  p2 <- ggplot(comb, aes(x = fine_state, y = gene, fill = fill_level)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_tile(
      data = comb[comb$is_boundary, , drop = FALSE],
      aes(fill = fill_level),
      color = "white",
      linewidth = 1.1
    ) +
    facet_grid(rowtype ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = fill_values,
      breaks = heat_fill_breaks,
      labels = heat_fill_breaks,
      name = "% cells\ndetected",
      drop = FALSE
    ) +
    scale_y_discrete(labels = function(l) {
      ifelse(
        grepl("^spacer\\|", l),
        sub("^spacer\\|", "", l),
        ifelse(l == "lineage", "", l)
      )
    }) +
    labs(x = "Fine cell state", y = NULL, title = "Gene x state detection") +
    style_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, color = "black"),
      axis.text.y = element_text(size = 7.5),
      panel.spacing = unit(4, "mm"),
      strip.text.y = element_blank(),
      strip.background.y = element_blank(),
      legend.key.height = unit(5, "mm"),
      plot.margin = margin(8, 8, 20, 8)
    )
  save_pdf_png(p2, "fig02_gene_x_state_heatmap", width = 12, height = 9)
  
  # Fig 3
  donor_long <- do.call(rbind, lapply(donor_cols, function(dn) {
    data.frame(
      signature_gene = gd$signature_gene,
      donor = dn,
      umi = gd[[dn]],
      stringsAsFactors = FALSE
    )
  }))
  
  donor_long$module <- det$module[match(donor_long$signature_gene, det$signature_gene)]
  donor_long$class <- det$class[match(donor_long$signature_gene, det$signature_gene)]
  donor_long$det_f <- factor(
    ifelse(is.na(donor_long$umi), "NA",
           ifelse(donor_long$umi > 0, "detected", "0 UMI")),
    levels = c("detected", "0 UMI", "NA")
  )
  donor_long$module_f <- factor(
    unname(module_short[as.character(donor_long$module)]),
    levels = unname(module_short)
  )
  donor_long$gene_f <- factor(donor_long$signature_gene, levels = rev(det$signature_gene))
  donor_long$donor_f <- factor(donor_long$donor, levels = donor_cols)
  
  p3 <- ggplot(donor_long, aes(x = donor_f, y = gene_f, fill = det_f)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_grid(module_f ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = raster_colors, name = NULL) +
    labs(x = "Donor", y = NULL, title = "Donor-level detection") +
    style_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.5),
      axis.text.y = element_text(size = 7.5),
      panel.spacing = unit(3, "mm"),
      strip.text.y = element_text(angle = 0, size = 8)
    )
  save_pdf_png(p3, "fig03_donor_detection_raster", width = 12, height = 8.5)
  
  # Fig 4
  d4 <- det[!is.na(det$global_total_umi), , drop = FALSE]
  d4$umi_plot <- d4$global_total_umi + 1
  
  p4 <- ggplot(d4, aes(x = umi_plot, y = n_cells_detected, color = class, shape = class)) +
    geom_point(size = 2.4, alpha = 0.9) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    geom_hline(
      yintercept = STABLE_MIN_DETECTED_CELLS,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = 3e3,
      y = STABLE_MIN_DETECTED_CELLS,
      label = paste0(">= ", STABLE_MIN_DETECTED_CELLS, " cells"),
      hjust = 0,
      vjust = 1.4,
      size = 3,
      color = "grey25"
    ) +
    scale_color_manual(values = det_colors, labels = class_labels, name = "Detection class") +
    scale_shape_manual(values = class_shapes, labels = class_labels, name = "Detection class") +
    labs(
      x = "Total UMI + 1 (log10)",
      y = "Cells detected",
      title = "UMI vs detected cells, thresholds"
    ) +
    style_theme
  save_pdf_png(p4, "fig04_umi_vs_cells_scatter", width = 7, height = 5.5)
  
  # Fig 5
  kdf <- data.frame(
    metric = factor(k_summary$metric, levels = k_summary$metric),
    k = k_summary$k,
    label = k_summary$k_of_32,
    stringsAsFactors = FALSE
  )
  
  p5 <- ggplot(kdf, aes(x = metric, y = k, fill = metric)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = k_colors, guide = "none") +
    scale_y_continuous(limits = c(0, 35), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Genes", title = "k / 32 summary") +
    style_theme +
    theme(axis.text.x = element_text(size = 10))
  save_pdf_png(p5, "fig05_k_of_32_summary", width = 5.5, height = 4.5)
  
  # Fig 6
  pal_df <- data.frame(
    cell = names(cell_colors),
    lineage = unname(lineage_of),
    stringsAsFactors = FALSE
  )
  pal_df$cell <- factor(pal_df$cell, levels = rev(names(cell_colors)))
  pal_df$lineage <- factor(
    pal_df$lineage,
    levels = lineage_order,
    labels = unname(lineage_short[lineage_order])
  )
  
  p6 <- ggplot(pal_df, aes(x = cell, y = 1, fill = cell)) +
    geom_tile(height = 0.95, color = "white", linewidth = 0.4) +
    facet_grid(. ~ lineage, scales = "free_x", space = "free") +
    scale_fill_manual(values = cell_colors, guide = "none") +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = NULL, title = "Cell-type palette (lineage-grouped)") +
    style_theme +
    theme(
      axis.text.y = element_text(size = 8, hjust = 1, margin = margin(r = 5)),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.spacing = unit(5, "mm"),
      strip.text.x = element_text(size = 7, face = "bold")
    )
  save_pdf_png(p6, "fig06_cell_palette", width = 10.5, height = 7)
  
  # ---------------------------------------------------------------------------
  # 9. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  write_csv_atomic(
    data.frame(
      key = c(
        "script",
        "script_version",
        "input_rds",
        "results_meta",
        "n_cells",
        "n_samples",
        "n_donors",
        "n_paired",
        "donor_source",
        "n_fine_states",
        "expected_fine_state_count",
        "module_source",
        "count_assay",
        "count_layer",
        "max_abs_diff_ncount",
        "stable_min_detected_cells",
        "stable_min_donor_rate",
        "mapped_k",
        "detected_k",
        "stable_k",
        "R_version"
      ),
      value = c(
        "01_gene_mapping_detectability_audit.R",
        "v4",
        INPUT_RDS,
        RESULTS_META,
        as.character(ncol(obj)),
        as.character(n_samples),
        as.character(n_donors),
        as.character(n_paired),
        donor_source,
        as.character(length(fine_states)),
        as.character(EXPECTED_FINE_STATE_COUNT),
        "Giannini E3 tables primary; MARS tmod fallback/cross-check",
        COUNT_ASSAY,
        COUNT_LAYER,
        as.character(max_diff),
        as.character(STABLE_MIN_DETECTED_CELLS),
        as.character(STABLE_MIN_DONOR_RATE),
        as.character(k_summary$k[k_summary$metric == "mapped"]),
        as.character(k_summary$k[k_summary$metric == "detected"]),
        as.character(k_summary$k[k_summary$metric == "stable"]),
        as.character(getRversion())
      ),
      stringsAsFactors = FALSE
    ),
    file.path(PROV_DIR, "run_provenance.csv")
  )
  
  writeLines(capture.output(sessionInfo()), file.path(PROV_DIR, "sessionInfo.txt"))
  
  cat("Script 01 v4 complete.\n")
  cat("Run dir:", OUT_ROOT, "\n")
  cat(sprintf(
    "k/32: mapped %d, detected %d, stable %d (main k/32 = %d)\n",
    k_summary$k[k_summary$metric == "mapped"],
    k_summary$k[k_summary$metric == "detected"],
    k_summary$k[k_summary$metric == "stable"],
    k_summary$k[k_summary$metric == "stable"]
  ))
  cat("Stable genes:",
      paste(det$signature_gene[det$class == "stable_detected"], collapse = ", "),
      "\n")
})