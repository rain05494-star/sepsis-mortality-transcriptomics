# ============================================================================
# SCRIPT 05a v1.5 -- Signature ablation: source-retention robustness
#   Purpose: leave-one-gene-out ablation on the full / overlap / residual gene
#   sets (cell-cycle and neutrophil modules) to test whether bulk-like source
#   localization is stable when individual genes are removed.
#   ROBUSTNESS DIAGNOSTIC, not causal gene-importance evidence.
#   rho = Spearman between the leave-one-out profile and the ORIGINAL same-set
#   profile. Guards: >=3 genes after removal, non-constant profile; <5 genes
#   flagged descriptive.
#   Outputs:
#     ablation_source_retention_summary.csv
#     ablation_set_summary.csv
#     figures/fig_ablation_residual_cellcycle.(pdf/png)
#     figures/fig_ablation_residual_neutrophil.(pdf/png)
# RUN: Rscript 05_signature_ablation_source_retention.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  PROJECT_ROOT <- "/home/sunshine/predicate/singlecell"
  OUTPUTS_ROOT <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  SCRIPT02_ROOT <- file.path(OUTPUTS_ROOT, "02_gene_cell_distribution")
  
  STAGE_NAME  <- "05_signature_ablation_source_retention"
  RUN_PURPOSE <- "signature_ablation"
  
  MODULE_FULL  <- c("Cell cycle", "Neutrophil")
  OVERLAP_TAG  <- "E3_overlap_15"
  RESIDUAL_TAG <- "E3_nonoverlap_residual_17"
  
  MIN_SET_FORMAL <- 5L
  MIN_RHO_GENES  <- 3L
  RHO_EVALUABLE_SD <- 1e-6
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  require_cols <- function(x, cols, label) {
    miss <- setdiff(cols, names(x))
    if (length(miss)) stop_msg(label, " missing required columns: ", paste(miss, collapse = ", "))
    invisible(TRUE)
  }
  safe_median <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else median(x)}
  safe_mean   <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else mean(x)}
  safe_min    <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else min(x)}
  safe_sd     <- function(x){x<-x[is.finite(x)]; if(length(x)<2) NA_real_ else sd(x)}
  fmt_num     <- function(x, digits=2) ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "NA")
  write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
    write.csv(x, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) { unlink(tmp); stop_msg("Failed to move temp CSV to final path: ", path) }
    invisible(path)
  }
  find_latest_file <- function(root, pattern) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    fs <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (!length(fs)) stop_msg("No file found under ", root, " matching: ", pattern)
    fs[which.max(file.info(fs)$mtime)]
  }
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
  module_short <- function(x) {
    y <- as.character(x)
    y[y %in% c("Cell cycle / proliferation","Cell cycle")] <- "Cell cycle"
    y[y %in% c("Neutrophil degranulation","Neutrophil")] <- "Neutrophil"
    y[y %in% c("Inflammatory / Down","Inflammatory")] <- "Inflammatory"
    y[y %in% c("Other")] <- "Other"
    y
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. INPUTS (single locked table dir per stage)
  # ---------------------------------------------------------------------------
  FROZEN_FILE <- find_latest_file(SCRIPT01_ROOT, "^frozen_signature_input\\.csv$")
  STATE_FILE  <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  LONG_FILE   <- find_latest_file(SCRIPT02_ROOT, "^donor_state_gene_counts_long\\.csv$")
  
  frozen <- read.csv(FROZEN_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  long <- read.csv(LONG_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  
  require_cols(frozen, c("gene","module","direction","E3_membership"), "frozen_signature_input.csv")
  require_cols(state01, c("display_order","fine_state"), "fine_state_display_order_audit.csv")
  require_cols(long, c("gene","donor","fine_state","umi","donor_gene_total_umi"), "donor_state_gene_counts_long.csv")
  
  long$umi <- suppressWarnings(as.numeric(long$umi))
  long$donor_gene_total_umi <- suppressWarnings(as.numeric(long$donor_gene_total_umi))
  if (any(!is.finite(long$umi))) stop_msg("Non-finite umi in donor_state_gene_counts_long.csv")
  if (any(!is.finite(long$donor_gene_total_umi))) {
    stop_msg("Non-finite donor_gene_total_umi in donor_state_gene_counts_long.csv")
  }
  
  if (!all(frozen$E3_membership %in% c(OVERLAP_TAG, RESIDUAL_TAG))) stop_msg("Unexpected E3 labels")
  stopifnot(nrow(frozen) == 32L, nrow(state01) == 24L)
  stopifnot(sum(frozen$E3_membership == OVERLAP_TAG) == 15L)
  stopifnot(sum(frozen$E3_membership == RESIDUAL_TAG) == 17L)
  stopifnot(all(frozen$gene %in% unique(long$gene)))
  
  fine_states <- state01$fine_state[order(state01$display_order)]
  frozen$module_s <- module_short(frozen$module)
  donors <- sort(unique(long$donor))
  stopifnot(length(donors) == 26L)
  
  sets <- list()
  for (mod in MODULE_FULL) {
    in_mod <- frozen$module_s == mod
    in_overlap <- frozen$E3_membership == OVERLAP_TAG
    in_residual <- frozen$E3_membership == RESIDUAL_TAG
    
    sets[[paste0("full_", mod)]] <- sort(frozen$gene[in_mod])
    sets[[paste0("overlap_", mod)]] <- sort(frozen$gene[in_mod & in_overlap])
    sets[[paste0("residual_", mod)]] <- sort(frozen$gene[in_mod & in_residual])
  }
  
  set_sizes <- vapply(sets, length, integer(1))
  if (any(set_sizes == 0L)) {
    stop_msg(
      "Empty ablation gene set(s): ",
      paste(names(set_sizes)[set_sizes == 0L], collapse = ", ")
    )
  }
  
  # ---------------------------------------------------------------------------
  # 3. SET-LEVEL PROFILES (donor-aware)
  # ---------------------------------------------------------------------------
  build_set_cf <- function(genes) {
    sub <- long[long$gene %in% genes,
                c("donor","gene","fine_state","umi","donor_gene_total_umi"), drop=FALSE]
    mat <- matrix(NA_real_, nrow=length(donors), ncol=length(fine_states),
                  dimnames=list(donors, fine_states))
    if (!nrow(sub)) return(mat)
    gd <- unique(sub[, c("donor","gene","donor_gene_total_umi"), drop=FALSE])
    tot_d <- tapply(gd$donor_gene_total_umi, gd$donor, sum, na.rm=TRUE)
    valid_donors <- intersect(names(tot_d)[is.finite(tot_d) & tot_d > 0], donors)
    mat[valid_donors, ] <- 0
    sub <- sub[sub$donor %in% valid_donors & sub$fine_state %in% fine_states, , drop=FALSE]
    umi_ds <- tapply(sub$umi, paste(sub$donor, sub$fine_state, sep="||"), sum, na.rm=TRUE)
    for (i in seq_along(umi_ds)) {
      p <- strsplit(names(umi_ds)[i], "\\|\\|")[[1]]
      d <- p[1]; s <- p[2]
      mat[d, s] <- umi_ds[i] / tot_d[[d]]
    }
    mat
  }
  set_profile <- function(cf) {
    vapply(fine_states, function(s) safe_median(cf[, s]), numeric(1))
  }
  top_state <- function(prof) {
    ok <- which(is.finite(prof))
    if (!length(ok)) return(NA_character_)
    fine_states[ok[which.max(prof[ok])]]
  }
  delta12 <- function(prof) {
    ok <- which(is.finite(prof))
    if (length(ok) < 2L) return(NA_real_)
    v <- sort(prof[ok], decreasing=TRUE)
    v[1] - v[2]
  }
  spearman_pair <- function(p1, p2) {
    ok <- is.finite(p1) & is.finite(p2)
    if (sum(ok) < 3L) return(NA_real_)
    if (safe_sd(p1[ok]) < RHO_EVALUABLE_SD || safe_sd(p2[ok]) < RHO_EVALUABLE_SD) return(NA_real_)
    cor(p1[ok], p2[ok], method="spearman")
  }
  
  set_profiles <- lapply(sets, function(gs) set_profile(build_set_cf(gs)))
  
  # ---------------------------------------------------------------------------
  # 4. LEAVE-ONE-GENE-OUT ABLATION
  # ---------------------------------------------------------------------------
  rows <- list()
  for (mod in MODULE_FULL) for (tag in c("full","overlap","residual")) {
    k <- paste0(tag, "_", mod)
    genes <- sets[[k]]
    if (length(genes) < 2L) next
    prof_orig <- set_profiles[[k]]
    top1_orig <- top_state(prof_orig)
    del_orig  <- delta12(prof_orig)
    for (g in genes) {
      genes_loo <- setdiff(genes, g)
      n_after <- length(genes_loo)
      if (n_after >= MIN_RHO_GENES) {
        prof_loo <- set_profile(build_set_cf(genes_loo))
        rho <- spearman_pair(prof_orig, prof_loo)
      } else {
        prof_loo <- NULL
        rho <- NA_real_
      }
      top1_loo <- if (!is.null(prof_loo)) top_state(prof_loo) else NA_character_
      flipped <- !is.na(top1_orig) && !is.na(top1_loo) && !identical(top1_loo, top1_orig)
      rows[[length(rows)+1L]] <- data.frame(
        module = mod,
        set = tag,
        removed_gene = g,
        set_n_genes = length(genes),
        n_genes_after = n_after,
        set_scope = if (n_after >= MIN_SET_FORMAL) "formal_set_level" else "descriptive_small_gene_set",
        descriptive_small_set = n_after < MIN_SET_FORMAL,
        top1_original = top1_orig,
        top1_leave_one_out = top1_loo,
        top1_flipped = flipped,
        rho_vs_original_set = rho,
        rho_evaluable = is.finite(rho),
        original_top1_minus_top2_delta = del_orig,
        stringsAsFactors=FALSE)
    }
  }
  abl <- do.call(rbind, rows)
  write_csv_atomic(abl, file.path(TABLE_DIR, "ablation_source_retention_summary.csv"))
  
  set_sum <- do.call(rbind, lapply(split(abl, paste(abl$module, abl$set, sep="|")), function(z) {
    n_eval <- sum(z$rho_evaluable, na.rm = TRUE)
    n_flip <- sum(z$top1_flipped, na.rm = TRUE)
    
    data.frame(
      module = z$module[1],
      set = z$set[1],
      set_n_genes = z$set_n_genes[1],
      n_ablations = nrow(z),
      n_top1_flips = n_flip,
      top1_retention_rate = 1 - n_flip / nrow(z),
      n_rho_evaluable = n_eval,
      rho_mean = safe_mean(z$rho_vs_original_set),
      rho_min = safe_min(z$rho_vs_original_set),
      stringsAsFactors = FALSE
    )
  }))
  write_csv_atomic(set_sum, file.path(TABLE_DIR, "ablation_set_summary.csv"))
  
  # ---------------------------------------------------------------------------
  # 5. FIGURES (residual sets)
  # ---------------------------------------------------------------------------
  if (requireNamespace("ggplot2", quietly=TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    for (mod in c("Cell cycle", "Neutrophil")) {
      z <- abl[abl$set == "residual" & abl$module == mod, , drop=FALSE]
      if (!nrow(z)) next
      z$gene_f <- factor(z$removed_gene, levels=rev(z$removed_gene))
      
      full_rho_ref <- spearman_pair(
        set_profiles[[paste0("full_", mod)]],
        set_profiles[[paste0("residual_", mod)]]
      )
      
      title_txt <- paste0(
        "Residual ", mod, ": leave-one-out source retention",
        if (is.finite(full_rho_ref)) paste0(" | full-vs-residual rho=", fmt_num(full_rho_ref, 2)) else ""
      )
      
      p <- ggplot(z, aes(x = gene_f, y = rho_vs_original_set, fill = top1_flipped)) +
        geom_col(width = 0.7, colour = "grey25", size = 0.25) +
        scale_fill_manual(
          values = c("TRUE" = "#E07A5F", "FALSE" = "#3D8DAE"),
          name = "Top-1 flipped"
        ) +
        coord_flip() +
        theme_classic(base_size = 8) +
        theme(
          axis.text.y = element_text(size = 7),
          legend.position = "top",
          panel.grid.major.x = element_line(colour = "#E8E8E8", size = 0.2),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold")
        ) +
        labs(
          x = NULL,
          y = "Spearman rho vs original set",
          title = title_txt,
          caption = paste(
            "Bars compare the leave-one-out profile with the original residual-set profile.",
            "Dashed line, if shown, marks full-module vs residual-set rho.",
            "NA = <3 genes or near-constant profile."
          )
        )
      
      if (is.finite(full_rho_ref)) {
        p <- p + geom_hline(
          yintercept = full_rho_ref,
          linetype = "dashed",
          colour = "#3D8DAE",
          linewidth = 0.35
        )
      }
      
      fn <- if (mod == "Cell cycle") "fig_ablation_residual_cellcycle" else "fig_ablation_residual_neutrophil"
      ggsave(file.path(FIG_DIR, paste0(fn, ".pdf")), p, width=6.5, height=4, device=cairo_pdf)
      ggsave(file.path(FIG_DIR, paste0(fn, ".png")), p, width=6.5, height=4, dpi=300)
    }
  } else {
    warning("Skipping figures: ggplot2 not available.")
  }
  
  # ---------------------------------------------------------------------------
  # 6. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  write_csv_atomic(data.frame(
    key=c("script","script_version","frozen_file","state_file","long_file",
          "min_rho_genes","rho_evaluable_sd","n_donors","n_ablation_rows","R_version"),
    value=c("05_signature_ablation_source_retention.R","v1.3",
            FROZEN_FILE,STATE_FILE,LONG_FILE,
            as.character(MIN_RHO_GENES), as.character(RHO_EVALUABLE_SD),
            as.character(length(donors)), as.character(nrow(abl)),
            as.character(getRversion())),
    stringsAsFactors=FALSE), file.path(TABLE_DIR,"run_provenance.csv"))
  
  cat("Script 05a v1.3 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("\n--- Ablation set summary ---\n")
  print(set_sum, row.names=FALSE)
  cat("\n--- Ablation detail (top-1 flips only) ---\n")
  flips <- abl[!is.na(abl$top1_flipped) & abl$top1_flipped,
               c("module","set","removed_gene","top1_original","top1_leave_one_out","rho_vs_original_set"),
               drop=FALSE]
  if (nrow(flips)) print(flips, row.names=FALSE) else cat("No top-1 source flips under leave-one-out.\n")
  cat("\nInterpretation: ablation is a robustness diagnostic for source\n")
  cat("localization, not evidence of causal gene importance. rho_vs_original_set\n")
  cat("compares the leave-one-out profile to the same set's original profile.\n")
})




# ============================================================================
# SCRIPT 05b-1 v1.3 -- Candidate TF donor-aware detectability precheck
#   Purpose: before any in silico TF perturbation, check whether candidate
#   granulopoietic TFs are detectable in the neutrophil-axis states where the
#   perturbation would apply. Whole-blood scRNA frequently under-detects
#   bone-marrow TFs; this gate decides whether CellOracle is feasible at all.
#   Acute primary cells (26 donors, Script 02 cohort). Donor-aware, >=20 cells.
#   logCPM = donor-state pseudobulk log2(sum_UMI / lib * 1e6 + 1).
#   Detection rate = per-cell fraction (median across eligible donors).
#   Entry threshold: eligible_donor_n >= 5 AND median_detection_rate >= 0.10
#   AND median pseudobulk logCPM > 0 in a target state.
#   Outputs:
#     tf_detectability_screen.csv
#     tf_perturbation_eligibility.csv
#     tf_detectability_passed_pairs.csv
# RUN: Rscript 05b_tf_detectability_precheck.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  PROJECT_ROOT <- "/home/sunshine/predicate/singlecell"
  INPUT_RDS    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  OUTPUTS_ROOT <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  
  STAGE_NAME  <- "05b_tf_detectability"
  RUN_PURPOSE <- "tf_detectability_precheck"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  MIN_CELLS_PER_DONOR_STATE <- 20L
  
  EXPECTED_N_ACUTE_DONORS <- 26L
  EXPECTED_N_PRIMARY_CELLS <- 151837L
  
  CANDIDATE_TFS <- c("CEBPA","CEBPB","CEBPE","GFI1","STAT3","SPI1","IRF8")
  
  TARGET_STATES <- c(
    "Cycling_neutrophil_progenitors",
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "S100A8-9_hi_neutrophils",
    "Mature_neutrophils",
    "Classical_monocytes"
  )
  
  MIN_ELIGIBLE_DONOR <- 5L
  MIN_DETECTION_RATE <- 0.10
  MIN_MEDIAN_LOGCPM  <- 0.0
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  safe_median <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else median(x)}
  write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
    write.csv(x, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) { unlink(tmp); stop_msg("Failed to move temp CSV to final path: ", path) }
    invisible(path)
  }
  find_latest_file <- function(root, pattern) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    fs <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (!length(fs)) stop_msg("No file found under ", root, " matching: ", pattern)
    fs[which.max(file.info(fs)$mtime)]
  }
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
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  needed <- c("SeuratObject","Matrix")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
  
  # ---------------------------------------------------------------------------
  # 2. INPUTS
  # ---------------------------------------------------------------------------
  STATE_FILE <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L)
  missing_states <- setdiff(TARGET_STATES, fine_states)
  if (length(missing_states)) stop_msg("Target states not in fine_states: ", paste(missing_states, collapse=", "))
  
  if (!file.exists(INPUT_RDS)) stop_msg("Input RDS not found: ", INPUT_RDS)
  obj <- readRDS(INPUT_RDS)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay=COUNT_ASSAY, layer=COUNT_LAYER, fast=FALSE)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    if ("layer" %in% names(formals(GetAssayData))) cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, layer="counts")
    else cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, slot="counts")
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  if (!identical(colnames(cnt), colnames(obj))) stop_msg("Count colnames != object cells")
  
  # ---------------------------------------------------------------------------
  # 3. ACUTE PRIMARY COHORT (Script 02 logic)
  # ---------------------------------------------------------------------------
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  
  sample_tab <- unique(data.frame(sample_id=cell_sample, diagnosis=cell_diag, stringsAsFactors=FALSE))
  sample_diag_n <- tapply(sample_tab$diagnosis, sample_tab$sample_id,
                          function(x) length(unique(x)))
  if (any(sample_diag_n != 1L)) {
    stop_msg("Some samples have multiple diagnosis values: ",
             paste(names(sample_diag_n)[sample_diag_n != 1L], collapse = ", "))
  }
  
  sample_tab$donor <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$is_acute <- sample_tab$diagnosis %in% ACUTE_DIAG
  acute_donors <- unique(sample_tab$donor[sample_tab$is_acute])
  stopifnot(length(acute_donors) == EXPECTED_N_ACUTE_DONORS)
  
  meta2 <- meta
  meta2$donor <- sub("_CONV$", "", cell_sample)
  meta2$is_acute <- meta2$donor %in% acute_donors & cell_diag %in% ACUTE_DIAG
  primary_idx <- which(meta2$is_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  
  cnt_p <- cnt[, primary_idx, drop=FALSE]
  meta_p <- meta2[primary_idx, , drop=FALSE]
  cell_state_p <- as.character(meta_p[[FINE_STATE_FIELD]])
  donor_p <- factor(as.character(meta_p$donor), levels=acute_donors)
  state_p <- factor(cell_state_p, levels=fine_states)
  stopifnot(!anyNA(donor_p), !anyNA(state_p))
  
  n_donor <- length(acute_donors); n_state <- length(fine_states)
  ds_idx <- (as.integer(donor_p)-1L)*n_state + as.integer(state_p)
  D <- Matrix::sparseMatrix(i=ds_idx, j=seq_along(ds_idx), x=1,
                            dims=c(n_donor*n_state, length(ds_idx)))
  
  n_cells_ds <- as.integer(Matrix::rowSums(D))
  n_cells_mat <- matrix(
    n_cells_ds,
    nrow = n_donor,
    ncol = n_state,
    byrow = TRUE,
    dimnames = list(acute_donors, fine_states)
  )
  
  lib_cell <- as.numeric(Matrix::colSums(cnt_p))
  if (any(!is.finite(lib_cell) | lib_cell <= 0)) {
    stop_msg("Some primary cells have non-positive library size.")
  }
  
  lib_ds <- as.numeric(D %*% lib_cell)
  lib_mat <- matrix(
    lib_ds,
    nrow = n_donor,
    ncol = n_state,
    byrow = TRUE,
    dimnames = list(acute_donors, fine_states)
  )
  
  eligible_mat <- n_cells_mat >= MIN_CELLS_PER_DONOR_STATE
  
  # ---------------------------------------------------------------------------
  # 4. TF DETECTABILITY
  # ---------------------------------------------------------------------------
  feat <- rownames(cnt)
  tf_missing <- setdiff(CANDIDATE_TFS, feat)
  if (length(tf_missing)) warning("TFs absent from count matrix: ", paste(tf_missing, collapse=", "))
  tf_present <- intersect(CANDIDATE_TFS, feat)
  if (!length(tf_present)) {
    stop_msg("None of the candidate TFs are present in the count matrix.")
  }
  gi <- match(tf_present, feat)
  
  nfeat_cell <- as.numeric(Matrix::colSums(cnt_p > 0))
  
  rows <- list()
  for (i in seq_along(tf_present)) {
    g <- tf_present[i]
    tf_umi <- as.numeric(cnt_p[gi[i], ])
    det_cell <- as.numeric(tf_umi > 0)
    
    tf_umi_ds <- as.numeric(D %*% tf_umi)
    tf_logcpm_ds <- log2(tf_umi_ds / lib_ds * 1e6 + 1)
    tf_logcpm_ds[!is.finite(tf_logcpm_ds) | lib_ds <= 0] <- NA_real_
    
    tf_logcpm_mat <- matrix(
      tf_logcpm_ds,
      nrow = n_donor,
      ncol = n_state,
      byrow = TRUE,
      dimnames = list(acute_donors, fine_states)
    )
    
    det_ds <- as.numeric(D %*% det_cell)
    det_rate_ds <- det_ds / n_cells_ds
    det_rate_ds[n_cells_ds == 0L] <- NA_real_
    
    det_rate_mat <- matrix(
      det_rate_ds,
      nrow = n_donor,
      ncol = n_state,
      byrow = TRUE,
      dimnames = list(acute_donors, fine_states)
    )
    
    med_nf_ds <- rep(NA_real_, n_donor * n_state)
    tf_nf <- tapply(nfeat_cell, ds_idx, safe_median)
    med_nf_ds[as.integer(names(tf_nf))] <- as.numeric(tf_nf)
    
    med_numi_ds <- rep(NA_real_, n_donor * n_state)
    tf_numi <- tapply(lib_cell, ds_idx, safe_median)
    med_numi_ds[as.integer(names(tf_numi))] <- as.numeric(tf_numi)
    
    med_nf_mat <- matrix(
      med_nf_ds,
      nrow = n_donor,
      ncol = n_state,
      byrow = TRUE,
      dimnames = list(acute_donors, fine_states)
    )
    
    med_numi_mat <- matrix(
      med_numi_ds,
      nrow = n_donor,
      ncol = n_state,
      byrow = TRUE,
      dimnames = list(acute_donors, fine_states)
    )
    
    for (s in TARGET_STATES) {
      si <- match(s, fine_states)
      elig <- eligible_mat[, si]
      
      n_elig <- sum(elig)
      med_det <- safe_median(det_rate_mat[elig, si])
      med_lc  <- safe_median(tf_logcpm_mat[elig, si])
      med_nf  <- safe_median(med_nf_mat[elig, si])
      med_nu  <- safe_median(med_numi_mat[elig, si])
      med_nc  <- safe_median(n_cells_mat[elig, si])
      
      pass <- n_elig >= MIN_ELIGIBLE_DONOR &&
        is.finite(med_det) && med_det >= MIN_DETECTION_RATE &&
        is.finite(med_lc) && med_lc > MIN_MEDIAN_LOGCPM
      
      rows[[length(rows) + 1L]] <- data.frame(
        tf = g,
        fine_state = s,
        eligible_donor_n = n_elig,
        eligible_donor_rate = n_elig / n_donor,
        median_detection_rate = med_det,
        median_logCPM = med_lc,
        median_n_cells = med_nc,
        median_nUMI = med_nu,
        median_nFeature = med_nf,
        eligibility = if (pass) "pass" else "insufficient_detectability",
        stringsAsFactors = FALSE
      )
    }
  }
  screen <- do.call(rbind, rows)
  write_csv_atomic(screen, file.path(TABLE_DIR, "tf_detectability_screen.csv"))
  
  primary_states <- c("MPO+_immature_neutrophils_or_progenitors",
                      "PADI4+_immature_neutrophils")
  tf_elig <- do.call(rbind, lapply(tf_present, function(g) {
    z <- screen[screen$tf == g, , drop=FALSE]
    pass_states <- z$fine_state[z$eligibility == "pass"]
    n_pass <- length(pass_states)
    n_pass_primary <- sum(pass_states %in% primary_states)
    
    overall <- if (n_pass_primary >= 1L) {
      "primary_immature_detectable"
    } else if (n_pass >= 1L) {
      "detectable_nonprimary_only"
    } else {
      "insufficient"
    }
    
    data.frame(
      tf = g,
      n_target_states = nrow(z),
      n_pass_states = n_pass,
      n_pass_primary_immature = n_pass_primary,
      pass_states = if (n_pass) paste(pass_states, collapse = ";") else "",
      overall = overall,
      stringsAsFactors=FALSE
    )
  }))
  write_csv_atomic(tf_elig, file.path(TABLE_DIR, "tf_perturbation_eligibility.csv"))
  
  write_csv_atomic(data.frame(
    key=c("script","script_version","input_rds","state_file",
          "min_eligible_donor","min_detection_rate","min_logcpm",
          "n_acute_donors","n_primary_cells","n_tf_present","n_tf_missing",
          "R_version","SeuratObject_version","Matrix_version"),
    value=c("05b_tf_detectability_precheck.R","v1.3",INPUT_RDS,STATE_FILE,
            as.character(MIN_ELIGIBLE_DONOR), as.character(MIN_DETECTION_RATE),
            as.character(MIN_MEDIAN_LOGCPM),
            as.character(length(acute_donors)), as.character(length(primary_idx)),
            as.character(length(tf_present)), as.character(length(tf_missing)),
            as.character(getRversion()),
            as.character(utils::packageVersion("SeuratObject")),
            as.character(utils::packageVersion("Matrix"))),
    stringsAsFactors=FALSE), file.path(TABLE_DIR,"run_provenance.csv"))
  
  cat("Script 05b-1 v1.3 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("TFs checked:", paste(tf_present, collapse=", "),
      if (length(tf_missing)) paste0(" | ABSENT: ", paste(tf_missing, collapse=", ")) else "", "\n")
  cat("\n--- Per-TF eligibility ---\n")
  print(tf_elig, row.names=FALSE)
  cat("\n--- Detectability screen (key immature states) ---\n")
  print(screen[screen$fine_state %in% primary_states,
               c("tf","fine_state","eligible_donor_n","median_detection_rate",
                 "median_logCPM","median_nUMI","median_nFeature","eligibility")],
        row.names = FALSE)
  cat("\n--- Passed TF-state pairs (all target states) ---\n")
  passed <- screen[screen$eligibility == "pass",
                   c("tf","fine_state","eligible_donor_n","median_detection_rate",
                     "median_logCPM","median_nUMI","median_nFeature","eligibility")]
  
  write_csv_atomic(
    passed,
    file.path(TABLE_DIR, "tf_detectability_passed_pairs.csv")
  )
  
  if (nrow(passed)) print(passed, row.names = FALSE) else cat("No TF-state pair passed the detectability gate.\n")
  cat("\nGate: eligible_donor_n >= ", MIN_ELIGIBLE_DONOR,
      " & median_detection_rate >= ", MIN_DETECTION_RATE,
      " & median pseudobulk logCPM > ", MIN_MEDIAN_LOGCPM, ".\n", sep="")
})




# ============================================================================
# SCRIPT 05b-2 v1.4 -- Donor-level TF-program co-expression structure
#   Purpose: weak prioritization screen BEFORE formal GRN perturbation.
#   For each TF x state x readout, correlate donor-level TF pseudobulk logCPM
#   against readout gene-set pseudobulk logCPM across eligible donors.
#
#   THIS IS NOT A GRN EDGE / REGULATION TEST. Donor-level pseudobulk
#   correlations share donor-state library-size / capture confounders, so the
#   output is a descriptive co-expression STRUCTURE, used only to prioritize
#   which TF-program axes to carry into CellOracle / motif-informed GRN.
#
#   Acute primary cells (26 donors, Script 02 cohort). Donor-aware, >=20 cells.
#   logCPM = donor-state pseudobulk log2(sum_UMI / lib * 1e6 + 1), consistent
#   with Script 02/03/05b-1. Readout logCPM = summed-UMI pseudobulk over the
#   member genes (module treated as one super-gene), same formula. NOTE: sums
#   scale with gene count, so readout absolute logCPM is NOT comparable across
#   readout sets of different sizes.
#
#   Readout gene sets follow Script 03b definitions (verified against source):
#     M2_full               = primary_azurophilic (10 genes)
#     primary_core          = primary_azurophilic_core (MPO/ELANE/CTSG/DEFA4)
#     secondary_granule     = secondary_specific_granule (LTF/LCN2)
#     tertiary_granule      = tertiary_gelatinase_granule (MMP8/MMP9)
#     mature_surface        = surface_migration (FCGR3B/CXCR2/CXCR1/FPR1/MME)
#     inflammatory_S100     = S100 (S100A8/S100A9/S100A12)
#     negative_control_hk   = housekeeping (NOT a biological readout; used to
#                             calibrate the descriptive band baseline against
#                             donor-level technical confounders)
#   Missing member genes are dropped with a warning; a readout with zero
#   present genes fails closed.
#
#   v1.4 changes (vs v1.3): TF median-detection floor
#   (MIN_TF_MED_DET_FOR_COR = 0.05). TFs whose median detection rate in a
#   state is below the floor (e.g. IRF8 in the immature states) become
#   not_evaluable_low_tf_detection instead of generating detection-floor
#   noise correlations. New column tf_detection_floor_pass. No other
#   thresholds changed.
#
#   Outputs (tables/):
#     tf_program_coexpression_structure.csv   (long table, all cells)
#     structure_matrix_rho_<state>.csv        (tf x readout rho, key states)
#     structure_matrix_label_<state>.csv      (tf x readout descriptive label)
#     readout_gene_membership.csv             (every readout gene, presence)
#     readout_gene_membership_summary.csv     (n present / total / missing)
#     run_provenance.csv
#   RUN: Rscript 05b_tf_program_coexpression.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  PROJECT_ROOT <- "/home/sunshine/predicate/singlecell"
  INPUT_RDS    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  OUTPUTS_ROOT <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  
  STAGE_NAME  <- "05b_tf_program_coexpression"
  RUN_PURPOSE <- "coexpression_structure"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  MIN_CELLS_PER_DONOR_STATE <- 20L
  
  EXPECTED_N_ACUTE_DONORS  <- 26L
  EXPECTED_N_PRIMARY_CELLS <- 151837L
  
  MIN_ELIGIBLE_DONOR <- 5L
  MIN_RANGE_LOGCPM   <- 0.5    # near-constant guard: logCPM range < this -> not_evaluable
  MIN_TF_MED_DET_FOR_COR <- 0.05   # TF detection floor: below -> not_evaluable_low_tf_detection
  
  # TF grouping (biological priors, NOT derived from 05b-1 detectability data)
  TF_GROUPS <- c(
    CEBPA = "primary_axis",
    CEBPE = "granulocytic_differentiation",
    GFI1  = "granulocytic_differentiation",
    CEBPB = "emergency_inflammatory",
    STAT3 = "emergency_inflammatory",
    SPI1  = "pan_myeloid_background",
    IRF8  = "monocyte_comparator"
  )
  CANDIDATE_TFS <- names(TF_GROUPS)
  
  # Readout gene sets (Script 03b definitions)
  READOUT_SETS <- list(
    M2_full             = c("CEACAM6","CEACAM8","CTSG","DEFA4","ELANE","MPO",
                            "MS4A3","OLFM4","RNASE3","TCN1"),
    primary_core        = c("MPO","ELANE","CTSG","DEFA4"),
    secondary_granule   = c("LTF","LCN2"),
    tertiary_granule    = c("MMP8","MMP9"),
    mature_surface      = c("FCGR3B","CXCR2","CXCR1","FPR1","MME"),
    inflammatory_S100   = c("S100A8","S100A9","S100A12"),
    negative_control_hk = c("ACTB","GAPDH","B2M","TPT1","FTL","FTH1")
  )
  READOUT_GROUPS <- c(
    M2_full             = "m2_overall",
    primary_core        = "primary_granule",
    secondary_granule   = "granule_maturation",
    tertiary_granule    = "granule_maturation",
    mature_surface      = "surface_migration",
    inflammatory_S100   = "inflammatory",
    negative_control_hk = "control"
  )
  
  TARGET_STATES <- c(
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "S100A8-9_hi_neutrophils",
    "Mature_neutrophils",
    "Classical_monocytes"
  )
  PRIMARY_STATES <- c("MPO+_immature_neutrophils_or_progenitors",
                      "PADI4+_immature_neutrophils")
  
  # optional donor-label permutation (0 = disable). Descriptive only: report
  # observed rho percentile relative to the permutation null, no formal P.
  PERM_N    <- 1000L
  PERM_SEED <- 20260804L
  
  # descriptive bands (NOT significance thresholds)
  RHO_POSITIVE <- 0.4
  RHO_WEAK     <- 0.2
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  safe_median <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else median(x)}
  write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
    write.csv(x, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) { unlink(tmp); stop_msg("Failed to move temp CSV to final path: ", path) }
    invisible(path)
  }
  find_latest_file <- function(root, pattern) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    fs <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (!length(fs)) stop_msg("No file found under ", root, " matching: ", pattern)
    fs[which.max(file.info(fs)$mtime)]
  }
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
  near_constant <- function(x) {
    x <- x[is.finite(x)]
    length(x) < 2L || (max(x) - min(x)) < MIN_RANGE_LOGCPM
  }
  structure_label <- function(r) {
    if (!is.finite(r)) return("not_evaluable")
    if (r >= RHO_POSITIVE) return("positive_structure")
    if (r >= RHO_WEAK)     return("weak_positive_structure")
    if (r > -RHO_WEAK)     return("near_null_structure")
    "inverse_structure"
  }
  sanitize_fn <- function(s) {
    s <- gsub("\\+", "plus", s)
    s <- gsub("[^A-Za-z0-9_\\-]+", "_", s)
    s <- gsub("_+", "_", s)
    s
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  needed <- c("SeuratObject","Matrix")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
  
  # ---------------------------------------------------------------------------
  # 2. INPUTS
  # ---------------------------------------------------------------------------
  STATE_FILE <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L)
  missing_states <- setdiff(TARGET_STATES, fine_states)
  if (length(missing_states)) stop_msg("Target states not in fine_states: ", paste(missing_states, collapse=", "))
  
  if (!file.exists(INPUT_RDS)) stop_msg("Input RDS not found: ", INPUT_RDS)
  obj <- readRDS(INPUT_RDS)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay=COUNT_ASSAY, layer=COUNT_LAYER, fast=FALSE)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    if ("layer" %in% names(formals(GetAssayData))) cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, layer="counts")
    else cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, slot="counts")
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  if (!identical(colnames(cnt), colnames(obj))) stop_msg("Count colnames != object cells")
  
  # ---------------------------------------------------------------------------
  # 3. ACUTE PRIMARY COHORT (Script 02 logic)
  # ---------------------------------------------------------------------------
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  
  sample_tab <- unique(data.frame(sample_id=cell_sample, diagnosis=cell_diag, stringsAsFactors=FALSE))
  sample_diag_n <- tapply(sample_tab$diagnosis, sample_tab$sample_id,
                          function(x) length(unique(x)))
  if (any(sample_diag_n != 1L)) {
    stop_msg("Some samples have multiple diagnosis values: ",
             paste(names(sample_diag_n)[sample_diag_n != 1L], collapse = ", "))
  }
  
  sample_tab$donor <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$is_acute <- sample_tab$diagnosis %in% ACUTE_DIAG
  acute_donors <- unique(sample_tab$donor[sample_tab$is_acute])
  stopifnot(length(acute_donors) == EXPECTED_N_ACUTE_DONORS)
  
  meta2 <- meta
  meta2$donor <- sub("_CONV$", "", cell_sample)
  meta2$is_acute <- meta2$donor %in% acute_donors & cell_diag %in% ACUTE_DIAG
  primary_idx <- which(meta2$is_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  stopifnot(all(!grepl("_CONV$", cell_sample[primary_idx])))
  
  cnt_p <- cnt[, primary_idx, drop=FALSE]
  meta_p <- meta2[primary_idx, , drop=FALSE]
  cell_state_p <- as.character(meta_p[[FINE_STATE_FIELD]])
  donor_p <- factor(as.character(meta_p$donor), levels=acute_donors)
  state_p <- factor(cell_state_p, levels=fine_states)
  stopifnot(!anyNA(donor_p), !anyNA(state_p))
  
  n_donor <- length(acute_donors); n_state <- length(fine_states)
  ds_idx <- (as.integer(donor_p)-1L)*n_state + as.integer(state_p)
  D <- Matrix::sparseMatrix(i=ds_idx, j=seq_along(ds_idx), x=1,
                            dims=c(n_donor*n_state, length(ds_idx)))
  
  n_cells_ds <- as.integer(Matrix::rowSums(D))
  n_cells_mat <- matrix(
    n_cells_ds,
    nrow = n_donor,
    ncol = n_state,
    byrow = TRUE,
    dimnames = list(acute_donors, fine_states)
  )
  
  lib_cell <- as.numeric(Matrix::colSums(cnt_p))
  if (any(!is.finite(lib_cell) | lib_cell <= 0)) {
    stop_msg("Some primary cells have non-positive library size.")
  }
  
  lib_ds <- as.numeric(D %*% lib_cell)
  lib_mat <- matrix(
    lib_ds,
    nrow = n_donor,
    ncol = n_state,
    byrow = TRUE,
    dimnames = list(acute_donors, fine_states)
  )
  
  eligible_mat <- n_cells_mat >= MIN_CELLS_PER_DONOR_STATE
  
  # ---------------------------------------------------------------------------
  # 4. GENE-LEVEL donor-state pseudobulk UMI / logCPM / detection rate
  #    (all candidate TFs + all readout member genes, in one pass; ds x gene)
  # ---------------------------------------------------------------------------
  feat <- rownames(cnt)
  tf_missing <- setdiff(CANDIDATE_TFS, feat)
  readout_member_genes <- unique(unlist(READOUT_SETS, use.names = FALSE))
  readout_missing <- setdiff(readout_member_genes, feat)
  all_missing <- unique(c(tf_missing, readout_missing))
  if (length(all_missing)) {
    warning("Genes absent from count matrix (dropped): ", paste(all_missing, collapse=", "))
  }
  tf_present <- intersect(CANDIDATE_TFS, feat)
  if (!length(tf_present)) {
    stop_msg("None of the candidate TFs are present in the count matrix.")
  }
  for (ro in names(READOUT_SETS)) {
    miss <- setdiff(READOUT_SETS[[ro]], feat)
    if (length(miss)) {
      warning("Readout set ", ro, " has genes absent from count matrix: ",
              paste(miss, collapse = ", "))
    }
    if (!length(intersect(READOUT_SETS[[ro]], feat))) {
      stop_msg("Readout '", ro, "' has no genes present in the count matrix.")
    }
  }
  
  # readout gene membership audit (for the manuscript / review)
  readout_membership <- do.call(rbind, lapply(names(READOUT_SETS), function(ro) {
    data.frame(
      readout = ro,
      readout_group = READOUT_GROUPS[[ro]],
      gene = READOUT_SETS[[ro]],
      present_in_counts = READOUT_SETS[[ro]] %in% feat,
      stringsAsFactors = FALSE
    )
  }))
  readout_summary <- aggregate(
    present_in_counts ~ readout + readout_group,
    data = readout_membership,
    FUN = function(x) sum(as.logical(x))
  )
  names(readout_summary)[names(readout_summary) == "present_in_counts"] <- "n_present_genes"
  readout_summary$n_total_genes <- vapply(
    readout_summary$readout,
    function(ro) length(READOUT_SETS[[ro]]),
    integer(1)
  )
  readout_summary$n_missing_genes <- readout_summary$n_total_genes - readout_summary$n_present_genes
  write_csv_atomic(readout_membership, file.path(TABLE_DIR, "readout_gene_membership.csv"))
  write_csv_atomic(readout_summary, file.path(TABLE_DIR, "readout_gene_membership_summary.csv"))
  
  present_genes <- intersect(unique(c(tf_present, readout_member_genes)), feat)
  gi <- match(present_genes, feat)
  gene_umi <- as.matrix(cnt_p[gi, , drop=FALSE])        # n_genes x n_cells
  gene_umi_ds <- as.matrix(D %*% t(gene_umi))           # n_ds x n_genes
  
  gene_logcpm_ds <- log2(sweep(gene_umi_ds, 1L, lib_ds, "/") * 1e6 + 1)
  gene_logcpm_ds[!is.finite(gene_logcpm_ds)] <- NA_real_
  gene_logcpm_ds[lib_ds <= 0, ] <- NA_real_
  colnames(gene_logcpm_ds) <- present_genes
  
  det_cell <- (gene_umi > 0) * 1                        # n_genes x n_cells (numeric)
  det_ds <- as.matrix(D %*% t(det_cell))                # n_ds x n_genes
  det_rate_ds <- sweep(det_ds, 1L, n_cells_ds, "/")
  det_rate_ds[!is.finite(det_rate_ds)] <- NA_real_
  det_rate_ds[n_cells_ds == 0L, ] <- NA_real_
  colnames(det_rate_ds) <- present_genes
  
  # ---------------------------------------------------------------------------
  # 5. TF x state x readout donor-level co-expression structure
  # ---------------------------------------------------------------------------
  if (PERM_N > 0L) set.seed(PERM_SEED)
  rows <- list()
  for (g in tf_present) {
    gcol <- match(g, present_genes)
    for (s in TARGET_STATES) {
      si <- match(s, fine_states)
      elig <- eligible_mat[, si]
      n_elig <- sum(elig)
      ds_keep <- seq.int(from = si, by = n_state, length.out = n_donor)[elig]
      
      tf_lc  <- gene_logcpm_ds[ds_keep, gcol]
      tf_det <- det_rate_ds[ds_keep, gcol]
      tf_med_det <- safe_median(tf_det)
      tf_med_lc  <- safe_median(tf_lc)
      det_floor_pass <- is.finite(tf_med_det) && tf_med_det >= MIN_TF_MED_DET_FOR_COR
      
      for (ro in names(READOUT_SETS)) {
        ro_genes <- intersect(READOUT_SETS[[ro]], present_genes)
        ro_cols <- match(ro_genes, present_genes)
        ro_umi <- rowSums(gene_umi_ds[ds_keep, ro_cols, drop = FALSE])
        ro_lc <- log2(ro_umi / lib_ds[ds_keep] * 1e6 + 1)
        ro_lc[!is.finite(ro_lc) | lib_ds[ds_keep] <= 0] <- NA_real_
        ro_med <- safe_median(ro_lc)
        
        evaluable <- n_elig >= MIN_ELIGIBLE_DONOR &&
          det_floor_pass &&
          !near_constant(tf_lc) && !near_constant(ro_lc)
        
        if (!evaluable) {
          ne_label <- if (!det_floor_pass) "not_evaluable_low_tf_detection" else "not_evaluable"
          rows[[length(rows) + 1L]] <- data.frame(
            tf = g, tf_group = TF_GROUPS[[g]], fine_state = s, readout = ro,
            readout_group = READOUT_GROUPS[[ro]], n_eligible_donors = n_elig,
            tf_median_detection_rate = tf_med_det, tf_median_logCPM = tf_med_lc,
            tf_detection_floor_pass = det_floor_pass,
            readout_median_logCPM = ro_med,
            rho = NA_real_, structure_label = ne_label,
            perm_n = if (PERM_N > 0L) PERM_N else 0L,
            perm_null_median = NA_real_, perm_null_pctl = NA_real_,
            stringsAsFactors = FALSE)
          next
        }
        
        rho <- suppressWarnings(cor(tf_lc, ro_lc, method = "spearman"))
        lab <- structure_label(rho)
        
        perm_med <- NA_real_; perm_pctl <- NA_real_
        if (PERM_N > 0L) {
          null <- vapply(seq_len(PERM_N), function(k) {
            suppressWarnings(cor(tf_lc, sample(ro_lc), method = "spearman"))
          }, numeric(1))
          null <- null[is.finite(null)]
          if (length(null)) {
            perm_med <- safe_median(null)
            perm_pctl <- 100 * mean(null <= rho)
          }
        }
        
        rows[[length(rows) + 1L]] <- data.frame(
          tf = g, tf_group = TF_GROUPS[[g]], fine_state = s, readout = ro,
          readout_group = READOUT_GROUPS[[ro]], n_eligible_donors = n_elig,
          tf_median_detection_rate = tf_med_det, tf_median_logCPM = tf_med_lc,
          tf_detection_floor_pass = det_floor_pass,
          readout_median_logCPM = ro_med,
          rho = rho, structure_label = lab,
          perm_n = if (PERM_N > 0L) PERM_N else 0L,
          perm_null_median = perm_med, perm_null_pctl = perm_pctl,
          stringsAsFactors = FALSE)
      }
    }
  }
  res <- do.call(rbind, rows)
  write_csv_atomic(res, file.path(TABLE_DIR, "tf_program_coexpression_structure.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. Structure matrices for the key immature states (tf x readout)
  # ---------------------------------------------------------------------------
  for (s in PRIMARY_STATES) {
    z <- res[res$fine_state == s, , drop = FALSE]
    mat <- matrix(NA_real_, nrow = length(CANDIDATE_TFS), ncol = length(READOUT_SETS),
                  dimnames = list(CANDIDATE_TFS, names(READOUT_SETS)))
    lab_mat <- matrix("not_evaluable", nrow = length(CANDIDATE_TFS), ncol = length(READOUT_SETS),
                      dimnames = dimnames(mat))
    for (g in CANDIDATE_TFS) for (ro in names(READOUT_SETS)) {
      w <- which(z$tf == g & z$readout == ro)
      if (length(w)) { mat[g, ro] <- z$rho[w[1L]]; lab_mat[g, ro] <- z$structure_label[w[1L]] }
    }
    fn <- sanitize_fn(s)
    rho_out <- data.frame(tf = rownames(mat), mat, check.names = FALSE)
    lab_out <- data.frame(tf = rownames(lab_mat), lab_mat, check.names = FALSE)
    write_csv_atomic(rho_out, file.path(TABLE_DIR, paste0("structure_matrix_rho_", fn, ".csv")))
    write_csv_atomic(lab_out, file.path(TABLE_DIR, paste0("structure_matrix_label_", fn, ".csv")))
  }
  
  write_csv_atomic(data.frame(
    key = c("script","script_version","input_rds","state_file",
            "min_eligible_donor","min_range_logcpm","min_tf_med_det_for_cor",
            "perm_n","perm_seed",
            "n_acute_donors","n_primary_cells",
            "n_tf_present","n_tf_missing","tf_missing",
            "readout_member_genes_missing",
            "readout_logcpm_note","R_version"),
    value = c("05b_tf_program_coexpression.R","v1.4",INPUT_RDS,STATE_FILE,
              as.character(MIN_ELIGIBLE_DONOR), as.character(MIN_RANGE_LOGCPM),
              as.character(MIN_TF_MED_DET_FOR_COR),
              as.character(PERM_N), as.character(PERM_SEED),
              as.character(length(acute_donors)), as.character(length(primary_idx)),
              as.character(length(tf_present)), as.character(length(tf_missing)),
              if (length(tf_missing)) paste(tf_missing, collapse=";") else "none",
              if (length(readout_missing)) paste(readout_missing, collapse=";") else "none",
              "Readout logCPM is summed-UMI pseudobulk; absolute values not directly comparable across readout sets
  with different gene counts.",
              as.character(getRversion())),
    stringsAsFactors = FALSE), file.path(TABLE_DIR, "run_provenance.csv"))
  
  # ---------------------------------------------------------------------------
  # 7. CONSOLE SUMMARY
  # ---------------------------------------------------------------------------
  cat("Script 05b-2 v1.4 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("TFs:", paste(tf_present, collapse=", "),
      if (length(tf_missing)) paste0(" | ABSENT: ", paste(tf_missing, collapse=", ")) else "", "\n")
  if (length(readout_missing)) cat("Readout genes absent:", paste(readout_missing, collapse=", "), "\n")
  
  cat("\n--- Readout gene membership summary ---\n")
  print(readout_summary, row.names = FALSE)
  
  cat("\n--- TF grouping ---\n")
  print(data.frame(tf = names(TF_GROUPS), group = unname(TF_GROUPS), stringsAsFactors = FALSE), row.names = FALSE)
  
  cat("\n--- Donor-level TF-program co-expression structure (Spearman rho; NE = not_evaluable) ---\n")
  for (s in PRIMARY_STATES) {
    z <- res[res$fine_state == s, , drop = FALSE]
    mat <- matrix(NA_real_, nrow = length(CANDIDATE_TFS), ncol = length(READOUT_SETS),
                  dimnames = list(CANDIDATE_TFS, names(READOUT_SETS)))
    for (g in CANDIDATE_TFS) for (ro in names(READOUT_SETS)) {
      w <- which(z$tf == g & z$readout == ro)
      if (length(w)) mat[g, ro] <- z$rho[w[1L]]
    }
    cat("\n", s, "\n", sep = "")
    print(round(mat, 2))
  }
  
  cat("\n--- Per-TF dominant readout in immature states (max |rho| among evaluable) ---\n")
  imm <- res[res$fine_state %in% PRIMARY_STATES & is.finite(res$rho), , drop = FALSE]
  if (nrow(imm)) {
    for (g in tf_present) {
      z <- imm[imm$tf == g, , drop = FALSE]
      if (!nrow(z)) { cat(g, ": no evaluable structure in immature states\n"); next }
      z <- z[order(-abs(z$rho)), , drop = FALSE]
      cat(sprintf("%-6s | top: %s (rho=%.2f, %s, n=%d)", g,
                  z$readout[1L], z$rho[1L], z$structure_label[1L], z$n_eligible_donors[1L]))
      if (nrow(z) >= 2L) cat(sprintf("; 2nd: %s (rho=%.2f, %s)", z$readout[2L], z$rho[2L], z$structure_label[2L]))
      cat("\n")
    }
  } else {
    cat("No evaluable structure in immature states.\n")
  }
  
  cat("\nNot-evaluable cells: ", sum(res$structure_label == "not_evaluable"),
      " (general), ", sum(res$structure_label == "not_evaluable_low_tf_detection"),
      " (low TF detection).\n", sep = "")
  cat("\nGate: eligible_donor_n >= ", MIN_ELIGIBLE_DONOR,
      " & TF median detection >= ", MIN_TF_MED_DET_FOR_COR,
      " & non-constant TF/readout (logCPM range >= ", MIN_RANGE_LOGCPM, ").\n", sep = "")
  cat("Note: rho bands are descriptive only; not a significance or regulation test.\n")
  cat("Note: readout logCPM is summed-UMI pseudobulk; absolute values are not directly\n")
  cat("comparable across readout sets with different gene counts.\n")
})



# ============================================================================
# SCRIPT 05b2b v1.3 -- TF-program partial-correlation stress test
#   Purpose: stress-test the most salvageable positive structure from 05b-2 --
#   CEBPB x primary_core in MPO+ immature (raw donor-level rho 0.52).
#   Question: does the positive coupling survive adjustment for donor-state
#   library size, or is it an artifact of shared X/library-size normalization?
#
#   Positioning: TECHNICAL ROBUSTNESS / REGULATORY PLAUSIBILITY CHECK.
#   NOT new evidence, NOT a regulation test, NOT a perturbation.
#
#   Method: rank-residual partial Spearman with ONE covariate:
#       cov = log10(donor_state_library_size + 1)
#   n_cells is deliberately excluded (highly collinear with library size;
#   n = 9-11 -> adding covariates destabilizes).
#
#   Interpretation (descriptive only, no significance):
#     non-finite raw           -> not_evaluable
#     raw > 0 & partial >= 0.30   -> retained_after_library_adjustment
#     raw > 0 & partial in [0,0.30) -> attenuated_library_sensitive
#     raw > 0 & partial < 0       -> do_not_cite_as_positive_structure
#     raw <= 0 (finite)           -> raw_not_positive (partial for completeness)
#     covariate constant / n small -> not_evaluable[_covariate_constant]
#     TF below detection floor   -> not_evaluable_low_tf_detection
#
#   v1.1: n_complete transparency; provenance split; robust crosscheck column;
#         key-pair interpretation.
#   v1.2: crosscheck via merge with NA-safe column selection (fixes
#         length(NA)==1 false-positive); provenance written after verdicts;
#         key filter + fail-closed column check.
#   v1.3: non-finite raw -> not_evaluable (not raw_not_positive); non-evaluable
#         rows retain n_complete_pairs; provenance splits raw-/partial-
#         evaluable counts.
#
#   Outputs (tables/):
#     tf_program_partial_correlation_stress_test.csv
#     key_tf_primary_core_partial_stress_test.csv
#     run_provenance.csv
#   RUN: Rscript 05b2b_tf_program_partial_correlation_stress_test.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  SCRIPT_NAME    <- "05b2b_tf_program_partial_correlation_stress_test.R"
  SCRIPT_VERSION <- "v1.3"
  
  PROJECT_ROOT  <- "/home/sunshine/predicate/singlecell"
  INPUT_FILE    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  OUTPUTS_ROOT  <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  SCRIPT05B2_ROOT <- file.path(OUTPUTS_ROOT, "05b_tf_program_coexpression")
  
  STAGE_NAME  <- "05b2b_tf_program_partial_correlation"
  RUN_PURPOSE <- "library_size_partial_correlation_stress_test"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  MIN_CELLS_PER_DONOR_STATE <- 20L
  
  EXPECTED_N_ACUTE_DONORS  <- 26L
  EXPECTED_N_PRIMARY_CELLS <- 151837L
  
  MIN_ELIGIBLE_DONOR  <- 5L
  MIN_RANGE_LOGCPM    <- 0.5      # near-constant guard (same as 05b-2)
  MIN_TF_MED_DET_FOR_COR <- 0.05  # TF detection floor (same as 05b-2)
  PARTIAL_RETAIN_THRESH <- 0.30   # descriptive retention band
  
  # Identical definitions to 05b-2 v1.4 (do not drift).
  TF_GROUPS <- c(
    CEBPA = "primary_axis",
    CEBPE = "granulocytic_differentiation",
    GFI1  = "granulocytic_differentiation",
    CEBPB = "emergency_inflammatory",
    STAT3 = "emergency_inflammatory",
    SPI1  = "pan_myeloid_background",
    IRF8  = "monocyte_comparator"
  )
  CANDIDATE_TFS <- names(TF_GROUPS)
  
  READOUT_SETS <- list(
    M2_full             = c("CEACAM6","CEACAM8","CTSG","DEFA4","ELANE","MPO",
                            "MS4A3","OLFM4","RNASE3","TCN1"),
    primary_core        = c("MPO","ELANE","CTSG","DEFA4"),
    secondary_granule   = c("LTF","LCN2"),
    tertiary_granule    = c("MMP8","MMP9"),
    mature_surface      = c("FCGR3B","CXCR2","CXCR1","FPR1","MME"),
    inflammatory_S100   = c("S100A8","S100A9","S100A12"),
    negative_control_hk = c("ACTB","GAPDH","B2M","TPT1","FTL","FTH1")
  )
  READOUT_GROUPS <- c(
    M2_full             = "m2_overall",
    primary_core        = "primary_granule",
    secondary_granule   = "granule_maturation",
    tertiary_granule    = "granule_maturation",
    mature_surface      = "surface_migration",
    inflammatory_S100   = "inflammatory",
    negative_control_hk = "control"
  )
  
  TARGET_STATES <- c(
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "S100A8-9_hi_neutrophils",
    "Mature_neutrophils",
    "Classical_monocytes"
  )
  KEY_TF      <- "CEBPB"
  KEY_READOUT <- "primary_core"
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  safe_median <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else median(x)}
  write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
    write.csv(x, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) { unlink(tmp); stop_msg("Failed to move temp CSV to final path: ", path) }
    invisible(path)
  }
  find_latest_file <- function(root, pattern) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    fs <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (!length(fs)) stop_msg("No file found under ", root, " matching: ", pattern)
    fs[which.max(file.info(fs)$mtime)]
  }
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
  near_constant <- function(x) {
    x <- x[is.finite(x)]
    length(x) < 2L || (max(x) - min(x)) < MIN_RANGE_LOGCPM
  }
  # Rank-residual partial Spearman with ONE covariate; always reports n used.
  partial_rho_library <- function(x, y, lib, min_n = MIN_ELIGIBLE_DONOR) {
    ok <- is.finite(x) & is.finite(y) & is.finite(lib)
    n_complete <- sum(ok)
    if (n_complete < min_n) {
      return(list(rho = NA_real_, covariate_varies = NA, n_complete = n_complete))
    }
    rx <- rank(x[ok], ties.method = "average")
    ry <- rank(y[ok], ties.method = "average")
    rc <- rank(lib[ok], ties.method = "average")
    if (length(unique(rc)) < 3L) {
      return(list(rho = NA_real_, covariate_varies = FALSE, n_complete = n_complete))
    }
    rxl <- resid(stats::lm(rx ~ rc))
    ryl <- resid(stats::lm(ry ~ rc))
    if (sum(!is.na(rxl) & !is.na(ryl)) < min_n) {
      return(list(rho = NA_real_, covariate_varies = TRUE, n_complete = n_complete))
    }
    list(rho = suppressWarnings(stats::cor(rxl, ryl, method = "pearson")),
         covariate_varies = TRUE, n_complete = n_complete)
  }
  verdict_partial <- function(raw, partial, cov_varies) {
    if (!is.finite(raw)) return("not_evaluable")
    if (is.na(partial)) {
      return(if (identical(cov_varies, FALSE)) "not_evaluable_covariate_constant" else "not_evaluable")
    }
    if (raw <= 0) return("raw_not_positive")
    if (partial >= PARTIAL_RETAIN_THRESH) return("retained_after_library_adjustment")
    if (partial >= 0) return("attenuated_library_sensitive")
    "do_not_cite_as_positive_structure"
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  needed <- c("SeuratObject","Matrix")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
  
  # ---------------------------------------------------------------------------
  # 2. INPUTS
  # ---------------------------------------------------------------------------
  STATE_FILE <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L)
  missing_states <- setdiff(TARGET_STATES, fine_states)
  if (length(missing_states)) stop_msg("Target states not in fine_states: ",
                                       paste(missing_states, collapse=", "))
  
  # Optional 05b-2 cross-check input (fail-open).
  CROSSCHECK_FILE <- NA_character_
  if (dir.exists(SCRIPT05B2_ROOT)) {
    CROSSCHECK_FILE <- tryCatch(
      find_latest_file(SCRIPT05B2_ROOT, "^tf_program_coexpression_structure\\.csv$"),
      error = function(e) NA_character_)
  }
  crosscheck <- NULL
  if (!is.na(CROSSCHECK_FILE) && file.exists(CROSSCHECK_FILE)) {
    crosscheck <- utils::read.csv(CROSSCHECK_FILE, stringsAsFactors = FALSE, check.names = FALSE)
    cat("Cross-check reference:", CROSSCHECK_FILE, "\n")
  } else {
    cat("No 05b-2 coexpression table found; running standalone (cross-check skipped).\n")
  }
  
  if (!file.exists(INPUT_FILE)) stop_msg("Input RDS not found: ", INPUT_FILE)
  obj <- readRDS(INPUT_FILE)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay=COUNT_ASSAY, layer=COUNT_LAYER, fast=FALSE)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    if ("layer" %in% names(formals(GetAssayData))) cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, layer="counts")
    else cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, slot="counts")
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  if (!identical(colnames(cnt), colnames(obj))) stop_msg("Count colnames != object cells")
  
  # ---------------------------------------------------------------------------
  # 3. ACUTE PRIMARY COHORT (identical to 05b-2 section 3)
  # ---------------------------------------------------------------------------
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  
  sample_tab <- unique(data.frame(sample_id=cell_sample, diagnosis=cell_diag, stringsAsFactors=FALSE))
  sample_diag_n <- tapply(sample_tab$diagnosis, sample_tab$sample_id, function(x) length(unique(x)))
  if (any(sample_diag_n != 1L)) {
    stop_msg("Some samples have multiple diagnosis values: ",
             paste(names(sample_diag_n)[sample_diag_n != 1L], collapse = ", "))
  }
  sample_tab$donor <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$is_acute <- sample_tab$diagnosis %in% ACUTE_DIAG
  acute_donors <- unique(sample_tab$donor[sample_tab$is_acute])
  stopifnot(length(acute_donors) == EXPECTED_N_ACUTE_DONORS)
  
  meta2 <- meta
  meta2$donor <- sub("_CONV$", "", cell_sample)
  meta2$is_acute <- meta2$donor %in% acute_donors & cell_diag %in% ACUTE_DIAG
  primary_idx <- which(meta2$is_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  
  cnt_p <- cnt[, primary_idx, drop=FALSE]
  donor_p <- factor(as.character(meta2$donor[primary_idx]), levels=acute_donors)
  state_p <- factor(as.character(meta2[[FINE_STATE_FIELD]][primary_idx]), levels=fine_states)
  stopifnot(!anyNA(donor_p), !anyNA(state_p))
  
  n_donor <- length(acute_donors); n_state <- length(fine_states)
  ds_idx <- (as.integer(donor_p)-1L)*n_state + as.integer(state_p)
  D <- Matrix::sparseMatrix(i=ds_idx, j=seq_along(ds_idx), x=1,
                            dims=c(n_donor*n_state, length(ds_idx)))
  
  n_cells_ds <- as.integer(Matrix::rowSums(D))
  n_cells_mat <- matrix(n_cells_ds, nrow=n_donor, ncol=n_state, byrow=TRUE,
                        dimnames=list(acute_donors, fine_states))
  
  lib_cell <- as.numeric(Matrix::colSums(cnt_p))
  if (any(!is.finite(lib_cell) | lib_cell <= 0)) stop_msg("Non-positive cell library sizes.")
  lib_ds <- as.numeric(D %*% lib_cell)
  lib_mat <- matrix(lib_ds, nrow=n_donor, ncol=n_state, byrow=TRUE,
                    dimnames=list(acute_donors, fine_states))
  
  eligible_mat <- n_cells_mat >= MIN_CELLS_PER_DONOR_STATE
  
  # ---------------------------------------------------------------------------
  # 4. GENE-LEVEL donor-state pseudobulk (TFs + readout members)
  # ---------------------------------------------------------------------------
  feat <- rownames(cnt)
  tf_missing <- setdiff(CANDIDATE_TFS, feat)
  readout_member_genes <- unique(unlist(READOUT_SETS, use.names=FALSE))
  readout_missing <- setdiff(readout_member_genes, feat)
  if (length(tf_missing)) warning("TFs absent (dropped): ", paste(tf_missing, collapse=", "))
  if (length(readout_missing)) warning("Readout genes absent (dropped): ", paste(readout_missing, collapse=", "))
  tf_present <- intersect(CANDIDATE_TFS, feat)
  if (!length(tf_present)) stop_msg("None of the candidate TFs present in the count matrix.")
  for (ro in names(READOUT_SETS)) {
    if (!length(intersect(READOUT_SETS[[ro]], feat))) {
      stop_msg("Readout '", ro, "' has no genes present in the count matrix.")
    }
  }
  
  present_genes <- intersect(unique(c(tf_present, readout_member_genes)), feat)
  gi <- match(present_genes, feat)
  gene_umi <- as.matrix(cnt_p[gi, , drop=FALSE])          # n_genes x n_cells
  gene_umi_ds <- as.matrix(D %*% t(gene_umi))             # n_ds x n_genes
  gene_logcpm_ds <- log2(sweep(gene_umi_ds, 1L, lib_ds, "/") * 1e6 + 1)
  gene_logcpm_ds[!is.finite(gene_logcpm_ds)] <- NA_real_
  gene_logcpm_ds[lib_ds <= 0, ] <- NA_real_
  colnames(gene_logcpm_ds) <- present_genes
  
  det_cell <- (gene_umi > 0) * 1
  det_ds <- as.matrix(D %*% t(det_cell))
  det_rate_ds <- sweep(det_ds, 1L, n_cells_ds, "/")
  det_rate_ds[!is.finite(det_rate_ds)] <- NA_real_
  det_rate_ds[n_cells_ds == 0L, ] <- NA_real_
  colnames(det_rate_ds) <- present_genes
  
  # ---------------------------------------------------------------------------
  # 5. PARTIAL-CORRELATION STRESS TEST (all TF x state x readout pairs)
  # ---------------------------------------------------------------------------
  rows <- list()
  for (g in tf_present) {
    gcol <- match(g, present_genes)
    for (s in TARGET_STATES) {
      si <- match(s, fine_states)
      elig <- eligible_mat[, si]
      n_elig <- sum(elig)
      ds_keep <- seq.int(from = si, by = n_state, length.out = n_donor)[elig]
      
      tf_lc  <- gene_logcpm_ds[ds_keep, gcol]
      tf_det <- det_rate_ds[ds_keep, gcol]
      tf_med_det <- safe_median(tf_det)
      det_floor_pass <- is.finite(tf_med_det) && tf_med_det >= MIN_TF_MED_DET_FOR_COR
      
      cov_lib <- log10(lib_ds[ds_keep] + 1)
      
      for (ro in names(READOUT_SETS)) {
        ro_genes <- intersect(READOUT_SETS[[ro]], present_genes)
        ro_cols <- match(ro_genes, present_genes)
        ro_umi <- rowSums(gene_umi_ds[ds_keep, ro_cols, drop=FALSE])
        ro_lc <- log2(ro_umi / lib_ds[ds_keep] * 1e6 + 1)
        ro_lc[!is.finite(ro_lc) | lib_ds[ds_keep] <= 0] <- NA_real_
        
        n_complete0 <- sum(is.finite(tf_lc) & is.finite(ro_lc) & is.finite(cov_lib))
        
        evaluable <- n_elig >= MIN_ELIGIBLE_DONOR &&
          det_floor_pass &&
          !near_constant(tf_lc) && !near_constant(ro_lc)
        
        if (!evaluable) {
          rows[[length(rows)+1L]] <- data.frame(
            tf = g, tf_group = TF_GROUPS[[g]], fine_state = s, readout = ro,
            readout_group = READOUT_GROUPS[[ro]],
            n_eligible_donors = n_elig,
            n_complete_pairs = n_complete0,
            tf_median_detection_rate = tf_med_det,
            tf_detection_floor_pass = det_floor_pass,
            raw_rho = NA_real_,
            partial_rho_library = NA_real_,
            rho_delta = NA_real_,
            covariate_log10_lib_varies = NA,
            verdict = if (!det_floor_pass) "not_evaluable_low_tf_detection" else "not_evaluable",
            stringsAsFactors = FALSE)
          next
        }
        
        raw_rho <- suppressWarnings(stats::cor(tf_lc, ro_lc, method = "spearman"))
        pr <- partial_rho_library(tf_lc, ro_lc, cov_lib)
        verdict <- verdict_partial(raw_rho, pr$rho, pr$covariate_varies)
        
        rows[[length(rows)+1L]] <- data.frame(
          tf = g, tf_group = TF_GROUPS[[g]], fine_state = s, readout = ro,
          readout_group = READOUT_GROUPS[[ro]],
          n_eligible_donors = n_elig,
          n_complete_pairs = pr$n_complete,
          tf_median_detection_rate = tf_med_det,
          tf_detection_floor_pass = det_floor_pass,
          raw_rho = raw_rho,
          partial_rho_library = pr$rho,
          rho_delta = pr$rho - raw_rho,
          covariate_log10_lib_varies = pr$covariate_varies,
          verdict = verdict,
          stringsAsFactors = FALSE)
      }
    }
  }
  res <- do.call(rbind, rows)
  
  # ---------------------------------------------------------------------------
  # 6. CROSS-CHECK vs 05b-2 raw rho (audit that this rebuild reproduces 05b-2)
  # ---------------------------------------------------------------------------
  rho_candidates <- intersect(c("rho","spearman_rho","raw_rho"), names(crosscheck))
  rho_col <- if (length(rho_candidates)) rho_candidates[1] else NA_character_
  cc_input_ok <- !is.null(crosscheck) &&
    !is.na(rho_col) &&
    all(c("tf","fine_state","readout") %in% names(crosscheck))
  
  if (cc_input_ok) {
    cc <- crosscheck[, c("tf","fine_state","readout", rho_col), drop = FALSE]
    names(cc)[names(cc) == rho_col] <- "rho_05b2"
    res <- merge(res, cc, by = c("tf","fine_state","readout"), all.x = TRUE, sort = FALSE)
    res$raw_rho_delta_vs_05b2 <- res$raw_rho - res$rho_05b2
    matched <- !is.na(res$rho_05b2) & is.finite(res$raw_rho)
    n_match <- sum(matched)
    n_bad <- sum(abs(res$raw_rho_delta_vs_05b2[matched]) >= 0.01)
    if (n_bad > 0L) {
      warning("Raw rho mismatch vs 05b-2 for ", n_bad, "/", n_match,
              " evaluable pairs. Cohort/logCPM reproduction diverged.")
    } else {
      cat("Cross-check: raw rho reproduced 05b-2 for", n_match, "pairs (col '", rho_col, "').\n", sep="")
    }
  } else {
    res$rho_05b2 <- NA_real_
    res$raw_rho_delta_vs_05b2 <- NA_real_
    if (!is.null(crosscheck)) {
      warning("05b-2 crosscheck table found, but rho column or key columns missing; cross-check skipped.")
    }
  }
  
  write_csv_atomic(res, file.path(TABLE_DIR, "tf_program_partial_correlation_stress_test.csv"))
  
  # ---------------------------------------------------------------------------
  # 7. KEY-PAIR OUTPUT (CEBPB x primary_core) with interpretation
  # ---------------------------------------------------------------------------
  interp_of <- c(
    retained_after_library_adjustment =
      "positive donor-level structure remains after adjusting for library size; prioritization cue only",
    attenuated_library_sensitive =
      "positive raw structure attenuates after library adjustment; do not use as independent support",
    do_not_cite_as_positive_structure =
      "positive raw structure reverses after adjustment; do not cite as positive structure",
    raw_not_positive =
      "raw donor-level structure was not positive; partial adjustment is not used as supportive evidence",
    not_evaluable =
      "not evaluable because donor coverage or variability gate failed",
    not_evaluable_covariate_constant =
      "library-size covariate constant among eligible donors; adjustment not meaningful",
    not_evaluable_low_tf_detection =
      "TF below detection floor in this state; not evaluable"
  )
  
  key <- res[res$tf %in% KEY_TF & res$fine_state %in% TARGET_STATES & res$readout %in% KEY_READOUT, , drop = FALSE]
  key$interpretation <- unname(interp_of[key$verdict])
  key$interpretation[is.na(key$interpretation)] <- "unmapped verdict; inspect manually"
  
  key <- key[order(match(key$tf, KEY_TF),
                   match(key$fine_state, TARGET_STATES),
                   match(key$readout, KEY_READOUT)), , drop = FALSE]
  
  key_cols <- c("tf","fine_state","readout","n_eligible_donors","n_complete_pairs",
                "raw_rho","partial_rho_library","rho_delta","verdict","interpretation")
  missing_key_cols <- setdiff(key_cols, names(key))
  if (length(missing_key_cols)) {
    stop_msg("Key-pair output missing columns: ", paste(missing_key_cols, collapse = ", "))
  }
  key_out <- key[, key_cols, drop = FALSE]
  write_csv_atomic(key_out, file.path(TABLE_DIR, "key_tf_primary_core_partial_stress_test.csv"))
  
  # ---------------------------------------------------------------------------
  # 8. PROVENANCE (written AFTER res$verdict exists)
  # ---------------------------------------------------------------------------
  provenance <- data.frame(
    field = c(
      "script", "version", "date", "input_file", "crosscheck_file",
      "covariate", "partial_method", "retain_threshold",
      "n_acute_donors", "n_primary_cells",
      "n_tf_present", "n_tf_missing",
      "n_total_pairs", "n_evaluable_pairs", "n_partial_evaluable_pairs",
      "n_positive_raw_pairs", "n_retained_positive_pairs",
      "R_version"
    ),
    value = c(
      SCRIPT_NAME, SCRIPT_VERSION, as.character(Sys.time()), INPUT_FILE,
      if (is.na(CROSSCHECK_FILE)) "none" else CROSSCHECK_FILE,
      "log10(donor_state_library_size + 1)",
      "rank-residual partial Spearman (single covariate)",
      as.character(PARTIAL_RETAIN_THRESH),
      as.character(length(acute_donors)),
      as.character(length(primary_idx)),
      as.character(length(tf_present)),
      as.character(length(tf_missing)),
      as.character(nrow(res)),
      as.character(sum(is.finite(res$raw_rho))),
      as.character(sum(is.finite(res$partial_rho_library))),
      as.character(sum(is.finite(res$raw_rho) & res$raw_rho > 0)),
      as.character(sum(res$verdict == "retained_after_library_adjustment", na.rm = TRUE)),
      as.character(getRversion())
    ),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(provenance, file.path(TABLE_DIR, "run_provenance.csv"))
  
  # ---------------------------------------------------------------------------
  # 9. CONSOLE SUMMARY
  # ---------------------------------------------------------------------------
  cat("Script 05b2b v1.3 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("Key pair:", KEY_TF, "x", KEY_READOUT, "(across target states)\n")
  print(key_out[, c("fine_state","n_eligible_donors","n_complete_pairs",
                    "raw_rho","partial_rho_library","rho_delta","verdict")],
        row.names = FALSE)
  
  kp <- key_out[key_out$fine_state == "MPO+_immature_neutrophils_or_progenitors", , drop = FALSE]
  if (nrow(kp) == 1L) {
    k <- kp[1, , drop = FALSE]
    if (identical(k$verdict, "retained_after_library_adjustment")) {
      cat("-> CEBPB-primary-core structure retained after library-size adjustment (prioritization cue only).\n")
    } else if (identical(k$verdict, "attenuated_library_sensitive")) {
      cat("-> Apparent CEBPB-primary-core coupling is library-size sensitive; not used as supporting evidence.\n")
    } else if (identical(k$verdict, "do_not_cite_as_positive_structure")) {
      cat("-> Partial rho negative: do not cite CEBPB-primary_core as positive structure.\n")
    } else {
      cat("-> Not evaluable; CEBPB-primary_core not used as supporting evidence.\n")
    }
  }
  
  cat("\n--- Verdict distribution (all evaluable, raw>0) ---\n")
  pos <- res[!is.na(res$raw_rho) & res$raw_rho > 0, , drop = FALSE]
  if (nrow(pos)) print(table(pos$verdict, useNA = "ifany")) else cat("No positive raw structures.\n")
  
  cat("\nNote: descriptive robustness check only. No significance, no regulatory claim.\n")
})



# ============================================================================
# SCRIPT 05d v1.5 -- doublet / ambient technical robustness (M1 sensitivity)
#   Purpose: exclude the two easiest reviewer attacks on the M1 cell-cycle
#   story and the residual granule signal, WITHOUT full re-processing.
#
#   Positioning: TECHNICAL ROBUSTNESS AUDIT, NOT a reprocessing or
#   decontamination pipeline. Tests (1) whether M1 source localization
#   survives plausible doublet-risk exclusion and (2) whether non-granulocyte
#   granule-marker signal tracks donor neutrophil burden as a descriptive
#   ambient-RNA proxy. It does NOT claim decontamination.
#
#   L1 doublet -- author score/prediction used if present in the object
#     metadata (the study ran Scrublet + Doublet Detect); otherwise
#     cross-lineage incompatible-pair risk ONLY (granule x T/NK,
#     granule x B/plasma, granule x platelet, granule x erythroid, mono x T).
#     NEVER granule x cell-cycle: that is the genuine cycling-granulocyte-
#     progenitor biology.
#   L2 ambient -- non-granulocyte (exclude 8 neutrophil-axis states + eos/mast)
#     granule score per donor vs donor neutrophil fraction, for BOTH the
#     primary-core and the residual-neutrophil marker panels.
#
#   Exclusion policy:
#     author score  -> top-fraction trimming (continuous risk metric).
#     author label  -> all author-labeled doublets.
#     cross-lineage -> ONLY cells passing the predeclared high-risk flag
#       (risk >= DOUBLE_RISK_QUANTILE); 5%/10% are CAPS, not deletion quotas.
#
#   v1.1: removed always-true placeholder; fixed M1 UMI aggregation; strict/
#     broad doublet-column detection; deterministic tie-break; cohort
#     assertions; require_cols; dual ambient panels; excluded-cell counts.
#   v1.2: EXPECTED_N_DONORS=39; sample_id matching; one-acute-sample assertion;
#     all-NA donors excluded from M1 class medians; audit marks selected
#     column; score-column priority; ambient scores hoisted; interpretation
#     column; signed+abs delta.
#   v1.3: fallback cross-lineage exclusion restricted to high-risk flagged
#     cells (no arbitrary top-risk quota); lineage/ambient marker presence
#     audits; expanded provenance (selected column/type/policy/flagged count);
#     cleaned unmapped-diagnosis message.
#   v1.4: author-score degeneracy guard -- a near-constant numeric doublet
#     column falls back to cross-lineage risk instead of arbitrary top-
#     fraction trimming.
#   v1.5: fix crash when a target state has zero complete paired donors
#     (fine_state/readout assigned with rep(, nrow) so empty pd rows are valid).
#
#   Outputs (tables/):
#     metadata_doublet_columns_audit.csv
#     lineage_marker_presence_audit.csv
#     ambient_marker_presence_audit.csv
#     cross_lineage_doublet_risk_by_cell.csv  (fallback mode only)
#     doublet_sensitivity_m1_source_summary.csv
#     ambient_granule_signal_non_neutrophil_by_donor.csv
#     ambient_neutrophil_fraction_correlation.csv
#     run_provenance.csv
#   RUN: Rscript 05d_technical_robustness_ambient_doublet.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  SCRIPT_NAME    <- "05d_technical_robustness_ambient_doublet.R"
  SCRIPT_VERSION <- "v1.5"
  
  PROJECT_ROOT  <- "/home/sunshine/predicate/singlecell"
  INPUT_FILE    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  OUTPUTS_ROOT  <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  
  STAGE_NAME  <- "05d_technical_robustness_ambient_doublet"
  RUN_PURPOSE <- "doublet_ambient_robustness"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  MIN_CELLS_PER_DONOR_STATE <- 20L
  MIN_NONGRAN_CELLS_PER_DONOR <- 20L
  
  EXPECTED_N_SAMPLES        <- 48L
  EXPECTED_N_DONORS         <- 39L
  EXPECTED_N_ACUTE_DONORS   <- 26L
  EXPECTED_N_HC_DONORS      <- 6L
  EXPECTED_N_SURGERY_DONORS <- 7L
  EXPECTED_N_CONV_SAMPLES   <- 9L
  EXPECTED_N_PAIRED_DONORS  <- 9L
  EXPECTED_N_PRIMARY_CELLS  <- 151837L
  
  EXCL_TOP_FRACS       <- c(0.05, 0.10)
  DOUBLE_RISK_QUANTILE <- 0.75   # both lineage percentile-ranks above this -> flagged
  RHO_STRONG           <- 0.50   # ambient interpretation band (descriptive)
  
  NEUTROPHIL_AXIS <- c(
    "Cycling_neutrophil_progenitors",
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "IL1R2+_immature_neutrophils",
    "S100A8-9_hi_neutrophils",
    "Mature_neutrophils",
    "Degranulating_neutrophils",
    "Apoptosing_neutrophils"
  )
  EOS_MAST <- c("Eosinophils", "Mast_cells/eosiniophils")
  
  # Cross-lineage marker panels (incompatible co-expression => doublet risk).
  LINEAGE_MARKERS <- list(
    granule   = c("MPO","ELANE","CTSG","DEFA4","CEACAM8","TCN1"),
    T_NK      = c("CD3D","CD3E","TRAC","NKG7","GNLY","KLRD1"),
    B_plasma  = c("MS4A1","CD79A","CD79B","MZB1","JCHAIN"),
    platelet  = c("PPBP","PF4","GP9"),
    erythroid = c("HBB","HBA1","HBA2","ALAS2"),
    mono      = c("LYZ","LST1","FCN1","MS4A7","FCGR3A")
  )
  INCOMPATIBLE_PAIRS <- list(
    c("granule","T_NK"),
    c("granule","B_plasma"),
    c("granule","platelet"),
    c("granule","erythroid"),
    c("mono","T_NK")
  )
  MIN_MARKERS_PER_LINEAGE <- 2L
  
  # Ambient marker panels (two layers, tied to Script 04 residual story).
  AMBIENT_MARKER_PANELS <- list(
    primary_core       = c("MPO","ELANE","CTSG","DEFA4"),
    residual_neutrophil = c("CEACAM8","OLFM4","RNASE3","TCN1")
  )
  MIN_MARKERS_PER_AMBIENT_PANEL <- 3L
  
  KEY_M1_STATES <- c(
    "Cycling_neutrophil_progenitors",
    "MPO+_immature_neutrophils_or_progenitors",
    "Cycling_TNK"
  )
  
  # Script 03/03c source classes (do not drift).
  M1_SOURCE_CLASS <- list(
    progenitor_granulopoiesis = c("HSPCs","Cycling_neutrophil_progenitors",
                                  "MPO+_immature_neutrophils_or_progenitors"),
    immature_neutrophil = c("PADI4+_immature_neutrophils","IL1R2+_immature_neutrophils",
                            "S100A8-9_hi_neutrophils"),
    cycling_lymphoid = c("Cycling_TNK")
  )
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  require_cols <- function(x, cols, name) {
    miss <- setdiff(cols, colnames(x))
    if (length(miss)) stop_msg(name, " missing columns: ", paste(miss, collapse = ", "))
    invisible(TRUE)
  }
  safe_median <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else median(x)}
  write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
    write.csv(x, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) { unlink(tmp); stop_msg("Failed to move temp CSV to final path: ", path) }
    invisible(path)
  }
  find_latest_file <- function(root, pattern) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    fs <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (!length(fs)) stop_msg("No file found under ", root, " matching: ", pattern)
    fs[which.max(file.info(fs)$mtime)]
  }
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
  module_short <- function(x) {
    y <- as.character(x)
    y[y %in% c("Cell cycle / proliferation","Cell cycle")] <- "Cell cycle"
    y[y %in% c("Neutrophil degranulation","Neutrophil")] <- "Neutrophil"
    y[y %in% c("Inflammatory / Down","Inflammatory")] <- "Inflammatory"
    y
  }
  is_doublet_label_col <- function(z) {
    zz <- tolower(as.character(z[!is.na(z)]))
    if (!length(zz)) return(FALSE)
    any(zz %in% c("doublet","multiplet","predicted_doublet","true","1","yes")) &&
      any(zz %in% c("singlet","false","0","negative","no"))
  }
  parse_doublet_label <- function(z) {
    if (is.logical(z)) return(ifelse(is.na(z), FALSE, z))
    tolower(as.character(z)) %in% c("doublet","multiplet","predicted_doublet","true","1","yes")
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  needed <- c("SeuratObject","Matrix")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
  
  # ---------------------------------------------------------------------------
  # 2. INPUTS (state order, M1 gene set, gene mapping)
  # ---------------------------------------------------------------------------
  STATE_FILE  <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  FROZEN_FILE <- find_latest_file(SCRIPT01_ROOT, "^frozen_signature_input\\.csv$")
  MAP_FILE    <- find_latest_file(SCRIPT01_ROOT, "^gene_mapping_audit\\.csv$")
  
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  frozen  <- read.csv(FROZEN_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  map01   <- read.csv(MAP_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  require_cols(state01, c("fine_state","display_order"), "fine_state_display_order_audit.csv")
  require_cols(frozen,  c("gene","module"), "frozen_signature_input.csv")
  require_cols(map01,   c("signature_gene","dataset_feature"), "gene_mapping_audit.csv")
  
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L, !anyDuplicated(fine_states))
  if (any(!c(NEUTROPHIL_AXIS, EOS_MAST, KEY_M1_STATES) %in% fine_states)) {
    stop_msg("Some configured states missing from fine_states.")
  }
  
  m1_genes <- frozen$gene[module_short(frozen$module) == "Cell cycle"]
  stopifnot(length(m1_genes) == 11L)
  m1_feat <- unique(unname(map01$dataset_feature[match(m1_genes, map01$signature_gene)]))
  if (anyNA(m1_feat) || length(m1_feat) != 11L) {
    stop_msg("M1 gene mapping incomplete.")
  }
  
  named_states <- unique(unlist(M1_SOURCE_CLASS, use.names = FALSE))
  M1_SOURCE_CLASS[["other"]] <- setdiff(fine_states, named_states)
  if (length(M1_SOURCE_CLASS[["other"]]) != 17L) {
    stop_msg("Unexpected 'other' class size: ", length(M1_SOURCE_CLASS[["other"]]))
  }
  
  # ---------------------------------------------------------------------------
  # 3. READ RDS + ACUTE PRIMARY COHORT (Script 02 logic + fail-closed)
  # ---------------------------------------------------------------------------
  if (!file.exists(INPUT_FILE)) stop_msg("Input RDS not found: ", INPUT_FILE)
  obj <- readRDS(INPUT_FILE)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay=COUNT_ASSAY, layer=COUNT_LAYER, fast=FALSE)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    if ("layer" %in% names(formals(GetAssayData))) cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, layer="counts")
    else cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, slot="counts")
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  if (!identical(colnames(cnt), colnames(obj))) stop_msg("Count colnames != object cells")
  m1_idx <- match(m1_feat, rownames(cnt))
  if (anyNA(m1_idx) || anyDuplicated(m1_idx)) stop_msg("M1 features not cleanly present in count matrix.")
  
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  
  sample_tab <- unique(data.frame(sample_id=cell_sample, diagnosis=cell_diag, stringsAsFactors=FALSE))
  sample_diag_n <- tapply(sample_tab$diagnosis, sample_tab$sample_id, function(x) length(unique(x)))
  if (any(sample_diag_n != 1L)) {
    stop_msg("Some samples have multiple diagnosis values: ",
             paste(names(sample_diag_n)[sample_diag_n != 1L], collapse = ", "))
  }
  sample_tab$donor <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort <- NA_character_
  sample_tab$cohort[sample_tab$diagnosis %in% ACUTE_DIAG] <- "Acute_sepsis"
  sample_tab$cohort[sample_tab$diagnosis == "Conv"] <- "Convalescent"
  sample_tab$cohort[sample_tab$diagnosis == "HV"] <- "Healthy_control"
  sample_tab$cohort[sample_tab$diagnosis == "CS"] <- "Surgery_control"
  if (anyNA(sample_tab$cohort)) {
    stop_msg("Unmapped diagnosis: ",
             paste(unique(sample_tab$diagnosis[is.na(sample_tab$cohort)]), collapse = ", "))
  }
  is_conv_suffix <- grepl("_CONV$", sample_tab$sample_id)
  is_conv_diag   <- sample_tab$diagnosis == "Conv"
  stopifnot(identical(is_conv_suffix, is_conv_diag))
  
  n_samples <- length(unique(sample_tab$sample_id))
  n_donors  <- length(unique(sample_tab$donor))
  acute_donors   <- unique(sample_tab$donor[sample_tab$cohort == "Acute_sepsis"])
  hc_donors      <- unique(sample_tab$donor[sample_tab$cohort == "Healthy_control"])
  surgery_donors <- unique(sample_tab$donor[sample_tab$cohort == "Surgery_control"])
  conv_donors    <- unique(sample_tab$donor[sample_tab$cohort == "Convalescent"])
  paired_donors  <- intersect(acute_donors, conv_donors)
  stopifnot(n_samples == EXPECTED_N_SAMPLES, n_donors == EXPECTED_N_DONORS)
  stopifnot(length(acute_donors)   == EXPECTED_N_ACUTE_DONORS)
  stopifnot(length(hc_donors)      == EXPECTED_N_HC_DONORS)
  stopifnot(length(surgery_donors) == EXPECTED_N_SURGERY_DONORS)
  stopifnot(sum(sample_tab$cohort == "Convalescent") == EXPECTED_N_CONV_SAMPLES)
  stopifnot(length(paired_donors) == EXPECTED_N_PAIRED_DONORS)
  stopifnot(all(conv_donors %in% acute_donors))
  acute_sample_n <- tapply(sample_tab$cohort == "Acute_sepsis", sample_tab$donor, sum)
  stopifnot(all(acute_sample_n[acute_donors] == 1L))
  
  # Cohort assignment: match by SAMPLE, never by donor alone (paired donors
  # have both an acute and a _CONV sample under the same donor id).
  meta2 <- meta
  m <- match(cell_sample, sample_tab$sample_id)
  if (anyNA(m)) stop_msg("Sample resolution failed for some cells.")
  meta2$sample_id_for_script05d <- sample_tab$sample_id[m]
  meta2$donor <- sample_tab$donor[m]
  meta2$cohort <- sample_tab$cohort[m]
  meta2$is_acute <- meta2$cohort == "Acute_sepsis"
  stopifnot(identical(grepl("_CONV$", cell_sample), meta2$cohort == "Convalescent"))
  stopifnot(all(!meta2$is_acute[grepl("_CONV$", cell_sample)]))
  primary_idx <- which(meta2$is_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  
  cnt_p  <- cnt[, primary_idx, drop=FALSE]
  donor_p <- factor(as.character(meta2$donor[primary_idx]), levels=acute_donors)
  state_p <- factor(as.character(meta2[[FINE_STATE_FIELD]][primary_idx]), levels=fine_states)
  stopifnot(!anyNA(donor_p), !anyNA(state_p))
  lib_cell <- as.numeric(Matrix::colSums(cnt_p))
  if (any(!is.finite(lib_cell) | lib_cell <= 0)) stop_msg("Non-positive cell library sizes.")
  n_primary <- length(primary_idx)
  
  # ---------------------------------------------------------------------------
  # 4. METADATA DOUBLET AUDIT (strict vs broad)
  # ---------------------------------------------------------------------------
  strict_doublet_cols <- grep("doublet|multiplet|scrub|scrublet|doubletfinder",
                              colnames(meta2), ignore.case = TRUE, value = TRUE)
  broad_qc_cols <- grep("doublet|multiplet|scrub|scrublet|doubletfinder|prediction|predicted",
                        colnames(meta2), ignore.case = TRUE, value = TRUE)
  audit <- data.frame(
    metadata_column = broad_qc_cols,
    strict_doublet_column = broad_qc_cols %in% strict_doublet_cols,
    class = vapply(broad_qc_cols, function(c) class(meta2[[c]])[1], character(1)),
    n_unique = vapply(broad_qc_cols, function(c) length(unique(meta2[[c]][!is.na(meta2[[c]])])), integer(1)),
    stringsAsFactors = FALSE)
  if (!nrow(audit)) {
    audit <- data.frame(metadata_column = "none_found", strict_doublet_column = FALSE,
                        class = NA_character_, n_unique = NA_integer_, stringsAsFactors = FALSE)
  }
  
  # Auto-select author doublet ONLY from strict columns.
  author_doublet <- NULL
  if (length(strict_doublet_cols)) {
    num_cols <- strict_doublet_cols[vapply(strict_doublet_cols,
                                           function(c) is.numeric(meta2[[c]]), logical(1))]
    label_cols <- strict_doublet_cols[vapply(strict_doublet_cols, function(c) {
      is.logical(meta2[[c]]) || is_doublet_label_col(meta2[[c]])
    }, logical(1))]
    if (length(num_cols)) {
      pri <- grep("score|scrublet|doubletfinder|prob|pann", num_cols,
                  ignore.case = TRUE, value = TRUE)
      author_doublet <- list(col = if (length(pri)) pri[1] else num_cols[1], type = "score")
    } else if (length(label_cols)) {
      author_doublet <- list(col = label_cols[1], type = "label")
    }
  }
  
  # Guard against a degenerate/near-constant numeric column being treated as a
  # usable author doublet score (would make top-fraction trimming arbitrary).
  if (!is.null(author_doublet) && identical(author_doublet$type, "score")) {
    v0 <- suppressWarnings(as.numeric(meta2[[author_doublet$col]][primary_idx]))
    v0f <- v0[is.finite(v0)]
    score_ok <- length(v0f) >= 10L &&
      length(unique(v0f)) >= 5L &&
      diff(range(v0f)) > 0
    if (!score_ok) {
      warning("Selected author doublet score column '", author_doublet$col,
              "' is degenerate or near-constant in primary cells; ",
              "falling back to cross-lineage incompatible-pair risk.")
      author_doublet <- NULL
    }
  }
  
  # Mark the actually-selected column in the audit.
  if (nrow(audit) && !identical(audit$metadata_column, "none_found")) {
    audit$selected_for_author_doublet <- FALSE
    if (!is.null(author_doublet)) {
      audit$selected_for_author_doublet[audit$metadata_column == author_doublet$col] <- TRUE
    }
  } else {
    audit$selected_for_author_doublet <- FALSE
  }
  write_csv_atomic(audit, file.path(TABLE_DIR, "metadata_doublet_columns_audit.csv"))
  
  # ---------------------------------------------------------------------------
  # 5. PER-CELL LINEAGE SCORES + MARKER PRESENCE AUDIT
  # ---------------------------------------------------------------------------
  lin_genes <- unique(unlist(LINEAGE_MARKERS, use.names = FALSE))
  lin_present <- intersect(lin_genes, rownames(cnt))
  lin_missing <- setdiff(lin_genes, lin_present)
  if (length(lin_missing)) warning("Lineage marker genes absent (dropped): ",
                                   paste(lin_missing, collapse = ", "))
  
  lineage_marker_audit <- do.call(rbind, lapply(names(LINEAGE_MARKERS), function(lg) {
    gs <- LINEAGE_MARKERS[[lg]]
    present <- intersect(gs, rownames(cnt))
    missing <- setdiff(gs, present)
    data.frame(
      lineage = lg,
      n_total = length(gs),
      n_present = length(present),
      n_missing = length(missing),
      usable = length(present) >= MIN_MARKERS_PER_LINEAGE,
      present_genes = paste(present, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE)
  }))
  write_csv_atomic(lineage_marker_audit, file.path(TABLE_DIR, "lineage_marker_presence_audit.csv"))
  
  X_lin <- as.matrix(cnt_p[match(lin_present, rownames(cnt)), , drop = FALSE])
  rownames(X_lin) <- lin_present
  logcpm_lin <- log2(sweep(X_lin, 2L, lib_cell, "/") * 1e6 + 1)
  
  lineage_score <- lapply(LINEAGE_MARKERS, function(gs) {
    g <- intersect(gs, lin_present)
    if (!length(g)) return(rep(NA_real_, n_primary))
    colMeans(logcpm_lin[match(g, rownames(logcpm_lin)), , drop = FALSE])
  })
  n_markers_lineage <- vapply(LINEAGE_MARKERS, function(gs) length(intersect(gs, lin_present)), integer(1))
  usable_lineages <- names(n_markers_lineage)[n_markers_lineage >= MIN_MARKERS_PER_LINEAGE]
  dropped_lineages <- setdiff(names(LINEAGE_MARKERS), usable_lineages)
  if (length(dropped_lineages)) warning("Lineages with <", MIN_MARKERS_PER_LINEAGE,
                                        " present markers excluded from pairs: ",
                                        paste(dropped_lineages, collapse = ", "))
  
  # ---------------------------------------------------------------------------
  # 6. DOUBLET-RISK DEFINITION + EXCLUSION SETS (deterministic, flag-gated)
  # ---------------------------------------------------------------------------
  doublet_source <- "none"
  risk <- rep(0, n_primary)
  
  if (!is.null(author_doublet) && author_doublet$type == "score") {
    v <- as.numeric(meta2[[author_doublet$col]][primary_idx])
    na_v <- is.na(v)
    if (any(na_v)) {
      warning("Author doublet score has ", sum(na_v), " NA primary cells; treated as lowest risk.")
      v[na_v] <- if (any(!na_v)) min(v[!na_v], na.rm = TRUE) else 0
    }
    risk <- v
    doublet_source <- paste0("author_score:", author_doublet$col)
    
  } else if (!is.null(author_doublet) && author_doublet$type == "label") {
    flagged <- parse_doublet_label(meta2[[author_doublet$col]][primary_idx])
    risk <- as.numeric(flagged)
    doublet_source <- paste0("author_label:", author_doublet$col)
    
  } else {
    doublet_source <- "cross_lineage_incompatible_pairs"
    z <- lapply(usable_lineages, function(lg) {
      s <- lineage_score[[lg]]
      rank(s, ties.method = "average") / length(s)
    })
    names(z) <- usable_lineages
    risk <- rep(0, n_primary)
    used_pairs <- character(0)
    for (pr in INCOMPATIBLE_PAIRS) {
      if (all(pr %in% usable_lineages)) {
        risk <- pmax(risk, pmin(z[[pr[1]]], z[[pr[2]]]))
        used_pairs <- c(used_pairs, paste(pr, collapse = " x "))
      }
    }
    if (!length(used_pairs)) stop_msg("No usable incompatible lineage pairs after marker filtering.")
    cat("Cross-lineage pairs used:", paste(used_pairs, collapse = "; "), "\n")
    
    risk_by_cell <- data.frame(
      cell_barcode = colnames(cnt_p),
      fine_state = as.character(state_p),
      donor = as.character(donor_p),
      stringsAsFactors = FALSE)
    for (lg in usable_lineages) risk_by_cell[[paste0("score_", lg)]] <- lineage_score[[lg]]
    risk_by_cell$doublet_risk <- risk
    risk_by_cell$flagged_doublet <- risk >= DOUBLE_RISK_QUANTILE
    write_csv_atomic(risk_by_cell, file.path(TABLE_DIR, "cross_lineage_doublet_risk_by_cell.csv"))
    cat("Cells flagged doublet-risk (both lineages > p", DOUBLE_RISK_QUANTILE * 100,
        "): ", sum(risk >= DOUBLE_RISK_QUANTILE), "\n", sep = "")
  }
  
  # Exclusion policy.
  #   author score  -> top-fraction trimming (continuous risk metric).
  #   author label  -> all author-labeled doublets.
  #   cross-lineage -> ONLY flag-passing cells; 5%/10% are CAPS, not quotas.
  n_crosslineage_flagged <- if (identical(doublet_source, "cross_lineage_incompatible_pairs")) {
    sum(risk >= DOUBLE_RISK_QUANTILE)
  } else {
    NA_integer_
  }
  excl_by_top_fraction <- function(f, require_crosslineage_flag = FALSE) {
    n_excl <- ceiling(n_primary * f)
    if (n_excl <= 0L) return(integer(0))
    ord <- order(-risk, colnames(cnt_p))
    if (require_crosslineage_flag) {
      ord <- ord[risk[ord] >= DOUBLE_RISK_QUANTILE]
    }
    if (!length(ord)) return(integer(0))
    ord[seq_len(min(n_excl, length(ord)))]
  }
  
  if (grepl("^author_score", doublet_source)) {
    excl_list <- lapply(EXCL_TOP_FRACS, function(f) excl_by_top_fraction(f, FALSE))
    names(excl_list) <- paste0("excl_author_score_top", EXCL_TOP_FRACS * 100, "pct")
  } else if (grepl("^author_label", doublet_source)) {
    flagged_cells <- which(risk == 1)
    excl_list <- list(excl_author_doublets = flagged_cells)
    cat("Author-labeled doublets in primary cells:", length(flagged_cells), "\n")
  } else {
    excl_list <- lapply(EXCL_TOP_FRACS, function(f) excl_by_top_fraction(f, TRUE))
    names(excl_list) <- paste0("excl_crosslineage_flagged_top", EXCL_TOP_FRACS * 100, "pct")
    cat("Cross-lineage flagged cells available for exclusion:", n_crosslineage_flagged, "\n")
  }
  
  # ---------------------------------------------------------------------------
  # 7. M1 SOURCE RECOMPUTATION (full vs excluded)
  # ---------------------------------------------------------------------------
  compute_m1_summary <- function(keep, label) {
    kp <- primary_idx[keep]
    sub_cnt <- cnt[, kp, drop = FALSE]
    sub_donor <- factor(as.character(meta2$donor[kp]), levels = acute_donors)
    sub_state <- factor(as.character(meta2[[FINE_STATE_FIELD]][kp]), levels = fine_states)
    stopifnot(!anyNA(sub_donor), !anyNA(sub_state))
    n_d <- nlevels(sub_donor); n_s <- length(fine_states)
    ds <- (as.integer(sub_donor) - 1L) * n_s + as.integer(sub_state)
    Ds <- Matrix::sparseMatrix(i = ds, j = seq_along(ds), x = 1, dims = c(n_d * n_s, length(ds)))
    nc_flat <- as.integer(Matrix::rowSums(Ds))
    ncell_mat <- matrix(nc_flat, nrow = n_d, ncol = n_s, byrow = TRUE,
                        dimnames = list(acute_donors, fine_states))
    elig <- ncell_mat >= MIN_CELLS_PER_DONOR_STATE
    lib_c <- as.numeric(Matrix::colSums(sub_cnt))
    if (any(!is.finite(lib_c) | lib_c <= 0)) stop_msg("Non-positive library size in subset.")
    lib_flat <- as.numeric(Ds %*% lib_c)
    # per-cell M1 total UMI, then donor-state aggregation
    m1_cell_umi <- as.numeric(Matrix::colSums(sub_cnt[m1_idx, , drop = FALSE]))
    umi_flat <- as.numeric(Ds %*% m1_cell_umi)
    umi_mat <- matrix(umi_flat, nrow = n_d, ncol = n_s, byrow = TRUE,
                      dimnames = list(acute_donors, fine_states))
    tot_d <- rowSums(umi_mat)
    cf <- sweep(umi_mat, 1, tot_d, "/")
    cf[tot_d == 0, ] <- NA_real_
    cls_med <- vapply(names(M1_SOURCE_CLASS), function(cl) {
      sts <- intersect(M1_SOURCE_CLASS[[cl]], fine_states)
      mat <- cf[, sts, drop = FALSE]
      rs <- rowSums(ifelse(is.finite(mat), mat, 0))
      rs[!is.finite(tot_d) | tot_d == 0] <- NA_real_   # no-M1-UMI donor is NA, not 0
      safe_median(rs)
    }, numeric(1))
    top_cls <- names(cls_med)[which.max(cls_med)]
    comp <- log2(umi_flat / lib_flat * 1e6 + 1)
    comp[!is.finite(comp) | lib_flat <= 0] <- NA_real_
    comp_mat <- matrix(comp, nrow = n_d, ncol = n_s, byrow = TRUE,
                       dimnames = list(acute_donors, fine_states))
    comp_mat[!elig] <- NA_real_
    key_med <- vapply(KEY_M1_STATES, function(s) safe_median(comp_mat[, s]), numeric(1))
    data.frame(
      set = label,
      n_cells = length(keep),
      n_cells_excluded = n_primary - length(keep),
      excluded_fraction = (n_primary - length(keep)) / n_primary,
      med_progenitor = cls_med[["progenitor_granulopoiesis"]],
      med_immature_neutrophil = cls_med[["immature_neutrophil"]],
      med_cycling_lymphoid = cls_med[["cycling_lymphoid"]],
      med_other = cls_med[["other"]],
      top_source_class = top_cls,
      progenitor_minus_other = cls_med[["progenitor_granulopoiesis"]] - cls_med[["other"]],
      M1_comp_Cycling_neut_prog = key_med[["Cycling_neutrophil_progenitors"]],
      M1_comp_MPOimm = key_med[["MPO+_immature_neutrophils_or_progenitors"]],
      M1_comp_Cycling_TNK = key_med[["Cycling_TNK"]],
      stringsAsFactors = FALSE)
  }
  
  full_set <- seq_len(n_primary)
  m1_sets <- list(full = full_set)
  for (nm in names(excl_list)) m1_sets[[nm]] <- setdiff(full_set, excl_list[[nm]])
  if (all(vapply(excl_list, length, integer(1)) == 0L)) {
    cat("All exclusion sets empty; only full-set M1 summary produced.\n")
    m1_sets <- list(full = full_set)
  }
  
  m1_sens <- do.call(rbind, lapply(names(m1_sets), function(nm)
    compute_m1_summary(m1_sets[[nm]], nm)))
  write_csv_atomic(m1_sens, file.path(TABLE_DIR, "doublet_sensitivity_m1_source_summary.csv"))
  
  # ---------------------------------------------------------------------------
  # 8. AMBIENT RNA CHECK (descriptive proxy; two marker panels)
  # ---------------------------------------------------------------------------
  ambient_marker_audit <- do.call(rbind, lapply(names(AMBIENT_MARKER_PANELS), function(pn) {
    gs <- AMBIENT_MARKER_PANELS[[pn]]
    present <- intersect(gs, rownames(cnt))
    missing <- setdiff(gs, present)
    data.frame(
      panel = pn,
      n_total = length(gs),
      n_present = length(present),
      n_missing = length(missing),
      usable = length(present) >= MIN_MARKERS_PER_AMBIENT_PANEL,
      present_genes = paste(present, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE)
  }))
  write_csv_atomic(ambient_marker_audit, file.path(TABLE_DIR, "ambient_marker_presence_audit.csv"))
  
  non_gran_states <- setdiff(fine_states, c(NEUTROPHIL_AXIS, EOS_MAST))
  non_gran_idx <- which(as.character(state_p) %in% non_gran_states)
  donor_non_gran <- split(non_gran_idx, as.character(donor_p)[non_gran_idx])
  donor_nf <- as.numeric(tapply(as.character(state_p) %in% NEUTROPHIL_AXIS, donor_p, mean))
  names(donor_nf) <- acute_donors
  
  panel_score <- function(gs) {
    g <- intersect(gs, rownames(cnt))
    if (length(g) < MIN_MARKERS_PER_AMBIENT_PANEL) {
      stop_msg("Ambient panel '", paste(gs, collapse = ","), "' has too few present markers (", length(g), ").")
    }
    X <- as.matrix(cnt_p[match(g, rownames(cnt)), , drop = FALSE])
    colMeans(log2(sweep(X, 2L, lib_cell, "/") * 1e6 + 1))
  }
  ambient_scores <- lapply(AMBIENT_MARKER_PANELS, panel_score)
  
  ambient_rows <- lapply(acute_donors, function(d) {
    idx <- donor_non_gran[[d]]
    n_non_gran <- if (is.null(idx)) 0L else length(idx)
    row <- data.frame(donor = d, neutrophil_fraction = donor_nf[[d]],
                      n_non_granulocyte_cells = n_non_gran, stringsAsFactors = FALSE)
    for (pn in names(AMBIENT_MARKER_PANELS)) {
      row[[paste0("median_", pn)]] <- if (n_non_gran >= MIN_NONGRAN_CELLS_PER_DONOR) {
        safe_median(ambient_scores[[pn]][idx])
      } else {
        NA_real_
      }
    }
    row
  })
  ambient_df <- do.call(rbind, ambient_rows)
  write_csv_atomic(ambient_df, file.path(TABLE_DIR, "ambient_granule_signal_non_neutrophil_by_donor.csv"))
  
  ambient_corr <- lapply(names(AMBIENT_MARKER_PANELS), function(pn) {
    col <- paste0("median_", pn)
    ok <- is.finite(ambient_df[[col]])
    sd_v <- ambient_df[[col]][ok]
    nf_v <- ambient_df$neutrophil_fraction[ok]
    n_ok <- sum(ok)
    score_median <- if (n_ok) median(sd_v) else NA_real_
    score_range  <- if (n_ok) diff(range(sd_v)) else NA_real_
    reason <- if (n_ok < 5L) {
      "not_evaluable_too_few_donors"
    } else if (all(abs(sd_v) < 1e-8)) {
      "not_evaluable_zero_signal"
    } else if (score_range < 1e-8) {
      "not_evaluable_constant_nonzero"
    } else {
      "evaluable"
    }
    rho <- if (reason == "evaluable") {
      suppressWarnings(cor(sd_v, nf_v, method = "spearman"))
    } else {
      NA_real_
    }
    interpretation <- if (reason != "evaluable") {
      reason
    } else if (rho >= RHO_STRONG) {
      "possible_ambient_contribution"
    } else {
      "ambient_not_supported_by_neutrophil_burden_proxy"
    }
    data.frame(panel = pn, n_donors_evaluable = n_ok,
               donor_score_median = score_median, donor_score_range = score_range,
               spearman_rho = rho, evaluation_reason = reason,
               interpretation = interpretation,
               stringsAsFactors = FALSE)
  })
  ambient_corr_df <- do.call(rbind, ambient_corr)
  write_csv_atomic(ambient_corr_df, file.path(TABLE_DIR, "ambient_neutrophil_fraction_correlation.csv"))
  
  # ---------------------------------------------------------------------------
  # 9. PROVENANCE + CONSOLE SUMMARY
  # ---------------------------------------------------------------------------
  selected_doublet_col  <- if (!is.null(author_doublet)) author_doublet$col else "none"
  selected_doublet_type <- if (!is.null(author_doublet)) author_doublet$type else "none"
  exclusion_policy <- if (grepl("^author_score", doublet_source)) {
    "author_score_top_fraction"
  } else if (grepl("^author_label", doublet_source)) {
    "author_label_all_flagged"
  } else {
    "cross_lineage_flagged_only_capped_by_top_fraction"
  }
  
  provenance <- data.frame(
    field = c("script","version","date","input_file","doublet_source",
              "selected_doublet_column","selected_doublet_type",
              "exclusion_policy","n_crosslineage_flagged",
              "exclusion_levels","doublet_risk_quantile",
              "ambient_min_cells_per_donor","ambient_panels",
              "n_acute_donors","n_primary_cells","lineages_dropped","R_version"),
    value = c(SCRIPT_NAME, SCRIPT_VERSION, as.character(Sys.time()), INPUT_FILE,
              doublet_source,
              selected_doublet_col, selected_doublet_type,
              exclusion_policy,
              as.character(n_crosslineage_flagged),
              paste(names(m1_sets), collapse = ";"),
              as.character(DOUBLE_RISK_QUANTILE),
              as.character(MIN_NONGRAN_CELLS_PER_DONOR),
              paste(names(AMBIENT_MARKER_PANELS), collapse = ";"),
              as.character(length(acute_donors)), as.character(n_primary),
              if (length(dropped_lineages)) paste(dropped_lineages, collapse = ";") else "none",
              as.character(getRversion())),
    stringsAsFactors = FALSE)
  write_csv_atomic(provenance, file.path(TABLE_DIR, "run_provenance.csv"))
  
  cat("Script 05d v1.4 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("Doublet source:", doublet_source, "| policy:", exclusion_policy, "\n")
  cat("\n--- M1 source sensitivity (full vs doublet-excluded) ---\n")
  print(m1_sens, row.names = FALSE)
  cat("\n--- M1 doublet interpretation ---\n")
  base <- m1_sens[m1_sens$set == "full", , drop = FALSE]
  if (nrow(m1_sens) > 1L) {
    for (nm in setdiff(m1_sens$set, "full")) {
      z <- m1_sens[m1_sens$set == nm, , drop = FALSE]
      top_same <- identical(z$top_source_class, base$top_source_class)
      prog_delta <- z$med_progenitor - base$med_progenitor
      cat(nm, ": top source", if (top_same) "unchanged" else paste0("FLIPPED to ", z$top_source_class),
          "| progenitor signed delta", sprintf("%.3f", prog_delta),
          "| abs delta", sprintf("%.3f", abs(prog_delta)), "\n")
    }
  } else {
    cat("No exclusion performed; M1 unchanged by construction.\n")
  }
  cat("\n--- Ambient (descriptive proxy, not decontamination) ---\n")
  for (i in seq_len(nrow(ambient_corr_df))) {
    cat(ambient_corr_df$panel[i], "| donors:", ambient_corr_df$n_donors_evaluable[i],
        "| spearman vs neutrophil fraction:",
        if (is.finite(ambient_corr_df$spearman_rho[i])) sprintf("%.3f", ambient_corr_df$spearman_rho[i]) else "NA",
        "|", ambient_corr_df$interpretation[i], "\n")
  }
  cat("\nNote: M1 stability under doublet exclusion means 'not obviously doublet-driven', not proof of biology.\n")
  cat("Ambient check is a descriptive proxy only; it does not exclude contamination.\n")
})



# ============================================================================
# SCRIPT 05e v1.5 -- internal pseudobulk sanity + state-intrinsic direction
#   Purpose: two-layer internal consistency check for the 32-gene signature
#   direction and for whether the M2/primary-core programs are intrinsically
#   higher in acute sepsis within the immature-neutrophil carrying states.
#
#   Positioning:
#     P3A = internal WHOLE-SAMPLE pseudobulk sanity check (composition-
#           sensitive, NOT a validation; only tests that scRNA does not
#           invert the bulk signature direction).
#     P3B = state-intrinsic module direction (THE main result): acute vs
#           convalescent, paired by donor, within immature neutrophil states.
#           HV/CS are descriptive only.
#   NOT a perturbation, NOT a regulation test, NOT matched-bulk validation.
#   total32_signed is an internal non-inversion sanity score, NOT a validated
#   quantitative MARS score.
#
#   Module scores use SUMMED-UMI composite (sum UMI / library -> log2+1),
#   consistent with Scripts 03/03b. Readout genes are resolved to dataset
#   features via the Script-01 mapping. Each readout has a per-set minimum
#   number of present genes (READOUT_MIN_PRESENT); below it the script fails
#   closed. The full 32/32 frozen signature must be present in counts.
#
#   P3B primary rows are MPO+/PADI4+ immature x {M2_full, primary_core}
#   (result_layer = primary_M2_state_intrinsic). Other state/readout rows are
#   background and are labelled secondary/supplementary. P3B direction calls
#   are DESCRIPTIVE (acute_higher_descriptive / acute_lower_descriptive);
#   they are not strong reverse-direction evidence.
#
#   P3B interpretation:
#     - if acute M2/primary_core > recovery within MPO+/PADI4+ immature,
#       sepsis up-regulates the early granulocyte program per cell;
#     - if state expression is unchanged but cell_fraction is higher in
#       acute, the bulk signal is largely composition-driven;
#     - if the state collapses in convalescence (<20 cells), not_evaluable,
#       no forced direction. not_evaluable_reason distinguishes collapse of
#       the acute side vs convalescent side vs insufficient paired donors.
#     Direction calls require BOTH a median effect AND donor support rate
#     (median >= 0.3 log2 AND >= 2/3 paired donors agreeing); otherwise
#     no_clear_direction. P3A uses symmetric descriptive thresholds,
#     with a separate +/-0.5 sensitivity classification.
#
#   v1.1: fixed sample_lib ordering; summed-UMI composite; dataset-feature
#     resolution; readout membership audit; P3A support rate; P3B HV/CS counts.
#   v1.2: support-rate-aware P3B paired direction calls; stricter P3A
#     opposite-direction flag; sample-state dimnames; provenance opposite
#     thresholds; total32_signed sanity-score note.
#   v1.3: readout-specific membership guard (READOUT_MIN_PRESENT, fail-closed);
#     P3B paired evaluability diagnostics + not_evaluable_reason; target-state
#     n-cell medians for sparse-state interpretation.
#   v1.4: P3B rows labelled primary/secondary/supplementary; fail-closed
#     32/32 signature-feature guard; module_composite dimension guard; P3B
#     direction calls clarified as descriptive acute_higher/acute_lower flags.
#   v1.5: symmetric P3A direction classification at +/-0.3, with a common
#     minimum of 4 paired donors and positive-pair fractions >=2/3 or <=1/3;
#     added +/-0.5 classification sensitivity. P3B calculations unchanged.
#
#   Outputs (tables/):
#     readout_gene_membership.csv
#     p3a_sample_module_scores.csv
#     p3a_gene_direction_consistency.csv
#     p3a_module_summary.csv
#     p3b_paired_delta_long.csv
#     p3b_state_readout_condition_summary.csv
#     p3b_direction_calls.csv
#     run_provenance.csv
#   RUN: Rscript 05e_internal_pseudobulk_sanity_and_state_direction.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  SCRIPT_NAME    <- "05e_internal_pseudobulk_sanity_and_state_direction.R"
  SCRIPT_VERSION <- "v1.4"
  
  PROJECT_ROOT  <- "/home/sunshine/predicate/singlecell"
  INPUT_FILE    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  OUTPUTS_ROOT  <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  
  STAGE_NAME  <- "05e_pseudobulk_sanity_state_direction"
  RUN_PURPOSE <- "internal_pseudobulk_sanity_and_state_direction"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  MIN_CELLS_PER_SAMPLE_STATE <- 20L
  
  EXPECTED_N_SAMPLES        <- 48L
  EXPECTED_N_DONORS         <- 39L
  EXPECTED_N_ACUTE_DONORS   <- 26L
  EXPECTED_N_HC_DONORS      <- 6L
  EXPECTED_N_SURGERY_DONORS <- 7L
  EXPECTED_N_CONV_SAMPLES   <- 9L
  EXPECTED_N_PAIRED_DONORS  <- 9L
  
  # Descriptive call thresholds (flags, not significance).
  # These two existing settings also serve P3B and remain unchanged.
  PAIRED_MIN_N    <- 4L
  DIRECTION_DELTA <- 0.3
  
  # P3A: symmetric primary classification and threshold sensitivity.
  P3A_ALIGNED_FRAC_MIN <- 2 / 3
  P3A_OPPOSITE_FRAC_MAX <- 1 / 3
  P3A_SENSITIVITY_DELTA <- 0.5
  
  # Script 03b readout modules (do not drift).
  READOUT_SETS <- list(
    M2_full            = c("CEACAM6","CEACAM8","CTSG","DEFA4","ELANE","MPO",
                           "MS4A3","OLFM4","RNASE3","TCN1"),
    primary_core       = c("MPO","ELANE","CTSG","DEFA4"),
    secondary_granule  = c("LTF","LCN2"),
    tertiary_granule   = c("MMP8","MMP9"),
    mature_surface     = c("FCGR3B","CXCR2","CXCR1","FPR1","MME"),
    inflammatory_S100  = c("S100A8","S100A9","S100A12"),
    M1_cell_cycle      = NULL    # filled from frozen Script-01 input
  )
  # Per-readout minimum number of present genes (fail-closed guard).
  READOUT_MIN_PRESENT <- c(
    M2_full = 8L,
    primary_core = 3L,
    secondary_granule = 2L,
    tertiary_granule = 2L,
    mature_surface = 3L,
    inflammatory_S100 = 2L,
    M1_cell_cycle = 8L
  )
  
  # P3B target states (immature neutrophil axis; the M2 carrying states).
  TARGET_STATES <- c(
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "Cycling_neutrophil_progenitors",
    "S100A8-9_hi_neutrophils"
  )
  # P3B primary rows (the main conclusion). Everything else is background.
  P3B_PRIMARY_READOUTS <- c("M2_full", "primary_core")
  P3B_PRIMARY_STATES <- c(
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils"
  )
  p3b_layer <- function(state, readout) {
    if (readout %in% P3B_PRIMARY_READOUTS && state %in% P3B_PRIMARY_STATES) {
      "primary_M2_state_intrinsic"
    } else if (readout %in% P3B_PRIMARY_READOUTS) {
      "secondary_M2_state_context"
    } else {
      "supplementary_readout_context"
    }
  }
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  require_cols <- function(x, cols, name) {
    miss <- setdiff(cols, colnames(x))
    if (length(miss)) stop_msg(name, " missing columns: ", paste(miss, collapse = ", "))
    invisible(TRUE)
  }
  safe_median <- function(x){x<-x[is.finite(x)]; if(!length(x)) NA_real_ else median(x)}
  classify_p3a_direction <- function(n_pairs, med_delta, positive_fraction,
                                     delta_cut) {
    if (!is.finite(n_pairs) ||
        n_pairs < PAIRED_MIN_N ||
        !is.finite(med_delta) ||
        !is.finite(positive_fraction)) {
      return("not_evaluable")
    }
    
    if (med_delta >= delta_cut &&
        positive_fraction >= P3A_ALIGNED_FRAC_MIN) {
      return("aligned")
    }
    
    if (med_delta <= -delta_cut &&
        positive_fraction <= P3A_OPPOSITE_FRAC_MAX) {
      return("opposite_direction")
    }
    
    "neither_criterion_met"
  }
   write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
    write.csv(x, tmp, row.names = FALSE, fileEncoding = "UTF-8")
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) { unlink(tmp); stop_msg("Failed to move temp CSV to final path: ", path) }
    invisible(path)
  }
  find_latest_file <- function(root, pattern) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    fs <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (!length(fs)) stop_msg("No file found under ", root, " matching: ", pattern)
    fs[which.max(file.info(fs)$mtime)]
  }
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
  module_short <- function(x) {
    y <- as.character(x)
    y[y %in% c("Cell cycle / proliferation","Cell cycle")] <- "Cell cycle"
    y[y %in% c("Neutrophil degranulation","Neutrophil")] <- "Neutrophil"
    y[y %in% c("Inflammatory / Down","Inflammatory")] <- "Inflammatory"
    y
  }
  # Summed-UMI composite (sum member-gene UMI / library -> log2+1),
  # consistent with Script 03/03b module composite. Dimension-guarded.
  module_composite <- function(count_mat, lib_vec, genes, min_genes) {
    if (length(lib_vec) != ncol(count_mat)) {
      stop_msg("module_composite: lib_vec length (", length(lib_vec),
               ") != ncol(count_mat) (", ncol(count_mat), ").")
    }
    g <- intersect(genes, rownames(count_mat))
    if (length(g) < min_genes) return(rep(NA_real_, ncol(count_mat)))
    x <- Matrix::colSums(count_mat[g, , drop = FALSE])
    y <- log2(x / lib_vec * 1e6 + 1)
    y[!is.finite(y) | lib_vec <= 0] <- NA_real_
    y
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  needed <- c("SeuratObject","Matrix")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
  
  # ---------------------------------------------------------------------------
  # 2. INPUTS + GENE DEFINITIONS (resolve to dataset features)
  # ---------------------------------------------------------------------------
  STATE_FILE  <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  FROZEN_FILE <- find_latest_file(SCRIPT01_ROOT, "^frozen_signature_input\\.csv$")
  MAP_FILE    <- find_latest_file(SCRIPT01_ROOT, "^gene_mapping_audit\\.csv$")
  
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  frozen  <- read.csv(FROZEN_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  map01   <- read.csv(MAP_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  require_cols(state01, c("fine_state","display_order"), "fine_state_display_order_audit.csv")
  require_cols(frozen,  c("gene","module","direction"), "frozen_signature_input.csv")
  require_cols(map01,   c("signature_gene","dataset_feature"), "gene_mapping_audit.csv")
  
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L, !anyDuplicated(fine_states))
  stopifnot(all(TARGET_STATES %in% fine_states))
  
  if (anyNA(frozen$direction) || !all(frozen$direction %in% c("Up","Down"))) {
    stop_msg("Invalid direction values in frozen input.")
  }
  stopifnot(sum(frozen$direction == "Up") == 29L, sum(frozen$direction == "Down") == 3L)
  
  m1_genes <- frozen$gene[module_short(frozen$module) == "Cell cycle"]
  stopifnot(length(m1_genes) == 11L)
  READOUT_SETS[["M1_cell_cycle"]] <- m1_genes
  
  sig_feat <- map01$dataset_feature[match(frozen$gene, map01$signature_gene)]
  if (anyNA(sig_feat)) {
    stop_msg("Some frozen signature genes lack dataset_feature mapping: ",
             paste(frozen$gene[is.na(sig_feat)], collapse = ", "))
  }
  if (anyDuplicated(sig_feat)) {
    stop_msg("Duplicated dataset_feature mappings in frozen signature: ",
             paste(unique(sig_feat[duplicated(sig_feat)]), collapse = ", "))
  }
  gene_to_feature <- setNames(map01$dataset_feature, map01$signature_gene)
  to_feature <- function(gs) {
    out <- gs
    hit <- gs %in% names(gene_to_feature)
    out[hit] <- unname(gene_to_feature[gs[hit]])
    out
  }
  READOUT_FEATURE_SETS <- lapply(READOUT_SETS, to_feature)
  if (!setequal(names(READOUT_FEATURE_SETS), names(READOUT_MIN_PRESENT))) {
    stop_msg("READOUT_MIN_PRESENT names do not match READOUT_FEATURE_SETS.")
  }
  READOUT_MIN_PRESENT <- READOUT_MIN_PRESENT[names(READOUT_FEATURE_SETS)]
  up_feat   <- sig_feat[frozen$direction == "Up"]
  down_feat <- sig_feat[frozen$direction == "Down"]
  
  # ---------------------------------------------------------------------------
  # 3. READ RDS + COUNT LAYER
  # ---------------------------------------------------------------------------
  if (!file.exists(INPUT_FILE)) stop_msg("Input RDS not found: ", INPUT_FILE)
  obj <- readRDS(INPUT_FILE)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay=COUNT_ASSAY, layer=COUNT_LAYER, fast=FALSE)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    if ("layer" %in% names(formals(GetAssayData))) cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, layer="counts")
    else cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, slot="counts")
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  if (!identical(colnames(cnt), colnames(obj))) stop_msg("Count colnames != object cells")
  
  all_feats <- unique(c(sig_feat, unlist(READOUT_FEATURE_SETS, use.names = FALSE)))
  present_feats <- intersect(all_feats, rownames(cnt))
  missing_feats <- setdiff(all_feats, present_feats)
  if (length(missing_feats)) warning("Requested features absent from counts (dropped): ",
                                     paste(missing_feats, collapse = ", "))
  if (!length(present_feats)) stop_msg("No requested features present in count matrix.")
  
  n32_present <- sum(sig_feat %in% rownames(cnt))
  if (n32_present != 32L) {
    stop_msg("Expected 32/32 frozen signature features present in counts, got ",
             n32_present, ".")
  }
  cat("32-gene features present in counts: 32 / 32\n")
  
  readout_membership <- do.call(rbind, lapply(names(READOUT_FEATURE_SETS), function(ro) {
    feats <- READOUT_FEATURE_SETS[[ro]]
    data.frame(
      readout = ro,
      requested_gene = paste(READOUT_SETS[[ro]], collapse = ";"),
      requested_feature = paste(feats, collapse = ";"),
      n_requested = length(feats),
      n_present = sum(feats %in% rownames(cnt)),
      present_features = paste(intersect(feats, rownames(cnt)), collapse = ";"),
      missing_features = paste(setdiff(feats, rownames(cnt)), collapse = ";"),
      stringsAsFactors = FALSE)
  }))
  readout_membership$min_present_required <- unname(READOUT_MIN_PRESENT[readout_membership$readout])
  readout_membership$membership_status <- ifelse(
    readout_membership$n_present >= readout_membership$min_present_required,
    "pass", "fail")
  if (any(readout_membership$membership_status == "fail")) {
    bad <- readout_membership[readout_membership$membership_status == "fail", , drop = FALSE]
    stop_msg("Readout gene membership below required minimum: ",
             paste(paste0(bad$readout, " ", bad$n_present, "/", bad$min_present_required),
                   collapse = "; "))
  }
  write_csv_atomic(readout_membership, file.path(TABLE_DIR, "readout_gene_membership.csv"))
  
  # ---------------------------------------------------------------------------
  # 4. COHORT (all 48 samples; fail-closed, same as Scripts 04b/05d)
  # ---------------------------------------------------------------------------
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  
  sample_tab <- unique(data.frame(sample_id=cell_sample, diagnosis=cell_diag, stringsAsFactors=FALSE))
  sample_diag_n <- tapply(sample_tab$diagnosis, sample_tab$sample_id, function(x) length(unique(x)))
  if (any(sample_diag_n != 1L)) {
    stop_msg("Some samples have multiple diagnosis values: ",
             paste(names(sample_diag_n)[sample_diag_n != 1L], collapse = ", "))
  }
  sample_tab$donor <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort <- NA_character_
  sample_tab$cohort[sample_tab$diagnosis %in% ACUTE_DIAG] <- "Acute"
  sample_tab$cohort[sample_tab$diagnosis == "Conv"] <- "Conv"
  sample_tab$cohort[sample_tab$diagnosis == "HV"] <- "HV"
  sample_tab$cohort[sample_tab$diagnosis == "CS"] <- "CS"
  if (anyNA(sample_tab$cohort)) {
    stop_msg("Unmapped diagnosis: ",
             paste(unique(sample_tab$diagnosis[is.na(sample_tab$cohort)]), collapse = ", "))
  }
  is_conv_suffix <- grepl("_CONV$", sample_tab$sample_id)
  is_conv_diag   <- sample_tab$diagnosis == "Conv"
  stopifnot(identical(is_conv_suffix, is_conv_diag))
  
  samples <- sample_tab$sample_id
  n_samples <- length(samples); n_donors <- length(unique(sample_tab$donor))
  acute_donors   <- unique(sample_tab$donor[sample_tab$cohort == "Acute"])
  conv_donors    <- unique(sample_tab$donor[sample_tab$cohort == "Conv"])
  hc_donors      <- unique(sample_tab$donor[sample_tab$cohort == "HV"])
  surgery_donors <- unique(sample_tab$donor[sample_tab$cohort == "CS"])
  paired_donors  <- intersect(acute_donors, conv_donors)
  stopifnot(n_samples == EXPECTED_N_SAMPLES, n_donors == EXPECTED_N_DONORS)
  stopifnot(length(acute_donors)   == EXPECTED_N_ACUTE_DONORS)
  stopifnot(length(hc_donors)      == EXPECTED_N_HC_DONORS)
  stopifnot(length(surgery_donors) == EXPECTED_N_SURGERY_DONORS)
  stopifnot(sum(sample_tab$cohort == "Conv") == EXPECTED_N_CONV_SAMPLES)
  stopifnot(length(paired_donors) == EXPECTED_N_PAIRED_DONORS)
  stopifnot(all(conv_donors %in% acute_donors))
  
  acute_samp <- setNames(sample_tab$sample_id[sample_tab$cohort == "Acute"],
                         sample_tab$donor[sample_tab$cohort == "Acute"])
  conv_samp  <- setNames(sample_tab$sample_id[sample_tab$cohort == "Conv"],
                         sample_tab$donor[sample_tab$cohort == "Conv"])
  
  # ---------------------------------------------------------------------------
  # 5. SAMPLE-STATE AGGREGATION (all cells, 48 x 24; named ds rows)
  # ---------------------------------------------------------------------------
  n_states <- length(fine_states)
  si_idx <- match(cell_sample, samples); st_idx <- match(cell_state, fine_states)
  stopifnot(!anyNA(si_idx), !anyNA(st_idx))
  ds_idx <- (si_idx - 1L) * n_states + st_idx
  D <- Matrix::sparseMatrix(i = ds_idx, j = seq_along(ds_idx), x = 1,
                            dims = c(n_samples * n_states, length(ds_idx)))
  
  n_cells_ds <- as.integer(Matrix::rowSums(D))
  lib_cell <- as.numeric(Matrix::colSums(cnt))
  lib_ds <- as.numeric(D %*% lib_cell)
  
  ds_names <- paste(rep(samples, each = n_states),
                    rep(fine_states, times = n_samples),
                    sep = "||")
  stopifnot(length(ds_names) == length(n_cells_ds),
            length(ds_names) == length(lib_ds))
  names(n_cells_ds) <- ds_names
  names(lib_ds) <- ds_names
  eligible_ds <- n_cells_ds >= MIN_CELLS_PER_SAMPLE_STATE & lib_ds > 0
  names(eligible_ds) <- ds_names
  
  n_cells_mat <- matrix(n_cells_ds, nrow=n_samples, ncol=n_states, byrow=TRUE,
                        dimnames=list(samples, fine_states))
  elig_mat <- matrix(eligible_ds, nrow=n_samples, ncol=n_states, byrow=TRUE,
                     dimnames=list(samples, fine_states))
  cellfrac_mat <- sweep(n_cells_mat, 1, rowSums(n_cells_mat), "/")
  
  gi <- match(present_feats, rownames(cnt))
  X <- cnt[gi, , drop=FALSE]; rownames(X) <- present_feats
  C_ds <- as.matrix(X %*% Matrix::t(D))              # n_genes x (48*24)
  colnames(C_ds) <- ds_names
  
  # Sample-level (whole-sample) gene UMI for P3A.
  sample_umi <- matrix(0, nrow = length(present_feats), ncol = n_samples,
                       dimnames = list(present_feats, samples))
  for (s in seq_len(n_samples)) {
    cols <- (s - 1L) * n_states + seq_len(n_states)
    sample_umi[, s] <- rowSums(C_ds[, cols, drop = FALSE])
  }
  sample_lib_named <- tapply(lib_ds, rep(samples, each = n_states), sum)
  sample_lib <- as.numeric(sample_lib_named[samples])
  names(sample_lib) <- samples
  stopifnot(!anyNA(sample_lib),
            identical(colnames(sample_umi), names(sample_lib)))
  
  # ---------------------------------------------------------------------------
  # 6. P3A -- whole-sample pseudobulk sanity check
  # ---------------------------------------------------------------------------
  score_up <- module_composite(sample_umi, sample_lib, up_feat, min_genes = 1L)
  score_dn <- module_composite(sample_umi, sample_lib, down_feat, min_genes = 1L)
  m_scores <- list(
    Up29_composite = score_up,
    Down3_composite = score_dn,
    total32_signed = score_up - score_dn,
    M1_cell_cycle = module_composite(sample_umi, sample_lib, READOUT_FEATURE_SETS$M1_cell_cycle,
                                     min_genes = READOUT_MIN_PRESENT[["M1_cell_cycle"]]),
    M2_full = module_composite(sample_umi, sample_lib, READOUT_FEATURE_SETS$M2_full,
                               min_genes = READOUT_MIN_PRESENT[["M2_full"]]),
    primary_core = module_composite(sample_umi, sample_lib, READOUT_FEATURE_SETS$primary_core,
                                    min_genes = READOUT_MIN_PRESENT[["primary_core"]])
  )
  p3a_scores <- do.call(rbind, lapply(names(m_scores), function(mn) {
    data.frame(score = mn, sample_id = samples,
               cohort = sample_tab$cohort[match(samples, sample_tab$sample_id)],
               donor = sample_tab$donor[match(samples, sample_tab$sample_id)],
               value = unname(m_scores[[mn]]), stringsAsFactors = FALSE)
  }))
  write_csv_atomic(p3a_scores, file.path(TABLE_DIR, "p3a_sample_module_scores.csv"))
  
  # Per-gene paired direction (acute vs conv, 9 donors), support-rate aware.
  sig_lc <- log2(sweep(sample_umi, 2, sample_lib, "/") * 1e6 + 1)
  sig_lc[!is.finite(sig_lc)] <- NA_real_
  sig_present <- intersect(sig_feat, present_feats)
  sig_lc <- sig_lc[sig_present, , drop=FALSE]
  sig_gene <- frozen$gene[match(rownames(sig_lc), sig_feat)]
  
  gene_dir_rows <- lapply(seq_len(nrow(sig_lc)), function(i) {
    g <- sig_gene[i]; feat <- rownames(sig_lc)[i]
    dir <- frozen$direction[match(g, frozen$gene)]
    exp_sign <- if (dir == "Up") 1 else -1
    deltas <- vapply(paired_donors, function(d) {
      a <- match(acute_samp[[d]], samples); c <- match(conv_samp[[d]], samples)
      exp_sign * (sig_lc[i, a] - sig_lc[i, c])
    }, numeric(1))
    deltas <- deltas[is.finite(deltas)]
    support_rate <- if (length(deltas)) mean(deltas > 0) else NA_real_
    med <- if (length(deltas)) median(deltas) else NA_real_
    call <- classify_p3a_direction(
      n_pairs = length(deltas),
      med_delta = med,
      positive_fraction = support_rate,
      delta_cut = DIRECTION_DELTA
    )
    
    call_abs0p5 <- classify_p3a_direction(
      n_pairs = length(deltas),
      med_delta = med,
      positive_fraction = support_rate,
      delta_cut = P3A_SENSITIVITY_DELTA
    )
    data.frame(gene = g, dataset_feature = feat, direction = dir,
               n_paired_donors = length(deltas),
               median_signed_delta = med,
               support_rate_signed_gt0 = support_rate,
               direction_call = call,
               direction_call_abs0p5 = call_abs0p5,
               stringsAsFactors = FALSE)
  })
  p3a_gene_dir <- do.call(rbind, gene_dir_rows)
  write_csv_atomic(p3a_gene_dir, file.path(TABLE_DIR, "p3a_gene_direction_consistency.csv"))
  
  # Module summary across conditions.
  p3a_module_summary <- do.call(rbind, lapply(names(m_scores), function(mn) {
    v <- setNames(unname(m_scores[[mn]]), samples)
    cond <- sample_tab$cohort[match(samples, sample_tab$sample_id)]
    paired_d <- vapply(paired_donors, function(d) {
      v[[acute_samp[[d]]]] - v[[conv_samp[[d]]]]
    }, numeric(1))
    paired_d <- paired_d[is.finite(paired_d)]
    data.frame(
      score = mn,
      acute_median = safe_median(v[cond == "Acute"]),
      conv_median  = safe_median(v[cond == "Conv"]),
      hv_median    = safe_median(v[cond == "HV"]),
      cs_median    = safe_median(v[cond == "CS"]),
      paired_delta_median = if (length(paired_d)) median(paired_d) else NA_real_,
      paired_delta_n = length(paired_d),
      acute_minus_hv = safe_median(v[cond == "Acute"]) - safe_median(v[cond == "HV"]),
      acute_minus_cs = safe_median(v[cond == "Acute"]) - safe_median(v[cond == "CS"]),
      stringsAsFactors = FALSE)
  }))
  write_csv_atomic(p3a_module_summary, file.path(TABLE_DIR, "p3a_module_summary.csv"))
  
  # ---------------------------------------------------------------------------
  # 7. P3B -- state-intrinsic module direction (paired, within target states)
  # ---------------------------------------------------------------------------
  p3b_long <- list()
  p3b_cond <- list()
  for (s in TARGET_STATES) {
    msi <- match(s, fine_states)
    col_state <- (seq_len(n_samples) - 1L) * n_states + msi
    elig_s <- elig_mat[, msi]
    cf_s <- cellfrac_mat[, msi]
    for (ro in names(READOUT_FEATURE_SETS)) {
      gs <- READOUT_FEATURE_SETS[[ro]]
      min_ro <- READOUT_MIN_PRESENT[[ro]]
      score_state <- module_composite(C_ds[, col_state, drop = FALSE],
                                      lib_ds[col_state], gs, min_genes = min_ro)
      score_state[!elig_s] <- NA_real_
      names(score_state) <- samples
      
      # Paired-side evaluability diagnostics (before building deltas).
      paired_acute_eval <- vapply(paired_donors, function(d) {
        a <- match(acute_samp[[d]], samples)
        is.finite(score_state[a])
      }, logical(1))
      paired_conv_eval <- vapply(paired_donors, function(d) {
        c <- match(conv_samp[[d]], samples)
        is.finite(score_state[c])
      }, logical(1))
      paired_acute_eval_n <- sum(paired_acute_eval)
      paired_conv_eval_n <- sum(paired_conv_eval)
      
      pd <- lapply(paired_donors, function(d) {
        a <- match(acute_samp[[d]], samples); c <- match(conv_samp[[d]], samples)
        if (is.finite(score_state[a]) && is.finite(score_state[c])) {
          data.frame(donor = d, delta = score_state[a] - score_state[c],
                     acute_value = score_state[a], conv_value = score_state[c],
                     stringsAsFactors = FALSE)
        } else NULL
      })
      pd <- do.call(rbind, pd)
      if (is.null(pd)) {
        pd <- data.frame(donor=character(), delta=numeric(),
                         acute_value=numeric(), conv_value=numeric(),
                         stringsAsFactors=FALSE)
      }
      pd$fine_state <- rep(s, nrow(pd))
      pd$readout <- rep(ro, nrow(pd))
      p3b_long[[length(p3b_long)+1L]] <- pd
      
      cond <- sample_tab$cohort[match(samples, sample_tab$sample_id)]
      med_delta <- if (nrow(pd) >= PAIRED_MIN_N) median(pd$delta) else NA_real_
      support_rate <- if (nrow(pd) >= PAIRED_MIN_N) mean(pd$delta > 0) else NA_real_
      call <- if (nrow(pd) < PAIRED_MIN_N || !is.finite(med_delta) || !is.finite(support_rate)) {
        "not_evaluable"
      } else if (med_delta >= DIRECTION_DELTA && support_rate >= 2/3) {
        "acute_higher_descriptive"
      } else if (med_delta <= -DIRECTION_DELTA && support_rate <= 1/3) {
        "acute_lower_descriptive"
      } else {
        "no_clear_direction"
      }
      not_evaluable_reason <- "evaluable"
      if (identical(call, "not_evaluable")) {
        if (nrow(pd) < PAIRED_MIN_N) {
          not_evaluable_reason <- paste0(
            "paired_complete<", PAIRED_MIN_N,
            "; acute_eval=", paired_acute_eval_n,
            "; conv_eval=", paired_conv_eval_n)
        } else {
          not_evaluable_reason <- "nonfinite_delta_or_support"
        }
      }
      p3b_cond[[length(p3b_cond)+1L]] <- data.frame(
        fine_state = s, readout = ro,
        result_layer = p3b_layer(s, ro),
        acute_median = safe_median(score_state[cond == "Acute"]),
        conv_median  = safe_median(score_state[cond == "Conv"]),
        hv_median    = safe_median(score_state[cond == "HV"]),
        cs_median    = safe_median(score_state[cond == "CS"]),
        n_acute_eval = sum(cond == "Acute" & is.finite(score_state)),
        n_conv_eval  = sum(cond == "Conv" & is.finite(score_state)),
        n_hv_eval    = sum(cond == "HV" & is.finite(score_state)),
        n_cs_eval    = sum(cond == "CS" & is.finite(score_state)),
        paired_delta_n = nrow(pd),
        paired_acute_eval_n = paired_acute_eval_n,
        paired_conv_eval_n = paired_conv_eval_n,
        paired_delta_median = med_delta,
        support_rate_delta_gt0 = support_rate,
        direction_call = call,
        not_evaluable_reason = not_evaluable_reason,
        cell_fraction_acute_median = safe_median(cf_s[cond == "Acute"]),
        cell_fraction_conv_median  = safe_median(cf_s[cond == "Conv"]),
        cell_fraction_hv_median    = safe_median(cf_s[cond == "HV"]),
        cell_fraction_cs_median    = safe_median(cf_s[cond == "CS"]),
        n_cells_acute_median = safe_median(n_cells_mat[cond == "Acute", msi]),
        n_cells_conv_median  = safe_median(n_cells_mat[cond == "Conv", msi]),
        n_cells_hv_median    = safe_median(n_cells_mat[cond == "HV", msi]),
        n_cells_cs_median    = safe_median(n_cells_mat[cond == "CS", msi]),
        stringsAsFactors = FALSE)
    }
  }
  p3b_paired <- do.call(rbind, p3b_long)
  p3b_cond_df <- do.call(rbind, p3b_cond)
  write_csv_atomic(p3b_paired, file.path(TABLE_DIR, "p3b_paired_delta_long.csv"))
  write_csv_atomic(p3b_cond_df, file.path(TABLE_DIR, "p3b_state_readout_condition_summary.csv"))
  
  p3b_calls <- p3b_cond_df[, c(
    "fine_state","readout","result_layer",
    "paired_delta_n",
    "paired_acute_eval_n","paired_conv_eval_n",
    "paired_delta_median","support_rate_delta_gt0",
    "direction_call","not_evaluable_reason",
    "cell_fraction_acute_median","cell_fraction_conv_median",
    "n_cells_acute_median","n_cells_conv_median")]
  write_csv_atomic(p3b_calls, file.path(TABLE_DIR, "p3b_direction_calls.csv"))
  
  # ---------------------------------------------------------------------------
  # 8. PROVENANCE + CONSOLE SUMMARY
  # ---------------------------------------------------------------------------
  provenance <- data.frame(
    field = c("script","version","date","input_file",
              "min_cells_per_sample_state","paired_min_n","direction_delta",
              "p3a_aligned_fraction_min","p3a_opposite_fraction_max",
              "p3a_sensitivity_delta",
              "readout_min_present",
              "n_samples","n_acute","n_conv","n_hv","n_cs","n_paired",
              "n32_genes_present","target_states","readouts","R_version"),
    value = c(SCRIPT_NAME, SCRIPT_VERSION, as.character(Sys.time()), INPUT_FILE,
              as.character(MIN_CELLS_PER_SAMPLE_STATE),
              as.character(PAIRED_MIN_N), as.character(DIRECTION_DELTA),
              as.character(P3A_ALIGNED_FRAC_MIN),
              as.character(P3A_OPPOSITE_FRAC_MAX),
              as.character(P3A_SENSITIVITY_DELTA),
              paste(names(READOUT_MIN_PRESENT), READOUT_MIN_PRESENT,
                    sep = "=", collapse = ";"),
              as.character(n_samples), as.character(length(acute_donors)),
              as.character(length(conv_donors)), as.character(length(hc_donors)),
              as.character(length(surgery_donors)), as.character(length(paired_donors)),
              as.character(n32_present),
              paste(TARGET_STATES, collapse = ";"),
              paste(names(READOUT_FEATURE_SETS), collapse = ";"),
              as.character(getRversion())),
    stringsAsFactors = FALSE)
  provenance <- rbind(
    provenance,
    data.frame(
      field = c("state_file", "frozen_file", "map_file",
                "p3a_classification_note"),
      value = c(
        STATE_FILE, FROZEN_FILE, MAP_FILE,
        "Revised descriptive classification; symmetric +/-0.3 primary rule and +/-0.5 sensitivity."
      ),
      stringsAsFactors = FALSE
    )
  )
  write_csv_atomic(provenance, file.path(TABLE_DIR, "run_provenance.csv"))
  
  cat("Script 05e v1.5 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("32-gene features present: 32 / 32\n")
  cat("\n--- Readout membership ---\n")
  print(readout_membership[, c("readout","n_requested","n_present",
                               "min_present_required","membership_status",
                               "missing_features")], row.names = FALSE)
  
  cat("\n--- P3A module summary (whole-sample composite; composition-sensitive sanity) ---\n")
  print(p3a_module_summary, row.names = FALSE)
  
  cat("\n--- P3A gene direction consistency (acute vs conv, paired) ---\n")
  up_ok <- p3a_gene_dir$direction_call == "aligned" &
    p3a_gene_dir$direction == "Up"
  
  dn_ok <- p3a_gene_dir$direction_call == "aligned" &
    p3a_gene_dir$direction == "Down"
  cat("Up genes aligned with the bulk-defined direction:", sum(up_ok, na.rm=TRUE), "/", sum(p3a_gene_dir$direction == "Up"), "\n")
  cat("Down genes aligned with the bulk-defined direction:", sum(dn_ok, na.rm=TRUE), "/", sum(p3a_gene_dir$direction == "Down"),
      "\n")
  cat("P3A direction-call distribution (+/-0.3):\n")
  print(table(p3a_gene_dir$direction_call, useNA = "ifany"))
  
  cat("P3A sensitivity distribution (+/-0.5):\n")
  print(table(p3a_gene_dir$direction_call_abs0p5, useNA = "ifany"))
  print(p3a_gene_dir[, c(
    "gene", "direction", "n_paired_donors",
    "median_signed_delta", "support_rate_signed_gt0",
    "direction_call", "direction_call_abs0p5"
  )], row.names = FALSE)
  
  cat("\n--- P3B direction calls (state-intrinsic; primary rows first) ---\n")
  layer_order <- c("primary_M2_state_intrinsic",
                   "secondary_M2_state_context",
                   "supplementary_readout_context")
  p3b_calls_print <- p3b_calls[order(
    factor(p3b_calls$result_layer, levels = layer_order),
    p3b_calls$fine_state, p3b_calls$readout), , drop = FALSE]
  print(p3b_calls_print, row.names = FALSE)
  
  cat("\n--- P3B primary M2/primary_core rows ---\n")
  foc <- p3b_cond_df[p3b_cond_df$result_layer == "primary_M2_state_intrinsic", , drop = FALSE]
  print(foc[, c("fine_state","readout",
                "acute_median","conv_median",
                "paired_delta_median","support_rate_delta_gt0",
                "direction_call","not_evaluable_reason",
                "cell_fraction_acute_median","cell_fraction_conv_median",
                "n_cells_acute_median","n_cells_conv_median")], row.names = FALSE)
  
  cat("\nNote: P3A is an internal sanity check only (composition-sensitive); not matched-bulk validation.\n")
  cat("P3B direction calls are descriptive flags (median + support rate), not significance tests.\n")
  cat("P3A total32_signed is an internal non-inversion sanity score, not a validated quantitative MARS score.\n")
})
