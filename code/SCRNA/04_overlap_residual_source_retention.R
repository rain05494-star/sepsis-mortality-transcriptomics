
# ============================================================================
# SCRIPT 04 v2.7 -- Overlap-exclusion single-cell source attribution and
#                   residual-gene prioritization
#   Purpose: single-cell extension of the bulk Giannini-E3 overlap-exclusion
#   sensitivity analysis. Within functional modules (cell-cycle, neutrophil),
#   compares full / E3-overlap / E3-residual gene sets for retained cell-state
#   localization; classifies residual genes by information novelty (referenced
#   against the OVERLAP set, with STRICT reference-reliability gating: if the
#   overlap reference is not formal, no retention/refinement/novelty is called);
#   assigns single-cell evidence candidate tiers (require >=8 eligible donors).
#   Other/Inflammatory are per-gene. All inputs for each stage read from a
#   SINGLE locked table directory. Cross-source consistency (frozen vs pgs,
#   evaluability) and data completeness (gene coverage) are fail-closed.
#   Outputs:
#     table0_signature_overlap_residual_membership.csv
#     table0b_module_overlap_residual_set_sizes.csv
#     table1_module_set_source_summary.csv
#     table2_all_gene_candidate.csv          (single-cell evidence scope)
#     table2_residual_gene_candidate.csv
#     fig04_profile_source_data.csv
#     fig04_residual_evidence_source_data.csv
#             (.svg if svglite available; .tiff if ragg available)
# RUN: Rscript 04_overlap_residual_source_retention.R
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
  
  STAGE_NAME  <- "04_overlap_residual_source_retention"
  RUN_PURPOSE <- "overlap_residual_source_attribution"
  
  MODULE_FULL <- c("Cell cycle", "Neutrophil")   # set-level modules
  OVERLAP_TAG  <- "E3_overlap_15"
  RESIDUAL_TAG <- "E3_nonoverlap_residual_17"
  MIN_SET_FORMAL <- 5L
  
  THRESH_DONOR_RATE      <- 0.50
  THRESH_ELIGIBLE_DONOR_N <- 8L
  THRESH_AMBIGUOUS_DELTA <- 0.05
  RHO_EVALUABLE_SD       <- 1e-6
  
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
  safe_sd <- function(x){x<-x[is.finite(x)]; if(length(x)<2) NA_real_ else sd(x)}
  one_chr <- function(x) {
    if (length(x) < 1L || is.na(x[1])) return(NA_character_)
    as.character(x[1])
  }
  is_true <- function(x) {
    z <- one_chr(x)
    !is.na(z) && tolower(z) %in% c("true", "t", "1", "yes")
  }
  is_stable_detected <- function(x) {
    identical(one_chr(x), "stable_detected")
  }
  is_low_or_missing_conf <- function(x) {
    z <- one_chr(x)
    is.na(z) || tolower(z) == "low"
  }
  coerce_numeric_cols <- function(df, cols, label) {
    for (cn in cols) {
      old <- df[[cn]]
      new <- suppressWarnings(as.numeric(old))
      bad <- is.na(new) & !is.na(old) & trimws(as.character(old)) != ""
      if (any(bad)) stop_msg(label, " has non-numeric values in column: ", cn)
      df[[cn]] <- new
    }
    df
  }
  safe_min_int <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) NA_integer_ else as.integer(min(x))
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
  # Lock a SINGLE table directory under `root` that contains ALL required files,
  # ranked by the newest mtime among the required files (avoids cross-run mixing
  # and is robust to directory-level touch operations).
  find_latest_table_dir <- function(root, required_files) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    anchors <- list.files(root, pattern = paste0("^", required_files[1], "$"),
                          recursive = TRUE, full.names = TRUE)
    dirs <- unique(dirname(anchors))
    if (!length(dirs)) stop_msg("No candidate table directory found under ", root)
    ok <- vapply(dirs, function(d) all(file.exists(file.path(d, required_files))), logical(1))
    dirs <- dirs[ok]
    if (!length(dirs)) {
      stop_msg("No complete table directory found under ", root,
               " with required files: ", paste(required_files, collapse = ", "))
    }
    score <- vapply(dirs, function(d) {
      max(as.numeric(file.info(file.path(d, required_files))$mtime), na.rm = TRUE)
    }, numeric(1))
    dirs[which.max(score)]
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
  short_state <- function(x) {
    z <- as.character(x)
    map <- c(
      "HSPCs"="HSPCs",
      "Cycling_neutrophil_progenitors"="Cycling neut prog",
      "MPO+_immature_neutrophils_or_progenitors"="MPO+ imm neut/prog",
      "PADI4+_immature_neutrophils"="PADI4+ imm neut",
      "IL1R2+_immature_neutrophils"="IL1R2+ imm neut",
      "S100A8-9_hi_neutrophils"="S100A8/9 hi neut",
      "Mature_neutrophils"="Mature neut",
      "Degranulating_neutrophils"="Degran neut",
      "Apoptosing_neutrophils"="Apoptosing neut",
      "Eosinophils"="Eosinophils","Mast_cells/eosiniophils"="Mast/eos",
      "Classical_monocytes"="Classical mono","Non-classical_monocytes"="Non-classical mono",
      "cDCs"="cDCs","pDCs"="pDCs",
      "Naive_CD4_T_cells"="Naive CD4 T","Memory_CD4_T_cells"="Memory CD4 T",
      "Naive_CD8_T_cells"="Naive CD8 T","CD8_T_cells"="CD8 T",
      "Cycling_TNK"="Cycling TNK","NK cells"="NK cells","B_cells"="B cells",
      "Plasmablasts"="Plasmablasts","Platelets"="Platelets")
    ifelse(z %in% names(map), unname(map[z]), z)
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. READ FROZEN INPUTS (single locked table dir per stage)
  # ---------------------------------------------------------------------------
  S01_TABLE_DIR <- find_latest_table_dir(
    SCRIPT01_ROOT,
    c("frozen_signature_input.csv","fine_state_display_order_audit.csv",
      "final_gene_evaluability_32.csv"))
  S02_TABLE_DIR <- find_latest_table_dir(
    SCRIPT02_ROOT,
    c("per_gene_cell_source_summary.csv","median_state_logCPM_matrix.csv",
      "donor_state_gene_counts_long.csv"))
  
  FROZEN_FILE <- file.path(S01_TABLE_DIR, "frozen_signature_input.csv")
  STATE_FILE  <- file.path(S01_TABLE_DIR, "fine_state_display_order_audit.csv")
  EVAL_FILE   <- file.path(S01_TABLE_DIR, "final_gene_evaluability_32.csv")
  PGS_FILE    <- file.path(S02_TABLE_DIR, "per_gene_cell_source_summary.csv")
  LCPM_FILE   <- file.path(S02_TABLE_DIR, "median_state_logCPM_matrix.csv")
  LONG_FILE   <- file.path(S02_TABLE_DIR, "donor_state_gene_counts_long.csv")
  
  frozen <- read.csv(FROZEN_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  eva01 <- read.csv(EVAL_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  pgs <- read.csv(PGS_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  lcpm <- read.csv(LCPM_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  long <- read.csv(LONG_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  
  require_cols(frozen, c("gene","module","direction","canonical_hgnc","E3_membership"), "frozen_signature_input.csv")
  require_cols(state01, c("display_order","fine_state","lineage","color"), "fine_state_display_order_audit.csv")
  require_cols(eva01, c("signature_gene","class"), "final_gene_evaluability_32.csv")
  require_cols(pgs, c("signature_gene","module","direction","E3_membership","evaluable_class",
                      "main_contributor_state","main_contributor_lineage",
                      "main_contribution_fraction_median_all26","second_contributor_state",
                      "contribution_delta_main_minus_second","source_call_ambiguous",
                      "source_call_confidence","main_state_donor_detection_rate_all26",
                      "main_state_eligible_donor_n"), "per_gene_cell_source_summary.csv")
  require_cols(long, c("gene","donor","fine_state","umi","donor_gene_total_umi"), "donor_state_gene_counts_long.csv")
  
  # ---- gene identity fail-closed (dup/missing would silently bias 15/17 split)
  if (anyDuplicated(frozen$gene)) {
    stop_msg("Duplicated genes in frozen_signature_input.csv: ",
             paste(unique(frozen$gene[duplicated(frozen$gene)]), collapse=", "))
  }
  if (anyDuplicated(pgs$signature_gene)) {
    stop_msg("Duplicated signature_gene in per_gene_cell_source_summary.csv: ",
             paste(unique(pgs$signature_gene[duplicated(pgs$signature_gene)]), collapse=", "))
  }
  missing_long_genes <- setdiff(frozen$gene, unique(long$gene))
  if (length(missing_long_genes)) {
    stop_msg("Genes missing from donor_state_gene_counts_long.csv: ",
             paste(missing_long_genes, collapse=", "))
  }
  
  # explicit numeric coercion on pgs key columns
  pgs_num_cols <- c("main_contribution_fraction_median_all26",
                    "contribution_delta_main_minus_second",
                    "main_state_donor_detection_rate_all26",
                    "main_state_eligible_donor_n")
  pgs <- coerce_numeric_cols(pgs, pgs_num_cols, "per_gene_cell_source_summary.csv")
  
  # force numeric on count columns
  long$umi <- suppressWarnings(as.numeric(long$umi))
  long$donor_gene_total_umi <- suppressWarnings(as.numeric(long$donor_gene_total_umi))
  if (any(!is.finite(long$umi))) stop_msg("Non-finite umi values in donor_state_gene_counts_long.csv")
  if (any(!is.finite(long$donor_gene_total_umi))) {
    stop_msg("Non-finite donor_gene_total_umi values in donor_state_gene_counts_long.csv")
  }
  
  # donor_gene_total_umi must be unique per donor-gene across states
  tot_audit <- unique(long[, c("donor","gene","donor_gene_total_umi"), drop=FALSE])
  tot_key <- paste(tot_audit$donor, tot_audit$gene, sep="||")
  tot_n <- tapply(tot_audit$donor_gene_total_umi, tot_key, function(z) length(unique(z[is.finite(z)])))
  bad_tot <- names(tot_n)[tot_n > 1L]
  if (length(bad_tot)) {
    stop_msg("Conflicting donor_gene_total_umi for donor-gene pairs: ",
             paste(head(bad_tot, 20), collapse=", "),
             if (length(bad_tot) > 20) " ..." else "")
  }
  
  stopifnot(nrow(frozen) == 32L, nrow(pgs) == 32L, nrow(state01) == 24L)
  fine_states <- state01$fine_state[order(state01$display_order)]
  state_lineage <- setNames(state01$lineage, state01$fine_state)
  state_colors  <- setNames(state01$color, state01$fine_state)
  donors <- sort(unique(long$donor))
  stopifnot(length(donors) == 26L)
  
  frozen$module_s <- module_short(frozen$module)
  pgs$module_s <- module_short(pgs$module)
  frozen$E3 <- frozen$E3_membership
  pgs$E3 <- pgs$E3_membership
  
  # ---- E3 membership fail-closed (lock the 15/17 split) ----
  if (!all(frozen$E3 %in% c(OVERLAP_TAG, RESIDUAL_TAG))) {
    stop_msg("Unexpected E3_membership labels: ",
             paste(setdiff(unique(frozen$E3), c(OVERLAP_TAG, RESIDUAL_TAG)), collapse = ", "))
  }
  stopifnot(sum(frozen$E3 == OVERLAP_TAG) == 15L)
  stopifnot(sum(frozen$E3 == RESIDUAL_TAG) == 17L)
  stopifnot(length(unique(frozen$gene)) == 32L)
  
  # frozen vs pgs agree on gene/module/E3/direction
  m_fp <- match(pgs$signature_gene, frozen$gene)
  stopifnot(all(!is.na(m_fp)))
  stopifnot(identical(pgs$module_s, frozen$module_s[m_fp]))
  stopifnot(identical(pgs$E3, frozen$E3[m_fp]))
  stopifnot(identical(pgs$direction, frozen$direction[m_fp]))
  
  # evaluability class consistency: Script01 frozen audit vs Script02 summary
  norm_class <- function(x) {
    gsub("[[:space:]-]+", "_", tolower(trimws(as.character(x))))
  }
  m_ep <- match(pgs$signature_gene, eva01$signature_gene)
  stopifnot(all(!is.na(m_ep)))
  eva_class_norm <- norm_class(eva01$class[m_ep])
  pgs_class_norm <- norm_class(pgs$evaluable_class)
  if (!identical(eva_class_norm, pgs_class_norm)) {
    bad <- pgs$signature_gene[eva_class_norm != pgs_class_norm]
    stop_msg("Evaluability class disagrees between Script01 and Script02 for: ",
             paste(bad, collapse=", "))
  }
  
  # Table 0: locked 15/17 gene identity (Methods + supplement)
  membership_table <- frozen[, c("gene","canonical_hgnc","module_s","direction","E3"), drop=FALSE]
  names(membership_table) <- c("gene","canonical_hgnc","module","direction","E3_membership")
  write_csv_atomic(membership_table, file.path(TABLE_DIR, "table0_signature_overlap_residual_membership.csv"))
  
  # set gene lists
  sets <- list()
  for (mod in MODULE_FULL) {
    sets[[paste0("full_", mod)]]     <- sort(frozen$gene[frozen$module_s == mod])
    sets[[paste0("overlap_", mod)]]  <- sort(frozen$gene[frozen$module_s == mod & frozen$E3 == OVERLAP_TAG])
    sets[[paste0("residual_", mod)]] <- sort(frozen$gene[frozen$module_s == mod & frozen$E3 == RESIDUAL_TAG])
  }
  stopifnot(length(sets[[paste0("full_","Cell cycle")]]) == 11L)
  stopifnot(length(sets[[paste0("full_","Neutrophil")]]) == 10L)
  stopifnot(all(vapply(sets, length, integer(1)) > 0L))
  
  # Table 0b: explicit set sizes per module (supplement)
  set_size_table <- do.call(rbind, lapply(MODULE_FULL, function(mod) {
    data.frame(
      module = mod,
      full_n = length(sets[[paste0("full_", mod)]]),
      overlap_n = length(sets[[paste0("overlap_", mod)]]),
      residual_n = length(sets[[paste0("residual_", mod)]]),
      overlap_scope = if (length(sets[[paste0("overlap_", mod)]]) >= MIN_SET_FORMAL) {
        "formal_or_check_delta"
      } else {
        "descriptive_small_gene_set"
      },
      residual_scope = if (length(sets[[paste0("residual_", mod)]]) >= MIN_SET_FORMAL) {
        "formal_set_level"
      } else {
        "descriptive_small_gene_set"
      },
      stringsAsFactors = FALSE
    )
  }))
  write_csv_atomic(set_size_table, file.path(TABLE_DIR, "table0b_module_overlap_residual_set_sizes.csv"))
  
  # ---------------------------------------------------------------------------
  # 3. SET-LEVEL DONOR-STATE CONTRIBUTION PROFILES
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
  set_cf <- lapply(sets, build_set_cf)
  
  set_main_state_donor_rate <- function(cf, main_s) {
    if (is.na(main_s) || !main_s %in% colnames(cf)) return(NA_real_)
    x <- cf[, main_s]
    mean(is.finite(x) & x > 0)
  }
  set_main_state_signal_donor_n <- function(cf, main_s) {
    if (is.na(main_s) || !main_s %in% colnames(cf)) return(NA_integer_)
    x <- cf[, main_s]
    sum(is.finite(x) & x > 0)
  }
  
  set_profile <- function(cf) {
    p <- vapply(fine_states, function(s) safe_median(cf[, s]), numeric(1))
    names(p) <- fine_states
    p
  }
  set_profiles <- lapply(set_cf, set_profile)
  
  # logCPM matrix: numeric coercion + full 32-gene coverage (fail-closed)
  require_cols(lcpm, c("gene", fine_states), "median_state_logCPM_matrix.csv")
  lcpm <- coerce_numeric_cols(lcpm, fine_states, "median_state_logCPM_matrix.csv")
  if (anyDuplicated(lcpm$gene)) {
    stop_msg("Duplicated genes in median_state_logCPM_matrix.csv: ",
             paste(unique(lcpm$gene[duplicated(lcpm$gene)]), collapse=", "))
  }
  missing_lcpm_genes <- setdiff(frozen$gene, lcpm$gene)
  if (length(missing_lcpm_genes)) {
    stop_msg("Genes missing from median_state_logCPM_matrix.csv: ",
             paste(missing_lcpm_genes, collapse=", "))
  }
  rownames(lcpm) <- lcpm$gene
  lcpm_mat <- as.matrix(lcpm[, fine_states, drop=FALSE])
  rownames(lcpm_mat) <- lcpm$gene
  set_logcpm_profile <- function(genes) {
    gs <- intersect(genes, rownames(lcpm_mat))
    if (!length(gs)) return(rep(NA_real_, length(fine_states)))
    vapply(fine_states, function(s) safe_median(lcpm_mat[gs, s]), numeric(1))
  }
  set_logcpm <- lapply(sets, set_logcpm_profile)
  
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
    if (safe_sd(p1[ok]) < RHO_EVALUABLE_SD || safe_sd(p2[ok]) < RHO_EVALUABLE_SD)
      return(NA_real_)  # correlation not evaluable (near-constant profile)
    cor(p1[ok], p2[ok], method="spearman")
  }
  
  # ---------------------------------------------------------------------------
  # 4. TABLE 1: module-set source summary (full / overlap / residual)
  # ---------------------------------------------------------------------------
  t1_rows <- list()
  for (mod in MODULE_FULL) {
    full_k <- paste0("full_", mod); ov_k <- paste0("overlap_", mod); res_k <- paste0("residual_", mod)
    prof_full <- set_profiles[[full_k]]
    prof_ov   <- set_profiles[[ov_k]]
    lc_full   <- set_logcpm[[full_k]]
    lc_ov     <- set_logcpm[[ov_k]]
    ov_delta  <- delta12(prof_ov)
    for (k in c(full_k, ov_k, res_k)) {
      prof <- set_profiles[[k]]
      lc <- set_logcpm[[k]]
      main_s <- top_state(prof)
      main_lg <- if (is.na(main_s)) NA_character_ else unname(state_lineage[main_s])
      top1_cf <- if (is.na(main_s)) NA_real_ else prof[main_s]
      del <- delta12(prof)
      gs <- sets[[k]]
      g_pgs <- pgs[pgs$signature_gene %in% gs, , drop=FALSE]
      t1_rows[[length(t1_rows)+1L]] <- data.frame(
        module = mod,
        set = k,
        n_genes = length(gs),
        set_analysis_scope = if (length(gs) >= MIN_SET_FORMAL) "formal_set_level" else "descriptive_small_gene_set",
        main_state = main_s,
        main_lineage = main_lg,
        top1_contribution_median = top1_cf,
        top1_minus_top2_delta = del,
        set_source_ambiguous = if (is.finite(del)) del < THRESH_AMBIGUOUS_DELTA else NA,
        rho_contribution_vs_full = if (k == full_k) 1 else spearman_pair(prof_full, prof),
        rho_logCPM_vs_full = if (k == full_k) 1 else spearman_pair(lc_full, lc),
        rho_contribution_vs_overlap = if (k == ov_k) 1 else spearman_pair(prof_ov, prof),
        rho_logCPM_vs_overlap = if (k == ov_k) 1 else spearman_pair(lc_ov, lc),
        set_main_state_donor_detection_rate = set_main_state_donor_rate(set_cf[[k]], main_s),
        set_main_state_signal_donor_n = set_main_state_signal_donor_n(set_cf[[k]], main_s),
        member_gene_main_state_donor_detection_median = safe_median(g_pgs$main_state_donor_detection_rate_all26),
        member_gene_main_state_eligible_donor_n_min = if (nrow(g_pgs)) {
          safe_min_int(g_pgs$main_state_eligible_donor_n)
        } else {
          NA_integer_
        },
        retained_full_source = if (k == full_k) NA else (main_s == top_state(prof_full)),
        retained_overlap_source = if (k %in% c(full_k, ov_k)) NA else (main_s == top_state(prof_ov)),
        overlap_reference_n_genes = length(sets[[ov_k]]),
        overlap_reference_delta = ov_delta,
        overlap_reference_scope = if (
          length(sets[[ov_k]]) >= MIN_SET_FORMAL &&
          is.finite(ov_delta) && ov_delta >= THRESH_AMBIGUOUS_DELTA
        ) "formal_reference" else "descriptive_reference",
        stringsAsFactors=FALSE)
    }
  }
  t1 <- do.call(rbind, t1_rows)
  write_csv_atomic(t1, file.path(TABLE_DIR, "table1_module_set_source_summary.csv"))
  
  # ---------------------------------------------------------------------------
  # 5. TABLE 2 (all 32 + residual subset): novelty (vs OVERLAP) + candidate tier
  # ---------------------------------------------------------------------------
  module_overlap_main_state <- setNames(vapply(MODULE_FULL, function(mod)
    top_state(set_profiles[[paste0("overlap_", mod)]]), character(1)), MODULE_FULL)
  module_overlap_main_lineage <- setNames(vapply(MODULE_FULL, function(mod) {
    ms <- module_overlap_main_state[[mod]]
    if (is.na(ms)) NA_character_ else unname(state_lineage[ms])
  }, character(1)), MODULE_FULL)
  module_overlap_delta <- setNames(vapply(MODULE_FULL, function(mod) {
    delta12(set_profiles[[paste0("overlap_", mod)]])
  }, numeric(1)), MODULE_FULL)
  module_overlap_scope <- setNames(vapply(MODULE_FULL, function(mod) {
    ov_k <- paste0("overlap_", mod)
    ov_delta <- delta12(set_profiles[[ov_k]])
    if (length(sets[[ov_k]]) >= MIN_SET_FORMAL &&
        is.finite(ov_delta) && ov_delta >= THRESH_AMBIGUOUS_DELTA) {
      "formal_reference"
    } else {
      "descriptive_reference"
    }
  }, character(1)), MODULE_FULL)
  
  classify_novelty <- function(g_pgs) {
    if (!is_stable_detected(g_pgs$evaluable_class)) return("low_detected_descriptive_only")
    if (is_true(g_pgs$source_call_ambiguous)) return("ambiguous")
    
    mod <- one_chr(g_pgs$module_s)
    if (!mod %in% MODULE_FULL) return("per_gene_module_not_set_level")
    
    ref_st <- module_overlap_main_state[[mod]]
    ref_lg <- module_overlap_main_lineage[[mod]]
    ref_formal <- identical(module_overlap_scope[[mod]], "formal_reference")
    
    g_st <- one_chr(g_pgs$main_contributor_state)
    g_lg <- one_chr(g_pgs$main_contributor_lineage)
    
    if (is.na(g_st) || is.na(g_lg) || is.na(ref_st) || is.na(ref_lg)) {
      return("ambiguous")
    }
    
    # STRICT reference gating: if the overlap reference is not formal, do not
    # call retention/refinement/novelty against it.
    if (!ref_formal) {
      return("overlap_reference_descriptive")
    }
    
    if (identical(g_st, ref_st)) return("retains_overlap_source")
    if (identical(g_lg, ref_lg)) return("refines_same_lineage")
    
    rate  <- as.numeric(g_pgs$main_state_donor_detection_rate_all26)
    delta <- as.numeric(g_pgs$contribution_delta_main_minus_second)
    elig  <- as.numeric(g_pgs$main_state_eligible_donor_n)
    conf  <- g_pgs$source_call_confidence
    
    if (is.finite(rate) && rate >= THRESH_DONOR_RATE &&
        is.finite(delta) && delta >= THRESH_AMBIGUOUS_DELTA &&
        is.finite(elig) && elig >= THRESH_ELIGIBLE_DONOR_N &&
        !is_low_or_missing_conf(conf)) {
      return("adds_new_cell_state")
    }
    
    "ambiguous"
  }
  
  novelty_for_all <- function(g_pgs) {
    if (identical(one_chr(g_pgs$E3), OVERLAP_TAG)) return("overlap_member")
    classify_novelty(g_pgs)
  }
  
  candidate_tier <- function(g_pgs, novelty) {
    stable <- is_stable_detected(g_pgs$evaluable_class)
    rate   <- as.numeric(g_pgs$main_state_donor_detection_rate_all26)
    delta  <- as.numeric(g_pgs$contribution_delta_main_minus_second)
    elig   <- as.numeric(g_pgs$main_state_eligible_donor_n)
    
    conf_bad <- is_low_or_missing_conf(g_pgs$source_call_confidence)
    amb      <- is_true(g_pgs$source_call_ambiguous)
    elig_ok  <- is.finite(elig) && elig >= THRESH_ELIGIBLE_DONOR_N
    
    if (!stable) return("Tier3_low_detected_descriptive")
    if (amb || conf_bad || !is.finite(rate) || rate < 0.30 || !elig_ok) return("Tier3_weak_source")
    if (identical(novelty, "overlap_reference_descriptive")) return("Tier2_supportive")
    
    if (rate >= THRESH_DONOR_RATE &&
        is.finite(delta) && delta >= THRESH_AMBIGUOUS_DELTA &&
        novelty %in% c("retains_overlap_source","adds_new_cell_state","refines_same_lineage",
                       "per_gene_module_not_set_level","overlap_member"))
      return("Tier1_priority")
    "Tier2_supportive"
  }
  
  build_candidate_row <- function(gp) {
    nov <- novelty_for_all(gp)
    tier <- candidate_tier(gp, nov)
    data.frame(
      gene = gp$signature_gene,
      module = gp$module_s,
      direction = gp$direction,
      E3_membership = gp$E3,
      evaluable_class = gp$evaluable_class,
      main_contributor_state = gp$main_contributor_state,
      main_contributor_lineage = gp$main_contributor_lineage,
      second_contributor_state = gp$second_contributor_state,
      top1_contribution_median = gp$main_contribution_fraction_median_all26,
      top1_minus_top2_delta = gp$contribution_delta_main_minus_second,
      main_state_donor_detection_rate = gp$main_state_donor_detection_rate_all26,
      main_state_eligible_donor_n = gp$main_state_eligible_donor_n,
      source_call_confidence = gp$source_call_confidence,
      source_call_ambiguous = gp$source_call_ambiguous,
      overlap_set_main_state = if (gp$module_s %in% MODULE_FULL) {
        module_overlap_main_state[[gp$module_s]]
      } else {
        NA_character_
      },
      overlap_reference_delta = if (gp$module_s %in% MODULE_FULL) {
        module_overlap_delta[[gp$module_s]]
      } else {
        NA_real_
      },
      overlap_reference_scope = if (gp$module_s %in% MODULE_FULL) {
        module_overlap_scope[[gp$module_s]]
      } else {
        NA_character_
      },
      novelty_class = nov,
      candidate_tier = tier,
      candidate_evidence_scope = "single_cell_source_only",
      stringsAsFactors=FALSE)
  }
  
  t2_all <- do.call(rbind, lapply(seq_len(nrow(pgs)), function(i) build_candidate_row(pgs[i, , drop=FALSE])))
  write_csv_atomic(t2_all, file.path(TABLE_DIR, "table2_all_gene_candidate.csv"))
  
  t2_res <- t2_all[t2_all$E3_membership == RESIDUAL_TAG, , drop=FALSE]
  write_csv_atomic(t2_res, file.path(TABLE_DIR, "table2_residual_gene_candidate.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. FIGURE SOURCE DATA (always exported, independent of plotting packages)
  # ---------------------------------------------------------------------------
  profile_fig_df <- do.call(rbind, lapply(MODULE_FULL, function(mod) {
    do.call(rbind, lapply(c("full","overlap","residual"), function(tag) {
      k <- paste0(tag, "_", mod)
      data.frame(module=mod, set=tag,
                 fine_state=fine_states, fine_state_short=short_state(fine_states),
                 contribution=unname(set_profiles[[k]]), stringsAsFactors=FALSE)
    }))
  }))
  write_csv_atomic(profile_fig_df, file.path(TABLE_DIR, "fig04_profile_source_data.csv"))
  
  nov_labels <- c(
    "retains_overlap_source"        = "Retains overlap source",
    "refines_same_lineage"          = "Same lineage, new state",
    "adds_new_cell_state"           = "New cell state",
    "ambiguous"                     = "Ambiguous",
    "low_detected_descriptive_only" = "Low-detected",
    "per_gene_module_not_set_level" = "Per-gene module",
    "overlap_reference_descriptive" = "Reference not formal"
  )
  tier_labels <- c(
    "Tier1_priority"                 = "Tier 1",
    "Tier2_supportive"               = "Tier 2",
    "Tier3_low_detected_descriptive" = "Low-detected",
    "Tier3_weak_source"              = "Weak source"
  )
  
  nov_df <- t2_res[!is.na(t2_res$novelty_class), , drop=FALSE]
  evidence_long <- rbind(
    data.frame(gene=nov_df$gene, module=nov_df$module,
               panel=factor("Novelty", levels=c("Novelty","Tier")),
               value=nov_df$novelty_class, stringsAsFactors=FALSE),
    data.frame(gene=nov_df$gene, module=nov_df$module,
               panel=factor("Tier", levels=c("Novelty","Tier")),
               value=nov_df$candidate_tier, stringsAsFactors=FALSE)
  )
  evidence_long$value_label <- ifelse(
    evidence_long$panel == "Novelty",
    nov_labels[evidence_long$value],
    tier_labels[evidence_long$value]
  )
  if (any(is.na(evidence_long$value_label))) {
    stop_msg("Unmapped evidence labels: ",
             paste(unique(evidence_long$value[is.na(evidence_long$value_label)]), collapse=", "))
  }
  write_csv_atomic(evidence_long, file.path(TABLE_DIR, "fig04_residual_evidence_source_data.csv"))
  
  # ============================================================================
  # 7. Figure: overlap/residual source retention
  #    Fresh / natural / bright publication-style layout
  # ============================================================================
  
  PLOT_PKGS <- c("ggplot2", "grid")
  HAS_PLOT <- all(vapply(PLOT_PKGS, requireNamespace, logical(1), quietly = TRUE))
  
  if (HAS_PLOT) {
    suppressPackageStartupMessages({
      library(ggplot2)
      library(grid)
    })
    
    # --------------------------------------------------------------------------
    # 7.1 Global plotting style
    # --------------------------------------------------------------------------
    
    nature_theme <- theme_classic(base_size = 8, base_family = "Arial") +
      theme(
        axis.line = element_line(size = 0.35, colour = "#2A2A2A"),
        axis.ticks = element_line(size = 0.30, colour = "#2A2A2A"),
        axis.text = element_text(colour = "#222222"),
        axis.title = element_text(colour = "#222222"),
        plot.title = element_text(
          face = "bold", size = 10.5, colour = "#111111",
          margin = margin(b = 1.5)
        ),
        plot.subtitle = element_text(
          size = 7.2, colour = "#555555",
          margin = margin(b = 5)
        ),
        legend.title = element_text(size = 7.0, colour = "#222222"),
        legend.text = element_text(size = 6.5, colour = "#222222"),
        legend.key.height = unit(4.0, "mm"),
        legend.key.width = unit(5.0, "mm"),
        panel.grid = element_blank(),
        plot.margin = margin(6, 8, 6, 6)
      )
    
    # Fresh-natural palettes selected for separate single-panel figures
    # A: coastal mint -> blue; B: warm coral; C: botanical evidence palette
    pal_cellcycle <- c(
      "#EEF8F2",
      "#C7EBDD",
      "#86D6C9",
      "#4CB9C7",
      "#3E82C4",
      "#2C3E8F"
    )
    
    pal_neutrophil <- c(
      "#FFF1E6",
      "#FFD0B8",
      "#FFA37F",
      "#FF6F59",
      "#EF3E36",
      "#B91D3A"
    )
    
    evidence_cols <- c(
      "Retains overlap source"  = "#6EC6FF",
      "Same lineage, new state" = "#A8E6B5",
      "New cell state"          = "#FFD166",
      "Per-gene module"         = "#FDBB84",
      "Reference not formal"    = "#D6CDC2",
      "Tier 1"                  = "#8ED9A8",
      "Tier 2"                  = "#F9D77E",
      "Weak source"             = "#DADDE2",
      "Low-detected"            = "#EEF0F4",
      "Ambiguous"               = "#CDB4DB"
    )
    
    # --------------------------------------------------------------------------
    # 7.2 Prepare profile heatmap data
    # --------------------------------------------------------------------------
    
    profile_fig_df$state_s <- factor(
      profile_fig_df$fine_state_short,
      levels = short_state(fine_states)
    )
    
    profile_fig_df$set_y <- factor(
      profile_fig_df$set,
      levels = rev(c("full", "overlap", "residual")),
      labels = rev(c("Full", "Overlap", "Residual"))
    )
    
    top_mark <- do.call(
      rbind,
      lapply(
        split(profile_fig_df, paste(profile_fig_df$module, profile_fig_df$set)),
        function(df) {
          df2 <- df[is.finite(df$contribution), , drop = FALSE]
          if (!nrow(df2)) return(NULL)
          mx <- max(df2$contribution, na.rm = TRUE)
          df2[df2$contribution == mx, , drop = FALSE][1, , drop = FALSE]
        }
      )
    )
    
    if (!is.null(top_mark) && nrow(top_mark)) {
      top_mark$state_s <- factor(
        top_mark$fine_state_short,
        levels = short_state(fine_states)
      )
      top_mark$set_y <- factor(
        top_mark$set,
        levels = rev(c("full", "overlap", "residual")),
        labels = rev(c("Full", "Overlap", "Residual"))
      )
    }
    
    make_profile_heatmap <- function(module_name, fill_cols, panel_title) {
      df <- profile_fig_df[profile_fig_df$module == module_name, , drop = FALSE]
      tm <- top_mark[top_mark$module == module_name, , drop = FALSE]
      
      upper <- max(df$contribution[is.finite(df$contribution)], na.rm = TRUE)
      if (!is.finite(upper) || upper <= 0) upper <- 1
      
      ggplot(df, aes(x = state_s, y = set_y, fill = contribution)) +
        geom_tile(
          width = 0.96,
          height = 0.80,
          colour = "white",
          size = 0.28
        ) +
        geom_point(
          data = tm,
          aes(x = state_s, y = set_y),
          inherit.aes = FALSE,
          shape = 21,
          size = 2.2,
          stroke = 0.50,
          colour = "#111111",
          fill = "white"
        ) +
        scale_fill_gradientn(
          colours = fill_cols,
          limits = c(0, upper),
          trans = "sqrt",
          name = "Median\ncontribution",
          na.value = "#EEF2F4"
        ) +
        labs(
          title = panel_title,
          subtitle = "Rows = full / overlap / residual; white dot = top contributing state",
          x = NULL,
          y = NULL
        ) +
        nature_theme +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1,
            size = 6.4,
            colour = "#222222"
          ),
          axis.text.y = element_text(
            size = 7.8,
            face = "bold",
            colour = "#222222"
          ),
          legend.position = "bottom",
          legend.justification = "left"
        ) +
        guides(
          fill = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5,
            barwidth = unit(36, "mm"),
            barheight = unit(3.5, "mm"),
            ticks = TRUE
          )
        )
    }
    
    pA <- make_profile_heatmap(
      module_name = "Cell cycle",
      fill_cols = pal_cellcycle,
      panel_title = "Cell-cycle module: full / overlap / residual"
    )
    
    pB <- make_profile_heatmap(
      module_name = "Neutrophil",
      fill_cols = pal_neutrophil,
      panel_title = "Neutrophil module: full / overlap / residual"
    )
    
    # --------------------------------------------------------------------------
    # 7.3 Residual-gene evidence panel
    #    改成长条为固定大小彩色方块，避免“彩色横线”
    # --------------------------------------------------------------------------
    
    evidence_long$gene_f <- factor(
      evidence_long$gene,
      levels = rev(unique(evidence_long$gene))
    )
    
    evidence_long$panel_f <- factor(
      as.character(evidence_long$panel),
      levels = c("Novelty", "Tier")
    )
    
    evidence_long$module_f <- factor(
      evidence_long$module,
      levels = c("Cell cycle", "Inflammatory", "Neutrophil", "Other")
    )
    
    pC <- ggplot(evidence_long, aes(x = panel_f, y = gene_f)) +
      geom_point(
        aes(fill = value_label),
        shape = 22,
        size = 4.7,
        stroke = 0.35,
        colour = "white"
      ) +
      facet_grid(
        module_f ~ .,
        scales = "free_y",
        space = "free_y",
        switch = "y"
      ) +
      scale_fill_manual(
        values = evidence_cols,
        name = NULL,
        drop = FALSE
      ) +
      scale_x_discrete(
        expand = expansion(add = c(0.55, 0.55))
      ) +
      labs(
        title = "Residual-gene evidence",
        subtitle = "Single-cell evidence only",
        x = NULL,
        y = NULL
      ) +
      nature_theme +
      theme(
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(
          angle = 0,
          face = "bold",
          size = 7.8,
          colour = "#222222",
          margin = margin(r = 6)
        ),
        axis.text.y = element_text(
          size = 7.3,
          colour = "#222222"
        ),
        axis.text.x = element_text(
          size = 7.8,
          face = "bold",
          colour = "#222222"
        ),
        panel.spacing.y = unit(6.2, "mm"),
        legend.position = "right",
        legend.justification = "center",
        legend.key.size = unit(4.4, "mm"),
        legend.text = element_text(size = 6.9),
        plot.margin = margin(8, 10, 8, 6)
      ) +
      guides(
        fill = guide_legend(
          ncol = 1,
          byrow = TRUE,
          override.aes = list(
            shape = 22,
            size = 5.2,
            colour = "white"
          )
        )
      )
    
    # --------------------------------------------------------------------------
    # 7.4 Export independent single-panel figures only
    # --------------------------------------------------------------------------
    # No patchwork assembly here.
    # Final multi-panel layout will be done manually in Adobe/Illustrator.
    
    save_single_plot <- function(plot, filename, width_mm, height_mm, dpi = 600) {
      width_in  <- width_mm / 25.4
      height_in <- height_mm / 25.4
      
      ggsave(
        filename = file.path(FIG_DIR, paste0(filename, ".pdf")),
        plot = plot,
        width = width_in,
        height = height_in,
        device = cairo_pdf
      )
      
      ggsave(
        filename = file.path(FIG_DIR, paste0(filename, ".png")),
        plot = plot,
        width = width_in,
        height = height_in,
        dpi = 450
      )
      
      if (requireNamespace("svglite", quietly = TRUE)) {
        ggsave(
          filename = file.path(FIG_DIR, paste0(filename, ".svg")),
          plot = plot,
          width = width_in,
          height = height_in,
          device = svglite::svglite
        )
      }
      
      if (requireNamespace("ragg", quietly = TRUE)) {
        ggsave(
          filename = file.path(FIG_DIR, paste0(filename, ".tiff")),
          plot = plot,
          width = width_in,
          height = height_in,
          dpi = dpi,
          device = ragg::agg_tiff,
          compression = "lzw"
        )
      }
    }
    
    save_single_plot(
      pA,
      "fig04a_cellcycle_full_overlap_residual_profile",
      width_mm = 190,
      height_mm = 82
    )
    
    save_single_plot(
      pB,
      "fig04b_neutrophil_full_overlap_residual_profile",
      width_mm = 190,
      height_mm = 82
    )
    
    save_single_plot(
      pC,
      "fig04c_residual_gene_evidence",
      width_mm = 125,
      height_mm = 150
    )
    
    message("Saved independent single-panel figures:")
    message("  - fig04a_cellcycle_full_overlap_residual_profile")
    message("  - fig04b_neutrophil_full_overlap_residual_profile")
    message("  - fig04c_residual_gene_evidence")
    
  } else {
    warning("Skipping figure: ggplot2 / grid not available.")
  }
  
  # ---------------------------------------------------------------------------
  # 8. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  has_svglite <- requireNamespace("svglite", quietly=TRUE)
  has_ragg    <- requireNamespace("ragg", quietly=TRUE)
  write_csv_atomic(data.frame(
    key=c("script","script_version","s01_table_dir","s02_table_dir",
          "frozen_file","state_file","eval_file",
          "pgs_file","lcpm_file","long_file","overlap_tag","residual_tag",
          "set_modules","min_set_formal","thresh_eligible_donor_n","n_donors",
          "has_svglite","has_ragg","candidate_tier_scope","R_version"),
    value=c("04_overlap_residual_source_retention.R","v2.7",
            S01_TABLE_DIR,S02_TABLE_DIR,
            FROZEN_FILE,STATE_FILE,EVAL_FILE,PGS_FILE,LCPM_FILE,LONG_FILE,
            OVERLAP_TAG,RESIDUAL_TAG,paste(MODULE_FULL,collapse=";"),
            as.character(MIN_SET_FORMAL),as.character(THRESH_ELIGIBLE_DONOR_N),
            as.character(length(donors)),
            as.character(has_svglite),as.character(has_ragg),
            "single-cell source attribution only; bulk and condition-direction evidence not included",
            as.character(getRversion())),
    stringsAsFactors=FALSE), file.path(TABLE_DIR,"run_provenance.csv"))
  
  cat("Script 04 v2.7 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("S01 table dir:", S01_TABLE_DIR, "\n")
  cat("S02 table dir:", S02_TABLE_DIR, "\n")
  cat("\n--- Table 0b: set sizes ---\n")
  print(set_size_table, row.names=FALSE)
  cat("\n--- Table 1: module-set source summary ---\n")
  print(t1, row.names=FALSE)
  cat("\n--- Table 2 (all 32): candidate tier distribution ---\n")
  print(table(t2_all$candidate_tier))
  cat("\n--- Table 2 (residual 17): novelty / tier ---\n")
  print(t2_res[, c("gene","module","novelty_class","candidate_tier",
                   "main_contributor_state","main_state_donor_detection_rate",
                   "main_state_eligible_donor_n","overlap_set_main_state",
                   "overlap_reference_scope")], row.names=FALSE)
  cat("\nNote: Table 2 is SINGLE-CELL candidate prioritization only; bulk and\n")
  cat("condition-direction evidence are NOT included (see candidate_evidence_scope).\n")
  cat("Novelty calls require a FORMAL overlap reference; otherwise genes are\n")
  cat("flagged overlap_reference_descriptive and capped at Tier2.\n")
  cat("Tier1/novelty require >= ", THRESH_ELIGIBLE_DONOR_N, " eligible donors.\n", sep="")
})




# ============================================================================
# SCRIPT 04b v1.4 -- Candidate-level condition-context direction
#   Purpose: for Step-4 candidate genes (Tier1/2), test whether expression in
#   the gene's main contributor state changes in the bulk-expected direction
#   (Up/Down) from acute sepsis to convalescence, with HV/CS as descriptive
#   secondary. Source-state abundance change is reported as a SEPARATE layer
#   and is NOT filtered by the >=20-cell expression eligibility.
#   Full cohort (48 samples) is re-aggregated from the RDS.
#   Per-paired-donor raw and signed deltas are exported (transparent audit).
#   v1.4: symmetric descriptive expression-direction classification at
#   +/-0.3, with >=4 paired donors and positive-pair fractions >=2/3 or
#   <=1/3; added +/-0.5 sensitivity. Cell-abundance rules unchanged.
#   Outputs:
#     table04b_candidate_condition_direction.csv
#     table04b_paired_expression_delta_long.csv
#     table04b_paired_cell_fraction_delta_long.csv
#     figures/fig04b_1_condition_logcpm_dotmap.(pdf/png)
#     figures/fig04b_2_paired_acute_conv_slope.(pdf/png)
#     figures/fig04b_3_main_state_cell_fraction_paired.(pdf/png)
# RUN: Rscript 04b_candidate_condition_direction.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG + PACKAGES
  # ---------------------------------------------------------------------------
  PROJECT_ROOT <- "/home/sunshine/predicate/singlecell"
  INPUT_RDS    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  OUTPUTS_ROOT <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  SCRIPT04_ROOT <- file.path(OUTPUTS_ROOT, "04_overlap_residual_source_retention")
  
  STAGE_NAME  <- "04b_candidate_condition_direction"
  RUN_PURPOSE <- "candidate_condition_direction"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  
  MIN_CELLS_PER_SAMPLE_STATE <- 20L
  # Descriptive expression-direction classification.
  CONDITION_PAIRED_MIN_N <- 4L
  CONDITION_DELTA <- 0.3
  CONDITION_ALIGNED_FRAC_MIN <- 2 / 3
  CONDITION_OPPOSITE_FRAC_MAX <- 1 / 3
  CONDITION_SENSITIVITY_DELTA <- 0.5
  
  EXPECTED_N_SAMPLES <- 48L
  EXPECTED_N_ACUTE   <- 26L
  EXPECTED_N_CONV    <- 9L
  EXPECTED_N_HV      <- 6L
  EXPECTED_N_CS      <- 7L
  EXPECTED_N_PAIRED  <- 9L
  
  CANDIDATE_TIERS <- c("Tier1_priority", "Tier2_supportive")
  
  needed <- c("SeuratObject", "Matrix")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  require_cols <- function(x, cols, label) {
    miss <- setdiff(cols, names(x))
    if (length(miss)) stop_msg(label, " missing required columns: ", paste(miss, collapse = ", "))
    invisible(TRUE)
  }
  one_chr <- function(x) {
    if (length(x) < 1L || is.na(x[1])) return(NA_character_)
    as.character(x[1])
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
    y[y %in% c("Other")] <- "Other"
    y
  }
  short_state <- function(x) {
    z <- as.character(x)
    map <- c(
      "HSPCs"="HSPCs",
      "Cycling_neutrophil_progenitors"="Cycling neut prog",
      "MPO+_immature_neutrophils_or_progenitors"="MPO+ imm neut/prog",
      "PADI4+_immature_neutrophils"="PADI4+ imm neut",
      "IL1R2+_immature_neutrophils"="IL1R2+ imm neut",
      "S100A8-9_hi_neutrophils"="S100A8/9 hi neut",
      "Mature_neutrophils"="Mature neut",
      "Degranulating_neutrophils"="Degran neut",
      "Apoptosing_neutrophils"="Apoptosing neut",
      "Eosinophils"="Eosinophils","Mast_cells/eosiniophils"="Mast/eos",
      "Classical_monocytes"="Classical mono","Non-classical_monocytes"="Non-classical mono",
      "cDCs"="cDCs","pDCs"="pDCs",
      "Naive_CD4_T_cells"="Naive CD4 T","Memory_CD4_T_cells"="Memory CD4 T",
      "Naive_CD8_T_cells"="Naive CD8 T","CD8_T_cells"="CD8 T",
      "Cycling_TNK"="Cycling TNK","NK cells"="NK cells","B_cells"="B cells",
      "Plasmablasts"="Plasmablasts","Platelets"="Platelets")
    ifelse(z %in% names(map), unname(map[z]), z)
  }
  
  call_condition_direction <- function(
    paired_n, signed_delta_median, support_rate,
    delta_cut = CONDITION_DELTA) {
    
    if (!is.finite(paired_n) ||
        paired_n < CONDITION_PAIRED_MIN_N ||
        !is.finite(signed_delta_median) ||
        !is.finite(support_rate)) {
      return("not_evaluable")
    }
    
    if (signed_delta_median >= delta_cut &&
        support_rate >= CONDITION_ALIGNED_FRAC_MIN) {
      return("aligned")
    }
    
    if (signed_delta_median <= -delta_cut &&
        support_rate <= CONDITION_OPPOSITE_FRAC_MAX) {
      return("opposite_direction")
    }
    
    "neither_criterion_met"
  }
  
  call_state_abundance <- function(paired_n, cell_fraction_delta_median, support_rate) {
    if (!is.finite(paired_n) || paired_n < 5) return("not_evaluable")
    if (is.finite(cell_fraction_delta_median) && cell_fraction_delta_median > 0 &&
        is.finite(support_rate) && support_rate >= 2/3) return("acute_expanded")
    if (is.finite(cell_fraction_delta_median) && cell_fraction_delta_median < 0 &&
        is.finite(support_rate) && support_rate <= 1/3) return("acute_depleted")
    "not_consistent"
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. CANDIDATE INPUT (from Script 04) + GENE MAPPING
  # ---------------------------------------------------------------------------
  CAND_FILE <- find_latest_file(SCRIPT04_ROOT, "^table2_all_gene_candidate\\.csv$")
  STATE_FILE <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  MAP_FILE   <- find_latest_file(SCRIPT01_ROOT, "^gene_mapping_audit\\.csv$")
  
  cands <- read.csv(CAND_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  state01 <- read.csv(STATE_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  map01 <- read.csv(MAP_FILE, check.names=FALSE, stringsAsFactors=FALSE)
  
  require_cols(cands, c("gene","module","direction","candidate_tier","novelty_class",
                        "main_contributor_state","main_contributor_lineage","E3_membership"),
               "table2_all_gene_candidate.csv")
  require_cols(state01, c("display_order","fine_state","lineage"), "fine_state_display_order_audit.csv")
  require_cols(map01, c("signature_gene","dataset_feature"), "gene_mapping_audit.csv")
  
  if (!all(cands$direction %in% c("Up","Down"))) {
    stop_msg("Unexpected direction values in candidate table: ",
             paste(unique(cands$direction[!cands$direction %in% c("Up","Down")]), collapse=", "))
  }
  
  cands$candidate_tier_raw <- cands$candidate_tier
  cands$candidate_tier <- ifelse(
    cands$candidate_tier %in% c("Tier 1", "Tier1", "Tier1_priority"),
    "Tier1_priority",
    ifelse(
      cands$candidate_tier %in% c("Tier 2", "Tier2", "Tier2_supportive"),
      "Tier2_supportive",
      cands$candidate_tier
    )
  )
  
  cands <- cands[cands$candidate_tier %in% CANDIDATE_TIERS, , drop=FALSE]
  if (!nrow(cands)) stop_msg("No candidate genes after tier filter: ", paste(CANDIDATE_TIERS, collapse="; "))
  if (anyDuplicated(cands$gene)) stop_msg("Duplicated candidate genes")
  
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L)
  if (any(!cands$main_contributor_state %in% fine_states)) {
    stop_msg("Candidate main_contributor_state not in fine_states: ",
             paste(unique(cands$main_contributor_state[!cands$main_contributor_state %in% fine_states]), collapse=",
  "))
  }
  
  feat_map <- setNames(map01$dataset_feature, map01$signature_gene)
  cand_feat <- feat_map[cands$gene]
  if (anyNA(cand_feat) || any(!nzchar(cand_feat))) stop_msg("Missing dataset_feature for some candidates")
  cands$dataset_feature <- cand_feat
  
  # ---------------------------------------------------------------------------
  # 3. READ RDS + FULL COHORT (all 48 samples, 4 conditions)
  # ---------------------------------------------------------------------------
  if (!file.exists(INPUT_RDS)) stop_msg("Input RDS not found: ", INPUT_RDS)
  obj <- readRDS(INPUT_RDS)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  require_cols(meta, c(SAMPLE_FIELD, CONDITION_FIELD, FINE_STATE_FIELD), "obj@meta.data")
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay=COUNT_ASSAY, layer=COUNT_LAYER, fast=FALSE)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    if ("layer" %in% names(formals(GetAssayData))) cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, layer="counts")
    else cnt <- GetAssayData(object=obj, assay=COUNT_ASSAY, slot="counts")
  }
  if (!identical(colnames(cnt), colnames(obj))) {
    stop_msg("Count matrix column names do not match object cells")
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  stopifnot(!anyNA(cell_sample), !anyNA(cell_diag), !anyNA(cell_state))
  
  samples <- sort(unique(cell_sample))
  sample_diag_list <- lapply(samples, function(s) unique(cell_diag[cell_sample == s]))
  bad_diag_samples <- samples[lengths(sample_diag_list) != 1L]
  if (length(bad_diag_samples)) {
    stop_msg("Samples with non-unique diagnosis: ", paste(bad_diag_samples, collapse=", "))
  }
  sample_diag <- vapply(sample_diag_list, `[`, character(1), 1L)
  
  sample_tab <- data.frame(sample_id=samples, diagnosis=unname(sample_diag), stringsAsFactors=FALSE)
  sample_tab$donor <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort <- NA_character_
  sample_tab$cohort[sample_tab$diagnosis %in% ACUTE_DIAG] <- "Acute"
  sample_tab$cohort[sample_tab$diagnosis == "Conv"] <- "Conv"
  sample_tab$cohort[sample_tab$diagnosis == "HV"] <- "HV"
  sample_tab$cohort[sample_tab$diagnosis == "CS"] <- "CS"
  if (anyNA(sample_tab$cohort)) stop_msg("Unmapped diagnosis")
  if (length(samples) != EXPECTED_N_SAMPLES) stop_msg("Expected 48 samples, found ", length(samples))
  
  stopifnot(sum(sample_tab$cohort=="Acute") == EXPECTED_N_ACUTE)
  stopifnot(sum(sample_tab$cohort=="Conv")  == EXPECTED_N_CONV)
  stopifnot(sum(sample_tab$cohort=="HV")    == EXPECTED_N_HV)
  stopifnot(sum(sample_tab$cohort=="CS")    == EXPECTED_N_CS)
  
  is_conv_suffix <- grepl("_CONV$", sample_tab$sample_id)
  is_conv_diag   <- sample_tab$diagnosis == "Conv"
  if (!identical(is_conv_suffix, is_conv_diag)) {
    stop_msg("_CONV suffix is not identical to diagnosis == 'Conv'")
  }
  
  acute_donors <- sample_tab$donor[sample_tab$cohort == "Acute"]
  conv_donors  <- sample_tab$donor[sample_tab$cohort == "Conv"]
  if (!all(conv_donors %in% acute_donors)) {
    stop_msg("Some convalescent donors are not present in acute donors: ",
             paste(setdiff(conv_donors, acute_donors), collapse=", "))
  }
  paired_donors <- intersect(acute_donors, conv_donors)
  stopifnot(length(paired_donors) == EXPECTED_N_PAIRED)
  
  # ---------------------------------------------------------------------------
  # 4. SAMPLE-STATE AGGREGATION (all cells)
  # ---------------------------------------------------------------------------
  n_samples <- length(samples); n_states <- length(fine_states)
  si_idx <- match(cell_sample, samples); st_idx <- match(cell_state, fine_states)
  stopifnot(!anyNA(si_idx), !anyNA(st_idx))
  ds_all <- (si_idx - 1L) * n_states + st_idx
  D <- Matrix::sparseMatrix(i=ds_all, j=seq_along(ds_all), x=1,
                            dims=c(n_samples*n_states, length(ds_all)))
  
  n_cells_ds <- as.integer(Matrix::rowSums(D))
  lib_cell <- as.numeric(Matrix::colSums(cnt))
  lib_ds <- as.numeric(D %*% lib_cell)
  eligible_ds <- n_cells_ds >= MIN_CELLS_PER_SAMPLE_STATE & lib_ds > 0
  
  n_cells_mat <- matrix(n_cells_ds, nrow=n_samples, ncol=n_states, byrow=TRUE,
                        dimnames=list(samples, fine_states))
  sample_total <- rowSums(n_cells_mat)
  cf_mat <- sweep(n_cells_mat, 1, sample_total, "/")
  
  gidx <- match(cand_feat, rownames(cnt))
  if (anyNA(gidx)) stop_msg("Candidate feature(s) absent from count matrix")
  X_c <- cnt[gidx, , drop=FALSE]; rownames(X_c) <- cands$gene
  C_ds <- as.matrix(X_c %*% Matrix::t(D))            # n_cand x (n_samples*n_states)
  stopifnot(ncol(C_ds) == n_samples*n_states)
  
  lc_ds <- log2(sweep(C_ds, 2, lib_ds, "/") * 1e6 + 1)
  lc_ds[, !eligible_ds] <- NA_real_
  
  acute_samp <- setNames(sample_tab$sample_id[sample_tab$cohort=="Acute"],
                         sample_tab$donor[sample_tab$cohort=="Acute"])
  conv_samp  <- setNames(sample_tab$sample_id[sample_tab$cohort=="Conv"],
                         sample_tab$donor[sample_tab$cohort=="Conv"])
  
  # ---------------------------------------------------------------------------
  # 5. PER-CANDIDATE CONDITION DIRECTION
  # ---------------------------------------------------------------------------
  rows <- list()
  expr_long_list <- list()   # per-paired-donor expression deltas (transparent audit)
  cf_long_list <- list()     # per-paired-donor cell-fraction deltas
  for (i in seq_len(nrow(cands))) {
    g <- cands$gene[i]
    ms <- cands$main_contributor_state[i]
    dir <- cands$direction[i]
    exp_sign <- if (dir == "Up") 1 else -1
    msi <- match(ms, fine_states)
    col_main <- (seq_len(n_samples) - 1L) * n_states + msi
    lg_g <- lc_ds[i, col_main]
    elig_s <- eligible_ds[col_main]
    cf_g <- cf_mat[, msi]
    n_cells_ms <- n_cells_mat[, msi]
    grp <- sample_tab$cohort
    
    acute_vals <- lg_g[grp=="Acute" & elig_s]
    conv_vals  <- lg_g[grp=="Conv"  & elig_s]
    hv_vals    <- lg_g[grp=="HV"    & elig_s]
    cs_vals    <- lg_g[grp=="CS"    & elig_s]
    
    delta_raw_expr <- c()
    cfrac_delta_all <- c()
    for (d in paired_donors) {
      a_s <- acute_samp[[d]]; c_s <- conv_samp[[d]]
      ai <- match(a_s, samples); ci <- match(c_s, samples)
      if (is.finite(cf_g[ai]) && is.finite(cf_g[ci])) {
        cfrac_delta_all <- c(cfrac_delta_all, cf_g[ai] - cf_g[ci])
        cf_long_list[[length(cf_long_list)+1L]] <- data.frame(
          gene = g, module = cands$module[i], direction = dir, donor = d,
          acute_cell_fraction = cf_g[ai], conv_cell_fraction = cf_g[ci],
          paired_cell_fraction_delta = cf_g[ai] - cf_g[ci],
          acute_main_state_cells = n_cells_ms[ai], conv_main_state_cells = n_cells_ms[ci],
          stringsAsFactors=FALSE)
      }
      if (elig_s[ai] && elig_s[ci]) {
        delta_raw_expr <- c(delta_raw_expr, lg_g[ai] - lg_g[ci])
        expr_long_list[[length(expr_long_list)+1L]] <- data.frame(
          gene = g, module = cands$module[i], direction = dir,
          expected_sign = exp_sign, donor = d,
          acute_logCPM = lg_g[ai], conv_logCPM = lg_g[ci],
          paired_delta_raw = lg_g[ai] - lg_g[ci],
          paired_delta_signed = exp_sign * (lg_g[ai] - lg_g[ci]),
          stringsAsFactors=FALSE)
      }
    }
    paired_n_expr <- length(delta_raw_expr)
    paired_n_abun <- length(cfrac_delta_all)
    signed_delta <- exp_sign * delta_raw_expr
    
    paired_delta_raw_med <- if (paired_n_expr) median(delta_raw_expr) else NA_real_
    paired_delta_raw_iqr <- if (paired_n_expr) IQR(delta_raw_expr) else NA_real_
    paired_delta_sig_med <- if (paired_n_expr) median(signed_delta) else NA_real_
    paired_delta_sig_iqr <- if (paired_n_expr) IQR(signed_delta) else NA_real_
    support_rate         <- if (paired_n_expr) mean(signed_delta > 0) else NA_real_
    wilcox_p_signed <- if (paired_n_expr >= 3) {
      tryCatch(wilcox.test(signed_delta, mu=0, exact=FALSE)$p.value, error=function(e) NA_real_)
    } else NA_real_
    
    cf_acute_med <- median(cf_g[grp=="Acute"])
    cf_conv_med  <- median(cf_g[grp=="Conv"])
    cf_hv_med    <- median(cf_g[grp=="HV"])
    cf_cs_med    <- median(cf_g[grp=="CS"])
    cf_delta_med <- if (paired_n_abun) median(cfrac_delta_all) else NA_real_
    cf_support_rate <- if (paired_n_abun) mean(cfrac_delta_all > 0) else NA_real_
    
    dir_call <- call_condition_direction(
      paired_n = paired_n_expr,
      signed_delta_median = paired_delta_sig_med,
      support_rate = support_rate,
      delta_cut = CONDITION_DELTA
    )
    
    dir_call_abs0p5 <- call_condition_direction(
      paired_n = paired_n_expr,
      signed_delta_median = paired_delta_sig_med,
      support_rate = support_rate,
      delta_cut = CONDITION_SENSITIVITY_DELTA
    )
    
    abun_call <- call_state_abundance(
      paired_n_abun, cf_delta_med, cf_support_rate
    )
    
    rows[[length(rows)+1L]] <- data.frame(
      gene = g,
      module = cands$module[i],
      direction = dir,
      expected_sign = exp_sign,
      overlap_status = cands$E3_membership[i],
      candidate_tier = cands$candidate_tier[i],
      novelty_class = cands$novelty_class[i],
      main_contributor_state = ms,
      main_lineage = cands$main_contributor_lineage[i],
      acute_n_evaluable = length(acute_vals),
      conv_n_evaluable = length(conv_vals),
      paired_n_evaluable = paired_n_expr,
      paired_cell_fraction_paired_n = paired_n_abun,
      HV_n_evaluable = length(hv_vals),
      CS_n_evaluable = length(cs_vals),
      acute_median_logCPM = median(acute_vals, na.rm=TRUE),
      conv_median_logCPM = median(conv_vals, na.rm=TRUE),
      HV_median_logCPM = median(hv_vals, na.rm=TRUE),
      CS_median_logCPM = median(cs_vals, na.rm=TRUE),
      paired_delta_raw_median = paired_delta_raw_med,
      paired_delta_raw_IQR = paired_delta_raw_iqr,
      paired_delta_signed_median = paired_delta_sig_med,
      paired_delta_signed_IQR = paired_delta_sig_iqr,
      paired_sign_support_rate = support_rate,
      paired_wilcox_p_signed_reference_only = wilcox_p_signed,
      main_state_cell_fraction_acute_median = cf_acute_med,
      main_state_cell_fraction_conv_median = cf_conv_med,
      main_state_cell_fraction_HV_median = cf_hv_med,
      main_state_cell_fraction_CS_median = cf_cs_med,
      paired_cell_fraction_delta_median = cf_delta_med,
      paired_cell_fraction_support_rate = cf_support_rate,
      condition_direction_call = dir_call,
      state_abundance_call = abun_call,
      condition_support_final = dir_call,
      condition_direction_call_abs0p5 = dir_call_abs0p5,
      stringsAsFactors=FALSE)
  }
  t <- do.call(rbind, rows)
  write_csv_atomic(t, file.path(TABLE_DIR, "table04b_candidate_condition_direction.csv"))
  
  expr_long <- if (length(expr_long_list)) do.call(rbind, expr_long_list) else data.frame()
  cf_long   <- if (length(cf_long_list))   do.call(rbind, cf_long_list)   else data.frame()
  write_csv_atomic(expr_long, file.path(TABLE_DIR, "table04b_paired_expression_delta_long.csv"))
  write_csv_atomic(cf_long,   file.path(TABLE_DIR, "table04b_paired_cell_fraction_delta_long.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. FIGURES (empty-data guarded)
  # ---------------------------------------------------------------------------
  if (requireNamespace("ggplot2", quietly=TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    
    cond_prefix <- c(Acute = "acute", Conv = "conv", HV = "HV", CS = "CS")
    
    cond_long <- do.call(rbind, lapply(c("Acute","Conv","HV","CS"), function(cond) {
      px <- cond_prefix[[cond]]
      med_col <- paste0(px, "_median_logCPM")
      n_col   <- paste0(px, "_n_evaluable")
      data.frame(gene=t$gene, module=t$module, direction=t$direction,
                 condition=cond,
                 logcpm=t[[med_col]],
                 n_eval=t[[n_col]],
                 stringsAsFactors=FALSE)
    }))
    cond_long$gene_f <- factor(cond_long$gene, levels=rev(t$gene))
    cond_long$condition <- factor(cond_long$condition, levels=c("Acute","Conv","HV","CS"))
    p1 <- ggplot(cond_long, aes(x=condition, y=gene_f, fill=logcpm)) +
      geom_tile(colour="white", width=0.85, height=0.85) +
      scale_fill_gradientn(colours=c("#E8EEF4","#9DBED1","#4E88A8","#2B5D7A"),
                           na.value="#EDEDED", name="median logCPM") +
      facet_grid(module ~ ., scales="free_y", space="free_y") +
      theme_classic(base_size=8) +
      theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
            axis.text.y=element_text(size=6.5),
            panel.spacing=grid::unit(4,"mm")) +
      labs(x=NULL, y=NULL, title="Candidate gene logCPM in main source state by condition",
           caption="Grey = not evaluable (<20 cells in state)")
    ggsave(file.path(FIG_DIR,"fig04b_1_condition_logcpm_dotmap.pdf"), p1, width=6.5, height=9, device=cairo_pdf)
    ggsave(file.path(FIG_DIR,"fig04b_1_condition_logcpm_dotmap.png"), p1, width=6.5, height=9, dpi=300)
    
    slope_list <- list()
    for (i in seq_len(nrow(cands))) {
      g <- cands$gene[i]; ms <- cands$main_contributor_state[i]
      msi <- match(ms, fine_states)
      col_main <- (seq_len(n_samples)-1L)*n_states + msi
      lg_g <- lc_ds[i, col_main]; elig_s <- eligible_ds[col_main]
      for (d in paired_donors) {
        a_s <- acute_samp[[d]]; c_s <- conv_samp[[d]]
        ai <- match(a_s, samples); ci <- match(c_s, samples)
        if (elig_s[ai] && elig_s[ci]) {
          slope_list[[length(slope_list)+1L]] <- data.frame(
            gene=g, module=cands$module[i], donor=d,
            timepoint=c("Acute","Conv"), logcpm=c(lg_g[ai], lg_g[ci]),
            stringsAsFactors=FALSE)
        }
      }
    }
    if (!length(slope_list)) {
      warning("No paired evaluable expression rows; skipping fig04b_2.")
    } else {
      slope_df <- do.call(rbind, slope_list)
      slope_df$timepoint <- factor(slope_df$timepoint, levels=c("Acute","Conv"))
      slope_df$gene_f <- factor(slope_df$gene, levels=rev(t$gene))
      p2 <- ggplot(slope_df, aes(x=timepoint, y=logcpm, group=donor)) +
        geom_line(colour="#9DBED1", linewidth=0.35) +
        geom_point(size=1.0, colour="#2B5D7A") +
        facet_wrap(~ gene_f, scales="free_y", ncol=4) +
        theme_classic(base_size=7) +
        theme(strip.text=element_text(size=6.5),
              axis.text.x=element_text(size=6.5),
              panel.spacing=grid::unit(2,"mm")) +
        labs(x=NULL, y="logCPM in main state",
             title="Paired acute-conv expression (main source state)",
             caption="Up genes: decline = support; Down genes: rise = support")
      ggsave(file.path(FIG_DIR,"fig04b_2_paired_acute_conv_slope.pdf"), p2, width=11, height=7, device=cairo_pdf)
      ggsave(file.path(FIG_DIR,"fig04b_2_paired_acute_conv_slope.png"), p2, width=11, height=7, dpi=300)
    }
    
    cf_list <- list()
    for (i in seq_len(nrow(cands))) {
      g <- cands$gene[i]; ms <- cands$main_contributor_state[i]
      msi <- match(ms, fine_states)
      cf_g <- cf_mat[, msi]
      n_cells_ms <- n_cells_mat[, msi]
      for (d in paired_donors) {
        a_s <- acute_samp[[d]]; c_s <- conv_samp[[d]]
        ai <- match(a_s, samples); ci <- match(c_s, samples)
        cf_list[[length(cf_list)+1L]] <- data.frame(
          gene=g, module=cands$module[i], donor=d,
          timepoint=c("Acute","Conv"), cell_fraction=c(cf_g[ai], cf_g[ci]),
          acute_main_state_cells=n_cells_ms[ai], conv_main_state_cells=n_cells_ms[ci],
          stringsAsFactors=FALSE)
      }
    }
    if (!length(cf_list)) {
      warning("No paired cell-fraction rows; skipping fig04b_3.")
    } else {
      cf_df <- do.call(rbind, cf_list)
      cf_df$timepoint <- factor(cf_df$timepoint, levels=c("Acute","Conv"))
      cf_df$gene_f <- factor(cf_df$gene, levels=rev(t$gene))
      p3 <- ggplot(cf_df, aes(x=timepoint, y=cell_fraction, group=donor)) +
        geom_line(colour="#86A873", linewidth=0.35) +
        geom_point(size=1.0, colour="#3F5C4A") +
        facet_wrap(~ gene_f, scales="free_y", ncol=4) +
        theme_classic(base_size=7) +
        theme(strip.text=element_text(size=6.5),
              axis.text.x=element_text(size=6.5),
              panel.spacing=grid::unit(2,"mm")) +
        labs(x=NULL, y="Main-state cell fraction",
             title="Paired acute-conv source-state abundance",
             caption="Decline = immature source state waning in convalescence")
      ggsave(file.path(FIG_DIR,"fig04b_3_main_state_cell_fraction_paired.pdf"), p3, width=11, height=7,
             device=cairo_pdf)
      ggsave(file.path(FIG_DIR,"fig04b_3_main_state_cell_fraction_paired.png"), p3, width=11, height=7, dpi=300)
    }
  } else {
    warning("Skipping figures: ggplot2 not available.")
  }
  
  # ---------------------------------------------------------------------------
  # 7. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  write_csv_atomic(data.frame(
    key=c("script","script_version","input_rds","candidate_file","state_file","map_file",
          "min_cells_per_sample_state","n_samples","n_acute","n_conv","n_hv","n_cs","n_paired",
          "candidate_tiers","n_candidates",
          "condition_paired_min_n","condition_delta",
          "condition_aligned_fraction_min","condition_opposite_fraction_max",
          "condition_sensitivity_delta","R_version"),
    value=c("04b_candidate_condition_direction.R","v1.4",INPUT_RDS,CAND_FILE,STATE_FILE,MAP_FILE,
            as.character(MIN_CELLS_PER_SAMPLE_STATE),
            as.character(length(samples)),as.character(sum(sample_tab$cohort=="Acute")),
            as.character(sum(sample_tab$cohort=="Conv")),as.character(sum(sample_tab$cohort=="HV")),
            as.character(sum(sample_tab$cohort=="CS")),as.character(length(paired_donors)),
            paste(CANDIDATE_TIERS,collapse=";"),as.character(nrow(cands)),
            as.character(CONDITION_PAIRED_MIN_N),
            as.character(CONDITION_DELTA),
            as.character(CONDITION_ALIGNED_FRAC_MIN),
            as.character(CONDITION_OPPOSITE_FRAC_MAX),
            as.character(CONDITION_SENSITIVITY_DELTA),
            as.character(getRversion())),
    stringsAsFactors=FALSE), file.path(TABLE_DIR,"run_provenance.csv"))
  
  cat("Script 04b v1.4 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("Candidates analyzed:", nrow(cands), "\n")
  cat("\n--- Expression-direction classification (+/-0.3) ---\n")
  print(table(t$condition_direction_call, useNA="ifany"))
  
  cat("\n--- Expression-direction sensitivity (+/-0.5) ---\n")
  print(table(t$condition_direction_call_abs0p5, useNA="ifany"))
  cat("\n--- state_abundance_call distribution ---\n")
  print(table(t$state_abundance_call, useNA="ifany"))
  cat("\n--- Condition direction table (key columns) ---\n")
  print(t[, c("gene","module","direction","main_contributor_state",
              "paired_n_evaluable","paired_delta_signed_median","paired_sign_support_rate",
              "condition_direction_call","condition_direction_call_abs0p5",
              "state_abundance_call")], row.names=FALSE)
  cat("\nPaired donor detail rows: expression =", nrow(expr_long),
      "| cell fraction =", nrow(cf_long), "\n")
  cat("\nNote: condition_support_final = condition_direction_call.\n")
  cat("Expression layer requires >=20 cells in both paired donor-states; abundance\n")
  cat("layer uses all paired donors (state collapse is evidence). HV/CS descriptive\n")
  cat("only. P values are reference-only. Per-paired-donor raw+signed deltas are in\n")
  cat("table04b_paired_*_delta_long.csv.\n")
})