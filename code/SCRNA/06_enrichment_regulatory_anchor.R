# ============================================================================
# SCRIPT 06a v1.4 -- Motif enrichment: CEBP-family sequence-level plausibility
#   Purpose: an expression-independent layer of regulatory plausibility for the
#   granule program. Asks: are the promoters of the M2/primary-core/residual
#   neutrophil genes enriched for CEBP-family (and GFI1/SPI1/STAT3) motifs,
#   relative to a matched background? IRF8 is a monocyte-TF comparator.
#
#   Positioning: SEQUENCE-LEVEL PLAUSIBILITY ONLY. M2_full (10 genes) is the
#   main enrichment set; primary_core/residual_neutrophil (4 genes each) are
#   direction-support only.
#
#   Background: matched on expression AND promoter GC in the carrying states,
#   target-excluded. R matched backgrounds per gene set; Fisher both-sided per
#   draw with per-draw BH over features; aggregated as median OR / median p /
#   fraction of draws with significant FDR (the call criterion).
#
#   Statistical honesty notes:
#   - Draws are NOT independent (all sample from the same matched pools);
#     prop_fdr_signif measures stability of the call under background
#     resampling, NOT a replication probability.
#   - median_p is a descriptive summary, NOT an inferential p-value.
#   - or_q2.5/or_q97.5 are across-draw quantiles of the resampling
#     distribution, NOT a confidence interval (draws are non-independent).
#   - Per-draw BH includes the aggregate rows (CEBP_family_any, per-TF _any)
#     which are deterministically correlated with their constituent motifs;
#     a significant aggregate + significant constituent is ONE piece of
#     evidence, not two. Gene sets are nested the same way.
#   - frac_2d (background drawn from the full expression x GC stratum) is
#     reported so the GC-match quality is a number, not a caveat.
#
#   v1.2: both-sided Fisher; R-draw resampling; expression x GC matching; trim;
#     per-TF-any; NA assertions; sessionInfo.
#   v1.3: target_rate/background_rate columns; fail-loud universe-dropout
#     guards; is_cebp_family flags CEBP_family_any; OR 2.5-97.5%; 2D-stratum
#     audit; res-NULL guard.
#   v1.4: identical(names(bg_bin), names(gc_bin)) alignment assertion;
#     per-set dropped-target warning; OR quantiles explicitly disclaimed as
#     resampling stability intervals, not confidence intervals.
#
#   Requires (Bioconductor): TFBSTools, motifmatchr, JASPAR2022,
#     BSgenome.Hsapiens.UCSC.hg38, TxDb.Hsapiens.UCSC.hg38.knownGene,
#     org.Hs.eg.db, AnnotationDbi, GenomicRanges, Biostrings,
#     SummarizedExperiment
#   RUN: Rscript 06a_motif_enrichment_regulatory_anchor.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  SCRIPT_NAME    <- "06a_motif_enrichment_regulatory_anchor.R"
  SCRIPT_VERSION <- "v1.4"
  
  PROJECT_ROOT  <- "/home/sunshine/predicate/singlecell"
  INPUT_FILE    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  OUTPUTS_ROOT  <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  
  STAGE_NAME  <- "06a_motif_enrichment"
  RUN_PURPOSE <- "cebp_family_motif_enrichment"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  
  EXPECTED_N_ACUTE_DONORS  <- 26L
  EXPECTED_N_PRIMARY_CELLS <- 151837L
  
  GENE_SETS <- list(
    primary_core       = c("MPO","ELANE","CTSG","DEFA4"),
    M2_full            = c("CEACAM6","CEACAM8","CTSG","DEFA4","ELANE","MPO",
                           "MS4A3","OLFM4","RNASE3","TCN1"),
    residual_neutrophil = c("CEACAM8","OLFM4","RNASE3","TCN1")
  )
  CARRYING_STATES <- c(
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils"
  )
  
  MOTIF_NAMES <- c("CEBPA","CEBPB","CEBPE","CEBPD","GFI1","SPI1","STAT3","IRF8")
  CEBP_FAMILY <- c("CEBPA","CEBPB","CEBPE","CEBPD")
  
  PROM_UP        <- 2000L
  PROM_DOWN      <- 2000L
  BG_PER_TARGET  <- 20L
  BG_MIN         <- 60L
  EXPR_BINS      <- 20L
  GC_BINS        <- 6L
  R_DRAWS        <- 100L
  MOTIF_P_CUTOFF <- 1e-6
  FDR_THRESH     <- 0.05
  PROP_SIGNIF_CALL <- 0.5
  GC_WARN_DIFF   <- 0.02
  BG_MIN_UNIVERSE <- 1000L
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS + PACKAGES (fail-closed)
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
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
  map_symbols_to_entrez <- function(symbols) {
    symbols <- unique(as.character(symbols))
    entrez <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
      column = "ENTREZID", multiVals = "first"))
    out <- data.frame(SYMBOL = names(entrez), ENTREZID = as.character(entrez),
                      stringsAsFactors = FALSE)
    out <- out[!is.na(out$ENTREZID) & nzchar(out$ENTREZID), , drop = FALSE]
    out <- out[!duplicated(out$SYMBOL) & !duplicated(out$ENTREZID), , drop = FALSE]
    out
  }
  
  REQUIRED_PKGS <- c("TFBSTools","motifmatchr","JASPAR2022",
                     "BSgenome.Hsapiens.UCSC.hg38",
                     "TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db",
                     "AnnotationDbi","GenomicRanges","Biostrings",
                     "SummarizedExperiment","GenomicFeatures")
  missing_pkgs <- REQUIRED_PKGS[!vapply(REQUIRED_PKGS, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop_msg("Missing Bioconductor packages: ", paste(missing_pkgs, collapse = ", "), "\n",
             "Install once with:\n",
             "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
             "  BiocManager::install(c('", paste(missing_pkgs, collapse = "','"), "'))")
  }
  suppressPackageStartupMessages({
    library(SeuratObject); library(Matrix)
    library(TFBSTools); library(motifmatchr); library(JASPAR2022)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicFeatures)
  })
  GENOME <- BSgenome.Hsapiens.UCSC.hg38
  TXDB   <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. STATE ORDER
  # ---------------------------------------------------------------------------
  STATE_FILE <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  state01 <- read.csv(STATE_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L, all(CARRYING_STATES %in% fine_states))
  
  # ---------------------------------------------------------------------------
  # 3. EXPRESSION-MATCHED BACKGROUND: pooled logCPM in carrying states
  # ---------------------------------------------------------------------------
  if (!file.exists(INPUT_FILE)) stop_msg("Input RDS not found: ", INPUT_FILE)
  obj <- readRDS(INPUT_FILE)
  meta <- obj@meta.data
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay = COUNT_ASSAY, layer = COUNT_LAYER, fast = FALSE)
  } else {
    if ("layer" %in% names(formals(GetAssayData))) cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, layer =
                                                                         "counts")
    else cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, slot = "counts")
  }
  if (!inherits(cnt, "dgCMatrix")) cnt <- methods::as(cnt, "dgCMatrix")
  
  cell_sample <- as.character(meta[[SAMPLE_FIELD]])
  cell_diag   <- as.character(meta[[CONDITION_FIELD]])
  cell_state  <- as.character(meta[[FINE_STATE_FIELD]])
  sample_tab <- unique(data.frame(sample_id = cell_sample, diagnosis = cell_diag, stringsAsFactors = FALSE))
  sample_tab$donor <- sub("_CONV$", "", sample_tab$sample_id)
  acute_donors <- unique(sample_tab$donor[sample_tab$diagnosis %in% ACUTE_DIAG])
  stopifnot(length(acute_donors) == EXPECTED_N_ACUTE_DONORS)
  is_acute <- sub("_CONV$", "", cell_sample) %in% acute_donors & cell_diag %in% ACUTE_DIAG
  primary_idx <- which(is_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  
  carry_idx <- which(cell_state[primary_idx] %in% CARRYING_STATES)
  if (!length(carry_idx)) stop_msg("No carrying-state cells found in acute primary cells.")
  sub_cnt <- cnt[, primary_idx[carry_idx], drop = FALSE]
  lib_sub <- as.numeric(Matrix::colSums(sub_cnt))
  if (any(lib_sub <= 0)) stop_msg("Non-positive library size in carrying-state cells.")
  gene_umi <- Matrix::rowSums(sub_cnt)
  lib_total <- sum(lib_sub)
  expr_lc <- log2(gene_umi / lib_total * 1e6 + 1)
  names(expr_lc) <- rownames(cnt)
  
  expressed_genes <- names(expr_lc)[is.finite(expr_lc) & expr_lc > 0]
  cat("Expressed genes in carrying states:", length(expressed_genes), "\n")
  
  # ---------------------------------------------------------------------------
  # 4. GENE SETS -> TSS -> PROMOTERS; BACKGROUND UNIVERSE (symbol-aware)
  # ---------------------------------------------------------------------------
  all_target <- unique(unlist(GENE_SETS, use.names = FALSE))
  sym_map <- map_symbols_to_entrez(all_target)
  genes_gr <- GenomicFeatures::genes(TXDB, columns = "GENEID",
                                     single.strand.genes.only = FALSE)
  genes_gr <- unlist(genes_gr)
  genes_gr <- genes_gr[order(width(genes_gr), decreasing = TRUE)]
  genes_gr <- genes_gr[!duplicated(names(genes_gr))]
  sym_map <- sym_map[sym_map$ENTREZID %in% names(genes_gr), , drop = FALSE]
  missing_target <- setdiff(all_target, sym_map$SYMBOL)
  if (length(missing_target)) warning("Target genes with no mapped TSS (dropped): ",
                                      paste(missing_target, collapse = ", "))
  if (!nrow(sym_map)) stop_msg("No target genes mapped to a TSS.")
  
  target_tss <- genes_gr[sym_map$ENTREZID]
  names(target_tss) <- sym_map$SYMBOL
  target_prom <- GenomicRanges::promoters(target_tss, upstream = PROM_UP, downstream = PROM_DOWN)
  GenomeInfoDb::seqinfo(target_prom) <- GenomeInfoDb::seqinfo(GENOME)[GenomeInfoDb::seqlevels(target_prom)]
  target_prom <- IRanges::trim(target_prom)
  target_prom <- target_prom[width(target_prom) > 0]
  
  expr_symbols <- names(expr_lc)[is.finite(expr_lc) & expr_lc > 0]
  bg_syms_tab <- map_symbols_to_entrez(expr_symbols)
  bg_syms_tab <- bg_syms_tab[bg_syms_tab$ENTREZID %in% names(genes_gr), , drop = FALSE]
  bg_syms_tab <- bg_syms_tab[bg_syms_tab$SYMBOL %in% names(expr_lc), , drop = FALSE]
  if (nrow(bg_syms_tab) < BG_MIN_UNIVERSE) {
    stop_msg("Too few expressed genes mapped to Entrez/promoters for motif background: ",
             nrow(bg_syms_tab))
  }
  bg_universe <- bg_syms_tab$SYMBOL
  bg_expr <- expr_lc[bg_universe]
  cat("Background universe (expressed + mapped):", length(bg_universe), "\n")
  
  qs <- stats::quantile(bg_expr, probs = seq(0, 1, length.out = EXPR_BINS + 1L), na.rm = TRUE)
  qs <- unique(qs)
  bg_bin <- cut(bg_expr, breaks = qs, include.lowest = TRUE, labels = FALSE)
  names(bg_bin) <- bg_universe
  
  # ---------------------------------------------------------------------------
  # 5. PROMOTERS + MOTIFS + SCAN (assay extraction; GC bins)
  # ---------------------------------------------------------------------------
  bg_prom <- GenomicRanges::promoters(genes_gr[bg_syms_tab$ENTREZID],
                                      upstream = PROM_UP, downstream = PROM_DOWN)
  names(bg_prom) <- bg_syms_tab$SYMBOL
  GenomeInfoDb::seqinfo(bg_prom) <- GenomeInfoDb::seqinfo(GENOME)[GenomeInfoDb::seqlevels(bg_prom)]
  bg_prom <- IRanges::trim(bg_prom)
  bg_prom <- bg_prom[width(bg_prom) > 0]
  bg_syms_tab <- bg_syms_tab[bg_syms_tab$SYMBOL %in% names(bg_prom), , drop = FALSE]
  bg_universe <- bg_syms_tab$SYMBOL
  bg_expr <- expr_lc[bg_universe]
  bg_bin <- bg_bin[bg_universe]
  
  pwmlist_all <- TFBSTools::getMatrixSet(
    JASPAR2022::JASPAR2022,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
  motif_tf_all <- TFBSTools::name(pwmlist_all)
  pwmlist <- pwmlist_all[motif_tf_all %in% MOTIF_NAMES]
  if (!length(pwmlist)) stop_msg("No requested motifs found in JASPAR2022.")
  motif_tfs <- TFBSTools::name(pwmlist)
  motif_ids <- TFBSTools::ID(pwmlist)
  motif_labels <- make.unique(paste0(motif_tfs, "|", motif_ids))
  cat("Resolved motifs:", paste(motif_labels, collapse = ", "), "\n")
  
  target_seqs <- BSgenome::getSeq(GENOME, target_prom)
  names(target_seqs) <- names(target_prom)
  tm <- motifmatchr::matchMotifs(
    pwmlist, target_seqs,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    out = "matches", p.cutoff = MOTIF_P_CUTOFF)
  tm_mat <- as.matrix(SummarizedExperiment::assay(tm))
  rownames(tm_mat) <- names(target_seqs)
  colnames(tm_mat) <- motif_labels
  stopifnot(!any(is.na(tm_mat)),
            nrow(tm_mat) == length(target_seqs), ncol(tm_mat) == length(pwmlist))
  
  bg_seqs <- BSgenome::getSeq(GENOME, bg_prom)
  names(bg_seqs) <- names(bg_prom)
  bm <- motifmatchr::matchMotifs(
    pwmlist, bg_seqs,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    out = "matches", p.cutoff = MOTIF_P_CUTOFF)
  bm_mat <- as.matrix(SummarizedExperiment::assay(bm))
  rownames(bm_mat) <- names(bg_seqs)
  colnames(bm_mat) <- motif_labels
  stopifnot(!any(is.na(bm_mat)),
            nrow(bm_mat) == length(bg_seqs), ncol(bm_mat) == length(pwmlist))
  
  gc_all_bg <- Biostrings::letterFrequency(bg_seqs, "GC", as.prob = TRUE)[, 1]
  names(gc_all_bg) <- names(bg_seqs)
  gc_qs <- stats::quantile(gc_all_bg, probs = seq(0, 1, length.out = GC_BINS + 1L), na.rm = TRUE)
  gc_qs <- unique(gc_qs)
  gc_bin <- cut(gc_all_bg, breaks = gc_qs, include.lowest = TRUE, labels = FALSE)
  names(gc_bin) <- names(bg_seqs)
  stopifnot(identical(names(bg_bin), names(gc_bin)))
  
  # ---------------------------------------------------------------------------
  # 6. EXTENDED FEATURES + RESAMPLED ENRICHMENT (both-sided, per-draw BH)
  # ---------------------------------------------------------------------------
  tm2 <- tm_mat; bm2 <- bm_mat
  feat_labels <- motif_labels; feat_tfs <- motif_tfs; feat_ids <- motif_ids
  feat_cols <- seq_along(motif_labels)
  for (tf in unique(motif_tfs)) {
    cols <- which(motif_tfs == tf)
    if (length(cols) > 1L) {
      tm2 <- cbind(tm2, rowSums(tm_mat[, cols, drop = FALSE]) > 0)
      bm2 <- cbind(bm2, rowSums(bm_mat[, cols, drop = FALSE]) > 0)
      feat_labels <- c(feat_labels, paste0(tf, "_any"))
      feat_tfs <- c(feat_tfs, tf)
      feat_ids <- c(feat_ids, paste0(tf, "_any"))
      feat_cols <- c(feat_cols, ncol(tm2))
    }
  }
  cebp_orig <- which(motif_tfs %in% CEBP_FAMILY)
  if (length(cebp_orig) > 1L) {
    tm2 <- cbind(tm2, rowSums(tm_mat[, cebp_orig, drop = FALSE]) > 0)
    bm2 <- cbind(bm2, rowSums(bm_mat[, cebp_orig, drop = FALSE]) > 0)
    feat_labels <- c(feat_labels, "CEBP_family_any")
    feat_tfs <- c(feat_tfs, "CEBP_family_any")
    feat_ids <- c(feat_ids, "CEBP_family_any")
    feat_cols <- c(feat_cols, ncol(tm2))
  }
  is_cebp_family <- feat_tfs %in% CEBP_FAMILY | feat_tfs == "CEBP_family_any"
  
  exclude_syms <- unique(unlist(GENE_SETS, use.names = FALSE))
  set.seed(20260805)
  nf <- length(feat_labels)
  
  set_res <- lapply(names(GENE_SETS), function(gs_name) {
    targets <- intersect(GENE_SETS[[gs_name]], rownames(tm2))
    dropped <- setdiff(GENE_SETS[[gs_name]], targets)
    if (length(dropped)) warning("Targets dropped in ", gs_name, ": ",
                                 paste(dropped, collapse = ", "))
    missing_bg <- setdiff(targets, names(bg_bin))
    if (length(missing_bg)) {
      stop_msg("Targets absent from background universe (expression/mapping/promoter filter): ",
               paste(missing_bg, collapse = ", "))
    }
    stopifnot(!any(is.na(bg_bin[targets])), !any(is.na(gc_bin[targets])))
    if (length(targets) < 2L) return(NULL)
    t_expr_bin <- bg_bin[targets]
    t_gc_bin <- gc_bin[targets]
    cand <- setdiff(bg_universe, exclude_syms)
    
    draw_bg <- function() {
      chosen <- character(0)
      sourced_2d <- character(0)
      for (ti in seq_along(targets)) {
        pool2d <- intersect(cand, names(bg_bin)[bg_bin == t_expr_bin[ti] & gc_bin == t_gc_bin[ti]])
        if (length(pool2d)) {
          take <- sample(pool2d, min(BG_PER_TARGET, length(pool2d)))
          chosen <- c(chosen, take); sourced_2d <- c(sourced_2d, take)
        } else {
          pool1d <- intersect(cand, names(bg_bin)[bg_bin == t_expr_bin[ti]])
          if (length(pool1d)) chosen <- c(chosen, sample(pool1d, min(BG_PER_TARGET, length(pool1d))))
        }
      }
      chosen <- unique(chosen)
      if (length(chosen) < BG_MIN) {
        n_need <- BG_MIN - length(chosen)
        cand_rest <- setdiff(cand, chosen)
        cand_rest <- cand_rest[is.finite(bg_bin[cand_rest])]
        mid <- round(stats::median(t_expr_bin, na.rm = TRUE))
        cand_rest <- cand_rest[order(abs(bg_bin[cand_rest] - mid))]
        chosen <- unique(c(chosen, head(cand_rest, n_need)))
      }
      frac_2d <- if (length(chosen)) mean(chosen %in% sourced_2d) else NA_real_
      list(bg = chosen, frac_2d = frac_2d)
    }
    
    bg_draws <- lapply(seq_len(R_DRAWS), function(r) draw_bg())
    nb <- vapply(bg_draws, function(z) length(z$bg), integer(1))
    bg_2d_frac <- vapply(bg_draws, function(z) z$frac_2d, numeric(1))
    
    P_g <- matrix(NA_real_, nf, R_DRAWS)
    P_l <- matrix(NA_real_, nf, R_DRAWS)
    OR <- matrix(NA_real_, nf, R_DRAWS)
    for (fi in seq_len(nf)) {
      ci <- feat_cols[fi]
      a <- sum(tm2[targets, ci]); b <- length(targets) - a
      for (r in seq_len(R_DRAWS)) {
        bg <- bg_draws[[r]]$bg
        c <- sum(bm2[bg, ci]); d <- nb[r] - c
        tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        P_g[fi, r] <- suppressWarnings(fisher.test(tab, alternative = "greater")$p.value)
        P_l[fi, r] <- suppressWarnings(fisher.test(tab, alternative = "less")$p.value)
        OR[fi, r] <- (a + 0.5) / (b + 0.5) / ((c + 0.5) / (d + 0.5))
      }
    }
    FDR_g <- apply(P_g, 2, function(col) stats::p.adjust(col, "BH"))
    FDR_l <- apply(P_l, 2, function(col) stats::p.adjust(col, "BH"))
    
    gc_target <- gc_all_bg[targets]
    gc_bg_draws <- vapply(bg_draws, function(z) mean(gc_all_bg[z$bg]), numeric(1))
    bg_rate_median <- vapply(seq_len(nf), function(fi) {
      ci <- feat_cols[fi]
      median(vapply(bg_draws, function(z) mean(bm2[z$bg, ci]), numeric(1)))
    }, numeric(1))
    or_q <- apply(OR, 1, function(col) stats::quantile(col, c(0.025, 0.975), na.rm = TRUE))
    
    data.frame(
      gene_set = gs_name,
      motif_label = feat_labels, motif_tf = feat_tfs, motif_id = feat_ids,
      is_cebp_family = is_cebp_family,
      n_target = length(targets),
      target_motif_pos = colSums(tm2[targets, feat_cols, drop = FALSE]),
      target_rate = colSums(tm2[targets, feat_cols, drop = FALSE]) / length(targets),
      n_background_median = median(nb),
      background_rate_median = bg_rate_median,
      background_2d_stratum_frac = median(bg_2d_frac, na.rm = TRUE),
      median_odds_ratio = apply(OR, 1, median),
      or_q2.5 = or_q[1, ], or_q97.5 = or_q[2, ],
      median_p_greater = apply(P_g, 1, median),
      median_p_less = apply(P_l, 1, median),
      prop_fdr_greater_signif = rowMeans(FDR_g < FDR_THRESH),
      prop_fdr_less_signif = rowMeans(FDR_l < FDR_THRESH),
      target_mean_gc = mean(gc_target),
      background_mean_gc_median = median(gc_bg_draws),
      stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, set_res[!vapply(set_res, is.null, logical(1))])
  if (is.null(res) || !nrow(res)) stop_msg("No enrichment results; check gene sets.")
  res$call <- ifelse(res$median_odds_ratio > 1 & res$prop_fdr_greater_signif >= PROP_SIGNIF_CALL, "enriched",
                     ifelse(res$median_odds_ratio < 1 & res$prop_fdr_less_signif >= PROP_SIGNIF_CALL, "depleted", "ns"))
  write_csv_atomic(res, file.path(TABLE_DIR, "motif_enrichment_results.csv"))
  
  bg_audit_df <- res[!duplicated(res$gene_set),
                     c("gene_set","n_target","n_background_median",
                       "background_2d_stratum_frac",
                       "target_mean_gc","background_mean_gc_median")]
  write_csv_atomic(bg_audit_df, file.path(TABLE_DIR, "background_audit.csv"))
  gc_bad <- abs(bg_audit_df$target_mean_gc - bg_audit_df$background_mean_gc_median) > GC_WARN_DIFF
  if (any(gc_bad)) {
    warning("Promoter GC imbalance > ", GC_WARN_DIFF, " for: ",
            paste(bg_audit_df$gene_set[gc_bad], collapse = ", "),
            " - interpret enrichment with GC caution.")
  }
  
  gene_motif <- data.frame(symbol = rownames(tm_mat), tm_mat, check.names = FALSE)
  write_csv_atomic(gene_motif, file.path(TABLE_DIR, "target_gene_motif_presence.csv"))
  
  # ---------------------------------------------------------------------------
  # 7. PROVENANCE + CONSOLE SUMMARY
  # ---------------------------------------------------------------------------
  provenance <- data.frame(
    field = c("script","version","date","input_file",
              "promoter_upstream","promoter_downstream","motif_p_cutoff",
              "bg_per_target","bg_min","expr_bins","gc_bins","n_draws",
              "fdr_thresh","prop_signif_call","carrying_states",
              "tss_definition","motifs","n_target_genes_mapped",
              "n_background_universe","R_version"),
    value = c(SCRIPT_NAME, SCRIPT_VERSION, as.character(Sys.time()), INPUT_FILE,
              as.character(PROM_UP), as.character(PROM_DOWN), as.character(MOTIF_P_CUTOFF),
              as.character(BG_PER_TARGET), as.character(BG_MIN), as.character(EXPR_BINS),
              as.character(GC_BINS), as.character(R_DRAWS),
              as.character(FDR_THRESH), as.character(PROP_SIGNIF_CALL),
              paste(CARRYING_STATES, collapse = ";"),
              "genes(TxDb) extreme 5' extent; MANE anchoring not applied",
              paste(motif_labels, collapse = ";"),
              as.character(nrow(sym_map)),
              as.character(length(bg_universe)),
              as.character(getRversion())),
    stringsAsFactors = FALSE)
  write_csv_atomic(provenance, file.path(TABLE_DIR, "run_provenance.csv"))
  writeLines(capture.output(sessionInfo()), file.path(TABLE_DIR, "sessionInfo.txt"))
  
  cat("Script 06a v1.4 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("\n--- Motif enrichment (expression x GC matched background, ", R_DRAWS, " draws) ---\n", sep = "")
  print(res[, c("gene_set","motif_label","is_cebp_family",
                "target_motif_pos","n_target","target_rate",
                "background_rate_median","median_odds_ratio","or_q2.5","or_q97.5",
                "median_p_greater","prop_fdr_greater_signif","prop_fdr_less_signif","call")],
        row.names = FALSE)
  cat("\n--- Background audit ---\n")
  print(bg_audit_df, row.names = FALSE)
  cat("\nNotes:\n")
  cat(" - Sequence-level plausibility only; not regulation evidence.\n")
  cat(" - M2_full is the main enrichment set; primary_core/residual_neutrophil (4 genes) are direction-support
  only.\n")
  cat(" - prop_fdr_signif = stability of the call under background resampling, not a replication probability.\n")
  cat(" - median_p is descriptive, not an inferential p-value.\n")
  cat(" - or_q2.5/or_q97.5 are resampling stability quantiles, not a confidence interval.\n")
  cat(" - CEBP_family_any may saturate (target_rate ~ background_rate ~ 0.8 at p.cutoff 1e-4 / 4kb);\n")
  cat("   that is a ceiling/power artifact, not evidence against CEBP plausibility.\n")
  cat(" - Gene sets are nested; a hit in two nested sets is one piece of evidence, not two.\n")
})