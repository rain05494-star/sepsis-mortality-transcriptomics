# ============================================================================
# SCRIPT 02 v6 -- GSE216009 32-gene major cell-state source attribution
#   FINAL figures: labeled UMAP, facet+shared-legend gene UMAPs, unified blue
#   logCPM / white-red contribution palettes, "not evaluable" grey, short_state
#   labels, low-detected gene asterisks, Fig04 primary-first ordering.
#   Analysis logic unchanged from v4.
# RUN: Rscript 02_gene_cell_distribution.R
# ============================================================================

options(stringsAsFactors = FALSE)

local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  PROJECT_ROOT <- "/home/sunshine/predicate/singlecell"
  INPUT_RDS    <- file.path(PROJECT_ROOT, "GSE216009_rhapsody_wholeblood_sobj.rds.gz")
  
  SCRIPT01_RUN_DIR <- file.path(
    PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis",
    "01_gene_mapping_detectability", "20260803_gene_mapping_detectability"
  )
  
  OUTPUTS_ROOT <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  STAGE_NAME   <- "02_gene_cell_distribution"
  RUN_PURPOSE  <- "gene_cell_distribution"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG       <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  
  MIN_CELLS_PER_DONOR_STATE <- 20L
  LOW_UMI_FOR_CONTRIBUTION  <- 3L
  CONF_AMBIGUOUS_DELTA      <- 0.05
  CONF_HIGH_DELTA           <- 0.10
  CONF_HIGH_ALL26           <- 0.50
  CONF_HIGH_ELIGIBLE_N      <- 8L
  CONF_MODERATE_ALL26       <- 0.30
  
  EXPECTED_N_SAMPLES        <- 48L
  EXPECTED_N_DONORS         <- 39L
  EXPECTED_N_ACUTE_DONORS   <- 26L
  EXPECTED_N_HC_DONORS      <- 6L
  EXPECTED_N_SURGERY_DONORS <- 7L
  EXPECTED_N_CONV_SAMPLES   <- 9L
  EXPECTED_N_PAIRED_DONORS  <- 9L
  EXPECTED_N_PRIMARY_CELLS  <- 151837L
  
  UMAP_GENES <- c("MKI67","TOP2A","ELANE","MPO","CHI3L1","SERPINB10","CX3CR1","IL1B")
  # NOTE: 2 stable genes per module, chosen for display only (spatially resolved);
  #       module-level inference always uses the full gene sets (Script 03 UCell).
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS + PACKAGES
  # ---------------------------------------------------------------------------
  needed <- c("SeuratObject","Matrix","ggplot2","scales","gridExtra","grid")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({
    library(SeuratObject); library(Matrix); library(ggplot2)
    library(scales); library(gridExtra); library(grid)
  })
  
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  require_cols <- function(x, cols, label) {
    m <- setdiff(cols, names(x))
    if (length(m)) stop_msg("Missing required columns in ", label, ": ", paste(m, collapse = ", "))
    invisible(TRUE)
  }
  write_csv_atomic <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    
    x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    
    # Excel / Power Query friendly:
    # write missing numeric values as blank cells, not literal "NA".
    num_cols <- vapply(x, is.numeric, logical(1))
    for (nm in names(x)[num_cols]) {
      bad <- !is.finite(x[[nm]])
      if (any(bad, na.rm = TRUE)) {
        x[[nm]][bad] <- NA_real_
      }
    }
    
    tmp <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = dirname(path),
      fileext = ".tmp"
    )
    
    write.csv(
      x,
      tmp,
      row.names = FALSE,
      fileEncoding = "UTF-8",
      na = ""
    )
    
    if (file.exists(path)) unlink(path)
    ok <- file.rename(tmp, path)
    if (!ok) {
      unlink(tmp)
      stop_msg("Failed to move temp CSV to final path: ", path)
    }
    
    invisible(path)
  }
  matrix_to_csv_df <- function(mat, row_id = "id") {
    data.frame(setNames(list(rownames(mat)), row_id),
               as.data.frame(mat, check.names = FALSE),
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  save_pdf_png <- function(p,
                           base,
                           dir = FIG_DIR,
                           width = 7,
                           height = 6,
                           dpi = 600,
                           make_pdf = TRUE,
                           make_png = TRUE,
                           make_tiff = TRUE,
                           make_svg = TRUE) {
    
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    
    pdf_file  <- file.path(dir, paste0(base, ".pdf"))
    png_file  <- file.path(dir, paste0(base, ".png"))
    tiff_file <- file.path(dir, paste0(base, ".tiff"))
    svg_file  <- file.path(dir, paste0(base, ".svg"))
    
    draw_obj <- function(x) {
      if (inherits(x, "ggplot")) {
        print(x)
      } else {
        grid::grid.newpage()
        grid::grid.draw(x)
      }
    }
    
    if (make_pdf) {
      grDevices::pdf(
        file = pdf_file,
        width = width,
        height = height,
        onefile = FALSE,
        family = "Helvetica",
        useDingbats = FALSE
      )
      draw_obj(p)
      grDevices::dev.off()
    }
    
    if (make_png) {
      ggsave(
        filename = png_file,
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white",
        limitsize = FALSE
      )
    }
    
    if (make_tiff) {
      ggsave(
        filename = tiff_file,
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        device = "tiff",
        compression = "lzw",
        bg = "white",
        limitsize = FALSE
      )
    }
    
    if (make_svg && requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(
        file = svg_file,
        width = width,
        height = height,
        bg = "white"
      )
      draw_obj(p)
      grDevices::dev.off()
    }
    
    invisible(TRUE)
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
  module_label <- function(x) {
    y <- as.character(x)
    y[y %in% c("Cell cycle / proliferation","Cell cycle")] <- "Cell cycle"
    y[y %in% c("Neutrophil degranulation","Neutrophil")] <- "Neutrophil"
    y[y %in% c("Inflammatory / Down","Inflammatory")] <- "Inflammatory"
    y[y %in% c("Other")] <- "Other"
    y
  }
  NOT_EVAL <- "not evaluable"
  short_state <- function(x) {
    z <- as.character(x)
    map <- c(
      "HSPCs" = "HSPCs",
      "Cycling_neutrophil_progenitors" = "Cycling neut prog",
      "MPO+_immature_neutrophils_or_progenitors" = "MPO+ imm neut/prog",
      "PADI4+_immature_neutrophils" = "PADI4+ imm neut",
      "IL1R2+_immature_neutrophils" = "IL1R2+ imm neut",
      "S100A8-9_hi_neutrophils" = "S100A8/A9 hi neut",
      "Mature_neutrophils" = "Mature neut",
      "Degranulating_neutrophils" = "Degran neut",
      "Apoptosing_neutrophils" = "Apoptosing neut",
      "Eosinophils" = "Eosinophils",
      "Mast_cells/eosiniophils" = "Mast/eos",
      "Classical_monocytes" = "Classical mono",
      "Non-classical_monocytes" = "Non-classical mono",
      "cDCs" = "cDCs",
      "pDCs" = "pDCs",
      "Naive_CD4_T_cells" = "Naive CD4+ T",
      "Memory_CD4_T_cells" = "Memory CD4+ T",
      "Naive_CD8_T_cells" = "Naive CD8+ T",
      "CD8_T_cells" = "CD8+ T",
      "Cycling_TNK" = "Cycling T/NK",
      "NK cells" = "NK cells",
      "B_cells" = "B cells",
      "Plasmablasts" = "Plasmablasts",
      "Platelets" = "Platelets"
    )
    out <- ifelse(z %in% names(map), unname(map[z]), z)
    out
  }
  assign_bin <- function(v, limits, labels) {
    if (!is.finite(v)) return(NOT_EVAL)
    idx <- which(v <= limits)[1L]
    if (is.na(idx)) idx <- length(labels)
    labels[idx]
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  PROV_DIR  <- file.path(OUT_ROOT, "provenance")
  for (d in c(TABLE_DIR, FIG_DIR, PROV_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. READ INPUTS
  # ---------------------------------------------------------------------------
  s01_tables <- c("frozen_signature_input.csv","gene_mapping_audit.csv",
                  "final_gene_evaluability_32.csv","sample_to_donor_mapping_audit.csv",
                  "fine_state_display_order_audit.csv")
  missing_t <- s01_tables[!file.exists(file.path(SCRIPT01_RUN_DIR, "tables", s01_tables))]
  prov_ok <- file.exists(file.path(SCRIPT01_RUN_DIR, "provenance", "run_provenance.csv"))
  if (length(missing_t) || !prov_ok)
    stop_msg("Script-01 run incomplete: missing tables ",
             paste(missing_t, collapse = "; "),
             if (!prov_ok) "; provenance/run_provenance.csv" else "",
             "\nFix SCRIPT01_RUN_DIR.")
  
  sig01   <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "frozen_signature_input.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  map01   <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "gene_mapping_audit.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  eva01   <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "final_gene_evaluability_32.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  state01 <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "fine_state_display_order_audit.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  
  require_cols(sig01, c("gene","module","direction","E3_membership","canonical_hgnc"), "frozen_signature_input.csv")
  require_cols(map01, c("signature_gene","dataset_feature"), "gene_mapping_audit.csv")
  require_cols(eva01, c("signature_gene","class"), "final_gene_evaluability_32.csv")
  require_cols(state01, c("display_order","fine_state","lineage","color"), "fine_state_display_order_audit.csv")
  stopifnot(nrow(sig01) == 32L, nrow(eva01) == 32L, nrow(state01) == 24L, !anyDuplicated(state01$fine_state))
  
  fine_states   <- state01$fine_state[order(state01$display_order)]
  state_lineage <- setNames(state01$lineage, state01$fine_state)
  state_colors  <- setNames(state01$color, state01$fine_state)
  
  # Vivid display palette for cell-state visualization.
  # This changes plotting colors only; annotations and source-attribution results are unchanged.
  USE_VIVID_STATE_COLORS <- TRUE
  
  if (USE_VIVID_STATE_COLORS) {
    vivid_state_colors <- c(
      "HSPCs" = "#7A5CCF",
      "Cycling_neutrophil_progenitors" = "#00B8FF",
      "MPO+_immature_neutrophils_or_progenitors" = "#00A5E3",
      "PADI4+_immature_neutrophils" = "#0086D1",
      "IL1R2+_immature_neutrophils" = "#006CB8",
      "S100A8-9_hi_neutrophils" = "#00549E",
      "Mature_neutrophils" = "#0072B2",
      "Degranulating_neutrophils" = "#004C99",
      "Apoptosing_neutrophils" = "#003B73",
      
      "Eosinophils" = "#E85D9A",
      "Mast_cells/eosiniophils" = "#B83B7D",
      
      "Classical_monocytes" = "#2FBF71",
      "Non-classical_monocytes" = "#009E4D",
      "cDCs" = "#007A3D",
      "pDCs" = "#006B4F",
      
      "Naive_CD4_T_cells" = "#F28E2B",
      "Memory_CD4_T_cells" = "#D55E00",
      "Naive_CD8_T_cells" = "#E64B35",
      "CD8_T_cells" = "#C0003B",
      "Cycling_TNK" = "#FF6F69",
      "NK cells" = "#A63603",
      
      "B_cells" = "#9B59B6",
      "Plasmablasts" = "#8E007E",
      "Platelets" = "#D9A300"
    )
    
    missing_cols <- setdiff(fine_states, names(vivid_state_colors))
    extra_cols <- setdiff(names(vivid_state_colors), fine_states)
    
    if (length(missing_cols)) {
      stop_msg("Vivid palette missing fine states: ", paste(missing_cols, collapse = ", "))
    }
    if (length(extra_cols)) {
      stop_msg("Vivid palette has unused fine states: ", paste(extra_cols, collapse = ", "))
    }
    
    state_colors <- vivid_state_colors[fine_states]
  }
  
  lineage_order <- unique(state01$lineage[order(state01$display_order)])
  lineage_colors <- setNames(
    vapply(lineage_order, function(lin) {
      fs <- state01$fine_state[state01$lineage == lin][1L]
      unname(state_colors[fs])
    }, character(1)), lineage_order)
  
  sig01 <- sig01[order(match(module_label(sig01$module), c("Cell cycle","Neutrophil","Other","Inflammatory")),
                       sig01$gene), ]
  genes_full  <- sig01$gene
  feat_map    <- setNames(map01$dataset_feature, map01$signature_gene)
  gene_class  <- setNames(eva01$class, eva01$signature_gene)
  genes_primary <- eva01$signature_gene[eva01$class == "stable_detected"]
  stopifnot(length(genes_primary) == 29L)
  stopifnot(setequal(setdiff(genes_full, genes_primary), c("HIST1H2BM","HIST1H3B","RHAG")))
  
  if (!file.exists(INPUT_RDS)) stop_msg("Input RDS not found: ", INPUT_RDS)
  obj <- readRDS(INPUT_RDS)
  meta <- obj@meta.data
  if (!is.data.frame(meta)) stop_msg("obj@meta.data is not a data.frame.")
  stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
  
  stopifnot(COUNT_ASSAY %in% Assays(obj))
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    stopifnot(COUNT_LAYER %in% Layers(obj[[COUNT_ASSAY]]))
    cnt <- LayerData(obj, assay = COUNT_ASSAY, layer = COUNT_LAYER, fast = FALSE)
  } else {
    stopifnot(COUNT_LAYER == "counts")
    cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, slot = "counts")
  }
  stopifnot(is(cnt, "dgCMatrix"), identical(colnames(cnt), colnames(obj)))
  stopifnot(c("UMAP_1","UMAP_2") %in% colnames(meta),
            all(is.finite(meta$UMAP_1)), all(is.finite(meta$UMAP_2)))
  
  # ---------------------------------------------------------------------------
  # 3. COHORT DEFINITION (fail-closed)
  # ---------------------------------------------------------------------------
  require_cols(meta, c(SAMPLE_FIELD, CONDITION_FIELD, FINE_STATE_FIELD), "obj@meta.data")
  meta$sample_id <- as.character(meta[[SAMPLE_FIELD]])
  meta$diagnosis_for_script02 <- as.character(meta[[CONDITION_FIELD]])
  
  sample_diag <- unique(meta[, c("sample_id","diagnosis_for_script02"), drop = FALSE])
  n_diag_per_sample <- tapply(sample_diag$diagnosis_for_script02, sample_diag$sample_id,
                              function(z) length(unique(z)))
  if (any(n_diag_per_sample != 1L)) stop_msg("sample_id with multiple diagnosis: ",
                                             paste(names(n_diag_per_sample)[n_diag_per_sample != 1L], collapse = ", "))
  
  sample_tab <- sample_diag
  sample_tab$donor_id <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort_group <- NA_character_
  sample_tab$cohort_group[sample_tab$diagnosis_for_script02 %in% ACUTE_DIAG] <- "Acute_sepsis"
  sample_tab$cohort_group[sample_tab$diagnosis_for_script02 == "Conv"] <- "Convalescent"
  sample_tab$cohort_group[sample_tab$diagnosis_for_script02 == "HV"] <- "Healthy_control"
  sample_tab$cohort_group[sample_tab$diagnosis_for_script02 == "CS"] <- "Surgery_control"
  if (anyNA(sample_tab$cohort_group)) stop_msg("Unmapped diagnosis: ",
                                               paste(unique(sample_tab$diagnosis_for_script02[is.na(sample_tab$cohort_group)]), collapse = ", "))
  
  is_conv_suffix <- grepl("_CONV$", sample_tab$sample_id)
  is_conv_diag   <- sample_tab$diagnosis_for_script02 == "Conv"
  stopifnot(identical(is_conv_suffix, is_conv_diag))
  
  n_samples <- length(unique(sample_tab$sample_id))
  n_donors  <- length(unique(sample_tab$donor_id))
  acute_samples  <- sample_tab[sample_tab$cohort_group == "Acute_sepsis", , drop = FALSE]
  hc_samples     <- sample_tab[sample_tab$cohort_group == "Healthy_control", , drop = FALSE]
  surgery_samples<- sample_tab[sample_tab$cohort_group == "Surgery_control", , drop = FALSE]
  conv_samples   <- sample_tab[sample_tab$cohort_group == "Convalescent", , drop = FALSE]
  acute_donors   <- unique(acute_samples$donor_id)
  hc_donors      <- unique(hc_samples$donor_id)
  surgery_donors <- unique(surgery_samples$donor_id)
  conv_donors    <- unique(conv_samples$donor_id)
  paired_donors  <- intersect(acute_donors, conv_donors)
  
  stopifnot(n_samples == EXPECTED_N_SAMPLES, n_donors == EXPECTED_N_DONORS)
  stopifnot(length(acute_donors) == EXPECTED_N_ACUTE_DONORS)
  stopifnot(length(hc_donors) == EXPECTED_N_HC_DONORS)
  stopifnot(length(surgery_donors) == EXPECTED_N_SURGERY_DONORS)
  stopifnot(nrow(conv_samples) == EXPECTED_N_CONV_SAMPLES)
  stopifnot(length(paired_donors) == EXPECTED_N_PAIRED_DONORS)
  stopifnot(all(conv_donors %in% acute_donors))
  
  acute_n_per_donor <- table(acute_samples$donor_id)
  if (any(acute_n_per_donor != 1L)) stop_msg("Acute donors without exactly one acute sample: ",
                                             paste(names(acute_n_per_donor)[acute_n_per_donor != 1L], collapse = ", "))
  if (!all(grepl("_CONV$", conv_samples$sample_id))) stop_msg("Convalescent sample without _CONV suffix.")
  
  m <- match(meta$sample_id, sample_tab$sample_id)
  if (anyNA(m)) stop_msg("meta$sample_id absent from sample_tab.")
  meta$donor_id    <- sample_tab$donor_id[m]
  meta$cohort_group<- sample_tab$cohort_group[m]
  meta$is_primary_acute <- meta$cohort_group == "Acute_sepsis"
  stopifnot(all(!meta$is_primary_acute[grepl("_CONV$", meta$sample_id)]))
  
  primary_idx <- which(meta$is_primary_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  
  write_csv_atomic(data.frame(
    metric = c("n_samples","n_donors","n_acute_donors","n_hc_donors","n_surgery_donors",
               "n_conv_samples","n_paired_donors","n_primary_cells"),
    observed = c(n_samples, n_donors, length(acute_donors), length(hc_donors),
                 length(surgery_donors), nrow(conv_samples), length(paired_donors),
                 length(primary_idx)),
    expected = c(EXPECTED_N_SAMPLES, EXPECTED_N_DONORS, EXPECTED_N_ACUTE_DONORS,
                 EXPECTED_N_HC_DONORS, EXPECTED_N_SURGERY_DONORS,
                 EXPECTED_N_CONV_SAMPLES, EXPECTED_N_PAIRED_DONORS, EXPECTED_N_PRIMARY_CELLS)),
    file.path(TABLE_DIR, "cohort_condition_audit.csv"))
  cat("Cohort passed. Acute primary cells:", length(primary_idx), "\n")
  
  # ---------------------------------------------------------------------------
  # 4. PRIMARY CELLS + DESIGN MATRICES
  # ---------------------------------------------------------------------------
  meta_p <- meta[primary_idx, , drop = FALSE]
  cnt_p  <- cnt[, primary_idx, drop = FALSE]
  donor_p <- factor(meta_p$donor_id, levels = acute_donors)
  state_p <- factor(meta_p[[FINE_STATE_FIELD]], levels = fine_states)
  stopifnot(!anyNA(donor_p), !anyNA(state_p),
            identical(levels(donor_p), acute_donors),
            setequal(unique(meta_p[[FINE_STATE_FIELD]]), fine_states))
  
  n_donor <- nlevels(donor_p); n_state <- length(fine_states)
  
  ds_idx <- (as.integer(donor_p) - 1L) * n_state + as.integer(state_p)
  Ds <- sparseMatrix(i = ds_idx, j = seq_along(ds_idx), x = 1,
                     dims = c(n_donor * n_state, length(ds_idx)))
  
  ncell_ds_flat <- Matrix::rowSums(Ds)
  n_cells_ds <- matrix(as.integer(ncell_ds_flat), nrow = n_donor, ncol = n_state, byrow = TRUE)
  rownames(n_cells_ds) <- levels(donor_p); colnames(n_cells_ds) <- fine_states
  eligible_ds <- n_cells_ds >= MIN_CELLS_PER_DONOR_STATE
  
  lib_cell <- as.numeric(Matrix::colSums(cnt_p))
  lib_ds_flat <- as.numeric(Ds %*% lib_cell)
  lib_ds <- matrix(lib_ds_flat, nrow = n_donor, ncol = n_state, byrow = TRUE)
  rownames(lib_ds) <- levels(donor_p); colnames(lib_ds) <- fine_states
  stopifnot(all(lib_ds[eligible_ds] > 0))
  
  feat <- rownames(cnt)
  idx <- match(feat_map[genes_full], feat)
  stopifnot(!anyNA(idx), !anyDuplicated(idx))
  X_p <- cnt_p[idx, , drop = FALSE]; rownames(X_p) <- genes_full
  C_flat <- as.matrix(X_p %*% Matrix::t(Ds))
  stopifnot(nrow(C_flat) == length(genes_full), ncol(C_flat) == n_donor * n_state)
  
  T_gd <- matrix(0, nrow = length(genes_full), ncol = n_donor,
                 dimnames = list(genes_full, levels(donor_p)))
  for (gi in seq_along(genes_full)) {
    Cg <- matrix(C_flat[gi, ], nrow = n_donor, ncol = n_state, byrow = TRUE)
    T_gd[gi, ] <- rowSums(Cg)
  }
  low_umi_gd <- T_gd < LOW_UMI_FOR_CONTRIBUTION
  
  # ---------------------------------------------------------------------------
  # 5. THREE-LAYER MATRICES
  # ---------------------------------------------------------------------------
  logcpm_gds <- array(NA_real_, dim = c(length(genes_full), n_donor, n_state),
                      dimnames = list(genes_full, levels(donor_p), fine_states))
  cf_gds     <- array(NA_real_, dim = c(length(genes_full), n_donor, n_state),
                      dimnames = list(genes_full, levels(donor_p), fine_states))
  cellfrac_ds <- sweep(n_cells_ds, 1, rowSums(n_cells_ds), "/")
  
  X_det_p <- X_p > 0
  det_cells_flat <- as.matrix(X_det_p %*% Matrix::t(Ds))
  
  for (gi in seq_along(genes_full)) {
    Cg <- matrix(C_flat[gi, ], nrow = n_donor, ncol = n_state, byrow = TRUE)
    cf <- Cg / T_gd[gi, ]; cf[T_gd[gi, ] == 0, ] <- NA_real_
    cf_gds[gi, , ] <- cf
    cpm <- Cg / lib_ds * 1e6
    lc <- log2(cpm + 1); lc[!eligible_ds] <- NA_real_
    logcpm_gds[gi, , ] <- lc
  }
  
  long_rows <- do.call(rbind, lapply(seq_along(genes_full), function(gi) {
    g <- unname(genes_full[gi])
    g_feature <- unname(feat_map[g])
    g_module <- unname(module_label(sig01$module[match(g, sig01$gene)]))
    g_class <- unname(gene_class[g])
    
    do.call(rbind, lapply(seq_len(n_donor), function(di) {
      d <- unname(levels(donor_p)[di])
      j <- ((di - 1L) * n_state + 1L):(di * n_state)
      
      data.frame(
        gene = rep(g, n_state),
        dataset_feature = rep(g_feature, n_state),
        module = rep(g_module, n_state),
        evaluable_class = rep(g_class, n_state),
        donor = rep(d, n_state),
        fine_state = unname(fine_states),
        lineage = unname(state_lineage[fine_states]),
        n_cells = unname(as.integer(n_cells_ds[di, ])),
        eligible_donor_state = unname(as.logical(eligible_ds[di, ])),
        library_size = unname(as.numeric(lib_ds[di, ])),
        umi = unname(as.numeric(C_flat[gi, j])),
        donor_gene_total_umi = rep(unname(as.numeric(T_gd[gi, di])), n_state),
        low_total_umi_for_contribution = rep(unname(as.logical(low_umi_gd[gi, di])), n_state),
        logCPM = unname(as.numeric(logcpm_gds[gi, di, ])),
        contribution_fraction = unname(as.numeric(cf_gds[gi, di, ])),
        cell_fraction = unname(as.numeric(cellfrac_ds[di, ])),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }))
  write_csv_atomic(long_rows, file.path(TABLE_DIR, "donor_state_gene_counts_long.csv"))
  write_csv_atomic(
    long_rows[, c("gene","dataset_feature","module","evaluable_class","donor",
                  "fine_state","lineage","n_cells","eligible_donor_state",
                  "library_size","logCPM"), drop = FALSE],
    file.path(TABLE_DIR, "donor_state_logCPM_long.csv"))
  write_csv_atomic(
    long_rows[, c("gene","dataset_feature","module","evaluable_class","donor",
                  "fine_state","lineage","n_cells","eligible_donor_state",
                  "umi","donor_gene_total_umi","low_total_umi_for_contribution",
                  "contribution_fraction","cell_fraction"), drop = FALSE],
    file.path(TABLE_DIR, "donor_state_contribution_fraction_long.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. PER-GENE STATE SUMMARIES
  # ---------------------------------------------------------------------------
  state_rank_df <- list()
  for (gi in seq_along(genes_full)) {
    g <- genes_full[gi]
    Cg <- matrix(C_flat[gi, ], nrow = n_donor, ncol = n_state, byrow = TRUE)
    Dg <- matrix(det_cells_flat[gi, ], nrow = n_donor, ncol = n_state, byrow = TRUE)
    Tg <- T_gd[gi, ]
    useT3 <- which(Tg >= LOW_UMI_FOR_CONTRIBUTION)
    
    med_cf_all26 <- apply(cf_gds[gi, , ], 2, median, na.rm = TRUE)
    mean_cf_all26 <- apply(cf_gds[gi, , ], 2, mean, na.rm = TRUE)
    med_cf_T3 <- vapply(seq_len(n_state), function(si)
      if (length(useT3)) median(cf_gds[gi, useT3, si], na.rm = TRUE) else NA_real_, numeric(1))
    med_logcpm <- apply(logcpm_gds[gi, , ], 2, median, na.rm = TRUE)
    
    det_all26 <- colSums(Cg >= 1) / n_donor
    det_elig  <- ifelse(colSums(eligible_ds) > 0,
                        colSums(Cg >= 1 & eligible_ds) / colSums(eligible_ds), NA_real_)
    
    cell_det_frac <- Dg / n_cells_ds
    cell_det_frac[n_cells_ds == 0] <- NA_real_
    med_cell_det <- apply(ifelse(eligible_ds, cell_det_frac, NA_real_), 2, median, na.rm = TRUE)
    
    cv_logcpm <- vapply(seq_len(n_state), function(si) {
      use <- eligible_ds[, si] & Cg[, si] >= 1
      if (sum(use) < 2) return(NA_real_)
      v <- logcpm_gds[gi, use, si]
      mu <- mean(v)
      if (!is.finite(mu) || mu < 0.1) return(NA_real_)
      sd(v) / mu
    }, numeric(1))
    
    state_rank_df[[gi]] <- data.frame(
      signature_gene = g, fine_state = fine_states,
      lineage = unname(state_lineage[fine_states]),
      median_contribution_fraction_all26 = med_cf_all26,
      mean_contribution_fraction_all26 = mean_cf_all26,
      median_contribution_fraction_Tge3 = med_cf_T3,
      median_logCPM = med_logcpm,
      donor_detection_rate_all26 = det_all26,
      donor_detection_rate_eligible = det_elig,
      cell_detection_fraction_median = med_cell_det,
      logCPM_CV = cv_logcpm, stringsAsFactors = FALSE)
  }
  state_rank <- do.call(rbind, state_rank_df)
  write_csv_atomic(state_rank, file.path(TABLE_DIR, "per_gene_all_state_rankings.csv"))
  
  summary_rows <- lapply(seq_along(genes_full), function(gi) {
    g <- genes_full[gi]
    Cg <- matrix(C_flat[gi, ], nrow = n_donor, ncol = n_state, byrow = TRUE)
    Tg <- T_gd[gi, ]
    r <- state_rank[state_rank$signature_gene == g, , drop = FALSE]
    r <- r[order(-r$median_contribution_fraction_all26, -r$mean_contribution_fraction_all26,
                 -r$donor_detection_rate_all26, match(r$fine_state, fine_states)), ]
    c1 <- r[1, ]; c2 <- r[2, ]
    e <- r[order(-r$median_logCPM, -r$donor_detection_rate_all26,
                 match(r$fine_state, fine_states)), ]
    e1 <- e[1, ]; e2 <- e[2, ]
    
    contrib_delta <- c1$median_contribution_fraction_all26 - c2$median_contribution_fraction_all26
    logcpm_delta  <- e1$median_logCPM - e2$median_logCPM
    
    si1 <- match(c1$fine_state, fine_states)
    if (is.na(si1)) stop_msg("main contributor state not found in fine_states: ", c1$fine_state)
    eligible_vec <- eligible_ds[, si1]
    eligible_n <- sum(eligible_vec)
    detected_n <- sum(eligible_vec & Cg[, si1] >= 1)
    usable_n   <- sum(Tg >= LOW_UMI_FOR_CONTRIBUTION)
    low_gene   <- gene_class[g] != "stable_detected"
    
    ambiguous  <- !is.finite(contrib_delta) || contrib_delta < CONF_AMBIGUOUS_DELTA
    same_lineage <- c1$lineage == e1$lineage
    if (low_gene || ambiguous || c1$donor_detection_rate_all26 < CONF_MODERATE_ALL26) {
      conf <- "low"
    } else if (c1$donor_detection_rate_all26 >= CONF_HIGH_ALL26 &&
               contrib_delta >= CONF_HIGH_DELTA &&
               eligible_n >= CONF_HIGH_ELIGIBLE_N &&
               (c1$fine_state == e1$fine_state || same_lineage)) {
      conf <- "high"
    } else {
      conf <- "moderate"
    }
    
    data.frame(
      signature_gene = g, dataset_feature = feat_map[g],
      canonical_hgnc = sig01$canonical_hgnc[match(g, sig01$gene)],
      module = module_label(sig01$module[match(g, sig01$gene)]),
      direction = sig01$direction[match(g, sig01$gene)],
      E3_membership = sig01$E3_membership[match(g, sig01$gene)],
      evaluable_class = gene_class[g],
      primary_analysis_flag = gene_class[g] == "stable_detected",
      main_contributor_state = c1$fine_state,
      main_contributor_lineage = c1$lineage,
      main_contribution_fraction_median_all26 = c1$median_contribution_fraction_all26,
      main_contribution_fraction_median_Tge3 = c1$median_contribution_fraction_Tge3,
      second_contributor_state = c2$fine_state,
      second_contribution_fraction_median = c2$median_contribution_fraction_all26,
      contribution_delta_main_minus_second = contrib_delta,
      usable_contribution_donor_n_Tge3 = usable_n,
      usable_contribution_donor_rate_Tge3 = usable_n / n_donor,
      main_enriched_state = e1$fine_state,
      main_enriched_lineage = e1$lineage,
      main_logCPM_median = e1$median_logCPM,
      second_enriched_state = e2$fine_state,
      second_logCPM_median = e2$median_logCPM,
      logCPM_delta_main_minus_second = logcpm_delta,
      main_state_donor_detection_rate_all26 = c1$donor_detection_rate_all26,
      main_state_donor_detection_rate_eligible = c1$donor_detection_rate_eligible,
      main_state_detected_donor_n = detected_n,
      main_state_eligible_donor_n = eligible_n,
      main_state_cell_detection_fraction_median = c1$cell_detection_fraction_median,
      main_state_logCPM_CV = c1$logCPM_CV,
      main_state_cell_fraction_median = median(cellfrac_ds[, si1]),
      main_state_n_cells_median = median(n_cells_ds[, si1]),
      source_call_ambiguous = ambiguous,
      source_call_confidence = conf,
      source_call_note = ifelse(low_gene, "low_detected; descriptive only",
                                ifelse(ambiguous, "contributor/second delta < 0.05",
                                       ifelse(c1$fine_state == e1$fine_state,
                                              "contributor == enriched",
                                              "contributor vs enriched discordant"))),
      stringsAsFactors = FALSE)
  })
  summary_df <- do.call(rbind, summary_rows)
  write_csv_atomic(summary_df, file.path(TABLE_DIR, "per_gene_cell_source_summary.csv"))
  
  median_logcpm_mat <- do.call(rbind, lapply(seq_along(genes_full), function(gi)
    apply(logcpm_gds[gi, , ], 2, median, na.rm = TRUE)))
  rownames(median_logcpm_mat) <- genes_full; colnames(median_logcpm_mat) <- fine_states
  median_cf_mat <- do.call(rbind, lapply(seq_along(genes_full), function(gi)
    apply(cf_gds[gi, , ], 2, median, na.rm = TRUE)))
  rownames(median_cf_mat) <- genes_full; colnames(median_cf_mat) <- fine_states
  
  write_csv_atomic(matrix_to_csv_df(median_logcpm_mat, "gene"), file.path(TABLE_DIR,
                                                                          "median_state_logCPM_matrix.csv"))
  write_csv_atomic(matrix_to_csv_df(median_cf_mat, "gene"), file.path(TABLE_DIR,
                                                                      "median_state_contribution_fraction_matrix.csv"))
  write_csv_atomic(matrix_to_csv_df(n_cells_ds, "donor"), file.path(TABLE_DIR, "donor_state_n_cells.csv"))
  write_csv_atomic(matrix_to_csv_df(lib_ds, "donor"), file.path(TABLE_DIR, "donor_state_library_size.csv"))
  write_csv_atomic(matrix_to_csv_df(cellfrac_ds, "donor"), file.path(TABLE_DIR, "donor_state_cell_fraction.csv"))
  
  # ---------------------------------------------------------------------------
  # 7. STYLE + FIGURES  (v6 final)
  # ---------------------------------------------------------------------------
  mod_pal <- c("Cell cycle" = "#2E6FD8","Neutrophil" = "#4cbef7",
               "Other" = "#62bb78","Inflammatory" = "#E86518")
  
  style_theme <- theme_classic(base_size = 11) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          axis.line = element_line(color = "grey40", linewidth = 0.4),
          axis.ticks = element_line(color = "grey30", linewidth = 0.3),
          panel.grid = element_blank(),
          plot.title = element_text(size = 11, hjust = 0, color = "black", margin = margin(b = 6)),
          plot.caption = element_text(size = 8, hjust = 0, color = "grey35"),
          plot.margin = margin(8, 8, 8, 8), legend.position = "right",
          legend.title = element_text(size = 9), legend.text = element_text(size = 8),
          legend.key = element_blank(), legend.key.size = unit(3.5, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 9, color = "black", face = "bold"))
  
  # ---- Fig 00: labeled UMAP landscape ----
  # Original all-cell atlas landscape is retained for audit.
  # A separate acute-infection main-panel version is exported for Figure 7A.
  
  LABEL_MIN_CELLS <- 30L
  
  make_umap_df <- function(plot_idx) {
    z <- data.frame(
      UMAP_1 = meta$UMAP_1[plot_idx],
      UMAP_2 = meta$UMAP_2[plot_idx],
      fine_state = factor(meta[[FINE_STATE_FIELD]][plot_idx], levels = fine_states),
      stringsAsFactors = FALSE
    )
    z[!is.na(z$fine_state), , drop = FALSE]
  }
  
  make_label_df <- function(umap_df, label_states = NULL) {
    lab <- aggregate(cbind(UMAP_1, UMAP_2) ~ fine_state, data = umap_df, FUN = median)
    label_n <- as.data.frame(table(umap_df$fine_state), stringsAsFactors = FALSE)
    colnames(label_n) <- c("fine_state", "n_cells_plot")
    lab <- merge(lab, label_n, by = "fine_state", all.x = TRUE)
    lab <- lab[lab$n_cells_plot >= LABEL_MIN_CELLS, , drop = FALSE]
    if (!is.null(label_states)) {
      lab <- lab[as.character(lab$fine_state) %in% label_states, , drop = FALSE]
    }
    lab$label <- short_state(lab$fine_state)
    lab
  }
  
  make_umap_landscape <- function(umap_df,
                                  label_df,
                                  plot_title,
                                  label_size = 5.2,
                                  point_size = 0.050,
                                  point_alpha = 0.42,
                                  base_size = 16,
                                  title_size = 24,
                                  axis_title_size = 18,
                                  axis_text_size = 15,
                                  legend_text_size = 12,
                                  legend_point_size = 4.2,
                                  with_legend = TRUE,
                                  boxed_labels = FALSE,
                                  use_repel = FALSE,
                                  label_fontface = "bold") {
    
    p <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = fine_state)) +
      geom_point(size = point_size, alpha = point_alpha, stroke = 0)
    
    if (boxed_labels) {
      p <- p +
        geom_label(
          data = label_df,
          aes(x = UMAP_1, y = UMAP_2, label = label),
          inherit.aes = FALSE,
          color = "black",
          fill = NA,
          label.size = NA,
          size = label_size,
          fontface = label_fontface
        )
    } else if (use_repel && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(x = UMAP_1, y = UMAP_2, label = label),
          inherit.aes = FALSE,
          color = "black",
          size = label_size,
          fontface = label_fontface,
          box.padding = 0.45,
          point.padding = 0.12,
          force = 1.8,
          force_pull = 0.25,
          min.segment.length = Inf,
          segment.color = NA,
          max.overlaps = Inf,
          seed = 20260803
        )
    } else {
      p <- p +
        geom_text(
          data = label_df,
          aes(x = UMAP_1, y = UMAP_2, label = label),
          inherit.aes = FALSE,
          color = "black",
          size = label_size,
          fontface = label_fontface,
          check_overlap = FALSE
        )
    }
    
    p +
      scale_color_manual(
        values = state_colors,
        labels = short_state,
        name = NULL,
        drop = FALSE
      ) +
      scale_x_continuous(expand = expansion(mult = 0.12)) +
      scale_y_continuous(expand = expansion(mult = 0.10)) +
      coord_fixed(clip = "off") +
      labs(
        x = "umap_1",
        y = "umap_2",
        title = plot_title
      ) +
      theme_classic(base_size = base_size) +
      theme(
        plot.title = element_text(
          size = title_size,
          face = "bold",
          hjust = 0.5,
          color = "black"
        ),
        axis.title = element_text(
          size = axis_title_size,
          face = "bold",
          color = "black"
        ),
        axis.text = element_text(
          size = axis_text_size,
          color = "black"
        ),
        axis.line = element_line(
          color = "black",
          linewidth = 0.65
        ),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.55
        ),
        legend.position = if (with_legend) "right" else "none",
        legend.text = element_text(
          size = legend_text_size,
          color = "black"
        ),
        legend.key = element_blank(),
        legend.key.size = unit(5.0, "mm"),
        legend.box.margin = margin(0, 0, 0, 18),
        legend.margin = margin(0, 0, 0, 4),
        plot.margin = margin(18, 36, 30, 36)
      ) +
      guides(
        color = guide_legend(
          override.aes = list(size = legend_point_size, alpha = 1),
          ncol = 1
        )
      )
  }
  
  # ---- Fig00 audit: all cells, all eligible labels, with legend -------------
  umap_df_all <- make_umap_df(seq_len(nrow(meta)))
  label_df_all <- make_label_df(umap_df_all, label_states = NULL)
  
  p00 <- make_umap_landscape(
    umap_df = umap_df_all,
    label_df = label_df_all,
    plot_title = "Cell types - GSE216009",
    label_size = 4.4,
    point_size = 0.120,
    point_alpha = 0.95,
    base_size = 16,
    title_size = 24,
    axis_title_size = 18,
    axis_text_size = 15,
    legend_text_size = 11.5,
    legend_point_size = 4.2,
    with_legend = TRUE,
    boxed_labels = FALSE,
    use_repel = TRUE,
    label_fontface = "bold"
  )
  
  save_pdf_png(p00, "fig00_umap_fine_states", width = 14.2, height = 9.8)
  
  # ---- Fig07A main-panel-ready: acute-infection subset, key labels only -----
  FIG07A_LABEL_STATES <- c(
    "HSPCs",
    "Cycling_neutrophil_progenitors",
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "IL1R2+_immature_neutrophils",
    "S100A8-9_hi_neutrophils",
    "Mature_neutrophils",
    "Degranulating_neutrophils",
    "Apoptosing_neutrophils",
    "Classical_monocytes",
    "Cycling_TNK",
    "NK cells",
    "B_cells",
    "Plasmablasts",
    "Platelets"
  )
  
  umap_df_main <- make_umap_df(primary_idx)
  label_df_main <- make_label_df(umap_df_main, label_states = FIG07A_LABEL_STATES)
  
  FIG07A_DIR <- file.path(OUT_ROOT, "fig07a_umap_mainpanel")
  dir.create(FIG07A_DIR, recursive = TRUE, showWarnings = FALSE)
  
  p00_main <- make_umap_landscape(
    umap_df = umap_df_main,
    label_df = label_df_main,
    plot_title = "Cell types",
    label_size = 4.4,
    point_size = 0.130,
    point_alpha = 0.95,
    base_size = 16,
    title_size = 24,
    axis_title_size = 18,
    axis_text_size = 15,
    legend_text_size = 11.5,
    legend_point_size = 4.2,
    with_legend = FALSE,
    boxed_labels = FALSE,
    use_repel = TRUE,
    label_fontface = "bold"
  )
  
  save_pdf_png(
    p00_main,
    "fig07a_umap_fine_states_mainpanel",
    dir = FIG07A_DIR,
    width = 7.2,
    height = 6.4
  )
  
  
  
  # ---- Fig 00b: representative gene UMAPs, per-gene panels + per-gene legends ----
  # Desired style: each gene is an independent ggplot panel with its own colorbar.
  # This improves visual clarity, but note that colors are no longer strictly comparable
  # across genes because each gene has its own robust 99.5% expression cap.
  
  expr_genes <- intersect(UMAP_GENES, genes_full)
  UMAP_EXPR_MAX_CELLS <- 60000L
  
  set.seed(20260803)
  expr_pos <- if (ncol(cnt_p) > UMAP_EXPR_MAX_CELLS) {
    sample(seq_len(ncol(cnt_p)), UMAP_EXPR_MAX_CELLS)
  } else {
    seq_len(ncol(cnt_p))
  }
  
  xlim_expr <- range(meta_p$UMAP_1[expr_pos], finite = TRUE)
  ylim_expr <- range(meta_p$UMAP_2[expr_pos], finite = TRUE)
  
  make_gene_umap_panel <- function(g) {
    gi <- match(g, genes_full)
    if (is.na(gi)) stop_msg("UMAP gene absent from genes_full: ", g)
    
    expr <- log2(as.numeric(X_p[gi, expr_pos]) / lib_cell[expr_pos] * 1e6 + 1)
    
    df <- data.frame(
      UMAP_1 = meta_p$UMAP_1[expr_pos],
      UMAP_2 = meta_p$UMAP_2[expr_pos],
      expr = expr,
      stringsAsFactors = FALSE
    )
    
    # Put high-expression cells on top
    df <- df[order(df$expr), , drop = FALSE]
    fg <- df[is.finite(df$expr) & df$expr > 0, , drop = FALSE]
    
    pos_expr <- df$expr[is.finite(df$expr) & df$expr > 0]
    gene_cap <- if (length(pos_expr)) {
      as.numeric(stats::quantile(pos_expr, 0.995, na.rm = TRUE))
    } else {
      1
    }
    gene_cap <- max(gene_cap, 1)
    
    gene_breaks <- pretty(c(0, gene_cap), n = 4)
    gene_breaks <- gene_breaks[gene_breaks >= 0 & gene_breaks <= gene_cap]
    if (length(gene_breaks) < 2L) gene_breaks <- c(0, gene_cap)
    
    ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(color = "#EDEDED", size = 0.035, alpha = 0.85, stroke = 0) +
      geom_point(
        data = fg,
        aes(color = pmin(expr, gene_cap)),
        size = 0.065,
        alpha = 0.90,
        stroke = 0
      ) +
      scale_color_gradientn(
        colors = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D"),
        limits = c(0, gene_cap),
        breaks = gene_breaks,
        oob = scales::squish,
        name = "log2 CPM"
      ) +
      coord_fixed(xlim = xlim_expr, ylim = ylim_expr, expand = FALSE) +
      labs(title = g, x = NULL, y = NULL) +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(
          size = 10.5,
          face = "bold",
          hjust = 0,
          color = "black",
          margin = margin(b = 4)
        ),
        axis.title = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.ticks = element_line(color = "grey40", linewidth = 0.35),
        axis.line = element_line(color = "grey45", linewidth = 0.45),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(6, 8, 6, 6),
        legend.position = "right",
        legend.title = element_text(size = 8.5, color = "black"),
        legend.text = element_text(size = 7.5, color = "black"),
        legend.key.width = unit(3.2, "mm"),
        legend.key.height = unit(13, "mm")
      ) +
      guides(
        color = guide_colorbar(
          title.position = "top",
          title.hjust = 0,
          barwidth = unit(3.2, "mm"),
          barheight = unit(18, "mm"),
          ticks = TRUE
        )
      )
  }
  
  p00b_list <- lapply(expr_genes, make_gene_umap_panel)
  
  p00b <- gridExtra::arrangeGrob(
    grobs = p00b_list,
    ncol = 4,
    top = NULL
  )
  
  save_pdf_png(
    p00b,
    "fig00b_umap_gene_expression",
    width = 13.5,
    height = 6.5
  )
  
  # ---- heatmap helper (bins + clamp; short x labels; not evaluable; * low-detected) ----
  make_heat <- function(mat, limits, labels, fill_values, legend_title, title) {
    gs <- rownames(mat)
    low_detected_genes <- setdiff(genes_full, genes_primary)
    gene_label_map <- setNames(ifelse(gs %in% low_detected_genes, paste0(gs, "*"), gs), gs)
    det_mod <- factor(module_label(sig01$module[match(gs, sig01$gene)]),
                      levels = c("Cell cycle","Neutrophil","Other","Inflammatory"))
    mod_list <- split(gs, det_mod, drop = TRUE)
    mods_rev <- rev(names(mod_list))
    gene_levels <- character(0)
    for (i in seq_along(mods_rev)) {
      gene_levels <- c(gene_levels, mod_list[[mods_rev[i]]])
      gene_levels <- c(gene_levels, paste0("spacer|", mods_rev[i]))
    }
    hd <- do.call(rbind, lapply(gene_levels, function(g)
      data.frame(gene = g, fine_state = fine_states, stringsAsFactors = FALSE)))
    hd$fill_level <- vapply(seq_len(nrow(hd)), function(i) {
      g <- hd$gene[i]; s <- hd$fine_state[i]
      if (grepl("^spacer", g)) "spacer" else assign_bin(mat[g, s], limits, labels)
    }, character(1))
    hd$rowtype <- "heat"
    strip_df <- data.frame(gene = "lineage", fine_state = fine_states,
                           fill_level = fine_states, rowtype = "strip", stringsAsFactors = FALSE)
    comb <- rbind(hd[, c("gene","fine_state","fill_level","rowtype")],
                  strip_df[, c("gene","fine_state","fill_level","rowtype")])
    comb$gene <- factor(comb$gene, levels = c(gene_levels, "lineage"))
    comb$fine_state <- factor(comb$fine_state, levels = fine_states)
    comb$rowtype <- factor(comb$rowtype, levels = c("heat","strip"))
    ggplot(comb, aes(x = fine_state, y = gene, fill = fill_level)) +
      geom_tile(color = "white", linewidth = 0.3) +
      facet_grid(rowtype ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(values = fill_values, breaks = c(labels, NOT_EVAL),
                        labels = c(labels, NOT_EVAL), name = legend_title, drop = FALSE) +
      scale_y_discrete(labels = function(l) {
        ifelse(grepl("^spacer\\|", l), sub("^spacer\\|", "", l),
               ifelse(l == "lineage", "",
                      ifelse(l %in% names(gene_label_map), gene_label_map[l], l)))
      }) +
      scale_x_discrete(labels = short_state) +
      labs(x = "Fine cell state", y = NULL, title = title,
           caption = "* low-detected genes; descriptive only") +
      style_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5, color = "black"),
            axis.text.y = element_text(size = 6.5), panel.spacing = unit(4, "mm"),
            strip.text.y = element_blank(), strip.background.y = element_blank(),
            legend.key.height = unit(5, "mm"), plot.margin = margin(8, 8, 20, 8))
  }
  
  lc_limits <- c(1, 2, 4, 8, 16)
  lc_labels <- c("[0,1]","(1,2]","(2,4]","(4,8]","(8,16]")
  lc_vals <- c("[0,1]" = "#F7FBFF","(1,2]" = "#C6DBEF","(2,4]" = "#6BAED6",
               "(4,8]" = "#2171B5","(8,16]" = "#08306B",
               "not evaluable" = "#E6E6E6","spacer" = "#FFFFFF")
  cf_limits <- c(0.05, 0.15, 0.30, 0.50, 1)
  cf_labels <- c("[0,0.05]","(0.05,0.15]","(0.15,0.3]","(0.3,0.5]","(0.5,1]")
  cf_vals <- c("[0,0.05]" = "#FFFFFF","(0.05,0.15]" = "#FEE5D9","(0.15,0.3]" = "#FCAE91",
               "(0.3,0.5]" = "#FB6A4A","(0.5,1]" = "#CB181D",
               "not evaluable" = "#E6E6E6","spacer" = "#FFFFFF")
  
  p02 <- make_heat(median_logcpm_mat, lc_limits, lc_labels, c(lc_vals, state_colors),
                   "median logCPM", "32-gene expression enrichment (logCPM)")
  save_pdf_png(p02, "fig02_logcpm_heatmap", width = 13.5, height = 9.2)
  p03 <- make_heat(median_cf_mat, cf_limits, cf_labels, c(cf_vals, state_colors),
                   "median contribution", "32-gene bulk-like contribution (fraction)")
  save_pdf_png(p03, "fig03_contribution_heatmap", width = 13.5, height = 9.2)
  
  # ---- Fig 04: main-source summary (primary first; * low-detected; short labels) ----
  s04 <- summary_df[
    order(!summary_df$primary_analysis_flag,
          -summary_df$main_contribution_fraction_median_all26),
  ]
  s04$state_label <- short_state(s04$main_contributor_state)
  s04$gene_label <- ifelse(s04$primary_analysis_flag, s04$signature_gene,
                           paste0(s04$signature_gene, "*"))
  s04$gene_f <- factor(s04$gene_label, levels = rev(s04$gene_label))
  s04$lineage_f <- factor(s04$main_contributor_lineage, levels = lineage_order)
  p04 <- ggplot(s04, aes(x = main_contribution_fraction_median_all26, y = gene_f)) +
    geom_col(aes(fill = lineage_f, alpha = primary_analysis_flag), width = 0.72) +
    geom_text(aes(x = pmax(main_contribution_fraction_median_all26 + 0.015, 0.025),
                  label = state_label), hjust = 0, size = 2.65) +
    scale_fill_manual(values = lineage_colors, name = "Main contributor lineage") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45), guide = "none") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    coord_cartesian(xlim = c(0, 1.18), clip = "off") +
    labs(x = "Median contribution fraction", y = NULL,
         title = "Per-gene main bulk contributor state",
         caption = "* low-detected genes; descriptive only") +
    style_theme +
    theme(axis.text.y = element_text(size = 7), plot.margin = margin(8, 135, 8, 8))
  save_pdf_png(p04, "fig04_main_source_summary", width = 11.5, height = 8.7)
  
  # ---- Fig 05: stability scatter with key-gene labels ----
  s05 <- summary_df
  s05$module_f <- factor(module_label(sig01$module[match(s05$signature_gene, sig01$gene)]),
                         levels = c("Cell cycle","Neutrophil","Other","Inflammatory"))
  label_genes <- c("RHAG","HIST1H2BM","HIST1H3B","CX3CR1","IL1B","SECTM1")
  s05$label <- ifelse(s05$signature_gene %in% label_genes | s05$source_call_confidence == "low",
                      s05$signature_gene, "")
  p05 <- ggplot(s05, aes(x = main_contribution_fraction_median_all26,
                         y = main_state_donor_detection_rate_all26,
                         color = module_f, shape = source_call_confidence)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    scale_color_manual(values = mod_pal, name = "Module") +
    scale_shape_manual(values = c(high = 16, moderate = 17, low = 4), name = "Confidence") +
    scale_x_continuous(limits = c(-0.03, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    scale_y_continuous(limits = c(0, 1.10), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(x = "Main contribution fraction (median, all26)",
         y = "Main-state donor detection (all26)",
         title = "Source stability vs contribution") +
    style_theme + theme(plot.margin = margin(8, 60, 8, 8))
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p05 <- p05 + ggrepel::geom_text_repel(aes(label = label), size = 2.8,
                                          max.overlaps = Inf, show.legend = FALSE, seed = 20260803)
  } else {
    p05 <- p05 + geom_text(aes(label = label), hjust = -0.1, vjust = 0.3,
                           size = 2.8, check_overlap = TRUE, show.legend = FALSE)
  }
  save_pdf_png(p05, "fig05_stability_scatter", width = 8.5, height = 6.5)
  
  # ---------------------------------------------------------------------------
  # 8. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  write_csv_atomic(data.frame(
    key = c("script","script_version","input_rds","script01_run_dir",
            "condition_field","acute_diag","n_acute_donors","n_primary_cells",
            "min_cells_per_donor_state","low_umi_for_contribution",
            "n_primary_genes","n_full_genes","R_version"),
    value = c("02_gene_cell_distribution.R","v6",INPUT_RDS,SCRIPT01_RUN_DIR,
              CONDITION_FIELD, paste(ACUTE_DIAG, collapse=";"),
              as.character(length(acute_donors)), as.character(length(primary_idx)),
              as.character(MIN_CELLS_PER_DONOR_STATE), as.character(LOW_UMI_FOR_CONTRIBUTION),
              as.character(length(genes_primary)), as.character(length(genes_full)),
              as.character(getRversion())),
    stringsAsFactors = FALSE), file.path(PROV_DIR, "run_provenance.csv"))
  writeLines(capture.output(sessionInfo()), file.path(PROV_DIR, "sessionInfo.txt"))
  
  cat("Script 02 v6 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("Primary cells (26 acute donors):", length(primary_idx), "\n")
  cat("Primary source-attribution genes:", nrow(summary_df), "\n")
  cat("Ambiguous source calls:", sum(summary_df$source_call_ambiguous, na.rm = TRUE), "\n")
})