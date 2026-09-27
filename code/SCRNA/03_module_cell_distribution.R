# ============================================================================
# SCRIPT 03 v16 (CLEAN WARNINGS + M2 AUDIT FLAG FIX) -- module-level cell-state attribution (M1-M4)
#   v16: based on v15; keeps warning fixes, fixes audit NA flagging, and adds
#        M2 UCell-vs-logCPM audit table for biological interpretation.
#       (1) visible plot text is ASCII-safe: >=20 cells; hyphen instead of middle dot;
#       (2) not-evaluable donor-states are pre-filtered before the M2 boxplot;
#       (3) module_neutrophil_axis_ucell_logcpm_audit.csv checks whether
#           zero/near-zero UCell scores are consistent with M2 gene logCPM.
#   Statistical analysis and figure logic otherwise unchanged from v13/v14.
# RUN: Rscript 03_module_cell_distribution_readability_v16_final_clean.R
# ============================================================================

options(stringsAsFactors = FALSE)
# ============================================================================
# GLOBAL FIGURE FONT CONTROL
# Required for Adobe Illustrator editable text.
# Recommended: Arial, installed both on Linux server and Windows/Adobe.
# ============================================================================

if (!requireNamespace("systemfonts", quietly = TRUE)) {
  stop("Package 'systemfonts' is required for font checking.", call. = FALSE)
}

PLOT_FONT_FAMILY <- Sys.getenv("FIG_FONT_FAMILY", unset = "Arial")

.font_info <- systemfonts::match_fonts(PLOT_FONT_FAMILY)
.font_path <- .font_info$path[1]

if (is.na(.font_path) || !nzchar(.font_path) || !file.exists(.font_path)) {
  stop(
    "Requested figure font is not available on this system: ",
    PLOT_FONT_FAMILY,
    "\nInstall this font on the Linux server, or run with FIG_FONT_FAMILY='Liberation Sans'.",
    call. = FALSE
  )
}

if (tolower(PLOT_FONT_FAMILY) == "arial" &&
    !grepl("arial", basename(.font_path), ignore.case = TRUE)) {
  stop(
    "FIG_FONT_FAMILY='Arial', but systemfonts matched a non-Arial file:\n",
    .font_path,
    "\nInstall real Arial on the Linux server, or explicitly use FIG_FONT_FAMILY='Liberation Sans' and install Liberation Sans in Windows/Adobe.",
    call. = FALSE
  )
}

message("Figure font: ", PLOT_FONT_FAMILY, " | ", .font_path)

fontify_plot <- function(p) {
  p + ggplot2::theme(
    text = ggplot2::element_text(family = PLOT_FONT_FAMILY)
  )
}
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
  SCRIPT02_RUN_DIR <- file.path(
    PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis",
    "02_gene_cell_distribution", "20260810_gene_cell_distribution_12"
  )
  
  OUTPUTS_ROOT <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  STAGE_NAME   <- "03_module_cell_distribution"
  RUN_PURPOSE  <- "module_cell_distribution"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  ACUTE_DIAG       <- c("Bacteraemia","Bili","CAP","CNS","IAS","IE","NF","Uro")
  MIN_CELLS_PER_DONOR_STATE <- 20L
  
  EXPECTED_N_SAMPLES        <- 48L
  EXPECTED_N_DONORS         <- 39L
  EXPECTED_N_ACUTE_DONORS   <- 26L
  EXPECTED_N_HC_DONORS      <- 6L
  EXPECTED_N_SURGERY_DONORS <- 7L
  EXPECTED_N_CONV_SAMPLES   <- 9L
  EXPECTED_N_PAIRED_DONORS  <- 9L
  EXPECTED_N_PRIMARY_CELLS  <- 151837L
  
  UMAP_MAX_CELLS   <- 80000L
  UMAP_MIN_PER_STATE <- 1000L
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS + PACKAGES (UCell required, fail-closed)
  # ---------------------------------------------------------------------------
  needed <- c(
    "SeuratObject","Matrix","ggplot2","scales","gridExtra","grid",
    "UCell","systemfonts","svglite"
  )
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!requireNamespace("UCell", quietly = TRUE))
    stop("UCell is required for Script 03. Install it before running.", call. = FALSE)
  UCELL_VERSION <- as.character(utils::packageVersion("UCell"))
  USE_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
  suppressPackageStartupMessages({
    library(SeuratObject); library(Matrix); library(ggplot2)
    library(scales); library(gridExtra); library(grid); library(UCell)
    if (USE_PATCHWORK) library(patchwork)
  })
  
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  require_cols <- function(x, cols, label) {
    m <- setdiff(cols, names(x))
    if (length(m)) stop_msg("Missing required columns in ", label, ": ", paste(m, collapse = ", "))
    invisible(TRUE)
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
  save_pdf_png <- function(p, base, dir = FIG_DIR, width = 7, height = 6) {
    p <- fontify_plot(p)
    
    pdf_file <- file.path(dir, paste0(base, ".pdf"))
    png_file <- file.path(dir, paste0(base, ".png"))
    svg_file <- file.path(dir, paste0(base, ".svg"))
    
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = PLOT_FONT_FAMILY,
      onefile = FALSE,
      bg = "white"
    )
    print(p)
    grDevices::dev.off()
    
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(
        png_file,
        width = width,
        height = height,
        units = "in",
        res = 450,
        background = "white"
      )
      print(p)
      grDevices::dev.off()
    } else {
      ggplot2::ggsave(
        png_file,
        p,
        width = width,
        height = height,
        dpi = 450,
        bg = "white",
        limitsize = FALSE
      )
    }
    
    svglite::svglite(
      svg_file,
      width = width,
      height = height,
      bg = "white",
      system_fonts = list(
        sans = PLOT_FONT_FAMILY,
        serif = PLOT_FONT_FAMILY,
        mono = PLOT_FONT_FAMILY
      )
    )
    print(p)
    grDevices::dev.off()
    
    invisible(TRUE)
  }
  save_grob_pdf_png <- function(grob, base, dir = FIG_DIR, width = 7, height = 6) {
    pdf_file <- file.path(dir, paste0(base, ".pdf"))
    png_file <- file.path(dir, paste0(base, ".png"))
    svg_file <- file.path(dir, paste0(base, ".svg"))
    
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = PLOT_FONT_FAMILY,
      onefile = FALSE,
      bg = "white"
    )
    grid::grid.newpage()
    grid::grid.draw(grob)
    grDevices::dev.off()
    
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(
        png_file,
        width = width,
        height = height,
        units = "in",
        res = 450,
        background = "white"
      )
    } else {
      grDevices::png(
        png_file,
        width = width,
        height = height,
        units = "in",
        res = 450,
        bg = "white"
      )
    }
    grid::grid.newpage()
    grid::grid.draw(grob)
    grDevices::dev.off()
    
    svglite::svglite(
      svg_file,
      width = width,
      height = height,
      bg = "white",
      system_fonts = list(
        sans = PLOT_FONT_FAMILY,
        serif = PLOT_FONT_FAMILY,
        mono = PLOT_FONT_FAMILY
      )
    )
    grid::grid.newpage()
    grid::grid.draw(grob)
    grDevices::dev.off()
    
    invisible(TRUE)
  }
  save_any_pdf_png <- function(p, base, dir = FIG_DIR, width = 7, height = 6) {
    if (inherits(p, c("grob", "gtable", "gTree"))) {
      save_grob_pdf_png(p, base, dir = dir, width = width, height = height)
    } else {
      save_pdf_png(p, base, dir = dir, width = width, height = height)
    }
    invisible(TRUE)
  }
  EXPORT_FIG8_SOURCE_PANELS <- TRUE
  EXPORT_FIG8_COMPOSITES    <- FALSE
  
  save_fig8_source_panel <- function(p, base, width = 7, height = 6,
                                     dpi_png = 450, dpi_tiff = 600) {
    fig8_dir <- file.path(dirname(dirname(OUT_ROOT)), "Figure8_source_panels")
    dir.create(fig8_dir, recursive = TRUE, showWarnings = FALSE)
    
    pdf_file  <- file.path(fig8_dir, paste0(base, ".pdf"))
    png_file  <- file.path(fig8_dir, paste0(base, ".png"))
    tiff_file <- file.path(fig8_dir, paste0(base, ".tiff"))
    svg_file  <- file.path(fig8_dir, paste0(base, ".svg"))
    
    p <- fontify_plot(p)
    
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = PLOT_FONT_FAMILY,
      onefile = FALSE,
      bg = "white"
    )
    print(p)
    grDevices::dev.off()
    
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(
        png_file,
        width = width,
        height = height,
        units = "in",
        res = dpi_png,
        background = "white"
      )
      print(p)
      grDevices::dev.off()
      
      ragg::agg_tiff(
        tiff_file,
        width = width,
        height = height,
        units = "in",
        res = dpi_tiff,
        compression = "lzw",
        background = "white"
      )
      print(p)
      grDevices::dev.off()
    } else {
      grDevices::png(
        png_file,
        width = width,
        height = height,
        units = "in",
        res = dpi_png,
        bg = "white"
      )
      print(p)
      grDevices::dev.off()
      
      grDevices::tiff(
        tiff_file,
        width = width,
        height = height,
        units = "in",
        res = dpi_tiff,
        compression = "lzw",
        bg = "white"
      )
      print(p)
      grDevices::dev.off()
    }
    
    svglite::svglite(
      svg_file,
      width = width,
      height = height,
      bg = "white",
      system_fonts = list(
        sans = PLOT_FONT_FAMILY,
        serif = PLOT_FONT_FAMILY,
        mono = PLOT_FONT_FAMILY
      )
    )
    print(p)
    grDevices::dev.off()
    
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
      "cDCs" = "cDCs", "pDCs" = "pDCs",
      "Naive_CD4_T_cells" = "Naive CD4+ T",
      "Memory_CD4_T_cells" = "Memory CD4+ T",
      "Naive_CD8_T_cells" = "Naive CD8+ T",
      "CD8_T_cells" = "CD8+ T",
      "Cycling_TNK" = "Cycling T/NK",
      "NK cells" = "NK cells", "B_cells" = "B cells",
      "Plasmablasts" = "Plasmablasts", "Platelets" = "Platelets"
    )
    ifelse(z %in% names(map), unname(map[z]), z)
  }
  
  wrap_caption <- function(x, width = 120) {
    paste(strwrap(x, width = width), collapse = "\n")
  }
  safe_median <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else median(x) }
  safe_mean   <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
  safe_max    <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }
  row_sum_preserve_all_na <- function(m) {
    mm <- as.matrix(m)
    finite <- is.finite(mm)
    out <- rowSums(ifelse(finite, mm, 0))
    out[rowSums(finite) == 0] <- NA_real_
    out
  }
  matrix_to_long <- function(mat, value_col) {
    out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(out) <- c("donor", "fine_state", value_col)
    out
  }
  pick_ucell_col <- function(scores, nm) {
    cand <- c(paste0(nm, "_UCell"), nm)
    hit <- cand[cand %in% colnames(scores)]
    if (!length(hit)) stop_msg("Cannot find UCell score column for ", nm,
                               ". Available: ", paste(colnames(scores), collapse = ", "))
    hit[1]
  }
  stratified_umap_sample <- function(meta_df, group_col, max_total = 40000L,
                                     min_per_group = 500L, seed = 20260803) {
    set.seed(seed)
    groups <- split(seq_len(nrow(meta_df)), meta_df[[group_col]])
    keep <- unlist(lapply(groups, function(idx) {
      if (length(idx) <= min_per_group) idx else sample(idx, min_per_group)
    }), use.names = FALSE)
    remaining <- setdiff(seq_len(nrow(meta_df)), keep)
    n_extra <- max_total - length(keep)
    if (n_extra > 0L && length(remaining) > 0L)
      keep <- c(keep, sample(remaining, min(n_extra, length(remaining))))
    sort(unique(keep))
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  PROV_DIR  <- file.path(OUT_ROOT, "provenance")
  for (d in c(TABLE_DIR, FIG_DIR, PROV_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. READ INPUTS
  # ---------------------------------------------------------------------------
  s01_t <- c("frozen_signature_input.csv","final_gene_evaluability_32.csv",
             "fine_state_display_order_audit.csv","gene_mapping_audit.csv")
  missing_t <- s01_t[!file.exists(file.path(SCRIPT01_RUN_DIR, "tables", s01_t))]
  if (length(missing_t)) stop_msg("Script-01 run missing tables: ", paste(missing_t, collapse="; "),
                                  "\nFix SCRIPT01_RUN_DIR.")
  s02_t <- c("per_gene_all_state_rankings.csv","per_gene_cell_source_summary.csv")
  missing_2 <- s02_t[!file.exists(file.path(SCRIPT02_RUN_DIR, "tables", s02_t))]
  if (length(missing_2)) stop_msg("Script-02 run missing tables: ", paste(missing_2, collapse="; "),
                                  "\nFix SCRIPT02_RUN_DIR.")
  
  sig01   <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "frozen_signature_input.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  eva01   <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "final_gene_evaluability_32.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  state01 <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "fine_state_display_order_audit.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  map01   <- read.csv(file.path(SCRIPT01_RUN_DIR, "tables", "gene_mapping_audit.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  rank02  <- read.csv(file.path(SCRIPT02_RUN_DIR, "tables", "per_gene_all_state_rankings.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  sum02   <- read.csv(file.path(SCRIPT02_RUN_DIR, "tables", "per_gene_cell_source_summary.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  
  require_cols(sig01, c("gene","module","direction"), "frozen_signature_input.csv")
  require_cols(eva01, c("signature_gene","class"), "final_gene_evaluability_32.csv")
  require_cols(state01, c("display_order","fine_state","lineage","color"), "fine_state_display_order_audit.csv")
  require_cols(map01, c("signature_gene","dataset_feature"), "gene_mapping_audit.csv")
  require_cols(rank02, c("signature_gene","fine_state","median_logCPM"), "per_gene_all_state_rankings.csv")
  require_cols(sum02, c("signature_gene","main_contributor_state","main_enriched_state",
                        "source_call_confidence","source_call_ambiguous","evaluable_class","module"),
               "per_gene_cell_source_summary.csv")
  
  stopifnot(nrow(sig01) == 32L, nrow(state01) == 24L)
  fine_states   <- state01$fine_state[order(state01$display_order)]
  state_lineage <- setNames(state01$lineage, state01$fine_state)
  state_colors  <- setNames(state01$color, state01$fine_state)
  lineage_order <- unique(state01$lineage[order(state01$display_order)])
  
  genes_full   <- sig01$gene
  gene_class   <- setNames(eva01$class, eva01$signature_gene)
  genes_primary<- eva01$signature_gene[eva01$class == "stable_detected"]
  feat_map     <- setNames(map01$dataset_feature, map01$signature_gene)
  
  m1_genes <- sig01$gene[module_label(sig01$module) == "Cell cycle"]
  m2_genes <- sig01$gene[module_label(sig01$module) == "Neutrophil"]
  m3_genes <- sig01$gene[module_label(sig01$module) == "Other"]
  m4_genes <- sig01$gene[module_label(sig01$module) == "Inflammatory"]
  stopifnot(length(m1_genes) == 11L, length(m2_genes) == 10L,
            length(m3_genes) == 8L, length(m4_genes) == 3L)
  
  M1_SOURCE_CLASS <- list(
    progenitor_granulopoiesis = c("HSPCs","Cycling_neutrophil_progenitors",
                                  "MPO+_immature_neutrophils_or_progenitors"),
    immature_neutrophil = c("PADI4+_immature_neutrophils","IL1R2+_immature_neutrophils",
                            "S100A8-9_hi_neutrophils"),
    cycling_lymphoid = c("Cycling_TNK"),
    other = setdiff(fine_states, c("HSPCs","Cycling_neutrophil_progenitors",
                                   "MPO+_immature_neutrophils_or_progenitors",
                                   "PADI4+_immature_neutrophils","IL1R2+_immature_neutrophils",
                                   "S100A8-9_hi_neutrophils","Cycling_TNK"))
  )
  M2_AXIS <- c("Cycling_neutrophil_progenitors","MPO+_immature_neutrophils_or_progenitors",
               "PADI4+_immature_neutrophils","IL1R2+_immature_neutrophils",
               "S100A8-9_hi_neutrophils","Mature_neutrophils",
               "Degranulating_neutrophils","Apoptosing_neutrophils")
  
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
  # 3. COHORT (fail-closed)
  # ---------------------------------------------------------------------------
  require_cols(meta, c(SAMPLE_FIELD, CONDITION_FIELD, FINE_STATE_FIELD), "obj@meta.data")
  meta$sample_id <- as.character(meta[[SAMPLE_FIELD]])
  meta$diagnosis_for_script03 <- as.character(meta[[CONDITION_FIELD]])
  sample_diag <- unique(meta[, c("sample_id","diagnosis_for_script03"), drop = FALSE])
  n_diag_per_sample <- tapply(sample_diag$diagnosis_for_script03, sample_diag$sample_id,
                              function(z) length(unique(z)))
  if (any(n_diag_per_sample != 1L)) stop_msg("sample_id with multiple diagnosis: ",
                                             paste(names(n_diag_per_sample)[n_diag_per_sample != 1L], collapse = ", "))
  sample_tab <- sample_diag
  sample_tab$donor_id <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort_group <- NA_character_
  sample_tab$cohort_group[sample_tab$diagnosis_for_script03 %in% ACUTE_DIAG] <- "Acute_sepsis"
  sample_tab$cohort_group[sample_tab$diagnosis_for_script03 == "Conv"] <- "Convalescent"
  sample_tab$cohort_group[sample_tab$diagnosis_for_script03 == "HV"] <- "Healthy_control"
  sample_tab$cohort_group[sample_tab$diagnosis_for_script03 == "CS"] <- "Surgery_control"
  if (anyNA(sample_tab$cohort_group)) stop_msg("Unmapped diagnosis: ",
                                               
                                               
                                               paste(unique(sample_tab$diagnosis_for_script03[is.na(sample_tab$cohort_group)]), collapse = ", "))
  is_conv_suffix <- grepl("_CONV$", sample_tab$sample_id)
  is_conv_diag   <- sample_tab$diagnosis_for_script03 == "Conv"
  stopifnot(identical(is_conv_suffix, is_conv_diag))
  n_samples <- length(unique(sample_tab$sample_id)); n_donors <- length(unique(sample_tab$donor_id))
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
  stopifnot(length(acute_donors) == EXPECTED_N_ACUTE_DONORS,
            length(hc_donors) == EXPECTED_N_HC_DONORS,
            length(surgery_donors) == EXPECTED_N_SURGERY_DONORS,
            nrow(conv_samples) == EXPECTED_N_CONV_SAMPLES,
            length(paired_donors) == EXPECTED_N_PAIRED_DONORS,
            all(conv_donors %in% acute_donors))
  m <- match(meta$sample_id, sample_tab$sample_id)
  if (anyNA(m)) stop_msg("meta$sample_id absent from sample_tab.")
  meta$donor_id <- sample_tab$donor_id[m]
  meta$cohort_group <- sample_tab$cohort_group[m]
  meta$is_primary_acute <- meta$cohort_group == "Acute_sepsis"
  stopifnot(all(!meta$is_primary_acute[grepl("_CONV$", meta$sample_id)]))
  primary_idx <- which(meta$is_primary_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  
  # ---------------------------------------------------------------------------
  # 4. PRIMARY CELLS + DESIGN MATRICES
  # ---------------------------------------------------------------------------
  meta_p <- meta[primary_idx, , drop = FALSE]
  cnt_p  <- cnt[, primary_idx, drop = FALSE]
  donor_p <- factor(meta_p$donor_id, levels = acute_donors)
  state_p <- factor(meta_p[[FINE_STATE_FIELD]], levels = fine_states)
  stopifnot(!anyNA(donor_p), !anyNA(state_p))
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
  idx <- match(feat_map[genes_full], rownames(cnt))
  stopifnot(!anyNA(idx), !anyDuplicated(idx))
  X_p <- cnt_p[idx, , drop = FALSE]; rownames(X_p) <- genes_full
  C_flat <- as.matrix(X_p %*% Matrix::t(Ds))
  stopifnot(nrow(C_flat) == length(genes_full), ncol(C_flat) == n_donor * n_state)
  
  # ---------------------------------------------------------------------------
  # 5. NORMALIZED EXPRESSION MATRIX for UCell
  # ---------------------------------------------------------------------------
  m1_features <- unique(unname(feat_map[m1_genes]))
  m2_features <- unique(unname(feat_map[m2_genes]))
  if (anyNA(m1_features) || anyNA(m2_features)) stop_msg("Missing mapped features for M1/M2.")
  data_layer <- NULL
  if (inherits(obj[[COUNT_ASSAY]], "Assay5")) {
    if ("data" %in% Layers(obj[[COUNT_ASSAY]])) {
      data_layer <- LayerData(obj, assay = COUNT_ASSAY, layer = "data", fast = FALSE)
    }
  } else {
    data_layer <- tryCatch(GetAssayData(object = obj, assay = COUNT_ASSAY, slot = "data"),
                           error = function(e) NULL)
  }
  data_ok <- !is.null(data_layer) &&
    (inherits(data_layer, "Matrix") || is.matrix(data_layer)) &&
    nrow(data_layer) > 0 && ncol(data_layer) == ncol(obj) &&
    all(c(m1_features, m2_features) %in% rownames(data_layer))
  if (data_ok) {
    expr_mat <- data_layer[, primary_idx, drop = FALSE]
    norm_note <- "RNA data layer"
  } else {
    if (any(!is.finite(lib_cell) | lib_cell <= 0)) stop_msg("Non-positive library size.")
    expr_mat <- cnt_p %*% Matrix::Diagonal(n = length(lib_cell), x = 1e6 / lib_cell)
    expr_mat@x <- log1p(expr_mat@x)
    norm_note <- "computed sparse log1p(CPM) from counts"
  }
  stopifnot(identical(colnames(expr_mat), rownames(meta_p)))
  missing_ucell_features <- setdiff(c(m1_features, m2_features), rownames(expr_mat))
  if (length(missing_ucell_features)) stop_msg("UCell features absent: ",
                                               paste(missing_ucell_features, collapse = ", "))
  cat("UCell input:", norm_note, "\n")
  
  # ---------------------------------------------------------------------------
  # 6. UCell scores + donor-state MEDIAN
  # ---------------------------------------------------------------------------
  ucell_scores <- UCell::ScoreSignatures_UCell(
    expr_mat, features = list(M1 = m1_features, M2 = m2_features))
  ucell_scores <- as.data.frame(ucell_scores, check.names = FALSE)
  if (!is.null(rownames(ucell_scores)) && all(colnames(expr_mat) %in% rownames(ucell_scores))) {
    ucell_scores <- ucell_scores[colnames(expr_mat), , drop = FALSE]
  }
  stopifnot(nrow(ucell_scores) == ncol(expr_mat))
  m1_col <- pick_ucell_col(ucell_scores, "M1"); m2_col <- pick_ucell_col(ucell_scores, "M2")
  m1_ucell_cell <- as.numeric(ucell_scores[[m1_col]])
  m2_ucell_cell <- as.numeric(ucell_scores[[m2_col]])
  names(m1_ucell_cell) <- colnames(expr_mat); names(m2_ucell_cell) <- colnames(expr_mat)
  stopifnot(length(m1_ucell_cell) == length(primary_idx),
            length(m2_ucell_cell) == length(primary_idx),
            all(is.finite(m1_ucell_cell)), all(is.finite(m2_ucell_cell)))
  agg_ds_median <- function(x) {
    tmp <- tapply(x, ds_idx, safe_median)
    v <- rep(NA_real_, n_donor * n_state)
    v[as.integer(names(tmp))] <- as.numeric(tmp)
    m <- matrix(v, nrow = n_donor, ncol = n_state, byrow = TRUE)
    rownames(m) <- levels(donor_p); colnames(m) <- fine_states
    m[!eligible_ds] <- NA_real_
    m
  }
  m1_ucell_ds <- agg_ds_median(m1_ucell_cell)
  m2_ucell_ds <- agg_ds_median(m2_ucell_cell)
  ucell_long <- rbind(
    transform(matrix_to_long(m1_ucell_ds, "ucell_median"), module = "M1"),
    transform(matrix_to_long(m2_ucell_ds, "ucell_median"), module = "M2"))
  ucell_long <- ucell_long[, c("module","donor","fine_state","ucell_median")]
  ucell_long$lineage <- unname(state_lineage[ucell_long$fine_state])
  ucell_long <- ucell_long[!is.na(ucell_long$ucell_median), , drop = FALSE]
  write_csv_atomic(ucell_long, file.path(TABLE_DIR, "module_ucell_donor_state_median.csv"))
  
  # ---------------------------------------------------------------------------
  # 7. MODULE COMPOSITE + MODULE CONTRIBUTION FRACTION
  # ---------------------------------------------------------------------------
  composite_list <- list(); module_cf_list <- list()
  for (mod in c("M1","M2")) {
    gs <- if (mod == "M1") m1_genes else m2_genes
    gidx <- which(genes_full %in% gs)
    mod_umi_ds_vec <- colSums(C_flat[gidx, , drop = FALSE])
    mod_umi_ds <- matrix(mod_umi_ds_vec, nrow = n_donor, ncol = n_state, byrow = TRUE)
    rownames(mod_umi_ds) <- levels(donor_p); colnames(mod_umi_ds) <- fine_states
    comp <- log2(mod_umi_ds_vec / lib_ds_flat * 1e6 + 1)
    comp_m <- matrix(comp, nrow = n_donor, ncol = n_state, byrow = TRUE)
    rownames(comp_m) <- levels(donor_p); colnames(comp_m) <- fine_states
    comp_m[!eligible_ds] <- NA_real_
    composite_list[[mod]] <- comp_m
    mod_total_d <- rowSums(mod_umi_ds)
    cf <- sweep(mod_umi_ds, 1, mod_total_d, "/")
    cf[mod_total_d == 0, ] <- NA_real_
    module_cf_list[[mod]] <- cf
  }
  composite_long <- rbind(
    transform(matrix_to_long(composite_list[["M1"]], "composite_logCPM"), module = "M1"),
    transform(matrix_to_long(composite_list[["M2"]], "composite_logCPM"), module = "M2"))
  composite_long <- composite_long[, c("module","donor","fine_state","composite_logCPM")]
  composite_long <- composite_long[!is.na(composite_long$composite_logCPM), , drop = FALSE]
  write_csv_atomic(composite_long, file.path(TABLE_DIR, "module_composite_donor_state.csv"))
  
  module_cf_long <- rbind(
    transform(matrix_to_long(module_cf_list[["M1"]], "module_contribution_fraction"), module = "M1"),
    transform(matrix_to_long(module_cf_list[["M2"]], "module_contribution_fraction"), module = "M2"))
  module_cf_long <- module_cf_long[, c("module","donor","fine_state","module_contribution_fraction")]
  module_cf_long$lineage <- unname(state_lineage[module_cf_long$fine_state])
  write_csv_atomic(module_cf_long, file.path(TABLE_DIR, "module_contribution_fraction_donor_state.csv"))
  
  # ---------------------------------------------------------------------------
  # 8. M1 SOURCE CLASSIFICATION (contribution-driven; all-NA preserved)
  # ---------------------------------------------------------------------------
  state_ucell_m1 <- apply(m1_ucell_ds, 2, safe_median)
  state_comp_m1  <- apply(composite_list[["M1"]], 2, safe_median)
  cellfrac_ds <- sweep(n_cells_ds, 1, rowSums(n_cells_ds), "/")
  
  class_rows <- lapply(names(M1_SOURCE_CLASS), function(cl) {
    sts <- intersect(M1_SOURCE_CLASS[[cl]], fine_states)
    ucell_sub <- m1_ucell_ds[, sts, drop = FALSE]
    comp_sub  <- composite_list[["M1"]][, sts, drop = FALSE]
    cf_sub    <- module_cf_list[["M1"]][, sts, drop = FALSE]
    donor_ucell_mean <- apply(ucell_sub, 1, safe_mean)
    donor_ucell_max  <- apply(ucell_sub, 1, safe_max)
    donor_comp_mean  <- apply(comp_sub, 1, safe_mean)
    donor_comp_max   <- apply(comp_sub, 1, safe_max)
    donor_cf_sum     <- row_sum_preserve_all_na(cf_sub)
    donor_cell_fraction <- rowSums(cellfrac_ds[, sts, drop = FALSE], na.rm = TRUE)
    donor_has_eligible_state <- rowSums(eligible_ds[, sts, drop = FALSE]) > 0
    data.frame(
      source_class = cl, n_states = length(sts), states = paste(sts, collapse = ";"),
      ucell_state_median_mean = safe_mean(state_ucell_m1[sts]),
      ucell_state_median_max = safe_max(state_ucell_m1[sts]),
      ucell_donor_mean_median = safe_median(donor_ucell_mean),
      ucell_donor_max_median = safe_median(donor_ucell_max),
      composite_state_median_mean = safe_mean(state_comp_m1[sts]),
      composite_state_median_max = safe_max(state_comp_m1[sts]),
      composite_donor_mean_median = safe_median(donor_comp_mean),
      composite_donor_max_median = safe_median(donor_comp_max),
      contribution_fraction_median = safe_median(donor_cf_sum),
      contribution_fraction_mean = safe_mean(donor_cf_sum),
      cell_fraction_median = safe_median(donor_cell_fraction),
      cell_fraction_mean = safe_mean(donor_cell_fraction),
      eligible_donor_n = sum(donor_has_eligible_state),
      eligible_donor_rate = mean(donor_has_eligible_state),
      stringsAsFactors = FALSE)
  })
  m1_source <- do.call(rbind, class_rows)
  m1_source$enrichment_ratio <- m1_source$contribution_fraction_median / m1_source$cell_fraction_median
  m1_source$enrichment_ratio[!is.finite(m1_source$enrichment_ratio)] <- NA_real_
  m1_source <- m1_source[order(-m1_source$contribution_fraction_median,
                               -m1_source$ucell_donor_max_median,
                               -m1_source$composite_donor_max_median), ]
  write_csv_atomic(m1_source, file.path(TABLE_DIR, "module_ucell_source_classification.csv"))
  
  # ---------------------------------------------------------------------------
  # 9. M2 NEUTROPHIL MATURATION AXIS
  # ---------------------------------------------------------------------------
  axis_states <- M2_AXIS[M2_AXIS %in% fine_states]
  m2_axis_ucell <- vapply(axis_states, function(s) safe_median(m2_ucell_ds[, s]), numeric(1))
  m2_axis_comp  <- vapply(axis_states, function(s) safe_median(composite_list[["M2"]][, s]), numeric(1))
  g_logcpm_axis <- rank02[rank02$signature_gene %in% m2_genes & rank02$fine_state %in% axis_states,
                          c("signature_gene","fine_state","median_logCPM"), drop = FALSE]
  write_csv_atomic(
    data.frame(fine_state = axis_states, m2_ucell_median = unname(m2_axis_ucell),
               m2_composite_median = unname(m2_axis_comp), stringsAsFactors = FALSE),
    file.path(TABLE_DIR, "module_neutrophil_axis_ucell.csv"))
  write_csv_atomic(g_logcpm_axis, file.path(TABLE_DIR, "module_neutrophil_axis_gene_logcpm.csv"))
  
  # M2 audit table: cross-check rank-based UCell against expression-based logCPM.
  # This is interpretation-critical when later/mature neutrophil states have
  # zero/near-zero UCell despite possible nonzero granule-gene expression.
  m2_axis_logcpm_state <- do.call(rbind, lapply(axis_states, function(st) {
    z <- g_logcpm_axis[g_logcpm_axis$fine_state == st, , drop = FALSE]
    data.frame(
      fine_state = st,
      m2_gene_logcpm_median = safe_median(z$median_logCPM),
      m2_gene_logcpm_mean = safe_mean(z$median_logCPM),
      m2_gene_logcpm_max = safe_max(z$median_logCPM),
      any_gene_logCPM_ge2 = any(is.finite(z$median_logCPM) & z$median_logCPM >= 2),
      m2_gene_n = nrow(z),
      stringsAsFactors = FALSE
    )
  }))
  m2_axis_audit <- merge(
    data.frame(
      fine_state = axis_states,
      m2_ucell_median = unname(m2_axis_ucell),
      m2_composite_median = unname(m2_axis_comp),
      stringsAsFactors = FALSE
    ),
    m2_axis_logcpm_state, by = "fine_state", all.x = TRUE, sort = FALSE
  )
  m2_axis_audit$fine_state <- factor(m2_axis_audit$fine_state, levels = axis_states)
  m2_axis_audit <- m2_axis_audit[order(m2_axis_audit$fine_state), , drop = FALSE]
  m2_axis_audit$fine_state <- as.character(m2_axis_audit$fine_state)
  # Flag logic:
  # - UCell is rank-based and may be zero even when absolute module-gene logCPM is nonzero.
  # - log2(CPM + 1) >= 2 corresponds to CPM >= 3; this is a deliberately permissive
  #   audit threshold used to flag expression present despite zero UCell, not a high-expression cutoff.
  # - not_evaluable is handled first; otherwise NA UCell would be incorrectly labelled positive.
  m2_axis_audit$interpretation_flag <- ifelse(
    !is.finite(m2_axis_audit$m2_ucell_median),
    "not_evaluable",
    ifelse(
      m2_axis_audit$m2_ucell_median > 0,
      "UCell_positive",
      ifelse(
        is.finite(m2_axis_audit$m2_gene_logcpm_median) &
          m2_axis_audit$m2_gene_logcpm_median >= 2,
        "UCell_zero_but_logCPM_present",
        "UCell_zero_or_near_zero"
      )
    )
  )
  write_csv_atomic(
    m2_axis_audit,
    file.path(TABLE_DIR, "module_neutrophil_axis_ucell_logcpm_audit.csv")
  )
  
  # ---------------------------------------------------------------------------
  # 10. WITHIN-MODULE HETEROGENEITY + M3/M4
  # ---------------------------------------------------------------------------
  het <- sum02[, c("signature_gene","module","main_contributor_state","main_enriched_state",
                   "source_call_confidence","source_call_ambiguous","evaluable_class"), drop = FALSE]
  het$module_of <- module_label(het$module)
  het$contributor_equals_enriched <- het$main_contributor_state == het$main_enriched_state
  write_csv_atomic(het, file.path(TABLE_DIR, "within_module_gene_heterogeneity.csv"))
  m3m4_from02 <- sum02[module_label(sum02$module) %in% c("Other","Inflammatory"), , drop = FALSE]
  write_csv_atomic(m3m4_from02, file.path(TABLE_DIR, "m3_m4_gene_level_from_script02.csv"))
  
  # ---------------------------------------------------------------------------
  # 11. STYLE + FIGURES (v16: v14/v15 layout + clean warnings + M2 audit)
  #     Figure encoding unchanged from v14; added M2 interpretation audit table.
  # ---------------------------------------------------------------------------
  
  # ---- colour systems --------------------------------------------------------
  # M1 (cell-cycle) uses a cool blue ramp; M2 (neutrophil) a warm orange ramp.
  # The two modules deliberately carry distinct scales (localisation contrast),
  # so no design meta-notes are needed inside the panels.
  ucell_m1_pal <- c("#F7FBFF", "#D6EAF8", "#9DCAE6", "#4F9ED1", "#1764AB", "#08306B")
  ucell_m2_pal <- c("#FFF7EC", "#FDD49E", "#FDBB84", "#FC8D59", "#D7301F", "#7F0000")
  m2_logcpm_pal <- c("#FFF7EC", "#FDD49E", "#FDBB84", "#EF6548", "#B2182B", "#67001F")
  
  # Categorical palette for the M1 competing-source classification (fig05):
  # four hue families with separated lightness -> deutan/protan distinguishable.
  m1_source_pal <- c(
    "progenitor_granulopoiesis" = "#3C5488",  # deep blue
    "immature_neutrophil"       = "#E64B35",  # vermilion
    "cycling_lymphoid"          = "#7E6148",  # brown-purple
    "other"                     = "#7F7F7F"   # neutral residual class
  )
  
  NA_GREY   <- "#D9D9D9"          # not evaluable / absent
  PANEL_BG  <- "#FFFFFF"
  BG_POINT  <- "#BFC7CE"          # non-highlighted background on UMAPs
  
  # Agreement column (fig06): hue + explicit text, colour-safe.
  agree_fill <- c("TRUE" = "#DCEFE0", "FALSE" = "#FBE8E0")
  agree_text <- c("TRUE" = "#1B6B3F", "FALSE" = "#A63F0E")
  
  theme_paper <- theme_classic(base_size = 11) +
    theme(
      axis.text   = element_text(color = "grey15", size = 8),
      axis.title  = element_text(color = "grey10", size = 9),
      axis.line   = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks  = element_line(color = "grey30", linewidth = 0.3),
      panel.grid  = element_blank(),
      plot.title    = element_text(size = 11, face = "bold", hjust = 0,
                                   color = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 8.5, color = "grey35", margin = margin(b = 6)),
      plot.caption  = element_text(size = 7.2, color = "grey40", hjust = 0,
                                   lineheight = 1.05, margin = margin(t = 5)),
      plot.margin = margin(7, 9, 7, 7),
      legend.position = "right",
      legend.title = element_text(size = 8.5, color = "grey10"),
      legend.text  = element_text(size = 7.5, color = "grey15"),
      legend.key   = element_blank(),
      strip.background = element_blank(),
      strip.text   = element_text(size = 8, face = "bold", color = "grey15")
    )
  
  # Optional patchwork layout. The script does not install packages at runtime;
  # if patchwork is absent, gridExtra is used and the analysis remains runnable.
  combine_panels <- function(plots, ncol = 2, widths = NULL, heights = NULL,
                             guides = "keep", tag = TRUE) {
    if (USE_PATCHWORK) {
      p <- Reduce(`+`, plots)
      layout_args <- list(ncol = ncol, guides = guides)
      if (!is.null(widths))  layout_args$widths  <- widths
      if (!is.null(heights)) layout_args$heights <- heights
      p <- p + do.call(patchwork::plot_layout, layout_args)
      if (isTRUE(tag)) p <- p + patchwork::plot_annotation(tag_levels = "A")
      p & theme(plot.tag = element_text(face = "bold", size = 12, color = "grey10"))
    } else {
      gridExtra::arrangeGrob(grobs = plots, ncol = ncol,
                             widths = widths, heights = heights)
    }
  }
  
  # ---- Fig 01: combined UMAP comparison (distinct module scales) ------------
  top_score_idx <- function(score_vec, n = 8000L) {
    ok <- which(is.finite(score_vec) & score_vec > 0)
    if (!length(ok)) return(integer(0))
    ok[order(score_vec[ok], decreasing = TRUE)][seq_len(min(n, length(ok)))]
  }
  
  umap_base <- stratified_umap_sample(
    meta_p,
    FINE_STATE_FIELD,
    max_total = UMAP_MAX_CELLS,
    min_per_group = UMAP_MIN_PER_STATE
  )
  
  # Retain high-score cells in addition to the stratified background sample.
  # v10 uses top 50% of positive UCell scores for display, rather than top 25%,
  # because top-quartile-only UMAPs looked like white boards for sparse modules.
  umap_pos <- sort(unique(c(
    umap_base,
    top_score_idx(m1_ucell_cell, 12000L),
    top_score_idx(m2_ucell_cell, 12000L)
  )))
  
  umap_xlim <- range(meta_p$UMAP_1[umap_pos], finite = TRUE)
  umap_ylim <- range(meta_p$UMAP_2[umap_pos], finite = TRUE)
  
  # Manual labels: avoid low-coverage states being selected only because of a
  # high median in very few cells/donors. These labels are explanatory, not used
  # for source calling.
  m1_label_states <- c(
    "MPO+_immature_neutrophils_or_progenitors",
    "Cycling_neutrophil_progenitors",
    "Cycling_TNK",
    "Plasmablasts"
  )
  m2_label_states <- c(
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "Degranulating_neutrophils"
  )
  
  umap_focus <- function(score_vec, title, pal, label_states,
                         show_quantile = 0.50) {
    df <- data.frame(
      UMAP_1 = meta_p$UMAP_1[umap_pos],
      UMAP_2 = meta_p$UMAP_2[umap_pos],
      score = as.numeric(score_vec[umap_pos]),
      fine_state = as.character(meta_p[[FINE_STATE_FIELD]][umap_pos]),
      stringsAsFactors = FALSE
    )
    df$score[!is.finite(df$score)] <- 0
    pos <- df$score[df$score > 0]
    
    if (length(pos)) {
      show_cut <- as.numeric(stats::quantile(pos, show_quantile, na.rm = TRUE))
      score_cap <- as.numeric(stats::quantile(pos, 0.995, na.rm = TRUE))
    } else {
      show_cut <- Inf
      score_cap <- 1
    }
    if (!is.finite(score_cap) || score_cap <= 0) score_cap <- 1
    
    hi <- df[df$score >= show_cut & df$score > 0, , drop = FALSE]
    hi$score_plot <- pmin(hi$score, score_cap)
    hi <- hi[order(hi$score_plot), , drop = FALSE]
    
    label_states <- intersect(label_states, unique(df$fine_state))
    lab_df <- do.call(rbind, lapply(label_states, function(st) {
      z <- hi[hi$fine_state == st, c("UMAP_1", "UMAP_2"), drop = FALSE]
      if (!nrow(z)) z <- df[df$fine_state == st, c("UMAP_1", "UMAP_2"), drop = FALSE]
      if (!nrow(z)) return(NULL)
      data.frame(
        UMAP_1 = stats::median(z$UMAP_1, na.rm = TRUE),
        UMAP_2 = stats::median(z$UMAP_2, na.rm = TRUE),
        state_label = short_state(st),
        stringsAsFactors = FALSE
      )
    }))
    if (is.null(lab_df)) {
      lab_df <- data.frame(UMAP_1 = numeric(), UMAP_2 = numeric(),
                           state_label = character())
    }
    
    ggplot() +
      geom_point(
        data = df,
        aes(UMAP_1, UMAP_2),
        color = BG_POINT, size = 0.09, alpha = 0.42, stroke = 0
      ) +
      geom_point(
        data = hi,
        aes(UMAP_1, UMAP_2, color = score_plot),
        size = 0.42, alpha = 1, stroke = 0
      ) +
      geom_text(
        data = lab_df,
        aes(UMAP_1, UMAP_2, label = state_label),
        inherit.aes = FALSE,
        size = 2.7, fontface = "bold", color = "grey10",
        check_overlap = TRUE
      ) +
      scale_color_gradientn(
        colors = pal,
        values = scales::rescale(c(0, 0.12, 0.35, 0.68, 1)),
        limits = c(0, score_cap),
        oob = scales::squish,
        name = "UCell",
        guide = guide_colorbar(
          barheight = unit(24, "mm"), barwidth = unit(3.4, "mm"),
          frame.colour = "grey35", ticks.colour = "white"
        )
      ) +
      coord_fixed(xlim = umap_xlim, ylim = umap_ylim, expand = FALSE, clip = "off") +
      labs(
        title = title,
        subtitle = paste0(
          "Grey = sampled primary cells; coloured = top ",
          scales::percent(1 - show_quantile), " of positive module scores"
        )
      ) +
      theme_void(base_size = 10) +
      theme(
        plot.background = element_rect(fill = PANEL_BG, color = NA),
        panel.border = element_rect(color = "grey78", fill = NA, linewidth = 0.35),
        plot.title = element_text(size = 11, face = "bold", hjust = 0,
                                  color = "grey10"),
        plot.subtitle = element_text(size = 7.8, color = "grey35", hjust = 0,
                                     margin = margin(b = 4)),
        legend.title = element_text(size = 8.5, color = "grey10"),
        legend.text = element_text(size = 7.5, color = "grey15"),
        plot.margin = margin(6, 18, 6, 7)
      )
  }
  
  p01a <- umap_focus(
    m1_ucell_cell,
    "M1 - cell-cycle module", ucell_m1_pal, m1_label_states, 0.50
  )
  p01b <- umap_focus(
    m2_ucell_cell,
    "M2 - neutrophil module", ucell_m2_pal, m2_label_states, 0.50
  )
  p01 <- combine_panels(list(p01a, p01b), ncol = 2, guides = "keep", tag = TRUE)
  save_any_pdf_png(p01, "fig01_module_ucell_umap_comparison", width = 11.8, height = 5.3)
  save_pdf_png(p01a, "fig01a_m1_ucell_umap", width = 6.2, height = 5.8)
  save_pdf_png(p01b, "fig01b_m2_ucell_umap", width = 6.2, height = 5.8)
  
  # ---- Fig 02: state summary dot plot, one panel per module (own scale) -----
  h3_all <- rbind(
    transform(matrix_to_long(m1_ucell_ds, "score"), module = "M1"),
    transform(matrix_to_long(m2_ucell_ds, "score"), module = "M2")
  )
  
  h3_state <- do.call(rbind, lapply(
    split(h3_all, list(h3_all$module, h3_all$fine_state), drop = TRUE),
    function(z) {
      data.frame(
        module = unique(z$module),
        fine_state = unique(z$fine_state),
        median_ucell = safe_median(z$score),
        eligible_donor_rate = mean(is.finite(z$score)),
        eligible_donor_n = sum(is.finite(z$score)),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  h3_state$module_f <- factor(h3_state$module, levels = c("M1", "M2"))
  h3_state$state_f <- factor(h3_state$fine_state, levels = rev(fine_states))
  h3_state$lineage_f <- factor(
    unname(state_lineage[h3_state$fine_state]),
    levels = lineage_order
  )
  
  make_state_dot <- function(mod, pal, title) {
    z <- h3_state[h3_state$module_f == mod, , drop = FALSE]
    cap <- as.numeric(stats::quantile(
      z$median_ucell[is.finite(z$median_ucell)], 0.99, na.rm = TRUE
    ))
    if (!is.finite(cap) || cap <= 0) cap <- 1
    
    z$top_rank <- ave(
      -z$median_ucell,
      z$module_f,
      FUN = function(x) rank(x, ties.method = "first", na.last = "keep")
    )
    z$value_label <- ifelse(
      is.finite(z$median_ucell) & z$top_rank <= 2,
      sprintf("%.2f", z$median_ucell),
      ""
    )
    
    ggplot(z, aes(factor(module), state_f)) +
      geom_point(
        aes(size = eligible_donor_rate, fill = median_ucell),
        shape = 21, color = "grey20", stroke = 0.30
      ) +
      geom_text(
        aes(label = value_label),
        nudge_x = 0.10, hjust = 0,
        size = 2.35, color = "grey15", na.rm = TRUE
      ) +
      facet_grid(lineage_f ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradientn(
        colors = pal,
        limits = c(0, cap), oob = scales::squish,
        na.value = NA_GREY, name = "Median\nUCell"
      ) +
      scale_size_area(
        max_size = 6.5, limits = c(0, 1),
        breaks = c(0.25, 0.50, 0.75, 1),
        labels = scales::percent_format(accuracy = 1),
        name = "Donors with\n>=20 cells"
      ) +
      scale_y_discrete(labels = short_state) +
      coord_cartesian(clip = "off") +
      labs(
        x = NULL, y = NULL,
        title = title,
        subtitle = "Fill = median UCell; dot size = donor coverage; labels = top two states"
      ) +
      theme_paper +
      theme(
        axis.text.x = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 7.6),
        panel.grid.major = element_blank(),
        panel.spacing.y = unit(1.4, "mm"),
        strip.text.y = element_text(angle = 0, size = 7.5),
        plot.margin = margin(7, 18, 7, 7)
      )
  }
  
  p02a <- make_state_dot("M1", ucell_m1_pal, "M1 - cell-cycle activity")
  p02b <- make_state_dot("M2", ucell_m2_pal, "M2 - neutrophil activity")
  
  if (TRUE) {
    save_fig8_source_panel(
      p02a,
      "Fig8A1_M1_UCell_state_summary",
      width = 5.8,
      height = 8.0
    )
    
    save_fig8_source_panel(
      p02b,
      "Fig8A2_M2_UCell_state_summary",
      width = 5.8,
      height = 8.0
    )
  }
  
  if (FALSE) {
    p02 <- combine_panels(list(p02a, p02b), ncol = 2, widths = c(1, 1), guides = "keep", tag = TRUE)
    save_any_pdf_png(p02, "fig02_ucell_state_summary", width = 11.6, height = 8.0)
  }
  
  # ---- Fig 03: detailed donor-state heatmaps -------------------------------
  # Use binned colours instead of a continuous near-white gradient. This makes
  # weak but evaluable donor-state scores visible while keeping grey exclusively
  # for not-evaluable donor-state combinations.
  ucell_bins <- c(-Inf, 0.02, 0.05, 0.10, 0.20, 0.40, Inf)
  ucell_bin_labels <- c("0-0.02", "0.02-0.05", "0.05-0.10",
                        "0.10-0.20", "0.20-0.40", ">0.40")
  ucell_m1_bin_pal <- c(
    "0-0.02" = "#E6F0F8", "0.02-0.05" = "#BDD7EA",
    "0.05-0.10" = "#6BAED6", "0.10-0.20" = "#3182BD",
    "0.20-0.40" = "#08519C", ">0.40" = "#08306B"
  )
  ucell_m2_bin_pal <- c(
    "0-0.02" = "#FEE8C8", "0.02-0.05" = "#FDBB84",
    "0.05-0.10" = "#FC8D59", "0.10-0.20" = "#EF6548",
    "0.20-0.40" = "#B30000", ">0.40" = "#7F0000"
  )
  
  make_module_heatmap <- function(mat, module_name, bin_pal, show_x = TRUE) {
    donor_score <- apply(mat, 1, safe_mean)
    donor_ord <- names(sort(donor_score, decreasing = TRUE, na.last = TRUE))
    
    d <- matrix_to_long(mat, "score")
    d$state_f <- factor(d$fine_state, levels = fine_states)
    d$donor_f <- factor(d$donor, levels = rev(donor_ord))
    d$score_bin <- cut(
      d$score,
      breaks = ucell_bins,
      labels = ucell_bin_labels,
      include.lowest = TRUE,
      right = TRUE
    )
    d$score_bin <- as.character(d$score_bin)
    d$score_bin[is.na(d$score_bin)] <- "not evaluable"
    d$score_bin <- factor(d$score_bin, levels = c(ucell_bin_labels, "not evaluable"))
    
    ggplot(d, aes(state_f, donor_f, fill = score_bin)) +
      geom_tile(width = 0.98, height = 0.98, color = "white", linewidth = 0.12) +
      scale_fill_manual(
        values = c(bin_pal, "not evaluable" = NA_GREY),
        breaks = c(ucell_bin_labels, "not evaluable"),
        drop = FALSE,
        name = paste0(module_name, "\nUCell")
      ) +
      scale_x_discrete(labels = short_state) +
      labs(
        x = if (show_x) "Fine cell state" else NULL,
        y = "Donor",
        title = paste0(module_name, " donor-state medians")
      ) +
      theme_paper +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 5.8),
        axis.text.x = if (show_x) {
          element_text(angle = 45, hjust = 1, size = 6.5)
        } else {
          element_blank()
        },
        axis.ticks.x = if (show_x) element_line() else element_blank(),
        plot.title = element_text(size = 9.5),
        legend.key.height = unit(3.4, "mm"),
        legend.key.width = unit(4.5, "mm")
      )
  }
  
  p03a <- make_module_heatmap(m1_ucell_ds, "M1", ucell_m1_bin_pal, show_x = FALSE)
  p03b <- make_module_heatmap(m2_ucell_ds, "M2", ucell_m2_bin_pal, show_x = TRUE) +
    labs(caption = paste0(
      "Grey = donor-state not evaluable (<20 cells or absent), not biological zero. ",
      "Blue/orange bins show evaluated UCell scores; the lowest nonzero bins are deliberately visible. ",
      "Donors are ordered separately within each module."
    ))
  p03 <- combine_panels(list(p03a, p03b), ncol = 1, heights = c(1, 1.14), guides = "keep", tag = TRUE)
  save_any_pdf_png(p03, "fig03_ucell_donor_state_heatmaps", width = 11.8, height = 8.7)
  
  # ---- Fig 04: M2 maturation axis ------------------------------------------
  m2_axis_long <- matrix_to_long(m2_ucell_ds[, axis_states, drop = FALSE], "ucell")
  m2_axis_long <- m2_axis_long[is.finite(m2_axis_long$ucell), , drop = FALSE]
  m2_axis_long$state_f <- factor(m2_axis_long$fine_state, levels = axis_states)
  
  p04a <- ggplot(m2_axis_long, aes(state_f, ucell)) +
    geom_boxplot(
      width = 0.55, outlier.shape = NA,
      fill = "white", color = "grey35", linewidth = 0.35
    ) +
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      size = 0.85, alpha = 0.35, color = "#D15F21"
    ) +
    stat_summary(
      fun = median, geom = "line", aes(group = 1),
      color = "#A63F0E", linewidth = 0.75
    ) +
    stat_summary(
      fun = median, geom = "point",
      shape = 21, size = 2.7, fill = "#A63F0E", color = "white", stroke = 0.35
    ) +
    scale_x_discrete(labels = short_state) +
    labs(
      x = NULL, y = "M2 UCell",
      title = "M2 activity across neutrophil maturation",
      subtitle = "Points = donors; boxes = IQR; line connects across-donor medians"
    ) +
    theme_paper +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(7, 8, 0, 8)
    )
  
  # Fixed M2 gene order: preserves the original M2 signature order and avoids
  # reordering rows in a way that makes the module hard to compare with Script 02.
  m2_gene_order <- c(
    "CEACAM6", "CEACAM8", "CTSG", "DEFA4", "ELANE",
    "MPO", "MS4A3", "OLFM4", "RNASE3", "TCN1"
  )
  m2_gene_order <- intersect(m2_gene_order, unique(g_logcpm_axis$signature_gene))
  
  g_logcpm_axis$state_f <- factor(g_logcpm_axis$fine_state, levels = axis_states)
  g_logcpm_axis$gene_f <- factor(g_logcpm_axis$signature_gene, levels = rev(m2_gene_order))
  logcpm_cap <- as.numeric(stats::quantile(
    g_logcpm_axis$median_logCPM[is.finite(g_logcpm_axis$median_logCPM)],
    0.99, na.rm = TRUE
  ))
  if (!is.finite(logcpm_cap) || logcpm_cap <= 0) logcpm_cap <- 1
  
  p04b <- ggplot(g_logcpm_axis, aes(state_f, gene_f, fill = median_logCPM)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradientn(
      colors = m2_logcpm_pal, limits = c(0, logcpm_cap),
      oob = scales::squish, na.value = "grey85", name = "Median\nlogCPM"
    ) +
    scale_x_discrete(labels = short_state) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
    labs(x = "Neutrophil maturation state", y = NULL) +
    theme_paper +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 7.2),
      axis.text.y = element_text(size = 7.8),
      plot.margin = margin(0, 8, 8, 8)
    )
  
  if (TRUE) {
    save_fig8_source_panel(
      p04a,
      "Fig8C1_M2_maturation_UCell_axis",
      width = 9.2,
      height = 2.35
    )
    
    save_fig8_source_panel(
      p04b,
      "Fig8C2_M2_maturation_gene_logCPM_heatmap",
      width = 9.2,
      height = 4.85
    )
  }
  
  if (FALSE) {
    p04 <- combine_panels(list(p04a, p04b), ncol = 1, heights = c(1.55, 3.2), guides = "keep", tag = TRUE)
    save_any_pdf_png(p04, "fig04_m2_maturation_axis", width = 9.2, height = 7.2)
  }
  
  # ---- Fig 05: M1 source classification with donor distributions -----------
  m1_cf_mat <- module_cf_list[["M1"]]
  m1_source_donor <- do.call(rbind, lapply(names(M1_SOURCE_CLASS), function(cl) {
    sts <- intersect(M1_SOURCE_CLASS[[cl]], colnames(m1_cf_mat))
    vals <- row_sum_preserve_all_na(m1_cf_mat[, sts, drop = FALSE])
    data.frame(
      donor = rownames(m1_cf_mat),
      source_class = cl,
      contribution_fraction = vals,
      stringsAsFactors = FALSE
    )
  }))
  
  # Fixed biological order. Keep the broad residual "other" class last even if
  # its median is high, because it is not one coherent cell population.
  source_levels <- c(
    "progenitor_granulopoiesis",
    "immature_neutrophil",
    "cycling_lymphoid",
    "other"
  )
  m1_source_donor$source_f <- factor(m1_source_donor$source_class, levels = source_levels)
  s05 <- m1_source
  s05$source_f <- factor(s05$source_class, levels = source_levels)
  
  source_axis_labels <- c(
    "progenitor_granulopoiesis" = "progenitor /\ngranulopoiesis",
    "immature_neutrophil" = "immature\nneutrophil",
    "cycling_lymphoid" = "cycling\nlymphoid",
    "other" = "other"
  )
  for (cl in source_levels) {
    source_axis_labels[cl] <- paste0(
      source_axis_labels[cl], "\n(", length(M1_SOURCE_CLASS[[cl]]), " states)"
    )
  }
  
  set.seed(20260803)
  p05 <- ggplot(m1_source_donor, aes(source_f, contribution_fraction)) +
    geom_boxplot(
      aes(fill = source_class), width = 0.52, outlier.shape = NA,
      alpha = 0.32, color = "black", linewidth = 0.45
    ) +
    geom_point(
      aes(color = source_class),
      position = position_jitter(width = 0.12, height = 0),
      size = 1.45, alpha = 0.78
    ) +
    stat_summary(
      fun = median, geom = "point",
      shape = 23, size = 3.6, fill = "white", color = "black", stroke = 0.55
    ) +
    geom_text(
      data = s05,
      aes(
        x = source_f,
        y = contribution_fraction_median,
        label = scales::percent(contribution_fraction_median, accuracy = 1)
      ),
      inherit.aes = FALSE, vjust = -1.15, size = 3.1, color = "black"
    ) +
    scale_fill_manual(values = m1_source_pal, guide = "none") +
    scale_color_manual(values = m1_source_pal, guide = "none") +
    scale_x_discrete(labels = source_axis_labels) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.02), expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "M1 contribution fraction per donor",
      title = "Competing cellular sources of the M1 cell-cycle module",
      subtitle = "Each point is one donor; diamond and label show the median",
      caption = wrap_caption(paste0(
        "The broad 'other' category contains all remaining fine states and should not be interpreted as one coherent cell population. ",
        "Medians are calculated across donors for each source class separately and therefore do not necessarily sum to 100%."
      ), width = 112)
    ) +
    theme_paper +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 7.7, lineheight = 0.95),
      axis.text.y = element_text(size = 8.5),
      plot.caption = element_text(size = 7.0, color = "grey40", hjust = 0, lineheight = 1.05),
      plot.margin = margin(8, 18, 14, 8)
    )
  
  save_pdf_png(p05, "fig05_m1_source_classification", width = 8.2, height = 5.9)
  
  # ---- Fig 06: within-module contributor-enrichment state map ---------------
  # Main Fig06 restores the per-gene heterogeneity map. The compact agreement
  # barplot is kept as Fig06b / supplement, not as a replacement.
  s06 <- het[het$module_of %in% c("Cell cycle", "Neutrophil"), , drop = FALSE]
  s06$contributor_equals_enriched[is.na(s06$contributor_equals_enriched)] <- FALSE
  
  agreement_detail <- data.frame(
    signature_gene = s06$signature_gene,
    module = s06$module_of,
    contributor_state = s06$main_contributor_state,
    enriched_state = s06$main_enriched_state,
    same_state = s06$contributor_equals_enriched,
    stringsAsFactors = FALSE
  )
  write_csv_atomic(
    agreement_detail,
    file.path(TABLE_DIR, "module_contributor_enriched_agreement_detail.csv")
  )
  
  module_levels <- c("Cell cycle", "Neutrophil")
  agreement_summary <- do.call(rbind, lapply(
    module_levels,
    function(mod) {
      z <- s06[s06$module_of == mod, , drop = FALSE]
      data.frame(
        module = mod,
        agree_n = sum(z$contributor_equals_enriched),
        total_n = nrow(z),
        stringsAsFactors = FALSE
      )
    }
  ))
  agreement_summary$agree_rate <- with(agreement_summary, agree_n / total_n)
  agreement_summary$module_strip <- paste0(
    agreement_summary$module, "  |  ",
    agreement_summary$agree_n, "/", agreement_summary$total_n,
    " concordant"
  )
  module_strip_map <- setNames(agreement_summary$module_strip, agreement_summary$module)
  
  fig06_gene_levels <- c(m1_genes, m2_genes)
  s06$gene_f <- factor(s06$signature_gene, levels = rev(fig06_gene_levels))
  s06$module_strip <- factor(
    unname(module_strip_map[s06$module_of]),
    levels = unname(module_strip_map[module_levels])
  )
  
  fig06_long <- rbind(
    data.frame(
      signature_gene = s06$signature_gene,
      module_strip = s06$module_strip,
      gene_f = s06$gene_f,
      role = "Bulk-like contributor",
      x = 1,
      state = s06$main_contributor_state,
      stringsAsFactors = FALSE
    ),
    data.frame(
      signature_gene = s06$signature_gene,
      module_strip = s06$module_strip,
      gene_f = s06$gene_f,
      role = "Per-cell enriched state",
      x = 2,
      state = s06$main_enriched_state,
      stringsAsFactors = FALSE
    )
  )
  fig06_long$role_f <- factor(
    fig06_long$role,
    levels = c("Bulk-like contributor", "Per-cell enriched state")
  )
  
  fig06_state_pal <- c(
    "MPO+_immature_neutrophils_or_progenitors" = "#0072B2",
    "Cycling_neutrophil_progenitors"           = "#6A3D9A",
    "Cycling_TNK"                              = "#D7191C",
    "Classical_monocytes"                      = "#1B9E77",
    "PADI4+_immature_neutrophils"              = "#56B4E9",
    "Plasmablasts"                             = "#CC79A7",
    "S100A8-9_hi_neutrophils"                  = "#E69F00"
  )
  fig06_states <- unique(fig06_long$state)
  missing_fig06_cols <- setdiff(fig06_states, names(fig06_state_pal))
  if (length(missing_fig06_cols)) {
    extra_cols <- grDevices::hcl.colors(length(missing_fig06_cols), palette = "Dark 3")
    names(extra_cols) <- missing_fig06_cols
    fig06_state_pal <- c(fig06_state_pal, extra_cols)
  }
  fig06_state_pal <- fig06_state_pal[names(fig06_state_pal) %in% fig06_states]
  fig06_state_labels <- setNames(short_state(names(fig06_state_pal)), names(fig06_state_pal))
  
  seg_same <- s06[s06$contributor_equals_enriched, , drop = FALSE]
  seg_diff <- s06[!s06$contributor_equals_enriched, , drop = FALSE]
  
  p06 <- ggplot() +
    geom_segment(
      data = seg_diff,
      aes(x = 1, xend = 2, y = gene_f, yend = gene_f),
      color = "grey78", linewidth = 0.35
    ) +
    geom_segment(
      data = seg_same,
      aes(x = 1, xend = 2, y = gene_f, yend = gene_f),
      color = "#66A61E", linewidth = 0.55
    ) +
    geom_point(
      data = fig06_long,
      aes(x = x, y = gene_f, fill = state, shape = role_f),
      size = 3.15, color = "grey15", stroke = 0.22, alpha = 0.98
    ) +
    facet_grid(module_strip ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Bulk-like\ncontributor", "Per-cell\nenriched state"),
      limits = c(0.72, 2.34), expand = expansion(mult = c(0.02, 0.08))
    ) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
    scale_fill_manual(values = fig06_state_pal, labels = fig06_state_labels, name = "Cell state") +
    scale_shape_manual(
      values = c("Bulk-like contributor" = 21, "Per-cell enriched state" = 24),
      name = "Role"
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(fill = "grey70", size = 3.2)),
      fill = guide_legend(order = 2, override.aes = list(shape = 21, size = 3.2))
    ) +
    labs(
      x = NULL, y = NULL,
      title = "Within-module heterogeneity: contributor versus enriched state",
      subtitle = "Discordance separates abundance-driven bulk-like contribution from per-cell enrichment",
      caption = wrap_caption(paste0(
        "Grey segments = discordant contributor/enriched states; green segments = concordant states. ",
        "M1 discordance reflects per-cell cycling enrichment but bulk-like contribution dominated by MPO+ immature neutrophil/progenitor states."
      ), width = 118)
    ) +
    theme_paper +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.y = unit(5.0, "mm"),
      axis.text.x = element_text(size = 8.7, face = "bold", color = "grey10", lineheight = 0.92),
      axis.text.y = element_text(size = 8.2, color = "grey10"),
      strip.background = element_blank(),
      strip.text.y = element_text(size = 9.0, face = "bold", angle = 0, color = "grey10",
                                  margin = margin(l = 4, r = 4)),
      legend.position = "right",
      legend.title = element_text(size = 8.5, face = "bold"),
      legend.text = element_text(size = 7.8),
      legend.key.height = unit(4.2, "mm"), legend.key.width = unit(4.2, "mm"),
      legend.spacing.y = unit(1.2, "mm"),
      plot.title = element_text(size = 11.2, face = "bold", color = "grey10"),
      plot.subtitle = element_text(size = 8.4, color = "grey35"),
      plot.caption = element_text(size = 7.1, color = "grey40", hjust = 0, lineheight = 1.05),
      plot.margin = margin(8, 18, 12, 8)
    )
  
  save_pdf_png(p06, "fig06_within_module_heterogeneity", width = 9.6, height = 6.5)
  
  # Fig06b: compact supplementary agreement summary.
  agreement_summary$module_f <- factor(agreement_summary$module, levels = module_levels)
  agreement_summary$bar_label <- paste0(
    agreement_summary$agree_n, "/", agreement_summary$total_n,
    " genes\n", scales::percent(agreement_summary$agree_rate, accuracy = 1)
  )
  agreement_pal <- c("Cell cycle" = "#3C5488", "Neutrophil" = "#E64B35")
  
  p06b <- ggplot(agreement_summary, aes(module_f, agree_rate, fill = module_f)) +
    geom_col(width = 0.56, color = "grey25", linewidth = 0.35, alpha = 0.92) +
    geom_text(
      aes(label = bar_label),
      vjust = -0.35, size = 3.4, fontface = "bold", color = "grey10", lineheight = 0.92
    ) +
    scale_fill_manual(values = agreement_pal, guide = "none") +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.25), limits = c(0, 1.08),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "Genes with same contributor and enriched state",
      title = "Contributor-enrichment concordance",
      subtitle = wrap_caption(
        "Supplementary summary; concordance is descriptive, not a quality score.",
        width = 70
      )
    ) +
    theme_paper +
    theme(
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 9.5, face = "bold"),
      axis.text.y = element_text(size = 8.5),
      plot.margin = margin(8, 16, 10, 8)
    )
  
  save_pdf_png(p06b, "fig06b_contributor_enriched_concordance_summary", width = 6.4, height = 4.4)
  
  # ---------------------------------------------------------------------------
  # 12. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  write_csv_atomic(data.frame(
    key = c("script","script_version","input_rds","script01_run_dir","script02_run_dir",
            "ucell_version","patchwork_available","normalization_input","n_primary_cells","n_donors",
            "m1_genes","m2_genes","m3_genes","m4_genes","R_version"),
    value = c("03_module_cell_distribution_readability_v16_final_clean.R","v16",INPUT_RDS,SCRIPT01_RUN_DIR,SCRIPT02_RUN_DIR,
              UCELL_VERSION, as.character(USE_PATCHWORK), norm_note, as.character(length(primary_idx)),
              as.character(n_donors),
              paste(m1_genes, collapse=";"), paste(m2_genes, collapse=";"),
              paste(m3_genes, collapse=";"), paste(m4_genes, collapse=";"),
              as.character(getRversion())),
    stringsAsFactors = FALSE), file.path(PROV_DIR, "run_provenance.csv"))
  writeLines(capture.output(sessionInfo()), file.path(PROV_DIR, "sessionInfo.txt"))
  
  cat("Script 03 v16 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("UCell:", UCELL_VERSION, "| patchwork:", USE_PATCHWORK, "| normalization:", norm_note, "\n")
  cat("M1 source classification (by contribution):\n")
  print(m1_source[, c("source_class","contribution_fraction_median","enrichment_ratio",
                      "ucell_donor_max_median","cell_fraction_median","eligible_donor_rate")])
  cat("M2 axis UCell medians:\n")
  print(data.frame(fine_state = axis_states, ucell = round(unname(m2_axis_ucell), 3)))
  cat("M2 axis UCell/logCPM audit:\n")
  print(m2_axis_audit[, c("fine_state", "m2_ucell_median", "m2_gene_logcpm_median",
                          "m2_gene_logcpm_max", "any_gene_logCPM_ge2", "interpretation_flag")])
})



# ============================================================================
# SCRIPT 03b v2.3 -- M2 cell-type program verification (Layer 1, must-do)
#   Question: why is the M2 neutrophil/granule module low in mature/downstream
#   neutrophil states in GSE216009?
#
#   Competing explanations:
#     (i)   developmental transcriptional decline along granulopoiesis,
#     (ii)  mature-granulocyte capture/QC limitation,
#     (iii) fine_annot boundary problem.
#
#   This script directly reads the Seurat RDS, reconstructs the same primary
#   acute-donor cohort as Scripts 02/03, performs donor-aware donor-state
#   aggregation, and outputs marker-set, per-gene, QC, and verdict tables.
#
#   RUN:
#     Rscript 03b_m2_celltype_program_verification_v2_3.R
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
  
  STAGE_NAME  <- "03b_m2_program_verification"
  RUN_PURPOSE <- "m2_celltype_program_audit"
  
  COUNT_ASSAY      <- "RNA"
  COUNT_LAYER      <- "counts"
  FINE_STATE_FIELD <- "fine_annot"
  SAMPLE_FIELD     <- "sample_id"
  CONDITION_FIELD  <- "diagnosis"
  
  ACUTE_DIAG <- c("Bacteraemia", "Bili", "CAP", "CNS", "IAS", "IE", "NF", "Uro")
  MIN_CELLS_PER_DONOR_STATE <- 20L
  
  EXPECTED_N_SAMPLES        <- 48L
  EXPECTED_N_DONORS         <- 39L
  EXPECTED_N_ACUTE_DONORS   <- 26L
  EXPECTED_N_HC_DONORS      <- 6L
  EXPECTED_N_SURGERY_DONORS <- 7L
  EXPECTED_N_CONV_SAMPLES   <- 9L
  EXPECTED_N_PAIRED_DONORS  <- 9L
  EXPECTED_N_PRIMARY_CELLS  <- 151837L
  
  # Marker-panel design.
  # MARS_M2_signature is the actual 10-gene M2 module from Scripts 01-03.
  # primary_azurophilic_core is separated because not all M2 genes are pure
  # canonical primary/azurophilic granule genes.
  marker_sets_raw <- list(
    MARS_M2_signature = c(
      "CEACAM6", "CEACAM8", "CTSG", "DEFA4", "ELANE",
      "MPO", "MS4A3", "OLFM4", "RNASE3", "TCN1"
    ),
    primary_azurophilic_core = c("MPO", "ELANE", "CTSG", "DEFA4"),
    secondary_specific_granule = c("LTF", "LCN2"),
    tertiary_gelatinase_granule = c("MMP8", "MMP9"),
    mature_surface_migration = c("FCGR3B", "CXCR2", "CXCR1", "FPR1", "MME"),
    inflammatory_S100 = c("S100A8", "S100A9", "S100A12"),
    capture_controls = c("B2M", "ACTB", "GAPDH", "RPLP0", "RPS18"),
    auxiliary_mature_markers_not_load_bearing = c("CD177", "CSF3R", "G0S2")
  )
  
  CRITICAL_MARKER_SETS <- c(
    "MARS_M2_signature",
    "primary_azurophilic_core",
    "mature_surface_migration",
    "capture_controls"
  )
  
  AXIS_STATES <- c(
    "Cycling_neutrophil_progenitors",
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils",
    "IL1R2+_immature_neutrophils",
    "S100A8-9_hi_neutrophils",
    "Mature_neutrophils",
    "Degranulating_neutrophils",
    "Apoptosing_neutrophils"
  )
  REF_STATE <- "MPO+_immature_neutrophils_or_progenitors"
  
  early_states <- c(
    "Cycling_neutrophil_progenitors",
    "MPO+_immature_neutrophils_or_progenitors",
    "PADI4+_immature_neutrophils"
  )
  transition_states <- c(
    "IL1R2+_immature_neutrophils",
    "S100A8-9_hi_neutrophils"
  )
  downstream_states <- c(
    "Mature_neutrophils",
    "Degranulating_neutrophils",
    "Apoptosing_neutrophils"
  )
  
  # Interpretive thresholds; these are flags, not hard biological truths.
  THRESH_PRESENT <- 0.20
  THRESH_REL_LOW <- 0.30
  THRESH_QC_LOW  <- 0.30
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS + PACKAGES
  # ---------------------------------------------------------------------------
  needed <- c("SeuratObject", "Matrix")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({
    library(SeuratObject)
    library(Matrix)
  })
  
  HAS_GGPLOT2 <- requireNamespace("ggplot2", quietly = TRUE)
  HAS_GRID    <- requireNamespace("grid", quietly = TRUE)
  
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  
  require_cols <- function(x, cols, label) {
    m <- setdiff(cols, names(x))
    if (length(m)) stop_msg("Missing required columns in ", label, ": ", paste(m, collapse = ", "))
    invisible(TRUE)
  }
  
  safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else median(x)
  }
  safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else max(x)
  }
  safe_iqr <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else unname(diff(quantile(x, c(0.25, 0.75))))
  }
  safe_div <- function(a, b) {
    if (!is.finite(a) || !is.finite(b) || b <= 0) return(NA_real_)
    a / b
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
  
  make_run_dir <- function(root, stage, purpose) {
    stage_dir <- file.path(root, stage)
    dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
    base <- paste0(format(Sys.Date(), "%Y%m%d"), "_", purpose)
    run_dir <- file.path(stage_dir, base)
    if (dir.exists(run_dir)) {
      i <- 2L
      repeat {
        cand <- file.path(stage_dir, paste0(base, "_", i))
        if (!dir.exists(cand)) {
          run_dir <- cand
          break
        }
        i <- i + 1L
      }
    }
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    run_dir
  }
  
  find_latest_file <- function(root, pattern) {
    if (!dir.exists(root)) stop_msg("Search root not found: ", root)
    fs <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (!length(fs)) stop_msg("No file found under ", root, " matching: ", pattern)
    fs[which.max(file.info(fs)$mtime)]
  }
  
  short_state <- function(x) {
    y <- x
    y <- gsub("MPO\\+_immature_neutrophils_or_progenitors", "MPO+ imm neut/prog", y)
    y <- gsub("Cycling_neutrophil_progenitors", "Cycling neut prog", y)
    y <- gsub("PADI4\\+_immature_neutrophils", "PADI4+ imm neut", y)
    y <- gsub("IL1R2\\+_immature_neutrophils", "IL1R2+ imm neut", y)
    y <- gsub("S100A8-9_hi_neutrophils", "S100A8/9 hi neut", y)
    y <- gsub("Mature_neutrophils", "Mature neut", y)
    y <- gsub("Degranulating_neutrophils", "Degran neut", y)
    y <- gsub("Apoptosing_neutrophils", "Apoptosing neut", y)
    y
  }
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  FIG_DIR   <- file.path(OUT_ROOT, "figures")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. READ INPUTS
  # ---------------------------------------------------------------------------
  STATE_ORDER_FILE <- find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  state01 <- read.csv(STATE_ORDER_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  require_cols(state01, c("display_order", "fine_state"), "fine_state_display_order_audit.csv")
  fine_states <- state01$fine_state[order(state01$display_order)]
  stopifnot(length(fine_states) == 24L, !anyDuplicated(fine_states))
  
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
    if ("layer" %in% names(formals(GetAssayData))) {
      cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, layer = "counts")
    } else {
      cnt <- GetAssayData(object = obj, assay = COUNT_ASSAY, slot = "counts")
    }
  }
  stopifnot(is(cnt, "dgCMatrix"), identical(colnames(cnt), colnames(obj)))
  
  # ---------------------------------------------------------------------------
  # 3. COHORT
  # ---------------------------------------------------------------------------
  require_cols(meta, c(SAMPLE_FIELD, CONDITION_FIELD, FINE_STATE_FIELD), "obj@meta.data")
  meta$sample_id <- as.character(meta[[SAMPLE_FIELD]])
  meta$diagnosis_for_verify <- as.character(meta[[CONDITION_FIELD]])
  
  sample_diag <- unique(meta[, c("sample_id", "diagnosis_for_verify"), drop = FALSE])
  n_diag_per_sample <- tapply(
    sample_diag$diagnosis_for_verify,
    sample_diag$sample_id,
    function(z) length(unique(z))
  )
  if (any(n_diag_per_sample != 1L)) stop_msg("sample_id with multiple diagnosis")
  
  sample_tab <- sample_diag
  sample_tab$donor_id <- sub("_CONV$", "", sample_tab$sample_id)
  sample_tab$cohort_group <- NA_character_
  sample_tab$cohort_group[sample_tab$diagnosis_for_verify %in% ACUTE_DIAG] <- "Acute_sepsis"
  sample_tab$cohort_group[sample_tab$diagnosis_for_verify == "Conv"] <- "Convalescent"
  sample_tab$cohort_group[sample_tab$diagnosis_for_verify == "HV"] <- "Healthy_control"
  sample_tab$cohort_group[sample_tab$diagnosis_for_verify == "CS"] <- "Surgery_control"
  if (anyNA(sample_tab$cohort_group)) {
    stop_msg("Unmapped diagnosis: ", paste(unique(sample_tab$diagnosis_for_verify[is.na(sample_tab$cohort_group)]), collapse = ", "))
  }
  
  is_conv_suffix <- grepl("_CONV$", sample_tab$sample_id)
  is_conv_diag   <- sample_tab$diagnosis_for_verify == "Conv"
  stopifnot(identical(is_conv_suffix, is_conv_diag))
  
  n_samples <- length(unique(sample_tab$sample_id))
  n_donors  <- length(unique(sample_tab$donor_id))
  acute_samples   <- sample_tab[sample_tab$cohort_group == "Acute_sepsis", , drop = FALSE]
  hc_samples      <- sample_tab[sample_tab$cohort_group == "Healthy_control", , drop = FALSE]
  surgery_samples <- sample_tab[sample_tab$cohort_group == "Surgery_control", , drop = FALSE]
  conv_samples    <- sample_tab[sample_tab$cohort_group == "Convalescent", , drop = FALSE]
  
  acute_donors   <- unique(acute_samples$donor_id)
  hc_donors      <- unique(hc_samples$donor_id)
  surgery_donors <- unique(surgery_samples$donor_id)
  conv_donors    <- unique(conv_samples$donor_id)
  paired_donors  <- intersect(acute_donors, conv_donors)
  
  stopifnot(n_samples == EXPECTED_N_SAMPLES, n_donors == EXPECTED_N_DONORS)
  stopifnot(
    length(acute_donors) == EXPECTED_N_ACUTE_DONORS,
    length(hc_donors) == EXPECTED_N_HC_DONORS,
    length(surgery_donors) == EXPECTED_N_SURGERY_DONORS,
    nrow(conv_samples) == EXPECTED_N_CONV_SAMPLES,
    length(paired_donors) == EXPECTED_N_PAIRED_DONORS,
    all(conv_donors %in% acute_donors)
  )
  
  m <- match(meta$sample_id, sample_tab$sample_id)
  if (anyNA(m)) stop_msg("meta$sample_id absent from sample_tab.")
  meta$donor_id <- sample_tab$donor_id[m]
  meta$cohort_group <- sample_tab$cohort_group[m]
  meta$is_primary_acute <- meta$cohort_group == "Acute_sepsis"
  stopifnot(all(!meta$is_primary_acute[grepl("_CONV$", meta$sample_id)]))
  
  primary_idx <- which(meta$is_primary_acute)
  stopifnot(length(primary_idx) == EXPECTED_N_PRIMARY_CELLS)
  
  # ---------------------------------------------------------------------------
  # 4. DONOR-STATE MATRICES
  # ---------------------------------------------------------------------------
  meta_p <- meta[primary_idx, , drop = FALSE]
  cnt_p  <- cnt[, primary_idx, drop = FALSE]
  
  donor_p <- factor(meta_p$donor_id, levels = acute_donors)
  state_p <- factor(meta_p[[FINE_STATE_FIELD]], levels = fine_states)
  stopifnot(!anyNA(donor_p), !anyNA(state_p))
  
  n_donor <- nlevels(donor_p)
  n_state <- length(fine_states)
  ds_idx <- (as.integer(donor_p) - 1L) * n_state + as.integer(state_p)
  
  Ds <- sparseMatrix(
    i = ds_idx,
    j = seq_along(ds_idx),
    x = 1,
    dims = c(n_donor * n_state, length(ds_idx))
  )
  
  n_cells_ds_flat <- as.integer(Matrix::rowSums(Ds))
  n_cells_ds <- matrix(n_cells_ds_flat, nrow = n_donor, ncol = n_state, byrow = TRUE)
  rownames(n_cells_ds) <- levels(donor_p)
  colnames(n_cells_ds) <- fine_states
  eligible_ds <- n_cells_ds >= MIN_CELLS_PER_DONOR_STATE
  
  lib_cell <- as.numeric(Matrix::colSums(cnt_p))
  lib_flat <- as.numeric(Ds %*% lib_cell)
  lib_ds <- matrix(lib_flat, nrow = n_donor, ncol = n_state, byrow = TRUE)
  rownames(lib_ds) <- levels(donor_p)
  colnames(lib_ds) <- fine_states
  
  nfeat_cell <- as.numeric(Matrix::colSums(cnt_p > 0))
  med_nfeat_ds <- rep(NA_real_, n_donor * n_state)
  med_numi_ds  <- rep(NA_real_, n_donor * n_state)
  tf1 <- tapply(nfeat_cell, ds_idx, safe_median)
  tf2 <- tapply(lib_cell,   ds_idx, safe_median)
  med_nfeat_ds[as.integer(names(tf1))] <- as.numeric(tf1)
  med_numi_ds[as.integer(names(tf2))]  <- as.numeric(tf2)
  
  states_report <- AXIS_STATES[AXIS_STATES %in% fine_states]
  ref_si <- match(REF_STATE, fine_states)
  if (is.na(ref_si)) stop_msg("REF_STATE not in fine_states: ", REF_STATE)
  ref_cols <- (seq_len(n_donor) - 1L) * n_state + ref_si
  ref_elig <- eligible_ds[, ref_si]
  ref_nUMI_med <- safe_median(med_numi_ds[ref_cols[ref_elig]])
  ref_nFeature_med <- safe_median(med_nfeat_ds[ref_cols[ref_elig]])
  
  qc_state_rows <- list()
  for (s in states_report) {
    si <- match(s, fine_states)
    cols <- (seq_len(n_donor) - 1L) * n_state + si
    elig <- eligible_ds[, si]
    qc_state_rows[[length(qc_state_rows) + 1L]] <- data.frame(
      fine_state = s,
      eligible_donor_n = sum(elig),
      n_cells_median = safe_median(n_cells_ds[elig, si]),
      library_size_median = safe_median(lib_ds[elig, si]),
      cell_nUMI_median = safe_median(med_numi_ds[cols[elig]]),
      cell_nFeature_median = safe_median(med_nfeat_ds[cols[elig]]),
      nUMI_ratio_vs_MPOimm = safe_div(safe_median(med_numi_ds[cols[elig]]), ref_nUMI_med),
      nFeature_ratio_vs_MPOimm = safe_div(safe_median(med_nfeat_ds[cols[elig]]), ref_nFeature_med),
      stringsAsFactors = FALSE
    )
  }
  qc_state_report <- do.call(rbind, qc_state_rows)
  write_csv_atomic(qc_state_report, file.path(TABLE_DIR, "axis_state_capture_qc_summary.csv"))
  
  # ---------------------------------------------------------------------------
  # 5. MARKER GENE MAPPING
  # ---------------------------------------------------------------------------
  feat <- rownames(cnt)
  all_requested_markers <- unique(unlist(marker_sets_raw, use.names = FALSE))
  
  marker_map <- data.frame(
    requested_gene = all_requested_markers,
    mapped_gene = ifelse(all_requested_markers %in% feat, all_requested_markers, NA_character_),
    mapped = all_requested_markers %in% feat,
    stringsAsFactors = FALSE
  )
  write_csv_atomic(marker_map, file.path(TABLE_DIR, "marker_gene_mapping_audit.csv"))
  
  missing_markers <- marker_map$requested_gene[!marker_map$mapped]
  if (length(missing_markers)) {
    warning("Some marker genes absent and will be skipped: ", paste(missing_markers, collapse = ", "))
  }
  
  marker_sets <- lapply(marker_sets_raw, function(gs) intersect(gs, feat))
  empty_sets <- names(marker_sets)[vapply(marker_sets, length, integer(1)) == 0L]
  critical_empty <- intersect(empty_sets, CRITICAL_MARKER_SETS)
  if (length(critical_empty)) {
    stop_msg("Critical marker sets with zero mapped genes: ", paste(critical_empty, collapse = ", "))
  }
  optional_empty <- setdiff(empty_sets, CRITICAL_MARKER_SETS)
  if (length(optional_empty)) {
    warning("Optional marker sets with zero mapped genes and will be reported as NA: ",
            paste(optional_empty, collapse = ", "))
  }
  
  all_markers <- unique(unlist(marker_sets[vapply(marker_sets, length, integer(1)) > 0L], use.names = FALSE))
  
  set_mapping_report <- do.call(rbind, lapply(names(marker_sets_raw), function(sn) {
    data.frame(
      marker_set = sn,
      requested_gene_n = length(marker_sets_raw[[sn]]),
      mapped_gene_n = length(marker_sets[[sn]]),
      requested_genes = paste(marker_sets_raw[[sn]], collapse = ";"),
      mapped_genes = paste(marker_sets[[sn]], collapse = ";"),
      missing_genes = paste(setdiff(marker_sets_raw[[sn]], marker_sets[[sn]]), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  write_csv_atomic(set_mapping_report, file.path(TABLE_DIR, "marker_set_mapping_audit.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. PER-GENE DONOR-STATE METRICS
  # ---------------------------------------------------------------------------
  gidx <- match(all_markers, feat)
  Xg <- cnt_p[gidx, , drop = FALSE]
  rownames(Xg) <- all_markers
  
  Xg_det <- Xg > 0
  det_cells_flat <- as.matrix(Xg_det %*% Matrix::t(Ds))
  det_rate_flat <- sweep(det_cells_flat, 2, n_cells_ds_flat, "/")
  det_rate_flat[, n_cells_ds_flat == 0] <- NA_real_
  
  C_gds <- as.matrix(Xg %*% Matrix::t(Ds))
  cpm_gds <- sweep(C_gds, 2, lib_flat, "/") * 1e6
  lc_gds <- log2(cpm_gds + 1)
  lc_gds[!is.finite(lc_gds)] <- NA_real_
  
  per_gene_rows <- list()
  for (g in all_markers) {
    gi <- match(g, all_markers)
    gene_sets <- names(marker_sets)[vapply(marker_sets, function(z) g %in% z, logical(1))]
    for (s in fine_states) {
      si <- match(s, fine_states)
      elig <- eligible_ds[, si]
      cols <- (seq_len(n_donor) - 1L) * n_state + si
      dr <- det_rate_flat[gi, cols]
      lc <- lc_gds[gi, cols]
      pooled_cells <- sum(n_cells_ds_flat[cols[elig]])
      pooled_det <- if (pooled_cells > 0) {
        sum(dr[elig] * n_cells_ds_flat[cols[elig]], na.rm = TRUE) / pooled_cells
      } else {
        NA_real_
      }
      per_gene_rows[[length(per_gene_rows) + 1L]] <- data.frame(
        gene = g,
        marker_sets = paste(gene_sets, collapse = ";"),
        fine_state = s,
        eligible_donor_n = sum(elig),
        median_detection_rate_per_donor = safe_median(dr[elig]),
        IQR_detection_rate_per_donor = safe_iqr(dr[elig]),
        pooled_detection_rate = pooled_det,
        median_logCPM = safe_median(lc[elig]),
        max_logCPM = safe_max(lc[elig]),
        stringsAsFactors = FALSE
      )
    }
  }
  per_gene_report <- do.call(rbind, per_gene_rows)
  write_csv_atomic(per_gene_report, file.path(TABLE_DIR, "per_gene_detection_logcpm_report.csv"))
  
  # ---------------------------------------------------------------------------
  # 7. MARKER-SET x STATE REPORT
  # ---------------------------------------------------------------------------
  set_names <- names(marker_sets)
  set_det <- matrix(NA_real_, nrow = length(set_names), ncol = n_donor * n_state,
                    dimnames = list(set_names, NULL))
  set_lc <- set_det
  
  for (sn in set_names) {
    gs <- intersect(marker_sets[[sn]], all_markers)
    if (length(gs)) {
      set_det[sn, ] <- colMeans(det_rate_flat[gs, , drop = FALSE], na.rm = FALSE)
      set_lc[sn, ]  <- colMeans(lc_gds[gs, , drop = FALSE], na.rm = FALSE)
    }
  }
  
  set_state_rows <- list()
  for (sn in set_names) {
    ref_det <- safe_median(set_det[sn, ref_cols[ref_elig]])
    ref_lc  <- safe_median(set_lc[sn,  ref_cols[ref_elig]])
    
    axis_det <- vapply(states_report, function(s2) {
      si2 <- match(s2, fine_states)
      e2 <- eligible_ds[, si2]
      cols2 <- (seq_len(n_donor) - 1L) * n_state + si2
      safe_median(set_det[sn, cols2[e2]])
    }, numeric(1))
    
    for (s in states_report) {
      si <- match(s, fine_states)
      elig <- eligible_ds[, si]
      cols <- (seq_len(n_donor) - 1L) * n_state + si
      dr <- set_det[sn, cols]
      lc <- set_lc[sn, cols]
      
      med_det <- safe_median(dr[elig])
      med_lc <- safe_median(lc[elig])
      pooled_cells <- sum(n_cells_ds_flat[cols[elig]])
      pooled_det <- if (pooled_cells > 0) {
        sum(dr[elig] * n_cells_ds_flat[cols[elig]], na.rm = TRUE) / pooled_cells
      } else {
        NA_real_
      }
      
      state_phase <- ifelse(
        s %in% early_states, "early",
        ifelse(s %in% transition_states, "transition",
               ifelse(s %in% downstream_states, "downstream", "other"))
      )
      
      rd <- safe_div(med_det, ref_det)
      rl <- safe_div(med_lc, ref_lc)
      
      set_state_rows[[length(set_state_rows) + 1L]] <- data.frame(
        marker_set = sn,
        mapped_gene_n = length(marker_sets[[sn]]),
        fine_state = s,
        state_phase = state_phase,
        eligible_donor_n = sum(elig),
        median_mean_gene_detection_rate = med_det,
        IQR_mean_gene_detection_rate = safe_iqr(dr[elig]),
        pooled_mean_gene_detection_rate = pooled_det,
        median_mean_gene_logCPM = med_lc,
        max_mean_gene_logCPM = safe_max(lc[elig]),
        relative_detection_vs_MPOimm = rd,
        relative_logCPM_vs_MPOimm = rl,
        relative_detection_vs_axis_max = safe_div(med_det, safe_max(axis_det)),
        median_library_size = safe_median(lib_ds[elig, si]),
        median_nUMI = safe_median(med_numi_ds[cols[elig]]),
        median_nFeature = safe_median(med_nfeat_ds[cols[elig]]),
        nUMI_ratio_vs_MPOimm = safe_div(safe_median(med_numi_ds[cols[elig]]), ref_nUMI_med),
        nFeature_ratio_vs_MPOimm = safe_div(safe_median(med_nfeat_ds[cols[elig]]), ref_nFeature_med),
        present_flag = is.finite(med_det) && med_det >= THRESH_PRESENT,
        relative_low_flag = is.finite(rd) && rd < THRESH_REL_LOW,
        stringsAsFactors = FALSE
      )
    }
  }
  set_state_report <- do.call(rbind, set_state_rows)
  write_csv_atomic(set_state_report, file.path(TABLE_DIR, "marker_set_state_report.csv"))
  
  # ---------------------------------------------------------------------------
  # 8. STATE-LEVEL VERDICT MATRIX
  # ---------------------------------------------------------------------------
  get_set_med_det <- function(sn, s) {
    if (!sn %in% rownames(set_det)) return(NA_real_)
    si <- match(s, fine_states)
    if (is.na(si)) return(NA_real_)
    e <- eligible_ds[, si]
    cols <- (seq_len(n_donor) - 1L) * n_state + si
    safe_median(set_det[sn, cols[e]])
  }
  get_set_med_lc <- function(sn, s) {
    if (!sn %in% rownames(set_lc)) return(NA_real_)
    si <- match(s, fine_states)
    if (is.na(si)) return(NA_real_)
    e <- eligible_ds[, si]
    cols <- (seq_len(n_donor) - 1L) * n_state + si
    safe_median(set_lc[sn, cols[e]])
  }
  
  imm_set  <- "MARS_M2_signature"
  pri_set  <- "primary_azurophilic_core"
  sec_set  <- "secondary_specific_granule"
  ter_set  <- "tertiary_gelatinase_granule"
  sur_set  <- "mature_surface_migration"
  s100_set <- "inflammatory_S100"
  cap_set  <- "capture_controls"
  
  imm_ref <- get_set_med_det(imm_set, REF_STATE)
  pri_ref <- get_set_med_det(pri_set, REF_STATE)
  
  verdict_rows <- list()
  for (s in states_report) {
    si <- match(s, fine_states)
    e <- eligible_ds[, si]
    cols <- (seq_len(n_donor) - 1L) * n_state + si
    
    qc_numi <- safe_div(safe_median(med_numi_ds[cols[e]]), ref_nUMI_med)
    qc_nfeat <- safe_div(safe_median(med_nfeat_ds[cols[e]]), ref_nFeature_med)
    qc_ok <- is.finite(qc_numi) && is.finite(qc_nfeat) &&
      qc_numi >= THRESH_QC_LOW && qc_nfeat >= THRESH_QC_LOW
    
    imm_det  <- get_set_med_det(imm_set, s)
    pri_det  <- get_set_med_det(pri_set, s)
    sec_det  <- get_set_med_det(sec_set, s)
    ter_det  <- get_set_med_det(ter_set, s)
    sur_det  <- get_set_med_det(sur_set, s)
    s100_det <- get_set_med_det(s100_set, s)
    cap_det  <- get_set_med_det(cap_set, s)
    
    imm_lc  <- get_set_med_lc(imm_set, s)
    pri_lc  <- get_set_med_lc(pri_set, s)
    sec_lc  <- get_set_med_lc(sec_set, s)
    ter_lc  <- get_set_med_lc(ter_set, s)
    sur_lc  <- get_set_med_lc(sur_set, s)
    s100_lc <- get_set_med_lc(s100_set, s)
    cap_lc  <- get_set_med_lc(cap_set, s)
    
    imm_present  <- is.finite(imm_det)  && imm_det  >= THRESH_PRESENT
    pri_present  <- is.finite(pri_det)  && pri_det  >= THRESH_PRESENT
    sec_present  <- is.finite(sec_det)  && sec_det  >= THRESH_PRESENT
    ter_present  <- is.finite(ter_det)  && ter_det  >= THRESH_PRESENT
    sur_present  <- is.finite(sur_det)  && sur_det  >= THRESH_PRESENT
    s100_present <- is.finite(s100_det) && s100_det >= THRESH_PRESENT
    cap_present  <- is.finite(cap_det)  && cap_det  >= THRESH_PRESENT
    
    imm_rel <- safe_div(imm_det, imm_ref)
    pri_rel <- safe_div(pri_det, pri_ref)
    imm_relative_low <- is.finite(imm_rel) && imm_rel < THRESH_REL_LOW
    pri_relative_low <- is.finite(pri_rel) && pri_rel < THRESH_REL_LOW
    
    state_phase <- ifelse(
      s %in% early_states, "early",
      ifelse(s %in% transition_states, "transition",
             ifelse(s %in% downstream_states, "downstream", "other"))
    )
    
    if (s %in% early_states && imm_present) {
      verdict <- "expected_early_granulopoiesis_program"
    } else if (s %in% early_states &&
               (sec_present || ter_present || sur_present || cap_present) &&
               qc_ok) {
      verdict <- "early_nonprimary_neutrophil_program_M2_low"
    } else if (s %in% early_states && !imm_present) {
      verdict <- "early_state_low_M2_unexpected"
    } else if (identical(s, "S100A8-9_hi_neutrophils") &&
               s100_present && sur_present && imm_relative_low && cap_present) {
      verdict <- "S100_inflammatory_mature_surface_program_primary_M2_low"
    } else if (s %in% transition_states &&
               imm_relative_low &&
               (sur_present || sec_present || ter_present || s100_present) &&
               qc_ok) {
      verdict <- "transition_or_inflammatory_neutrophil_program"
    } else if (s %in% downstream_states &&
               sur_present && (sec_present || ter_present) &&
               imm_relative_low && qc_ok) {
      verdict <- "developmental_decline"
    } else if (s %in% downstream_states &&
               !qc_ok && !sur_present && !sec_present && !ter_present && !cap_present) {
      verdict <- "capture_or_QC_limitation"
    } else if (s %in% downstream_states &&
               sur_present && imm_relative_low && !qc_ok && cap_present) {
      verdict <- "granule_low_surface_preserved_QC_reduced"
    } else if (s %in% downstream_states &&
               imm_present && !sur_present) {
      verdict <- "annotation_boundary_suspect"
    } else if (s %in% downstream_states &&
               sur_present && !sec_present && !ter_present && imm_relative_low) {
      verdict <- "granule_programs_low_surface_preserved"
    } else if (s %in% downstream_states &&
               imm_relative_low && qc_ok && cap_present) {
      verdict <- "low_M2_with_preserved_capture_controls"
    } else {
      verdict <- "needs_manual_review"
    }
    
    verdict_rows[[length(verdict_rows) + 1L]] <- data.frame(
      fine_state = s,
      short_state = short_state(s),
      state_phase = state_phase,
      eligible_donor_n = sum(e),
      MARS_M2_detection = imm_det,
      primary_core_detection = pri_det,
      secondary_granule_detection = sec_det,
      tertiary_granule_detection = ter_det,
      mature_surface_detection = sur_det,
      inflammatory_S100_detection = s100_det,
      capture_control_detection = cap_det,
      MARS_M2_logCPM = imm_lc,
      primary_core_logCPM = pri_lc,
      secondary_granule_logCPM = sec_lc,
      tertiary_granule_logCPM = ter_lc,
      mature_surface_logCPM = sur_lc,
      inflammatory_S100_logCPM = s100_lc,
      capture_control_logCPM = cap_lc,
      MARS_M2_relative_detection_vs_MPOimm = safe_div(imm_det, imm_ref),
      primary_core_relative_detection_vs_MPOimm = safe_div(pri_det, pri_ref),
      qc_nUMI_ratio_vs_MPOimm = qc_numi,
      qc_nFeature_ratio_vs_MPOimm = qc_nfeat,
      qc_ok = qc_ok,
      MARS_M2_present = imm_present,
      primary_core_present = pri_present,
      secondary_present = sec_present,
      tertiary_present = ter_present,
      surface_present = sur_present,
      S100_present = s100_present,
      capture_present = cap_present,
      MARS_M2_relative_low = imm_relative_low,
      primary_core_relative_low = pri_relative_low,
      verdict = verdict,
      stringsAsFactors = FALSE
    )
  }
  verdict_df <- do.call(rbind, verdict_rows)
  write_csv_atomic(verdict_df, file.path(TABLE_DIR, "state_verdict_matrix.csv"))
  
  # ---------------------------------------------------------------------------
  # 9. OPTIONAL FIGURES
  # ---------------------------------------------------------------------------
  if (HAS_GGPLOT2 && HAS_GRID) {
    suppressPackageStartupMessages({
      library(ggplot2)
      library(grid)
    })
    
    theme_audit <- function(base_size = 10) {
      theme_classic(base_size = base_size) %+replace%
        theme(
          plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0,
                                    margin = margin(b = 3)),
          plot.subtitle = element_text(size = base_size, colour = "grey30", hjust = 0,
                                       lineheight = 1.05, margin = margin(b = 8)),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(face = "plain"),
          legend.title = element_text(face = "plain"),
          panel.grid = element_blank(),
          plot.margin = margin(8, 12, 8, 8)
        )
    }
    
    set_state_axis <- set_state_report[set_state_report$fine_state %in% states_report, , drop = FALSE]
    set_state_axis$short_state <- factor(short_state(set_state_axis$fine_state), levels = short_state(states_report))
    set_order <- c(
      "MARS_M2_signature",
      "primary_azurophilic_core",
      "secondary_specific_granule",
      "tertiary_gelatinase_granule",
      "mature_surface_migration",
      "inflammatory_S100",
      "capture_controls"
    )
    set_state_axis <- set_state_axis[set_state_axis$marker_set %in% set_order, , drop = FALSE]
    set_state_axis$marker_set <- factor(set_state_axis$marker_set, levels = rev(set_order))
    
    max_fill <- safe_max(set_state_axis$median_mean_gene_detection_rate)
    if (!is.finite(max_fill)) max_fill <- 0.5
    
    p_marker <- ggplot(
      set_state_axis,
      aes(x = short_state, y = marker_set, fill = median_mean_gene_detection_rate)
    ) +
      geom_tile(colour = "white", linewidth = 0.25) +
      scale_fill_gradientn(
        colours = c("#F6F1E8", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
        limits = c(0, max(0.5, max_fill)),
        name = "Donor-aware\nmean detection"
      ) +
      labs(
        title = "M2 program verification across neutrophil states",
        subtitle = "Rows = marker programs; fill = median donor-aware mean-gene detection rate",
        x = "Neutrophil-axis state", y = NULL
      ) +
      theme_audit(10)
    
    ggsave(file.path(FIG_DIR, "fig_m2_marker_set_detection_heatmap.pdf"), p_marker,
           width = 11.5, height = 4.2, useDingbats = FALSE)
    ggsave(file.path(FIG_DIR, "fig_m2_marker_set_detection_heatmap.png"), p_marker,
           width = 11.5, height = 4.2, dpi = 320)
    
    verdict_plot <- verdict_df
    verdict_plot$short_state <- factor(verdict_plot$short_state, levels = short_state(states_report))
    verdict_long <- rbind(
      data.frame(fine_state = verdict_plot$short_state, metric = "MARS M2", value = verdict_plot$MARS_M2_detection),
      data.frame(fine_state = verdict_plot$short_state, metric = "Primary core", value = verdict_plot$primary_core_detection),
      data.frame(fine_state = verdict_plot$short_state, metric = "Secondary", value = verdict_plot$secondary_granule_detection),
      data.frame(fine_state = verdict_plot$short_state, metric = "Tertiary", value = verdict_plot$tertiary_granule_detection),
      data.frame(fine_state = verdict_plot$short_state, metric = "Mature surface", value = verdict_plot$mature_surface_detection),
      data.frame(fine_state = verdict_plot$short_state, metric = "Capture ctrl", value = verdict_plot$capture_control_detection)
    )
    verdict_long$metric <- factor(
      verdict_long$metric,
      levels = c("MARS M2", "Primary core", "Secondary", "Tertiary", "Mature surface", "Capture ctrl")
    )
    
    p_profile <- ggplot(verdict_long, aes(x = fine_state, y = value, group = metric, colour = metric)) +
      geom_line(linewidth = 0.7, alpha = 0.85) +
      geom_point(size = 2.2) +
      scale_colour_manual(values = c(
        "MARS M2" = "#B2182B",
        "Primary core" = "#D6604D",
        "Secondary" = "#2166AC",
        "Tertiary" = "#4393C3",
        "Mature surface" = "#1B7837",
        "Capture ctrl" = "#4D4D4D"
      )) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        title = "Marker-program detection profile",
        subtitle = "Developmental decline is supported when primary/M2 falls while surface or later-granule markers remain detectable",
        x = "Neutrophil-axis state",
        y = "Median donor-aware mean-gene detection rate",
        colour = "Program"
      ) +
      theme_audit(10)
    
    ggsave(file.path(FIG_DIR, "fig_m2_marker_program_detection_profile.pdf"), p_profile,
           width = 11.5, height = 4.5, useDingbats = FALSE)
    ggsave(file.path(FIG_DIR, "fig_m2_marker_program_detection_profile.png"), p_profile,
           width = 11.5, height = 4.5, dpi = 320)
    
    qc_long <- rbind(
      data.frame(short_state = verdict_df$short_state, metric = "nUMI",
                 ratio = verdict_df$qc_nUMI_ratio_vs_MPOimm, stringsAsFactors = FALSE),
      data.frame(short_state = verdict_df$short_state, metric = "nFeature",
                 ratio = verdict_df$qc_nFeature_ratio_vs_MPOimm, stringsAsFactors = FALSE)
    )
    qc_long$short_state <- factor(qc_long$short_state, levels = short_state(states_report))
    qc_long$metric <- factor(qc_long$metric, levels = c("nUMI", "nFeature"))
    
    p_qc <- ggplot(qc_long, aes(x = short_state, y = ratio, fill = metric)) +
      geom_col(position = position_dodge(width = 0.68), width = 0.30) +
      geom_hline(yintercept = THRESH_QC_LOW, linetype = "dashed", colour = "grey40") +
      scale_fill_manual(values = c("nUMI" = "#7570B3", "nFeature" = "#1B9E77")) +
      labs(
        title = "Capture/QC ratios relative to MPO+ immature neutrophils",
        subtitle = "Dashed line marks the low-QC interpretive threshold",
        x = "Neutrophil-axis state",
        y = "Ratio vs MPO+ immature state",
        fill = NULL
      ) +
      theme_audit(10)
    
    ggsave(file.path(FIG_DIR, "fig_m2_axis_capture_qc_ratio.pdf"), p_qc,
           width = 10.5, height = 4.0, useDingbats = FALSE)
    ggsave(file.path(FIG_DIR, "fig_m2_axis_capture_qc_ratio.png"), p_qc,
           width = 10.5, height = 4.0, dpi = 320)
  }
  
  # ---------------------------------------------------------------------------
  # 10. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  write_csv_atomic(data.frame(
    key = c(
      "script", "script_version", "input_rds", "script01_root", "state_order_file",
      "ref_state", "threshold_present", "threshold_rel_low", "threshold_qc_low",
      "n_acute_donors", "n_primary_cells", "min_cells_per_donor_state",
      "ggplot2_available", "R_version"
    ),
    value = c(
      "03b_m2_celltype_program_verification_v2_3.R", "v2.3", INPUT_RDS, SCRIPT01_ROOT, STATE_ORDER_FILE,
      REF_STATE, as.character(THRESH_PRESENT), as.character(THRESH_REL_LOW), as.character(THRESH_QC_LOW),
      as.character(length(acute_donors)), as.character(length(primary_idx)),
      as.character(MIN_CELLS_PER_DONOR_STATE), as.character(HAS_GGPLOT2), as.character(getRversion())
    ),
    stringsAsFactors = FALSE
  ), file.path(TABLE_DIR, "run_provenance.csv"))
  
  cat("Script 03b v2.3 complete.\n")
  cat("Run dir:", OUT_ROOT, "\n")
  cat("State order file:", STATE_ORDER_FILE, "\n")
  cat("Primary cells:", length(primary_idx), " | acute donors:", length(acute_donors), "\n")
  cat("Requested marker genes:", length(all_requested_markers),
      " | mapped:", length(all_markers),
      " | missing:", if (length(missing_markers)) paste(missing_markers, collapse = ", ") else "none", "\n")
  
  cat("\n--- Marker set mapping ---\n")
  print(set_mapping_report, row.names = FALSE)
  
  cat("\n--- State verdict matrix ---\n")
  print(verdict_df[, c(
    "fine_state", "state_phase", "MARS_M2_detection", "primary_core_detection",
    "secondary_granule_detection", "tertiary_granule_detection",
    "mature_surface_detection", "capture_control_detection",
    "qc_nUMI_ratio_vs_MPOimm", "qc_nFeature_ratio_vs_MPOimm", "verdict"
  )], row.names = FALSE)
  
  cat("\n--- MARS M2 signature: donor-aware detection by neutrophil-axis state ---\n")
  sub <- set_state_report[set_state_report$marker_set == imm_set,
                          c("fine_state", "eligible_donor_n",
                            "median_mean_gene_detection_rate",
                            "median_mean_gene_logCPM",
                            "relative_detection_vs_MPOimm",
                            "nUMI_ratio_vs_MPOimm",
                            "nFeature_ratio_vs_MPOimm")]
  print(sub, row.names = FALSE)
})


# ============================================================================
# SCRIPT 03c v1.7 -- M1 source-contribution donor bootstrap
#   Purpose: donor-level bootstrap support for whether M1 bulk-like contribution
#   is progenitor-dominated or shared with the broad residual compartment.
#   Reuses Script 03's donor-state module contribution fractions (no UCell/Seurat
#   rerun).
#
#   Current Script 03 contribution table is expected to be a complete finite
#   donor x fine-state grid. Under this NA-free input, fig05_matched and
#   all26_zero_absent are expected to be identical; the sensitivity branch is
#   retained only for backward compatibility and will emit a figure only if
#   future inputs differ.
#
#   Fail-closed (audit-grade): override path existence, state-order schema,
#   fine-state NA/empty/duplicate checks, fine-state completeness AND donor x
#   fine-state grid completeness, numeric/[0,1] validation, duplicate rows,
#   26-donor count, donor-internal contribution closure, source-class partition.
#   Tie-aware pairwise probability; main interpretation uses the fractional-tie
#   probability. A diagnostic records whether the two analyses are identical.
#   "other" is displayed as "broad residual".
# RUN: Rscript 03c_M1_source_bootstrap.R
# ============================================================================

options(stringsAsFactors = FALSE)
local({
  
  # ---------------------------------------------------------------------------
  # 0. CONFIG
  # ---------------------------------------------------------------------------
  PROJECT_ROOT  <- "/home/sunshine/predicate/singlecell"
  OUTPUTS_ROOT  <- file.path(PROJECT_ROOT, "outputs", "GSE216009_32gene_analysis")
  SCRIPT01_ROOT <- file.path(OUTPUTS_ROOT, "01_gene_mapping_detectability")
  SCRIPT03_ROOT <- file.path(OUTPUTS_ROOT, "03_module_cell_distribution")
  
  STAGE_NAME  <- "03c_m1_source_bootstrap"
  RUN_PURPOSE <- "m1_source_bootstrap"
  
  N_BOOT <- 10000L
  SEED   <- 20260804
  EXPECTED_N_ACUTE_DONORS <- 26L
  EXPECTED_N_FINE_STATES  <- 24L
  
  STATE_ORDER_FILE_OVERRIDE <- NA_character_
  CF_FILE_OVERRIDE          <- NA_character_
  
  M1_SOURCE_CLASS <- list(
    progenitor_granulopoiesis = c("HSPCs","Cycling_neutrophil_progenitors",
                                  "MPO+_immature_neutrophils_or_progenitors"),
    immature_neutrophil = c("PADI4+_immature_neutrophils","IL1R2+_immature_neutrophils",
                            "S100A8-9_hi_neutrophils"),
    cycling_lymphoid = c("Cycling_TNK")
  )
  NAMED_STATES <- unique(unlist(M1_SOURCE_CLASS, use.names = FALSE))
  
  # ---------------------------------------------------------------------------
  # 1. HELPERS
  # ---------------------------------------------------------------------------
  stop_msg <- function(...) stop(paste0(...), call. = FALSE)
  require_cols <- function(x, cols, table_name = "table") {
    miss <- setdiff(cols, names(x))
    if (length(miss)) stop_msg(table_name, " missing required columns: ", paste(miss, collapse = ", "))
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
  
  OUT_ROOT  <- make_run_dir(OUTPUTS_ROOT, STAGE_NAME, RUN_PURPOSE)
  TABLE_DIR <- file.path(OUT_ROOT, "tables")
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 2. INPUTS (fail-closed)
  # ---------------------------------------------------------------------------
  STATE_ORDER_FILE <- if (is.na(STATE_ORDER_FILE_OVERRIDE)) {
    find_latest_file(SCRIPT01_ROOT, "^fine_state_display_order_audit\\.csv$")
  } else STATE_ORDER_FILE_OVERRIDE
  
  CF_FILE <- if (is.na(CF_FILE_OVERRIDE)) {
    find_latest_file(SCRIPT03_ROOT, "^module_contribution_fraction_donor_state\\.csv$")
  } else CF_FILE_OVERRIDE
  
  if (!file.exists(STATE_ORDER_FILE)) stop_msg("STATE_ORDER_FILE not found: ", STATE_ORDER_FILE)
  if (!file.exists(CF_FILE))          stop_msg("CF_FILE not found: ", CF_FILE)
  
  state01 <- read.csv(STATE_ORDER_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  require_cols(state01, c("fine_state", "display_order"), "fine_state_display_order_audit.csv")
  
  state01$fine_state <- trimws(state01$fine_state)
  if (anyNA(state01$fine_state) || any(!nzchar(state01$fine_state))) {
    stop_msg("fine_state contains NA or empty values in state order file.")
  }
  state01$display_order <- suppressWarnings(as.integer(state01$display_order))
  if (anyNA(state01$display_order)) stop_msg("display_order contains NA/non-integer values in state order file.")
  if (anyDuplicated(state01$display_order)) stop_msg("Duplicated display_order values in state order file.")
  
  fine_states <- state01$fine_state[order(state01$display_order)]
  if (length(fine_states) != EXPECTED_N_FINE_STATES) {
    stop_msg("Expected ", EXPECTED_N_FINE_STATES, " fine states, found: ", length(fine_states))
  }
  if (anyDuplicated(fine_states)) {
    stop_msg("Duplicated fine_state values in state order file: ",
             paste(unique(fine_states[duplicated(fine_states)]), collapse = ", "))
  }
  
  M1_SOURCE_CLASS[["other"]] <- setdiff(fine_states, NAMED_STATES)
  stopifnot(length(M1_SOURCE_CLASS[["other"]]) == 17L)
  
  class_state_map <- do.call(rbind, lapply(names(M1_SOURCE_CLASS), function(cl) {
    data.frame(source_class = cl, fine_state = M1_SOURCE_CLASS[[cl]], stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(class_state_map$fine_state)) {
    dup_states <- unique(class_state_map$fine_state[duplicated(class_state_map$fine_state)])
    stop_msg("Fine states assigned to multiple source classes: ",
             paste(dup_states, collapse = ", "))
  }
  if (!setequal(class_state_map$fine_state, fine_states)) {
    stop_msg(
      "Source-class partition does not exactly cover all fine states.\n",
      "Missing from source classes: ",
      paste(setdiff(fine_states, class_state_map$fine_state), collapse = ", "),
      "\nExtra in source classes: ",
      paste(setdiff(class_state_map$fine_state, fine_states), collapse = ", ")
    )
  }
  write_csv_atomic(class_state_map, file.path(TABLE_DIR, "m1_source_class_state_map.csv"))
  
  dat <- read.csv(CF_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  require_cols(dat, c("module", "donor", "fine_state", "module_contribution_fraction"),
               "module_contribution_fraction_donor_state.csv")
  dat <- dat[dat$module == "M1", , drop = FALSE]
  if (!nrow(dat)) stop_msg("No M1 rows found in module_contribution_fraction_donor_state.csv: ", CF_FILE)
  
  dat$module_contribution_fraction <- suppressWarnings(as.numeric(dat$module_contribution_fraction))
  bad_value <- !is.finite(dat$module_contribution_fraction) |
    dat$module_contribution_fraction < -1e-8 |
    dat$module_contribution_fraction > 1 + 1e-8
  if (any(bad_value)) {
    stop_msg(
      "Invalid module_contribution_fraction values in M1 table: ",
      paste(
        paste0(dat$donor[bad_value], "||", dat$fine_state[bad_value], "=",
               dat$module_contribution_fraction[bad_value]),
        collapse = "; "
      )
    )
  }
  dat$module_contribution_fraction[dat$module_contribution_fraction < 0] <- 0
  dat$module_contribution_fraction[dat$module_contribution_fraction > 1] <- 1
  
  dup_key <- paste(dat$donor, dat$fine_state, sep = "||")
  if (anyDuplicated(dup_key)) {
    bad <- unique(dup_key[duplicated(dup_key)])
    stop_msg("Duplicate donor-fine_state rows in M1 contribution table: ",
             paste(head(bad, 10), collapse = ", "),
             if (length(bad) > 10) " ..." else "")
  }
  
  donors <- sort(unique(dat$donor))
  if (length(donors) != EXPECTED_N_ACUTE_DONORS)
    stop_msg("Expected ", EXPECTED_N_ACUTE_DONORS,
             " acute donors for Script 03c, found: ", length(donors))
  
  dat_states <- unique(dat$fine_state)
  missing_states <- setdiff(fine_states, dat_states)
  extra_states   <- setdiff(dat_states, fine_states)
  if (length(missing_states) || length(extra_states)) {
    stop_msg(
      "Fine-state mismatch between Script 01 state order and Script 03 contribution table.\n",
      "Missing from Script 03: ",
      if (length(missing_states)) paste(missing_states, collapse = ", ") else "none",
      "\nExtra in Script 03: ",
      if (length(extra_states)) paste(extra_states, collapse = ", ") else "none"
    )
  }
  
  expected_grid <- expand.grid(donor = donors, fine_state = fine_states, stringsAsFactors = FALSE)
  expected_key <- paste(expected_grid$donor, expected_grid$fine_state, sep = "||")
  observed_key <- paste(dat$donor, dat$fine_state, sep = "||")
  missing_pairs <- setdiff(expected_key, observed_key)
  extra_pairs   <- setdiff(observed_key, expected_key)
  if (length(missing_pairs) || length(extra_pairs)) {
    stop_msg(
      "Donor-fine_state grid is incomplete or contains unexpected pairs.\n",
      "Missing pairs: ",
      if (length(missing_pairs)) paste(head(missing_pairs, 20), collapse = ", ") else "none",
      if (length(missing_pairs) > 20) " ..." else "",
      "\nExtra pairs: ",
      if (length(extra_pairs)) paste(head(extra_pairs, 20), collapse = ", ") else "none",
      if (length(extra_pairs) > 20) " ..." else ""
    )
  }
  states <- fine_states
  
  # ---------------------------------------------------------------------------
  # 3. DONOR x STATE MATRIX + CLOSURE CHECK + CLASS SUMS
  # ---------------------------------------------------------------------------
  m <- matrix(NA_real_, nrow = length(donors), ncol = length(states),
              dimnames = list(donors, states))
  for (i in seq_len(nrow(dat))) m[dat$donor[i], dat$fine_state[i]] <- dat$module_contribution_fraction[i]
  
  donor_state_total <- rowSums(ifelse(is.finite(m), m, 0))
  bad_total <- abs(donor_state_total - 1) > 1e-4
  if (any(bad_total)) {
    stop_msg(
      "M1 module contribution fractions do not sum to 1 for some donors: ",
      paste(
        paste0(names(donor_state_total)[bad_total], "=", signif(donor_state_total[bad_total], 4)),
        collapse = "; "
      )
    )
  }
  
  class_sum <- function(sts, all_na_as = c("NA", "zero")) {
    all_na_as <- match.arg(all_na_as)
    sts_present <- intersect(sts, colnames(m))
    out <- rep(NA_real_, nrow(m))
    names(out) <- rownames(m)
    if (!length(sts_present)) {
      if (all_na_as == "zero") out[] <- 0
      return(out)
    }
    sub <- m[, sts_present, drop = FALSE]
    fin <- is.finite(sub)
    out <- rowSums(ifelse(fin, sub, 0))
    no_info <- rowSums(fin) == 0L
    if (all_na_as == "NA") out[no_info] <- NA_real_ else out[no_info] <- 0
    out
  }
  
  dcf_fig05 <- lapply(M1_SOURCE_CLASS, function(sts) class_sum(sts, all_na_as = "NA"))
  dcf_all26 <- lapply(M1_SOURCE_CLASS, function(sts) class_sum(sts, all_na_as = "zero"))
  
  dcf_identical <- isTRUE(all.equal(
    as.data.frame(dcf_fig05), as.data.frame(dcf_all26),
    tolerance = 1e-12, check.attributes = FALSE
  ))
  if (dcf_identical) {
    message("Note: fig05_matched and all26_zero_absent are identical because no NA donor-class values were present.")
  }
  
  write_donor_long <- function(dcf, path) {
    df <- do.call(rbind, lapply(names(dcf), function(cl)
      data.frame(donor = donors, source_class = cl,
                 contribution_fraction = unname(dcf[[cl]]),
                 stringsAsFactors = FALSE)))
    write_csv_atomic(df, path)
  }
  write_donor_long(dcf_fig05, file.path(TABLE_DIR, "m1_source_donor_class_contribution_fig05_matched.csv"))
  write_donor_long(dcf_all26, file.path(TABLE_DIR, "m1_source_donor_class_contribution_all26_zero_absent.csv"))
  
  # ---------------------------------------------------------------------------
  # 4. DONOR BOOTSTRAP (tie-aware; raw distributions incl. per-iteration medians)
  # ---------------------------------------------------------------------------
  run_bootstrap <- function(dcf, n_boot = N_BOOT, seed = SEED) {
    class_names <- names(dcf)
    n_don <- length(dcf[[1L]])
    set.seed(seed)
    diff_vec <- rep(NA_real_, n_boot)
    top_weight <- matrix(0, nrow = n_boot, ncol = length(class_names),
                         dimnames = list(NULL, class_names))
    boot_medians <- matrix(NA_real_, nrow = n_boot, ncol = length(class_names),
                           dimnames = list(NULL, paste0("median_", class_names)))
    diff_valid <- rep(FALSE, n_boot)
    top_valid  <- rep(FALSE, n_boot)
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n_don, n_don, replace = TRUE)
      meds <- vapply(dcf, function(v) safe_median(v[idx]), numeric(1))
      boot_medians[b, ] <- meds
      finite_meds <- is.finite(meds)
      if (sum(finite_meds) > 0L) {
        mx <- max(meds[finite_meds])
        is_top <- finite_meds & abs(meds - mx) <= sqrt(.Machine$double.eps)
        top_weight[b, is_top] <- 1 / sum(is_top)
        top_valid[b] <- TRUE
      }
      if (is.finite(meds[["progenitor_granulopoiesis"]]) &&
          is.finite(meds[["other"]])) {
        diff_vec[b] <- meds[["progenitor_granulopoiesis"]] - meds[["other"]]
        diff_valid[b] <- TRUE
      }
    }
    n_diff <- sum(diff_valid); n_top <- sum(top_valid)
    diff_use <- diff_vec[diff_valid]
    eps <- sqrt(.Machine$double.eps)
    P_prog_gt_other <- if (n_diff) mean(diff_use > eps) else NA_real_
    P_prog_eq_other <- if (n_diff) mean(abs(diff_use) <= eps) else NA_real_
    P_prog_gt_other_fractional_tie <- if (n_diff) {
      mean(diff_use > eps) + 0.5 * mean(abs(diff_use) <= eps)
    } else NA_real_
    ci <- if (n_diff) unname(quantile(diff_use, c(0.025, 0.975))) else c(NA_real_, NA_real_)
    P_class_top <- if (n_top) colSums(top_weight[top_valid, , drop = FALSE]) / n_top
    else rep(NA_real_, length(class_names))
    names(P_class_top) <- class_names
    list(n_diff = n_diff, n_top = n_top,
         P_prog_gt_other = P_prog_gt_other,
         P_prog_eq_other = P_prog_eq_other,
         P_prog_gt_other_fractional_tie = P_prog_gt_other_fractional_tie,
         ci = ci, P_class_top = P_class_top,
         diff_vec = diff_vec, diff_valid = diff_valid,
         top_weight = top_weight, top_valid = top_valid,
         boot_medians = boot_medians)
  }
  
  write_boot_distribution <- function(analysis, boot, path) {
    top_weight_df <- as.data.frame(boot$top_weight, check.names = FALSE)
    names(top_weight_df) <- paste0("top_weight_", names(top_weight_df))
    boot_median_df <- as.data.frame(boot$boot_medians, check.names = FALSE)
    df <- data.frame(
      analysis = analysis,
      bootstrap_id = seq_along(boot$diff_vec),
      diff_progenitor_minus_residual = boot$diff_vec,
      diff_valid = boot$diff_valid,
      top_valid = boot$top_valid,
      boot_median_df,
      top_weight_df,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    write_csv_atomic(df, path)
  }
  
  boot_fig05 <- run_bootstrap(dcf_fig05)
  boot_all26 <- run_bootstrap(dcf_all26)
  
  write_boot_distribution("fig05_matched", boot_fig05,
                          file.path(TABLE_DIR, "m1_source_bootstrap_distribution_fig05_matched.csv"))
  write_boot_distribution("all26_zero_absent", boot_all26,
                          file.path(TABLE_DIR, "m1_source_bootstrap_distribution_all26_zero_absent.csv"))
  
  # ---------------------------------------------------------------------------
  # 5. INTERPRETATION + SUMMARY (human long + machine wide)
  #   Main interpretation uses the fractional-tie probability (ties share mass).
  # ---------------------------------------------------------------------------
  interpret <- function(P, ci) {
    prob <- if (is.na(P)) "NA" else if (P > 0.8) "progenitor tends to exceed broad residual"
    else if (P >= 0.5) "progenitor and residual contributions comparable"
    else "bulk-like contribution not dominated by progenitor"
    ci_txt <- if (anyNA(ci)) "NA"
    else if (ci[1] <= 0 && ci[2] >= 0) "CI includes 0: report shared/comparable, not dominant"
    else if (ci[1] > 0) "CI excludes 0 above: progenitor exceeds residual"
    else "CI excludes 0 below: residual exceeds progenitor"
    c(prob, ci_txt)
  }
  
  build_summary <- function(analysis, dcf, boot) {
    obs_med <- vapply(dcf, safe_median, numeric(1))
    obs_diff <- obs_med[["progenitor_granulopoiesis"]] - obs_med[["other"]]
    n_fin <- vapply(dcf, function(v) sum(is.finite(v)), integer(1))
    n_nz  <- vapply(dcf, function(v) sum(is.finite(v) & v > 0), integer(1))
    interp <- interpret(boot$P_prog_gt_other_fractional_tie, boot$ci)
    data.frame(
      analysis = analysis,
      metric = c(
        "n_donors","n_bootstrap","n_bootstrap_diff_usable","n_bootstrap_top_usable",
        paste0("obs_median_", names(obs_med)),
        paste0("obs_n_finite_", names(n_fin)),
        paste0("obs_n_nonzero_", names(n_nz)),
        "obs_median_diff_progenitor_minus_other",
        "P_progenitor_gt_other",
        "P_progenitor_eq_other",
        "P_progenitor_gt_other_fractional_tie",
        paste0("P_", names(boot$P_class_top), "_top"),
        "diff_CI_2.5","diff_CI_97.5",
        "interpretation_probability","interpretation_CI"),
      value = c(
        length(dcf[[1L]]), N_BOOT, boot$n_diff, boot$n_top,
        format(obs_med, digits = 4),
        as.character(n_fin), as.character(n_nz),
        format(obs_diff, digits = 4),
        format(boot$P_prog_gt_other, digits = 3),
        format(boot$P_prog_eq_other, digits = 3),
        format(boot$P_prog_gt_other_fractional_tie, digits = 3),
        format(boot$P_class_top, digits = 3),
        format(boot$ci[1], digits = 4), format(boot$ci[2], digits = 4),
        interp[1], interp[2]),
      stringsAsFactors = FALSE)
  }
  
  build_summary_wide <- function(analysis, dcf, boot) {
    obs_med <- vapply(dcf, safe_median, numeric(1))
    obs_diff <- obs_med[["progenitor_granulopoiesis"]] - obs_med[["other"]]
    n_fin <- vapply(dcf, function(v) sum(is.finite(v)), integer(1))
    n_nz  <- vapply(dcf, function(v) sum(is.finite(v) & v > 0), integer(1))
    data.frame(
      analysis = analysis,
      n_donors = length(dcf[[1L]]),
      n_bootstrap = N_BOOT,
      obs_median_progenitor_granulopoiesis = obs_med[["progenitor_granulopoiesis"]],
      obs_median_other = obs_med[["other"]],
      obs_median_immature_neutrophil = obs_med[["immature_neutrophil"]],
      obs_median_cycling_lymphoid = obs_med[["cycling_lymphoid"]],
      obs_median_diff_progenitor_minus_other = obs_diff,
      P_progenitor_gt_other = boot$P_prog_gt_other,
      P_progenitor_eq_other = boot$P_prog_eq_other,
      P_progenitor_gt_other_fractional_tie = boot$P_prog_gt_other_fractional_tie,
      P_progenitor_top = boot$P_class_top[["progenitor_granulopoiesis"]],
      P_other_top = boot$P_class_top[["other"]],
      P_immature_neutrophil_top = boot$P_class_top[["immature_neutrophil"]],
      P_cycling_lymphoid_top = boot$P_class_top[["cycling_lymphoid"]],
      diff_CI_2.5 = boot$ci[1],
      diff_CI_97.5 = boot$ci[2],
      n_finite_progenitor_granulopoiesis = n_fin[["progenitor_granulopoiesis"]],
      n_finite_other = n_fin[["other"]],
      n_finite_immature_neutrophil = n_fin[["immature_neutrophil"]],
      n_finite_cycling_lymphoid = n_fin[["cycling_lymphoid"]],
      n_nonzero_progenitor_granulopoiesis = n_nz[["progenitor_granulopoiesis"]],
      n_nonzero_other = n_nz[["other"]],
      n_nonzero_immature_neutrophil = n_nz[["immature_neutrophil"]],
      n_nonzero_cycling_lymphoid = n_nz[["cycling_lymphoid"]],
      residual_definition = "other = broad residual compartment of 17 fine states",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  summary_df <- rbind(
    build_summary("fig05_matched",    dcf_fig05, boot_fig05),
    build_summary("all26_zero_absent", dcf_all26, boot_all26))
  write_csv_atomic(summary_df, file.path(TABLE_DIR, "m1_source_bootstrap_summary.csv"))
  
  summary_wide <- rbind(
    build_summary_wide("fig05_matched",    dcf_fig05, boot_fig05),
    build_summary_wide("all26_zero_absent", dcf_all26, boot_all26))
  write_csv_atomic(summary_wide, file.path(TABLE_DIR, "m1_source_bootstrap_summary_wide.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. NATURE-STYLE VISUALIZATION
  # ---------------------------------------------------------------------------
  FIG_DIR <- file.path(OUT_ROOT, "figures")
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  
  PLOT_PKGS <- c("ggplot2", "patchwork")
  HAS_PLOT <- all(vapply(PLOT_PKGS, requireNamespace, logical(1), quietly = TRUE))
  
  if (!HAS_PLOT) {
    warning("Skipping figures because required plotting packages are missing: ",
            paste(PLOT_PKGS[!vapply(PLOT_PKGS, requireNamespace, logical(1), quietly = TRUE)],
                  collapse = ", "))
  } else {
    
    suppressPackageStartupMessages({
      library(ggplot2)
      library(patchwork)
    })
    
    source_order <- c("progenitor_granulopoiesis", "other",
                      "immature_neutrophil", "cycling_lymphoid")
    
    source_labels <- c(
      progenitor_granulopoiesis = "progenitor /\ngranulopoiesis",
      other                     = "broad residual\n(17 states)",
      immature_neutrophil       = "immature\nneutrophil",
      cycling_lymphoid          = "cycling\nlymphoid"
    )
    
    source_cols <- c(
      progenitor_granulopoiesis = "#5B8CC0",
      other                     = "#9A938A",
      immature_neutrophil       = "#58A87D",
      cycling_lymphoid          = "#B47AA5"
    )
    
    source_fill <- c(
      progenitor_granulopoiesis = "#C9D6E8",
      other                     = "#D9D4CE",
      immature_neutrophil       = "#CFE6D8",
      cycling_lymphoid          = "#E3CDDC"
    )
    
    theme_nature_03c <- function(base_size = 7, base_family = PLOT_FONT_FAMILY) {
      theme_classic(base_size = base_size, base_family = base_family) +
        theme(
          axis.line = element_line(linewidth = 0.35, colour = "#202020"),
          axis.ticks = element_line(linewidth = 0.35, colour = "#202020"),
          axis.title = element_text(size = base_size + 0.5, colour = "#202020"),
          axis.text = element_text(size = base_size, colour = "#202020"),
          plot.title = element_text(size = base_size + 1.2, face = "bold", colour = "#202020",
                                    margin = margin(b = 2)),
          plot.subtitle = element_text(size = base_size, colour = "#5F5F5F",
                                       margin = margin(b = 5)),
          plot.caption = element_text(size = base_size - 0.5, colour = "#666666", hjust = 0),
          legend.position = "none",
          panel.grid.major.y = element_line(linewidth = 0.20, colour = "#E8E8E8"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(5, 6, 5, 6)
        )
    }
    
    pct_lab <- function(x) paste0(round(100 * x), "%")
    
    save_pub_plot <- function(plot, filename, width_mm = 183, height_mm = 105, dpi = 600) {
      w <- width_mm / 25.4
      h <- height_mm / 25.4
      
      plot <- fontify_plot(plot)
      
      svglite::svglite(
        paste0(filename, ".svg"),
        width = w,
        height = h,
        bg = "white",
        system_fonts = list(
          sans = PLOT_FONT_FAMILY,
          serif = PLOT_FONT_FAMILY,
          mono = PLOT_FONT_FAMILY
        )
      )
      print(plot)
      grDevices::dev.off()
      
      grDevices::cairo_pdf(
        filename = paste0(filename, ".pdf"),
        width = w,
        height = h,
        family = PLOT_FONT_FAMILY,
        onefile = FALSE,
        bg = "white"
      )
      print(plot)
      grDevices::dev.off()
      
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_tiff(
          paste0(filename, ".tiff"),
          width = w,
          height = h,
          units = "in",
          res = dpi,
          compression = "lzw",
          background = "white"
        )
        print(plot)
        grDevices::dev.off()
        
        ragg::agg_png(
          paste0(filename, ".png"),
          width = w,
          height = h,
          units = "in",
          res = 450,
          background = "white"
        )
        print(plot)
        grDevices::dev.off()
      } else {
        grDevices::tiff(
          paste0(filename, ".tiff"),
          width = w,
          height = h,
          units = "in",
          res = dpi,
          compression = "lzw",
          bg = "white"
        )
        print(plot)
        grDevices::dev.off()
        
        grDevices::png(
          paste0(filename, ".png"),
          width = w,
          height = h,
          units = "in",
          res = 450,
          bg = "white"
        )
        print(plot)
        grDevices::dev.off()
      }
      
      invisible(filename)
    }
    
    donor_long_df <- function(dcf, analysis_label) {
      df <- do.call(rbind, lapply(names(dcf), function(cl) {
        data.frame(
          donor = donors,
          source_class = cl,
          contribution_fraction = unname(dcf[[cl]]),
          analysis = analysis_label,
          stringsAsFactors = FALSE
        )
      }))
      df$source_class <- factor(df$source_class, levels = source_order)
      df$source_label <- source_labels[as.character(df$source_class)]
      df
    }
    
    make_donor_panel <- function(dcf, title, subtitle) {
      df <- donor_long_df(dcf, title)
      df_plot <- df[is.finite(df$contribution_fraction), , drop = FALSE]
      
      med_df <- aggregate(contribution_fraction ~ source_class, df_plot, median)
      med_df$label <- pct_lab(med_df$contribution_fraction)
      
      y_upper <- max(df_plot$contribution_fraction, med_df$contribution_fraction, na.rm = TRUE)
      y_upper <- min(1.05, max(0.20, y_upper + 0.12))
      med_df$y_lab <- pmin(med_df$contribution_fraction + 0.055, y_upper * 0.96)
      
      ggplot(df_plot, aes(x = source_class, y = contribution_fraction)) +
        geom_boxplot(aes(fill = source_class),
                     width = 0.55, outlier.shape = NA,
                     linewidth = 0.40, colour = "#202020", alpha = 0.85) +
        geom_jitter(aes(colour = source_class),
                    width = 0.12, height = 0,
                    size = 1.35, alpha = 0.82) +
        stat_summary(fun = median, geom = "point",
                     shape = 23, size = 2.3,
                     fill = "white", colour = "#202020",
                     linewidth = 0.35) +
        geom_text(data = med_df,
                  aes(x = source_class, y = y_lab, label = label),
                  inherit.aes = FALSE,
                  size = 2.35, colour = "#202020") +
        scale_fill_manual(values = source_fill, drop = FALSE) +
        scale_colour_manual(values = source_cols, drop = FALSE) +
        scale_x_discrete(labels = source_labels, drop = FALSE) +
        scale_y_continuous(labels = pct_lab,
                           limits = c(0, y_upper),
                           expand = expansion(mult = c(0, 0.03))) +
        labs(title = title, subtitle = subtitle,
             x = NULL, y = "M1 contribution fraction per donor") +
        theme_nature_03c() +
        theme(
          axis.text.x = element_text(size = 6.6, lineheight = 0.90),
          plot.title = element_text(size = 8.5, face = "bold")
        )
    }
    
    make_top_prob_panel <- function(boot, title = "Top-contributor probability") {
      df <- data.frame(
        source_class = names(boot$P_class_top),
        probability = unname(boot$P_class_top),
        stringsAsFactors = FALSE
      )
      
      # Horizontal display: first source class appears at the top.
      df$source_class <- factor(df$source_class, levels = rev(source_order))
      df <- df[order(df$source_class), , drop = FALSE]
      
      df$label <- ifelse(is.finite(df$probability), pct_lab(df$probability), "NA")
      
      x_upper <- min(
        1,
        max(0.60, max(df$probability, na.rm = TRUE) + 0.12)
      )
      df$x_lab <- pmin(df$probability + 0.025, x_upper * 0.96)
      
      ggplot(df, aes(y = source_class, x = probability)) +
        geom_col(
          aes(fill = source_class),
          width = 0.62,
          colour = "#202020",
          linewidth = 0.30,
          alpha = 0.95
        ) +
        geom_text(
          aes(x = x_lab, label = label),
          hjust = 0,
          size = 2.8,
          colour = "#202020"
        ) +
        scale_fill_manual(values = source_cols, guide = "none", drop = FALSE) +
        scale_y_discrete(labels = source_labels, drop = FALSE) +
        scale_x_continuous(
          labels = pct_lab,
          limits = c(0, x_upper),
          breaks = c(0, 0.25, 0.50, 0.75, 1.00),
          expand = expansion(mult = c(0, 0.03))
        ) +
        labs(
          title = title,
          x = "Probability",
          y = NULL
        ) +
        theme_nature_03c() +
        theme(
          axis.text.y = element_text(size = 7.2, lineheight = 0.92),
          axis.text.x = element_text(size = 7.0),
          axis.title.x = element_text(size = 7.4, face = "bold"),
          plot.title = element_text(size = 8.2, face = "bold"),
          plot.margin = margin(4, 10, 4, 4)
        )
    }
    
    make_diff_ci_panel <- function(dcf, boot, title = "Progenitor minus residual contribution") {
      obs_med <- vapply(dcf, safe_median, numeric(1))
      est <- obs_med[["progenitor_granulopoiesis"]] - obs_med[["other"]]
      lo <- boot$ci[1]
      hi <- boot$ci[2]
      
      lim <- max(abs(c(est, lo, hi)), na.rm = TRUE)
      if (!is.finite(lim) || lim <= 0) lim <- 0.05
      lim <- lim * 1.35
      
      df <- data.frame(
        y = 1,
        estimate = est,
        low = lo,
        high = hi,
        label = paste0("Δ = ", sprintf("%.2f", est),
                       "\n95% CI ", sprintf("%.2f", lo), " to ", sprintf("%.2f", hi)),
        stringsAsFactors = FALSE
      )
      
      ggplot(df, aes(y = y)) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.35, colour = "#777777") +
        geom_segment(aes(x = low, xend = high, yend = y),
                     linewidth = 0.65, colour = "#202020") +
        geom_point(aes(x = estimate),
                   shape = 23, size = 2.9,
                   fill = "white", colour = "#202020",
                   linewidth = 0.45) +
        geom_text(aes(x = estimate, y = 1.18, label = label),
                  size = 2.25, lineheight = 0.88, colour = "#202020") +
        scale_x_continuous(limits = c(-lim, lim),
                           labels = function(x) sprintf("%.2f", x)) +
        scale_y_continuous(limits = c(0.82, 1.28), breaks = NULL) +
        labs(title = title,
             x = "Median difference in contribution fraction",
             y = NULL) +
        theme_nature_03c() +
        theme(
          panel.grid = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 8.0, face = "bold")
        )
    }
    
    # Figure 8 source-panel export switches for this 03c local block.
    export_fig8_source_panels <- TRUE
    export_fig8_composites    <- FALSE
    # Main manuscript-facing bootstrap figure (all26_zero_absent).
    p_main_a <- make_donor_panel(
      dcf_all26,
      title = "M1 bulk-like contribution is shared across source classes",
      subtitle = "Each point is one acute donor; broad residual denotes 17 remaining fine states"
    )
    p_main_b <- make_top_prob_panel(boot_all26, title = "Top-contributor probability")
    p_main_c <- make_diff_ci_panel(dcf_all26, boot_all26,
                                   title = "Progenitor minus residual contribution")
    
    fig8_dir <- file.path(dirname(dirname(OUT_ROOT)), "Figure8_source_panels")
    dir.create(fig8_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (TRUE) {
      save_pub_plot(
        p_main_a,
        file.path(fig8_dir, "Fig8B1_M1_source_donor_contribution"),
        width_mm = 110,
        height_mm = 105,
        dpi = 600
      )
      
      save_pub_plot(
        p_main_b,
        file.path(fig8_dir, "Fig8B2_top_contributor_probability"),
        width_mm = 95,
        height_mm = 58,
        dpi = 600
      )
      
      save_pub_plot(
        p_main_c,
        file.path(fig8_dir, "Fig8B3_progenitor_minus_residual_CI"),
        width_mm = 75,
        height_mm = 50,
        dpi = 600
      )
    }
    
    if (FALSE) {
      right_main <- (p_main_b / p_main_c) +
        plot_layout(heights = c(1, 0.70))
      
      fig_main <- (p_main_a | right_main) +
        plot_layout(widths = c(1.75, 1.00)) +
        plot_annotation(
          tag_levels = "A",
          caption = paste(
            "Each point is one acute donor.",
            "Bootstrap resamples donors with replacement;",
            "ties in top-source calls share probability mass."
          )
        ) &
        theme(
          plot.tag = element_text(size = 9, face = "bold", colour = "#202020"),
          plot.caption = element_text(size = 6.3, colour = "#666666", hjust = 0),
          plot.margin = margin(6, 7, 6, 7)
        )
      
      save_pub_plot(
        fig_main,
        file.path(FIG_DIR, "fig_m1_source_bootstrap_all26_nature"),
        width_mm = 183,
        height_mm = 115,
        dpi = 600
      )
    }
    
    # Sensitivity supplement: only emitted when the two NA-handling analyses
    # actually differ, so a reader is never shown two identical figures.
    if (!dcf_identical) {
      p_sens_a <- make_donor_panel(
        dcf_fig05,
        title = "Sensitivity: Fig. 05-matched handling",
        subtitle = "Donor-class values with no evaluable state remain NA"
      )
      p_sens_b <- make_top_prob_panel(boot_fig05, title = "Top-contributor probability")
      p_sens_c <- make_diff_ci_panel(dcf_fig05, boot_fig05,
                                     title = "Progenitor minus residual contribution")
      
      right_sens <- (p_sens_b / p_sens_c) +
        plot_layout(heights = c(1, 0.70))
      
      fig_sens <- (p_sens_a | right_sens) +
        plot_layout(widths = c(1.75, 1.00)) +
        plot_annotation(
          tag_levels = "A",
          caption = "Sensitivity analysis uses the same NA handling as Script 03 Fig. 05."
        ) &
        theme(
          plot.tag = element_text(size = 9, face = "bold", colour = "#202020"),
          plot.caption = element_text(size = 6.3, colour = "#666666", hjust = 0),
          plot.margin = margin(6, 7, 6, 7)
        )
      
      save_pub_plot(
        fig_sens,
        file.path(FIG_DIR, "fig_m1_source_bootstrap_fig05_matched_nature"),
        width_mm = 183, height_mm = 115, dpi = 600
      )
    } else {
      message("Skipping sensitivity figure because fig05_matched and all26_zero_absent are identical.")
    }
  }
  
  # ---------------------------------------------------------------------------
  # 7. PROVENANCE + SUMMARY
  # ---------------------------------------------------------------------------
  write_csv_atomic(data.frame(
    key = c("script","script_version","source_table","state_order_file",
            "n_boot","seed","n_donors","analyses","primary_analysis",
            "grid_check","closure_check","partition_check","fig05_all26_identical",
            "sensitivity_figure","R_version"),
    value = c("03c_M1_source_bootstrap.R","v1.7",CF_FILE,STATE_ORDER_FILE,
              as.character(N_BOOT), as.character(SEED), as.character(length(donors)),
              "fig05_matched;all26_zero_absent",
              "all26_zero_absent; fig05_matched is sensitivity analysis",
              "pass","pass","pass",
              as.character(dcf_identical),
              if (dcf_identical) "skipped (identical)" else "emitted",
              as.character(getRversion())),
    stringsAsFactors = FALSE), file.path(TABLE_DIR, "run_provenance.csv"))
  
  cat("Script 03c v1.7 complete.\nRun dir:", OUT_ROOT, "\n")
  cat("Source table:", CF_FILE, "\n")
  cat("Grid completeness: all ", length(donors), " x ", length(states),
      " donor-state pairs present | donor closure: all donors sum to 1\n", sep = "")
  cat("fig05 vs all26 analyses identical:", dcf_identical,
      " | sensitivity figure:", if (dcf_identical) "skipped (identical)" else "emitted", "\n")
  for (an in c("fig05_matched", "all26_zero_absent")) {
    dcf <- if (an == "fig05_matched") dcf_fig05 else dcf_all26
    boot <- if (an == "fig05_matched") boot_fig05 else boot_all26
    obs_med <- vapply(dcf, safe_median, numeric(1))
    n_fin <- vapply(dcf, function(v) sum(is.finite(v)), integer(1))
    cat("\n--- ", an, " ---\n", sep = "")
    cat("Observed medians:", paste(names(obs_med), format(obs_med, digits = 3), collapse = "; "), "\n")
    cat("n finite donors per class:", paste(names(n_fin), n_fin, collapse = "; "), "\n")
    cat("P(progenitor > other):", format(boot$P_prog_gt_other, digits = 3),
        " | P(=):", format(boot$P_prog_eq_other, digits = 3),
        " | P(>, ties=0.5):", format(boot$P_prog_gt_other_fractional_tie, digits = 3), "\n")
    cat("P(progenitor top):", format(boot$P_class_top[["progenitor_granulopoiesis"]], digits = 3),
        " | P(other top):", format(boot$P_class_top[["other"]], digits = 3),
        " | P(cycling top):", format(boot$P_class_top[["cycling_lymphoid"]], digits = 3), "\n")
    cat("Median diff (progenitor - other):", format(obs_med[["progenitor_granulopoiesis"]] - obs_med[["other"]],
                                                    digits = 4),
        " 95% CI [", format(boot$ci[1], digits = 4), ", ", format(boot$ci[2], digits = 4), "]\n")
    cat("Interpretation:", interpret(boot$P_prog_gt_other_fractional_tie, boot$ci), "\n")
  }
  cat("\nNote: 'broad residual' is a union of 17 fine states, not one coherent\n")
  cat("cell population; if P(residual top) is high, write 'the broad residual\n")
  cat("compartment contributed comparably', never 'residual cells dominate'.\n")
})
