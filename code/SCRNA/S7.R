# ============================================================
# Supplementary Figures S6/S7
# Panel S6A: E3-layer reference-state contribution profiles
#
# Standalone plotting script
# No Seurat object required
# No panel letter embedded in the figure
# ============================================================

options(
  width = 220,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 1. Package audit
# ------------------------------------------------------------

required_packages_s6a <- c(
  "ggplot2",
  "dplyr",
  "patchwork",
  "svglite",
  "ragg"
)

missing_packages_s6a <- required_packages_s6a[
  !vapply(
    required_packages_s6a,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_s6a) > 0L) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_packages_s6a, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ------------------------------------------------------------
# 2. Fixed project and output paths
# ------------------------------------------------------------

SC_ROOT_S6A <- "/home/sunshine/predicate/singlecell"

PROFILE_FILE_S6A <- file.path(
  SC_ROOT_S6A,
  "outputs",
  "GSE216009_32gene_analysis",
  "04_overlap_residual_source_retention",
  "20260804_overlap_residual_source_attribution_6",
  "tables",
  "fig04_profile_source_data.csv"
)

SUMMARY_FILE_S6A <- file.path(
  SC_ROOT_S6A,
  "outputs",
  "GSE216009_32gene_analysis",
  "04_overlap_residual_source_retention",
  "20260804_overlap_residual_source_attribution_6",
  "tables",
  "table1_module_set_source_summary.csv"
)

STATE_ORDER_FILE_S6A <- file.path(
  SC_ROOT_S6A,
  "outputs",
  "GSE216009_32gene_analysis",
  "01_gene_mapping_detectability",
  "20260803_gene_mapping_detectability",
  "tables",
  "fine_state_display_order_audit.csv"
)

OUTPUT_ROOT_S6A <- file.path(
  SC_ROOT_S6A,
  "outputs",
  "GSE216009_32gene_analysis",
  "manuscript_revision",
  "Supplementary_Figures_S6_S7"
)

PANEL_DIR_S6A <- file.path(
  OUTPUT_ROOT_S6A,
  "S6A_reference_state_profiles"
)

dir.create(
  PANEL_DIR_S6A,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(PROFILE_FILE_S6A),
  file.exists(SUMMARY_FILE_S6A),
  file.exists(STATE_ORDER_FILE_S6A)
)

# ------------------------------------------------------------
# 3. Read source data
# ------------------------------------------------------------

profile_raw_s6a <- read.csv(
  PROFILE_FILE_S6A,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

summary_raw_s6a <- read.csv(
  SUMMARY_FILE_S6A,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

state_order_s6a <- read.csv(
  STATE_ORDER_FILE_S6A,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_profile_columns_s6a <- c(
  "module",
  "set",
  "fine_state",
  "fine_state_short",
  "contribution"
)

required_summary_columns_s6a <- c(
  "module",
  "set",
  "n_genes",
  "main_state",
  "rho_contribution_vs_full",
  "rho_logCPM_vs_full"
)

required_state_columns_s6a <- c(
  "display_order",
  "fine_state",
  "lineage"
)

stopifnot(
  all(required_profile_columns_s6a %in% names(profile_raw_s6a)),
  all(required_summary_columns_s6a %in% names(summary_raw_s6a)),
  all(required_state_columns_s6a %in% names(state_order_s6a))
)

state_order_s6a <- state_order_s6a[
  order(state_order_s6a$display_order),
  ,
  drop = FALSE
]

stopifnot(
  nrow(state_order_s6a) == 24L,
  identical(
    state_order_s6a$display_order,
    seq_len(24L)
  ),
  !anyDuplicated(state_order_s6a$fine_state)
)

# ------------------------------------------------------------
# 4. Lock the six biological display groups
# ------------------------------------------------------------

expected_progenitor_states_s6a <- c(
  "HSPCs",
  "Cycling_neutrophil_progenitors",
  "MPO+_immature_neutrophils_or_progenitors",
  "PADI4+_immature_neutrophils"
)

stopifnot(
  identical(
    state_order_s6a$fine_state[1:4],
    expected_progenitor_states_s6a
  )
)

state_order_s6a$display_group <- c(
  rep("Progenitor/granulopoiesis", 4L),
  rep("Neutrophil", 5L),
  rep("Eosinophil/mast", 2L),
  rep("Monocyte/DC", 4L),
  rep("Lymphoid", 8L),
  "Platelets"
)

group_boundaries_s6a <- c(
  4.5,
  9.5,
  11.5,
  15.5,
  23.5
)

# ------------------------------------------------------------
# 5. Standardize module and gene-set names
# ------------------------------------------------------------

module_key_s6a <- function(x) {
  
  result <- rep(NA_character_, length(x))
  
  result[
    grepl(
      "cell",
      x,
      ignore.case = TRUE
    )
  ] <- "M1"
  
  result[
    grepl(
      "neutrophil",
      x,
      ignore.case = TRUE
    )
  ] <- "M2"
  
  result
}

set_key_s6a <- function(x) {
  
  result <- rep(NA_character_, length(x))
  
  result[
    grepl(
      "^full",
      x,
      ignore.case = TRUE
    )
  ] <- "full"
  
  result[
    grepl(
      "^overlap",
      x,
      ignore.case = TRUE
    )
  ] <- "overlap"
  
  result[
    grepl(
      "^residual",
      x,
      ignore.case = TRUE
    )
  ] <- "residual"
  
  result
}

profile_s6a <- profile_raw_s6a %>%
  mutate(
    module_key = module_key_s6a(module),
    set_key = set_key_s6a(set)
  )

summary_s6a <- summary_raw_s6a %>%
  mutate(
    module_key = module_key_s6a(module),
    set_key = set_key_s6a(set)
  )

stopifnot(
  !anyNA(profile_s6a$module_key),
  !anyNA(profile_s6a$set_key),
  !anyNA(summary_s6a$module_key),
  !anyNA(summary_s6a$set_key)
)

# ------------------------------------------------------------
# 6. State-label lookup and plotting data
# ------------------------------------------------------------

short_state_lookup_s6a <- unique(
  profile_s6a[
    ,
    c(
      "fine_state",
      "fine_state_short"
    )
  ]
)

stopifnot(
  nrow(short_state_lookup_s6a) == 24L,
  !anyDuplicated(short_state_lookup_s6a$fine_state)
)

state_order_s6a <- state_order_s6a %>%
  left_join(
    short_state_lookup_s6a,
    by = "fine_state"
  )

stopifnot(
  !anyNA(state_order_s6a$fine_state_short)
)

set_y_positions_s6a <- c(
  full = 3,
  overlap = 2,
  residual = 1
)

set_display_names_s6a <- c(
  full = "Full",
  overlap = "E3-overlap",
  residual = "E3-residual"
)

profile_s6a <- profile_s6a %>%
  left_join(
    state_order_s6a[
      ,
      c(
        "fine_state",
        "fine_state_short",
        "display_order",
        "display_group",
        "lineage"
      )
    ],
    by = c(
      "fine_state",
      "fine_state_short"
    )
  ) %>%
  mutate(
    set_y = unname(
      set_y_positions_s6a[set_key]
    )
  )

summary_s6a <- summary_s6a %>%
  mutate(
    set_y = unname(
      set_y_positions_s6a[set_key]
    ),
    set_display = unname(
      set_display_names_s6a[set_key]
    ),
    row_label = paste0(
      set_display,
      ", n=",
      n_genes
    ),
    contribution_rho_label = ifelse(
      set_key == "full",
      "\u2014",
      sprintf("%.3f", rho_contribution_vs_full)
    ),
    expression_rho_label = ifelse(
      set_key == "full",
      "\u2014",
      sprintf("%.3f", rho_logCPM_vs_full)
    )
  )

profile_s6a <- profile_s6a %>%
  left_join(
    summary_s6a[
      ,
      c(
        "module_key",
        "set_key",
        "n_genes",
        "main_state",
        "rho_contribution_vs_full",
        "rho_logCPM_vs_full",
        "row_label"
      )
    ],
    by = c(
      "module_key",
      "set_key"
    )
  )

# ------------------------------------------------------------
# 7. Numerical integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(profile_s6a) == 144L,
  nrow(summary_s6a) == 6L,
  all(table(profile_s6a$module_key, profile_s6a$set_key) == 24L),
  all(is.finite(profile_s6a$contribution)),
  all(profile_s6a$contribution >= 0)
)

expected_sizes_s6a <- data.frame(
  module_key = c(
    "M1", "M1", "M1",
    "M2", "M2", "M2"
  ),
  set_key = rep(
    c("full", "overlap", "residual"),
    2L
  ),
  expected_n = c(
    11L, 5L, 6L,
    10L, 6L, 4L
  ),
  stringsAsFactors = FALSE
)

size_audit_s6a <- expected_sizes_s6a %>%
  left_join(
    summary_s6a[
      ,
      c(
        "module_key",
        "set_key",
        "n_genes"
      )
    ],
    by = c(
      "module_key",
      "set_key"
    )
  )

stopifnot(
  all(size_audit_s6a$expected_n == size_audit_s6a$n_genes)
)

expected_main_state_s6a <-
  "MPO+_immature_neutrophils_or_progenitors"

stopifnot(
  all(
    summary_s6a$main_state ==
      expected_main_state_s6a
  )
)

m1_residual_s6a <- summary_s6a %>%
  filter(
    module_key == "M1",
    set_key == "residual"
  )

m2_residual_s6a <- summary_s6a %>%
  filter(
    module_key == "M2",
    set_key == "residual"
  )

stopifnot(
  abs(
    m1_residual_s6a$rho_contribution_vs_full -
      0.8738898
  ) < 1e-6,
  abs(
    m1_residual_s6a$rho_logCPM_vs_full -
      0.9979296
  ) < 1e-6,
  abs(
    m2_residual_s6a$rho_contribution_vs_full -
      0.957
  ) < 0.002,
  abs(
    m2_residual_s6a$rho_logCPM_vs_full -
      0.929
  ) < 0.002
)

peak_data_s6a <- profile_s6a %>%
  group_by(
    module_key,
    set_key
  ) %>%
  slice_max(
    order_by = contribution,
    n = 1L,
    with_ties = FALSE
  ) %>%
  ungroup()

stopifnot(
  nrow(peak_data_s6a) == 6L,
  all(
    peak_data_s6a$fine_state ==
      peak_data_s6a$main_state
  )
)

# ------------------------------------------------------------
# 8. Palettes
# ------------------------------------------------------------

pal_cellcycle_s6a <- c(
  "#EEF8F2",
  "#C7EBDD",
  "#86D6C9",
  "#4CB9C7",
  "#3E82C4",
  "#2C3E8F"
)

pal_neutrophil_s6a <- c(
  "#FFF1E6",
  "#FFD0B8",
  "#FFA37F",
  "#FF6F59",
  "#EF3E36",
  "#B91D3A"
)

BASE_FAMILY_S6A <- "Helvetica"

# ------------------------------------------------------------
# 9. Heatmap-block function
# ------------------------------------------------------------

make_s6a_block <- function(
    requested_module,
    palette_values,
    module_title,
    show_x_labels = TRUE
) {
  
  block_data <- profile_s6a %>%
    filter(
      module_key == requested_module
    )
  
  block_peaks <- peak_data_s6a %>%
    filter(
      module_key == requested_module
    )
  
  block_summary <- summary_s6a %>%
    filter(
      module_key == requested_module
    ) %>%
    arrange(desc(set_y))
  
  fill_upper <- max(
    block_data$contribution,
    na.rm = TRUE
  )
  
  y_breaks <- c(3, 2, 1)
  
  y_labels <- block_summary$row_label[
    match(
      y_breaks,
      block_summary$set_y
    )
  ]
  
  p_heat <- ggplot(
    block_data,
    aes(
      x = display_order,
      y = set_y,
      fill = contribution
    )
  ) +
    geom_tile(
      width = 0.96,
      height = 0.78,
      colour = "white",
      linewidth = 0.22
    ) +
    geom_vline(
      xintercept = group_boundaries_s6a,
      colour = "#666666",
      linewidth = 0.32
    ) +
    geom_point(
      data = block_peaks,
      aes(
        x = display_order,
        y = set_y
      ),
      inherit.aes = FALSE,
      shape = 21,
      fill = NA,
      colour = "#111111",
      size = 2.5,
      stroke = 0.65
    ) +
    scale_fill_gradientn(
      colours = palette_values,
      limits = c(0, fill_upper),
      trans = "sqrt",
      name = "Median donor-normalized\ncontribution"
    ) +
    scale_x_continuous(
      breaks = state_order_s6a$display_order,
      labels = state_order_s6a$fine_state_short,
      limits = c(0.5, 24.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(0.5, 3.5),
      expand = c(0, 0)
    ) +
    labs(
      x = if (show_x_labels) {
        "Reference fine state"
      } else {
        NULL
      },
      y = NULL,
      title = module_title
    ) +
    guides(
      fill = guide_colourbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0,
        barwidth = grid::unit(34, "mm"),
        barheight = grid::unit(2.3, "mm"),
        ticks = FALSE
      )
    ) +
    theme_minimal(
      base_size = 6.5,
      base_family = BASE_FAMILY_S6A
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(
        size = 6.8,
        margin = margin(t = 4)
      ),
      axis.text.y = element_text(
        size = 6.8,
        colour = "#222222",
        hjust = 1
      ),
      axis.ticks.x = element_line(
        colour = "#333333",
        linewidth = 0.3
      ),
      axis.ticks.length.x = grid::unit(1.2, "mm"),
      plot.title = element_text(
        size = 8.1,
        face = "bold",
        colour = "#111111",
        margin = margin(b = 1)
      ),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(
        size = 6.0,
        colour = "#333333"
      ),
      legend.text = element_text(
        size = 5.7,
        colour = "#333333"
      ),
      plot.margin = margin(
        t = 4,
        r = 3,
        b = 3,
        l = 4
      )
    )
  
  if (!show_x_labels) {
    p_heat <- p_heat +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  } else {
    p_heat <- p_heat +
      theme(
        axis.text.x = element_text(
          size = 5.4,
          angle = 55,
          hjust = 1,
          vjust = 1,
          colour = "#333333"
        )
      )
  }
  
  side_long <- rbind(
    data.frame(
      metric = "Contribution\nrho vs full",
      set_y = block_summary$set_y,
      display_value =
        block_summary$contribution_rho_label,
      stringsAsFactors = FALSE
    ),
    data.frame(
      metric = "Expression\nrho vs full",
      set_y = block_summary$set_y,
      display_value =
        block_summary$expression_rho_label,
      stringsAsFactors = FALSE
    )
  )
  
  side_long$metric <- factor(
    side_long$metric,
    levels = c(
      "Contribution\nrho vs full",
      "Expression\nrho vs full"
    )
  )
  
  p_side <- ggplot(
    side_long,
    aes(
      x = metric,
      y = set_y,
      label = display_value
    )
  ) +
    geom_text(
      size = 2.55,
      family = BASE_FAMILY_S6A,
      colour = "#222222"
    ) +
    scale_x_discrete(
      position = "top",
      expand = expansion(
        mult = c(0.10, 0.10)
      )
    ) +
    scale_y_continuous(
      limits = c(0.5, 3.5),
      expand = c(0, 0)
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_void(
      base_family = BASE_FAMILY_S6A
    ) +
    theme(
      axis.text.x.top = element_text(
        size = 5.8,
        face = "bold",
        colour = "#333333",
        lineheight = 0.95,
        margin = margin(b = 5)
      ),
      plot.margin = margin(
        t = 4,
        r = 4,
        b = 3,
        l = 1,
        unit = "pt"
      )
    )
  
  (p_heat | p_side) +
    plot_layout(
      widths = c(5.7, 1.35)
    )
}

# ------------------------------------------------------------
# 10. Build S6A
# ------------------------------------------------------------

block_m1_s6a <- make_s6a_block(
  requested_module = "M1",
  palette_values = pal_cellcycle_s6a,
  module_title = "M1 cell-cycle/proliferation",
  show_x_labels = FALSE
)

block_m2_s6a <- make_s6a_block(
  requested_module = "M2",
  palette_values = pal_neutrophil_s6a,
  module_title = "M2 neutrophil-degranulation",
  show_x_labels = TRUE
)

plot_s6a <- (
  block_m1_s6a /
    block_m2_s6a
) +
  plot_layout(
    heights = c(1.00, 1.23)
  ) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = "white",
        colour = NA
      )
    )
  )

# ------------------------------------------------------------
# 11. Source-data export
# ------------------------------------------------------------

source_data_s6a <- profile_s6a %>%
  mutate(
    peak_reference_state =
      fine_state == main_state
  ) %>%
  arrange(
    factor(
      module_key,
      levels = c("M1", "M2")
    ),
    desc(set_y),
    display_order
  ) %>%
  dplyr::select(
    module_key,
    module,
    set_key,
    n_genes,
    fine_state,
    fine_state_short,
    display_order,
    display_group,
    lineage,
    contribution,
    peak_reference_state,
    rho_contribution_vs_full,
    rho_logCPM_vs_full
  )

SOURCE_DATA_FILE_S6A <- file.path(
  PANEL_DIR_S6A,
  "Supplementary_Figure_S6A_source_data.csv"
)

write.csv(
  source_data_s6a,
  SOURCE_DATA_FILE_S6A,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 12. Export PDF, SVG, TIFF and PNG
# ------------------------------------------------------------

FILE_STEM_S6A <- file.path(
  PANEL_DIR_S6A,
  "Supplementary_Figure_S6A_reference_state_profiles"
)

WIDTH_MM_S6A <- 183
HEIGHT_MM_S6A <- 115
DPI_S6A <- 600

width_in_s6a <- WIDTH_MM_S6A / 25.4
height_in_s6a <- HEIGHT_MM_S6A / 25.4

# PDF: standard Helvetica-compatible device
grDevices::pdf(
  file = paste0(FILE_STEM_S6A, ".pdf"),
  width = width_in_s6a,
  height = height_in_s6a,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE
)

print(plot_s6a)

grDevices::dev.off()

# SVG
svglite::svglite(
  file = paste0(FILE_STEM_S6A, ".svg"),
  width = width_in_s6a,
  height = height_in_s6a
)

print(plot_s6a)

grDevices::dev.off()

# TIFF
ragg::agg_tiff(
  filename = paste0(FILE_STEM_S6A, ".tiff"),
  width = width_in_s6a,
  height = height_in_s6a,
  units = "in",
  res = DPI_S6A,
  background = "white",
  compression = "lzw"
)

print(plot_s6a)

grDevices::dev.off()

# PNG preview
ragg::agg_png(
  filename = paste0(FILE_STEM_S6A, ".png"),
  width = width_in_s6a,
  height = height_in_s6a,
  units = "in",
  res = 300,
  background = "white"
)

print(plot_s6a)

grDevices::dev.off()

# ------------------------------------------------------------
# 13. Final output audit
# ------------------------------------------------------------

expected_outputs_s6a <- c(
  paste0(FILE_STEM_S6A, ".pdf"),
  paste0(FILE_STEM_S6A, ".svg"),
  paste0(FILE_STEM_S6A, ".tiff"),
  paste0(FILE_STEM_S6A, ".png"),
  SOURCE_DATA_FILE_S6A
)

stopifnot(
  all(file.exists(expected_outputs_s6a)),
  all(file.info(expected_outputs_s6a)$size > 0)
)

cat(
  "\n====================================================\n",
  "S6A COMPLETED SUCCESSFULLY\n",
  "====================================================\n",
  sep = ""
)

cat(
  "Output directory:\n",
  normalizePath(
    PANEL_DIR_S6A,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

print(
  data.frame(
    file = basename(expected_outputs_s6a),
    size_bytes = file.info(expected_outputs_s6a)$size,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

cat(
  "\nPanel letter embedded: no\n",
  "Seurat object loaded: no\n",
  "Online query executed: no\n",
  sep = ""
)

# ============================================================
# Panel S6B: Leave-one-gene-out reference-state stability
#
# Each point represents removal of one gene.
# No panel letter embedded in the figure.
# ============================================================

# ------------------------------------------------------------
# 14. Fixed input and output paths
# ------------------------------------------------------------

ABLATION_SET_FILE_S6B <- file.path(
  SC_ROOT_S6A,
  "outputs",
  "GSE216009_32gene_analysis",
  "05_signature_ablation_source_retention",
  "20260804_signature_ablation_2",
  "tables",
  "ablation_set_summary.csv"
)

ABLATION_DETAIL_FILE_S6B <- file.path(
  SC_ROOT_S6A,
  "outputs",
  "GSE216009_32gene_analysis",
  "05_signature_ablation_source_retention",
  "20260804_signature_ablation_2",
  "tables",
  "ablation_source_retention_summary.csv"
)

PANEL_DIR_S6B <- file.path(
  OUTPUT_ROOT_S6A,
  "S6B_leave_one_gene_out_stability"
)

dir.create(
  PANEL_DIR_S6B,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(ABLATION_SET_FILE_S6B),
  file.exists(ABLATION_DETAIL_FILE_S6B)
)

# ------------------------------------------------------------
# 15. Read source tables
# ------------------------------------------------------------

ablation_set_raw_s6b <- read.csv(
  ABLATION_SET_FILE_S6B,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ablation_raw_s6b <- read.csv(
  ABLATION_DETAIL_FILE_S6B,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_set_columns_s6b <- c(
  "module",
  "set",
  "set_n_genes",
  "n_ablations",
  "n_top1_flips",
  "top1_retention_rate",
  "n_rho_evaluable",
  "rho_mean",
  "rho_min"
)

required_detail_columns_s6b <- c(
  "module",
  "set",
  "removed_gene",
  "set_n_genes",
  "n_genes_after",
  "top1_original",
  "top1_leave_one_out",
  "top1_flipped",
  "rho_vs_original_set",
  "rho_evaluable"
)

stopifnot(
  all(required_set_columns_s6b %in%
        names(ablation_set_raw_s6b)),
  all(required_detail_columns_s6b %in%
        names(ablation_raw_s6b))
)

# ------------------------------------------------------------
# 16. Standardize fields
# ------------------------------------------------------------

module_key_s6b <- function(x) {
  result <- rep(NA_character_, length(x))
  
  result[
    grepl("cell", x, ignore.case = TRUE)
  ] <- "M1"
  
  result[
    grepl("neutrophil", x, ignore.case = TRUE)
  ] <- "M2"
  
  result
}

set_key_s6b <- function(x) {
  result <- rep(NA_character_, length(x))
  
  result[
    grepl("^full", x, ignore.case = TRUE)
  ] <- "full"
  
  result[
    grepl("^overlap", x, ignore.case = TRUE)
  ] <- "overlap"
  
  result[
    grepl("^residual", x, ignore.case = TRUE)
  ] <- "residual"
  
  result
}

logical_s6b <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  value <- toupper(trimws(as.character(x)))
  result <- rep(NA, length(value))
  
  result[value %in% c("TRUE", "T", "1")] <- TRUE
  result[value %in% c("FALSE", "F", "0")] <- FALSE
  
  result
}

ablation_set_s6b <- ablation_set_raw_s6b %>%
  dplyr::mutate(
    module_key = module_key_s6b(module),
    set_key = set_key_s6b(set)
  )

ablation_s6b <- ablation_raw_s6b %>%
  dplyr::mutate(
    module_key = module_key_s6b(module),
    set_key = set_key_s6b(set),
    top1_flipped = logical_s6b(top1_flipped),
    rho_evaluable = logical_s6b(rho_evaluable),
    rho_vs_original_set =
      as.numeric(rho_vs_original_set)
  )

stopifnot(
  !anyNA(ablation_set_s6b$module_key),
  !anyNA(ablation_set_s6b$set_key),
  !anyNA(ablation_s6b$module_key),
  !anyNA(ablation_s6b$set_key)
)

# ------------------------------------------------------------
# 17. Fixed row order and expected results
# ------------------------------------------------------------

display_map_s6b <- data.frame(
  module_key = c(
    "M1", "M1", "M1",
    "M2", "M2", "M2"
  ),
  set_key = rep(
    c("full", "overlap", "residual"),
    2L
  ),
  row_y = c(6, 5, 4, 3, 2, 1),
  row_label = c(
    "M1 full",
    "M1 E3-overlap",
    "M1 E3-residual",
    "M2 full",
    "M2 E3-overlap",
    "M2 E3-residual"
  ),
  expected_n = c(
    11L, 5L, 6L,
    10L, 6L, 4L
  ),
  expected_min_gene = c(
    "CKS2",
    "CKS2",
    "MKI67",
    "MS4A3",
    "MS4A3",
    "RNASE3"
  ),
  expected_min_rho = c(
    0.940,
    0.896,
    0.924,
    0.957,
    0.956,
    0.992
  ),
  stringsAsFactors = FALSE
)

ablation_s6b <- ablation_s6b %>%
  dplyr::left_join(
    display_map_s6b,
    by = c(
      "module_key",
      "set_key"
    )
  )

stopifnot(
  !anyNA(ablation_s6b$row_y),
  !anyNA(ablation_s6b$row_label)
)

# ------------------------------------------------------------
# 18. Numerical integrity checks
# ------------------------------------------------------------

count_audit_s6b <- ablation_s6b %>%
  dplyr::count(
    module_key,
    set_key,
    name = "observed_n"
  ) %>%
  dplyr::left_join(
    display_map_s6b[
      ,
      c(
        "module_key",
        "set_key",
        "expected_n"
      )
    ],
    by = c(
      "module_key",
      "set_key"
    )
  )

stopifnot(
  nrow(ablation_s6b) == 42L,
  nrow(count_audit_s6b) == 6L,
  all(
    count_audit_s6b$observed_n ==
      count_audit_s6b$expected_n
  ),
  !anyNA(ablation_s6b$rho_vs_original_set),
  !anyNA(ablation_s6b$rho_evaluable),
  all(ablation_s6b$rho_evaluable),
  !anyNA(ablation_s6b$top1_flipped),
  all(!ablation_s6b$top1_flipped),
  !anyNA(ablation_s6b$top1_original),
  !anyNA(ablation_s6b$top1_leave_one_out),
  all(
    ablation_s6b$top1_original ==
      ablation_s6b$top1_leave_one_out
  )
)

summary_audit_s6b <- display_map_s6b %>%
  dplyr::left_join(
    ablation_set_s6b,
    by = c(
      "module_key",
      "set_key"
    )
  )

stopifnot(
  nrow(summary_audit_s6b) == 6L,
  !anyNA(summary_audit_s6b$n_ablations),
  all(
    summary_audit_s6b$n_ablations ==
      summary_audit_s6b$expected_n
  ),
  all(summary_audit_s6b$n_top1_flips == 0L),
  all(
    abs(
      summary_audit_s6b$top1_retention_rate - 1
    ) < 1e-12
  ),
  all(
    summary_audit_s6b$n_rho_evaluable ==
      summary_audit_s6b$expected_n
  )
)

# ------------------------------------------------------------
# 19. Plotting data
# ------------------------------------------------------------

# Deterministic vertical offsets reveal overlapping rho values.
ablation_s6b <- ablation_s6b %>%
  dplyr::group_by(
    module_key,
    set_key,
    row_y
  ) %>%
  dplyr::arrange(
    removed_gene,
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    point_y = row_y +
      seq(
        from = -0.12,
        to = 0.12,
        length.out = dplyr::n()
      )
  ) %>%
  dplyr::ungroup()

range_data_s6b <- ablation_s6b %>%
  dplyr::group_by(
    module_key,
    set_key,
    row_y,
    row_label
  ) %>%
  dplyr::summarise(
    rho_min = min(rho_vs_original_set),
    rho_max = max(rho_vs_original_set),
    .groups = "drop"
  )

minimum_data_s6b <- ablation_s6b %>%
  dplyr::group_by(
    module_key,
    set_key,
    row_y,
    row_label
  ) %>%
  dplyr::slice_min(
    order_by = rho_vs_original_set,
    n = 1L,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    minimum_label = sprintf(
      "%s  %.3f",
      removed_gene,
      rho_vs_original_set
    ),
    label_x = rho_vs_original_set +
      ifelse(
        rho_vs_original_set >= 0.95,
        -0.004,
        0.004
      ),
    label_hjust = ifelse(
      rho_vs_original_set >= 0.95,
      1,
      0
    )
  )

minimum_check_s6b <- minimum_data_s6b %>%
  dplyr::arrange(
    dplyr::desc(row_y)
  )

expected_minimum_check_s6b <- display_map_s6b %>%
  dplyr::arrange(
    dplyr::desc(row_y)
  )

stopifnot(
  identical(
    minimum_check_s6b$removed_gene,
    expected_minimum_check_s6b$expected_min_gene
  ),
  all(
    abs(
      minimum_check_s6b$rho_vs_original_set -
        expected_minimum_check_s6b$expected_min_rho
    ) < 0.0015
  )
)

retention_data_s6b <- summary_audit_s6b %>%
  dplyr::mutate(
    retention_label = paste0(
      n_ablations,
      "/",
      n_ablations,
      " retained"
    )
  )

# ------------------------------------------------------------
# 20. Plot S6B
# ------------------------------------------------------------

module_colours_s6b <- c(
  M1 = "#4E79A7",
  M2 = "#E45756"
)

set_shapes_s6b <- c(
  full = 21,
  overlap = 22,
  residual = 24
)

plot_s6b <- ggplot() +
  geom_vline(
    xintercept = 1,
    colour = "#666666",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_segment(
    data = range_data_s6b,
    aes(
      x = rho_min,
      xend = rho_max,
      y = row_y,
      yend = row_y,
      colour = module_key
    ),
    linewidth = 1.2,
    alpha = 0.42,
    show.legend = FALSE
  ) +
  geom_point(
    data = ablation_s6b,
    aes(
      x = rho_vs_original_set,
      y = point_y,
      fill = module_key,
      shape = set_key
    ),
    colour = "#202020",
    size = 2.35,
    stroke = 0.35,
    alpha = 0.92
  ) +
  geom_text(
    data = minimum_data_s6b,
    aes(
      x = label_x,
      y = row_y + 0.24,
      label = minimum_label,
      hjust = label_hjust,
      colour = module_key
    ),
    family = BASE_FAMILY_S6A,
    fontface = "bold",
    size = 2.35,
    show.legend = FALSE
  ) +
  geom_text(
    data = retention_data_s6b,
    aes(
      x = 1.022,
      y = row_y,
      label = retention_label
    ),
    family = BASE_FAMILY_S6A,
    fontface = "bold",
    size = 2.35,
    hjust = 0,
    colour = "#333333"
  ) +
  scale_fill_manual(
    values = module_colours_s6b,
    breaks = c("M1", "M2"),
    labels = c(
      "M1 cell-cycle/proliferation",
      "M2 neutrophil-degranulation"
    ),
    name = "Module"
  ) +
  scale_colour_manual(
    values = module_colours_s6b,
    guide = "none"
  ) +
  scale_shape_manual(
    values = set_shapes_s6b,
    breaks = c(
      "full",
      "overlap",
      "residual"
    ),
    labels = c(
      "Full",
      "E3-overlap",
      "E3-residual"
    ),
    name = "Gene-set definition"
  ) +
  scale_x_continuous(
    breaks = c(
      0.85,
      0.90,
      0.95,
      1.00
    ),
    labels = sprintf(
      "%.2f",
      c(
        0.85,
        0.90,
        0.95,
        1.00
      )
    )
  ) +
  scale_y_continuous(
    breaks = display_map_s6b$row_y,
    labels = display_map_s6b$row_label
  ) +
  coord_cartesian(
    xlim = c(0.85, 1.07),
    ylim = c(0.55, 6.45),
    clip = "off"
  ) +
  labs(
    title =
      "Leave-one-gene-out reference-state profile stability",
    subtitle =
      "42/42 analyses retained the original highest-contributor state",
    x = "Spearman rho versus original gene set",
    y = NULL
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 2.5
      )
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        fill = "white",
        colour = "#202020",
        size = 2.5
      )
    )
  ) +
  theme_classic(
    base_size = 6.5,
    base_family = BASE_FAMILY_S6A
  ) +
  theme(
    panel.grid.major.x = element_line(
      colour = "#E2E2E2",
      linewidth = 0.28
    ),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(
      colour = "#222222",
      linewidth = 0.4
    ),
    axis.ticks.x = element_line(
      colour = "#222222",
      linewidth = 0.35
    ),
    axis.text.y = element_text(
      size = 6.4,
      colour = "#222222",
      hjust = 1
    ),
    axis.text.x = element_text(
      size = 6.0,
      colour = "#333333"
    ),
    axis.title.x = element_text(
      size = 6.7,
      margin = margin(t = 5)
    ),
    plot.title = element_text(
      size = 8.1,
      face = "bold",
      colour = "#111111",
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 6.4,
      colour = "#444444",
      margin = margin(b = 8)
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(
      size = 6.0,
      face = "bold"
    ),
    legend.text = element_text(
      size = 5.8
    ),
    legend.key.width = grid::unit(4.5, "mm"),
    legend.spacing.x = grid::unit(2, "mm"),
    plot.margin = margin(
      t = 5,
      r = 8,
      b = 4,
      l = 5,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# 21. Export source data
# ------------------------------------------------------------

SOURCE_DATA_FILE_S6B <- file.path(
  PANEL_DIR_S6B,
  "Supplementary_Figure_S6B_source_data.csv"
)

SET_SUMMARY_FILE_S6B <- file.path(
  PANEL_DIR_S6B,
  "Supplementary_Figure_S6B_set_summary.csv"
)

write.csv(
  ablation_s6b,
  SOURCE_DATA_FILE_S6B,
  row.names = FALSE,
  na = ""
)

write.csv(
  retention_data_s6b,
  SET_SUMMARY_FILE_S6B,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 22. Export figure
# ------------------------------------------------------------

FILE_STEM_S6B <- file.path(
  PANEL_DIR_S6B,
  "Supplementary_Figure_S6B_leave_one_gene_out_stability"
)

WIDTH_MM_S6B <- 160
HEIGHT_MM_S6B <- 95
DPI_S6B <- 600

width_in_s6b <- WIDTH_MM_S6B / 25.4
height_in_s6b <- HEIGHT_MM_S6B / 25.4

grDevices::pdf(
  file = paste0(FILE_STEM_S6B, ".pdf"),
  width = width_in_s6b,
  height = height_in_s6b,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE
)

print(plot_s6b)
grDevices::dev.off()

svglite::svglite(
  file = paste0(FILE_STEM_S6B, ".svg"),
  width = width_in_s6b,
  height = height_in_s6b
)

print(plot_s6b)
grDevices::dev.off()

ragg::agg_tiff(
  filename = paste0(FILE_STEM_S6B, ".tiff"),
  width = width_in_s6b,
  height = height_in_s6b,
  units = "in",
  res = DPI_S6B,
  background = "white",
  compression = "lzw"
)

print(plot_s6b)
grDevices::dev.off()

ragg::agg_png(
  filename = paste0(FILE_STEM_S6B, ".png"),
  width = width_in_s6b,
  height = height_in_s6b,
  units = "in",
  res = 300,
  background = "white"
)

print(plot_s6b)
grDevices::dev.off()

# ------------------------------------------------------------
# 23. Final output audit
# ------------------------------------------------------------

expected_outputs_s6b <- c(
  paste0(FILE_STEM_S6B, ".pdf"),
  paste0(FILE_STEM_S6B, ".svg"),
  paste0(FILE_STEM_S6B, ".tiff"),
  paste0(FILE_STEM_S6B, ".png"),
  SOURCE_DATA_FILE_S6B,
  SET_SUMMARY_FILE_S6B
)

stopifnot(
  all(file.exists(expected_outputs_s6b)),
  all(file.info(expected_outputs_s6b)$size > 0)
)

cat(
  "\n====================================================\n",
  "S6B COMPLETED SUCCESSFULLY\n",
  "====================================================\n",
  sep = ""
)

cat(
  "Output directory:\n",
  normalizePath(
    PANEL_DIR_S6B,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

print(
  data.frame(
    file = basename(expected_outputs_s6b),
    size_bytes = file.info(expected_outputs_s6b)$size,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

cat(
  "\nPanel letter embedded: no\n",
  "Upstream single-cell analysis rerun: no\n",
  "Online query executed: no\n",
  sep = ""
)


# ============================================================
# Panel S6C: High doublet-score cell exclusion sensitivity
#
# Left: exclusion-induced change in the M1 contribution contrast
# Right: top-contributor bootstrap probability
# No panel letter embedded
# ============================================================

# ------------------------------------------------------------
# 24. Fixed input and output paths
# ------------------------------------------------------------

QC_TABLE_DIR_S6C <- file.path(
  SC_ROOT_S6A,
  "outputs",
  "GSE216009_32gene_analysis",
  "qc_ambient_followup",
  "20260813_ambient_followup_TEST_2",
  "tables"
)

CONTRIBUTION_FILE_S6C <- file.path(
  QC_TABLE_DIR_S6C,
  "m1_contribution_groups.csv"
)

PROBABILITY_FILE_S6C <- file.path(
  QC_TABLE_DIR_S6C,
  "m1_top_compartment_support.csv"
)

PANEL_DIR_S6C <- file.path(
  OUTPUT_ROOT_S6A,
  "S6C_doublet_score_exclusion_sensitivity"
)

dir.create(
  PANEL_DIR_S6C,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(CONTRIBUTION_FILE_S6C),
  file.exists(PROBABILITY_FILE_S6C)
)

# ------------------------------------------------------------
# 25. Read source data
# ------------------------------------------------------------

contribution_raw_s6c <- read.csv(
  CONTRIBUTION_FILE_S6C,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

probability_raw_s6c <- read.csv(
  PROBABILITY_FILE_S6C,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_contribution_columns_s6c <- c(
  "group",
  "n_donors_evaluable_M1",
  "diff_prog_minus_resid_common_donors",
  "paired_diff_median",
  "paired_diff_ci_lo",
  "paired_diff_ci_hi",
  "paired_diff_contains_zero",
  "rank_changed_vs_full"
)

required_probability_columns_s6c <- c(
  "group",
  "progenitor_granulopoiesis",
  "broad_residual"
)

stopifnot(
  all(required_contribution_columns_s6c %in%
        names(contribution_raw_s6c)),
  all(required_probability_columns_s6c %in%
        names(probability_raw_s6c))
)

# ------------------------------------------------------------
# 26. Fixed group order
# ------------------------------------------------------------

group_order_s6c <- c(
  "full",
  "per_sample_top5",
  "per_sample_top10",
  "global_acute_primary_top5",
  "global_acute_primary_top10"
)

group_display_s6c <- c(
  full = "Full data",
  per_sample_top5 = "Per-sample top 5%",
  per_sample_top10 = "Per-sample top 10%",
  global_acute_primary_top5 = "Global acute top 5%",
  global_acute_primary_top10 = "Global acute top 10%"
)

strategy_colours_s6c <- c(
  full = "#777777",
  per_sample_top5 = "#B8CCE4",
  per_sample_top10 = "#84A9D2",
  global_acute_primary_top5 = "#527FAF",
  global_acute_primary_top10 = "#2F5F8F"
)

logical_s6c <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  value <- toupper(trimws(as.character(x)))
  result <- rep(NA, length(value))
  
  result[value %in% c("TRUE", "T", "1")] <- TRUE
  result[value %in% c("FALSE", "F", "0")] <- FALSE
  
  result
}

contribution_s6c <- contribution_raw_s6c %>%
  dplyr::mutate(
    paired_diff_contains_zero =
      logical_s6c(paired_diff_contains_zero),
    rank_changed_vs_full =
      logical_s6c(rank_changed_vs_full)
  )

merged_s6c <- contribution_s6c %>%
  dplyr::left_join(
    probability_raw_s6c,
    by = "group"
  ) %>%
  dplyr::mutate(
    group_display =
      unname(group_display_s6c[group]),
    full_contrast_pp =
      100 * diff_prog_minus_resid_common_donors,
    median_change_pp =
      100 * paired_diff_median,
    ci_low_pp =
      100 * paired_diff_ci_lo,
    ci_high_pp =
      100 * paired_diff_ci_hi,
    progenitor_probability_pct =
      100 * progenitor_granulopoiesis
  )

merged_s6c <- merged_s6c[
  match(
    group_order_s6c,
    merged_s6c$group
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 27. Numerical integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(contribution_s6c) == 5L,
  nrow(probability_raw_s6c) == 5L,
  identical(
    merged_s6c$group,
    group_order_s6c
  ),
  !anyNA(merged_s6c$group_display),
  all(merged_s6c$n_donors_evaluable_M1 == 26L),
  all(
    abs(
      merged_s6c$progenitor_granulopoiesis +
        merged_s6c$broad_residual - 1
    ) < 1e-10
  )
)

full_row_s6c <- merged_s6c[
  merged_s6c$group == "full",
  ,
  drop = FALSE
]

exclusion_rows_s6c <- merged_s6c[
  merged_s6c$group != "full",
  ,
  drop = FALSE
]

stopifnot(
  nrow(full_row_s6c) == 1L,
  nrow(exclusion_rows_s6c) == 4L,
  abs(
    full_row_s6c$full_contrast_pp -
      (-0.24213929)
  ) < 1e-5,
  all(
    is.finite(
      exclusion_rows_s6c$median_change_pp
    )
  ),
  all(
    is.finite(
      exclusion_rows_s6c$ci_low_pp
    )
  ),
  all(
    is.finite(
      exclusion_rows_s6c$ci_high_pp
    )
  ),
  all(
    exclusion_rows_s6c$paired_diff_contains_zero
  ),
  all(
    exclusion_rows_s6c$rank_changed_vs_full
  )
)

expected_probability_s6c <- c(
  46.44,
  48.07,
  52.09,
  48.53,
  49.83
)

stopifnot(
  all(
    abs(
      merged_s6c$progenitor_probability_pct -
        expected_probability_s6c
    ) < 0.001
  )
)

# ------------------------------------------------------------
# 28. Prepare panel data
# ------------------------------------------------------------

left_order_s6c <- c(
  "per_sample_top5",
  "per_sample_top10",
  "global_acute_primary_top5",
  "global_acute_primary_top10"
)

left_map_s6c <- data.frame(
  group = left_order_s6c,
  row_y = c(4, 3, 2, 1),
  stringsAsFactors = FALSE
)

left_data_s6c <- exclusion_rows_s6c %>%
  dplyr::left_join(
    left_map_s6c,
    by = "group"
  )

right_map_s6c <- data.frame(
  group = group_order_s6c,
  row_y = c(5, 4, 3, 2, 1),
  stringsAsFactors = FALSE
)

right_data_s6c <- merged_s6c %>%
  dplyr::left_join(
    right_map_s6c,
    by = "group"
  ) %>%
  dplyr::mutate(
    probability_label = sprintf(
      "%.1f%%",
      progenitor_probability_pct
    ),
    probability_label_x =
      progenitor_probability_pct +
      ifelse(
        progenitor_probability_pct >= 50,
        3,
        -3
      ),
    probability_label_hjust =
      ifelse(
        progenitor_probability_pct >= 50,
        0,
        1
      )
  )

stopifnot(
  !anyNA(left_data_s6c$row_y),
  !anyNA(right_data_s6c$row_y)
)

# ------------------------------------------------------------
# 29. Left panel: paired bootstrap change
# ------------------------------------------------------------

plot_left_s6c <- ggplot(
  left_data_s6c,
  aes(
    y = row_y,
    colour = group
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "#666666",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      x = ci_low_pp,
      xend = ci_high_pp,
      yend = row_y
    ),
    linewidth = 0.9
  ) +
  geom_segment(
    aes(
      x = ci_low_pp,
      xend = ci_low_pp,
      y = row_y - 0.09,
      yend = row_y + 0.09
    ),
    linewidth = 0.6
  ) +
  geom_segment(
    aes(
      x = ci_high_pp,
      xend = ci_high_pp,
      y = row_y - 0.09,
      yend = row_y + 0.09
    ),
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      x = median_change_pp
    ),
    size = 2.8,
    shape = 21,
    fill = "white",
    stroke = 0.9
  ) +
  scale_colour_manual(
    values = strategy_colours_s6c,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(
      -8,
      -4,
      0,
      4,
      8
    ),
    labels = c(
      "-8",
      "-4",
      "0",
      "4",
      "8"
    )
  ) +
  scale_y_continuous(
    breaks = left_map_s6c$row_y,
    labels = unname(
      group_display_s6c[
        left_map_s6c$group
      ]
    )
  ) +
  coord_cartesian(
    xlim = c(-8, 8),
    ylim = c(0.55, 4.45),
    clip = "off"
  ) +
  labs(
    title =
      "Exclusion-induced change in the M1 contribution contrast",
    subtitle = sprintf(
      "Full-data contrast: %.2f percentage points",
      full_row_s6c$full_contrast_pp
    ),
    x = paste0(
      "Change in progenitor-minus-residual contribution ",
      "versus full\n(percentage points)"
    ),
    y = NULL
  ) +
  theme_classic(
    base_size = 6.5,
    base_family = BASE_FAMILY_S6A
  ) +
  theme(
    panel.grid.major.x = element_line(
      colour = "#E4E4E4",
      linewidth = 0.28
    ),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      size = 6.2,
      colour = "#222222"
    ),
    axis.text.x = element_text(
      size = 5.9,
      colour = "#333333"
    ),
    axis.title.x = element_text(
      size = 6.3,
      lineheight = 1.05,
      margin = margin(t = 5)
    ),
    plot.title = element_text(
      size = 7.1,
      face = "bold",
      colour = "#111111",
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 6.0,
      colour = "#444444",
      margin = margin(b = 7)
    ),
    plot.margin = margin(
      t = 4,
      r = 5,
      b = 4,
      l = 3,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# 30. Right panel: top-contributor probability
# ------------------------------------------------------------

plot_right_s6c <- ggplot(
  right_data_s6c,
  aes(
    y = row_y,
    colour = group
  )
) +
  geom_vline(
    xintercept = 50,
    colour = "#666666",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      x = progenitor_probability_pct
    ),
    size = 3.0
  ) +
  geom_text(
    aes(
      x = probability_label_x,
      label = probability_label,
      hjust = probability_label_hjust
    ),
    family = BASE_FAMILY_S6A,
    fontface = "bold",
    size = 2.35,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = strategy_colours_s6c,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    labels = c(
      "0",
      "25",
      "50",
      "75",
      "100"
    )
  ) +
  scale_y_continuous(
    breaks = right_map_s6c$row_y,
    labels = unname(
      group_display_s6c[
        right_map_s6c$group
      ]
    )
  ) +
  coord_cartesian(
    xlim = c(0, 100),
    ylim = c(0.55, 5.45),
    clip = "off"
  ) +
  labs(
    title = "Top-contributor probability",
    subtitle = paste0(
      "P(progenitor/granulopoiesis ",
      "> broad residual)"
    ),
    x = "Bootstrap probability (%)",
    y = NULL
  ) +
  theme_classic(
    base_size = 6.5,
    base_family = BASE_FAMILY_S6A
  ) +
  theme(
    panel.grid.major.x = element_line(
      colour = "#E4E4E4",
      linewidth = 0.28
    ),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      size = 6.2,
      colour = "#222222"
    ),
    axis.text.x = element_text(
      size = 5.9,
      colour = "#333333"
    ),
    axis.title.x = element_text(
      size = 6.3,
      margin = margin(t = 5)
    ),
    plot.title = element_text(
      size = 7.1,
      face = "bold",
      colour = "#111111",
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 6.0,
      colour = "#444444",
      margin = margin(b = 7)
    ),
    plot.margin = margin(
      t = 4,
      r = 6,
      b = 4,
      l = 5,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# 31. Assemble S6C
# ------------------------------------------------------------

plot_s6c <- (
  plot_left_s6c |
    plot_right_s6c
) +
  patchwork::plot_layout(
    widths = c(1.15, 1)
  ) +
  patchwork::plot_annotation(
    title = "M1 contribution-group sensitivity",
    subtitle =
      "Sensitivity analyses: high doublet-score cells removed",
    theme = theme(
      plot.title = element_text(
        family = BASE_FAMILY_S6A,
        size = 8.3,
        face = "bold",
        colour = "#111111"
      ),
      plot.subtitle = element_text(
        family = BASE_FAMILY_S6A,
        size = 6.4,
        colour = "#444444",
        margin = margin(b = 6)
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      )
    )
  )

# ------------------------------------------------------------
# 32. Export source data
# ------------------------------------------------------------

SOURCE_DATA_FILE_S6C <- file.path(
  PANEL_DIR_S6C,
  "Supplementary_Figure_S6C_source_data.csv"
)

write.csv(
  merged_s6c,
  SOURCE_DATA_FILE_S6C,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 33. Export PDF, SVG, TIFF and PNG
# ------------------------------------------------------------

FILE_STEM_S6C <- file.path(
  PANEL_DIR_S6C,
  "Supplementary_Figure_S6C_doublet_score_exclusion_sensitivity"
)

WIDTH_MM_S6C <- 183
HEIGHT_MM_S6C <- 88
DPI_S6C <- 600

width_in_s6c <- WIDTH_MM_S6C / 25.4
height_in_s6c <- HEIGHT_MM_S6C / 25.4

grDevices::pdf(
  file = paste0(FILE_STEM_S6C, ".pdf"),
  width = width_in_s6c,
  height = height_in_s6c,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE
)

print(plot_s6c)
grDevices::dev.off()

svglite::svglite(
  file = paste0(FILE_STEM_S6C, ".svg"),
  width = width_in_s6c,
  height = height_in_s6c
)

print(plot_s6c)
grDevices::dev.off()

ragg::agg_tiff(
  filename = paste0(FILE_STEM_S6C, ".tiff"),
  width = width_in_s6c,
  height = height_in_s6c,
  units = "in",
  res = DPI_S6C,
  background = "white",
  compression = "lzw"
)

print(plot_s6c)
grDevices::dev.off()

ragg::agg_png(
  filename = paste0(FILE_STEM_S6C, ".png"),
  width = width_in_s6c,
  height = height_in_s6c,
  units = "in",
  res = 300,
  background = "white"
)

print(plot_s6c)
grDevices::dev.off()

# ------------------------------------------------------------
# 34. Final output audit
# ------------------------------------------------------------

expected_outputs_s6c <- c(
  paste0(FILE_STEM_S6C, ".pdf"),
  paste0(FILE_STEM_S6C, ".svg"),
  paste0(FILE_STEM_S6C, ".tiff"),
  paste0(FILE_STEM_S6C, ".png"),
  SOURCE_DATA_FILE_S6C
)

stopifnot(
  all(file.exists(expected_outputs_s6c)),
  all(file.info(expected_outputs_s6c)$size > 0)
)

cat(
  "\n====================================================\n",
  "S6C COMPLETED SUCCESSFULLY\n",
  "====================================================\n",
  sep = ""
)

cat(
  "Output directory:\n",
  normalizePath(
    PANEL_DIR_S6C,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

print(
  data.frame(
    file = basename(expected_outputs_s6c),
    size_bytes = file.info(expected_outputs_s6c)$size,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

cat(
  "\nPanel letter embedded: no\n",
  "High doublet-score cells, not confirmed doublets: yes\n",
  "Upstream single-cell analysis rerun: no\n",
  "Online query executed: no\n",
  sep = ""
)



# ============================================================
# Supplementary Figure S6D
# Donor-level granule-transcript correlations
# No embedded panel letter
# ============================================================

cat(
  "\n========== S6D: GRANULE CORRELATION FOREST ==========\n"
)

# ------------------------------------------------------------
# 1. Runtime and package checks
# ------------------------------------------------------------

required_packages_s6d <- c(
  "ggplot2",
  "svglite",
  "ragg"
)

missing_packages_s6d <- required_packages_s6d[
  !vapply(
    required_packages_s6d,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_s6d) > 0L) {
  stop(
    paste0(
      "Missing required packages for S6D: ",
      paste(missing_packages_s6d, collapse = ", ")
    )
  )
}

sc_root_s6d <- if (exists("SC_ROOT_S6A", inherits = TRUE)) {
  SC_ROOT_S6A
} else {
  "/home/sunshine/predicate/singlecell"
}

output_root_s6d <- if (
  exists("OUTPUT_ROOT_S6A", inherits = TRUE)
) {
  OUTPUT_ROOT_S6A
} else {
  file.path(
    sc_root_s6d,
    "manuscript_revision",
    "Supplementary_Figures_S6_S7"
  )
}

base_family_s6d <- if (
  exists("BASE_FAMILY_S6A", inherits = TRUE)
) {
  BASE_FAMILY_S6A
} else {
  "Helvetica"
}

# ------------------------------------------------------------
# 2. Fixed input files
# ------------------------------------------------------------

input_dir_s6d <- file.path(
  sc_root_s6d,
  "outputs",
  "GSE216009_32gene_analysis",
  "qc_ambient_followup",
  "20260813_ambient_followup_TEST_2",
  "tables"
)

donor_signal_file_s6d <- file.path(
  input_dir_s6d,
  "ambient_doublet_sensitivity_donor_signals_long.csv"
)

jackknife_file_s6d <- file.path(
  input_dir_s6d,
  "donor_jackknife_summary.csv"
)

required_files_s6d <- c(
  donor_signal_file_s6d,
  jackknife_file_s6d
)

file_audit_s6d <- data.frame(
  file = required_files_s6d,
  exists = file.exists(required_files_s6d),
  stringsAsFactors = FALSE
)

print(file_audit_s6d, row.names = FALSE)

if (!all(file_audit_s6d$exists)) {
  stop(
    paste0(
      "S6D input file(s) missing:\n",
      paste(
        file_audit_s6d$file[
          !file_audit_s6d$exists
        ],
        collapse = "\n"
      )
    )
  )
}

# ------------------------------------------------------------
# 3. Output directory
# ------------------------------------------------------------

output_dir_s6d <- file.path(
  output_root_s6d,
  "S6D_granule_correlation_forest"
)

dir.create(
  output_dir_s6d,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 4. Read source tables
# ------------------------------------------------------------

donor_signal_s6d <- read.csv(
  donor_signal_file_s6d,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

jackknife_s6d <- read.csv(
  jackknife_file_s6d,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_donor_columns_s6d <- c(
  "donor",
  "panel",
  "scope",
  "neutrophil_fraction_full",
  "ambient_signal",
  "in_correlation"
)

missing_donor_columns_s6d <- setdiff(
  required_donor_columns_s6d,
  colnames(donor_signal_s6d)
)

if (length(missing_donor_columns_s6d) > 0L) {
  stop(
    paste0(
      "Missing columns in donor-level signal table: ",
      paste(missing_donor_columns_s6d, collapse = ", ")
    )
  )
}

required_jackknife_columns_s6d <- c(
  "panel",
  "n_donors",
  "rho_full",
  "rho_min_loo",
  "rho_max_loo"
)

missing_jackknife_columns_s6d <- setdiff(
  required_jackknife_columns_s6d,
  colnames(jackknife_s6d)
)

if (length(missing_jackknife_columns_s6d) > 0L) {
  stop(
    paste0(
      "Missing columns in jackknife table: ",
      paste(missing_jackknife_columns_s6d, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 5. Normalize logical and numeric fields
# ------------------------------------------------------------

as_logical_s6d <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  x_character <- tolower(trimws(as.character(x)))
  
  output <- rep(
    NA,
    length(x_character)
  )
  
  output[
    x_character %in% c(
      "true",
      "t",
      "1",
      "yes"
    )
  ] <- TRUE
  
  output[
    x_character %in% c(
      "false",
      "f",
      "0",
      "no"
    )
  ] <- FALSE
  
  output
}

donor_signal_s6d$in_correlation <- as_logical_s6d(
  donor_signal_s6d$in_correlation
)

donor_signal_s6d$neutrophil_fraction_full <- suppressWarnings(
  as.numeric(
    donor_signal_s6d$neutrophil_fraction_full
  )
)

donor_signal_s6d$ambient_signal <- suppressWarnings(
  as.numeric(
    donor_signal_s6d$ambient_signal
  )
)

# ------------------------------------------------------------
# 6. Define the formal four-panel family
#
# This order is retained during bootstrap reconstruction
# because the original random seed was set once before
# iterating through the panel family.
# ------------------------------------------------------------

formal_panel_order_s6d <- c(
  "primary_azurophilic",
  "secondary_specific",
  "tertiary_gelatinase",
  "residual_E3_granulocyte_panel"
)

analysis_data_s6d <- donor_signal_s6d[
  donor_signal_s6d$scope == "full" &
    donor_signal_s6d$panel %in%
    formal_panel_order_s6d &
    donor_signal_s6d$in_correlation %in% TRUE &
    is.finite(
      donor_signal_s6d$neutrophil_fraction_full
    ) &
    is.finite(
      donor_signal_s6d$ambient_signal
    ),
  ,
  drop = FALSE
]

panel_counts_s6d <- table(
  factor(
    analysis_data_s6d$panel,
    levels = formal_panel_order_s6d
  )
)

cat("\nEligible donor counts by panel:\n")
print(panel_counts_s6d)

stopifnot(
  identical(
    as.integer(panel_counts_s6d),
    rep(21L, 4L)
  )
)

unique_donor_check_s6d <- vapply(
  formal_panel_order_s6d,
  function(current_panel) {
    
    current_donors <- analysis_data_s6d$donor[
      analysis_data_s6d$panel == current_panel
    ]
    
    length(current_donors) ==
      length(unique(current_donors))
  },
  logical(1)
)

stopifnot(all(unique_donor_check_s6d))

# ------------------------------------------------------------
# 7. Load frozen full-data correlation statistics
#    and verify them against donor-level data
# ------------------------------------------------------------

correlation_source_file_s6d <- file.path(
  sc_root_s6d,
  "outputs",
  "GSE216009_32gene_analysis",
  "qc_sample_summary_ambient",
  "20260813_singlecell_qc_TEST_5",
  "tables",
  "ambient_vs_neutrophil_fraction_by_panel.csv"
)

if (!file.exists(correlation_source_file_s6d)) {
  stop(
    paste0(
      "Frozen S6D correlation source file was not found:\n",
      correlation_source_file_s6d
    )
  )
}

formal_correlation_source_s6d <- read.csv(
  correlation_source_file_s6d,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_correlation_columns_s6d <- c(
  "panel",
  "n_donors_with_data",
  "rho",
  "ci_low",
  "ci_high",
  "p_value",
  "padj"
)

missing_correlation_columns_s6d <- setdiff(
  required_correlation_columns_s6d,
  colnames(formal_correlation_source_s6d)
)

if (length(missing_correlation_columns_s6d) > 0L) {
  stop(
    paste0(
      "Missing columns in frozen S6D correlation source: ",
      paste(
        missing_correlation_columns_s6d,
        collapse = ", "
      )
    )
  )
}

formal_correlation_source_s6d$panel[
  formal_correlation_source_s6d$panel ==
    "residual_neutrophil_signature"
] <- "residual_E3_granulocyte_panel"

formal_row_match_s6d <- match(
  formal_panel_order_s6d,
  formal_correlation_source_s6d$panel
)

if (anyNA(formal_row_match_s6d)) {
  stop(
    paste0(
      "The frozen correlation source does not contain ",
      "all four formal panels."
    )
  )
}

formal_correlation_source_s6d <-
  formal_correlation_source_s6d[
    formal_row_match_s6d,
    required_correlation_columns_s6d,
    drop = FALSE
  ]

numeric_source_columns_s6d <- c(
  "n_donors_with_data",
  "rho",
  "ci_low",
  "ci_high",
  "p_value",
  "padj"
)

formal_correlation_source_s6d[
  numeric_source_columns_s6d
] <- lapply(
  formal_correlation_source_s6d[
    numeric_source_columns_s6d
  ],
  as.numeric
)

stopifnot(
  identical(
    formal_correlation_source_s6d$panel,
    formal_panel_order_s6d
  ),
  all(
    formal_correlation_source_s6d$
      n_donors_with_data == 21L
  ),
  all(
    is.finite(
      as.matrix(
        formal_correlation_source_s6d[
          ,
          numeric_source_columns_s6d,
          drop = FALSE
        ]
      )
    )
  )
)

correlation_summary_s6d <- data.frame(
  panel =
    formal_correlation_source_s6d$panel,
  n_donors = as.integer(
    formal_correlation_source_s6d$
      n_donors_with_data
  ),
  rho_full =
    formal_correlation_source_s6d$rho,
  P_value =
    formal_correlation_source_s6d$p_value,
  bootstrap_CI_low =
    formal_correlation_source_s6d$ci_low,
  bootstrap_CI_high =
    formal_correlation_source_s6d$ci_high,
  FDR =
    formal_correlation_source_s6d$padj,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 7.1 Donor-level verification of rho, P and FDR
#
# Bootstrap confidence intervals are not recalculated here.
# They are read from the frozen main-QC summary above.
# ------------------------------------------------------------

donor_audit_rows_s6d <- lapply(
  formal_panel_order_s6d,
  function(current_panel_s6d) {
    
    current_data_s6d <- analysis_data_s6d[
      analysis_data_s6d$panel ==
        current_panel_s6d,
      ,
      drop = FALSE
    ]
    
    current_test_s6d <- suppressWarnings(
      stats::cor.test(
        x = current_data_s6d$
          neutrophil_fraction_full,
        y = current_data_s6d$
          ambient_signal,
        method = "spearman",
        exact = FALSE
      )
    )
    
    data.frame(
      panel = current_panel_s6d,
      n_donors_audit =
        nrow(current_data_s6d),
      rho_audit =
        unname(current_test_s6d$estimate),
      P_value_audit =
        current_test_s6d$p.value,
      stringsAsFactors = FALSE
    )
  }
)

donor_audit_s6d <- do.call(
  rbind,
  donor_audit_rows_s6d
)

rownames(donor_audit_s6d) <- NULL

donor_audit_s6d$FDR_audit <-
  stats::p.adjust(
    donor_audit_s6d$P_value_audit,
    method = "BH",
    n = 4L
  )

donor_audit_match_s6d <- match(
  correlation_summary_s6d$panel,
  donor_audit_s6d$panel
)

stopifnot(
  !anyNA(donor_audit_match_s6d),
  identical(
    correlation_summary_s6d$n_donors,
    as.integer(
      donor_audit_s6d$n_donors_audit[
        donor_audit_match_s6d
      ]
    )
  ),
  isTRUE(
    all.equal(
      correlation_summary_s6d$rho_full,
      donor_audit_s6d$rho_audit[
        donor_audit_match_s6d
      ],
      tolerance = 1e-10,
      check.attributes = FALSE
    )
  ),
  isTRUE(
    all.equal(
      correlation_summary_s6d$P_value,
      donor_audit_s6d$P_value_audit[
        donor_audit_match_s6d
      ],
      tolerance = 1e-10,
      check.attributes = FALSE
    )
  ),
  isTRUE(
    all.equal(
      correlation_summary_s6d$FDR,
      donor_audit_s6d$FDR_audit[
        donor_audit_match_s6d
      ],
      tolerance = 1e-10,
      check.attributes = FALSE
    )
  )
)

cat(
  "\nFrozen S6D correlation statistics:\n"
)

print(
  correlation_summary_s6d,
  row.names = FALSE,
  digits = 9
)

cat(
  "\nPASS: donor-level rho, P value and FDR ",
  "reproduced the frozen main-QC summary.\n",
  sep = ""
)

cat(
  "Bootstrap confidence intervals were read from:\n",
  normalizePath(
    correlation_source_file_s6d,
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
# ------------------------------------------------------------
# 8. Add leave-one-donor-out ranges
# ------------------------------------------------------------

jackknife_match_s6d <- match(
  correlation_summary_s6d$panel,
  jackknife_s6d$panel
)

if (anyNA(jackknife_match_s6d)) {
  stop(
    "One or more formal panels were absent from donor_jackknife_summary.csv."
  )
}

correlation_summary_s6d$rho_min_LODO <-
  as.numeric(
    jackknife_s6d$rho_min_loo[
      jackknife_match_s6d
    ]
  )

correlation_summary_s6d$rho_max_LODO <-
  as.numeric(
    jackknife_s6d$rho_max_loo[
      jackknife_match_s6d
    ]
  )

correlation_summary_s6d$
  jackknife_n_donors <-
  as.integer(
    jackknife_s6d$n_donors[
      jackknife_match_s6d
    ]
  )

stopifnot(
  identical(
    as.integer(
      correlation_summary_s6d$n_donors
    ),
    correlation_summary_s6d$
      jackknife_n_donors
  )
)

# ------------------------------------------------------------
# 9. Exact-value audit
# ------------------------------------------------------------

expected_statistics_s6d <- data.frame(
  panel = formal_panel_order_s6d,
  rho_full_expected = c(
    0.523376623,
    0.024707433,
    0.293506494,
    0.677922078
  ),
  FDR_expected = c(
    0.029800013,
    0.915339285,
    0.262123356,
    0.002930127
  ),
  bootstrap_CI_low_expected = c(
    0.090,
    -0.430,
    -0.234,
    0.345
  ),
  bootstrap_CI_high_expected = c(
    0.815,
    0.433,
    0.716,
    0.851
  ),
  rho_min_LODO_expected = c(
    0.448120301,
    -0.078313342,
    0.181954887,
    0.627067669
  ),
  rho_max_LODO_expected = c(
    0.607518800,
    0.118976040,
    0.434586470,
    0.723308270
  ),
  stringsAsFactors = FALSE
)

expected_match_s6d <- match(
  correlation_summary_s6d$panel,
  expected_statistics_s6d$panel
)

stopifnot(
  max(
    abs(
      correlation_summary_s6d$rho_full -
        expected_statistics_s6d$
        rho_full_expected[
          expected_match_s6d
        ]
    )
  ) < 1e-6,
  max(
    abs(
      correlation_summary_s6d$FDR -
        expected_statistics_s6d$
        FDR_expected[
          expected_match_s6d
        ]
    )
  ) < 1e-6,
  max(
    abs(
      correlation_summary_s6d$rho_min_LODO -
        expected_statistics_s6d$
        rho_min_LODO_expected[
          expected_match_s6d
        ]
    )
  ) < 1e-6,
  max(
    abs(
      correlation_summary_s6d$rho_max_LODO -
        expected_statistics_s6d$
        rho_max_LODO_expected[
          expected_match_s6d
        ]
    )
  ) < 1e-6
)

bootstrap_difference_s6d <- max(
  c(
    abs(
      correlation_summary_s6d$
        bootstrap_CI_low -
        expected_statistics_s6d$
        bootstrap_CI_low_expected[
          expected_match_s6d
        ]
    ),
    abs(
      correlation_summary_s6d$
        bootstrap_CI_high -
        expected_statistics_s6d$
        bootstrap_CI_high_expected[
          expected_match_s6d
        ]
    )
  )
)

if (bootstrap_difference_s6d > 0.006) {
  
  bootstrap_audit_s6d <- data.frame(
    panel = correlation_summary_s6d$panel,
    observed_low =
      correlation_summary_s6d$
      bootstrap_CI_low,
    expected_low =
      expected_statistics_s6d$
      bootstrap_CI_low_expected[
        expected_match_s6d
      ],
    observed_high =
      correlation_summary_s6d$
      bootstrap_CI_high,
    expected_high =
      expected_statistics_s6d$
      bootstrap_CI_high_expected[
        expected_match_s6d
      ]
  )
  
  print(
    bootstrap_audit_s6d,
    row.names = FALSE,
    digits = 6
  )
  
  stop(
    paste0(
      "Reconstructed bootstrap confidence intervals ",
      "did not reproduce the validated S6D values. ",
      "Do not draw the figure until this discrepancy ",
      "is resolved."
    )
  )
}

cat("\nValidated S6D statistics:\n")

print(
  correlation_summary_s6d,
  row.names = FALSE,
  digits = 6
)

# ------------------------------------------------------------
# 10. Display structure
# ------------------------------------------------------------

display_order_s6d <- c(
  "residual_E3_granulocyte_panel",
  "primary_azurophilic",
  "tertiary_gelatinase",
  "secondary_specific"
)

display_labels_s6d <- c(
  "Residual-E3 granulocyte panel",
  "Primary azurophilic",
  "Tertiary gelatinase",
  "Secondary specific"
)

display_match_s6d <- match(
  display_order_s6d,
  correlation_summary_s6d$panel
)

stopifnot(!anyNA(display_match_s6d))

plot_data_s6d <-
  correlation_summary_s6d[
    display_match_s6d,
    ,
    drop = FALSE
  ]

plot_data_s6d$display_label <-
  display_labels_s6d

plot_data_s6d$row_y <- c(
  4,
  3,
  2,
  1
)

plot_data_s6d$bootstrap_y <-
  plot_data_s6d$row_y + 0.10

plot_data_s6d$lodo_y <-
  plot_data_s6d$row_y - 0.10

plot_data_s6d$FDR_status <- ifelse(
  plot_data_s6d$FDR < 0.05,
  "FDR < 0.05",
  "FDR >= 0.05"
)

format_fdr_s6d <- function(x) {
  
  if (x < 0.01) {
    return(sprintf("%.5f", x))
  }
  
  if (x < 0.10) {
    return(sprintf("%.4f", x))
  }
  
  sprintf("%.3f", x)
}

plot_data_s6d$FDR_label <- vapply(
  plot_data_s6d$FDR,
  format_fdr_s6d,
  character(1)
)

# ------------------------------------------------------------
# 11. Colours
# ------------------------------------------------------------

fdr_colours_s6d <- c(
  "FDR < 0.05" = "#E45756",
  "FDR >= 0.05" = "#8A8A8A"
)

# ------------------------------------------------------------
# 12. Construct the forest plot
# ------------------------------------------------------------

plot_s6d <- ggplot2::ggplot() +
  
  ggplot2::geom_hline(
    yintercept = plot_data_s6d$row_y,
    colour = "#EEEEEE",
    linewidth = 0.35
  ) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "#555555",
    linewidth = 0.45
  ) +
  
  # Leave-one-donor-out range
  ggplot2::geom_segment(
    data = plot_data_s6d,
    ggplot2::aes(
      x = rho_min_LODO,
      xend = rho_max_LODO,
      y = lodo_y,
      yend = lodo_y,
      linewidth =
        "Leave-one-donor-out range"
    ),
    colour = "#3F3F3F",
    lineend = "butt"
  ) +
  
  ggplot2::geom_segment(
    data = plot_data_s6d,
    ggplot2::aes(
      x = rho_min_LODO,
      xend = rho_min_LODO,
      y = lodo_y - 0.045,
      yend = lodo_y + 0.045
    ),
    colour = "#3F3F3F",
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_segment(
    data = plot_data_s6d,
    ggplot2::aes(
      x = rho_max_LODO,
      xend = rho_max_LODO,
      y = lodo_y - 0.045,
      yend = lodo_y + 0.045
    ),
    colour = "#3F3F3F",
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  
  # Bootstrap 95% confidence interval
  ggplot2::geom_segment(
    data = plot_data_s6d,
    ggplot2::aes(
      x = bootstrap_CI_low,
      xend = bootstrap_CI_high,
      y = bootstrap_y,
      yend = bootstrap_y,
      linewidth = "Bootstrap 95% CI",
      colour = FDR_status
    ),
    lineend = "round"
  ) +
  
  ggplot2::geom_segment(
    data = plot_data_s6d,
    ggplot2::aes(
      x = bootstrap_CI_low,
      xend = bootstrap_CI_low,
      y = bootstrap_y - 0.055,
      yend = bootstrap_y + 0.055,
      colour = FDR_status
    ),
    linewidth = 0.70,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_segment(
    data = plot_data_s6d,
    ggplot2::aes(
      x = bootstrap_CI_high,
      xend = bootstrap_CI_high,
      y = bootstrap_y - 0.055,
      yend = bootstrap_y + 0.055,
      colour = FDR_status
    ),
    linewidth = 0.70,
    show.legend = FALSE
  ) +
  
  # Full-data Spearman rho
  ggplot2::geom_point(
    data = plot_data_s6d,
    ggplot2::aes(
      x = rho_full,
      y = bootstrap_y,
      fill = FDR_status
    ),
    shape = 21,
    size = 3.6,
    stroke = 0.65,
    colour = "black"
  ) +
  
  # Right-side FDR column
  ggplot2::geom_text(
    data = plot_data_s6d,
    ggplot2::aes(
      x = 1.16,
      y = row_y,
      label = FDR_label
    ),
    hjust = 1,
    size = 2.75,
    family = base_family_s6d,
    colour = "#222222"
  ) +
  
  ggplot2::annotate(
    geom = "text",
    x = 1.16,
    y = 4.48,
    label = "FDR",
    hjust = 1,
    vjust = 1,
    size = 2.75,
    fontface = "bold",
    family = base_family_s6d,
    colour = "#222222"
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = c(
      -0.5,
      -0.25,
      0,
      0.25,
      0.5,
      0.75,
      1.0
    ),
    labels = function(x) {
      sprintf("%.2f", x)
    },
    expand = ggplot2::expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  ggplot2::scale_fill_manual(
    name = "FDR status",
    values = fdr_colours_s6d,
    breaks = c(
      "FDR < 0.05",
      "FDR >= 0.05"
    )
  ) +
  
  ggplot2::scale_colour_manual(
    values = fdr_colours_s6d,
    guide = "none"
  ) +
  
  ggplot2::scale_linewidth_manual(
    name = "Interval",
    values = c(
      "Bootstrap 95% CI" = 1.45,
      "Leave-one-donor-out range" = 0.50
    ),
    breaks = c(
      "Bootstrap 95% CI",
      "Leave-one-donor-out range"
    )
  ) +
  
  ggplot2::coord_cartesian(
    xlim = c(-0.50, 1.20),
    ylim = c(0.55, 4.50),
    clip = "off"
  ) +
  
  ggplot2::labs(
    title = paste0(
      "Residual-E3 and primary-azurophilic signals\n",
      "track neutrophil fraction across donors"
    ),
    subtitle = paste0(
      "Strict non-granulocyte compartment; ",
      "n = 21 donors per panel"
    ),
    x = expression(
      "Spearman " * rho *
        " with neutrophil fraction"
    ),
    y = NULL
  ) +
  
  ggplot2::guides(
    linewidth = ggplot2::guide_legend(
      order = 1,
      override.aes = list(
        colour = "#333333"
      )
    ),
    fill = ggplot2::guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        size = 3.2,
        colour = "black"
      )
    )
  ) +
  
  ggplot2::theme_classic(
    base_size = 8,
    base_family = base_family_s6d
  ) +
  
  ggplot2::theme(
    pplot.title = ggplot2::element_text(
      size = 9.2,
      face = "bold",
      colour = "black",
      lineheight = 0.98,
      margin = ggplot2::margin(
        b = 3
      )
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.3,
      colour = "#444444",
      margin = ggplot2::margin(
        b = 9
      )
    ),
    axis.title.x = ggplot2::element_text(
      size = 7.8,
      margin = ggplot2::margin(
        t = 7
      )
    ),
    axis.text.x = ggplot2::element_text(
      size = 7.0,
      colour = "#222222"
    ),
    axis.text.y = ggplot2::element_text(
      size = 7.4,
      colour = "#222222",
      margin = ggplot2::margin(
        r = 6
      )
    ),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    panel.grid.major.x =
      ggplot2::element_line(
        colour = "#E6E6E6",
        linewidth = 0.30
      ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_text(
      size = 7.0,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(
      size = 6.7
    ),
    legend.key.width = grid::unit(
      8,
      "mm"
    ),
    legend.spacing.x = grid::unit(
      2,
      "mm"
    ),
    plot.margin = ggplot2::margin(
      t = 8,
      r = 10,
      b = 5,
      l = 5
    )
  )

print(plot_s6d)

# ------------------------------------------------------------
# 13. Export source data
# ------------------------------------------------------------

source_data_s6d <- plot_data_s6d[
  ,
  c(
    "panel",
    "display_label",
    "n_donors",
    "rho_full",
    "P_value",
    "FDR",
    "bootstrap_CI_low",
    "bootstrap_CI_high",
    "rho_min_LODO",
    "rho_max_LODO",
    "FDR_status"
  ),
  drop = FALSE
]

source_data_file_s6d <- file.path(
  output_dir_s6d,
  paste0(
    "Supplementary_Figure_S6D_",
    "source_data.csv"
  )
)

write.csv(
  source_data_s6d,
  source_data_file_s6d,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 14. Publication-format export
# ------------------------------------------------------------

width_mm_s6d <- 180
height_mm_s6d <- 105
width_in_s6d <- width_mm_s6d / 25.4
height_in_s6d <- height_mm_s6d / 25.4

output_stem_s6d <- file.path(
  output_dir_s6d,
  paste0(
    "Supplementary_Figure_S6D_",
    "granule_correlation_forest"
  )
)

pdf_file_s6d <- paste0(
  output_stem_s6d,
  ".pdf"
)

svg_file_s6d <- paste0(
  output_stem_s6d,
  ".svg"
)

tiff_file_s6d <- paste0(
  output_stem_s6d,
  ".tiff"
)

png_file_s6d <- paste0(
  output_stem_s6d,
  ".png"
)

save_pdf_s6d <- function(
    filename,
    plot_object
) {
  
  grDevices::pdf(
    file = filename,
    width = width_in_s6d,
    height = height_in_s6d,
    family = "Helvetica",
    useDingbats = FALSE,
    onefile = TRUE
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  print(plot_object)
}

save_svg_s6d <- function(
    filename,
    plot_object
) {
  
  svglite::svglite(
    file = filename,
    width = width_in_s6d,
    height = height_in_s6d,
    bg = "white"
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  print(plot_object)
}

save_tiff_s6d <- function(
    filename,
    plot_object
) {
  
  ragg::agg_tiff(
    filename = filename,
    width = width_in_s6d,
    height = height_in_s6d,
    units = "in",
    res = 600,
    background = "white"
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  print(plot_object)
}

save_png_s6d <- function(
    filename,
    plot_object
) {
  
  ragg::agg_png(
    filename = filename,
    width = width_in_s6d,
    height = height_in_s6d,
    units = "in",
    res = 300,
    background = "white"
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  print(plot_object)
}

save_pdf_s6d(
  pdf_file_s6d,
  plot_s6d
)

save_svg_s6d(
  svg_file_s6d,
  plot_s6d
)

save_tiff_s6d(
  tiff_file_s6d,
  plot_s6d
)

save_png_s6d(
  png_file_s6d,
  plot_s6d
)

# ------------------------------------------------------------
# 15. Output audit
# ------------------------------------------------------------

expected_output_files_s6d <- c(
  pdf_file_s6d,
  svg_file_s6d,
  tiff_file_s6d,
  png_file_s6d,
  source_data_file_s6d
)

output_audit_s6d <- data.frame(
  file = expected_output_files_s6d,
  exists = file.exists(
    expected_output_files_s6d
  ),
  size_bytes = file.info(
    expected_output_files_s6d
  )$size,
  stringsAsFactors = FALSE
)

print(
  output_audit_s6d,
  row.names = FALSE
)

stopifnot(
  all(output_audit_s6d$exists),
  all(
    output_audit_s6d$size_bytes > 0
  )
)

cat(
  "\nS6D completed successfully.\n"
)

cat(
  "Output directory:\n",
  normalizePath(
    output_dir_s6d,
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "No panel letter was embedded in the figure.\n"
)

cat(
  "====================================================\n"
)


# ============================================================
# Supplementary Figure S7A
# Whole-sample condition contrasts
#
# No embedded panel letter
# Standalone from frozen CSV source data
# ============================================================

cat(
  "\n========== S7A: WHOLE-SAMPLE CONDITION CONTRASTS ==========\n"
)

# ------------------------------------------------------------
# 1. Package audit
# ------------------------------------------------------------

required_packages_s7a <- c(
  "ggplot2",
  "patchwork",
  "svglite",
  "ragg"
)

missing_packages_s7a <- required_packages_s7a[
  !vapply(
    required_packages_s7a,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_s7a) > 0L) {
  stop(
    paste0(
      "Missing packages required for S7A: ",
      paste(missing_packages_s7a, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# 2. Project, input and output paths
# ------------------------------------------------------------

sc_root_s7a <- if (
  exists("SC_ROOT_S6A", inherits = TRUE)
) {
  SC_ROOT_S6A
} else {
  "/home/sunshine/predicate/singlecell"
}

output_root_s7a <- if (
  exists("OUTPUT_ROOT_S6A", inherits = TRUE)
) {
  OUTPUT_ROOT_S6A
} else {
  file.path(
    sc_root_s7a,
    "outputs",
    "GSE216009_32gene_analysis",
    "manuscript_revision",
    "Supplementary_Figures_S6_S7"
  )
}

base_family_s7a <- if (
  exists("BASE_FAMILY_S6A", inherits = TRUE)
) {
  BASE_FAMILY_S6A
} else {
  "Helvetica"
}

input_dir_s7a <- file.path(
  sc_root_s7a,
  "outputs",
  "GSE216009_32gene_analysis",
  "05e_pseudobulk_sanity_state_direction",
  "20260805_internal_pseudobulk_sanity_and_state_direction_3",
  "tables"
)

summary_file_s7a <- file.path(
  input_dir_s7a,
  "p3a_module_summary.csv"
)

sample_file_s7a <- file.path(
  input_dir_s7a,
  "p3a_sample_module_scores.csv"
)

output_dir_s7a <- file.path(
  output_root_s7a,
  "S7A_whole_sample_condition_contrasts"
)

dir.create(
  output_dir_s7a,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files_s7a <- c(
  summary_file_s7a,
  sample_file_s7a
)

file_audit_s7a <- data.frame(
  file = required_files_s7a,
  exists = file.exists(required_files_s7a),
  stringsAsFactors = FALSE
)

print(
  file_audit_s7a,
  row.names = FALSE
)

if (!all(file_audit_s7a$exists)) {
  stop(
    paste0(
      "Missing S7A input file(s):\n",
      paste(
        file_audit_s7a$file[
          !file_audit_s7a$exists
        ],
        collapse = "\n"
      )
    )
  )
}

# ------------------------------------------------------------
# 3. Read frozen source data
# ------------------------------------------------------------

module_summary_s7a <- read.csv(
  summary_file_s7a,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_scores_s7a <- read.csv(
  sample_file_s7a,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_summary_columns_s7a <- c(
  "score",
  "acute_median",
  "conv_median",
  "hv_median",
  "cs_median",
  "paired_delta_median",
  "paired_delta_n",
  "acute_minus_hv",
  "acute_minus_cs"
)

required_sample_columns_s7a <- c(
  "score",
  "sample_id",
  "cohort",
  "donor",
  "value"
)

missing_summary_columns_s7a <- setdiff(
  required_summary_columns_s7a,
  colnames(module_summary_s7a)
)

missing_sample_columns_s7a <- setdiff(
  required_sample_columns_s7a,
  colnames(sample_scores_s7a)
)

if (length(missing_summary_columns_s7a) > 0L) {
  stop(
    paste0(
      "Missing summary columns: ",
      paste(
        missing_summary_columns_s7a,
        collapse = ", "
      )
    )
  )
}

if (length(missing_sample_columns_s7a) > 0L) {
  stop(
    paste0(
      "Missing sample-level columns: ",
      paste(
        missing_sample_columns_s7a,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# 4. Freeze readout and condition definitions
# ------------------------------------------------------------

score_order_s7a <- c(
  "Up29_composite",
  "Down3_composite",
  "total32_signed",
  "M1_cell_cycle",
  "M2_full",
  "primary_core"
)

condition_order_s7a <- c(
  "Acute",
  "Conv",
  "HV",
  "CS"
)

condition_axis_labels_s7a <- c(
  Acute = "Acute (n=26)",
  Conv = "Convalescent (n=9)",
  HV = "Healthy volunteers (n=6)",
  CS = "Surgery controls (n=7)"
)

readout_labels_s7a <- c(
  Up29_composite =
    "29-gene NS-higher composite",
  Down3_composite =
    "3-gene NS-lower composite",
  total32_signed =
    "32-gene signed score",
  M1_cell_cycle =
    "M1 cell cycle/proliferation",
  M2_full =
    "M2 neutrophil degranulation",
  primary_core =
    "Primary azurophilic/core"
)

stopifnot(
  nrow(module_summary_s7a) == 6L,
  nrow(sample_scores_s7a) == 288L,
  identical(
    as.character(module_summary_s7a$score),
    score_order_s7a
  ),
  setequal(
    unique(as.character(sample_scores_s7a$score)),
    score_order_s7a
  ),
  setequal(
    unique(as.character(sample_scores_s7a$cohort)),
    condition_order_s7a
  )
)

sample_scores_s7a$score <- as.character(
  sample_scores_s7a$score
)

sample_scores_s7a$sample_id <- as.character(
  sample_scores_s7a$sample_id
)

sample_scores_s7a$cohort <- as.character(
  sample_scores_s7a$cohort
)

sample_scores_s7a$donor <- as.character(
  sample_scores_s7a$donor
)

sample_scores_s7a$value <- suppressWarnings(
  as.numeric(sample_scores_s7a$value)
)

numeric_summary_columns_s7a <- setdiff(
  required_summary_columns_s7a,
  "score"
)

for (current_column_s7a in numeric_summary_columns_s7a) {
  
  module_summary_s7a[[current_column_s7a]] <-
    suppressWarnings(
      as.numeric(
        module_summary_s7a[[current_column_s7a]]
      )
    )
}

stopifnot(
  all(is.finite(sample_scores_s7a$value)),
  all(
    vapply(
      module_summary_s7a[
        numeric_summary_columns_s7a
      ],
      function(x) all(is.finite(x)),
      logical(1)
    )
  ),
  all(module_summary_s7a$paired_delta_n == 9L),
  !anyDuplicated(
    sample_scores_s7a[
      c("score", "sample_id")
    ]
  )
)

# ------------------------------------------------------------
# 5. Verify sample counts
# ------------------------------------------------------------

observed_counts_s7a <- table(
  factor(
    sample_scores_s7a$score,
    levels = score_order_s7a
  ),
  factor(
    sample_scores_s7a$cohort,
    levels = condition_order_s7a
  )
)

expected_counts_s7a <- matrix(
  rep(
    c(26L, 9L, 6L, 7L),
    times = length(score_order_s7a)
  ),
  nrow = length(score_order_s7a),
  byrow = TRUE,
  dimnames = list(
    score = score_order_s7a,
    cohort = condition_order_s7a
  )
)

cat("\nObserved score-by-condition counts:\n")
print(observed_counts_s7a)

stopifnot(
  identical(
    as.integer(observed_counts_s7a),
    as.integer(expected_counts_s7a)
  )
)

# ------------------------------------------------------------
# 6. Identify the nine paired donors
# ------------------------------------------------------------

paired_donors_by_score_s7a <- lapply(
  score_order_s7a,
  function(current_score_s7a) {
    
    current_data_s7a <- sample_scores_s7a[
      sample_scores_s7a$score ==
        current_score_s7a,
      ,
      drop = FALSE
    ]
    
    acute_donors_s7a <- unique(
      current_data_s7a$donor[
        current_data_s7a$cohort == "Acute"
      ]
    )
    
    conv_donors_s7a <- unique(
      current_data_s7a$donor[
        current_data_s7a$cohort == "Conv"
      ]
    )
    
    sort(
      intersect(
        acute_donors_s7a,
        conv_donors_s7a
      )
    )
  }
)

names(paired_donors_by_score_s7a) <-
  score_order_s7a

stopifnot(
  all(
    vapply(
      paired_donors_by_score_s7a,
      length,
      integer(1)
    ) == 9L
  )
)

paired_donors_s7a <-
  paired_donors_by_score_s7a[[1L]]

same_paired_set_s7a <- vapply(
  paired_donors_by_score_s7a,
  function(x) {
    identical(
      x,
      paired_donors_s7a
    )
  },
  logical(1)
)

stopifnot(
  all(same_paired_set_s7a)
)

# ------------------------------------------------------------
# 7. Recalculate paired medians from sample-level data
# ------------------------------------------------------------

recalculated_delta_s7a <- vapply(
  score_order_s7a,
  function(current_score_s7a) {
    
    current_data_s7a <- sample_scores_s7a[
      sample_scores_s7a$score ==
        current_score_s7a &
        sample_scores_s7a$cohort %in%
        c("Acute", "Conv") &
        sample_scores_s7a$donor %in%
        paired_donors_s7a,
      ,
      drop = FALSE
    ]
    
    acute_data_s7a <- current_data_s7a[
      current_data_s7a$cohort == "Acute",
      c("donor", "value"),
      drop = FALSE
    ]
    
    conv_data_s7a <- current_data_s7a[
      current_data_s7a$cohort == "Conv",
      c("donor", "value"),
      drop = FALSE
    ]
    
    stopifnot(
      !anyDuplicated(acute_data_s7a$donor),
      !anyDuplicated(conv_data_s7a$donor)
    )
    
    acute_data_s7a <- acute_data_s7a[
      match(
        paired_donors_s7a,
        acute_data_s7a$donor
      ),
      ,
      drop = FALSE
    ]
    
    conv_data_s7a <- conv_data_s7a[
      match(
        paired_donors_s7a,
        conv_data_s7a$donor
      ),
      ,
      drop = FALSE
    ]
    
    stopifnot(
      identical(
        acute_data_s7a$donor,
        conv_data_s7a$donor
      )
    )
    
    stats::median(
      acute_data_s7a$value -
        conv_data_s7a$value
    )
  },
  numeric(1)
)

summary_delta_s7a <-
  module_summary_s7a$paired_delta_median[
    match(
      score_order_s7a,
      module_summary_s7a$score
    )
  ]

stopifnot(
  max(
    abs(
      recalculated_delta_s7a -
        summary_delta_s7a
    )
  ) < 1e-8
)

module_summary_s7a$paired_delta_recalculated <-
  recalculated_delta_s7a[
    match(
      module_summary_s7a$score,
      names(recalculated_delta_s7a)
    )
  ]

m1_delta_s7a <- module_summary_s7a$
  paired_delta_median[
    module_summary_s7a$score ==
      "M1_cell_cycle"
  ]

m2_delta_s7a <- module_summary_s7a$
  paired_delta_median[
    module_summary_s7a$score ==
      "M2_full"
  ]

stopifnot(
  abs(
    m1_delta_s7a -
      0.3515559137323
  ) < 1e-8,
  abs(
    m2_delta_s7a -
      1.0906944112015
  ) < 1e-8,
  module_summary_s7a$score[
    which.max(
      module_summary_s7a$paired_delta_median
    )
  ] == "M2_full"
)

# ------------------------------------------------------------
# 8. Plotting palettes
# ------------------------------------------------------------

m1_condition_palette_s7a <- c(
  Acute = "#3E82C4",
  Conv = "#86D6C9",
  HV = "#B8BDC3",
  CS = "#D8DBDE"
)

m2_condition_palette_s7a <- c(
  Acute = "#EF3E36",
  Conv = "#FFA37F",
  HV = "#B8BDC3",
  CS = "#D8DBDE"
)

delta_palette_s7a <- c(
  M1 = "#3E82C4",
  M2 = "#EF3E36",
  Other = "#8C8C8C"
)

# ------------------------------------------------------------
# 9. Prepare deterministic paired-point offsets
# ------------------------------------------------------------

paired_offset_s7a <- seq(
  from = -0.12,
  to = 0.12,
  length.out = length(paired_donors_s7a)
)

names(paired_offset_s7a) <-
  paired_donors_s7a

sample_scores_s7a$condition_position <- match(
  sample_scores_s7a$cohort,
  condition_order_s7a
)

sample_scores_s7a$paired_ac_conv <-
  sample_scores_s7a$donor %in%
  paired_donors_s7a &
  sample_scores_s7a$cohort %in%
  c("Acute", "Conv")

sample_scores_s7a$paired_position <-
  sample_scores_s7a$condition_position

paired_rows_s7a <-
  sample_scores_s7a$paired_ac_conv

sample_scores_s7a$paired_position[
  paired_rows_s7a
] <- sample_scores_s7a$condition_position[
  paired_rows_s7a
] +
  unname(
    paired_offset_s7a[
      sample_scores_s7a$donor[
        paired_rows_s7a
      ]
    ]
  )

sample_scores_s7a$readout_display <-
  unname(
    readout_labels_s7a[
      sample_scores_s7a$score
    ]
  )

sample_scores_s7a$condition_display <-
  unname(
    condition_axis_labels_s7a[
      sample_scores_s7a$cohort
    ]
  )

# ------------------------------------------------------------
# 10. Function for M1 and M2 condition distributions
# ------------------------------------------------------------

make_condition_plot_s7a <- function(
    requested_score,
    plot_title,
    condition_palette,
    show_y_title = TRUE,
    jitter_seed = 1L
) {
  
  current_data_s7a <- sample_scores_s7a[
    sample_scores_s7a$score ==
      requested_score,
    ,
    drop = FALSE
  ]
  
  paired_data_s7a <- current_data_s7a[
    current_data_s7a$paired_ac_conv,
    ,
    drop = FALSE
  ]
  
  unpaired_data_s7a <- current_data_s7a[
    !current_data_s7a$paired_ac_conv,
    ,
    drop = FALSE
  ]
  
  current_y_title_s7a <- if (
    show_y_title
  ) {
    expression(
      atop(
        "Whole-sample composite expression",
        log[2] * "(CPM + 1)"
      )
    )
  } else {
    NULL
  }
  
  ggplot(
    current_data_s7a,
    aes(
      x = condition_position,
      y = value
    )
  ) +
    geom_boxplot(
      aes(
        group = cohort,
        fill = cohort
      ),
      width = 0.56,
      linewidth = 0.38,
      colour = "#333333",
      alpha = 0.28,
      outlier.shape = NA
    ) +
    geom_line(
      data = paired_data_s7a,
      aes(
        x = paired_position,
        y = value,
        group = donor
      ),
      inherit.aes = FALSE,
      linewidth = 0.34,
      colour = "#9A9A9A",
      alpha = 0.65
    ) +
    geom_point(
      data = unpaired_data_s7a,
      aes(
        x = condition_position,
        y = value,
        fill = cohort
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.65,
      stroke = 0.27,
      colour = "#333333",
      position = position_jitter(
        width = 0.09,
        height = 0,
        seed = jitter_seed
      )
    ) +
    geom_point(
      data = paired_data_s7a,
      aes(
        x = paired_position,
        y = value,
        fill = cohort
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.75,
      stroke = 0.30,
      colour = "#222222"
    ) +
    scale_fill_manual(
      values = condition_palette,
      breaks = condition_order_s7a,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(condition_order_s7a),
      labels = condition_axis_labels_s7a[
        condition_order_s7a
      ],
      limits = c(0.55, 4.45),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.06, 0.12)
      )
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = current_y_title_s7a
    ) +
    theme_classic(
      base_size = 8.0,
      base_family = base_family_s7a
    ) +
    theme(
      plot.title = element_text(
        size = 9.2,
        face = "bold",
        colour = "#111111",
        margin = margin(b = 6)
      ),
      axis.title.y = element_text(
        size = 8.4,
        margin = margin(r = 6)
      ),
      axis.text.y = element_text(
        size = 7.5,
        colour = "#222222"
      ),
      axis.text.x = element_text(
        size = 6.7,
        angle = 60,
        hjust = 1,
        vjust = 1,
        colour = "#222222",
        lineheight = 1,
        margin = margin(t = 6)
      ),
      axis.ticks = element_line(
        linewidth = 0.35,
        colour = "#222222"
      ),
      axis.line = element_line(
        linewidth = 0.40,
        colour = "#222222"
      ),
      legend.position = "none",
      plot.margin = margin(
        t = 4,
        r = 5,
        b = 3,
        l = 4
      )
    )
}

plot_m1_s7a <- make_condition_plot_s7a(
  requested_score = "M1_cell_cycle",
  plot_title = "M1 cell-cycle/\nproliferation",
  condition_palette = m1_condition_palette_s7a,
  show_y_title = TRUE,
  jitter_seed = 20260821L
)

plot_m2_s7a <- make_condition_plot_s7a(
  requested_score = "M2_full",
  plot_title = "M2 neutrophil-\ndegranulation",
  condition_palette = m2_condition_palette_s7a,
  show_y_title = FALSE,
  jitter_seed = 20260822L
)

# ------------------------------------------------------------
# 11. Six-readout paired-change plot
# ------------------------------------------------------------

readout_top_to_bottom_s7a <- unname(
  readout_labels_s7a[
    score_order_s7a
  ]
)

module_summary_s7a$readout_display <-
  unname(
    readout_labels_s7a[
      module_summary_s7a$score
    ]
  )

module_summary_s7a$readout_factor <- factor(
  module_summary_s7a$readout_display,
  levels = rev(
    readout_top_to_bottom_s7a
  )
)

module_summary_s7a$display_group <- "Other"

module_summary_s7a$display_group[
  module_summary_s7a$score ==
    "M1_cell_cycle"
] <- "M1"

module_summary_s7a$display_group[
  module_summary_s7a$score ==
    "M2_full"
] <- "M2"

highlight_delta_s7a <- module_summary_s7a[
  module_summary_s7a$score %in%
    c(
      "M1_cell_cycle",
      "M2_full"
    ),
  ,
  drop = FALSE
]

highlight_delta_s7a$delta_label <- ifelse(
  highlight_delta_s7a$score ==
    "M2_full",
  sprintf(
    "+%.2f",
    highlight_delta_s7a$paired_delta_median
  ),
  sprintf(
    "+%.3f",
    highlight_delta_s7a$paired_delta_median
  )
)

plot_delta_s7a <- ggplot(
  module_summary_s7a,
  aes(
    y = readout_factor,
    x = paired_delta_median,
    colour = display_group
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.42,
    linetype = "dashed",
    colour = "#555555"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = paired_delta_median,
      yend = readout_factor
    ),
    linewidth = 1.05,
    alpha = 0.42,
    lineend = "round"
  ) +
  geom_point(
    size = 2.75,
    stroke = 0.25
  ) +
  geom_text(
    data = highlight_delta_s7a,
    aes(
      label = delta_label
    ),
    hjust = 0,
    nudge_x = 0.055,
    size = 3.05,
    family = base_family_s7a,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = delta_palette_s7a,
    breaks = c(
      "M1",
      "M2",
      "Other"
    )
  ) +
  scale_x_continuous(
    breaks = seq(
      from = -0.5,
      to = 1.5,
      by = 0.25
    ),
    limits = c(-0.55, 1.55),
    expand = c(0, 0)
  ) +
  labs(
    title = "Paired change across six readouts",
    subtitle = "Median within-donor difference; n = 9 pairs",
    x = "Median paired change (acute - convalescent)",
    y = NULL
  ) +
  theme_classic(
    base_size = 8.0,
    base_family = base_family_s7a
  ) +
  theme(
    plot.title = element_text(
      size = 9.2,
      face = "bold",
      colour = "#111111",
      margin = margin(b = 3)
    ),
    plot.subtitle = element_text(
      size = 7.2,
      colour = "#555555",
      margin = margin(b = 6)
    ),
    axis.title.x = element_text(
      size = 8.2,
      margin = margin(t = 6)
    ),
    axis.text.x = element_text(
      size = 7.2,
      colour = "#222222"
    ),
    axis.text.y = element_text(
      size = 7.2,
      colour = "#222222",
      lineheight = 0.95
    ),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(
      linewidth = 0.40
    ),
    axis.line.y = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.45
    ),
    legend.position = "none",
    plot.margin = margin(
      t = 4,
      r = 12,
      b = 3,
      l = 8
    )
  )

# ------------------------------------------------------------
# 12. Assemble figure
#
# wrap_elements() prevents patchwork from flattening the
# two left-hand plots into the outer two-column layout.
# ------------------------------------------------------------

left_condition_block_s7a <- patchwork::wrap_plots(
  plot_m1_s7a,
  plot_m2_s7a,
  nrow = 1,
  ncol = 2,
  widths = c(1, 1)
)

left_condition_container_s7a <-
  patchwork::wrap_elements(
    full = left_condition_block_s7a,
    clip = FALSE
  )

plot_s7a <- patchwork::wrap_plots(
  left_condition_container_s7a,
  plot_delta_s7a,
  nrow = 1,
  ncol = 2,
  widths = c(
    2.05,
    1.45
  )
) +
  patchwork::plot_annotation(
    title = paste0(
      "Whole-sample condition contrasts identified ",
      "the largest paired acute shift in M2"
    ),
    subtitle = paste0(
      "Points represent samples; lines connect the nine ",
      "acute-convalescent donor pairs."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = base_family_s7a,
        size = 11.0,
        face = "bold",
        colour = "#111111",
        hjust = 0,
        margin = ggplot2::margin(
          b = 4
        )
      ),
      plot.subtitle = ggplot2::element_text(
        family = base_family_s7a,
        size = 7.8,
        colour = "#555555",
        hjust = 0,
        margin = ggplot2::margin(
          b = 7
        )
      ),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 9,
        b = 9,
        l = 9
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      )
    )
  )

print(plot_s7a)

# ------------------------------------------------------------
# 13. Export source data
# ------------------------------------------------------------

summary_source_file_s7a <- file.path(
  output_dir_s7a,
  "Supplementary_Figure_S7A_summary_source_data.csv"
)

sample_source_file_s7a <- file.path(
  output_dir_s7a,
  "Supplementary_Figure_S7A_sample_level_source_data.csv"
)

write.csv(
  module_summary_s7a,
  summary_source_file_s7a,
  row.names = FALSE,
  na = ""
)

write.csv(
  sample_scores_s7a,
  sample_source_file_s7a,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 14. Export PDF, SVG, TIFF and PNG
# ------------------------------------------------------------

file_stem_s7a <- file.path(
  output_dir_s7a,
  "Supplementary_Figure_S7A_whole_sample_condition_contrasts"
)

width_mm_s7a <- 183
height_mm_s7a <- 118
dpi_s7a <- 600

width_in_s7a <- width_mm_s7a / 25.4
height_in_s7a <- height_mm_s7a / 25.4

# PDF
grDevices::pdf(
  file = paste0(file_stem_s7a, ".pdf"),
  width = width_in_s7a,
  height = height_in_s7a,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE
)

print(plot_s7a)
grDevices::dev.off()

# SVG
svglite::svglite(
  file = paste0(file_stem_s7a, ".svg"),
  width = width_in_s7a,
  height = height_in_s7a
)

print(plot_s7a)
grDevices::dev.off()

# TIFF
ragg::agg_tiff(
  filename = paste0(file_stem_s7a, ".tiff"),
  width = width_in_s7a,
  height = height_in_s7a,
  units = "in",
  res = dpi_s7a,
  background = "white",
  compression = "lzw"
)

print(plot_s7a)
grDevices::dev.off()

# PNG preview
ragg::agg_png(
  filename = paste0(file_stem_s7a, ".png"),
  width = width_in_s7a,
  height = height_in_s7a,
  units = "in",
  res = 300,
  background = "white"
)

print(plot_s7a)
grDevices::dev.off()

# ------------------------------------------------------------
# 15. Final output audit
# ------------------------------------------------------------

figure_files_s7a <- c(
  paste0(file_stem_s7a, ".pdf"),
  paste0(file_stem_s7a, ".svg"),
  paste0(file_stem_s7a, ".tiff"),
  paste0(file_stem_s7a, ".png")
)

output_audit_file_s7a <- file.path(
  output_dir_s7a,
  "Supplementary_Figure_S7A_output_audit.csv"
)

output_audit_s7a <- data.frame(
  output_type = c(
    "PDF figure",
    "SVG figure",
    "TIFF figure",
    "PNG preview",
    "Summary source data",
    "Sample-level source data"
  ),
  file = c(
    figure_files_s7a,
    summary_source_file_s7a,
    sample_source_file_s7a
  ),
  stringsAsFactors = FALSE
)

output_audit_s7a$exists <- file.exists(
  output_audit_s7a$file
)

output_audit_s7a$size_bytes <- file.info(
  output_audit_s7a$file
)$size

write.csv(
  output_audit_s7a,
  output_audit_file_s7a,
  row.names = FALSE
)

expected_outputs_s7a <- c(
  output_audit_s7a$file,
  output_audit_file_s7a
)

stopifnot(
  all(file.exists(expected_outputs_s7a)),
  all(file.info(expected_outputs_s7a)$size > 0)
)

cat(
  "\n====================================================\n",
  "S7A COMPLETED SUCCESSFULLY\n",
  "====================================================\n",
  sep = ""
)

cat(
  "Output directory:\n",
  normalizePath(
    output_dir_s7a,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

print(
  data.frame(
    file = expected_outputs_s7a,
    size_bytes = file.info(
      expected_outputs_s7a
    )$size,
    row.names = NULL
  ),
  row.names = FALSE
)

cat(
  "\nVerified paired changes:\n",
  sprintf(
    "M1 acute - convalescent median = %.6f\n",
    m1_delta_s7a
  ),
  sprintf(
    "M2 acute - convalescent median = %.6f\n",
    m2_delta_s7a
  ),
  sep = ""
)

cat(
  "No panel letter was embedded in the figure.\n"
)



# ============================================================
# Supplementary Figure S7B
# Paired pseudobulk direction across 32 signature genes
#
# No embedded panel letter
# ============================================================

cat(
  "\n========== S7B: PAIRED GENE DIRECTION ==========\n"
)

# ------------------------------------------------------------
# 1. Package audit
# ------------------------------------------------------------

required_packages_s7b <- c(
  "ggplot2",
  "patchwork",
  "svglite",
  "ragg"
)

missing_packages_s7b <- required_packages_s7b[
  !vapply(
    required_packages_s7b,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_s7b) > 0L) {
  stop(
    paste0(
      "Missing packages required for S7B: ",
      paste(missing_packages_s7b, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# 2. Input and output paths
# ------------------------------------------------------------

sc_root_s7b <- "/home/sunshine/predicate/singlecell"

output_root_s7b <- file.path(
  sc_root_s7b,
  "outputs",
  "GSE216009_32gene_analysis",
  "manuscript_revision",
  "Supplementary_Figures_S6_S7"
)

base_family_s7b <- "Helvetica"

input_dir_s7b <- file.path(
  sc_root_s7b,
  "outputs",
  "GSE216009_32gene_analysis",
  "05e_pseudobulk_sanity_state_direction",
  "20260923_internal_pseudobulk_sanity_and_state_direction",
  "tables"
)

direction_file_s7b <- file.path(
  input_dir_s7b,
  "p3a_gene_direction_consistency.csv"
)

output_dir_s7b <- file.path(
  output_root_s7b,
  paste0(
    "S9B_paired_gene_direction_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

if (dir.exists(output_dir_s7b)) {
  stop("Output directory already exists; use a new output directory.")
}

dir.create(
  output_dir_s7b,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(direction_file_s7b)) {
  stop(
    paste0(
      "S7B source file was not found:\n",
      direction_file_s7b
    )
  )
}

# ------------------------------------------------------------
# 3. Read and validate frozen source table
# ------------------------------------------------------------

direction_raw_s7b <- read.csv(
  direction_file_s7b,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns_s7b <- c(
  "gene",
  "dataset_feature",
  "direction",
  "n_paired_donors",
  "median_signed_delta",
  "support_rate_signed_gt0",
  "direction_call",
  "direction_call_abs0p5"
)

missing_columns_s7b <- setdiff(
  required_columns_s7b,
  colnames(direction_raw_s7b)
)

if (length(missing_columns_s7b) > 0L) {
  stop(
    paste0(
      "Missing S7B columns: ",
      paste(
        missing_columns_s7b,
        collapse = ", "
      )
    )
  )
}

direction_raw_s7b$gene <- as.character(
  direction_raw_s7b$gene
)

direction_raw_s7b$dataset_feature <- as.character(
  direction_raw_s7b$dataset_feature
)

direction_raw_s7b$direction <- as.character(
  direction_raw_s7b$direction
)

direction_raw_s7b$direction_call <- as.character(
  direction_raw_s7b$direction_call
)

numeric_columns_s7b <- c(
  "n_paired_donors",
  "median_signed_delta",
  "support_rate_signed_gt0"
)

for (current_column_s7b in numeric_columns_s7b) {
  direction_raw_s7b[[current_column_s7b]] <-
    suppressWarnings(
      as.numeric(
        direction_raw_s7b[[current_column_s7b]]
      )
    )
}

stopifnot(
  nrow(direction_raw_s7b) == 32L,
  !anyDuplicated(direction_raw_s7b$gene),
  all(direction_raw_s7b$n_paired_donors == 9L),
  all(is.finite(direction_raw_s7b$median_signed_delta)),
  all(is.finite(direction_raw_s7b$support_rate_signed_gt0)),
  all(
    direction_raw_s7b$support_rate_signed_gt0 >= 0 &
      direction_raw_s7b$support_rate_signed_gt0 <= 1
  )
)

# ------------------------------------------------------------
# 4. Freeze the signature-module order
# ------------------------------------------------------------

m1_genes_s7b <- c(
  "ANLN",
  "ASPM",
  "BUB1",
  "CCNA2",
  "CKS2",
  "KIF11",
  "MKI67",
  "RRM2",
  "TOP2A",
  "TPX2",
  "TYMS"
)

m2_genes_s7b <- c(
  "CEACAM6",
  "CEACAM8",
  "CTSG",
  "DEFA4",
  "ELANE",
  "MPO",
  "MS4A3",
  "OLFM4",
  "RNASE3",
  "TCN1"
)

inflammatory_genes_s7b <- c(
  "CX3CR1",
  "IL1B",
  "SECTM1"
)

other_genes_s7b <- c(
  "ABCA13",
  "CA2",
  "CD24",
  "CHI3L1",
  "HIST1H2BM",
  "HIST1H3B",
  "RHAG",
  "SERPINB10"
)

signature_gene_order_s7b <- c(
  m1_genes_s7b,
  m2_genes_s7b,
  inflammatory_genes_s7b,
  other_genes_s7b
)

module_order_s7b <- c(
  "M1",
  "M2",
  "Inflammatory/down",
  "Other"
)

module_vector_s7b <- rep(
  module_order_s7b,
  times = c(
    length(m1_genes_s7b),
    length(m2_genes_s7b),
    length(inflammatory_genes_s7b),
    length(other_genes_s7b)
  )
)

module_map_s7b <- data.frame(
  gene = signature_gene_order_s7b,
  module = module_vector_s7b,
  stringsAsFactors = FALSE
)

stopifnot(
  length(signature_gene_order_s7b) == 32L,
  !anyDuplicated(signature_gene_order_s7b),
  setequal(
    direction_raw_s7b$gene,
    signature_gene_order_s7b
  )
)

plot_data_s7b <- direction_raw_s7b[
  match(
    signature_gene_order_s7b,
    direction_raw_s7b$gene
  ),
  ,
  drop = FALSE
]

plot_data_s7b$module <- module_map_s7b$module

stopifnot(
  identical(
    plot_data_s7b$gene,
    signature_gene_order_s7b
  )
)

# ------------------------------------------------------------
# 5. Normalize the frozen direction-call labels
# ------------------------------------------------------------

direction_call_text_s7b <- tolower(
  trimws(plot_data_s7b$direction_call)
)

direction_label_map_s7b <- c(
  aligned = "Aligned",
  neither_criterion_met = "Neither criterion met",
  opposite_direction = "Opposite direction"
)

plot_data_s7b$direction_call_display <- unname(
  direction_label_map_s7b[direction_call_text_s7b]
)

if (anyNA(plot_data_s7b$direction_call_display)) {
  stop(
    paste0(
      "Unexpected direction-call labels: ",
      paste(
        unique(
          plot_data_s7b$direction_call[
            is.na(plot_data_s7b$direction_call_display)
          ]
        ),
        collapse = ", "
      ),
      ". Check the revised 05e input."
    )
  )
}

direction_call_order_s7b <- c(
  "Aligned",
  "Neither criterion met",
  "Opposite direction"
)

observed_call_counts_s7b <- table(
  factor(
    plot_data_s7b$direction_call_display,
    levels = direction_call_order_s7b
  )
)

expected_call_counts_s7b <- c(
  Aligned = 14L,
  `Neither criterion met` = 13L,
  `Opposite direction` = 5L
)

cat("\nDirection-call counts:\n")
print(observed_call_counts_s7b)

stopifnot(
  identical(
    as.integer(observed_call_counts_s7b),
    as.integer(expected_call_counts_s7b)
  )
)

observed_opposite_genes_s7b <- sort(
  plot_data_s7b$gene[
    plot_data_s7b$direction_call_display ==
      "Opposite direction"
  ]
)

expected_opposite_genes_s7b <- sort(
  c("CTSG", "CHI3L1", "IL1B", "CKS2", "CA2")
)

stopifnot(
  identical(
    observed_opposite_genes_s7b,
    expected_opposite_genes_s7b
  )
)
# ------------------------------------------------------------
# 6. Plotting variables
# ------------------------------------------------------------

plot_data_s7b$support_percentage <-
  100 *
  plot_data_s7b$support_rate_signed_gt0

plot_data_s7b$plot_order <-
  seq_len(nrow(plot_data_s7b))

plot_data_s7b$gene_factor <- factor(
  plot_data_s7b$gene,
  levels = rev(
    signature_gene_order_s7b
  )
)

plot_data_s7b$module <- factor(
  plot_data_s7b$module,
  levels = module_order_s7b
)

plot_data_s7b$direction_call_display <- factor(
  plot_data_s7b$direction_call_display,
  levels = direction_call_order_s7b
)

group_sizes_s7b <- c(
  length(m1_genes_s7b),
  length(m2_genes_s7b),
  length(inflammatory_genes_s7b),
  length(other_genes_s7b)
)

separator_positions_s7b <-
  length(signature_gene_order_s7b) -
  cumsum(group_sizes_s7b)[
    -length(group_sizes_s7b)
  ] +
  0.5

stopifnot(
  identical(
    as.numeric(separator_positions_s7b),
    c(21.5, 11.5, 8.5)
  )
)

# ------------------------------------------------------------
# 7. Palettes
# ------------------------------------------------------------

module_palette_s7b <- c(
  M1 = "#3E82C4",
  M2 = "#EF3E36",
  `Inflammatory/down` = "#59A14F",
  Other = "#B07AA1"
)

direction_palette_s7b <- c(
  Aligned = "#176D82",
  `Neither criterion met` = "#B7BBC0",
  `Opposite direction` = "#A33B72"
)

direction_scale_s7b <- function() {
  scale_colour_manual(
    values = direction_palette_s7b,
    breaks = direction_call_order_s7b,
    drop = FALSE,
    name = "Direction call"
  )
}

# ------------------------------------------------------------
# 8. Determine a symmetric left-axis range
# ------------------------------------------------------------

maximum_absolute_delta_s7b <- max(
  abs(
    c(
      plot_data_s7b$median_signed_delta,
      -0.5,
      0.3
    )
  )
)

left_axis_limit_s7b <- ceiling(
  (
    maximum_absolute_delta_s7b +
      0.10
  ) /
    0.25
) *
  0.25

left_axis_limits_s7b <- c(
  -left_axis_limit_s7b,
  left_axis_limit_s7b
)

left_axis_breaks_s7b <- pretty(
  left_axis_limits_s7b,
  n = 5
)

left_axis_breaks_s7b <- left_axis_breaks_s7b[
  left_axis_breaks_s7b >=
    left_axis_limits_s7b[1] &
    left_axis_breaks_s7b <=
    left_axis_limits_s7b[2]
]

# ------------------------------------------------------------
# 9. Common theme
# ------------------------------------------------------------

common_theme_s7b <- theme_classic(
  base_size = 9.0,
  base_family = base_family_s7b
) +
  theme(
    plot.title = element_text(
      size = 10.2,
      face = "bold",
      colour = "#111111",
      lineheight = 0.95,
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 7.8,
      colour = "#555555",
      lineheight = 0.95,
      margin = margin(b = 7)
    ),
    axis.title = element_text(
      size = 8.8,
      colour = "#111111"
    ),
    axis.text = element_text(
      size = 8.0,
      colour = "#222222"
    ),
    axis.line = element_line(
      linewidth = 0.42,
      colour = "#222222"
    ),
    axis.ticks = element_line(
      linewidth = 0.36,
      colour = "#222222"
    ),
    plot.margin = margin(
      t = 4,
      r = 5,
      b = 4,
      l = 5
    )
  )

# ------------------------------------------------------------
# 10. Module-identity strip and gene labels
# ------------------------------------------------------------

plot_module_s7b <- ggplot(
  plot_data_s7b,
  aes(
    x = 1,
    y = gene_factor,
    fill = module
  )
) +
  geom_tile(
    width = 0.30,
    height = 0.91
  ) +
  geom_hline(
    yintercept = separator_positions_s7b,
    linewidth = 0.45,
    colour = "#555555"
  ) +
  scale_fill_manual(
    values = module_palette_s7b,
    breaks = module_order_s7b,
    labels = c(
      "P-program",
      "G-program",
      "NS-lower",
      "Other NS-higher"
    ),
    drop = FALSE,
    name = "Functional group"
  ) +
  scale_x_continuous(
    limits = c(0.80, 1.20),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(
      add = 0.48
    )
  ) +
  labs(
    title = "Gene",
    subtitle = " ",
    x = NULL,
    y = NULL
  ) +
  common_theme_s7b +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(
      size = 8.2,
      face = "italic",
      colour = "#222222",
      margin = margin(r = 5)
    ),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom"
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

# ------------------------------------------------------------
# 11. Direction-oriented paired-change plot
# ------------------------------------------------------------

plot_delta_gene_s7b <- ggplot(
  plot_data_s7b,
  aes(
    x = median_signed_delta,
    y = gene_factor,
    colour = direction_call_display
  )
) +
  geom_hline(
    yintercept = separator_positions_s7b,
    linewidth = 0.45,
    colour = "#555555"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.48,
    colour = "#333333"
  ) +
  geom_vline(
    xintercept = c(-0.3, 0.3),
    linewidth = 0.40,
    linetype = "dashed",
    colour = "#7A7A7A"
  ) +
  geom_point(
    size = 3.05,
    stroke = 0.30
  ) +
  direction_scale_s7b() +
  scale_x_continuous(
    limits = left_axis_limits_s7b,
    breaks = left_axis_breaks_s7b,
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(
      add = 0.48
    )
  ) +
  labs(
    title = "Direction-oriented\npaired change",
    subtitle = "Oriented by the bulk-defined direction",
    x = "Median direction-oriented change",
    y = NULL
  ) +
  common_theme_s7b +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom"
  ) +
  guides(
    colour = "none"
  )

# ------------------------------------------------------------
# 12. Directional-support plot
# ------------------------------------------------------------

plot_support_gene_s7b <- ggplot(
  plot_data_s7b,
  aes(
    x = support_percentage,
    y = gene_factor,
    colour = direction_call_display
  )
) +
  geom_hline(
    yintercept = separator_positions_s7b,
    linewidth = 0.45,
    colour = "#555555"
  ) +
  geom_vline(
    xintercept = c(
      100 / 3,
      100 * 2 / 3
    ),
    linewidth = 0.40,
    linetype = "dashed",
    colour = "#7A7A7A"
  ) +
  geom_point(
    size = 2.65,
    stroke = 0.25
  ) +
  direction_scale_s7b() +
  scale_x_continuous(
    limits = c(-2, 102),
    breaks = c(
      0,
      20,
      40,
      60,
      80,
      100
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(
      add = 0.48
    )
  ) +
  labs(
    title = "Positive change across\npaired donors",
    subtitle = "Oriented change > 0; n = 9",
    x = "Pairs with positive oriented change",
    y = NULL
  ) +
  common_theme_s7b +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom"
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 3.0
      )
    )
  )

# ------------------------------------------------------------
# 13. Assemble S7B
# ------------------------------------------------------------

plot_s7b <- patchwork::wrap_plots(
  plot_module_s7b,
  plot_delta_gene_s7b,
  plot_support_gene_s7b,
  nrow = 1,
  ncol = 3,
  widths = c(
    0.18,
    1.62,
    1.30
  ),
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Direction-oriented acute-convalescent changes in signature genes",
    subtitle = paste0(
      observed_call_counts_s7b[["Aligned"]], " aligned | ",
      observed_call_counts_s7b[["Neither criterion met"]],
      " neither criterion met | ",
      observed_call_counts_s7b[["Opposite direction"]],
      " opposite direction\n",
      "n = 9 acute-convalescent donor pairs"
    ),
    theme = theme(
      plot.title = element_text(
        family = base_family_s7b,
        size = 12.0,
        face = "bold",
        colour = "#111111",
        hjust = 0,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        family = base_family_s7b,
        size = 8.6,
        colour = "#555555",
        hjust = 0,
        margin = margin(b = 9)
      ),
      plot.margin = margin(
        t = 10,
        r = 10,
        b = 8,
        l = 10
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(
      family = base_family_s7b,
      size = 8.3,
      face = "bold"
    ),
    legend.text = element_text(
      family = base_family_s7b,
      size = 7.8
    )
  )

print(plot_s7b)

# ------------------------------------------------------------
# 14. Export source data
# ------------------------------------------------------------

source_data_file_s7b <- file.path(
  output_dir_s7b,
  "Supplementary_Figure_S9B_source_data.csv"
)

summary_file_s7b <- file.path(
  output_dir_s7b,
  "Supplementary_Figure_S9B_direction_summary.csv"
)

source_data_s7b <- plot_data_s7b

source_data_s7b$gene_factor <- as.character(
  source_data_s7b$gene_factor
)

source_data_s7b$module <- as.character(
  source_data_s7b$module
)

source_data_s7b$direction_call_display <-
  as.character(
    source_data_s7b$direction_call_display
  )

summary_data_s7b <- data.frame(
  direction_call = direction_call_order_s7b,
  gene_n = as.integer(observed_call_counts_s7b),
  stringsAsFactors = FALSE
)

write.csv(
  source_data_s7b,
  source_data_file_s7b,
  row.names = FALSE,
  na = ""
)

write.csv(
  summary_data_s7b,
  summary_file_s7b,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 15. Export PDF, SVG, TIFF and PNG
# ------------------------------------------------------------

file_stem_s7b <- file.path(
  output_dir_s7b,
  "Supplementary_Figure_S9B_paired_gene_direction"
)

width_mm_s7b <- 183
height_mm_s7b <- 180
dpi_s7b <- 600

width_in_s7b <- width_mm_s7b / 25.4
height_in_s7b <- height_mm_s7b / 25.4

# PDF
grDevices::pdf(
  file = paste0(file_stem_s7b, ".pdf"),
  width = width_in_s7b,
  height = height_in_s7b,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE
)

print(plot_s7b)
grDevices::dev.off()

# SVG
svglite::svglite(
  file = paste0(file_stem_s7b, ".svg"),
  width = width_in_s7b,
  height = height_in_s7b
)

print(plot_s7b)
grDevices::dev.off()

# TIFF
ragg::agg_tiff(
  filename = paste0(file_stem_s7b, ".tiff"),
  width = width_in_s7b,
  height = height_in_s7b,
  units = "in",
  res = dpi_s7b,
  background = "white",
  compression = "lzw"
)

print(plot_s7b)
grDevices::dev.off()

# PNG preview
ragg::agg_png(
  filename = paste0(file_stem_s7b, ".png"),
  width = width_in_s7b,
  height = height_in_s7b,
  units = "in",
  res = 300,
  background = "white"
)

print(plot_s7b)
grDevices::dev.off()

# ------------------------------------------------------------
# 16. Final output audit
# ------------------------------------------------------------

figure_files_s7b <- c(
  paste0(file_stem_s7b, ".pdf"),
  paste0(file_stem_s7b, ".svg"),
  paste0(file_stem_s7b, ".tiff"),
  paste0(file_stem_s7b, ".png")
)

output_audit_file_s7b <- file.path(
  output_dir_s7b,
  "Supplementary_Figure_S9B_output_audit.csv"
)

output_audit_s7b <- data.frame(
  output_type = c(
    "PDF figure",
    "SVG figure",
    "TIFF figure",
    "PNG preview",
    "Gene-level source data",
    "Direction summary"
  ),
  file = c(
    figure_files_s7b,
    source_data_file_s7b,
    summary_file_s7b
  ),
  stringsAsFactors = FALSE
)

output_audit_s7b$exists <- file.exists(
  output_audit_s7b$file
)

output_audit_s7b$size_bytes <- file.info(
  output_audit_s7b$file
)$size

write.csv(
  output_audit_s7b,
  output_audit_file_s7b,
  row.names = FALSE
)

expected_outputs_s7b <- c(
  output_audit_s7b$file,
  output_audit_file_s7b
)

stopifnot(
  all(file.exists(expected_outputs_s7b)),
  all(file.info(expected_outputs_s7b)$size > 0)
)

cat(
  "\n====================================================\n",
  "S7B COMPLETED SUCCESSFULLY\n",
  "====================================================\n",
  sep = ""
)

cat(
  "Output directory:\n",
  normalizePath(
    output_dir_s7b,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

cat("Direction-call counts:\n")
print(observed_call_counts_s7b)

cat(
  "\nOpposite-direction genes:\n",
  paste(
    observed_opposite_genes_s7b,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "\nNo panel letter was embedded in the figure.\n"
)





# ============================================================
# Supplementary Figure S7C
# Gene-state condition-context classification matrix
#
# No embedded panel letter
# Standalone from frozen CSV source data
# ============================================================

cat(
  "\n========== S7C: GENE-STATE CONDITION CONTEXT ==========\n"
)

# ------------------------------------------------------------
# 1. Package audit
# ------------------------------------------------------------

required_packages_s7c <- c(
  "ggplot2",
  "svglite",
  "ragg"
)

missing_packages_s7c <- required_packages_s7c[
  !vapply(
    required_packages_s7c,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_s7c) > 0L) {
  stop(
    paste0(
      "Missing packages required for S7C: ",
      paste(missing_packages_s7c, collapse = ", ")
    )
  )
}

library(ggplot2)

# ------------------------------------------------------------
# 2. Input and output paths
# ------------------------------------------------------------

source_file_s7c <- file.path(
  "/home/sunshine/predicate/singlecell",
  "outputs",
  "GSE216009_32gene_analysis",
  "04b_candidate_condition_direction",
  "20260804_candidate_condition_direction_2",
  "tables",
  "table04b_candidate_condition_direction.csv"
)

supplementary_root_s7c <- file.path(
  "/home/sunshine/predicate/singlecell",
  "outputs",
  "GSE216009_32gene_analysis",
  "manuscript_revision",
  "Supplementary_Figures_S6_S7"
)

figure_dir_s7c <- file.path(
  supplementary_root_s7c,
  "Supplementary_Figure_S7C_gene_state_condition_context"
)

dir.create(
  figure_dir_s7c,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(source_file_s7c)) {
  stop(
    paste0(
      "S7C source file was not found:\n",
      source_file_s7c
    )
  )
}

# ------------------------------------------------------------
# 3. Read frozen source data
# ------------------------------------------------------------

data_s7c <- utils::read.csv(
  source_file_s7c,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

required_columns_s7c <- c(
  "gene",
  "module",
  "main_contributor_state",
  "paired_n_evaluable",
  "paired_delta_signed_median",
  "condition_direction_call",
  "state_abundance_call",
  "condition_support_final"
)

missing_columns_s7c <- setdiff(
  required_columns_s7c,
  names(data_s7c)
)

if (length(missing_columns_s7c) > 0L) {
  stop(
    paste0(
      "Missing required S7C columns: ",
      paste(missing_columns_s7c, collapse = ", ")
    )
  )
}

stopifnot(
  nrow(data_s7c) == 24L,
  anyDuplicated(data_s7c$gene) == 0L
)

# ------------------------------------------------------------
# 4. Clean source fields
# ------------------------------------------------------------

text_columns_s7c <- c(
  "gene",
  "module",
  "main_contributor_state",
  "condition_direction_call",
  "state_abundance_call",
  "condition_support_final"
)

for (current_column_s7c in text_columns_s7c) {
  data_s7c[[current_column_s7c]] <- trimws(
    as.character(data_s7c[[current_column_s7c]])
  )
}

data_s7c$source_order_s7c <- seq_len(nrow(data_s7c))

# Fixed module order, consistent with the bulk and S7B figures.
module_order_s7c <- c(
  "Cell cycle",
  "Neutrophil",
  "Inflammatory",
  "Other"
)

if (any(!data_s7c$module %in% module_order_s7c)) {
  stop(
    paste0(
      "Unexpected module labels: ",
      paste(
        setdiff(
          unique(data_s7c$module),
          module_order_s7c
        ),
        collapse = ", "
      )
    )
  )
}

data_s7c <- data_s7c[
  order(
    match(data_s7c$module, module_order_s7c),
    data_s7c$source_order_s7c
  ),
  ,
  drop = FALSE
]

row.names(data_s7c) <- NULL

expected_gene_order_s7c <- c(
  # M1 cell cycle/proliferation
  "ANLN",
  "CCNA2",
  "CKS2",
  "KIF11",
  "MKI67",
  "RRM2",
  "TOP2A",
  "TPX2",
  "TYMS",
  
  # M2 neutrophil degranulation
  "CEACAM6",
  "CTSG",
  "DEFA4",
  "ELANE",
  "MPO",
  "MS4A3",
  "RNASE3",
  "TCN1",
  
  # Inflammatory/downregulated
  "CX3CR1",
  "IL1B",
  "SECTM1",
  
  # Other
  "ABCA13",
  "CD24",
  "CHI3L1",
  "SERPINB10"
)

if (!identical(data_s7c$gene, expected_gene_order_s7c)) {
  stop(
    paste0(
      "The reconstructed S7C gene order does not match ",
      "the prespecified module-grouped order."
    )
  )
}

# Top row receives the highest y coordinate.
data_s7c$plot_y_s7c <- rev(seq_len(nrow(data_s7c)))

# ------------------------------------------------------------
# 5. Display-label mappings
# ------------------------------------------------------------

state_label_map_s7c <- c(
  "MPO+_immature_neutrophils_or_progenitors" =
    "MPO+ imm neut/prog",
  "Classical_monocytes" =
    "Classical mono",
  "Mature_neutrophils" =
    "Mature neut",
  "S100A8-9_hi_neutrophils" =
    "S100A8/A9-high neut"
)

support_label_map_s7c <- c(
  "moderate_support" =
    "Moderate support",
  "not_supported" =
    "Unsupported",
  "opposite_direction" =
    "Opposite direction",
  "not_evaluable" =
    "NE"
)

abundance_label_map_s7c <- c(
  "acute_expanded" =
    "Acute-expanded",
  "acute_depleted" =
    "Acute-depleted",
  "not_consistent" =
    "Not consistent"
)

data_s7c$state_display_s7c <- unname(
  state_label_map_s7c[
    data_s7c$main_contributor_state
  ]
)

data_s7c$support_display_s7c <- unname(
  support_label_map_s7c[
    data_s7c$condition_support_final
  ]
)

data_s7c$abundance_display_s7c <- unname(
  abundance_label_map_s7c[
    data_s7c$state_abundance_call
  ]
)

if (
  any(is.na(data_s7c$state_display_s7c)) ||
  any(is.na(data_s7c$support_display_s7c)) ||
  any(is.na(data_s7c$abundance_display_s7c))
) {
  stop(
    paste0(
      "At least one S7C source category was not mapped ",
      "to a display label."
    )
  )
}

# ------------------------------------------------------------
# 6. Numerical display fields
# ------------------------------------------------------------

data_s7c$paired_n_label_s7c <- ifelse(
  is.na(data_s7c$paired_n_evaluable),
  "0",
  as.character(
    as.integer(data_s7c$paired_n_evaluable)
  )
)

data_s7c$delta_label_s7c <- "NE"

finite_delta_s7c <- is.finite(
  data_s7c$paired_delta_signed_median
)

data_s7c$delta_label_s7c[
  finite_delta_s7c
] <- sprintf(
  "%+.3f",
  data_s7c$paired_delta_signed_median[
    finite_delta_s7c
  ]
)

# ------------------------------------------------------------
# 7. Formal source-data assertions
# ------------------------------------------------------------

support_levels_source_s7c <- c(
  "moderate_support",
  "not_supported",
  "opposite_direction",
  "not_evaluable"
)

expected_support_counts_s7c <- c(
  3L,
  3L,
  1L,
  17L
)

observed_support_counts_s7c <- as.integer(
  table(
    factor(
      data_s7c$condition_support_final,
      levels = support_levels_source_s7c
    )
  )
)

abundance_levels_source_s7c <- c(
  "acute_expanded",
  "acute_depleted",
  "not_consistent"
)

expected_abundance_counts_s7c <- c(
  4L,
  3L,
  17L
)

observed_abundance_counts_s7c <- as.integer(
  table(
    factor(
      data_s7c$state_abundance_call,
      levels = abundance_levels_source_s7c
    )
  )
)

stopifnot(
  identical(
    observed_support_counts_s7c,
    expected_support_counts_s7c
  ),
  identical(
    observed_abundance_counts_s7c,
    expected_abundance_counts_s7c
  ),
  setequal(
    data_s7c$gene[
      data_s7c$condition_support_final ==
        "moderate_support"
    ],
    c("SERPINB10", "CX3CR1", "IL1B")
  ),
  setequal(
    data_s7c$gene[
      data_s7c$condition_support_final ==
        "not_supported"
    ],
    c("CKS2", "TCN1", "SECTM1")
  ),
  setequal(
    data_s7c$gene[
      data_s7c$condition_support_final ==
        "opposite_direction"
    ],
    "CHI3L1"
  ),
  sum(data_s7c$paired_n_evaluable == 9L) == 7L,
  sum(data_s7c$paired_n_evaluable == 0L) == 17L
)

# ------------------------------------------------------------
# 8. Colour contract
# ------------------------------------------------------------

# Module colours match the established S7B family.
module_colour_s7c <- c(
  "Cell cycle" = "#3F86C5",
  "Neutrophil" = "#F23B35",
  "Inflammatory" = "#59A14F",
  "Other" = "#AF7AA1"
)

module_block_label_s7c <- c(
  "Cell cycle" = "M1",
  "Neutrophil" = "M2",
  "Inflammatory" = "Inflam./down",
  "Other" = "Other"
)

# Within-state expression classification.
support_fill_s7c <- c(
  "moderate_support" = "#18788C",
  "not_supported" = "#C6CBD1",
  "opposite_direction" = "#A93670",
  "not_evaluable" = "#F5F5F5"
)

support_text_colour_s7c <- c(
  "moderate_support" = "#FFFFFF",
  "not_supported" = "#252525",
  "opposite_direction" = "#FFFFFF",
  "not_evaluable" = "#777777"
)

# Source-state abundance classification.
abundance_fill_s7c <- c(
  "acute_expanded" = "#E68A45",
  "acute_depleted" = "#4E7FAF",
  "not_consistent" = "#E2E2E2"
)

abundance_text_colour_s7c <- c(
  "acute_expanded" = "#FFFFFF",
  "acute_depleted" = "#FFFFFF",
  "not_consistent" = "#4F4F4F"
)

delta_text_colour_s7c <- c(
  "moderate_support" = "#12677A",
  "not_supported" = "#4F555B",
  "opposite_direction" = "#A12F68",
  "not_evaluable" = "#888888"
)

data_s7c$module_fill_s7c <- unname(
  module_colour_s7c[data_s7c$module]
)

data_s7c$support_fill_s7c <- unname(
  support_fill_s7c[
    data_s7c$condition_support_final
  ]
)

data_s7c$support_text_colour_s7c <- unname(
  support_text_colour_s7c[
    data_s7c$condition_support_final
  ]
)

data_s7c$abundance_fill_s7c <- unname(
  abundance_fill_s7c[
    data_s7c$state_abundance_call
  ]
)

data_s7c$abundance_text_colour_s7c <- unname(
  abundance_text_colour_s7c[
    data_s7c$state_abundance_call
  ]
)

data_s7c$delta_text_colour_s7c <- unname(
  delta_text_colour_s7c[
    data_s7c$condition_support_final
  ]
)

# ------------------------------------------------------------
# 9. Module-block geometry
# ------------------------------------------------------------

module_block_data_s7c <- do.call(
  rbind,
  lapply(
    module_order_s7c,
    function(current_module_s7c) {
      current_y_s7c <- data_s7c$plot_y_s7c[
        data_s7c$module == current_module_s7c
      ]
      
      data.frame(
        module = current_module_s7c,
        ymin = min(current_y_s7c) - 0.46,
        ymax = max(current_y_s7c) + 0.46,
        label_y = mean(current_y_s7c),
        stringsAsFactors = FALSE
      )
    }
  )
)

row.names(module_block_data_s7c) <- NULL

module_block_data_s7c$fill_colour <- unname(
  module_colour_s7c[
    module_block_data_s7c$module
  ]
)

module_block_data_s7c$block_label <- unname(
  module_block_label_s7c[
    module_block_data_s7c$module
  ]
)

module_counts_s7c <- as.integer(
  table(
    factor(
      data_s7c$module,
      levels = module_order_s7c
    )
  )
)

module_separator_y_s7c <- (
  nrow(data_s7c) -
    cumsum(
      module_counts_s7c[
        -length(module_counts_s7c)
      ]
    ) +
    0.5
)

# ------------------------------------------------------------
# 10. Column geometry
# ------------------------------------------------------------

x_module_s7c <- 0.22
x_gene_s7c <- 0.55
x_state_s7c <- 2.15
x_pairs_s7c <- 6.00
x_expression_s7c <- 7.85
x_abundance_s7c <- 10.65
x_delta_s7c <- 13.40

header_data_s7c <- data.frame(
  x = c(
    x_gene_s7c,
    x_state_s7c,
    x_pairs_s7c,
    x_expression_s7c,
    x_abundance_s7c,
    x_delta_s7c
  ),
  y = rep(25.32, 6L),
  label = c(
    "Gene / module",
    "Highest-contributor\nstate",
    "Paired\ndonors, n",
    "Within-state\nexpression call",
    "Source-state\nabundance call",
    "Median\noriented delta"
  ),
  hjust = c(
    0,
    0,
    0.5,
    0.5,
    0.5,
    0.5
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 11. Construct the classification matrix
# ------------------------------------------------------------

plot_s7c <- ggplot() +
  
  # Fine row guides.
  geom_hline(
    yintercept = seq(
      0.5,
      nrow(data_s7c) + 0.5,
      by = 1
    ),
    linewidth = 0.28,
    colour = "#E5E5E5"
  ) +
  
  # Module separators.
  geom_hline(
    yintercept = module_separator_y_s7c,
    linewidth = 0.72,
    colour = "#585858"
  ) +
  
  # Major column separators.
  geom_vline(
    xintercept = c(
      1.93,
      5.38,
      6.63,
      9.22,
      12.08
    ),
    linewidth = 0.30,
    colour = "#D1D1D1"
  ) +
  
  # Module blocks.
  geom_rect(
    data = module_block_data_s7c,
    aes(
      xmin = 0.03,
      xmax = 0.40,
      ymin = ymin,
      ymax = ymax,
      fill = fill_colour
    ),
    colour = NA
  ) +
  
  geom_text(
    data = module_block_data_s7c,
    aes(
      x = x_module_s7c,
      y = label_y,
      label = block_label
    ),
    angle = 90,
    colour = "white",
    fontface = "bold",
    family = "sans",
    size = 3.55
  ) +
  
  # Neutral cells for paired-donor n.
  geom_tile(
    data = data_s7c,
    aes(
      x = x_pairs_s7c,
      y = plot_y_s7c
    ),
    width = 0.94,
    height = 0.82,
    fill = "#FAFAFA",
    colour = "#D8D8D8",
    linewidth = 0.30
  ) +
  
  # Within-state expression cells.
  geom_tile(
    data = data_s7c,
    aes(
      x = x_expression_s7c,
      y = plot_y_s7c,
      fill = support_fill_s7c
    ),
    width = 2.36,
    height = 0.82,
    colour = "#D0D0D0",
    linewidth = 0.30
  ) +
  
  # Source-state abundance cells.
  geom_tile(
    data = data_s7c,
    aes(
      x = x_abundance_s7c,
      y = plot_y_s7c,
      fill = abundance_fill_s7c
    ),
    width = 2.42,
    height = 0.82,
    colour = "#D0D0D0",
    linewidth = 0.30
  ) +
  
  # Neutral cells for numerical oriented delta.
  geom_tile(
    data = data_s7c,
    aes(
      x = x_delta_s7c,
      y = plot_y_s7c
    ),
    width = 1.52,
    height = 0.82,
    fill = "#FAFAFA",
    colour = "#D8D8D8",
    linewidth = 0.30
  ) +
  
  # Gene labels.
  geom_text(
    data = data_s7c,
    aes(
      x = x_gene_s7c,
      y = plot_y_s7c,
      label = gene
    ),
    hjust = 0,
    family = "sans",
    fontface = "italic",
    colour = "#202020",
    size = 4.10
  ) +
  
  # Short highest-contributor state labels.
  geom_text(
    data = data_s7c,
    aes(
      x = x_state_s7c,
      y = plot_y_s7c,
      label = state_display_s7c
    ),
    hjust = 0,
    family = "sans",
    colour = "#303030",
    size = 3.80
  ) +
  
  # Paired-donor counts.
  geom_text(
    data = data_s7c,
    aes(
      x = x_pairs_s7c,
      y = plot_y_s7c,
      label = paired_n_label_s7c
    ),
    family = "sans",
    fontface = "bold",
    colour = "#303030",
    size = 3.90
  ) +
  
  # Within-state expression labels.
  geom_text(
    data = data_s7c,
    aes(
      x = x_expression_s7c,
      y = plot_y_s7c,
      label = support_display_s7c,
      colour = support_text_colour_s7c
    ),
    family = "sans",
    fontface = "bold",
    size = 3.70
  ) +
  
  # Source-state abundance labels.
  geom_text(
    data = data_s7c,
    aes(
      x = x_abundance_s7c,
      y = plot_y_s7c,
      label = abundance_display_s7c,
      colour = abundance_text_colour_s7c
    ),
    family = "sans",
    fontface = "bold",
    size = 3.55
  ) +
  
  # Median oriented delta.
  geom_text(
    data = data_s7c,
    aes(
      x = x_delta_s7c,
      y = plot_y_s7c,
      label = delta_label_s7c,
      colour = delta_text_colour_s7c
    ),
    family = "sans",
    fontface = "bold",
    size = 3.75
  ) +
  
  # Column headings.
  geom_text(
    data = header_data_s7c,
    aes(
      x = x,
      y = y,
      label = label,
      hjust = hjust
    ),
    family = "sans",
    fontface = "bold",
    colour = "#171717",
    lineheight = 0.96,
    size = 4.45
  ) +
  
  # Header rule.
  annotate(
    geom = "segment",
    x = 0.03,
    xend = 14.18,
    y = 24.68,
    yend = 24.68,
    linewidth = 0.75,
    colour = "#303030"
  ) +
  
  scale_fill_identity() +
  scale_colour_identity() +
  
  scale_x_continuous(
    limits = c(0, 14.22),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0.35, 25.75),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    title = paste0(
      "Condition-context analysis resolved seven ",
      "evaluable gene–state comparisons"
    ),
    subtitle = paste0(
      "Within-state expression: 3 moderate support | ",
      "3 unsupported | 1 opposite direction | ",
      "17 not evaluable\n",
      "Source-state abundance: 4 acute-expanded | ",
      "3 acute-depleted | 17 not consistently changed"
    ),
    caption = paste0(
      "Positive oriented delta supports the bulk-defined direction. ",
      "NE, not evaluable. Module strip: M1, cell cycle/proliferation; ",
      "M2, neutrophil degranulation; Inflam./down, ",
      "inflammatory/downregulated."
    )
  ) +
  
  theme_void(
    base_size = 12.5,
    base_family = "sans"
  ) +
  
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.title = element_text(
      size = 18.0,
      face = "bold",
      colour = "#111111",
      hjust = 0,
      margin = margin(
        b = 7
      )
    ),
    plot.subtitle = element_text(
      size = 12.2,
      colour = "#4E4E4E",
      hjust = 0,
      lineheight = 1.12,
      margin = margin(
        b = 14
      )
    ),
    plot.caption = element_text(
      size = 10.1,
      colour = "#565656",
      hjust = 0,
      lineheight = 1.08,
      margin = margin(
        t = 12
      )
    ),
    plot.margin = margin(
      t = 22,
      r = 24,
      b = 20,
      l = 22
    )
  )

stopifnot(
  inherits(plot_s7c, "ggplot")
)

print(plot_s7c)

# ------------------------------------------------------------
# 12. Export source data
# ------------------------------------------------------------

source_export_s7c <- data.frame(
  display_order = seq_len(nrow(data_s7c)),
  gene = data_s7c$gene,
  module = data_s7c$module,
  highest_contributor_state =
    data_s7c$main_contributor_state,
  highest_contributor_state_display =
    data_s7c$state_display_s7c,
  paired_n_evaluable =
    data_s7c$paired_n_evaluable,
  median_direction_oriented_delta =
    data_s7c$paired_delta_signed_median,
  within_state_expression_call =
    data_s7c$condition_support_final,
  source_state_abundance_call =
    data_s7c$state_abundance_call,
  stringsAsFactors = FALSE
)

source_data_file_s7c <- file.path(
  figure_dir_s7c,
  "Supplementary_Figure_S7C_source_data.csv"
)

utils::write.csv(
  source_export_s7c,
  source_data_file_s7c,
  row.names = FALSE,
  na = ""
)

summary_export_s7c <- data.frame(
  metric = c(
    "Prespecified gene-state comparisons",
    "Moderate support",
    "Unsupported",
    "Opposite direction",
    "Not evaluable",
    "Acute-expanded source state",
    "Acute-depleted source state",
    "Source state not consistently changed"
  ),
  value = c(
    24L,
    3L,
    3L,
    1L,
    17L,
    4L,
    3L,
    17L
  ),
  stringsAsFactors = FALSE
)

summary_file_s7c <- file.path(
  figure_dir_s7c,
  "Supplementary_Figure_S7C_summary.csv"
)

utils::write.csv(
  summary_export_s7c,
  summary_file_s7c,
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Export figure
# ------------------------------------------------------------

figure_width_s7c <- 14.4
figure_height_s7c <- 12.8

pdf_file_s7c <- file.path(
  figure_dir_s7c,
  "Supplementary_Figure_S7C_gene_state_condition_context.pdf"
)

svg_file_s7c <- file.path(
  figure_dir_s7c,
  "Supplementary_Figure_S7C_gene_state_condition_context.svg"
)

tiff_file_s7c <- file.path(
  figure_dir_s7c,
  "Supplementary_Figure_S7C_gene_state_condition_context.tiff"
)

png_file_s7c <- file.path(
  figure_dir_s7c,
  "Supplementary_Figure_S7C_gene_state_condition_context.png"
)

# PDF: do not pass compression = "lzw" to grDevices::pdf().
grDevices::pdf(
  file = pdf_file_s7c,
  width = figure_width_s7c,
  height = figure_height_s7c,
  family = "sans",
  onefile = TRUE,
  useDingbats = FALSE
)

print(plot_s7c)

grDevices::dev.off()

svglite::svglite(
  file = svg_file_s7c,
  width = figure_width_s7c,
  height = figure_height_s7c,
  bg = "white"
)

print(plot_s7c)

grDevices::dev.off()

ragg::agg_tiff(
  filename = tiff_file_s7c,
  width = figure_width_s7c,
  height = figure_height_s7c,
  units = "in",
  res = 600,
  scaling = 1,
  compression = "lzw",
  background = "white"
)

print(plot_s7c)

grDevices::dev.off()

ragg::agg_png(
  filename = png_file_s7c,
  width = figure_width_s7c,
  height = figure_height_s7c,
  units = "in",
  res = 300,
  scaling = 1,
  background = "white"
)

print(plot_s7c)

grDevices::dev.off()

# ------------------------------------------------------------
# 14. Output audit
# ------------------------------------------------------------

expected_output_files_s7c <- c(
  pdf_file_s7c,
  svg_file_s7c,
  tiff_file_s7c,
  png_file_s7c,
  source_data_file_s7c,
  summary_file_s7c
)

if (!all(file.exists(expected_output_files_s7c))) {
  stop(
    paste0(
      "One or more expected S7C output files ",
      "were not created."
    )
  )
}

output_info_s7c <- file.info(
  expected_output_files_s7c
)

output_audit_s7c <- data.frame(
  file = normalizePath(
    expected_output_files_s7c,
    mustWork = TRUE
  ),
  size_bytes = output_info_s7c$size,
  modified = output_info_s7c$mtime,
  row.names = NULL,
  stringsAsFactors = FALSE
)

cat("\nS7C output directory:\n")
cat(
  normalizePath(
    figure_dir_s7c,
    mustWork = TRUE
  ),
  "\n"
)

cat("\nS7C output files:\n")
print(
  output_audit_s7c,
  row.names = FALSE
)

cat(
  paste0(
    "\nPASS: S7C retained all 24 prespecified ",
    "gene-state comparisons.\n"
  )
)

cat(
  paste0(
    "PASS: the figure displays 7 evaluable and ",
    "17 non-evaluable comparisons explicitly.\n"
  )
)

cat(
  paste0(
    "PASS: PDF, SVG, TIFF, PNG and source-data ",
    "files were created.\n"
  )
)

cat(
  "====================================================\n"
)





# ============================================================
# Supplementary Figure S7D
# Immature-neutrophil state-intrinsic paired analysis
#
# No embedded panel letter
# Standalone from frozen CSV source data
# ============================================================

cat(
  "\n========== S7D: STATE-INTRINSIC PAIRED ANALYSIS ==========\n"
)

# ------------------------------------------------------------
# 1. Package audit
# ------------------------------------------------------------

required_packages_s7d <- c(
  "ggplot2",
  "patchwork",
  "svglite",
  "ragg"
)

missing_packages_s7d <- required_packages_s7d[
  !vapply(
    required_packages_s7d,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_s7d) > 0L) {
  stop(
    paste0(
      "Missing packages required for S7D: ",
      paste(missing_packages_s7d, collapse = ", ")
    )
  )
}

library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 2. Input paths
# ------------------------------------------------------------

table_dir_s7d <- file.path(
  "/home/sunshine/predicate/singlecell",
  "outputs",
  "GSE216009_32gene_analysis",
  "05e_pseudobulk_sanity_state_direction",
  "20260805_internal_pseudobulk_sanity_and_state_direction_3",
  "tables"
)

direction_file_s7d <- file.path(
  table_dir_s7d,
  "p3b_direction_calls.csv"
)

paired_file_s7d <- file.path(
  table_dir_s7d,
  "p3b_paired_delta_long.csv"
)

condition_file_s7d <- file.path(
  table_dir_s7d,
  "p3b_state_readout_condition_summary.csv"
)

required_input_files_s7d <- c(
  direction_file_s7d,
  paired_file_s7d,
  condition_file_s7d
)

if (!all(file.exists(required_input_files_s7d))) {
  stop(
    paste0(
      "One or more required S7D input files were not found:\n",
      paste(
        required_input_files_s7d[
          !file.exists(required_input_files_s7d)
        ],
        collapse = "\n"
      )
    )
  )
}

# ------------------------------------------------------------
# 3. Output directory
# ------------------------------------------------------------

supplementary_root_s7d <- file.path(
  "/home/sunshine/predicate/singlecell",
  "outputs",
  "GSE216009_32gene_analysis",
  "manuscript_revision",
  "Supplementary_Figures_S6_S7"
)

figure_dir_s7d <- file.path(
  supplementary_root_s7d,
  "Supplementary_Figure_S7D_state_intrinsic_analysis"
)

dir.create(
  figure_dir_s7d,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 4. Read frozen source data
# ------------------------------------------------------------

direction_data_s7d <- utils::read.csv(
  direction_file_s7d,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

paired_data_s7d <- utils::read.csv(
  paired_file_s7d,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

condition_data_s7d <- utils::read.csv(
  condition_file_s7d,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

stopifnot(
  nrow(direction_data_s7d) == 28L,
  nrow(paired_data_s7d) == 63L,
  nrow(condition_data_s7d) == 28L
)

required_direction_columns_s7d <- c(
  "fine_state",
  "readout",
  "result_layer",
  "paired_delta_n",
  "paired_acute_eval_n",
  "paired_conv_eval_n",
  "paired_delta_median",
  "support_rate_delta_gt0",
  "direction_call",
  "not_evaluable_reason"
)

required_paired_columns_s7d <- c(
  "donor",
  "delta",
  "acute_value",
  "conv_value",
  "fine_state",
  "readout"
)

required_condition_columns_s7d <- c(
  "fine_state",
  "readout",
  "acute_median",
  "conv_median",
  "hv_median",
  "cs_median",
  "paired_delta_n",
  "paired_delta_median",
  "direction_call"
)

if (
  length(
    setdiff(
      required_direction_columns_s7d,
      names(direction_data_s7d)
    )
  ) > 0L ||
  length(
    setdiff(
      required_paired_columns_s7d,
      names(paired_data_s7d)
    )
  ) > 0L ||
  length(
    setdiff(
      required_condition_columns_s7d,
      names(condition_data_s7d)
    )
  ) > 0L
) {
  stop("One or more required S7D source columns are missing.")
}

# ------------------------------------------------------------
# 5. Prespecified state and readout order
# ------------------------------------------------------------

state_order_s7d <- c(
  "MPO+_immature_neutrophils_or_progenitors",
  "PADI4+_immature_neutrophils",
  "Cycling_neutrophil_progenitors",
  "S100A8-9_hi_neutrophils"
)

readout_order_s7d <- c(
  "M2_full",
  "primary_core",
  "secondary_granule",
  "tertiary_granule",
  "mature_surface",
  "inflammatory_S100",
  "M1_cell_cycle"
)

state_display_map_s7d <- c(
  "MPO+_immature_neutrophils_or_progenitors" =
    "MPO+ immature neutrophils/\nprogenitors",
  "PADI4+_immature_neutrophils" =
    "PADI4+ immature\nneutrophils",
  "Cycling_neutrophil_progenitors" =
    "Cycling neutrophil\nprogenitors",
  "S100A8-9_hi_neutrophils" =
    "S100A8/A9-high\nneutrophils"
)

readout_display_map_s7d <- c(
  "M2_full" =
    "M2 full",
  "primary_core" =
    "Primary\ncore",
  "secondary_granule" =
    "Secondary\ngranule",
  "tertiary_granule" =
    "Tertiary\ngranule",
  "mature_surface" =
    "Mature\nsurface",
  "inflammatory_S100" =
    "Inflammatory\nS100",
  "M1_cell_cycle" =
    "M1 cell\ncycle"
)

stopifnot(
  setequal(
    unique(direction_data_s7d$fine_state),
    state_order_s7d
  ),
  setequal(
    unique(direction_data_s7d$readout),
    readout_order_s7d
  ),
  setequal(
    unique(paired_data_s7d$readout),
    readout_order_s7d
  ),
  identical(
    unique(paired_data_s7d$fine_state),
    "S100A8-9_hi_neutrophils"
  )
)

# ------------------------------------------------------------
# 6. Prepare the 4-state x 7-readout matrix
# ------------------------------------------------------------

matrix_data_s7d <- direction_data_s7d

matrix_data_s7d$state_order_index_s7d <- match(
  matrix_data_s7d$fine_state,
  state_order_s7d
)

matrix_data_s7d$readout_order_index_s7d <- match(
  matrix_data_s7d$readout,
  readout_order_s7d
)

matrix_data_s7d <- matrix_data_s7d[
  order(
    matrix_data_s7d$state_order_index_s7d,
    matrix_data_s7d$readout_order_index_s7d
  ),
  ,
  drop = FALSE
]

row.names(matrix_data_s7d) <- NULL

matrix_data_s7d$state_display_s7d <- unname(
  state_display_map_s7d[
    matrix_data_s7d$fine_state
  ]
)

matrix_data_s7d$readout_display_s7d <- unname(
  readout_display_map_s7d[
    matrix_data_s7d$readout
  ]
)

matrix_data_s7d$state_display_s7d <- factor(
  matrix_data_s7d$state_display_s7d,
  levels = rev(
    unname(
      state_display_map_s7d[state_order_s7d]
    )
  )
)

matrix_data_s7d$readout_display_s7d <- factor(
  matrix_data_s7d$readout_display_s7d,
  levels = unname(
    readout_display_map_s7d[
      readout_order_s7d
    ]
  )
)

matrix_evaluable_s7d <- matrix_data_s7d[
  matrix_data_s7d$paired_delta_n > 0L,
  ,
  drop = FALSE
]

matrix_not_evaluable_s7d <- matrix_data_s7d[
  matrix_data_s7d$paired_delta_n == 0L,
  ,
  drop = FALSE
]

matrix_evaluable_s7d$delta_label_s7d <- sprintf(
  "%+.3f",
  matrix_evaluable_s7d$paired_delta_median
)

matrix_evaluable_s7d$text_colour_s7d <- ifelse(
  abs(matrix_evaluable_s7d$paired_delta_median) >= 0.8,
  "#FFFFFF",
  "#202020"
)

matrix_m2_outline_s7d <- matrix_data_s7d[
  matrix_data_s7d$readout == "M2_full",
  ,
  drop = FALSE
]

stopifnot(
  nrow(matrix_data_s7d) == 28L,
  nrow(matrix_evaluable_s7d) == 7L,
  nrow(matrix_not_evaluable_s7d) == 21L,
  all(
    matrix_not_evaluable_s7d$direction_call ==
      "not_evaluable"
  )
)

# ------------------------------------------------------------
# 7. Reconstruct donor-level median and IQR summaries
# ------------------------------------------------------------

paired_summary_s7d <- do.call(
  rbind,
  lapply(
    readout_order_s7d,
    function(current_readout_s7d) {
      current_data_s7d <- paired_data_s7d[
        paired_data_s7d$readout ==
          current_readout_s7d,
        ,
        drop = FALSE
      ]
      
      data.frame(
        readout = current_readout_s7d,
        n_pairs = nrow(current_data_s7d),
        median_delta = stats::median(
          current_data_s7d$delta
        ),
        Q1_delta = unname(
          stats::quantile(
            current_data_s7d$delta,
            probs = 0.25,
            type = 7
          )
        ),
        Q3_delta = unname(
          stats::quantile(
            current_data_s7d$delta,
            probs = 0.75,
            type = 7
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

formal_s100_rows_s7d <- matrix_data_s7d[
  matrix_data_s7d$fine_state ==
    "S100A8-9_hi_neutrophils",
  ,
  drop = FALSE
]

formal_s100_rows_s7d <- formal_s100_rows_s7d[
  match(
    readout_order_s7d,
    formal_s100_rows_s7d$readout
  ),
  ,
  drop = FALSE
]

stopifnot(
  all(paired_summary_s7d$n_pairs == 9L),
  all(
    abs(
      paired_summary_s7d$median_delta -
        formal_s100_rows_s7d$paired_delta_median
    ) < 1e-10
  )
)

expected_medians_s7d <- c(
  -0.166,
  -0.663,
  0.790,
  1.459,
  -0.532,
  1.065,
  -0.633
)

stopifnot(
  all(
    abs(
      paired_summary_s7d$median_delta -
        expected_medians_s7d
    ) < 0.005
  )
)

paired_summary_s7d$readout_display_s7d <- factor(
  unname(
    readout_display_map_s7d[
      paired_summary_s7d$readout
    ]
  ),
  levels = rev(
    unname(
      readout_display_map_s7d[
        readout_order_s7d
      ]
    )
  )
)

paired_summary_s7d$summary_colour_s7d <- ifelse(
  paired_summary_s7d$readout == "M2_full",
  "#EF3B33",
  "#687078"
)

paired_summary_s7d$summary_fill_s7d <- ifelse(
  paired_summary_s7d$readout == "M2_full",
  "#EF3B33",
  "#7E858C"
)

paired_summary_s7d$median_label_s7d <- sprintf(
  "%+.3f",
  paired_summary_s7d$median_delta
)


paired_data_s7d$readout_display_s7d <- factor(
  unname(
    readout_display_map_s7d[
      paired_data_s7d$readout
    ]
  ),
  levels = rev(
    unname(
      readout_display_map_s7d[
        readout_order_s7d
      ]
    )
  )
)

paired_data_s7d$point_fill_s7d <- ifelse(
  paired_data_s7d$readout == "M2_full",
  "#EF3B33",
  "#AEB4BA"
)

# ------------------------------------------------------------
# 8. Left plot: evaluability and median-delta matrix
# ------------------------------------------------------------

matrix_plot_s7d <- ggplot() +
  
  geom_tile(
    data = matrix_not_evaluable_s7d,
    aes(
      x = readout_display_s7d,
      y = state_display_s7d
    ),
    width = 0.92,
    height = 0.84,
    fill = "#F3F3F3",
    colour = "#D0D0D0",
    linewidth = 0.45
  ) +
  
  geom_tile(
    data = matrix_evaluable_s7d,
    aes(
      x = readout_display_s7d,
      y = state_display_s7d,
      fill = paired_delta_median
    ),
    width = 0.92,
    height = 0.84,
    colour = "#BCBCBC",
    linewidth = 0.45
  ) +
  
  # Highlight the focal M2-full column without changing sign coding.
  geom_tile(
    data = matrix_m2_outline_s7d,
    aes(
      x = readout_display_s7d,
      y = state_display_s7d
    ),
    width = 0.94,
    height = 0.86,
    fill = NA,
    colour = "#EF3B33",
    linewidth = 1.05
  ) +
  
  geom_text(
    data = matrix_not_evaluable_s7d,
    aes(
      x = readout_display_s7d,
      y = state_display_s7d
    ),
    label = "NE",
    family = "sans",
    fontface = "bold",
    colour = "#858585",
    size = 4.20
  ) +
  
  geom_text(
    data = matrix_evaluable_s7d,
    aes(
      x = readout_display_s7d,
      y = state_display_s7d,
      label = delta_label_s7d,
      colour = text_colour_s7d
    ),
    family = "sans",
    fontface = "bold",
    size = 4.05
  ) +
  
  scale_fill_gradient2(
    name = "Median paired delta",
    low = "#4878A8",
    mid = "#F7F7F7",
    high = "#E76F51",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    breaks = c(-1.5, 0, 1.5),
    labels = c(
      "−1.5",
      "0",
      "+1.5"
    ),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(6.0, "cm"),
      barheight = grid::unit(0.42, "cm")
    )
  ) +
  
  scale_colour_identity() +
  
  scale_x_discrete(
    drop = FALSE,
    expand = expansion(
      add = 0.08
    )
  ) +
  
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(
      add = 0.10
    )
  ) +
  
  labs(
    title = "State evaluability and median paired change",
    subtitle = paste0(
      "Cell values are median acute − convalescent delta; ",
      "the focal M2 composite is outlined in red"
    ),
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 12.5,
    base_family = "sans"
  ) +
  
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    axis.text.x = element_text(
      size = 10.6,
      colour = "#202020",
      face = "bold",
      lineheight = 0.96,
      margin = margin(
        t = 8
      )
    ),
    axis.text.y = element_text(
      size = 10.7,
      colour = "#202020",
      lineheight = 0.98,
      margin = margin(
        r = 8
      )
    ),
    plot.title = element_text(
      size = 14.5,
      face = "bold",
      colour = "#111111",
      hjust = 0,
      margin = margin(
        b = 5
      )
    ),
    plot.subtitle = element_text(
      size = 10.9,
      colour = "#555555",
      hjust = 0,
      lineheight = 1.05,
      margin = margin(
        b = 10
      )
    ),
    legend.position = "bottom",
    legend.title = element_text(
      size = 10.8,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10.1
    ),
    plot.margin = margin(
      t = 8,
      r = 14,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------
# 9. Right plot: donor-level paired deltas
# ------------------------------------------------------------

delta_plot_s7d <- ggplot() +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.75,
    linetype = "dashed",
    colour = "#565656"
  ) +
  
  # IQR for each readout.
  geom_segment(
    data = paired_summary_s7d,
    aes(
      x = Q1_delta,
      xend = Q3_delta,
      y = readout_display_s7d,
      yend = readout_display_s7d,
      colour = summary_colour_s7d
    ),
    linewidth = 2.5,
    lineend = "round"
  ) +
  
  # Nine paired donors per readout.
  geom_point(
    data = paired_data_s7d,
    aes(
      x = delta,
      y = readout_display_s7d,
      fill = point_fill_s7d
    ),
    position = position_jitter(
      width = 0,
      height = 0.095,
      seed = 20260905
    ),
    shape = 21,
    size = 3.05,
    stroke = 0.38,
    colour = "#3E3E3E",
    alpha = 0.78
  ) +
  
  # Median diamond.
  geom_point(
    data = paired_summary_s7d,
    aes(
      x = median_delta,
      y = readout_display_s7d,
      fill = summary_fill_s7d
    ),
    shape = 23,
    size = 5.0,
    stroke = 0.85,
    colour = "#202020"
  ) +
  
  geom_text(
    data = paired_summary_s7d,
    aes(
      x = median_delta,
      y = readout_display_s7d,
      label = median_label_s7d,
      colour = summary_colour_s7d
    ),
    position = position_nudge(
      x = 0,
      y = 0.27
    ),
    hjust = 0.5,
    vjust = 0.5,
    family = "sans",
    fontface = "bold",
    size = 4.05
  ) +
  
  scale_colour_identity() +
  scale_fill_identity() +
  
  scale_x_continuous(
    limits = c(-6, 6),
    breaks = c(
      -6,
      -3,
      0,
      3,
      6
    ),
    labels = c(
      "−6",
      "−3",
      "0",
      "+3",
      "+6"
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(
      add = c(
        0.45,
        0.75
      )
    )
  ) +
  
  labs(
    title = "S100A8/A9-high donor-level paired changes",
    subtitle = paste0(
      "Points are paired donors; diamonds are medians; ",
      "horizontal bars are IQRs; n = 9"
    ),
    x = expression(
      paste(
        "Acute ",
        "\u2212",
        " convalescent ",
        Delta,
        " log"[2],
        "(CPM + 1)"
      )
    ),
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12.5,
    base_family = "sans"
  ) +
  
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid.major.x = element_line(
      linewidth = 0.35,
      colour = "#E4E4E4"
    ),
    panel.grid.major.y = element_line(
      linewidth = 0.35,
      colour = "#EEEEEE"
    ),
    axis.line = element_line(
      linewidth = 0.65,
      colour = "#252525"
    ),
    axis.ticks = element_line(
      linewidth = 0.55,
      colour = "#252525"
    ),
    axis.text.x = element_text(
      size = 10.8,
      colour = "#202020"
    ),
    axis.text.y = element_text(
      size = 11.0,
      colour = "#202020",
      lineheight = 0.98,
      margin = margin(
        r = 7
      )
    ),
    axis.title.x = element_text(
      size = 11.8,
      colour = "#202020",
      margin = margin(
        t = 10
      )
    ),
    plot.title = element_text(
      size = 14.5,
      face = "bold",
      colour = "#111111",
      hjust = 0,
      margin = margin(
        b = 5
      )
    ),
    plot.subtitle = element_text(
      size = 10.9,
      colour = "#555555",
      hjust = 0,
      lineheight = 1.05,
      margin = margin(
        b = 10
      )
    ),
    legend.position = "none",
    plot.margin = margin(
      t = 8,
      r = 22,
      b = 8,
      l = 14
    )
  )

# ------------------------------------------------------------
# 10. Assemble final S7D figure
# ------------------------------------------------------------

plot_s7d <- (
  matrix_plot_s7d |
    delta_plot_s7d
) +
  
  plot_layout(
    widths = c(
      1.16,
      1.08
    ),
    guides = "collect"
  ) +
  
  plot_annotation(
    title = paste0(
      "S100A8/A9-high neutrophils showed mixed ",
      "state-intrinsic acute–convalescent changes"
    ),
    subtitle = paste0(
      "M2 full showed no clear paired direction; ",
      "the other three prespecified immature-neutrophil ",
      "states were not evaluable"
    ),
    caption = paste0(
      "Delta represents acute minus convalescent expression. ",
      "NE, not evaluable because no complete paired comparison ",
      "met the prespecified evaluability requirement."
    ),
    theme = theme(
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      plot.title = element_text(
        family = "sans",
        size = 19.0,
        face = "bold",
        colour = "#111111",
        hjust = 0,
        margin = margin(
          b = 7
        )
      ),
      plot.subtitle = element_text(
        family = "sans",
        size = 12.5,
        colour = "#505050",
        hjust = 0,
        lineheight = 1.10,
        margin = margin(
          b = 14
        )
      ),
      plot.caption = element_text(
        family = "sans",
        size = 10.5,
        colour = "#555555",
        hjust = 0,
        lineheight = 1.08,
        margin = margin(
          t = 12
        )
      ),
      plot.margin = margin(
        t = 20,
        r = 22,
        b = 18,
        l = 20
      )
    )
  )

plot_s7d <- plot_s7d &
  theme(
    legend.position = "bottom"
  )

stopifnot(
  inherits(plot_s7d, "patchwork")
)

print(plot_s7d)

# ------------------------------------------------------------
# 11. Export source-data tables
# ------------------------------------------------------------

matrix_source_export_s7d <- matrix_data_s7d[
  ,
  c(
    "fine_state",
    "readout",
    "result_layer",
    "paired_delta_n",
    "paired_acute_eval_n",
    "paired_conv_eval_n",
    "paired_delta_median",
    "support_rate_delta_gt0",
    "direction_call",
    "not_evaluable_reason"
  ),
  drop = FALSE
]

paired_source_export_s7d <- paired_data_s7d[
  ,
  c(
    "donor",
    "fine_state",
    "readout",
    "acute_value",
    "conv_value",
    "delta"
  ),
  drop = FALSE
]

summary_source_export_s7d <- paired_summary_s7d[
  ,
  c(
    "readout",
    "n_pairs",
    "median_delta",
    "Q1_delta",
    "Q3_delta"
  ),
  drop = FALSE
]

condition_source_export_s7d <- condition_data_s7d[
  condition_data_s7d$fine_state %in%
    state_order_s7d,
  required_condition_columns_s7d,
  drop = FALSE
]

matrix_source_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_matrix_source_data.csv"
)

paired_source_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_paired_donor_source_data.csv"
)

summary_source_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_paired_summary.csv"
)

condition_source_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_condition_summary.csv"
)

utils::write.csv(
  matrix_source_export_s7d,
  matrix_source_file_s7d,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  paired_source_export_s7d,
  paired_source_file_s7d,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  summary_source_export_s7d,
  summary_source_file_s7d,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  condition_source_export_s7d,
  condition_source_file_s7d,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 12. Export figure
# ------------------------------------------------------------

figure_width_s7d <- 16.2
figure_height_s7d <- 9.4

pdf_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_state_intrinsic_analysis.pdf"
)

svg_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_state_intrinsic_analysis.svg"
)

tiff_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_state_intrinsic_analysis.tiff"
)

png_file_s7d <- file.path(
  figure_dir_s7d,
  "Supplementary_Figure_S7D_state_intrinsic_analysis.png"
)

# PDF: compression must not be passed to grDevices::pdf().
grDevices::pdf(
  file = pdf_file_s7d,
  width = figure_width_s7d,
  height = figure_height_s7d,
  family = "sans",
  onefile = TRUE,
  useDingbats = FALSE
)

print(plot_s7d)

grDevices::dev.off()

svglite::svglite(
  file = svg_file_s7d,
  width = figure_width_s7d,
  height = figure_height_s7d,
  bg = "white"
)

print(plot_s7d)

grDevices::dev.off()

ragg::agg_tiff(
  filename = tiff_file_s7d,
  width = figure_width_s7d,
  height = figure_height_s7d,
  units = "in",
  res = 600,
  scaling = 1,
  compression = "lzw",
  background = "white"
)

print(plot_s7d)

grDevices::dev.off()

ragg::agg_png(
  filename = png_file_s7d,
  width = figure_width_s7d,
  height = figure_height_s7d,
  units = "in",
  res = 300,
  scaling = 1,
  background = "white"
)

print(plot_s7d)

grDevices::dev.off()

# ------------------------------------------------------------
# 13. Output audit
# ------------------------------------------------------------

expected_output_files_s7d <- c(
  pdf_file_s7d,
  svg_file_s7d,
  tiff_file_s7d,
  png_file_s7d,
  matrix_source_file_s7d,
  paired_source_file_s7d,
  summary_source_file_s7d,
  condition_source_file_s7d
)

if (!all(file.exists(expected_output_files_s7d))) {
  stop(
    "One or more expected S7D output files were not created."
  )
}

output_info_s7d <- file.info(
  expected_output_files_s7d
)

output_audit_s7d <- data.frame(
  file = normalizePath(
    expected_output_files_s7d,
    mustWork = TRUE
  ),
  size_bytes = output_info_s7d$size,
  modified = output_info_s7d$mtime,
  row.names = NULL,
  stringsAsFactors = FALSE
)

cat("\nS7D output directory:\n")
cat(
  normalizePath(
    figure_dir_s7d,
    mustWork = TRUE
  ),
  "\n"
)

cat("\nS7D output files:\n")
print(
  output_audit_s7d,
  row.names = FALSE
)

cat(
  paste0(
    "\nPASS: the final matrix contains all ",
    "28 state-readout comparisons.\n"
  )
)

cat(
  paste0(
    "PASS: the donor-level plot contains all ",
    "63 paired delta observations.\n"
  )
)

cat(
  paste0(
    "PASS: PDF, SVG, TIFF, PNG and four ",
    "source-data tables were created.\n"
  )
)

cat(
  "====================================================\n"
)