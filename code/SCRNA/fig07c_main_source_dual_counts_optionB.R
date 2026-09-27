#!/usr/bin/env Rscript
# =============================================================================
# fig07c_main_source_dual_counts_optionB.R
#
# Figure 7C — compact dual mini-panel for gene-level source attribution.
#   left  : main enriched state    (logCPM peak state)
#   right : main contributor state (bulk-like contribution state)
#
# Reads only:
#   per_gene_cell_source_summary.csv
#
# No Seurat / RDS dependency.
# Uses Option B story-focused source-state palette.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- config ---------------------------------------------------------------
RUN_DIR <- Sys.getenv(
  "FIG07C_RUN_DIR",
  unset = "/home/sunshine/predicate/singlecell/outputs/GSE216009_32gene_analysis/02_gene_cell_distribution/20260803_gene_cell_distribution_5"
)

CSV_IN  <- file.path(RUN_DIR, "tables", "per_gene_cell_source_summary.csv")
OUT_DIR <- file.path(RUN_DIR, "fig07c_main_source_dual")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PNG  <- file.path(OUT_DIR, "fig07c_main_enriched_contributor_counts_optionB.png")
OUT_PDF  <- file.path(OUT_DIR, "fig07c_main_enriched_contributor_counts_optionB.pdf")
OUT_TIFF <- file.path(OUT_DIR, "fig07c_main_enriched_contributor_counts_optionB.tiff")
OUT_DATA <- file.path(OUT_DIR, "fig07c_main_enriched_contributor_counts_data.csv")

stopifnot(file.exists(CSV_IN))

# ---- helpers --------------------------------------------------------------
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

module_label <- function(x) {
  y <- as.character(x)
  y[y %in% c("Cell cycle / proliferation", "Cell cycle")] <- "Cell cycle"
  y[y %in% c("Neutrophil degranulation", "Neutrophil")] <- "Neutrophil"
  y[y %in% c("Inflammatory / Down", "Inflammatory")] <- "Inflammatory"
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
    "S100A8-9_hi_neutrophils" = "S100A8/A9-hi neut",
    "Mature_neutrophils" = "Mature neut",
    "Degranulating_neutrophils" = "Degran neut",
    "Apoptosing_neutrophils" = "Apoptosing neut",
    "Eosinophils" = "Eosinophils",
    "Mast_cells/eosiniophils" = "Mast/eos",
    "Classical_monocytes" = "Classical mono",
    "Non-classical_monocytes" = "Non-classical mono",
    "cDCs" = "cDCs",
    "pDCs" = "pDCs",
    "Naive_CD4_T_cells" = "Naive CD4 T",
    "Memory_CD4_T_cells" = "Memory CD4 T",
    "Naive_CD8_T_cells" = "Naive CD8 T",
    "CD8_T_cells" = "CD8 T",
    "Cycling_TNK" = "Cycling T/NK",
    "NK cells" = "NK cells",
    "B_cells" = "B cells",
    "Plasmablasts" = "Plasmablasts",
    "Platelets" = "Platelets"
  )
  ifelse(z %in% names(map), unname(map[z]), z)
}

# Automatic text-color choice using WCAG relative luminance and contrast ratio.
hex_to_rgb01 <- function(hex) {
  h <- gsub("#", "", hex)
  c(
    r = strtoi(substr(h, 1, 2), 16L) / 255,
    g = strtoi(substr(h, 3, 4), 16L) / 255,
    b = strtoi(substr(h, 5, 6), 16L) / 255
  )
}

srgb_to_linear <- function(x) {
  ifelse(x <= 0.03928, x / 12.92, ((x + 0.055) / 1.055)^2.4)
}

relative_luminance <- function(hex) {
  rgb <- hex_to_rgb01(hex)
  lin <- srgb_to_linear(rgb)
  0.2126 * lin[["r"]] + 0.7152 * lin[["g"]] + 0.0722 * lin[["b"]]
}

contrast_ratio <- function(hex1, hex2) {
  l1 <- relative_luminance(hex1)
  l2 <- relative_luminance(hex2)
  lighter <- max(l1, l2)
  darker  <- min(l1, l2)
  (lighter + 0.05) / (darker + 0.05)
}

# Option B: story-focused source-state palette.
state_pal <- c(
  "Cycling_neutrophil_progenitors" = "#00B8FF",
  "MPO+_immature_neutrophils_or_progenitors" = "#00A5E3",
  "PADI4+_immature_neutrophils" = "#0086D1",
  "Cycling_TNK" = "#FF6F69",
  "S100A8-9_hi_neutrophils" = "#00549E",
  "Mature_neutrophils" = "#0072B2",
  "Apoptosing_neutrophils" = "#003B73",
  "NK cells" = "#A63603",
  "Classical_monocytes" = "#2FBF71",
  "Plasmablasts" = "#8E007E",
  "Mast_cells/eosiniophils" = "#B83B7D"
)

state_order <- names(state_pal)

TYPE_ENR <- "Highest-expression reference state\n(logCPM peak)"
TYPE_CON <- "Highest-contributor reference state\n(donor-normalized UMI contribution)"
TYPE_LEVELS <- c(TYPE_ENR, TYPE_CON)

module_order <- c(
  "Cell cycle",
  "Neutrophil",
  "Inflammatory",
  "Other"
)

module_labs <- c(
  "Cell cycle" = "Cell-cycle/\nproliferation\n(n=11)",
  "Neutrophil" = "Neutrophil\ndegranulation\n(n=10)",
  "Inflammatory" = "Inflammatory/\ndownregulated\n(n=3)",
  "Other" = "Other\n(n=5)"
)

# ---- read ----------------------------------------------------------------
d0 <- read.csv(CSV_IN, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c(
  "signature_gene",
  "module",
  "main_enriched_state",
  "main_contributor_state",
  "primary_analysis_flag"
)
stopifnot(all(required_cols %in% names(d0)))
stopifnot(nrow(d0) == 32L)

d0$primary_analysis_flag <- as_bool(d0$primary_analysis_flag)
d0$module_f <- module_label(d0$module)

# 7C shows only stable / primary-analysis genes.
d <- d0[d0$primary_analysis_flag, , drop = FALSE]

# ---- fail-closed checks ---------------------------------------------------
stopifnot(nrow(d) == 29L)

low_genes <- setdiff(d0$signature_gene, d$signature_gene)
stopifnot(setequal(low_genes, c("HIST1H2BM", "HIST1H3B", "RHAG")))

module_counts <- table(factor(d$module_f, levels = module_order))
stopifnot(
  module_counts[["Cell cycle"]] == 11L,
  module_counts[["Neutrophil"]] == 10L,
  module_counts[["Inflammatory"]] == 3L,
  module_counts[["Other"]] == 5L
)

# Ensure all source states are explicitly colored.
all_states <- unique(c(d$main_enriched_state, d$main_contributor_state))
stopifnot(all(all_states %in% names(state_pal)))

# Headline enriched-state checks used in §3.7 text.
stopifnot(
  sum(d$module_f == "Cell cycle" &
        d$main_enriched_state == "Cycling_neutrophil_progenitors") == 9L,
  sum(d$module_f == "Cell cycle" &
        d$main_enriched_state == "Cycling_TNK") == 2L,
  sum(d$module_f == "Neutrophil" &
        d$main_enriched_state == "MPO+_immature_neutrophils_or_progenitors") == 7L,
  sum(d$module_f == "Neutrophil" &
        d$main_enriched_state == "PADI4+_immature_neutrophils") == 3L
)

# Contributor-state checks.
stopifnot(
  sum(d$module_f == "Cell cycle" &
        d$main_contributor_state == "MPO+_immature_neutrophils_or_progenitors") == 9L,
  sum(d$module_f == "Neutrophil" &
        d$main_contributor_state == "MPO+_immature_neutrophils_or_progenitors") == 8L,
  sum(d$module_f == "Neutrophil" &
        d$main_contributor_state == "PADI4+_immature_neutrophils") == 1L,
  sum(d$module_f == "Neutrophil" &
        d$main_contributor_state == "S100A8-9_hi_neutrophils") == 1L
)

# ---- build compact count table -------------------------------------------
enr <- data.frame(
  source_type = factor(TYPE_ENR, levels = TYPE_LEVELS),
  module_f = d$module_f,
  source_state = d$main_enriched_state,
  stringsAsFactors = FALSE
)

con <- data.frame(
  source_type = factor(TYPE_CON, levels = TYPE_LEVELS),
  module_f = d$module_f,
  source_state = d$main_contributor_state,
  stringsAsFactors = FALSE
)

plot_d <- rbind(enr, con)
plot_d$source_type <- factor(plot_d$source_type, levels = TYPE_LEVELS)
plot_d$module_f <- factor(plot_d$module_f, levels = module_order)
plot_d$source_state <- factor(plot_d$source_state, levels = state_order)

plot_d <- as.data.frame(
  table(
    source_type = plot_d$source_type,
    module_f = plot_d$module_f,
    source_state = plot_d$source_state
  ),
  stringsAsFactors = FALSE
)
names(plot_d)[names(plot_d) == "Freq"] <- "n"
plot_d <- plot_d[plot_d$n > 0, , drop = FALSE]

plot_d$source_type <- factor(plot_d$source_type, levels = TYPE_LEVELS)
plot_d$module_y <- factor(
  as.character(plot_d$module_f),
  levels = rev(module_order)
)
plot_d$state_short <- short_state(plot_d$source_state)
# Explicitly lock the states used in the plot to the named palette.
used_states <- state_order[
  state_order %in% unique(as.character(plot_d$source_state))
]

plot_d$source_state <- factor(
  as.character(plot_d$source_state),
  levels = used_states
)

state_pal_used <- state_pal[used_states]

stopifnot(
  identical(names(state_pal_used), used_states),
  all(as.character(plot_d$source_state) %in% names(state_pal_used))
)

legend_audit <- data.frame(
  source_state = used_states,
  display_label = short_state(used_states),
  fill_color = unname(state_pal_used),
  stringsAsFactors = FALSE
)

write.csv(
  legend_audit,
  file.path(
    OUT_DIR,
    "fig07c_legend_state_color_audit.csv"
  ),
  row.names = FALSE
)
plot_d$label <- as.character(plot_d$n)

# Dynamic label color from fill contrast.
fill_hex <- state_pal[as.character(plot_d$source_state)]

contrast_with_black <- vapply(
  fill_hex,
  contrast_ratio,
  numeric(1),
  hex2 = "#000000"
)

contrast_with_white <- vapply(
  fill_hex,
  contrast_ratio,
  numeric(1),
  hex2 = "#FFFFFF"
)

plot_d$label_col <- ifelse(
  contrast_with_black >= contrast_with_white,
  "black",
  "white"
)

totals <- aggregate(
  n ~ source_type + module_y,
  data = plot_d,
  FUN = sum
)
totals$source_type <- factor(totals$source_type, levels = TYPE_LEVELS)

write.csv(plot_d, OUT_DATA, row.names = FALSE)

# ---- plot ----------------------------------------------------------------
p <- ggplot(
  plot_d,
  aes(
    x = n,
    y = module_y,
    fill = source_state,
    group = source_state
  )
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = label,
      color = label_col,
      group = source_state
    ),
    position = position_stack(vjust = 0.5),
    size = 3.3,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ source_type,
    nrow = 1
  ) +
  scale_fill_manual(
    values = state_pal_used,
    limits = used_states,
    breaks = used_states,
    labels = short_state(used_states),
    name = "Reference state",
    drop = FALSE
  ) +
  scale_color_identity() +
  scale_y_discrete(
    labels = module_labs
  ) +
  scale_x_continuous(
    breaks = seq(0, 12, by = 2),
    expand = expansion(mult = c(0, 0.13))
  ) +
  coord_cartesian(
    xlim = c(0, 12.4),
    clip = "off"
  ) +
  labs(
    x = "Genes",
    y = NULL,
    title = "Reference-state mapping of 29 stably detected signature genes"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0
    ),
    axis.text.x = element_text(
      size = 9,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.line = element_line(
      color = "grey35",
      linewidth = 0.4
    ),
    axis.ticks = element_line(
      color = "grey35",
      linewidth = 0.3
    ),
    panel.spacing.x = unit(10, "mm"),
    legend.position = "bottom",
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    legend.text = element_text(size = 8),
    legend.key.size = unit(4, "mm"),
    plot.margin = margin(8, 18, 8, 8)
  ) +
  guides(
    fill = guide_legend(
      nrow = 3,
      byrow = TRUE,
      title.position = "top"
    )
  )

# ---- save ----------------------------------------------------------------
ggsave(OUT_PNG, p, width = 9.2, height = 5.2, dpi = 600, bg = "white")
ggsave(OUT_PDF, p, width = 9.2, height = 5.2, bg = "white")

ggsave(
  OUT_TIFF,
  p,
  width = 9.2,
  height = 5.2,
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

cat("Wrote:\n")
cat("  ", OUT_PNG, "\n", sep = "")
cat("  ", OUT_PDF, "\n", sep = "")
cat("  ", OUT_TIFF, "\n", sep = "")
cat("  ", OUT_DATA, "\n", sep = "")