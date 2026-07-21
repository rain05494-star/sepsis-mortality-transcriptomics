# Run this script in a clean R session.
# ============================================================
# GSE272769 + GSE95233 联合 Meta 分析
# 7 个子文件夹，每节独立 checkpoint
# ============================================================
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(VennDiagram)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(metafor)
  library(grid)
  library(dplyr)
  library(enrichplot)
  library(pheatmap)
  library(enrichR)
  library(fgsea)
  library(msigdbr)
  library(xCell)
  library(patchwork)
  library(httr)
  library(jsonlite)
  library(reshape2)
})

# ============================================================
# 0. 全局参数
# ============================================================
LFC_CUT <- 0.5
FDR_CUT <- 0.05

# Nominal within-cohort support threshold.
# Used for descriptive annotation only, not signature selection.
INDIVIDUAL_P_CUT <- 0.05

I2_BREAKS <- c(-Inf, 25, 50, 75, Inf)

I2_LABELS <- c(
  "<25%",
  "25-<50%",
  "50-<75%",
  ">=75%"
)
# ---- Project directories ----

command_args <- commandArgs(trailingOnly = TRUE)

if (length(command_args) >= 1L && nzchar(command_args[1])) {
  PROJECT_DIR <- normalizePath(
    command_args[1],
    winslash = "/",
    mustWork = TRUE
  )
} else {
  PROJECT_DIR <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}
DATA_DIR <- file.path(PROJECT_DIR, "data")
ROOT_DIR <- file.path(PROJECT_DIR, "results_meta")
if (!dir.exists(DATA_DIR)) {
  stop(sprintf(
    paste0(
      "Input directory was not found: %s\n",
      "Run the script from the project root or provide the ",
      "project directory as the first command-line argument."
    ),
    DATA_DIR
  ))
}

dir.create(
  ROOT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

S1_DIR <- file.path(ROOT_DIR, "01_data_prep")
S2_DIR <- file.path(ROOT_DIR, "02_meta_analysis")
S3_DIR <- file.path(ROOT_DIR, "03_discovery_figs")
S4_DIR <- file.path(ROOT_DIR, "04_GSE65682")
S5_DIR <- file.path(ROOT_DIR, "05_triple_cohort")
S6_DIR <- file.path(ROOT_DIR, "06_drug")
S7_DIR <- file.path(ROOT_DIR, "07_supp")

output_directories <- c(
  S1_DIR,
  S2_DIR,
  S3_DIR,
  S4_DIR,
  S5_DIR,
  S6_DIR,
  S7_DIR
)

for (output_directory in output_directories) {
  dir.create(
    output_directory,
    showWarnings = FALSE,
    recursive = TRUE
  )
}

# ---- Local GEO input files ----

gse272769_file <- file.path(
  DATA_DIR,
  "GSE272769_series_matrix.txt.gz"
)

gse95233_file <- file.path(
  DATA_DIR,
  "GSE95233_series_matrix.txt.gz"
)

gse65682_file <- file.path(
  DATA_DIR,
  "GSE65682_series_matrix.txt.gz"
)

gpl17692_file <- file.path(
  DATA_DIR,
  "GPL17692.soft.gz"
)

gpl570_file <- file.path(
  DATA_DIR,
  "GPL570.soft.gz"
)

gpl13667_file <- file.path(
  DATA_DIR,
  "GPL13667.soft.gz"
)

required_input_files <- c(
  gse272769_file,
  gse95233_file,
  gse65682_file,
  gpl17692_file,
  gpl570_file,
  gpl13667_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(sprintf(
    "Missing required input files:\n%s",
    paste(
      paste0("  - ", missing_input_files),
      collapse = "\n"
    )
  ))
}

cat(sprintf("Project directory: %s\n", PROJECT_DIR))
cat(sprintf("Input directory:   %s\n", DATA_DIR))
cat(sprintf("Output directory:  %s\n", ROOT_DIR))

# ---- 颜色 ----
cols_v  <- c("NS"="grey80", "Up"="#E15759", "Down"="#4E79A7")
cols_sc <- c("Not significant"="grey85", "Up"="#E15759", "Down"="#4E79A7")
cols_i2 <- setNames(c("#4E79A7","#B07AA1","#E15759","#B6992D","grey70"), c(I2_LABELS,"NA"))
support_cols <- c(
  "Both cohorts" = "#4E79A7",
  "GSE272769 only" = "#E15759",
  "GSE95233 only" = "#B07AA1",
  "Neither cohort" = "grey50"
)
logFC_COL_95233 <- "logFC_95233 (D01)"
pval_COL_95233  <- "pval_95233 (D01)"
fdr_COL_95233   <- "fdr_95233 (D01)"
t_COL_95233     <- "t_95233 (D01)"
ave_COL_95233   <- "AveExpr_95233 (D01)"

# ---- 工具 ----
save_plot <- function(p, name, d = S3_DIR, w = 7.5, h = 7, res = 300) {
  pdf(file.path(d, paste0(name,".pdf")), w, h); print(p); dev.off()
  png(file.path(d, paste0(name,".png")), w, h, units="in", res=res); print(p); dev.off()
  cat(sprintf("  saved: %s/%s\n", basename(d), name))
}

check_no_na <- function(x, label) {
  if (anyNA(x)) stop(sprintf("%s: %d NA values", label, sum(is.na(x))))
}

check_group <- function(raw, name_map, label) {
  cat(sprintf("\n[%s] raw group table:\n", label))
  raw <- trimws(as.character(raw))
  raw[raw %in% c("","NA","N/A","na","n/a")] <- NA_character_
  print(table(raw, useNA = "ifany"))
  if (anyNA(raw)) stop(sprintf("%s: raw group contains NA", label))
  mapped <- rep(NA_character_, length(raw))
  for (tgt in names(name_map)) {
    vals <- trimws(as.character(name_map[[tgt]]))
    mapped[raw %in% vals] <- tgt
  }
  unrec <- is.na(mapped)
  if (any(unrec)) { cat(sprintf("  Unrecognized (%d):\n",sum(unrec))); print(unique(raw[unrec])); stop(sprintf("%s: unrecognized",label)) }
  f <- factor(mapped, levels=c("Survivor","NonSurvivor"))
  cat(sprintf("  Survivor=%d | NonSurvivor=%d\n",sum(f=="Survivor"),sum(f=="NonSurvivor")))
  check_no_na(f, label); f
}

run_limma <- function(emat, grp, name) {
  stopifnot(is.factor(grp), !anyNA(grp), ncol(emat)==length(grp))
  tbl <- table(grp)
  if (any(tbl < 2)) stop(sprintf("%s: need >=2 per group, got: %s", name, paste(tbl,collapse=",")))
  design <- model.matrix(~0+grp); colnames(design) <- levels(grp)
  contrast <- makeContrasts(NonSurvivor-Survivor, levels=design)
  fit <- lmFit(emat,design); fit2 <- contrasts.fit(fit,contrast); fit2 <- eBayes(fit2,trend=TRUE)
  res <- topTable(fit2,coef=1,number=Inf,sort.by="none"); res$gene <- rownames(res)
  cat(sprintf("%s: %d genes | Up=%d Down=%d\n", name, nrow(res),
              sum(res$logFC>LFC_CUT&res$adj.P.Val<FDR_CUT,na.rm=TRUE),
              sum(res$logFC<(-LFC_CUT)&res$adj.P.Val<FDR_CUT,na.rm=TRUE)))
  res
}

invalid_gene <- function(x) {
  is.na(x) | x=="" | x=="---" | x=="NA" | x=="N/A" | nchar(x)>30
}

parse_gene_17692 <- function(x) {
  sapply(strsplit(trimws(as.character(x)), " /// "), function(block) {
    if (length(block)==0) return(NA_character_)
    fields <- strsplit(trimws(block[1])," // ")[[1]]
    sym <- if (length(fields)>=2) trimws(fields[2]) else NA_character_
    if (invalid_gene(sym)) NA_character_ else sym
  }, USE.NAMES=FALSE)
}

parse_gene_570 <- function(x) {
  sapply(strsplit(trimws(as.character(x))," /// "), function(block) {
    sym <- trimws(block[1])
    if (invalid_gene(sym)) NA_character_ else sym
  }, USE.NAMES=FALSE)
}

parse_gene_symbol <- function(x) {
  sapply(trimws(as.character(x)), function(z) {
    if (is.na(z) || !nzchar(z)) return(NA_character_)
    parts <- unlist(strsplit(z, "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*"))
    parts <- trimws(parts); parts <- parts[!invalid_gene(parts)]
    if (length(parts)==0) NA_character_ else parts[1]
  }, USE.NAMES=FALSE)
}

probes_to_genes <- function(emat, pmap) {
  pmap$probe <- trimws(as.character(pmap$probe))
  pmap$gene  <- trimws(as.character(pmap$gene))
  pmap <- pmap[!invalid_gene(pmap$gene), ]
  common <- intersect(rownames(emat), pmap$probe)
  if (length(common)==0) stop("No probes matched")
  cat(sprintf("  probes: %d -> matched: %d\n", nrow(emat), length(common)))
  esub <- emat[common,,drop=FALSE]; msub <- pmap[match(common,pmap$probe),]
  gs <- split(seq_len(nrow(esub)), msub$gene)
  gmat <- matrix(NA_real_, length(gs), ncol(esub))
  rownames(gmat) <- names(gs); colnames(gmat) <- colnames(esub)
  for (i in seq_along(gs)) {
    idx <- gs[[i]]
    vars <- apply(esub[idx,,drop=FALSE],1,var,na.rm=TRUE)
    gmat[i,] <- if (all(is.na(vars))) esub[idx[1],] else esub[idx[which.max(vars)],]
  }
  if (anyNA(gmat)) {
    na_by_gene <- rowSums(
      is.na(gmat)
    )
    
    affected_genes <- names(
      na_by_gene[na_by_gene > 0]
    )
    
    na_by_sample <- colSums(
      is.na(gmat)
    )
    
    affected_samples <- names(
      na_by_sample[na_by_sample > 0]
    )
    
    stop(
      paste0(
        "Gene-level expression matrix contains ",
        sum(is.na(gmat)),
        " missing values across ",
        length(affected_genes),
        " genes and ",
        length(affected_samples),
        " samples. ",
        "Affected genes include: ",
        paste(
          head(affected_genes, 10),
          collapse = ", "
        ),
        if (
          length(affected_genes) > 10
        ) {
          ", ..."
        } else {
          ""
        },
        ". Missing values must be investigated before analysis."
      )
    )
  }
  
  cat(sprintf(
    "  genes: %d\n",
    nrow(gmat)
  ))
  
  check_no_na(
    gmat,
    "gene matrix"
  )
  
  gmat
}

get_direction_vec <- function(lfc1, lfc2, lfc_cut=LFC_CUT) {
  res <- rep("NS", length(lfc1))
  ok  <- !is.na(lfc1) & !is.na(lfc2)
  res[ok & lfc1 >  lfc_cut & lfc2 >  lfc_cut] <- "Up"
  res[ok & lfc1 < -lfc_cut & lfc2 < -lfc_cut] <- "Down"
  res
}

run_re_meta_one <- function(lfc1, se1, lfc2, se2) {
  yi <- c(lfc1,lfc2); sei <- c(se1,se2); ok <- is.finite(yi)&is.finite(sei)&sei>0
  if (sum(ok)<2) return(list(converged=FALSE,re_p=NA_real_,tau2=NA_real_,I2=NA_real_,re_estimate=NA_real_,re_se=NA_real_))
  fit <- tryCatch(suppressWarnings(rma(yi=yi[ok],sei=sei[ok],method="REML")),error=function(e)NULL)
  if (is.null(fit)) return(list(converged=FALSE,re_p=NA_real_,tau2=NA_real_,I2=NA_real_,re_estimate=NA_real_,re_se=NA_real_))
  list(converged=TRUE,re_p=as.numeric(fit$pval),tau2=as.numeric(fit$tau2),I2=as.numeric(fit$I2),re_estimate=as.numeric(fit$b[1]),re_se=as.numeric(fit$se[1]))
}

write_lines_safe <- function(x, path) {
  if (length(x)>0) writeLines(x, path) else { cat(sprintf("  Empty -> skip %s\n",basename(path))); file.create(path) }
}

run_pca <- function(emat, genes_use, group_vec, title) {
  idx <- intersect(genes_use, rownames(emat))
  if (length(idx)<100) stop(sprintf("%s: too few genes (%d)",title,length(idx)))
  pca <- prcomp(t(emat[idx,,drop=FALSE]), scale.=TRUE)
  cat(sprintf("  %s: %d genes, %d samples | PC1=%.1f%% PC2=%.1f%%\n",title,length(idx),ncol(emat),
              summary(pca)$importance[2,1]*100, summary(pca)$importance[2,2]*100))
  pca
}

draw_heatmap_triple <- function(mat, title, fname, d, w=8, h=7) {
  max_abs <- max(abs(mat),na.rm=TRUE)
  pdf(file.path(d,paste0(fname,".pdf")), w, h)
  pheatmap(mat, color=colorRampPalette(c("#4E79A7","white","#E15759"))(100),
           breaks=seq(-max_abs,max_abs,length.out=101),
           cluster_rows=FALSE,cluster_cols=FALSE,display_numbers=TRUE,number_format="%.2f",
           fontsize_number=7,fontsize_row=9,fontsize_col=11,angle_col=0,border_color=NA,main=title)
  dev.off()
  png(file.path(d,paste0(fname,".png")),w,h,units="in",res=300)
  pheatmap(mat, color=colorRampPalette(c("#4E79A7","white","#E15759"))(100),
           breaks=seq(-max_abs,max_abs,length.out=101),
           cluster_rows=FALSE,cluster_cols=FALSE,display_numbers=TRUE,number_format="%.2f",
           fontsize_number=7,fontsize_row=9,fontsize_col=11,angle_col=0,border_color=NA,main=title)
  dev.off()
  cat(sprintf("  %s saved.\n", fname))
}

enrichr_retry_v2 <- function(genes, dbs, max_retry=5, base_wait=3) {
  for (i in seq_len(max_retry)) {
    cat(sprintf("  Attempt %d/%d ... ", i, max_retry))
    probe <- tryCatch(httr::GET("https://maayanlab.cloud/Enrichr/", httr::timeout(8)), error=function(e)NULL)
    if (is.null(probe) || probe$status_code!=200) { wait<-base_wait*2^(i-1); cat(sprintf("probe failed (waiting %ds)\n",wait)); Sys.sleep(wait); next }
    res <- tryCatch(suppressMessages(enrichR::enrichr(genes, dbs)), error=function(e){cat(sprintf("error: %s ",e$message));NULL})
    has_data <- !is.null(res)&&length(res)>0&&any(sapply(res,function(x) is.data.frame(x)&&nrow(x)>0))
    if (has_data) { cat("SUCCESS\n"); return(res) }
    wait<-base_wait*2^(i-1); cat(sprintf("empty (waiting %ds)\n",wait)); Sys.sleep(wait)
  }
  cat("  ALL ATTEMPTS FAILED\n"); NULL
}

safe_enrichr_df <- function(res, db) {
  if (is.null(res)) return(data.frame())
  if (!db %in% names(res)) return(data.frame())
  df <- res[[db]]
  if (is.null(df)||!is.data.frame(df)||nrow(df)==0) return(data.frame())
  df
}

extract_tf_symbol <- function(term_vec) {
  sapply(term_vec, function(x) {
    x <- trimws(as.character(x)); if (is.na(x)||!nzchar(x)) return(NA_character_)
    parts <- strsplit(x,"[\\s_]+",perl=TRUE)[[1]]; parts<-parts[parts!=""]
    if (length(parts)==0) return(NA_character_)
    for (p in parts) if (grepl("^[A-Z][A-Za-z0-9._-]{1,14}$",p)) return(p)
    parts[1]
  }, USE.NAMES=FALSE)
}

# Extract compound names from Enrichr perturbation terms.
# The original Term column is retained in all exported result tables.
extract_drug <- function(x) {
  x <- as.character(x)
  x <- trimws(gsub("\\s+", " ", x))
  
  species_regex <- paste(
    c(
      "homo sapiens",
      "mus musculus",
      "rattus norvegicus",
      "sus scrofa",
      "danio rerio",
      "bos taurus",
      "macaca mulatta",
      "drosophila melanogaster",
      "caenorhabditis elegans",
      "saccharomyces cerevisiae"
    ),
    collapse = "|"
  )
  
  # Remove species and all following GEO/platform/direction descriptors.
  x <- sub(
    paste0("\\s+(?:", species_regex, ")\\b.*$"),
    "",
    x,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  # Fallback for terms lacking an explicit species descriptor.
  x <- sub(
    "\\s+(?:gpl\\d+|gse\\d+|gds\\d+|chdir\\s+(?:up|down))\\b.*$",
    "",
    x,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  x <- trimws(gsub("\\s+", " ", x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

# Conservative compound key used consistently across Enrichr and DGIdb.
# Punctuation is retained to avoid unintentionally merging distinct database
# records, salts, formulations, or source-specific identifiers.
norm_drug <- function(x) {
  x <- iconv(
    as.character(x),
    from = "",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  x <- tolower(x)
  x <- trimws(gsub("\\s+", " ", x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

# Display formatting only; this function does not change the matching key.
format_drug_label <- function(x) {
  x <- trimws(as.character(x))
  out <- x
  
  valid <- !is.na(x) & nzchar(x)
  source_identifier <- valid & grepl(
    "^[A-Za-z][A-Za-z0-9._-]*:",
    x,
    perl = TRUE
  )
  
  ordinary_name <- valid & !source_identifier
  out[ordinary_name] <- tools::toTitleCase(tolower(x[ordinary_name]))
  out
}

# Retain one Enrichr record per normalized compound key.
# If several perturbation terms map to the same compound, retain the record
# with the highest Combined Score; adjusted P value and Term provide
# deterministic tie-breaking.
top_enrichr_by_drug <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  if (!"drug_name" %in% colnames(df)) {
    stop("top_enrichr_by_drug requires a drug_name column.")
  }
  
  df$Combined.Score <- suppressWarnings(as.numeric(df$Combined.Score))
  df$drug_key <- norm_drug(df$drug_name)
  
  valid <- !is.na(df$drug_key) &
    nzchar(df$drug_key) &
    !is.na(df$Combined.Score)
  
  df <- df[valid, , drop = FALSE]
  if (nrow(df) == 0) return(df)
  
  adjusted_p <- if ("Adjusted.P.value" %in% colnames(df)) {
    suppressWarnings(as.numeric(df$Adjusted.P.value))
  } else {
    rep(Inf, nrow(df))
  }
  
  term_text <- if ("Term" %in% colnames(df)) {
    as.character(df$Term)
  } else {
    rep("", nrow(df))
  }
  
  ord <- order(
    -df$Combined.Score,
    adjusted_p,
    term_text,
    na.last = TRUE
  )
  
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df$drug_key), , drop = FALSE]
  rownames(df) <- NULL
  df
}
# ============================================================
# §1 — 数据准备
# ============================================================
cat("\n========== §1 Data Preparation ==========\n")

gse272769 <- getGEO(filename=gse272769_file,getGPL=FALSE)
gse95233  <- getGEO(filename=gse95233_file, getGPL=FALSE)
pdata_272769 <- pData(gse272769); pdata_95233 <- pData(gse95233)
cat(sprintf("GSE272769: %d probes x %d\n", nrow(exprs(gse272769)), ncol(exprs(gse272769))))
cat(sprintf("GSE95233:  %d probes x %d\n", nrow(exprs(gse95233)),  ncol(exprs(gse95233))))

cat("\n--- Group assignment ---\n")
group_272769 <- check_group(pdata_272769$`mort30:ch1`, name_map=list(NonSurvivor="Yes",Survivor="No"), label="GSE272769")

pdata_95233$tp <- as.character(pdata_95233$`time point:ch1`)
pdata_95233$sv <- trimws(as.character(pdata_95233$`survival:ch1`))
pdata_95233$sv[pdata_95233$sv %in% c("","NA","N/A","na","n/a")] <- NA_character_
cat("\nGSE95233 timepoint x survival:\n"); print(table(pdata_95233$tp,pdata_95233$sv,useNA="ifany"))
keep_idx <- which(!is.na(pdata_95233$tp)&pdata_95233$tp=="D01"&!is.na(pdata_95233$sv)&pdata_95233$sv%in%c("Non Survivor","Survivor"))
stopifnot(length(keep_idx)>0)
gse95233_D01 <- gse95233[,keep_idx]; pdata_95233_D01 <- pData(gse95233_D01)
group_95233 <- check_group(pdata_95233_D01$`survival:ch1`, name_map=list(NonSurvivor="Non Survivor",Survivor="Survivor"), label="GSE95233 D01")

cat("\n--- Probe annotation ---\n")
gpl17692 <- getGEO(filename=gpl17692_file,GSEMatrix=FALSE)
gpl570   <- getGEO(filename=gpl570_file,  GSEMatrix=FALSE)
tab17692 <- Table(gpl17692); tab570 <- Table(gpl570)
stopifnot("gene_assignment"%in%colnames(tab17692),"Gene Symbol"%in%colnames(tab570))
map17692 <- data.frame(probe=as.character(tab17692$ID),gene=parse_gene_17692(tab17692$gene_assignment),stringsAsFactors=FALSE)
map17692 <- map17692[!invalid_gene(map17692$gene),]
map570   <- data.frame(probe=as.character(tab570$ID),gene=parse_gene_570(tab570$`Gene Symbol`),stringsAsFactors=FALSE)
map570   <- map570[!invalid_gene(map570$gene),]

cat("\n--- Probe to gene ---\n")
expr_272769 <- exprs(gse272769); expr_95233 <- exprs(gse95233_D01)
gmat_272769 <- probes_to_genes(expr_272769,map17692)
gmat_95233  <- probes_to_genes(expr_95233, map570)
common_genes <- intersect(rownames(gmat_272769),rownames(gmat_95233))
stopifnot(length(common_genes)>0)
cat(sprintf("Common genes: %d\n",length(common_genes)))
gmat_272769 <- gmat_272769[common_genes,,drop=FALSE]
gmat_95233  <- gmat_95233[common_genes,,drop=FALSE]

cat("\n--- Data space check ---\n")

is_log2_like <- function(x) {
  finite_values <- x[is.finite(x)]
  
  if (length(finite_values) == 0) {
    return(FALSE)
  }
  
  max(finite_values) < 30 &&
    min(finite_values) > -2
}

range_272769 <- range(
  gmat_272769,
  na.rm = TRUE
)

range_95233 <- range(
  gmat_95233,
  na.rm = TRUE
)

is_log2_272769 <- is_log2_like(
  gmat_272769
)

is_log2_95233 <- is_log2_like(
  gmat_95233
)

cat(sprintf(
  paste0(
    "GSE272769 range: %.3f to %.3f | log2-like: %s\n",
    "GSE95233 D01 range: %.3f to %.3f | log2-like: %s\n"
  ),
  range_272769[1],
  range_272769[2],
  is_log2_272769,
  range_95233[1],
  range_95233[2],
  is_log2_95233
))

if (!is_log2_272769 || !is_log2_95233) {
  stop(
    paste0(
      "Both discovery cohorts were expected to contain ",
      "log2-scale processed expression values. ",
      "Inspect and transform each cohort explicitly before continuing."
    )
  )
}

# xCell requires non-negative approximately linear-scale expression.
gmat_272769_linear <- pmax(
  2^gmat_272769 - 1,
  0
)

gmat_95233_linear <- pmax(
  2^gmat_95233 - 1,
  0
)
saveRDS(list(
  gmat_272769=gmat_272769, gmat_95233=gmat_95233,
  gmat_272769_linear=gmat_272769_linear, gmat_95233_linear=gmat_95233_linear,
  group_272769=group_272769, group_95233=group_95233,
  common_genes=common_genes
), file.path(S1_DIR,"checkpoint_s1.rds"))
cat("§1 complete.\n")

# ============================================================
# §2 — Meta-analysis + PCA + Tables S1/S2
# ============================================================
cat("\n========== §2 Meta-analysis ==========\n")

res_272769 <- run_limma(gmat_272769,group_272769,"GSE272769")
res_95233  <- run_limma(gmat_95233, group_95233,"GSE95233 D01")

df1 <- res_272769[,c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")]
colnames(df1) <- c("gene","logFC_272769","AveExpr_272769","t_272769","pval_272769","fdr_272769")
df2 <- res_95233[,c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")]
colnames(df2) <- c("gene",logFC_COL_95233,ave_COL_95233,t_COL_95233,pval_COL_95233,fdr_COL_95233)
meta_df <- merge(df1,df2,by="gene",all=FALSE)

meta_df$se_272769 <- pmax(abs(meta_df$logFC_272769/meta_df$t_272769),1e-8)
meta_df$se_95233  <- pmax(abs(meta_df[[logFC_COL_95233]]/meta_df[[t_COL_95233]]),1e-8)
meta_df$direction <- get_direction_vec(meta_df$logFC_272769, meta_df[[logFC_COL_95233]])
meta_df$sign_agree <- sign(meta_df$logFC_272769)==sign(meta_df[[logFC_COL_95233]])
dir_rate <- mean(meta_df$sign_agree,na.rm=TRUE)*100

# Fisher
fisher_p_272769 <- pmax(meta_df$pval_272769,.Machine$double.xmin)
fisher_p_95233 <- pmax(meta_df[[pval_COL_95233]],.Machine$double.xmin)
meta_df$fisher_stat <- -2*(log(fisher_p_272769)+log(fisher_p_95233))
meta_df$fisher_pval <- pchisq(meta_df$fisher_stat,df=4,lower.tail=FALSE)
meta_df$fisher_fdr  <- p.adjust(meta_df$fisher_pval,"BH")

# RE meta
cat("Running RE meta...\n")
re_list <- lapply(seq_len(nrow(meta_df)),function(i) {
  run_re_meta_one(meta_df$logFC_272769[i],meta_df$se_272769[i],meta_df[[logFC_COL_95233]][i],meta_df$se_95233[i])
})
re_map <- c("converged"="logical","re_p"="numeric","tau2"="numeric","I2"="numeric","re_estimate"="numeric","re_se"="numeric")
re_cols <- names(re_map)
re_df <- as.data.frame(setNames(lapply(re_cols,function(col) {
  vapply(re_list,function(x) if(is.list(x)) x[[col]] else NA, FUN.VALUE=vector(re_map[col],1L))
}),re_cols), stringsAsFactors=FALSE)
meta_df <- cbind(meta_df,re_df); meta_df$re_fdr <- p.adjust(meta_df$re_p,"BH")

fisher_out <- meta_df[meta_df$direction%in%c("Up","Down")&meta_df$fisher_fdr<FDR_CUT,,drop=FALSE]
if (nrow(fisher_out)>0) {
  fisher_out$mean_abs_logFC <- (abs(fisher_out$logFC_272769)+abs(fisher_out[[logFC_COL_95233]]))/2
  fisher_out <- fisher_out[order(fisher_out$fisher_fdr,-fisher_out$mean_abs_logFC,na.last=TRUE),]
}
n_up<-sum(fisher_out$direction=="Up",na.rm=TRUE); n_down<-sum(fisher_out$direction=="Down",na.rm=TRUE)
cat(sprintf("Fisher: Up=%d Down=%d | Dir=%.1f%%\n",n_up,n_down,dir_rate))

# Individual-cohort nominal support
fisher_out$individual_support_272769 <-
  !is.na(fisher_out$pval_272769) &
  abs(fisher_out$logFC_272769) > LFC_CUT &
  fisher_out$pval_272769 < INDIVIDUAL_P_CUT

fisher_out$individual_support_95233 <-
  !is.na(fisher_out[[pval_COL_95233]]) &
  abs(fisher_out[[logFC_COL_95233]]) > LFC_CUT &
  fisher_out[[pval_COL_95233]] < INDIVIDUAL_P_CUT

fisher_out$individual_cohort_support <- ifelse(
  fisher_out$individual_support_272769 &
    fisher_out$individual_support_95233,
  "Both cohorts",
  ifelse(
    fisher_out$individual_support_272769,
    "GSE272769 only",
    ifelse(
      fisher_out$individual_support_95233,
      "GSE95233 only",
      "Neither cohort"
    )
  )
)

fisher_out$individual_cohort_support <- factor(
  fisher_out$individual_cohort_support,
  levels = names(support_cols)
)

n_both <- sum(
  fisher_out$individual_cohort_support == "Both cohorts",
  na.rm = TRUE
)

cat(sprintf(
  paste0(
    "Individual-cohort nominal support: %d/%d genes met ",
    "|log2FC| > %.2f and nominal P < %.2f in both cohorts\n"
  ),
  n_both,
  nrow(fisher_out),
  LFC_CUT,
  INDIVIDUAL_P_CUT
))

# ---- Functional annotation labels for reporting and visualization ----
# Defined after Fisher-based gene selection; these labels are not selection criteria.
MODULE_NEUTRO      <- c("OLFM4","DEFA4","CEACAM8","CEACAM6","ELANE","MPO","CTSG","RNASE3","TCN1","MS4A3")
MODULE_CELL_CYCLE   <- c("TOP2A","MKI67","BUB1","CCNA2","TPX2","KIF11","ASPM","CKS2","ANLN","TYMS","RRM2")
MODULE_DOWN_DEFENSE <- c("IL1B","CX3CR1","SECTM1")
MODULE_OTHER        <- c("CD24","CA2","SERPINB10","CHI3L1","ABCA13","RHAG","HIST1H3B","HIST1H2BM")
module_map <- c(
  setNames(rep("Neutrophil degranulation", length(MODULE_NEUTRO)),       MODULE_NEUTRO),
  setNames(rep("Cell cycle / proliferation", length(MODULE_CELL_CYCLE)),  MODULE_CELL_CYCLE),
  setNames(rep("Inflammatory / Down",        length(MODULE_DOWN_DEFENSE)), MODULE_DOWN_DEFENSE),
  setNames(rep("Other",                      length(MODULE_OTHER)),       MODULE_OTHER)
)
MODULE_LEVELS <- c("Neutrophil degranulation","Cell cycle / proliferation","Inflammatory / Down","Other","Unclassified")
mod_cols <- c("Neutrophil degranulation"="#E15759","Cell cycle / proliferation"="#B07AA1",
              "Inflammatory / Down"="#4E79A7","Other"="#B6992D","Unclassified"="grey50")

# ---- Table S1: Genome-wide direction cross-classification ----
dir1 <- ifelse(meta_df$logFC_272769>LFC_CUT,"Up",ifelse(meta_df$logFC_272769<(-LFC_CUT),"Down","NS"))
dir2 <- ifelse(meta_df[[logFC_COL_95233]]>LFC_CUT,"Up",ifelse(meta_df[[logFC_COL_95233]]<(-LFC_CUT),"Down","NS"))
tab_s1 <- table(dir1, dir2)
tab_s1 <- tab_s1[c("Up","NS","Down"), c("Up","NS","Down")]
cat("\nTable S1: Genome-wide direction cross-classification\n")
print(tab_s1)
cat(sprintf("Direction agreement: %.1f%%\n",dir_rate))
write.csv(as.data.frame.matrix(tab_s1), file.path(S2_DIR,"Table_S1_direction_cross.csv"))

# ---- Table S2: Individual-cohort nominal support ----
tab_s2 <- data.frame(
  gene = fisher_out$gene,
  logFC_272769 = fisher_out$logFC_272769,
  nominal_P_272769 = fisher_out$pval_272769,
  support_272769 = fisher_out$individual_support_272769,
  logFC_95233_D01 = fisher_out[[logFC_COL_95233]],
  nominal_P_95233_D01 = fisher_out[[pval_COL_95233]],
  support_95233_D01 = fisher_out$individual_support_95233,
  individual_cohort_support =
    as.character(fisher_out$individual_cohort_support),
  stringsAsFactors = FALSE
)

tab_s2 <- tab_s2[
  order(
    match(tab_s2$individual_cohort_support, names(support_cols)),
    tab_s2$gene
  ),
  ,
  drop = FALSE
]

rownames(tab_s2) <- NULL

write.csv(
  tab_s2,
  file.path(S2_DIR, "Table_S2_individual_cohort_support.csv"),
  row.names = FALSE
)

cat(sprintf(
  "Table S2: individual-cohort support reported for %d genes\n",
  nrow(tab_s2)
))
# ============================================================
# STRING web and Cytoscape input files
# ============================================================

# Plain gene list submitted to the STRING web interface
string_input_genes <- unique(as.character(fisher_out$gene))
string_input_genes <- string_input_genes[
  !is.na(string_input_genes) &
    nzchar(string_input_genes)
]

stopifnot(length(string_input_genes) == 32)

writeLines(
  string_input_genes,
  file.path(
    S3_DIR,
    "STRING_web_input_32gene_signature.txt"
  )
)

# Node attributes imported into Cytoscape after exporting the
# STRING interaction table from the STRING web interface
string_cytoscape_nodes <- fisher_out[, c(
  "gene",
  "direction",
  "logFC_272769",
  logFC_COL_95233,
  "fisher_fdr",
  "re_estimate",
  "I2"
)]

colnames(string_cytoscape_nodes) <- c(
  "gene",
  "direction",
  "logFC_GSE272769",
  "logFC_GSE95233_D01",
  "Fisher_FDR",
  "random_effect_estimate",
  "I2_percent"
)

string_cytoscape_nodes <- string_cytoscape_nodes[
  order(string_cytoscape_nodes$direction,
        string_cytoscape_nodes$gene),
]

rownames(string_cytoscape_nodes) <- NULL

stopifnot(
  nrow(string_cytoscape_nodes) == 32,
  !anyDuplicated(string_cytoscape_nodes$gene),
  !anyNA(string_cytoscape_nodes$gene)
)

write.csv(
  string_cytoscape_nodes,
  file.path(
    S3_DIR,
    "STRING_Cytoscape_node_attributes_32gene_signature.csv"
  ),
  row.names = FALSE
)

cat(sprintf(
  "STRING/Cytoscape input exported: %d genes\n",
  nrow(string_cytoscape_nodes)
))

# ---- PCA 发现集 ----
cat("\n--- Discovery PCA ---\n")
pca_272769 <- run_pca(gmat_272769, common_genes, group_272769, "GSE272769")
pca_272769_df <- data.frame(PC1=pca_272769$x[,1],PC2=pca_272769$x[,2],Outcome=group_272769)
p_pca_d1 <- ggplot(pca_272769_df,aes(PC1,PC2,color=Outcome))+
  geom_point(alpha=0.5,size=1)+stat_ellipse(level=0.95,lwd=0.5)+
  scale_color_manual(values=c(Survivor="#4E79A7",NonSurvivor="#E15759"))+
  labs(title="GSE272769 - PCA by Outcome",subtitle=sprintf("%d genes, %d samples | PC1=%.1f%% PC2=%.1f%%",
                                                           length(common_genes),ncol(gmat_272769),summary(pca_272769)$importance[2,1]*100,summary(pca_272769)$importance[2,2]*100))+
  theme_bw(11)+theme(legend.position="bottom")
save_plot(p_pca_d1,"Figure_PCA_GSE272769",S2_DIR,7,6)

pca_95233 <- run_pca(gmat_95233, common_genes, group_95233, "GSE95233")
pca_95233_df <- data.frame(PC1=pca_95233$x[,1],PC2=pca_95233$x[,2],Outcome=group_95233)
p_pca_d2 <- ggplot(pca_95233_df,aes(PC1,PC2,color=Outcome))+
  geom_point(alpha=0.5,size=1)+stat_ellipse(level=0.95,lwd=0.5)+
  scale_color_manual(values=c(Survivor="#4E79A7",NonSurvivor="#E15759"))+
  labs(title="GSE95233 D01 - PCA by Outcome",subtitle=sprintf("%d genes, %d samples | PC1=%.1f%% PC2=%.1f%%",
                                                              length(common_genes),ncol(gmat_95233),summary(pca_95233)$importance[2,1]*100,summary(pca_95233)$importance[2,2]*100))+
  theme_bw(11)+theme(legend.position="bottom")
save_plot(p_pca_d2,"Figure_PCA_GSE95233",S2_DIR,7,6)

# ---- Save ----
write.csv(meta_df, file.path(S2_DIR,sprintf("meta_all_%d_genes.csv",nrow(meta_df))),row.names=FALSE)
write.csv(fisher_out, file.path(S2_DIR,sprintf("meta_significant_%d_genes.csv",nrow(fisher_out))),row.names=FALSE)

saveRDS(list(
  meta_df=meta_df, fisher_out=fisher_out,
  gmat_272769=gmat_272769, gmat_95233=gmat_95233,
  gmat_272769_linear=gmat_272769_linear, gmat_95233_linear=gmat_95233_linear,
  group_272769=group_272769, group_95233=group_95233,
  common_genes=common_genes, logFC_COL_95233=logFC_COL_95233,
  pval_COL_95233=pval_COL_95233, module_map=module_map,
  res_272769=res_272769, res_95233=res_95233, dir_rate=dir_rate,
  n_up=n_up, n_down=n_down, n_both=n_both
), file.path(S2_DIR,"checkpoint_s2.rds"))
cat("§2 complete.\n")

# ============================================================
# §3 — Discovery figures + GO/KEGG
# ============================================================
cat("\n========== §3 Discovery Figures + GO/KEGG ==========\n")

# Scatter
plot_data <- meta_df
plot_data$sig <- "Not significant"
plot_data$sig[meta_df$sign_agree&abs(meta_df$logFC_272769)>LFC_CUT&abs(meta_df[[logFC_COL_95233]])>LFC_CUT&meta_df$fisher_fdr<FDR_CUT&meta_df$logFC_272769>0]<-"Up"
plot_data$sig[meta_df$sign_agree&abs(meta_df$logFC_272769)>LFC_CUT&abs(meta_df[[logFC_COL_95233]])>LFC_CUT&meta_df$fisher_fdr<FDR_CUT&meta_df$logFC_272769<0]<-"Down"
plot_data$sig <- factor(plot_data$sig,c("Not significant","Up","Down"))
plot_data <- plot_data[!is.na(plot_data$logFC_272769)&!is.na(plot_data[[logFC_COL_95233]]),]
p_scatter <- ggplot(plot_data,aes(logFC_272769,.data[[logFC_COL_95233]],color=sig))+
  geom_point(alpha=0.4,size=0.6)+scale_color_manual(values=cols_sc)+
  geom_hline(yintercept=0,lwd=0.25)+geom_vline(xintercept=0,lwd=0.25)+geom_abline(slope=1,lty="dotted",lwd=0.25)+
  labs(x=expression(log[2]~FC~GSE272769),y=expression(log[2]~FC~GSE95233),title="Cross-cohort logFC concordance",
       subtitle=sprintf("%d genes | Dir=%.1f%% | Fisher: %d Up, %d Down",nrow(meta_df),dir_rate,n_up,n_down))+
  theme_bw(11)+theme(legend.position="bottom",legend.title=element_blank())
save_plot(p_scatter,"Figure1_logFC_scatter")

# Volcano
plot_volcano <- function(
    res,
    ttl,
    fname
) {
  required_columns <- c(
    "logFC",
    "adj.P.Val"
  )
  
  if (!all(
    required_columns %in% colnames(res)
  )) {
    stop(
      "Volcano-plot input is missing required columns: ",
      paste(
        setdiff(
          required_columns,
          colnames(res)
        ),
        collapse = ", "
      )
    )
  }
  
  res$sig <- "Not significant"
  
  res$sig[
    !is.na(res$logFC) &
      !is.na(res$adj.P.Val) &
      res$logFC > LFC_CUT &
      res$adj.P.Val < FDR_CUT
  ] <- "Up"
  
  res$sig[
    !is.na(res$logFC) &
      !is.na(res$adj.P.Val) &
      res$logFC < -LFC_CUT &
      res$adj.P.Val < FDR_CUT
  ] <- "Down"
  
  res$sig <- factor(
    res$sig,
    levels = c(
      "Up",
      "Down",
      "Not significant"
    )
  )
  
  res$plot_FDR <- pmax(
    res$adj.P.Val,
    .Machine$double.xmin
  )
  
  n_up <- sum(
    res$sig == "Up",
    na.rm = TRUE
  )
  
  n_down <- sum(
    res$sig == "Down",
    na.rm = TRUE
  )
  
  p <- ggplot(
    res,
    aes(
      x = logFC,
      y = -log10(plot_FDR),
      color = sig
    )
  ) +
    geom_point(
      alpha = 0.55,
      size = 0.75
    ) +
    scale_color_manual(
      values = c(
        "Up" = "#E15759",
        "Down" = "#4E79A7",
        "Not significant" = "grey80"
      )
    ) +
    geom_hline(
      yintercept = -log10(FDR_CUT),
      linetype = "dashed",
      linewidth = 0.3,
      color = "grey50"
    ) +
    geom_vline(
      xintercept = c(
        -LFC_CUT,
        LFC_CUT
      ),
      linetype = "dashed",
      linewidth = 0.3,
      color = "grey50"
    ) +
    labs(
      title = ttl,
      subtitle = sprintf(
        "FDR < %.2f and |log2FC| > %.1f | Up: %d | Down: %d",
        FDR_CUT,
        LFC_CUT,
        n_up,
        n_down
      ),
      x = expression(log[2]~FC),
      y = expression(-log[10]~FDR),
      color = NULL
    ) +
    theme_bw(11) +
    theme(
      legend.position = "none"
    )
  
  save_plot(
    p,
    fname
  )
}
plot_volcano(res_272769,"GSE272769 (n=161, 60 NS vs 101 S)","FigureS1_volcano_GSE272769")
plot_volcano(res_95233,"GSE95233 D01 (n=51, 17 NS vs 34 S)","FigureS2_volcano_GSE95233")

# Venn
up1<-meta_df$gene[meta_df$logFC_272769>LFC_CUT&meta_df$fdr_272769<FDR_CUT]
dn1<-meta_df$gene[meta_df$logFC_272769<(-LFC_CUT)&meta_df$fdr_272769<FDR_CUT]
up2<-meta_df$gene[meta_df[[logFC_COL_95233]]>LFC_CUT&meta_df[[fdr_COL_95233]]<FDR_CUT]
dn2<-meta_df$gene[meta_df[[logFC_COL_95233]]<(-LFC_CUT)&meta_df[[fdr_COL_95233]]<FDR_CUT]
draw_venn_one<-function(set1,set2,title,file_prefix,fill_col="#E15759"){
  p<-venn.diagram(x=list(GSE272769=unique(set1),GSE95233=unique(set2)),filename=NULL,fill=c(fill_col,fill_col),alpha=c(0.4,0.7),cat.cex=1.3,cex=1.5,main=title,main.cex=1.3,scaled=FALSE,euler.d=FALSE)
  pdf(file.path(S3_DIR,paste0(file_prefix,".pdf")),7.5,7.5); grid.newpage(); grid.draw(p); dev.off()
  png(file.path(S3_DIR,paste0(file_prefix,".png")),2400,2400,res=300); grid.newpage(); grid.draw(p); dev.off()
}
draw_venn_one(
  up1,
  up2,
  sprintf(
    "Cohort-level upregulated genes\nFDR < %.2f; log2FC > %.1f",
    FDR_CUT,
    LFC_CUT
  ),
  "FigureS3a_Venn_Up",
  "#E15759"
)

draw_venn_one(
  dn1,
  dn2,
  sprintf(
    "Cohort-level downregulated genes\nFDR < %.2f; log2FC < -%.1f",
    FDR_CUT,
    LFC_CUT
  ),
  "FigureS3b_Venn_Down",
  "#4E79A7"
)

# I2 barplot
if (nrow(fisher_out) > 0) {
  fisher_out$I2_group <- cut(fisher_out$I2, I2_BREAKS, I2_LABELS, right = FALSE)
  fisher_out$I2_group <- as.character(fisher_out$I2_group)
  fisher_out$I2_group[is.na(fisher_out$I2_group)] <- "NA"
  
  i2d <- fisher_out[order(fisher_out$I2), ]
  i2d$gene <- factor(i2d$gene, levels = i2d$gene)
  
  p_i2 <- ggplot(i2d, aes(gene, I2, fill = I2_group)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = cols_i2) +
    geom_hline(yintercept = c(25, 50, 75), lty = 2, lwd = 0.3, color = "grey40") +
    coord_flip() +
    labs(
      y = "I2 (%)",
      x = NULL,
      title = "32-gene signature: I2 heterogeneity",
      subtitle = sprintf("k=2 | %d/%d I2<25%%", sum(fisher_out$I2 < 25, na.rm = TRUE), nrow(fisher_out))
    ) +
    theme_bw(10) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major.y = element_blank()
    )
  
  save_plot(p_i2, "FigureS4_I2_barplot", w = 6, h = 8)
}

# Module heatmap (2-cohort)
if (nrow(fisher_out)>=3) {
  plot_genes<-fisher_out$gene; plot_module<-ifelse(plot_genes%in%names(module_map),module_map[plot_genes],"Unclassified")
  plot_module<-factor(plot_module,levels=MODULE_LEVELS)
  ord<-order(plot_module,plot_genes); plot_genes<-plot_genes[ord]; plot_module<-plot_module[ord]
  mod_counts<-table(plot_module); n_seen<-0; gaps_row<-c()
  for(m in MODULE_LEVELS[1:4]){n_seen<-n_seen+ifelse(is.na(mod_counts[m]),0L,mod_counts[m]); if(!is.na(mod_counts[m])&&mod_counts[m]>0&&n_seen<length(plot_genes)) gaps_row<-c(gaps_row,n_seen)}
  if(length(gaps_row)==0) gaps_row<-NULL
  annot<-fisher_out[match(plot_genes,fisher_out$gene),]; stopifnot(identical(as.character(annot$gene),plot_genes))
  mat<-as.matrix(annot[,c("logFC_272769",logFC_COL_95233)]); rownames(mat)<-plot_genes; colnames(mat)<-c("GSE272769","GSE95233 D01")
  i2_group <- as.character(cut(
    annot$I2,
    breaks = I2_BREAKS,
    labels = I2_LABELS,
    right = FALSE
  ))
  i2_group[is.na(i2_group)] <- "NA"
  
  row_annot <- data.frame(
    Module = plot_module,
    `Individual-cohort support` = factor(
      as.character(annot$individual_cohort_support),
      levels = names(support_cols)
    ),
    I2_group = factor(
      i2_group,
      levels = c(I2_LABELS, "NA")
    ),
    row.names = plot_genes,
    check.names = FALSE
  )
  
  ann_colors <- list(
    Module = mod_cols,
    `Individual-cohort support` = support_cols,
    I2_group = cols_i2
  )
  max_abs<-max(abs(mat),na.rm=TRUE); breaks<-seq(-max_abs,max_abs,length.out=101); palette<-colorRampPalette(c("#4E79A7","white","#E15759"))(100)
  draw_hm<-function(file,dev=c("pdf","png"),w=10,h=12,r=300){
    if(dev=="pdf") pdf(file,w,h) else png(file,w,h,units="in",res=r)
    args<-list(mat,color=palette,breaks=breaks,cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=row_annot,annotation_colors=ann_colors,display_numbers=TRUE,number_format="%.2f",number_color="grey30",fontsize_number=6,fontsize_row=7,fontsize_col=10,angle_col=0,border_color=NA,annotation_legend=TRUE,legend_breaks=seq(-max_abs,max_abs,length.out=5),legend_labels=sprintf("%.1f",seq(-max_abs,max_abs,length.out=5)),main=sprintf("log2 FC (NS vs S) | %d-gene signature",nrow(mat)))
    if(!is.null(gaps_row)) args$gaps_row<-gaps_row
    do.call(pheatmap,args); dev.off()
  }
  draw_hm(file.path(S3_DIR,"Figure_Module_heatmap.pdf"),"pdf"); draw_hm(file.path(S3_DIR,"Figure_Module_heatmap.png"),"png")
  export <- data.frame(
    gene = plot_genes,
    module = as.character(plot_module),
    individual_cohort_support =
      as.character(row_annot[["Individual-cohort support"]]),
    I2 = annot$I2,
    I2_group = as.character(row_annot$I2_group),
    logFC_272769 = annot$logFC_272769,
    logFC_95233_D01 = annot[[logFC_COL_95233]],
    stringsAsFactors = FALSE
  )
  write.csv(export,file.path(S3_DIR,"Figure_Module_heatmap_source_data.csv"),row.names=FALSE)
}

# Forest plot
if (nrow(fisher_out)>=3) {
  forest_df<-fisher_out[order(fisher_out$direction,-fisher_out$mean_abs_logFC),]; forest_df$gene<-as.character(forest_df$gene)
  plot_data<-data.frame(gene=c(forest_df$gene,forest_df$gene,forest_df$gene),logFC=c(forest_df$logFC_272769,forest_df[[logFC_COL_95233]],forest_df$re_estimate),se=c(forest_df$se_272769,forest_df$se_95233,forest_df$re_se),study=rep(c("GSE272769","GSE95233 D01","RE summary"),each=nrow(forest_df)),direction=rep(forest_df$direction,3),stringsAsFactors=FALSE)
  plot_data$lower<-plot_data$logFC-1.96*plot_data$se; plot_data$upper<-plot_data$logFC+1.96*plot_data$se
  plot_data<-plot_data[is.finite(plot_data$logFC)&is.finite(plot_data$se)&is.finite(plot_data$lower)&is.finite(plot_data$upper),,drop=FALSE]
  plot_data$gene<-factor(plot_data$gene,levels=unique(forest_df$gene))
  study_cols<-c("GSE272769"="#E15759","GSE95233 D01"="#4E79A7","RE summary"="black")
  study_shapes<-c("GSE272769"=16,"GSE95233 D01"=17,"RE summary"=15)
  for(dir in c("Up","Down")) {
    sub<-plot_data[plot_data$direction==dir,]; if(nrow(sub)==0) next
    n_genes<-length(unique(sub$gene))
    p<-ggplot(sub,aes(logFC,gene,color=study,shape=study))+geom_vline(xintercept=0,lty="dashed",lwd=0.3,col="grey50")+
      geom_errorbar(
        aes(
          xmin = lower,
          xmax = upper
        ),
        width = 0.2,
        orientation = "y",
        position = position_dodge(width = 0.55)
      ) +
      geom_point(
        size = 2.5,
        position = position_dodge(width = 0.55)
      ) +
      scale_color_manual(values = study_cols) +scale_shape_manual(values=study_shapes)+
      labs(x="log2 Fold Change (NS vs S)",y=NULL,title=sprintf("Directional effect-size plot | %s-regulated (n=%d)",dir,n_genes),
           subtitle=sprintf("GSE272769 (n=%d) | GSE95233 D01 (n=%d) | RE",length(group_272769),length(group_95233)))+theme_bw(11)+theme(legend.position="bottom")
    ggsave(file.path(S3_DIR,sprintf("Figure_Forest_%s.pdf",tolower(dir))),p,width=9,height=n_genes*0.35+2)
    ggsave(file.path(S3_DIR,sprintf("Figure_Forest_%s.png",tolower(dir))),p,width=9,height=n_genes*0.35+2,dpi=300)
  }
}

# ---- GO/KEGG ----
cat("\n--- GO/KEGG enrichment ---\n")
# Legacy symbols were retained in expression and signature outputs.
# The following one-to-one updates are used only for GO/KEGG annotation.
annotation_symbol_updates <- c(
  "HIST1H2BM" = "H2BC14",
  "HIST1H3B" = "H3C2"
)

update_annotation_symbols <- function(symbols) {
  symbols <- as.character(symbols)
  
  replacement <- unname(
    annotation_symbol_updates[symbols]
  )
  
  has_replacement <- !is.na(replacement)
  
  symbols[has_replacement] <- replacement[has_replacement]
  
  symbols
}

annotation_symbol_update_table <- data.frame(
  original_symbol = names(annotation_symbol_updates),
  annotation_symbol = unname(annotation_symbol_updates),
  stringsAsFactors = FALSE
)

write.csv(
  annotation_symbol_update_table,
  file.path(S3_DIR, "GO_KEGG_symbol_updates.csv"),
  row.names = FALSE
)
common_universe <- unique(
  update_annotation_symbols(common_genes)
)

universe_eg <- unique(
  bitr(
    common_universe,
    "SYMBOL",
    "ENTREZID",
    org.Hs.eg.db
  )$ENTREZID
)

enrich_one <- function(genes, direction) {
  genes <- unique(
    update_annotation_symbols(genes)
  )
  if(length(genes)==0){cat(sprintf("  %s: 0 genes\n",direction));return(list(GO=NULL,KEGG=NULL))}
  eg<-unique(bitr(genes,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID)
  if(length(eg)==0){cat(sprintf("  %s: no ENTREZ\n",direction));return(list(GO=NULL,KEGG=NULL))}
  go_obj<-tryCatch(enrichGO(gene=eg,universe=universe_eg,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="BP",pAdjustMethod="BH",qvalueCutoff=FDR_CUT,readable=TRUE),error=function(e){cat(sprintf("  GO %s failed\n",direction));NULL})
  kegg_obj <- tryCatch(
    enrichKEGG(
      gene = eg,
      universe = universe_eg,
      organism = "hsa",
      pAdjustMethod = "BH",
      qvalueCutoff = FDR_CUT
    ),
    error = function(e) {
      cat(sprintf(
        "  KEGG %s failed: %s\n",
        direction,
        conditionMessage(e)
      ))
      NULL
    }
  )
  list(GO=go_obj,KEGG=kegg_obj)
}
up_annotation_symbols <- update_annotation_symbols(
  fisher_out$gene[
    fisher_out$direction == "Up"
  ]
)

up_annotation_mapping <- suppressWarnings(
  bitr(
    up_annotation_symbols,
    "SYMBOL",
    "ENTREZID",
    org.Hs.eg.db
  )
)

stopifnot(
  length(unique(up_annotation_mapping$SYMBOL)) == 29
)
enrich_up<-enrich_one(fisher_out$gene[fisher_out$direction=="Up"],"Up")
enrich_down<-enrich_one(fisher_out$gene[fisher_out$direction=="Down"],"Down")

# GO visualization
for (d in c("Up","Down")) {
  go_obj <- if (d=="Up") enrich_up$GO else enrich_down$GO
  if (!is.null(go_obj) && nrow(as.data.frame(go_obj))>0) {
    go_df <- as.data.frame(go_obj)
    write.csv(go_df, file.path(S3_DIR,sprintf("GO_BP_%s.csv",d)), row.names=FALSE)
    pg <- ggplot(go_df[1:min(10,nrow(go_df)),], aes(Count, reorder(Description,Count), size=Count, color=p.adjust))+
      geom_point()+scale_color_gradient(low="#E15759",high="#4E79A7",trans="log10")+scale_size(range=c(3,7))+
      labs(title=sprintf("GO BP - %s-regulated",d),x="Gene count",y="")+theme_bw(11)
    save_plot(pg, sprintf("Figure_GO_BP_%s",d), w=10, h=min(6, nrow(go_df)*0.3+2))
    cat(sprintf("GO %s: %d terms\n",d,nrow(go_df)))
  }
}

# KEGG visualization
for (d in c("Up","Down")) {
  kegg_obj <- if (d=="Up") enrich_up$KEGG else enrich_down$KEGG
  if (!is.null(kegg_obj) && nrow(as.data.frame(kegg_obj))>0) {
    kegg_df <- as.data.frame(kegg_obj)
    write.csv(kegg_df, file.path(S3_DIR,sprintf("KEGG_%s.csv",d)), row.names=FALSE)
    pk <- ggplot(kegg_df, aes(Count, reorder(Description,Count), size=Count, color=p.adjust))+
      geom_point()+scale_color_gradient(low="#E15759",high="#4E79A7",trans="log10")+scale_size(range=c(3,7))+
      labs(title=sprintf("KEGG - %s-regulated",d),x="Gene count",y="")+theme_bw(11)
    save_plot(pk, sprintf("Figure_KEGG_%s",d), w=10, h=min(4, nrow(kegg_df)*0.4+2))
    cat(sprintf("KEGG %s: %d pathways\n",d,nrow(kegg_df)))
  }
}

saveRDS(list(enrich_up=enrich_up, enrich_down=enrich_down), file.path(S3_DIR,"checkpoint_s3.rds"))
cat("§3 complete.\n")

# ============================================================
# §4 — GSE65682 external validation
# ============================================================
cat("\n========== §4 GSE65682 Validation ==========\n")

if (!file.exists(gse65682_file)) {
  stop(sprintf("Missing local GSE65682 Series Matrix file: %s", gse65682_file))
}
gset <- getGEO(filename=gse65682_file, getGPL=FALSE)
pdata <- pData(gset)
common_samples <- intersect(colnames(exprs(gset)), rownames(pdata))
stopifnot(length(common_samples)>0)
gset <- gset[,common_samples]; pdata <- pdata[common_samples,,drop=FALSE]
stopifnot(identical(colnames(exprs(gset)),rownames(pdata)))
cat(sprintf("GSE65682: %d probes x %d samples\n", nrow(exprs(gset)), ncol(exprs(gset))))

if (!file.exists(gpl13667_file)) {
  stop(sprintf("Missing local GPL13667 platform annotation file: %s", gpl13667_file))
}
gpl13667 <- getGEO(filename=gpl13667_file,GSEMatrix=FALSE)
tab13667 <- Table(gpl13667)
gene_col <- intersect(c("Gene Symbol","GENE_SYMBOL","Gene_Symbol","Symbol"),colnames(tab13667))[1]
stopifnot(!is.na(gene_col))
map13667 <- data.frame(probe=as.character(tab13667$ID),gene=parse_gene_symbol(tab13667[[gene_col]]),stringsAsFactors=FALSE)
map13667 <- map13667[!invalid_gene(map13667$gene),]

mort <- trimws(
  as.character(
    pdata$`mortality_event_28days:ch1`
  )
)

endo <- trimws(
  as.character(
    pdata$`endotype_class:ch1`
  )
)

valid_mortality <- mort %in% c(
  "0",
  "1"
)

valid_endotype <- endo %in% c(
  "Mars1",
  "Mars2",
  "Mars3",
  "Mars4"
)

include_65682 <- valid_mortality & valid_endotype

selection_summary_65682 <- data.frame(
  category = c(
    "All GSE65682 series-matrix samples",
    "Valid 28-day mortality information",
    "Valid MARS endotype assignment",
    "Included: valid mortality and MARS endotype",
    "Excluded: invalid or missing mortality only",
    "Excluded: invalid or missing endotype only",
    "Excluded: invalid or missing mortality and endotype"
  ),
  n_samples = c(
    length(mort),
    sum(valid_mortality),
    sum(valid_endotype),
    sum(include_65682),
    sum(!valid_mortality & valid_endotype),
    sum(valid_mortality & !valid_endotype),
    sum(!valid_mortality & !valid_endotype)
  ),
  stringsAsFactors = FALSE
)

sample_audit_65682 <- data.frame(
  sample_id = colnames(exprs(gset)),
  mortality_event_28days = mort,
  endotype_class = endo,
  valid_mortality = valid_mortality,
  valid_endotype = valid_endotype,
  included = include_65682,
  stringsAsFactors = FALSE
)

write.csv(
  selection_summary_65682,
  file.path(
    S4_DIR,
    "GSE65682_sample_selection_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_audit_65682,
  file.path(
    S4_DIR,
    "GSE65682_sample_inclusion_audit.csv"
  ),
  row.names = FALSE
)

cat(
  "\nGSE65682 sample selection summary:\n"
)

print(
  selection_summary_65682,
  row.names = FALSE
)

keep_idx <- which(include_65682)

stopifnot(
  length(keep_idx) > 0
)
gset_sepsis <- gset[,keep_idx]; pdata_sepsis <- pdata[keep_idx,]
group_val <- check_group(pdata_sepsis$`mortality_event_28days:ch1`,name_map=list(NonSurvivor="1",Survivor="0"),label="GSE65682")
endo_val  <- trimws(as.character(pdata_sepsis$`endotype_class:ch1`))

exprs_raw <- exprs(gset_sepsis)

range_65682 <- range(
  exprs_raw,
  na.rm = TRUE
)

# Reuse the is_log2_like() function defined during
# discovery-cohort preprocessing.
is_log2_65682 <- is_log2_like(
  exprs_raw
)

cat(sprintf(
  "GSE65682 range: %.3f to %.3f | log2-like: %s\n",
  range_65682[1],
  range_65682[2],
  is_log2_65682
))

if (!is_log2_65682) {
  stop(
    paste0(
      "GSE65682 was expected to contain log2-scale ",
      "processed expression values. Explicit inspection ",
      "and transformation are required before limma analysis."
    )
  )
}

exprs_log <- exprs_raw

# Select one representative probe per gene on the log2 scale.
gmat_val <- probes_to_genes(
  exprs_log,
  map13667
)

# Derive the xCell input from the same gene-level matrix.
# This ensures that limma and xCell use the same selected probe.
gmat_val_linear <- pmax(
  2^gmat_val - 1,
  0
)

stopifnot(
  identical(
    rownames(gmat_val),
    rownames(gmat_val_linear)
  ),
  identical(
    colnames(gmat_val),
    colnames(gmat_val_linear)
  ),
  !anyNA(gmat_val),
  !anyNA(gmat_val_linear),
  all(gmat_val_linear >= 0)
)
res_val <- run_limma(gmat_val, group_val, "GSE65682")

# Concordance
val_genes <- fisher_out$gene
genes_found <- intersect(val_genes, res_val$gene)
val_fc <- res_val[match(genes_found, res_val$gene),]
val_fc$disc_dir <- fisher_out$direction[match(genes_found, fisher_out$gene)]
val_fc$disc_sign <- ifelse(val_fc$disc_dir=="Up",1L,ifelse(val_fc$disc_dir=="Down",-1L,NA_integer_))
val_fc$agree_sign <- sign(val_fc$logFC)==val_fc$disc_sign
val_fc$agree_threshold <- val_fc$agree_sign & abs(val_fc$logFC)>LFC_CUT
n_valid <- sum(!is.na(val_fc$agree_sign)); val_fc$agree_sign[is.na(val_fc$agree_sign)]<-FALSE; val_fc$agree_threshold[is.na(val_fc$agree_threshold)]<-FALSE
cat(sprintf("\nGSE65682: %d genes (%d valid)\n",length(genes_found),n_valid))
cat(sprintf("  Sign: %d/%d (%.1f%%)\n",sum(val_fc$agree_sign),n_valid,100*sum(val_fc$agree_sign)/max(1,n_valid)))
cat(sprintf("  Threshold: %d/%d (%.1f%%)\n",sum(val_fc$agree_threshold),n_valid,100*sum(val_fc$agree_threshold)/max(1,n_valid)))

# Per-endotype
cat("\nPer-endotype:\n")
endo_results <- data.frame()
for (endo_i in c("Mars1","Mars2","Mars3","Mars4")) {
  ie<-which(endo_val==endo_i); ge<-group_val[ie]
  if(length(ie)<10||length(unique(ge))<2){ next }
  ee<-gmat_val[,ie]; de<-model.matrix(~0+ge); colnames(de)<-levels(ge)
  fe<-lmFit(ee,de); f2e<-contrasts.fit(fe,makeContrasts(NonSurvivor-Survivor,levels=de)); f2e<-eBayes(f2e,trend=TRUE)
  re<-topTable(f2e,coef=1,number=Inf,sort.by="none"); re$gene<-rownames(re); fce<-re[match(genes_found,re$gene),]
  matched<-!is.na(fce$logFC)&!is.na(val_fc$disc_sign); ae<-rep(NA,length(genes_found)); ae[matched]<-sign(fce$logFC[matched])==val_fc$disc_sign[matched]
  endo_results<-rbind(endo_results,data.frame(Endotype=endo_i,S=sum(ge=="Survivor"),NS=sum(ge=="NonSurvivor"),agree=sum(ae,na.rm=TRUE),total=sum(!is.na(ae)),pct=round(mean(ae,na.rm=TRUE)*100,1)))
}
print(endo_results)
write.csv(endo_results,file.path(S4_DIR,"GSE65682_per_endotype.csv"),row.names=FALSE)

if(nrow(endo_results)>=2) {
  ep <- endo_results[!is.na(endo_results$pct), ]
  ep$Endotype <- factor(ep$Endotype, levels = c("Mars1", "Mars2", "Mars3", "Mars4"))
  ep$label <- sprintf("%d/%d\n(%.0f%%)", ep$agree, ep$total, ep$pct)
  
  p_endo <- ggplot(ep, aes(Endotype, pct, fill = Endotype)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_hline(yintercept = 50, lty = "dotted", color = "grey50") +
    geom_text(aes(label = label, y = pmin(pct + 4, 106)), size = 3.4) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(0, 110), expand = expansion(mult = c(0, 0.03))) +
    labs(
      x = NULL,
      y = "Sign concordance (%)",
      title = "Discovery signature concordance by MARS endotype",
      subtitle = sprintf("%d detectable genes tested", n_valid)
    ) +
    theme_bw(11) +
    theme(legend.position = "none")
  
  save_plot(p_endo, "Figure_Val_per_endotype", S4_DIR, 7.5, 5.2)
}

# Volcano GSE65682
res_val$sig<-ifelse(res_val$adj.P.Val<FDR_CUT&abs(res_val$logFC)>LFC_CUT,ifelse(res_val$logFC>0,"Up","Down"),"NS")
res_val$sig<-factor(res_val$sig,c("Up","Down","NS"))
top_lab<-res_val[res_val$sig!="NS",]; top_lab<-top_lab[order(top_lab$adj.P.Val),]; top_lab<-rbind(head(top_lab[top_lab$sig=="Up",],15),head(top_lab[top_lab$sig=="Down",],15))
p_volc<-ggplot(res_val,aes(logFC,-log10(pmax(adj.P.Val,.Machine$double.xmin)),color=sig))+geom_point(alpha=0.5,size=0.8)+
  geom_vline(xintercept=c(-LFC_CUT,LFC_CUT),lty="dashed",color="grey50",lwd=0.3)+geom_hline(yintercept=-log10(FDR_CUT),lty="dashed",color="grey50",lwd=0.3)+
  geom_text_repel(data=top_lab,aes(label=gene),size=2.5,max.overlaps=25)+scale_color_manual(values=c(Up="#E15759",Down="#4E79A7",NS="grey70"))+
  labs(title="GSE65682 External Validation - Volcano",subtitle=sprintf("Up=%d Down=%d",sum(res_val$sig=="Up"),sum(res_val$sig=="Down")),x="log2 FC",y=expression(-log[10](FDR)))+theme_bw(11)+theme(legend.position="bottom")
save_plot(p_volc,"Figure_Val_GSE65682_volcano",S4_DIR,7.5,6)

# ============================================================
# §4 GSE65682 — 完整 PCA 套件 (A-L 共12张)
# ============================================================

genes_32 <- fisher_out$gene
genes_32_val <- intersect(genes_32, rownames(gmat_val))
common_genes_val <- intersect(common_genes, rownames(gmat_val))
endo_levels <- c("Mars1","Mars2","Mars3","Mars4")
outcome_cols <- c(Survivor="#4E79A7", NonSurvivor="#E15759")

# ---- A. All genes × All samples — PCA by Outcome ----
pca_val <- run_pca(gmat_val, common_genes_val, group_val, "GSE65682 all-genes")
pca_val_df <- data.frame(PC1=pca_val$x[,1], PC2=pca_val$x[,2], Outcome=group_val, Endotype=endo_val)

pA <- ggplot(pca_val_df, aes(PC1, PC2, color=Outcome)) +
  geom_point(alpha=0.5, size=0.8) + stat_ellipse(level=0.95, lwd=0.5) +
  scale_color_manual(values=outcome_cols) +
  labs(title="A. All genes x All - Outcome",
       subtitle=sprintf("%d genes | PC1=%.1f%% PC2=%.1f%% | S=%d NS=%d",
                        length(common_genes_val), summary(pca_val)$importance[2,1]*100,
                        summary(pca_val)$importance[2,2]*100,
                        sum(group_val=="Survivor"), sum(group_val=="NonSurvivor"))) +
  theme_bw(10) + theme(legend.position="bottom")

# ---- B. All genes × All samples — PCA by Endotype ----
pB <- ggplot(pca_val_df, aes(PC1, PC2, color=Endotype)) +
  geom_point(alpha=0.5, size=0.8) + scale_color_brewer(palette="Set2") +
  labs(title="B. All genes x All - Endotype",
       subtitle=sprintf("M1=%d M2=%d M3=%d M4=%d",
                        sum(endo_val=="Mars1"), sum(endo_val=="Mars2"),
                        sum(endo_val=="Mars3"), sum(endo_val=="Mars4"))) +
  theme_bw(10) + theme(legend.position="bottom")

# ---- C+D. Detectable signature genes × all samples ----
pca_32 <- prcomp(t(gmat_val[genes_32_val, , drop=FALSE]), scale.=TRUE)
pca_32_df <- data.frame(PC1=pca_32$x[,1], PC2=pca_32$x[,2], Outcome=group_val, Endotype=endo_val)

pC <- ggplot(pca_32_df, aes(PC1, PC2, color=Outcome)) +
  geom_point(alpha=0.6, size=1.2) + stat_ellipse(level=0.95, lwd=0.6) +
  scale_color_manual(values=outcome_cols) +
  labs(
    title = sprintf(
      "C. Detectable signature genes (%d/%d) - Outcome",
      length(genes_32_val),
      length(genes_32)
    ),
       subtitle=sprintf("%d/%d genes | PC1=%.1f%% PC2=%.1f%% | S=%d NS=%d",
                        length(genes_32_val), length(genes_32),
                        summary(pca_32)$importance[2,1]*100, summary(pca_32)$importance[2,2]*100,
                        sum(group_val=="Survivor"), sum(group_val=="NonSurvivor"))) +
  theme_bw(10) + theme(legend.position="bottom")
pD <- ggplot(
  pca_32_df,
  aes(
    PC1,
    PC2,
    color = Endotype
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.2
  ) +
  stat_ellipse(
    level = 0.95,
    lwd = 0.6
  ) +
  scale_color_brewer(
    palette = "Set2"
  ) +
  labs(
    title = sprintf(
      "D. Detectable signature genes (%d/%d) - Endotype",
      length(genes_32_val),
      length(genes_32)
    ),
    subtitle = sprintf(
      "PC1=%.1f%% PC2=%.1f%%",
      summary(pca_32)$importance[2, 1] * 100,
      summary(pca_32)$importance[2, 2] * 100
    )
  ) +
  theme_bw(10) +
  theme(
    legend.position = "bottom"
  )
# Save A-D
p_AD <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = sprintf(
      paste0(
        "GSE65682 PCA: all genes (top) versus ",
        "%d detectable genes from the %d-gene signature (bottom)"
      ),
      length(genes_32_val),
      length(genes_32)
    )
  )
pdf(file.path(S4_DIR, "Figure_PCA_GSE65682_A_D.pdf"), 14, 12); print(p_AD); dev.off()
png(file.path(S4_DIR, "Figure_PCA_GSE65682_A_D.png"), 14, 12, units="in", res=300); print(p_AD); dev.off()
cat("A-D saved.\n")

# ---- E-H. Per-endotype PCA (ALL genes, NS vs S) ----
plots_EH <- list()
for (e in endo_levels) {
  idx <- which(endo_val == e); grp <- group_val[idx]; tbl <- table(grp)
  if (any(tbl < 2)) next
  pca <- prcomp(t(gmat_val[common_genes_val, idx, drop=FALSE]), scale.=TRUE)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Outcome=grp)
  plots_EH[[e]] <- ggplot(df, aes(PC1, PC2, color=Outcome)) +
    geom_point(alpha=0.6, size=1) + stat_ellipse(level=0.95, lwd=0.5) +
    scale_color_manual(values=outcome_cols) +
    labs(title=sprintf("%s - All genes", e),
         subtitle=sprintf("S=%d NS=%d | PC1=%.1f%% PC2=%.1f%%",
                          tbl["Survivor"], tbl["NonSurvivor"],
                          summary(pca)$importance[2,1]*100, summary(pca)$importance[2,2]*100)) +
    theme_bw(10) + theme(legend.position="bottom")
}
p_EH <- wrap_plots(plots_EH, ncol=2) +
  plot_annotation(title="E-H. Per-Endotype PCA: All Common Genes",
                  subtitle=sprintf("%d genes", length(common_genes_val)))
pdf(file.path(S4_DIR, "Figure_PCA_per_endotype_all_genes.pdf"), 10, 9); print(p_EH); dev.off()
png(file.path(S4_DIR, "Figure_PCA_per_endotype_all_genes.png"), 10, 9, units="in", res=300); print(p_EH); dev.off()
cat("E-H saved.\n")

# ---- I-L. Per-endotype PCA using detectable signature genes ----
plots_IL <- list()
for (e in endo_levels) {
  idx <- which(endo_val == e); grp <- group_val[idx]; tbl <- table(grp)
  if (any(tbl < 2)) next
  pca <- prcomp(t(gmat_val[genes_32_val, idx, drop=FALSE]), scale.=TRUE)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Outcome=grp)
  plots_IL[[e]] <- ggplot(df, aes(PC1, PC2, color=Outcome)) +
    geom_point(alpha=0.7, size=1.5) + stat_ellipse(level=0.95, lwd=0.6) +
    scale_color_manual(values=outcome_cols) +
    labs(
      title = sprintf(
        "%s - %d detectable signature genes",
        e,
        length(genes_32_val)
      ),
         subtitle=sprintf("S=%d NS=%d | PC1=%.1f%% PC2=%.1f%%",
                          tbl["Survivor"], tbl["NonSurvivor"],
                          summary(pca)$importance[2,1]*100, summary(pca)$importance[2,2]*100)) +
    theme_bw(10) + theme(legend.position="bottom")
}
p_IL <- wrap_plots(plots_IL, ncol=2) +
  plot_annotation(
    title = sprintf(
      "I-L. Per-endotype PCA: %d detectable genes from the %d-gene signature",
      length(genes_32_val),
      length(genes_32)
    )
  )
pdf(file.path(S4_DIR, "Figure_PCA_per_endotype_detectable_signature_genes.pdf"), 10, 9); print(p_IL); dev.off()
png(file.path(S4_DIR, "Figure_PCA_per_endotype_detectable_signature_genes.png"), 10, 9, units="in", res=300); print(p_IL); dev.off()
cat("I-L saved.\n")
# ========== 保存全部PCA数据 ==========

pca_dir <- file.path(S4_DIR, "PCA_data")
dir.create(pca_dir, showWarnings=FALSE)

# --- A+C: all-sample PCA ---
save_pca_data <- function(pca_obj, scores_df, name) {
  # scores
  write.csv(scores_df, file.path(pca_dir, sprintf("%s_scores.csv", name)), row.names=FALSE)
  # variance
  var_df <- data.frame(
    PC = paste0("PC", seq_along(pca_obj$sdev)),
    SD = pca_obj$sdev,
    Variance = pca_obj$sdev^2 / sum(pca_obj$sdev^2) * 100,
    Cumulative = cumsum(pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100
  )
  write.csv(var_df, file.path(pca_dir, sprintf("%s_variance.csv", name)), row.names=FALSE)
  cat(sprintf("  %s: %d PCs saved\n", name, nrow(var_df)))
}

# A+B shared pca_val
save_pca_data(pca_val, pca_val_df, "A_B_allgenes_Outcome_Endotype")
# C+D shared pca_32
save_pca_data(pca_32, pca_32_df, "C_D_32genes_Outcome_Endotype")

# --- E-H: per-endotype all-genes ---
for (e in endo_levels) {
  idx <- which(endo_val == e); grp <- group_val[idx]
  if (any(table(grp) < 2)) next
  pca <- prcomp(t(gmat_val[common_genes_val, idx, drop=FALSE]), scale.=TRUE)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Outcome=grp)
  save_pca_data(pca, df, sprintf("%s_allgenes", e))
}

# --- I-L: per-endotype 32-genes ---
for (e in endo_levels) {
  idx <- which(endo_val == e); grp <- group_val[idx]
  if (any(table(grp) < 2)) next
  pca <- prcomp(t(gmat_val[genes_32_val, idx, drop=FALSE]), scale.=TRUE)
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Outcome=grp)
  save_pca_data(pca, df, sprintf("%s_32genes", e))
}

cat(sprintf("All PCA data saved to %s\n", pca_dir))

# ---- 32-gene triple-cohort heatmap ----
mod_df <- data.frame(
  gene = c(MODULE_NEUTRO, MODULE_CELL_CYCLE, MODULE_DOWN_DEFENSE, MODULE_OTHER),
  module = c(
    rep("Neutrophil degranulation",  length(MODULE_NEUTRO)),
    rep("Cell cycle / proliferation", length(MODULE_CELL_CYCLE)),
    rep("Inflammatory / Down",        length(MODULE_DOWN_DEFENSE)),
    rep("Other",                      length(MODULE_OTHER))
  ),
  stringsAsFactors = FALSE
)
module_colors <- c("Cell cycle / proliferation"="#4E79A7","Neutrophil degranulation"="#E15759","Inflammatory / Down"="#59A14F","Other"="#B07AA1")
meta_sub  <- meta_df[meta_df$gene%in%mod_df$gene,c("gene","logFC_272769")]; names(meta_sub)[2]<-"GSE272769"
fc_95233  <- meta_df[meta_df$gene%in%mod_df$gene,c("gene",logFC_COL_95233)]; names(fc_95233)[2]<-"GSE95233"
fc_65682  <- res_val[res_val$gene%in%mod_df$gene,c("gene","logFC")]; names(fc_65682)[2]<-"GSE65682"
fc_mat <- merge(meta_sub, fc_95233, by="gene")
fc_mat <- merge(fc_mat, fc_65682, by="gene", all.x=TRUE)
fc_mat <- merge(fc_mat, mod_df, by="gene")
fc_mat$module <- factor(fc_mat$module,levels=names(module_colors))
fc_mat$mean_abs <- rowMeans(abs(fc_mat[,c("GSE272769","GSE95233","GSE65682")]),na.rm=TRUE)
fc_mat <- fc_mat[order(fc_mat$module,-fc_mat$mean_abs),]; fc_mat$gene<-factor(fc_mat$gene,levels=fc_mat$gene)
fc_long <- reshape(fc_mat,varying=c("GSE272769","GSE95233","GSE65682"),v.names="logFC",timevar="cohort",times=c("GSE272769","GSE95233","GSE65682"),direction="long")
fc_long$cohort <- factor(fc_long$cohort,levels=c("GSE272769","GSE95233","GSE65682"))
fc_lim <- max(abs(fc_long$logFC),na.rm=TRUE)*1.05
color_bar <- fc_mat[,c("gene","module")]; color_bar$gene<-factor(color_bar$gene,levels=levels(fc_mat$gene))
p_bar <- ggplot(color_bar,aes(x=1,y=gene,fill=module))+geom_tile()+scale_fill_manual(values=module_colors,guide="none")+scale_x_continuous(expand=c(0,0))+theme_void()
p_heat <- ggplot(fc_long,aes(cohort,gene,fill=logFC))+geom_tile(color="white",linewidth=0.5)+
  geom_text(aes(label=ifelse(is.na(logFC),"-",sprintf("%.2f",logFC))),size=3.0)+
  scale_fill_gradient2(low="#4E79A7",mid="white",high="#E15759",limits=c(-fc_lim,fc_lim),na.value="grey90")+
  facet_grid(module~.,scales="free_y",space="free_y",switch="y")+scale_x_discrete(expand=c(0,0))+
  labs(title="32-gene signature: logFC across three cohorts",subtitle="Grouped by functional module | - = gene not detected on platform",x=NULL,y=NULL)+
  theme_minimal(12)+theme(panel.grid=element_blank(),axis.text.y=element_text(size=9,face="italic"),axis.text.x=element_text(size=10,face="bold"),strip.text.y.left=element_text(size=10,face="bold",angle=0,hjust=1),strip.placement="outside",legend.position="right",plot.title=element_text(face="bold",size=13))
p_final <- p_bar+p_heat+plot_layout(widths=c(0.03,1))
save_plot(p_final,"Figure_32gene_triple_heatmap",S4_DIR,14,10)
fc_out <- fc_mat[,c("gene","module","GSE272769","GSE95233","GSE65682","mean_abs")]; fc_out<-fc_out[order(fc_out$module,-fc_out$mean_abs),]; rownames(fc_out)<-NULL
write.csv(fc_out,file.path(S4_DIR,"Table_S4_32genes_triple_logFC.csv"),row.names=FALSE)

write.csv(res_val,file.path(S4_DIR,"DE_GSE65682.csv"),row.names=FALSE)
saveRDS(list(res_val=res_val,gmat_val=gmat_val,gmat_val_linear=gmat_val_linear,group_val=group_val,endo_val=endo_val,genes_found=genes_found,val_fc=val_fc,endo_results=endo_results),file.path(S4_DIR,"checkpoint_s4.rds"))
cat("§4 complete.\n")

# ============================================================
# §5 — Triple-cohort xCell + GSEA + TF
# ============================================================
cat("\n========== §5 Triple-cohort ==========\n")

# xCell
run_xcell<-function(emat_linear,group_vec,name){
  group_vec<-factor(group_vec,levels=c("Survivor","NonSurvivor"))
  if(min(emat_linear,na.rm=TRUE)<0){emat_linear[emat_linear<0]<-0}
  xc<-tryCatch(xCellAnalysis(emat_linear,parallel.sz=1),error=function(e)NULL)
  if(is.null(xc)) return(NULL)
  results<-do.call(rbind,lapply(rownames(xc),function(ct){
    scores<-as.numeric(xc[ct,]); ns_idx<-which(group_vec=="NonSurvivor"); s_idx<-which(group_vec=="Survivor")
    if(length(ns_idx)<2||length(s_idx)<2) return(NULL)
    wt<-tryCatch(wilcox.test(scores[ns_idx],scores[s_idx],exact=FALSE),error=function(e)NULL)
    data.frame(cell_type=ct,mean_NS=mean(scores[ns_idx]),mean_S=mean(scores[s_idx]),delta=mean(scores[ns_idx])-mean(scores[s_idx]),pval=if(is.null(wt))NA else wt$p.value,stringsAsFactors=FALSE)
  }))
  if(is.null(results)||nrow(results)==0) return(NULL)
  results$fdr<-p.adjust(results$pval,"BH"); results$cohort<-name; results[order(results$fdr),]
}
xc_272769 <- run_xcell(
  gmat_272769_linear,
  group_272769,
  "GSE272769"
)

xc_95233 <- run_xcell(
  gmat_95233_linear,
  group_95233,
  "GSE95233"
)

xc_65682 <- run_xcell(
  gmat_val_linear,
  group_val,
  "GSE65682"
)

xc_list <- list(
  GSE272769 = xc_272769,
  GSE95233 = xc_95233,
  GSE65682 = xc_65682
)

xc_ok <- vapply(
  xc_list,
  function(x) {
    !is.null(x) &&
      is.data.frame(x) &&
      nrow(x) > 0
  },
  logical(1)
)

cat("\nxCell cohort status:\n")
print(xc_ok)

if (!all(xc_ok)) {
  stop(
    paste0(
      "xCell failed or returned no results in: ",
      paste(
        names(xc_list)[!xc_ok],
        collapse = ", "
      ),
      ". Three-cohort xCell integration was not performed."
    )
  )
}

# Export cohort-specific results only after all three analyses succeed.
write.csv(
  xc_272769,
  file.path(
    S5_DIR,
    "xCell_GSE272769.csv"
  ),
  row.names = FALSE
)

write.csv(
  xc_95233,
  file.path(
    S5_DIR,
    "xCell_GSE95233.csv"
  ),
  row.names = FALSE
)

write.csv(
  xc_65682,
  file.path(
    S5_DIR,
    "xCell_GSE65682.csv"
  ),
  row.names = FALSE
)

# Merge the three cohort-specific results.
xc_all <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "cell_type",
      all = FALSE
    )
  },
  lapply(
    names(xc_list),
    function(nm) {
      df <- xc_list[[nm]][, c(
        "cell_type",
        "delta",
        "fdr"
      )]
      
      names(df)[2:3] <- paste0(
        c("delta_", "fdr_"),
        nm
      )
      
      df
    }
  )
)

xc_all <- as.data.frame(
  xc_all
)

expected_delta_cols <- c(
  "delta_GSE272769",
  "delta_GSE95233",
  "delta_GSE65682"
)

expected_fdr_cols <- c(
  "fdr_GSE272769",
  "fdr_GSE95233",
  "fdr_GSE65682"
)

stopifnot(
  all(expected_delta_cols %in% colnames(xc_all)),
  all(expected_fdr_cols %in% colnames(xc_all))
)

for (col in expected_fdr_cols) {
  xc_all[[col]][
    is.na(xc_all[[col]])
  ] <- 1
}

delta_cols <- expected_delta_cols

xc_all$mean_abs_delta <- rowMeans(
  abs(
    xc_all[, delta_cols, drop = FALSE]
  ),
  na.rm = TRUE
)

direction_matrix <- sign(
  as.matrix(
    xc_all[, delta_cols, drop = FALSE]
  )
)

xc_all$all_same_dir <- apply(
  direction_matrix,
  1,
  function(x) {
    all(is.finite(x)) &&
      length(unique(x)) == 1
  }
)

xc_all$n_same_pairs <- apply(
  direction_matrix,
  1,
  function(x) {
    sum(
      outer(x, x, "==")[
        upper.tri(
          outer(x, x, "==")
        )
      ]
    )
  }
)

xc_all$total_pairs <- 3L

cat(sprintf(
  paste0(
    "Triple-cohort xCell: %d cell types | ",
    "All three cohorts directionally concordant: %d (%.1f%%)\n"
  ),
  nrow(xc_all),
  sum(xc_all$all_same_dir),
  100 * mean(xc_all$all_same_dir)
))

write.csv(
  xc_all,
  file.path(
    S5_DIR,
    "xCell_triple_cohort.csv"
  ),
  row.names = FALSE
)

# xCell 9-type bar
selected<-c("CMP","GMP","HSC","Erythrocytes","iDC","CD4+ Tcm","CD4+ Tem","Th1 cells","NKT")
df_xc<-xc_all[xc_all$cell_type%in%selected,c("cell_type",delta_cols)]
df_long<-tidyr::pivot_longer(df_xc,cols=starts_with("delta_"),names_to="cohort",values_to="delta")
df_long$cohort<-dplyr::recode(df_long$cohort,"delta_GSE272769"="GSE272769","delta_GSE95233"="GSE95233","delta_GSE65682"="GSE65682")
df_long$cell_type<-factor(df_long$cell_type,levels=selected)
cohort_cols<-c("GSE272769"="#4E79A7","GSE95233"="#E15759","GSE65682"="#B07AA1")
p_xc<-ggplot(df_long,aes(cell_type,delta,fill=cohort))+geom_bar(stat="identity",position=position_dodge(0.8),width=0.7,color="white",linewidth=0.3)+
  geom_hline(yintercept=0,linewidth=0.5,color="gray40")+scale_fill_manual(values=cohort_cols)+
  labs(title="Triple-cohort concordant immune cell types",subtitle="9 of 24 concordant cell types",x=NULL,y=expression(Delta~"xCell score (NS - S)"))+
  theme_bw(12)+theme(legend.position="bottom")
save_plot(p_xc,"Figure_xCell_9types_bar",S5_DIR,10,6)

# GSEA Hallmark
msig_h <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
hallmark_list <- split(msig_h$gene_symbol, msig_h$gs_name)

build_ranked <- function(df, lfc_col, p_col = NULL) {
  if (!is.null(p_col) && p_col %in% names(df)) {
    pvec <- pmax(
      suppressWarnings(as.numeric(df[[p_col]])),
      1e-300
    )
    ranks <- -log10(pvec) *
      sign(suppressWarnings(as.numeric(df[[lfc_col]])))
  } else {
    ranks <- suppressWarnings(as.numeric(df[[lfc_col]]))
  }
  
  names(ranks) <- as.character(df$gene)
  ranks <- ranks[
    is.finite(ranks) &
      !is.na(names(ranks)) &
      names(ranks) != ""
  ]
  ranks <- ranks[!duplicated(names(ranks))]
  sort(ranks, decreasing = TRUE)
}

rank_272769 <- build_ranked(
  meta_df,
  "logFC_272769",
  "pval_272769"
)

rank_95233 <- build_ranked(
  meta_df,
  logFC_COL_95233,
  pval_COL_95233
)

rank_65682 <- build_ranked(
  res_val,
  "logFC"
)

run_fgsea <- function(ranked, pathways, name) {
  fg <- tryCatch(
    fgseaMultilevel(
      pathways = pathways,
      stats = ranked,
      minSize = 10,
      maxSize = 500
    ),
    error = function(e) {
      cat(sprintf("%s: FAILED — %s\n", name, conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(fg)) {
    return(NULL)
  }
  
  cat(sprintf("%s: OK\n", name))
  fg
}

set.seed(42)
fg_h_272769 <- run_fgsea(
  rank_272769,
  hallmark_list,
  "GSE272769"
)

set.seed(42)
fg_h_95233 <- run_fgsea(
  rank_95233,
  hallmark_list,
  "GSE95233"
)

set.seed(42)
fg_h_65682 <- run_fgsea(
  rank_65682,
  hallmark_list,
  "GSE65682"
)

fg_list <- list(
  GSE272769 = fg_h_272769,
  GSE95233 = fg_h_95233,
  GSE65682 = fg_h_65682
)

fg_ok <- vapply(
  fg_list,
  function(x) {
    !is.null(x) &&
      is.data.frame(x) &&
      nrow(x) > 0
  },
  logical(1)
)

cat("\nHallmark GSEA cohort status:\n")
print(fg_ok)

if (!all(fg_ok)) {
  stop(
    paste0(
      "Hallmark GSEA failed or returned no results in: ",
      paste(
        names(fg_list)[!fg_ok],
        collapse = ", "
      ),
      ". Three-cohort Hallmark integration was not performed."
    )
  )
}

# Retain only pathway, NES, and adjusted P value
# for the cross-cohort summary.
hall_3 <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "pathway",
      all = FALSE
    )
  },
  lapply(
    names(fg_list),
    function(nm) {
      df <- as.data.frame(
        fg_list[[nm]]
      )
      
      df <- df[, c(
        "pathway",
        "NES",
        "padj"
      )]
      
      names(df)[2:3] <- paste0(
        c("NES_", "padj_"),
        nm
      )
      
      df
    }
  )
)

hall_3 <- as.data.frame(
  hall_3
)

expected_nes_cols <- c(
  "NES_GSE272769",
  "NES_GSE95233",
  "NES_GSE65682"
)

expected_padj_cols <- c(
  "padj_GSE272769",
  "padj_GSE95233",
  "padj_GSE65682"
)

stopifnot(
  all(expected_nes_cols %in% colnames(hall_3)),
  all(expected_padj_cols %in% colnames(hall_3))
)

for (col in expected_padj_cols) {
  hall_3[[col]][
    is.na(hall_3[[col]])
  ] <- 1
}

nes_cols <- expected_nes_cols

# Used only to order pathways in the visualization.
hall_3$mean_abs_NES <- rowMeans(
  abs(
    hall_3[, nes_cols, drop = FALSE]
  ),
  na.rm = TRUE
)

nes_direction_matrix <- sign(
  as.matrix(
    hall_3[, nes_cols, drop = FALSE]
  )
)

hall_3$all_same_dir <- apply(
  nes_direction_matrix,
  1,
  function(x) {
    all(is.finite(x)) &&
      length(unique(x)) == 1
  }
)

hall_3$n_same_pairs <- apply(
  nes_direction_matrix,
  1,
  function(x) {
    agreement_matrix <- outer(
      x,
      x,
      "=="
    )
    
    sum(
      agreement_matrix[
        upper.tri(
          agreement_matrix
        )
      ]
    )
  }
)

hall_3$total_pairs <- 3L

cat(sprintf(
  paste0(
    "Triple-cohort Hallmark GSEA: %d pathways | ",
    "All three cohorts directionally concordant: %d (%.1f%%)\n"
  ),
  nrow(hall_3),
  sum(hall_3$all_same_dir),
  100 * mean(hall_3$all_same_dir)
))

write.csv(
  hall_3,
  file.path(
    S5_DIR,
    "GSEA_Hallmark_triple_cohort.csv"
  ),
  row.names = FALSE
)

# Display the 15 pathways with the highest mean absolute NES.
# This statistic is used only for display ordering.
if (nrow(hall_3) >= 3) {
  hall_top <- head(
    hall_3[
      order(
        -hall_3$mean_abs_NES
      ),
    ],
    min(
      15,
      nrow(hall_3)
    )
  )
  
  hm_h <- as.matrix(
    hall_top[
      ,
      nes_cols,
      drop = FALSE
    ]
  )
  
  rownames(hm_h) <- gsub(
    "^HALLMARK_",
    "",
    hall_top$pathway
  )
  
  hm_h <- hm_h[
    nrow(hm_h):1,
    ,
    drop = FALSE
  ]
  
  draw_heatmap_triple(
    hm_h,
    "GSEA Hallmark - Triple-cohort",
    "Figure_GSEA_triple_heatmap",
    S5_DIR,
    8,
    min(
      8,
      nrow(hm_h) * 0.35 + 2
    )
  )
}

# TF enrichment
cat("\n--- TF Enrichment ---\n")
httr::set_config(httr::config(http_version=1))
genes_up<-unique(na.omit(fisher_out$gene[fisher_out$direction=="Up"])); genes_down<-unique(na.omit(fisher_out$gene[fisher_out$direction=="Down"]))
TF_DBS_DEFAULT<-c("ChEA_2022","ENCODE_TF_ChIP-seq_2015","TRANSFAC_and_JASPAR_PWMs","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

run_tf_online<-function(genes,direction){
  if(length(genes)==0) return(NULL)
  res<-enrichr_retry_v2(genes,TF_DBS_DEFAULT); if(is.null(res)) return(NULL)
  for(db in names(res)){df<-safe_enrichr_df(res,db); if(nrow(df)>0) write.csv(df,file.path(S5_DIR,sprintf("TF_%s_%s.csv",db,direction)),row.names=FALSE)}
  chea<-safe_enrichr_df(res,"ChEA_2022")
  if(nrow(chea)>=3){chea$TF_symbol<-extract_tf_symbol(chea$Term); chea$Combined.Score<-suppressWarnings(as.numeric(chea$Combined.Score)); chea<-chea[!is.na(chea$TF_symbol)&chea$TF_symbol!=""&!is.na(chea$Combined.Score),]
  # 在 chea 过滤后、by() 之前加 Overlap_n
  chea$Overlap_n <- suppressWarnings(as.numeric(sub("/.*", "", chea$Overlap)))
  chea$Overlap_n[is.na(chea$Overlap_n)] <- 1
  if(nrow(chea)>=3){chea_dd<-do.call(rbind,by(chea,chea$TF_symbol,function(x)x[which.max(x$Combined.Score),,drop=FALSE])); chea_dd<-head(chea_dd[order(-chea_dd$Combined.Score),],10)
  p<-ggplot(chea_dd,aes(Combined.Score,reorder(TF_symbol,Combined.Score),size=Overlap_n,color=-log10(pmax(as.numeric(Adjusted.P.value),.Machine$double.xmin))))+
    geom_point()+scale_color_gradient(low="#4E79A7",high="#E15759")+scale_size(range=c(3,7))+
    labs(title=sprintf("TF Enrichment - ChEA 2022 (%s)",direction),x="Combined Score",y="")+theme_bw(11)
  save_plot(p,sprintf("Figure_TF_ChEA2022_%s",direction),S5_DIR,10,5)}}; res}
tf_up<-run_tf_online(genes_up,"Up"); tf_down<-run_tf_online(genes_down,"Down")

saveRDS(list(xc_272769=xc_272769,xc_95233=xc_95233,xc_65682=xc_65682,xc_all=xc_all,hall_3=hall_3),file.path(S5_DIR,"checkpoint_s5.rds"))
cat("§5 complete.\n")

# ============================================================
# §6 — Drug Repurposing 
# 产出: 06_drug/
# ============================================================
cat("\n========== §6 Drug Repurposing ==========\n")
library(rentrez)
library(digest)

drug_dbs <- c("Drug_Perturbations_from_GEO_2014")

# =============================================================
# 0. Enrichr v3（缓存 + 不 probe 主页）
# =============================================================
enrichr_cache_dir <- file.path(ROOT_DIR, "enrichr_cache")
dir.create(enrichr_cache_dir, showWarnings = FALSE, recursive = TRUE)

enrichr_retry_v3 <- function(genes, dbs, max_retry = 5, base_wait = 3,
                             use_cache = TRUE, label = NULL) {
  genes_sorted <- sort(genes)
  cache_key <- digest::digest(list(genes_sorted, dbs), algo = "md5")
  cache_file <- file.path(enrichr_cache_dir, paste0(cache_key, ".rds"))
  
  if (use_cache && file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (!is.null(cached) && length(cached) > 0 &&
        any(sapply(cached, function(x) is.data.frame(x) && nrow(x) > 0))) {
      cat(sprintf("  Enrichr CACHE HIT (%s)\n",
                  ifelse(is.null(label), paste(length(genes), "genes"), label)))
      return(cached)
    } else {
      file.remove(cache_file)
    }
  }
  
  last_error <- NULL
  for (i in seq_len(max_retry)) {
    cat(sprintf("  Enrichr attempt %d/%d (%s) ... ",
                i, max_retry, ifelse(is.null(label), paste(length(genes), "genes"), label)))
    
    res <- tryCatch(
      enrichR::enrichr(genes, dbs),
      error = function(e) { last_error <<- e$message; NULL }
    )
    
    if (!is.null(res) && length(res) > 0 &&
        any(sapply(res, function(x) is.data.frame(x) && nrow(x) > 0))) {
      cat("SUCCESS\n")
      if (use_cache) saveRDS(res, cache_file)
      return(res)
    }
    
    if (is.null(res)) {
      cat(sprintf("error (%s) ", ifelse(is.null(last_error), "unknown", substr(last_error, 1, 40))))
    } else {
      cat("empty result ")
    }
    
    if (i < max_retry) {
      wait <- base_wait * 2^(i - 1)
      cat(sprintf("| waiting %ds...\n", wait))
      Sys.sleep(wait)
    } else {
      cat("| FINAL FAILURE\n")
    }
  }
  if (!is.null(last_error)) {
    cat(sprintf("  Last error: %s\n", substr(last_error, 1, 200)))
  }
  return(NULL)
}

# =============================================================
# 1. 全签名药物扰动 (32 genes)
# =============================================================
cat("\n--- 1. Full-signature drug perturbation ---\n")
drug_up <- enrichr_retry_v3(
  genes_up,
  drug_dbs,
  max_retry = 5,
  label = "Upregulated-signature query"
)

Sys.sleep(2)

drug_down <- enrichr_retry_v3(
  genes_down,
  drug_dbs,
  max_retry = 5,
  label = "Downregulated-signature query"
)

process_drug <- function(res, direction, n_input) {
  if (is.null(res)) return(NULL)
  
  df <- safe_enrichr_df(res, drug_dbs[1])
  if (nrow(df) == 0) return(NULL)
  
  df$direction <- direction
  df$n_input_genes <- n_input
  
  # Preserve the original Term and derive both display name and matching key.
  df$drug_name <- extract_drug(df$Term)
  df$drug_key <- norm_drug(df$drug_name)
  
  df$Overlap_n <- suppressWarnings(as.numeric(sub("/.*", "", df$Overlap)))
  df$Combined.Score <- suppressWarnings(as.numeric(df$Combined.Score))
  df$Adjusted.P.value <- suppressWarnings(as.numeric(df$Adjusted.P.value))
  
  valid <- !is.na(df$drug_key) &
    nzchar(df$drug_key) &
    !is.na(df$Combined.Score) &
    !is.na(df$Adjusted.P.value)
  
  df <- df[valid, , drop = FALSE]
  df[order(-df$Combined.Score), , drop = FALSE]
}

drug_all <- rbind(
  process_drug(
    drug_up,
    "Upregulated-signature query",
    length(genes_up)
  ),
  process_drug(
    drug_down,
    "Downregulated-signature query",
    length(genes_down)
  )
)
if (!is.null(drug_all) && nrow(drug_all) > 0) {
  write.csv(drug_all, file.path(S6_DIR, "Drug_perturbation_all.csv"), row.names = FALSE)
  cat(sprintf("  Full signature: %d drugs\n", length(unique(na.omit(drug_all$drug_key)))))
}

# ---- Figure: Enrichr full-signature perturbation bar (top 25) ----
enrichr_top <- top_enrichr_by_drug(drug_all)
enrichr_top <- enrichr_top[
  !is.na(enrichr_top$drug_key) &
    enrichr_top$drug_key != "doxorubicin",
  ,
  drop = FALSE
]
enrichr_top <- head(enrichr_top[order(-enrichr_top$Combined.Score), ], 25)
enrichr_top$drug_name <- factor(enrichr_top$drug_name, levels = rev(enrichr_top$drug_name))
p_enrichr_full <- ggplot(enrichr_top, aes(Combined.Score, drug_name)) +
  geom_bar(stat = "identity", fill = "#E15759", alpha = 0.85, width = 0.7) +
  labs(title = "Enrichr: Full-Signature Perturbation (32 genes)",
       subtitle = "Top 25 displayed after excluding doxorubicin",
       x = "Combined Score", y = NULL) + theme_bw(11)
save_plot(p_enrichr_full, "Figure_Enrichr_perturbation_full", S6_DIR, 10, 10)

# =============================================================
# 2. 非周期药物扰动 (21 genes)
# =============================================================
cat("\n--- 2. Non-cycle drug perturbation ---\n")
cycle_genes_A <- MODULE_CELL_CYCLE
non_cycle_up_A <- setdiff(genes_up, cycle_genes_A)
non_cycle_dn_A <- setdiff(genes_down, cycle_genes_A)

# DGIdb query sets: 18 non-cell-cycle upregulated genes and
# three downregulated signature genes.
noncycle_upregulated_genes_dgidb <- non_cycle_up_A
downregulated_signature_genes_dgidb <- intersect(MODULE_DOWN_DEFENSE, genes_down)
stopifnot(length(noncycle_upregulated_genes_dgidb) == 18)
stopifnot(length(downregulated_signature_genes_dgidb) == 3)
arm1 <- NULL
arm2 <- NULL

drug_nc_all <- rbind(
  process_drug(
    enrichr_retry_v3(
      non_cycle_up_A,
      drug_dbs,
      max_retry = 5,
      label = "Non-cycle upregulated query"
    ),
    "Non-cycle upregulated query",
    length(non_cycle_up_A)
  ),
  process_drug(
    enrichr_retry_v3(
      non_cycle_dn_A,
      drug_dbs,
      max_retry = 5,
      label = "Downregulated-signature query"
    ),
    "Downregulated-signature query",
    length(non_cycle_dn_A)
  )
)
if (!is.null(drug_nc_all) && nrow(drug_nc_all) > 0) {
  write.csv(drug_nc_all, file.path(S6_DIR, "Drug_perturbation_non_cycle.csv"), row.names = FALSE)
  cat(sprintf("  Non-cycle: %d drugs\n", length(unique(na.omit(drug_nc_all$drug_key)))))
}

# ---- Figure: Enrichr non-cycle perturbation bar (top 25) ----
enrichr_nc_top <- top_enrichr_by_drug(drug_nc_all)
enrichr_nc_top <- enrichr_nc_top[
  !is.na(enrichr_nc_top$drug_key) &
    enrichr_nc_top$drug_key != "doxorubicin",
  ,
  drop = FALSE
]
enrichr_nc_top <- head(enrichr_nc_top[order(-enrichr_nc_top$Combined.Score), ], 25)
enrichr_nc_top$drug_name <- factor(enrichr_nc_top$drug_name, levels = rev(enrichr_nc_top$drug_name))
p_enrichr_nc <- ggplot(enrichr_nc_top, aes(Combined.Score, drug_name)) +
  geom_bar(stat = "identity", fill = "#B07AA1", alpha = 0.85, width = 0.7) +
  labs(title = "Enrichr: Non-Cycle Perturbation (21 genes)",
       subtitle = "Top 25 displayed after excluding doxorubicin",
       x = "Combined Score", y = NULL) + theme_bw(11)
save_plot(p_enrichr_nc, "Figure_Enrichr_perturbation_non_cycle", S6_DIR, 10, 10)

# =============================================================
# 3. Dual-query Enrichr perturbational overlap
# =============================================================
cat("\n--- 3. Dual-query Enrichr perturbational overlap ---\n")
downregulated_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

gene_list <- lapply(
  downregulated_pathways,
  function(pw) {
    idx <- which(fg_h_272769$pathway == pw)
    if (length(idx) == 0) return(NULL)
    
    genes <- fg_h_272769$leadingEdge[[idx[1]]]
    if (is.null(genes) || length(genes) == 0) return(NULL)
    
    genes
  }
)

gene_list <- gene_list[!vapply(gene_list, is.null, logical(1))]
gene_freq <- table(unlist(gene_list))

downregulated_pathway_genes <- head(
  names(sort(gene_freq, decreasing = TRUE)),
  150
)

drug_downregulated_pathway <- enrichr_retry_v3(
  downregulated_pathway_genes,
  drug_dbs,
  max_retry = 5,
  label = "Downregulated-pathway query (150 genes)"
)

drug_noncycle_upregulated <- enrichr_retry_v3(
  non_cycle_up_A,
  drug_dbs,
  max_retry = 5,
  label = "Non-cycle upregulated query (18 genes)"
)
bi_summary <- NULL
if (!is.null(drug_downregulated_pathway) && !is.null(drug_noncycle_upregulated)) {
  df_downregulated_pathway <- safe_enrichr_df(drug_downregulated_pathway, drug_dbs[1])
  df_noncycle_upregulated <- safe_enrichr_df(drug_noncycle_upregulated, drug_dbs[1])
  df_downregulated_pathway$drug_name <- extract_drug(df_downregulated_pathway$Term)
  df_noncycle_upregulated$drug_name <- extract_drug(df_noncycle_upregulated$Term)
  
  df_downregulated_pathway$drug_key <- norm_drug(df_downregulated_pathway$drug_name)
  df_noncycle_upregulated$drug_key <- norm_drug(df_noncycle_upregulated$drug_name)
  
  df_downregulated_pathway$Combined.Score <- suppressWarnings(
    as.numeric(df_downregulated_pathway$Combined.Score)
  )
  df_noncycle_upregulated$Combined.Score <- suppressWarnings(
    as.numeric(df_noncycle_upregulated$Combined.Score)
  )
  
  df_downregulated_pathway <- df_downregulated_pathway[
    !is.na(df_downregulated_pathway$drug_key) &
      nzchar(df_downregulated_pathway$drug_key) &
      !is.na(df_downregulated_pathway$Combined.Score),
    ,
    drop = FALSE
  ]
  
  df_noncycle_upregulated <- df_noncycle_upregulated[
    !is.na(df_noncycle_upregulated$drug_key) &
      nzchar(df_noncycle_upregulated$drug_key) &
      !is.na(df_noncycle_upregulated$Combined.Score),
    ,
    drop = FALSE
  ]
  
  immune_drug_keys <- unique(df_downregulated_pathway$drug_key)
  neutro_drug_keys <- unique(df_noncycle_upregulated$drug_key)
  bi_drug_keys <- intersect(immune_drug_keys, neutro_drug_keys)
  
  cat(sprintf(
    "  Cross-query compound overlap: %d\n",
    length(bi_drug_keys)
  ))
  
  if (length(bi_drug_keys) > 0) {
    
    representative_name <- vapply(
      bi_drug_keys,
      function(k) {
        candidates <- rbind(
          data.frame(
            drug_name = df_downregulated_pathway$drug_name[df_downregulated_pathway$drug_key == k],
            score = df_downregulated_pathway$Combined.Score[df_downregulated_pathway$drug_key == k],
            stringsAsFactors = FALSE
          ),
          data.frame(
            drug_name = df_noncycle_upregulated$drug_name[df_noncycle_upregulated$drug_key == k],
            score = df_noncycle_upregulated$Combined.Score[df_noncycle_upregulated$drug_key == k],
            stringsAsFactors = FALSE
          )
        )
        
        candidates <- candidates[
          order(-candidates$score, candidates$drug_name),
          ,
          drop = FALSE
        ]
        
        candidates$drug_name[1]
      },
      character(1)
    )
    
    bi_summary <- data.frame(
      drug_key = bi_drug_keys,
      drug = representative_name,
      downregulated_pathway_score = vapply(
        bi_drug_keys,
        function(k) {
          max(
            df_downregulated_pathway$Combined.Score[df_downregulated_pathway$drug_key == k],
            na.rm = TRUE
          )
        },
        numeric(1)
      ),
      noncycle_upregulated_score = vapply(
        bi_drug_keys,
        function(k) {
          max(
            df_noncycle_upregulated$Combined.Score[df_noncycle_upregulated$drug_key == k],
            na.rm = TRUE
          )
        },
        numeric(1)
      ),
      stringsAsFactors = FALSE
    )
    
    bi_summary$bidirectional_score <- sqrt(
      bi_summary$downregulated_pathway_score *
        bi_summary$noncycle_upregulated_score
    )
    
    bi_summary$pass_threshold <-
      bi_summary$downregulated_pathway_score > 30 &
      bi_summary$noncycle_upregulated_score > 30
    
    bi_summary <- bi_summary[
      order(-bi_summary$bidirectional_score),
      ,
      drop = FALSE
    ]
    
    write.csv(
      bi_summary,
      file.path(S6_DIR, "Drug_bidirectional_candidates.csv"),
      row.names = FALSE
    )
    
    label_df <- subset(bi_summary, pass_threshold)
    
    p_bi <- ggplot(
      bi_summary,
      aes(
        x = downregulated_pathway_score,
        y = noncycle_upregulated_score,
        size = bidirectional_score
      )
    ) +
      geom_point(alpha = 0.55, color = "#E15759") +
      geom_vline(xintercept = 30, lty = "dashed", linewidth = 0.3, color = "grey55") +
      geom_hline(yintercept = 30, lty = "dashed", linewidth = 0.3, color = "grey55") +
      geom_text_repel(
        data = label_df,
        aes(
          x = downregulated_pathway_score,
          y = noncycle_upregulated_score,
          label = drug
        ),
        inherit.aes = FALSE,
        size = 2.7,
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.18,
        min.segment.length = 0,
        segment.size = 0.25,
        segment.color = "grey55"
      ) +
      scale_x_continuous(
        trans = scales::pseudo_log_trans(sigma = 10, base = 10),
        limits = c(0, 1250),
        breaks = c(0, 30, 250, 500, 750, 1000, 1250),
        minor_breaks = NULL
      ) +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 10, base = 10),
        limits = c(0, 20000),
        breaks = c(0, 5000, 10000, 15000, 20000),
        minor_breaks = NULL
      ) +
      scale_size_continuous(
        name = "Geometric-mean score",
        range = c(1.1, 4.2),
        breaks = c(100, 200, 300, 400, 500)
      ) +
      labs(
        title = "Enrichr: Dual-Query Perturbational Matches",
        subtitle = sprintf(
          "%d compound matches | %d had CS > 30 in both query arms",
          nrow(bi_summary),
          sum(bi_summary$pass_threshold)
        ),
        x = "Downregulated-pathway arm score",
        y = "Non-cycle upregulated arm score"
      ) +
      theme_bw(9) +
      theme(
        plot.title = element_text(size = 10.5, face = "plain"),
        plot.subtitle = element_text(size = 8.8),
        axis.title = element_text(size = 8.8),
        axis.text = element_text(size = 7.8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
        plot.margin = margin(6, 18, 8, 6)
      )
    
    save_plot(p_bi, "Figure_Enrichr_bidirectional", S6_DIR, 8, 6)
  }
}

# =============================================================
# 4. DGIdb v5
# =============================================================
cat("\n--- 4. DGIdb v5 ---\n")
query_dgidb_v5 <- function(genes, label, retry = 3) {
  gene_str <- paste0('"', genes, '"', collapse = ",")
  query <- sprintf('{genes(names:[%s]){nodes{name interactions{drug{name approved}interactionTypes{type}interactionScore}}}}', gene_str)
  for (i in 1:retry) {
    resp <- tryCatch(httr::POST("https://dgidb.org/api/graphql",
                                body = toJSON(list(query = query), auto_unbox = TRUE), encode = "raw",
                                httr::add_headers("Content-Type" = "application/json"), httr::timeout(60)), error = function(e) NULL)
    if (is.null(resp)) { Sys.sleep(3); next }
    if (httr::status_code(resp) != 200) { Sys.sleep(5); next }
    parsed <- tryCatch(fromJSON(httr::content(resp, "text", encoding = "UTF-8")), error = function(e) NULL)
    if (is.null(parsed) || !is.null(parsed$errors)) break
    gd <- parsed$data$genes$nodes; if (is.null(gd)) gd <- parsed$data$genes
    if (is.null(gd) || length(gd) == 0) break
    row_list <- list()
    n_genes <- if (is.data.frame(gd)) nrow(gd) else length(gd)
    for (j in seq_len(n_genes)) {
      gene_name <- tryCatch(as.character(gd$name[j]), error = function(e) NA_character_)
      ints <- tryCatch(gd$interactions[[j]], error = function(e) NULL)
      if (is.null(ints) || length(ints) == 0 || !is.data.frame(ints)) next
      n <- nrow(ints)
      drug_names <- if (!is.null(ints$drug) && !is.null(ints$drug$name)) ints$drug$name else rep(NA_character_, n)
      approved <- if (!is.null(ints$drug) && !is.null(ints$drug$approved)) ints$drug$approved else rep(NA, n)
      scores <- if (!is.null(ints$interactionScore)) ints$interactionScore else rep(NA_real_, n)
      types <- tryCatch({sapply(seq_len(n), function(k) {
        t <- ints$interactionTypes[[k]]
        if (is.data.frame(t) && nrow(t) > 0) paste(t$type, collapse = "; ") else NA_character_
      })}, error = function(e) rep(NA_character_, n))
      row_list[[length(row_list) + 1]] <- data.frame(
        gene = gene_name, drug = as.character(drug_names), approved = as.logical(approved),
        interaction_type = as.character(types), score = as.numeric(scores), stringsAsFactors = FALSE)
    }
    results <- if (length(row_list) > 0) do.call(rbind, row_list) else NULL
    if (!is.null(results) && nrow(results) > 0) {
      cat(sprintf("  %s: %d interactions, %d unique drugs\n", label, nrow(results), length(unique(results$drug))))
      return(results)
    }
    break
  }
  return(NULL)
}

# ---- 4a. Single-gene DGIdb ----
dgidb_all <- query_dgidb_v5(c(genes_up, genes_down), "32_genes")
dgidb_query_time_utc <- format(
  Sys.time(),
  tz = "UTC",
  usetz = TRUE
)

if (!is.null(dgidb_all) && nrow(dgidb_all) > 0) {
  write.csv(
    dgidb_all,
    file.path(S6_DIR, "DGIdb_32genes_interactions_raw.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      database = "DGIdb",
      api_endpoint = "https://dgidb.org/api/graphql",
      retrieval_time_utc = dgidb_query_time_utc,
      n_query_genes = length(unique(c(genes_up, genes_down))),
      n_raw_interaction_records = nrow(dgidb_all),
      n_raw_unique_drugs = length(unique(dgidb_all$drug)),
      stringsAsFactors = FALSE
    ),
    file.path(S6_DIR, "DGIdb_query_provenance.csv"),
    row.names = FALSE
  )
}
gene_summary <- NULL
if (!is.null(dgidb_all) && nrow(dgidb_all) > 0) {
  dgidb_all$drug_name <- trimws(as.character(dgidb_all$drug))
  dgidb_all$drug_key <- norm_drug(dgidb_all$drug_name)
  
  dgidb_all <- dgidb_all[
    !is.na(dgidb_all$drug_key) & nzchar(dgidb_all$drug_key),
    ,
    drop = FALSE
  ]
  
  # Deduplicate at the gene-compound-key level.
  dgidb_all <- dgidb_all[
    !duplicated(dgidb_all[, c("gene", "drug_key")]),
    ,
    drop = FALSE
  ]
  write.csv(
    dgidb_all,
    file.path(S6_DIR, "DGIdb_32genes_interactions_deduplicated.csv"),
    row.names = FALSE
  )
  arm1 <- dgidb_all[dgidb_all$gene %in% noncycle_upregulated_genes_dgidb, , drop = FALSE]
  arm2 <- dgidb_all[dgidb_all$gene %in% downregulated_signature_genes_dgidb, , drop = FALSE]
  cat(sprintf("  DGIdb arms: neutrophil=%d genes with records | immune=%d genes with records\n",
              length(unique(arm1$gene)), length(unique(arm2$gene))))
  gene_summary <- aggregate(drug ~ gene, data = dgidb_all,
                            FUN = function(x) paste(sort(unique(x)), collapse = "; "))
  gene_summary$n_drugs <- lengths(strsplit(gene_summary$drug, "; "))
  gene_summary <- gene_summary[order(-gene_summary$n_drugs), ]
  write.csv(gene_summary, file.path(S6_DIR, "DGIdb_32genes_summary.csv"), row.names = FALSE)
  cat(sprintf("  DGIdb single-gene: %d/%d genes with records\n", length(unique(dgidb_all$gene)), 32))
}

# ---- Figure: DGIdb per-gene drug count ----
if (!is.null(gene_summary) && nrow(gene_summary) > 0) {
  gene_summary$gene <- factor(gene_summary$gene, levels = rev(gene_summary$gene))
  gene_mod <- ifelse(as.character(gene_summary$gene) %in% names(module_map),
                     as.character(module_map[as.character(gene_summary$gene)]), "Unclassified")
  gene_summary$module <- factor(gene_mod, levels = MODULE_LEVELS)
  mod_fill <- c("Neutrophil degranulation" = "#E15759", "Cell cycle / proliferation" = "#4E79A7",
                "Inflammatory / Down" = "#59A14F", "Other" = "#B07AA1", "Unclassified" = "grey70")
  p_dgidb_gene <- ggplot(gene_summary, aes(n_drugs, gene, fill = module)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.7) + scale_fill_manual(values = mod_fill) +
    labs(title = "DGIdb: Known Drug Interactions per Signature Gene",
         subtitle = sprintf("%d/%d genes with records | %d unique drugs", nrow(gene_summary), 32, length(unique(dgidb_all$drug))),
         x = "Number of interacting drugs", y = NULL, fill = "Module") + theme_bw(11) + theme(legend.position = "bottom")
  save_plot(p_dgidb_gene, "Figure_DGIdb_per_gene", S6_DIR, 9, nrow(gene_summary) * 0.35 + 2)
}

# ---- 4b. DGIdb cross-query interaction overlap ----
bi_df <- NULL
if (!is.null(arm1) && !is.null(arm2) && nrow(arm1) > 0 && nrow(arm2) > 0) {
  arm1$drug_key <- norm_drug(arm1$drug_name)
  arm2$drug_key <- norm_drug(arm2$drug_name)
  
  bi_keys <- intersect(
    unique(na.omit(arm1$drug_key)),
    unique(na.omit(arm2$drug_key))
  )
  
  cat(sprintf(
    "  DGIdb cross-query overlap: %d\n",
    length(bi_keys)
  ))
  
  if (length(bi_keys) > 0) {
    
    representative_dgidb_name <- vapply(
      bi_keys,
      function(k) {
        values <- c(
          arm1$drug_name[arm1$drug_key == k],
          arm2$drug_name[arm2$drug_key == k]
        )
        values <- unique(values[!is.na(values) & nzchar(values)])
        values[1]
      },
      character(1)
    )
    
    bi_df <- data.frame(
      drug_key = bi_keys,
      drug = representative_dgidb_name,
      stringsAsFactors = FALSE
    )
    
    bi_df$noncycle_upregulated_genes <- vapply(
      bi_keys,
      function(k) {
        paste(
          sort(unique(arm1$gene[arm1$drug_key == k])),
          collapse = "; "
        )
      },
      character(1)
    )
    
    bi_df$downregulated_signature_genes <- vapply(
      bi_keys,
      function(k) {
        paste(
          sort(unique(arm2$gene[arm2$drug_key == k])),
          collapse = "; "
        )
      },
      character(1)
    )
    
    bi_df$n_noncycle_upregulated <- lengths(
      strsplit(bi_df$noncycle_upregulated_genes, "; ", fixed = TRUE)
    )
    
    bi_df$n_downregulated_signature <- lengths(
      strsplit(bi_df$downregulated_signature_genes, "; ", fixed = TRUE)
    )
    
    bi_df$combined_n <-
      bi_df$n_noncycle_upregulated +
      bi_df$n_downregulated_signature
    
    bi_df$approved <- vapply(
      bi_keys,
      function(k) {
        any(arm1$approved[arm1$drug_key == k] %in% TRUE) |
          any(arm2$approved[arm2$drug_key == k] %in% TRUE)
      },
      logical(1)
    )
    
    bi_df$clean_name <- format_drug_label(bi_df$drug)
    bi_df <- bi_df[order(-bi_df$approved, -bi_df$combined_n), ]
    write.csv(bi_df, file.path(S6_DIR, "DGIdb_bidirectional_candidates.csv"), row.names = FALSE)
    
    if (nrow(bi_df) >= 1) {
      bi_df$drug_label <- ifelse(bi_df$approved, paste0(bi_df$clean_name, " *"), bi_df$clean_name)
      bi_df$drug_label <- factor(bi_df$drug_label, levels = rev(bi_df$drug_label))
      p_dgidb_bi <- ggplot(bi_df, aes(x = drug_label)) +
        geom_col(aes(y = n_noncycle_upregulated, fill = "Non-cycle upregulated targets"), alpha = 0.85, width = 0.5) +
        geom_col(aes(y = -n_downregulated_signature, fill = "Downregulated-signature targets"), alpha = 0.85, width = 0.5) +
        geom_hline(yintercept = 0, linewidth = 0.5) +
        scale_fill_manual(values = c(
          "Non-cycle upregulated targets" = "#E15759",
          "Downregulated-signature targets" = "#4E79A7"
        )) +
        scale_y_continuous(labels = abs, breaks = seq(-5, 5, 1)) + coord_flip() +
        labs(title = "DGIdb: Curated Interaction Overlap Across Query Sets",
             subtitle = sprintf(   "%d compound records with curated interactions in both query sets",   nrow(bi_df) ), x = "", y = "Target genes") +
        theme_bw(11) + theme(legend.position = "bottom")
      save_plot(p_dgidb_bi, "Figure_DGIdb_bidirectional_bar", S6_DIR, 7, nrow(bi_df) * 0.8 + 2)
    }
  }
} else {
  cat("  DGIdb cross-query overlap skipped: one or both arms had no DGIdb records.\n")
}

# ---- Figure: DGIdb records by query gene set ----
if (!is.null(dgidb_all) && nrow(dgidb_all) > 0) {
  
  noncycle_upregulated_dgidb_genes <- intersect(
    noncycle_upregulated_genes_dgidb,
    unique(dgidb_all$gene)
  )
  
  downregulated_signature_dgidb_genes <- intersect(
    downregulated_signature_genes_dgidb,
    unique(dgidb_all$gene)
  )
  
  dgidb_query_genes <- unique(
    c(
      noncycle_upregulated_dgidb_genes,
      downregulated_signature_dgidb_genes
    )
  )
  
  if (length(dgidb_query_genes) > 0) {
    dgidb_nc <- dgidb_all[
      dgidb_all$gene %in% dgidb_query_genes,
      ,
      drop = FALSE
    ]
    
    dgidb_nc_summary <- aggregate(
      drug ~ gene,
      data = dgidb_nc,
      FUN = function(x) paste(sort(unique(x)), collapse = "; ")
    )
    
    dgidb_nc_summary$n_drugs <- lengths(
      strsplit(dgidb_nc_summary$drug, "; ", fixed = TRUE)
    )
    
    dgidb_nc_summary <- dgidb_nc_summary[
      order(-dgidb_nc_summary$n_drugs),
      ,
      drop = FALSE
    ]
    
    dgidb_nc_summary$query_set <- ifelse(
      dgidb_nc_summary$gene %in% noncycle_upregulated_dgidb_genes,
      "Non-cycle upregulated",
      "Downregulated signature"
    )
    
    dgidb_nc_summary$gene <- factor(
      dgidb_nc_summary$gene,
      levels = rev(dgidb_nc_summary$gene)
    )
    
    p_dgidb_nc <- ggplot(
      dgidb_nc_summary,
      aes(n_drugs, gene, fill = query_set)
    ) +
      geom_col(
        alpha = 0.85,
        width = 0.7
      ) +
      scale_fill_manual(
        values = c(
          "Non-cycle upregulated" = "#E15759",
          "Downregulated signature" = "#4E79A7"
        )
      ) +
      labs(
        title = "DGIdb: Curated Interaction Records by Query Gene Set",
        subtitle = sprintf(
          "%d non-cycle upregulated genes | %d downregulated signature genes",
          length(noncycle_upregulated_dgidb_genes),
          length(downregulated_signature_dgidb_genes)
        ),
        x = "Number of compound records",
        y = NULL,
        fill = "Query set"
      ) +
      theme_bw(11) +
      theme(legend.position = "bottom")
    
    # Retain the existing basename for downstream compatibility.
    save_plot(
      p_dgidb_nc,
      "Figure_DGIdb_per_gene_non_cycle",
      S6_DIR,
      8,
      nrow(dgidb_nc_summary) * 0.4 + 2
    )
  }
}
# =============================================================
# 5. Descriptive DGIdb compound-level interaction summary
# =============================================================
cat("\n--- 5. Descriptive DGIdb compound-level summary ---\n")

dgidb_per_drug <- NULL

if (!is.null(dgidb_all) && nrow(dgidb_all) > 0) {
  
  required_columns <- c(
    "gene",
    "drug_name",
    "drug_key",
    "approved"
  )
  
  if (!all(required_columns %in% colnames(dgidb_all))) {
    stop(
      "DGIdb compound-level summary requires gene, drug_name, ",
      "drug_key, and approved columns."
    )
  }
  
  # Count unique signature genes with curated interaction records
  # for each standardized compound key.
  dgidb_per_drug <- aggregate(
    gene ~ drug_key,
    data = dgidb_all,
    FUN = function(x) length(unique(x))
  )
  
  colnames(dgidb_per_drug)[
    colnames(dgidb_per_drug) == "gene"
  ] <- "n_signature_genes"
  
  display_lookup <- tapply(
    dgidb_all$drug_name,
    dgidb_all$drug_key,
    function(x) {
      x <- unique(x[!is.na(x) & nzchar(x)])
      if (length(x) == 0) return(NA_character_)
      format_drug_label(x[1])
    }
  )
  
  approval_lookup <- tapply(
    dgidb_all$approved,
    dgidb_all$drug_key,
    function(x) any(x %in% TRUE)
  )
  
  dgidb_per_drug$drug_name <- unname(
    display_lookup[dgidb_per_drug$drug_key]
  )
  
  dgidb_per_drug$approved <- unname(
    approval_lookup[dgidb_per_drug$drug_key]
  )
  
  dgidb_per_drug <- dgidb_per_drug[
    order(
      -dgidb_per_drug$n_signature_genes,
      dgidb_per_drug$drug_key
    ),
    ,
    drop = FALSE
  ]
  
  write.csv(
    dgidb_per_drug,
    file.path(S6_DIR, "DGIdb_per_compound_summary.csv"),
    row.names = FALSE
  )
  
  # Descriptive display of compounds linked to at least two
  # signature genes.
  dgidb_multi <- dgidb_per_drug[
    dgidb_per_drug$n_signature_genes >= 2,
    ,
    drop = FALSE
  ]
  
  if (nrow(dgidb_multi) > 0) {
    n_show <- min(20L, nrow(dgidb_multi))
    
    dgidb_multi_plot <- head(
      dgidb_multi,
      n_show
    )
    
    # Standardize display labels only; do not alter records or ranking.
    dgidb_multi_plot$display_name <- as.character(
      dgidb_multi_plot$drug_name
    )
    
    dgidb_multi_plot$display_name[
      tolower(dgidb_multi_plot$drug_key) == "drugsatfda.nda:008378"
    ] <- "Prednisone (NDA 008378)"
    
    dgidb_multi_plot$display_name[
      tolower(dgidb_multi_plot$drug_key) == "lithium cation"
    ] <- "Lithium"
    
    dgidb_multi_plot$display_name <- gsub(
      "\\[pmid:\\s*([0-9]+)\\]",
      "(PMID \\1)",
      dgidb_multi_plot$display_name,
      ignore.case = TRUE
    )
    
    dgidb_multi_plot$display_name <- factor(
      dgidb_multi_plot$display_name,
      levels = rev(dgidb_multi_plot$display_name)
    )
    
    subtitle_text <- if (nrow(dgidb_multi) <= 20L) {
      sprintf(
        "All %d DGIdb entities linked to at least two signature genes",
        nrow(dgidb_multi)
      )
    } else {
      sprintf(
        "Top %d of %d DGIdb entities linked to at least two signature genes",
        n_show,
        nrow(dgidb_multi)
      )
    }
    
    p_dgidb_multi <- ggplot(
      dgidb_multi_plot,
      aes(n_signature_genes, display_name)
    ) +
      geom_col(
        fill = "#4E79A7",
        alpha = 0.85,
        width = 0.7
      ) +
      labs(
        title = "DGIdb: Curated Multi-gene Interaction Records",
        subtitle = subtitle_text,
        x = "Number of linked signature genes",
        y = NULL
      ) +
      theme_bw(11)
    
    # Retain the existing file basename for downstream compatibility.
    save_plot(
      p_dgidb_multi,
      "Figure_DGIdb_multitarget",
      S6_DIR,
      8,
      nrow(dgidb_multi_plot) * 0.35 + 1.5
    )
  }
}

# ============================================================
# Descriptive Enrichr Stage 1 and Stage 2 top-20 panel
# ============================================================
suppressPackageStartupMessages(library(ggplot2))

drug_dir <- S6_DIR

# ---- Stage 1: Full signature ----
enrichr_full_df <- read.csv(
  file.path(drug_dir, "Drug_perturbation_all.csv"),
  stringsAsFactors = FALSE
)

enrichr_full_df$drug_label <- if ("Term" %in% colnames(enrichr_full_df)) {
  extract_drug(enrichr_full_df$Term)
} else {
  extract_drug(enrichr_full_df$drug_name)
}

enrichr_full_df$drug_name <- enrichr_full_df$drug_label
enrichr_full_df <- top_enrichr_by_drug(enrichr_full_df)

# top_enrichr_by_drug() has created the standardized drug_key column.
enrichr_full_df <- enrichr_full_df[
  !is.na(enrichr_full_df$drug_key) &
    enrichr_full_df$drug_key != "doxorubicin",
  ,
  drop = FALSE
]

enrichr_full_top <- head(enrichr_full_df, 20)
enrichr_full_top <- enrichr_full_top[order(enrichr_full_top$Combined.Score), ]
enrichr_full_top$drug_label <- factor(enrichr_full_top$drug_label, levels = enrichr_full_top$drug_label)

fig6A_full_signature_plot <- ggplot(enrichr_full_top, aes(Combined.Score, drug_label)) +
  geom_col(fill = "#4E79A7", width = 0.7) +
  geom_text(aes(label = round(Combined.Score)), hjust = -0.2, size = 3) +
  labs(title = "Stage 1: Full signature (32 genes)", x = "Combined Score", y = "") +
  theme_bw(11) +
  theme(panel.grid.major.y = element_blank())

# ---- Stage 2: Non-cycle ----
enrichr_noncycle_df <- read.csv(
  file.path(drug_dir, "Drug_perturbation_non_cycle.csv"),
  stringsAsFactors = FALSE
)

enrichr_noncycle_df$drug_label <- if ("Term" %in% colnames(enrichr_noncycle_df)) {
  extract_drug(enrichr_noncycle_df$Term)
} else {
  extract_drug(enrichr_noncycle_df$drug_name)
}

enrichr_noncycle_df$drug_name <- enrichr_noncycle_df$drug_label
enrichr_noncycle_df <- top_enrichr_by_drug(enrichr_noncycle_df)

# top_enrichr_by_drug() has created the standardized drug_key column.
enrichr_noncycle_df <- enrichr_noncycle_df[
  !is.na(enrichr_noncycle_df$drug_key) &
    enrichr_noncycle_df$drug_key != "doxorubicin",
  ,
  drop = FALSE
]

enrichr_noncycle_top <- head(enrichr_noncycle_df, 20)
enrichr_noncycle_top <- enrichr_noncycle_top[order(enrichr_noncycle_top$Combined.Score), ]
enrichr_noncycle_top$drug_label <- factor(enrichr_noncycle_top$drug_label, levels = enrichr_noncycle_top$drug_label)

fig6B_noncycle_plot <- ggplot(enrichr_noncycle_top, aes(Combined.Score, drug_label)) +
  geom_col(fill = "#E15759", width = 0.7) +
  geom_text(aes(label = round(Combined.Score)), hjust = -0.2, size = 3) +
  labs(title = "Stage 2: Non-cycle (21 genes)", x = "Combined Score", y = "") +
  theme_bw(11) +
  theme(panel.grid.major.y = element_blank())

# ---- Combine & save ----
pdf(file.path(drug_dir, "Figure_Enrichr_stage1_stage2.pdf"), 14, 6)
gridExtra::grid.arrange(fig6A_full_signature_plot, fig6B_noncycle_plot, nrow = 1)
dev.off()

png(file.path(drug_dir, "Figure_Enrichr_stage1_stage2.png"), 14, 6, units = "in", res = 300)
gridExtra::grid.arrange(fig6A_full_signature_plot, fig6B_noncycle_plot, nrow = 1)
dev.off()

cat("Enrichr Stage 1/Stage 2 panel written.\n")
cat(sprintf("Stage 1: %d unique drugs, top 3 = %s\n",
            nrow(enrichr_full_df),
            paste(head(enrichr_full_df$drug_label, 3), collapse = ", ")))
cat(sprintf("Stage 2: %d unique drugs, top 3 = %s\n",
            nrow(enrichr_noncycle_df),
            paste(head(enrichr_noncycle_df$drug_label, 3), collapse = ", ")))
# ============================================================
# Figure 7D — Compound ranks across three exploratory
# Enrichr analyses
# ============================================================
cat("\n--- Three-analysis Enrichr rank comparison ---\n")

TRAJECTORY_TOP_N <- 20L

tri_stage_keys <- character(0)
traj_df <- NULL
three_stage_overlap_n <- NA_integer_

if (
  !is.null(drug_all) &&
  nrow(drug_all) > 0 &&
  !is.null(drug_nc_all) &&
  nrow(drug_nc_all) > 0 &&
  !is.null(bi_summary) &&
  nrow(bi_summary) > 0
) {
  
  # Stage 1: top 20 standardized compounds from the
  # full-signature query.
  #
  # Doxorubicin is retained here because its exclusion from the
  # Stage 1/2 bar plots was a display-only decision. Retaining it
  # reproduces the existing Table S17 ranking.
  trajectory_stage1 <- top_enrichr_by_drug(drug_all)
  trajectory_stage1 <- trajectory_stage1[
    order(
      -trajectory_stage1$Combined.Score,
      trajectory_stage1$drug_key
    ),
    ,
    drop = FALSE
  ]
  trajectory_stage1 <- head(
    trajectory_stage1,
    TRAJECTORY_TOP_N
  )
  trajectory_stage1$rank_in_analysis <- seq_len(
    nrow(trajectory_stage1)
  )
  
  # Stage 2: top 20 standardized compounds from the
  # 21-gene non-cycle query.
  trajectory_stage2 <- top_enrichr_by_drug(drug_nc_all)
  trajectory_stage2 <- trajectory_stage2[
    order(
      -trajectory_stage2$Combined.Score,
      trajectory_stage2$drug_key
    ),
    ,
    drop = FALSE
  ]
  trajectory_stage2 <- head(
    trajectory_stage2,
    TRAJECTORY_TOP_N
  )
  trajectory_stage2$rank_in_analysis <- seq_len(
    nrow(trajectory_stage2)
  )
  
  # Stage 3: all cross-query compound overlaps ranked by the
  # custom geometric-mean score.
  trajectory_stage3 <- bi_summary[
    !is.na(bi_summary$drug_key) &
      nzchar(bi_summary$drug_key) &
      is.finite(bi_summary$bidirectional_score),
    ,
    drop = FALSE
  ]
  
  trajectory_stage3 <- trajectory_stage3[
    order(
      -trajectory_stage3$bidirectional_score,
      trajectory_stage3$drug_key
    ),
    ,
    drop = FALSE
  ]
  
  trajectory_stage3 <- trajectory_stage3[
    !duplicated(trajectory_stage3$drug_key),
    ,
    drop = FALSE
  ]
  
  trajectory_stage3$rank_in_analysis <- seq_len(
    nrow(trajectory_stage3)
  )
  
  # Standardized compounds represented in all three analyses.
  tri_stage_keys <- Reduce(
    intersect,
    list(
      trajectory_stage1$drug_key,
      trajectory_stage2$drug_key,
      trajectory_stage3$drug_key
    )
  )
  
  three_stage_overlap_n <- length(tri_stage_keys)
  
  cat(sprintf(
    "  Compounds shared by Stage 1 top 20, Stage 2 top 20, and Stage 3: %d\n",
    three_stage_overlap_n
  ))
  
  if (three_stage_overlap_n > 0) {
    
    traj_df <- data.frame(
      drug = tri_stage_keys,
      Stage1_Full = trajectory_stage1$rank_in_analysis[
        match(
          tri_stage_keys,
          trajectory_stage1$drug_key
        )
      ],
      Stage2_NonCycle = trajectory_stage2$rank_in_analysis[
        match(
          tri_stage_keys,
          trajectory_stage2$drug_key
        )
      ],
      Stage3_Bidirectional = trajectory_stage3$rank_in_analysis[
        match(
          tri_stage_keys,
          trajectory_stage3$drug_key
        )
      ],
      Bidirectional_Score =
        trajectory_stage3$bidirectional_score[
          match(
            tri_stage_keys,
            trajectory_stage3$drug_key
          )
        ],
      stringsAsFactors = FALSE
    )
    
    traj_df <- traj_df[
      order(
        -traj_df$Bidirectional_Score,
        traj_df$drug
      ),
      ,
      drop = FALSE
    ]
    rownames(traj_df) <- NULL
    
    write.csv(
      traj_df,
      file.path(
        S6_DIR,
        "Drug_trajectory_3stage.csv"
      ),
      row.names = FALSE
    )
    
    traj_long <- reshape2::melt(
      traj_df,
      id.vars = c(
        "drug",
        "Bidirectional_Score"
      ),
      measure.vars = c(
        "Stage1_Full",
        "Stage2_NonCycle",
        "Stage3_Bidirectional"
      ),
      variable.name = "Analysis",
      value.name = "Rank"
    )
    
    traj_long$Analysis <- factor(
      traj_long$Analysis,
      levels = c(
        "Stage1_Full",
        "Stage2_NonCycle",
        "Stage3_Bidirectional"
      ),
      labels = c(
        "Full signature\n(top 20)",
        "Non-cycle set\n(top 20)",
        "Dual-query\noverlap"
      )
    )
    
    trajectory_labels <- traj_long[
      traj_long$Analysis == "Dual-query\noverlap",
      ,
      drop = FALSE
    ]
    
    p_traj <- ggplot(
      traj_long,
      aes(
        x = Analysis,
        y = Rank,
        group = drug,
        color = Bidirectional_Score
      )
    ) +
      geom_line(
        aes(linewidth = Bidirectional_Score),
        alpha = 0.65
      ) +
      geom_point(size = 2.4) +
      scale_color_gradient(
        low = "#4E79A7",
        high = "#E15759",
        trans = "log10",
        name = "Geometric-mean\nscore"
      ) +
      scale_linewidth(
        range = c(0.3, 1.4),
        guide = "none"
      ) +
      scale_y_reverse() +
      geom_text_repel(
        data = trajectory_labels,
        aes(
          label = format_drug_label(drug)
        ),
        seed = 2024,
        size = 3,
        direction = "y",
        hjust = -0.1,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        show.legend = FALSE
      ) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste(
          "Compound Ranks Across Three",
          "Exploratory Enrichr Analyses"
        ),
        subtitle = sprintf(
          "%d compounds represented in all three analyses",
          nrow(traj_df)
        ),
        x = NULL,
        y = "Rank within analysis"
      ) +
      theme_bw(11) +
      theme(
        legend.position = "right",
        plot.margin = margin(
          t = 6,
          r = 55,
          b = 6,
          l = 6
        )
      )
    
    save_plot(
      p_traj,
      "Figure_Drug_trajectory_Enrichr",
      S6_DIR,
      10,
      7
    )
    
  } else {
    warning(
      paste(
        "No standardized compounds were shared by",
        "the Stage 1 top 20, Stage 2 top 20,",
        "and Stage 3 dual-query results."
      )
    )
  }
}
# ============================================================
# Drug-analysis quality-control and provenance summary
# ============================================================

safe_unique_n <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_integer_)
  }
  
  length(unique(na.omit(x)))
}

safe_nrow <- function(x) {
  if (is.null(x)) {
    return(NA_integer_)
  }
  
  nrow(x)
}

stage1_keys <- if (!is.null(drug_all) && nrow(drug_all) > 0) {
  unique(na.omit(drug_all$drug_key))
} else {
  character(0)
}

stage2_keys <- if (!is.null(drug_nc_all) && nrow(drug_nc_all) > 0) {
  unique(na.omit(drug_nc_all$drug_key))
} else {
  character(0)
}

stage3_keys <- if (!is.null(bi_summary) && nrow(bi_summary) > 0) {
  unique(na.omit(bi_summary$drug_key))
} else {
  character(0)
}

# Defined above from the Stage 1 top 20, Stage 2 top 20,
# and all Stage 3 dual-query overlaps.
if (!exists("three_stage_overlap_n", inherits = FALSE)) {
  three_stage_overlap_n <- NA_integer_
}

enrichr_dual_threshold_n <- if (
  !is.null(bi_summary) &&
  nrow(bi_summary) > 0 &&
  "pass_threshold" %in% colnames(bi_summary)
) {
  sum(bi_summary$pass_threshold, na.rm = TRUE)
} else {
  NA_integer_
}

dgidb_overlap_names <- if (!is.null(bi_df) && nrow(bi_df) > 0) {
  paste(sort(unique(bi_df$drug)), collapse = "; ")
} else {
  NA_character_
}

drug_qc_summary <- data.frame(
  metric = c(
    "Enrichr Stage 1 unique standardized compounds",
    "Enrichr Stage 2 unique standardized compounds",
    "Enrichr Stage 3 cross-query compound overlaps",
    "Enrichr Stage 3 compounds with Combined Score >30 in both arms",
    "Enrichr compounds shared by the Stage 1 top 20, Stage 2 top 20, and Stage 3 dual-query analysis",
    "DGIdb deduplicated gene-compound pairs",
    "DGIdb unique standardized compounds",
    "DGIdb signature genes with curated records",
    "DGIdb cross-query compound overlaps",
    "DGIdb cross-query compound names"
  ),
  value = c(
    length(stage1_keys),
    length(stage2_keys),
    length(stage3_keys),
    enrichr_dual_threshold_n,
    three_stage_overlap_n,
    safe_nrow(dgidb_all),
    if (!is.null(dgidb_all) && nrow(dgidb_all) > 0) {
      safe_unique_n(dgidb_all$drug_key)
    } else {
      NA_integer_
    },
    if (!is.null(dgidb_all) && nrow(dgidb_all) > 0) {
      safe_unique_n(dgidb_all$gene)
    } else {
      NA_integer_
    },
    safe_nrow(bi_df),
    dgidb_overlap_names
  ),
  stringsAsFactors = FALSE
)

write.csv(
  drug_qc_summary,
  file.path(S6_DIR, "Drug_analysis_QC_summary.csv"),
  row.names = FALSE
)

cat("\nDrug-analysis QC summary:\n")
print(drug_qc_summary)


# =============================================================
# 7. PubMed 文献计数
# =============================================================
cat("\n--- 7. PubMed literature ---\n")
PUBMED_SEARCH_DATE <- "2026.07.19"
PUBMED_SEARCH_DATE_NCBI <- gsub("\\.", "/", PUBMED_SEARCH_DATE)
PUBMED_DATE_FILTER <- sprintf('("1800/01/01"[Date - Publication] : "%s"[Date - Publication])',
                              PUBMED_SEARCH_DATE_NCBI)
pubmed_drugs <- unique(c(
  "tocilizumab", "sorafenib", "cholecalciferol", "dexamethasone", "prednisolone",
  "cyclophosphamide", "sirolimus", "doxorubicin", "decitabine", "etoposide",
  "sivelestat", "alpha-1 antitrypsin", "canakinumab", "hydrocortisone",
  "methylprednisolone", "lithium", "prednisone",
  "vitamin d", "celecoxib", "imatinib", "valproic acid"))

cat(sprintf("Querying %d drugs against PubMed...\n", length(pubmed_drugs)))
pubmed_res_full <- do.call(rbind, lapply(seq_along(pubmed_drugs), function(i) {
  d <- pubmed_drugs[i]
  query_term <- sprintf('("%s"[Title/Abstract]) AND sepsis[Title/Abstract] AND %s',
                        d, PUBMED_DATE_FILTER)
  cnt <- tryCatch(rentrez::entrez_search(db = "pubmed",
                                         term = query_term,
                                         retmax = 0)$count, error = function(e) NA)
  Sys.sleep(0.35)
  cat(sprintf("  %-25s %5d\n", d, cnt))
  data.frame(drug = d, sepsis_pmids = cnt, search_date = PUBMED_SEARCH_DATE,
             pubmed_query = query_term, stringsAsFactors = FALSE)
}))
write.csv(pubmed_res_full, file.path(S6_DIR, "Drug_pubmed_sepsis_full.csv"), row.names = FALSE)

pubmed_res_full <- pubmed_res_full[order(-pubmed_res_full$sepsis_pmids), ]
pubmed_res_full$drug <- factor(pubmed_res_full$drug, levels = rev(pubmed_res_full$drug))
p_pubmed <- ggplot(pubmed_res_full, aes(sepsis_pmids, drug)) +
  geom_bar(stat = "identity", fill = "#4E79A7", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sepsis_pmids), hjust = -0.2, size = 3, color = "grey30") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "PubMed Literature: Drug Mentions in Sepsis",
       subtitle = sprintf('Query: ("drug"[Title/Abstract]) AND sepsis[Title/Abstract]; publication date <= %s',
                          PUBMED_SEARCH_DATE),
       x = "Number of PubMed articles", y = NULL,
       caption = sprintf("PubMed publication-date cutoff: %s | rentrez::entrez_search", PUBMED_SEARCH_DATE)) + theme_bw(11)
pdf(file.path(S6_DIR, "Figure_PubMed_literature.pdf"), 8, nrow(pubmed_res_full) * 0.35 + 1.5)
print(p_pubmed); dev.off()
png(file.path(S6_DIR, "Figure_PubMed_literature.png"), 8, nrow(pubmed_res_full) * 0.35 + 1.5,
    units = "in", res = 300); print(p_pubmed); dev.off()
cat("PubMed literature saved.\n")

cat(   "\n========== §6 PRIMARY DATABASE QUERIES COMPLETE ==========\n" )

# ============================================================
# Figure 6D — PubMed sepsis literature (full 21 drugs)
# ============================================================
suppressPackageStartupMessages(library(ggplot2))

drug_dir <- S6_DIR
if (!exists("PUBMED_SEARCH_DATE")) PUBMED_SEARCH_DATE <- "2026.07.19"

pub <- read.csv(file.path(drug_dir, "Drug_pubmed_sepsis_full.csv"), stringsAsFactors = FALSE)
if (!"search_date" %in% colnames(pub)) pub$search_date <- PUBMED_SEARCH_DATE
pub <- pub[order(-pub$sepsis_pmids), ]
pub$drug <- factor(pub$drug, levels = rev(pub$drug))

# Color: glucocorticoids in one color, others in another
pub$class <- ifelse(grepl("dexamethasone|hydrocortisone|prednisolone|prednisone|methylprednisolone",
                          pub$drug, ignore.case = TRUE),
                    "Glucocorticoid", "Other")
class_cols <- c("Glucocorticoid" = "#E15759", "Other" = "#4E79A7")

p <- ggplot(pub, aes(sepsis_pmids, drug, fill = class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sepsis_pmids), hjust = -0.2, size = 3, color = "grey30") +
  scale_fill_manual(values = class_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Figure 6D. PubMed sepsis literature",
       subtitle = sprintf('Query: ("drug"[Title/Abstract]) AND sepsis[Title/Abstract]; publication date <= %s',
                          PUBMED_SEARCH_DATE),
       x = "Publications", y = "",
       caption = sprintf("PubMed publication-date cutoff: %s | rentrez::entrez_search", PUBMED_SEARCH_DATE)) +
  theme_bw(11) +
  theme(panel.grid.major.y = element_blank())

pdf(file.path(drug_dir, "Figure_6D_pubmed_sepsis_full.pdf"), 8, nrow(pub) * 0.32 + 1.5)
print(p)
dev.off()

png(file.path(drug_dir, "Figure_6D_pubmed_sepsis_full.png"), 8, nrow(pub) * 0.32 + 1.5,
    units = "in", res = 300)
print(p)
dev.off()

cat(sprintf("Figure 6D (full) written. %d drugs.\n", nrow(pub)))
print(pub[, c("drug", "sepsis_pmids")])


# ============================================================
# §7 — Supplementary tables
# ============================================================
cat("\n========== §7 Supplementary ==========\n")

# Table 1
meta <- read.csv(file.path(S2_DIR,sprintf("meta_significant_%d_genes.csv",nrow(fisher_out))),stringsAsFactors=FALSE)
triple <- read.csv(file.path(S4_DIR,"Table_S4_32genes_triple_logFC.csv"),stringsAsFactors=FALSE)
colnames(meta)[colnames(meta)=="logFC_95233..D01."] <- "logFC_95233"
colnames(meta)[colnames(meta)=="fdr_95233..D01."] <- "fdr_95233"
tab <- merge(meta,triple[,c("gene","GSE65682","module","mean_abs")],by="gene",all.x=TRUE)
tab <- tab[order(-tab$mean_abs), , drop = FALSE]
tab$logFC_272769 <- round(tab$logFC_272769,3)
tab$logFC_95233 <- round(tab$logFC_95233,3)
tab$GSE65682 <- ifelse(is.na(tab$GSE65682),NA,round(tab$GSE65682,3))
fmt_fdr <- function(x) ifelse(x<0.001,sprintf("%.1e",x),sprintf("%.3f",x))
tab_out <- tab[, c("gene","module","direction","logFC_272769","fdr_272769",
                   "logFC_95233","fdr_95233","GSE65682","fisher_fdr",
                   "I2","individual_cohort_support"), drop = FALSE]
tab_out$fdr_272769 <- fmt_fdr(tab_out$fdr_272769)
tab_out$fdr_95233 <- fmt_fdr(tab_out$fdr_95233)
tab_out$fisher_fdr <- fmt_fdr(tab_out$fisher_fdr)
tab_out$I2 <- round(tab_out$I2,1)
tab_out$pass_05 <- ifelse(is.na(tab_out$GSE65682),"NA",
                          ifelse(abs(tab_out$GSE65682)>0.5,"YES","NO"))
colnames(tab_out) <- c("Gene","Module","Direction","logFC_272769","FDR_272769",
                       "logFC_95233","FDR_95233","logFC_65682","Fisher_FDR",
                       "I2_pct","Individual_cohort_support","|FC|>0.5_65682")
rownames(tab_out)<-NULL
write.csv(tab_out,file.path(S7_DIR,"Table1_32gene_signature.csv"),row.names=FALSE)
cat(sprintf("Table 1: %d genes\n",nrow(tab_out)))

# ============================================================
# 改进 #1: 三种筛选标准对比表

# 产出: 02_meta_analysis/Table_S3_screening_comparison.csv
# ============================================================

# B. 三种标准
# 标准A: Fisher FDR<0.05 + |logFC|>0.5 both cohorts + direction concordance (current, = fisher_out)
genes_A <- fisher_out$gene
N_A     <- length(genes_A)
sub_A   <- fisher_out
n_up_A  <- sum(sub_A$direction == "Up")
n_down_A<- sum(sub_A$direction == "Down")

# 标准B: 标准A + both cohorts nominal P<0.05
genes_B <- fisher_out$gene[!is.na(fisher_out$pval_272769) &
                             fisher_out$pval_272769 < 0.05 &
                             !is.na(fisher_out[[pval_COL_95233]]) &
                             fisher_out[[pval_COL_95233]] < 0.05]
N_B   <- length(genes_B)
sub_B <- fisher_out[fisher_out$gene %in% genes_B, ]

# 标准C: 标准A + both cohorts FDR<0.05
genes_C <- fisher_out$gene[!is.na(fisher_out$fdr_272769) &
                             fisher_out$fdr_272769 < 0.05 &
                             !is.na(fisher_out[[fdr_COL_95233]]) &
                             fisher_out[[fdr_COL_95233]] < 0.05]
N_C   <- length(genes_C)
sub_C <- fisher_out[fisher_out$gene %in% genes_C, ]

# C. 汇总表
screening_tab <- data.frame(
  Criterion = c(
    "A: Fisher FDR<0.05 + |logFC|>0.5 both cohorts + direction concordance (current 32-gene signature)",
    "B: Criterion A + nominal P<0.05 in both cohorts individually",
    "C: Criterion A + FDR<0.05 in both cohorts individually"
  ),
  N_genes = c(N_A, N_B, N_C),
  N_Up    = c(n_up_A, sum(sub_B$direction == "Up"), sum(sub_C$direction == "Up")),
  N_Down  = c(n_down_A, sum(sub_B$direction == "Down"), sum(sub_C$direction == "Down")),
  Genes   = c(
    paste(sort(genes_A), collapse = "; "),
    paste(sort(genes_B), collapse = "; "),
    paste(sort(genes_C), collapse = "; ")
  ),
  stringsAsFactors = FALSE
)
write.csv(screening_tab, file.path(S2_DIR, "Table_S3_screening_comparison.csv"), row.names = FALSE)

cat(sprintf("Screening comparison: A=%d  B=%d  C=%d genes\n", N_A, N_B, N_C))



# ============================================================
# Block B — §4 补充（置换检验 + PCA null + Mars3 深度解析 + 签名评分/AUC）
# ============================================================

suppressPackageStartupMessages(library(pROC))

genes_32 <- fisher_out$gene
genes_32_val <- intersect(genes_32, rownames(gmat_val))
n_val <- length(genes_32_val) 

# =============================================================
# #2: Outcome-label permutation test for sign concordance
# =============================================================

val_fc <- res_val[match(genes_32_val, res_val$gene), , drop = FALSE]

disc_dir <- fisher_out$direction[
  match(genes_32_val, fisher_out$gene)
]

disc_sign <- ifelse(
  disc_dir == "Up",
  1L,
  ifelse(disc_dir == "Down", -1L, NA_integer_)
)

valid_test <- is.finite(val_fc$logFC) & !is.na(disc_sign)

if (!any(valid_test)) {
  stop("No signature genes were evaluable for the sign-concordance permutation test.")
}

genes_tested <- genes_32_val[valid_test]
disc_sign_tested <- disc_sign[valid_test]
observed_sign <- sign(val_fc$logFC[valid_test])

n_comparable <- length(genes_tested)
actual_concordant <- sum(observed_sign == disc_sign_tested)
actual_concord_rate <- actual_concordant / n_comparable

up_idx <- disc_sign_tested == 1L
down_idx <- disc_sign_tested == -1L

n_up_disc <- sum(up_idx)
n_down_disc <- sum(down_idx)

actual_concord_up <- sum(observed_sign[up_idx] == 1L)
actual_concord_down <- sum(observed_sign[down_idx] == -1L)

# The group-only limma coefficient is equivalent in direction to
# the difference between the NS and S group means.
expr_tested <- gmat_val[genes_tested, , drop = FALSE]
group_chr <- as.character(group_val)

stopifnot(
  ncol(expr_tested) == length(group_chr),
  identical(rownames(expr_tested), genes_tested),
  !anyNA(expr_tested),
  all(group_chr %in% c("Survivor", "NonSurvivor"))
)

set.seed(2024)
n_perm <- 10000L

perm_concordant_count <- replicate(n_perm, {
  perm_group <- sample(group_chr, replace = FALSE)
  
  perm_ns <- perm_group == "NonSurvivor"
  perm_s <- perm_group == "Survivor"
  
  perm_logFC <- rowMeans(
    expr_tested[, perm_ns, drop = FALSE]
  ) - rowMeans(
    expr_tested[, perm_s, drop = FALSE]
  )
  
  sum(sign(perm_logFC) == disc_sign_tested)
})

perm_concord_rate <- perm_concordant_count / n_comparable

# Plus-one correction prevents an impossible empirical P value of zero.
n_exceed <- sum(perm_concordant_count >= actual_concordant)

perm_pval <- (n_exceed + 1) / (n_perm + 1)

p_text <- if (perm_pval < 0.001) {
  format(perm_pval, scientific = TRUE, digits = 3)
} else {
  sprintf("%.4f", perm_pval)
}

cat(sprintf(
  paste0(
    "Outcome-label permutation: Up=%d/%d, Down=%d/%d, ",
    "overall=%d/%d (%.1f%%), empirical P=%s\n"
  ),
  actual_concord_up,
  n_up_disc,
  actual_concord_down,
  n_down_disc,
  actual_concordant,
  n_comparable,
  100 * actual_concord_rate,
  p_text
))

permutation_summary <- data.frame(
  evaluable_genes = n_comparable,
  concordant_genes = actual_concordant,
  observed_concordance_rate = actual_concord_rate,
  permutations = n_perm,
  permutations_at_least_as_extreme = n_exceed,
  empirical_P = perm_pval,
  stringsAsFactors = FALSE
)

write.csv(
  permutation_summary,
  file.path(S4_DIR, "Sign_concordance_permutation_summary.csv"),
  row.names = FALSE
)

permutation_source <- data.frame(
  permutation = seq_len(n_perm),
  concordant_genes = perm_concordant_count,
  concordance_rate = perm_concord_rate
)

write.csv(
  permutation_source,
  file.path(
    S4_DIR,
    "Figure_outcome_label_permutation_source_data.csv"
  ),
  row.names = FALSE
)

bin_width <- 1 / n_comparable

p_perm <- ggplot(
  permutation_source,
  aes(x = concordance_rate)
) +
  geom_histogram(
    binwidth = bin_width,
    boundary = -bin_width / 2,
    fill = "grey80",
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = actual_concord_rate,
    color = "#E15759",
    linewidth = 1
  ) +
  annotate(
    "label",
    x = actual_concord_rate,
    y = Inf,
    label = sprintf(
      "Observed: %d/%d\nEmpirical P = %s",
      actual_concordant,
      n_comparable,
      p_text
    ),
    hjust = 1.05,
    vjust = 1.2,
    size = 3.2,
    color = "#E15759"
  ) +
 scale_x_continuous(
  limits = c(
    -bin_width / 2,
    1 + bin_width / 2
  ),
  breaks = seq(0, 1, 0.2),
  labels = scales::percent_format(accuracy = 1),
  expand = expansion(mult = 0)
) +
  coord_cartesian(
    xlim = c(0, 1)
  ) +
  labs(
    title = "Outcome-label permutation null distribution",
    subtitle = sprintf(
      "%s permutations of 28-day outcome labels",
      format(n_perm, big.mark = ",")
    ),
    x = "Directional concordance",
    y = "Number of permutations"
  ) +
  theme_bw(11)

save_plot(
  p_perm,
  "Figure_outcome_label_permutation_sign_concordance",
  S4_DIR,
  7,
  5
)
# =============================================================
# #3: Random-gene reference distributions for PC1 coherence
# =============================================================

required_pc1_objects <- c(
  "genes_32",
  "genes_32_val",
  "pca_32",
  "gmat_val",
  "common_genes",
  "S4_DIR"
)

missing_pc1_objects <- required_pc1_objects[
  !vapply(required_pc1_objects, exists, logical(1))
]

if (length(missing_pc1_objects) > 0) {
  stop(sprintf(
    "Missing objects required for the PC1 reference analysis: %s",
    paste(missing_pc1_objects, collapse = ", ")
  ))
}

n_signature_total <- length(genes_32)
n_signature_detectable <- length(genes_32_val)

if (n_signature_total != 32L) {
  warning(sprintf(
    "Expected a 32-gene discovery signature but found %d genes.",
    n_signature_total
  ))
}

if (n_signature_detectable < 3L) {
  stop(sprintf(
    paste0(
      "Only %d signature genes were detectable in GSE65682; ",
      "PCA reference analysis was not performed."
    ),
    n_signature_detectable
  ))
}

# -------------------------------------------------------------
# A. Check the observed signature-expression matrix
# -------------------------------------------------------------

signature_expr <- gmat_val[
  genes_32_val,
  ,
  drop = FALSE
]

stopifnot(
  nrow(signature_expr) == n_signature_detectable,
  identical(rownames(signature_expr), genes_32_val),
  ncol(signature_expr) == length(group_val)
)

if (anyNA(signature_expr)) {
  stop("The detectable signature-expression matrix contains missing values.")
}

signature_sd <- apply(
  signature_expr,
  1,
  sd,
  na.rm = TRUE
)

invalid_signature_variance <- (
  !is.finite(signature_sd) |
    signature_sd <= 0
)

if (any(invalid_signature_variance)) {
  invalid_genes <- names(signature_sd)[
    invalid_signature_variance
  ]
  
  stop(sprintf(
    paste0(
      "The detectable signature contains genes with ",
      "zero or non-finite variance: %s"
    ),
    paste(invalid_genes, collapse = ", ")
  ))
}

actual_pc1 <- summary(pca_32)$importance[2, 1] * 100

if (!is.finite(actual_pc1)) {
  stop("The observed signature PC1 variance was not finite.")
}

# -------------------------------------------------------------
# B. Define eligible random-gene pools
# -------------------------------------------------------------

eligible_gene_pool <- function(
    candidate_genes,
    expression_matrix,
    excluded_genes = character(0)
) {
  
  candidate_genes <- unique(
    intersect(
      candidate_genes,
      rownames(expression_matrix)
    )
  )
  
  candidate_genes <- setdiff(
    candidate_genes,
    excluded_genes
  )
  
  if (length(candidate_genes) == 0) {
    return(character(0))
  }
  
  candidate_expression <- expression_matrix[
    candidate_genes,
    ,
    drop = FALSE
  ]
  
  candidate_sd <- apply(
    candidate_expression,
    1,
    sd,
    na.rm = TRUE
  )
  
  candidate_genes[
    is.finite(candidate_sd) &
      candidate_sd > 0
  ]
}

# Pool 1: genes represented in both discovery cohorts and GSE65682
three_cohort_candidates <- intersect(
  common_genes,
  rownames(gmat_val)
)

three_cohort_pool <- eligible_gene_pool(
  candidate_genes = three_cohort_candidates,
  expression_matrix = gmat_val,
  excluded_genes = genes_32_val
)

# Pool 2: all genes represented in GSE65682 after gene summarization
full_platform_pool <- eligible_gene_pool(
  candidate_genes = rownames(gmat_val),
  expression_matrix = gmat_val,
  excluded_genes = genes_32_val
)

if (length(three_cohort_pool) < n_signature_detectable) {
  stop(sprintf(
    paste0(
      "The eligible three-cohort common-gene pool contained only ",
      "%d genes, fewer than the %d genes required per random set."
    ),
    length(three_cohort_pool),
    n_signature_detectable
  ))
}

if (length(full_platform_pool) < n_signature_detectable) {
  stop(sprintf(
    paste0(
      "The eligible full-platform pool contained only %d genes, ",
      "fewer than the %d genes required per random set."
    ),
    length(full_platform_pool),
    n_signature_detectable
  ))
}

cat(sprintf(
  paste0(
    "PC1 reference pools: three-cohort common=%d genes; ",
    "full GPL13667=%d genes; random-set size=%d genes\n"
  ),
  length(three_cohort_pool),
  length(full_platform_pool),
  n_signature_detectable
))

# -------------------------------------------------------------
# C. Function for random-gene PC1 reference simulation
# -------------------------------------------------------------

simulate_pc1_reference <- function(
    gene_pool,
    expression_matrix,
    n_genes,
    n_iterations,
    seed
) {
  
  if (length(gene_pool) < n_genes) {
    stop("The random-gene pool was smaller than the requested set size.")
  }
  
  set.seed(seed)
  
  random_pc1 <- numeric(n_iterations)
  
  for (iteration in seq_len(n_iterations)) {
    
    random_genes <- sample(
      gene_pool,
      size = n_genes,
      replace = FALSE
    )
    
    random_expression <- expression_matrix[
      random_genes,
      ,
      drop = FALSE
    ]
    
    random_pca <- prcomp(
      t(random_expression),
      center = TRUE,
      scale. = TRUE
    )
    
    random_pc1[iteration] <- (
      summary(random_pca)$importance[2, 1] * 100
    )
  }
  
  random_pc1
}

# -------------------------------------------------------------
# D. Generate the two reference distributions
# -------------------------------------------------------------

n_random_sets <- 10000L

random_pc1_common <- simulate_pc1_reference(
  gene_pool = three_cohort_pool,
  expression_matrix = gmat_val,
  n_genes = n_signature_detectable,
  n_iterations = n_random_sets,
  seed = 2024
)

random_pc1_full <- simulate_pc1_reference(
  gene_pool = full_platform_pool,
  expression_matrix = gmat_val,
  n_genes = n_signature_detectable,
  n_iterations = n_random_sets,
  seed = 2025
)

if (
  any(!is.finite(random_pc1_common)) ||
  any(!is.finite(random_pc1_full))
) {
  stop("Non-finite values were generated in a PC1 reference distribution.")
}

# -------------------------------------------------------------
# E. Calculate descriptive upper-tail proportions
# -------------------------------------------------------------

common_exceedances <- sum(
  random_pc1_common >= actual_pc1
)

full_exceedances <- sum(
  random_pc1_full >= actual_pc1
)

# Plus-one correction avoids a value of exactly zero.
common_upper_tail <- (
  common_exceedances + 1
) / (
  n_random_sets + 1
)

full_upper_tail <- (
  full_exceedances + 1
) / (
  n_random_sets + 1
)

# -------------------------------------------------------------
# F. Save the summary table
# -------------------------------------------------------------

pc1_reference_summary <- data.frame(
  reference_pool = c(
    "Three-cohort common-gene pool",
    "Full GPL13667 gene pool"
  ),
  pool_size = c(
    length(three_cohort_pool),
    length(full_platform_pool)
  ),
  excluded_signature_genes = c(
    n_signature_detectable,
    n_signature_detectable
  ),
  genes_per_random_set = c(
    n_signature_detectable,
    n_signature_detectable
  ),
  random_sets = c(
    n_random_sets,
    n_random_sets
  ),
  observed_PC1_variance_pct = c(
    actual_pc1,
    actual_pc1
  ),
  random_mean_PC1_variance_pct = c(
    mean(random_pc1_common),
    mean(random_pc1_full)
  ),
  random_SD_PC1_variance_pct = c(
    sd(random_pc1_common),
    sd(random_pc1_full)
  ),
  random_median_PC1_variance_pct = c(
    median(random_pc1_common),
    median(random_pc1_full)
  ),
  random_sets_at_least_observed = c(
    common_exceedances,
    full_exceedances
  ),
  adjusted_upper_tail_proportion = c(
    common_upper_tail,
    full_upper_tail
  ),
  stringsAsFactors = FALSE
)

write.csv(
  pc1_reference_summary,
  file.path(
    S4_DIR,
    "Table_PCA_random_gene_reference_summary.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------
# G. Save the complete source distributions
# -------------------------------------------------------------

pc1_reference_source <- rbind(
  data.frame(
    random_set = seq_len(n_random_sets),
    reference_pool = "Three-cohort common-gene pool",
    PC1_variance_pct = random_pc1_common,
    stringsAsFactors = FALSE
  ),
  data.frame(
    random_set = seq_len(n_random_sets),
    reference_pool = "Full GPL13667 gene pool",
    PC1_variance_pct = random_pc1_full,
    stringsAsFactors = FALSE
  )
)

pc1_reference_source$reference_pool <- factor(
  pc1_reference_source$reference_pool,
  levels = c(
    "Three-cohort common-gene pool",
    "Full GPL13667 gene pool"
  )
)

write.csv(
  pc1_reference_source,
  file.path(
    S4_DIR,
    "Figure_PCA_random_gene_reference_source_data.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------
# H. Prepare labels for the two-panel figure
# -------------------------------------------------------------

pc1_label_data <- data.frame(
  reference_pool = factor(
    c(
      "Three-cohort common-gene pool",
      "Full GPL13667 gene pool"
    ),
    levels = levels(
      pc1_reference_source$reference_pool
    )
  ),
  x = c(
    actual_pc1,
    actual_pc1
  ),
  y = c(
    Inf,
    Inf
  ),
  label = c(
    sprintf(
      paste0(
        "Observed PC1 = %.1f%%\n",
        "Upper-tail proportion = %.3g"
      ),
      actual_pc1,
      common_upper_tail
    ),
    sprintf(
      paste0(
        "Observed PC1 = %.1f%%\n",
        "Upper-tail proportion = %.3g"
      ),
      actual_pc1,
      full_upper_tail
    )
  ),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------
# I. Plot the reference distributions
# -------------------------------------------------------------

p_pc1_reference <- ggplot(
  pc1_reference_source,
  aes(x = PC1_variance_pct)
) +
  geom_histogram(
    bins = 60,
    fill = "grey80",
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = actual_pc1,
    color = "#E15759",
    linewidth = 1
  ) +
  geom_label(
    data = pc1_label_data,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.2,
    size = 3,
    color = "#E15759"
  ) +
  facet_wrap(
    ~reference_pool,
    nrow = 1,
    scales = "fixed"
  ) +
  labs(
    title = "PC1 coherence relative to random-gene reference sets",
    subtitle = sprintf(
      paste0(
        "Observed PCA: %d detectable genes from the ",
        "%d-gene signature; %s random sets per reference pool."
      ),
      n_signature_detectable,
      n_signature_total,
      format(n_random_sets, big.mark = ",")
    ),
    x = "PC1 variance explained (%)",
    y = "Number of random gene sets"
  ) +
  theme_bw(11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

save_plot(
  p_pc1_reference,
  "Figure_PCA_random_gene_reference",
  S4_DIR,
  10,
  5.5
)

cat(sprintf(
  paste0(
    "PC1 random-gene reference: observed=%.1f%%; ",
    "three-cohort upper tail=%.4g; ",
    "full-platform upper tail=%.4g\n"
  ),
  actual_pc1,
  common_upper_tail,
  full_upper_tail
))

# =============================================================
# #4a: Per-endotype per-gene logFC heatmap (Mars3 deep-dive)
# =============================================================
endo_levels <- c("Mars1", "Mars2", "Mars3", "Mars4")
per_endo_fc <- matrix(NA_real_, nrow = n_val, ncol = 4,
                      dimnames = list(genes_32_val, endo_levels))

for (j in seq_along(endo_levels)) {
  e <- endo_levels[j]
  ie <- which(endo_val == e)
  ge <- group_val[ie]
  if (length(unique(ge)) < 2) next
  de <- model.matrix(~0 + ge); colnames(de) <- levels(ge)
  fe <- lmFit(gmat_val[genes_32_val, ie], de)
  f2e <- contrasts.fit(fe, makeContrasts(NonSurvivor - Survivor, levels = de))
  f2e <- eBayes(f2e, trend = TRUE)
  per_endo_fc[, j] <- f2e$coefficients[, 1]
}

# Direction-discordance flag: sign in this endotype differs from discovery direction
REVERSAL_STAR_CUT <- 0.2
reversal_flag <- sign(per_endo_fc) != matrix(disc_sign, nrow = n_val, ncol = 4)
colnames(reversal_flag) <- endo_levels

# Display star only for substantive direction reversal
substantive_reversal_flag <- reversal_flag & abs(per_endo_fc) > REVERSAL_STAR_CUT
colnames(substantive_reversal_flag) <- endo_levels


# ============================================================
# Per-endotype per-gene logFC heatmap
# ============================================================

# Row annotations
mod_vec <- ifelse(
  genes_32_val %in% names(module_map),
  as.character(module_map[genes_32_val]),
  "Unclassified"
)

mod_vec <- factor(mod_vec, levels = MODULE_LEVELS)

mod_colors_heat <- c(
  "Neutrophil degranulation"   = "#E15759",
  "Cell cycle / proliferation" = "#4E79A7",
  "Inflammatory / Down"        = "#59A14F",
  "Other"                      = "#B07AA1",
  "Unclassified"               = "grey70"
)

row_annot <- data.frame(Module = mod_vec, row.names = genes_32_val)

# Sort by module then gene
ord <- order(mod_vec, genes_32_val)
hm_mat <- per_endo_fc[ord, , drop = FALSE]
row_annot <- row_annot[ord, , drop = FALSE]

# Ensure annotation color mapping is complete
row_annot <- as.data.frame(row_annot)
row_annot$Module <- as.character(row_annot$Module)

missing_modules <- setdiff(unique(row_annot$Module), names(mod_colors_heat))
if (length(missing_modules) > 0) {
  mod_colors_heat[missing_modules] <- "grey70"
}

row_annot$Module <- factor(row_annot$Module, levels = names(mod_colors_heat))

annotation_colors_heat <- list(
  Module = mod_colors_heat[levels(row_annot$Module)]
)

# Heatmap color scale
max_abs <- max(abs(hm_mat), na.rm = TRUE)
palette <- colorRampPalette(c("#4E79A7", "white", "#E15759"))(100)
breaks  <- seq(-max_abs, max_abs, length.out = 101)

# Mark reversal cells with *
rev_mat <- substantive_reversal_flag[ord, , drop = FALSE]
hm_numbers <- matrix(sprintf("%.2f", hm_mat), nrow = nrow(hm_mat), ncol = ncol(hm_mat))
hm_numbers[rev_mat] <- paste0(hm_numbers[rev_mat], "*")
dimnames(hm_numbers) <- dimnames(hm_mat)

# PDF output
pdf(file.path(S4_DIR, "Figure_Mars3_per_endo_logFC_heatmap.pdf"), 8.6, 8.8)

pheatmap(
  hm_mat,
  color = palette,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annot,
  annotation_colors = annotation_colors_heat,
  annotation_names_row = FALSE,
  display_numbers = hm_numbers,
  number_format = "%.2f",
  number_color = "grey20",
  fontsize_number = 6.8,
  fontsize_row = 8.5,
  fontsize_col = 10,
  cellheight = 13,
  cellwidth = 42,
  angle_col = 0,
  border_color = NA,
  main = "Per-endotype log2FC by MARS endotype",
  na_col = "grey90"
)

dev.off()

png(file.path(S4_DIR, "Figure_Mars3_per_endo_logFC_heatmap.png"), 8.6, 8.8, units = "in", res = 300)

pheatmap(
  hm_mat,
  color = palette,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annot,
  annotation_colors = annotation_colors_heat,
  annotation_names_row = FALSE,
  display_numbers = hm_numbers,
  number_format = "%.2f",
  number_color = "grey20",
  fontsize_number = 6.8,
  fontsize_row = 8.5,
  fontsize_col = 10,
  cellheight = 13,
  cellwidth = 42,
  angle_col = 0,
  border_color = NA,
  main = "Per-endotype log2FC by MARS endotype",
  na_col = "grey90"
)

dev.off()
cat("Per-endotype heatmap saved.\n")

# Save source data
per_endo_out <- data.frame(
  gene = genes_32_val[ord],
  module = as.character(mod_vec[ord]),
  Mars1 = hm_mat[, 1], Mars2 = hm_mat[, 2], Mars3 = hm_mat[, 3], Mars4 = hm_mat[, 4],
  N_substantive_reversals = rowSums(substantive_reversal_flag[ord, , drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(per_endo_out, file.path(S4_DIR, "Table_S5_per_endotype_logFC.csv"), row.names = FALSE)
cat("Per-endotype heatmap saved.\n")

# =============================================================
# #4b: Mars3 reversal gene module attribution
# =============================================================
mars3_rev <- genes_32_val[substantive_reversal_flag[, "Mars3"]]
mars3_rev_mod <- ifelse(mars3_rev %in% names(module_map),
                        as.character(module_map[mars3_rev]), "Unclassified")
rev_summary <- data.frame(
  gene = mars3_rev,
  module = mars3_rev_mod,
  logFC_Mars3  = per_endo_fc[substantive_reversal_flag[, "Mars3"], "Mars3"],
  discovery_dir = disc_dir[substantive_reversal_flag[, "Mars3"]],
  stringsAsFactors = FALSE
)
rev_summary <- rev_summary[order(rev_summary$module, rev_summary$gene), ]
write.csv(rev_summary, file.path(S4_DIR, "Table_S6_Mars3_reversal_genes.csv"), row.names = FALSE)

cat(sprintf("Mars3 reversals: %d/%d genes (%.0f%%)\n",
            nrow(rev_summary), n_val, nrow(rev_summary) / n_val * 100))
cat("Module attribution:\n")
print(table(rev_summary$module))

# =============================================================
# #7: 32-gene signature score + AUC in GSE65682
# =============================================================
# 断言：sig_score/gmat_val列序/group_val/endo_val 四者严格对齐
stopifnot(length(group_val) == ncol(gmat_val),
          length(endo_val)  == ncol(gmat_val),
          all(rownames(pca_32$x) == colnames(gmat_val)))

# PC1 as signature score, direction: higher = NonSurvivor
sig_score <- as.numeric(pca_32$x[, 1])
if (mean(sig_score[group_val == "NonSurvivor"]) < mean(sig_score[group_val == "Survivor"])) {
  sig_score <- -sig_score
}

roc_obj <- roc(group_val, sig_score, levels = c("Survivor", "NonSurvivor"), direction = "<")
auc_val <- as.numeric(auc(roc_obj))
auc_ci  <- as.numeric(ci(roc_obj, method = "delong"))

pdf(file.path(S4_DIR, "Figure_32gene_AUC_GSE65682.pdf"), 6, 6)
plot(roc_obj, main = "32-gene Signature Score in GSE65682",
     col = "#E15759", lwd = 2.5, print.auc = FALSE)
legend("bottomright", sprintf("AUC = %.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]),
       bty = "n", cex = 1.1)
dev.off()
png(file.path(S4_DIR, "Figure_32gene_AUC_GSE65682.png"), 6, 6, units = "in", res = 300)
plot(roc_obj, main = "32-gene Signature Score in GSE65682",
     col = "#E15759", lwd = 2.5, print.auc = FALSE)
legend("bottomright", sprintf("AUC = %.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]),
       bty = "n", cex = 1.1)
dev.off()

# Save scores — 用 PCA 行名，不做位置索引
auc_out <- data.frame(
  Sample   = rownames(pca_32$x),
  Signature_Score = sig_score,
  Outcome  = as.character(group_val),
  Endotype = as.character(endo_val),
  stringsAsFactors = FALSE
)
write.csv(auc_out, file.path(S4_DIR, "Table_S7_signature_score_AUC.csv"), row.names = FALSE)

cat(sprintf("GSE65682 AUC: %.3f (%.3f-%.3f)\n", auc_val, auc_ci[1], auc_ci[3]))

# Per-endotype AUC — 固定 direction="<"，Mars3 低 AUC 即反转信号
cat("\nPer-endotype AUC:\n")
for (e in endo_levels) {
  idx <- which(endo_val == e)
  grp <- group_val[idx]
  if (length(unique(grp)) < 2) next
  sc  <- sig_score[idx]
  roc_e <- tryCatch(roc(grp, sc, levels = c("Survivor", "NonSurvivor"), direction = "<"),
                    error = function(e) NULL)
  if (!is.null(roc_e)) {
    auc_e <- as.numeric(auc(roc_e))
    cat(sprintf("  %s: AUC=%.3f\n", e, auc_e))
  }
}

suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

genes_32      <- fisher_out$gene
genes_32_val  <- intersect(genes_32, rownames(gmat_val))
n_val         <- length(genes_32_val)
endo_levels   <- c("Mars1", "Mars2", "Mars3", "Mars4")
common_genes_vec <- intersect(common_genes, rownames(gmat_val))
all_genes_gpl13667 <- rownames(gmat_val)

# ============================================================
# 1. AUC 补充: Youden + Sensitivity/Specificity + Score-Outcome
# ============================================================
cat("\n===== 1. AUC + Youden + Score-Outcome =====\n")

sig_score <- as.numeric(pca_32$x[, 1])
if (mean(sig_score[group_val == "NonSurvivor"]) < mean(sig_score[group_val == "Survivor"])) {
  sig_score <- -sig_score
}

roc_obj <- roc(group_val, sig_score, levels = c("Survivor", "NonSurvivor"), direction = "<")
auc_val <- as.numeric(auc(roc_obj))
auc_ci  <- as.numeric(ci(roc_obj, method = "delong"))

# Youden index（防多并列）
roc_coords    <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))
roc_coords    <- roc_coords[1, ]
youden_thresh <- roc_coords$threshold
youden_sens   <- roc_coords$sensitivity
youden_spec   <- roc_coords$specificity
cat(sprintf("Youden: threshold=%.3f | Sens=%.2f | Spec=%.2f\n",
            youden_thresh, youden_sens, youden_spec))

# Score-outcome relationship（非 calibration）
sig_binary <- ifelse(group_val == "NonSurvivor", 1, 0)
calib_df   <- data.frame(predicted = sig_score, observed = sig_binary)

p_score_rel <- ggplot(calib_df, aes(predicted, observed)) +
  geom_point(alpha = 0.03, size = 0.5, position = position_jitter(height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.7,
              color = "#E15759", fill = "#E15759", alpha = 0.15) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  labs(
    title = "PC1 score and 28-day outcome",
    subtitle = "Descriptive LOESS smoother (span = 0.7)",
    x = "PC1 signature score",
    y = "Observed 28-day mortality"
  ) +
  theme_bw(12)

pdf(file.path(S4_DIR, "Figure_score_outcome_GSE65682.pdf"), 6, 6); print(p_score_rel); dev.off()
png(file.path(S4_DIR, "Figure_score_outcome_GSE65682.png"), 6, 6, units = "in", res = 300); print(p_score_rel); dev.off()

auc_metrics <- data.frame(
  Metric = c("AUC", "AUC_95CI_low", "AUC_95CI_high",
             "Youden_threshold", "Sensitivity", "Specificity",
             "N_total", "N_Survivor", "N_NonSurvivor"),
  Value = c(auc_val, auc_ci[1], auc_ci[3],
            youden_thresh, youden_sens, youden_spec,
            length(group_val), sum(group_val == "Survivor"), sum(group_val == "NonSurvivor")),
  stringsAsFactors = FALSE
)
write.csv(auc_metrics, file.path(S4_DIR, "Table_AUC_metrics.csv"), row.names = FALSE)
cat("AUC + score-outcome saved.\n")

# ============================================================
# 3. Mars3 效量阈值分类 + 瀑布图
# ============================================================
cat("\n===== 3. Mars3 effect-size waterfall =====\n")

disc_dir  <- fisher_out$direction[match(genes_32_val, fisher_out$gene)]
disc_sign <- ifelse(disc_dir == "Up", 1L, -1L)

# 重算 per-endotype logFC
per_endo_fc <- matrix(NA_real_, nrow = n_val, ncol = 4,
                      dimnames = list(genes_32_val, endo_levels))
for (j in seq_along(endo_levels)) {
  e <- endo_levels[j]; ie <- which(endo_val == e); ge <- group_val[ie]
  if (length(unique(ge)) < 2) next
  de <- model.matrix(~0 + ge); colnames(de) <- levels(ge)
  fe <- lmFit(gmat_val[genes_32_val, ie], de)
  f2e <- contrasts.fit(fe, makeContrasts(NonSurvivor - Survivor, levels = de))
  f2e <- eBayes(f2e, trend = TRUE)
  per_endo_fc[, j] <- f2e$coefficients[, 1]
}

EFFECT_CUT <- 0.2
mars3_logfc <- per_endo_fc[, "Mars3"]

mars3_strong_rev   <- (abs(mars3_logfc) > EFFECT_CUT) & (sign(mars3_logfc) != disc_sign)
mars3_weak_rev     <- (abs(mars3_logfc) <= EFFECT_CUT) & (sign(mars3_logfc) != disc_sign)
mars3_strong_norev <- (abs(mars3_logfc) > EFFECT_CUT) & (sign(mars3_logfc) == disc_sign)

cat(sprintf("Mars3 reversal breakdown (|logFC|>%.1f):\n", EFFECT_CUT))
cat(sprintf("  Strong reversal:          %d/%d (%.0f%%)\n",
            sum(mars3_strong_rev, na.rm = TRUE), n_val, sum(mars3_strong_rev, na.rm = TRUE)/n_val*100))
cat(sprintf("  Weak reversal (near-zero): %d/%d (%.0f%%)\n",
            sum(mars3_weak_rev, na.rm = TRUE), n_val, sum(mars3_weak_rev, na.rm = TRUE)/n_val*100))
cat(sprintf("  Strong non-reversal:       %d/%d (%.0f%%)\n",
            sum(mars3_strong_norev, na.rm = TRUE), n_val, sum(mars3_strong_norev, na.rm = TRUE)/n_val*100))

# 分类表
rev_class <- data.frame(
  gene = genes_32_val,
  logFC_Mars3 = round(mars3_logfc, 3),
  discovery_dir = disc_dir,
  class = ifelse(mars3_strong_rev, "Strong_reversal",
                 ifelse(mars3_weak_rev, "Near_zero",
                        ifelse(abs(mars3_logfc) > EFFECT_CUT, "Consistent", "Near_zero_consistent"))),
  stringsAsFactors = FALSE
)
rev_class <- rev_class[order(-abs(rev_class$logFC_Mars3)), ]
write.csv(rev_class, file.path(S4_DIR, "Table_S6b_Mars3_reversal_classified.csv"), row.names = FALSE)

# Waterfall plot
rev_class$gene <- factor(rev_class$gene, levels = rev_class$gene)

rev_class$class <- factor(
  rev_class$class,
  levels = c("Strong_reversal", "Near_zero", "Near_zero_consistent", "Consistent")
)

class_cols <- c(
  "Strong_reversal" = "#E15759",
  "Near_zero" = "grey60",
  "Near_zero_consistent" = "grey85",
  "Consistent" = "#4E79A7"
)

near_zero_total <- sum(abs(mars3_logfc) <= EFFECT_CUT, na.rm = TRUE)
consistent_total <- sum(abs(mars3_logfc) > EFFECT_CUT & !mars3_strong_rev, na.rm = TRUE)

p_waterfall <- ggplot(rev_class, aes(logFC_Mars3, gene, fill = class)) +
  geom_col(width = 0.82) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_vline(xintercept = c(-EFFECT_CUT, EFFECT_CUT), lty = "dashed", linewidth = 0.3) +
  scale_fill_manual(
    values = class_cols,
    breaks = c("Strong_reversal", "Near_zero", "Near_zero_consistent", "Consistent"),
    labels = c("Reversal", "Near-zero discordant", "Near-zero concordant", "Concordant"),
    name = NULL
  ) +
  labs(
    title = "Mars3 gene-level log2FC",
    subtitle = sprintf(
      "Reversal=%d | Near-zero=%d | Concordant >%.1f=%d",
      sum(mars3_strong_rev, na.rm = TRUE),
      near_zero_total,
      EFFECT_CUT,
      consistent_total
    ),
    x = "log2FC in Mars3",
    y = NULL
  ) +
  theme_bw(10) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )
pdf(file.path(S4_DIR, "Figure_Mars3_waterfall_classified.pdf"), 4.2, 8.8)
print(p_waterfall)
dev.off()

png(file.path(S4_DIR, "Figure_Mars3_waterfall_classified.png"), 4.2, 8.8, units = "in", res = 300)
print(p_waterfall)
dev.off()
cat("Mars3 waterfall saved.\n")

# ============================================================
# 4. Per-endotype AUC + 95% CI + Forest plot
# ============================================================
cat("\n===== 4. Per-endotype AUC with CI =====\n")

endo_auc <- data.frame(
  Endotype = character(), AUC = numeric(), CI_low = numeric(), CI_high = numeric(),
  N = numeric(), N_NS = numeric(), N_S = numeric(), stringsAsFactors = FALSE
)
for (e in endo_levels) {
  idx <- which(endo_val == e); grp <- group_val[idx]; sc <- sig_score[idx]
  if (length(unique(grp)) < 2) next
  roc_e <- tryCatch(roc(grp, sc, levels = c("Survivor", "NonSurvivor"), direction = "<"),
                    error = function(e) NULL)
  if (is.null(roc_e)) next
  ci_e <- tryCatch(as.numeric(ci(roc_e, method = "delong")), error = function(e) rep(NA, 3))
  endo_auc <- rbind(endo_auc, data.frame(
    Endotype = e, AUC = as.numeric(auc(roc_e)),
    CI_low = ci_e[1], CI_high = ci_e[3],
    N = length(grp), N_NS = sum(grp == "NonSurvivor"), N_S = sum(grp == "Survivor"),
    stringsAsFactors = FALSE
  ))
}
write.csv(endo_auc, file.path(S4_DIR, "Table_per_endotype_AUC.csv"), row.names = FALSE)
print(endo_auc)

# Forest plot
p_auc_forest <- ggplot(endo_auc, aes(AUC, Endotype)) +
  geom_point(size = 3, color = "#E15759") +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    width = 0.2,
    orientation = "y",
    linewidth = 1,
    color = "#E15759"
  ) +
  geom_vline(xintercept = 0.5, lty = "dashed", lwd = 0.5, color = "grey50") +
  geom_text(aes(label = sprintf("AUC=%.3f\n(N=%d, NS=%d)", AUC, N, N_NS)),
            hjust = -0.1, size = 3, color = "grey30") +
  scale_x_continuous(limits = c(0.3, 1.05)) +
  labs(
    title = "Per-Endotype AUC (30 detectable signature genes)",
    subtitle = "Fixed global direction (higher = NS) | Error bars: 95% CI (DeLong)",
    x = "AUC",
    y = NULL
  ) +
  theme_bw(11)

pdf(file.path(S4_DIR, "Figure_per_endotype_AUC_forest.pdf"), 8, 4); print(p_auc_forest); dev.off()
png(file.path(S4_DIR, "Figure_per_endotype_AUC_forest.png"), 8, 4, units = "in", res = 300); print(p_auc_forest); dev.off()
cat("Per-endotype AUC forest saved.\n")

# ============================================================
# 5. 平台探针数统计
# ============================================================
cat("\n===== 5. Platform probe count =====\n")

probe_counts <- data.frame(
  gene = genes_32_val,
  n_probes_272769 = sapply(genes_32_val, function(g) sum(map17692$gene == g)),
  n_probes_95233  = sapply(genes_32_val, function(g) sum(map570$gene == g)),
  stringsAsFactors = FALSE
)
probe_counts$disparity <- abs(probe_counts$n_probes_272769 - probe_counts$n_probes_95233)
probe_counts$I2 <- fisher_out$I2[match(genes_32_val, fisher_out$gene)]

cat(sprintf("Mean probes: 272769=%.1f  95233=%.1f\n",
            mean(probe_counts$n_probes_272769), mean(probe_counts$n_probes_95233)))
cat(sprintf("Correlation(I2, probe_disparity): r=%.3f\n",
            cor(probe_counts$disparity, probe_counts$I2, use = "complete.obs")))
write.csv(probe_counts, file.path(S4_DIR, "Table_platform_probe_counts.csv"), row.names = FALSE)

# ============================================================
# 6. Tocilizumab 间接性 + 化疗药占比
# ============================================================
cat("\n===== 6. Drug evidence quantitative check =====\n")

il6_related <- intersect(c("IL6R", "IL6ST", "JAK1", "JAK2", "STAT3", "SOCS3"), genes_32_val)
cat(sprintf("IL6-pathway genes in 32-gene signature: %d\n", length(il6_related)))
cat(sprintf("  Present: %s\n", if(length(il6_related) > 0) paste(il6_related, collapse = ", ") else "NONE"))
cat(sprintf("  IL6R itself: %s\n", ifelse("IL6R" %in% genes_32_val, "YES (direct target)", "NO (indirect only)")))

chemo_kw <- c("doxorubicin", "etoposide", "cyclophosphamide", "decitabine", "bleomycin",
              "cisplatin", "methotrexate", "fluorouracil", "gemcitabine", "vincristine",
              "paclitaxel", "docetaxel", "irinotecan", "topotecan", "mitoxantrone",
              "daunorubicin", "epirubicin", "idarubicin", "melphalan", "busulfan",
              "carboplatin", "oxaliplatin", "ifosfamide", "tamoxifen", "vorinostat")
if (exists("bi_summary")) {
  bi_drugs <- bi_summary$drug
  chemo_hits <- bi_drugs[tolower(bi_drugs) %in% chemo_kw]
  cat(sprintf(
    paste0(
      "Chemotherapy-associated compounds among %d ",
      "dual-query Enrichr overlaps: %d (%.0f%%)\n"
    ),
    length(bi_drugs),
    length(chemo_hits),
    length(chemo_hits) / max(1, length(bi_drugs)) * 100
  ))
  if (length(chemo_hits) > 0) cat(sprintf("  Names: %s\n", paste(chemo_hits, collapse = ", ")))
}

if (exists("dgidb_all")) {
  top2a_count <- sum(dgidb_all$gene == "TOP2A", na.rm = TRUE)
  total_ints  <- nrow(dgidb_all)
  cat(sprintf("TOP2A accounts for %d/%d DGIdb interactions (%.0f%%)\n",
              top2a_count, total_ints, top2a_count/max(1, total_ints)*100))
}

# ============================================================
# 安全 concordance 计算（防除以 0 + 防 NA）
# ============================================================
safe_concordance_pct <- function(val_logFC, disc_sign) {
  if (length(val_logFC) == 0) return(NA_real_)
  ok <- !is.na(val_logFC) & !is.na(disc_sign)
  if (!any(ok)) return(NA_real_)
  mean(sign(val_logFC[ok]) == disc_sign[ok]) * 100
}

# ============================================================
# 阈值敏感性分析
# ============================================================
lfc_thresholds <- c(0.3, 0.5, 0.7, 1.0)

threshold_summary <- do.call(rbind, lapply(lfc_thresholds, function(th) {
  lfc1 <- meta_df$logFC_272769
  lfc2 <- meta_df[[logFC_COL_95233]]
  
  keep <- !is.na(meta_df$fisher_fdr) &
    meta_df$fisher_fdr < 0.05 &
    !is.na(lfc1) & !is.na(lfc2) &
    abs(lfc1) > th &
    abs(lfc2) > th &
    sign(lfc1) == sign(lfc2) &
    sign(lfc1) != 0
  
  selected_genes <- meta_df$gene[keep]
  discovery_sign <- sign(lfc1[keep])
  
  genes_found <- intersect(selected_genes, res_val$gene)
  val_idx <- match(genes_found, res_val$gene)
  disc_idx <- match(genes_found, selected_genes)
  
  val_logfc <- res_val$logFC[val_idx]
  disc_sign_found <- discovery_sign[disc_idx]
  
  evaluable <- !is.na(val_logfc) & !is.na(disc_sign_found) & disc_sign_found != 0
  concordant <- sign(val_logfc[evaluable]) == disc_sign_found[evaluable]
  
  data.frame(
    log2FC_threshold = th,
    selected_genes = length(selected_genes),
    upregulated_genes = sum(lfc1[keep] > 0, na.rm = TRUE),
    downregulated_genes = sum(lfc1[keep] < 0, na.rm = TRUE),
    evaluable_in_GSE65682 = sum(evaluable),
    directionally_concordant_in_GSE65682 = sum(concordant, na.rm = TRUE),
    directional_concordance_rate = ifelse(
      sum(evaluable) > 0,
      round(100 * sum(concordant, na.rm = TRUE) / sum(evaluable), 1),
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  threshold_summary,
  file.path(S2_DIR, "Table_SX_logFC_threshold_sensitivity.csv"),
  row.names = FALSE
)

print(threshold_summary)


# ============================================================
# Figure: Pairwise Hallmark NES-direction concordance
# ============================================================

if (!exists("hall_3") ||
    !is.data.frame(hall_3) ||
    nrow(hall_3) == 0) {
  stop("Three-cohort Hallmark GSEA results were unavailable.")
}

gsea_plot_data <- hall_3

gsea_pairs <- list(
  list(
    x = "NES_GSE272769",
    y = "NES_GSE95233",
    xlab = "GSE272769 NES",
    ylab = "GSE95233 D01 NES"
  ),
  list(
    x = "NES_GSE272769",
    y = "NES_GSE65682",
    xlab = "GSE272769 NES",
    ylab = "GSE65682 NES"
  ),
  list(
    x = "NES_GSE95233",
    y = "NES_GSE65682",
    xlab = "GSE95233 D01 NES",
    ylab = "GSE65682 NES"
  )
)

direction_cols <- c(
  "Concordant" = "#4E79A7",
  "Discordant" = "#B07AA1"
)

plot_list <- lapply(gsea_pairs, function(pair) {
  
  pair_df <- data.frame(
    pathway = gsea_plot_data$pathway,
    x_NES = suppressWarnings(
      as.numeric(gsea_plot_data[[pair$x]])
    ),
    y_NES = suppressWarnings(
      as.numeric(gsea_plot_data[[pair$y]])
    ),
    stringsAsFactors = FALSE
  )
  
  pair_df <- pair_df[
    is.finite(pair_df$x_NES) &
      is.finite(pair_df$y_NES),
    ,
    drop = FALSE
  ]
  
  if (nrow(pair_df) == 0) {
    stop(sprintf(
      "No finite NES values were available for %s versus %s.",
      pair$x,
      pair$y
    ))
  }
  
  pair_df$direction_status <- ifelse(
    sign(pair_df$x_NES) == sign(pair_df$y_NES),
    "Concordant",
    "Discordant"
  )
  
  pair_df$direction_status <- factor(
    pair_df$direction_status,
    levels = c("Concordant", "Discordant")
  )
  
  n_total <- nrow(pair_df)
  n_concordant <- sum(
    pair_df$direction_status == "Concordant",
    na.rm = TRUE
  )
  concordance_pct <- 100 * n_concordant / n_total
  
  ggplot(
    pair_df,
    aes(
      x = x_NES,
      y = y_NES,
      color = direction_status
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      color = "grey70"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      color = "grey70"
    ) +
    geom_point(
      size = 2.4,
      alpha = 0.8
    ) +
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = sprintf(
        "Directionally concordant:\n%d/%d (%.0f%%)",
        n_concordant,
        n_total,
        concordance_pct
      ),
      hjust = -0.05,
      vjust = 1.1,
      size = 3.1,
      color = "grey20"
    ) +
    scale_color_manual(
      values = direction_cols,
      drop = FALSE,
      name = NULL
    ) +
    labs(
      x = pair$xlab,
      y = pair$ylab
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_bw(11) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(0.6, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    )
})

combined_gsea <- patchwork::wrap_plots(
  plot_list,
  ncol = 3,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Pairwise Hallmark NES-direction concordance",
    subtitle = paste0(
      "NES magnitudes are displayed descriptively; ",
      "cross-cohort interpretation focuses on direction because ",
      "ranking statistics differed between discovery and external cohorts."
    )
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    S5_DIR,
    "Figure_GSEA_pairwise_direction_concordance.pdf"
  ),
  combined_gsea,
  width = 15,
  height = 5.8,
  device = "pdf"
)

ggsave(
  file.path(
    S5_DIR,
    "Figure_GSEA_pairwise_direction_concordance.png"
  ),
  combined_gsea,
  width = 15,
  height = 5.8,
  dpi = 300
)

cat(
  paste0(
    "Saved: Figure_GSEA_pairwise_direction_concordance.pdf / .png\n"
  )
)



# ============================================================
# Figure: Cross-database TF target-set enrichment — FOXM1 and E2F4
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# ---- Load data ----
chea_2022_up <- read.csv(
  file.path(
    S5_DIR,
    "TF_ChEA_2022_Up.csv"
  ),
  stringsAsFactors = FALSE
)
encode_consensus <- read.csv(
  file.path(
    S5_DIR,
    "TF_ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_Up.csv"
  ),
  stringsAsFactors = FALSE
)

encode_2015 <- read.csv(
  file.path(
    S5_DIR,
    "TF_ENCODE_TF_ChIP-seq_2015_Up.csv"
  ),
  stringsAsFactors = FALSE
)

# ---- Extract TF symbol from Term ----
extract_tf <- function(term) {
  trimws(sub(" .*", "", term))
}

# ---- Build cross-database table for FOXM1 and E2F4 ----
build_tf_table <- function(tf_name, chea3, encode_c, encode_5) {
  rows <- list()
  
  chea3$TF <- extract_tf(chea3$Term)
  chea3_tf <- chea3[chea3$TF == tf_name, ]
  if (nrow(chea3_tf) > 0) {
    best <- chea3_tf[which.max(chea3_tf$Combined.Score), ]
    rows[["ChEA_2022"]] <- data.frame(
      TF = tf_name, Database = "ChEA_2022",
      Combined_Score = best$Combined.Score,
      Source = gsub(".*ChIP-\\w+ ", "", best$Term),
      stringsAsFactors = FALSE
    )
  }
  
  encode_c$TF <- extract_tf(encode_c$Term)
  encode_c_tf <- encode_c[encode_c$TF == tf_name, ]
  if (nrow(encode_c_tf) > 0) {
    best <- encode_c_tf[which.max(encode_c_tf$Combined.Score), ]
    rows[["ENCODE Consensus"]] <- data.frame(
      TF = tf_name, Database = "ENCODE Consensus",
      Combined_Score = best$Combined.Score,
      Source = "ENCODE Consensus",
      stringsAsFactors = FALSE
    )
  }
  
  encode_5$TF <- extract_tf(encode_5$Term)
  encode_5_tf <- encode_5[encode_5$TF == tf_name, ]
  if (nrow(encode_5_tf) > 0) {
    best <- encode_5_tf[which.max(encode_5_tf$Combined.Score), ]
    rows[["ENCODE 2015"]] <- data.frame(
      TF = tf_name, Database = "ENCODE 2015",
      Combined_Score = best$Combined.Score,
      Source = gsub(".*ChIP-seq ", "", best$Term),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, rows)
}

foxm1 <- build_tf_table("FOXM1", chea_2022_up, encode_consensus, encode_2015)
e2f4  <- build_tf_table("E2F4",  chea_2022_up, encode_consensus, encode_2015)

# ---- Colors ----
db_cols <- c(
  "ChEA_2022"         = "#4E79A7",
  "ENCODE Consensus" = "#E15759",
  "ENCODE 2015"      = "#B07AA1"
)
# ---- Panel A: FOXM1 ----
p_foxm1 <- ggplot(foxm1, aes(x = Database, y = Combined_Score, fill = Database)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f", Combined_Score)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = db_cols, guide = "none") +
  labs(title = "FOXM1", x = NULL, y = "Combined Score") +
  expand_limits(y = max(foxm1$Combined_Score) * 1.3) +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 10))

# ---- Panel B: E2F4 ----
p_e2f4 <- ggplot(e2f4, aes(x = Database, y = Combined_Score, fill = Database)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f", Combined_Score)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = db_cols, guide = "none") +
  labs(title = "E2F4", x = NULL, y = "Combined Score") +
  expand_limits(y = max(e2f4$Combined_Score) * 1.3) +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 10))

# ---- Combine ----
combined <- ggarrange(p_foxm1, p_e2f4, ncol = 2, nrow = 1)
combined <- annotate_figure(
  combined,
  top = text_grob(
    "Cross-database TF Target-Set Enrichment",
    face = "bold",
    size = 13
  )
)

# ---- Save ----
out_dir <- S5_DIR
ggsave(file.path(out_dir, "Figure_TF_cross_database.pdf"),
       combined, width = 10, height = 5, device = "pdf")
ggsave(file.path(out_dir, "Figure_TF_cross_database.png"),
       combined, width = 10, height = 5, dpi = 300)

cat("Saved: Figure_TF_cross_database.pdf / .png\n")
# ============================================================
# Four TF enrichment visualizations
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)

out_dir <- S5_DIR

# ---- Fig 1: ChEA_2022 Up (top 15 TFs) ----
chea_2022_up <- read.csv(file.path(out_dir, "TF_ChEA_2022_Up.csv"), stringsAsFactors = FALSE)
chea_2022_up$TF <- trimws(sub(" .*", "", chea_2022_up$Term))
chea_2022_up$neg_log10_p <- -log10(chea_2022_up$P.value)
chea_2022_up$Overlap_n <- as.numeric(sub("/.*", "", chea_2022_up$Overlap))

# Keep top entry per TF, take top 15 by Combined Score
chea_2022_up_top <- chea_2022_up %>%
  group_by(TF) %>%
  slice_max(Combined.Score, n = 1) %>%
  ungroup() %>%
  slice_max(Combined.Score, n = 15) %>%
  mutate(TF = factor(TF, levels = rev(TF)))

chea_2022_up_plot <- ggplot(chea_2022_up_top, aes(x = Combined.Score, y = TF)) +
  geom_point(aes(size = Overlap_n, color = neg_log10_p), alpha = 0.85) +
  scale_color_gradient(low = "#B07AA1", high = "#E15759", name = "-log10(P)") +
  scale_size_continuous(range = c(3, 10), name = "n Genes") +
  labs(title = "ChEA_2022: Up-regulated Genes (29 genes)",
       x = "Combined Score", y = NULL) +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(out_dir, "Figure_TF_ChEA2022_Up_dotplot.pdf"), chea_2022_up_plot, width = 8, height = 6, device = "pdf")
ggsave(file.path(out_dir, "Figure_TF_ChEA2022_Up_dotplot.png"), chea_2022_up_plot, width = 8, height = 6, dpi = 300)

# ---- Fig 2: ChEA_2022 Down (top 15 TFs) ----
chea_2022_down <- read.csv(
  file.path(out_dir, "TF_ChEA_2022_Down.csv"),
  stringsAsFactors = FALSE
)

chea_2022_down$TF <- trimws(
  sub(" .*", "", chea_2022_down$Term)
)

chea_2022_down$neg_log10_p <- -log10(
  pmax(
    chea_2022_down$P.value,
    .Machine$double.xmin
  )
)

chea_2022_down$Overlap_n <- as.numeric(
  sub("/.*", "", chea_2022_down$Overlap)
)

# Retain the highest-scoring term for each TF,
# followed by the 15 highest-scoring TFs.
chea_2022_down_top <- chea_2022_down %>%
  group_by(TF) %>%
  slice_max(
    order_by = Combined.Score,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  slice_max(
    order_by = Combined.Score,
    n = 15,
    with_ties = FALSE
  ) %>%
  arrange(Combined.Score) %>%
  mutate(
    TF = factor(TF, levels = TF)
  )

chea_2022_down_top$log10_CombinedScore <- log10(
  chea_2022_down_top$Combined.Score + 1
)

stopifnot(
  nrow(chea_2022_down_top) <= 15,
  !anyNA(chea_2022_down_top$TF),
  !anyNA(chea_2022_down_top$Combined.Score)
)

chea_2022_down_plot <- ggplot(chea_2022_down_top, aes(x = log10_CombinedScore, y = TF)) +
  geom_point(aes(size = Overlap_n, color = neg_log10_p), alpha = 0.85) +
  scale_color_gradient(low = "#4E79A7", high = "#E15759", name = "-log10(P)") +
  scale_size_continuous(range = c(3, 10), name = "n Genes") +
  labs(title = "ChEA_2022: Down-regulated Genes (3 genes)",
       x = "log10(Combined Score + 1)", y = NULL) +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold", size = 12))
ggsave(file.path(out_dir, "Figure_TF_ChEA2022_Down_dotplot.pdf"), chea_2022_down_plot, width = 8, height = 6, device = "pdf")
ggsave(file.path(out_dir, "Figure_TF_ChEA2022_Down_dotplot.png"), chea_2022_down_plot, width = 8, height = 6, dpi = 300)

# ---- Fig 3: ENCODE Consensus (top 15 TFs) ----
encode_consensus <- read.csv(
  file.path(out_dir, "TF_ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_Up.csv"),
  stringsAsFactors = FALSE)
encode_consensus$TF <- trimws(sub(" .*", "", encode_consensus$Term))
encode_consensus$neg_log10_p <- -log10(encode_consensus$P.value)
encode_consensus$Overlap_n <- as.numeric(sub("/.*", "", encode_consensus$Overlap))

encode_consensus_top <- encode_consensus %>%
  group_by(TF) %>%
  slice_max(Combined.Score, n = 1) %>%
  ungroup() %>%
  slice_max(Combined.Score, n = 15) %>%
  mutate(TF = factor(TF, levels = rev(TF)))

p3 <- ggplot(encode_consensus_top, aes(x = Combined.Score, y = TF)) +
  geom_point(aes(size = Overlap_n, color = neg_log10_p), alpha = 0.85) +
  scale_color_gradient(low = "#B07AA1", high = "#E15759", name = "-log10(P)") +
  scale_size_continuous(range = c(3, 10), name = "n Genes") +
  labs(title = "ENCODE & ChEA Consensus: Up-regulated Genes",
       x = "Combined Score", y = NULL) +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(out_dir, "Figure_TF_ENCODE_Consensus_dotplot.pdf"), p3, width = 8, height = 6, device = "pdf")
ggsave(file.path(out_dir, "Figure_TF_ENCODE_Consensus_dotplot.png"), p3, width = 8, height = 6, dpi = 300)

# ---- Fig 4: ENCODE 2015 (top 15 TFs) ----
encode_2015 <- read.csv(
  file.path(out_dir, "TF_ENCODE_TF_ChIP-seq_2015_Up.csv"),
  stringsAsFactors = FALSE)
encode_2015$TF <- trimws(sub(" .*", "", encode_2015$Term))
encode_2015$neg_log10_p <- -log10(encode_2015$P.value)
encode_2015$Overlap_n <- as.numeric(sub("/.*", "", encode_2015$Overlap))

encode_2015_top <- encode_2015 %>%
  group_by(TF) %>%
  slice_max(Combined.Score, n = 1) %>%
  ungroup() %>%
  slice_max(Combined.Score, n = 15) %>%
  mutate(TF = factor(TF, levels = rev(TF)))

p4 <- ggplot(encode_2015_top, aes(x = Combined.Score, y = TF)) +
  geom_point(aes(size = Overlap_n, color = neg_log10_p), alpha = 0.85) +
  scale_color_gradient(low = "#B07AA1", high = "#E15759", name = "-log10(P)") +
  scale_size_continuous(range = c(3, 10), name = "n Genes") +
  labs(title = "ENCODE TF ChIP-seq 2015: Up-regulated Genes",
       x = "Combined Score", y = NULL) +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(out_dir, "Figure_TF_ENCODE_2015_dotplot.pdf"), p4, width = 8, height = 6, device = "pdf")
ggsave(file.path(out_dir, "Figure_TF_ENCODE_2015_dotplot.png"), p4, width = 8, height = 6, dpi = 300)

# ---- Done ----
cat("Saved 8 files (4 PDF + 4 PNG):\n",
    "  Figure_TF_ChEA2022_Up_dotplot\n",
    "  Figure_TF_ChEA2022_Down_dotplot\n",
    "  Figure_TF_ENCODE_Consensus_dotplot\n",
    "  Figure_TF_ENCODE_2015_dotplot\n")


# ============================================================
# 验证队列 logFC 散点图（最终版）
# 30/30 方向一致 | 按阈值分类着色 | 标 top 10
# ============================================================
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

# ============================================================
# Figure 4B — External effect-direction assessment
# ============================================================
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

signature_genes <- unique(as.character(fisher_out$gene))
genes_evaluable <- intersect(signature_genes, res_val$gene)
genes_missing <- setdiff(signature_genes, genes_evaluable)

plot_df <- data.frame(
  gene = genes_evaluable,
  logFC_disc = fisher_out$re_estimate[
    match(genes_evaluable, fisher_out$gene)
  ],
  logFC_val = res_val$logFC[
    match(genes_evaluable, res_val$gene)
  ],
  external_t = res_val$t[
    match(genes_evaluable, res_val$gene)
  ],
  external_P = res_val$P.Value[
    match(genes_evaluable, res_val$gene)
  ],
  external_FDR = res_val$adj.P.Val[
    match(genes_evaluable, res_val$gene)
  ],
  discovery_direction = fisher_out$direction[
    match(genes_evaluable, fisher_out$gene)
  ],
  stringsAsFactors = FALSE
)

# Reproducibility checks
stopifnot(
  length(signature_genes) == 32L,
  nrow(plot_df) == 30L,
  setequal(
    genes_missing,
    c("HIST1H3B", "HIST1H2BM")
  ),
  !anyNA(plot_df$logFC_disc),
  !anyNA(plot_df$logFC_val)
)

# Functional modules
plot_df$module <- ifelse(
  plot_df$gene %in% names(module_map),
  as.character(module_map[plot_df$gene]),
  "Unclassified"
)

plot_df$module <- factor(
  plot_df$module,
  levels = MODULE_LEVELS
)

# Directional concordance
plot_df$directionally_concordant <-
  sign(plot_df$logFC_disc) == sign(plot_df$logFC_val)

n_total <- length(signature_genes)
n_evaluable <- nrow(plot_df)
n_concordant <- sum(
  plot_df$directionally_concordant,
  na.rm = TRUE
)

# External effect-size display threshold
EXTERNAL_LFC_DISPLAY_CUT <- 0.5

plot_df$external_above_threshold <-
  abs(plot_df$logFC_val) > EXTERNAL_LFC_DISPLAY_CUT

n_above_threshold <- sum(
  plot_df$external_above_threshold,
  na.rm = TRUE
)

class_levels <- c(
  sprintf(
    "External |log2FC| > 0.5 (n=%d)",
    n_above_threshold
  ),
  sprintf(
    "External |log2FC| <= 0.5 (n=%d)",
    n_evaluable - n_above_threshold
  )
)

plot_df$class <- ifelse(
  plot_df$external_above_threshold,
  class_levels[1],
  class_levels[2]
)

plot_df$class <- factor(
  plot_df$class,
  levels = class_levels
)

class_colors <- setNames(
  c("#E15759", "#4E79A7"),
  class_levels
)

# Label the 10 genes with the largest external absolute log2FC
# Label all genes exceeding the external effect-size display threshold.
plot_df$label <- ifelse(
  plot_df$external_above_threshold,
  plot_df$gene,
  ""
)
# Descriptive association statistics
r_pearson <- cor(
  plot_df$logFC_disc,
  plot_df$logFC_val,
  method = "pearson",
  use = "complete.obs"
)

rho_spearman <- cor(
  plot_df$logFC_disc,
  plot_df$logFC_val,
  method = "spearman",
  use = "complete.obs"
)

ols_fit <- lm(
  logFC_val ~ logFC_disc,
  data = plot_df
)

ols_intercept <- unname(coef(ols_fit)[1])
ols_slope <- unname(coef(ols_fit)[2])
# Descriptive regression confidence intervals
ols_ci <- confint(
  ols_fit,
  level = 0.95
)

# Correlation sensitivity analyses
pearson_without_cx3cr1 <- with(
  subset(plot_df, gene != "CX3CR1"),
  cor(
    logFC_disc,
    logFC_val,
    method = "pearson",
    use = "complete.obs"
  )
)

pearson_upregulated_only <- with(
  subset(
    plot_df,
    discovery_direction == "Up"
  ),
  cor(
    logFC_disc,
    logFC_val,
    method = "pearson",
    use = "complete.obs"
  )
)

loo_pearson <- vapply(
  seq_len(nrow(plot_df)),
  function(i) {
    cor(
      plot_df$logFC_disc[-i],
      plot_df$logFC_val[-i],
      method = "pearson",
      use = "complete.obs"
    )
  },
  numeric(1)
)

correlation_sensitivity <- data.frame(
  analysis = c(
    "Pearson correlation, all evaluable genes",
    "Spearman correlation, all evaluable genes",
    "Pearson correlation excluding CX3CR1",
    "Pearson correlation, upregulated genes only",
    "Leave-one-gene-out Pearson minimum",
    "Leave-one-gene-out Pearson maximum",
    "Descriptive OLS slope",
    "Descriptive OLS slope lower 95% CI",
    "Descriptive OLS slope upper 95% CI",
    "Descriptive OLS intercept"
  ),
  n_genes = c(
    nrow(plot_df),
    nrow(plot_df),
    nrow(plot_df) - 1L,
    sum(
      plot_df$discovery_direction == "Up"
    ),
    nrow(plot_df) - 1L,
    nrow(plot_df) - 1L,
    nrow(plot_df),
    nrow(plot_df),
    nrow(plot_df),
    nrow(plot_df)
  ),
  estimate = c(
    r_pearson,
    rho_spearman,
    pearson_without_cx3cr1,
    pearson_upregulated_only,
    min(loo_pearson, na.rm = TRUE),
    max(loo_pearson, na.rm = TRUE),
    ols_slope,
    ols_ci["logFC_disc", 1],
    ols_ci["logFC_disc", 2],
    ols_intercept
  ),
  stringsAsFactors = FALSE
)

write.csv(
  correlation_sensitivity,
  file.path(
    S4_DIR,
    "External_effect_correlation_sensitivity.csv"
  ),
  row.names = FALSE
)

# Identical x- and y-axis limits are required for meaningful
# visual comparison with the y = x identity line.
common_limits <- range(
  c(
    plot_df$logFC_disc,
    plot_df$logFC_val
  ),
  finite = TRUE
)

limit_padding <- max(
  diff(common_limits) * 0.08,
  0.05
)

common_limits <- c(
  common_limits[1] - limit_padding,
  common_limits[2] + limit_padding
)

p_validation_scatter <- ggplot(
  plot_df,
  aes(
    x = logFC_disc,
    y = logFC_val
  )
) +
  # Zero-effect reference lines
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    color = "grey55"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey55"
  ) +
  # External magnitude threshold
  geom_hline(
    yintercept = c(
      -EXTERNAL_LFC_DISPLAY_CUT,
      EXTERNAL_LFC_DISPLAY_CUT
    ),
    linetype = "dotted",
    linewidth = 0.35,
    color = "grey72"
  ) +
  # Absolute-agreement reference
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "solid",
    linewidth = 0.5,
    color = "grey60"
  ) +
  # Descriptive OLS fit and confidence band
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey20",
    fill = "grey65",
    linewidth = 0.65,
    linetype = "dashed",
    alpha = 0.18
  ) +
  # Fixed point size avoids redundant encoding of the x-axis
  geom_point(
    aes(fill = class),
    shape = 21,
    size = 3.8,
    alpha = 0.9,
    color = "white",
    stroke = 0.3
  ) +
  geom_text_repel(
    aes(label = label),
    seed = 2024,
    size = 3.0,
    color = "grey15",
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.color = "grey50",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = class_colors,
    drop = FALSE
  ) +
  coord_equal(
    xlim = common_limits,
    ylim = common_limits,
    clip = "off"
  ) +
  labs(
    title = paste(
      "External assessment of 30 evaluable genes",
      "from the 32-gene signature"
    ),
    subtitle = sprintf(
      "%d of %d genes evaluable; all %d directionally concordant | Pearson r = %.3f",
      n_evaluable,
      n_total,
      n_concordant,
      r_pearson
    ),
    x = paste0(
      "Discovery random-effects estimate\n",
      "(log\u2082FC, NS vs S)"
    ),
    y = paste0(
      "GSE65682 estimate\n",
      "(log\u2082FC, NS vs S)"
    ),
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    plot.subtitle = element_text(size = 9.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(
      t = 7,
      r = 20,
      b = 7,
      l = 7
    )
  )

pdf(
  file.path(
    S4_DIR,
    "Figure_validation_logFC_scatter.pdf"
  ),
  width = 8,
  height = 7
)
print(p_validation_scatter)
dev.off()

png(
  file.path(
    S4_DIR,
    "Figure_validation_logFC_scatter.png"
  ),
  width = 8,
  height = 7,
  units = "in",
  res = 300
)
print(p_validation_scatter)
dev.off()

# Export source data for the panel and supplementary reporting
write.csv(
  plot_df,
  file.path(
    S4_DIR,
    "Figure_validation_logFC_scatter_source_data.csv"
  ),
  row.names = FALSE
)

cat(sprintf(
  paste0(
    "\nExternal logFC scatter:\n",
    "  Signature genes: %d\n",
    "  Evaluable genes: %d\n",
    "  Missing genes: %s\n",
    "  Directionally concordant: %d/%d\n",
    "  External |log2FC| > 0.5: %d/%d\n",
    "  Pearson r: %.3f\n",
    "  Spearman rho: %.3f\n",
    "  OLS intercept: %.3f\n",
    "  OLS slope: %.3f\n",
    "  External nominal P < 0.05: %d/%d\n",
    "  External FDR < 0.05: %d/%d\n"
  ),
  n_total,
  n_evaluable,
  paste(genes_missing, collapse = ", "),
  n_concordant,
  n_evaluable,
  n_above_threshold,
  n_evaluable,
  r_pearson,
  rho_spearman,
  ols_intercept,
  ols_slope,
  sum(plot_df$external_P < 0.05, na.rm = TRUE),
  n_evaluable,
  sum(plot_df$external_FDR < 0.05, na.rm = TRUE),
  n_evaluable
))




# ============================================================
# Figure 1: xCell 三队列全细胞类型 delta 热图
# ============================================================
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))

xc_all <- read.csv(file.path(S5_DIR, "xCell_triple_cohort.csv"), stringsAsFactors = FALSE)

# 取 delta 列构建矩阵
delta_cols <- grep("^delta_", names(xc_all), value = TRUE)
xc_mat <- as.matrix(xc_all[, delta_cols])
rownames(xc_mat) <- xc_all$cell_type
colnames(xc_mat) <- gsub("delta_", "", delta_cols)

# 按 mean_abs_delta 排序
ord <- order(-xc_all$mean_abs_delta)
xc_mat <- xc_mat[ord, , drop = FALSE]

# 标记三队列方向一致性
all_same <- xc_all$all_same_dir[ord]
rownames(xc_mat) <- ifelse(all_same,
                           paste0(rownames(xc_mat), "  yes"),
                           rownames(xc_mat))

# 列注释：一致性
row_annot <- data.frame(
  Concordance = ifelse(all_same, "All 3 same direction", "Discordant"),
  row.names = rownames(xc_mat)
)
ann_colors <- list(
  Concordance = c("All 3 same direction" = "#4E79A7", "Discordant" = "grey80")
)

max_abs <- max(abs(xc_mat), na.rm = TRUE)
palette <- colorRampPalette(c("#4E79A7", "white", "#E15759"))(100)
breaks  <- seq(-max_abs, max_abs, length.out = 101)

pdf(file.path(S5_DIR, "Figure_xCell_all_types_delta_heatmap.pdf"), 10, 16)
pheatmap(xc_mat, color = palette, breaks = breaks,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_annot,
         annotation_colors = ann_colors,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 5, fontsize_row = 6, fontsize_col = 10,
         angle_col = 0, border_color = NA,
         main = "xCell: All 67 Cell Types - Delta (NS - S) Across 3 Cohorts",
         na_col = "grey90")
dev.off()
png(file.path(S5_DIR, "Figure_xCell_all_types_delta_heatmap.png"), 10, 16, units = "in", res = 300)
pheatmap(xc_mat, color = palette, breaks = breaks,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_annot,
         annotation_colors = ann_colors,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 5, fontsize_row = 6, fontsize_col = 10,
         angle_col = 0, border_color = NA,
         main = "xCell: All 67 Cell Types - Delta (NS - S) Across 3 Cohorts",
         na_col = "grey90")
dev.off()
cat("xCell all-types heatmap saved.\n")

# ============================================================
# Figure: Exploratory overlap across three Enrichr query analyses
# This visualization describes cross-analysis overlap and does
# not represent sequential experimental validation.
# ============================================================
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggalluvial))


# ============================================================
# Nature 配色方案
# ============================================================
NATURE_RED    <- "#D43F3A"   # 明亮砖红
NATURE_BLUE   <- "#3A6EA5"   # 明亮钢蓝
NATURE_PURPLE <- "#8B5E9E"   # 明亮紫
NATURE_GREY   <- "#E8E3DC"   # 暖灰
NATURE_BG     <- "#FAFAF8"

# DGIdb 暖色系（浅黄 / 浅绿）
NATURE_GREEN       <- "#7A9B3C"   # 明亮鼠尾草绿
NATURE_GREEN_LIGHT <- "#A3C45A"   # 浅绿
NATURE_GOLD        <- "#D4A830"   # 明亮琥珀金
NATURE_GOLD_LIGHT  <- "#E5D28A"   # 浅麦秆黄

theme_nature <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      panel.background  = element_rect(fill = NATURE_BG, color = NA),
      plot.background   = element_rect(fill = NATURE_BG, color = NA),
      plot.title        = element_text(size = 13, face = "bold", hjust = 0.5,
                                       color = "#222222", margin = margin(b = 6)),
      plot.subtitle     = element_text(size = 10, hjust = 0.5,
                                       color = "#666666", margin = margin(b = 12)),
      axis.text.x       = element_text(size = 10, color = "#333333"),
      axis.text.y       = element_blank(),
      axis.title.y      = element_text(size = 10, color = "#555555"),
      legend.position   = "none"
    )
}

out_dir <- S6_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Enrichr alluvial plot: overlap across three exploratory analyses
# ============================================================

# ---- 构建排名表 ----
# S1: Full signature
s1_uniq <- top_enrichr_by_drug(drug_all)
s1_uniq <- s1_uniq[order(-s1_uniq$Combined.Score), ]
s1_uniq$rank <- seq_len(nrow(s1_uniq))
s1_lookup_rank  <- setNames(s1_uniq$rank, norm_drug(s1_uniq$drug_name))

# S2: Non-cycle
s2_uniq <- top_enrichr_by_drug(drug_nc_all)
s2_uniq <- s2_uniq[order(-s2_uniq$Combined.Score), ]
s2_uniq$rank <- seq_len(nrow(s2_uniq))
s2_lookup_rank  <- setNames(s2_uniq$rank, norm_drug(s2_uniq$drug_name))

# S3: Bidirectional
s3_lookup_pass  <- setNames(bi_summary$pass_threshold, norm_drug(bi_summary$drug))

# 全药物
all_drugs <- unique(c(names(s1_lookup_rank), names(s2_lookup_rank), names(s3_lookup_pass)))

# 构建 flow
flow <- data.frame(
  drug = all_drugs,
  s1_rank = s1_lookup_rank[all_drugs],
  s2_rank = s2_lookup_rank[all_drugs],
  s3_pass = s3_lookup_pass[all_drugs],
  stringsAsFactors = FALSE
)

# 三阶段分级
flow$S1 <- ifelse(
  !is.na(flow$s1_rank) & flow$s1_rank <= 20,
  "Top 20",
  ifelse(
    !is.na(flow$s1_rank),
    "Other full-signature match",
    "Not identified"
  )
)

flow$S2 <- ifelse(
  !is.na(flow$s2_rank) & flow$s2_rank <= 15,
  "Top 15",
  ifelse(
    !is.na(flow$s2_rank),
    "Other non-cycle-set match",
    "Not identified"
  )
)
flow$S3 <- ifelse(
  !is.na(flow$s3_pass) & flow$s3_pass,
  "CS > 30 in both arms",
  ifelse(
    !is.na(flow$s3_pass),
    "CS <= 30 in at least one arm",
    "Not identified"
  )
)
# 因子
status_levels <- c(
  "Top 20",
  "Top 15",
  "CS > 30 in both arms",
  "Other full-signature match",
  "Other non-cycle-set match",
  "CS <= 30 in at least one arm",
  "Not identified"
)

flow$S1 <- factor(
  flow$S1,
  levels = c(
    "Top 20",
    "Other full-signature match",
    "Not identified"
  )
)

flow$S2 <- factor(
  flow$S2,
  levels = c(
    "Top 15",
    "Other non-cycle-set match",
    "Not identified"
  )
)

flow$S3 <- factor(
  flow$S3,
  levels = c(
    "CS > 30 in both arms",
    "CS <= 30 in at least one arm",
    "Not identified"
  )
)
# 长格式
flow_long <- data.frame(
  drug = rep(flow$drug, 3),
  Stage = rep(
    c(
      "Stage 1: Full signature\n(32 genes)",
      "Stage 2: Non-cycle set\n(21 genes)",
      "Stage 3: Dual-query overlap"
    ),
    each = nrow(flow)
  ),
  Status = factor(
    c(
      as.character(flow$S1),
      as.character(flow$S2),
      as.character(flow$S3)
    ),
    levels = status_levels
  ),
  stringsAsFactors = FALSE
)

flow_long$Stage <- factor(
  flow_long$Stage,
  levels = c(
    "Stage 1: Full signature\n(32 genes)",
    "Stage 2: Non-cycle set\n(21 genes)",
    "Stage 3: Dual-query overlap"
  )
)

# Nature 分层配色
status_cols_enrichr <- c(
  "Top 20" = NATURE_RED,
  "Top 15" = "#E0736E",
  "CS > 30 in both arms" = NATURE_RED,
  "Other full-signature match" = NATURE_BLUE,
  "Other non-cycle-set match" = "#6B9CC8",
  "CS <= 30 in at least one arm" = NATURE_PURPLE,
  "Not identified" = NATURE_GREY
)

n_s1 <- sum(!is.na(flow$s1_rank))
n_s2 <- sum(!is.na(flow$s2_rank))
n_s3 <- sum(!is.na(flow$s3_pass))
n_pass <- sum(flow$s3_pass, na.rm = TRUE)

p_enrichr <- ggplot(flow_long,
                    aes(x = Stage, stratum = Status, alluvium = drug, fill = Status)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "white", width = 0.35, alpha = 0.55) +
  geom_stratum(width = 0.30, color = "white", linewidth = 0.4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.2, color = "#333333", fontface = "bold") +
  scale_fill_manual(values = status_cols_enrichr) +
  labs(
    title = "Overlap Across Exploratory Enrichr Analyses",
    subtitle = sprintf(
      "%d full-signature matches | %d non-cycle-set matches | %d dual-query overlaps; %d met CS > 30 in both arms",
      n_s1,
      n_s2,
      n_s3,
      n_pass
    ),
    x = NULL,
    y = "Number of compounds"
  ) +
  theme_nature()

pdf(file.path(out_dir, "Figure_Sankey_Enrichr.pdf"), 11, 7)
print(p_enrichr)
dev.off()
png(file.path(out_dir, "Figure_Sankey_Enrichr.png"), 11, 7, units = "in", res = 300)
print(p_enrichr)
dev.off()


if (!is.null(gene_summary) && nrow(gene_summary) > 0) {
  s2_rows <- do.call(rbind, lapply(seq_len(nrow(gene_summary)), function(i) {
    drugs <- trimws(unlist(strsplit(gene_summary$drug[i], ";")))
    drugs <- drugs[drugs != ""]
    if (length(drugs) == 0) return(NULL)
    data.frame(gene = gene_summary$gene[i], drug = drugs, stringsAsFactors = FALSE)
  }))
  
  write.csv(
    s2_rows,
    file.path(
      S6_DIR,
      "DGIdb_32genes_gene_drug_interactions.csv"
    ),
    row.names = FALSE
  )
  
  cat(sprintf(
    "DGIdb interaction export: %d gene–compound pairs, %d genes, %d compounds\n",
    nrow(s2_rows),
    length(unique(s2_rows$gene)),
    length(unique(s2_rows$drug))
  ))
} else {
  cat("Table S2 skipped: no DGIdb gene-level records were available.\n")
}

# ============================================================
# Final run metadata and session information
# ============================================================

RUN_FINISHED_AT <- Sys.time()

# ------------------------------------------------------------
# A. Input-file manifest
# ------------------------------------------------------------

if (
  exists("required_input_files") &&
  length(required_input_files) > 0
) {
  
  input_file_info <- file.info(
    required_input_files
  )
  
  input_file_manifest <- data.frame(
    file = basename(required_input_files),
    relative_location = file.path(
      "data",
      basename(required_input_files)
    ),
    size_bytes = input_file_info$size,
    last_modified = format(
      input_file_info$mtime,
      "%Y-%m-%d %H:%M:%S"
    ),
    md5 = unname(
      tools::md5sum(required_input_files)
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    input_file_manifest,
    file.path(
      S7_DIR,
      "input_file_manifest.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# B. Analysis-run manifest
# ------------------------------------------------------------

analysis_run_manifest <- data.frame(
  item = c(
    "run_finished_at",
    "R_version",
    "platform",
    "log2FC_threshold",
    "FDR_threshold",
    "individual_cohort_nominal_P_threshold",
    "random_seed_sign_concordance",
    "outcome_label_permutations",
    "random_seed_PC1_common_pool",
    "random_seed_PC1_full_platform_pool",
    "random_gene_sets_per_pool"
  ),
  value = c(
    format(
      RUN_FINISHED_AT,
      "%Y-%m-%d %H:%M:%S %Z"
    ),
    R.version.string,
    R.version$platform,
    as.character(LFC_CUT),
    as.character(FDR_CUT),
    if (
      exists("INDIVIDUAL_P_CUT")
    ) {
      as.character(INDIVIDUAL_P_CUT)
    } else {
      NA_character_
    },
    "2024",
    if (
      exists("n_perm")
    ) {
      as.character(n_perm)
    } else {
      NA_character_
    },
    "2024",
    "2025",
    if (
      exists("n_random_sets")
    ) {
      as.character(n_random_sets)
    } else {
      NA_character_
    }
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_run_manifest,
  file.path(
    S7_DIR,
    "analysis_run_manifest.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# C. Complete R session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    S7_DIR,
    "sessionInfo.txt"
  )
)

# ------------------------------------------------------------
# D. Warnings generated during the run
# ------------------------------------------------------------

run_warnings <- warnings()

if (!is.null(run_warnings)) {
  
  writeLines(
    capture.output(run_warnings),
    con = file.path(
      S7_DIR,
      "warnings.txt"
    )
  )
  
} else {
  
  writeLines(
    "No warnings were retained at the end of the R session.",
    con = file.path(
      S7_DIR,
      "warnings.txt"
    )
  )
}

cat(sprintf(
  paste0(
    "\n====== ANALYSIS COMPLETE ======\n",
    "Finished: %s\n",
    "Session information: %s\n",
    "Input manifest: %s\n",
    "Run manifest: %s\n"
  ),
  format(
    RUN_FINISHED_AT,
    "%Y-%m-%d %H:%M:%S %Z"
  ),
  file.path(S7_DIR, "sessionInfo.txt"),
  file.path(S7_DIR, "input_file_manifest.csv"),
  file.path(S7_DIR, "analysis_run_manifest.csv")
))