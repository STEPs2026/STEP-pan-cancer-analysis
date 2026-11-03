# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

BASE_DIR <- getwd()

EXP_DIR <- file.path(BASE_DIR, "TCGA_201906cBioP", "exp_gct")
STEP_GMT <- file.path(BASE_DIR, "gene set", "STEP_all_human_gene_sets_merged.gmt")

HALLMARK_GMT <- NA_character_

STEP_SET_NAME <- "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"
CAF3_SET_NAME  <- "CAF_3gene_score"

HALLMARK_SCORE_METHOD <- "ssGSEA"

RUN_HALLMARK_ZMEAN_SENSITIVITY <- FALSE

OPTIONAL_ENRICHMENT_GMTS <- c(
  REACTOME = file.path(BASE_DIR, "gene set", "c2.cp.reactome.Hs.symbols.gmt"),
  GOBP     = file.path(BASE_DIR, "gene set", "c5.go.bp.Hs.symbols.gmt")
)

RUN_GENE_LEVEL_CORRELATION <- TRUE
RUN_PER_CANCER_GSEA        <- TRUE
RUN_PAN_CANCER_GSEA        <- TRUE
RUN_ORA                     <- TRUE

MIN_MATCHED_STEP_GENES    <- 1L
MIN_MATCHED_CAF3_GENES     <- 1L
MIN_MATCHED_HALLMARK_GENES <- 5L
MIN_SAMPLES_PER_CANCER     <- 10L
MIN_COR_PAIRS              <- 10L
MIN_GROUP_N                <- 5L
GSEA_MIN_SIZE              <- 10L
GSEA_MAX_SIZE              <- 500L

MIN_CANCERS_FOR_GENE_META      <- 20L
MIN_CANCER_COVERAGE_FRACTION   <- 0.80

PAN_GSEA_PRIMARY_META_MODEL  <- "random"    
PAN_GSEA_PRIMARY_RANK_METRIC <- "meta_zstat" 

RUN_FIXED_EFFECT_META_Z_SENSITIVITY <- TRUE
RUN_META_RHO_RANK_SENSITIVITY       <- TRUE
RUN_DEOVERLAP_SENSITIVITY           <- TRUE
REMOVE_CAF3_GENES_IN_DEOVERLAP      <- TRUE
PLOT_SENSITIVITY_CURVES             <- FALSE

RUN_LEAVE_ONE_CANCER_OUT <- FALSE
LOCO_COLLECTIONS <- c("HALLMARK")
LOCO_CORE_PATHWAYS <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING"
)

GSEA_CURVE_FDR_CUTOFF    <- 0.05
GSEA_DOTPLOT_FDR_CUTOFF  <- 0.05
GSEA_CURVE_MIN_ABS_NES   <- 0.00  
TOP_CURVES_POSITIVE      <- 6L
TOP_CURVES_NEGATIVE      <- 6L
TOP_DOTPLOT_PER_DIRECTION <- 12L
SAVE_INDIVIDUAL_CLASSIC_CURVES <- TRUE
REACTOME_COLLAPSE_FOR_PLOTS     <- TRUE

ORA_META_FDR_CUTOFF <- 0.05
ORA_ABS_RHO_CUTOFF  <- 0.10
ORA_PLOT_FDR_CUTOFF <- 0.05

RUN_GENE_CORE_INTEGRATION <- TRUE
GENE_CORE_META_FDR_CUTOFF <- 0.05
GENE_CORE_MIN_COVERAGE_FRACTION <- 0.80
GENE_CORE_MIN_DIRECTION_CONSISTENCY <- 0.80

GENE_CORE_STRICT_MAX_I2 <- 75
GENE_CORE_MIN_HALLMARK_LE_CANCERS <- 5L
GENE_CORE_MIN_EMT_LE_CANCERS <- 3L
GENE_CORE_TOP_PER_DIRECTION <- 20L
GENE_CORE_SHOW_HEATMAP_STARS <- TRUE
GENE_CORE_MAKE_FOREST_PLOTS <- TRUE
GENE_CORE_FOREST_PER_DIRECTION <- 5L
EMT_LEADING_EDGE_TOP_N <- 20L
PATHWAY_LEADING_EDGE_TOP_N <- 30L

EMT_HALLMARK_PATHWAY <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
CORE_HALLMARK_PATHWAYS <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING"
)

GENE_FUNCTIONAL_MODULES <- list(
  EMT_CYTOSKELETON_REMODELING = c(
    "LGALS1", "PFN1", "CFL1", "TMSB10", "TIMP1", "PLAUR",
    "SERPINE1", "EMP3", "NNMT", "LOXL1", "IL32", "TNFRSF12A"
  ),
  METABOLIC_OXIDATIVE_STRESS = c(
    "GPX4", "HMOX1", "LDHA", "SOD1", "ECH1", "ATOX1", "MRPL52"
  ),
  IMMUNE_INTERFERON = c(
    "IRF1", "IRF7", "CXCL10", "IFI30", "ISG20", "SOCS1", "IL6", "CCL2"
  ),
  STRESS_PROTEOSTASIS_CELL_CYCLE = c(
    "CDKN1A", "PPIA", "PPIB", "NEDD8", "RIPK2", "CCND3"
  )
)
USER_FUNCTIONAL_CORE_GENES <- unique(unlist(GENE_FUNCTIONAL_MODULES, use.names = FALSE))

LOW_COLOR  <- "#4DBBD5"   
HIGH_COLOR <- "#E64B35"   
NEG_COLOR  <- "#3C5488"   
ZERO_COLOR <- "#FFFFFF"
POS_COLOR  <- "#DC0000"   
GROUP_COLORS <- c("Low" = LOW_COLOR, "High" = HIGH_COLOR)

HEATMAP_CLUSTER_ROWS <- TRUE
HEATMAP_CLUSTER_COLUMNS <- TRUE
HEATMAP_CLUSTER_DISTANCE <- "euclidean"   
HEATMAP_CLUSTER_METHOD <- "complete"      
HEATMAP_CLUSTER_NA_FILL <- 0               
HEATMAP_EXPORT_CLUSTER_ORDER <- TRUE

OUT_DIR <- file.path(
  BASE_DIR,
  paste0("TCGA_STEP_Hallmark_CAF3_Zmean_ssGSEA_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

cran_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "stringr", "purrr",
  "tibble", "matrixStats", "scales", "ggrepel", "patchwork"
)
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c("GSVA", "BiocParallel", "fgsea", "AnnotationDbi", "org.Hs.eg.db")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(tibble)
  library(matrixStats)
  library(scales)
  library(ggrepel)
  library(patchwork)
})

subdirs <- c(
  "00_logs",
  "01_gene_set_QC",
  "02_scores_by_cancer",
  "03_CAF3_residualization",
  "04_Hallmark_correlations",
  "05_Hallmark_heatmaps",
  "06_Hallmark_meta_summary",
  "07_Hallmark_high_low",
  "08_method_concordance",
  "09_gene_correlations",
  "10_gene_meta_analysis",
  "11_GSEA",
  "12_ORA",
  "13_leading_edge",
  "14_figure_ready",
  "15_sensitivity",
  "16_gene_core_integration"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (d in subdirs) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "00_logs", "run_log.txt")
cat("Run started: ", as.character(Sys.time()), "\n", file = LOG_FILE)

msg <- function(...) {
  txt <- paste0(...)
  message(txt)
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
}

normalize_tcga_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  x <- gsub("_", "-", x)
  trimws(x)
}

sample15 <- function(x) substr(normalize_tcga_id(x), 1, 15)
patient12 <- function(x) substr(normalize_tcga_id(x), 1, 12)
sample_type_code <- function(x) substr(normalize_tcga_id(x), 14, 15)

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

fmt_p <- function(p, prefix = "p = ") {
  ifelse(
    is.na(p) | !is.finite(p), paste0(prefix, "NA"),
    ifelse(p < 1e-4, paste0(prefix, formatC(p, format = "e", digits = 2)),
           paste0(prefix, signif(p, 3)))
  )
}

p_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

clean_pathway_name <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  stringr::str_to_title(tolower(x))
}

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8, color = "black"),
      strip.text = element_text(size = 8, face = "bold"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
}

save_plot <- function(p, stem, width = 10, height = 7, dpi = 600) {
  print(p)
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

cluster_long_heatmap <- function(df, row_var, col_var, value_var,
                                 default_row_order = NULL,
                                 default_col_order = NULL,
                                 stem = NULL) {
  x <- df %>%
    transmute(
      .row = as.character(.data[[row_var]]),
      .col = as.character(.data[[col_var]]),
      .value = suppressWarnings(as.numeric(.data[[value_var]]))
    ) %>%
    filter(!is.na(.row), nzchar(.row), !is.na(.col), nzchar(.col)) %>%
    group_by(.row, .col) %>%
    summarise(
      .value = if (all(!is.finite(.value))) NA_real_ else mean(.value[is.finite(.value)], na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(x) == 0) {
    return(list(
      data = df,
      row_order = character(0),
      col_order = character(0)
    ))
  }

  observed_rows <- unique(x$.row)
  observed_cols <- unique(x$.col)

  if (is.null(default_row_order)) default_row_order <- sort(observed_rows)
  if (is.null(default_col_order)) default_col_order <- sort(observed_cols)

  default_row_order <- unique(c(
    intersect(as.character(default_row_order), observed_rows),
    setdiff(observed_rows, as.character(default_row_order))
  ))
  default_col_order <- unique(c(
    intersect(as.character(default_col_order), observed_cols),
    setdiff(observed_cols, as.character(default_col_order))
  ))

  mat <- matrix(
    NA_real_,
    nrow = length(default_row_order),
    ncol = length(default_col_order),
    dimnames = list(default_row_order, default_col_order)
  )
  rr <- match(x$.row, rownames(mat))
  cc <- match(x$.col, colnames(mat))
  mat[cbind(rr, cc)] <- x$.value

  mat_cluster <- mat
  mat_cluster[!is.finite(mat_cluster)] <- HEATMAP_CLUSTER_NA_FILL

  row_order <- default_row_order
  col_order <- default_col_order

  if (isTRUE(HEATMAP_CLUSTER_ROWS) && nrow(mat_cluster) > 1) {
    row_ok <- tryCatch({
      d <- stats::dist(mat_cluster, method = HEATMAP_CLUSTER_DISTANCE)
      hc <- stats::hclust(d, method = HEATMAP_CLUSTER_METHOD)
      rownames(mat_cluster)[hc$order]
    }, error = function(e) {
      warning("Row clustering failed; original order retained. Error: ", conditionMessage(e))
      default_row_order
    })
    row_order <- row_ok
  }

  if (isTRUE(HEATMAP_CLUSTER_COLUMNS) && ncol(mat_cluster) > 1) {
    col_ok <- tryCatch({
      d <- stats::dist(t(mat_cluster), method = HEATMAP_CLUSTER_DISTANCE)
      hc <- stats::hclust(d, method = HEATMAP_CLUSTER_METHOD)
      colnames(mat_cluster)[hc$order]
    }, error = function(e) {
      warning("Column clustering failed; original order retained. Error: ", conditionMessage(e))
      default_col_order
    })
    col_order <- col_ok
  }

  out_df <- df

  out_df[[row_var]] <- factor(as.character(out_df[[row_var]]), levels = rev(row_order))
  out_df[[col_var]] <- factor(as.character(out_df[[col_var]]), levels = col_order)

  if (isTRUE(HEATMAP_EXPORT_CLUSTER_ORDER) && !is.null(stem) && nzchar(stem)) {
    data.table::fwrite(
      data.frame(position_top_to_bottom = seq_along(row_order), item = row_order),
      paste0(stem, "_clustered_row_order.csv")
    )
    data.table::fwrite(
      data.frame(position_left_to_right = seq_along(col_order), item = col_order),
      paste0(stem, "_clustered_column_order.csv")
    )
  }

  list(
    data = out_df,
    row_order = row_order,
    col_order = col_order,
    cluster_matrix = mat_cluster
  )
}

safe_cor_test <- function(x, y, min_n = MIN_COR_PAIRS) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < min_n || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(data.frame(n = n, rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(n = n, rho = unname(ct$estimate), p = ct$p.value)
}

safe_wilcox <- function(value, group, min_group_n = MIN_GROUP_N) {
  ok <- is.finite(value) & !is.na(group)
  value <- value[ok]
  group <- droplevels(factor(group[ok], levels = c("Low", "High")))
  low <- value[group == "Low"]
  high <- value[group == "High"]
  if (length(low) < min_group_n || length(high) < min_group_n) {
    return(data.frame(
      n_low = length(low), n_high = length(high),
      median_low = NA_real_, median_high = NA_real_, median_diff = NA_real_,
      rank_biserial = NA_real_, p = NA_real_
    ))
  }
  wt <- suppressWarnings(wilcox.test(high, low, exact = FALSE))
  
  U <- as.numeric(wt$statistic)
  rbc <- 2 * U / (length(high) * length(low)) - 1
  data.frame(
    n_low = length(low), n_high = length(high),
    median_low = median(low, na.rm = TRUE),
    median_high = median(high, na.rm = TRUE),
    median_diff = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
    rank_biserial = rbc,
    p = wt$p.value
  )
}

read_gmt <- function(file) {
  if (!file.exists(file)) stop("GMT file not found: ", file)
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  out <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) next
    nm <- parts[1]
    genes <- unique(toupper(trimws(parts[-c(1, 2)])))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    out[[nm]] <- genes
  }
  if (length(out) == 0) stop("No gene sets parsed from: ", file)
  out
}

resolve_hallmark_gmt <- function() {
  if (length(HALLMARK_GMT) == 1 && !is.na(HALLMARK_GMT) && file.exists(HALLMARK_GMT)) {
    return(HALLMARK_GMT)
  }
  candidates <- list.files(
    file.path(BASE_DIR, "gene set"),
    pattern = "^h\\.all.*Hs\\.symbols\\.gmt$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(candidates) == 0) {
    stop(
      "Cannot find Hallmark GMT. Download the HUMAN Gene Symbols file named like ",
      "h.all.vXX.X.Hs.symbols.gmt and place it in: ", file.path(BASE_DIR, "gene set")
    )
  }
  candidates[order(file.info(candidates)$mtime, decreasing = TRUE)][1]
}

resolve_set_name <- function(gene_sets, preferred, fallback_patterns, label) {
  if (preferred %in% names(gene_sets)) return(preferred)
  hits <- unique(unlist(lapply(fallback_patterns, function(pt) {
    grep(pt, names(gene_sets), value = TRUE, ignore.case = TRUE)
  })))
  if (length(hits) == 1) {
    warning(label, " exact name not found; using: ", hits)
    return(hits)
  }
  stop(
    "Cannot uniquely identify ", label, " gene set. Preferred name: ", preferred,
    "\nPossible matches: ", paste(hits, collapse = ", ")
  )
}

extract_symbol_tokens <- function(x, symbol_reference) {
  ref <- unique(toupper(symbol_reference))
  x <- toupper(trimws(as.character(x)))
  vapply(x, function(v) {
    if (is.na(v) || !nzchar(v) || v == "NA") return(NA_character_)
    if (v %in% ref) return(v)
    toks <- unlist(strsplit(v, "\\s*(?:\\|+|///|;|,)\\s*", perl = TRUE))
    toks <- trimws(toks)
    toks <- toks[nzchar(toks)]
    if (!length(toks)) return(v)
    
    toks_match <- sub("^(ENSG[0-9]+)\\.[0-9]+$", "\\1", toks, perl = TRUE)
    hit <- which(toks %in% ref | toks_match %in% ref)
    if (length(hit)) {
      z <- toks[hit[1]]
      if (sub("^(ENSG[0-9]+)\\.[0-9]+$", "\\1", z, perl = TRUE) %in% ref) {
        z <- sub("^(ENSG[0-9]+)\\.[0-9]+$", "\\1", z, perl = TRUE)
      }
      return(z)
    }
    
    sym_like <- toks[grepl("^[A-Z][A-Z0-9._-]*$", toks) &
                       !grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", toks) &
                       !grepl("^[0-9]+$", toks)]
    if (length(sym_like)) return(sym_like[1])
    toks[1]
  }, character(1))
}

map_identifier_vector <- function(x, keytype) {
  keys <- toupper(trimws(as.character(x)))
  if (keytype == "ENSEMBL") keys <- sub("\\.[0-9]+$", "", keys)
  ans <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(keys[!is.na(keys) & nzchar(keys)]),
    column = "SYMBOL", keytype = keytype, multiVals = "first"
  ))
  unname(ans[keys])
}

read_gct_expression <- function(file, symbol_reference) {
  first <- readLines(file, n = 3, warn = FALSE)
  skip_n <- if (length(first) >= 1 && grepl("^#1\\.", first[1])) 2L else 0L
  dt <- data.table::fread(file, skip = skip_n, data.table = TRUE, check.names = FALSE)
  if (ncol(dt) < 3) stop("Expression file has too few columns: ", file)

  candidate_cols <- intersect(
    c("Name", "Description", "GeneSymbol", "gene_symbol", "Gene", "gene", "SYMBOL", "symbol"),
    colnames(dt)
  )
  if (length(candidate_cols) == 0) candidate_cols <- colnames(dt)[1:min(2, ncol(dt))]

  parsed_candidates <- lapply(candidate_cols, function(cc) {
    extract_symbol_tokens(dt[[cc]], symbol_reference)
  })
  names(parsed_candidates) <- candidate_cols
  overlap_n <- vapply(parsed_candidates, function(vals) {
    length(intersect(unique(vals[!is.na(vals)]), symbol_reference))
  }, numeric(1))

  gene_col <- candidate_cols[which.max(overlap_n)]
  gene <- parsed_candidates[[gene_col]]
  mapping_mode <- "parsed_symbol"

  if (max(overlap_n, na.rm = TRUE) == 0) {
    best_gene <- gene
    best_overlap <- 0L
    best_col <- gene_col
    best_mode <- mapping_mode

    for (cc in candidate_cols) {
      raw <- toupper(trimws(as.character(dt[[cc]])))
      raw_non_na <- raw[!is.na(raw) & nzchar(raw)]
      if (!length(raw_non_na)) next

      numeric_frac <- mean(grepl("^[0-9]+$", raw_non_na))
      ensembl_frac <- mean(grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", raw_non_na))

      mapped_list <- list()
      if (is.finite(numeric_frac) && numeric_frac >= 0.5) {
        mapped_list$ENTREZID <- map_identifier_vector(raw, "ENTREZID")
      }
      if (is.finite(ensembl_frac) && ensembl_frac >= 0.5) {
        mapped_list$ENSEMBL <- map_identifier_vector(raw, "ENSEMBL")
      }

      for (kt in names(mapped_list)) {
        mapped <- toupper(mapped_list[[kt]])
        ov <- length(intersect(unique(mapped[!is.na(mapped)]), symbol_reference))
        if (ov > best_overlap) {
          best_overlap <- ov
          best_gene <- mapped
          best_col <- cc
          best_mode <- paste0("mapped_", kt)
        }
      }
    }

    gene <- best_gene
    gene_col <- best_col
    mapping_mode <- best_mode
    overlap_n <- c(overlap_n, FALLBACK_MAPPING = best_overlap)
  }

  tcga_cols <- grep("^TCGA[-._]", colnames(dt), value = TRUE, ignore.case = TRUE)
  if (length(tcga_cols) > 0) {
    sample_cols <- tcga_cols
  } else {
    non_gene <- setdiff(colnames(dt), candidate_cols)
    numeric_fraction <- vapply(non_gene, function(cc) {
      v <- suppressWarnings(as.numeric(as.character(dt[[cc]])))
      mean(is.finite(v))
    }, numeric(1))
    sample_cols <- non_gene[numeric_fraction >= 0.8]
  }
  if (length(sample_cols) < MIN_SAMPLES_PER_CANCER) {
    stop("Too few expression sample columns in: ", file)
  }

  gene <- toupper(trimws(as.character(gene)))
  keep_gene <- !is.na(gene) & nzchar(gene) & gene != "NA"
  dt2 <- dt[keep_gene, sample_cols, with = FALSE]
  dt2[, GeneSymbol := gene[keep_gene]]
  data.table::setcolorder(dt2, c("GeneSymbol", sample_cols))
  for (cc in sample_cols) data.table::set(dt2, j = cc, value = suppressWarnings(as.numeric(dt2[[cc]])))

  dt2 <- dt2[, lapply(.SD, function(v) mean(v, na.rm = TRUE)),
             by = GeneSymbol, .SDcols = sample_cols]
  mat <- as.matrix(dt2[, -1, with = FALSE])
  rownames(mat) <- dt2$GeneSymbol
  storage.mode(mat) <- "double"

  if (any(!is.finite(mat))) {
    meds <- matrixStats::rowMedians(mat, na.rm = TRUE)
    meds[!is.finite(meds)] <- 0
    bad <- which(!is.finite(mat), arr.ind = TRUE)
    mat[bad] <- meds[bad[, 1]]
  }

  norm_ids <- normalize_tcga_id(colnames(mat))
  type_code <- sample_type_code(norm_ids)
  primary <- type_code == "01"
  if (sum(primary) >= MIN_SAMPLES_PER_CANCER) {
    mat <- mat[, primary, drop = FALSE]
    norm_ids <- norm_ids[primary]
  } else {
    warning("No reliable primary-tumor filter for ", basename(file), "; retaining all samples.")
  }

  s15 <- sample15(norm_ids)
  groups <- split(seq_along(s15), s15)
  mat2 <- do.call(cbind, lapply(groups, function(idx) {
    if (length(idx) == 1) mat[, idx] else rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
  }))
  if (is.null(dim(mat2))) mat2 <- matrix(mat2, ncol = 1)
  rownames(mat2) <- rownames(mat)
  colnames(mat2) <- names(groups)

  finite_vals <- mat2[is.finite(mat2)]
  q99 <- if (length(finite_vals)) unname(quantile(finite_vals, 0.99, na.rm = TRUE)) else NA_real_
  min_v <- if (length(finite_vals)) min(finite_vals, na.rm = TRUE) else NA_real_
  transformed <- FALSE
  if (is.finite(q99) && is.finite(min_v) && min_v >= 0 && q99 > 50) {
    mat2 <- log2(mat2 + 1)
    transformed <- TRUE
  }

  sds <- matrixStats::rowSds(mat2, na.rm = TRUE)
  mat2 <- mat2[is.finite(sds) & sds > 0, , drop = FALSE]

  post_overlap <- length(intersect(rownames(mat2), symbol_reference))

  list(
    expr = mat2,
    gene_col = gene_col,
    mapping_mode = mapping_mode,
    overlap = post_overlap,
    candidate_overlaps = paste(names(overlap_n), overlap_n, sep = "=", collapse = ";"),
    gct_skip = skip_n,
    log2_transformed = transformed,
    q99_before_transform = q99
  )
}

make_gene_set_qc <- function(expr, gene_sets, cancer, collection, method) {
  bind_rows(lapply(names(gene_sets), function(nm) {
    matched <- intersect(gene_sets[[nm]], rownames(expr))
    data.frame(
      cancer = cancer,
      collection = collection,
      method = method,
      gene_set = nm,
      n_input = length(gene_sets[[nm]]),
      n_matched = length(matched),
      matched_genes = paste(matched, collapse = ";"),
      missing_genes = paste(setdiff(gene_sets[[nm]], rownames(expr)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
}

zmean_score <- function(expr, gene_sets, min_matched = 1L) {
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  out <- matrix(NA_real_, nrow = length(gene_sets), ncol = ncol(expr),
                dimnames = list(names(gene_sets), colnames(expr)))
  for (nm in names(gene_sets)) {
    matched <- intersect(gene_sets[[nm]], rownames(z))
    if (length(matched) >= min_matched) {
      out[nm, ] <- colMeans(z[matched, , drop = FALSE], na.rm = TRUE)
    }
  }
  out
}

ssgsea_score <- function(expr, gene_sets, min_matched = 1L) {
  filtered <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr)))
  filtered <- filtered[lengths(filtered) >= min_matched]
  if (length(filtered) == 0) stop("No gene set passed ssGSEA matched-gene threshold.")

  exports <- getNamespaceExports("GSVA")
  ans <- NULL
  if ("ssgseaParam" %in% exports) {
    param <- tryCatch(
      GSVA::ssgseaParam(
        exprData = expr,
        geneSets = filtered,
        minSize = min_matched,
        maxSize = Inf,
        alpha = 0.25,
        normalize = TRUE
      ),
      error = function(e) {
        msg("Could not construct ssgseaParam: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(param)) {
      ans <- tryCatch(
        GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()),
        error = function(e1) {
          msg("GSVA ssGSEA with SerialParam failed; retrying without BPPARAM. Error: ", conditionMessage(e1))
          tryCatch(
            GSVA::gsva(param, verbose = FALSE),
            error = function(e2) {
              msg("New GSVA ssGSEA API failed; trying legacy API. Error: ", conditionMessage(e2))
              NULL
            }
          )
        }
      )
    }
  }

  if (is.null(ans)) {
    ans <- GSVA::gsva(
      expr = expr,
      gset.idx.list = filtered,
      method = "ssgsea",
      min.sz = min_matched,
      max.sz = Inf,
      parallel.sz = 1,
      verbose = FALSE,
      ssgsea.norm = TRUE
    )
  }
  ans <- as.matrix(ans)

  out <- matrix(NA_real_, nrow = length(gene_sets), ncol = ncol(expr),
                dimnames = list(names(gene_sets), colnames(expr)))
  common <- intersect(rownames(out), rownames(ans))
  out[common, ] <- ans[common, colnames(expr), drop = FALSE]
  out
}

residualize_caf3 <- function(step, caf3, cancer, method) {
  ok <- is.finite(step) & is.finite(caf3)
  result <- data.frame(
    residual = rep(NA_real_, length(step)),
    residual_z = rep(NA_real_, length(step))
  )
  if (sum(ok) < MIN_COR_PAIRS || length(unique(caf3[ok])) < 3 || length(unique(step[ok])) < 3) {
    params <- data.frame(
      cancer = cancer, method = method, n = sum(ok), alpha = NA_real_, beta_CAF3 = NA_real_,
      beta_p = NA_real_, r_squared = NA_real_, stringsAsFactors = FALSE
    )
    return(list(values = result, params = params))
  }
  fit <- lm(step[ok] ~ caf3[ok])
  res <- residuals(fit)
  result$residual[ok] <- res
  result$residual_z[ok] <- safe_z(res)
  sm <- summary(fit)
  params <- data.frame(
    cancer = cancer,
    method = method,
    n = sum(ok),
    alpha = unname(coef(fit)[1]),
    beta_CAF3 = unname(coef(fit)[2]),
    beta_p = unname(sm$coefficients[2, 4]),
    r_squared = unname(sm$r.squared),
    stringsAsFactors = FALSE
  )
  list(values = result, params = params)
}

random_effects_meta_cor <- function(rho, n) {
  ok <- is.finite(rho) & is.finite(n) & n > 3 & abs(rho) < 1
  rho <- rho[ok]
  n <- n[ok]
  k <- length(rho)
  if (k == 0) {
    return(data.frame(k = 0, meta_rho = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p = NA_real_, tau2 = NA_real_, I2 = NA_real_))
  }
  yi <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  vi <- 1 / (n - 3)
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - fixed)^2)
  C <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (k > 1 && C > 0) max(0, (Q - (k - 1)) / C) else 0
  wr <- 1 / (vi + tau2)
  est <- sum(wr * yi) / sum(wr)
  se <- sqrt(1 / sum(wr))
  zstat <- est / se
  p <- 2 * pnorm(-abs(zstat))
  ci <- est + c(-1, 1) * qnorm(0.975) * se
  I2 <- if (k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  data.frame(
    k = k,
    meta_rho = tanh(est),
    ci_low = tanh(ci[1]),
    ci_high = tanh(ci[2]),
    p = p,
    tau2 = tau2,
    I2 = I2
  )
}

vectorized_gene_spearman <- function(expr, score) {
  ok <- is.finite(score)
  expr <- expr[, ok, drop = FALSE]
  score <- score[ok]
  n <- length(score)
  if (n < MIN_COR_PAIRS || length(unique(score)) < 3) {
    return(data.frame(gene = rownames(expr), n = n, rho = NA_real_, p = NA_real_))
  }
  yr <- rank(score, ties.method = "average")
  yr <- yr - mean(yr)
  xr <- matrixStats::rowRanks(expr, ties.method = "average")
  xr <- xr - rowMeans(xr)
  denom <- sqrt(rowSums(xr^2) * sum(yr^2))
  rho <- as.numeric(xr %*% yr) / denom
  rho[!is.finite(rho)] <- NA_real_
  rho <- pmin(pmax(rho, -0.999999), 0.999999)
  tstat <- rho * sqrt((n - 2) / pmax(1 - rho^2, .Machine$double.eps))
  p <- 2 * pt(-abs(tstat), df = n - 2)
  data.frame(gene = rownames(expr), n = n, rho = rho, p = p)
}

run_fgsea_one <- function(stats_vec, pathways, collection, cancer, method) {
  stats_vec <- stats_vec[is.finite(stats_vec)]
  stats_vec <- stats_vec[!duplicated(names(stats_vec))]
  stats_vec <- stats_vec[order(-stats_vec, names(stats_vec))]
  if (length(stats_vec) < 100) return(NULL)

  pathways2 <- lapply(pathways, function(gs) intersect(unique(gs), names(stats_vec)))
  pathways2 <- pathways2[
    lengths(pathways2) >= GSEA_MIN_SIZE &
      lengths(pathways2) <= GSEA_MAX_SIZE
  ]
  if (length(pathways2) == 0) return(NULL)

  fg_raw <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = pathways2,
    stats = stats_vec,
    minSize = GSEA_MIN_SIZE,
    maxSize = GSEA_MAX_SIZE,
    eps = 0,
    scoreType = "std"
  ))
  if (nrow(fg_raw) == 0) return(NULL)

  main_pathways <- as.character(fg_raw$pathway)
  parent_map <- setNames(rep(NA_character_, length(main_pathways)), main_pathways)
  if (identical(collection, "REACTOME") && identical(cancer, "PAN_CANCER_META") &&
      isTRUE(REACTOME_COLLAPSE_FOR_PLOTS)) {
    cp <- tryCatch(
      fgsea::collapsePathways(
        fgseaRes = fg_raw,
        pathways = pathways2,
        stats = stats_vec
      ),
      error = function(e) NULL
    )
    if (!is.null(cp)) {
      main_pathways <- as.character(cp$mainPathways)
      if (!is.null(cp$parentPathways)) {
        parent_map <- as.character(cp$parentPathways)
        names(parent_map) <- names(cp$parentPathways)
      }
    }
  }

  fg <- as.data.frame(fg_raw)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg$is_main_pathway <- fg$pathway %in% main_pathways
  fg$parent_pathway <- unname(parent_map[fg$pathway])
  fg$collection <- collection
  fg$cancer <- cancer
  fg$method <- method
  fg
}

format_pvalue_short <- function(p, label = "P") {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return(paste0(label, " = NA"))
  if (p < 1e-4) paste0(label, " = ", formatC(p, format = "e", digits = 2))
  else paste0(label, " = ", formatC(p, format = "fg", digits = 3))
}

split_gene_string <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  unique(trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE))))
}

safe_file_token <- function(x, max_n = 120L) {
  y <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  substr(y, 1, max_n)
}

calculate_running_es <- function(stats_vec, gene_set, exponent = 1) {
  stats_vec <- stats_vec[is.finite(stats_vec)]
  stats_vec <- stats_vec[!duplicated(names(stats_vec))]
  stats_vec <- stats_vec[order(-stats_vec, names(stats_vec))]

  genes <- names(stats_vec)
  hits <- genes %in% intersect(gene_set, genes)
  n_total <- length(stats_vec)
  n_hit <- sum(hits)
  n_miss <- n_total - n_hit
  if (n_hit == 0 || n_miss == 0) return(NULL)

  weights <- abs(stats_vec)^exponent
  hit_denom <- sum(weights[hits])
  if (!is.finite(hit_denom) || hit_denom <= 0) {
    weights <- rep(1, length(stats_vec))
    hit_denom <- sum(weights[hits])
  }

  increments <- numeric(n_total)
  increments[hits] <- weights[hits] / hit_denom
  increments[!hits] <- -1 / n_miss
  running_es <- cumsum(increments)
  max_i <- which.max(running_es)
  min_i <- which.min(running_es)
  peak_i <- if (abs(running_es[max_i]) >= abs(running_es[min_i])) max_i else min_i

  data.frame(
    rank = seq_len(n_total),
    gene = genes,
    rank_metric = as.numeric(stats_vec),
    hit = hits,
    running_es = running_es,
    peak = seq_len(n_total) == peak_i,
    stringsAsFactors = FALSE
  )
}

make_gsea_rank_gradient <- function(prof) {
  dd <- prof %>%
    dplyr::transmute(
      rank = rank,
      rank_metric = as.numeric(rank_metric)
    )

  lim <- max(abs(dd$rank_metric), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1

  dd$rank_metric_clipped <- pmax(pmin(dd$rank_metric, lim), -lim)
  dd
}

plot_gsea_hit_strip <- function(prof) {
  ggplot(prof %>% dplyr::filter(hit), aes(rank)) +
    geom_segment(
      aes(xend = rank, y = 0, yend = 1),
      linewidth = 0.28, color = "black", alpha = 0.78
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
}

plot_gsea_rank_gradient_strip <- function(prof, rank_label = NULL, show_x_axis = TRUE) {
  gd <- make_gsea_rank_gradient(prof)
  lim <- max(abs(gd$rank_metric_clipped), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1

  p <- ggplot(gd, aes(rank, y = 1, fill = rank_metric_clipped)) +
    geom_tile(height = 1, width = 1.05) +
    scale_fill_gradient2(
      low = NEG_COLOR,
      mid = "#F4B6C2",
      high = POS_COLOR,
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      guide = "none"
    ) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 5.5, if (show_x_axis) 5.5 else 0, 5.5))

  if (isTRUE(show_x_axis)) {
    p <- p +
      theme(
        axis.text.x = element_text(size = 7, color = "black"),
        axis.ticks.x = element_line(linewidth = 0.3, color = "black"),
        axis.title.x = element_text(size = 7.2, color = "black", margin = margin(t = 3)),
        axis.line.x = element_line(linewidth = 0.3, color = "black")
      ) +
      labs(x = rank_label)
  }

  p
}

curve_stat_label <- function(gsea_row) {
  leading_n <- length(split_gene_string(gsea_row$leadingEdge[1]))
  paste0(
    "NES = ", sprintf("%.3f", as.numeric(gsea_row$NES[1])), "\n",
    format_pvalue_short(as.numeric(gsea_row$pval[1]), "Nominal P"), "\n",
    format_pvalue_short(as.numeric(gsea_row$padj[1]), "BH FDR"), "\n",
    "Set size = ", gsea_row$size[1], "\n",
    "Leading edge = ", leading_n
  )
}

plot_compact_pan_enrichment <- function(stats_vec, gene_set, gsea_row, rank_label) {
  prof <- calculate_running_es(stats_vec, gene_set)
  if (is.null(prof)) return(NULL)

  nes <- as.numeric(gsea_row$NES[1])
  curve_color <- ifelse(nes >= 0, POS_COLOR, NEG_COLOR)
  peak_rank <- prof$rank[which(prof$peak)[1]]

  p1 <- ggplot(prof, aes(rank, running_es)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, color = "grey55") +
    geom_line(linewidth = 0.85, color = curve_color) +
    geom_vline(xintercept = peak_rank, linetype = 3, linewidth = 0.35, color = curve_color) +
    annotate(
      "label", x = Inf, y = Inf, label = curve_stat_label(gsea_row),
      hjust = 1.03, vjust = 1.05, size = 2.35,
      label.size = 0.2, fill = scales::alpha("white", 0.90)
    ) +
    theme_pub(9) +
    theme(
      plot.title = element_text(size = 9.0, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7.0, hjust = 0.5),
      axis.text.y = element_text(size = 6.8),
      axis.title.y = element_text(size = 7.2),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    ) +
    labs(
      title = clean_pathway_name(gsea_row$pathway[1]),
      subtitle = paste0(
        "CAF3-adjusted STEP ", gsea_row$method[1],
        " | ", gsea_row$analysis_variant[1]
      ),
      y = "Running enrichment score"
    )

  p2 <- plot_gsea_hit_strip(prof)

  p3 <- plot_gsea_rank_gradient_strip(prof, rank_label = rank_label, show_x_axis = TRUE)

  p1 / p2 / p3 +
    patchwork::plot_layout(heights = c(3.15, 0.48, 0.42))
}

plot_classic_pan_enrichment <- function(stats_vec, gene_set, gsea_row, rank_label) {
  prof <- calculate_running_es(stats_vec, gene_set)
  if (is.null(prof)) return(NULL)

  nes <- as.numeric(gsea_row$NES[1])
  curve_color <- ifelse(nes >= 0, POS_COLOR, NEG_COLOR)
  peak_rank <- prof$rank[which(prof$peak)[1]]

  p1 <- ggplot(prof, aes(rank, running_es)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, color = "grey55") +
    geom_line(linewidth = 0.95, color = curve_color) +
    geom_vline(xintercept = peak_rank, linetype = 3, linewidth = 0.4, color = curve_color) +
    annotate(
      "label", x = Inf, y = Inf, label = curve_stat_label(gsea_row),
      hjust = 1.03, vjust = 1.05, size = 2.75,
      label.size = 0.22, fill = scales::alpha("white", 0.92)
    ) +
    theme_pub(10) +
    theme(
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), plot.margin = margin(5.5, 5.5, 0, 5.5)
    ) +
    labs(
      title = clean_pathway_name(gsea_row$pathway[1]),
      subtitle = paste0(
        "CAF3-adjusted STEP ", gsea_row$method[1], " | ",
        gsea_row$collection[1], " | ", gsea_row$analysis_variant[1]
      ),
      y = "Running ES"
    )

  p2 <- plot_gsea_hit_strip(prof)

  p3 <- plot_gsea_rank_gradient_strip(prof, rank_label = rank_label, show_x_axis = TRUE)

  p1 / p2 / p3 +
    patchwork::plot_layout(heights = c(3.2, 0.55, 0.50))
}

select_top_gsea_curves <- function(fg) {
  sig <- fg %>%
    filter(
      is.finite(padj), padj < GSEA_CURVE_FDR_CUTOFF,
      is.finite(NES), abs(NES) >= GSEA_CURVE_MIN_ABS_NES
    )
  if ("is_main_pathway" %in% colnames(sig) &&
      any(sig$collection == "REACTOME") && isTRUE(REACTOME_COLLAPSE_FOR_PLOTS)) {
    sig <- sig %>% filter(collection != "REACTOME" | is_main_pathway)
  }

  pos <- sig %>% filter(NES > 0) %>%
    arrange(padj, desc(abs(NES)), desc(NES)) %>%
    slice_head(n = TOP_CURVES_POSITIVE) %>% mutate(direction = "Positive")
  neg <- sig %>% filter(NES < 0) %>%
    arrange(padj, desc(abs(NES)), NES) %>%
    slice_head(n = TOP_CURVES_NEGATIVE) %>% mutate(direction = "Negative")
  bind_rows(pos, neg)
}

save_pan_curve_outputs <- function(fg, stats_vec, pathways, step_method,
                                   collection, analysis_variant, rank_label) {
  selected <- select_top_gsea_curves(fg)
  stem_prefix <- paste0(
    "PAN_CANCER_META_", step_method, "_", collection, "_", analysis_variant
  )
  fwrite(
    selected,
    file.path(OUT_DIR, "11_GSEA", paste0(stem_prefix, "_pathways_selected_for_curves.csv"))
  )
  if (nrow(selected) == 0) {
    msg("No FDR-significant pathways available for curves: ", stem_prefix)
    return(invisible(NULL))
  }

  compact_plots <- lapply(seq_len(nrow(selected)), function(i) {
    rr <- selected[i, , drop = FALSE]
    gs <- intersect(pathways[[rr$pathway[1]]], names(stats_vec))
    plot_compact_pan_enrichment(stats_vec, gs, rr, rank_label)
  })
  keep <- !vapply(compact_plots, is.null, logical(1))
  compact_plots <- compact_plots[keep]
  selected2 <- selected[keep, , drop = FALSE]

  pos_plots <- compact_plots[selected2$direction == "Positive"]
  neg_plots <- compact_plots[selected2$direction == "Negative"]

  if (length(pos_plots)) {
    ppos <- patchwork::wrap_plots(pos_plots, ncol = 2) +
      patchwork::plot_annotation(
        title = paste0(stem_prefix, ": Top significant positive pathways"),
        subtitle = paste0("Selection: BH FDR < ", GSEA_CURVE_FDR_CUTOFF,
                          "; up to ", TOP_CURVES_POSITIVE, " pathways")
      )
    save_plot(ppos, file.path(OUT_DIR, "11_GSEA", paste0(stem_prefix, "_top6_positive_enrichment_curves")),
              width = 13, height = max(7, ceiling(length(pos_plots) / 2) * 4.2))
  }
  if (length(neg_plots)) {
    pneg <- patchwork::wrap_plots(neg_plots, ncol = 2) +
      patchwork::plot_annotation(
        title = paste0(stem_prefix, ": Top significant negative pathways"),
        subtitle = paste0("Selection: BH FDR < ", GSEA_CURVE_FDR_CUTOFF,
                          "; up to ", TOP_CURVES_NEGATIVE, " pathways")
      )
    save_plot(pneg, file.path(OUT_DIR, "11_GSEA", paste0(stem_prefix, "_top6_negative_enrichment_curves")),
              width = 13, height = max(7, ceiling(length(neg_plots) / 2) * 4.2))
  }
  if (length(compact_plots)) {
    pall <- patchwork::wrap_plots(compact_plots, ncol = 3) +
      patchwork::plot_annotation(
        title = paste0(stem_prefix, ": Top 6 positive + Top 6 negative"),
        subtitle = paste0(
          "All plotted pathways satisfy BH FDR < ", GSEA_CURVE_FDR_CUTOFF,
          "; positive=red, negative=blue"
        )
      )
    save_plot(pall, file.path(OUT_DIR, "11_GSEA", paste0(stem_prefix, "_top6_positive_plus_top6_negative_enrichment_curves")),
              width = 17, height = max(9, ceiling(length(compact_plots) / 3) * 4.0))
  }

  if (isTRUE(SAVE_INDIVIDUAL_CLASSIC_CURVES)) {
    indiv_dir <- file.path(OUT_DIR, "11_GSEA", "individual_classic_curves", stem_prefix)
    dir.create(indiv_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_len(nrow(selected2))) {
      rr <- selected2[i, , drop = FALSE]
      gs <- intersect(pathways[[rr$pathway[1]]], names(stats_vec))
      pc <- plot_classic_pan_enrichment(stats_vec, gs, rr, rank_label)
      if (is.null(pc)) next
      direction_rank <- sum(selected2$direction[seq_len(i)] == rr$direction[1])
      stem <- file.path(
        indiv_dir,
        paste0(sprintf("%02d", direction_rank), "_", rr$direction[1], "_", safe_file_token(rr$pathway[1]))
      )
      save_plot(pc, stem, width = 10, height = 8)
    }
  }
  invisible(selected)
}

compute_gene_meta_table <- function(gene_dt_input) {
  dt <- data.table::as.data.table(data.table::copy(gene_dt_input))
  dt <- dt[is.finite(rho) & is.finite(n) & n > 3]
  if (nrow(dt) == 0) return(data.frame())

  totals <- dt[, .(total_tested_cancers = uniqueN(cancer)), by = method]
  ans <- dt[, {
    rho_use <- pmin(pmax(rho, -0.999999), 0.999999)
    yi <- atanh(rho_use)
    vi <- 1 / (n - 3)
    wi <- 1 / vi

    fixed_est <- sum(wi * yi) / sum(wi)
    fixed_se <- sqrt(1 / sum(wi))
    fixed_z <- fixed_est / fixed_se

    Q <- sum(wi * (yi - fixed_est)^2)
    C <- sum(wi) - sum(wi^2) / sum(wi)
    tau2 <- if (.N > 1 && is.finite(C) && C > 0) max(0, (Q - (.N - 1)) / C) else 0
    wr <- 1 / (vi + tau2)
    random_est <- sum(wr * yi) / sum(wr)
    random_se <- sqrt(1 / sum(wr))
    random_z <- random_est / random_se
    I2 <- if (.N > 1 && Q > 0) max(0, (Q - (.N - 1)) / Q) * 100 else 0

    list(
      k_cancers = .N,
      total_n = sum(n, na.rm = TRUE),
      fixed_meta_rho = tanh(fixed_est),
      fixed_meta_zstat = fixed_z,
      fixed_p = 2 * pnorm(-abs(fixed_z)),
      random_meta_rho = tanh(random_est),
      random_meta_zstat = random_z,
      random_p = 2 * pnorm(-abs(random_z)),
      Q = Q,
      Q_p = if (.N > 1) pchisq(Q, df = .N - 1, lower.tail = FALSE) else NA_real_,
      tau2 = tau2,
      I2 = I2,
      positive_sig_cancers = sum(FDR < 0.05 & rho > 0, na.rm = TRUE),
      negative_sig_cancers = sum(FDR < 0.05 & rho < 0, na.rm = TRUE)
    )
  }, by = .(method, gene)]

  ans <- merge(ans, totals, by = "method", all.x = TRUE)
  ans[, coverage_fraction := k_cancers / total_tested_cancers]
  ans[, min_required_cancers := pmax(
    MIN_CANCERS_FOR_GENE_META,
    ceiling(MIN_CANCER_COVERAGE_FRACTION * total_tested_cancers)
  )]
  ans[, eligible_for_pan_gsea := k_cancers >= min_required_cancers]
  ans[, fixed_FDR := p.adjust(fixed_p, method = "BH"), by = method]
  ans[, random_FDR := p.adjust(random_p, method = "BH"), by = method]

  if (identical(PAN_GSEA_PRIMARY_META_MODEL, "random")) {
    ans[, `:=`(
      meta_rho = random_meta_rho,
      meta_zstat = random_meta_zstat,
      p = random_p,
      FDR = random_FDR
    )]
  } else {
    ans[, `:=`(
      meta_rho = fixed_meta_rho,
      meta_zstat = fixed_meta_zstat,
      p = fixed_p,
      FDR = fixed_FDR
    )]
  }
  as.data.frame(ans)
}

build_pan_variant_specs <- function() {
  primary_variant <- paste0(
    "PRIMARY_", toupper(PAN_GSEA_PRIMARY_META_MODEL), "_",
    toupper(gsub("meta_", "", PAN_GSEA_PRIMARY_RANK_METRIC))
  )
  out <- data.frame(
    analysis_variant = primary_variant,
    meta_model = PAN_GSEA_PRIMARY_META_MODEL,
    rank_metric = PAN_GSEA_PRIMARY_RANK_METRIC,
    deoverlap = FALSE,
    plot_curves = TRUE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(RUN_FIXED_EFFECT_META_Z_SENSITIVITY) &&
      !(PAN_GSEA_PRIMARY_META_MODEL == "fixed" && PAN_GSEA_PRIMARY_RANK_METRIC == "meta_zstat")) {
    out <- bind_rows(out, data.frame(
      analysis_variant = "SENSITIVITY_FIXED_META_Z",
      meta_model = "fixed", rank_metric = "meta_zstat",
      deoverlap = FALSE, plot_curves = PLOT_SENSITIVITY_CURVES
    ))
  }
  if (isTRUE(RUN_META_RHO_RANK_SENSITIVITY) &&
      PAN_GSEA_PRIMARY_RANK_METRIC != "meta_rho") {
    out <- bind_rows(out, data.frame(
      analysis_variant = paste0("SENSITIVITY_", toupper(PAN_GSEA_PRIMARY_META_MODEL), "_META_RHO"),
      meta_model = PAN_GSEA_PRIMARY_META_MODEL, rank_metric = "meta_rho",
      deoverlap = FALSE, plot_curves = PLOT_SENSITIVITY_CURVES
    ))
  }
  if (isTRUE(RUN_DEOVERLAP_SENSITIVITY)) {
    out <- bind_rows(out, data.frame(
      analysis_variant = paste0(primary_variant, "_DEOVERLAP"),
      meta_model = PAN_GSEA_PRIMARY_META_MODEL,
      rank_metric = PAN_GSEA_PRIMARY_RANK_METRIC,
      deoverlap = TRUE, plot_curves = PLOT_SENSITIVITY_CURVES
    ))
  }
  out %>% distinct(analysis_variant, .keep_all = TRUE)
}

run_ora_one <- function(selected, universe, pathways, collection, direction, method) {
  universe <- unique(universe)
  selected <- intersect(unique(selected), universe)
  if (length(selected) == 0) return(NULL)
  out <- bind_rows(lapply(names(pathways), function(nm) {
    gs <- intersect(unique(pathways[[nm]]), universe)
    if (length(gs) < GSEA_MIN_SIZE || length(gs) > GSEA_MAX_SIZE) return(NULL)
    a <- length(intersect(selected, gs))
    b <- length(setdiff(selected, gs))
    c <- length(setdiff(gs, selected))
    d <- length(setdiff(universe, union(selected, gs)))
    if (a == 0) return(NULL)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    data.frame(
      method = method,
      collection = collection,
      direction = direction,
      pathway = nm,
      overlap = a,
      selected_n = length(selected),
      pathway_n = length(gs),
      gene_ratio = a / length(selected),
      odds_ratio = unname(ft$estimate),
      p = ft$p.value,
      overlap_genes = paste(intersect(selected, gs), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(out)) out$FDR <- p.adjust(out$p, method = "BH")
  out
}

if (!dir.exists(EXP_DIR)) stop("Expression folder not found: ", EXP_DIR)
if (!file.exists(STEP_GMT)) stop("STEP GMT not found: ", STEP_GMT)
HALLMARK_GMT <- resolve_hallmark_gmt()

step_all_sets <- read_gmt(STEP_GMT)
hallmark_sets <- read_gmt(HALLMARK_GMT)
step_name <- resolve_set_name(
  step_all_sets, STEP_SET_NAME,
  c("STEP.*UP.*84", "STEP_UP_HYBRID_LEADING_EDGE", "^STEP_UP$"),
  "STEP_UP"
)
caf3_name <- resolve_set_name(
  step_all_sets, CAF3_SET_NAME,
  c("CAF.*3.*gene", "CAF_3gene", "CAF3"),
  "CAF3"
)
core_sets <- list(
  STEP_UP = step_all_sets[[step_name]],
  CAF3 = step_all_sets[[caf3_name]]
)

enrichment_collections <- list(HALLMARK = hallmark_sets)
for (nm in names(OPTIONAL_ENRICHMENT_GMTS)) {
  f <- OPTIONAL_ENRICHMENT_GMTS[[nm]]
  
  if (!file.exists(f)) {
    prefix <- if (nm == "REACTOME") "^c2\\.cp\\.reactome.*Hs\\.symbols\\.gmt$" else
      if (nm == "GOBP") "^c5\\.go\\.bp.*Hs\\.symbols\\.gmt$" else NULL
    if (!is.null(prefix)) {
      cand <- list.files(file.path(BASE_DIR, "gene set"), pattern = prefix,
                         full.names = TRUE, ignore.case = TRUE)
      if (length(cand)) f <- cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
    }
  }
  if (file.exists(f)) {
    enrichment_collections[[nm]] <- read_gmt(f)
    msg("Loaded optional enrichment collection ", nm, ": ", f)
  } else {
    msg("Optional enrichment collection skipped (file absent): ", nm)
  }
}

all_symbol_reference <- unique(c(
  unlist(core_sets, use.names = FALSE),
  unlist(hallmark_sets, use.names = FALSE),
  unlist(enrichment_collections, recursive = TRUE, use.names = FALSE)
))

expr_files <- list.files(EXP_DIR, pattern = "_tcga_expression\\.gct$",
                         full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (length(expr_files) == 0) stop("No *_tcga_expression.gct files found under: ", EXP_DIR)

input_summary <- data.frame(
  item = c("BASE_DIR", "EXP_DIR", "STEP_GMT", "HALLMARK_GMT", "STEP_SET", "CAF3_SET", "N_EXPRESSION_FILES"),
  value = c(BASE_DIR, EXP_DIR, STEP_GMT, HALLMARK_GMT, step_name, caf3_name, length(expr_files)),
  stringsAsFactors = FALSE
)
write.csv(input_summary, file.path(OUT_DIR, "00_logs", "input_summary.csv"), row.names = FALSE)

all_scores <- list()
all_qc <- list()
all_resid_params <- list()
all_hallmark_cor <- list()
all_hallmark_highlow <- list()
all_concordance <- list()
all_gene_cor <- list()
all_gsea <- list()

for (fi in seq_along(expr_files)) {
  f <- expr_files[fi]
  cancer <- toupper(sub("_tcga_expression\\.gct$", "", basename(f), ignore.case = TRUE))
  msg("[", fi, "/", length(expr_files), "] Processing ", cancer, ": ", basename(f))

  one <- tryCatch({
    read_info <- read_gct_expression(f, all_symbol_reference)
    expr <- read_info$expr
    if (ncol(expr) < MIN_SAMPLES_PER_CANCER) {
      stop("Only ", ncol(expr), " usable samples after filtering.")
    }

    qc_core_z <- make_gene_set_qc(expr, core_sets, cancer, "CORE", "Zmean")
    qc_core_s <- make_gene_set_qc(expr, core_sets, cancer, "CORE", "ssGSEA")
    qc_hall_s <- make_gene_set_qc(expr, hallmark_sets, cancer, "HALLMARK", "ssGSEA")
    all_qc[[paste0(cancer, "_core_z")]] <- qc_core_z
    all_qc[[paste0(cancer, "_core_s")]] <- qc_core_s
    all_qc[[paste0(cancer, "_hall_s")]] <- qc_hall_s

    combined_sets <- c(core_sets, hallmark_sets)
    nmatch <- vapply(combined_sets, function(gs) length(intersect(gs, rownames(expr))), numeric(1))
    msg(
      cancer, " matched genes: STEP=", nmatch["STEP_UP"],
      ", CAF3=", nmatch["CAF3"],
      ", Hallmark sets with >=1 gene=", sum(nmatch[names(hallmark_sets)] >= 1), "/", length(hallmark_sets),
      ", Hallmark sets with >=", MIN_MATCHED_HALLMARK_GENES, " genes=",
      sum(nmatch[names(hallmark_sets)] >= MIN_MATCHED_HALLMARK_GENES), "/", length(hallmark_sets)
    )
    if (sum(nmatch >= 1) == 0) {
      stop(
        "Zero genes matched between expression row names and all STEP/CAF3/Hallmark gene sets. ",
        "This is an identifier/annotation parsing problem, not a threshold problem. ",
        "Selected gene column=", read_info$gene_col,
        "; mapping_mode=", read_info$mapping_mode,
        "; candidate overlaps=", read_info$candidate_overlaps
      )
    }
    z_scores <- zmean_score(expr, combined_sets, min_matched = 1L)
    s_scores <- ssgsea_score(expr, combined_sets, min_matched = 1L)

    if (nmatch["STEP_UP"] < MIN_MATCHED_STEP_GENES) {
      z_scores["STEP_UP", ] <- NA_real_
      s_scores["STEP_UP", ] <- NA_real_
    }
    if (nmatch["CAF3"] < MIN_MATCHED_CAF3_GENES) {
      z_scores["CAF3", ] <- NA_real_
      s_scores["CAF3", ] <- NA_real_
    }
    low_hallmark <- names(hallmark_sets)[nmatch[names(hallmark_sets)] < MIN_MATCHED_HALLMARK_GENES]
    if (length(low_hallmark)) {
      z_scores[low_hallmark, ] <- NA_real_
      s_scores[low_hallmark, ] <- NA_real_
    }

    ids <- colnames(expr)
    score_df <- data.frame(
      sample15 = ids,
      patient12 = patient12(ids),
      cancer = cancer,
      STEP_Zmean_raw = as.numeric(z_scores["STEP_UP", ids]),
      CAF3_Zmean = as.numeric(z_scores["CAF3", ids]),
      STEP_ssGSEA_raw = as.numeric(s_scores["STEP_UP", ids]),
      CAF3_ssGSEA = as.numeric(s_scores["CAF3", ids]),
      stringsAsFactors = FALSE
    )

    rz <- residualize_caf3(score_df$STEP_Zmean_raw, score_df$CAF3_Zmean, cancer, "Zmean")
    rs <- residualize_caf3(score_df$STEP_ssGSEA_raw, score_df$CAF3_ssGSEA, cancer, "ssGSEA")
    score_df$STEP_Zmean_residual_CAF3 <- rz$values$residual
    score_df$STEP_Zmean_residual_CAF3_z <- rz$values$residual_z
    score_df$STEP_ssGSEA_residual_CAF3 <- rs$values$residual
    score_df$STEP_ssGSEA_residual_CAF3_z <- rs$values$residual_z
    all_resid_params[[paste0(cancer, "_Z")]] <- rz$params
    all_resid_params[[paste0(cancer, "_S")]] <- rs$params

    hall_s <- t(s_scores[names(hallmark_sets), ids, drop = FALSE])
    colnames(hall_s) <- paste0("HALLMARK_ssGSEA__", colnames(hall_s))
    score_df <- cbind(score_df, hall_s)
    if (RUN_HALLMARK_ZMEAN_SENSITIVITY) {
      hall_z <- t(z_scores[names(hallmark_sets), ids, drop = FALSE])
      colnames(hall_z) <- paste0("HALLMARK_Zmean__", colnames(hall_z))
      score_df <- cbind(score_df, hall_z)
    }

    fwrite(score_df, file.path(OUT_DIR, "02_scores_by_cancer", paste0(cancer, "_STEP_CAF3_Hallmark_scores.csv")))
    all_scores[[cancer]] <- score_df

    resid_detail <- score_df[, c(
      "sample15", "patient12", "cancer",
      "STEP_Zmean_raw", "CAF3_Zmean", "STEP_Zmean_residual_CAF3", "STEP_Zmean_residual_CAF3_z",
      "STEP_ssGSEA_raw", "CAF3_ssGSEA", "STEP_ssGSEA_residual_CAF3", "STEP_ssGSEA_residual_CAF3_z"
    )]
    fwrite(resid_detail, file.path(OUT_DIR, "03_CAF3_residualization", paste0(cancer, "_CAF3_residualization_sample_details.csv")))

    conc <- safe_cor_test(
      score_df$STEP_Zmean_residual_CAF3,
      score_df$STEP_ssGSEA_residual_CAF3,
      min_n = MIN_COR_PAIRS
    )
    conc$cancer <- cancer
    all_concordance[[cancer]] <- conc

    hallmark_method_list <- list(ssGSEA = s_scores[names(hallmark_sets), ids, drop = FALSE])
    if (RUN_HALLMARK_ZMEAN_SENSITIVITY) {
      hallmark_method_list$Zmean <- z_scores[names(hallmark_sets), ids, drop = FALSE]
    }

    step_specs <- list(
      Zmean_raw = score_df$STEP_Zmean_raw,
      Zmean_residual = score_df$STEP_Zmean_residual_CAF3,
      ssGSEA_raw = score_df$STEP_ssGSEA_raw,
      ssGSEA_residual = score_df$STEP_ssGSEA_residual_CAF3
    )

    cor_rows <- list()
    kk <- 1L
    for (hall_method in names(hallmark_method_list)) {
      hm <- hallmark_method_list[[hall_method]]
      for (spec_name in names(step_specs)) {
        parts <- strsplit(spec_name, "_", fixed = TRUE)[[1]]
        step_method <- parts[1]
        adjustment <- parts[2]
        x <- step_specs[[spec_name]]
        for (pathway in rownames(hm)) {
          ct <- safe_cor_test(x, as.numeric(hm[pathway, ids]), min_n = MIN_COR_PAIRS)
          cor_rows[[kk]] <- data.frame(
            cancer = cancer,
            STEP_method = step_method,
            adjustment = adjustment,
            Hallmark_method = hall_method,
            pathway = pathway,
            n = ct$n,
            rho = ct$rho,
            p = ct$p,
            stringsAsFactors = FALSE
          )
          kk <- kk + 1L
        }
      }
    }
    cor_df <- bind_rows(cor_rows) %>%
      group_by(cancer, STEP_method, adjustment, Hallmark_method) %>%
      mutate(FDR = p.adjust(p, method = "BH"), sig = p_star(FDR)) %>%
      ungroup()
    fwrite(cor_df, file.path(OUT_DIR, "04_Hallmark_correlations", paste0(cancer, "_STEP_Hallmark_correlations.csv")))
    all_hallmark_cor[[cancer]] <- cor_df

    hl_rows <- list()
    box_data_list <- list()
    kk <- 1L
    for (step_method in c("Zmean", "ssGSEA")) {
      resid_col <- paste0("STEP_", step_method, "_residual_CAF3")
      residual <- score_df[[resid_col]]
      cutoff <- median(residual, na.rm = TRUE)
      group <- ifelse(residual <= cutoff, "Low", "High")
      group[!is.finite(residual)] <- NA
      group <- factor(group, levels = c("Low", "High"))

      hm <- hallmark_method_list[[HALLMARK_SCORE_METHOD]]
      if (is.null(hm)) stop("Requested HALLMARK_SCORE_METHOD was not calculated: ", HALLMARK_SCORE_METHOD)
      for (pathway in rownames(hm)) {
        y <- as.numeric(hm[pathway, ids])
        wt <- safe_wilcox(y, group, MIN_GROUP_N)
        hl_rows[[kk]] <- data.frame(
          cancer = cancer,
          STEP_method = step_method,
          Hallmark_method = HALLMARK_SCORE_METHOD,
          pathway = pathway,
          cutoff = cutoff,
          wt,
          stringsAsFactors = FALSE
        )
        box_data_list[[kk]] <- data.frame(
          cancer = cancer,
          STEP_method = step_method,
          pathway = pathway,
          score = y,
          group = group,
          stringsAsFactors = FALSE
        )
        kk <- kk + 1L
      }
    }
    hl_df <- bind_rows(hl_rows) %>%
      group_by(cancer, STEP_method, Hallmark_method) %>%
      mutate(FDR = p.adjust(p, method = "BH"), sig = p_star(FDR)) %>%
      ungroup()
    fwrite(hl_df, file.path(OUT_DIR, "07_Hallmark_high_low", paste0(cancer, "_STEP_HighLow_Hallmark_Wilcoxon.csv")))
    all_hallmark_highlow[[cancer]] <- hl_df

    box_all <- bind_rows(box_data_list) %>%
      left_join(hl_df %>% select(cancer, STEP_method, pathway, p, FDR),
                by = c("cancer", "STEP_method", "pathway")) %>%
      mutate(pathway_label = clean_pathway_name(pathway))

    for (step_method in c("Zmean", "ssGSEA")) {
      dd <- box_all %>% filter(STEP_method == step_method, is.finite(score))
      if (nrow(dd) == 0) next
      ann <- dd %>%
        distinct(pathway, pathway_label, p, FDR) %>%
        group_by(pathway, pathway_label) %>%
        summarise(p = first(p), FDR = first(FDR), .groups = "drop")
      yr <- dd %>% group_by(pathway, pathway_label) %>%
        summarise(y = max(score, na.rm = TRUE), ymin = min(score, na.rm = TRUE), .groups = "drop")
      ann <- left_join(ann, yr, by = c("pathway", "pathway_label")) %>%
        mutate(y = y + pmax((y - ymin) * 0.10, 0.05),
               label = paste0(fmt_p(p), "\nFDR = ", ifelse(is.na(FDR), "NA", signif(FDR, 2))))
      pbox <- ggplot(dd, aes(group, score, fill = group)) +
        geom_boxplot(width = 0.62, outlier.size = 0.35, linewidth = 0.35, na.rm = TRUE) +
        geom_jitter(width = 0.12, size = 0.25, alpha = 0.35, na.rm = TRUE) +
        geom_text(data = ann, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 2.1) +
        scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
        facet_wrap(~pathway_label, scales = "free_y", ncol = 5) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
        theme_pub(9) +
        theme(legend.position = "none", strip.text = element_text(size = 6.5),
              axis.text.x = element_text(size = 6.5)) +
        labs(
          title = paste0(cancer, ": Hallmark ", HALLMARK_SCORE_METHOD,
                         " by CAF3-adjusted STEP ", step_method, " High/Low"),
          x = NULL, y = paste0("Hallmark ", HALLMARK_SCORE_METHOD, " score")
        )
      save_plot(
        pbox,
        file.path(OUT_DIR, "07_Hallmark_high_low",
                  paste0(cancer, "_", step_method, "_all_Hallmark_HighLow_boxplots")),
        width = 15, height = 18
      )
    }

    if (RUN_GENE_LEVEL_CORRELATION) {
      gene_methods <- list(
        Zmean = score_df$STEP_Zmean_residual_CAF3,
        ssGSEA = score_df$STEP_ssGSEA_residual_CAF3
      )
      for (step_method in names(gene_methods)) {
        gc <- vectorized_gene_spearman(expr, gene_methods[[step_method]])
        gc$FDR <- p.adjust(gc$p, method = "BH")
        gc$cancer <- cancer
        gc$method <- step_method
        fwrite(gc, file.path(OUT_DIR, "09_gene_correlations",
                             paste0(cancer, "_", step_method, "_CAF3adjusted_STEP_gene_correlations.csv")))
        all_gene_cor[[paste(cancer, step_method, sep = "__")]] <- gc[, c("gene", "n", "rho", "p", "FDR", "cancer", "method")]

        top_lab <- gc %>%
          filter(is.finite(rho), is.finite(FDR)) %>%
          arrange(FDR, desc(abs(rho))) %>%
          slice_head(n = 12)
        pv <- ggplot(gc, aes(rho, -log10(pmax(FDR, 1e-300)))) +
          geom_point(aes(color = case_when(
            FDR < 0.05 & rho > 0 ~ "Positive",
            FDR < 0.05 & rho < 0 ~ "Negative",
            TRUE ~ "NS"
          )), size = 0.65, alpha = 0.65) +
          geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
          geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35) +
          ggrepel::geom_text_repel(data = top_lab, aes(label = gene), size = 2.4,
                                   max.overlaps = 20, min.segment.length = 0) +
          scale_color_manual(values = c("Negative" = NEG_COLOR, "NS" = "grey75", "Positive" = POS_COLOR)) +
          theme_pub(10) +
          labs(
            title = paste0(cancer, ": genes correlated with CAF3-adjusted STEP ", step_method),
            x = "Spearman rho", y = expression(-log[10](FDR)), color = NULL
          )
        save_plot(pv, file.path(OUT_DIR, "09_gene_correlations",
                                paste0(cancer, "_", step_method, "_gene_correlation_volcano")),
                  width = 7.5, height = 6)

        if (RUN_PER_CANCER_GSEA) {
          rank_vec <- gc$rho
          names(rank_vec) <- gc$gene
          for (collection in names(enrichment_collections)) {
            fg <- run_fgsea_one(rank_vec, enrichment_collections[[collection]],
                                collection, cancer, step_method)
            if (!is.null(fg)) {
              fwrite(fg, file.path(OUT_DIR, "11_GSEA",
                                   paste0(cancer, "_", step_method, "_", collection, "_preranked_GSEA.csv")))
              all_gsea[[paste(cancer, step_method, collection, sep = "__")]] <- fg
            }
          }
        }
      }
    }

    diag <- data.frame(
      cancer = cancer,
      expression_file = f,
      n_genes = nrow(expr),
      n_primary_samples = ncol(expr),
      gene_column = read_info$gene_col,
      mapping_mode = read_info$mapping_mode,
      candidate_annotation_overlaps = read_info$candidate_overlaps,
      gene_column_reference_overlap = read_info$overlap,
      gct_skip = read_info$gct_skip,
      log2_transformed = read_info$log2_transformed,
      q99_before_transform = read_info$q99_before_transform,
      stringsAsFactors = FALSE
    )
    fwrite(diag, file.path(OUT_DIR, "00_logs", paste0(cancer, "_expression_diagnostics.csv")))
    TRUE
  }, error = function(e) {
    msg("ERROR in ", cancer, ": ", conditionMessage(e))
    fwrite(data.frame(cancer = cancer, file = f, error = conditionMessage(e)),
           file.path(OUT_DIR, "00_logs", paste0(cancer, "_ERROR.csv")))
    FALSE
  })
}

if (length(all_scores) == 0) stop("No cancer completed successfully. Check 00_logs.")

scores_all <- bind_rows(all_scores)
qc_all <- bind_rows(all_qc)
resid_params_all <- bind_rows(all_resid_params)
cor_all <- bind_rows(all_hallmark_cor)
hl_all <- bind_rows(all_hallmark_highlow)
concordance_all <- bind_rows(all_concordance)

fwrite(scores_all, file.path(OUT_DIR, "14_figure_ready", "ALL_cancers_STEP_CAF3_Hallmark_scores.csv"))
fwrite(qc_all, file.path(OUT_DIR, "01_gene_set_QC", "ALL_cancers_gene_set_QC.csv"))
fwrite(resid_params_all, file.path(OUT_DIR, "03_CAF3_residualization", "ALL_cancers_CAF3_residualization_alpha_beta.csv"))
fwrite(cor_all, file.path(OUT_DIR, "04_Hallmark_correlations", "ALL_cancers_STEP_Hallmark_correlations.csv"))
fwrite(hl_all, file.path(OUT_DIR, "07_Hallmark_high_low", "ALL_cancers_STEP_HighLow_Hallmark_Wilcoxon.csv"))
fwrite(concordance_all, file.path(OUT_DIR, "08_method_concordance", "ALL_cancers_Zmean_vs_ssGSEA_corrected_score_concordance.csv"))

plot_cor_heatmap <- function(df, step_method, adjustment, hall_method, stem) {
  dd <- df %>%
    filter(.data$STEP_method == .env$step_method,
           .data$adjustment == .env$adjustment,
           .data$Hallmark_method == .env$hall_method) %>%
    mutate(
      pathway_label = clean_pathway_name(pathway),
      cancer = as.character(cancer)
    )
  if (nrow(dd) == 0) return(NULL)

  default_pathway_order <- dd %>%
    group_by(pathway_label) %>%
    summarise(mean_abs = mean(abs(rho), na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_abs) %>%
    pull(pathway_label)
  default_cancer_order <- sort(unique(dd$cancer))

  cl <- cluster_long_heatmap(
    dd, row_var = "pathway_label", col_var = "cancer", value_var = "rho",
    default_row_order = default_pathway_order,
    default_col_order = default_cancer_order,
    stem = stem
  )
  dd <- cl$data

  p <- ggplot(dd, aes(cancer, pathway_label, fill = rho)) +
    geom_tile(color = "grey90", linewidth = 0.15) +
    geom_text(aes(label = sig), size = 2.0, color = "black") +
    scale_fill_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR,
                         midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                         na.value = "grey92") +
    theme_pub(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 6.8), panel.grid = element_blank()) +
    labs(
      title = paste0("TCGA Hallmark correlations: ", step_method, " STEP (", adjustment, ")"),
      subtitle = paste0("Hallmark score = ", hall_method, "; stars use within-cancer BH FDR"),
      x = "Cancer type", y = NULL, fill = "Spearman\nrho"
    )
  save_plot(p, stem, width = 15, height = 12)
}

for (step_method in c("Zmean", "ssGSEA")) {
  for (adjustment in c("raw", "residual")) {
    plot_cor_heatmap(
      cor_all, step_method, adjustment, HALLMARK_SCORE_METHOD,
      file.path(OUT_DIR, "05_Hallmark_heatmaps",
                paste0("PanCancer_", step_method, "_", adjustment, "_STEP_vs_Hallmark_heatmap"))
    )
  }
}

delta_cor <- cor_all %>%
  filter(Hallmark_method == HALLMARK_SCORE_METHOD) %>%
  select(cancer, STEP_method, adjustment, pathway, rho, p, FDR) %>%
  pivot_wider(names_from = adjustment, values_from = c(rho, p, FDR)) %>%
  mutate(delta_rho = rho_residual - rho_raw)
fwrite(delta_cor, file.path(OUT_DIR, "04_Hallmark_correlations", "ALL_cancers_Hallmark_delta_rho_corrected_minus_raw.csv"))

for (step_method in c("Zmean", "ssGSEA")) {
  dd <- delta_cor %>% filter(STEP_method == step_method) %>%
    mutate(pathway_label = clean_pathway_name(pathway),
           cancer = as.character(cancer))
  if (nrow(dd)) {
    order_y <- dd %>% group_by(pathway_label) %>%
      summarise(mean_abs = mean(abs(delta_rho), na.rm = TRUE), .groups = "drop") %>%
      arrange(mean_abs) %>% pull(pathway_label)
    delta_stem <- file.path(OUT_DIR, "05_Hallmark_heatmaps",
                            paste0("PanCancer_", step_method, "_Hallmark_delta_rho_heatmap"))
    cl <- cluster_long_heatmap(
      dd, "pathway_label", "cancer", "delta_rho",
      default_row_order = order_y,
      default_col_order = sort(unique(dd$cancer)),
      stem = delta_stem
    )
    dd <- cl$data
    lim <- max(abs(dd$delta_rho), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) lim <- 0.5
    p <- ggplot(dd, aes(cancer, pathway_label, fill = delta_rho)) +
      geom_tile(color = "grey90", linewidth = 0.15) +
      scale_fill_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR,
                           midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                           na.value = "grey92") +
      theme_pub(9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.text.y = element_text(size = 6.8), panel.grid = element_blank()) +
      labs(
        title = paste0("CAF3 adjustment effect on Hallmark correlations: ", step_method),
        subtitle = expression(Delta*rho == rho[adjusted] - rho[raw]),
        x = "Cancer type", y = NULL, fill = expression(Delta*rho)
      )
    save_plot(p, delta_stem, width = 15, height = 12)
  }
}

meta_cor <- cor_all %>%
  group_by(STEP_method, adjustment, Hallmark_method, pathway) %>%
  group_modify(~ random_effects_meta_cor(.x$rho, .x$n)) %>%
  ungroup() %>%
  group_by(STEP_method, adjustment, Hallmark_method) %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  ungroup()

sig_counts <- cor_all %>%
  group_by(STEP_method, adjustment, Hallmark_method, pathway) %>%
  summarise(
    positive_sig_cancers = sum(FDR < 0.05 & rho > 0, na.rm = TRUE),
    negative_sig_cancers = sum(FDR < 0.05 & rho < 0, na.rm = TRUE),
    total_sig_cancers = positive_sig_cancers + negative_sig_cancers,
    tested_cancers = sum(is.finite(rho)),
    .groups = "drop"
  )
meta_cor <- left_join(meta_cor, sig_counts,
                      by = c("STEP_method", "adjustment", "Hallmark_method", "pathway"))
fwrite(meta_cor, file.path(OUT_DIR, "06_Hallmark_meta_summary", "Hallmark_random_effects_meta_analysis.csv"))

meta_res <- meta_cor %>%
  filter(adjustment == "residual", Hallmark_method == HALLMARK_SCORE_METHOD) %>%
  mutate(pathway_label = clean_pathway_name(pathway))
order_meta <- meta_res %>% group_by(pathway_label) %>%
  summarise(max_abs = max(abs(meta_rho), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs) %>% pull(pathway_label)
meta_res$pathway_label <- factor(meta_res$pathway_label, levels = order_meta)

pmeta <- ggplot(meta_res, aes(meta_rho, pathway_label, color = meta_rho, size = total_sig_cancers)) +
  geom_segment(aes(x = 0, xend = meta_rho, yend = pathway_label), color = "grey75", linewidth = 0.45) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
  scale_color_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0) +
  facet_wrap(~STEP_method, ncol = 2) +
  theme_pub(9) +
  theme(axis.text.y = element_text(size = 6.7)) +
  labs(
    title = "Pan-cancer Hallmark associations with CAF3-adjusted STEP",
    subtitle = "Random-effects meta-analysis; point size = cancers with within-cancer FDR < 0.05",
    x = "Meta-analysis Spearman rho", y = NULL,
    color = "Meta rho", size = "Significant\ncancers"
  )
save_plot(pmeta, file.path(OUT_DIR, "06_Hallmark_meta_summary", "CAF3adjusted_STEP_Hallmark_meta_lollipop"),
          width = 13, height = 12)

meta_wide <- meta_cor %>%
  filter(Hallmark_method == HALLMARK_SCORE_METHOD) %>%
  select(STEP_method, adjustment, pathway, meta_rho, FDR) %>%
  pivot_wider(names_from = adjustment, values_from = c(meta_rho, FDR)) %>%
  mutate(pathway_label = clean_pathway_name(pathway),
         delta_meta_rho = meta_rho_residual - meta_rho_raw)
fwrite(meta_wide, file.path(OUT_DIR, "06_Hallmark_meta_summary", "Hallmark_meta_raw_vs_CAF3adjusted.csv"))

for (step_method in c("Zmean", "ssGSEA")) {
  dd <- meta_wide %>% filter(STEP_method == step_method) %>%
    arrange(delta_meta_rho) %>%
    mutate(pathway_label = factor(pathway_label, levels = pathway_label))
  pd <- ggplot(dd) +
    geom_segment(aes(x = meta_rho_raw, xend = meta_rho_residual,
                     y = pathway_label, yend = pathway_label), color = "grey65", linewidth = 0.65) +
    geom_point(aes(meta_rho_raw, pathway_label, color = "Raw"), size = 2.1) +
    geom_point(aes(meta_rho_residual, pathway_label, color = "CAF3-adjusted"), size = 2.1) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
    scale_color_manual(values = c("Raw" = "grey45", "CAF3-adjusted" = POS_COLOR)) +
    theme_pub(9) +
    theme(axis.text.y = element_text(size = 6.7), legend.position = "top") +
    labs(
      title = paste0("Hallmark meta-correlations before and after CAF3 adjustment: ", step_method),
      x = "Meta-analysis Spearman rho", y = NULL, color = NULL
    )
  save_plot(pd, file.path(OUT_DIR, "06_Hallmark_meta_summary",
                          paste0(step_method, "_Hallmark_raw_vs_CAF3adjusted_dumbbell")),
            width = 10, height = 12)
}

count_df <- sig_counts %>%
  filter(adjustment == "residual", Hallmark_method == HALLMARK_SCORE_METHOD) %>%
  mutate(pathway_label = clean_pathway_name(pathway),
         negative_plot = -negative_sig_cancers)
for (step_method in c("Zmean", "ssGSEA")) {
  dd <- count_df %>% filter(STEP_method == step_method) %>%
    arrange(positive_sig_cancers + negative_sig_cancers) %>%
    mutate(pathway_label = factor(pathway_label, levels = pathway_label))
  pc <- ggplot(dd) +
    geom_col(aes(positive_sig_cancers, pathway_label, fill = "Positive"), width = 0.72) +
    geom_col(aes(negative_plot, pathway_label, fill = "Negative"), width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    scale_fill_manual(values = c("Negative" = NEG_COLOR, "Positive" = POS_COLOR)) +
    scale_x_continuous(labels = abs) +
    theme_pub(9) +
    theme(axis.text.y = element_text(size = 6.7), legend.position = "top") +
    labs(
      title = paste0("Cancer types with significant Hallmark correlation: ", step_method, " CAF3-adjusted STEP"),
      subtitle = "Within-cancer BH FDR < 0.05",
      x = "Number of cancer types", y = NULL, fill = NULL
    )
  save_plot(pc, file.path(OUT_DIR, "06_Hallmark_meta_summary",
                          paste0(step_method, "_Hallmark_significant_cancer_counts")),
            width = 10, height = 12)
}

for (step_method in c("Zmean", "ssGSEA")) {
  dd <- hl_all %>%
    filter(STEP_method == step_method, Hallmark_method == HALLMARK_SCORE_METHOD) %>%
    mutate(pathway_label = clean_pathway_name(pathway),
           cancer = as.character(cancer))
  order_y <- dd %>% group_by(pathway_label) %>%
    summarise(mean_abs = mean(abs(median_diff), na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_abs) %>% pull(pathway_label)
  hl_stem <- file.path(OUT_DIR, "07_Hallmark_high_low",
                       paste0("PanCancer_", step_method, "_Hallmark_HighLow_median_difference_heatmap"))
  cl <- cluster_long_heatmap(
    dd, "pathway_label", "cancer", "median_diff",
    default_row_order = order_y,
    default_col_order = sort(unique(dd$cancer)),
    stem = hl_stem
  )
  dd <- cl$data
  lim <- max(abs(dd$median_diff), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  ph <- ggplot(dd, aes(cancer, pathway_label, fill = median_diff)) +
    geom_tile(color = "grey90", linewidth = 0.15) +
    geom_text(aes(label = sig), size = 2.0) +
    scale_fill_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR,
                         midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                         na.value = "grey92") +
    theme_pub(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 6.8), panel.grid = element_blank()) +
    labs(
      title = paste0("Hallmark ", HALLMARK_SCORE_METHOD, ": STEP High minus Low (", step_method, ")"),
      subtitle = "Low=blue, High=red in boxplots; stars use within-cancer BH FDR",
      x = "Cancer type", y = NULL, fill = "Median\ndifference"
    )
  save_plot(ph, hl_stem, width = 15, height = 12)
}

hl_meta <- hl_all %>%
  group_by(STEP_method, Hallmark_method, pathway) %>%
  summarise(
    mean_rank_biserial = weighted.mean(rank_biserial, w = pmax(n_low + n_high, 1), na.rm = TRUE),
    positive_sig_cancers = sum(FDR < 0.05 & rank_biserial > 0, na.rm = TRUE),
    negative_sig_cancers = sum(FDR < 0.05 & rank_biserial < 0, na.rm = TRUE),
    tested_cancers = sum(is.finite(rank_biserial)),
    .groups = "drop"
  )
fwrite(hl_meta, file.path(OUT_DIR, "07_Hallmark_high_low", "PanCancer_Hallmark_HighLow_effect_summary.csv"))

concordance_all$FDR <- p.adjust(concordance_all$p, method = "BH")
concordance_all <- concordance_all %>% arrange(rho) %>%
  mutate(cancer = factor(cancer, levels = cancer))
fwrite(concordance_all, file.path(OUT_DIR, "08_method_concordance", "ALL_cancers_Zmean_vs_ssGSEA_corrected_score_concordance_with_FDR.csv"))

pconc <- ggplot(concordance_all, aes(rho, cancer)) +
  geom_segment(aes(x = 0, xend = rho, yend = cancer), color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = rho), size = 2.8) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
  scale_color_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0,
                        limits = c(-1, 1), oob = scales::squish) +
  theme_pub(10) +
  labs(
    title = "Concordance of CAF3-adjusted STEP scores",
    subtitle = "Spearman correlation: Z-mean residual versus ssGSEA residual within each cancer",
    x = "Spearman rho", y = NULL, color = "rho"
  )
save_plot(pconc, file.path(OUT_DIR, "08_method_concordance", "Zmean_vs_ssGSEA_CAF3adjusted_STEP_concordance"),
          width = 8, height = 9)

meta_gene_all <- NULL
gene_dt <- NULL
if (RUN_GENE_LEVEL_CORRELATION && length(all_gene_cor) > 0) {
  gene_dt <- rbindlist(all_gene_cor, fill = TRUE)
  meta_gene_all <- compute_gene_meta_table(gene_dt)
  fwrite(meta_gene_all, file.path(OUT_DIR, "10_gene_meta_analysis", "PanCancer_gene_correlation_meta_analysis_fixed_and_random.csv"))
  
  fwrite(meta_gene_all, file.path(OUT_DIR, "10_gene_meta_analysis", "PanCancer_gene_correlation_meta_analysis.csv"))

  for (step_method in c("Zmean", "ssGSEA")) {
    dd <- meta_gene_all %>% filter(method == step_method)
    dd_eligible <- dd %>% filter(eligible_for_pan_gsea)
    top_lab <- dd_eligible %>% arrange(FDR, desc(abs(meta_rho))) %>% slice_head(n = 20)
    pgv <- ggplot(dd, aes(meta_rho, -log10(pmax(FDR, 1e-300)))) +
      geom_point(aes(color = case_when(
        eligible_for_pan_gsea & FDR < 0.05 & meta_rho > 0 ~ "Positive",
        eligible_for_pan_gsea & FDR < 0.05 & meta_rho < 0 ~ "Negative",
        TRUE ~ "NS/ineligible"
      )), size = 0.65, alpha = 0.65) +
      geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
      geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35) +
      ggrepel::geom_text_repel(data = top_lab, aes(label = gene), size = 2.5,
                               max.overlaps = 30, min.segment.length = 0) +
      scale_color_manual(values = c(
        "Negative" = NEG_COLOR,
        "NS/ineligible" = "grey75",
        "Positive" = POS_COLOR
      )) +
      theme_pub(10) +
      labs(
        title = paste0("Pan-cancer genes associated with CAF3-adjusted STEP ", step_method),
        subtitle = paste0(
          stringr::str_to_title(PAN_GSEA_PRIMARY_META_MODEL),
          "-effects Fisher-z meta-analysis; pan-GSEA-eligible genes satisfy cancer-coverage threshold"
        ),
        x = "Primary meta-analysis Spearman rho",
        y = expression(-log[10](FDR)), color = NULL
      )
    save_plot(pgv, file.path(OUT_DIR, "10_gene_meta_analysis",
                             paste0(step_method, "_PanCancer_gene_meta_volcano")),
              width = 8, height = 6.5)

    hi <- dd_eligible %>% arrange(FDR, desc(abs(meta_rho))) %>% slice_head(n = 50) %>%
      mutate(gene = factor(gene, levels = rev(gene)))
    if (nrow(hi)) {
      phi <- ggplot(hi, aes(I2, gene, color = meta_rho, size = -log10(pmax(FDR, 1e-300)))) +
        geom_point(alpha = 0.9) +
        scale_color_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0) +
        theme_pub(9) +
        theme(axis.text.y = element_text(size = 6.5)) +
        labs(
          title = paste0("Heterogeneity of top STEP-associated genes: ", step_method),
          x = expression(I^2~"(%)"), y = NULL, color = "Meta rho",
          size = expression(-log[10](FDR))
        )
      save_plot(phi, file.path(OUT_DIR, "10_gene_meta_analysis",
                               paste0(step_method, "_top50_gene_meta_I2_plot")),
                width = 9, height = 10)
    }

    top_genes <- dd_eligible %>% arrange(FDR, desc(abs(meta_rho))) %>% slice_head(n = 50) %>% pull(gene)
    gh <- as.data.frame(gene_dt[method == step_method & gene %in% top_genes])
    gh$gene <- toupper(as.character(gh$gene))
    gh$cancer <- as.character(gh$cancer)
    gh$sig <- p_star(gh$FDR)
    if (nrow(gh)) {
      gene_hm_stem <- file.path(OUT_DIR, "10_gene_meta_analysis",
                                paste0(step_method, "_top50_gene_by_cancer_heatmap"))
      cl <- cluster_long_heatmap(
        gh, "gene", "cancer", "rho",
        default_row_order = top_genes,
        default_col_order = sort(unique(gh$cancer)),
        stem = gene_hm_stem
      )
      gh <- cl$data
      pg <- ggplot(gh, aes(cancer, gene, fill = rho)) +
        geom_tile(color = "grey90", linewidth = 0.15) +
        geom_text(aes(label = sig), size = 1.7) +
        scale_fill_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR,
                             midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                             na.value = "grey92") +
        theme_pub(9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
              axis.text.y = element_text(size = 6.2), panel.grid = element_blank()) +
        labs(
          title = paste0("Top pan-cancer genes: CAF3-adjusted STEP ", step_method),
          subtitle = "Stars use per-cancer gene-correlation BH FDR",
          x = "Cancer type", y = NULL, fill = "Spearman\nrho"
        )
      save_plot(pg, gene_hm_stem, width = 15, height = 12)
    }
  }
}

gsea_all <- if (length(all_gsea)) bind_rows(all_gsea) else data.frame()
if (nrow(gsea_all)) {
  fwrite(gsea_all, file.path(OUT_DIR, "11_GSEA", "ALL_cancers_per_cancer_preranked_GSEA.csv"))

  for (collection in unique(gsea_all$collection)) {
    for (step_method in c("Zmean", "ssGSEA")) {
      dd <- gsea_all %>% filter(.data$collection == .env$collection, .data$method == .env$step_method)
      if (nrow(dd) == 0) next
      top_paths <- dd %>% group_by(pathway) %>%
        summarise(sig_n = sum(padj < 0.05, na.rm = TRUE),
                  mean_abs_NES = mean(abs(NES), na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(sig_n), desc(mean_abs_NES)) %>% slice_head(n = 30) %>% pull(pathway)
      pp <- dd %>% filter(pathway %in% top_paths) %>%
        mutate(pathway_label = clean_pathway_name(pathway),
               pathway_label = factor(pathway_label, levels = rev(clean_pathway_name(top_paths))),
               cancer = factor(cancer, levels = sort(unique(cancer))),
               neglogFDR = -log10(pmax(padj, 1e-300)))
      pb <- ggplot(pp, aes(cancer, pathway_label)) +
        geom_point(aes(size = neglogFDR, color = NES), alpha = 0.9) +
        scale_color_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0) +
        scale_size_continuous(range = c(0.5, 5)) +
        theme_pub(9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
              axis.text.y = element_text(size = 6.5)) +
        labs(
          title = paste0(collection, " GSEA across TCGA: CAF3-adjusted STEP ", step_method),
          x = "Cancer type", y = NULL, color = "NES", size = expression(-log[10](FDR))
        )
      save_plot(pb, file.path(OUT_DIR, "11_GSEA",
                              paste0(collection, "_", step_method, "_perCancer_GSEA_bubble")),
                width = 15, height = 10)
    }
  }
}

pan_gsea_all <- list()
if (RUN_PAN_CANCER_GSEA && !is.null(meta_gene_all)) {
  variant_specs <- build_pan_variant_specs()
  fwrite(variant_specs, file.path(OUT_DIR, "15_sensitivity", "pan_cancer_GSEA_variant_specifications.csv"))

  excluded_score_genes <- unique(core_sets$STEP_UP)
  if (isTRUE(REMOVE_CAF3_GENES_IN_DEOVERLAP)) {
    excluded_score_genes <- unique(c(excluded_score_genes, core_sets$CAF3))
  }
  fwrite(
    data.frame(gene = excluded_score_genes, reason = "Removed in de-overlap sensitivity"),
    file.path(OUT_DIR, "15_sensitivity", "genes_removed_in_deoverlap_sensitivity.csv")
  )

  for (step_method in c("Zmean", "ssGSEA")) {
    md0 <- meta_gene_all %>% filter(method == step_method, eligible_for_pan_gsea)
    if (nrow(md0) == 0) {
      msg("No pan-GSEA-eligible genes for method: ", step_method)
      next
    }

    for (vi in seq_len(nrow(variant_specs))) {
      spec <- variant_specs[vi, , drop = FALSE]
      rank_col <- paste0(spec$meta_model, "_", spec$rank_metric)
      if (!rank_col %in% colnames(md0)) {
        msg("Skipping unavailable rank column: ", rank_col)
        next
      }

      md <- md0 %>% filter(is.finite(.data[[rank_col]]))
      stats_vec <- md[[rank_col]]
      names(stats_vec) <- md$gene
      if (isTRUE(spec$deoverlap)) {
        stats_vec <- stats_vec[!names(stats_vec) %in% excluded_score_genes]
      }
      stats_vec <- stats_vec[is.finite(stats_vec)]
      stats_vec <- stats_vec[!duplicated(names(stats_vec))]
      stats_vec <- stats_vec[order(-stats_vec, names(stats_vec))]
      if (length(stats_vec) < 100) next

      rank_label <- paste0(
        "Gene rank by ", spec$meta_model, "-effects ",
        gsub("_", "-", spec$rank_metric),
        ": positive STEP association  →  negative association"
      )

      for (collection in names(enrichment_collections)) {
        pathways_use <- enrichment_collections[[collection]]
        if (isTRUE(spec$deoverlap)) {
          pathways_use <- lapply(pathways_use, function(gs) setdiff(gs, excluded_score_genes))
        }

        fg <- run_fgsea_one(stats_vec, pathways_use, collection, "PAN_CANCER_META", step_method)
        if (is.null(fg)) next
        fg$analysis_variant <- spec$analysis_variant
        fg$meta_model <- spec$meta_model
        fg$rank_metric <- spec$rank_metric
        fg$deoverlap <- spec$deoverlap
        fg$ranked_gene_n <- length(stats_vec)
        fg$minimum_cancers_required <- min(md$min_required_cancers, na.rm = TRUE)

        file_stem <- paste0(
          "PAN_CANCER_META_", step_method, "_", collection, "_", spec$analysis_variant
        )
        fwrite(fg, file.path(OUT_DIR, "11_GSEA", paste0(file_stem, "_GSEA.csv")))
        if (vi == 1L) {
          
          fwrite(fg, file.path(OUT_DIR, "11_GSEA",
                               paste0("PAN_CANCER_META_", step_method, "_", collection, "_GSEA.csv")))
        }
        pan_gsea_all[[paste(step_method, collection, spec$analysis_variant, sep = "__")]] <- fg

        fg_plot <- fg
        if (identical(collection, "REACTOME") && isTRUE(REACTOME_COLLAPSE_FOR_PLOTS)) {
          fg_plot <- fg_plot %>% filter(is_main_pathway)
        }

        sig <- fg_plot %>% filter(is.finite(padj), padj < GSEA_DOTPLOT_FDR_CUTOFF, is.finite(NES))
        top_pos <- sig %>% filter(NES > 0) %>% arrange(padj, desc(abs(NES))) %>%
          slice_head(n = TOP_DOTPLOT_PER_DIRECTION)
        top_neg <- sig %>% filter(NES < 0) %>% arrange(padj, desc(abs(NES))) %>%
          slice_head(n = TOP_DOTPLOT_PER_DIRECTION)
        top <- bind_rows(top_pos, top_neg) %>%
          mutate(
            pathway_label = clean_pathway_name(pathway),
            pathway_label = factor(pathway_label, levels = pathway_label[order(NES)]),
            fdr_star = p_star(padj)
          )
        if (nrow(top)) {
          pgsea <- ggplot(top, aes(NES, pathway_label, color = NES,
                                   size = -log10(pmax(padj, 1e-300)))) +
            geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
            geom_point() +
            geom_text(aes(label = fdr_star), color = "black", size = 2.4, vjust = -1.1) +
            scale_color_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0) +
            theme_pub(9) +
            theme(axis.text.y = element_text(size = 6.8)) +
            labs(
              title = paste0("Pan-cancer meta-gene GSEA: ", collection, " / STEP ", step_method),
              subtitle = paste0(
                spec$analysis_variant, " | significant pathways only: BH FDR < ",
                GSEA_DOTPLOT_FDR_CUTOFF
              ),
              x = "NES", y = NULL, color = "NES", size = expression(-log[10](FDR))
            )
          save_plot(pgsea, file.path(OUT_DIR, "11_GSEA", paste0(file_stem, "_significant_GSEA_dotplot")),
                    width = 11, height = 9)
        }

        if (isTRUE(spec$plot_curves)) {
          save_pan_curve_outputs(
            fg = fg,
            stats_vec = stats_vec,
            pathways = pathways_use,
            step_method = step_method,
            collection = collection,
            analysis_variant = spec$analysis_variant,
            rank_label = rank_label
          )
        }
      }
    }
  }
}

if (length(pan_gsea_all)) {
  pan_gsea_df <- bind_rows(pan_gsea_all)
  fwrite(pan_gsea_df, file.path(OUT_DIR, "11_GSEA", "ALL_pan_cancer_meta_GSEA_all_variants.csv"))
  primary_variant_name <- build_pan_variant_specs()$analysis_variant[1]
  pan_primary <- pan_gsea_df %>% filter(analysis_variant == primary_variant_name)
  fwrite(pan_primary, file.path(OUT_DIR, "11_GSEA", "ALL_pan_cancer_meta_GSEA.csv"))
}

if (isTRUE(RUN_LEAVE_ONE_CANCER_OUT) && !is.null(gene_dt) && nrow(gene_dt)) {
  loco_all <- list()
  for (step_method in c("Zmean", "ssGSEA")) {
    cancers_method <- sort(unique(gene_dt$`cancer`[gene_dt$`method` == step_method]))
    for (left_out in cancers_method) {
      dt_sub <- gene_dt[method == step_method & cancer != left_out]
      meta_sub <- compute_gene_meta_table(dt_sub)
      md <- meta_sub %>% filter(method == step_method, eligible_for_pan_gsea)
      rank_col <- paste0(PAN_GSEA_PRIMARY_META_MODEL, "_", PAN_GSEA_PRIMARY_RANK_METRIC)
      md <- md %>% filter(is.finite(.data[[rank_col]]))
      stats_vec <- md[[rank_col]]
      names(stats_vec) <- md$gene
      stats_vec <- stats_vec[!duplicated(names(stats_vec))]
      stats_vec <- stats_vec[order(-stats_vec, names(stats_vec))]

      for (collection in intersect(LOCO_COLLECTIONS, names(enrichment_collections))) {
        fg <- run_fgsea_one(stats_vec, enrichment_collections[[collection]],
                            collection, paste0("LEAVE_OUT_", left_out), step_method)
        if (is.null(fg)) next
        fg$left_out_cancer <- left_out
        loco_all[[paste(step_method, left_out, collection, sep = "__")]] <- fg
      }
    }
  }
  if (length(loco_all)) {
    loco_df <- bind_rows(loco_all)
    fwrite(loco_df, file.path(OUT_DIR, "15_sensitivity", "leave_one_cancer_out_GSEA_full.csv"))
    loco_core <- loco_df %>% filter(pathway %in% LOCO_CORE_PATHWAYS) %>%
      mutate(pathway_label = clean_pathway_name(pathway),
             left_out_cancer = as.character(left_out_cancer),
             sig = p_star(padj))
    fwrite(loco_core, file.path(OUT_DIR, "15_sensitivity", "leave_one_cancer_out_core_pathways.csv"))
    if (nrow(loco_core)) {
      loco_stem <- file.path(OUT_DIR, "15_sensitivity", "leave_one_cancer_out_core_pathway_NES_heatmap")

      cl <- cluster_long_heatmap(
        loco_core, "pathway_label", "left_out_cancer", "NES",
        default_row_order = clean_pathway_name(LOCO_CORE_PATHWAYS),
        default_col_order = sort(unique(loco_core$left_out_cancer)),
        stem = loco_stem
      )
      loco_core <- cl$data
      ploco <- ggplot(loco_core, aes(left_out_cancer, pathway_label, fill = NES)) +
        geom_tile(color = "grey90", linewidth = 0.15) +
        geom_text(aes(label = sig), size = 2.0) +
        facet_wrap(~method, scales = "free_x") +
        scale_fill_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0) +
        theme_pub(9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
              axis.text.y = element_text(size = 7)) +
        labs(
          title = "Leave-one-cancer-out robustness of core pan-cancer GSEA pathways",
          subtitle = "Each column shows the result after removing the indicated cancer; stars use BH FDR",
          x = "Cancer removed from meta-analysis", y = NULL, fill = "NES"
        )
      save_plot(ploco, loco_stem, width = 16, height = 8)
    }
  }
}

ora_all <- list()
if (RUN_ORA && !is.null(meta_gene_all)) {
  for (step_method in c("Zmean", "ssGSEA")) {
    md <- meta_gene_all %>% filter(method == step_method, eligible_for_pan_gsea, is.finite(meta_rho), is.finite(FDR))
    universe <- md$gene
    selected_pos <- md %>% filter(FDR < ORA_META_FDR_CUTOFF, meta_rho >= ORA_ABS_RHO_CUTOFF) %>% pull(gene)
    selected_neg <- md %>% filter(FDR < ORA_META_FDR_CUTOFF, meta_rho <= -ORA_ABS_RHO_CUTOFF) %>% pull(gene)
    for (collection in names(enrichment_collections)) {
      op <- run_ora_one(selected_pos, universe, enrichment_collections[[collection]],
                        collection, "Positive", step_method)
      on <- run_ora_one(selected_neg, universe, enrichment_collections[[collection]],
                        collection, "Negative", step_method)
      oo <- bind_rows(op, on)
      if (nrow(oo) == 0) next
      fwrite(oo, file.path(OUT_DIR, "12_ORA",
                           paste0(step_method, "_", collection, "_panCancer_meta_gene_ORA.csv")))
      ora_all[[paste(step_method, collection, sep = "__")]] <- oo

      top <- oo %>% filter(FDR < ORA_PLOT_FDR_CUTOFF) %>% group_by(direction) %>%
        arrange(FDR, desc(overlap)) %>% slice_head(n = 15) %>% ungroup() %>%
        mutate(pathway_label = clean_pathway_name(pathway),
               signed_ratio = ifelse(direction == "Positive", gene_ratio, -gene_ratio))
      if (nrow(top)) {
        top$pathway_label <- factor(top$pathway_label, levels = rev(unique(top$pathway_label)))
        po <- ggplot(top, aes(signed_ratio, pathway_label, size = overlap, color = -log10(pmax(FDR, 1e-300)))) +
          geom_vline(xintercept = 0, linewidth = 0.35) +
          geom_point() +
          scale_x_continuous(labels = abs) +
          scale_color_gradient(low = "grey70", high = POS_COLOR) +
          theme_pub(9) +
          theme(axis.text.y = element_text(size = 6.8)) +
          labs(
            title = paste0("ORA of pan-cancer STEP-associated genes: ", step_method, " / ", collection),
            subtitle = paste0("BH FDR < ", ORA_PLOT_FDR_CUTOFF, "; left=negative-associated genes; right=positive-associated genes"),
            x = "Gene ratio", y = NULL, size = "Overlap", color = expression(-log[10](FDR))
          )
        save_plot(po, file.path(OUT_DIR, "12_ORA",
                                paste0(step_method, "_", collection, "_panCancer_ORA_dotplot")),
                  width = 11, height = 9)
      }
    }
  }
}
if (length(ora_all)) fwrite(bind_rows(ora_all), file.path(OUT_DIR, "12_ORA", "ALL_panCancer_meta_gene_ORA.csv"))

make_leading_long <- function(df, source_label) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  if (!"analysis_variant" %in% colnames(df)) df$analysis_variant <- source_label
  df %>%
    filter(padj < 0.05, nzchar(leadingEdge)) %>%
    select(method, collection, cancer, pathway, NES, padj, analysis_variant, leadingEdge) %>%
    separate_rows(leadingEdge, sep = ";") %>%
    rename(gene = leadingEdge) %>%
    filter(nzchar(gene)) %>%
    mutate(source = source_label)
}

summarise_leading <- function(lead_df) {
  if (nrow(lead_df) == 0) return(data.frame())
  lead_df %>%
    group_by(source, analysis_variant, method, collection, gene) %>%
    summarise(
      n_cancers = n_distinct(cancer[cancer != "PAN_CANCER_META"]),
      n_pathways = n_distinct(pathway),
      total_occurrences = n(),
      positive_NES_occurrences = sum(NES > 0),
      negative_NES_occurrences = sum(NES < 0),
      pathways = paste(sort(unique(pathway)), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(source, analysis_variant, method, collection,
            desc(total_occurrences), desc(n_cancers))
}

per_lead_long <- make_leading_long(gsea_all, "PER_CANCER")
pan_all_df <- if (length(pan_gsea_all)) bind_rows(pan_gsea_all) else data.frame()
primary_variant_name <- build_pan_variant_specs()$analysis_variant[1]
pan_primary_df <- if (nrow(pan_all_df)) pan_all_df %>% filter(analysis_variant == primary_variant_name) else data.frame()
pan_primary_lead_long <- make_leading_long(pan_primary_df, "PAN_CANCER_PRIMARY")
pan_sensitivity_lead_long <- if (nrow(pan_all_df)) {
  make_leading_long(pan_all_df %>% filter(analysis_variant != primary_variant_name), "PAN_CANCER_SENSITIVITY")
} else data.frame()

if (nrow(per_lead_long)) {
  fwrite(per_lead_long, file.path(OUT_DIR, "13_leading_edge", "PER_CANCER_significant_GSEA_leading_edge_long.csv"))
  fwrite(summarise_leading(per_lead_long),
         file.path(OUT_DIR, "13_leading_edge", "PER_CANCER_leading_edge_gene_frequency_summary.csv"))
}
if (nrow(pan_primary_lead_long)) {
  fwrite(pan_primary_lead_long,
         file.path(OUT_DIR, "13_leading_edge", "PAN_CANCER_PRIMARY_significant_GSEA_leading_edge_long.csv"))
  fwrite(summarise_leading(pan_primary_lead_long),
         file.path(OUT_DIR, "13_leading_edge", "PAN_CANCER_PRIMARY_leading_edge_gene_frequency_summary.csv"))
}
if (nrow(pan_sensitivity_lead_long)) {
  fwrite(pan_sensitivity_lead_long,
         file.path(OUT_DIR, "13_leading_edge", "PAN_CANCER_SENSITIVITY_significant_GSEA_leading_edge_long.csv"))
  fwrite(summarise_leading(pan_sensitivity_lead_long),
         file.path(OUT_DIR, "13_leading_edge", "PAN_CANCER_SENSITIVITY_leading_edge_gene_frequency_summary.csv"))
}

lead_sig <- bind_rows(per_lead_long, pan_primary_lead_long)
if (nrow(lead_sig)) {
  lead_summary <- summarise_leading(lead_sig)
  fwrite(lead_sig, file.path(OUT_DIR, "13_leading_edge", "ALL_significant_GSEA_leading_edge_long.csv"))
  fwrite(lead_summary, file.path(OUT_DIR, "13_leading_edge", "leading_edge_gene_frequency_summary.csv"))

  plot_summary <- lead_summary %>% filter(source == "PER_CANCER")
  for (collection in unique(plot_summary$collection)) {
    for (step_method in c("Zmean", "ssGSEA")) {
      dd <- plot_summary %>%
        filter(.data$collection == .env$collection, .data$method == .env$step_method) %>%
        slice_head(n = 30) %>% arrange(total_occurrences) %>%
        mutate(gene = factor(gene, levels = gene))
      if (nrow(dd) == 0) next
      pl <- ggplot(dd, aes(total_occurrences, gene, fill = n_cancers)) +
        geom_col(width = 0.72) +
        scale_fill_gradient(low = "grey80", high = POS_COLOR) +
        theme_pub(9) +
        labs(
          title = paste0("Per-cancer leading-edge genes: ", collection, " / STEP ", step_method),
          x = "Occurrences in significant per-cancer GSEA results", y = NULL,
          fill = "Cancer\ntypes"
        )
      save_plot(pl, file.path(OUT_DIR, "13_leading_edge",
                              paste0(collection, "_", step_method, "_top30_perCancer_leading_edge_genes")),
                width = 9, height = 8)
    }
  }
}

if (isTRUE(RUN_GENE_CORE_INTEGRATION) &&
    !is.null(gene_dt) && nrow(gene_dt) > 0 &&
    exists("meta_gene_all") && nrow(meta_gene_all) > 0) {

  msg("Running integrated gene prioritization and article-ready heatmaps ...")

  core_out <- file.path(OUT_DIR, "16_gene_core_integration")
  dir.create(core_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "01_integrated_tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "02_top_positive_negative_heatmaps"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "03_functional_core_heatmaps"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "04_EMT_leading_edge_heatmaps"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "05_method_concordance"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "06_pathway_specific_leading_edge"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "07_gene_forest_plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(core_out, "08_evidence_summary_plots"), recursive = TRUE, showWarnings = FALSE)

  core_pstar <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", ""))))
  }

  core_clean_pathway <- function(x) {
    x <- gsub("^HALLMARK_", "", x)
    x <- gsub("^REACTOME_", "", x)
    stringr::str_to_title(gsub("_", " ", x))
  }

  core_assign_modules <- function(gene_vector) {
    module_map <- bind_rows(lapply(names(GENE_FUNCTIONAL_MODULES), function(module_nm) {
      data.frame(
        gene = toupper(GENE_FUNCTIONAL_MODULES[[module_nm]]),
        functional_module = module_nm,
        stringsAsFactors = FALSE
      )
    }))
    module_map <- module_map %>%
      group_by(gene) %>%
      summarise(
        functional_module = paste(unique(functional_module), collapse = ";"),
        .groups = "drop"
      )

    data.frame(
      gene = unique(toupper(as.character(gene_vector))),
      stringsAsFactors = FALSE
    ) %>%
      left_join(module_map, by = "gene", relationship = "many-to-one") %>%
      mutate(functional_module = ifelse(is.na(functional_module), "Other", functional_module)) %>%
      distinct(gene, .keep_all = TRUE)
  }

  core_meta_ci <- function(meta_rho, meta_zstat) {
    est_z <- atanh(pmin(pmax(meta_rho, -0.999999), 0.999999))
    se_z <- ifelse(is.finite(meta_zstat) & abs(meta_zstat) > 1e-12,
                   abs(est_z / meta_zstat), NA_real_)
    data.frame(
      meta_ci_low = tanh(est_z - 1.96 * se_z),
      meta_ci_high = tanh(est_z + 1.96 * se_z)
    )
  }

  core_heatmap_order <- function(gene_list, method_names = c("Zmean", "ssGSEA")) {
    gene_list <- unique(toupper(gene_list))
    dd0 <- as.data.frame(
      gene_dt[method %in% method_names & toupper(gene) %in% gene_list]
    )
    if (nrow(dd0) == 0) return(NULL)
    dd0$gene <- toupper(as.character(dd0$gene))
    dd0$cancer <- as.character(dd0$cancer)
    cluster_long_heatmap(
      dd0, "gene", "cancer", "rho",
      default_row_order = gene_list,
      default_col_order = sort(unique(dd0$cancer)),
      stem = NULL
    )
  }

  core_basic_heatmap <- function(gene_list, method_name, title_text,
                                 show_stars = FALSE, row_modules = NULL,
                                 order_obj = NULL) {
    gene_list <- unique(toupper(gene_list))
    dd <- as.data.frame(gene_dt[method == method_name & toupper(gene) %in% gene_list])
    if (nrow(dd) == 0) return(NULL)
    dd$gene <- toupper(as.character(dd$gene))
    dd$cancer <- as.character(dd$cancer)

    if (is.null(order_obj)) {
      order_obj <- core_heatmap_order(gene_list, method_names = method_name)
    }
    if (is.null(order_obj)) return(NULL)

    cancer_levels <- order_obj$col_order
    gene_levels_top_to_bottom <- order_obj$row_order
    dd <- tidyr::complete(dd, gene = gene_list, cancer = cancer_levels)
    dd$gene <- factor(dd$gene, levels = rev(gene_levels_top_to_bottom))
    dd$cancer <- factor(dd$cancer, levels = cancer_levels)
    dd$sig <- if (show_stars) core_pstar(dd$FDR) else ""

    p <- ggplot(dd, aes(cancer, gene, fill = rho)) +
      geom_tile(color = "grey90", linewidth = 0.15) +
      scale_fill_gradient2(
        low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR,
        midpoint = 0, limits = c(-1, 1), oob = scales::squish,
        na.value = "grey92"
      ) +
      theme_pub(9) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        plot.title = element_text(size = 12, face = "bold")
      ) +
      labs(title = title_text, x = "Cancer type", y = NULL, fill = "Spearman\nrho")
    if (show_stars) p <- p + geom_text(aes(label = sig), size = 1.7)
    p
  }

  core_annotated_heatmap <- function(gene_list, method_name, title_text,
                                     show_stars = GENE_CORE_SHOW_HEATMAP_STARS) {
    gene_list <- unique(toupper(gene_list))
    order_obj <- core_heatmap_order(gene_list, method_names = method_name)
    if (is.null(order_obj)) return(NULL)
    evidence_method <- integrated_gene_evidence %>%
      filter(method == method_name, gene %in% gene_list) %>%
      mutate(gene = factor(gene, levels = rev(order_obj$row_order)))
    if (nrow(evidence_method) == 0) return(NULL)

    p_heat <- core_basic_heatmap(
      gene_list, method_name, title_text, show_stars, order_obj = order_obj
    )
    if (is.null(p_heat)) return(NULL)

    y_theme <- theme_pub(8) +
      theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5)
      )

    p_meta <- ggplot(evidence_method, aes(meta_rho, gene, color = meta_rho)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
      geom_segment(aes(x = 0, xend = meta_rho, yend = gene), linewidth = 0.45) +
      geom_point(size = 2.1) +
      scale_color_gradient2(low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0) +
      y_theme + theme(legend.position = "none") +
      labs(title = "Meta rho", x = NULL)

    p_dir <- ggplot(evidence_method, aes(direction_consistency, gene)) +
      geom_vline(xintercept = GENE_CORE_MIN_DIRECTION_CONSISTENCY,
                 linetype = 2, color = "grey55") +
      geom_point(aes(color = dominant_direction), size = 2.1) +
      scale_color_manual(values = c("Positive" = POS_COLOR, "Negative" = NEG_COLOR, "Tie" = "grey60")) +
      scale_x_continuous(limits = c(0.45, 1), breaks = c(0.5, 0.75, 1)) +
      y_theme + theme(legend.position = "none") +
      labs(title = "Direction", x = "Consistency")

    p_i2 <- ggplot(evidence_method, aes(I2, gene, color = I2)) +
      geom_vline(xintercept = c(50, 75), linetype = c(3, 2), color = "grey55") +
      geom_point(size = 2.1) +
      scale_color_gradient(low = "#2E8B57", high = "#B22222", limits = c(0, 100), oob = scales::squish) +
      scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
      y_theme + theme(legend.position = "none") +
      labs(title = expression(I^2), x = "%")

    max_le <- max(evidence_method$hallmark_le_cancers, na.rm = TRUE)
    if (!is.finite(max_le) || max_le < 1) max_le <- 1
    p_le <- ggplot(evidence_method, aes(hallmark_le_cancers, gene)) +
      geom_point(aes(size = hallmark_le_occurrences, color = emt_le_cancers), alpha = 0.9) +
      scale_color_gradient(low = "grey75", high = POS_COLOR) +
      scale_size_continuous(range = c(1.8, 4.5)) +
      scale_x_continuous(limits = c(0, max_le), breaks = pretty(c(0, max_le), n = 3)) +
      y_theme + theme(legend.position = "none") +
      labs(title = "Leading edge", x = "Hallmark cancers")

    p_heat + p_meta + p_dir + p_i2 + p_le +
      patchwork::plot_layout(widths = c(8.2, 1.5, 1.4, 1.25, 1.5), guides = "collect")
  }

  core_make_forest_plot <- function(gene_name, method_name) {
    dd <- as.data.frame(gene_dt[method == method_name & toupper(gene) == toupper(gene_name)])
    if (nrow(dd) == 0) return(NULL)
    dd <- dd %>%
      mutate(
        rho_clip = pmin(pmax(rho, -0.999999), 0.999999),
        z = atanh(rho_clip),
        se = 1 / sqrt(pmax(n - 3, 1)),
        ci_low = tanh(z - 1.96 * se),
        ci_high = tanh(z + 1.96 * se)
      ) %>%
      arrange(rho)

    mm <- integrated_gene_evidence %>%
      filter(method == method_name, gene == toupper(gene_name)) %>%
      slice_head(n = 1)
    if (nrow(mm) == 0) return(NULL)
    meta_ci <- core_meta_ci(mm$meta_rho, mm$meta_zstat)
    meta_row <- data.frame(
      cancer = "Random-effects meta",
      n = sum(dd$n, na.rm = TRUE),
      rho = mm$meta_rho,
      ci_low = meta_ci$meta_ci_low,
      ci_high = meta_ci$meta_ci_high,
      FDR = mm$FDR,
      stringsAsFactors = FALSE
    )
    plot_dt <- bind_rows(dd %>% select(cancer, n, rho, ci_low, ci_high, FDR), meta_row)
    plot_dt$cancer <- factor(plot_dt$cancer, levels = c(as.character(dd$cancer), "Random-effects meta"))
    plot_dt$is_meta <- as.character(plot_dt$cancer) == "Random-effects meta"

    subtitle_text <- paste0(
      "Meta rho = ", sprintf("%.3f", mm$meta_rho),
      " | BH FDR = ", format.pval(mm$FDR, digits = 3),
      " | I2 = ", sprintf("%.1f%%", mm$I2),
      " | direction consistency = ", sprintf("%.1f%%", 100 * mm$direction_consistency)
    )

    ggplot(plot_dt, aes(rho, cancer)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
      geom_segment(aes(x = ci_low, xend = ci_high, y = cancer, yend = cancer, color = is_meta), linewidth = 0.55) +
      geom_point(aes(size = ifelse(is_meta, max(n, na.rm = TRUE), n), color = is_meta), alpha = 0.95) +
      scale_color_manual(values = c("FALSE" = "grey35", "TRUE" = ifelse(mm$meta_rho >= 0, POS_COLOR, NEG_COLOR))) +
      scale_size_continuous(range = c(1.8, 4.5), guide = "none") +
      coord_cartesian(xlim = c(-1, 1)) +
      theme_pub(9) +
      theme(legend.position = "none", axis.text.y = element_text(size = 7)) +
      labs(
        title = paste0(gene_name, " vs CAF3-adjusted STEP ", method_name),
        subtitle = subtitle_text,
        x = "Spearman rho (95% CI)", y = NULL
      )
  }

  direction_stats <- as.data.frame(gene_dt)[, c("method", "gene", "cancer", "rho", "FDR")] %>%
    mutate(gene = toupper(gene)) %>%
    filter(is.finite(rho)) %>%
    group_by(method, gene) %>%
    summarise(
      positive_cancers = sum(rho > 0, na.rm = TRUE),
      negative_cancers = sum(rho < 0, na.rm = TRUE),
      zero_cancers = sum(rho == 0, na.rm = TRUE),
      significant_cancers = sum(FDR < 0.05, na.rm = TRUE),
      positive_significant_cancers = sum(FDR < 0.05 & rho > 0, na.rm = TRUE),
      negative_significant_cancers = sum(FDR < 0.05 & rho < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      direction_consistency = pmax(positive_cancers, negative_cancers) /
        pmax(positive_cancers + negative_cancers, 1),
      significant_direction_consistency = pmax(positive_significant_cancers, negative_significant_cancers) /
        pmax(positive_significant_cancers + negative_significant_cancers, 1),
      dominant_direction = case_when(
        positive_cancers > negative_cancers ~ "Positive",
        negative_cancers > positive_cancers ~ "Negative",
        TRUE ~ "Tie"
      )
    )

  hallmark_le_summary <- data.frame()
  emt_le_summary <- data.frame()
  core_pathway_le_summary <- data.frame()
  if (exists("per_lead_long") && nrow(per_lead_long) > 0) {
    lead_use <- per_lead_long %>%
      mutate(gene = toupper(gene)) %>%
      filter(cancer != "PAN_CANCER_META")

    hallmark_le_summary <- lead_use %>%
      filter(collection == "HALLMARK") %>%
      group_by(method, gene) %>%
      summarise(
        hallmark_le_cancers = n_distinct(cancer),
        hallmark_le_pathways = n_distinct(pathway),
        hallmark_le_occurrences = n(),
        .groups = "drop"
      )

    emt_le_summary <- lead_use %>%
      filter(collection == "HALLMARK", pathway == EMT_HALLMARK_PATHWAY) %>%
      group_by(method, gene) %>%
      summarise(
        emt_le_cancers = n_distinct(cancer),
        emt_le_occurrences = n(),
        emt_mean_NES = mean(NES, na.rm = TRUE),
        .groups = "drop"
      )

    core_pathway_le_summary <- lead_use %>%
      filter(collection == "HALLMARK", pathway %in% CORE_HALLMARK_PATHWAYS) %>%
      group_by(method, gene) %>%
      summarise(
        core_pathway_le_cancers = n_distinct(cancer),
        core_pathway_le_pathways = n_distinct(pathway),
        core_pathway_le_occurrences = n(),
        .groups = "drop"
      )
  }

  meta_use <- meta_gene_all %>%
    mutate(
      method = as.character(method),
      gene = toupper(as.character(gene))
    ) %>%
    arrange(method, gene, FDR) %>%
    group_by(method, gene) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    left_join(direction_stats, by = c("method", "gene"), relationship = "one-to-one") %>%
    left_join(hallmark_le_summary, by = c("method", "gene"), relationship = "one-to-one") %>%
    left_join(emt_le_summary, by = c("method", "gene"), relationship = "one-to-one") %>%
    left_join(core_pathway_le_summary, by = c("method", "gene"), relationship = "one-to-one")

  numeric_zero_cols <- c(
    "hallmark_le_cancers", "hallmark_le_pathways", "hallmark_le_occurrences",
    "emt_le_cancers", "emt_le_occurrences", "core_pathway_le_cancers",
    "core_pathway_le_pathways", "core_pathway_le_occurrences"
  )
  for (cc in numeric_zero_cols) {
    if (!cc %in% colnames(meta_use)) meta_use[[cc]] <- 0
    meta_use[[cc]][is.na(meta_use[[cc]])] <- 0
  }

  step_members <- toupper(core_sets$STEP_UP)
  caf3_members <- toupper(core_sets$CAF3)
  functional_modules <- core_assign_modules(meta_use$gene)

  integrated_gene_evidence <- meta_use %>%
    left_join(functional_modules, by = "gene", relationship = "many-to-one") %>%
    distinct(method, gene, .keep_all = TRUE) %>%
    mutate(
      is_STEP_UP_member = gene %in% step_members,
      is_CAF3_member = gene %in% caf3_members,
      is_user_functional_core = gene %in% toupper(USER_FUNCTIONAL_CORE_GENES),
      hallmark_le_cancer_fraction = hallmark_le_cancers / pmax(total_tested_cancers, 1),
      emt_le_cancer_fraction = emt_le_cancers / pmax(total_tested_cancers, 1),
      I2_available = is.finite(I2),
      low_heterogeneity_score = ifelse(is.finite(I2), pmax(0, 1 - I2 / 100), 0),
      eligible_by_meta = is.finite(FDR) & FDR < GENE_CORE_META_FDR_CUTOFF &
        coverage_fraction >= GENE_CORE_MIN_COVERAGE_FRACTION &
        direction_consistency >= GENE_CORE_MIN_DIRECTION_CONSISTENCY,
      strict_core_candidate = eligible_by_meta & is.finite(I2) &
        I2 <= GENE_CORE_STRICT_MAX_I2 &
        hallmark_le_cancers >= GENE_CORE_MIN_HALLMARK_LE_CANCERS,
      evidence_count =
        as.integer(is.finite(FDR) & FDR < GENE_CORE_META_FDR_CUTOFF) +
        as.integer(coverage_fraction >= GENE_CORE_MIN_COVERAGE_FRACTION) +
        as.integer(direction_consistency >= GENE_CORE_MIN_DIRECTION_CONSISTENCY) +
        as.integer(is.finite(I2) & I2 <= GENE_CORE_STRICT_MAX_I2) +
        as.integer(hallmark_le_cancers >= GENE_CORE_MIN_HALLMARK_LE_CANCERS) +
        as.integer(emt_le_cancers >= GENE_CORE_MIN_EMT_LE_CANCERS) +
        as.integer(is_STEP_UP_member) +
        as.integer(is_user_functional_core),
      figure_priority_score =
        3 * abs(meta_rho) + direction_consistency + coverage_fraction +
        low_heterogeneity_score + hallmark_le_cancer_fraction +
        0.5 * emt_le_cancer_fraction + 0.25 * as.numeric(is_STEP_UP_member) +
        0.25 * as.numeric(is_user_functional_core)
    )

  fwrite(
    integrated_gene_evidence,
    file.path(core_out, "01_integrated_tables", "Integrated_gene_evidence_all_methods.csv")
  )
  for (method_name in c("Zmean", "ssGSEA")) {
    fwrite(
      integrated_gene_evidence %>% filter(method == method_name),
      file.path(core_out, "01_integrated_tables", paste0(method_name, "_Integrated_gene_evidence.csv"))
    )
    fwrite(
      integrated_gene_evidence %>%
        filter(method == method_name, strict_core_candidate) %>%
        arrange(desc(evidence_count), desc(abs(meta_rho)), I2),
      file.path(core_out, "01_integrated_tables", paste0(method_name, "_Strict_core_candidates.csv"))
    )
  }

  selected_by_method <- list()
  for (method_name in c("Zmean", "ssGSEA")) {
    ev <- integrated_gene_evidence %>%
      filter(method == method_name, eligible_by_meta, is.finite(meta_rho))

    pos <- ev %>%
      filter(meta_rho > 0) %>%
      arrange(desc(evidence_count), desc(hallmark_le_cancers), desc(abs(meta_rho)), I2) %>%
      slice_head(n = GENE_CORE_TOP_PER_DIRECTION) %>%
      mutate(selection_direction = "Positive")
    neg <- ev %>%
      filter(meta_rho < 0) %>%
      arrange(desc(evidence_count), desc(hallmark_le_cancers), desc(abs(meta_rho)), I2) %>%
      slice_head(n = GENE_CORE_TOP_PER_DIRECTION) %>%
      mutate(selection_direction = "Negative")
    selected <- bind_rows(pos, neg)
    selected_by_method[[method_name]] <- selected
    fwrite(
      selected,
      file.path(core_out, "01_integrated_tables",
                paste0(method_name, "_Top", GENE_CORE_TOP_PER_DIRECTION,
                       "_positive_and_negative_multiEvidence_genes.csv"))
    )

    gene_order <- c(pos$gene, neg$gene)
    if (length(gene_order)) {
      p_clean <- core_annotated_heatmap(
        gene_order, method_name,
        paste0("Multi-evidence top genes vs CAF3-adjusted STEP ", method_name),
        show_stars = FALSE
      )
      if (!is.null(p_clean)) {
        save_plot(
          p_clean,
          file.path(core_out, "02_top_positive_negative_heatmaps",
                    paste0(method_name, "_Top_positive_negative_multiEvidence_annotated_heatmap_clean")),
          width = 19, height = max(8, 0.28 * length(gene_order) + 3)
        )
      }
      if (isTRUE(GENE_CORE_SHOW_HEATMAP_STARS)) {
        p_star <- core_annotated_heatmap(
          gene_order, method_name,
          paste0("Multi-evidence top genes vs CAF3-adjusted STEP ", method_name,
                 " (per-cancer FDR stars)"),
          show_stars = TRUE
        )
        if (!is.null(p_star)) {
          save_plot(
            p_star,
            file.path(core_out, "02_top_positive_negative_heatmaps",
                      paste0(method_name, "_Top_positive_negative_multiEvidence_annotated_heatmap_withStars")),
            width = 19, height = max(8, 0.28 * length(gene_order) + 3)
          )
        }
      }
    }
  }

  functional_gene_order <- unique(unlist(GENE_FUNCTIONAL_MODULES, use.names = FALSE))
  functional_gene_order <- toupper(functional_gene_order)
  functional_gene_order <- functional_gene_order[
    functional_gene_order %in% unique(toupper(as.character(gene_dt$gene)))
  ]
  if (length(functional_gene_order)) {
    functional_table <- integrated_gene_evidence %>%
      filter(gene %in% functional_gene_order) %>%
      arrange(match(gene, functional_gene_order), method)
    fwrite(functional_table,
           file.path(core_out, "01_integrated_tables", "Functional_core_gene_evidence.csv"))

    method_heatmaps <- list()
    for (method_name in c("Zmean", "ssGSEA")) {
      pp <- core_annotated_heatmap(
        functional_gene_order, method_name,
        paste0("STEP functional core genes: ", method_name),
        show_stars = FALSE
      )
      if (!is.null(pp)) {
        method_heatmaps[[method_name]] <- pp
        save_plot(
          pp,
          file.path(core_out, "03_functional_core_heatmaps",
                    paste0(method_name, "_functional_core_gene_annotated_heatmap")),
          width = 19, height = max(7, 0.30 * length(functional_gene_order) + 3)
        )
      }
    }

    shared_functional_order <- core_heatmap_order(
      functional_gene_order, method_names = c("Zmean", "ssGSEA")
    )
    p_z_basic <- core_basic_heatmap(
      functional_gene_order, "Zmean", "Z-mean", FALSE,
      order_obj = shared_functional_order
    )
    p_s_basic <- core_basic_heatmap(
      functional_gene_order, "ssGSEA", "ssGSEA", FALSE,
      order_obj = shared_functional_order
    )
    if (!is.null(p_z_basic) && !is.null(p_s_basic)) {
      p_pair <- p_z_basic + p_s_basic +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(
          title = "Functional-core gene correlations with CAF3-adjusted STEP",
          subtitle = "Identical genes and cancer order; left=Z-mean, right=ssGSEA"
        ) & theme(legend.position = "right")
      save_plot(
        p_pair,
        file.path(core_out, "03_functional_core_heatmaps",
                  "Zmean_vs_ssGSEA_functional_core_gene_heatmaps"),
        width = 22, height = max(7, 0.30 * length(functional_gene_order) + 3)
      )
    }
  }

  if (exists("per_lead_long") && nrow(per_lead_long) > 0) {
    emt_le_gene_table <- per_lead_long %>%
      mutate(gene = toupper(gene)) %>%
      filter(collection == "HALLMARK", pathway == EMT_HALLMARK_PATHWAY,
             cancer != "PAN_CANCER_META") %>%
      group_by(method, gene) %>%
      summarise(
        emt_le_cancers = n_distinct(cancer),
        emt_le_occurrences = n(),
        emt_mean_NES = mean(NES, na.rm = TRUE),
        emt_min_FDR = min(padj, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(
        integrated_gene_evidence %>%
          select(method, gene, meta_rho, FDR, I2, coverage_fraction,
                 direction_consistency, hallmark_le_cancers,
                 is_STEP_UP_member, functional_module),
        by = c("method", "gene")
      ) %>%
      arrange(method, desc(emt_le_cancers), desc(emt_le_occurrences), desc(abs(meta_rho)))

    fwrite(emt_le_gene_table,
           file.path(core_out, "01_integrated_tables", "EMT_specific_leading_edge_gene_evidence.csv"))

    emt_top_by_method <- list()
    for (method_name in c("Zmean", "ssGSEA")) {
      emt_top <- emt_le_gene_table %>%
        filter(method == method_name) %>%
        arrange(desc(emt_le_cancers), desc(emt_le_occurrences), desc(abs(meta_rho))) %>%
        slice_head(n = EMT_LEADING_EDGE_TOP_N)
      emt_top_by_method[[method_name]] <- emt_top$gene
      if (nrow(emt_top)) {
        pp <- core_annotated_heatmap(
          emt_top$gene, method_name,
          paste0("Hallmark EMT leading-edge genes vs CAF3-adjusted STEP ", method_name),
          show_stars = FALSE
        )
        if (!is.null(pp)) {
          save_plot(
            pp,
            file.path(core_out, "04_EMT_leading_edge_heatmaps",
                      paste0(method_name, "_EMT_leading_edge_gene_annotated_heatmap")),
            width = 19, height = max(7, 0.32 * nrow(emt_top) + 3)
          )
        }
      }
    }

    emt_union <- unique(c(emt_top_by_method$Zmean, emt_top_by_method$ssGSEA))
    if (length(emt_union)) {
      emt_union_order <- emt_le_gene_table %>%
        filter(gene %in% emt_union) %>%
        group_by(gene) %>%
        summarise(max_cancers = max(emt_le_cancers, na.rm = TRUE),
                  total_occurrences = sum(emt_le_occurrences, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_cancers), desc(total_occurrences)) %>% pull(gene)
      shared_emt_order <- core_heatmap_order(
        emt_union_order, method_names = c("Zmean", "ssGSEA")
      )
      pz <- core_basic_heatmap(
        emt_union_order, "Zmean", "Z-mean", FALSE,
        order_obj = shared_emt_order
      )
      ps <- core_basic_heatmap(
        emt_union_order, "ssGSEA", "ssGSEA", FALSE,
        order_obj = shared_emt_order
      )
      if (!is.null(pz) && !is.null(ps)) {
        ppair <- pz + ps +
          patchwork::plot_layout(guides = "collect") +
          patchwork::plot_annotation(
            title = "Hallmark EMT leading-edge gene correlations",
            subtitle = "Union of the top EMT leading-edge genes from both STEP scoring methods"
          ) & theme(legend.position = "right")
        save_plot(
          ppair,
          file.path(core_out, "04_EMT_leading_edge_heatmaps",
                    "Zmean_vs_ssGSEA_EMT_leading_edge_gene_heatmaps"),
          width = 22, height = max(7, 0.30 * length(emt_union_order) + 3)
        )
      }
    }
  }

  if (exists("per_lead_long") && nrow(per_lead_long) > 0) {
    pathway_le_counts <- per_lead_long %>%
      mutate(gene = toupper(gene)) %>%
      filter(collection == "HALLMARK", pathway %in% CORE_HALLMARK_PATHWAYS,
             cancer != "PAN_CANCER_META") %>%
      group_by(method, pathway, gene) %>%
      summarise(
        n_cancers = n_distinct(cancer),
        mean_NES = mean(NES, na.rm = TRUE),
        signed_cancer_count = n_cancers * sign(mean_NES),
        .groups = "drop"
      )

    fwrite(pathway_le_counts,
           file.path(core_out, "01_integrated_tables", "Core_Hallmark_pathway_leading_edge_cancer_counts.csv"))

    for (method_name in c("Zmean", "ssGSEA")) {
      dd <- pathway_le_counts %>% filter(method == method_name)
      if (nrow(dd) == 0) next
      top_path_genes <- dd %>%
        group_by(gene) %>%
        summarise(
          total_cancers = sum(n_cancers, na.rm = TRUE),
          n_pathways = n_distinct(pathway),
          .groups = "drop"
        ) %>%
        arrange(desc(n_pathways), desc(total_cancers)) %>%
        slice_head(n = PATHWAY_LEADING_EDGE_TOP_N) %>%
        pull(gene)
      ppd <- dd %>%
        filter(gene %in% top_path_genes) %>%
        mutate(
          pathway_label = core_clean_pathway(pathway),
          pathway_label = as.character(pathway_label),
          gene = as.character(gene)
        )
      le_matrix_stem <- file.path(
        core_out, "06_pathway_specific_leading_edge",
        paste0(method_name, "_core_Hallmark_pathway_leading_edge_matrix")
      )
      cl <- cluster_long_heatmap(
        ppd, "gene", "pathway_label", "signed_cancer_count",
        default_row_order = top_path_genes,
        default_col_order = core_clean_pathway(CORE_HALLMARK_PATHWAYS),
        stem = le_matrix_stem
      )
      ppd <- cl$data
      max_count <- max(abs(ppd$signed_cancer_count), na.rm = TRUE)
      if (!is.finite(max_count) || max_count < 1) max_count <- 1
      pmat <- ggplot(ppd, aes(pathway_label, gene, fill = signed_cancer_count)) +
        geom_tile(color = "white", linewidth = 0.25) +
        geom_text(aes(label = n_cancers), size = 2.3) +
        scale_fill_gradient2(
          low = NEG_COLOR, mid = ZERO_COLOR, high = POS_COLOR, midpoint = 0,
          limits = c(-max_count, max_count), oob = scales::squish
        ) +
        theme_pub(9) +
        theme(
          axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 7), panel.grid = element_blank()
        ) +
        labs(
          title = paste0("Core Hallmark pathway-specific leading-edge genes: ", method_name),
          subtitle = "Cell values = number of cancer types; red=positive NES, blue=negative NES",
          x = NULL, y = NULL, fill = "Signed\ncancer count"
        )
      save_plot(
        pmat, le_matrix_stem,
        width = 12, height = max(7, 0.27 * length(top_path_genes) + 3)
      )
    }
  }

  concordance_value_cols <- c(
    "meta_rho", "FDR", "I2", "coverage_fraction",
    "direction_consistency", "hallmark_le_cancers",
    "emt_le_cancers", "evidence_count"
  )

  concordance_long <- integrated_gene_evidence %>%
    select(method, gene, all_of(concordance_value_cols)) %>%
    mutate(
      method = as.character(method),
      gene = toupper(as.character(gene)),
      across(
        all_of(concordance_value_cols),
        ~ suppressWarnings(as.numeric(as.character(.x)))
      )
    ) %>%
    group_by(method, gene) %>%
    summarise(
      across(
        all_of(concordance_value_cols),
        ~ {
          zz <- .x[is.finite(.x)]
          if (length(zz)) zz[[1]] else NA_real_
        }
      ),
      .groups = "drop"
    )

  duplicate_qc <- integrated_gene_evidence %>%
    count(method, gene, name = "n_rows") %>%
    filter(n_rows > 1)
  fwrite(
    duplicate_qc,
    file.path(
      core_out, "01_integrated_tables",
      "QC_duplicate_method_gene_rows_before_concordance.csv"
    )
  )

  concordance <- concordance_long %>%
    tidyr::pivot_wider(
      id_cols = gene,
      names_from = method,
      values_from = all_of(concordance_value_cols),
      names_sep = "__",
      values_fill = NA_real_
    ) %>%
    mutate(
      same_meta_direction = case_when(
        is.finite(meta_rho__Zmean) & is.finite(meta_rho__ssGSEA) ~
          sign(meta_rho__Zmean) == sign(meta_rho__ssGSEA),
        TRUE ~ NA
      ),
      mean_abs_meta_rho = rowMeans(
        cbind(abs(meta_rho__Zmean), abs(meta_rho__ssGSEA)),
        na.rm = TRUE
      ),
      mean_abs_meta_rho = ifelse(is.finite(mean_abs_meta_rho), mean_abs_meta_rho, NA_real_)
    )
  fwrite(concordance,
         file.path(core_out, "01_integrated_tables", "Zmean_vs_ssGSEA_gene_meta_concordance.csv"))

  concord_plot_dt <- concordance %>%
    filter(is.finite(meta_rho__Zmean), is.finite(meta_rho__ssGSEA))
  if (nrow(concord_plot_dt)) {
    ct <- suppressWarnings(cor.test(concord_plot_dt$meta_rho__Zmean,
                                    concord_plot_dt$meta_rho__ssGSEA,
                                    method = "spearman", exact = FALSE))
    pc <- ggplot(concord_plot_dt, aes(meta_rho__Zmean, meta_rho__ssGSEA)) +
      geom_hline(yintercept = 0, color = "grey70") +
      geom_vline(xintercept = 0, color = "grey70") +
      geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
      geom_point(aes(color = same_meta_direction), alpha = 0.45, size = 1.2) +
      scale_color_manual(values = c("TRUE" = POS_COLOR, "FALSE" = "black")) +
      theme_pub(10) +
      labs(
        title = "Gene-level meta-correlation concordance",
        subtitle = paste0("Spearman rho = ", sprintf("%.3f", unname(ct$estimate)),
                          "; P = ", format.pval(ct$p.value, digits = 3)),
        x = "Z-mean random-effects meta rho",
        y = "ssGSEA random-effects meta rho",
        color = "Same direction"
      )
    save_plot(pc,
              file.path(core_out, "05_method_concordance", "Zmean_vs_ssGSEA_gene_meta_rho_scatter"),
              width = 7.5, height = 6.5)
  }

  for (method_name in c("Zmean", "ssGSEA")) {
    selected <- selected_by_method[[method_name]]
    if (is.null(selected) || nrow(selected) == 0) next
    gene_order <- c(selected$gene[selected$meta_rho > 0], selected$gene[selected$meta_rho < 0])
    ev <- integrated_gene_evidence %>%
      filter(method == method_name, gene %in% gene_order) %>%
      mutate(gene = factor(gene, levels = rev(gene_order)))
    ev_long <- ev %>%
      transmute(
        gene,
        dominant_direction,
        `Absolute meta rho` = pmin(abs(meta_rho) / max(abs(meta_rho), na.rm = TRUE), 1),
        `Direction consistency` = direction_consistency,
        `Cancer coverage` = coverage_fraction,
        `Low heterogeneity` = low_heterogeneity_score,
        `Hallmark leading-edge coverage` = hallmark_le_cancer_fraction,
        `EMT leading-edge coverage` = emt_le_cancer_fraction
      ) %>%
      pivot_longer(-c(gene, dominant_direction), names_to = "evidence", values_to = "value")
    ev_long$value[!is.finite(ev_long$value)] <- 0
    pe <- ggplot(ev_long, aes(evidence, gene)) +
      geom_point(aes(size = value, fill = dominant_direction), shape = 21, color = "grey25") +
      scale_fill_manual(values = c("Positive" = POS_COLOR, "Negative" = NEG_COLOR, "Tie" = "grey65")) +
      scale_size_continuous(range = c(0.5, 6), limits = c(0, 1)) +
      theme_pub(9) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7), panel.grid.major = element_line(color = "grey92")
      ) +
      labs(
        title = paste0("Multi-evidence profile of prioritized genes: ", method_name),
        subtitle = "Bubble size is normalized evidence strength; color is dominant correlation direction",
        x = NULL, y = NULL, size = "Evidence", fill = "Direction"
      )
    save_plot(
      pe,
      file.path(core_out, "08_evidence_summary_plots",
                paste0(method_name, "_prioritized_gene_multiEvidence_dotplot")),
      width = 11, height = max(8, 0.28 * length(gene_order) + 3)
    )
  }

  if (isTRUE(GENE_CORE_MAKE_FOREST_PLOTS)) {
    for (method_name in c("Zmean", "ssGSEA")) {
      selected <- selected_by_method[[method_name]]
      if (is.null(selected) || nrow(selected) == 0) next
      forest_genes <- bind_rows(
        selected %>% filter(meta_rho > 0) %>% slice_head(n = GENE_CORE_FOREST_PER_DIRECTION),
        selected %>% filter(meta_rho < 0) %>% slice_head(n = GENE_CORE_FOREST_PER_DIRECTION)
      ) %>% pull(gene) %>% unique()
      method_dir <- file.path(core_out, "07_gene_forest_plots", method_name)
      dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
      for (gg in forest_genes) {
        pf <- core_make_forest_plot(gg, method_name)
        if (!is.null(pf)) {
          save_plot(
            pf,
            file.path(method_dir, paste0(method_name, "_", gsub("[^A-Za-z0-9_.-]", "_", gg), "_perCancer_correlation_forest")),
            width = 8.5, height = 9.5
          )
        }
      }
    }
  }

  msg("Integrated gene prioritization completed: ", core_out)
} else {
  msg("Integrated gene prioritization skipped: gene-level correlation/meta-analysis results are unavailable or RUN_GENE_CORE_INTEGRATION=FALSE.")
}

compact_hallmark <- cor_all %>%
  filter(adjustment == "residual", Hallmark_method == HALLMARK_SCORE_METHOD) %>%
  left_join(
    meta_cor %>% filter(adjustment == "residual", Hallmark_method == HALLMARK_SCORE_METHOD) %>%
      select(STEP_method, pathway, meta_rho, ci_low, ci_high, meta_p = p, meta_FDR = FDR,
             I2, positive_sig_cancers, negative_sig_cancers),
    by = c("STEP_method", "pathway")
  )
fwrite(compact_hallmark, file.path(OUT_DIR, "14_figure_ready", "FigureReady_CAF3adjusted_STEP_Hallmark_correlations.csv"))

settings <- data.frame(
  parameter = c(
    "STEP_SET_NAME", "CAF3_SET_NAME", "HALLMARK_SCORE_METHOD",
    "MIN_MATCHED_STEP_GENES", "MIN_MATCHED_CAF3_GENES", "MIN_MATCHED_HALLMARK_GENES",
    "MIN_SAMPLES_PER_CANCER", "MIN_COR_PAIRS", "MIN_GROUP_N",
    "MIN_CANCERS_FOR_GENE_META", "MIN_CANCER_COVERAGE_FRACTION",
    "PAN_GSEA_PRIMARY_META_MODEL", "PAN_GSEA_PRIMARY_RANK_METRIC",
    "GSEA_CURVE_FDR_CUTOFF", "GSEA_CURVE_MIN_ABS_NES",
    "TOP_CURVES_POSITIVE", "TOP_CURVES_NEGATIVE",
    "RUN_FIXED_EFFECT_META_Z_SENSITIVITY", "RUN_META_RHO_RANK_SENSITIVITY",
    "RUN_DEOVERLAP_SENSITIVITY", "RUN_LEAVE_ONE_CANCER_OUT",
    "RUN_GENE_LEVEL_CORRELATION", "RUN_PER_CANCER_GSEA", "RUN_PAN_CANCER_GSEA", "RUN_ORA",
    "RUN_GENE_CORE_INTEGRATION", "GENE_CORE_META_FDR_CUTOFF",
    "GENE_CORE_MIN_COVERAGE_FRACTION", "GENE_CORE_MIN_DIRECTION_CONSISTENCY",
    "GENE_CORE_STRICT_MAX_I2", "GENE_CORE_TOP_PER_DIRECTION",
    "GENE_CORE_MIN_HALLMARK_LE_CANCERS", "GENE_CORE_MIN_EMT_LE_CANCERS",
    "EMT_LEADING_EDGE_TOP_N", "PATHWAY_LEADING_EDGE_TOP_N",
    "HEATMAP_CLUSTER_ROWS", "HEATMAP_CLUSTER_COLUMNS",
    "HEATMAP_CLUSTER_DISTANCE", "HEATMAP_CLUSTER_METHOD",
    "HEATMAP_CLUSTER_NA_FILL", "HEATMAP_EXPORT_CLUSTER_ORDER",
    "LOW_COLOR", "HIGH_COLOR"
  ),
  value = c(
    step_name, caf3_name, HALLMARK_SCORE_METHOD,
    MIN_MATCHED_STEP_GENES, MIN_MATCHED_CAF3_GENES, MIN_MATCHED_HALLMARK_GENES,
    MIN_SAMPLES_PER_CANCER, MIN_COR_PAIRS, MIN_GROUP_N,
    MIN_CANCERS_FOR_GENE_META, MIN_CANCER_COVERAGE_FRACTION,
    PAN_GSEA_PRIMARY_META_MODEL, PAN_GSEA_PRIMARY_RANK_METRIC,
    GSEA_CURVE_FDR_CUTOFF, GSEA_CURVE_MIN_ABS_NES,
    TOP_CURVES_POSITIVE, TOP_CURVES_NEGATIVE,
    RUN_FIXED_EFFECT_META_Z_SENSITIVITY, RUN_META_RHO_RANK_SENSITIVITY,
    RUN_DEOVERLAP_SENSITIVITY, RUN_LEAVE_ONE_CANCER_OUT,
    RUN_GENE_LEVEL_CORRELATION, RUN_PER_CANCER_GSEA, RUN_PAN_CANCER_GSEA, RUN_ORA,
    RUN_GENE_CORE_INTEGRATION, GENE_CORE_META_FDR_CUTOFF,
    GENE_CORE_MIN_COVERAGE_FRACTION, GENE_CORE_MIN_DIRECTION_CONSISTENCY,
    GENE_CORE_STRICT_MAX_I2, GENE_CORE_TOP_PER_DIRECTION,
    GENE_CORE_MIN_HALLMARK_LE_CANCERS, GENE_CORE_MIN_EMT_LE_CANCERS,
    EMT_LEADING_EDGE_TOP_N, PATHWAY_LEADING_EDGE_TOP_N,
    HEATMAP_CLUSTER_ROWS, HEATMAP_CLUSTER_COLUMNS,
    HEATMAP_CLUSTER_DISTANCE, HEATMAP_CLUSTER_METHOD,
    HEATMAP_CLUSTER_NA_FILL, HEATMAP_EXPORT_CLUSTER_ORDER,
    LOW_COLOR, HIGH_COLOR
  ),
  stringsAsFactors = FALSE
)
fwrite(settings, file.path(OUT_DIR, "00_logs", "analysis_settings.csv"))

readme <- c(
  "TCGA STEP + Hallmark CAF3-adjusted analysis",
  "",
  paste0("Output directory: ", OUT_DIR),
  paste0("Expression directory: ", EXP_DIR),
  paste0("STEP GMT: ", STEP_GMT),
  paste0("Hallmark GMT: ", HALLMARK_GMT),
  paste0("STEP set: ", step_name),
  paste0("CAF3 set: ", caf3_name),
  "",
  "Score definitions:",
  "1. STEP_Zmean_raw: mean gene-wise z-scored expression of matched STEP genes.",
  "2. CAF3_Zmean: mean gene-wise z-scored expression of matched CAF3 genes.",
  "3. STEP_Zmean_residual_CAF3: residual from STEP_Zmean_raw ~ CAF3_Zmean within each cancer.",
  "4. STEP_ssGSEA_raw and CAF3_ssGSEA: GSVA/ssGSEA scores calculated from the same cancer expression matrix.",
  "5. STEP_ssGSEA_residual_CAF3: residual from STEP_ssGSEA_raw ~ CAF3_ssGSEA within each cancer.",
  "6. High/Low groups use the within-cancer median of the corresponding corrected STEP score.",
  "",
  "Main outputs:",
  "14_figure_ready/ALL_cancers_STEP_CAF3_Hallmark_scores.csv",
  "03_CAF3_residualization/ALL_cancers_CAF3_residualization_alpha_beta.csv",
  "04_Hallmark_correlations/ALL_cancers_STEP_Hallmark_correlations.csv",
  "05_Hallmark_heatmaps/*raw*, *residual* and *delta_rho* heatmaps",
  "06_Hallmark_meta_summary/Hallmark_random_effects_meta_analysis.csv",
  "07_Hallmark_high_low/ALL_cancers_STEP_HighLow_Hallmark_Wilcoxon.csv",
  "08_method_concordance/ALL_cancers_Zmean_vs_ssGSEA_corrected_score_concordance_with_FDR.csv",
  "09_gene_correlations/*gene_correlations.csv and volcano plots",
  "10_gene_meta_analysis/PanCancer_gene_correlation_meta_analysis_fixed_and_random.csv",
  "11_GSEA/*per-cancer and pan-cancer pre-ranked GSEA",
  "11_GSEA/*top6_positive*, *top6_negative* and *top6_positive_plus_top6_negative* curves",
  "15_sensitivity/*meta-model, meta-rho, de-overlap and optional leave-one-cancer-out results",
  "12_ORA/*pan-cancer ORA",
  "13_leading_edge/leading_edge_gene_frequency_summary.csv",
  "16_gene_core_integration/01_integrated_tables/Integrated_gene_evidence_all_methods.csv",
  "16_gene_core_integration/02_top_positive_negative_heatmaps/*multiEvidence* heatmaps",
  "16_gene_core_integration/03_functional_core_heatmaps/*functional_core* heatmaps",
  "16_gene_core_integration/04_EMT_leading_edge_heatmaps/*EMT_leading_edge* heatmaps",
  "16_gene_core_integration/06_pathway_specific_leading_edge/*leading_edge_matrix*",
  "16_gene_core_integration/07_gene_forest_plots/*perCancer_correlation_forest*",
  "",
  "Plot conventions:",
  paste0("Low group = blue ", LOW_COLOR),
  paste0("High group = red ", HIGH_COLOR),
  "Correlation/enrichment: negative=blue, zero=white, positive=red.",
  "P values and BH-adjusted FDR are exported; boxplots and pan-cancer GSEA curves display both nominal P and BH FDR.",
  paste0("Pan-cancer curves require BH FDR < ", GSEA_CURVE_FDR_CUTOFF,
         " and are selected separately as Top ", TOP_CURVES_POSITIVE,
         " positive + Top ", TOP_CURVES_NEGATIVE, " negative pathways."),
  "Heatmap stars represent BH FDR: * <0.05, ** <0.01, *** <0.001.",
  paste0("Heatmaps use hierarchical clustering: rows=", HEATMAP_CLUSTER_ROWS,
         ", columns=", HEATMAP_CLUSTER_COLUMNS,
         ", distance=", HEATMAP_CLUSTER_DISTANCE,
         ", linkage=", HEATMAP_CLUSTER_METHOD, "."),
  "NA values are replaced by HEATMAP_CLUSTER_NA_FILL only for clustering distances; plotted NA values remain unchanged.",
  "",
  "Important interpretation:",
  "Use Hallmark correlations with continuous corrected STEP as the primary statistical analysis.",
  "Use High/Low boxplots mainly for visualization because median splitting loses information.",
  "Compare raw versus corrected and delta-rho results to determine which associations are CAF-dependent or retained after CAF3 correction.",
  "Primary pan-cancer ranking uses the configured fixed/random-effects model and excludes genes failing the cancer-coverage threshold.",
  "The de-overlap sensitivity removes STEP_UP (and optionally CAF3) genes from both the rank and tested pathway gene sets.",
  "Reactome collapsePathways affects figure selection only; full uncollapsed GSEA tables remain available.",
  "Integrated gene prioritization uses meta FDR, cancer coverage, cross-cancer direction consistency, I2, Hallmark/EMT leading-edge recurrence and STEP membership as transparent evidence columns.",
  "The Top positive/negative heatmaps do not use I2 as a hard exclusion; I2 is shown as an annotation because biologically strong pan-cancer genes can remain heterogeneous."
)
writeLines(readme, file.path(OUT_DIR, "README_results.txt"))

msg("Analysis completed successfully.")
msg("Results saved to: ", OUT_DIR)
cat("Run finished: ", as.character(Sys.time()), "\n", file = LOG_FILE, append = TRUE)
