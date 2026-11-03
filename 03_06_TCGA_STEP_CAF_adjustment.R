# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
RUN_CLINICAL_MODULE <- TRUE
RUN_ADVANCED_CLINICAL_MODULE <- TRUE
RUN_ESTIMATE_EXPANSION <- TRUE
RUN_ALL_CANCER_KM_OVERVIEW <- TRUE
RUN_SURVIVAL_VOLCANO <- TRUE
RUN_TIDE_MAIN <- TRUE
RUN_TIDE_FIG7_EXTRA <- TRUE
RUN_IMMUNE_ACTIVITY_EXTRA <- TRUE
RUN_CHECKPOINT_MODULE <- TRUE
RUN_SUMMARY_BUBBLE <- TRUE
RUN_CHECKPOINT_SYMBOL_RESCUE <- TRUE

TME_MIN_MATCHED_GENES <- 1L
TME_MIN_COR_PAIRS <- 3L
TME_MIN_UNIQUE_VALUES <- 2L
TME_MIN_GROUP_N <- 2L
TME_MIN_TOTAL_SAMPLES <- 4L

IMMUNE_EXTRA_MIN_MATCHED_GENES <- 1L
IMMUNE_EXTRA_MIN_COR_PAIRS <- 3L
IMMUNE_EXTRA_MIN_GROUP_N <- 2L

CHECKPOINT_MIN_MATCHED_GENES <- 1L
CHECKPOINT_MIN_COR_PAIRS <- 3L
CHECKPOINT_MIN_GROUP_N <- 2L

INPUT_OVERRIDE_EXP_GCT_DIR <- NA_character_
INPUT_OVERRIDE_TME_GMT <- NA_character_
INPUT_OVERRIDE_IMMUNE_ACTIVITY_GMT <- NA_character_
INPUT_OVERRIDE_CHECKPOINT_GMT <- NA_character_
INPUT_OVERRIDE_TIDE_BATCH_FILE <- NA_character_
INPUT_OVERRIDE_TIDE_RESULT_DIR <- NA_character_

is_existing_file_v112 <- function(x) {
  length(x) == 1 && !is.na(x) && nzchar(x) && file.exists(x)
}
is_existing_dir_v112 <- function(x) {
  length(x) == 1 && !is.na(x) && nzchar(x) && dir.exists(x)
}

run_optional_module_v11 <- function(label, expr) {
  started <- Sys.time()
  status <- "SUCCESS"
  detail <- ""
  ans <- tryCatch(
    force(expr),
    error = function(e) {
      status <<- "ERROR"
      detail <<- conditionMessage(e)
      message("[", label, "] stopped after error: ", detail)
      invisible(NULL)
    }
  )

  log_root <- if (exists("MASTER_RESULT_DIR", inherits = TRUE)) {
    get("MASTER_RESULT_DIR", inherits = TRUE)
  } else {
    getwd()
  }
  log_dir <- file.path(log_root, "00_logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  status_file <- file.path(log_dir, "v11_2_optional_module_status.csv")
  row <- data.frame(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    module = label,
    status = status,
    detail = detail,
    elapsed_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    stringsAsFactors = FALSE
  )
  write.table(
    row, status_file, sep = ",", row.names = FALSE,
    col.names = !file.exists(status_file), append = file.exists(status_file),
    qmethod = "double"
  )
  ans
}

options(stringsAsFactors = FALSE)
set.seed(123)

cran_pkgs <- c("ggplot2", "pheatmap", "RColorBrewer", "survival", "survminer",
               "scales", "ggrepel", "patchwork", "dplyr", "tidyr", "readr", "stringr",
               "gridExtra", "data.table", "purrr", "forcats", "tibble")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(survival)
  library(survminer)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(gridExtra)
})

STEP_LOW_COLOR  <- "#00BFC4"   
STEP_HIGH_COLOR <- "#F8766D"   
STEP_GROUP_COLORS <- c("Low" = STEP_LOW_COLOR, "High" = STEP_HIGH_COLOR)
STEP_GROUP_COLORS_REVERSED <- c("High" = STEP_HIGH_COLOR, "Low" = STEP_LOW_COLOR)
STEP_LIGHT_BLUE_LOW <- "#F7FBFF"
STEP_LIGHT_BLUE_HIGH <- "#9ECAE1"

format_p_step <- function(p, prefix = "p = ") {
  if (length(p) == 0) return(character(0))
  vapply(p, function(pi) {
    if (is.na(pi) || !is.finite(pi)) return(paste0(prefix, "NA"))
    if (pi < 1e-4) return(paste0(prefix, formatC(pi, format = "e", digits = 2)))
    paste0(prefix, signif(pi, 3))
  }, character(1))
}

format_p_short_step <- function(p) {
  if (length(p) == 0) return(character(0))
  vapply(p, function(pi) {
    if (is.na(pi) || !is.finite(pi)) return("NA")
    if (pi < 1e-4) return(formatC(pi, format = "e", digits = 1))
    format(signif(pi, 2), trim = TRUE)
  }, character(1))
}

wrap_title_step <- function(x, width = 58) {
  stringr::str_wrap(as.character(x), width = width)
}

theme_step <- function(base_size = 10, title_size = 10, axis_title_size = 9,
                        axis_text_size = 8, strip_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5,
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = max(7, title_size - 1), hjust = 0.5,
                                   margin = margin(b = 4)),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size, color = "black"),
      strip.text = element_text(size = strip_size, face = "bold"),
      legend.title = element_text(size = max(7, axis_title_size - 1)),
      legend.text = element_text(size = max(7, axis_text_size - 1)),
      panel.grid.minor = element_blank()
    )
}

base_dir <- getwd()

exp_dir <- if (is_existing_dir_v112(INPUT_OVERRIDE_EXP_GCT_DIR)) INPUT_OVERRIDE_EXP_GCT_DIR else file.path(base_dir, "TCGA_201906cBioP", "exp_gct")
gmt_file <- file.path(base_dir, "gene set", "STEP_all_human_gene_sets_merged.gmt")

tme_gmt_candidates <- c(
  INPUT_OVERRIDE_TME_GMT,
  file.path(base_dir, "gene set", "TME_CIBERSORTConsensusMCP.gmt"),
  file.path(base_dir, "gene set", "TME_CIBERSORTConsensusMCP(1).gmt"),
  file.path(base_dir, "TME_CIBERSORTConsensusMCP.gmt"),
  file.path(base_dir, "TME_CIBERSORTConsensusMCP(1).gmt")
)
tme_gmt_file <- tme_gmt_candidates[file.exists(tme_gmt_candidates)][1]
purity_file <- file.path(base_dir, "TCGA_Tumor_purity_TCGAbiolinks.csv")

RUN_OFFICIAL_ESTIMATE <- TRUE
REQUIRE_OFFICIAL_ESTIMATE <- TRUE

STOP_IF_ESTIMATE_FAILS_PER_CANCER <- FALSE

out_dir <- file.path(base_dir, paste0("TCGA_STEP_full_results_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subdirs <- c(
  "00_logs",
  "01_scores_by_cancer",
  "02_gene_set_QC",
  "03_correlation_tables",
  "04_correlation_heatmaps",
  "05_scatter_plots",
  "06_purity_merged",
  "07_clinical_merged",
  "08_survival",
  "08b_survival_residual_CPE_AGE_models",
  "08c_survival_all_endpoints",
  "09_pan_cancer_summary",
  "10_figure_ready_CAF_adjusted",
  "11_CAF_residualization_details",
  "12_TME_scores",
  "13_TME_correlations",
  "14_TME_high_low_plots"
)
for (d in subdirs) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(
    parameter = c("TME_MIN_MATCHED_GENES", "TME_MIN_COR_PAIRS", "TME_MIN_UNIQUE_VALUES",
                  "TME_MIN_GROUP_N", "TME_MIN_TOTAL_SAMPLES"),
    value = c(TME_MIN_MATCHED_GENES, TME_MIN_COR_PAIRS, TME_MIN_UNIQUE_VALUES,
              TME_MIN_GROUP_N, TME_MIN_TOTAL_SAMPLES),
    interpretation = c(
      "Score a signature when at least this many genes are matched",
      "Minimum paired observations for Spearman correlation",
      "Minimum distinct values in each variable for correlation",
      "Minimum samples in each High/Low group for Wilcoxon testing",
      "Minimum samples required to retain a cancer expression matrix"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "12_TME_scores", "TME_minimum_threshold_settings_v11_2.csv"),
  row.names = FALSE
)

log_file <- file.path(out_dir, "00_logs", "run_log.txt")
cat("Run started:", as.character(Sys.time()), "\n", file = log_file)

msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

clean_na <- function(x) {
  x <- as.character(x)
  x[x %in% c("[Not Available]", "[Not Applicable]", "[Unknown]", "[Not Evaluated]",
             "", "NA", "N/A", "NULL", "null")] <- NA
  x
}

num_clean <- function(x) {
  x <- clean_na(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

normalize_tcga_id <- function(x) {
  
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  x <- gsub("_", "-", x)
  x <- trimws(x)
  x
}

get_sample15 <- function(x) {
  x <- normalize_tcga_id(x)
  substr(x, 1, 15)
}

get_patient12 <- function(x) {
  x <- normalize_tcga_id(x)
  substr(x, 1, 12)
}

read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  gs <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) >= 3) {
      nm <- parts[1]
      genes <- unique(toupper(trimws(parts[-c(1, 2)])))
      genes <- genes[!is.na(genes) & genes != ""]
      gs[[nm]] <- genes
    }
  }
  gs
}

add_filtered_step_gene_sets <- function(gene_sets, out_dir = NULL) {
  step_up_name <- "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"
  if (!(step_up_name %in% names(gene_sets))) {
    warning("Cannot find ", step_up_name, "; skip filtered STEP_UP construction.")
    return(gene_sets)
  }

  step_up <- unique(toupper(gene_sets[[step_up_name]]))

  conservative_blacklist <- unique(toupper(c(
    "COL1A1", "COL1A2", "COL3A1",
    "FN1", "TGFB2", "TGFBR2", "TGFBI",
    "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL5A3",
    "COL6A1", "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COL11A1", "COL12A1", "COL16A1",
    "DCN", "LUM", "SPARC", "THBS1", "THBS2", "POSTN", "BGN", "VCAN", "TNC",
    "FAP", "ACTA2", "TAGLN", "THY1", "PDGFRA", "PDGFRB",
    "MMP2", "MMP11", "MMP14", "LOX", "LOXL1", "LOXL2",
    "CTGF", "CCN2", "SERPINE1", "SERPINE2",
    "PLOD1", "PLOD2", "PLOD3", "P3H1", "PCOLCE", "PCOLCE2",
    "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FMOD", "MFAP5",
    "LRRC15", "MXRA5", "ECM1", "ECM2", "COMP", "CTHRC1"
  )))

  conservative_removed <- intersect(step_up, conservative_blacklist)
  conservative_kept <- setdiff(step_up, conservative_blacklist)

  tme_control_sets <- intersect(c("CAF_3gene_score", "CAF_4gene_index", "ECM_score", "Fibroblast_score"), names(gene_sets))
  tme_overlap_blacklist <- unique(toupper(unlist(gene_sets[tme_control_sets])))
  tme_removed <- intersect(step_up, tme_overlap_blacklist)
  tme_kept <- setdiff(step_up, tme_overlap_blacklist)

  gene_sets[["STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"]] <- conservative_kept
  gene_sets[["STEP_UP_TME_OVERLAP_FILTERED"]] <- tme_kept

  if (!is.null(out_dir)) {
    filter_dir <- file.path(out_dir, "02_gene_set_QC")
    dir.create(filter_dir, recursive = TRUE, showWarnings = FALSE)

    filter_summary <- data.frame(
      filtered_set = c("STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED", "STEP_UP_TME_OVERLAP_FILTERED"),
      original_set = step_up_name,
      original_n = length(step_up),
      removed_n = c(length(conservative_removed), length(tme_removed)),
      kept_n = c(length(conservative_kept), length(tme_kept)),
      removed_genes = c(paste(conservative_removed, collapse = ";"), paste(tme_removed, collapse = ";")),
      kept_genes = c(paste(conservative_kept, collapse = ";"), paste(tme_kept, collapse = ";")),
      filtering_basis = c(
        "Conservative published/canonical CAF, ECM, collagen, fibroblast activation, and matrix-remodeling blacklist",
        paste0("Overlap with project TME-control gene sets: ", paste(tme_control_sets, collapse = ";"))
      ),
      stringsAsFactors = FALSE
    )
    write.csv(filter_summary, file.path(filter_dir, "STEP_UP_filtered_gene_set_construction_summary.csv"), row.names = FALSE)

    make_filter_df <- function(filtered_set, action, genes) {
      if (length(genes) == 0) {
        return(data.frame(
          filtered_set = character(),
          action = character(),
          gene = character(),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        filtered_set = rep(filtered_set, length(genes)),
        action = rep(action, length(genes)),
        gene = genes,
        stringsAsFactors = FALSE
      )
    }

    filter_long <- rbind(
      make_filter_df("STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED", "removed", conservative_removed),
      make_filter_df("STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED", "kept", conservative_kept),
      make_filter_df("STEP_UP_TME_OVERLAP_FILTERED", "removed", tme_removed),
      make_filter_df("STEP_UP_TME_OVERLAP_FILTERED", "kept", tme_kept)
    )
    write.csv(filter_long, file.path(filter_dir, "STEP_UP_filtered_gene_level_details.csv"), row.names = FALSE)
  }

  gene_sets
}

detect_gct_header_skip <- function(file, max_lines = 80) {
  lines <- readLines(file, n = max_lines, warn = FALSE)
  for (i in seq_along(lines)) {
    fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    has_sample <- any(grepl("^TCGA[\\.-]", fields, ignore.case = TRUE))
    has_gene_col <- any(fields %in% c("Symbol", "Name", "NAME", "Gene", "gene",
                                      "Hugo_Symbol", "Hugo.Symbol", "Description"))
    if (length(fields) >= 3 && has_sample && has_gene_col) return(i - 1)
  }
  return(2)
}

write_gct_diagnostic <- function(file, cancer, out_dir) {
  diag_dir <- file.path(out_dir, "00_logs", "gct_diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  lines <- readLines(file, n = 30, warn = FALSE)
  writeLines(lines, con = file.path(diag_dir, paste0(cancer, "_first30_lines.txt")))
  diag <- data.frame()
  for (sk in c(0, 1, 2, 3, 4, 5, 10)) {
    tmp <- tryCatch(read.delim(file, skip = sk, nrows = 5, check.names = FALSE,
                               stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(tmp)) {
      diag <- rbind(diag, data.frame(
        cancer = cancer,
        skip = sk,
        ncol = ncol(tmp),
        colnames_preview = paste(head(colnames(tmp), 15), collapse = ";"),
        stringsAsFactors = FALSE
      ))
    }
  }
  write.csv(diag, file.path(diag_dir, paste0(cancer, "_skip_column_diagnostic.csv")), row.names = FALSE)
}

read_gct <- function(file, gene_sets_ref = NULL, cancer_name = "UNKNOWN", out_dir_for_diag = NULL) {
  skip_n <- detect_gct_header_skip(file)
  dat <- read.delim(file, skip = skip_n, check.names = FALSE, stringsAsFactors = FALSE)

  candidate_cols <- intersect(c("Symbol", "Name", "NAME", "Gene", "gene",
                                "Hugo_Symbol", "Hugo.Symbol", "Description",
                                "gene_name", "geneSymbol"), colnames(dat))
  if (length(candidate_cols) == 0) candidate_cols <- colnames(dat)[1]

  if (!is.null(gene_sets_ref)) {
    gene_universe <- unique(toupper(unlist(gene_sets_ref)))
    overlaps <- sapply(candidate_cols, function(cc) {
      vals <- toupper(trimws(as.character(dat[[cc]])))
      vals <- vals[!is.na(vals) & vals != ""]
      sum(unique(vals) %in% gene_universe)
    })
    gene_col <- candidate_cols[which.max(overlaps)]
  } else {
    gene_col <- candidate_cols[1]
    overlaps <- NA
  }

  genes_raw <- as.character(dat[[gene_col]])
  genes <- toupper(trimws(genes_raw))

  bad_gene <- is.na(genes) | genes == "" | genes == "NA" | genes == "?" | genes == "---" | genes == "NULL"
  if (any(bad_gene)) {
    dat <- dat[!bad_gene, , drop = FALSE]
    genes <- genes[!bad_gene]
  }

  non_expr <- intersect(c("Symbol", "Description", "Name", "NAME", "Gene", "gene",
                          "Hugo_Symbol", "Hugo.Symbol", "gene_name", "geneSymbol",
                          "Entrez_Gene_Id", "Entrez.Gene.Id"), colnames(dat))
  expr <- dat[, setdiff(colnames(dat), non_expr), drop = FALSE]

  tcga_cols <- grep("^TCGA[\\.-]", colnames(expr), value = TRUE, ignore.case = TRUE)
  if (length(tcga_cols) >= 3) expr <- expr[, tcga_cols, drop = FALSE]

  colnames(expr) <- normalize_tcga_id(colnames(expr))
  expr <- as.data.frame(lapply(expr, function(x) suppressWarnings(as.numeric(as.character(x)))))

  all_na <- apply(expr, 1, function(x) all(is.na(x)))
  if (any(all_na)) {
    expr <- expr[!all_na, , drop = FALSE]
    genes <- genes[!all_na]
  }

  if (nrow(expr) == 0 || length(genes) == 0) {
    if (!is.null(out_dir_for_diag)) write_gct_diagnostic(file, cancer_name, out_dir_for_diag)
    stop("No valid gene rows remain after filtering. Gene column used: ", gene_col,
         "; skip used: ", skip_n,
         "; candidate columns: ", paste(candidate_cols, collapse = ";"))
  }

  rownames(expr) <- make.unique(genes)
  expr <- as.matrix(expr)

  gene_symbol <- sub("\\.\\d+$", "", rownames(expr))
  expr_df <- data.frame(gene = gene_symbol, expr, check.names = FALSE)
  if (nrow(expr_df) == 0) {
    if (!is.null(out_dir_for_diag)) write_gct_diagnostic(file, cancer_name, out_dir_for_diag)
    stop("No rows to aggregate after gene-symbol processing. Gene column used: ", gene_col)
  }
  expr_agg <- aggregate(. ~ gene, data = expr_df, FUN = function(x) mean(x, na.rm = TRUE))
  expr_agg <- expr_agg[!is.na(expr_agg$gene) & expr_agg$gene != "", , drop = FALSE]
  rownames(expr_agg) <- expr_agg$gene
  expr_agg$gene <- NULL

  expr_mat <- as.matrix(expr_agg)
  mode(expr_mat) <- "numeric"
  attr(expr_mat, "gene_col_used") <- gene_col
  attr(expr_mat, "gene_col_overlap") <- if (length(candidate_cols) > 0 && !is.null(gene_sets_ref)) max(overlaps) else NA
  attr(expr_mat, "candidate_gene_cols") <- paste(candidate_cols, collapse = ";")
  attr(expr_mat, "gct_skip_used") <- skip_n
  attr(expr_mat, "n_bad_gene_removed") <- sum(bad_gene)
  attr(expr_mat, "n_all_na_removed") <- sum(all_na)
  expr_mat
}

auto_log2 <- function(expr) {
  q95 <- suppressWarnings(quantile(expr, 0.95, na.rm = TRUE))
  mx <- suppressWarnings(max(expr, na.rm = TRUE))
  if (is.finite(mx) && (mx > 100 || q95 > 50)) {
    expr2 <- log2(expr + 1)
    attr(expr2, "log2_transformed") <- TRUE
  } else {
    expr2 <- expr
    attr(expr2, "log2_transformed") <- FALSE
  }
  expr2
}

add_gene_alias_rows <- function(expr) {
  alias_pairs <- list(
    c("CYR61", "CCN1"),
    c("CTGF", "CCN2"),
    c("NOV", "CCN3"),
    c("WISP1", "CCN4"),
    c("WISP2", "CCN5"),
    c("WISP3", "CCN6"),
    c("SEPT11", "SEPTIN11"),
    c("FAM129B", "NIBAN2"),
    c("CASC5", "KNL1"),
    c("MARCHF5", "MARCH5")
  )

  for (pair in alias_pairs) {
    old <- pair[1]
    new <- pair[2]
    if (old %in% rownames(expr) && !(new %in% rownames(expr))) {
      add <- expr[old, , drop = FALSE]
      rownames(add) <- new
      expr <- rbind(expr, add)
    }
    if (new %in% rownames(expr) && !(old %in% rownames(expr))) {
      add <- expr[new, , drop = FALSE]
      rownames(add) <- old
      expr <- rbind(expr, add)
    }
  }
  expr
}

zscore_rows <- function(expr) {
  mu <- rowMeans(expr, na.rm = TRUE)
  sdv <- apply(expr, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- NA
  z <- sweep(expr, 1, mu, "-")
  z <- sweep(z, 1, sdv, "/")
  z[is.na(z)] <- 0
  z
}

score_gene_sets_zmean <- function(expr_log, gene_sets, min_matched_genes = 1L) {
  min_matched_genes <- max(1L, as.integer(min_matched_genes))
  expr_log <- add_gene_alias_rows(expr_log)
  z <- zscore_rows(expr_log)

  scores <- matrix(NA_real_, nrow = ncol(expr_log), ncol = length(gene_sets))
  rownames(scores) <- colnames(expr_log)
  colnames(scores) <- names(gene_sets)

  qc <- data.frame()

  for (nm in names(gene_sets)) {
    genes <- unique(toupper(gene_sets[[nm]]))
    matched <- intersect(genes, rownames(z))
    missing <- setdiff(genes, rownames(z))

    if (length(matched) >= min_matched_genes) {
      scores[, nm] <- colMeans(z[matched, , drop = FALSE], na.rm = TRUE)
    }

    qc <- rbind(qc, data.frame(
      gene_set = nm,
      n_genes = length(genes),
      n_matched = length(matched),
      n_missing = length(missing),
      matched_fraction = ifelse(length(genes) > 0, length(matched) / length(genes), NA_real_),
      min_matched_genes_required = min_matched_genes,
      score_included = length(matched) >= min_matched_genes,
      score_exclusion_reason = ifelse(length(matched) >= min_matched_genes,
                                      "included",
                                      "no_matched_gene"),
      matched_genes = paste(matched, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }

  list(scores = as.data.frame(scores), qc = qc)
}

make_feature_by_cancer_matrix <- function(cor_df, target_col, variable_order = NULL) {
  if (is.null(cor_df) || nrow(cor_df) == 0) return(NULL)
  dat <- cor_df[cor_df$target == target_col, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)
  dat <- dat[!duplicated(dat[, c("cancer", "variable")]), , drop = FALSE]
  wide <- reshape(dat[, c("variable", "cancer", "cor")],
                  idvar = "variable", timevar = "cancer", direction = "wide")
  rownames(wide) <- wide$variable
  wide$variable <- NULL
  colnames(wide) <- sub("^cor\\.", "", colnames(wide))
  mat <- as.matrix(wide)
  mode(mat) <- "numeric"
  if (!is.null(variable_order)) {
    keep <- intersect(variable_order, rownames(mat))
    extras <- setdiff(rownames(mat), keep)
    mat <- mat[c(keep, extras), , drop = FALSE]
  }
  mat
}

write_tme_cor_heatmap <- function(cor_df, target_col, variables, variable_order,
                                  out_prefix, title, width = 12, height = 9) {
  dat <- cor_df[cor_df$target == target_col & cor_df$variable %in% variables, , drop = FALSE]
  if (nrow(dat) == 0) return(invisible(FALSE))
  write.csv(dat, paste0(out_prefix, "_long.csv"), row.names = FALSE)
  mat <- make_feature_by_cancer_matrix(dat, target_col = target_col, variable_order = variable_order)
  if (is.null(mat)) return(invisible(FALSE))
  write.csv(mat, paste0(out_prefix, "_matrix.csv"))
  safe_cor_heatmap(
    mat,
    out_pdf = paste0(out_prefix, ".pdf"),
    out_png = paste0(out_prefix, ".png"),
    title = title,
    width = width,
    height = height,
    png_width = round(width * 300),
    png_height = round(height * 300)
  )
}

sem_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

make_high_low_stats <- function(df, score_col, features, cancer, min_n_per_group = TME_MIN_GROUP_N) {
  if (!(score_col %in% colnames(df))) return(data.frame())
  features <- unique(features[features %in% colnames(df)])
  if (length(features) == 0) return(data.frame())

  score <- df[[score_col]]
  ok_score <- is.finite(score)
  min_n_per_group <- max(2L, as.integer(min_n_per_group))
  if (sum(ok_score) < 2L * min_n_per_group || length(unique(score[ok_score])) < 2) return(data.frame())
  cutoff <- stats::median(score[ok_score], na.rm = TRUE)

  group <- rep(NA_character_, nrow(df))
  group[ok_score & score <= cutoff] <- "Low"
  group[ok_score & score > cutoff] <- "High"

  out <- data.frame()
  for (feat in features) {
    y <- df[[feat]]
    dat <- data.frame(value = y, group = group, stringsAsFactors = FALSE)
    dat <- dat[is.finite(dat$value) & !is.na(dat$group), , drop = FALSE]
    n_low <- sum(dat$group == "Low")
    n_high <- sum(dat$group == "High")
    if (n_low < min_n_per_group || n_high < min_n_per_group) next

    pval <- tryCatch(stats::wilcox.test(value ~ group, data = dat)$p.value, error = function(e) NA_real_)
    low_values <- dat$value[dat$group == "Low"]
    high_values <- dat$value[dat$group == "High"]
    out <- rbind(out, data.frame(
      cancer = cancer,
      score = score_col,
      feature = feat,
      cutoff = cutoff,
      n_low = n_low,
      n_high = n_high,
      mean_low = mean(low_values, na.rm = TRUE),
      mean_high = mean(high_values, na.rm = TRUE),
      median_low = stats::median(low_values, na.rm = TRUE),
      median_high = stats::median(high_values, na.rm = TRUE),
      sem_low = sem_safe(low_values),
      sem_high = sem_safe(high_values),
      mean_diff_high_minus_low = mean(high_values, na.rm = TRUE) - mean(low_values, na.rm = TRUE),
      median_diff_high_minus_low = stats::median(high_values, na.rm = TRUE) - stats::median(low_values, na.rm = TRUE),
      p = pval,
      low_n_warning = n_low < 5 || n_high < 5,
      stringsAsFactors = FALSE
    ))
  }
  if (nrow(out) > 0) {
    out$FDR <- p.adjust(out$p, method = "BH")
  }
  out
}

plot_high_low_bar_features <- function(stats_df, feature_order, out_prefix, title, width = 12, height = 8) {
  if (is.null(stats_df) || nrow(stats_df) == 0) return(invisible(FALSE))
  feature_order <- unique(feature_order[feature_order %in% stats_df$feature])
  if (length(feature_order) == 0) feature_order <- unique(stats_df$feature)
  stats_df$feature <- factor(stats_df$feature, levels = rev(feature_order))

  plot_df <- rbind(
    data.frame(feature = stats_df$feature, group = "Low",
               mean = stats_df$mean_low, sem = stats_df$sem_low, stringsAsFactors = FALSE),
    data.frame(feature = stats_df$feature, group = "High",
               mean = stats_df$mean_high, sem = stats_df$sem_high, stringsAsFactors = FALSE)
  )
  plot_df$group <- factor(plot_df$group, levels = c("Low", "High"))

  value_range <- range(c(plot_df$mean - plot_df$sem, plot_df$mean + plot_df$sem),
                       na.rm = TRUE)
  span <- diff(value_range)
  if (!is.finite(span) || span == 0) span <- max(abs(value_range), na.rm = TRUE)
  if (!is.finite(span) || span == 0) span <- 1

  ann <- stats_df
  ann$label_y <- pmax(ann$mean_low + ifelse(is.finite(ann$sem_low), ann$sem_low, 0),
                      ann$mean_high + ifelse(is.finite(ann$sem_high), ann$sem_high, 0),
                      na.rm = TRUE) + 0.05 * span
  ann$p_label <- vapply(ann$p, format_p_step, character(1), prefix = "p = ")

  p <- ggplot(plot_df, aes(x = feature, y = mean, fill = group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             color = "grey35", linewidth = 0.15) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  position = position_dodge(width = 0.75), width = 0.25,
                  linewidth = 0.25, na.rm = TRUE) +
    geom_text(data = ann, aes(x = feature, y = label_y, label = p_label),
              inherit.aes = FALSE, hjust = 0, size = 2.7) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(title = wrap_title_step(title, 70), x = NULL,
         y = "Mean Z-mean signature score",
         fill = "CAF-corrected STEP_UP") +
    theme_step(base_size = 10, title_size = 10, axis_title_size = 9,
                axis_text_size = 8, strip_size = 8) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      plot.margin = margin(6, 34, 6, 6)
    )

  ggsave(paste0(out_prefix, ".pdf"), p, width = width, height = height,
         limitsize = FALSE, useDingbats = FALSE)
  ggsave(paste0(out_prefix, ".png"), p, width = width, height = height,
         dpi = 300, limitsize = FALSE)
  print(p)
  invisible(TRUE)
}

plot_single_high_low_feature <- function(df, score_col, feature, cancer, out_prefix, title) {
  if (!(score_col %in% colnames(df)) || !(feature %in% colnames(df))) return(invisible(FALSE))
  score <- suppressWarnings(as.numeric(df[[score_col]]))
  ok_score <- is.finite(score)
  if (sum(ok_score) < 2L * TME_MIN_GROUP_N || length(unique(score[ok_score])) < 2) return(invisible(FALSE))
  cutoff <- stats::median(score[ok_score], na.rm = TRUE)
  dat <- data.frame(
    group = ifelse(ok_score & score > cutoff, "High",
                   ifelse(ok_score, "Low", NA_character_)),
    value = suppressWarnings(as.numeric(df[[feature]])),
    stringsAsFactors = FALSE
  )
  dat <- dat[is.finite(dat$value) & !is.na(dat$group), , drop = FALSE]
  if (sum(dat$group == "Low") < TME_MIN_GROUP_N || sum(dat$group == "High") < TME_MIN_GROUP_N) return(invisible(FALSE))
  dat$group <- factor(dat$group, levels = c("Low", "High"))
  pval <- tryCatch(stats::wilcox.test(value ~ group, data = dat, exact = FALSE)$p.value,
                   error = function(e) NA_real_)

  yr <- range(dat$value, na.rm = TRUE)
  span <- diff(yr)
  if (!is.finite(span) || span == 0) span <- 1
  y_ann <- yr[2] + 0.10 * span

  p <- ggplot(dat, aes(x = group, y = value, fill = group, color = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.82,
                 linewidth = 0.35) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.45,
                show.legend = FALSE) +
    annotate("text", x = 1.5, y = y_ann,
             label = format_p_step(pval, "Wilcoxon p = "), size = 3.4) +
    scale_fill_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
    scale_color_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    coord_cartesian(clip = "off") +
    labs(title = wrap_title_step(title, 52), x = NULL, y = feature,
         fill = "CAF-corrected STEP_UP") +
    theme_step(base_size = 10.5, title_size = 10, axis_title_size = 9,
                axis_text_size = 8.5) +
    theme(legend.position = "none",
          plot.margin = margin(6, 6, 12, 6))
  ggsave(paste0(out_prefix, ".pdf"), p, width = 4.5, height = 4.5,
         useDingbats = FALSE)
  ggsave(paste0(out_prefix, ".png"), p, width = 4.5, height = 4.5, dpi = 300)
  print(p)
  invisible(TRUE)
}

read_purity <- function(file) {
  pur <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(pur) <- make.names(colnames(pur))

  if (!("Sample.ID" %in% colnames(pur))) {
    stop("Purity file does not contain Sample.ID column.")
  }

  pur$Sample.ID.norm <- normalize_tcga_id(pur$Sample.ID)
  pur$sample15 <- get_sample15(pur$Sample.ID.norm)
  pur$patient_id <- get_patient12(pur$Sample.ID.norm)

  for (cc in intersect(c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"), colnames(pur))) {
    pur[[cc]] <- num_clean(pur[[cc]])
  }
  pur
}

dedup_purity_by_sample15 <- function(pur) {
  purity_cols <- intersect(c("CPE", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC"), colnames(pur))
  if (length(purity_cols) == 0) {
    pur <- pur[!duplicated(pur$sample15), , drop = FALSE]
    return(pur)
  }
  pur$non_missing_purity_n <- rowSums(!is.na(pur[, purity_cols, drop = FALSE]))
  pur <- pur[order(pur$sample15, -pur$non_missing_purity_n), , drop = FALSE]
  pur <- pur[!duplicated(pur$sample15), , drop = FALSE]
  pur$non_missing_purity_n <- NULL
  pur
}

dedup_clinical_by_patient <- function(clin_std) {
  clin_cols <- intersect(c("OS_MONTHS", "OS_EVENT", "DFS_MONTHS", "DFS_EVENT", "AGE", "STAGE"), colnames(clin_std))
  if (length(clin_cols) == 0) {
    clin_std <- clin_std[!duplicated(clin_std$patient_id), , drop = FALSE]
    return(clin_std)
  }
  clin_std$non_missing_clin_n <- rowSums(!is.na(clin_std[, clin_cols, drop = FALSE]))
  clin_std <- clin_std[order(clin_std$patient_id, -clin_std$non_missing_clin_n), , drop = FALSE]
  clin_std <- clin_std[!duplicated(clin_std$patient_id), , drop = FALSE]
  clin_std$non_missing_clin_n <- NULL
  clin_std
}

find_clinical_file <- function(cancer_lower, exp_dir) {
  files <- list.files(exp_dir, pattern = "\\.txt$", full.names = TRUE, ignore.case = TRUE)
  cand <- files[
    grepl(paste0("^", cancer_lower, "_"), basename(files), ignore.case = TRUE) &
      grepl("clinical", basename(files), ignore.case = TRUE)
  ]
  if (length(cand) == 0) return(NA_character_)
  cand[1]
}

read_cbio_clinical <- function(file) {
  lines <- readLines(file, warn = FALSE)
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  header_idx <- which(sapply(split_lines, function(x) any(x == "PATIENT_ID")))

  if (length(header_idx) > 0) {
    dat <- read.delim(file, skip = header_idx[1] - 1, check.names = FALSE,
                      stringsAsFactors = FALSE, comment.char = "")
  } else {
    dat <- read.delim(file, check.names = FALSE,
                      stringsAsFactors = FALSE, comment.char = "")
  }

  colnames(dat) <- make.names(colnames(dat), unique = TRUE)
  dat[] <- lapply(dat, clean_na)
  dat
}

standardize_clinical <- function(clin) {
  pick <- function(possible) {
    hit <- intersect(possible, colnames(clin))
    if (length(hit) > 0) hit[1] else NA_character_
  }

  patient_col <- pick(c("PATIENT_ID", "Patient.Identifier"))
  if (is.na(patient_col)) stop("No PATIENT_ID column found.")

  out <- data.frame(patient_id = get_patient12(clin[[patient_col]]), stringsAsFactors = FALSE)

  os_status_col <- pick(c("OS_STATUS", "Overall.Survival.Status"))
  os_month_col  <- pick(c("OS_MONTHS", "Overall.Survival..Months.", "Overall.Survival.Months"))
  dfs_status_col <- pick(c("DFS_STATUS", "Disease.Free.Status"))
  dfs_month_col  <- pick(c("DFS_MONTHS", "Disease.Free..Months.", "Disease.Free.Months"))
  age_col <- pick(c("AGE", "Diagnosis.Age", "Age"))
  sex_col <- pick(c("SEX", "Sex"))
  stage_col <- pick(c("AJCC_PATHOLOGIC_TUMOR_STAGE", "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code",
                      "CLINICAL_STAGE", "PATH_STAGE", "STAGE"))

  out$OS_STATUS_raw <- if (!is.na(os_status_col)) clin[[os_status_col]] else NA
  out$OS_MONTHS <- if (!is.na(os_month_col)) num_clean(clin[[os_month_col]]) else NA_real_
  out$DFS_STATUS_raw <- if (!is.na(dfs_status_col)) clin[[dfs_status_col]] else NA
  out$DFS_MONTHS <- if (!is.na(dfs_month_col)) num_clean(clin[[dfs_month_col]]) else NA_real_
  out$AGE <- if (!is.na(age_col)) num_clean(clin[[age_col]]) else NA_real_
  out$SEX <- if (!is.na(sex_col)) clin[[sex_col]] else NA
  out$STAGE <- if (!is.na(stage_col)) clin[[stage_col]] else NA

  out$OS_EVENT <- ifelse(grepl("DECEASED|Dead|1:", out$OS_STATUS_raw, ignore.case = TRUE), 1,
                         ifelse(grepl("LIVING|Alive|0:", out$OS_STATUS_raw, ignore.case = TRUE), 0, NA))

  out$DFS_EVENT <- ifelse(grepl("Recurred|Progressed|1:", out$DFS_STATUS_raw, ignore.case = TRUE), 1,
                          ifelse(grepl("DiseaseFree|Disease Free|0:", out$DFS_STATUS_raw, ignore.case = TRUE), 0, NA))

  out
}

cor_one_vs_all <- function(df, target_cols, method = "spearman",
                           min_pairs = TME_MIN_COR_PAIRS,
                           min_unique = TME_MIN_UNIQUE_VALUES) {
  min_pairs <- max(3L, as.integer(min_pairs))
  min_unique <- max(2L, as.integer(min_unique))
  numeric_cols <- colnames(df)[sapply(df, is.numeric)]
  res <- data.frame()

  for (target in target_cols) {
    if (!(target %in% numeric_cols)) next
    for (v in numeric_cols) {
      if (v == target) next
      x <- suppressWarnings(as.numeric(df[[target]]))
      y <- suppressWarnings(as.numeric(df[[v]]))
      ok <- is.finite(x) & is.finite(y)
      n_ok <- sum(ok)
      if (n_ok < min_pairs ||
          length(unique(x[ok])) < min_unique ||
          length(unique(y[ok])) < min_unique) next

      ct <- tryCatch(
        suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE)),
        error = function(e) NULL
      )
      if (is.null(ct)) next
      res <- rbind(res, data.frame(
        target = target,
        variable = v,
        n = n_ok,
        cor = unname(ct$estimate),
        p = ct$p.value,
        method = method,
        low_n_warning = n_ok < 10,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(res) > 0) res$FDR <- p.adjust(res$p, method = "BH")
  res
}

zscore_vec <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

safe_residualize_one <- function(df, y_col, caf_col, out_col, cancer = NA_character_) {
  qc_empty <- data.frame(
    cancer = cancer,
    score = y_col,
    caf_covariate = caf_col,
    residual_score = out_col,
    n = NA_integer_,
    alpha_intercept = NA_real_,
    beta_CAF = NA_real_,
    p_CAF = NA_real_,
    r2 = NA_real_,
    formula = paste0(y_col, " = alpha + beta * ", caf_col, " + residual"),
    status = "not_run",
    stringsAsFactors = FALSE
  )

  if (!(y_col %in% colnames(df)) || !(caf_col %in% colnames(df))) {
    qc_empty$status <- "missing_score_or_caf"
    df[[paste0(out_col, "_predicted_by_", caf_col)]] <- NA_real_
    df[[paste0(out_col, "_predicted_by_", caf_col)]] <- NA_real_
    df[[out_col]] <- NA_real_
    df[[paste0(out_col, "_z")]] <- NA_real_
    return(list(df = df, qc = qc_empty))
  }

  ok <- complete.cases(df[[y_col]], df[[caf_col]])
  if (sum(ok) < 10 || sd(df[[y_col]][ok], na.rm = TRUE) == 0 || sd(df[[caf_col]][ok], na.rm = TRUE) == 0) {
    qc_empty$n <- sum(ok)
    qc_empty$status <- "too_few_samples_or_zero_variance"
    df[[out_col]] <- NA_real_
    df[[paste0(out_col, "_z")]] <- NA_real_
    return(list(df = df, qc = qc_empty))
  }

  dat <- data.frame(
    y = as.numeric(df[[y_col]][ok]),
    caf = as.numeric(df[[caf_col]][ok])
  )
  fit <- lm(y ~ caf, data = dat)
  sm <- summary(fit)

  pred_vec <- rep(NA_real_, nrow(df))
  resid_vec <- rep(NA_real_, nrow(df))
  pred_vec[ok] <- fitted(fit)
  resid_vec[ok] <- residuals(fit)

  pred_col <- paste0(out_col, "_predicted_by_", caf_col)
  df[[pred_col]] <- pred_vec
  df[[out_col]] <- resid_vec
  df[[paste0(out_col, "_z")]] <- zscore_vec(resid_vec)

  qc <- data.frame(
    cancer = cancer,
    score = y_col,
    caf_covariate = caf_col,
    residual_score = out_col,
    n = sum(ok),
    alpha_intercept = unname(coef(fit)["(Intercept)"]),
    beta_CAF = unname(coef(fit)["caf"]),
    p_CAF = sm$coefficients["caf", "Pr(>|t|)"],
    r2 = sm$r.squared,
    formula = paste0(y_col, " = ", signif(unname(coef(fit)["(Intercept)"]), 5),
                     " + ", signif(unname(coef(fit)["caf"]), 5), " * ", caf_col,
                     " + residual"),
    status = "ok",
    stringsAsFactors = FALSE
  )
  list(df = df, qc = qc)
}

add_caf_adjusted_step_scores <- function(df, cancer = NA_character_) {
  up_col <- "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"
  down_col <- "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN"

  qc_all <- data.frame()

  if (all(c(up_col, down_col) %in% colnames(df))) {
    df$STEP_UP_MINUS_DOWN_RAW <- as.numeric(df[[up_col]]) - as.numeric(df[[down_col]])
    df$STEP_UP_MINUS_DOWN_RAW_z <- zscore_vec(df$STEP_UP_MINUS_DOWN_RAW)
  } else {
    df$STEP_UP_MINUS_DOWN_RAW <- NA_real_
    df$STEP_UP_MINUS_DOWN_RAW_z <- NA_real_
  }

  for (one in list(
    c(up_col, "CAF_3gene_score", "STEP_UP_residual_CAF3"),
    c(down_col, "CAF_3gene_score", "STEP_DOWN_residual_CAF3"),
    c("STEP_UP_MINUS_DOWN_RAW", "CAF_3gene_score", "STEP_BALANCE_residual_CAF3"),
    c(up_col, "CAF_4gene_index", "STEP_UP_residual_CAF4"),
    c(down_col, "CAF_4gene_index", "STEP_DOWN_residual_CAF4"),
    c("STEP_UP_MINUS_DOWN_RAW", "CAF_4gene_index", "STEP_BALANCE_residual_CAF4")
  )) {
    rr <- safe_residualize_one(
      df = df,
      y_col = one[1],
      caf_col = one[2],
      out_col = one[3],
      cancer = cancer
    )
    df <- rr$df
    qc_all <- rbind(qc_all, rr$qc)
  }

  list(df = df, qc = qc_all)
}

save_gg_both <- function(p, file_base, width = 6.5, height = 5.2) {
  print(p)
  ggsave(paste0(file_base, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(file_base, ".png"), p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

plot_scatter <- function(df, x, y, cancer, outbase) {
  if (!(x %in% colnames(df)) || !(y %in% colnames(df))) return(NULL)
  ok <- complete.cases(df[[x]], df[[y]])
  if (sum(ok) < 5) return(NULL)

  ct <- suppressWarnings(cor.test(df[[x]][ok], df[[y]][ok], method = "spearman", exact = FALSE))
  subtitle <- paste0("Spearman rho = ", round(unname(ct$estimate), 3),
                     ", p = ", signif(ct$p.value, 3), ", n = ", sum(ok))

  p <- ggplot(df[ok, ], aes_string(x = x, y = y)) +
    geom_point(size = 1.8, alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(title = paste0(cancer, ": ", y, " vs ", x),
         subtitle = subtitle, x = x, y = y)

  save_gg_both(p, outbase, width = 6.3, height = 5.2)
}

save_ggsurvplot_both <- function(g, file_base, width = 6.8, height = 6.3) {
  grDevices::pdf(paste0(file_base, ".pdf"), width = width, height = height,
                 useDingbats = FALSE)
  print(g, newpage = FALSE)
  grDevices::dev.off()

  grDevices::png(paste0(file_base, ".png"), width = round(width * 300),
                 height = round(height * 300), res = 300)
  print(g, newpage = FALSE)
  grDevices::dev.off()

  print(g)
  invisible(TRUE)
}

compact_score_label_step <- function(score_col) {
  x <- score_col
  x <- gsub("_HUMAN_TO_HUMAN", "", x)
  x <- gsub("STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", "STEP_UP", x)
  x <- gsub("STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN", "STEP_DOWN", x)
  x <- gsub("_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", "_UP", x)
  x <- gsub("_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN", "_DOWN", x)
  x
}

run_survival <- function(df, score_col, cancer, out_prefix) {
  needed <- c(score_col, "OS_MONTHS", "OS_EVENT")
  if (!all(needed %in% colnames(df))) return(NULL)

  surv_df <- df
  surv_df$OS_MONTHS <- suppressWarnings(as.numeric(surv_df$OS_MONTHS))
  surv_df$OS_EVENT <- normalize_event_vector(surv_df$OS_EVENT)
  surv_df[[score_col]] <- suppressWarnings(as.numeric(surv_df[[score_col]]))
  surv_df <- surv_df[complete.cases(surv_df[, needed]), , drop = FALSE]
  surv_df <- surv_df[is.finite(surv_df$OS_MONTHS) & surv_df$OS_MONTHS > 0 &
                       surv_df$OS_EVENT %in% c(0, 1) &
                       is.finite(surv_df[[score_col]]), , drop = FALSE]
  if (nrow(surv_df) < 20 || length(unique(surv_df$OS_EVENT)) < 2) return(NULL)

  med <- median(surv_df[[score_col]], na.rm = TRUE)
  surv_df$group <- ifelse(surv_df[[score_col]] >= med, "High", "Low")
  surv_df$group <- factor(surv_df$group, levels = c("Low", "High"))
  if (length(unique(surv_df$group)) < 2) return(NULL)

  fit <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ group, data = surv_df)
  lr <- survdiff(Surv(OS_MONTHS, OS_EVENT) ~ group, data = surv_df)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  p_label <- format_p_step(lr_p, "log-rank p = ")
  score_short <- compact_score_label_step(score_col)

  g <- survminer::ggsurvplot(
    fit,
    data = surv_df,
    risk.table = "absolute",
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.col = "strata",
    pval = p_label,
    pval.size = 3.5,
    conf.int = FALSE,
    censor = TRUE,
    censor.shape = 124,
    censor.size = 3,
    palette = unname(STEP_GROUP_COLORS[c("Low", "High")]),
    legend.title = "CAF3-corrected STEP_UP",
    legend.labs = c("Low", "High"),
    xlab = "Months",
    ylab = "Overall survival probability",
    title = wrap_title_step(paste0(cancer, " OS by ", score_short), 55),
    font.title = c(10, "bold"),
    font.x = c(9),
    font.y = c(9),
    font.tickslab = c(8),
    ggtheme = theme_step(base_size = 10, title_size = 10,
                          axis_title_size = 9, axis_text_size = 8)
  )
  g$plot <- g$plot +
    theme(legend.position = "top",
          plot.margin = margin(5.5, 8, 5.5, 5.5))
  g$table <- g$table +
    theme_step(base_size = 8, title_size = 8,
                axis_title_size = 8, axis_text_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, face = "plain", hjust = 0))
  save_ggsurvplot_both(g, paste0(out_prefix, "_KM"), width = 6.8, height = 6.2)

  surv_df$score_scaled <- as.numeric(scale(surv_df[[score_col]]))
  cox_uni <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ score_scaled, data = surv_df)
  su <- summary(cox_uni)

  out <- data.frame(
    cancer = cancer,
    score = score_col,
    model = "univariate_score",
    n = nrow(surv_df),
    logrank_p = lr_p,
    HR = su$coefficients[1, "exp(coef)"],
    p = su$coefficients[1, "Pr(>|z|)"],
    lower95 = su$conf.int[1, "lower .95"],
    upper95 = su$conf.int[1, "upper .95"],
    stringsAsFactors = FALSE
  )

  if (all(c("CPE", "CAF_3gene_score", "AGE") %in% colnames(surv_df))) {
    mv_df <- surv_df
    mv_df$CPE_scaled <- as.numeric(scale(suppressWarnings(as.numeric(mv_df$CPE))))
    mv_df$CAF_scaled <- as.numeric(scale(suppressWarnings(as.numeric(mv_df$CAF_3gene_score))))
    mv_df$AGE <- suppressWarnings(as.numeric(mv_df$AGE))
    mv_df <- mv_df[complete.cases(mv_df[, c("OS_MONTHS", "OS_EVENT",
                                             "score_scaled", "CPE_scaled",
                                             "CAF_scaled", "AGE")]), , drop = FALSE]

    if (nrow(mv_df) >= 20 && length(unique(mv_df$OS_EVENT)) >= 2) {
      cox_mv <- coxph(Surv(OS_MONTHS, OS_EVENT) ~
                        score_scaled + CPE_scaled + CAF_scaled + AGE,
                      data = mv_df)
      sm <- summary(cox_mv)
      mv_out <- data.frame(
        cancer = cancer,
        score = score_col,
        model = rownames(sm$coefficients),
        n = nrow(mv_df),
        logrank_p = NA_real_,
        HR = sm$coefficients[, "exp(coef)"],
        p = sm$coefficients[, "Pr(>|z|)"],
        lower95 = sm$conf.int[, "lower .95"],
        upper95 = sm$conf.int[, "upper .95"],
        stringsAsFactors = FALSE
      )
      out <- rbind(out, mv_out)
    }
  }

  write.csv(out, paste0(out_prefix, "_cox_summary.csv"), row.names = FALSE)
  out
}

run_survival_residual_cpe_age <- function(df, score_col, cancer, out_prefix) {

  needed <- c(score_col, "OS_MONTHS", "OS_EVENT")
  if (!all(needed %in% colnames(df))) return(NULL)

  surv_df <- df[complete.cases(df[, needed]), , drop = FALSE]
  if (nrow(surv_df) < 20 || length(unique(surv_df$OS_EVENT)) < 2) return(NULL)

  med <- median(surv_df[[score_col]], na.rm = TRUE)
  surv_df$group <- ifelse(surv_df[[score_col]] >= med, "High", "Low")
  surv_df$group <- factor(surv_df$group, levels = c("Low", "High"))

  lr_p <- NA_real_
  if (length(unique(surv_df$group)) == 2) {
    lr <- survdiff(Surv(OS_MONTHS, OS_EVENT) ~ group, data = surv_df)
    lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  }

  surv_df$score_scaled <- as.numeric(scale(surv_df[[score_col]]))

  cox_uni <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ score_scaled, data = surv_df)
  su <- summary(cox_uni)

  out <- data.frame(
    cancer = cancer,
    score = score_col,
    model_set = "residual_CPE_AGE",
    model = "univariate_score",
    covariates = "score_only",
    n = nrow(surv_df),
    logrank_p = lr_p,
    HR = su$coefficients[1, "exp(coef)"],
    p = su$coefficients[1, "Pr(>|z|)"],
    lower95 = su$conf.int[1, "lower .95"],
    upper95 = su$conf.int[1, "upper .95"],
    stringsAsFactors = FALSE
  )

  if (all(c("CPE", "AGE") %in% colnames(surv_df))) {
    mv_df <- surv_df
    mv_df$CPE_scaled <- as.numeric(scale(mv_df$CPE))
    mv_df <- mv_df[complete.cases(mv_df[, c("OS_MONTHS", "OS_EVENT", "score_scaled", "CPE_scaled", "AGE")]), ]

    if (nrow(mv_df) >= 20 && length(unique(mv_df$OS_EVENT)) >= 2) {
      cox_mv <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ score_scaled + CPE_scaled + AGE, data = mv_df)
      sm <- summary(cox_mv)
      mv_out <- data.frame(
        cancer = cancer,
        score = score_col,
        model_set = "residual_CPE_AGE",
        model = rownames(sm$coefficients),
        covariates = "score + CPE + AGE",
        n = nrow(mv_df),
        logrank_p = NA_real_,
        HR = sm$coefficients[, "exp(coef)"],
        p = sm$coefficients[, "Pr(>|z|)"],
        lower95 = sm$conf.int[, "lower .95"],
        upper95 = sm$conf.int[, "upper .95"],
        stringsAsFactors = FALSE
      )
      out <- rbind(out, mv_out)
    }
  }

  write.csv(out, paste0(out_prefix, "_cox_residual_CPE_AGE_summary.csv"), row.names = FALSE)
  out
}

export_caf_residualization_sample_details <- function(df, residual_qc, cancer, out_dir) {

  if (is.null(residual_qc) || nrow(residual_qc) == 0) return(NULL)

  all_detail <- data.frame()

  get_col_or_na <- function(dat, col) {
    if (!is.na(col) && col %in% colnames(dat)) dat[[col]] else rep(NA_real_, nrow(dat))
  }

  get_chr_col_or_na <- function(dat, col) {
    if (!is.na(col) && col %in% colnames(dat)) as.character(dat[[col]]) else rep(NA_character_, nrow(dat))
  }

  sample_raw <- get_chr_col_or_na(df, "sample_raw")
  sample <- get_chr_col_or_na(df, "sample")
  sample15 <- get_chr_col_or_na(df, "sample15")
  patient_id <- get_chr_col_or_na(df, "patient_id")

  for (i in seq_len(nrow(residual_qc))) {
    yy <- as.character(residual_qc$score[i])
    cc <- as.character(residual_qc$caf_covariate[i])
    rr <- as.character(residual_qc$residual_score[i])
    pred_col <- paste0(rr, "_predicted_by_", cc)
    rz <- paste0(rr, "_z")

    tmp <- data.frame(
      cancer = cancer,
      sample_raw = sample_raw,
      sample = sample,
      sample15 = sample15,
      patient_id = patient_id,

      residual_score = rr,
      original_score_name = yy,
      caf_covariate = cc,

      original_score_value = get_col_or_na(df, yy),
      caf_score_value = get_col_or_na(df, cc),
      predicted_by_CAF = get_col_or_na(df, pred_col),
      residual_value = get_col_or_na(df, rr),
      residual_z = get_col_or_na(df, rz),

      alpha_intercept = residual_qc$alpha_intercept[i],
      beta_CAF = residual_qc$beta_CAF[i],
      model_formula = residual_qc$formula[i],
      model_r2 = residual_qc$r2[i],
      model_p_CAF = residual_qc$p_CAF[i],
      model_status = residual_qc$status[i],

      stringsAsFactors = FALSE
    )

    all_detail <- rbind(all_detail, tmp)
  }

  if (nrow(all_detail) > 0) {
    write.csv(all_detail,
              file.path(out_dir, "11_CAF_residualization_details",
                        paste0(cancer, "_CAF_residualization_sample_details.csv")),
              row.names = FALSE)
  }
  invisible(all_detail)
}

plot_all_cancer_endpoint_forest <- function(cox_df, out_dir, endpoint = "OS",
                                            model_type = c("adjusted", "unadjusted")) {
  model_type <- match.arg(model_type)
  if (is.null(cox_df) || nrow(cox_df) == 0) return(NULL)

  show_scores <- c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_UP_residual_CAF3",
    "STEP_UP_residual_CAF4"
  )

  model_to_use <- ifelse(model_type == "adjusted", "score_scaled", "univariate_score")
  model_label <- ifelse(model_type == "adjusted", "Adjusted: score + CPE purity + age", "Unadjusted: score only")
  file_tag <- ifelse(model_type == "adjusted", "CPE_AGE_adjusted", "unadjusted")

  dat <- cox_df[
    cox_df$endpoint == endpoint &
      cox_df$model == model_to_use &
      cox_df$score %in% show_scores,
    , drop = FALSE
  ]
  if (nrow(dat) == 0) return(NULL)

  dat$score_label <- dat$score
  dat$score_label[dat$score == "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"] <- "Before CAF adjustment\nRaw STEP_UP"
  dat$score_label[dat$score == "STEP_UP_residual_CAF3"] <- "After CAF3 adjustment\nSTEP_UP residual"
  dat$score_label[dat$score == "STEP_UP_residual_CAF4"] <- "After CAF4 adjustment\nSTEP_UP residual"

  dat$neglog10p <- -log10(dat$p)
  dat$significance <- ifelse(!is.na(dat$FDR_by_endpoint_score) & dat$FDR_by_endpoint_score < 0.05,
                             "FDR < 0.05",
                             ifelse(dat$p < 0.05, "p < 0.05", "NS"))

  raw <- dat[dat$score == "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", c("cancer", "HR")]
  raw <- raw[order(raw$HR, decreasing = TRUE), , drop = FALSE]
  cancer_levels <- unique(c(raw$cancer, sort(unique(dat$cancer))))
  dat$cancer <- factor(dat$cancer, levels = rev(cancer_levels))

  write.csv(dat,
            file.path(out_dir, "10_figure_ready_CAF_adjusted",
                      paste0("FigureReady_", endpoint, "_all_cancers_RAW_vs_CAF_adjusted_forest_data_", file_tag, ".csv")),
            row.names = FALSE)

  h <- max(7, 0.23 * length(unique(dat$cancer)) + 2.2)

  p <- ggplot(dat, aes(x = HR, y = cancer)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.16, linewidth = 0.35) +
    geom_point(aes(shape = significance), size = 1.9) +
    facet_wrap(~ score_label, nrow = 1) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 8),
      legend.position = "bottom",
      plot.title = element_text(size = 11)
    ) +
    labs(title = paste0(endpoint, " forest plot across all cancers: before vs after CAF adjustment"),
         subtitle = model_label,
         x = "Hazard ratio per 1 SD increase",
         y = NULL,
         shape = "Significance")

  ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted",
                   paste0("FigureReady_", endpoint, "_all_cancers_RAW_vs_CAF_adjusted_forest_", file_tag, ".pdf")),
         p, width = 12.5, height = h)
  ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted",
                   paste0("FigureReady_", endpoint, "_all_cancers_RAW_vs_CAF_adjusted_forest_", file_tag, ".png")),
         p, width = 12.5, height = h, dpi = 300)
  print(p)

  for (sc in unique(dat$score_label)) {
    dd <- dat[dat$score_label == sc, , drop = FALSE]
    if (nrow(dd) == 0) next
    hh <- max(7, 0.25 * nrow(dd) + 2)
    p2 <- ggplot(dd, aes(x = HR, y = cancer)) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
      geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.16, linewidth = 0.35) +
      geom_point(aes(shape = significance), size = 2) +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 7),
            legend.position = "bottom",
            plot.title = element_text(size = 11)) +
      labs(title = paste0(endpoint, " forest plot: ", gsub("\\n", " ", sc)),
           subtitle = model_label,
           x = "Hazard ratio per 1 SD increase",
           y = NULL,
           shape = "Significance")

    clean_sc <- gsub("[^A-Za-z0-9]+", "_", sc)
    ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted",
                     paste0("FigureReady_", endpoint, "_all_cancers_forest_", clean_sc, "_", file_tag, ".pdf")),
           p2, width = 6.5, height = hh)
    ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted",
                     paste0("FigureReady_", endpoint, "_all_cancers_forest_", clean_sc, "_", file_tag, ".png")),
           p2, width = 6.5, height = hh, dpi = 300)
    print(p2)
  }

  invisible(dat)
}

detect_survival_endpoints <- function(df) {
  cn <- colnames(df)
  endpoints <- list()

  add_if_present <- function(name, time_candidates, event_candidates) {
    time_col <- time_candidates[time_candidates %in% cn][1]
    event_col <- event_candidates[event_candidates %in% cn][1]
    if (!is.na(time_col) && !is.na(event_col)) {
      endpoints[[name]] <<- list(time = time_col, event = event_col)
    }
  }

  add_if_present("OS",
                 c("OS_MONTHS", "OS_MONTH", "OS.time", "OS_TIME", "OS_DAYS"),
                 c("OS_EVENT", "OS_STATUS", "OS", "OS.event"))
  add_if_present("DFS",
                 c("DFS_MONTHS", "DFS_MONTH", "DFS.time", "DFS_TIME", "DFS_DAYS"),
                 c("DFS_EVENT", "DFS_STATUS", "DFS", "DFS.event"))
  add_if_present("PFS",
                 c("PFS_MONTHS", "PFS_MONTH", "PFS.time", "PFS_TIME", "PFS_DAYS"),
                 c("PFS_EVENT", "PFS_STATUS", "PFS", "PFS.event"))
  add_if_present("RFS",
                 c("RFS_MONTHS", "RFS_MONTH", "RFS.time", "RFS_TIME", "RFS_DAYS"),
                 c("RFS_EVENT", "RFS_STATUS", "RFS", "RFS.event"))

  endpoints
}

normalize_event_vector <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  xx <- toupper(trimws(as.character(x)))
  out <- rep(NA_real_, length(xx))
  out[xx %in% c("1", "TRUE", "DECEASED", "DEAD", "EVENT", "YES",
                "PROGRESSED", "RECURRENCE", "RECURRED", "RELAPSED")] <- 1
  out[xx %in% c("0", "FALSE", "LIVING", "ALIVE", "CENSORED", "NO",
                "DISEASEFREE", "DISEASE_FREE", "NO_RECURRENCE")] <- 0
  out[grepl("^1", xx)] <- 1
  out[grepl("^0", xx)] <- 0
  out
}

run_survival_endpoint <- function(df, score_col, cancer, endpoint_name, time_col, event_col,
                                  out_prefix, covariate_mode = c("score_CPE_AGE")) {
  covariate_mode <- match.arg(covariate_mode)
  needed <- c(score_col, time_col, event_col)
  if (!all(needed %in% colnames(df))) return(NULL)

  surv_df <- df
  surv_df$.time <- suppressWarnings(as.numeric(surv_df[[time_col]]))
  surv_df$.event <- normalize_event_vector(surv_df[[event_col]])
  surv_df$.score <- suppressWarnings(as.numeric(surv_df[[score_col]]))
  surv_df <- surv_df[complete.cases(surv_df[, c(".time", ".event", ".score")]), , drop = FALSE]
  surv_df <- surv_df[is.finite(surv_df$.time) & surv_df$.time > 0 &
                       surv_df$.event %in% c(0, 1) &
                       is.finite(surv_df$.score), , drop = FALSE]
  if (nrow(surv_df) < 20 || length(unique(surv_df$.event)) < 2) return(NULL)

  med <- median(surv_df$.score, na.rm = TRUE)
  surv_df$group <- ifelse(surv_df$.score >= med, "High", "Low")
  surv_df$group <- factor(surv_df$group, levels = c("Low", "High"))
  if (length(unique(surv_df$group)) < 2) return(NULL)

  fit <- survfit(Surv(.time, .event) ~ group, data = surv_df)
  lr <- survdiff(Surv(.time, .event) ~ group, data = surv_df)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  p_label <- format_p_step(lr_p, "log-rank p = ")
  score_short <- compact_score_label_step(score_col)

  g <- survminer::ggsurvplot(
    fit,
    data = surv_df,
    risk.table = "absolute",
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.col = "strata",
    pval = p_label,
    pval.size = 3.5,
    conf.int = FALSE,
    censor = TRUE,
    censor.shape = 124,
    censor.size = 3,
    palette = unname(STEP_GROUP_COLORS[c("Low", "High")]),
    legend.title = "CAF3-corrected STEP_UP",
    legend.labs = c("Low", "High"),
    xlab = "Months",
    ylab = paste0(endpoint_name, " probability"),
    title = wrap_title_step(paste0(cancer, " ", endpoint_name, " by ", score_short), 55),
    font.title = c(10, "bold"),
    font.x = c(9),
    font.y = c(9),
    font.tickslab = c(8),
    ggtheme = theme_step(base_size = 10, title_size = 10,
                          axis_title_size = 9, axis_text_size = 8)
  )
  g$plot <- g$plot + theme(legend.position = "top")
  g$table <- g$table +
    theme_step(base_size = 8, title_size = 8,
                axis_title_size = 8, axis_text_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, face = "plain", hjust = 0))
  save_ggsurvplot_both(g, paste0(out_prefix, "_KM"), width = 6.8, height = 6.2)

  surv_df$score_scaled <- as.numeric(scale(surv_df$.score))
  cox_uni <- coxph(Surv(.time, .event) ~ score_scaled, data = surv_df)
  su <- summary(cox_uni)

  out <- data.frame(
    cancer = cancer, endpoint = endpoint_name, score = score_col,
    model_set = covariate_mode, model = "univariate_score",
    covariates = "score_only", time_col = time_col, event_col = event_col,
    n = nrow(surv_df), n_events = sum(surv_df$.event == 1, na.rm = TRUE),
    logrank_p = lr_p,
    HR = su$coefficients[1, "exp(coef)"],
    p = su$coefficients[1, "Pr(>|z|)"],
    lower95 = su$conf.int[1, "lower .95"],
    upper95 = su$conf.int[1, "upper .95"],
    stringsAsFactors = FALSE
  )

  if (all(c("CPE", "AGE") %in% colnames(surv_df))) {
    mv_df <- surv_df
    mv_df$CPE_scaled <- as.numeric(scale(suppressWarnings(as.numeric(mv_df$CPE))))
    mv_df$AGE <- suppressWarnings(as.numeric(mv_df$AGE))
    mv_df <- mv_df[complete.cases(mv_df[, c(".time", ".event",
                                             "score_scaled", "CPE_scaled", "AGE")]), , drop = FALSE]
    if (nrow(mv_df) >= 20 && length(unique(mv_df$.event)) >= 2) {
      cox_mv <- coxph(Surv(.time, .event) ~ score_scaled + CPE_scaled + AGE,
                      data = mv_df)
      sm <- summary(cox_mv)
      mv_out <- data.frame(
        cancer = cancer, endpoint = endpoint_name, score = score_col,
        model_set = covariate_mode, model = rownames(sm$coefficients),
        covariates = "score + CPE + AGE",
        time_col = time_col, event_col = event_col,
        n = nrow(mv_df), n_events = sum(mv_df$.event == 1, na.rm = TRUE),
        logrank_p = NA_real_,
        HR = sm$coefficients[, "exp(coef)"],
        p = sm$coefficients[, "Pr(>|z|)"],
        lower95 = sm$conf.int[, "lower .95"],
        upper95 = sm$conf.int[, "upper .95"],
        stringsAsFactors = FALSE
      )
      out <- rbind(out, mv_out)
    }
  }

  write.csv(out, paste0(out_prefix, "_cox_summary.csv"), row.names = FALSE)
  out
}

plot_caf_adjustment_validation <- function(cor_long, out_dir) {
  if (is.null(cor_long) || nrow(cor_long) == 0) return(NULL)
  dat <- cor_long[cor_long$target %in% c(
      "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
      "STEP_UP_residual_CAF3",
      "STEP_UP_residual_CAF4"
    ) &
    cor_long$variable %in% c("CAF_3gene_score", "CAF_4gene_index",
                             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                             "HUMAN_MESENCHYMAL_SCORE_18",
                             "CPE", "ESTIMATE_ImmuneScore"), , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  dat$target_label <- dat$target
  dat$target_label[dat$target == "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"] <- "Raw STEP_UP"
  dat$target_label[dat$target == "STEP_UP_residual_CAF3"] <- "CAF3-adjusted STEP_UP"
  dat$target_label[dat$target == "STEP_UP_residual_CAF4"] <- "CAF4-adjusted STEP_UP"

  dat$variable_label <- dat$variable
  dat$variable_label[dat$variable == "CAF_3gene_score"] <- "CAF3"
  dat$variable_label[dat$variable == "CAF_4gene_index"] <- "CAF4"
  dat$variable_label[dat$variable == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "Hallmark EMT"
  dat$variable_label[dat$variable == "HUMAN_MESENCHYMAL_SCORE_18"] <- "M score"
  dat$variable_label[dat$variable == "CPE"] <- "CPE purity"
  dat$variable_label[dat$variable == "ESTIMATE_ImmuneScore"] <- "ImmuneScore"

  summary_dat <- aggregate(cor ~ target_label + variable_label, data = dat, FUN = function(x) median(x, na.rm = TRUE))
  colnames(summary_dat)[3] <- "median_cor"
  write.csv(summary_dat,
            file.path(out_dir, "10_figure_ready_CAF_adjusted", "FigureReady_CAF_adjustment_median_correlations.csv"),
            row.names = FALSE)

  ann <- dat %>%
    group_by(variable_label) %>%
    summarise(
      p = tryCatch(kruskal.test(cor ~ target_label)$p.value,
                   error = function(e) NA_real_),
      ymax = max(cor, na.rm = TRUE),
      ymin = min(cor, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(x = 2,
           y = ymax + 0.12 * ifelse(ymax > ymin, ymax - ymin, 1),
           label = vapply(p, format_p_step, character(1),
                          prefix = "Kruskal-Wallis p = "))

  target_colors <- c("Raw STEP_UP" = "#BDBDBD",
                     "CAF3-adjusted STEP_UP" = STEP_HIGH_COLOR,
                     "CAF4-adjusted STEP_UP" = STEP_LOW_COLOR)

  p <- ggplot(dat, aes(x = target_label, y = cor,
                       fill = target_label, color = target_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.5, width = 0.65, alpha = 0.82) +
    geom_text(data = ann, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, color = "black", size = 2.6) +
    scale_fill_manual(values = target_colors, drop = FALSE) +
    scale_color_manual(values = target_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    facet_wrap(~ variable_label, nrow = 2, scales = "free_y") +
    theme_step(base_size = 10, title_size = 10,
                axis_title_size = 9, axis_text_size = 8,
                strip_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          legend.position = "none") +
    labs(title = "CAF adjustment validation and retained biological associations",
         x = NULL, y = "Spearman correlation across cancers")

  ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted", "FigureReady_CAF_adjustment_validation_boxplot.pdf"),
         p, width = 10, height = 6)
  ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted", "FigureReady_CAF_adjustment_validation_boxplot.png"),
         p, width = 10, height = 6, dpi = 300)
  print(p)
}

plot_residual_cox_forest <- function(cox_df, out_dir, endpoint = "OS") {
  if (is.null(cox_df) || nrow(cox_df) == 0) return(NULL)
  show_scores <- c("STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
                   "STEP_UP_residual_CAF3", "STEP_UP_residual_CAF4")
  show_cancers <- c("KIRC", "HNSC", "LUAD", "UCEC", "LIHC", "LGG", "SKCM")
  dat <- cox_df[cox_df$endpoint == endpoint & cox_df$model == "score_scaled" &
                  cox_df$score %in% show_scores & cox_df$cancer %in% show_cancers, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  dat$score_label <- dat$score
  dat$score_label[dat$score == "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"] <- "Raw STEP_UP"
  dat$score_label[dat$score == "STEP_UP_residual_CAF3"] <- "CAF3-adjusted"
  dat$score_label[dat$score == "STEP_UP_residual_CAF4"] <- "CAF4-adjusted"

  write.csv(dat,
            file.path(out_dir, "10_figure_ready_CAF_adjusted", paste0("FigureReady_", endpoint, "_CAF_adjusted_Cox_forest_data.csv")),
            row.names = FALSE)

  p <- ggplot(dat, aes(x = HR, y = cancer)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.18) +
    geom_point(size = 2) +
    facet_wrap(~ score_label, nrow = 1) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 9),
          plot.title = element_text(size = 11)) +
    labs(title = paste0(endpoint, " Cox forest plot after CAF adjustment"),
         x = "Hazard ratio per 1 SD increase", y = NULL)

  ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted", paste0("FigureReady_", endpoint, "_CAF_adjusted_Cox_forest.pdf")),
         p, width = 9.5, height = 4.5)
  ggsave(file.path(out_dir, "10_figure_ready_CAF_adjusted", paste0("FigureReady_", endpoint, "_CAF_adjusted_Cox_forest.png")),
         p, width = 9.5, height = 4.5, dpi = 300)
  print(p)
}

safe_cor_heatmap <- function(mat, out_pdf, out_png, title = "Spearman correlation",
                             width = 8.5, height = 7.5,
                             png_width = 2550, png_height = 2250) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"

  keep_r <- apply(mat, 1, function(x) any(is.finite(x)))
  keep_c <- apply(mat, 2, function(x) any(is.finite(x)))
  mat <- mat[keep_r, keep_c, drop = FALSE]

  if (nrow(mat) < 2 || ncol(mat) < 2) {
    message("Skip heatmap: too few finite rows/cols for ", title)
    return(invisible(FALSE))
  }

  mat[!is.finite(mat)] <- 0

  row_unique <- nrow(unique(as.data.frame(mat))) > 1
  col_unique <- nrow(unique(as.data.frame(t(mat)))) > 1

  fs_col <- ifelse(ncol(mat) > 40, 4.2, ifelse(ncol(mat) > 25, 5.2, ifelse(ncol(mat) > 15, 6.2, 8)))
  fs_row <- ifelse(nrow(mat) > 35, 5.5, ifelse(nrow(mat) > 25, 6.5, 8))
  fs_num <- ifelse(ncol(mat) > 40, 3.5, ifelse(ncol(mat) > 25, 4.5, ifelse(ncol(mat) > 15, 5.5, 7)))

  pdf(out_pdf, width = width, height = height)
  pheatmap(mat,
           cluster_rows = row_unique,
           cluster_cols = col_unique,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_col = fs_col,
           fontsize_row = fs_row,
           fontsize_number = fs_num,
           angle_col = 90,
           main = title)
  dev.off()

  png(out_png, width = png_width, height = png_height, res = 300)
  pheatmap(mat,
           cluster_rows = row_unique,
           cluster_cols = col_unique,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_col = fs_col,
           fontsize_row = fs_row,
           fontsize_number = fs_num,
           angle_col = 90,
           main = title)
  dev.off()

  pheatmap(mat,
           cluster_rows = row_unique,
           cluster_cols = col_unique,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_col = fs_col,
           fontsize_row = fs_row,
           fontsize_number = fs_num,
           angle_col = 90,
           main = title)

  invisible(TRUE)
}

estimate_package_available <- function() {
  requireNamespace("estimate", quietly = TRUE)
}

check_required_estimate_or_stop <- function() {
  if (!RUN_OFFICIAL_ESTIMATE) return(invisible(TRUE))
  if (!estimate_package_available()) {
    stop(
      "Package 'estimate' is not installed, but this v7 script requires it.\n",
      "Please run first:\n",
      "install.packages('estimate', repos='http://R-Forge.R-project.org', dependencies=TRUE)\n",
      "Then test with: library(estimate)"
    )
  }
  suppressPackageStartupMessages(library(estimate))
  required_funs <- c("filterCommonGenes", "estimateScore")
  missing_funs <- required_funs[!sapply(required_funs, function(f) exists(f, where = asNamespace("estimate"), inherits = FALSE))]
  if (length(missing_funs) > 0) {
    stop("Package 'estimate' is installed, but missing required functions: ", paste(missing_funs, collapse = ", "))
  }
  invisible(TRUE)
}

write_estimate_input_gct <- function(expr_log, file, format = c("gct", "plain")) {
  format <- match.arg(format)
  mat <- expr_log
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]

  original_samples <- colnames(mat)
  safe_samples <- paste0("S", seq_len(ncol(mat)))
  colnames(mat) <- safe_samples

  map <- data.frame(
    safe_sample = safe_samples,
    original_sample = original_samples,
    sample15 = get_sample15(original_samples),
    stringsAsFactors = FALSE
  )

  mat_df <- as.data.frame(mat, check.names = FALSE)
  mat_df <- data.frame(NAME = rownames(mat_df), Description = rownames(mat_df), mat_df, check.names = FALSE)

  if (format == "gct") {
    con <- file(file, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines("#1.2", con)
    writeLines(paste(nrow(mat_df), ncol(mat_df) - 2, sep = "\t"), con)
    write.table(mat_df, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    write.table(mat_df, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }

  map
}

read_estimate_score_file <- function(score_file, sample_map, cancer) {
  
  sc <- tryCatch(read.delim(score_file, skip = 2, check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(sc) || nrow(sc) == 0 || ncol(sc) < 3) {
    sc <- read.delim(score_file, skip = 0, check.names = FALSE, stringsAsFactors = FALSE)
  }

  name_col <- intersect(c("NAME", "Name", "name"), colnames(sc))[1]
  if (is.na(name_col)) name_col <- colnames(sc)[1]

  row_names <- as.character(sc[[name_col]])
  sample_cols <- setdiff(colnames(sc), c("NAME", "Name", "name", "Description", "description"))
  mat <- sc[, sample_cols, drop = FALSE]
  mat <- as.data.frame(lapply(mat, num_clean), check.names = FALSE)
  rownames(mat) <- row_names

  wanted <- intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"), rownames(mat))
  if (length(wanted) == 0) {
    stop("Official ESTIMATE score file for ", cancer,
         " does not contain expected rows. Observed rows: ",
         paste(head(rownames(mat), 20), collapse = ", "))
  }

  out <- as.data.frame(t(as.matrix(mat[wanted, , drop = FALSE])), stringsAsFactors = FALSE)
  colnames(out) <- paste0("ESTIMATE_", colnames(out))
  out$safe_sample <- rownames(out)
  out <- merge(out, sample_map, by = "safe_sample", all.x = TRUE)
  out$sample <- normalize_tcga_id(out$original_sample)
  out$sample15 <- get_sample15(out$sample)
  out <- out[!duplicated(out$sample15), , drop = FALSE]
  out
}

calculate_official_estimate_scores <- function(expr_log, cancer, out_dir) {
  if (!RUN_OFFICIAL_ESTIMATE) return(NULL)
  check_required_estimate_or_stop()

  est_dir <- file.path(out_dir, "10_official_ESTIMATE_scores")
  dir.create(est_dir, recursive = TRUE, showWarnings = FALSE)

  formats_to_try <- c("gct", "plain")
  errors <- c()

  for (fmt in formats_to_try) {
    input_file <- file.path(est_dir, paste0(cancer, "_estimate_input_", fmt, ifelse(fmt == "gct", ".gct", ".txt")))
    filtered_file <- file.path(est_dir, paste0(cancer, "_estimate_filtered_", fmt, ".gct"))
    score_file <- file.path(est_dir, paste0(cancer, "_estimate_scores_", fmt, ".gct"))
    map_file <- file.path(est_dir, paste0(cancer, "_estimate_sample_map_", fmt, ".csv"))

    sample_map <- write_estimate_input_gct(expr_log, input_file, format = fmt)
    write.csv(sample_map, map_file, row.names = FALSE)

    ok <- tryCatch({
      estimate::filterCommonGenes(input.f = input_file, output.f = filtered_file, id = "GeneSymbol")
      estimate::estimateScore(filtered_file, score_file, platform = "illumina")
      TRUE
    }, error = function(e) {
      errors <<- c(errors, paste0("[", fmt, "] ", conditionMessage(e)))
      FALSE
    })

    if (isTRUE(ok) && file.exists(score_file)) {
      out <- tryCatch(read_estimate_score_file(score_file, sample_map, cancer),
                      error = function(e) {
                        errors <<- c(errors, paste0("[", fmt, " read_score] ", conditionMessage(e)))
                        NULL
                      })
      if (!is.null(out)) {
        msg("Official ESTIMATE succeeded for ", cancer, " using input format: ", fmt)
        return(out)
      }
    }
  }

  stop("Official ESTIMATE failed for ", cancer, " using all input formats. Errors: ",
       paste(errors, collapse = " | "))
}

if (!file.exists(gmt_file)) stop("Cannot find GMT file: ", gmt_file)
if (!file.exists(purity_file)) stop("Cannot find purity file: ", purity_file)
if (!dir.exists(exp_dir)) stop("Cannot find expression folder: ", exp_dir)

check_required_estimate_or_stop()
msg("Official ESTIMATE package is installed and required functions are available.")
msg("v7.4 mode: full 01-09 summary outputs will be generated; official ESTIMATE will be merged when successful.")

gene_sets <- read_gmt(gmt_file)
gene_sets <- add_filtered_step_gene_sets(gene_sets, out_dir = out_dir)

if (length(tme_gmt_file) == 0 || is.na(tme_gmt_file) || !file.exists(tme_gmt_file)) {
  msg("TME GMT not found. TME/immune cell module will be skipped. Expected one of: ", paste(tme_gmt_candidates, collapse = " | "))
  tme_gene_sets <- list()
} else {
  tme_gene_sets <- read_gmt(tme_gmt_file)
  msg("TME GMT loaded: ", tme_gmt_file, "; TME signatures=", length(tme_gene_sets))
}

tme_signature_names <- names(tme_gene_sets)
epithelial_compare_vars <- intersect(c("HUMAN_EPITHELIAL_SCORE_14", "KOHN_EMT_EPITHELIAL"), names(gene_sets))
tme_module_all_vars <- c(tme_signature_names, epithelial_compare_vars)
tme_module_core_vars <- c(
  "CD8 T cells",
  "T cells_M",
  "Cytotoxic lymphocytes_M",
  "Cytotoxic lymphocytes_co",
  "NK cells_M",
  "NK activated",
  "B lineage_M",
  "Plasma cells",
  "Tregs",
  "Tfh",
  "Monocytic lineage_M",
  "Monocytes",
  "Macrophage M0",
  "Macrophage M1",
  "Macrophage M2",
  "Macrophages_co",
  "Dendritic activated",
  "mDendritic cells_M",
  "Neutrophils_M",
  "Neutrophils",
  "Endothelial cells_M",
  "Fibroblasts_M",
  epithelial_compare_vars
)

cd8_tme_signature_names <- grep("CD8", tme_signature_names, value = TRUE, ignore.case = TRUE)
tme_module_core_vars <- unique(c(cd8_tme_signature_names, tme_module_core_vars))
tme_module_core_vars <- unique(tme_module_core_vars[tme_module_core_vars %in% tme_module_all_vars])

purity <- read_purity(purity_file)
purity_raw_n <- nrow(purity)
purity <- dedup_purity_by_sample15(purity)
msg("Purity deduplicated by sample15: raw=", purity_raw_n, "; unique=", nrow(purity))

gct_files <- list.files(exp_dir, pattern = "_tcga_expression\\.gct$", full.names = TRUE, ignore.case = TRUE)

msg("Gene sets loaded: ", length(gene_sets))
msg("TME gene sets loaded: ", length(tme_gene_sets))
msg("Purity rows loaded: ", nrow(purity))
msg("Expression GCT files found: ", length(gct_files))

if (length(gct_files) == 0) stop("No *_tcga_expression.gct files found.")

all_qc <- list()
all_cor <- list()
all_surv <- list()
all_surv_residual_cpe_age <- list()
all_surv_all_endpoints <- list()
all_merged <- list()
all_residual_qc <- list()
all_tme_qc <- list()
all_tme_highlow <- list()
skipped_cancers <- data.frame(cancer = character(), file = character(), reason = character(), stringsAsFactors = FALSE)
estimate_failed_cancers <- data.frame(cancer = character(), reason = character(), stringsAsFactors = FALSE)

for (gct in gct_files) {

  cancer <- sub("_tcga_expression\\.gct$", "", basename(gct), ignore.case = TRUE)
  cancer_upper <- toupper(cancer)
  cancer_lower <- tolower(cancer)

  msg("\n==============================")
  msg("Processing cancer: ", cancer_upper)
  msg("==============================")

  expr_raw <- tryCatch(
    read_gct(gct, gene_sets_ref = gene_sets, cancer_name = cancer_upper, out_dir_for_diag = out_dir),
    error = function(e) {
      msg("SKIP ", cancer_upper, ": failed to read/aggregate GCT: ", conditionMessage(e))
      skipped_cancers <<- rbind(skipped_cancers, data.frame(
        cancer = cancer_upper,
        file = basename(gct),
        reason = paste0("read_gct failed: ", conditionMessage(e)),
        stringsAsFactors = FALSE
      ))
      return(NULL)
    }
  )
  if (is.null(expr_raw)) next

  n_bad_gene_removed <- attr(expr_raw, "n_bad_gene_removed")
  n_all_na_removed <- attr(expr_raw, "n_all_na_removed")
  gene_col_used <- attr(expr_raw, "gene_col_used")
  gene_col_overlap <- attr(expr_raw, "gene_col_overlap")
  candidate_gene_cols <- attr(expr_raw, "candidate_gene_cols")
  gct_skip_used <- attr(expr_raw, "gct_skip_used")

  if (nrow(expr_raw) < 50 || ncol(expr_raw) < TME_MIN_TOTAL_SAMPLES) {
    reason <- paste0("too few genes/samples after reading: genes=", nrow(expr_raw), ", samples=", ncol(expr_raw))
    msg("SKIP ", cancer_upper, ": ", reason)
    skipped_cancers <<- rbind(skipped_cancers, data.frame(
      cancer = cancer_upper,
      file = basename(gct),
      reason = reason,
      stringsAsFactors = FALSE
    ))
    next
  }

  expr <- auto_log2(expr_raw)

  msg("Genes: ", nrow(expr), "; samples: ", ncol(expr),
      "; log2 transformed: ", attr(expr, "log2_transformed"),
      "; gene column used: ", gene_col_used,
      "; GMT overlap: ", gene_col_overlap,
      "; candidate gene columns: ", candidate_gene_cols,
      "; gct skip used: ", gct_skip_used,
      "; empty gene rows removed: ", n_bad_gene_removed,
      "; all-NA rows removed: ", n_all_na_removed)

  scored <- score_gene_sets_zmean(expr, gene_sets)
  scores <- scored$scores
  scores$sample_raw <- rownames(scores)
  scores$sample <- normalize_tcga_id(scores$sample_raw)
  scores$sample15 <- get_sample15(scores$sample)
  scores$patient_id <- get_patient12(scores$sample)
  scores$cancer_from_file <- cancer_upper

  scores <- scores[, c("sample_raw", "sample", "sample15", "patient_id", "cancer_from_file",
                       setdiff(colnames(scores), c("sample_raw", "sample", "sample15", "patient_id", "cancer_from_file")))]

  qc <- scored$qc
  qc$cancer <- cancer_upper
  all_qc[[cancer_upper]] <- qc

  write.csv(scores, file.path(out_dir, "01_scores_by_cancer", paste0(cancer_upper, "_gene_set_scores_zmean.csv")), row.names = FALSE)
  write.csv(qc, file.path(out_dir, "02_gene_set_QC", paste0(cancer_upper, "_gene_set_matched_missing.csv")), row.names = FALSE)

  if (length(tme_gene_sets) > 0) {
    tme_scored <- score_gene_sets_zmean(expr, tme_gene_sets, min_matched_genes = TME_MIN_MATCHED_GENES)
    tme_scores <- tme_scored$scores
    tme_scores$sample_raw <- rownames(tme_scores)
    tme_scores$sample <- normalize_tcga_id(tme_scores$sample_raw)
    tme_scores$sample15 <- get_sample15(tme_scores$sample)
    tme_scores <- tme_scores[, c("sample_raw", "sample", "sample15", setdiff(colnames(tme_scores), c("sample_raw", "sample", "sample15")))]

    tme_qc <- tme_scored$qc
    tme_score_diagnostics <- do.call(rbind, lapply(intersect(names(tme_gene_sets), colnames(tme_scores)), function(sig) {
      vv <- suppressWarnings(as.numeric(tme_scores[[sig]]))
      vv <- vv[is.finite(vv)]
      n_non_na <- length(vv)
      n_unique <- length(unique(vv))
      data.frame(
        gene_set = sig,
        n_score_non_na = n_non_na,
        n_score_unique = n_unique,
        score_sd = if (n_non_na >= 2) stats::sd(vv) else NA_real_,
        correlation_eligible_minimum = n_non_na >= TME_MIN_COR_PAIRS && n_unique >= TME_MIN_UNIQUE_VALUES,
        high_low_eligible_minimum = n_non_na >= 2L * TME_MIN_GROUP_N && n_unique >= 2,
        stringsAsFactors = FALSE
      )
    }))
    if (!is.null(tme_score_diagnostics) && nrow(tme_score_diagnostics) > 0) {
      tme_qc <- merge(tme_qc, tme_score_diagnostics, by = "gene_set", all.x = TRUE, sort = FALSE)
    }
    tme_qc$analysis_status <- ifelse(
      tme_qc$n_matched < TME_MIN_MATCHED_GENES,
      "excluded_no_matched_gene",
      ifelse(is.na(tme_qc$n_score_non_na) | tme_qc$n_score_non_na == 0,
             "excluded_all_scores_NA",
             ifelse(tme_qc$n_score_unique < 2,
                    "score_available_but_no_variation",
                    "score_available"))
    )
    tme_qc$cancer <- cancer_upper
    all_tme_qc[[cancer_upper]] <- tme_qc

    write.csv(tme_scores,
              file.path(out_dir, "12_TME_scores", paste0(cancer_upper, "_TME_signature_scores_zmean.csv")),
              row.names = FALSE)
    write.csv(tme_qc,
              file.path(out_dir, "12_TME_scores", paste0(cancer_upper, "_TME_signature_matched_missing.csv")),
              row.names = FALSE)

    tme_keep_cols <- c("sample15", intersect(tme_signature_names, colnames(tme_scores)))
    scores <- merge(scores, tme_scores[, tme_keep_cols, drop = FALSE], by = "sample15", all.x = TRUE)
    msg("TME signatures scored and merged for ", cancer_upper, ": ", length(intersect(tme_signature_names, colnames(scores))), " signatures")
  }

  pur_cols <- intersect(c("Sample.ID", "sample15", "Cancer.type", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"), colnames(purity))
  merged <- merge(scores, purity[, pur_cols, drop = FALSE], by = "sample15", all.x = TRUE)

  if (RUN_OFFICIAL_ESTIMATE) {
    est_scores <- tryCatch(
      calculate_official_estimate_scores(expr, cancer_upper, out_dir),
      error = function(e) {
        msg("Official ESTIMATE failed for ", cancer_upper, ": ", conditionMessage(e))
        estimate_failed_cancers <<- rbind(estimate_failed_cancers, data.frame(
          cancer = cancer_upper,
          reason = conditionMessage(e),
          stringsAsFactors = FALSE
        ))
        if (isTRUE(STOP_IF_ESTIMATE_FAILS_PER_CANCER)) stop(e)
        NULL
      }
    )
    if (!is.null(est_scores)) {
      merged <- merge(merged, est_scores, by = "sample15", all.x = TRUE)
      msg("Official ESTIMATE merged for ", cancer_upper,
          ": StromalScore=", if ("ESTIMATE_StromalScore" %in% colnames(merged)) sum(!is.na(merged$ESTIMATE_StromalScore)) else 0,
          "; ImmuneScore=", if ("ESTIMATE_ImmuneScore" %in% colnames(merged)) sum(!is.na(merged$ESTIMATE_ImmuneScore)) else 0,
          "; ESTIMATEScore=", if ("ESTIMATE_ESTIMATEScore" %in% colnames(merged)) sum(!is.na(merged$ESTIMATE_ESTIMATEScore)) else 0,
          "; TumorPurity=", if ("ESTIMATE_TumorPurity" %in% colnames(merged)) sum(!is.na(merged$ESTIMATE_TumorPurity)) else 0)
    }
  }

  write.csv(merged, file.path(out_dir, "06_purity_merged", paste0(cancer_upper, "_scores_plus_purity.csv")), row.names = FALSE)

  clin_file <- find_clinical_file(cancer_lower, exp_dir)
  if (!is.na(clin_file)) {
    clin <- read_cbio_clinical(clin_file)
    clin_std <- standardize_clinical(clin)
    clin_std <- dedup_clinical_by_patient(clin_std)
    merged2 <- merge(merged, clin_std, by = "patient_id", all.x = TRUE)
  } else {
    msg("Clinical file not found for ", cancer_upper)
    merged2 <- merged
  }

  residual_res <- add_caf_adjusted_step_scores(merged2, cancer = cancer_upper)
  merged2 <- residual_res$df
  residual_qc <- residual_res$qc
  all_residual_qc[[cancer_upper]] <- residual_qc
  write.csv(residual_qc,
            file.path(out_dir, "02_gene_set_QC", paste0(cancer_upper, "_CAF_adjusted_residual_score_QC.csv")),
            row.names = FALSE)

  export_caf_residualization_sample_details(merged2, residual_qc, cancer_upper, out_dir)

  match_summary <- data.frame(
    cancer = cancer_upper,
    n_expression_samples = nrow(scores),
    n_purity_CPE_matched = if ("CPE" %in% colnames(merged2)) sum(!is.na(merged2$CPE)) else NA,
    n_purity_ESTIMATE_matched = if ("ESTIMATE" %in% colnames(merged2)) sum(!is.na(merged2$ESTIMATE)) else NA,
    n_purity_ABSOLUTE_matched = if ("ABSOLUTE" %in% colnames(merged2)) sum(!is.na(merged2$ABSOLUTE)) else NA,
    n_OS_months_matched = if ("OS_MONTHS" %in% colnames(merged2)) sum(!is.na(merged2$OS_MONTHS)) else NA,
    n_OS_event_matched = if ("OS_EVENT" %in% colnames(merged2)) sum(!is.na(merged2$OS_EVENT)) else NA,
    n_AGE_matched = if ("AGE" %in% colnames(merged2)) sum(!is.na(merged2$AGE)) else NA,
    n_ESTIMATE_StromalScore_matched = if ("ESTIMATE_StromalScore" %in% colnames(merged2)) sum(!is.na(merged2$ESTIMATE_StromalScore)) else NA,
    n_ESTIMATE_ImmuneScore_matched = if ("ESTIMATE_ImmuneScore" %in% colnames(merged2)) sum(!is.na(merged2$ESTIMATE_ImmuneScore)) else NA,
    n_ESTIMATE_ESTIMATEScore_matched = if ("ESTIMATE_ESTIMATEScore" %in% colnames(merged2)) sum(!is.na(merged2$ESTIMATE_ESTIMATEScore)) else NA,
    gene_col_used = gene_col_used,
    gene_col_overlap = gene_col_overlap,
    gct_skip_used = gct_skip_used,
    stringsAsFactors = FALSE
  )
  write.csv(match_summary, file.path(out_dir, "00_logs", paste0(cancer_upper, "_sample_matching_summary.csv")), row.names = FALSE)
  msg("Matching summary for ", cancer_upper,
      ": CPE=", match_summary$n_purity_CPE_matched,
      "; ESTIMATE purity=", match_summary$n_purity_ESTIMATE_matched,
      "; OS_MONTHS=", match_summary$n_OS_months_matched,
      "; OS_EVENT=", match_summary$n_OS_event_matched,
      "; ESTIMATE_StromalScore=", match_summary$n_ESTIMATE_StromalScore_matched,
      "; ESTIMATE_ImmuneScore=", match_summary$n_ESTIMATE_ImmuneScore_matched)

  write.csv(merged2, file.path(out_dir, "07_clinical_merged", paste0(cancer_upper, "_scores_purity_clinical.csv")), row.names = FALSE)

  highlow_targets <- intersect(c("STEP_UP_residual_CAF3", "STEP_UP_residual_CAF4"), colnames(merged2))
  for (score_col in highlow_targets) {
    score_safe <- gsub("[^A-Za-z0-9]+", "_", score_col)

    all_features <- unique(tme_module_all_vars[tme_module_all_vars %in% colnames(merged2)])
    if (length(all_features) > 0) {
      stats_all <- make_high_low_stats(merged2, score_col, all_features, cancer_upper)
      if (nrow(stats_all) > 0) {
        all_tme_highlow[[paste(cancer_upper, score_col, "ALL", sep = "_")]] <- stats_all
        all_prefix <- file.path(out_dir, "14_TME_high_low_plots", paste0(cancer_upper, "_", score_safe, "_ALL_TME_plus_epithelial_high_low_bar"))
        write.csv(stats_all, paste0(all_prefix, "_stats.csv"), row.names = FALSE)
        plot_high_low_bar_features(
          stats_all, feature_order = all_features, out_prefix = all_prefix,
          title = paste0(cancer_upper, ": ", score_col, " high vs low - all TME + epithelial"),
          width = 10, height = max(7, 0.25 * length(unique(stats_all$feature)) + 2)
        )
      }
    }

    core_features <- unique(tme_module_core_vars[tme_module_core_vars %in% colnames(merged2)])
    if (length(core_features) > 0) {
      stats_core <- make_high_low_stats(merged2, score_col, core_features, cancer_upper)
      if (nrow(stats_core) > 0) {
        all_tme_highlow[[paste(cancer_upper, score_col, "CORE", sep = "_")]] <- stats_core
        core_prefix <- file.path(out_dir, "14_TME_high_low_plots", paste0(cancer_upper, "_", score_safe, "_CORE_TME_plus_epithelial_high_low_bar"))
        write.csv(stats_core, paste0(core_prefix, "_stats.csv"), row.names = FALSE)
        plot_high_low_bar_features(
          stats_core, feature_order = core_features, out_prefix = core_prefix,
          title = paste0(cancer_upper, ": ", score_col, " high vs low - core TME + epithelial"),
          width = 9, height = max(6, 0.28 * length(unique(stats_core$feature)) + 2)
        )
      }
    }

    estimate_features <- unique(intersect(c(
      "ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore", "ESTIMATE_ESTIMATEScore", "ESTIMATE_TumorPurity",
      "HUMAN_EPITHELIAL_SCORE_14", "KOHN_EMT_EPITHELIAL",
      "HUMAN_MESENCHYMAL_SCORE_18", "KOHN_EMT_MESENCHYMAL",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
      "CAF_3gene_score", "CAF_4gene_index", "ECM_score"
    ), colnames(merged2)))
    if (length(estimate_features) > 0) {
      stats_est <- make_high_low_stats(merged2, score_col, estimate_features, cancer_upper)
      if (nrow(stats_est) > 0) {
        all_tme_highlow[[paste(cancer_upper, score_col, "ESTIMATE_EMT", sep = "_")]] <- stats_est
        est_prefix <- file.path(out_dir, "14_TME_high_low_plots", paste0(cancer_upper, "_", score_safe, "_ESTIMATE_EMT_epithelial_CAF_ECM_high_low_bar"))
        write.csv(stats_est, paste0(est_prefix, "_stats.csv"), row.names = FALSE)
        plot_high_low_bar_features(
          stats_est, feature_order = estimate_features, out_prefix = est_prefix,
          title = paste0(cancer_upper, ": ", score_col, " high vs low - ESTIMATE/EMT/epithelial/CAF/ECM"),
          width = 8, height = max(5, 0.35 * length(unique(stats_est$feature)) + 2)
        )
      }
    }

    if ("ESTIMATE_ImmuneScore" %in% colnames(merged2)) {
      immune_prefix <- file.path(out_dir, "14_TME_high_low_plots", paste0(cancer_upper, "_", score_safe, "_ESTIMATE_ImmuneScore_boxplot"))
      plot_single_high_low_feature(
        merged2, score_col = score_col, feature = "ESTIMATE_ImmuneScore", cancer = cancer_upper,
        out_prefix = immune_prefix,
        title = paste0(cancer_upper, ": ESTIMATE ImmuneScore by ", score_col, " high/low")
      )
    }
  }

  all_merged[[cancer_upper]] <- merged2

  numeric_df <- merged2[, sapply(merged2, is.numeric), drop = FALSE]
  if (ncol(numeric_df) >= 2 && nrow(numeric_df) >= TME_MIN_COR_PAIRS) {

    cor_mat <- suppressWarnings(cor(numeric_df, use = "pairwise.complete.obs", method = "spearman"))
    write.csv(cor_mat, file.path(out_dir, "03_correlation_tables", paste0(cancer_upper, "_all_numeric_spearman_matrix.csv")))

    targets <- intersect(c(
      "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
      "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
      "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED",
      "STEP_UP_TME_OVERLAP_FILTERED",
      "STEP_UP_MINUS_DOWN_RAW",
      "STEP_UP_residual_CAF3",
      "STEP_DOWN_residual_CAF3",
      "STEP_BALANCE_residual_CAF3",
      "STEP_UP_residual_CAF4",
      "STEP_DOWN_residual_CAF4",
      "STEP_BALANCE_residual_CAF4",
      "HUMAN_EPITHELIAL_SCORE_14",
      "HUMAN_MESENCHYMAL_SCORE_18",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
    ), colnames(merged2))

    cor_long <- cor_one_vs_all(merged2, targets, min_pairs = TME_MIN_COR_PAIRS, min_unique = TME_MIN_UNIQUE_VALUES)
    if (nrow(cor_long) > 0) {
      cor_long$cancer <- cancer_upper
      all_cor[[cancer_upper]] <- cor_long
      write.csv(cor_long, file.path(out_dir, "03_correlation_tables", paste0(cancer_upper, "_target_vs_all_spearman.csv")), row.names = FALSE)
    }

    key_vars <- intersect(c(
      "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
      "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
      "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED",
      "STEP_UP_TME_OVERLAP_FILTERED",
      "STEP_UP_MINUS_DOWN_RAW",
      "STEP_UP_residual_CAF3",
      "STEP_DOWN_residual_CAF3",
      "STEP_BALANCE_residual_CAF3",
      "STEP_UP_residual_CAF4",
      "STEP_DOWN_residual_CAF4",
      "STEP_BALANCE_residual_CAF4",
      "HUMAN_EPITHELIAL_SCORE_14",
      "HUMAN_MESENCHYMAL_SCORE_18",
      "KOHN_EMT_EPITHELIAL",
      "KOHN_EMT_MESENCHYMAL",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
      "CAF_3gene_score",
      "CAF_4gene_index",
      "ECM_score",
      "Fibroblast_score",
      "TGFb_score",
      "ESTIMATE",
      "ABSOLUTE",
      "LUMP",
      "IHC",
      "CPE",
      "ESTIMATE_StromalScore",
      "ESTIMATE_ImmuneScore",
      "ESTIMATE_ESTIMATEScore",
      "ESTIMATE_TumorPurity",
      tme_module_core_vars
    ), colnames(cor_mat))

    if (length(key_vars) >= 3) {
      hm <- cor_mat[key_vars, key_vars, drop = FALSE]

      safe_cor_heatmap(
        hm,
        out_pdf = file.path(out_dir, "04_correlation_heatmaps", paste0(cancer_upper, "_key_score_purity_spearman_heatmap.pdf")),
        out_png = file.path(out_dir, "04_correlation_heatmaps", paste0(cancer_upper, "_key_score_purity_spearman_heatmap.png")),
        title = paste0(cancer_upper, " key score/purity Spearman"),
        width = 8.5,
        height = 7.5,
        png_width = 2550,
        png_height = 2250
      )
    }
  }

  scatter_pairs <- list(
    c("STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"),
    c("STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", "STEP_UP_TME_OVERLAP_FILTERED"),
    c("CAF_3gene_score", "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"),
    c("ECM_score", "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"),
    c("ESTIMATE_StromalScore", "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"),
    c("CPE", "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"),
    c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"),
    c("HUMAN_MESENCHYMAL_SCORE_18", "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED"),
    c("CAF_3gene_score", "STEP_UP_TME_OVERLAP_FILTERED"),
    c("ECM_score", "STEP_UP_TME_OVERLAP_FILTERED"),
    c("ESTIMATE_StromalScore", "STEP_UP_TME_OVERLAP_FILTERED"),
    c("CPE", "STEP_UP_TME_OVERLAP_FILTERED"),
    c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "STEP_UP_TME_OVERLAP_FILTERED"),
    c("HUMAN_MESENCHYMAL_SCORE_18", "STEP_UP_TME_OVERLAP_FILTERED"),
    c("CAF_3gene_score", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("CAF_4gene_index", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("ECM_score", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("Fibroblast_score", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("TGFb_score", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("CPE", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("ESTIMATE", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("ESTIMATE_StromalScore", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("ESTIMATE_ImmuneScore", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("ESTIMATE_ESTIMATEScore", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("HUMAN_MESENCHYMAL_SCORE_18", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"),
    c("CAF_3gene_score", "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN"),
    c("CPE", "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN"),
    c("HUMAN_EPITHELIAL_SCORE_14", "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN"),
    c("CAF_3gene_score", "STEP_UP_residual_CAF3"),
    c("CAF_3gene_score", "STEP_DOWN_residual_CAF3"),
    c("CAF_3gene_score", "STEP_BALANCE_residual_CAF3"),
    c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "STEP_UP_residual_CAF3"),
    c("HUMAN_MESENCHYMAL_SCORE_18", "STEP_UP_residual_CAF3"),
    c("HUMAN_EPITHELIAL_SCORE_14", "STEP_DOWN_residual_CAF3"),
    c("CAF_4gene_index", "STEP_UP_residual_CAF4"),
    c("CAF_4gene_index", "STEP_DOWN_residual_CAF4"),
    c("CAF_4gene_index", "STEP_BALANCE_residual_CAF4")
  )

  for (pair in scatter_pairs) {
    x <- pair[1]
    y <- pair[2]
    safe <- paste0(cancer_upper, "_", gsub("[^A-Za-z0-9]+", "_", y), "_vs_", gsub("[^A-Za-z0-9]+", "_", x))
    plot_scatter(merged2, x, y, cancer_upper, file.path(out_dir, "05_scatter_plots", safe))
  }

  surv_scores <- intersect(c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
    "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED",
    "STEP_UP_TME_OVERLAP_FILTERED",
    "STEP_UP_MINUS_DOWN_RAW",
    "STEP_UP_residual_CAF3",
    "STEP_DOWN_residual_CAF3",
    "STEP_BALANCE_residual_CAF3",
    "STEP_UP_residual_CAF4",
    "STEP_DOWN_residual_CAF4",
    "STEP_BALANCE_residual_CAF4",
    "HUMAN_EPITHELIAL_SCORE_14",
    "HUMAN_MESENCHYMAL_SCORE_18",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "CAF_3gene_score",
    "CAF_4gene_index",
    "ECM_score",
    "Fibroblast_score",
    "TGFb_score"
  ), colnames(merged2))

  for (sc in surv_scores) {
    out_prefix <- file.path(out_dir, "08_survival", paste0(cancer_upper, "_", gsub("[^A-Za-z0-9]+", "_", sc)))
    tab <- tryCatch(
      run_survival(merged2, sc, cancer_upper, out_prefix),
      error = function(e) {
        msg("Survival failed for ", cancer_upper, " / ", sc, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(tab)) {
      all_surv[[paste(cancer_upper, sc, sep = "__")]] <- tab
    }
  }

  residual_surv_scores <- intersect(c(
    "STEP_UP_residual_CAF3",
    "STEP_DOWN_residual_CAF3",
    "STEP_BALANCE_residual_CAF3",
    "STEP_UP_residual_CAF4",
    "STEP_DOWN_residual_CAF4",
    "STEP_BALANCE_residual_CAF4"
  ), colnames(merged2))

  for (sc in residual_surv_scores) {
    out_prefix <- file.path(out_dir, "08b_survival_residual_CPE_AGE_models",
                            paste0(cancer_upper, "_", gsub("[^A-Za-z0-9]+", "_", sc)))
    tab2 <- tryCatch(
      run_survival_residual_cpe_age(merged2, sc, cancer_upper, out_prefix),
      error = function(e) {
        msg("Residual CPE+AGE survival failed for ", cancer_upper, " / ", sc, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(tab2)) {
      all_surv_residual_cpe_age[[paste(cancer_upper, sc, sep = "__")]] <- tab2
    }
  }

  endpoints <- detect_survival_endpoints(merged2)
  endpoint_scores <- intersect(c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_UP_residual_CAF3",
    "STEP_UP_residual_CAF4",
    "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
    "STEP_DOWN_residual_CAF3",
    "STEP_DOWN_residual_CAF4"
  ), colnames(merged2))

  if (length(endpoints) > 0 && length(endpoint_scores) > 0) {
    for (ep in names(endpoints)) {
      for (sc in endpoint_scores) {
        out_prefix <- file.path(out_dir, "08c_survival_all_endpoints",
                                paste0(cancer_upper, "_", ep, "_", gsub("[^A-Za-z0-9]+", "_", sc)))
        tab_ep <- tryCatch(
          run_survival_endpoint(merged2, sc, cancer_upper, ep,
                                time_col = endpoints[[ep]]$time,
                                event_col = endpoints[[ep]]$event,
                                out_prefix = out_prefix,
                                covariate_mode = "score_CPE_AGE"),
          error = function(e) {
            msg("Endpoint survival failed for ", cancer_upper, " / ", ep, " / ", sc, ": ", conditionMessage(e))
            NULL
          }
        )
        if (!is.null(tab_ep)) {
          all_surv_all_endpoints[[paste(cancer_upper, ep, sc, sep = "__")]] <- tab_ep
        }
      }
    }
  }
}

if (length(all_qc) > 0) {
  qc_all <- do.call(rbind, all_qc)
  write.csv(qc_all, file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_gene_set_QC_matched_missing.csv"), row.names = FALSE)
}

if (length(all_cor) > 0) {
  cor_all <- do.call(rbind, all_cor)
  write.csv(cor_all, file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_target_vs_all_spearman.csv"), row.names = FALSE)

  target_up <- "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"
  variables <- c(
    "HUMAN_MESENCHYMAL_SCORE_18",
    "HUMAN_EPITHELIAL_SCORE_14",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "CAF_3gene_score",
    "CAF_4gene_index",
    "ECM_score",
    "Fibroblast_score",
    "TGFb_score",
    "CPE",
    "ESTIMATE",
    "ABSOLUTE",
    "LUMP",
    "IHC",
    "ESTIMATE_StromalScore",
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_ESTIMATEScore",
    "ESTIMATE_TumorPurity"
  )

  heat_dat <- cor_all[cor_all$target == target_up & cor_all$variable %in% variables, ]
  if (nrow(heat_dat) > 0) {
    wide <- reshape(heat_dat[, c("cancer", "variable", "cor")],
                    idvar = "cancer", timevar = "variable", direction = "wide")
    rownames(wide) <- wide$cancer
    wide$cancer <- NULL
    colnames(wide) <- sub("^cor\\.", "", colnames(wide))
    wide <- as.matrix(wide)

    safe_cor_heatmap(
      wide,
      out_pdf = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_UP_key_correlations_heatmap.pdf"),
      out_png = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_UP_key_correlations_heatmap.png"),
      title = "Pan-cancer STEP_UP key correlations",
      width = 11,
      height = 8,
      png_width = 3300,
      png_height = 2400
    )
  }
}

  target_down <- "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN"
  variables_down <- c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED",
    "STEP_UP_TME_OVERLAP_FILTERED",
    "HUMAN_EPITHELIAL_SCORE_14",
    "HUMAN_MESENCHYMAL_SCORE_18",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "KOHN_EMT_EPITHELIAL",
    "KOHN_EMT_MESENCHYMAL",
    "CAF_3gene_score",
    "CAF_4gene_index",
    "ECM_score",
    "Fibroblast_score",
    "TGFb_score",
    "ESTIMATE_StromalScore",
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_ESTIMATEScore",
    "CPE",
    "ESTIMATE",
    "ABSOLUTE",
    "LUMP",
    "IHC"
  )

  heat_down <- cor_all[cor_all$target == target_down & cor_all$variable %in% variables_down, ]
  if (nrow(heat_down) > 0) {
    write.csv(heat_down, file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_DOWN_key_correlations_long.csv"), row.names = FALSE)

    wide_down <- reshape(heat_down[, c("cancer", "variable", "cor")],
                         idvar = "cancer", timevar = "variable", direction = "wide")
    rownames(wide_down) <- wide_down$cancer
    wide_down$cancer <- NULL
    colnames(wide_down) <- sub("^cor\\.", "", colnames(wide_down))
    wide_down <- as.matrix(wide_down)

    safe_cor_heatmap(
      wide_down,
      out_pdf = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_DOWN_key_correlations_heatmap.pdf"),
      out_png = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_DOWN_key_correlations_heatmap.png"),
      title = "Pan-cancer STEP_DOWN key correlations",
      width = 12,
      height = 8,
      png_width = 3600,
      png_height = 2400
    )
  }

  compact_targets <- c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_UP_residual_CAF3",
    "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
    "STEP_DOWN_residual_CAF3",
    "STEP_UP_MINUS_DOWN_RAW",
    "STEP_BALANCE_residual_CAF3"
  )
  compact_variables <- c(
    "CAF_3gene_score",
    "CAF_4gene_index",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HUMAN_MESENCHYMAL_SCORE_18",
    "HUMAN_EPITHELIAL_SCORE_14",
    "CPE",
    "ESTIMATE",
    "ESTIMATE_StromalScore",
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_ESTIMATEScore"
  )
  compact_cor <- cor_all[cor_all$target %in% compact_targets & cor_all$variable %in% compact_variables, ]
  if (nrow(compact_cor) > 0) {
    compact_cor$key <- paste(compact_cor$target, compact_cor$variable, sep = " | ")
    write.csv(compact_cor,
              file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_compact_correlations_long.csv"),
              row.names = FALSE)

    compact_wide <- reshape(compact_cor[, c("cancer", "key", "cor")],
                            idvar = "cancer", timevar = "key", direction = "wide")
    rownames(compact_wide) <- compact_wide$cancer
    compact_wide$cancer <- NULL
    colnames(compact_wide) <- sub("^cor\\.", "", colnames(compact_wide))
    compact_wide <- as.matrix(compact_wide)

    safe_cor_heatmap(
      compact_wide,
      out_pdf = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_compact_correlations_heatmap.pdf"),
      out_png = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_compact_correlations_heatmap.png"),
      title = "CAF-adjusted STEP compact correlations",
      width = 16,
      height = 8.5,
      png_width = 4800,
      png_height = 2550
    )
  }

  if (exists("compact_cor") && nrow(compact_cor) > 0) {
    compact_summary <- do.call(rbind, lapply(split(compact_cor, compact_cor$key), function(xx) {
      data.frame(
        key = unique(xx$key),
        target = unique(xx$target),
        variable = unique(xx$variable),
        n_cancers = sum(is.finite(xx$cor)),
        median_cor = median(xx$cor, na.rm = TRUE),
        mean_cor = mean(xx$cor, na.rm = TRUE),
        min_cor = min(xx$cor, na.rm = TRUE),
        max_cor = max(xx$cor, na.rm = TRUE),
        n_pos_0.2 = sum(xx$cor >= 0.2, na.rm = TRUE),
        n_neg_0.2 = sum(xx$cor <= -0.2, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    compact_summary <- compact_summary[order(compact_summary$target, compact_summary$variable), ]
    write.csv(compact_summary,
              file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_compact_correlations_summary.csv"),
              row.names = FALSE)
  }

  residual_targets <- c(
    "STEP_UP_MINUS_DOWN_RAW",
    "STEP_UP_residual_CAF3",
    "STEP_DOWN_residual_CAF3",
    "STEP_BALANCE_residual_CAF3",
    "STEP_UP_residual_CAF4",
    "STEP_DOWN_residual_CAF4",
    "STEP_BALANCE_residual_CAF4"
  )
  residual_variables <- c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
    "CAF_3gene_score",
    "CAF_4gene_index",
    "HUMAN_EPITHELIAL_SCORE_14",
    "HUMAN_MESENCHYMAL_SCORE_18",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "KOHN_EMT_EPITHELIAL",
    "KOHN_EMT_MESENCHYMAL",
    "ECM_score",
    "Fibroblast_score",
    "TGFb_score",
    "CPE",
    "ESTIMATE",
    "LUMP",
    "ABSOLUTE",
    "ESTIMATE_StromalScore",
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_ESTIMATEScore"
  )

  residual_cor <- cor_all[cor_all$target %in% residual_targets & cor_all$variable %in% residual_variables, ]
  if (nrow(residual_cor) > 0) {
    residual_cor$key <- paste(residual_cor$target, residual_cor$variable, sep = " | ")
    write.csv(residual_cor, file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_residual_correlations_long.csv"), row.names = FALSE)

    residual_wide <- reshape(residual_cor[, c("cancer", "key", "cor")],
                             idvar = "cancer", timevar = "key", direction = "wide")
    rownames(residual_wide) <- residual_wide$cancer
    residual_wide$cancer <- NULL
    colnames(residual_wide) <- sub("^cor\\.", "", colnames(residual_wide))
    residual_wide <- as.matrix(residual_wide)

    safe_cor_heatmap(
      residual_wide,
      out_pdf = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_residual_correlations_heatmap.pdf"),
      out_png = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_residual_correlations_heatmap.png"),
      title = "Pan-cancer CAF-adjusted STEP residual correlations",
      width = 24,
      height = 9,
      png_width = 7200,
      png_height = 2700
    )
  }

  up_down_targets <- c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
    "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED",
    "STEP_UP_TME_OVERLAP_FILTERED",
    "STEP_UP_MINUS_DOWN_RAW",
    "STEP_UP_residual_CAF3",
    "STEP_DOWN_residual_CAF3",
    "STEP_BALANCE_residual_CAF3",
    "STEP_UP_residual_CAF4",
    "STEP_DOWN_residual_CAF4",
    "STEP_BALANCE_residual_CAF4"
  )
  up_down_variables <- c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
    "STEP_UP_MINUS_DOWN_RAW",
    "STEP_UP_residual_CAF3",
    "STEP_DOWN_residual_CAF3",
    "STEP_BALANCE_residual_CAF3",
    "STEP_UP_residual_CAF4",
    "STEP_DOWN_residual_CAF4",
    "STEP_BALANCE_residual_CAF4",
    "HUMAN_EPITHELIAL_SCORE_14",
    "HUMAN_MESENCHYMAL_SCORE_18",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "CAF_3gene_score",
    "CAF_4gene_index",
    "ECM_score",
    "Fibroblast_score",
    "TGFb_score",
    "ESTIMATE_StromalScore",
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_ESTIMATEScore",
    "CPE",
    "ESTIMATE",
    "LUMP",
    "ABSOLUTE"
  )

  up_down <- cor_all[cor_all$target %in% up_down_targets & cor_all$variable %in% up_down_variables, ]
  if (nrow(up_down) > 0) {
    up_down$key <- paste(up_down$target, up_down$variable, sep = " | ")
    write.csv(up_down, file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_UP_DOWN_core_correlations_long.csv"), row.names = FALSE)

    up_down_wide <- reshape(up_down[, c("cancer", "key", "cor")],
                            idvar = "cancer", timevar = "key", direction = "wide")
    rownames(up_down_wide) <- up_down_wide$cancer
    up_down_wide$cancer <- NULL
    colnames(up_down_wide) <- sub("^cor\\.", "", colnames(up_down_wide))
    up_down_wide <- as.matrix(up_down_wide)

    safe_cor_heatmap(
      up_down_wide,
      out_pdf = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_UP_DOWN_core_correlations_heatmap.pdf"),
      out_png = file.path(out_dir, "09_pan_cancer_summary", "PanCancer_STEP_UP_DOWN_core_correlations_heatmap.png"),
      title = "Pan-cancer STEP_UP/DOWN core correlations",
      width = 22,
      height = 9,
      png_width = 6600,
      png_height = 2700
    )
  }

if (exists("cor_all") && nrow(cor_all) > 0) {
  compare_targets <- c(
    "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
    "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED",
    "STEP_UP_TME_OVERLAP_FILTERED"
  )
  compare_variables <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HUMAN_MESENCHYMAL_SCORE_18",
    "HUMAN_EPITHELIAL_SCORE_14",
    "CAF_3gene_score",
    "CAF_4gene_index",
    "ECM_score",
    "Fibroblast_score",
    "TGFb_score",
    "ESTIMATE_StromalScore",
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_ESTIMATEScore",
    "CPE",
    "ESTIMATE",
    "LUMP",
    "ABSOLUTE"
  )
  comp <- cor_all[cor_all$target %in% compare_targets & cor_all$variable %in% compare_variables, ]
  if (nrow(comp) > 0) {
    comp$key <- paste(comp$target, comp$variable, sep = " | ")
    write.csv(comp, file.path(out_dir, "09_pan_cancer_summary", "Original_vs_Filtered_STEP_UP_key_correlations_long.csv"), row.names = FALSE)

    comp_wide <- reshape(comp[, c("cancer", "key", "cor")], idvar = "cancer", timevar = "key", direction = "wide")
    rownames(comp_wide) <- comp_wide$cancer
    comp_wide$cancer <- NULL
    colnames(comp_wide) <- sub("^cor\\.", "", colnames(comp_wide))
    comp_wide <- as.matrix(comp_wide)

    safe_cor_heatmap(
      comp_wide,
      out_pdf = file.path(out_dir, "09_pan_cancer_summary", "Original_vs_Filtered_STEP_UP_key_correlations_heatmap.pdf"),
      out_png = file.path(out_dir, "09_pan_cancer_summary", "Original_vs_Filtered_STEP_UP_key_correlations_heatmap.png"),
      title = "Original vs filtered STEP_UP key correlations",
      width = 22,
      height = 9,
      png_width = 6600,
      png_height = 2700
    )
  }
}

if (exists("cor_all") && nrow(cor_all) > 0 && length(tme_module_all_vars) > 0) {
  tme_targets <- intersect(c("STEP_UP_residual_CAF3", "STEP_UP_residual_CAF4"), unique(cor_all$target))
  for (tt in tme_targets) {
    tt_safe <- gsub("[^A-Za-z0-9]+", "_", tt)

    all_vars_present <- unique(tme_module_all_vars[tme_module_all_vars %in% unique(cor_all$variable)])
    if (length(all_vars_present) > 0) {
      write_tme_cor_heatmap(
        cor_df = cor_all,
        target_col = tt,
        variables = all_vars_present,
        variable_order = all_vars_present,
        out_prefix = file.path(out_dir, "13_TME_correlations", paste0("PanCancer_", tt_safe, "_ALL_TME_plus_epithelial_correlation_heatmap")),
        title = paste0("Pan-cancer ", tt, " vs all TME signatures + epithelial"),
        width = 13,
        height = 12
      )
    }

    core_vars_present <- unique(tme_module_core_vars[tme_module_core_vars %in% unique(cor_all$variable)])
    if (length(core_vars_present) > 0) {
      write_tme_cor_heatmap(
        cor_df = cor_all,
        target_col = tt,
        variables = core_vars_present,
        variable_order = core_vars_present,
        out_prefix = file.path(out_dir, "13_TME_correlations", paste0("PanCancer_", tt_safe, "_CORE_TME_plus_epithelial_correlation_heatmap")),
        title = paste0("Pan-cancer ", tt, " vs core TME signatures + epithelial"),
        width = 12,
        height = 8.5
      )
    }
  }

  estimate_emt_vars <- c(
    "ESTIMATE_StromalScore", "ESTIMATE_ImmuneScore", "ESTIMATE_ESTIMATEScore", "ESTIMATE_TumorPurity",
    "HUMAN_EPITHELIAL_SCORE_14", "KOHN_EMT_EPITHELIAL",
    "HUMAN_MESENCHYMAL_SCORE_18", "KOHN_EMT_MESENCHYMAL",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "CAF_3gene_score", "CAF_4gene_index", "ECM_score"
  )
  estimate_emt_vars <- unique(estimate_emt_vars[estimate_emt_vars %in% unique(cor_all$variable)])
  for (tt in tme_targets) {
    tt_safe <- gsub("[^A-Za-z0-9]+", "_", tt)
    if (length(estimate_emt_vars) > 0) {
      write_tme_cor_heatmap(
        cor_df = cor_all,
        target_col = tt,
        variables = estimate_emt_vars,
        variable_order = estimate_emt_vars,
        out_prefix = file.path(out_dir, "13_TME_correlations", paste0("PanCancer_", tt_safe, "_ESTIMATE_EMT_epithelial_CAF_ECM_correlation_heatmap")),
        title = paste0("Pan-cancer ", tt, " vs ESTIMATE/EMT/epithelial/CAF/ECM"),
        width = 11,
        height = 7
      )
    }
  }
}

if (length(all_tme_qc) > 0) {
  tme_qc_all <- do.call(rbind, all_tme_qc)
  write.csv(tme_qc_all,
            file.path(out_dir, "13_TME_correlations", "ALL_cancers_TME_signature_QC_matched_missing.csv"),
            row.names = FALSE)
  write.csv(tme_qc_all,
            file.path(out_dir, "12_TME_scores", "ALL_cancers_TME_signature_threshold_audit_v11_2.csv"),
            row.names = FALSE)
  cd8_qc_all <- tme_qc_all[grepl("CD8", tme_qc_all$gene_set, ignore.case = TRUE), , drop = FALSE]
  if (nrow(cd8_qc_all) > 0) {
    write.csv(cd8_qc_all,
              file.path(out_dir, "12_TME_scores", "ALL_cancers_CD8_signature_threshold_audit_v11_2.csv"),
              row.names = FALSE)
  }
  tme_qc_summary <- aggregate(
    cbind(n_genes, n_matched, n_score_non_na, n_score_unique) ~ gene_set,
    data = tme_qc_all,
    FUN = function(x) stats::median(x, na.rm = TRUE)
  )
  names(tme_qc_summary)[-1] <- paste0("median_", names(tme_qc_summary)[-1])
  write.csv(tme_qc_summary,
            file.path(out_dir, "12_TME_scores", "PanCancer_TME_signature_median_match_and_score_summary_v11_2.csv"),
            row.names = FALSE)
}

if (length(all_tme_highlow) > 0) {
  highlow_all <- do.call(rbind, all_tme_highlow)
  highlow_all$FDR_by_score_feature_class <- ave(highlow_all$p, highlow_all$score, highlow_all$feature,
                                                FUN = function(x) p.adjust(x, method = "BH"))
  write.csv(highlow_all,
            file.path(out_dir, "14_TME_high_low_plots", "ALL_cancers_STEP_UP_residual_high_low_TME_stats.csv"),
            row.names = FALSE)

  for (score_col in unique(highlow_all$score)) {
    score_safe <- gsub("[^A-Za-z0-9]+", "_", score_col)
    for (feature_set_name in c("ALL_TME_plus_epithelial", "CORE_TME_plus_epithelial", "ESTIMATE_EMT_epithelial_CAF_ECM")) {
      vars <- switch(feature_set_name,
                     ALL_TME_plus_epithelial = tme_module_all_vars,
                     CORE_TME_plus_epithelial = tme_module_core_vars,
                     ESTIMATE_EMT_epithelial_CAF_ECM = c(
                       "ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore", "ESTIMATE_ESTIMATEScore", "ESTIMATE_TumorPurity",
                       "HUMAN_EPITHELIAL_SCORE_14", "KOHN_EMT_EPITHELIAL",
                       "HUMAN_MESENCHYMAL_SCORE_18", "KOHN_EMT_MESENCHYMAL",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                       "CAF_3gene_score", "CAF_4gene_index", "ECM_score"
                     ))
      dat <- highlow_all[highlow_all$score == score_col & highlow_all$feature %in% vars, , drop = FALSE]
      if (nrow(dat) == 0) next
      dat <- dat[!duplicated(dat[, c("cancer", "feature")]), , drop = FALSE]
      wide <- reshape(dat[, c("feature", "cancer", "median_diff_high_minus_low")],
                      idvar = "feature", timevar = "cancer", direction = "wide")
      rownames(wide) <- wide$feature
      wide$feature <- NULL
      colnames(wide) <- sub("^median_diff_high_minus_low\\.", "", colnames(wide))
      mat <- as.matrix(wide)
      mode(mat) <- "numeric"
      vars_present <- unique(vars[vars %in% rownames(mat)])
      mat <- mat[vars_present, , drop = FALSE]
      if (nrow(mat) > 0 && ncol(mat) > 0) {
        out_prefix <- file.path(out_dir, "14_TME_high_low_plots",
                                paste0("PanCancer_", score_safe, "_", feature_set_name, "_median_diff_high_minus_low_heatmap"))
        write.csv(mat, paste0(out_prefix, "_matrix.csv"))
        safe_cor_heatmap(
          mat,
          out_pdf = paste0(out_prefix, ".pdf"),
          out_png = paste0(out_prefix, ".png"),
          title = paste0(score_col, " high-low median difference: ", feature_set_name),
          width = ifelse(feature_set_name == "ALL_TME_plus_epithelial", 13, 11),
          height = ifelse(feature_set_name == "ALL_TME_plus_epithelial", 12, 8),
          png_width = ifelse(feature_set_name == "ALL_TME_plus_epithelial", 3900, 3300),
          png_height = ifelse(feature_set_name == "ALL_TME_plus_epithelial", 3600, 2400)
        )
      }
    }
  }
}

if (length(all_residual_qc) > 0) {
  residual_qc_all <- do.call(rbind, all_residual_qc)
  residual_qc_all$FDR_CAF_by_residual_score <- ave(residual_qc_all$p_CAF, residual_qc_all$residual_score,
                                                    FUN = function(x) p.adjust(x, method = "BH"))
  write.csv(residual_qc_all,
            file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_CAF_adjusted_residual_score_QC.csv"),
            row.names = FALSE)
}

if (length(all_surv) > 0) {
  surv_all <- do.call(rbind, all_surv)
  surv_all$FDR_by_score <- ave(surv_all$p, surv_all$score, FUN = function(x) p.adjust(x, method = "BH"))
  write.csv(surv_all, file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_OS_cox_summary.csv"), row.names = FALSE)
}

if (length(all_surv_residual_cpe_age) > 0) {
  surv_resid <- do.call(rbind, all_surv_residual_cpe_age)
  surv_resid$FDR_by_score <- ave(surv_resid$p, surv_resid$score, FUN = function(x) p.adjust(x, method = "BH"))
  write.csv(surv_resid,
            file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_OS_cox_summary_residual_CPE_AGE.csv"),
            row.names = FALSE)

  resid_sig <- surv_resid[
    surv_resid$model == "score_scaled" &
      is.finite(surv_resid$p) &
      surv_resid$p < 0.05,
    , drop = FALSE
  ]
  resid_sig <- resid_sig[order(resid_sig$score, resid_sig$FDR_by_score, resid_sig$p), ]
  write.csv(resid_sig,
            file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_OS_cox_summary_residual_CPE_AGE_nominal_significant.csv"),
            row.names = FALSE)

  resid_sig_fdr <- surv_resid[
    surv_resid$model == "score_scaled" &
      is.finite(surv_resid$FDR_by_score) &
      surv_resid$FDR_by_score < 0.05,
    , drop = FALSE
  ]
  resid_sig_fdr <- resid_sig_fdr[order(resid_sig_fdr$score, resid_sig_fdr$FDR_by_score, resid_sig_fdr$p), ]
  write.csv(resid_sig_fdr,
            file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_OS_cox_summary_residual_CPE_AGE_FDR_significant.csv"),
            row.names = FALSE)
}

if (length(all_surv_all_endpoints) > 0) {
  endpoint_surv_all <- do.call(rbind, all_surv_all_endpoints)
  endpoint_surv_all$FDR_by_endpoint_score <- ave(
    endpoint_surv_all$p,
    endpoint_surv_all$endpoint,
    endpoint_surv_all$score,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  write.csv(endpoint_surv_all,
            file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_all_endpoints_cox_summary.csv"),
            row.names = FALSE)

  endpoint_unadj <- endpoint_surv_all[endpoint_surv_all$model == "univariate_score", , drop = FALSE]
  if (nrow(endpoint_unadj) > 0) {
    write.csv(endpoint_unadj,
              file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_all_endpoints_cox_unadjusted_score_only.csv"),
              row.names = FALSE)

    endpoint_unadj_sig <- endpoint_unadj[
      is.finite(endpoint_unadj$p) & endpoint_unadj$p < 0.05,
      , drop = FALSE
    ]
    endpoint_unadj_sig <- endpoint_unadj_sig[order(endpoint_unadj_sig$endpoint, endpoint_unadj_sig$score,
                                                   endpoint_unadj_sig$FDR_by_endpoint_score, endpoint_unadj_sig$p), ]
    write.csv(endpoint_unadj_sig,
              file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_all_endpoints_cox_unadjusted_nominal_significant.csv"),
              row.names = FALSE)

    endpoint_unadj_sig_fdr <- endpoint_unadj[
      is.finite(endpoint_unadj$FDR_by_endpoint_score) & endpoint_unadj$FDR_by_endpoint_score < 0.05,
      , drop = FALSE
    ]
    endpoint_unadj_sig_fdr <- endpoint_unadj_sig_fdr[order(endpoint_unadj_sig_fdr$endpoint, endpoint_unadj_sig_fdr$score,
                                                           endpoint_unadj_sig_fdr$FDR_by_endpoint_score, endpoint_unadj_sig_fdr$p), ]
    write.csv(endpoint_unadj_sig_fdr,
              file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_all_endpoints_cox_unadjusted_FDR_significant.csv"),
              row.names = FALSE)
  }

  endpoint_sig <- endpoint_surv_all[
    endpoint_surv_all$model == "score_scaled" &
      is.finite(endpoint_surv_all$p) &
      endpoint_surv_all$p < 0.05,
    , drop = FALSE
  ]
  endpoint_sig <- endpoint_sig[order(endpoint_sig$endpoint, endpoint_sig$score,
                                     endpoint_sig$FDR_by_endpoint_score, endpoint_sig$p), ]
  write.csv(endpoint_sig,
            file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_all_endpoints_cox_nominal_significant.csv"),
            row.names = FALSE)

  endpoint_sig_fdr <- endpoint_surv_all[
    endpoint_surv_all$model == "score_scaled" &
      is.finite(endpoint_surv_all$FDR_by_endpoint_score) &
      endpoint_surv_all$FDR_by_endpoint_score < 0.05,
    , drop = FALSE
  ]
  endpoint_sig_fdr <- endpoint_sig_fdr[order(endpoint_sig_fdr$endpoint, endpoint_sig_fdr$score,
                                             endpoint_sig_fdr$FDR_by_endpoint_score, endpoint_sig_fdr$p), ]
  write.csv(endpoint_sig_fdr,
            file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_all_endpoints_cox_FDR_significant.csv"),
            row.names = FALSE)
}

if (length(all_merged) > 0) {
  all_cols <- unique(c("cancer_loop", unlist(lapply(all_merged, colnames))))
  all_merged2 <- do.call(rbind, lapply(names(all_merged), function(nm) {
    x <- all_merged[[nm]]
    x$cancer_loop <- nm
    miss <- setdiff(all_cols, colnames(x))
    for (cc in miss) x[[cc]] <- NA
    x[, all_cols, drop = FALSE]
  }))
  write.csv(all_merged2, file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_scores_purity_clinical_merged.csv"), row.names = FALSE)
}

if (nrow(estimate_failed_cancers) > 0) {
  write.csv(estimate_failed_cancers, file.path(out_dir, "09_pan_cancer_summary", "ESTIMATE_failed_cancers.csv"), row.names = FALSE)
  write.csv(estimate_failed_cancers, file.path(out_dir, "00_logs", "ESTIMATE_failed_cancers.csv"), row.names = FALSE)
}

if (nrow(skipped_cancers) > 0) {
  write.csv(skipped_cancers, file.path(out_dir, "09_pan_cancer_summary", "skipped_cancers.csv"), row.names = FALSE)
  write.csv(skipped_cancers, file.path(out_dir, "00_logs", "skipped_cancers.csv"), row.names = FALSE)
}

match_files <- list.files(file.path(out_dir, "00_logs"), pattern = "_sample_matching_summary\\.csv$", full.names = TRUE)
if (length(match_files) > 0) {
  match_all <- do.call(rbind, lapply(match_files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
  write.csv(match_all, file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_sample_matching_summary.csv"), row.names = FALSE)
}

tryCatch({
  compact_long_file <- file.path(out_dir, "09_pan_cancer_summary", "PanCancer_CAF_adjusted_STEP_compact_correlations_long.csv")
  if (file.exists(compact_long_file)) {
    compact_long_df <- read.csv(compact_long_file, check.names = FALSE)
    plot_caf_adjustment_validation(compact_long_df, out_dir)
  }
}, error = function(e) {
  msg("Figure-ready CAF adjustment validation plot failed: ", conditionMessage(e))
})

tryCatch({
  endpoint_cox_file <- file.path(out_dir, "09_pan_cancer_summary", "ALL_cancers_all_endpoints_cox_summary.csv")
  if (file.exists(endpoint_cox_file)) {
    endpoint_cox_df <- read.csv(endpoint_cox_file, check.names = FALSE)
    for (ep in intersect(c("OS", "DFS", "PFS", "RFS"), unique(endpoint_cox_df$endpoint))) {
      plot_residual_cox_forest(endpoint_cox_df, out_dir, endpoint = ep)
      plot_all_cancer_endpoint_forest(endpoint_cox_df, out_dir, endpoint = ep, model_type = "adjusted")
      plot_all_cancer_endpoint_forest(endpoint_cox_df, out_dir, endpoint = ep, model_type = "unadjusted")
    }
  }
}, error = function(e) {
  msg("Figure-ready Cox forest plot failed: ", conditionMessage(e))
})

input_manifest_v112 <- data.frame(
  module = c(
    "Main TME score/correlation/high-low",
    "Immune activity extra",
    "Checkpoint/immunotherapy",
    "TIDE batch preferred",
    "TIDE per-cancer fallback",
    "Shared TCGA expression",
    "Shared merged STEP table"
  ),
  required_input = c(
    "TME_CIBERSORTConsensusMCP.gmt",
    "immune_activity_extra.gmt",
    "checkpoint_immunotherapy_gene_sets.gmt",
    "ALL_TCGA_TIDE_results_merged.csv or .tsv",
    "*_TIDE_result.txt files",
    "*_tcga_expression.gct files",
    "ALL_cancers_scores_purity_clinical_merged.csv"
  ),
  automatic_expected_location = c(
    file.path(base_dir, "gene set"),
    file.path(base_dir, "gene set"),
    file.path(base_dir, "gene set"),
    file.path(base_dir, "18_TIDEpy_results", "tables"),
    paste(file.path(base_dir, "18_TIDEpy_results", "02_TIDE_results_by_cancer"),
          file.path(base_dir, "CIBERSORT_upload_inputs", "linear_normalized_for_CIBERSORT"), sep = " OR "),
    file.path(base_dir, "TCGA_201906cBioP", "exp_gct"),
    file.path(out_dir, "09_pan_cancer_summary")
  ),
  manual_override_variable = c(
    "INPUT_OVERRIDE_TME_GMT",
    "INPUT_OVERRIDE_IMMUNE_ACTIVITY_GMT",
    "INPUT_OVERRIDE_CHECKPOINT_GMT",
    "INPUT_OVERRIDE_TIDE_BATCH_FILE",
    "INPUT_OVERRIDE_TIDE_RESULT_DIR",
    "INPUT_OVERRIDE_EXP_GCT_DIR",
    "generated by main pipeline"
  ),
  stringsAsFactors = FALSE
)
write.csv(input_manifest_v112,
          file.path(out_dir, "00_logs", "v11_2_required_input_manifest.csv"),
          row.names = FALSE)

msg("\nMain v11.2 TME/high-low module finished: ", as.character(Sys.time()))
msg("Output directory: ", out_dir)

MASTER_RESULT_DIR <- normalizePath(out_dir, mustWork = FALSE)
MASTER_PROJECT_ROOT <- normalizePath(base_dir, mustWork = FALSE)
message("MASTER_RESULT_DIR: ", MASTER_RESULT_DIR)
message("MASTER_PROJECT_ROOT: ", MASTER_PROJECT_ROOT)

if (isTRUE(RUN_CLINICAL_MODULE)) {
  run_optional_module_v11("15 clinical", {

    local({
      message("========== 15_clinical_sex_stage_CAF3: unified module ==========")
      result_dir <- MASTER_RESULT_DIR
      input_file <- file.path(result_dir, "09_pan_cancer_summary",
                              "ALL_cancers_scores_purity_clinical_merged.csv")
      if (!file.exists(input_file)) {
        warning("Clinical module skipped: missing ", input_file)
      } else {
        clin <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
        out15 <- file.path(result_dir, "15_clinical_sex_stage_CAF3")
        dir.create(out15, recursive = TRUE, showWarnings = FALSE)
        per_dir <- file.path(out15, "per_cancer")
        sum_dir15 <- file.path(out15, "all_cancer_summary")
        dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(sum_dir15, recursive = TRUE, showWarnings = FALSE)

        pick_first <- function(patterns, nms) {
          for (pat in patterns) {
            hit <- grep(pat, nms, value = TRUE, ignore.case = TRUE)
            if (length(hit) > 0) return(hit[1])
          }
          NA_character_
        }

        cancer_col <- pick_first(c("^cancer_loop$", "^cancer_from_file$",
                                   "^cancer_type$", "^cancer$", "^CANCER$"), colnames(clin))
        score_col <- pick_first(c("^STEP_UP_residual_CAF3$",
                                  "STEP.*UP.*residual.*CAF3"), colnames(clin))
        sex_col <- pick_first(c("^SEX$", "^gender$", "^sex$"), colnames(clin))
        stage_col <- pick_first(c("^STAGE_STD$", "^STAGE$", "^STAGE_RAW$",
                                  "AJCC.*STAGE", "PATHOLOGIC.*STAGE",
                                  "CLINICAL.*STAGE", "tumor_stage"), colnames(clin))
        age_col <- pick_first(c("^AGE$", "AGE_AT_DIAGNOSIS", "AGE.*DIAG",
                                "diagnosis_age"), colnames(clin))

        if (is.na(cancer_col) || is.na(score_col)) {
          warning("Clinical module skipped: cannot detect cancer or STEP_UP_residual_CAF3 column.")
        } else {
          clean_chr <- function(x) {
            x <- trimws(as.character(x))
            x[x %in% c("", "NA", "N/A", "[Not Available]", "[Not Applicable]",
                       "[Unknown]", "Unknown", "Not Reported")] <- NA_character_
            x
          }

          standardize_sex_v11 <- function(x) {
            y <- toupper(clean_chr(x))
            out <- rep(NA_character_, length(y))
            out[grepl("^M$|MALE", y)] <- "Male"
            out[grepl("^F$|FEMALE", y)] <- "Female"
            factor(out, levels = c("Male", "Female"))
          }

          standardize_stage_v11 <- function(x) {
            y <- toupper(clean_chr(x))
            y <- gsub("STAGE", "", y)
            y <- gsub("[^A-Z0-9]", "", y)
            out <- rep(NA_character_, length(y))
            out[grepl("^0", y)] <- "0"
            out[grepl("^I($|A|B|C|1)", y) | grepl("^1", y)] <- "I"
            out[grepl("^II($|A|B|C|2)", y) | grepl("^2", y)] <- "II"
            out[grepl("^III($|A|B|C|3)", y) | grepl("^3", y)] <- "III"
            out[grepl("^IV($|A|B|C|4)", y) | grepl("^4", y)] <- "IV"
            factor(out, levels = c("0", "I", "II", "III", "IV"))
          }

          age_numeric_v11 <- function(x) {
            y <- suppressWarnings(as.numeric(as.character(x)))
            
            if (sum(is.finite(y)) > 0 && median(y, na.rm = TRUE) > 200) y <- y / 365.25
            y
          }

          safe_comp_test_v11 <- function(tab) {
            if (length(tab) == 0 || any(dim(tab) < 2) || sum(tab) == 0) {
              return(list(test = "NA", p = NA_real_))
            }
            cs <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
            use_fisher <- is.null(cs) || any(cs$expected < 5)
            if (use_fisher) {
              ft <- tryCatch(fisher.test(tab, simulate.p.value = any(dim(tab) > 2),
                                         B = 10000), error = function(e) NULL)
              return(list(test = "Fisher", p = if (is.null(ft)) NA_real_ else ft$p.value))
            }
            list(test = "Chi-square", p = cs$p.value)
          }

          save_plot_v11 <- function(p, base, width = 6, height = 4.8) {
            print(p)
            ggsave(paste0(base, ".pdf"), p, width = width, height = height,
                   limitsize = FALSE, useDingbats = FALSE)
            ggsave(paste0(base, ".png"), p, width = width, height = height,
                   dpi = 300, limitsize = FALSE)
          }

          composition_plot_v11 <- function(dat, category_col, category_colors,
                                           test_label, title, fill_title) {
            tmp <- dat %>%
              filter(!is.na(.data[[category_col]]), !is.na(STEP_group)) %>%
              count(STEP_group, category = .data[[category_col]], name = "n") %>%
              group_by(STEP_group) %>%
              mutate(prop = n / sum(n)) %>%
              ungroup()
            if (nrow(tmp) == 0) return(NULL)

            p <- ggplot(tmp, aes(x = STEP_group, y = prop, fill = category)) +
              geom_col(width = 0.70, color = "grey35", linewidth = 0.2) +
              geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
                        position = position_stack(vjust = 0.5), size = 3,
                        color = "black", check_overlap = TRUE) +
              annotate("text", x = 1.5, y = 1.08, label = test_label,
                       size = 3.2, fontface = "plain") +
              scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                 limits = c(0, 1.14),
                                 expand = expansion(mult = c(0, 0))) +
              labs(title = wrap_title_step(title, 54),
                   x = "CAF3-corrected STEP_UP group",
                   y = "Composition", fill = fill_title) +
              theme_step(base_size = 10, title_size = 10,
                          axis_title_size = 9, axis_text_size = 8) +
              theme(legend.position = "top")
            if (!is.null(category_colors)) {
              p <- p + scale_fill_manual(values = category_colors, drop = FALSE)
            }
            p
          }

          boxplot_v11 <- function(dat, x_col, y_col, fill_values, p_label,
                                  title, xlab, ylab, x_order = NULL) {
            d <- dat %>% filter(!is.na(.data[[x_col]]), is.finite(.data[[y_col]]))
            if (nrow(d) == 0) return(NULL)
            if (!is.null(x_order)) d[[x_col]] <- factor(d[[x_col]], levels = x_order)
            yr <- range(d[[y_col]], na.rm = TRUE)
            span <- diff(yr)
            if (!is.finite(span) || span == 0) span <- 1
            p <- ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]],
                               fill = .data[[x_col]], color = .data[[x_col]])) +
              geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82,
                           linewidth = 0.35) +
              geom_jitter(width = 0.15, size = 0.8, alpha = 0.42,
                          show.legend = FALSE) +
              annotate("text", x = (length(unique(na.omit(d[[x_col]]))) + 1) / 2,
                       y = yr[2] + 0.11 * span, label = p_label, size = 3.2) +
              scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
              labs(title = wrap_title_step(title, 54), x = xlab, y = ylab) +
              theme_step(base_size = 10, title_size = 10,
                          axis_title_size = 9, axis_text_size = 8) +
              theme(legend.position = "none")
            if (!is.null(fill_values)) {
              p <- p +
                scale_fill_manual(values = fill_values, drop = FALSE) +
                scale_color_manual(values = fill_values, drop = FALSE)
            }
            p
          }

          clin$cancer_v11 <- toupper(as.character(clin[[cancer_col]]))
          clin$score_v11 <- suppressWarnings(as.numeric(clin[[score_col]]))
          clin$sex_v11 <- if (!is.na(sex_col)) standardize_sex_v11(clin[[sex_col]]) else NA
          clin$stage_v11 <- if (!is.na(stage_col)) standardize_stage_v11(clin[[stage_col]]) else NA
          clin$age_v11 <- if (!is.na(age_col)) age_numeric_v11(clin[[age_col]]) else NA_real_

          clin <- clin %>%
            filter(!is.na(cancer_v11), cancer_v11 != "", is.finite(score_v11)) %>%
            group_by(cancer_v11) %>%
            mutate(
              cutoff_v11 = median(score_v11, na.rm = TRUE),
              STEP_group = factor(ifelse(score_v11 > cutoff_v11, "High", "Low"),
                                   levels = c("Low", "High"))
            ) %>%
            ungroup()

          clin$stage_binary_v11 <- ifelse(clin$stage_v11 %in% c("0", "I", "II"),
                                          "Early (0/I/II)",
                                          ifelse(clin$stage_v11 %in% c("III", "IV"),
                                                 "Advanced (III/IV)", NA_character_))
          clin$stage_binary_v11 <- factor(clin$stage_binary_v11,
                                          levels = c("Early (0/I/II)",
                                                     "Advanced (III/IV)"))

          stage_colors <- c("0" = "#D9D9D9", "I" = "#A6CEE3", "II" = "#6BAED6",
                            "III" = "#FB6A4A", "IV" = "#CB181D")
          sex_colors <- c("Male" = STEP_LOW_COLOR, "Female" = STEP_HIGH_COLOR)
          stage_binary_colors <- c("Early (0/I/II)" = STEP_LOW_COLOR,
                                   "Advanced (III/IV)" = STEP_HIGH_COLOR)

          all_stats15 <- list()
          sex_comp_all <- list()
          stage_comp_all <- list()
          stagebin_comp_all <- list()
          sex_box_all <- list()
          stage_box_all <- list()
          stagebin_box_all <- list()
          age_box_all <- list()
          kk <- 1

          for (ca in sort(unique(clin$cancer_v11))) {
            sub <- clin[clin$cancer_v11 == ca, , drop = FALSE]
            ca_dir <- file.path(per_dir, ca)
            dir.create(ca_dir, recursive = TRUE, showWarnings = FALSE)
            write.csv(sub, file.path(ca_dir, paste0(ca, "_clinical_CAF3_plot_data.csv")),
                      row.names = FALSE)

            sex_dat <- sub %>% filter(!is.na(sex_v11), !is.na(STEP_group))
            if (nrow(sex_dat) >= 20 && n_distinct(sex_dat$sex_v11) >= 2 &&
                n_distinct(sex_dat$STEP_group) >= 2) {
              comp <- safe_comp_test_v11(table(sex_dat$STEP_group, sex_dat$sex_v11))
              comp_lab <- paste0(comp$test, " ", format_p_step(comp$p))
              p1 <- composition_plot_v11(
                sex_dat, "sex_v11", sex_colors, comp_lab,
                paste0(ca, ": sex composition by CAF3-corrected STEP_UP"),
                "Sex"
              )
              save_plot_v11(p1, file.path(ca_dir,
                paste0(ca, "_CAF3_STEP_UP_high_low_sex_composition")), 5.7, 4.8)

              pw <- tryCatch(wilcox.test(score_v11 ~ sex_v11, data = sex_dat,
                                        exact = FALSE)$p.value,
                             error = function(e) NA_real_)
              p2 <- boxplot_v11(
                sex_dat, "sex_v11", "score_v11", sex_colors,
                format_p_step(pw, "Wilcoxon p = "),
                paste0(ca, ": CAF3-corrected STEP_UP by sex"),
                "Sex", score_col, c("Male", "Female")
              )
              save_plot_v11(p2, file.path(ca_dir,
                paste0(ca, "_CAF3_STEP_UP_score_by_sex_boxplot")), 5.2, 4.8)

              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "STEP_highlow_vs_sex_composition",
                test = comp$test, p = comp$p, n = nrow(sex_dat),
                effect = NA_real_, stringsAsFactors = FALSE); kk <- kk + 1
              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "score_by_sex", test = "Wilcoxon", p = pw,
                n = nrow(sex_dat),
                effect = median(sex_dat$score_v11[sex_dat$sex_v11 == "Female"],
                                na.rm = TRUE) -
                         median(sex_dat$score_v11[sex_dat$sex_v11 == "Male"],
                                na.rm = TRUE),
                stringsAsFactors = FALSE); kk <- kk + 1

              sex_comp_all[[ca]] <- sex_dat
              sex_box_all[[ca]] <- sex_dat
            }

            stage_dat <- sub %>% filter(!is.na(stage_v11), !is.na(STEP_group))
            if (nrow(stage_dat) >= 20 && n_distinct(stage_dat$stage_v11) >= 2 &&
                n_distinct(stage_dat$STEP_group) >= 2) {
              comp <- safe_comp_test_v11(table(stage_dat$STEP_group,
                                               stage_dat$stage_v11))
              comp_lab <- paste0(comp$test, " ", format_p_step(comp$p))
              p3 <- composition_plot_v11(
                stage_dat, "stage_v11", stage_colors, comp_lab,
                paste0(ca, ": stage composition by CAF3-corrected STEP_UP"),
                "Stage"
              )
              save_plot_v11(p3, file.path(ca_dir,
                paste0(ca, "_CAF3_STEP_UP_high_low_stage_composition")), 6.3, 4.8)

              pk <- tryCatch(kruskal.test(score_v11 ~ stage_v11,
                                         data = stage_dat)$p.value,
                             error = function(e) NA_real_)
              p4 <- boxplot_v11(
                stage_dat, "stage_v11", "score_v11", stage_colors,
                format_p_step(pk, "Kruskal-Wallis p = "),
                paste0(ca, ": CAF3-corrected STEP_UP by stage"),
                "Pathologic/clinical stage", score_col,
                c("0", "I", "II", "III", "IV")
              )
              save_plot_v11(p4, file.path(ca_dir,
                paste0(ca, "_CAF3_STEP_UP_score_by_stage_boxplot")), 6.3, 4.8)

              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "STEP_highlow_vs_stage_composition",
                test = comp$test, p = comp$p, n = nrow(stage_dat),
                effect = NA_real_, stringsAsFactors = FALSE); kk <- kk + 1
              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "score_by_stage", test = "Kruskal-Wallis", p = pk,
                n = nrow(stage_dat), effect = NA_real_,
                stringsAsFactors = FALSE); kk <- kk + 1

              stage_comp_all[[ca]] <- stage_dat
              stage_box_all[[ca]] <- stage_dat
            }

            sb_dat <- sub %>% filter(!is.na(stage_binary_v11), !is.na(STEP_group))
            if (nrow(sb_dat) >= 20 && n_distinct(sb_dat$stage_binary_v11) >= 2 &&
                n_distinct(sb_dat$STEP_group) >= 2) {
              comp <- safe_comp_test_v11(table(sb_dat$STEP_group,
                                               sb_dat$stage_binary_v11))
              comp_lab <- paste0(comp$test, " ", format_p_step(comp$p))
              p5 <- composition_plot_v11(
                sb_dat, "stage_binary_v11", stage_binary_colors, comp_lab,
                paste0(ca, ": early/advanced stage by CAF3-corrected STEP_UP"),
                "Stage group"
              )
              save_plot_v11(p5, file.path(ca_dir,
                paste0(ca, "_CAF3_STEP_UP_high_low_early_advanced_stage_composition")),
                6.3, 4.8)

              pw <- tryCatch(wilcox.test(score_v11 ~ stage_binary_v11,
                                        data = sb_dat, exact = FALSE)$p.value,
                             error = function(e) NA_real_)
              p6 <- boxplot_v11(
                sb_dat, "stage_binary_v11", "score_v11",
                stage_binary_colors,
                format_p_step(pw, "Wilcoxon p = "),
                paste0(ca, ": CAF3-corrected STEP_UP by early/advanced stage"),
                "Stage group", score_col,
                c("Early (0/I/II)", "Advanced (III/IV)")
              )
              save_plot_v11(p6, file.path(ca_dir,
                paste0(ca, "_CAF3_STEP_UP_score_by_early_advanced_stage_boxplot")),
                6.4, 4.8)

              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "STEP_highlow_vs_early_advanced_stage_composition",
                test = comp$test, p = comp$p, n = nrow(sb_dat),
                effect = NA_real_, stringsAsFactors = FALSE); kk <- kk + 1
              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "score_by_early_advanced_stage",
                test = "Wilcoxon", p = pw, n = nrow(sb_dat),
                effect = median(sb_dat$score_v11[
                           sb_dat$stage_binary_v11 == "Advanced (III/IV)"],
                           na.rm = TRUE) -
                         median(sb_dat$score_v11[
                           sb_dat$stage_binary_v11 == "Early (0/I/II)"],
                           na.rm = TRUE),
                stringsAsFactors = FALSE); kk <- kk + 1

              stagebin_comp_all[[ca]] <- sb_dat
              stagebin_box_all[[ca]] <- sb_dat
            }

            age_dat <- sub %>% filter(is.finite(age_v11), !is.na(STEP_group))
            if (nrow(age_dat) >= 20 && n_distinct(age_dat$STEP_group) >= 2) {
              pw <- tryCatch(wilcox.test(age_v11 ~ STEP_group,
                                        data = age_dat, exact = FALSE)$p.value,
                             error = function(e) NA_real_)
              ct <- tryCatch(cor.test(age_dat$score_v11, age_dat$age_v11,
                                     method = "spearman", exact = FALSE),
                             error = function(e) NULL)
              rho <- if (is.null(ct)) NA_real_ else unname(ct$estimate)
              pcor <- if (is.null(ct)) NA_real_ else ct$p.value
              p7 <- boxplot_v11(
                age_dat, "STEP_group", "age_v11", STEP_GROUP_COLORS,
                format_p_step(pw, "Wilcoxon p = "),
                paste0(ca, ": age by CAF3-corrected STEP_UP group"),
                "CAF3-corrected STEP_UP group", "Age (years)",
                c("Low", "High")
              )
              save_plot_v11(p7, file.path(ca_dir,
                paste0(ca, "_CAF3_STEP_UP_high_low_age_boxplot")), 5.2, 4.8)

              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "age_by_STEP_highlow", test = "Wilcoxon",
                p = pw, n = nrow(age_dat),
                effect = median(age_dat$age_v11[age_dat$STEP_group == "High"],
                                na.rm = TRUE) -
                         median(age_dat$age_v11[age_dat$STEP_group == "Low"],
                                na.rm = TRUE),
                stringsAsFactors = FALSE); kk <- kk + 1
              all_stats15[[kk]] <- data.frame(cancer = ca,
                analysis = "score_vs_age_correlation",
                test = paste0("Spearman rho=", signif(rho, 3)),
                p = pcor, n = nrow(age_dat), effect = rho,
                stringsAsFactors = FALSE); kk <- kk + 1
              age_box_all[[ca]] <- age_dat
            }
          }

          make_faceted_composition <- function(dat_list, category_col, colors,
                                               title, filename, legend_title,
                                               stats_analysis) {
            if (length(dat_list) == 0) return(invisible(NULL))
            d <- bind_rows(dat_list)
            comp <- d %>%
              filter(!is.na(.data[[category_col]]), !is.na(STEP_group)) %>%
              count(cancer_v11, STEP_group,
                    category = .data[[category_col]], name = "n") %>%
              group_by(cancer_v11, STEP_group) %>%
              mutate(prop = n / sum(n)) %>%
              ungroup()
            st <- if (length(all_stats15) > 0) bind_rows(all_stats15) else data.frame()
            ann <- st %>%
              filter(analysis == stats_analysis) %>%
              transmute(cancer_v11 = cancer, x = 1.5, y = 1.08,
                        label = format_p_step(p))
            p <- ggplot(comp, aes(x = STEP_group, y = prop, fill = category)) +
              geom_col(width = 0.72, color = "grey35", linewidth = 0.15) +
              geom_text(aes(label = scales::percent(prop, accuracy = 1)),
                        position = position_stack(vjust = 0.5), size = 1.9,
                        color = "black", check_overlap = TRUE) +
              geom_text(data = ann, aes(x = x, y = y, label = label),
                        inherit.aes = FALSE, size = 2.3) +
              facet_wrap(~ cancer_v11, ncol = 5) +
              scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                 limits = c(0, 1.14),
                                 expand = expansion(mult = c(0, 0))) +
              scale_fill_manual(values = colors, drop = FALSE) +
              labs(title = wrap_title_step(title, 80),
                   x = "CAF3-corrected STEP_UP group",
                   y = "Composition", fill = legend_title) +
              theme_step(base_size = 8.5, title_size = 10,
                          axis_title_size = 8, axis_text_size = 6.8,
                          strip_size = 7.2) +
              theme(legend.position = "top",
                    axis.text.x = element_text(angle = 35, hjust = 1))
            nca <- n_distinct(comp$cancer_v11)
            save_plot_v11(p, file.path(sum_dir15, filename),
                          width = 13, height = max(6, ceiling(nca / 5) * 2.6))
          }

          make_faceted_box <- function(dat_list, x_col, y_col, colors,
                                       title, filename, stats_analysis,
                                       x_order = NULL, ylab = score_col) {
            if (length(dat_list) == 0) return(invisible(NULL))
            d <- bind_rows(dat_list)
            if (!is.null(x_order)) d[[x_col]] <- factor(d[[x_col]], levels = x_order)
            st <- if (length(all_stats15) > 0) bind_rows(all_stats15) else data.frame()
            ann <- st %>%
              filter(analysis == stats_analysis) %>%
              transmute(cancer_v11 = cancer, label = format_p_step(p))
            ypos <- d %>%
              group_by(cancer_v11) %>%
              summarise(ymax = max(.data[[y_col]], na.rm = TRUE),
                        ymin = min(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
              mutate(y = ymax + 0.12 * ifelse(ymax > ymin, ymax - ymin, 1),
                     x = (length(x_order %||% unique(d[[x_col]])) + 1) / 2)
            ann <- left_join(ann, ypos, by = "cancer_v11")
            p <- ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]],
                               fill = .data[[x_col]], color = .data[[x_col]])) +
              geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82,
                           linewidth = 0.25) +
              geom_jitter(width = 0.12, size = 0.35, alpha = 0.28,
                          show.legend = FALSE) +
              geom_text(data = ann, aes(x = x, y = y, label = label),
                        inherit.aes = FALSE, size = 2.2, color = "black") +
              facet_wrap(~ cancer_v11, scales = "free_y", ncol = 5) +
              scale_fill_manual(values = colors, drop = FALSE) +
              scale_color_manual(values = colors, drop = FALSE) +
              scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
              labs(title = wrap_title_step(title, 80), x = NULL, y = ylab) +
              theme_step(base_size = 8.5, title_size = 10,
                          axis_title_size = 8, axis_text_size = 6.8,
                          strip_size = 7.2) +
              theme(legend.position = "none",
                    axis.text.x = element_text(angle = 35, hjust = 1))
            nca <- n_distinct(d$cancer_v11)
            save_plot_v11(p, file.path(sum_dir15, filename),
                          width = 13, height = max(6, ceiling(nca / 5) * 2.7))
          }

          `%||%` <- function(a, b) if (!is.null(a)) a else b

          make_faceted_composition(
            sex_comp_all, "sex_v11", sex_colors,
            "All cancers: sex composition by CAF3-corrected STEP_UP group",
            "ALL_cancers_sex_composition_faceted_with_percent_and_p",
            "Sex", "STEP_highlow_vs_sex_composition"
          )
          make_faceted_composition(
            stage_comp_all, "stage_v11", stage_colors,
            "All cancers: stage composition by CAF3-corrected STEP_UP group",
            "ALL_cancers_stage_composition_faceted_with_percent_and_p",
            "Stage", "STEP_highlow_vs_stage_composition"
          )
          make_faceted_composition(
            stagebin_comp_all, "stage_binary_v11", stage_binary_colors,
            "All cancers: early/advanced stage composition by CAF3-corrected STEP_UP group",
            "ALL_cancers_early_advanced_stage_composition_faceted_with_percent_and_p",
            "Stage group", "STEP_highlow_vs_early_advanced_stage_composition"
          )
          make_faceted_box(
            sex_box_all, "sex_v11", "score_v11", sex_colors,
            "All cancers: CAF3-corrected STEP_UP score by sex",
            "ALL_cancers_score_by_sex_boxplot_faceted_with_p",
            "score_by_sex", c("Male", "Female"), score_col
          )
          make_faceted_box(
            stage_box_all, "stage_v11", "score_v11", stage_colors,
            "All cancers: CAF3-corrected STEP_UP score by stage",
            "ALL_cancers_score_by_stage_boxplot_faceted_with_p",
            "score_by_stage", c("0", "I", "II", "III", "IV"), score_col
          )
          make_faceted_box(
            stagebin_box_all, "stage_binary_v11", "score_v11",
            stage_binary_colors,
            "All cancers: CAF3-corrected STEP_UP score by early/advanced stage",
            "ALL_cancers_score_by_early_advanced_stage_boxplot_faceted_with_p",
            "score_by_early_advanced_stage",
            c("Early (0/I/II)", "Advanced (III/IV)"), score_col
          )
          make_faceted_box(
            age_box_all, "STEP_group", "age_v11", STEP_GROUP_COLORS,
            "All cancers: age by CAF3-corrected STEP_UP group",
            "ALL_cancers_age_by_STEP_group_boxplot_faceted_with_p",
            "age_by_STEP_highlow", c("Low", "High"), "Age (years)"
          )

          if (length(all_stats15) > 0) {
            stats15 <- bind_rows(all_stats15)
            stats15$FDR_by_analysis <- ave(stats15$p, stats15$analysis,
                                          FUN = function(x) p.adjust(x, method = "BH"))
            stats15$p_label <- vapply(stats15$p, format_p_short_step, character(1))
            stats15$neglog10FDR <- -log10(pmax(stats15$FDR_by_analysis, 1e-300))
            write.csv(stats15,
                      file.path(out15,
                        "ALL_cancers_CAF3_STEP_UP_clinical_association_stats.csv"),
                      row.names = FALSE)

            p_sum <- ggplot(stats15,
                            aes(x = factor(cancer, levels = sort(unique(cancer))),
                                y = analysis, fill = neglog10FDR)) +
              geom_tile(color = "white", linewidth = 0.25) +
              geom_text(aes(label = paste0("p=", p_label)), size = 2.0,
                        color = "black") +
              scale_fill_gradient(low = "#FFF7EC", high = "#D7301F",
                                  name = expression(-log[10](FDR))) +
              labs(title = "CAF3-corrected STEP_UP clinical associations across cancers",
                   subtitle = "Numbers shown in each tile are raw p values",
                   x = "TCGA cancer type", y = NULL) +
              theme_step(base_size = 9, title_size = 10,
                          axis_title_size = 8, axis_text_size = 7.2,
                          strip_size = 7) +
              theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
            save_plot_v11(p_sum,
              file.path(sum_dir15,
                        "ALL_cancers_clinical_association_pvalue_labeled_heatmap"),
              width = max(12, length(unique(stats15$cancer)) * 0.45),
              height = max(5.5, length(unique(stats15$analysis)) * 0.65 + 2))
          }
          message("Clinical module finished: ", out15)
        }
      }
    })

  })
}

if (isTRUE(RUN_ADVANCED_CLINICAL_MODULE)) {
  run_optional_module_v11("15 advanced clinical", {
    local({
      message("========== 15_clinical_sex_stage_CAF3: advanced additions ==========")

      result_dir <- MASTER_RESULT_DIR
      in_file <- file.path(result_dir, "09_pan_cancer_summary",
                           "ALL_cancers_scores_purity_clinical_merged.csv")
      if (!file.exists(in_file)) stop("Missing clinical-score table: ", in_file)

      dat <- read.csv(in_file, check.names = FALSE, stringsAsFactors = FALSE)
      out_root <- file.path(result_dir, "15_clinical_sex_stage_CAF3",
                            "advanced_clinical_additions")
      per_dir <- file.path(out_root, "per_cancer")
      sum_dir <- file.path(out_root, "all_cancer_summary")
      tab_dir <- file.path(out_root, "tables")
      for (dd in c(out_root, per_dir, sum_dir, tab_dir)) {
        dir.create(dd, recursive = TRUE, showWarnings = FALSE)
      }

      pick_exact_or_pattern <- function(exact, patterns = character()) {
        hit <- exact[exact %in% colnames(dat)]
        if (length(hit) > 0) return(hit[1])
        for (pat in patterns) {
          hit <- grep(pat, colnames(dat), value = TRUE, ignore.case = TRUE)
          if (length(hit) > 0) return(hit[1])
        }
        NA_character_
      }

      cancer_col <- pick_exact_or_pattern(
        c("cancer_loop", "cancer_from_file", "cancer_type", "cancer", "CANCER")
      )
      score_col <- pick_exact_or_pattern(
        c("STEP_UP_residual_CAF3", "STEP_UP_residual_CAF3_z"),
        c("STEP.*UP.*residual.*CAF3")
      )
      stage_std_col <- pick_exact_or_pattern(
        c("STAGE_STD", "STAGE", "STAGE_RAW"),
        c("AJCC.*STAGE", "PATHOLOGIC.*STAGE", "CLINICAL.*STAGE")
      )
      stage_num_col <- pick_exact_or_pattern(c("STAGE_NUM"))
      stage_group_col <- pick_exact_or_pattern(c("STAGE_GROUP"))
      t_col <- pick_exact_or_pattern(c("T_STAGE_BEST_STD", "PATH_T_MAJOR",
                                       "PATH_T_STD", "CLIN_T_MAJOR", "CLIN_T_STD"))
      n_col <- pick_exact_or_pattern(c("N_STAGE_BEST_STD", "PATH_N_MAJOR",
                                       "PATH_N_STD", "CLIN_N_MAJOR", "CLIN_N_STD"))
      m_col <- pick_exact_or_pattern(c("M_STAGE_BEST_STD", "PATH_M_MAJOR",
                                       "PATH_M_STD", "CLIN_M_MAJOR", "CLIN_M_STD"))
      grade_col <- pick_exact_or_pattern(c("GRADE_GROUP", "GRADE_STD", "GRADE_RAW"))
      gleason_col <- pick_exact_or_pattern(c("GLEASON_SCORE"))
      primary_col <- pick_exact_or_pattern(c("IS_PRIMARY_TUMOR"))

      if (is.na(cancer_col) || is.na(score_col)) {
        stop("Advanced clinical module cannot detect cancer or CAF3-corrected STEP score.")
      }

      dat$cancer_adv <- toupper(as.character(dat[[cancer_col]]))
      dat$score_adv <- suppressWarnings(as.numeric(dat[[score_col]]))
      if (!is.na(primary_col)) {
        pflag <- tolower(as.character(dat[[primary_col]]))
        keep_primary <- pflag %in% c("true", "t", "1", "yes", "y")
        if (sum(keep_primary, na.rm = TRUE) > 0) dat <- dat[keep_primary, , drop = FALSE]
      }
      dat <- dat[!is.na(dat$cancer_adv) & dat$cancer_adv != "" &
                   is.finite(dat$score_adv), , drop = FALSE]

      clean_stage <- function(x) {
        y <- toupper(trimws(as.character(x)))
        y[y %in% c("", "NA", "N/A", "UNKNOWN", "[NOT AVAILABLE]",
                   "[UNKNOWN]", "NOT REPORTED")] <- NA_character_
        y <- gsub("^STAGE[ _-]*", "", y)
        y <- gsub("[^A-Z0-9]", "", y)
        out <- rep(NA_character_, length(y))
        out[grepl("^0", y)] <- "0"
        out[grepl("^I($|A|B|C|1)", y) | grepl("^1", y)] <- "I"
        out[grepl("^II($|A|B|C|2)", y) | grepl("^2", y)] <- "II"
        out[grepl("^III($|A|B|C|3)", y) | grepl("^3", y)] <- "III"
        out[grepl("^IV($|A|B|C|4)", y) | grepl("^4", y)] <- "IV"
        factor(out, levels = c("0", "I", "II", "III", "IV"), ordered = TRUE)
      }

      clean_tnm <- function(x, letter) {
        y <- toupper(trimws(as.character(x)))
        y[y %in% c("", "NA", "N/A", "NX", "MX", "TX", "UNKNOWN",
                   "[NOT AVAILABLE]", "[UNKNOWN]")] <- NA_character_
        y <- gsub(paste0("^P?", letter), "", y)
        y <- gsub("[^0-9A-Z]", "", y)
        major <- sub("^([0-9]+).*", "\\1", y)
        major[!grepl("^[0-9]+", y)] <- NA_character_
        paste0(letter, major)
      }

      clean_grade <- function(x) {
        y <- toupper(trimws(as.character(x)))
        y[y %in% c("", "NA", "N/A", "GX", "UNKNOWN", "[NOT AVAILABLE]",
                   "[UNKNOWN]")] <- NA_character_
        g <- sub(".*?([1-4]).*", "\\1", y)
        g[!grepl("[1-4]", y)] <- NA_character_
        factor(paste0("G", g), levels = c("G1", "G2", "G3", "G4"), ordered = TRUE)
      }

      dat$stage_adv <- if (!is.na(stage_std_col)) clean_stage(dat[[stage_std_col]]) else NA
      dat$stage_num_adv <- if (!is.na(stage_num_col)) {
        suppressWarnings(as.numeric(dat[[stage_num_col]]))
      } else {
        as.numeric(dat$stage_adv) - 1
      }
      dat$stage_group_adv <- if (!is.na(stage_group_col)) {
        z <- trimws(as.character(dat[[stage_group_col]]))
        z[grepl("early|I-II|I/II", z, ignore.case = TRUE)] <- "Early"
        z[grepl("advanced|III-IV|III/IV", z, ignore.case = TRUE)] <- "Advanced"
        z
      } else {
        ifelse(dat$stage_adv %in% c("0", "I", "II"), "Early",
               ifelse(dat$stage_adv %in% c("III", "IV"), "Advanced", NA_character_))
      }
      dat$stage_group_adv <- factor(dat$stage_group_adv,
                                    levels = c("Early", "Advanced"))
      dat$T_adv <- if (!is.na(t_col)) clean_tnm(dat[[t_col]], "T") else NA_character_
      dat$N_adv <- if (!is.na(n_col)) clean_tnm(dat[[n_col]], "N") else NA_character_
      dat$M_adv <- if (!is.na(m_col)) clean_tnm(dat[[m_col]], "M") else NA_character_
      dat$Grade_adv <- if (!is.na(grade_col)) clean_grade(dat[[grade_col]]) else NA
      dat$Gleason_adv <- if (!is.na(gleason_col)) suppressWarnings(as.numeric(dat[[gleason_col]])) else NA_real_

      dat <- dat %>%
        group_by(cancer_adv) %>%
        mutate(
          score_z_adv = as.numeric(scale(score_adv)),
          STEP_group_adv = factor(
            ifelse(score_adv > median(score_adv, na.rm = TRUE), "High", "Low"),
            levels = c("Low", "High")
          )
        ) %>%
        ungroup()

      save_both_adv <- function(p, base, width = 6, height = 4.8) {
        print(p)
        ggsave(paste0(base, ".pdf"), p, width = width, height = height,
               limitsize = FALSE, useDingbats = FALSE)
        ggsave(paste0(base, ".png"), p, width = width, height = height,
               dpi = 300, limitsize = FALSE)
      }

      make_box_adv <- function(d, xvar, title, p_text, colors = NULL,
                               x_order = NULL, ylab = score_col) {
        d <- d[!is.na(d[[xvar]]) & is.finite(d$score_adv), , drop = FALSE]
        if (nrow(d) == 0 || length(unique(d[[xvar]])) < 2) return(NULL)
        if (!is.null(x_order)) d[[xvar]] <- factor(d[[xvar]], levels = x_order)
        yr <- range(d$score_adv, na.rm = TRUE)
        span <- diff(yr)
        if (!is.finite(span) || span == 0) span <- 1
        p <- ggplot(d, aes(x = .data[[xvar]], y = score_adv,
                           fill = .data[[xvar]], color = .data[[xvar]])) +
          geom_boxplot(outlier.shape = NA, width = 0.58, alpha = 0.82,
                       linewidth = 0.3) +
          geom_jitter(width = 0.12, size = 0.55, alpha = 0.35,
                      show.legend = FALSE) +
          annotate("text", x = (length(unique(na.omit(d[[xvar]]))) + 1) / 2,
                   y = yr[2] + 0.12 * span, label = p_text, size = 3) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
          labs(title = wrap_title_step(title, 58), x = NULL, y = ylab) +
          theme_step(base_size = 10, title_size = 9.5,
                      axis_title_size = 8.5, axis_text_size = 7.5) +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 35, hjust = 1))
        if (!is.null(colors)) {
          p <- p + scale_fill_manual(values = colors, drop = FALSE) +
            scale_color_manual(values = colors, drop = FALSE)
        }
        p
      }

      make_composition_adv <- function(d, category, title, p_text, colors = NULL) {
        dd <- d %>%
          filter(!is.na(.data[[category]]), !is.na(STEP_group_adv)) %>%
          count(STEP_group_adv, category_value = .data[[category]], name = "n") %>%
          group_by(STEP_group_adv) %>%
          mutate(prop = n / sum(n)) %>%
          ungroup()
        if (nrow(dd) == 0) return(NULL)
        p <- ggplot(dd, aes(x = STEP_group_adv, y = prop, fill = category_value)) +
          geom_col(width = 0.70, color = "grey35", linewidth = 0.2) +
          geom_text(aes(label = scales::percent(prop, accuracy = 1)),
                    position = position_stack(vjust = 0.5),
                    size = 2.5, color = "black", check_overlap = TRUE) +
          annotate("text", x = 1.5, y = 1.08, label = p_text, size = 3) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(0, 1.14),
                             expand = expansion(mult = c(0, 0))) +
          labs(title = wrap_title_step(title, 58), x = "STEP group",
               y = "Composition", fill = NULL) +
          theme_step(base_size = 10, title_size = 9.5,
                      axis_title_size = 8.5, axis_text_size = 7.5) +
          theme(legend.position = "top")
        if (!is.null(colors)) p <- p + scale_fill_manual(values = colors, drop = FALSE)
        p
      }

      safe_comp_p_adv <- function(d, category) {
        dd <- d[!is.na(d[[category]]) & !is.na(d$STEP_group_adv), , drop = FALSE]
        if (nrow(dd) < 10 || length(unique(dd[[category]])) < 2 ||
            length(unique(dd$STEP_group_adv)) < 2) return(NA_real_)
        tab <- table(dd$STEP_group_adv, dd[[category]])
        cs <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
        if (is.null(cs) || any(cs$expected < 5)) {
          return(tryCatch(fisher.test(tab, simulate.p.value = any(dim(tab) > 2),
                                      B = 10000)$p.value,
                          error = function(e) NA_real_))
        }
        cs$p.value
      }

      stage_colors <- c("0" = "#D9D9D9", "I" = "#A6CEE3",
                        "II" = "#6BAED6", "III" = "#FB6A4A", "IV" = "#CB181D")
      t_colors <- c("T0" = "#D9D9D9", "T1" = "#A6CEE3", "T2" = "#6BAED6",
                    "T3" = "#FB6A4A", "T4" = "#CB181D")
      n_colors <- c("N0" = STEP_LOW_COLOR, "N1" = "#FDAE6B",
                    "N2" = "#F16913", "N3" = STEP_HIGH_COLOR)
      m_colors <- c("M0" = STEP_LOW_COLOR, "M1" = STEP_HIGH_COLOR)
      grade_colors <- c("G1" = "#A6CEE3", "G2" = "#6BAED6",
                        "G3" = "#FB6A4A", "G4" = "#CB181D")

      stats_list <- list()
      logistic_list <- list()
      facet_data <- list(stage = list(), stage_group = list(), T = list(),
                         N = list(), M = list(), Grade = list(),
                         Gleason = list())
      k <- 1
      j <- 1

      for (ca in sort(unique(dat$cancer_adv))) {
        d <- dat[dat$cancer_adv == ca, , drop = FALSE]
        ca_dir <- file.path(per_dir, ca)
        dir.create(ca_dir, recursive = TRUE, showWarnings = FALSE)

        ds <- d[!is.na(d$stage_adv) & is.finite(d$score_adv), , drop = FALSE]
        if (nrow(ds) >= 20 && length(unique(ds$stage_adv)) >= 2) {
          p_kw <- tryCatch(kruskal.test(score_adv ~ stage_adv, data = ds)$p.value,
                           error = function(e) NA_real_)
          ct <- tryCatch(cor.test(ds$score_adv, ds$stage_num_adv,
                                  method = "spearman", exact = FALSE),
                         error = function(e) NULL)
          rho <- if (is.null(ct)) NA_real_ else unname(ct$estimate)
          p_trend <- if (is.null(ct)) NA_real_ else ct$p.value
          label <- paste0("Kruskal-Wallis ", format_p_step(p_kw),
                          "\nTrend rho=", signif(rho, 3), ", ",
                          format_p_step(p_trend))
          p <- make_box_adv(ds, "stage_adv",
                            paste0(ca, ": STEP_UP residual across stage"),
                            label, stage_colors, c("0", "I", "II", "III", "IV"))
          if (!is.null(p)) save_both_adv(
            p, file.path(ca_dir, paste0(ca, "_score_by_stage_with_trend")),
            6.2, 5.0
          )
          stats_list[[k]] <- data.frame(
            cancer = ca, analysis = "overall_stage_Kruskal",
            effect = NA_real_, p = p_kw, n = nrow(ds),
            stringsAsFactors = FALSE); k <- k + 1
          stats_list[[k]] <- data.frame(
            cancer = ca, analysis = "overall_stage_ordered_trend",
            effect = rho, p = p_trend, n = nrow(ds),
            stringsAsFactors = FALSE); k <- k + 1
          facet_data$stage[[ca]] <- ds
        }

        div <- d[d$stage_adv %in% c("I", "IV") & is.finite(d$score_adv), , drop = FALSE]
        if (sum(div$stage_adv == "I") >= 5 && sum(div$stage_adv == "IV") >= 5) {
          div$stage_I_IV <- factor(as.character(div$stage_adv), levels = c("I", "IV"))
          p_iv <- tryCatch(wilcox.test(score_adv ~ stage_I_IV, data = div,
                                      exact = FALSE)$p.value,
                           error = function(e) NA_real_)
          p <- make_box_adv(div, "stage_I_IV",
                            paste0(ca, ": Stage I versus Stage IV"),
                            format_p_step(p_iv, "Wilcoxon p = "),
                            c("I" = STEP_LOW_COLOR, "IV" = STEP_HIGH_COLOR),
                            c("I", "IV"))
          if (!is.null(p)) save_both_adv(
            p, file.path(ca_dir, paste0(ca, "_Stage_I_vs_IV_boxplot")),
            5.2, 4.8
          )
          stats_list[[k]] <- data.frame(
            cancer = ca, analysis = "Stage_I_vs_IV",
            effect = median(div$score_adv[div$stage_I_IV == "IV"], na.rm = TRUE) -
              median(div$score_adv[div$stage_I_IV == "I"], na.rm = TRUE),
            p = p_iv, n = nrow(div), stringsAsFactors = FALSE); k <- k + 1
        }

        dg <- d[!is.na(d$stage_group_adv) & is.finite(d$score_adv), , drop = FALSE]
        if (sum(dg$stage_group_adv == "Early") >= 5 &&
            sum(dg$stage_group_adv == "Advanced") >= 5) {
          p_g <- tryCatch(wilcox.test(score_adv ~ stage_group_adv, data = dg,
                                     exact = FALSE)$p.value,
                          error = function(e) NA_real_)
          p <- make_box_adv(
            dg, "stage_group_adv",
            paste0(ca, ": Early versus advanced stage"),
            format_p_step(p_g, "Wilcoxon p = "),
            c("Early" = STEP_LOW_COLOR, "Advanced" = STEP_HIGH_COLOR),
            c("Early", "Advanced")
          )
          if (!is.null(p)) save_both_adv(
            p, file.path(ca_dir, paste0(ca, "_Early_vs_Advanced_boxplot")),
            5.3, 4.8
          )

          dg$advanced_binary <- as.integer(dg$stage_group_adv == "Advanced")
          fit <- tryCatch(glm(advanced_binary ~ score_z_adv,
                              family = binomial(), data = dg),
                          error = function(e) NULL)
          if (!is.null(fit)) {
            sm <- summary(fit)$coefficients
            if ("score_z_adv" %in% rownames(sm)) {
              beta <- sm["score_z_adv", "Estimate"]
              se <- sm["score_z_adv", "Std. Error"]
              p_log <- sm["score_z_adv", "Pr(>|z|)"]
              logistic_list[[j]] <- data.frame(
                cancer = ca, outcome = "Advanced_vs_Early",
                n = nrow(dg), OR_per_SD = exp(beta),
                CI_low = exp(beta - 1.96 * se),
                CI_high = exp(beta + 1.96 * se),
                p = p_log, stringsAsFactors = FALSE
              ); j <- j + 1
            }
          }
          stats_list[[k]] <- data.frame(
            cancer = ca, analysis = "Early_vs_Advanced",
            effect = median(dg$score_adv[dg$stage_group_adv == "Advanced"], na.rm = TRUE) -
              median(dg$score_adv[dg$stage_group_adv == "Early"], na.rm = TRUE),
            p = p_g, n = nrow(dg), stringsAsFactors = FALSE); k <- k + 1
          facet_data$stage_group[[ca]] <- dg
        }

        dgscore <- d[is.finite(d$Gleason_adv) & is.finite(d$score_adv), , drop = FALSE]
        if (nrow(dgscore) >= 20 && length(unique(dgscore$Gleason_adv)) >= 2) {
          dgscore$Gleason_factor <- factor(
            dgscore$Gleason_adv,
            levels = sort(unique(dgscore$Gleason_adv))
          )
          p_gl <- tryCatch(
            kruskal.test(score_adv ~ Gleason_factor, data = dgscore)$p.value,
            error = function(e) NA_real_
          )
          gleason_levels <- levels(dgscore$Gleason_factor)
          gleason_cols <- setNames(
            grDevices::colorRampPalette(c("#A6CEE3", "#CB181D"))(
              length(gleason_levels)
            ),
            gleason_levels
          )
          p <- make_box_adv(
            dgscore, "Gleason_factor",
            paste0(ca, ": STEP_UP residual by Gleason score"),
            format_p_step(p_gl, "Kruskal-Wallis p = "),
            gleason_cols, gleason_levels
          )
          if (!is.null(p)) save_both_adv(
            p, file.path(ca_dir, paste0(ca, "_score_by_Gleason")),
            6.0, 4.9
          )
          stats_list[[k]] <- data.frame(
            cancer = ca, analysis = "Gleason_score_Kruskal",
            effect = NA_real_, p = p_gl, n = nrow(dgscore),
            stringsAsFactors = FALSE
          ); k <- k + 1
          facet_data$Gleason[[ca]] <- dgscore
        }

        for (spec in list(
          list(var = "T_adv", label = "T stage", colors = t_colors, bucket = "T"),
          list(var = "N_adv", label = "N stage", colors = n_colors, bucket = "N"),
          list(var = "M_adv", label = "M stage", colors = m_colors, bucket = "M"),
          list(var = "Grade_adv", label = "Histologic grade",
               colors = grade_colors, bucket = "Grade")
        )) {
          vv <- spec$var
          dx <- d[!is.na(d[[vv]]) & is.finite(d$score_adv), , drop = FALSE]
          if (nrow(dx) < 20 || length(unique(dx[[vv]])) < 2) next
          p_kw <- tryCatch(kruskal.test(dx$score_adv, as.factor(dx[[vv]]))$p.value,
                           error = function(e) NA_real_)
          p <- make_box_adv(
            dx, vv, paste0(ca, ": STEP_UP residual by ", spec$label),
            format_p_step(p_kw, "Kruskal-Wallis p = "),
            spec$colors
          )
          if (!is.null(p)) save_both_adv(
            p, file.path(ca_dir, paste0(ca, "_score_by_", spec$bucket)),
            6.0, 4.9
          )

          p_comp <- safe_comp_p_adv(dx, vv)
          pc <- make_composition_adv(
            dx, vv, paste0(ca, ": ", spec$label,
                           " composition by STEP group"),
            format_p_step(p_comp, "Composition p = "),
            spec$colors
          )
          if (!is.null(pc)) save_both_adv(
            pc, file.path(ca_dir, paste0(ca, "_", spec$bucket,
                                        "_composition_by_STEP_group")),
            5.8, 4.8
          )

          stats_list[[k]] <- data.frame(
            cancer = ca, analysis = paste0(spec$bucket, "_stage_Kruskal"),
            effect = NA_real_, p = p_kw, n = nrow(dx),
            stringsAsFactors = FALSE); k <- k + 1
          stats_list[[k]] <- data.frame(
            cancer = ca, analysis = paste0(spec$bucket,
                                           "_composition_by_STEP_group"),
            effect = NA_real_, p = p_comp, n = nrow(dx),
            stringsAsFactors = FALSE); k <- k + 1
          facet_data[[spec$bucket]][[ca]] <- dx
        }
      }

      stats_adv <- if (length(stats_list) > 0) bind_rows(stats_list) else data.frame()
      logistic_adv <- if (length(logistic_list) > 0) bind_rows(logistic_list) else data.frame()
      if (nrow(stats_adv) > 0) {
        stats_adv$FDR_by_analysis <- ave(
          stats_adv$p, stats_adv$analysis,
          FUN = function(x) p.adjust(x, method = "BH")
        )
        write.csv(stats_adv,
                  file.path(tab_dir, "ALL_cancers_advanced_clinical_stats.csv"),
                  row.names = FALSE)

        stats_adv$label <- format_p_short_step(stats_adv$p)
        stats_adv$fill_value <- -log10(pmax(stats_adv$FDR_by_analysis, 1e-300))
        p_heat <- ggplot(stats_adv,
                         aes(x = factor(cancer,
                                        levels = sort(unique(cancer))),
                             y = analysis, fill = fill_value)) +
          geom_tile(color = "white", linewidth = 0.2) +
          geom_text(aes(label = paste0("p=", label)), size = 1.9,
                    color = "black") +
          scale_fill_gradient(low = "#FFF7EC", high = "#D7301F",
                              name = expression(-log[10](FDR))) +
          labs(title = "Advanced clinical associations of CAF3-corrected STEP_UP",
               subtitle = "Raw p values are printed in tiles",
               x = "TCGA cancer", y = NULL) +
          theme_step(base_size = 8.5, title_size = 9.5,
                      axis_title_size = 8, axis_text_size = 6.8) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1))
        save_both_adv(
          p_heat,
          file.path(sum_dir, "ALL_cancers_advanced_clinical_pvalue_heatmap"),
          max(12, 0.43 * length(unique(stats_adv$cancer))),
          max(6, 0.55 * length(unique(stats_adv$analysis)) + 2)
        )
      }

      if (nrow(logistic_adv) > 0) {
        logistic_adv$FDR <- p.adjust(logistic_adv$p, method = "BH")
        write.csv(
          logistic_adv,
          file.path(tab_dir,
                    "ALL_cancers_advanced_stage_logistic_OR_per_SD_STEP.csv"),
          row.names = FALSE
        )
        logistic_adv <- logistic_adv %>%
          filter(is.finite(OR_per_SD), is.finite(CI_low), is.finite(CI_high),
                 OR_per_SD > 0, CI_low > 0, CI_high > 0)
        if (nrow(logistic_adv) > 0) {
          logistic_adv$cancer <- factor(
            logistic_adv$cancer,
            levels = rev(logistic_adv$cancer[order(logistic_adv$OR_per_SD)])
          )
          logistic_adv$p_label <- format_p_step(logistic_adv$p)
          p_for <- ggplot(logistic_adv,
                          aes(y = cancer, x = OR_per_SD,
                              xmin = CI_low, xmax = CI_high)) +
            geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
            geom_errorbarh(height = 0.18, linewidth = 0.45) +
            geom_point(aes(fill = p < 0.05), shape = 21, size = 2.4) +
            geom_text(aes(x = CI_high * 1.08, label = p_label),
                      hjust = 0, size = 2.5) +
            scale_x_log10() +
            scale_fill_manual(values = c("FALSE" = "white",
                                         "TRUE" = STEP_HIGH_COLOR),
                              guide = "none") +
            coord_cartesian(clip = "off") +
            labs(
              title = "Association of CAF3-corrected STEP_UP with advanced stage",
              subtitle = "Logistic regression; odds ratio per 1-SD increase in STEP",
              x = "Odds ratio (95% CI, log scale)", y = NULL
            ) +
            theme_step(base_size = 9, title_size = 9.5,
                        axis_title_size = 8.5, axis_text_size = 7.2) +
            theme(plot.margin = margin(5.5, 60, 5.5, 5.5))
          save_both_adv(
            p_for,
            file.path(sum_dir,
                      "ALL_cancers_advanced_stage_logistic_OR_forest"),
            9, max(5, 0.28 * nrow(logistic_adv) + 2.2)
          )
        }
      }

      make_facet_adv <- function(lst, xvar, colors, title, filename) {
        if (length(lst) == 0) return(invisible(NULL))
        dd <- bind_rows(lst)
        if (nrow(dd) == 0) return(invisible(NULL))
        ann <- dd %>%
          group_by(cancer_adv) %>%
          summarise(
            p = tryCatch(kruskal.test(score_adv,
                                      as.factor(.data[[xvar]]))$p.value,
                         error = function(e) NA_real_),
            ymax = max(score_adv, na.rm = TRUE),
            ymin = min(score_adv, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(
            x = (vapply(cancer_adv, function(cc) {
              length(unique(na.omit(dd[[xvar]][dd$cancer_adv == cc])))
            }, numeric(1)) + 1) / 2,
            y = ymax + 0.12 * ifelse(ymax > ymin, ymax - ymin, 1),
            label = format_p_step(p)
          )
        p <- ggplot(dd, aes(x = .data[[xvar]], y = score_adv,
                            fill = .data[[xvar]], color = .data[[xvar]])) +
          geom_boxplot(outlier.shape = NA, width = 0.58, alpha = 0.82,
                       linewidth = 0.25) +
          geom_jitter(width = 0.10, size = 0.28, alpha = 0.25,
                      show.legend = FALSE) +
          geom_text(data = ann, aes(x = x, y = y, label = label),
                    inherit.aes = FALSE, size = 2.0, color = "black") +
          facet_wrap(~ cancer_adv, scales = "free_y", ncol = 5) +
          scale_fill_manual(values = colors, drop = FALSE) +
          scale_color_manual(values = colors, drop = FALSE) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
          labs(title = wrap_title_step(title, 85), x = NULL,
               y = score_col) +
          theme_step(base_size = 8, title_size = 9.5,
                      axis_title_size = 7.5, axis_text_size = 6.2,
                      strip_size = 6.8) +
          theme(legend.position = "top",
                axis.text.x = element_text(angle = 45, hjust = 1))
        nca <- n_distinct(dd$cancer_adv)
        save_both_adv(
          p, file.path(sum_dir, filename), 13,
          max(6, ceiling(nca / 5) * 2.6)
        )
      }

      make_facet_adv(facet_data$stage, "stage_adv", stage_colors,
                     "All cancers: CAF3-corrected STEP_UP across overall stage",
                     "ALL_cancers_score_by_overall_stage_faceted")
      make_facet_adv(facet_data$stage_group, "stage_group_adv",
                     c("Early" = STEP_LOW_COLOR,
                       "Advanced" = STEP_HIGH_COLOR),
                     "All cancers: CAF3-corrected STEP_UP in early versus advanced stage",
                     "ALL_cancers_score_early_vs_advanced_faceted")
      make_facet_adv(facet_data$T, "T_adv", t_colors,
                     "All cancers: CAF3-corrected STEP_UP across T stage",
                     "ALL_cancers_score_by_T_stage_faceted")
      make_facet_adv(facet_data$N, "N_adv", n_colors,
                     "All cancers: CAF3-corrected STEP_UP across N stage",
                     "ALL_cancers_score_by_N_stage_faceted")
      make_facet_adv(facet_data$M, "M_adv", m_colors,
                     "All cancers: CAF3-corrected STEP_UP across M stage",
                     "ALL_cancers_score_by_M_stage_faceted")
      make_facet_adv(facet_data$Grade, "Grade_adv", grade_colors,
                     "All cancers: CAF3-corrected STEP_UP across histologic grade",
                     "ALL_cancers_score_by_grade_faceted")
      if (length(facet_data$Gleason) > 0) {
        gleason_all <- bind_rows(facet_data$Gleason)
        gleason_all$Gleason_factor <- factor(
          gleason_all$Gleason_adv,
          levels = sort(unique(gleason_all$Gleason_adv))
        )
        gleason_levels <- levels(gleason_all$Gleason_factor)
        gleason_cols <- setNames(
          grDevices::colorRampPalette(c("#A6CEE3", "#CB181D"))(
            length(gleason_levels)
          ),
          gleason_levels
        )
        make_facet_adv(
          split(gleason_all, gleason_all$cancer_adv),
          "Gleason_factor", gleason_cols,
          "CAF3-corrected STEP_UP across Gleason score",
          "ALL_cancers_score_by_Gleason_faceted"
        )
      }

      write.csv(
        data.frame(
          variable = c("cancer", "score", "stage_std", "stage_num",
                       "stage_group", "T", "N", "M", "grade", "gleason",
                       "primary_tumor"),
          detected_column = c(cancer_col, score_col, stage_std_col,
                              stage_num_col, stage_group_col, t_col, n_col,
                              m_col, grade_col, gleason_col, primary_col),
          stringsAsFactors = FALSE
        ),
        file.path(tab_dir, "advanced_clinical_detected_columns.csv"),
        row.names = FALSE
      )

      message("Advanced clinical additions finished: ", out_root)
    })
  })
}

if (isTRUE(RUN_ESTIMATE_EXPANSION)) {
  run_optional_module_v11("13/14 ESTIMATE expansion", {

    local({
      message("========== 13/14 ESTIMATE expansion: unified module ==========")
      result_dir <- MASTER_RESULT_DIR
      input_file <- file.path(result_dir, "09_pan_cancer_summary",
                              "ALL_cancers_scores_purity_clinical_merged.csv")
      if (!file.exists(input_file)) {
        warning("ESTIMATE expansion skipped: missing ", input_file)
      } else {
        dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
        pick_col_v11 <- function(candidates, nms) {
          hit <- candidates[candidates %in% nms]
          if (length(hit) > 0) return(hit[1])
          for (pat in candidates) {
            hit2 <- grep(pat, nms, value = TRUE, ignore.case = TRUE)
            if (length(hit2) > 0) return(hit2[1])
          }
          NA_character_
        }
        cancer_col <- pick_col_v11(c("cancer_loop", "cancer_from_file",
                                     "cancer_type", "cancer", "CANCER"),
                                   colnames(dat))
        score_col <- pick_col_v11(c("STEP_UP_residual_CAF3",
                                    "STEP.*UP.*residual.*CAF3"),
                                  colnames(dat))
        est_map <- c(
          "ESTIMATE score" = pick_col_v11(c("ESTIMATE_ESTIMATEScore",
                                            "^ESTIMATEScore$", "^ESTIMATE_score$"),
                                          colnames(dat)),
          "Immune score" = pick_col_v11(c("ESTIMATE_ImmuneScore",
                                          "^ImmuneScore$", "^Immune_score$"),
                                        colnames(dat)),
          "Stromal score" = pick_col_v11(c("ESTIMATE_StromalScore",
                                           "^StromalScore$", "^Stromal_score$"),
                                         colnames(dat))
        )
        est_map <- est_map[!is.na(est_map)]

        if (is.na(cancer_col) || is.na(score_col) || length(est_map) == 0) {
          warning("ESTIMATE expansion skipped: required columns were not detected.")
        } else {
          out13 <- file.path(result_dir, "13_TME_correlations")
          out14 <- file.path(result_dir, "14_TME_high_low_plots")
          d13_pc <- file.path(out13, "ESTIMATE_boxplots_per_cancer")
          d13_sum <- file.path(out13, "ESTIMATE_boxplots_all_cancers")
          d14_pc <- file.path(out14, "ESTIMATE_boxplots_per_cancer")
          d14_sum <- file.path(out14, "ESTIMATE_boxplots_all_cancers")
          for (d in c(out13, out14, d13_pc, d13_sum, d14_pc, d14_sum)) {
            dir.create(d, recursive = TRUE, showWarnings = FALSE)
          }

          dat$cancer_v11 <- toupper(as.character(dat[[cancer_col]]))
          dat$score_v11 <- suppressWarnings(as.numeric(dat[[score_col]]))
          for (nm in names(est_map)) {
            dat[[nm]] <- suppressWarnings(as.numeric(dat[[est_map[[nm]]]]))
          }
          dat <- dat %>%
            filter(!is.na(cancer_v11), cancer_v11 != "", is.finite(score_v11)) %>%
            group_by(cancer_v11) %>%
            mutate(cutoff_v11 = median(score_v11, na.rm = TRUE),
                   STEP_group = factor(ifelse(score_v11 > cutoff_v11,
                                               "High", "Low"),
                                        levels = c("Low", "High"))) %>%
            ungroup()

          long <- dat %>%
            select(cancer_v11, score_v11, STEP_group,
                   all_of(names(est_map))) %>%
            pivot_longer(cols = all_of(names(est_map)),
                         names_to = "feature", values_to = "value") %>%
            filter(is.finite(value))

          save_plot_est <- function(p, base, width, height) {
            print(p)
            ggsave(paste0(base, ".pdf"), p, width = width, height = height,
                   useDingbats = FALSE, limitsize = FALSE)
            ggsave(paste0(base, ".png"), p, width = width, height = height,
                   dpi = 300, limitsize = FALSE)
          }

          highlow_stats <- long %>%
            group_by(cancer_v11, feature) %>%
            group_modify(~ {
              d <- .x
              nl <- sum(d$STEP_group == "Low")
              nh <- sum(d$STEP_group == "High")
              pv <- if (nl >= 3 && nh >= 3) {
                tryCatch(wilcox.test(value ~ STEP_group, data = d,
                                     exact = FALSE)$p.value,
                         error = function(e) NA_real_)
              } else NA_real_
              tibble(
                n_low = nl, n_high = nh,
                median_low = if (nl > 0) median(d$value[d$STEP_group == "Low"],
                                                na.rm = TRUE) else NA_real_,
                median_high = if (nh > 0) median(d$value[d$STEP_group == "High"],
                                                 na.rm = TRUE) else NA_real_,
                p = pv
              )
            }) %>%
            ungroup() %>%
            mutate(
              median_diff_high_minus_low = median_high - median_low,
              FDR_by_feature = ave(p, feature,
                                   FUN = function(x) p.adjust(x, method = "BH")),
              p_label = vapply(p, format_p_short_step, character(1))
            )
          write.csv(highlow_stats,
                    file.path(out14,
                      "ALL_cancers_ESTIMATE_high_low_boxplot_statistics.csv"),
                    row.names = FALSE)
          write.csv(highlow_stats,
                    file.path(out13,
                      "ALL_cancers_ESTIMATE_high_low_boxplot_statistics.csv"),
                    row.names = FALSE)

          cor_stats <- long %>%
            group_by(cancer_v11, feature) %>%
            group_modify(~ {
              d <- .x
              ok <- is.finite(d$score_v11) & is.finite(d$value)
              if (sum(ok) < 10 || length(unique(d$score_v11[ok])) < 2 ||
                  length(unique(d$value[ok])) < 2) {
                return(tibble(n = sum(ok), rho = NA_real_, p = NA_real_))
              }
              ct <- tryCatch(cor.test(d$score_v11[ok], d$value[ok],
                                      method = "spearman", exact = FALSE),
                             error = function(e) NULL)
              tibble(n = sum(ok),
                     rho = if (is.null(ct)) NA_real_ else unname(ct$estimate),
                     p = if (is.null(ct)) NA_real_ else ct$p.value)
            }) %>%
            ungroup() %>%
            mutate(
              FDR_by_feature = ave(p, feature,
                                   FUN = function(x) p.adjust(x, method = "BH")),
              p_label = vapply(p, format_p_short_step, character(1))
            )
          write.csv(cor_stats,
                    file.path(out13,
                      "ALL_cancers_ESTIMATE_Spearman_correlation_statistics.csv"),
                    row.names = FALSE)

          for (ca in sort(unique(long$cancer_v11))) {
            dca <- long %>% filter(cancer_v11 == ca)
            sca <- highlow_stats %>% filter(cancer_v11 == ca)
            if (nrow(dca) == 0) next
            ypos <- dca %>% group_by(feature) %>%
              summarise(ymax = max(value, na.rm = TRUE),
                        ymin = min(value, na.rm = TRUE), .groups = "drop") %>%
              mutate(y = ymax + 0.12 * ifelse(ymax > ymin, ymax - ymin, 1),
                     x = 1.5)
            ann <- left_join(sca, ypos, by = "feature") %>%
              mutate(label = paste0("p=", p_label))
            p <- ggplot(dca, aes(x = STEP_group, y = value,
                                 fill = STEP_group, color = STEP_group)) +
              geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82,
                           linewidth = 0.35) +
              geom_jitter(width = 0.13, size = 0.65, alpha = 0.38,
                          show.legend = FALSE) +
              geom_text(data = ann, aes(x = x, y = y, label = label),
                        inherit.aes = FALSE, color = "black", size = 3) +
              facet_wrap(~ feature, scales = "free_y", nrow = 1) +
              scale_fill_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
              scale_color_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
              scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
              labs(title = wrap_title_step(
                     paste0(ca, ": ESTIMATE scores by CAF3-corrected STEP_UP group"),
                     65),
                   x = "CAF3-corrected STEP_UP group",
                   y = "ESTIMATE score") +
              theme_step(base_size = 10, title_size = 10,
                          axis_title_size = 9, axis_text_size = 8,
                          strip_size = 8) +
              theme(legend.position = "none")
            for (dd in c(d13_pc, d14_pc)) {
              save_plot_est(p,
                file.path(dd, paste0(ca,
                  "_ESTIMATE_Immune_Stromal_high_low_boxplots_with_p")),
                10.5, 4.8)
            }
          }

          for (feat in unique(long$feature)) {
            d <- long %>% filter(feature == feat)
            s <- highlow_stats %>% filter(feature == feat)
            ypos <- d %>% group_by(cancer_v11) %>%
              summarise(ymax = max(value, na.rm = TRUE),
                        ymin = min(value, na.rm = TRUE), .groups = "drop") %>%
              mutate(y = ymax + 0.13 * ifelse(ymax > ymin, ymax - ymin, 1),
                     x = 1.5)
            ann <- left_join(s, ypos, by = "cancer_v11") %>%
              mutate(label = paste0("p=", p_label))
            p <- ggplot(d, aes(x = STEP_group, y = value,
                               fill = STEP_group, color = STEP_group)) +
              geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82,
                           linewidth = 0.25) +
              geom_jitter(width = 0.12, size = 0.35, alpha = 0.25,
                          show.legend = FALSE) +
              geom_text(data = ann, aes(x = x, y = y, label = label),
                        inherit.aes = FALSE, color = "black", size = 2.2) +
              facet_wrap(~ cancer_v11, scales = "free_y", ncol = 5) +
              scale_fill_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
              scale_color_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
              scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
              labs(title = wrap_title_step(
                     paste0("All cancers: ", feat,
                            " by CAF3-corrected STEP_UP group"), 80),
                   x = "CAF3-corrected STEP_UP group", y = feat) +
              theme_step(base_size = 8.5, title_size = 10,
                          axis_title_size = 8, axis_text_size = 6.8,
                          strip_size = 7.2) +
              theme(legend.position = "none",
                    axis.text.x = element_text(angle = 35, hjust = 1))
            nca <- n_distinct(d$cancer_v11)
            safe_feat <- gsub("[^A-Za-z0-9]+", "_", feat)
            for (dd in c(d13_sum, d14_sum)) {
              save_plot_est(p,
                file.path(dd, paste0("ALL_cancers_", safe_feat,
                  "_high_low_boxplot_faceted_with_p")),
                13, max(6, ceiling(nca / 5) * 2.7))
            }
          }

          ypos_all <- long %>% group_by(cancer_v11, feature) %>%
            summarise(ymax = max(value, na.rm = TRUE),
                      ymin = min(value, na.rm = TRUE), .groups = "drop") %>%
            mutate(y = ymax + 0.13 * ifelse(ymax > ymin, ymax - ymin, 1),
                   x = 1.5)
          ann_all <- left_join(highlow_stats, ypos_all,
                               by = c("cancer_v11", "feature")) %>%
            mutate(label = paste0("p=", p_label))
          p_all <- ggplot(long, aes(x = STEP_group, y = value,
                                    fill = STEP_group, color = STEP_group)) +
            geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82,
                         linewidth = 0.22) +
            geom_jitter(width = 0.11, size = 0.28, alpha = 0.20,
                        show.legend = FALSE) +
            geom_text(data = ann_all, aes(x = x, y = y, label = label),
                      inherit.aes = FALSE, color = "black", size = 1.8) +
            facet_grid(feature ~ cancer_v11, scales = "free_y") +
            scale_fill_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
            scale_color_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
            labs(title = "All cancers: ESTIMATE / Immune / Stromal scores by CAF3-corrected STEP_UP group",
                 x = "CAF3-corrected STEP_UP group", y = NULL) +
            theme_step(base_size = 7.5, title_size = 10,
                        axis_title_size = 8, axis_text_size = 5.8,
                        strip_size = 6.5) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 60, hjust = 1))
          nca_all <- n_distinct(long$cancer_v11)
          for (dd in c(d13_sum, d14_sum)) {
            save_plot_est(p_all,
              file.path(dd,
                "ALL_cancers_ALL_ESTIMATE_metrics_high_low_boxplot_grid_with_p"),
              max(14, nca_all * 0.62), 8.5)
          }

          global_effect_p <- highlow_stats %>%
            group_by(feature) %>%
            summarise(
              global_p = if (sum(is.finite(median_diff_high_minus_low)) >= 3)
                tryCatch(wilcox.test(median_diff_high_minus_low,
                                     mu = 0, exact = FALSE)$p.value,
                         error = function(e) NA_real_) else NA_real_,
              y = max(median_diff_high_minus_low, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            mutate(label = vapply(global_p, format_p_step, character(1),
                                  prefix = "Across-cancer p = "),
                   y = y + 0.12 * ifelse(is.finite(y) & y != 0, abs(y), 1))
          p_eff <- ggplot(highlow_stats,
                          aes(x = feature, y = median_diff_high_minus_low)) +
            geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
            geom_boxplot(width = 0.55, outlier.shape = NA, fill = "white",
                         color = "grey35") +
            geom_jitter(aes(fill = median_diff_high_minus_low > 0),
                        shape = 21, width = 0.12, size = 2.4,
                        alpha = 0.85, color = "black") +
            geom_text(data = global_effect_p,
                      aes(x = feature, y = y, label = label),
                      inherit.aes = FALSE, size = 3) +
            scale_fill_manual(values = c("FALSE" = STEP_LOW_COLOR,
                                         "TRUE" = STEP_HIGH_COLOR),
                              guide = "none") +
            scale_y_continuous(expand = expansion(mult = c(0.12, 0.20))) +
            labs(title = "ESTIMATE high-low effects across cancers",
                 subtitle = "Each point represents one TCGA cancer type",
                 x = NULL, y = "Median(High) - Median(Low)") +
            theme_step(base_size = 10, title_size = 10,
                        axis_title_size = 9, axis_text_size = 8) +
            theme(axis.text.x = element_text(angle = 20, hjust = 1))
          for (dd in c(d13_sum, d14_sum)) {
            save_plot_est(p_eff,
              file.path(dd,
                "ESTIMATE_high_low_effect_each_cancer_as_one_point_boxplot"),
              8.2, 5.6)
          }

          global_rho_p <- cor_stats %>%
            group_by(feature) %>%
            summarise(
              global_p = if (sum(is.finite(rho)) >= 3)
                tryCatch(wilcox.test(rho, mu = 0, exact = FALSE)$p.value,
                         error = function(e) NA_real_) else NA_real_,
              y = max(rho, na.rm = TRUE), .groups = "drop"
            ) %>%
            mutate(label = vapply(global_p, format_p_step, character(1),
                                  prefix = "Across-cancer p = "),
                   y = pmin(0.95, y + 0.12))
          p_rho <- ggplot(cor_stats, aes(x = feature, y = rho)) +
            geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
            geom_boxplot(width = 0.55, outlier.shape = NA, fill = "white",
                         color = "grey35") +
            geom_jitter(aes(fill = rho > 0), shape = 21, width = 0.12,
                        size = 2.4, alpha = 0.85, color = "black") +
            geom_text(data = global_rho_p,
                      aes(x = feature, y = y, label = label),
                      inherit.aes = FALSE, size = 3) +
            scale_fill_manual(values = c("FALSE" = STEP_LOW_COLOR,
                                         "TRUE" = STEP_HIGH_COLOR),
                              guide = "none") +
            coord_cartesian(ylim = c(-1, 1), clip = "off") +
            labs(title = "STEP_UP_residual_CAF3 vs ESTIMATE correlations across cancers",
                 subtitle = "Each point represents one TCGA cancer type",
                 x = NULL, y = "Spearman rho") +
            theme_step(base_size = 10, title_size = 10,
                        axis_title_size = 9, axis_text_size = 8) +
            theme(axis.text.x = element_text(angle = 20, hjust = 1))
          save_plot_est(p_rho,
            file.path(d13_sum,
              "ESTIMATE_correlation_each_cancer_as_one_point_boxplot"),
            8.2, 5.6)

          message("ESTIMATE expansion finished: ", out13, " and ", out14)
        }
      }
    })

  })
}

if (isTRUE(RUN_ALL_CANCER_KM_OVERVIEW)) {
  run_optional_module_v11("08 all-cancer KM overview", {

    local({
      message("========== 08_survival: all-cancer KM overview ==========")
      result_dir <- MASTER_RESULT_DIR
      input_file <- file.path(result_dir, "09_pan_cancer_summary",
                              "ALL_cancers_scores_purity_clinical_merged.csv")
      if (!file.exists(input_file)) {
        warning("All-cancer KM overview skipped: missing ", input_file)
      } else {
        sdat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
        find_col_surv <- function(candidates, nms) {
          hit <- candidates[candidates %in% nms]
          if (length(hit) > 0) return(hit[1])
          NA_character_
        }
        cancer_col <- find_col_surv(c("cancer_loop", "cancer_from_file",
                                      "cancer_type", "cancer", "CANCER"),
                                    colnames(sdat))
        if (is.na(cancer_col)) {
          warning("All-cancer KM overview skipped: cancer column missing.")
        } else {
          out08_over <- file.path(result_dir, "08_survival",
                                  "00_ALL_CANCERS_OVERVIEW")
          dir.create(out08_over, recursive = TRUE, showWarnings = FALSE)

          endpoint_defs <- list(
            OS = c(find_col_surv(c("OS_MONTHS", "OS_MONTH", "OS.time", "OS_TIME"),
                                 colnames(sdat)),
                   find_col_surv(c("OS_EVENT", "OS_STATUS", "OS", "OS.event"),
                                 colnames(sdat))),
            DFS = c(find_col_surv(c("DFS_MONTHS", "DFS_MONTH", "DFS.time", "DFS_TIME"),
                                  colnames(sdat)),
                    find_col_surv(c("DFS_EVENT", "DFS_STATUS", "DFS", "DFS.event"),
                                  colnames(sdat))),
            PFS = c(find_col_surv(c("PFS_MONTHS", "PFS_MONTH", "PFS.time", "PFS_TIME"),
                                  colnames(sdat)),
                    find_col_surv(c("PFS_EVENT", "PFS_STATUS", "PFS", "PFS.event"),
                                  colnames(sdat))),
            RFS = c(find_col_surv(c("RFS_MONTHS", "RFS_MONTH", "RFS.time", "RFS_TIME"),
                                  colnames(sdat)),
                    find_col_surv(c("RFS_EVENT", "RFS_STATUS", "RFS", "RFS.event"),
                                  colnames(sdat)))
          )
          endpoint_defs <- endpoint_defs[vapply(endpoint_defs,
            function(z) length(z) == 2 && all(!is.na(z)), logical(1))]

          overview_scores <- intersect(c(
            "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
            "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN",
            "STEP_UP_CONSERVATIVE_CAF_ECM_FILTERED",
            "STEP_UP_TME_OVERLAP_FILTERED",
            "STEP_UP_MINUS_DOWN_RAW",
            "STEP_UP_residual_CAF3",
            "STEP_DOWN_residual_CAF3",
            "STEP_BALANCE_residual_CAF3",
            "STEP_UP_residual_CAF4",
            "STEP_DOWN_residual_CAF4",
            "STEP_BALANCE_residual_CAF4",
            "HUMAN_EPITHELIAL_SCORE_14",
            "HUMAN_MESENCHYMAL_SCORE_18",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "CAF_3gene_score", "CAF_4gene_index",
            "ECM_score", "Fibroblast_score", "TGFb_score"
          ), colnames(sdat))

          make_facet_km <- function(df, score_col, endpoint, time_col, event_col) {
            curve_list <- list()
            p_list <- list()
            qc_list <- list()
            ii <- 1
            for (ca in sort(unique(toupper(as.character(df[[cancer_col]]))))) {
              d <- df[toupper(as.character(df[[cancer_col]])) == ca, , drop = FALSE]
              d$.time <- suppressWarnings(as.numeric(d[[time_col]]))
              d$.event <- normalize_event_vector(d[[event_col]])
              d$.score <- suppressWarnings(as.numeric(d[[score_col]]))
              d <- d[complete.cases(d[, c(".time", ".event", ".score")]) &
                       is.finite(d$.time) & d$.time > 0 &
                       d$.event %in% c(0, 1) & is.finite(d$.score), , drop = FALSE]
              if (nrow(d) < 20 || length(unique(d$.event)) < 2 ||
                  length(unique(d$.score)) < 2) next
              med <- median(d$.score, na.rm = TRUE)
              d$group <- factor(ifelse(d$.score > med, "High", "Low"),
                                levels = c("Low", "High"))
              if (length(unique(d$group)) < 2 ||
                  min(table(d$group)) < 5) next

              fit <- tryCatch(survfit(Surv(.time, .event) ~ group, data = d),
                              error = function(e) NULL)
              if (is.null(fit)) next
              ss <- summary(fit)
              if (length(ss$time) == 0) next
              strata_chr <- sub("^group=", "", as.character(ss$strata))
              curve <- data.frame(
                cancer = ca,
                time = ss$time,
                surv = ss$surv,
                group = factor(strata_chr, levels = c("Low", "High")),
                stringsAsFactors = FALSE
              )
              starts <- data.frame(cancer = ca, time = 0, surv = 1,
                                   group = factor(c("Low", "High"),
                                                  levels = c("Low", "High")))
              curve_list[[ii]] <- bind_rows(starts, curve)

              lr <- tryCatch(survdiff(Surv(.time, .event) ~ group, data = d),
                             error = function(e) NULL)
              lp <- if (is.null(lr)) NA_real_ else
                1 - pchisq(lr$chisq, df = length(lr$n) - 1)
              p_list[[ii]] <- data.frame(
                cancer = ca,
                x = max(d$.time, na.rm = TRUE) * 0.04,
                y = 0.08,
                label = format_p_step(lp, "p = "),
                stringsAsFactors = FALSE
              )
              qc_list[[ii]] <- data.frame(
                cancer = ca, endpoint = endpoint, score = score_col,
                n = nrow(d), n_event = sum(d$.event == 1),
                n_low = sum(d$group == "Low"),
                n_high = sum(d$group == "High"),
                cutoff = med, logrank_p = lp,
                stringsAsFactors = FALSE
              )
              ii <- ii + 1
            }
            if (length(curve_list) < 3) return(NULL)
            curves <- bind_rows(curve_list)
            anns <- bind_rows(p_list)
            qc <- bind_rows(qc_list)
            list(curves = curves, anns = anns, qc = qc)
          }

          all_qc_over <- list()
          qi <- 1
          for (ep in names(endpoint_defs)) {
            time_col <- endpoint_defs[[ep]][1]
            event_col <- endpoint_defs[[ep]][2]
            for (sc in overview_scores) {
              obj <- make_facet_km(sdat, sc, ep, time_col, event_col)
              if (is.null(obj)) next
              nca <- n_distinct(obj$curves$cancer)
              p <- ggplot(obj$curves,
                          aes(x = time, y = surv, color = group)) +
                geom_step(linewidth = 0.62, na.rm = TRUE) +
                geom_text(data = obj$anns,
                          aes(x = x, y = y, label = label),
                          inherit.aes = FALSE, size = 2.2, hjust = 0,
                          color = "black") +
                facet_wrap(~ cancer, scales = "free_x", ncol = 5) +
                scale_color_manual(values = STEP_GROUP_COLORS,
                                   drop = FALSE) +
                scale_y_continuous(limits = c(0, 1),
                                   breaks = c(0, 0.25, 0.5, 0.75, 1)) +
                labs(
                  title = wrap_title_step(
                    paste0("All cancers: ", ep, " by ",
                           compact_score_label_step(sc)), 85),
                  subtitle = "High/Low cutoff is the within-cancer median; p values are log-rank tests",
                  x = "Months", y = paste0(ep, " probability"),
                  color = "Group"
                ) +
                theme_step(base_size = 8.5, title_size = 10,
                            axis_title_size = 8, axis_text_size = 6.8,
                            strip_size = 7.2) +
                theme(legend.position = "top")
              safe_sc <- gsub("[^A-Za-z0-9]+", "_", sc)
              base_out <- file.path(out08_over,
                paste0("ALL_cancers_", ep, "_", safe_sc,
                       "_KM_faceted_with_logrank_p"))
              ggsave(paste0(base_out, ".pdf"), p, width = 13,
                     height = max(6, ceiling(nca / 5) * 2.6),
                     limitsize = FALSE, useDingbats = FALSE)
              ggsave(paste0(base_out, ".png"), p, width = 13,
                     height = max(6, ceiling(nca / 5) * 2.6),
                     dpi = 300, limitsize = FALSE)
              print(p)
              all_qc_over[[qi]] <- obj$qc
              qi <- qi + 1
            }
          }
          if (length(all_qc_over) > 0) {
            write.csv(bind_rows(all_qc_over),
                      file.path(out08_over,
                        "ALL_cancers_KM_overview_sample_size_and_logrank_stats.csv"),
                      row.names = FALSE)
          }
          message("All-cancer KM overview finished: ", out08_over)
        }
      }
    })

  })
}

if (isTRUE(RUN_SURVIVAL_VOLCANO)) {
  run_optional_module_v11("16 survival volcano", {

    local({

      suppressPackageStartupMessages({
        library(ggplot2)
        library(dplyr)
        library(readr)
      })

      msg <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)

      find_root_dir <- function(start_dir = getwd()) {
        d <- normalizePath(start_dir, mustWork = TRUE)
        for (i in 1:8) {
          if (dir.exists(file.path(d, "TCGA_201906cBioP")) ||
              dir.exists(file.path(d, "gene set")) ||
              length(list.files(d, pattern = "^TCGA_STEP_full_results_", full.names = TRUE)) > 0) {
            return(d)
          }
          parent <- dirname(d)
          if (identical(parent, d)) break
          d <- parent
        }
        stop("Cannot locate TCGA_TME_correlation root from: ", start_dir,
             "\nPlease setwd('<project_root>') and rerun.")
      }

      root_dir <- find_root_dir(getwd())
      msg("Using root_dir: ", root_dir)

      res_dirs <- list.files(root_dir, pattern = "^TCGA_STEP_full_results_", full.names = TRUE)
      if (length(res_dirs) == 0) stop("Cannot find TCGA_STEP_full_results_* under: ", root_dir)
      res_dirs <- res_dirs[dir.exists(res_dirs)]
      res_dir <- MASTER_RESULT_DIR
      msg("Using result folder: ", res_dir)

      sum_dir <- file.path(res_dir, "09_pan_cancer_summary")
      if (!dir.exists(sum_dir)) stop("Cannot find 09_pan_cancer_summary under result folder: ", res_dir)

      out_dir <- file.path(res_dir, "16_survival_volcano_CAF3_OS_DFS")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      all_endpoint_file <- file.path(sum_dir, "ALL_cancers_all_endpoints_cox_summary.csv")
      unadj_file <- file.path(sum_dir, "ALL_cancers_all_endpoints_cox_unadjusted_score_only.csv")
      os_file <- file.path(sum_dir, "ALL_cancers_OS_cox_summary.csv")

      if (!file.exists(all_endpoint_file) && !file.exists(unadj_file) && !file.exists(os_file)) {
        stop("Cannot find Cox summary tables in: ", sum_dir,
             "\nExpected at least one of:\n",
             all_endpoint_file, "\n", unadj_file, "\n", os_file)
      }

      read_csv_quiet <- function(f) {
        if (!file.exists(f)) return(NULL)
        msg("Reading: ", f)
        suppressMessages(readr::read_csv(f, show_col_types = FALSE))
      }

      cox_all <- read_csv_quiet(all_endpoint_file)
      cox_unadj <- read_csv_quiet(unadj_file)
      cox_os <- read_csv_quiet(os_file)

      clean_cox <- function(dat) {
        if (is.null(dat) || nrow(dat) == 0) return(NULL)
        nm <- names(dat)
        
        if (!"endpoint" %in% nm) dat$endpoint <- "OS"
        
        if (!"cancer" %in% names(dat)) {
          cand <- intersect(c("cancer_type", "CANCER", "cancer_loop"), names(dat))
          if (length(cand) > 0) names(dat)[match(cand[1], names(dat))] <- "cancer"
        }
        if (!"HR" %in% names(dat)) {
          cand <- intersect(c("hazard_ratio", "hr"), names(dat))
          if (length(cand) > 0) names(dat)[match(cand[1], names(dat))] <- "HR"
        }
        if (!"p" %in% names(dat)) {
          cand <- intersect(c("p.value", "pval", "p_value", "P"), names(dat))
          if (length(cand) > 0) names(dat)[match(cand[1], names(dat))] <- "p"
        }
        dat
      }

      cox_all <- clean_cox(cox_all)
      cox_unadj <- clean_cox(cox_unadj)
      cox_os <- clean_cox(cox_os)

      select_score_rows <- function(dat, mode = c("univariate", "multivariate")) {
        mode <- match.arg(mode)
        if (is.null(dat) || nrow(dat) == 0) return(NULL)
        if (!all(c("cancer", "endpoint", "HR", "p") %in% names(dat))) return(NULL)
        dat <- dat %>% mutate(endpoint = as.character(endpoint), cancer = as.character(cancer))
        if ("score" %in% names(dat)) {
          dat <- dat %>% filter(score == "STEP_UP_residual_CAF3" | grepl("STEP.*UP.*residual.*CAF3", score, ignore.case = TRUE))
        }
        if ("model" %in% names(dat)) {
          if (mode == "univariate") {
            dat <- dat %>% filter(model %in% c("univariate_score", "score_only", "unadjusted_score_only") |
                                    grepl("univariate|unadjusted|score_only", model, ignore.case = TRUE))
          } else {
            
            dat <- dat %>% filter(model %in% c("score_scaled", "multivariate_score", "score") |
                                    grepl("score_scaled|multivariate.*score", model, ignore.case = TRUE))
          }
        } else if (mode == "multivariate") {
          return(NULL)
        }
        dat <- dat %>%
          mutate(HR = as.numeric(HR), p = as.numeric(p),
                 loge_HR = log(HR), neg_log10_p = -log10(p)) %>%
          filter(is.finite(HR), is.finite(p), p > 0, is.finite(loge_HR), is.finite(neg_log10_p))
        if (nrow(dat) == 0) return(NULL)
        dat
      }

      univ <- bind_rows(
        select_score_rows(cox_unadj, "univariate"),
        
        if (is.null(cox_unadj)) select_score_rows(cox_os, "univariate") else NULL
      ) %>% distinct(endpoint, cancer, .keep_all = TRUE)

      multi <- select_score_rows(cox_all, "multivariate")
      if (!is.null(multi)) multi <- multi %>% distinct(endpoint, cancer, .keep_all = TRUE)

      write_csv(univ, file.path(out_dir, "CAF3_STEP_UP_univariate_OS_DFS_available.csv"))
      if (!is.null(multi)) write_csv(multi, file.path(out_dir, "CAF3_STEP_UP_multivariate_OS_DFS_available.csv"))

      plot_volcano <- function(dat, endpoint_name, model_label, out_prefix) {
        dat <- dat %>% filter(endpoint == endpoint_name)
        if (nrow(dat) == 0) return(invisible(NULL))
        dat <- dat %>% mutate(
          sig = p < 0.05,
          p_text = vapply(p, format_p_short_step, character(1)),
          label = ifelse(sig, paste0(cancer, "\np=", p_text), ""),
          direction = ifelse(loge_HR >= 0, "HR > 1", "HR < 1")
        )
        p <- ggplot(dat, aes(x = loge_HR, y = neg_log10_p)) +
          geom_point(aes(color = sig), size = 2.7, alpha = 0.9) +
          geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
          ggrepel::geom_text_repel(aes(label = label), size = 3.2, max.overlaps = Inf, box.padding = 0.4, show.legend = FALSE) +
          scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#C74A37"), labels = c("P >= 0.05", "P < 0.05")) +
          labs(title = wrap_title_step(paste0("CAF3-corrected STEP_UP - ", endpoint_name, " Cox volcano"), 52),
               subtitle = wrap_title_step(model_label, 60),
               x = expression(log[e](HR)), y = expression(-log[10](P~value)), color = NULL) +
          theme_step(base_size = 10, title_size = 9.5,
                      axis_title_size = 9, axis_text_size = 8) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 9.5),
                plot.subtitle = element_text(size = 8.2))
        ggsave(paste0(out_prefix, ".pdf"), p, width = 5.2, height = 4.4)
        ggsave(paste0(out_prefix, ".png"), p, width = 5.2, height = 4.4, dpi = 300)
      }

      plot_bubble <- function(dat, model_label, out_prefix) {
        if (is.null(dat) || nrow(dat) == 0) return(invisible(NULL))
        dat <- dat %>% filter(endpoint %in% c("OS", "DFS"))
        if (nrow(dat) == 0) return(invisible(NULL))
        dat <- dat %>% mutate(
          direction = ifelse(loge_HR >= 0, "Risk (HR>1)", "Protective (HR<1)"),
          sig_size = pmin(neg_log10_p, 15),
          p_text = vapply(p, format_p_short_step, character(1)),
          endpoint = factor(endpoint, levels = c("OS", "DFS")),
          cancer = factor(cancer, levels = sort(unique(cancer)))
        )
        p <- ggplot(dat, aes(x = cancer, y = endpoint)) +
          geom_point(aes(size = sig_size, fill = neg_log10_p, color = direction), shape = 21, stroke = 0.7) +
          geom_text(aes(label = ifelse(p < 0.05, paste0("p=", p_text), "")),
                    size = 2.0, vjust = -1.1, color = "black", check_overlap = TRUE) +
          scale_size_continuous(range = c(1.2, 9), name = expression(-log[10](P))) +
          scale_fill_gradient(low = "#FFFFCC", high = "#D7301F", name = expression(-log[10](P))) +
          scale_color_manual(values = c("Risk (HR>1)" = "black", "Protective (HR<1)" = "#2C7FB8"), name = "Direction") +
          labs(title = "CAF3-corrected STEP_UP survival summary",
               subtitle = wrap_title_step(paste0(model_label, "; significant bubbles are labeled with raw p values"), 90),
               x = NULL, y = NULL) +
          theme_step(base_size = 9.5, title_size = 9.5,
                      axis_title_size = 8.5, axis_text_size = 7.3) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
                panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 9.5),
                plot.subtitle = element_text(size = 8))
        w <- max(10, length(unique(dat$cancer)) * 0.33)
        ggsave(paste0(out_prefix, ".pdf"), p, width = w, height = 3.3)
        ggsave(paste0(out_prefix, ".png"), p, width = w, height = 3.3, dpi = 300)
      }

      if (!is.null(univ) && nrow(univ) > 0) {
        for (ep in intersect(c("OS", "DFS"), unique(univ$endpoint))) {
          plot_volcano(univ, ep, "Univariate Cox: Surv(time, event) ~ STEP_UP_residual_CAF3",
                       file.path(out_dir, paste0("CAF3_STEP_UP_", ep, "_univariate_score_cox_volcano")))
        }
        plot_bubble(univ, "Univariate Cox: score only",
                    file.path(out_dir, "CAF3_STEP_UP_survival_endpoint_bubble_univariate_OS_DFS"))
      }

      if (!is.null(multi) && nrow(multi) > 0) {
        for (ep in intersect(c("OS", "DFS"), unique(multi$endpoint))) {
          plot_volcano(multi, ep, "Multivariate Cox: score term adjusted by available covariates",
                       file.path(out_dir, paste0("CAF3_STEP_UP_", ep, "_multivariate_score_scaled_cox_volcano")))
        }
        plot_bubble(multi, "Multivariate Cox: score term adjusted by available covariates",
                    file.path(out_dir, "CAF3_STEP_UP_survival_endpoint_bubble_multivariate_OS_DFS"))
      }

      combined <- bind_rows(
        if (!is.null(univ)) mutate(univ, plot_model = "univariate_score") else NULL,
        if (!is.null(multi)) mutate(multi, plot_model = "multivariate_score_scaled") else NULL
      )
      write_csv(combined, file.path(out_dir, "CAF3_STEP_UP_survival_volcano_combined_OS_DFS_tables.csv"))

      msg("Done. Outputs saved to: ", out_dir)

    })

  })
}

if (isTRUE(RUN_TIDE_MAIN)) {
  run_optional_module_v11("18 TIDE main", {

    local({

      options(stringsAsFactors = FALSE)

      suppressPackageStartupMessages({
        library(readr)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(stringr)
        library(tibble)
      })

      TIDE_GROUP_COLORS <- STEP_GROUP_COLORS
      TIDE_RESPONSE_COLORS <- c("Predicted non-responder" = "#D9D9D9",
                                "Predicted responder" = "#F4A6A6")

      ROOT <- MASTER_PROJECT_ROOT
      STEP_DIR <- MASTER_RESULT_DIR
      SUMMARY_DIR <- file.path(STEP_DIR, "09_pan_cancer_summary")
      STEP_FILE <- file.path(SUMMARY_DIR, "ALL_cancers_scores_purity_clinical_merged.csv")

      tide_batch_candidates <- c(
        INPUT_OVERRIDE_TIDE_BATCH_FILE,
        file.path(ROOT, "18_TIDEpy_results", "tables", "ALL_TCGA_TIDE_results_merged.csv"),
        file.path(ROOT, "18_TIDEpy_results", "tables", "ALL_TCGA_TIDE_results_merged.tsv"),
        file.path(ROOT, "ALL_TCGA_TIDE_results_merged.csv")
      )
      tide_batch_candidates_valid <- tide_batch_candidates[!is.na(tide_batch_candidates) & nzchar(tide_batch_candidates)]
      TIDE_BATCH_FILE <- tide_batch_candidates_valid[file.exists(tide_batch_candidates_valid)][1]
      if (length(TIDE_BATCH_FILE) == 0 || is.na(TIDE_BATCH_FILE)) TIDE_BATCH_FILE <- tide_batch_candidates_valid[1]
      TIDE_RESULT_DIR_0 <- if (is_existing_dir_v112(INPUT_OVERRIDE_TIDE_RESULT_DIR)) INPUT_OVERRIDE_TIDE_RESULT_DIR else NA_character_
      TIDE_RESULT_DIR_1 <- file.path(ROOT, "18_TIDEpy_results", "02_TIDE_results_by_cancer")
      TIDE_RESULT_DIR_2 <- file.path(ROOT, "CIBERSORT_upload_inputs", "linear_normalized_for_CIBERSORT")

      OUT_DIR <- file.path(STEP_DIR, "18_TIDE_CAF3_STEP_UP", "v8_fixed_columns")
      TABLE_DIR <- file.path(OUT_DIR, "tables")
      PANCAN_DIR <- file.path(OUT_DIR, "pancancer_plots")
      PER_CANCER_DIR <- file.path(OUT_DIR, "per_cancer_plots")
      RESPONDER_DIR <- file.path(OUT_DIR, "responder_plots")
      LOG_DIR <- file.path(OUT_DIR, "logs")

      for (d in c(OUT_DIR, TABLE_DIR, PANCAN_DIR, PER_CANCER_DIR, RESPONDER_DIR, LOG_DIR)) {
        dir.create(d, showWarnings = FALSE, recursive = TRUE)
      }

      log_file <- file.path(LOG_DIR, paste0("v8_fixed_columns_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
      log_msg <- function(...) {
        msg <- paste0(...)
        message(msg)
        cat(msg, "\n", file = log_file, append = TRUE)
      }

      log_msg("============================================================")
      log_msg("v8 fixed-column STEP/CAF3-TIDE analysis started")
      log_msg("ROOT: ", ROOT)
      log_msg("STEP_FILE: ", STEP_FILE)
      log_msg("TIDE_BATCH_FILE: ", TIDE_BATCH_FILE)
      log_msg("OUT_DIR: ", OUT_DIR)
      log_msg("============================================================")

      clean_sample_id <- function(x) {
        x <- as.character(x)
        x <- gsub("\\.", "-", x)
        x <- trimws(x)
        x
      }

      as_numeric_safely <- function(x) suppressWarnings(as.numeric(as.character(x)))

      parse_responder_bool <- function(x) {
        y <- tolower(trimws(as.character(x)))
        dplyr::case_when(
          y %in% c("true", "t", "1", "yes", "y", "responder", "response", "respond") ~ TRUE,
          y %in% c("false", "f", "0", "no", "n", "nonresponder", "non-responder", "no response", "no_response") ~ FALSE,
          TRUE ~ NA
        )
      }

      safe_wilcox <- function(dat, group_col = "STEP_UP_CAF3_group", value_col = "value") {
        dat <- dat %>%
          filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]]), .data[[group_col]] %in% c("Low", "High"))

        n_low <- sum(dat[[group_col]] == "Low", na.rm = TRUE)
        n_high <- sum(dat[[group_col]] == "High", na.rm = TRUE)
        med_low <- ifelse(n_low > 0, median(dat[[value_col]][dat[[group_col]] == "Low"], na.rm = TRUE), NA_real_)
        med_high <- ifelse(n_high > 0, median(dat[[value_col]][dat[[group_col]] == "High"], na.rm = TRUE), NA_real_)

        if (n_low < 2 || n_high < 2 || length(unique(dat[[value_col]])) < 2) {
          return(tibble(n_low = n_low, n_high = n_high,
                        median_low = med_low, median_high = med_high,
                        delta_high_low = med_high - med_low,
                        statistic = NA_real_, p_value = NA_real_))
        }

        wt <- tryCatch(
          wilcox.test(dat[[value_col]] ~ dat[[group_col]], exact = FALSE),
          error = function(e) NULL
        )

        if (is.null(wt)) {
          return(tibble(n_low = n_low, n_high = n_high,
                        median_low = med_low, median_high = med_high,
                        delta_high_low = med_high - med_low,
                        statistic = NA_real_, p_value = NA_real_))
        }

        tibble(n_low = n_low, n_high = n_high,
               median_low = med_low, median_high = med_high,
               delta_high_low = med_high - med_low,
               statistic = unname(wt$statistic), p_value = wt$p.value)
      }

      safe_spearman <- function(dat, x_col = "STEP_UP_residual_CAF3", y_col = "value") {
        dat <- dat %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
        n_pair <- nrow(dat)

        if (n_pair < 3 || length(unique(dat[[x_col]])) < 2 || length(unique(dat[[y_col]])) < 2) {
          return(tibble(n = n_pair, spearman_r = NA_real_, p_value = NA_real_))
        }

        ct <- tryCatch(
          suppressWarnings(cor.test(dat[[x_col]], dat[[y_col]], method = "spearman", exact = FALSE)),
          error = function(e) NULL
        )

        if (is.null(ct)) {
          return(tibble(n = n_pair, spearman_r = NA_real_, p_value = NA_real_))
        }

        tibble(n = n_pair, spearman_r = unname(ct$estimate), p_value = ct$p.value)
      }

      p_label <- function(p) {
        ifelse(is.na(p), "p = NA",
               ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 3))))
      }

      save_plot_both <- function(p, filename_base, out_dir, width = 7, height = 5) {
        pdf_file <- file.path(out_dir, paste0(filename_base, ".pdf"))
        png_file <- file.path(out_dir, paste0(filename_base, ".png"))
        print(p)  
        ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE)
        ggsave(png_file, p, width = width, height = height, dpi = 300)
        log_msg("Saved plot: ", pdf_file)
        log_msg("Saved plot: ", png_file)
      }

      if (!file.exists(STEP_FILE)) {
        stop("Cannot find STEP file: ", STEP_FILE)
      }

      step_raw <- readr::read_csv(STEP_FILE, show_col_types = FALSE)

      first_existing_step_col <- function(candidates) {
        hit <- candidates[candidates %in% colnames(step_raw)]
        if (length(hit) > 0) return(hit[1])
        NA_character_
      }
      step_sample_col <- first_existing_step_col(c(
        "sample15", "sample15_for_merge", "sample", "sample_raw",
        "sample_id", "Tumor_Sample_Barcode", "patient_key"
      ))
      step_cancer_col <- first_existing_step_col(c(
        "cancer_loop", "cancer_from_file", "cancer_type", "cancer", "CANCER"
      ))
      step_score_col <- first_existing_step_col(c(
        "STEP_UP_residual_CAF3", "STEP_UP_residual_CAF3_z"
      ))
      if (is.na(step_sample_col) || is.na(step_cancer_col) || is.na(step_score_col)) {
        readr::write_csv(
          tibble::tibble(column_name = colnames(step_raw)),
          file.path(TABLE_DIR, "DEBUG_STEP_columns_detected.csv")
        )
        stop(
          "STEP table column detection failed. Need sample, cancer and CAF3-corrected STEP columns. ",
          "Detected sample=", step_sample_col, "; cancer=", step_cancer_col,
          "; score=", step_score_col
        )
      }

      step <- step_raw %>%
        transmute(
          sample_id_15 = substr(clean_sample_id(.data[[step_sample_col]]), 1, 15),
          patient_id = if ("patient_id" %in% colnames(step_raw)) clean_sample_id(patient_id) else substr(clean_sample_id(.data[[step_sample_col]]), 1, 12),
          cancer_type = toupper(as.character(.data[[step_cancer_col]])),
          STEP_UP_residual_CAF3 = as_numeric_safely(.data[[step_score_col]]),
          CAF_3gene_score = if ("CAF_3gene_score" %in% colnames(step_raw)) as_numeric_safely(CAF_3gene_score) else NA_real_,
          STEP_UP_raw = if ("STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN" %in% colnames(step_raw)) as_numeric_safely(STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN) else NA_real_,
          STEP_DOWN_raw = if ("STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN" %in% colnames(step_raw)) as_numeric_safely(STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN) else NA_real_
        ) %>%
        filter(!is.na(sample_id_15), sample_id_15 != "", !is.na(cancer_type), cancer_type != "") %>%
        distinct(cancer_type, sample_id_15, .keep_all = TRUE) %>%
        group_by(cancer_type) %>%
        mutate(
          STEP_UP_CAF3_median = median(STEP_UP_residual_CAF3, na.rm = TRUE),
          STEP_UP_CAF3_group = case_when(
            is.na(STEP_UP_residual_CAF3) ~ NA_character_,
            STEP_UP_residual_CAF3 > STEP_UP_CAF3_median ~ "High",
            STEP_UP_residual_CAF3 <= STEP_UP_CAF3_median ~ "Low",
            TRUE ~ NA_character_
          ),
          STEP_UP_CAF3_group = factor(STEP_UP_CAF3_group, levels = c("Low", "High"))
        ) %>%
        ungroup()

      readr::write_csv(step, file.path(TABLE_DIR, "plot_v8_STEP_fixed_columns_used.csv"))
      readr::write_csv(
        step %>% count(cancer_type, STEP_UP_CAF3_group, name = "n"),
        file.path(TABLE_DIR, "plot_v8_STEP_CAF3_group_counts_by_cancer.csv")
      )

      log_msg("Read STEP rows: ", nrow(step))

      read_one_tide_txt <- function(path) {
        cancer <- basename(path) %>% str_replace("_TIDE_result\\.txt$", "")
        df <- readr::read_tsv(path, show_col_types = FALSE)

        if (colnames(df)[1] %in% c("", "...1", "X1")) {
          colnames(df)[1] <- "sample_id"
        }
        if (!"sample_id" %in% colnames(df)) {
          colnames(df)[1] <- "sample_id"
        }

        df %>%
          mutate(
            cancer_type = toupper(cancer),
            sample_id = clean_sample_id(sample_id),
            sample_id_15 = substr(sample_id, 1, 15),
            sample_id_12 = substr(sample_id, 1, 12),
            sample_id_16 = substr(sample_id, 1, 16)
          )
      }

      if (file.exists(TIDE_BATCH_FILE)) {
        log_msg("Reading batch TIDE table: ", TIDE_BATCH_FILE)
        tide_raw <- if (grepl("\\.tsv$", TIDE_BATCH_FILE, ignore.case = TRUE)) {
          readr::read_tsv(TIDE_BATCH_FILE, show_col_types = FALSE)
        } else {
          readr::read_csv(TIDE_BATCH_FILE, show_col_types = FALSE)
        }

        first_existing_tide_col <- function(candidates) {
          hit <- candidates[candidates %in% colnames(tide_raw)]
          if (length(hit) > 0) return(hit[1])
          NA_character_
        }
        tide_sample_col <- first_existing_tide_col(c(
          "sample_id", "sample_id_15", "Sample", "sample", "Mixture",
          "Tumor_Sample_Barcode", "...1", "X1"
        ))
        if (is.na(tide_sample_col)) tide_sample_col <- colnames(tide_raw)[1]
        tide_cancer_col <- first_existing_tide_col(c(
          "cancer_type", "cancer_loop", "cancer_from_file", "cancer", "CANCER"
        ))

        tide_raw$sample_id_v11 <- tide_raw[[tide_sample_col]]
        tide_raw$cancer_type_v11 <- if (!is.na(tide_cancer_col)) {
          tide_raw[[tide_cancer_col]]
        } else {
          NA_character_
        }

        tide <- tide_raw %>%
          mutate(
            sample_id = clean_sample_id(sample_id_v11),
            sample_id_15 = substr(sample_id, 1, 15),
            sample_id_12 = substr(sample_id, 1, 12),
            sample_id_16 = substr(sample_id, 1, 16),
            cancer_type = toupper(as.character(cancer_type_v11))
          )
      } else {
        log_msg("Batch TIDE table not found. Falling back to *_TIDE_result.txt files.")
        tide_files <- c(
          if (is_existing_dir_v112(TIDE_RESULT_DIR_0)) list.files(TIDE_RESULT_DIR_0, pattern = "_TIDE_result\\.txt$", full.names = TRUE) else character(0),
          list.files(TIDE_RESULT_DIR_1, pattern = "_TIDE_result\\.txt$", full.names = TRUE),
          list.files(TIDE_RESULT_DIR_2, pattern = "_TIDE_result\\.txt$", full.names = TRUE)
        ) %>% unique()

        if (length(tide_files) == 0) {
          stop("No TIDE results found. Expected: ", TIDE_BATCH_FILE, " or *_TIDE_result.txt")
        }
        tide <- bind_rows(lapply(tide_files, read_one_tide_txt))
      }

      find_metric_alias_v11 <- function(canonical, patterns) {
        if (canonical %in% colnames(tide)) return(invisible(TRUE))
        hit <- unique(unlist(lapply(patterns, function(pat) {
          grep(pat, colnames(tide), value = TRUE, ignore.case = TRUE)
        })))
        if (length(hit) > 0) colnames(tide)[match(hit[1], colnames(tide))] <<- canonical
        invisible(length(hit) > 0)
      }
      find_metric_alias_v11("TIDE", c("^TIDE$", "TIDE.*score"))
      find_metric_alias_v11("Dysfunction", c("^Dysfunction$", "dysfunction"))
      find_metric_alias_v11("Exclusion", c("^Exclusion$", "exclusion"))
      if (!"Responder" %in% colnames(tide)) {
        resp_hit <- grep("responder|predicted.*response|response", colnames(tide),
                         value = TRUE, ignore.case = TRUE)
        if (length(resp_hit) > 0) colnames(tide)[match(resp_hit[1], colnames(tide))] <- "Responder"
      }

      required_tide_metrics <- c("TIDE", "Dysfunction", "Exclusion")
      missing_tide_metrics <- setdiff(required_tide_metrics, colnames(tide))
      if (length(missing_tide_metrics) > 0) {
        readr::write_csv(
          tibble::tibble(column_name = colnames(tide)),
          file.path(TABLE_DIR, "DEBUG_TIDE_columns_detected.csv")
        )
        stop("TIDE table is missing required metric columns: ", paste(missing_tide_metrics, collapse = ", "))
      }

      tide <- tide %>%
        mutate(
          TIDE = as_numeric_safely(TIDE),
          Dysfunction = as_numeric_safely(Dysfunction),
          Exclusion = as_numeric_safely(Exclusion),
          Responder_bool = if ("Responder" %in% colnames(.)) parse_responder_bool(Responder) else NA
        ) %>%
        distinct(cancer_type, sample_id_15, .keep_all = TRUE)

      readr::write_csv(tide, file.path(TABLE_DIR, "plot_v8_ALL_TIDE_table_used.csv"))
      readr::write_csv(
        tide %>% count(cancer_type, name = "n_tide_samples"),
        file.path(TABLE_DIR, "plot_v8_TIDE_sample_counts_by_cancer.csv")
      )

      log_msg("Read TIDE rows: ", nrow(tide))

      merged <- step %>%
        inner_join(
          tide,
          by = c("cancer_type", "sample_id_15"),
          suffix = c(".STEP", ".TIDE")
        )

      if (nrow(merged) == 0) {
        
        log_msg("Cancer+sample merge returned 0 rows. Trying sample_id_15-only merge.")
        merged <- step %>%
          inner_join(tide, by = "sample_id_15", suffix = c(".STEP", ".TIDE")) %>%
          mutate(
            cancer_type = ifelse(!is.na(cancer_type.STEP), cancer_type.STEP, cancer_type.TIDE)
          )
      }

      if (nrow(merged) == 0) {
        stop("STEP and TIDE did not merge. Check sample_id_15 in both tables.")
      }

      readr::write_csv(merged, file.path(TABLE_DIR, "plot_v8_ALL_TCGA_TIDE_STEP_CAF3_merged_fixed_columns.csv"))

      merge_qc <- tibble(
        n_step = nrow(step),
        n_tide = nrow(tide),
        n_merged = nrow(merged),
        n_cancers_merged = n_distinct(merged$cancer_type)
      )
      readr::write_csv(merge_qc, file.path(TABLE_DIR, "plot_v8_merge_QC_summary.csv"))
      readr::write_csv(
        merged %>% count(cancer_type, STEP_UP_CAF3_group, name = "n"),
        file.path(TABLE_DIR, "plot_v8_group_counts_after_merge.csv")
      )

      log_msg("Merged rows: ", nrow(merged))
      log_msg("Merged cancers: ", n_distinct(merged$cancer_type))

      metrics <- c("TIDE", "Dysfunction", "Exclusion")

      plot_long <- merged %>%
        select(cancer_type, sample_id_15, STEP_UP_residual_CAF3, STEP_UP_CAF3_group,
               all_of(metrics), Responder_bool) %>%
        pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
        mutate(value = as_numeric_safely(value)) %>%
        filter(!is.na(value), !is.na(STEP_UP_CAF3_group))

      readr::write_csv(plot_long, file.path(TABLE_DIR, "plot_v8_boxplot_input_samples_long.csv"))
      readr::write_csv(
        plot_long %>% count(cancer_type, metric, STEP_UP_CAF3_group, name = "n"),
        file.path(TABLE_DIR, "plot_v8_input_QC_by_cancer_metric_group.csv")
      )

      if (nrow(plot_long) == 0) {
        readr::write_csv(
          tibble(reason = "plot_long has zero rows after filtering non-missing metrics and High/Low groups"),
          file.path(TABLE_DIR, "plot_v8_EMPTY_PLOT_DEBUG.csv")
        )
        stop("No valid rows for plotting. See plot_v8_EMPTY_PLOT_DEBUG.csv")
      }

      pancan_wilcox <- plot_long %>%
        group_by(metric) %>%
        group_modify(~ safe_wilcox(.x)) %>%
        ungroup() %>%
        mutate(
          p_adj_BH = p.adjust(p_value, method = "BH"),
          p_label = p_label(p_value),
          p_adj_label = p_label(p_adj_BH)
        )
      readr::write_csv(pancan_wilcox, file.path(TABLE_DIR, "plot_v8_pancancer_TIDE_metrics_wilcox_summary.csv"))

      per_cancer_wilcox <- plot_long %>%
        group_by(cancer_type, metric) %>%
        group_modify(~ safe_wilcox(.x)) %>%
        ungroup() %>%
        group_by(metric) %>%
        mutate(
          p_adj_BH = p.adjust(p_value, method = "BH"),
          p_label = p_label(p_value),
          p_adj_label = p_label(p_adj_BH)
        ) %>%
        ungroup()
      readr::write_csv(per_cancer_wilcox, file.path(TABLE_DIR, "plot_v8_per_cancer_TIDE_metrics_wilcox_summary.csv"))

      pancan_spearman <- plot_long %>%
        group_by(metric) %>%
        group_modify(~ safe_spearman(.x)) %>%
        ungroup() %>%
        mutate(
          p_adj_BH = p.adjust(p_value, method = "BH"),
          p_label = p_label(p_value),
          p_adj_label = p_label(p_adj_BH)
        )
      readr::write_csv(pancan_spearman, file.path(TABLE_DIR, "plot_v8_pancancer_TIDE_metrics_spearman_with_STEP_UP_CAF3_residual.csv"))

      per_cancer_spearman <- plot_long %>%
        group_by(cancer_type, metric) %>%
        group_modify(~ safe_spearman(.x)) %>%
        ungroup() %>%
        group_by(metric) %>%
        mutate(
          p_adj_BH = p.adjust(p_value, method = "BH"),
          p_label = p_label(p_value),
          p_adj_label = p_label(p_adj_BH)
        ) %>%
        ungroup()
      readr::write_csv(per_cancer_spearman, file.path(TABLE_DIR, "plot_v8_per_cancer_TIDE_metrics_spearman_with_STEP_UP_CAF3_residual.csv"))

      pancan_ann <- plot_long %>%
        group_by(metric) %>%
        summarise(ymax = max(value, na.rm = TRUE),
                  ymin = min(value, na.rm = TRUE), .groups = "drop") %>%
        left_join(pancan_wilcox %>% select(metric, p_label), by = "metric") %>%
        mutate(x = 1.5,
               y = ymax + 0.10 * ifelse(ymax > ymin, ymax - ymin, 1))
      p_fac <- ggplot(plot_long, aes(x = STEP_UP_CAF3_group, y = value, fill = STEP_UP_CAF3_group, color = STEP_UP_CAF3_group)) +
        geom_boxplot(width = 0.55, outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 0.8, alpha = 0.45, show.legend = FALSE) +
        scale_fill_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
        scale_color_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
        geom_text(data = pancan_ann, aes(x = x, y = y, label = p_label),
                  inherit.aes = FALSE, color = "black", size = 3) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        labs(
          x = "CAF3-corrected STEP_UP group",
          y = "TIDE-derived score",
          title = wrap_title_step("Pan-cancer TIDE metrics by CAF3-corrected STEP_UP group", 65)
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
              axis.text = element_text(color = "black"))
      save_plot_both(p_fac, "plot_v8_PanCancer_TIDE_Dysfunction_Exclusion_CAF3_STEP_UP_high_low_facet", PANCAN_DIR, width = 9, height = 4.5)

      for (m in metrics) {
        dat_m <- plot_long %>% filter(metric == m)
        stat_m <- pancan_wilcox %>% filter(metric == m)
        y_top <- max(dat_m$value, na.rm = TRUE)
        y_rng <- diff(range(dat_m$value, na.rm = TRUE))
        if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1
        label <- if (nrow(stat_m) > 0) stat_m$p_label[1] else "p = NA"

        p <- ggplot(dat_m, aes(x = STEP_UP_CAF3_group, y = value, fill = STEP_UP_CAF3_group, color = STEP_UP_CAF3_group)) +
          geom_boxplot(width = 0.55, outlier.shape = NA) +
          geom_jitter(width = 0.15, size = 1.0, alpha = 0.55, show.legend = FALSE) +
          scale_fill_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
          scale_color_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
          annotate("text", x = 1.5, y = y_top + 0.08 * y_rng, label = label, size = 4.2) +
          labs(
            x = "CAF3-corrected STEP_UP group",
            y = m,
            title = paste0("Pan-cancer ", m)
          ) +
          theme_classic(base_size = 12) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                axis.text = element_text(color = "black"))

        save_plot_both(p, paste0("plot_v8_PanCancer_", m, "_CAF3_STEP_UP_high_low_boxplot"), PANCAN_DIR, width = 5.2, height = 4.6)
      }

      for (m in metrics) {
        dat_m <- plot_long %>% filter(metric == m)
        if (nrow(dat_m) == 0) next
        ann_m <- dat_m %>%
          group_by(cancer_type) %>%
          summarise(ymax = max(value, na.rm = TRUE),
                    ymin = min(value, na.rm = TRUE), .groups = "drop") %>%
          left_join(per_cancer_wilcox %>%
                      filter(metric == m) %>%
                      select(cancer_type, p_label),
                    by = "cancer_type") %>%
          mutate(x = 1.5,
                 y = ymax + 0.12 * ifelse(ymax > ymin, ymax - ymin, 1))

        p <- ggplot(dat_m, aes(x = STEP_UP_CAF3_group, y = value, fill = STEP_UP_CAF3_group, color = STEP_UP_CAF3_group)) +
          geom_boxplot(width = 0.55, outlier.shape = NA) +
          geom_jitter(width = 0.12, size = 0.45, alpha = 0.4, show.legend = FALSE) +
          scale_fill_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
          scale_color_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
          geom_text(data = ann_m, aes(x = x, y = y, label = p_label),
                    inherit.aes = FALSE, color = "black", size = 2.1) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
          facet_wrap(~ cancer_type, scales = "free_y") +
          labs(
            x = "CAF3-corrected STEP_UP group",
            y = m,
            title = wrap_title_step(paste0(m, " by CAF3-corrected STEP_UP group across cancers"), 70)
          ) +
          theme_classic(base_size = 9) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                strip.text = element_text(face = "bold", size = 7),
                axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
                axis.text.y = element_text(color = "black"))

        save_plot_both(p, paste0("plot_v8_PerCancer_facet_", m, "_CAF3_STEP_UP_high_low_boxplot"), PER_CANCER_DIR, width = 11, height = 8)
      }

      for (cancer in sort(unique(plot_long$cancer_type))) {
        for (m in metrics) {
          dat_cm <- plot_long %>% filter(cancer_type == cancer, metric == m)
          stat_cm <- per_cancer_wilcox %>% filter(cancer_type == cancer, metric == m)
          if (nrow(dat_cm) == 0) next

          n_low <- sum(dat_cm$STEP_UP_CAF3_group == "Low", na.rm = TRUE)
          n_high <- sum(dat_cm$STEP_UP_CAF3_group == "High", na.rm = TRUE)
          if (n_low < 2 || n_high < 2) next

          y_top <- max(dat_cm$value, na.rm = TRUE)
          y_rng <- diff(range(dat_cm$value, na.rm = TRUE))
          if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1
          label <- if (nrow(stat_cm) > 0) stat_cm$p_label[1] else "p = NA"

          p <- ggplot(dat_cm, aes(x = STEP_UP_CAF3_group, y = value, fill = STEP_UP_CAF3_group, color = STEP_UP_CAF3_group)) +
            geom_boxplot(width = 0.55, outlier.shape = NA) +
            geom_jitter(width = 0.15, size = 1.0, alpha = 0.55, show.legend = FALSE) +
            scale_fill_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
            scale_color_manual(values = TIDE_GROUP_COLORS, drop = FALSE) +
            annotate("text", x = 1.5, y = y_top + 0.08 * y_rng, label = label, size = 4.0) +
            labs(
              x = "CAF3-corrected STEP_UP group",
              y = m,
              title = paste0(cancer, " - ", m)
            ) +
            theme_classic(base_size = 12) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                  axis.text = element_text(color = "black"))

          save_plot_both(p, paste0("plot_v8_", cancer, "_", m, "_CAF3_STEP_UP_high_low_boxplot"), PER_CANCER_DIR, width = 5.0, height = 4.5)
        }
      }

      responder_data <- merged %>%
        select(cancer_type, sample_id_15, STEP_UP_residual_CAF3, STEP_UP_CAF3_group, Responder_bool) %>%
        filter(!is.na(Responder_bool), !is.na(STEP_UP_CAF3_group)) %>%
        mutate(Responder_label = ifelse(Responder_bool, "Predicted responder", "Predicted non-responder"))

      readr::write_csv(responder_data, file.path(TABLE_DIR, "plot_v8_responder_plot_input_samples.csv"))

      safe_fisher <- function(dat) {
        dat <- dat %>% filter(!is.na(STEP_UP_CAF3_group), !is.na(Responder_bool))
        n_low <- sum(dat$STEP_UP_CAF3_group == "Low")
        n_high <- sum(dat$STEP_UP_CAF3_group == "High")
        resp_low <- sum(dat$STEP_UP_CAF3_group == "Low" & dat$Responder_bool)
        resp_high <- sum(dat$STEP_UP_CAF3_group == "High" & dat$Responder_bool)
        prop_low <- ifelse(n_low > 0, resp_low / n_low, NA_real_)
        prop_high <- ifelse(n_high > 0, resp_high / n_high, NA_real_)

        if (n_low == 0 || n_high == 0 || length(unique(dat$Responder_bool)) < 2) {
          return(tibble(n_low = n_low, n_high = n_high,
                        responder_low = resp_low, responder_high = resp_high,
                        responder_prop_low = prop_low, responder_prop_high = prop_high,
                        delta_prop_high_low = prop_high - prop_low,
                        p_value = NA_real_))
        }

        tab <- table(dat$STEP_UP_CAF3_group, dat$Responder_bool)
        ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
        p <- if (is.null(ft)) NA_real_ else ft$p.value

        tibble(n_low = n_low, n_high = n_high,
               responder_low = resp_low, responder_high = resp_high,
               responder_prop_low = prop_low, responder_prop_high = prop_high,
               delta_prop_high_low = prop_high - prop_low,
               p_value = p)
      }

      if (nrow(responder_data) > 0) {
        pancan_resp_fisher <- safe_fisher(responder_data) %>%
          mutate(scope = "PanCancer", p_label = p_label(p_value))
        readr::write_csv(pancan_resp_fisher, file.path(TABLE_DIR, "plot_v8_pancancer_responder_fisher_summary.csv"))

        per_cancer_resp_fisher <- responder_data %>%
          group_by(cancer_type) %>%
          group_modify(~ safe_fisher(.x)) %>%
          ungroup() %>%
          mutate(
            p_adj_BH = p.adjust(p_value, method = "BH"),
            p_label = p_label(p_value),
            p_adj_label = p_label(p_adj_BH)
          )
        readr::write_csv(per_cancer_resp_fisher, file.path(TABLE_DIR, "plot_v8_per_cancer_responder_fisher_summary.csv"))

        resp_prop <- responder_data %>%
          group_by(STEP_UP_CAF3_group, Responder_label) %>%
          summarise(n = n(), .groups = "drop") %>%
          group_by(STEP_UP_CAF3_group) %>%
          mutate(prop = n / sum(n)) %>%
          ungroup()

        p_resp <- ggplot(resp_prop, aes(x = STEP_UP_CAF3_group, y = prop, fill = Responder_label)) +
          geom_col(position = "fill", width = 0.7, color = "grey35", linewidth = 0.2) +
          geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
                    position = position_stack(vjust = 0.5), size = 3.2, color = "black") +
          annotate("text", x = 1.5, y = 1.08,
                   label = pancan_resp_fisher$p_label[1], size = 3.3) +
          scale_fill_manual(values = TIDE_RESPONSE_COLORS, drop = FALSE) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(0, 1.14), expand = expansion(mult = c(0, 0))) +
          labs(
            x = "CAF3-corrected STEP_UP group",
            y = "Predicted response proportion",
            fill = "TIDE prediction",
            title = paste0("Pan-cancer predicted responder proportion (", pancan_resp_fisher$p_label[1], ")")
          ) +
          theme_classic(base_size = 12) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                axis.text = element_text(color = "black"))
        save_plot_both(p_resp, "plot_v8_PanCancer_TIDE_predicted_responder_proportion_CAF3_STEP_UP", RESPONDER_DIR, width = 6.2, height = 4.8)

        resp_prop_cancer <- responder_data %>%
          group_by(cancer_type, STEP_UP_CAF3_group, Responder_label) %>%
          summarise(n = n(), .groups = "drop") %>%
          group_by(cancer_type, STEP_UP_CAF3_group) %>%
          mutate(prop = n / sum(n)) %>%
          ungroup()

        if (nrow(resp_prop_cancer) > 0) {
          resp_ann <- per_cancer_resp_fisher %>%
            transmute(cancer_type, x = 1.5, y = 1.08, label = p_label)
          p_resp_cancer <- ggplot(resp_prop_cancer, aes(x = STEP_UP_CAF3_group, y = prop, fill = Responder_label)) +
            geom_col(position = "fill", width = 0.7, color = "grey35", linewidth = 0.15) +
            geom_text(aes(label = scales::percent(prop, accuracy = 1)),
                      position = position_stack(vjust = 0.5), size = 1.8, color = "black") +
            geom_text(data = resp_ann, aes(x = x, y = y, label = label),
                      inherit.aes = FALSE, size = 2.0, color = "black") +
            facet_wrap(~ cancer_type) +
            scale_fill_manual(values = TIDE_RESPONSE_COLORS, drop = FALSE) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                               limits = c(0, 1.14), expand = expansion(mult = c(0, 0))) +
            labs(
              x = "CAF3-corrected STEP_UP group",
              y = "Predicted response proportion",
              fill = "TIDE prediction",
              title = wrap_title_step("Predicted responder proportion across cancers", 70)
            ) +
            theme_classic(base_size = 9) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                  strip.text = element_text(face = "bold", size = 7),
                  axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
                  axis.text.y = element_text(color = "black"))
          save_plot_both(p_resp_cancer, "plot_v8_PerCancer_TIDE_predicted_responder_proportion_CAF3_STEP_UP", RESPONDER_DIR, width = 11, height = 8)
        }
      } else {
        readr::write_csv(tibble(reason = "No valid Responder_bool rows"), file.path(TABLE_DIR, "plot_v8_responder_EMPTY_DEBUG.csv"))
        log_msg("No valid responder rows. Responder plots skipped.")
      }

      final_summary <- tibble(
        root = ROOT,
        step_file = STEP_FILE,
        tide_batch_file = ifelse(file.exists(TIDE_BATCH_FILE), TIDE_BATCH_FILE, NA_character_),
        output_dir = OUT_DIR,
        n_step = nrow(step),
        n_tide = nrow(tide),
        n_merged = nrow(merged),
        n_plot_long = nrow(plot_long),
        n_responder_rows = nrow(responder_data),
        n_cancers_merged = n_distinct(merged$cancer_type),
        group_definition = "Within each cancer: High if STEP_UP_residual_CAF3 > median; Low otherwise"
      )
      readr::write_csv(final_summary, file.path(TABLE_DIR, "plot_v8_final_summary.csv"))

      log_msg("============================================================")
      log_msg("v8 finished successfully")
      log_msg("Output tables: ", TABLE_DIR)
      log_msg("Pan-cancer plots: ", PANCAN_DIR)
      log_msg("Per-cancer plots: ", PER_CANCER_DIR)
      log_msg("Responder plots: ", RESPONDER_DIR)
      log_msg("============================================================")

    })

  })
}

if (isTRUE(RUN_TIDE_FIG7_EXTRA)) {
  run_optional_module_v11("18 TIDE Fig7A/C extra", {

    local({

      suppressPackageStartupMessages({
        library(readr)
        library(dplyr)
        library(tidyr)
        library(stringr)
        library(ggplot2)
        library(forcats)
        library(purrr)
        library(tibble)
      })

      FIG7_GROUP_COLORS <- STEP_GROUP_COLORS_REVERSED
      FIG7_RESPONDER_COLORS <- c("FALSE" = "#D9D9D9", "TRUE" = "#F4A6A6")

      ROOT <- MASTER_PROJECT_ROOT
      STEP_DIR <- MASTER_RESULT_DIR
      V8_DIR <- file.path(STEP_DIR, "18_TIDE_CAF3_STEP_UP", "v8_fixed_columns")
      TABLE_DIR <- file.path(V8_DIR, "tables")

      OUT_DIR <- file.path(V8_DIR, "Fig7AC_extra_v2")
      OUT_TABLE_DIR <- file.path(OUT_DIR, "tables")
      OUT_PANCAN_DIR <- file.path(OUT_DIR, "pan_cancer")
      OUT_RANK_DIR <- file.path(OUT_DIR, "per_cancer_ranked_TIDE")
      OUT_HEAT_DIR <- file.path(OUT_DIR, "per_cancer_responder_heatmap")
      OUT_BAR_DIR <- file.path(OUT_DIR, "per_cancer_responder_barplot")

      for (dd in c(OUT_DIR, OUT_TABLE_DIR, OUT_PANCAN_DIR, OUT_RANK_DIR, OUT_HEAT_DIR, OUT_BAR_DIR)) {
        dir.create(dd, recursive = TRUE, showWarnings = FALSE)
      }

      message2 <- function(...) cat(paste0(..., "\n"))

      save_plot_both <- function(p, filename, outdir, width = 7, height = 5) {
        pdf_file <- file.path(outdir, paste0(filename, ".pdf"))
        png_file <- file.path(outdir, paste0(filename, ".png"))
        print(p)  
        ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE)
        ggsave(png_file, p, width = width, height = height, dpi = 300)
      }

      parse_responder <- function(x) {
        x_chr <- tolower(trimws(as.character(x)))
        dplyr::case_when(
          x_chr %in% c("true", "t", "1", "yes", "y", "responder", "response") ~ TRUE,
          x_chr %in% c("false", "f", "0", "no", "n", "nonresponder", "non-responder", "no response") ~ FALSE,
          TRUE ~ NA
        )
      }

      safe_fisher_2x2 <- function(dat) {
        dat <- dat %>%
          filter(!is.na(STEP_group), !is.na(Responder_bool)) %>%
          mutate(
            STEP_group = factor(STEP_group, levels = c("High", "Low")),
            Responder_chr = ifelse(Responder_bool, "TRUE", "FALSE"),
            Responder_chr = factor(Responder_chr, levels = c("FALSE", "TRUE"))
          )

        tab <- table(dat$STEP_group, dat$Responder_chr)

        full_tab <- matrix(0, nrow = 2, ncol = 2,
                           dimnames = list(c("High", "Low"), c("FALSE", "TRUE")))
        full_tab[rownames(tab), colnames(tab)] <- tab

        n_total <- sum(full_tab)
        n_high <- sum(full_tab["High", ])
        n_low <- sum(full_tab["Low", ])
        n_resp_true <- sum(full_tab[, "TRUE"])
        n_resp_false <- sum(full_tab[, "FALSE"])

        resp_rate_high <- ifelse(n_high > 0, full_tab["High", "TRUE"] / n_high, NA_real_)
        resp_rate_low <- ifelse(n_low > 0, full_tab["Low", "TRUE"] / n_low, NA_real_)

        if (n_total == 0 || n_high == 0 || n_low == 0 || (n_resp_true == 0 && n_resp_false == 0)) {
          return(tibble(
            n_total = n_total,
            n_high = n_high,
            n_low = n_low,
            high_false = full_tab["High", "FALSE"],
            high_true = full_tab["High", "TRUE"],
            low_false = full_tab["Low", "FALSE"],
            low_true = full_tab["Low", "TRUE"],
            responder_rate_high = resp_rate_high,
            responder_rate_low = resp_rate_low,
            fisher_p = NA_real_,
            odds_ratio = NA_real_
          ))
        }

        ft <- tryCatch(fisher.test(full_tab), error = function(e) NULL)
        tibble(
          n_total = n_total,
          n_high = n_high,
          n_low = n_low,
          high_false = full_tab["High", "FALSE"],
          high_true = full_tab["High", "TRUE"],
          low_false = full_tab["Low", "FALSE"],
          low_true = full_tab["Low", "TRUE"],
          responder_rate_high = resp_rate_high,
          responder_rate_low = resp_rate_low,
          fisher_p = ifelse(is.null(ft), NA_real_, ft$p.value),
          odds_ratio = ifelse(is.null(ft), NA_real_, unname(ft$estimate))
        )
      }

      p_label <- function(p) {
        ifelse(is.na(p), "P = NA",
               ifelse(p < 1e-4, paste0("P.fisher = ", formatC(p, format = "e", digits = 2)),
                      paste0("P.fisher = ", signif(p, 3))))
      }

      candidate_files <- c(
        file.path(TABLE_DIR, "plot_v8_ALL_TCGA_TIDE_STEP_CAF3_merged_fixed_columns.csv"),
        file.path(TABLE_DIR, "plot_v8_ALL_TCGA_TIDE_STEP_CAF3_merged.csv"),
        file.path(TABLE_DIR, "ALL_TCGA_TIDE_STEP_CAF3_merged.csv"),
        file.path(TABLE_DIR, "plot_v8_merged_for_plotting.csv"),
        file.path(TABLE_DIR, "plot_v8_final_merged.csv"),
        file.path(STEP_DIR, "18_TIDE_CAF3_STEP_UP", "tables", "ALL_TCGA_TIDE_STEP_CAF3_merged.csv")
      )

      existing_candidates <- candidate_files[file.exists(candidate_files)]

      if (length(existing_candidates) == 0) {
        all_csv <- list.files(TABLE_DIR, pattern = "\\.csv$", full.names = TRUE)
        
        candidate_auto <- all_csv[str_detect(basename(all_csv), regex("merged|TIDE", ignore_case = TRUE))]
        if (length(candidate_auto) == 0) {
          stop("No merged table found in: ", TABLE_DIR, "\nPlease run v8 first or check the merged table name.")
        }
        merged_file <- candidate_auto[1]
      } else {
        merged_file <- existing_candidates[1]
      }

      message2("Reading merged table: ", merged_file)
      merged <- suppressMessages(readr::read_csv(merged_file, show_col_types = FALSE))

      readr::write_csv(tibble(merged_file_used = merged_file), file.path(OUT_TABLE_DIR, "Fig7AC_v2_merged_file_used.csv"))
      readr::write_csv(tibble(column_names = names(merged)), file.path(OUT_TABLE_DIR, "Fig7AC_v2_input_columns.csv"))

      if (!"TIDE" %in% names(merged)) stop("Column TIDE not found in merged table.")
      if (!"Responder" %in% names(merged)) stop("Column Responder not found in merged table.")

      cancer_col <- dplyr::case_when(
        "cancer_type_final" %in% names(merged) ~ "cancer_type_final",
        "cancer_from_file" %in% names(merged) ~ "cancer_from_file",
        "cancer_type" %in% names(merged) ~ "cancer_type",
        "Cancer.type" %in% names(merged) ~ "Cancer.type",
        TRUE ~ NA_character_
      )
      if (is.na(cancer_col)) stop("No cancer type column found.")

      sample_col <- dplyr::case_when(
        "sample_id_plot" %in% names(merged) ~ "sample_id_plot",
        "sample15" %in% names(merged) ~ "sample15",
        "sample_id_15" %in% names(merged) ~ "sample_id_15",
        "sample_id" %in% names(merged) ~ "sample_id",
        "sample_raw" %in% names(merged) ~ "sample_raw",
        TRUE ~ NA_character_
      )
      if (is.na(sample_col)) stop("No sample ID column found.")

      group_col <- dplyr::case_when(
        "STEP_UP_CAF3_group_for_plot" %in% names(merged) ~ "STEP_UP_CAF3_group_for_plot",
        "STEP_UP_CAF3_group" %in% names(merged) ~ "STEP_UP_CAF3_group",
        "STEP_UP_residual_CAF3_group" %in% names(merged) ~ "STEP_UP_residual_CAF3_group",
        TRUE ~ NA_character_
      )

      resid_col <- dplyr::case_when(
        "STEP_UP_CAF3_residual_for_plot" %in% names(merged) ~ "STEP_UP_CAF3_residual_for_plot",
        "STEP_UP_residual_CAF3" %in% names(merged) ~ "STEP_UP_residual_CAF3",
        "STEP_UP_CAF3_residual" %in% names(merged) ~ "STEP_UP_CAF3_residual",
        TRUE ~ NA_character_
      )

      readr::write_csv(
        tibble(
          cancer_col = cancer_col,
          sample_col = sample_col,
          group_col = ifelse(is.na(group_col), "NA", group_col),
          resid_col = ifelse(is.na(resid_col), "NA", resid_col)
        ),
        file.path(OUT_TABLE_DIR, "Fig7AC_v2_detected_columns.csv")
      )

      plot_df <- merged %>%
        transmute(
          cancer_type = as.character(.data[[cancer_col]]),
          sample_id = as.character(.data[[sample_col]]),
          TIDE = suppressWarnings(as.numeric(.data[["TIDE"]])),
          Responder_bool = parse_responder(.data[["Responder"]]),
          group_raw = if (!is.na(group_col)) as.character(.data[[group_col]]) else NA_character_,
          STEP_UP_residual_CAF3 = if (!is.na(resid_col)) suppressWarnings(as.numeric(.data[[resid_col]])) else NA_real_
        ) %>%
        mutate(
          group_raw_lower = tolower(trimws(group_raw)),
          STEP_group = case_when(
            group_raw_lower %in% c("high", "h", "up", "1") ~ "High",
            group_raw_lower %in% c("low", "l", "down", "0") ~ "Low",
            TRUE ~ NA_character_
          )
        )

      if (all(is.na(plot_df$STEP_group)) && !all(is.na(plot_df$STEP_UP_residual_CAF3))) {
        message2("Group column unavailable/all NA; deriving High/Low by per-cancer median of STEP_UP_residual_CAF3.")
        plot_df <- plot_df %>%
          group_by(cancer_type) %>%
          mutate(
            med_resid = median(STEP_UP_residual_CAF3, na.rm = TRUE),
            STEP_group = case_when(
              is.na(STEP_UP_residual_CAF3) | is.na(med_resid) ~ NA_character_,
              STEP_UP_residual_CAF3 > med_resid ~ "High",
              TRUE ~ "Low"
            )
          ) %>%
          ungroup() %>%
          select(-med_resid)
      }

      plot_df <- plot_df %>%
        filter(!is.na(cancer_type), !is.na(sample_id)) %>%
        mutate(
          STEP_group = factor(STEP_group, levels = c("High", "Low")),
          Responder_chr = case_when(
            is.na(Responder_bool) ~ NA_character_,
            Responder_bool ~ "TRUE",
            TRUE ~ "FALSE"
          ),
          Responder_chr = factor(Responder_chr, levels = c("FALSE", "TRUE"))
        )

      readr::write_csv(plot_df, file.path(OUT_TABLE_DIR, "Fig7AC_v2_clean_plot_input_all_samples.csv"))

      qc_by_cancer <- plot_df %>%
        group_by(cancer_type) %>%
        summarise(
          n_total = n(),
          n_with_TIDE = sum(!is.na(TIDE)),
          n_with_group = sum(!is.na(STEP_group)),
          n_with_responder = sum(!is.na(Responder_bool)),
          n_high = sum(STEP_group == "High", na.rm = TRUE),
          n_low = sum(STEP_group == "Low", na.rm = TRUE),
          n_responder_true = sum(Responder_bool == TRUE, na.rm = TRUE),
          n_responder_false = sum(Responder_bool == FALSE, na.rm = TRUE),
          .groups = "drop"
        )
      readr::write_csv(qc_by_cancer, file.path(OUT_TABLE_DIR, "Fig7AC_v2_QC_by_cancer.csv"))

      plot_df_valid <- plot_df %>%
        filter(!is.na(TIDE), !is.na(STEP_group), !is.na(Responder_bool))

      if (nrow(plot_df_valid) == 0) {
        stop("No valid rows with TIDE + STEP High/Low group + Responder. Check Fig7AC_v2_QC_by_cancer.csv")
      }

      pancan_ranked <- plot_df_valid %>%
        arrange(desc(TIDE)) %>%
        mutate(rank = row_number())
      readr::write_csv(pancan_ranked, file.path(OUT_TABLE_DIR, "Fig7A_pancancer_ranked_TIDE_input.csv"))
      pan_rank_p <- tryCatch(wilcox.test(TIDE ~ STEP_group, data = pancan_ranked,
                                         exact = FALSE)$p.value,
                             error = function(e) NA_real_)
      pan_rank_y <- max(pancan_ranked$TIDE, na.rm = TRUE)
      pan_rank_span <- diff(range(pancan_ranked$TIDE, na.rm = TRUE))
      if (!is.finite(pan_rank_span) || pan_rank_span == 0) pan_rank_span <- 1

      pA_pan <- ggplot(pancan_ranked, aes(x = rank, y = TIDE, fill = STEP_group)) +
        geom_col(width = 1) +
        geom_hline(yintercept = 0, linewidth = 0.3) +
        annotate("text", x = max(pancan_ranked$rank, na.rm = TRUE) * 0.5,
                 y = pan_rank_y + 0.08 * pan_rank_span,
                 label = format_p_step(pan_rank_p, "Wilcoxon p = "),
                 size = 3.4) +
        scale_fill_manual(values = FIG7_GROUP_COLORS, drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
        labs(
          x = "Patients ranked by TIDE score",
          y = "TIDE",
          fill = "CAF3-corrected\nSTEP_UP",
          title = "Pan-cancer ranked TIDE score"
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(color = "black")
        )

      save_plot_both(pA_pan, "Fig7A_like_pancancer_ranked_TIDE_by_CAF3_STEP_UP_group", OUT_PANCAN_DIR, width = 8.5, height = 4.2)

      pancan_fisher <- safe_fisher_2x2(plot_df_valid) %>%
        mutate(cancer_type = "PanCancer", p_adj_BH = fisher_p)
      readr::write_csv(pancan_fisher, file.path(OUT_TABLE_DIR, "Fig7C_pancancer_fisher_summary.csv"))

      pancan_count <- plot_df_valid %>%
        count(STEP_group, Responder_chr, name = "n") %>%
        complete(STEP_group = factor(c("High", "Low"), levels = c("High", "Low")),
                 Responder_chr = factor(c("FALSE", "TRUE"), levels = c("FALSE", "TRUE")),
                 fill = list(n = 0))
      readr::write_csv(pancan_count, file.path(OUT_TABLE_DIR, "Fig7C_pancancer_count_matrix_long.csv"))

      pC_pan <- ggplot(pancan_count, aes(x = Responder_chr, y = STEP_group, fill = n)) +
        geom_tile(color = "white", linewidth = 0.8) +
        geom_text(aes(label = n), size = 5, color = "black") +
        scale_fill_gradient(low = STEP_LIGHT_BLUE_LOW,
                            high = STEP_LIGHT_BLUE_HIGH) +
        labs(
          x = "Predicted responder by TIDE",
          y = "CAF3-corrected STEP_UP group",
          fill = "Count",
          title = p_label(pancan_fisher$fisher_p[1])
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.text = element_text(color = "black")
        )

      save_plot_both(pC_pan, "Fig7C_like_pancancer_STEP_UP_group_by_TIDE_responder_heatmap", OUT_PANCAN_DIR, width = 5.2, height = 4.6)

      pB_pan_dat <- plot_df_valid %>%
        count(STEP_group, Responder_chr, name = "n") %>%
        group_by(STEP_group) %>%
        mutate(percent = n / sum(n) * 100) %>%
        ungroup()
      readr::write_csv(pB_pan_dat, file.path(OUT_TABLE_DIR, "Fig7B_pancancer_responder_percent_input.csv"))

      pB_pan <- ggplot(pB_pan_dat, aes(x = STEP_group, y = percent, fill = Responder_chr)) +
        geom_col(width = 0.6, color = "grey35", linewidth = 0.2) +
        geom_text(aes(label = paste0(round(percent, 1), "%")),
                  position = position_stack(vjust = 0.5), size = 4, color = "black") +
        annotate("text", x = 1.5, y = 106,
                 label = p_label(pancan_fisher$fisher_p[1]), size = 3.3) +
        scale_fill_manual(values = FIG7_RESPONDER_COLORS, drop = FALSE) +
        scale_y_continuous(limits = c(0, 112), expand = expansion(mult = c(0, 0))) +
        labs(x = "CAF3-corrected STEP_UP group", y = "Percent", fill = "Responder", title = "Predicted response rate") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10), axis.text = element_text(color = "black"))

      save_plot_both(pB_pan, "Fig7B_like_pancancer_TIDE_responder_percentage", OUT_PANCAN_DIR, width = 4.5, height = 4.8)

      min_rank_n <- 10
      eligible_rank_cancers <- qc_by_cancer %>%
        filter(n_with_TIDE >= min_rank_n, n_high > 0, n_low > 0) %>%
        pull(cancer_type)

      readr::write_csv(tibble(cancer_type = eligible_rank_cancers), file.path(OUT_TABLE_DIR, "Fig7A_per_cancer_eligible_cancers.csv"))

      for (cc in eligible_rank_cancers) {
        dat_cc <- plot_df %>%
          filter(cancer_type == cc, !is.na(TIDE), !is.na(STEP_group)) %>%
          arrange(desc(TIDE)) %>%
          mutate(rank = row_number())

        readr::write_csv(dat_cc, file.path(OUT_TABLE_DIR, paste0("Fig7A_ranked_TIDE_input_", cc, ".csv")))
        p_rank_cc <- tryCatch(wilcox.test(TIDE ~ STEP_group, data = dat_cc,
                                          exact = FALSE)$p.value,
                              error = function(e) NA_real_)
        y_rank_cc <- max(dat_cc$TIDE, na.rm = TRUE)
        span_rank_cc <- diff(range(dat_cc$TIDE, na.rm = TRUE))
        if (!is.finite(span_rank_cc) || span_rank_cc == 0) span_rank_cc <- 1

        p <- ggplot(dat_cc, aes(x = rank, y = TIDE, fill = STEP_group)) +
          geom_col(width = 1) +
          geom_hline(yintercept = 0, linewidth = 0.3) +
          annotate("text", x = max(dat_cc$rank, na.rm = TRUE) * 0.5,
                   y = y_rank_cc + 0.08 * span_rank_cc,
                   label = format_p_step(p_rank_cc, "Wilcoxon p = "),
                   size = 3.2) +
          scale_fill_manual(values = FIG7_GROUP_COLORS, drop = FALSE) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
          labs(
            x = "Patients ranked by TIDE score",
            y = "TIDE",
            fill = "CAF3-corrected\nSTEP_UP",
            title = paste0(cc, ": ranked TIDE score")
          ) +
          theme_classic(base_size = 12) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text = element_text(color = "black")
          )

        save_plot_both(p, paste0("Fig7A_like_ranked_TIDE_", cc), OUT_RANK_DIR, width = 7.2, height = 4.0)
      }

      facet_rank <- plot_df %>%
        filter(cancer_type %in% eligible_rank_cancers, !is.na(TIDE), !is.na(STEP_group)) %>%
        group_by(cancer_type) %>%
        arrange(desc(TIDE), .by_group = TRUE) %>%
        mutate(rank = row_number()) %>%
        ungroup()

      if (nrow(facet_rank) > 0) {
        facet_rank_p <- facet_rank %>%
          group_by(cancer_type) %>%
          summarise(
            p = tryCatch(wilcox.test(TIDE[STEP_group == "High"],
                                     TIDE[STEP_group == "Low"],
                                     exact = FALSE)$p.value,
                         error = function(e) NA_real_),
            x = max(rank, na.rm = TRUE) * 0.5,
            ymax = max(TIDE, na.rm = TRUE),
            ymin = min(TIDE, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(y = ymax + 0.10 * ifelse(ymax > ymin, ymax - ymin, 1),
                 label = vapply(p, format_p_step, character(1), prefix = "p = "))
        pA_facet <- ggplot(facet_rank, aes(x = rank, y = TIDE, fill = STEP_group)) +
          geom_col(width = 1) +
          geom_hline(yintercept = 0, linewidth = 0.2) +
          geom_text(data = facet_rank_p, aes(x = x, y = y, label = label),
                    inherit.aes = FALSE, size = 2.2, color = "black") +
          scale_fill_manual(values = FIG7_GROUP_COLORS, drop = FALSE) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
          facet_wrap(~ cancer_type, scales = "free_x", ncol = 4) +
          labs(x = "Rank within cancer", y = "TIDE", fill = "CAF3-corrected\nSTEP_UP", title = "Per-cancer ranked TIDE score") +
          theme_classic(base_size = 10) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.text = element_text(face = "bold"),
            axis.text = element_text(color = "black")
          )
        save_plot_both(pA_facet, "Fig7A_like_per_cancer_ranked_TIDE_faceted", OUT_PANCAN_DIR, width = 12, height = 8.5)
      }

      per_cancer_fisher <- plot_df_valid %>%
        group_by(cancer_type) %>%
        group_modify(~ safe_fisher_2x2(.x)) %>%
        ungroup() %>%
        mutate(
          p_adj_BH = p.adjust(fisher_p, method = "BH"),
          responder_rate_delta_high_minus_low = responder_rate_high - responder_rate_low,
          signif = case_when(
            is.na(p_adj_BH) ~ "NA",
            p_adj_BH <= 0.0001 ~ "****",
            p_adj_BH <= 0.001 ~ "***",
            p_adj_BH <= 0.01 ~ "**",
            p_adj_BH <= 0.05 ~ "*",
            TRUE ~ "ns"
          )
        ) %>%
        arrange(fisher_p)

      readr::write_csv(per_cancer_fisher, file.path(OUT_TABLE_DIR, "Fig7C_per_cancer_fisher_summary.csv"))

      min_heatmap_n <- 10
      eligible_heat_cancers <- per_cancer_fisher %>%
        filter(n_total >= min_heatmap_n, n_high > 0, n_low > 0) %>%
        pull(cancer_type)
      readr::write_csv(tibble(cancer_type = eligible_heat_cancers), file.path(OUT_TABLE_DIR, "Fig7C_per_cancer_eligible_cancers.csv"))

      for (cc in eligible_heat_cancers) {
        dat_cc <- plot_df_valid %>% filter(cancer_type == cc)
        fish_cc <- per_cancer_fisher %>% filter(cancer_type == cc)

        count_cc <- dat_cc %>%
          count(STEP_group, Responder_chr, name = "n") %>%
          complete(STEP_group = factor(c("High", "Low"), levels = c("High", "Low")),
                   Responder_chr = factor(c("FALSE", "TRUE"), levels = c("FALSE", "TRUE")),
                   fill = list(n = 0))

        readr::write_csv(count_cc, file.path(OUT_TABLE_DIR, paste0("Fig7C_count_matrix_", cc, ".csv")))

        p_heat <- ggplot(count_cc, aes(x = Responder_chr, y = STEP_group, fill = n)) +
          geom_tile(color = "white", linewidth = 0.8) +
          geom_text(aes(label = n), size = 5, color = "black") +
          scale_fill_gradient(low = STEP_LIGHT_BLUE_LOW,
                              high = STEP_LIGHT_BLUE_HIGH) +
          labs(
            x = "Predicted responder by TIDE",
            y = "CAF3-corrected STEP_UP group",
            fill = "Count",
            title = paste0(cc, ": ", p_label(fish_cc$fisher_p[1]))
          ) +
          theme_classic(base_size = 12) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10), axis.text = element_text(color = "black"))

        save_plot_both(p_heat, paste0("Fig7C_like_responder_heatmap_", cc), OUT_HEAT_DIR, width = 5.2, height = 4.6)

        bar_cc <- dat_cc %>%
          count(STEP_group, Responder_chr, name = "n") %>%
          group_by(STEP_group) %>%
          mutate(percent = n / sum(n) * 100) %>%
          ungroup()

        p_bar <- ggplot(bar_cc, aes(x = STEP_group, y = percent, fill = Responder_chr)) +
          geom_col(width = 0.6, color = "grey35", linewidth = 0.2) +
          geom_text(aes(label = paste0(round(percent, 1), "%")),
                    position = position_stack(vjust = 0.5), size = 4, color = "black") +
          annotate("text", x = 1.5, y = 106,
                   label = p_label(fish_cc$fisher_p[1]), size = 3.2) +
          scale_fill_manual(values = FIG7_RESPONDER_COLORS, drop = FALSE) +
          scale_y_continuous(limits = c(0, 112), expand = expansion(mult = c(0, 0))) +
          labs(x = "CAF3-corrected STEP_UP group", y = "Percent", fill = "Responder", title = paste0(cc, ": predicted response rate")) +
          theme_classic(base_size = 12) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10), axis.text = element_text(color = "black"))

        save_plot_both(p_bar, paste0("Fig7B_like_responder_percentage_", cc), OUT_BAR_DIR, width = 4.5, height = 4.8)
      }

      facet_heat <- plot_df_valid %>%
        filter(cancer_type %in% eligible_heat_cancers) %>%
        count(cancer_type, STEP_group, Responder_chr, name = "n") %>%
        complete(cancer_type, STEP_group = factor(c("High", "Low"), levels = c("High", "Low")),
                 Responder_chr = factor(c("FALSE", "TRUE"), levels = c("FALSE", "TRUE")),
                 fill = list(n = 0))

      if (nrow(facet_heat) > 0) {
        facet_heat_ann <- per_cancer_fisher %>%
          filter(cancer_type %in% eligible_heat_cancers) %>%
          transmute(cancer_type, Responder_chr = factor("FALSE", levels = c("FALSE", "TRUE")),
                    STEP_group = factor("High", levels = c("High", "Low")),
                    label = p_label(fisher_p))
        pC_facet <- ggplot(facet_heat, aes(x = Responder_chr, y = STEP_group, fill = n)) +
          geom_tile(color = "white", linewidth = 0.4) +
          geom_text(aes(label = n), size = 3, color = "black") +
          scale_fill_gradient(low = STEP_LIGHT_BLUE_LOW,
                              high = STEP_LIGHT_BLUE_HIGH) +
          geom_text(data = facet_heat_ann,
                    aes(x = Responder_chr, y = STEP_group, label = label),
                    inherit.aes = FALSE, size = 1.8, hjust = -0.05,
                    vjust = -1.2, color = "black") +
          facet_wrap(~ cancer_type, ncol = 4) +
          labs(x = "Responder", y = "CAF3-corrected STEP_UP", fill = "Count", title = "Per-cancer STEP_UP group by TIDE responder") +
          theme_classic(base_size = 10) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10), strip.text = element_text(face = "bold"), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1))
        save_plot_both(pC_facet, "Fig7C_like_per_cancer_responder_heatmap_faceted", OUT_PANCAN_DIR, width = 12, height = 8.5)
      }

      message2("Done. Outputs saved to: ", OUT_DIR)
      message2("Key tables:")
      message2("  ", file.path(OUT_TABLE_DIR, "Fig7AC_v2_QC_by_cancer.csv"))
      message2("  ", file.path(OUT_TABLE_DIR, "Fig7C_per_cancer_fisher_summary.csv"))
      message2("Key plot dirs:")
      message2("  ", OUT_PANCAN_DIR)
      message2("  ", OUT_RANK_DIR)
      message2("  ", OUT_HEAT_DIR)
      message2("  ", OUT_BAR_DIR)

    })

  })
}

if (isTRUE(RUN_IMMUNE_ACTIVITY_EXTRA)) {
  run_optional_module_v11("17 immune activity", {

    local({

      message("========== APPEND: immune_activity_extra signatures ==========")

      suppressPackageStartupMessages({
        if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install/load pheatmap.")
        if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install/load ggplot2.")
      })

      base_dir <- MASTER_PROJECT_ROOT

      find_latest_result_dir <- function(base_dir) {
        
        is_valid_result_dir <- function(x) {
          dir.exists(x) &&
            grepl("TCGA_STEP_full_results_", basename(x)) &&
            dir.exists(file.path(x, "09_pan_cancer_summary"))
        }

        if (exists("out_dir") && is_valid_result_dir(out_dir)) return(normalizePath(out_dir))

        starts <- unique(normalizePath(c(getwd(), base_dir), mustWork = FALSE))
        for (st in starts) {
          cur <- st
          repeat {
            if (is_valid_result_dir(cur)) return(normalizePath(cur))
            parent <- dirname(cur)
            if (identical(parent, cur)) break
            cur <- parent
          }
        }

        cand <- character(0)
        for (st in starts) {
          if (dir.exists(st)) {
            tmp <- list.dirs(st, recursive = FALSE, full.names = TRUE)
            cand <- c(cand, tmp[grepl("TCGA_STEP_full_results_", basename(tmp))])
          }
        }
        cand <- unique(cand[file.exists(file.path(cand, "09_pan_cancer_summary",
                                                 "ALL_cancers_scores_purity_clinical_merged.csv"))])
        if (length(cand) == 0) {
          stop("Cannot find a valid TCGA_STEP_full_results_* folder containing 09_pan_cancer_summary. Please setwd() to TCGA_TME_correlation or the main result folder.")
        }
        cand[which.max(file.info(cand)$mtime)]
      }

      res_dir <- MASTER_RESULT_DIR
      message("Using result folder: ", res_dir)

      sum_dir <- file.path(res_dir, "09_pan_cancer_summary")
      merged_file <- file.path(sum_dir, "ALL_cancers_scores_purity_clinical_merged.csv")
      if (!file.exists(merged_file)) {
        stop("Cannot find merged table: ", merged_file)
      }

      extra_out <- file.path(res_dir, "17_immune_activity_extra_CAF3")
      dir.create(extra_out, showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(extra_out, "per_cancer_high_low_barplots"), showWarnings = FALSE, recursive = TRUE)

      extra_gmt_candidates <- c(
        INPUT_OVERRIDE_IMMUNE_ACTIVITY_GMT,
        file.path(base_dir, "gene set", "immune_activity_extra.gmt"),
        file.path(base_dir, "gene set", "immune_activity_extra_(1).gmt"),
        file.path(dirname(res_dir), "gene set", "immune_activity_extra.gmt")
      )
      extra_gmt <- extra_gmt_candidates[file.exists(extra_gmt_candidates)][1]
      if (length(extra_gmt) == 0 || is.na(extra_gmt)) {
        stop("Cannot find immune activity GMT. Tried: ",
             paste(extra_gmt_candidates, collapse = "; "))
      }

      exp_dir <- if (is_existing_dir_v112(INPUT_OVERRIDE_EXP_GCT_DIR)) INPUT_OVERRIDE_EXP_GCT_DIR else file.path(base_dir, "TCGA_201906cBioP", "exp_gct")
      if (!dir.exists(exp_dir)) {
        exp_dir2 <- file.path(dirname(res_dir), "TCGA_201906cBioP", "exp_gct")
        if (dir.exists(exp_dir2)) exp_dir <- exp_dir2
      }
      if (!dir.exists(exp_dir)) {
        stop("Cannot find expression GCT folder: ", exp_dir)
      }

      read_gmt_simple <- function(gmt_file) {
        lines <- readLines(gmt_file, warn = FALSE)
        sets <- list()
        for (ln in lines) {
          parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
          if (length(parts) < 3) next
          nm <- parts[1]
          genes <- parts[-c(1,2)]
          genes <- genes[!is.na(genes) & genes != ""]
          genes <- unique(toupper(trimws(genes)))
          if (length(genes) > 0) sets[[nm]] <- genes
        }
        sets
      }

      read_gct_expr <- function(gct_file) {
        x <- tryCatch(
          read.delim(gct_file, skip = 2, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
          error = function(e) read.delim(gct_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
        )
        gene_col <- intersect(c("Name", "Gene", "Hugo_Symbol", "gene_symbol", "GeneSymbol"), colnames(x))[1]
        if (is.na(gene_col)) gene_col <- colnames(x)[1]
        genes <- toupper(trimws(as.character(x[[gene_col]])))

        drop_cols <- intersect(c(gene_col, "Description", "description", "Entrez_Gene_Id"), colnames(x))
        sample_cols <- setdiff(colnames(x), drop_cols)
        mat <- as.matrix(x[, sample_cols, drop = FALSE])
        suppressWarnings(storage.mode(mat) <- "numeric")
        rownames(mat) <- genes

        keep <- !is.na(rownames(mat)) & rownames(mat) != ""
        mat <- mat[keep, , drop = FALSE]
        if (any(duplicated(rownames(mat)))) {
          mat <- rowsum(mat, group = rownames(mat), reorder = FALSE) /
            as.vector(table(rownames(mat)))[match(rownames(rowsum(mat, group = rownames(mat), reorder = FALSE)),
                                            names(table(rownames(mat))))]
        }

        colnames(mat) <- substr(gsub("\\.", "-", colnames(mat)), 1, 15)
        mat
      }

      zmean_scores <- function(expr_mat, gene_sets) {
        z <- t(scale(t(expr_mat)))
        z[is.nan(z)] <- NA
        out <- matrix(NA_real_, nrow = ncol(expr_mat), ncol = length(gene_sets))
        rownames(out) <- colnames(expr_mat)
        colnames(out) <- names(gene_sets)
        qc <- data.frame(signature = names(gene_sets), n_genes = NA_integer_, n_matched = NA_integer_,
                         matched_genes = NA_character_, missing_genes = NA_character_, stringsAsFactors = FALSE)

        for (i in seq_along(gene_sets)) {
          gs <- unique(toupper(gene_sets[[i]]))
          hit <- intersect(gs, rownames(z))
          miss <- setdiff(gs, rownames(z))
          qc$n_genes[i] <- length(gs)
          qc$n_matched[i] <- length(hit)
          qc$matched_genes[i] <- paste(hit, collapse = ";")
          qc$missing_genes[i] <- paste(miss, collapse = ";")
          if (length(hit) >= IMMUNE_EXTRA_MIN_MATCHED_GENES) {
            out[, i] <- colMeans(z[hit, , drop = FALSE], na.rm = TRUE)
          }
        }
        list(score = as.data.frame(out), qc = qc)
      }

      detect_col <- function(df, candidates) {
        hit <- candidates[candidates %in% colnames(df)]
        if (length(hit) == 0) return(NA_character_)
        hit[1]
      }

      spearman_one <- function(x, y) {
        ok <- is.finite(x) & is.finite(y)
        n <- sum(ok)
        if (n < IMMUNE_EXTRA_MIN_COR_PAIRS || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
          return(c(rho = NA_real_, p = NA_real_, n = n))
        }
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        c(rho = unname(ct$estimate), p = ct$p.value, n = n)
      }

      plot_heatmap_pdf_png <- function(mat, out_prefix, title, limits = NULL, cluster_rows = TRUE, cluster_cols = TRUE) {
        mat <- as.matrix(mat)
        mat <- mat[rowSums(!is.na(mat)) > 0, colSums(!is.na(mat)) > 0, drop = FALSE]
        if (nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))

        if (is.null(limits)) {
          mx <- max(abs(mat), na.rm = TRUE)
          if (!is.finite(mx) || mx == 0) mx <- 1
        } else {
          mx <- limits
        }
        breaks <- seq(-mx, mx, length.out = 101)
        cols <- colorRampPalette(c("#2b6cb0", "#e0f3f8", "#ffffbf", "#fdae61", "#d7191c"))(100)
        num <- matrix(sprintf("%.2f", mat), nrow = nrow(mat), dimnames = dimnames(mat))

        w <- max(8, min(26, 0.35 * ncol(mat) + 4))
        h <- max(4, min(16, 0.35 * nrow(mat) + 3))

        pdf(paste0(out_prefix, ".pdf"), width = w, height = h)
        pheatmap::pheatmap(mat, color = cols, breaks = breaks, main = title,
                           cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                           display_numbers = num, fontsize_number = 6,
                           fontsize_row = 9, fontsize_col = 8,
                           na_col = "grey90", border_color = "grey70")
        dev.off()

        png(paste0(out_prefix, ".png"), width = w, height = h, units = "in", res = 300)
        pheatmap::pheatmap(mat, color = cols, breaks = breaks, main = title,
                           cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                           display_numbers = num, fontsize_number = 6,
                           fontsize_row = 9, fontsize_col = 8,
                           na_col = "grey90", border_color = "grey70")
        dev.off()
      }

      df_all <- read.csv(merged_file, check.names = FALSE, stringsAsFactors = FALSE)
      score_col <- "STEP_UP_residual_CAF3"
      if (!score_col %in% colnames(df_all)) {
        score_candidates <- grep("STEP.*UP.*residual.*CAF3", colnames(df_all), value = TRUE, ignore.case = TRUE)
        if (length(score_candidates) == 0) stop("Cannot find STEP_UP_residual_CAF3 in merged table.")
        score_col <- score_candidates[1]
      }
      message("Using score_col: ", score_col)

      cancer_col <- detect_col(df_all, c("cancer_loop", "cancer_from_file", "cancer", "cancer_type", "CANCER"))
      if (is.na(cancer_col)) stop("Cannot find cancer column in merged table.")

      sample_col <- detect_col(df_all, c("sample15", "sample", "SAMPLE", "Sample", "PATIENT_ID", "patient_id"))
      if (is.na(sample_col)) stop("Cannot find sample column in merged table.")
      df_all$sample15_for_merge <- substr(gsub("\\.", "-", as.character(df_all[[sample_col]])), 1, 15)
      df_all[[cancer_col]] <- toupper(as.character(df_all[[cancer_col]]))

      extra_sets <- read_gmt_simple(extra_gmt)
      message("Loaded immune activity extra gene sets: ", paste(names(extra_sets), collapse = ", "))

      gct_files <- list.files(exp_dir, pattern = "_tcga_expression\\.gct$", full.names = TRUE, ignore.case = TRUE)
      if (length(gct_files) == 0) stop("No *_tcga_expression.gct found in: ", exp_dir)

      all_scores <- list()
      all_qc <- list()
      all_cor <- list()
      all_hl <- list()

      for (gf in gct_files) {
        cancer <- toupper(sub("_tcga_expression\\.gct$", "", basename(gf), ignore.case = TRUE))
        message("Processing immune_activity_extra: ", cancer)

        expr <- tryCatch(read_gct_expr(gf), error = function(e) {
          message("  Skip ", cancer, " due to GCT read error: ", conditionMessage(e))
          NULL
        })
        if (is.null(expr)) next

        sc <- zmean_scores(expr, extra_sets)
        scores <- sc$score
        
        scores$sample15_for_merge <- substr(gsub("\\.", "-", rownames(scores)), 1, 15)

        qc <- sc$qc
        qc$cancer <- cancer
        all_qc[[cancer]] <- qc

        subdf <- df_all[df_all[[cancer_col]] == cancer, , drop = FALSE]
        if (nrow(subdf) == 0) {
          subdf <- df_all
        }
        dat <- merge(subdf, scores, by = "sample15_for_merge", all = FALSE)
        
        if (!cancer_col %in% colnames(dat)) {
          dat[[cancer_col]] <- cancer
        }
        if (nrow(dat) < max(IMMUNE_EXTRA_MIN_COR_PAIRS, 2L * IMMUNE_EXTRA_MIN_GROUP_N)) {
          message("  Skip ", cancer, ": too few merged samples: ", nrow(dat))
          next
        }

        keep_cols <- unique(c("sample15_for_merge", cancer_col, score_col, names(extra_sets)))
        keep_cols <- keep_cols[keep_cols %in% colnames(dat)]
        all_scores[[cancer]] <- dat[, keep_cols, drop = FALSE]

        cors <- lapply(names(extra_sets), function(feat) {
          r <- spearman_one(as.numeric(dat[[score_col]]), as.numeric(dat[[feat]]))
          data.frame(cancer = cancer, feature = feat, rho = r["rho"], p = r["p"], n = r["n"], stringsAsFactors = FALSE)
        })
        cors <- do.call(rbind, cors)
        cors$FDR_by_cancer <- p.adjust(cors$p, method = "BH")
        all_cor[[cancer]] <- cors

        score_vec <- as.numeric(dat[[score_col]])
        cutoff <- median(score_vec, na.rm = TRUE)
        dat$STEP_UP_CAF3_group <- ifelse(score_vec > cutoff, "High", "Low")
        hls <- lapply(names(extra_sets), function(feat) {
          xh <- as.numeric(dat[dat$STEP_UP_CAF3_group == "High", feat])
          xl <- as.numeric(dat[dat$STEP_UP_CAF3_group == "Low", feat])
          okh <- is.finite(xh); okl <- is.finite(xl)
          p <- NA_real_
          if (sum(okh) >= IMMUNE_EXTRA_MIN_GROUP_N && sum(okl) >= IMMUNE_EXTRA_MIN_GROUP_N && length(unique(c(xh[okh], xl[okl]))) >= 2) {
            p <- suppressWarnings(wilcox.test(xh[okh], xl[okl])$p.value)
          }
          data.frame(
            cancer = cancer, feature = feat, cutoff = cutoff,
            n_high = sum(okh), n_low = sum(okl),
            median_high = ifelse(sum(okh) > 0, median(xh[okh], na.rm = TRUE), NA_real_),
            median_low  = ifelse(sum(okl) > 0, median(xl[okl], na.rm = TRUE), NA_real_),
            median_diff_high_minus_low = ifelse(sum(okh) > 0 && sum(okl) > 0, median(xh[okh], na.rm = TRUE) - median(xl[okl], na.rm = TRUE), NA_real_),
            p = p,
            stringsAsFactors = FALSE
          )
        })
        hls <- do.call(rbind, hls)
        hls$FDR_by_cancer <- p.adjust(hls$p, method = "BH")
        all_hl[[cancer]] <- hls

        plot_dat <- hls
        plot_dat$feature <- factor(plot_dat$feature, levels = plot_dat$feature[order(plot_dat$median_diff_high_minus_low)])
        plot_dat$p_label <- vapply(plot_dat$p, format_p_step, character(1), prefix = "p = ")
        p_bar <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = feature, y = median_diff_high_minus_low,
                                                        fill = median_diff_high_minus_low > 0)) +
          ggplot2::geom_col(color = "grey35", linewidth = 0.2) +
          ggplot2::geom_text(ggplot2::aes(label = p_label,
                                         hjust = ifelse(median_diff_high_minus_low >= 0, -0.05, 1.05)),
                             size = 2.6, color = "black") +
          ggplot2::scale_fill_manual(values = c("FALSE" = STEP_LOW_COLOR,
                                                "TRUE" = STEP_HIGH_COLOR),
                                     guide = "none") +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.22, 0.28))) +
          ggplot2::coord_flip(clip = "off") +
          ggplot2::geom_hline(yintercept = 0, linetype = 2) +
          theme_step(base_size = 10, title_size = 9.5,
                      axis_title_size = 8.5, axis_text_size = 7.5) +
          ggplot2::labs(title = wrap_title_step(paste0(cancer, ": STEP_UP_residual_CAF3 high-low immune activity"), 62),
                        x = NULL, y = "median(High) - median(Low)")
        ggplot2::ggsave(file.path(extra_out, "per_cancer_high_low_barplots", paste0(cancer, "_immune_activity_extra_high_low_bar.pdf")),
                        p_bar, width = 8, height = 5)
        ggplot2::ggsave(file.path(extra_out, "per_cancer_high_low_barplots", paste0(cancer, "_immune_activity_extra_high_low_bar.png")),
                        p_bar, width = 8, height = 5, dpi = 300)
      }

      score_out <- if (length(all_scores) > 0) do.call(rbind, all_scores) else data.frame()
      qc_out <- if (length(all_qc) > 0) do.call(rbind, all_qc) else data.frame()
      cor_out <- if (length(all_cor) > 0) do.call(rbind, all_cor) else data.frame()
      hl_out <- if (length(all_hl) > 0) do.call(rbind, all_hl) else data.frame()

      if (nrow(cor_out) > 0) cor_out$FDR_global <- p.adjust(cor_out$p, method = "BH")
      if (nrow(hl_out) > 0) hl_out$FDR_global <- p.adjust(hl_out$p, method = "BH")

      write.csv(score_out, file.path(extra_out, "ALL_cancers_immune_activity_extra_scores.csv"), row.names = FALSE)
      write.csv(qc_out, file.path(extra_out, "ALL_cancers_immune_activity_extra_gene_set_QC.csv"), row.names = FALSE)
      write.csv(cor_out, file.path(extra_out, "ALL_cancers_STEP_UP_residual_CAF3_immune_activity_extra_correlations.csv"), row.names = FALSE)
      write.csv(hl_out, file.path(extra_out, "ALL_cancers_STEP_UP_residual_CAF3_immune_activity_extra_high_low_median_diff.csv"), row.names = FALSE)

      if (nrow(cor_out) > 0) {
        cor_wide <- reshape(cor_out[, c("feature", "cancer", "rho")],
                            idvar = "feature", timevar = "cancer", direction = "wide")
        rownames(cor_wide) <- cor_wide$feature
        cor_wide$feature <- NULL
        colnames(cor_wide) <- sub("^rho\\.", "", colnames(cor_wide))
        cor_mat <- as.matrix(cor_wide)
        write.csv(cor_mat, file.path(extra_out, "STEP_UP_residual_CAF3_immune_activity_extra_correlation_matrix.csv"))
        plot_heatmap_pdf_png(
          cor_mat,
          file.path(extra_out, "PanCancer_STEP_UP_residual_CAF3_immune_activity_extra_correlation_heatmap"),
          "Pan-cancer STEP_UP_residual_CAF3 vs immune activity extra signatures",
          limits = 0.6
        )
      }

      if (nrow(hl_out) > 0) {
        hl_wide <- reshape(hl_out[, c("feature", "cancer", "median_diff_high_minus_low")],
                           idvar = "feature", timevar = "cancer", direction = "wide")
        rownames(hl_wide) <- hl_wide$feature
        hl_wide$feature <- NULL
        colnames(hl_wide) <- sub("^median_diff_high_minus_low\\.", "", colnames(hl_wide))
        hl_mat <- as.matrix(hl_wide)
        write.csv(hl_mat, file.path(extra_out, "STEP_UP_residual_CAF3_immune_activity_extra_high_low_median_diff_matrix.csv"))
        mx <- max(abs(hl_mat), na.rm = TRUE)
        if (!is.finite(mx) || mx == 0) mx <- 1
        plot_heatmap_pdf_png(
          hl_mat,
          file.path(extra_out, "PanCancer_STEP_UP_residual_CAF3_immune_activity_extra_high_low_median_diff_heatmap"),
          "STEP_UP_residual_CAF3 high-low median difference: immune activity extra",
          limits = mx
        )
      }

      message("========== DONE: immune_activity_extra appended module ==========")
      message("Outputs saved to: ", extra_out)

    })

  })
}

if (isTRUE(RUN_CHECKPOINT_MODULE)) {
  run_optional_module_v11("18 checkpoint", {

    local({

      suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(pheatmap)
        library(grid)
      })

      base_dir <- MASTER_PROJECT_ROOT

      gene_set_dir <- file.path(base_dir, "gene set")
      gmt_candidates <- c(
        INPUT_OVERRIDE_CHECKPOINT_GMT,
        file.path(gene_set_dir, "checkpoint_immunotherapy_gene_sets.gmt"),
        file.path(gene_set_dir, "checkpoint_immunotherapy_gene_sets(1).gmt")
      )
      gmt_file <- gmt_candidates[file.exists(gmt_candidates)][1]
      if (length(gmt_file) == 0 || is.na(gmt_file)) {
        stop("Cannot find checkpoint GMT. Tried: ", paste(gmt_candidates, collapse = "; "))
      }

      exp_dir <- if (is_existing_dir_v112(INPUT_OVERRIDE_EXP_GCT_DIR)) INPUT_OVERRIDE_EXP_GCT_DIR else file.path(base_dir, "TCGA_201906cBioP", "exp_gct")

      addon_name <- "18_checkpoint_immunotherapy_CAF3"

      main_score_col <- "STEP_UP_residual_CAF3"
      sample_merge_col <- "sample15_for_merge"
      cancer_col_candidates <- c("cancer_loop", "cancer_from_file", "cancer", "cancer_type", "CANCER")

      heatmap_limit_cor <- 0.6
      heatmap_limit_diff <- 1.2
      min_n_for_test <- CHECKPOINT_MIN_COR_PAIRS
      top_single_genes_to_boxplot <- c(
        "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
        "TIGIT", "IDO1", "CD276", "VTCN1", "ENTPD1", "NT5E",
        "TGFB1", "CXCL13", "TOX"
      )

      find_latest_result_dir <- function(base_dir) {
        dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
        hits <- dirs[grepl("TCGA_STEP_full_results_", basename(dirs))]
        if (length(hits) == 0) {
          stop("Cannot find TCGA_STEP_full_results_* under: ", base_dir)
        }
        hits[order(file.info(hits)$mtime, decreasing = TRUE)][1]
      }

      detect_cancer_col <- function(df) {
        hit <- cancer_col_candidates[cancer_col_candidates %in% colnames(df)]
        if (length(hit) == 0) return(NA_character_)
        hit[1]
      }

      read_gmt <- function(file) {
        if (!file.exists(file)) stop("GMT file not found: ", file)
        lines <- readLines(file, warn = FALSE)
        lines <- lines[nchar(trimws(lines)) > 0]
        gs <- list()
        for (ln in lines) {
          parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
          if (length(parts) < 3) next
          nm <- parts[1]
          genes <- unique(trimws(parts[-c(1, 2)]))
          genes <- genes[!is.na(genes) & genes != ""]
          gs[[nm]] <- genes
        }
        if (length(gs) == 0) stop("No gene sets parsed from: ", file)
        gs
      }

      read_gct_expr <- function(gct_file) {
        x <- data.table::fread(gct_file, skip = 2, data.table = FALSE)
        colnames(x)[1] <- "GeneSymbol"
        if ("Description" %in% colnames(x)) x$Description <- NULL

        x$GeneSymbol <- trimws(as.character(x$GeneSymbol))
        x <- x[!is.na(x$GeneSymbol) & x$GeneSymbol != "" & x$GeneSymbol != "NA", , drop = FALSE]

        sample_cols <- setdiff(colnames(x), "GeneSymbol")
        for (cc in sample_cols) {
          x[[cc]] <- suppressWarnings(as.numeric(as.character(x[[cc]])))
        }

        x2 <- aggregate(
          x[, sample_cols, drop = FALSE],
          by = list(GeneSymbol = x$GeneSymbol),
          FUN = function(v) mean(v, na.rm = TRUE)
        )
        x2[is.na(x2)] <- 0

        mat <- as.matrix(x2[, sample_cols, drop = FALSE])
        rownames(mat) <- x2$GeneSymbol

        keep <- rowSums(mat, na.rm = TRUE) > 0
        mat <- mat[keep, , drop = FALSE]
        mat
      }

      zmean_score <- function(expr_mat, gene_sets) {

        expr_mat <- expr_mat[apply(expr_mat, 1, function(v) sd(v, na.rm = TRUE) > 0), , drop = FALSE]
        z <- t(scale(t(expr_mat)))
        z[is.na(z)] <- 0

        score_list <- list()
        qc <- data.frame()

        for (nm in names(gene_sets)) {
          genes <- gene_sets[[nm]]
          matched <- intersect(genes, rownames(z))
          missing <- setdiff(genes, rownames(z))

          qc <- rbind(qc, data.frame(
            gene_set = nm,
            n_input = length(genes),
            n_matched = length(matched),
            n_missing = length(missing),
            matched_genes = paste(matched, collapse = ";"),
            missing_genes = paste(missing, collapse = ";"),
            stringsAsFactors = FALSE
          ))

          if (length(matched) == 0) {
            score_list[[nm]] <- rep(NA_real_, ncol(z))
          } else if (length(matched) == 1) {
            score_list[[nm]] <- as.numeric(z[matched, ])
          } else {
            score_list[[nm]] <- colMeans(z[matched, , drop = FALSE], na.rm = TRUE)
          }
        }

        scores <- as.data.frame(score_list, check.names = FALSE)
        scores$sample_id <- colnames(expr_mat)
        scores$sample15_for_merge <- substr(scores$sample_id, 1, 15)
        scores <- scores[, c("sample_id", "sample15_for_merge", names(gene_sets)), drop = FALSE]
        list(scores = scores, qc = qc)
      }

      safe_cor <- function(x, y) {
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < min_n_for_test ||
            length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
          return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
        }
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
      }

      save_plot_pdf_png <- function(plot_obj, filename_base, width = 10, height = 6) {
        pdf_file <- paste0(filename_base, ".pdf")
        png_file <- paste0(filename_base, ".png")
        ggplot2::ggsave(pdf_file, plot_obj, width = width, height = height, useDingbats = FALSE)
        ggplot2::ggsave(png_file, plot_obj, width = width, height = height, dpi = 300)
        print(plot_obj)  
      }

      pheatmap_save <- function(mat, filename_base, main = "", cluster_rows = TRUE, cluster_cols = TRUE,
                                color_limit = 1, width = 14, height = 5) {
        mat2 <- as.matrix(mat)
        mat2[!is.finite(mat2)] <- NA

        brks <- seq(-color_limit, color_limit, length.out = 101)
        cols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#D73027"))(100)

        pdf_file <- paste0(filename_base, ".pdf")
        png_file <- paste0(filename_base, ".png")

        ph <- pheatmap::pheatmap(
          mat2,
          color = cols,
          breaks = brks,
          cluster_rows = cluster_rows,
          cluster_cols = cluster_cols,
          display_numbers = TRUE,
          number_format = "%.2f",
          fontsize_number = 6,
          fontsize_row = 10,
          fontsize_col = 8,
          angle_col = 90,
          border_color = "grey80",
          main = main,
          silent = TRUE
        )

        grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
        grDevices::dev.off()

        grDevices::png(png_file, width = width, height = height, units = "in", res = 300)
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
        grDevices::dev.off()

        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
      }

      get_primary_sample_cols <- function(sample_cols) {

        is_tcga <- grepl("^TCGA-", sample_cols)
        if (!any(is_tcga)) return(sample_cols)
        keep <- substr(sample_cols, 14, 15) == "01"
        sample_cols[keep]
      }

      result_dir <- MASTER_RESULT_DIR
      message("Latest result dir: ", result_dir)

      summary_dir <- file.path(result_dir, "09_pan_cancer_summary")
      merged_file <- file.path(summary_dir, "ALL_cancers_scores_purity_clinical_merged.csv")
      if (!file.exists(merged_file)) {
        stop("Cannot find merged score table: ", merged_file)
      }

      out_dir <- file.path(result_dir, addon_name)
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      dir_scores <- file.path(out_dir, "01_checkpoint_scores")
      dir_qc <- file.path(out_dir, "02_gene_set_QC")
      dir_cor <- file.path(out_dir, "03_correlations")
      dir_highlow <- file.path(out_dir, "04_high_low")
      dir_single <- file.path(out_dir, "05_single_gene_checkpoint")
      dir_plots <- file.path(out_dir, "06_per_cancer_boxplots")

      for (dd in c(dir_scores, dir_qc, dir_cor, dir_highlow, dir_single, dir_plots)) {
        dir.create(dd, showWarnings = FALSE, recursive = TRUE)
      }

      df_all <- data.table::fread(merged_file, data.table = FALSE)
      cancer_col <- detect_cancer_col(df_all)
      if (is.na(cancer_col)) {
        stop("Cannot find cancer column in merged file. Columns include: ", paste(colnames(df_all), collapse = ", "))
      }
      if (!(main_score_col %in% colnames(df_all))) {
        stop("Cannot find main score column: ", main_score_col)
      }

      find_tcga_sample_col <- function(df) {
        
        candidates <- c(
          "sample15_for_merge", "sample_id", "sample", "sample_barcode", "Sample",
          "Tumor_Sample_Barcode", "tumor_sample_barcode", "barcode", "Barcode",
          "sampleID", "SampleID", "sample_id_short", "sample15", "sample16",
          "patient_sample", "PATIENT_SAMPLE"
        )
        hit <- candidates[candidates %in% colnames(df)]
        if (length(hit) > 0) return(hit[1])

        fuzzy <- grep("sample|barcode|tumor", colnames(df), ignore.case = TRUE, value = TRUE)
        fuzzy <- fuzzy[sapply(fuzzy, function(cc) any(grepl("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", as.character(df[[cc]]))))]
        if (length(fuzzy) > 0) return(fuzzy[1])

        for (cc in colnames(df)) {
          v <- as.character(df[[cc]])
          if (any(grepl("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", v))) return(cc)
        }

        NA_character_
      }

      if (!(sample_merge_col %in% colnames(df_all))) {
        sample_hit <- find_tcga_sample_col(df_all)
        if (is.na(sample_hit)) {
          diag_file <- file.path(out_dir, "DEBUG_merged_table_columns.csv")
          data.table::fwrite(data.frame(column = colnames(df_all)), diag_file)
          stop(
            "Cannot find sample merge column or TCGA sample-id column. ",
            "I wrote available columns to: ", diag_file,
            "\nPlease open that file and check which column contains TCGA sample IDs."
          )
        }
        message("Creating ", sample_merge_col, " from column: ", sample_hit)
        df_all[[sample_merge_col]] <- substr(as.character(df_all[[sample_hit]]), 1, 15)
      }

      df_all[[sample_merge_col]] <- substr(as.character(df_all[[sample_merge_col]]), 1, 15)
      df_all[[cancer_col]] <- toupper(as.character(df_all[[cancer_col]]))

      gene_sets <- read_gmt(gmt_file)
      message("Loaded gene sets: ", paste(names(gene_sets), collapse = ", "))

      all_checkpoint_genes <- unique(unlist(gene_sets, use.names = FALSE))
      single_genes <- unique(c(top_single_genes_to_boxplot, all_checkpoint_genes))

      gct_files <- list.files(exp_dir, pattern = "_tcga_expression\\.gct$", full.names = TRUE)
      if (length(gct_files) == 0) stop("No *_tcga_expression.gct files found in: ", exp_dir)

      all_scores <- list()
      all_qc <- list()
      all_single_expr <- list()

      for (gct in gct_files) {
        cancer <- toupper(sub("_tcga_expression\\.gct$", "", basename(gct)))
        message("Processing cancer: ", cancer)

        expr <- read_gct_expr(gct)

        primary_cols <- get_primary_sample_cols(colnames(expr))
        expr <- expr[, primary_cols, drop = FALSE]

        if (ncol(expr) < TME_MIN_TOTAL_SAMPLES) {
          warning("Skipping ", cancer, " because too few samples after primary-tumor filtering: ", ncol(expr))
          next
        }

        sc <- zmean_score(expr, gene_sets)
        scores <- sc$scores
        scores$cancer_loop <- cancer
        all_scores[[cancer]] <- scores

        qc <- sc$qc
        qc$cancer_loop <- cancer
        all_qc[[cancer]] <- qc

        present_genes <- intersect(single_genes, rownames(expr))
        expr2 <- expr[present_genes, , drop = FALSE]
        expr2 <- expr2[apply(expr2, 1, function(v) sd(v, na.rm = TRUE) > 0), , drop = FALSE]
        zexpr <- t(scale(t(expr2)))
        zexpr[is.na(zexpr)] <- 0

        if (nrow(zexpr) > 0) {

          sge <- as.data.frame(t(zexpr), check.names = FALSE)
          gene_cols_sge <- intersect(rownames(zexpr), colnames(sge))
          sge$sample_id <- rownames(sge)
          sge$sample15_for_merge <- substr(sge$sample_id, 1, 15)
          sge$cancer_loop <- cancer
          sge <- sge[, c("sample_id", "sample15_for_merge", "cancer_loop", gene_cols_sge), drop = FALSE]
          all_single_expr[[cancer]] <- sge
        }
      }

      checkpoint_scores <- bind_rows(all_scores)
      checkpoint_qc <- bind_rows(all_qc)
      single_expr_all <- bind_rows(all_single_expr)

      data.table::fwrite(checkpoint_scores, file.path(dir_scores, "ALL_cancers_checkpoint_immunotherapy_scores.csv"))
      data.table::fwrite(checkpoint_qc, file.path(dir_qc, "ALL_cancers_checkpoint_immunotherapy_gene_set_QC.csv"))
      data.table::fwrite(single_expr_all, file.path(dir_single, "ALL_cancers_checkpoint_single_gene_z_expression.csv"))

      df_key <- df_all %>%
        mutate(
          sample15_for_merge = .data[[sample_merge_col]],
          cancer_loop_for_merge = .data[[cancer_col]]
        ) %>%
        select(sample15_for_merge, cancer_loop_for_merge, all_of(main_score_col), everything())

      merged <- checkpoint_scores %>%
        left_join(
          df_key %>% select(sample15_for_merge, cancer_loop_for_merge, all_of(main_score_col)),
          by = c("sample15_for_merge" = "sample15_for_merge", "cancer_loop" = "cancer_loop_for_merge")
        )

      data.table::fwrite(merged, file.path(dir_scores, "ALL_cancers_checkpoint_scores_merged_with_STEP.csv"))

      single_merged <- single_expr_all %>%
        left_join(
          df_key %>% select(sample15_for_merge, cancer_loop_for_merge, all_of(main_score_col)),
          by = c("sample15_for_merge" = "sample15_for_merge", "cancer_loop" = "cancer_loop_for_merge")
        )

      data.table::fwrite(single_merged, file.path(dir_single, "ALL_cancers_checkpoint_single_gene_z_expression_merged_with_STEP.csv"))

      gene_set_names <- names(gene_sets)
      cancers <- sort(unique(merged$cancer_loop))

      cor_res <- list()
      highlow_res <- list()

      for (cc in cancers) {
        dat <- merged %>% filter(cancer_loop == cc, is.finite(.data[[main_score_col]]))
        if (nrow(dat) < min_n_for_test) next

        med <- median(dat[[main_score_col]], na.rm = TRUE)
        dat$STEP_group <- ifelse(dat[[main_score_col]] >= med, "High", "Low")

        for (gs in gene_set_names) {
          if (!(gs %in% colnames(dat))) next

          ct <- safe_cor(dat[[main_score_col]], dat[[gs]])
          cor_res[[length(cor_res) + 1]] <- data.frame(
            cancer_loop = cc,
            feature = gs,
            spearman_rho = ct["rho"],
            p_value = ct["p"],
            n = ct["n"],
            stringsAsFactors = FALSE
          )

          xh <- dat[[gs]][dat$STEP_group == "High"]
          xl <- dat[[gs]][dat$STEP_group == "Low"]
          if (sum(is.finite(xh)) >= CHECKPOINT_MIN_GROUP_N && sum(is.finite(xl)) >= CHECKPOINT_MIN_GROUP_N) {
            wt <- suppressWarnings(wilcox.test(xh, xl))
            highlow_res[[length(highlow_res) + 1]] <- data.frame(
              cancer_loop = cc,
              feature = gs,
              median_high = median(xh, na.rm = TRUE),
              median_low = median(xl, na.rm = TRUE),
              median_diff_high_minus_low = median(xh, na.rm = TRUE) - median(xl, na.rm = TRUE),
              p_value = wt$p.value,
              n_high = sum(is.finite(xh)),
              n_low = sum(is.finite(xl)),
              stringsAsFactors = FALSE
            )
          }
        }
      }

      cor_df <- bind_rows(cor_res) %>%
        group_by(feature) %>%
        mutate(FDR_by_feature = p.adjust(p_value, method = "BH")) %>%
        ungroup()

      highlow_df <- bind_rows(highlow_res) %>%
        group_by(feature) %>%
        mutate(FDR_by_feature = p.adjust(p_value, method = "BH")) %>%
        ungroup()

      data.table::fwrite(cor_df, file.path(dir_cor, "ALL_cancers_STEP_UP_residual_CAF3_checkpoint_gene_set_correlations.csv"))
      data.table::fwrite(highlow_df, file.path(dir_highlow, "ALL_cancers_STEP_UP_residual_CAF3_checkpoint_gene_set_high_low_median_diff.csv"))

      cor_mat <- cor_df %>%
        select(feature, cancer_loop, spearman_rho) %>%
        pivot_wider(names_from = cancer_loop, values_from = spearman_rho) %>%
        as.data.frame()
      rownames(cor_mat) <- cor_mat$feature
      cor_mat$feature <- NULL

      pheatmap_save(
        cor_mat,
        file.path(dir_cor, "PanCancer_STEP_UP_residual_CAF3_checkpoint_gene_set_correlation_heatmap"),
        main = "Pan-cancer STEP_UP_residual_CAF3 vs checkpoint/immunotherapy signatures",
        color_limit = heatmap_limit_cor,
        width = 15,
        height = 5.5
      )

      diff_mat <- highlow_df %>%
        select(feature, cancer_loop, median_diff_high_minus_low) %>%
        pivot_wider(names_from = cancer_loop, values_from = median_diff_high_minus_low) %>%
        as.data.frame()
      rownames(diff_mat) <- diff_mat$feature
      diff_mat$feature <- NULL

      pheatmap_save(
        diff_mat,
        file.path(dir_highlow, "PanCancer_STEP_UP_residual_CAF3_checkpoint_gene_set_high_low_median_diff_heatmap"),
        main = "STEP_UP_residual_CAF3 high-low median difference: checkpoint/immunotherapy signatures",
        color_limit = heatmap_limit_diff,
        width = 15,
        height = 5.5
      )

      single_gene_cols <- setdiff(colnames(single_merged), c("sample_id", "sample15_for_merge", "cancer_loop", main_score_col))
      single_cor_res <- list()
      single_highlow_res <- list()

      for (cc in sort(unique(single_merged$cancer_loop))) {
        dat <- single_merged %>% filter(cancer_loop == cc, is.finite(.data[[main_score_col]]))
        if (nrow(dat) < min_n_for_test) next

        med <- median(dat[[main_score_col]], na.rm = TRUE)
        dat$STEP_group <- ifelse(dat[[main_score_col]] >= med, "High", "Low")

        for (gn in single_gene_cols) {
          if (!(gn %in% colnames(dat))) next
          ct <- safe_cor(dat[[main_score_col]], dat[[gn]])
          single_cor_res[[length(single_cor_res) + 1]] <- data.frame(
            cancer_loop = cc,
            gene = gn,
            spearman_rho = ct["rho"],
            p_value = ct["p"],
            n = ct["n"],
            stringsAsFactors = FALSE
          )

          xh <- dat[[gn]][dat$STEP_group == "High"]
          xl <- dat[[gn]][dat$STEP_group == "Low"]
          if (sum(is.finite(xh)) >= CHECKPOINT_MIN_GROUP_N && sum(is.finite(xl)) >= CHECKPOINT_MIN_GROUP_N) {
            wt <- suppressWarnings(wilcox.test(xh, xl))
            single_highlow_res[[length(single_highlow_res) + 1]] <- data.frame(
              cancer_loop = cc,
              gene = gn,
              median_high = median(xh, na.rm = TRUE),
              median_low = median(xl, na.rm = TRUE),
              median_diff_high_minus_low = median(xh, na.rm = TRUE) - median(xl, na.rm = TRUE),
              p_value = wt$p.value,
              n_high = sum(is.finite(xh)),
              n_low = sum(is.finite(xl)),
              stringsAsFactors = FALSE
            )
          }
        }
      }

      single_cor_df <- bind_rows(single_cor_res) %>%
        group_by(gene) %>%
        mutate(FDR_by_gene = p.adjust(p_value, method = "BH")) %>%
        ungroup()

      single_highlow_df <- bind_rows(single_highlow_res) %>%
        group_by(gene) %>%
        mutate(FDR_by_gene = p.adjust(p_value, method = "BH")) %>%
        ungroup()

      data.table::fwrite(single_cor_df, file.path(dir_single, "ALL_cancers_STEP_UP_residual_CAF3_single_checkpoint_gene_correlations.csv"))
      data.table::fwrite(single_highlow_df, file.path(dir_single, "ALL_cancers_STEP_UP_residual_CAF3_single_checkpoint_gene_high_low_median_diff.csv"))

      focus_genes <- intersect(top_single_genes_to_boxplot, unique(single_cor_df$gene))
      if (length(focus_genes) >= 2) {
        sg_cor_mat <- single_cor_df %>%
          filter(gene %in% focus_genes) %>%
          select(gene, cancer_loop, spearman_rho) %>%
          pivot_wider(names_from = cancer_loop, values_from = spearman_rho) %>%
          as.data.frame()
        rownames(sg_cor_mat) <- sg_cor_mat$gene
        sg_cor_mat$gene <- NULL

        pheatmap_save(
          sg_cor_mat,
          file.path(dir_single, "PanCancer_STEP_UP_residual_CAF3_focused_single_checkpoint_gene_correlation_heatmap"),
          main = "Pan-cancer STEP_UP_residual_CAF3 vs focused checkpoint genes",
          color_limit = heatmap_limit_cor,
          width = 15,
          height = 6
        )

        sg_diff_mat <- single_highlow_df %>%
          filter(gene %in% focus_genes) %>%
          select(gene, cancer_loop, median_diff_high_minus_low) %>%
          pivot_wider(names_from = cancer_loop, values_from = median_diff_high_minus_low) %>%
          as.data.frame()
        rownames(sg_diff_mat) <- sg_diff_mat$gene
        sg_diff_mat$gene <- NULL

        pheatmap_save(
          sg_diff_mat,
          file.path(dir_single, "PanCancer_STEP_UP_residual_CAF3_focused_single_checkpoint_gene_high_low_median_diff_heatmap"),
          main = "STEP_UP_residual_CAF3 high-low median difference: focused checkpoint genes",
          color_limit = heatmap_limit_diff,
          width = 15,
          height = 6
        )
      }

      key_gene_sets <- c(
        "IMMUNE_CHECKPOINT_CORE_INHIBITORY",
        "IMMUNE_CHECKPOINT_INHIBITORY",
        "IMMUNE_CHECKPOINT_STIMULATORY",
        "ANTIGEN_PRESENTATION_MHC",
        "T_CELL_EXHAUSTION_EXTENDED",
        "ADENOSINE_AXIS",
        "TGF_BETA_IMMUNOSUPPRESSION"
      )
      key_gene_sets <- intersect(key_gene_sets, gene_set_names)

      for (cc in cancers) {
        dat <- merged %>% filter(cancer_loop == cc, is.finite(.data[[main_score_col]]))
        if (nrow(dat) < min_n_for_test) next

        med <- median(dat[[main_score_col]], na.rm = TRUE)
        dat$STEP_group <- factor(ifelse(dat[[main_score_col]] >= med, "High", "Low"), levels = c("Low", "High"))

        long <- dat %>%
          select(STEP_group, all_of(key_gene_sets)) %>%
          pivot_longer(cols = all_of(key_gene_sets), names_to = "feature", values_to = "score") %>%
          filter(is.finite(score))

        if (nrow(long) > 0) {
          ann <- long %>%
            group_by(feature) %>%
            summarise(
              p = tryCatch(wilcox.test(score[STEP_group == "High"],
                                       score[STEP_group == "Low"],
                                       exact = FALSE)$p.value,
                           error = function(e) NA_real_),
              ymax = max(score, na.rm = TRUE),
              ymin = min(score, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            mutate(x = 1.5,
                   y = ymax + 0.12 * ifelse(ymax > ymin, ymax - ymin, 1),
                   label = vapply(p, format_p_step, character(1), prefix = "p = "))
          p <- ggplot(long, aes(x = STEP_group, y = score,
                                fill = STEP_group, color = STEP_group)) +
            geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.82) +
            geom_jitter(width = 0.12, size = 0.6, alpha = 0.5, show.legend = FALSE) +
            geom_text(data = ann, aes(x = x, y = y, label = label),
                      inherit.aes = FALSE, color = "black", size = 2.6) +
            scale_fill_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
            scale_color_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
            facet_wrap(~feature, scales = "free_y", ncol = 3) +
            theme_step(base_size = 10, title_size = 9.5, axis_title_size = 8.5, axis_text_size = 7.5) +
            labs(
              title = paste0(cc, ": checkpoint/immunotherapy signature scores by STEP_UP_residual_CAF3 group"),
              x = "STEP_UP_residual_CAF3 group",
              y = "Z-mean signature score"
            ) +
            theme(
              strip.text = element_text(size = 7.5),
              plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              legend.position = "none"
            )
          save_plot_pdf_png(p, file.path(dir_plots, paste0(cc, "_checkpoint_gene_set_high_low_boxplots")), width = 11, height = 7)
        }
      }

      for (cc in sort(unique(single_merged$cancer_loop))) {
        dat <- single_merged %>% filter(cancer_loop == cc, is.finite(.data[[main_score_col]]))
        if (nrow(dat) < min_n_for_test) next

        med <- median(dat[[main_score_col]], na.rm = TRUE)
        dat$STEP_group <- factor(ifelse(dat[[main_score_col]] >= med, "High", "Low"), levels = c("Low", "High"))

        use_genes <- intersect(top_single_genes_to_boxplot, colnames(dat))
        if (length(use_genes) < 2) next

        long <- dat %>%
          select(STEP_group, all_of(use_genes)) %>%
          pivot_longer(cols = all_of(use_genes), names_to = "gene", values_to = "z_expression") %>%
          filter(is.finite(z_expression))

        if (nrow(long) > 0) {
          ann <- long %>%
            group_by(gene) %>%
            summarise(
              p = tryCatch(wilcox.test(z_expression[STEP_group == "High"],
                                       z_expression[STEP_group == "Low"],
                                       exact = FALSE)$p.value,
                           error = function(e) NA_real_),
              ymax = max(z_expression, na.rm = TRUE),
              ymin = min(z_expression, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            mutate(x = 1.5,
                   y = ymax + 0.12 * ifelse(ymax > ymin, ymax - ymin, 1),
                   label = vapply(p, format_p_step, character(1), prefix = "p = "))
          p <- ggplot(long, aes(x = STEP_group, y = z_expression,
                                fill = STEP_group, color = STEP_group)) +
            geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.82) +
            geom_jitter(width = 0.12, size = 0.6, alpha = 0.5, show.legend = FALSE) +
            geom_text(data = ann, aes(x = x, y = y, label = label),
                      inherit.aes = FALSE, color = "black", size = 2.3) +
            scale_fill_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
            scale_color_manual(values = STEP_GROUP_COLORS, drop = FALSE) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
            facet_wrap(~gene, scales = "free_y", ncol = 5) +
            theme_step(base_size = 10, title_size = 9.5, axis_title_size = 8.5, axis_text_size = 7.5) +
            labs(
              title = paste0(cc, ": focused checkpoint genes by STEP_UP_residual_CAF3 group"),
              x = "STEP_UP_residual_CAF3 group",
              y = "Gene-wise z-expression"
            ) +
            theme(
              strip.text = element_text(size = 7.5),
              plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              legend.position = "none"
            )
          save_plot_pdf_png(p, file.path(dir_plots, paste0(cc, "_focused_single_checkpoint_gene_high_low_boxplots")), width = 13, height = 7)
        }
      }

      readme <- c(
        "Checkpoint / immunotherapy add-on for TCGA STEP CAF3 analysis",
        "",
        paste0("Main score: ", main_score_col),
        paste0("GMT file: ", gmt_file),
        paste0("Expression folder: ", exp_dir),
        paste0("Merged score table: ", merged_file),
        "",
        "Important outputs:",
        "01_checkpoint_scores/ALL_cancers_checkpoint_immunotherapy_scores.csv",
        "01_checkpoint_scores/ALL_cancers_checkpoint_scores_merged_with_STEP.csv",
        "02_gene_set_QC/ALL_cancers_checkpoint_immunotherapy_gene_set_QC.csv",
        "03_correlations/ALL_cancers_STEP_UP_residual_CAF3_checkpoint_gene_set_correlations.csv",
        "03_correlations/PanCancer_STEP_UP_residual_CAF3_checkpoint_gene_set_correlation_heatmap.pdf/png",
        "04_high_low/ALL_cancers_STEP_UP_residual_CAF3_checkpoint_gene_set_high_low_median_diff.csv",
        "04_high_low/PanCancer_STEP_UP_residual_CAF3_checkpoint_gene_set_high_low_median_diff_heatmap.pdf/png",
        "05_single_gene_checkpoint/ALL_cancers_STEP_UP_residual_CAF3_single_checkpoint_gene_correlations.csv",
        "05_single_gene_checkpoint/ALL_cancers_STEP_UP_residual_CAF3_single_checkpoint_gene_high_low_median_diff.csv",
        "05_single_gene_checkpoint/PanCancer_STEP_UP_residual_CAF3_focused_single_checkpoint_gene_correlation_heatmap.pdf/png",
        "05_single_gene_checkpoint/PanCancer_STEP_UP_residual_CAF3_focused_single_checkpoint_gene_high_low_median_diff_heatmap.pdf/png",
        "06_per_cancer_boxplots/*_checkpoint_gene_set_high_low_boxplots.pdf/png",
        "06_per_cancer_boxplots/*_focused_single_checkpoint_gene_high_low_boxplots.pdf/png",
        "",
        "How to interpret:",
        "1) Correlation heatmap: positive Spearman rho means higher STEP_UP_residual_CAF3 is associated with higher checkpoint/immunotherapy signature activity.",
        "2) High-low heatmap: positive median_diff_high_minus_low means STEP-high tumors have higher signature score than STEP-low tumors within the same cancer type.",
        "3) Single-gene heatmaps are useful for direct checkpoint targets such as CD274, PDCD1, CTLA4, LAG3, HAVCR2, TIGIT, IDO1, CD276, VTCN1, ENTPD1, NT5E, TGFB1.",
        "4) Check gene-set QC before interpretation, especially PVRL2/NECTIN2 and other alias-sensitive genes."
      )
      writeLines(readme, file.path(out_dir, "README_checkpoint_immunotherapy_addon.txt"))

      message("Done.")
      message("Output folder: ", out_dir)

    })

  })
}

if (isTRUE(RUN_SUMMARY_BUBBLE)) {
  run_optional_module_v11("19 summary bubble", {

    local({

      suppressPackageStartupMessages({
        library(data.table)
        library(ggplot2)
      })

      base_dir <- getwd()
      main_score_col <- "STEP_UP_residual_CAF3"
      out_subdir <- "19_summary_bubble_plot_CAF3"

      find_latest_results_dir <- function(base_dir) {
        dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
        hits <- dirs[grepl("TCGA_STEP_full_results_", basename(dirs))]
        if (length(hits) == 0 && grepl("TCGA_STEP_full_results_", basename(base_dir))) return(base_dir)
        if (length(hits) == 0) stop("Cannot find TCGA_STEP_full_results_* under: ", base_dir)
        hits[order(file.info(hits)$mtime, decreasing = TRUE)][1]
      }

      res_dir <- MASTER_RESULT_DIR
      out_dir <- file.path(res_dir, out_subdir)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      message("Using results dir: ", res_dir)
      message("Output dir: ", out_dir)

      read_any_csv <- function(f) data.table::fread(f, data.table = FALSE, check.names = FALSE)

      pick_col <- function(df, candidates) {
        hit <- candidates[candidates %in% colnames(df)]
        if (length(hit) == 0) return(NA_character_)
        hit[1]
      }

      standardize_cor_table <- function(df, source_name) {
        cancer_col <- pick_col(df, c("cancer", "cancer_loop", "cancer_type", "CANCER", "Cancer", "TCGA_cancer"))
        feature_col <- pick_col(df, c("feature", "gene_set", "signature", "module", "cell_type", "variable", "term", "GeneSet", "gene", "checkpoint_gene"))
        rho_col <- pick_col(df, c("rho", "spearman_rho", "cor", "correlation", "Spearman_rho", "estimate", "r"))
        p_col <- pick_col(df, c("p", "p_value", "pvalue", "P.Value", "pval"))
        fdr_col <- pick_col(df, c("FDR", "fdr", "padj", "p_adj", "q_value", "qvalue"))

        if (is.na(feature_col)) {
          stat_cols <- c(cancer_col, rho_col, p_col, fdr_col, "n", "n_samples", "method")
          char_cols <- colnames(df)[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
          char_cols <- setdiff(char_cols, stat_cols)
          if (length(char_cols) > 0) feature_col <- char_cols[1]
        }

        if (is.na(cancer_col) || is.na(feature_col) || is.na(rho_col)) {
          warning("Skipping table from ", source_name, ": cannot infer cancer/feature/rho columns. Columns: ", paste(colnames(df), collapse = ", "))
          return(NULL)
        }

        ans <- data.frame(
          cancer = as.character(df[[cancer_col]]),
          feature = as.character(df[[feature_col]]),
          rho = suppressWarnings(as.numeric(df[[rho_col]])),
          source = source_name,
          stringsAsFactors = FALSE
        )
        ans$p <- if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_
        ans$FDR <- if (!is.na(fdr_col)) suppressWarnings(as.numeric(df[[fdr_col]])) else NA_real_
        ans
      }

      cor_files <- c(
        list.files(file.path(res_dir, "13_TME_correlations"), pattern = "cor.*\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
        list.files(file.path(res_dir, "17_immune_activity_extra_CAF3"), pattern = "cor.*\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
        list.files(file.path(res_dir, "18_checkpoint_immunotherapy_CAF3"), pattern = "cor.*\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
      )
      cor_files <- unique(cor_files[file.exists(cor_files)])
      if (length(cor_files) == 0) stop("No correlation csv files found in 13/17/18 result folders. Please run TME, immune activity, and checkpoint modules first.")

      all_cor <- list()
      for (f in cor_files) {
        message("Reading correlation table: ", f)
        df <- tryCatch(read_any_csv(f), error = function(e) NULL)
        if (is.null(df)) next
        src <- if (grepl("17_immune_activity", f)) "immune_activity_extra" else if (grepl("18_checkpoint", f)) "checkpoint" else "TME"
        st <- standardize_cor_table(df, src)
        if (!is.null(st)) all_cor[[length(all_cor) + 1]] <- st
      }
      all_cor <- do.call(rbind, all_cor)
      all_cor <- all_cor[!is.na(all_cor$rho) & !is.na(all_cor$cancer) & !is.na(all_cor$feature), ]

      all_cor$key <- paste(all_cor$cancer, all_cor$feature, all_cor$source, sep = "||")
      all_cor <- all_cor[!duplicated(all_cor$key), ]
      all_cor$key <- NULL

      module_map <- data.frame(
        module_label = c(
          "ImmuneScore",
          "StromalScore",
          "Macrophage/M2",
          "Monocyte/Myeloid",
          "Treg",
          "NK activated",
          "Cytotoxic lymphocytes",
          "T cell core",
          "CD8 T extended",
          "CYT score",
          "Cytotoxic T/NK",
          "T cell exhaustion",
          "Checkpoint inhibitory",
          "Checkpoint core",
          "Checkpoint stimulatory",
          "Antigen presentation/MHC",
          "TIM3/HAVCR2 axis",
          "TIGIT axis",
          "Adenosine axis",
          "TGF-beta axis"
        ),
        pattern = c(
          "ESTIMATE_ImmuneScore|ImmuneScore$",
          "ESTIMATE_StromalScore|StromalScore$",
          "Macrophage M2|Macrophages_co|Macrophage|M2",
          "Monocytes|Monocytic lineage|Myeloid",
          "Tregs|T cells regulatory|Regulatory T",
          "NK activated|NK cells activated",
          "Cytotoxic lymphocytes_co|Cytotoxic lymphocytes_M",
          "^T_cell_core$",
          "^CD8_T_extended$",
          "^CYT_score$",
          "^Cytotoxic_T_NK$",
          "T_CELL_EXHAUSTION_EXTENDED|^T_cell_exhaustion$",
          "^IMMUNE_CHECKPOINT_INHIBITORY$",
          "^IMMUNE_CHECKPOINT_CORE_INHIBITORY$",
          "^IMMUNE_CHECKPOINT_STIMULATORY$",
          "^ANTIGEN_PRESENTATION_MHC$",
          "^TIM3_GALECTIN9_AXIS$",
          "^TIGIT_AXIS$",
          "^ADENOSINE_AXIS$",
          "^TGF_BETA_IMMUNOSUPPRESSION$"
        ),
        stringsAsFactors = FALSE
      )

      select_one_feature <- function(cor_df, pattern) {
        hit <- cor_df[grepl(pattern, cor_df$feature, ignore.case = TRUE), ]
        if (nrow(hit) == 0) return(NULL)
        
        hit$absrho <- abs(hit$rho)
        hit <- hit[order(hit$cancer, -hit$absrho), ]
        hit <- hit[!duplicated(hit$cancer), ]
        hit$absrho <- NULL
        hit
      }

      bubble_df <- list()
      for (i in seq_len(nrow(module_map))) {
        tmp <- select_one_feature(all_cor, module_map$pattern[i])
        if (!is.null(tmp) && nrow(tmp) > 0) {
          tmp$module <- module_map$module_label[i]
          bubble_df[[length(bubble_df) + 1]] <- tmp
        }
      }
      bubble_df <- do.call(rbind, bubble_df)
      if (is.null(bubble_df) || nrow(bubble_df) == 0) stop("No selected modules were found in correlation tables. Check feature names in output CSVs.")

      bubble_df$module <- factor(bubble_df$module, levels = rev(module_map$module_label[module_map$module_label %in% bubble_df$module]))
      
      cancer_order <- aggregate(rho ~ cancer, bubble_df, median, na.rm = TRUE)
      cancer_order <- cancer_order[order(cancer_order$rho, decreasing = TRUE), "cancer"]
      bubble_df$cancer <- factor(bubble_df$cancer, levels = cancer_order)

      bubble_df$neglog10FDR <- ifelse(!is.na(bubble_df$FDR) & bubble_df$FDR > 0, -log10(bubble_df$FDR), NA_real_)
      bubble_df$size_value <- ifelse(!is.na(bubble_df$neglog10FDR), pmin(bubble_df$neglog10FDR, 6), pmin(abs(bubble_df$rho) * 5, 6))
      bubble_df$significance <- ifelse(!is.na(bubble_df$FDR) & bubble_df$FDR < 0.05, "FDR<0.05",
                                       ifelse(!is.na(bubble_df$p) & bubble_df$p < 0.05, "p<0.05", "ns/NA"))
      bubble_df$p_text <- ifelse(!is.na(bubble_df$p) & bubble_df$p < 0.05,
                                 paste0("p=", vapply(bubble_df$p, format_p_short_step, character(1))), "")

      fwrite(bubble_df, file.path(out_dir, "summary_bubble_plot_input_long_table.csv"))
      fwrite(all_cor, file.path(out_dir, "all_collected_correlation_tables_standardized.csv"))

      p <- ggplot(bubble_df, aes(x = cancer, y = module)) +
        geom_point(aes(size = size_value, fill = rho), shape = 21, color = "grey25", stroke = 0.2) +
        geom_text(aes(label = p_text), size = 1.35, color = "black", check_overlap = TRUE) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                             limits = c(-0.6, 0.6), oob = scales::squish, name = "Spearman rho") +
        scale_size_continuous(range = c(1.2, 6.5), name = "significance\n/ |rho|") +
        labs(
          title = paste0(main_score_col, " summary immune-module associations"),
          x = "TCGA cancer type",
          y = "Immune/TME module"
        ) +
        theme_step(base_size = 9, title_size = 9.5, axis_title_size = 8.5, axis_text_size = 7) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )

      print(p)
      ggsave(file.path(out_dir, "STEP_UP_residual_CAF3_summary_immune_module_bubble_plot.pdf"), p, width = 15, height = 7.5)
      ggsave(file.path(out_dir, "STEP_UP_residual_CAF3_summary_immune_module_bubble_plot.png"), p, width = 15, height = 7.5, dpi = 300)

      keep_compact <- c("ImmuneScore", "Macrophage/M2", "T cell core", "CYT score", "T cell exhaustion",
                        "Checkpoint core", "Antigen presentation/MHC", "TIM3/HAVCR2 axis", "TIGIT axis", "Adenosine axis")
      compact_df <- bubble_df[as.character(bubble_df$module) %in% keep_compact, ]
      compact_df$module <- factor(as.character(compact_df$module), levels = rev(keep_compact[keep_compact %in% as.character(compact_df$module)]))

      p2 <- ggplot(compact_df, aes(x = cancer, y = module)) +
        geom_point(aes(size = size_value, fill = rho), shape = 21, color = "grey25", stroke = 0.2) +
        geom_text(aes(label = p_text), size = 1.45, color = "black", check_overlap = TRUE) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                             limits = c(-0.6, 0.6), oob = scales::squish, name = "Spearman rho") +
        scale_size_continuous(range = c(1.5, 7), name = "significance\n/ |rho|") +
        labs(title = paste0(main_score_col, " compact immune-module summary"), x = "TCGA cancer type", y = NULL) +
        theme_step(base_size = 9, title_size = 9.5, axis_title_size = 8.5, axis_text_size = 7) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              panel.grid.major = element_line(color = "grey90"),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold"))

      print(p2)
      ggsave(file.path(out_dir, "STEP_UP_residual_CAF3_compact_summary_bubble_plot.pdf"), p2, width = 14, height = 5.5)
      ggsave(file.path(out_dir, "STEP_UP_residual_CAF3_compact_summary_bubble_plot.png"), p2, width = 14, height = 5.5, dpi = 300)

      message("Done. Check: ", out_dir)

    })

  })
}

if (isTRUE(RUN_CHECKPOINT_SYMBOL_RESCUE)) {
  run_optional_module_v11("20 checkpoint symbol rescue", {

    local({

      suppressPackageStartupMessages({
        library(data.table)
      })

      base_dir <- getwd()
      gct_dir <- file.path(base_dir, "TCGA_201906cBioP", "exp_gct")
      gmt_file <- file.path(base_dir, "gene set", "checkpoint_immunotherapy_gene_sets.gmt")
      out_dir <- file.path(MASTER_RESULT_DIR, "20_checkpoint_gene_symbol_rescue_report")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      if (!dir.exists(gct_dir)) stop("Cannot find GCT dir: ", gct_dir)
      if (!file.exists(gmt_file)) stop("Cannot find GMT file: ", gmt_file)

      read_gmt <- function(gmt) {
        lines <- readLines(gmt, warn = FALSE)
        ans <- lapply(lines, function(z) {
          sp <- strsplit(z, "\t")[[1]]
          if (length(sp) < 3) return(NULL)
          data.frame(gene_set = sp[1], gene = unique(sp[-c(1, 2)]), stringsAsFactors = FALSE)
        })
        do.call(rbind, ans)
      }

      gmt_df <- read_gmt(gmt_file)
      gmt_df$gene <- trimws(gmt_df$gene)
      gmt_df <- gmt_df[gmt_df$gene != "" & !is.na(gmt_df$gene), ]

      target_genes <- unique(c(
        gmt_df$gene,
        "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "LGALS9", "TIGIT", "IDO1",
        "CD276", "VTCN1", "ENTPD1", "NT5E", "TGFB1", "CXCL13", "TOX", "TOX2",
        "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9", "B2M", "NLRC5", "CIITA"
      ))

      alias_map <- list(
        PDCD1 = c("PDCD1", "PD1", "PD-1", "CD279"),
        CD274 = c("CD274", "PDCD1L1", "PDL1", "PD-L1", "B7-H1"),
        PDCD1LG2 = c("PDCD1LG2", "PDCD1L2", "PDL2", "PD-L2", "CD273"),
        CTLA4 = c("CTLA4", "CTLA-4", "CD152"),
        HAVCR2 = c("HAVCR2", "TIM3", "TIM-3", "TIMD3"),
        LGALS9 = c("LGALS9", "GAL9", "Galectin-9"),
        TIGIT = c("TIGIT", "VSTM3"),
        IDO1 = c("IDO1", "IDO"),
        VTCN1 = c("VTCN1", "B7H4", "B7-H4", "B7S1"),
        CD276 = c("CD276", "B7H3", "B7-H3"),
        ENTPD1 = c("ENTPD1", "CD39"),
        NT5E = c("NT5E", "CD73"),
        PVR = c("PVR", "CD155"),
        NECTIN2 = c("NECTIN2", "PVRL2", "CD112"),
        VSIR = c("VSIR", "C10orf54", "VISTA", "B7H5"),
        TGFB1 = c("TGFB1", "TGF-BETA1", "TGFbeta1", "TGFB"),
        CXCL13 = c("CXCL13", "BCA1", "BLC"),
        TOX = c("TOX"),
        TOX2 = c("TOX2")
      )

      for (g in target_genes) if (is.null(alias_map[[g]])) alias_map[[g]] <- g

      files <- list.files(gct_dir, pattern = "_tcga_expression\\.gct$", full.names = TRUE)
      if (length(files) == 0) stop("No *_tcga_expression.gct files found in: ", gct_dir)

      scan_one_gct <- function(f) {
        cancer <- toupper(sub("_tcga_expression\\.gct$", "", basename(f)))
        x <- fread(f, skip = 2, data.table = FALSE, check.names = FALSE, nThread = 2)
        colnames(x)[1] <- "Symbol"
        if (!"Description" %in% colnames(x)) x$Description <- NA_character_
        symbols <- trimws(as.character(x$Symbol))
        desc <- trimws(as.character(x$Description))
        symbols_upper <- toupper(symbols)
        desc_upper <- toupper(desc)

        rows <- lapply(names(alias_map), function(canonical) {
          aliases <- unique(c(canonical, alias_map[[canonical]]))
          alias_upper <- toupper(aliases)
          exact_hit <- symbols[symbols_upper %in% alias_upper]
          desc_hit <- symbols[desc_upper %in% alias_upper | grepl(paste(alias_upper, collapse = "|"), desc_upper)]
          data.frame(
            cancer = cancer,
            canonical_gene = canonical,
            found_exact_symbol = length(exact_hit) > 0,
            exact_symbols_found = paste(unique(exact_hit), collapse = ";"),
            found_in_description = length(desc_hit) > 0,
            symbols_with_description_hit = paste(unique(desc_hit), collapse = ";"),
            aliases_checked = paste(aliases, collapse = ";"),
            stringsAsFactors = FALSE
          )
        })
        do.call(rbind, rows)
      }

      message("Scanning GCT files...")
      report <- do.call(rbind, lapply(files, scan_one_gct))
      fwrite(report, file.path(out_dir, "checkpoint_gene_presence_by_cancer_alias_report.csv"))

      summary_report <- aggregate(found_exact_symbol ~ canonical_gene, report, function(v) sum(v, na.rm = TRUE))
      colnames(summary_report)[2] <- "n_cancers_with_exact_symbol"
      summary_report$n_cancers_total <- length(unique(report$cancer))
      summary_report$presence_fraction <- summary_report$n_cancers_with_exact_symbol / summary_report$n_cancers_total
      summary_report <- summary_report[order(summary_report$presence_fraction, summary_report$canonical_gene), ]
      fwrite(summary_report, file.path(out_dir, "checkpoint_gene_presence_summary.csv"))

      gs_rows <- list()
      for (gs in unique(gmt_df$gene_set)) {
        genes <- unique(gmt_df$gene[gmt_df$gene_set == gs])
        for (ca in unique(report$cancer)) {
          sub <- report[report$cancer == ca & report$canonical_gene %in% genes, ]
          matched <- sub$canonical_gene[sub$found_exact_symbol]
          missing <- setdiff(genes, matched)
          gs_rows[[length(gs_rows) + 1]] <- data.frame(
            cancer = ca,
            gene_set = gs,
            n_input = length(genes),
            n_matched = length(matched),
            n_missing = length(missing),
            matched_genes = paste(matched, collapse = ";"),
            missing_genes = paste(missing, collapse = ";"),
            stringsAsFactors = FALSE
          )
        }
      }
      gs_qc <- do.call(rbind, gs_rows)
      fwrite(gs_qc, file.path(out_dir, "checkpoint_gmt_gene_set_QC_by_cancer_exact_symbol.csv"))

      replacement_suggestions <- data.frame(
        old_symbol = c("PVRL2", "C10orf54", "PDL1", "PD1", "TIM3", "CD39", "CD73", "B7H3", "B7H4"),
        preferred_symbol = c("NECTIN2", "VSIR", "CD274", "PDCD1", "HAVCR2", "ENTPD1", "NT5E", "CD276", "VTCN1"),
        note = c(
          "NECTIN2 is the current HGNC symbol commonly used in expression matrices.",
          "VSIR is the current HGNC symbol for VISTA/C10orf54.",
          "Use formal gene symbol CD274 for PD-L1.",
          "Use formal gene symbol PDCD1 for PD-1.",
          "Use formal gene symbol HAVCR2 for TIM-3.",
          "Use formal gene symbol ENTPD1 for CD39.",
          "Use formal gene symbol NT5E for CD73.",
          "Use formal gene symbol CD276 for B7-H3.",
          "Use formal gene symbol VTCN1 for B7-H4."
        ),
        stringsAsFactors = FALSE
      )
      fwrite(replacement_suggestions, file.path(out_dir, "gene_symbol_alias_replacement_suggestions.csv"))

      message("Done. Check output dir: ", out_dir)
      message("Most important files:")
      message("  checkpoint_gene_presence_summary.csv")
      message("  checkpoint_gene_presence_by_cancer_alias_report.csv")
      message("  checkpoint_gmt_gene_set_QC_by_cancer_exact_symbol.csv")

    })

  })
}

message("============================================================")
message("Unified v11 pipeline finished.")
message("Main output directory: ", MASTER_RESULT_DIR)
message("============================================================")
