# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(20260806)

AUTO_INSTALL <- TRUE
MIN_GENESET_SIZE <- 2
PNG_DPI <- 600
MIN_SUBGROUP_N <- 25
MIN_SUBGROUP_EVENTS <- 5

COL_LOW <- "#1F77B4"
COL_HIGH <- "#D62728"
COL_NONPCR <- "#1F77B4"
COL_PCR <- "#D62728"

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  getwd()
}

PROJECT_DIR <- get_script_dir()
EXPR_FILE <- file.path(PROJECT_DIR, "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz")
SERIES_20078 <- file.path(PROJECT_DIR, "GSE194040-GPL20078_series_matrix.txt.gz")
SERIES_30493 <- file.path(PROJECT_DIR, "GSE194040-GPL30493_series_matrix.txt.gz")
GMT_FILE <- file.path(PROJECT_DIR, "STEP_all_human_gene_sets_merged.gmt")
OUT_DIR <- file.path(PROJECT_DIR, "GSE194040_STEP_ISPY2_results_v1_1")

subdirs <- c(
  "00_QC", "01_scores", "02_global_pCR", "03_subtype_MP",
  "04_treatment_arm", "05_models", "06_matched_control", "07_summary",
  "08_source_data"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (d in subdirs) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUT_DIR, "analysis_log.txt")
sink(log_file, split = TRUE)
cat("GSE194040 STEP I-SPY2 full analysis v1.1\n")
cat("Started: ", format(Sys.time()), "\n", sep = "")
cat("Project directory: ", PROJECT_DIR, "\n\n", sep = "")

cran_pkgs <- c(
  "dplyr", "tidyr", "stringr", "ggplot2", "ggpubr", "broom",
  "pROC", "pheatmap", "data.table", "scales", "janitor"
)
bioc_pkgs <- c("GSVA")

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  if (!AUTO_INSTALL) stop("Missing packages: ", paste(missing, collapse = ", "))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing, dependencies = TRUE)
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(broom)
  library(pROC)
  library(pheatmap)
  library(data.table)
  library(scales)
  library(janitor)
  library(GSVA)
})

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
)

`%||%` <- function(x, y) if (is.null(x)) y else x
safe_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

save_gg <- function(p, out_base, width = 6.5, height = 5.2) {
  grDevices::pdf(
    paste0(out_base, ".pdf"), width = width, height = height,
    onefile = FALSE, useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()

  png_file <- paste0(out_base, ".png")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_file, width = width, height = height, units = "in", res = PNG_DPI)
    print(p)
    grDevices::dev.off()
  } else {
    grDevices::png(png_file, width = width, height = height, units = "in", res = PNG_DPI)
    print(p)
    grDevices::dev.off()
  }
}

save_pheatmap <- function(mat, out_base, width = 10, height = 8, ...) {
  grDevices::pdf(
    paste0(out_base, ".pdf"), width = width, height = height,
    onefile = FALSE, useDingbats = FALSE
  )
  pheatmap::pheatmap(mat, ...)
  grDevices::dev.off()

  png_file <- paste0(out_base, ".png")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_file, width = width, height = height, units = "in", res = PNG_DPI)
  } else {
    grDevices::png(png_file, width = width, height = height, units = "in", res = PNG_DPI)
  }
  pheatmap::pheatmap(mat, ...)
  grDevices::dev.off()
}

required_files <- c(EXPR_FILE, SERIES_20078, SERIES_30493, GMT_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required input files are missing:\n", paste(missing_files, collapse = "\n"))
}

input_inventory <- data.frame(
  file = basename(required_files),
  full_path = required_files,
  exists = file.exists(required_files),
  size_MB = round(file.info(required_files)$size / 1024^2, 3)
)
data.table::fwrite(input_inventory, file.path(OUT_DIR, "00_QC", "input_file_inventory.csv"))

readme_lines <- c(
  "GSE194040 I-SPY2 STEP analysis v1.0",
  "",
  "Primary endpoint: pathological complete response (pCR).",
  "All expression samples are pretreatment breast cancer biopsies.",
  "No OS/PFS/RFS data are present in the supplied GEO series matrices.",
  "No paired post-treatment expression matrix is present.",
  "The combined gene-level expression matrix is already platform/batch adjusted.",
  "Probe-level, RAW, annotation, and linear-adjustment files are not required.",
  "988 expression profiles correspond to 987 patients; the duplicated patient profile is averaged."
)
writeLines(readme_lines, file.path(OUT_DIR, "README_ANALYSIS_SCOPE.txt"))

expr_header <- strsplit(readLines(gzfile(EXPR_FILE), n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
expr_raw <- tryCatch(
  data.table::fread(EXPR_FILE, skip = 1, header = FALSE, data.table = FALSE, showProgress = TRUE),
  error = function(e) {
    message("fread failed; using read.delim fallback: ", conditionMessage(e))
    read.delim(gzfile(EXPR_FILE), skip = 1, header = FALSE, check.names = FALSE, quote = "", comment.char = "")
  }
)

if (ncol(expr_raw) != length(expr_header) + 1) {
  stop(
    "Unexpected expression structure: ", ncol(expr_raw), " columns read, but ",
    length(expr_header), " sample headers were detected."
  )
}
colnames(expr_raw) <- c("gene_symbol", expr_header)

expr_raw$gene_symbol <- toupper(trimws(as.character(expr_raw$gene_symbol)))
expr_raw <- expr_raw[!is.na(expr_raw$gene_symbol) & nzchar(expr_raw$gene_symbol), , drop = FALSE]
gene_symbols <- expr_raw$gene_symbol
expr_raw$gene_symbol <- NULL

non_numeric_cols <- which(!vapply(expr_raw, is.numeric, logical(1)))
if (length(non_numeric_cols) > 0) {
  for (j in non_numeric_cols) expr_raw[[j]] <- suppressWarnings(as.numeric(expr_raw[[j]]))
}
expr_mat <- as.matrix(expr_raw)
rownames(expr_mat) <- gene_symbols
storage.mode(expr_mat) <- "double"
rm(expr_raw, gene_symbols)
gc(verbose = FALSE)

if (anyDuplicated(rownames(expr_mat)) > 0) {
  gene_counts <- table(rownames(expr_mat))
  expr_sum <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE, na.rm = TRUE)
  expr_mat <- expr_sum / as.numeric(gene_counts[rownames(expr_sum)])
}

na_count_before <- sum(is.na(expr_mat))
na_rows <- which(rowSums(is.na(expr_mat)) > 0)
if (length(na_rows) > 0) {
  for (i in na_rows) {
    med_i <- stats::median(expr_mat[i, ], na.rm = TRUE)
    if (!is.finite(med_i)) med_i <- 0
    expr_mat[i, is.na(expr_mat[i, ])] <- med_i
  }
}
if (any(!is.finite(expr_mat))) stop("Expression matrix contains non-finite values after imputation.")

original_sample_ids <- colnames(expr_mat)
patient_ids <- sub("-(GPL16233|GPL20078|GPL30493)$", "", original_sample_ids, ignore.case = TRUE)
duplicate_patient_ids <- unique(patient_ids[duplicated(patient_ids)])

sample_duplicate_audit <- data.frame(
  original_sample_id = original_sample_ids,
  patient_id = patient_ids,
  duplicated_patient = patient_ids %in% duplicate_patient_ids
)
data.table::fwrite(sample_duplicate_audit, file.path(OUT_DIR, "00_QC", "expression_sample_duplicate_audit.csv"))

if (anyDuplicated(patient_ids) > 0) {
  unique_ids <- unique(patient_ids)
  collapsed <- matrix(
    NA_real_, nrow = nrow(expr_mat), ncol = length(unique_ids),
    dimnames = list(rownames(expr_mat), unique_ids)
  )
  for (j in seq_along(unique_ids)) {
    idx <- which(patient_ids == unique_ids[j])
    collapsed[, j] <- if (length(idx) == 1) expr_mat[, idx] else rowMeans(expr_mat[, idx, drop = FALSE])
  }
  expr_mat <- collapsed
  rm(collapsed)
} else {
  colnames(expr_mat) <- patient_ids
}

data.table::fwrite(
  data.frame(
    metric = c(
      "genes", "original_expression_columns", "unique_patients_after_collapse",
      "duplicated_patient_count", "NA_values_imputed"
    ),
    value = c(
      nrow(expr_mat), length(original_sample_ids), ncol(expr_mat),
      length(duplicate_patient_ids), na_count_before
    )
  ),
  file.path(OUT_DIR, "00_QC", "expression_matrix_QC.csv")
)

cat("Expression matrix after patient collapse: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " patients\n", sep = "")

strip_geo_quotes <- function(x) {
  x <- sub('^"', '', x)
  x <- sub('"$', '', x)
  x
}

read_geo_series_metadata <- function(file) {
  lines <- readLines(gzfile(file), warn = FALSE)
  sample_lines <- lines[grepl("^!Sample_", lines)]
  if (length(sample_lines) == 0) stop("No !Sample_ metadata found in ", basename(file))

  parsed <- lapply(sample_lines, function(z) {
    y <- strsplit(z, "\t", fixed = TRUE)[[1]]
    list(key = sub("^!", "", y[1]), values = strip_geo_quotes(y[-1]))
  })

  title_row <- parsed[vapply(parsed, function(z) z$key == "Sample_title", logical(1))]
  geo_row <- parsed[vapply(parsed, function(z) z$key == "Sample_geo_accession", logical(1))]
  platform_row <- parsed[vapply(parsed, function(z) z$key == "Sample_platform_id", logical(1))]
  if (length(title_row) == 0) stop("Sample_title missing in ", basename(file))

  n <- length(title_row[[1]]$values)
  out <- data.frame(
    sample_title = title_row[[1]]$values,
    geo_accession = if (length(geo_row) > 0) geo_row[[1]]$values else rep(NA_character_, n),
    platform = if (length(platform_row) > 0) platform_row[[1]]$values else rep(NA_character_, n),
    stringsAsFactors = FALSE
  )

  char_rows <- parsed[vapply(parsed, function(z) z$key == "Sample_characteristics_ch1", logical(1))]
  for (z in char_rows) {
    if (length(z$values) != n) next
    prefixes <- tolower(trimws(sub(":.*$", "", z$values)))
    if (length(unique(prefixes)) != 1) next
    nm <- janitor::make_clean_names(unique(prefixes))
    values <- trimws(sub("^[^:]+:\\s*", "", z$values))
    if (nm %in% names(out)) nm <- make.unique(c(names(out), nm))[length(names(out)) + 1]
    out[[nm]] <- values
  }
  out
}

meta_20078 <- read_geo_series_metadata(SERIES_20078)
meta_30493 <- read_geo_series_metadata(SERIES_30493)
meta_all <- dplyr::bind_rows(meta_20078, meta_30493)

required_meta <- c("patient_id", "hr", "her2", "mp", "pcr", "arm")
missing_meta <- setdiff(required_meta, names(meta_all))
if (length(missing_meta) > 0) {
  stop("GEO metadata are missing required columns: ", paste(missing_meta, collapse = ", "))
}

meta_all <- meta_all %>%
  dplyr::mutate(
    patient_id = as.character(patient_id),
    hr = suppressWarnings(as.integer(hr)),
    her2 = suppressWarnings(as.integer(her2)),
    mp = suppressWarnings(as.integer(mp)),
    pcr = suppressWarnings(as.integer(pcr)),
    arm = stringr::str_squish(as.character(arm)),
    arm = dplyr::recode(arm, "Paclitaxel + AMG-386" = "Paclitaxel + AMG 386")
  )

meta_duplicate_audit <- meta_all %>%
  dplyr::group_by(patient_id) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::ungroup()
data.table::fwrite(meta_duplicate_audit, file.path(OUT_DIR, "00_QC", "metadata_duplicate_patient_audit.csv"))

clin <- meta_all %>%
  dplyr::distinct(patient_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    pcr_label = factor(ifelse(pcr == 1, "pCR", "non-pCR"), levels = c("non-pCR", "pCR")),
    hr_label = factor(ifelse(hr == 1, "HR+", "HR-"), levels = c("HR-", "HR+")),
    her2_label = factor(ifelse(her2 == 1, "HER2+", "HER2-"), levels = c("HER2-", "HER2+")),
    mp_group = factor(ifelse(mp == 1, "MP2 ultra-high", "MP1 high"), levels = c("MP1 high", "MP2 ultra-high")),
    breast_subtype = dplyr::case_when(
      hr == 1 & her2 == 0 ~ "HR+/HER2-",
      hr == 1 & her2 == 1 ~ "HR+/HER2+",
      hr == 0 & her2 == 1 ~ "HR-/HER2+",
      hr == 0 & her2 == 0 ~ "TNBC",
      TRUE ~ "Unknown"
    ),
    breast_subtype = factor(
      breast_subtype,
      levels = c("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "TNBC", "Unknown")
    ),
    arm_class = dplyr::case_when(
      arm %in% c("Paclitaxel", "Paclitaxel + Trastuzumab") ~ "Control backbone",
      arm == "Paclitaxel + Pembrolizumab" ~ "Immune checkpoint",
      arm == "Paclitaxel + ABT 888 + Carboplatin" ~ "DNA damage",
      arm %in% c("Paclitaxel + Neratinib", "Paclitaxel + Pertuzumab + Trastuzumab", "T-DM1 + Pertuzumab") ~ "HER2 targeted",
      arm %in% c("Paclitaxel + MK-2206", "Paclitaxel + MK-2206 + Trastuzumab") ~ "PI3K/AKT targeted",
      arm %in% c("Paclitaxel + AMG 386", "Paclitaxel + AMG 386 + Trastuzumab") ~ "Anti-angiogenic",
      arm == "Paclitaxel + Ganitumab" ~ "IGF1R targeted",
      arm == "Paclitaxel + Ganetespib" ~ "HSP90 targeted",
      TRUE ~ "Other"
    )
  )

missing_in_clin <- setdiff(colnames(expr_mat), clin$patient_id)
missing_in_expr <- setdiff(clin$patient_id, colnames(expr_mat))
match_qc <- data.frame(
  metric = c("expression_unique_patients", "metadata_unique_patients", "expression_missing_in_metadata", "metadata_missing_in_expression"),
  value = c(ncol(expr_mat), nrow(clin), length(missing_in_clin), length(missing_in_expr))
)
data.table::fwrite(match_qc, file.path(OUT_DIR, "00_QC", "expression_metadata_matching_QC.csv"))
if (length(missing_in_clin) > 0 || length(missing_in_expr) > 0) {
  writeLines(missing_in_clin, file.path(OUT_DIR, "00_QC", "expression_ids_missing_in_metadata.txt"))
  writeLines(missing_in_expr, file.path(OUT_DIR, "00_QC", "metadata_ids_missing_in_expression.txt"))
  stop("Expression/metadata patient matching is incomplete; inspect 00_QC outputs.")
}

clin <- clin[match(colnames(expr_mat), clin$patient_id), , drop = FALSE]
stopifnot(identical(clin$patient_id, colnames(expr_mat)))

cohort_counts <- dplyr::bind_rows(
  clin %>% dplyr::count(pcr_label, name = "n") %>% dplyr::mutate(variable = "pCR", level = as.character(pcr_label)) %>% dplyr::select(variable, level, n),
  clin %>% dplyr::count(breast_subtype, name = "n") %>% dplyr::mutate(variable = "Breast subtype", level = as.character(breast_subtype)) %>% dplyr::select(variable, level, n),
  clin %>% dplyr::count(mp_group, name = "n") %>% dplyr::mutate(variable = "MammaPrint", level = as.character(mp_group)) %>% dplyr::select(variable, level, n),
  clin %>% dplyr::count(arm, name = "n") %>% dplyr::mutate(variable = "Treatment arm", level = arm) %>% dplyr::select(variable, level, n)
)
data.table::fwrite(cohort_counts, file.path(OUT_DIR, "00_QC", "cohort_counts.csv"))
data.table::fwrite(clin, file.path(OUT_DIR, "00_QC", "GSE194040_clinical_metadata_987_patients.csv"))

read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  out <- lapply(lines, function(z) {
    x <- strsplit(z, "\t", fixed = TRUE)[[1]]
    unique(toupper(trimws(x[-c(1, 2)])))
  })
  names(out) <- vapply(lines, function(z) strsplit(z, "\t", fixed = TRUE)[[1]][1], character(1))
  out
}

gene_sets_raw <- read_gmt(GMT_FILE)
short_map <- c(
  "CAF_3gene_score" = "CAF3",
  "CAF_4gene_index" = "CAF4",
  "ECM_score" = "ECM",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "HALLMARK_EMT",
  "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN" = "STEP_UP",
  "STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN" = "STEP_DOWN",
  "HUMAN_EPITHELIAL_SCORE_14" = "EPITHELIAL",
  "HUMAN_MESENCHYMAL_SCORE_18" = "MESENCHYMAL",
  "KOHN_EMT_EPITHELIAL" = "KOHN_EPITHELIAL",
  "KOHN_EMT_MESENCHYMAL" = "KOHN_MESENCHYMAL"
)

set_names <- names(gene_sets_raw)
renamed <- ifelse(set_names %in% names(short_map), unname(short_map[set_names]), safe_name(set_names))
names(gene_sets_raw) <- renamed

gene_sets_present <- lapply(gene_sets_raw, function(gs) intersect(gs, rownames(expr_mat)))
geneset_qc <- data.frame(
  gene_set = names(gene_sets_raw),
  requested_n = vapply(gene_sets_raw, length, integer(1)),
  present_n = vapply(gene_sets_present, length, integer(1)),
  missing_n = vapply(gene_sets_raw, length, integer(1)) - vapply(gene_sets_present, length, integer(1)),
  present_fraction = vapply(gene_sets_present, length, integer(1)) / vapply(gene_sets_raw, length, integer(1))
)
geneset_qc$missing_genes <- vapply(names(gene_sets_raw), function(nm) {
  paste(setdiff(gene_sets_raw[[nm]], gene_sets_present[[nm]]), collapse = ";")
}, character(1))
data.table::fwrite(geneset_qc, file.path(OUT_DIR, "00_QC", "gene_set_expression_overlap.csv"))

gene_sets_present <- gene_sets_present[vapply(gene_sets_present, length, integer(1)) >= MIN_GENESET_SIZE]
required_sets <- c("STEP_UP", "STEP_DOWN", "CAF3")
if (!all(required_sets %in% names(gene_sets_present))) {
  stop("Required STEP/CAF3 gene sets do not have sufficient expression overlap.")
}

score_gene_union <- unique(unlist(gene_sets_present, use.names = FALSE))
gene_z <- t(scale(t(expr_mat[score_gene_union, , drop = FALSE])))
gene_z[!is.finite(gene_z)] <- 0

z_score_matrix <- vapply(names(gene_sets_present), function(nm) {
  colMeans(gene_z[gene_sets_present[[nm]], , drop = FALSE], na.rm = TRUE)
}, numeric(ncol(expr_mat)))
rownames(z_score_matrix) <- colnames(expr_mat)
z_scores <- as.data.frame(z_score_matrix, check.names = FALSE)
z_scores$patient_id <- rownames(z_scores)
z_scores <- z_scores %>% dplyr::select(patient_id, dplyr::everything())
names(z_scores)[-1] <- paste0("Z_", names(z_scores)[-1])
rm(gene_z, z_score_matrix)

run_ssgsea_compatible <- function(expr, sets) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"

  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    ans <- tryCatch({
      available_args <- names(formals(GSVA::ssgseaParam))
      param_args <- list(
        exprData = expr,
        geneSets = sets,
        minSize = MIN_GENESET_SIZE,
        maxSize = Inf,
        alpha = 0.25,
        normalize = TRUE,
        checkNA = "auto",
        use = "na.rm"
      )
      param_args <- param_args[names(param_args) %in% available_args]
      param <- do.call(GSVA::ssgseaParam, param_args)
      GSVA::gsva(param, verbose = FALSE)
    }, error = function(e) {
      message("Current GSVA API attempt failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(ans)) return(ans)
  }

  GSVA::gsva(
    expr = expr,
    gset.idx.list = sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE,
    min.sz = MIN_GENESET_SIZE,
    max.sz = Inf,
    ssgsea.norm = TRUE,
    verbose = FALSE,
    parallel.sz = 1
  )
}

ssgsea_mat <- run_ssgsea_compatible(expr_mat, gene_sets_present)
ssgsea_scores <- as.data.frame(t(ssgsea_mat), check.names = FALSE)
ssgsea_scores$patient_id <- rownames(ssgsea_scores)
ssgsea_scores <- ssgsea_scores %>% dplyr::select(patient_id, dplyr::everything())
names(ssgsea_scores)[-1] <- paste0("SSGSEA_", names(ssgsea_scores)[-1])

analysis_df <- clin %>%
  dplyr::left_join(z_scores, by = "patient_id") %>%
  dplyr::left_join(ssgsea_scores, by = "patient_id") %>%
  dplyr::mutate(
    Z_STEP_BALANCE = Z_STEP_UP - Z_STEP_DOWN,
    SSGSEA_STEP_BALANCE = SSGSEA_STEP_UP - SSGSEA_STEP_DOWN
  )

fit_resid_z <- stats::lm(Z_STEP_UP ~ Z_CAF3, data = analysis_df)
analysis_df$Z_STEP_UP_residual_CAF3 <- stats::residuals(fit_resid_z)
fit_resid_ss <- stats::lm(SSGSEA_STEP_UP ~ SSGSEA_CAF3, data = analysis_df)
analysis_df$SSGSEA_STEP_UP_residual_CAF3 <- stats::residuals(fit_resid_ss)

resid_coeff <- dplyr::bind_rows(
  data.frame(
    method = "Z-mean",
    alpha_intercept = unname(stats::coef(fit_resid_z)[1]),
    beta_CAF3 = unname(stats::coef(fit_resid_z)[2]),
    r_squared = summary(fit_resid_z)$r.squared,
    n = stats::nobs(fit_resid_z)
  ),
  data.frame(
    method = "ssGSEA",
    alpha_intercept = unname(stats::coef(fit_resid_ss)[1]),
    beta_CAF3 = unname(stats::coef(fit_resid_ss)[2]),
    r_squared = summary(fit_resid_ss)$r.squared,
    n = stats::nobs(fit_resid_ss)
  )
)
data.table::fwrite(resid_coeff, file.path(OUT_DIR, "01_scores", "CAF3_residualization_coefficients.csv"))

safe_spearman <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  if (sum(ok) < 3) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

qc_cor_list <- list(
  "Z-mean|raw STEP_UP vs CAF3" = safe_spearman(analysis_df$Z_STEP_UP, analysis_df$Z_CAF3),
  "Z-mean|adjusted STEP_UP vs CAF3" = safe_spearman(analysis_df$Z_STEP_UP_residual_CAF3, analysis_df$Z_CAF3),
  "ssGSEA|raw STEP_UP vs CAF3" = safe_spearman(analysis_df$SSGSEA_STEP_UP, analysis_df$SSGSEA_CAF3),
  "ssGSEA|adjusted STEP_UP vs CAF3" = safe_spearman(analysis_df$SSGSEA_STEP_UP_residual_CAF3, analysis_df$SSGSEA_CAF3)
)
resid_qc <- dplyr::bind_rows(lapply(names(qc_cor_list), function(nm) {
  parts <- strsplit(nm, "\\|", fixed = FALSE)[[1]]
  x <- qc_cor_list[[nm]]
  data.frame(method = parts[1], correlation = parts[2], spearman_rho = x["rho"], p_value = x["p"], n = x["n"])
}))
data.table::fwrite(resid_qc, file.path(OUT_DIR, "01_scores", "CAF3_residualization_correlation_QC.csv"))

score_cols <- grep("^(Z_|SSGSEA_)", names(analysis_df), value = TRUE)
data.table::fwrite(analysis_df, file.path(OUT_DIR, "01_scores", "GSE194040_clinical_STEP_scores_merged_987_patients.csv"))

score_cor <- stats::cor(analysis_df[, score_cols, drop = FALSE], use = "pairwise.complete.obs", method = "spearman")
save_pheatmap(
  score_cor,
  file.path(OUT_DIR, "01_scores", "all_score_spearman_correlation_heatmap"),
  width = 12, height = 10, border_color = NA, fontsize = 7,
  main = "GSE194040 I-SPY2 score correlations"
)

make_caf_scatter <- function(dat, x, y, method_label, adjusted = FALSE) {
  ct <- safe_spearman(dat[[x]], dat[[y]])
  ggplot(dat, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
    labs(
      title = paste0(method_label, ifelse(adjusted, " CAF3-adjusted STEP", " raw STEP vs CAF3")),
      subtitle = paste0("Spearman rho = ", sprintf("%.3f", ct["rho"]), ", p = ", format.pval(ct["p"], digits = 3)),
      x = x, y = y
    )
}

p <- make_caf_scatter(analysis_df, "Z_CAF3", "Z_STEP_UP", "Z-mean", FALSE)
save_gg(p, file.path(OUT_DIR, "01_scores", "CAF3_vs_raw_STEP_UP_Zmean"), 6.2, 5.2)
p <- make_caf_scatter(analysis_df, "Z_CAF3", "Z_STEP_UP_residual_CAF3", "Z-mean", TRUE)
save_gg(p, file.path(OUT_DIR, "01_scores", "CAF3_vs_adjusted_STEP_UP_Zmean"), 6.2, 5.2)
p <- make_caf_scatter(analysis_df, "SSGSEA_CAF3", "SSGSEA_STEP_UP", "ssGSEA", FALSE)
save_gg(p, file.path(OUT_DIR, "01_scores", "CAF3_vs_raw_STEP_UP_ssGSEA"), 6.2, 5.2)
p <- make_caf_scatter(analysis_df, "SSGSEA_CAF3", "SSGSEA_STEP_UP_residual_CAF3", "ssGSEA", TRUE)
save_gg(p, file.path(OUT_DIR, "01_scores", "CAF3_vs_adjusted_STEP_UP_ssGSEA"), 6.2, 5.2)

primary_scores <- intersect(c(
  "Z_STEP_UP", "Z_STEP_DOWN", "Z_CAF3", "Z_STEP_UP_residual_CAF3",
  "SSGSEA_STEP_UP", "SSGSEA_STEP_DOWN", "SSGSEA_CAF3", "SSGSEA_STEP_UP_residual_CAF3"
), names(analysis_df))

secondary_scores <- intersect(c(
  "Z_STEP_BALANCE", "Z_EPITHELIAL", "Z_MESENCHYMAL", "Z_HALLMARK_EMT",
  "SSGSEA_STEP_BALANCE", "SSGSEA_EPITHELIAL", "SSGSEA_MESENCHYMAL", "SSGSEA_HALLMARK_EMT"
), names(analysis_df))
all_response_scores <- unique(c(primary_scores, secondary_scores))

response_tests <- dplyr::bind_rows(lapply(all_response_scores, function(sc) {
  d <- analysis_df[, c("pcr", "pcr_label", sc), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  x1 <- d[d$pcr == 1, sc]
  x0 <- d[d$pcr == 0, sc]
  wt <- suppressWarnings(stats::wilcox.test(x1, x0, exact = FALSE))
  data.frame(
    score = sc,
    n = nrow(d),
    pCR_n = length(x1),
    non_pCR_n = length(x0),
    pCR_median = stats::median(x1),
    non_pCR_median = stats::median(x0),
    median_difference_pCR_minus_non_pCR = stats::median(x1) - stats::median(x0),
    wilcoxon_p = wt$p.value
  )
})) %>%
  dplyr::mutate(FDR = p.adjust(wilcoxon_p, method = "BH"))
data.table::fwrite(response_tests, file.path(OUT_DIR, "02_global_pCR", "global_pCR_group_comparisons.csv"))

plot_two_group <- function(dat, yvar, title, subtitle = NULL) {
  ggplot(dat, aes(x = pcr_label, y = .data[[yvar]], fill = pcr_label)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
    geom_jitter(width = 0.14, size = 1.2, alpha = 0.45) +
    scale_fill_manual(values = c("non-pCR" = COL_NONPCR, "pCR" = COL_PCR)) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
    labs(title = title, subtitle = subtitle, x = NULL, y = yvar) +
    theme(legend.position = "none")
}

for (sc in primary_scores) {
  p <- plot_two_group(analysis_df, sc, paste0("GSE194040 pCR: ", sc))
  save_gg(p, file.path(OUT_DIR, "02_global_pCR", paste0("pCR_boxplot_", sc)), 5.8, 5.2)
}

make_roc_table_plot <- function(dat, outcome, scores, title, out_base) {
  roc_list <- list()
  auc_list <- list()
  for (sc in scores) {
    ok <- stats::complete.cases(dat[[outcome]], dat[[sc]])
    if (sum(ok) < 20 || length(unique(dat[[outcome]][ok])) < 2) next
    rr <- pROC::roc(
      response = dat[[outcome]][ok], predictor = dat[[sc]][ok],
      levels = c(0, 1), direction = "<", quiet = TRUE
    )
    ci <- suppressWarnings(pROC::ci.auc(rr))
    roc_list[[sc]] <- rr
    auc_list[[sc]] <- data.frame(
      score = sc,
      AUC_higher_score_predicts_pCR = as.numeric(pROC::auc(rr)),
      AUC_low = as.numeric(ci[1]),
      AUC_high = as.numeric(ci[3]),
      n = sum(ok)
    )
  }
  auc_df <- dplyr::bind_rows(auc_list)
  data.table::fwrite(auc_df, paste0(out_base, "_AUC.csv"))

  if (length(roc_list) > 0) {
    display_names <- paste0(
      names(roc_list), " (AUC=",
      sprintf("%.3f", auc_df$AUC_higher_score_predicts_pCR[match(names(roc_list), auc_df$score)]), ")"
    )
    names(roc_list) <- display_names
    p <- pROC::ggroc(roc_list, legacy.axes = TRUE, linewidth = 0.8) +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      coord_equal() +
      labs(title = title, x = "1 - Specificity", y = "Sensitivity", color = "Score") +
      theme(legend.position = "right", legend.text = element_text(size = 7))
    save_gg(p, out_base, 8.5, 6.6)
  }
  auc_df
}

z_roc_scores <- primary_scores[grepl("^Z_", primary_scores)]
ss_roc_scores <- primary_scores[grepl("^SSGSEA_", primary_scores)]
auc_z <- make_roc_table_plot(
  analysis_df, "pcr", z_roc_scores,
  "GSE194040 pCR ROC: Z-mean scores",
  file.path(OUT_DIR, "02_global_pCR", "ROC_pCR_Zmean_primary_scores")
)
auc_ss <- make_roc_table_plot(
  analysis_df, "pcr", ss_roc_scores,
  "GSE194040 pCR ROC: ssGSEA scores",
  file.path(OUT_DIR, "02_global_pCR", "ROC_pCR_ssGSEA_primary_scores")
)

high_low_results <- list()
for (sc in c("Z_STEP_UP_residual_CAF3", "SSGSEA_STEP_UP_residual_CAF3")) {
  if (!sc %in% names(analysis_df)) next
  d <- analysis_df %>%
    dplyr::filter(!is.na(.data[[sc]]), !is.na(pcr)) %>%
    dplyr::mutate(score_group = factor(ifelse(.data[[sc]] > stats::median(.data[[sc]]), "High", "Low"), levels = c("Low", "High")))
  tab <- table(d$score_group, d$pcr)
  p_fisher <- if (all(dim(tab) == c(2, 2))) stats::fisher.test(tab)$p.value else NA_real_
  rate <- d %>%
    dplyr::group_by(score_group) %>%
    dplyr::summarise(n = dplyr::n(), pCR_n = sum(pcr == 1), pCR_rate = mean(pcr == 1), .groups = "drop") %>%
    dplyr::mutate(score = sc, fisher_p = p_fisher)
  high_low_results[[sc]] <- rate
  p <- ggplot(rate, aes(x = score_group, y = pCR_rate, fill = score_group)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(pCR_n, "/", n, "\n", scales::percent(pCR_rate, accuracy = 0.1))), vjust = -0.25, size = 3.5) +
    scale_fill_manual(values = c("Low" = COL_LOW, "High" = COL_HIGH)) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, min(1, max(rate$pCR_rate) + 0.16))) +
    labs(title = paste0("pCR rate by ", sc), subtitle = paste0("Fisher p = ", format.pval(p_fisher, digits = 3)), x = NULL, y = "pCR rate") +
    theme(legend.position = "none")
  save_gg(p, file.path(OUT_DIR, "02_global_pCR", paste0("pCR_rate_High_Low_", sc)), 5.8, 5.2)
}
data.table::fwrite(dplyr::bind_rows(high_low_results), file.path(OUT_DIR, "02_global_pCR", "pCR_rate_by_CAF_adjusted_STEP_high_low.csv"))

fit_logistic_score <- function(dat, score, covars = character(), model_name = "") {
  vars <- unique(c("pcr", score, covars))
  d <- dat[, vars, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 20 || sum(d$pcr == 1) < 5 || sum(d$pcr == 0) < 5) return(NULL)
  if (stats::sd(d[[score]]) == 0) return(NULL)

  d$score_std <- as.numeric(scale(d[[score]]))
  f <- stats::reformulate(c("score_std", covars), response = "pcr")
  warnings_seen <- character()
  fit <- withCallingHandlers(
    tryCatch(stats::glm(f, data = d, family = stats::binomial()), error = function(e) NULL),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  if (!"score_std" %in% rownames(co)) return(NULL)
  beta <- co["score_std", "Estimate"]
  se <- co["score_std", "Std. Error"]
  pval <- co["score_std", "Pr(>|z|)"]
  data.frame(
    score = score,
    model = model_name,
    covariates = ifelse(length(covars) == 0, "None", paste(covars, collapse = "+")),
    n = nrow(d),
    pCR_events = sum(d$pcr == 1),
    OR_per_SD = exp(beta),
    CI_low = exp(beta - 1.96 * se),
    CI_high = exp(beta + 1.96 * se),
    p_value = pval,
    converged = isTRUE(fit$converged),
    warning = paste(unique(warnings_seen), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

model_definitions <- list(
  "Univariable" = character(),
  "HR_HER2_MP_adjusted" = c("hr_label", "her2_label", "mp_group"),
  "Treatment_arm_adjusted" = c("arm"),
  "Full_HR_HER2_MP_arm_adjusted" = c("hr_label", "her2_label", "mp_group", "arm")
)

logistic_results <- dplyr::bind_rows(lapply(names(model_definitions), function(mdl) {
  dplyr::bind_rows(lapply(all_response_scores, function(sc) {
    fit_logistic_score(analysis_df, sc, model_definitions[[mdl]], mdl)
  }))
})) %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()
data.table::fwrite(logistic_results, file.path(OUT_DIR, "05_models", "global_pCR_logistic_regression_all_scores.csv"))

plot_or_forest <- function(tbl, model_name, scores, title, out_base, height = 5.5) {
  dd <- tbl %>%
    dplyr::filter(model == model_name, score %in% scores) %>%
    dplyr::mutate(score = factor(score, levels = rev(scores)))
  if (nrow(dd) == 0) return(invisible(NULL))
  p <- ggplot(dd, aes(x = OR_per_SD, y = score)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.18, orientation = "y") +
    geom_point(size = 2.5) +
    scale_x_log10() +
    labs(title = title, x = "OR for pCR per 1-SD score increase (log scale)", y = NULL)
  save_gg(p, out_base, 7.8, height)
}

for (mdl in names(model_definitions)) {
  plot_or_forest(
    logistic_results, mdl, primary_scores,
    paste0("GSE194040 pCR logistic model: ", mdl),
    file.path(OUT_DIR, "05_models", paste0("forest_global_pCR_", safe_name(mdl))),
    height = 6.0
  )
}

pcr_summary_by <- function(dat, variable) {
  dat %>%
    dplyr::group_by(.data[[variable]]) %>%
    dplyr::summarise(
      n = dplyr::n(), pCR_n = sum(pcr == 1), non_pCR_n = sum(pcr == 0),
      pCR_rate = mean(pcr == 1), .groups = "drop"
    ) %>%
    dplyr::rename(level = dplyr::all_of(variable)) %>%
    dplyr::mutate(variable = variable, level = as.character(level)) %>%
    dplyr::select(variable, level, dplyr::everything())
}

subtype_pcr_summary <- pcr_summary_by(analysis_df, "breast_subtype")
mp_pcr_summary <- pcr_summary_by(analysis_df, "mp_group")
data.table::fwrite(subtype_pcr_summary, file.path(OUT_DIR, "03_subtype_MP", "pCR_rate_by_breast_subtype.csv"))
data.table::fwrite(mp_pcr_summary, file.path(OUT_DIR, "03_subtype_MP", "pCR_rate_by_MammaPrint.csv"))

p <- ggplot(subtype_pcr_summary, aes(x = level, y = pCR_rate, fill = level)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(pCR_n, "/", n, "\n", scales::percent(pCR_rate, accuracy = 0.1))), vjust = -0.2, size = 3.3) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, min(1, max(subtype_pcr_summary$pCR_rate) + 0.15))) +
  labs(title = "GSE194040 pCR rate by HR/HER2 subtype", x = NULL, y = "pCR rate") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_gg(p, file.path(OUT_DIR, "03_subtype_MP", "pCR_rate_by_breast_subtype"), 7.2, 5.4)

focus_scores <- intersect(c(
  "Z_STEP_UP", "Z_STEP_UP_residual_CAF3",
  "SSGSEA_STEP_UP", "SSGSEA_STEP_UP_residual_CAF3"
), names(analysis_df))

SUBTYPE_SOURCE_DIR <- file.path(OUT_DIR, "08_source_data", "03_subtype_pCR")
dir.create(SUBTYPE_SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

export_subtype_pcr_source <- function(dat, score_var, out_dir) {

  required_cols <- c("patient_id", "pcr", "pcr_label", "breast_subtype", score_var)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    warning(
      "Cannot export subtype pCR source data for ", score_var,
      "; missing columns: ", paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }

  keep_cols <- unique(c(
    "patient_id",
    "hr", "her2", "hr_label", "her2_label",
    "breast_subtype",
    "pcr", "pcr_label",
    "mp", "mp_group",
    "arm", "arm_class",
    if (grepl("^SSGSEA_", score_var)) c(
      "SSGSEA_STEP_UP",
      "SSGSEA_CAF3",
      "SSGSEA_STEP_UP_residual_CAF3"
    ) else c(
      "Z_STEP_UP",
      "Z_CAF3",
      "Z_STEP_UP_residual_CAF3"
    ),
    score_var
  ))
  keep_cols <- intersect(keep_cols, names(dat))

  source_dat <- dat[, keep_cols, drop = FALSE]

  use <- stats::complete.cases(
    source_dat[, intersect(c("breast_subtype", "pcr_label", score_var), names(source_dat)), drop = FALSE]
  )
  source_dat <- source_dat[use, , drop = FALSE]
  source_dat <- source_dat[is.finite(suppressWarnings(as.numeric(source_dat[[score_var]]))), , drop = FALSE]

  source_dat$score_name <- score_var
  source_dat$score_value <- suppressWarnings(as.numeric(source_dat[[score_var]]))

  data.table::fwrite(
    source_dat,
    file.path(out_dir, paste0("pCR_by_subtype_", safe_name(score_var), "_source_data.csv"))
  )

  subtype_levels <- unique(as.character(source_dat$breast_subtype))
  subtype_levels <- subtype_levels[!is.na(subtype_levels) & nzchar(subtype_levels)]

  stat_list <- lapply(subtype_levels, function(lv) {

    dd <- source_dat[as.character(source_dat$breast_subtype) == lv, , drop = FALSE]
    dd$pcr_label <- droplevels(factor(dd$pcr_label, levels = c("non-pCR", "pCR")))

    if (nrow(dd) == 0 || length(unique(as.character(dd$pcr_label))) < 2) {
      return(data.frame(
        score = score_var,
        breast_subtype = lv,
        n = nrow(dd),
        non_pCR_n = sum(as.character(dd$pcr_label) == "non-pCR", na.rm = TRUE),
        pCR_n = sum(as.character(dd$pcr_label) == "pCR", na.rm = TRUE),
        non_pCR_median = NA_real_,
        pCR_median = NA_real_,
        non_pCR_IQR = NA_real_,
        pCR_IQR = NA_real_,
        median_difference_pCR_minus_non_pCR = NA_real_,
        wilcoxon_W = NA_real_,
        wilcoxon_p = NA_real_,
        p_format_for_plot = NA_character_,
        p_signif = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    x_non <- dd$score_value[as.character(dd$pcr_label) == "non-pCR"]
    x_pcr <- dd$score_value[as.character(dd$pcr_label) == "pCR"]

    cmp <- ggpubr::compare_means(
      formula = stats::reformulate("pcr_label", response = "score_value"),
      data = dd,
      method = "wilcox.test"
    )

    wt <- suppressWarnings(
      stats::wilcox.test(score_value ~ pcr_label, data = dd)
    )

    data.frame(
      score = score_var,
      breast_subtype = lv,
      n = nrow(dd),
      non_pCR_n = length(x_non),
      pCR_n = length(x_pcr),
      non_pCR_median = stats::median(x_non, na.rm = TRUE),
      pCR_median = stats::median(x_pcr, na.rm = TRUE),
      non_pCR_IQR = stats::IQR(x_non, na.rm = TRUE),
      pCR_IQR = stats::IQR(x_pcr, na.rm = TRUE),
      median_difference_pCR_minus_non_pCR =
        stats::median(x_pcr, na.rm = TRUE) - stats::median(x_non, na.rm = TRUE),
      wilcoxon_W = unname(wt$statistic),
      wilcoxon_p = cmp$p[1],
      p_format_for_plot = cmp$p.format[1],
      p_signif = cmp$p.signif[1],
      stringsAsFactors = FALSE
    )
  })

  stat_df <- dplyr::bind_rows(stat_list) %>%
    dplyr::mutate(
      FDR_across_subtypes_for_this_score = stats::p.adjust(wilcoxon_p, method = "BH"),
      test = "Wilcoxon rank-sum test",
      figure_method = 'ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format")'
    )

  data.table::fwrite(
    stat_df,
    file.path(out_dir, paste0("pCR_by_subtype_", safe_name(score_var), "_wilcoxon_statistics.csv"))
  )

  stat_df
}

subtype_source_stats_list <- list()

for (sc in focus_scores) {

  subtype_source_stats_list[[sc]] <- export_subtype_pcr_source(
    analysis_df, sc, SUBTYPE_SOURCE_DIR
  )
  p <- ggplot(analysis_df, aes(x = pcr_label, y = .data[[sc]], fill = pcr_label)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
    geom_jitter(width = 0.13, size = 0.9, alpha = 0.4) +
    scale_fill_manual(values = c("non-pCR" = COL_NONPCR, "pCR" = COL_PCR)) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.2) +
    facet_wrap(~ breast_subtype, scales = "free_y", ncol = 2) +
    labs(title = paste0("pCR and ", sc, " within breast subtypes"), x = NULL, y = sc) +
    theme(legend.position = "none")
  save_gg(p, file.path(OUT_DIR, "03_subtype_MP", paste0("pCR_by_subtype_", sc)), 9.2, 7.4)

  p <- ggplot(analysis_df, aes(x = mp_group, y = .data[[sc]], fill = mp_group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
    geom_jitter(width = 0.13, size = 0.9, alpha = 0.35) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5) +
    labs(title = paste0(sc, " by MammaPrint group"), x = NULL, y = sc) +
    theme(legend.position = "none")
  save_gg(p, file.path(OUT_DIR, "03_subtype_MP", paste0("score_by_MammaPrint_", sc)), 6.3, 5.3)
}

subtype_source_stats_all <- dplyr::bind_rows(subtype_source_stats_list)
if (nrow(subtype_source_stats_all) > 0) {
  data.table::fwrite(
    subtype_source_stats_all,
    file.path(
      SUBTYPE_SOURCE_DIR,
      "pCR_by_subtype_ALL_focus_scores_wilcoxon_statistics.csv"
    )
  )
}

fit_subgroup_score <- function(dat, subgroup_var, subgroup_level, score) {
  d <- dat[as.character(dat[[subgroup_var]]) == subgroup_level, c("pcr", score), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < MIN_SUBGROUP_N || sum(d$pcr == 1) < MIN_SUBGROUP_EVENTS || sum(d$pcr == 0) < MIN_SUBGROUP_EVENTS) return(NULL)
  if (stats::sd(d[[score]]) == 0) return(NULL)
  d$score_std <- as.numeric(scale(d[[score]]))
  fit <- tryCatch(stats::glm(pcr ~ score_std, data = d, family = stats::binomial()), error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  co <- summary(fit)$coefficients["score_std", ]
  data.frame(
    subgroup_variable = subgroup_var,
    subgroup = subgroup_level,
    score = score,
    n = nrow(d), pCR_events = sum(d$pcr == 1),
    OR_per_SD = exp(co["Estimate"]),
    CI_low = exp(co["Estimate"] - 1.96 * co["Std. Error"]),
    CI_high = exp(co["Estimate"] + 1.96 * co["Std. Error"]),
    p_value = co["Pr(>|z|)"]
  )
}

subgroup_results <- dplyr::bind_rows(
  lapply(levels(droplevels(analysis_df$breast_subtype)), function(lv) {
    dplyr::bind_rows(lapply(focus_scores, function(sc) fit_subgroup_score(analysis_df, "breast_subtype", lv, sc)))
  }),
  lapply(levels(droplevels(analysis_df$mp_group)), function(lv) {
    dplyr::bind_rows(lapply(focus_scores, function(sc) fit_subgroup_score(analysis_df, "mp_group", lv, sc)))
  })
) %>%
  dplyr::group_by(subgroup_variable, score) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()
data.table::fwrite(subgroup_results, file.path(OUT_DIR, "03_subtype_MP", "subgroup_pCR_logistic_results.csv"))

for (sc in focus_scores) {
  dd <- subgroup_results %>%
    dplyr::filter(score == sc, subgroup_variable == "breast_subtype") %>%
    dplyr::mutate(subgroup = factor(subgroup, levels = rev(unique(subgroup))))
  if (nrow(dd) > 0) {
    p <- ggplot(dd, aes(x = OR_per_SD, y = subgroup)) +
      geom_vline(xintercept = 1, linetype = 2) +
      geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.18, orientation = "y") +
      geom_point(size = 2.5) +
      scale_x_log10() +
      labs(title = paste0("Subtype-specific pCR association: ", sc), x = "OR per 1-SD increase", y = NULL)
    save_gg(p, file.path(OUT_DIR, "03_subtype_MP", paste0("forest_subtype_", sc)), 7.2, 5.2)
  }
}

arm_pcr_summary <- analysis_df %>%
  dplyr::group_by(arm, arm_class) %>%
  dplyr::summarise(
    n = dplyr::n(), pCR_n = sum(pcr == 1), non_pCR_n = sum(pcr == 0),
    pCR_rate = mean(pcr == 1), .groups = "drop"
  ) %>%
  dplyr::arrange(pCR_rate)
data.table::fwrite(arm_pcr_summary, file.path(OUT_DIR, "04_treatment_arm", "pCR_rate_by_treatment_arm.csv"))

p <- ggplot(arm_pcr_summary, aes(x = pCR_rate, y = reorder(stringr::str_wrap(arm, 32), pCR_rate), fill = arm_class)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = paste0(pCR_n, "/", n)), hjust = -0.12, size = 3.0) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, min(1, max(arm_pcr_summary$pCR_rate) + 0.12))) +
  labs(title = "GSE194040 pCR rate by treatment arm", x = "pCR rate", y = NULL, fill = "Arm class")
save_gg(p, file.path(OUT_DIR, "04_treatment_arm", "pCR_rate_by_treatment_arm"), 10.5, 7.8)

for (sc in c("Z_STEP_UP_residual_CAF3", "SSGSEA_STEP_UP_residual_CAF3")) {
  if (!sc %in% names(analysis_df)) next
  plot_dat <- analysis_df %>%
    dplyr::mutate(arm_plot = stringr::str_wrap(arm, 28))
  p <- ggplot(plot_dat, aes(x = pcr_label, y = .data[[sc]], fill = pcr_label)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
    geom_jitter(width = 0.12, size = 0.65, alpha = 0.35) +
    scale_fill_manual(values = c("non-pCR" = COL_NONPCR, "pCR" = COL_PCR)) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.5) +
    facet_wrap(~ arm_plot, scales = "free_y", ncol = 4) +
    labs(title = paste0("pCR and ", sc, " within treatment arms"), x = NULL, y = sc) +
    theme(legend.position = "none", strip.text = element_text(size = 8))
  save_gg(p, file.path(OUT_DIR, "04_treatment_arm", paste0("pCR_by_arm_", sc)), 13.5, 11.0)
}

arm_subgroup_results <- dplyr::bind_rows(lapply(unique(analysis_df$arm), function(lv) {
  dplyr::bind_rows(lapply(focus_scores, function(sc) fit_subgroup_score(analysis_df, "arm", lv, sc)))
})) %>%
  dplyr::group_by(score) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()
data.table::fwrite(arm_subgroup_results, file.path(OUT_DIR, "04_treatment_arm", "treatment_arm_specific_pCR_logistic_results.csv"))

for (sc in focus_scores) {
  dd <- arm_subgroup_results %>%
    dplyr::filter(score == sc) %>%
    dplyr::arrange(OR_per_SD) %>%
    dplyr::mutate(arm_label = factor(stringr::str_wrap(subgroup, 34), levels = stringr::str_wrap(subgroup, 34)))
  if (nrow(dd) > 0) {
    p <- ggplot(dd, aes(x = OR_per_SD, y = arm_label)) +
      geom_vline(xintercept = 1, linetype = 2) +
      geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.17, orientation = "y") +
      geom_point(aes(size = n), alpha = 0.85) +
      scale_x_log10() +
      labs(title = paste0("Treatment-arm-specific pCR association: ", sc), x = "OR per 1-SD increase", y = NULL, size = "n")
    save_gg(p, file.path(OUT_DIR, "04_treatment_arm", paste0("forest_treatment_arm_", sc)), 9.2, max(5.5, 0.48 * nrow(dd) + 1.8))
  }
}

fit_interaction_lrt <- function(dat, score, modifier) {
  d <- dat[, c("pcr", score, modifier), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 40 || length(unique(d$pcr)) < 2 || length(unique(d[[modifier]])) < 2) return(NULL)
  if (stats::sd(d[[score]]) == 0) return(NULL)
  d$score_std <- as.numeric(scale(d[[score]]))
  f0 <- stats::reformulate(c("score_std", modifier), response = "pcr")
  f1 <- stats::as.formula(paste("pcr ~ score_std *", modifier))
  m0 <- tryCatch(stats::glm(f0, data = d, family = stats::binomial()), error = function(e) NULL)
  m1 <- tryCatch(stats::glm(f1, data = d, family = stats::binomial()), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(NULL)
  lrt <- tryCatch(stats::anova(m0, m1, test = "LRT"), error = function(e) NULL)
  if (is.null(lrt)) return(NULL)
  data.frame(
    score = score, modifier = modifier, n = nrow(d),
    modifier_levels = length(unique(d[[modifier]])),
    interaction_LRT_p = lrt$`Pr(>Chi)`[2]
  )
}

interaction_scores <- focus_scores
interaction_modifiers <- c("breast_subtype", "mp_group", "hr_label", "her2_label", "arm", "arm_class")
interaction_results <- dplyr::bind_rows(lapply(interaction_scores, function(sc) {
  dplyr::bind_rows(lapply(interaction_modifiers, function(mod) fit_interaction_lrt(analysis_df, sc, mod)))
})) %>%
  dplyr::mutate(FDR = p.adjust(interaction_LRT_p, method = "BH"))
data.table::fwrite(interaction_results, file.path(OUT_DIR, "05_models", "global_pCR_interaction_tests.csv"))

comparison_pairs <- data.frame(
  comparison = c(
    "ABT888_Carboplatin_vs_Paclitaxel_HER2neg",
    "AMG386_vs_Paclitaxel_HER2neg",
    "Ganitumab_vs_Paclitaxel_HER2neg",
    "Ganetespib_vs_Paclitaxel_HER2neg",
    "MK2206_vs_Paclitaxel_HER2neg",
    "Pembrolizumab_vs_Paclitaxel_HER2neg",
    "Neratinib_vs_Paclitaxel_HER2neg",
    "MK2206_Trastuzumab_vs_Paclitaxel_Trastuzumab_HER2pos",
    "AMG386_Trastuzumab_vs_Paclitaxel_Trastuzumab_HER2pos",
    "Pertuzumab_Trastuzumab_vs_Paclitaxel_Trastuzumab_HER2pos",
    "TDM1_Pertuzumab_vs_Paclitaxel_Trastuzumab_HER2pos",
    "Neratinib_vs_Paclitaxel_Trastuzumab_HER2pos"
  ),
  experimental_arm = c(
    "Paclitaxel + ABT 888 + Carboplatin",
    "Paclitaxel + AMG 386",
    "Paclitaxel + Ganitumab",
    "Paclitaxel + Ganetespib",
    "Paclitaxel + MK-2206",
    "Paclitaxel + Pembrolizumab",
    "Paclitaxel + Neratinib",
    "Paclitaxel + MK-2206 + Trastuzumab",
    "Paclitaxel + AMG 386 + Trastuzumab",
    "Paclitaxel + Pertuzumab + Trastuzumab",
    "T-DM1 + Pertuzumab",
    "Paclitaxel + Neratinib"
  ),
  control_arm = c(
    rep("Paclitaxel", 7),
    rep("Paclitaxel + Trastuzumab", 5)
  ),
  her2_filter = c(rep(0, 7), rep(1, 5)),
  stringsAsFactors = FALSE
)
data.table::fwrite(comparison_pairs, file.path(OUT_DIR, "06_matched_control", "matched_control_comparison_definitions.csv"))

fit_matched_control_interaction <- function(dat, pair_row, score) {
  d <- dat %>%
    dplyr::filter(
      arm %in% c(pair_row$experimental_arm, pair_row$control_arm),
      her2 == pair_row$her2_filter
    ) %>%
    dplyr::select(pcr, arm, hr_label, mp_group, dplyr::all_of(score)) %>%
    tidyr::drop_na()

  if (nrow(d) < 40 || sum(d$pcr == 1) < 8 || sum(d$pcr == 0) < 8) return(NULL)
  if (!all(c(pair_row$experimental_arm, pair_row$control_arm) %in% unique(d$arm))) return(NULL)
  if (stats::sd(d[[score]]) == 0) return(NULL)

  d$treatment <- factor(ifelse(d$arm == pair_row$experimental_arm, "Experimental", "Control"), levels = c("Control", "Experimental"))
  d$score_std <- as.numeric(scale(d[[score]]))

  covars <- character()
  if (length(unique(d$hr_label)) > 1) covars <- c(covars, "hr_label")
  if (length(unique(d$mp_group)) > 1) covars <- c(covars, "mp_group")
  base_terms <- c("score_std", "treatment", covars)
  int_terms <- c("score_std * treatment", covars)
  f0 <- stats::reformulate(base_terms, response = "pcr")
  f1 <- stats::as.formula(paste("pcr ~", paste(int_terms, collapse = " + ")))

  m0 <- tryCatch(stats::glm(f0, data = d, family = stats::binomial()), error = function(e) NULL)
  m1 <- tryCatch(stats::glm(f1, data = d, family = stats::binomial()), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(NULL)
  lrt <- tryCatch(stats::anova(m0, m1, test = "LRT"), error = function(e) NULL)
  if (is.null(lrt)) return(NULL)

  co <- summary(m1)$coefficients
  int_term <- grep("score_std:treatmentExperimental|treatmentExperimental:score_std", rownames(co), value = TRUE)
  if (length(int_term) != 1) return(NULL)
  beta <- co[int_term, "Estimate"]
  se <- co[int_term, "Std. Error"]

  data.frame(
    comparison = pair_row$comparison,
    experimental_arm = pair_row$experimental_arm,
    control_arm = pair_row$control_arm,
    HER2_subset = ifelse(pair_row$her2_filter == 1, "HER2+", "HER2-"),
    score = score,
    n = nrow(d),
    experimental_n = sum(d$treatment == "Experimental"),
    control_n = sum(d$treatment == "Control"),
    pCR_events = sum(d$pcr == 1),
    interaction_OR_ratio = exp(beta),
    CI_low = exp(beta - 1.96 * se),
    CI_high = exp(beta + 1.96 * se),
    interaction_Wald_p = co[int_term, "Pr(>|z|)"],
    interaction_LRT_p = lrt$`Pr(>Chi)`[2],
    adjusted_covariates = ifelse(length(covars) == 0, "None", paste(covars, collapse = "+"))
  )
}

matched_control_results <- dplyr::bind_rows(lapply(seq_len(nrow(comparison_pairs)), function(i) {
  pair_row <- comparison_pairs[i, , drop = FALSE]
  dplyr::bind_rows(lapply(focus_scores, function(sc) fit_matched_control_interaction(analysis_df, pair_row, sc)))
})) %>%
  dplyr::group_by(score) %>%
  dplyr::mutate(FDR_LRT = p.adjust(interaction_LRT_p, method = "BH")) %>%
  dplyr::ungroup()
data.table::fwrite(matched_control_results, file.path(OUT_DIR, "06_matched_control", "matched_control_score_treatment_interactions.csv"))

for (sc in focus_scores) {
  dd <- matched_control_results %>%
    dplyr::filter(score == sc) %>%
    dplyr::arrange(interaction_OR_ratio) %>%
    dplyr::mutate(comparison_label = factor(stringr::str_wrap(comparison, 38), levels = stringr::str_wrap(comparison, 38)))
  if (nrow(dd) > 0) {
    p <- ggplot(dd, aes(x = interaction_OR_ratio, y = comparison_label)) +
      geom_vline(xintercept = 1, linetype = 2) +
      geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.17, orientation = "y") +
      geom_point(aes(size = n), alpha = 0.85) +
      scale_x_log10() +
      labs(
        title = paste0("Score-by-treatment interaction: ", sc),
        subtitle = "OR ratio >1: score-pCR association is more positive in experimental vs control arm",
        x = "Interaction OR ratio (log scale)", y = NULL, size = "n"
      )
    save_gg(p, file.path(OUT_DIR, "06_matched_control", paste0("forest_matched_control_interaction_", sc)), 10.0, max(6.0, 0.52 * nrow(dd) + 2.0))
  }
}

key_summary <- dplyr::bind_rows(
  response_tests %>%
    dplyr::filter(score %in% primary_scores) %>%
    dplyr::transmute(
      analysis = "Global pCR Wilcoxon",
      score, estimate = median_difference_pCR_minus_non_pCR,
      CI_low = NA_real_, CI_high = NA_real_, p_value = wilcoxon_p, FDR
    ),
  logistic_results %>%
    dplyr::filter(score %in% primary_scores, model == "Full_HR_HER2_MP_arm_adjusted") %>%
    dplyr::transmute(
      analysis = "Global pCR full-adjusted logistic",
      score, estimate = OR_per_SD, CI_low, CI_high, p_value, FDR
    )
)
data.table::fwrite(key_summary, file.path(OUT_DIR, "07_summary", "key_primary_results_summary.csv"))

summary_text <- c(
  "GSE194040 I-SPY2 STEP analysis completed.",
  paste0("Patients analyzed: ", nrow(analysis_df)),
  paste0("pCR: ", sum(analysis_df$pcr == 1), "; non-pCR: ", sum(analysis_df$pcr == 0)),
  "Primary endpoint: pCR.",
  "No survival or paired treatment-change analysis was performed because those data are absent from the supplied files.",
  "Interpret matched-control and treatment-arm interaction results as exploratory, particularly where confidence intervals are wide.",
  "Continuous-score logistic models are primary; median-split results are visualization-only."
)
writeLines(summary_text, file.path(OUT_DIR, "07_summary", "analysis_completion_summary.txt"))

capture.output(sessionInfo(), file = file.path(OUT_DIR, "07_summary", "sessionInfo.txt"))
cat("\nCompleted: ", format(Sys.time()), "\n", sep = "")
cat("Results directory: ", OUT_DIR, "\n", sep = "")
sink()

message("GSE194040 STEP I-SPY2 analysis completed. Results: ", OUT_DIR)
