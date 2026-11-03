# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(123)

cran_pkgs <- c("ggplot2", "pheatmap", "RColorBrewer", "survival")
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
})

base_dir <- getwd()

exp_dir <- file.path(base_dir, "TCGA_201906cBioP", "exp_gct")
gmt_file <- file.path(base_dir, "gene set", "STEP_all_human_gene_sets_merged.gmt")
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
  "11_CAF_residualization_details"
)
for (d in subdirs) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)

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
      genes <- unique(toupper(parts[-c(1, 2)]))
      genes <- genes[genes != ""]
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

score_gene_sets_zmean <- function(expr_log, gene_sets) {
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

    if (length(matched) >= 1) {
      scores[, nm] <- colMeans(z[matched, , drop = FALSE], na.rm = TRUE)
    }

    qc <- rbind(qc, data.frame(
      gene_set = nm,
      n_genes = length(genes),
      n_matched = length(matched),
      n_missing = length(missing),
      matched_genes = paste(matched, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }

  list(scores = as.data.frame(scores), qc = qc)
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

cor_one_vs_all <- function(df, target_cols, method = "spearman") {
  numeric_cols <- colnames(df)[sapply(df, is.numeric)]
  res <- data.frame()

  for (target in target_cols) {
    if (!(target %in% numeric_cols)) next
    for (v in numeric_cols) {
      if (v == target) next
      x <- df[[target]]
      y <- df[[v]]
      ok <- complete.cases(x, y)
      if (sum(ok) >= 5) {
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
        res <- rbind(res, data.frame(
          target = target,
          variable = v,
          n = sum(ok),
          cor = unname(ct$estimate),
          p = ct$p.value,
          method = method,
          stringsAsFactors = FALSE
        ))
      }
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

run_survival <- function(df, score_col, cancer, out_prefix) {
  needed <- c(score_col, "OS_MONTHS", "OS_EVENT")
  if (!all(needed %in% colnames(df))) return(NULL)

  surv_df <- df[complete.cases(df[, needed]), , drop = FALSE]
  if (nrow(surv_df) < 20 || length(unique(surv_df$OS_EVENT)) < 2) return(NULL)

  med <- median(surv_df[[score_col]], na.rm = TRUE)
  surv_df$group <- ifelse(surv_df[[score_col]] >= med, "High", "Low")
  surv_df$group <- factor(surv_df$group, levels = c("Low", "High"))

  fit <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ group, data = surv_df)

  lr <- survdiff(Surv(OS_MONTHS, OS_EVENT) ~ group, data = surv_df)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  p_label <- ifelse(lr_p < 0.001,
                    paste0("log-rank p = ", format(lr_p, scientific = TRUE, digits = 2)),
                    paste0("log-rank p = ", signif(lr_p, 3)))

  score_short <- gsub("_HUMAN_TO_HUMAN", "", score_col)
  score_short <- gsub("_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", "_UP", score_short)
  score_short <- gsub("_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN", "_DOWN", score_short)
  score_short <- gsub("STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", "STEP_UP", score_short)
  score_short <- gsub("STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN", "STEP_DOWN", score_short)

  draw_km <- function() {
    op <- par(mar = c(4.2, 4.2, 3.2, 1.2), cex.axis = 0.82, cex.lab = 0.9, cex.main = 0.85)
    on.exit(par(op), add = TRUE)
    plot(fit, col = c("black", "gray50"), lwd = 2,
         xlab = "Months", ylab = "Overall survival probability",
         main = paste0(cancer, " OS by ", score_short))
    legend("bottomleft", legend = c("Low", "High"), col = c("black", "gray50"), lwd = 2, bty = "n", cex = 0.85)
    usr <- par("usr")
    text(x = usr[1] + 0.05 * (usr[2] - usr[1]),
         y = usr[3] + 0.13 * (usr[4] - usr[3]),
         labels = p_label,
         adj = c(0, 0),
         cex = 0.85)
  }

  pdf(paste0(out_prefix, "_KM.pdf"), width = 6, height = 5)
  draw_km()
  dev.off()

  png(paste0(out_prefix, "_KM.png"), width = 1800, height = 1500, res = 300)
  draw_km()
  dev.off()

  draw_km()

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
    mv_df$CPE_scaled <- as.numeric(scale(mv_df$CPE))
    mv_df$CAF_scaled <- as.numeric(scale(mv_df$CAF_3gene_score))
    mv_df <- mv_df[complete.cases(mv_df[, c("OS_MONTHS", "OS_EVENT", "score_scaled", "CPE_scaled", "CAF_scaled", "AGE")]), ]

    if (nrow(mv_df) >= 20 && length(unique(mv_df$OS_EVENT)) >= 2) {
      cox_mv <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ score_scaled + CPE_scaled + CAF_scaled + AGE, data = mv_df)
      sm <- summary(cox_mv)
      mv_out <- data.frame(
        cancer = cancer,
        score = score_col,
        model = rownames(sm$coefficients),
        n = nrow(mv_df),
        logrank_p = NA,
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
  surv_df$.time <- as.numeric(surv_df[[time_col]])
  surv_df$.event <- normalize_event_vector(surv_df[[event_col]])
  surv_df$.score <- as.numeric(surv_df[[score_col]])
  surv_df <- surv_df[complete.cases(surv_df[, c(".time", ".event", ".score")]), , drop = FALSE]
  surv_df <- surv_df[is.finite(surv_df$.time) & surv_df$.time > 0 & surv_df$.event %in% c(0, 1), , drop = FALSE]

  if (nrow(surv_df) < 20 || length(unique(surv_df$.event)) < 2) return(NULL)

  med <- median(surv_df$.score, na.rm = TRUE)
  surv_df$group <- ifelse(surv_df$.score >= med, "High", "Low")
  surv_df$group <- factor(surv_df$group, levels = c("Low", "High"))

  fit <- survfit(Surv(.time, .event) ~ group, data = surv_df)
  lr <- survdiff(Surv(.time, .event) ~ group, data = surv_df)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  p_label <- ifelse(lr_p < 0.001,
                    paste0("log-rank p = ", format(lr_p, scientific = TRUE, digits = 2)),
                    paste0("log-rank p = ", signif(lr_p, 3)))

  score_short <- score_col
  score_short <- gsub("STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN", "STEP_UP", score_short)
  score_short <- gsub("STEP_DOWN_HYBRID_LEADING_EDGE_170_MOUSE_TO_HUMAN", "STEP_DOWN", score_short)

  draw_km <- function() {
    op <- par(mar = c(4.2, 4.2, 3.2, 1.2), cex.axis = 0.82, cex.lab = 0.9, cex.main = 0.85)
    on.exit(par(op), add = TRUE)
    plot(fit, col = c("black", "gray50"), lwd = 2,
         xlab = "Months", ylab = paste0(endpoint_name, " probability"),
         main = paste0(cancer, " ", endpoint_name, " by ", score_short))
    legend("bottomleft", legend = c("Low", "High"), col = c("black", "gray50"), lwd = 2, bty = "n", cex = 0.85)
    usr <- par("usr")
    text(x = usr[1] + 0.05 * (usr[2] - usr[1]),
         y = usr[3] + 0.13 * (usr[4] - usr[3]),
         labels = p_label, adj = c(0, 0), cex = 0.85)
  }

  pdf(paste0(out_prefix, "_KM.pdf"), width = 6, height = 5)
  draw_km()
  dev.off()
  png(paste0(out_prefix, "_KM.png"), width = 1800, height = 1500, res = 300)
  draw_km()
  dev.off()
  draw_km()

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
    mv_df$CPE_scaled <- as.numeric(scale(mv_df$CPE))
    mv_df <- mv_df[complete.cases(mv_df[, c(".time", ".event", "score_scaled", "CPE_scaled", "AGE")]), , drop = FALSE]
    if (nrow(mv_df) >= 20 && length(unique(mv_df$.event)) >= 2) {
      cox_mv <- coxph(Surv(.time, .event) ~ score_scaled + CPE_scaled + AGE, data = mv_df)
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

  p <- ggplot(dat, aes(x = target_label, y = cor)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.5, width = 0.65) +
    facet_wrap(~ variable_label, nrow = 2) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 9),
          plot.title = element_text(size = 11)) +
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
purity <- read_purity(purity_file)
purity_raw_n <- nrow(purity)
purity <- dedup_purity_by_sample15(purity)
msg("Purity deduplicated by sample15: raw=", purity_raw_n, "; unique=", nrow(purity))

gct_files <- list.files(exp_dir, pattern = "_tcga_expression\\.gct$", full.names = TRUE, ignore.case = TRUE)

msg("Gene sets loaded: ", length(gene_sets))
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

  if (nrow(expr_raw) < 50 || ncol(expr_raw) < 5) {
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
  all_merged[[cancer_upper]] <- merged2

  numeric_df <- merged2[, sapply(merged2, is.numeric), drop = FALSE]
  if (ncol(numeric_df) >= 2 && nrow(numeric_df) >= 5) {

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

    cor_long <- cor_one_vs_all(merged2, targets)
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
      "ESTIMATE_TumorPurity"
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
  all_cols <- unique(unlist(lapply(all_merged, colnames)))
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

msg("\nRun finished: ", as.character(Sys.time()))
msg("Output directory: ", out_dir)
