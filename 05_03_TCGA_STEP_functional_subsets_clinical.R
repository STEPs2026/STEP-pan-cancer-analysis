# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(123)

RESULT_DIR <- NA_character_   
MIN_CLINICAL_N <- 10L
MIN_GROUP_N <- 5L
PNG_DPI <- 300

LOW_COLOR  <- "#006FB9"
HIGH_COLOR <- "#E60011"
GROUP_COLORS <- c("Low" = LOW_COLOR, "High" = HIGH_COLOR)

required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr",
  "scales", "ggrepel", "readxl"
)
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
  library(ggrepel)
})

get_script_dir_v13 <- function() {
  p <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(p) && length(p) == 1 && !is.na(p) && nzchar(p)) {
    return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p2 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p2)) return(dirname(normalizePath(p2, mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

is_functional_result_dir <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x) &&
    file.exists(file.path(
      x, "09_pan_cancer_summary",
      "ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv"
    ))
}

resolve_result_dir_v13 <- function(manual = RESULT_DIR) {
  if (is_functional_result_dir(manual)) {
    return(normalizePath(manual, mustWork = TRUE))
  }

  script_dir <- get_script_dir_v13()
  wd <- normalizePath(getwd(), mustWork = FALSE)

  direct <- unique(c(wd, script_dir, dirname(wd), dirname(script_dir)))
  for (x in direct) {
    if (is_functional_result_dir(x)) return(normalizePath(x, mustWork = TRUE))
  }

  roots <- unique(c(wd, script_dir, dirname(wd), dirname(script_dir)))
  candidates <- character(0)

  for (rr in roots) {
    if (!dir.exists(rr)) next
    subs <- list.dirs(rr, recursive = FALSE, full.names = TRUE)
    keep <- grepl("^TCGA_STEP_functional_subsets_MYC7_Exosome_v1_[23]_", basename(subs))
    subs <- subs[keep]
    if (length(subs) > 0) {
      ok <- vapply(subs, is_functional_result_dir, logical(1))
      candidates <- c(candidates, subs[ok])
    }
  }

  candidates <- unique(candidates)
  if (length(candidates) == 0) {
    stop(
      "Cannot auto-detect the completed v1.2 functional-subset result folder.\n",
      "Please set RESULT_DIR manually at the top of this script.\n",
      "It must contain:\n",
      "09_pan_cancer_summary/ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv"
    )
  }

  mt <- file.info(candidates)$mtime
  normalizePath(candidates[order(mt, decreasing = TRUE)][1], mustWork = TRUE)
}

RESULT_DIR <- resolve_result_dir_v13(RESULT_DIR)
BASE_DIR <- dirname(RESULT_DIR)

TCGA_ROOT_CANDIDATES <- unique(c(
  file.path(BASE_DIR, "TCGA_201906cBioP"),
  file.path(dirname(BASE_DIR), "TCGA_201906cBioP")
))
TCGA_ROOT_HITS <- TCGA_ROOT_CANDIDATES[dir.exists(TCGA_ROOT_CANDIDATES)]
TCGA_ROOT <- if (length(TCGA_ROOT_HITS) > 0) TCGA_ROOT_HITS[1] else NA_character_

INPUT_MERGED <- file.path(
  RESULT_DIR, "09_pan_cancer_summary",
  "ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv"
)
COX_FILE <- file.path(
  RESULT_DIR, "09_pan_cancer_summary",
  "ALL_cancers_all_endpoints_Cox_summary.csv"
)

OUT_CLIN <- file.path(RESULT_DIR, "15_clinical_sex_stage_functional_subsets")
OUT_CLIN_SUM <- file.path(OUT_CLIN, "00_ALL_CANCERS_SUMMARY_v1_3")
OUT_SUMMARY <- file.path(RESULT_DIR, "10_summary_figure_ready")
OUT_SURV <- file.path(RESULT_DIR, "16_survival_volcano_functional_subsets_OS_DFS")
OUT_LOG <- file.path(RESULT_DIR, "00_logs")

for (d in c(OUT_CLIN, OUT_CLIN_SUM, OUT_SUMMARY, OUT_SURV, OUT_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(OUT_LOG, "functional_subset_posthoc_fix_v1_3.log")
msg <- function(...) {
  z <- paste0(...)
  message(z)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      " | ", z, "\n",
      file = LOG_FILE, append = TRUE, sep = "")
}

msg("RESULT_DIR = ", RESULT_DIR)
msg("BASE_DIR = ", BASE_DIR)
msg("TCGA_ROOT = ", ifelse(is.na(TCGA_ROOT), "NOT FOUND", TCGA_ROOT))
msg("INPUT_MERGED = ", INPUT_MERGED)
msg("COX_FILE = ", COX_FILE)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

format_p <- function(p, prefix = "p = ") {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return(paste0(prefix, "NA"))
  if (p < 1e-4) return(paste0(prefix, formatC(p, format = "e", digits = 2)))
  paste0(prefix, signif(p, 3))
}

format_p_short <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 1))
  format(signif(p, 2), trim = TRUE)
}

theme_step <- function(base_size = 10, title_size = 10,
                        axis_title_size = 9, axis_text_size = 8,
                        strip_size = 8) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = max(7, title_size - 1), hjust = 0.5),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size, color = "black"),
      strip.text = element_text(size = strip_size, face = "bold"),
      legend.title = element_text(size = max(7, axis_title_size - 1)),
      legend.text = element_text(size = max(7, axis_text_size - 1)),
      panel.grid.minor = element_blank()
    )
}

save_gg_both <- function(p, file_base, width = 6.5, height = 5.2) {
  print(p)
  ggsave(paste0(file_base, ".pdf"), p, width = width, height = height,
         limitsize = FALSE, useDingbats = FALSE)
  ggsave(paste0(file_base, ".png"), p, width = width, height = height,
         dpi = PNG_DPI, limitsize = FALSE)
}

normalize_tcga_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  x <- gsub("_", "-", x)
  trimws(x)
}
get_patient12 <- function(x) substr(normalize_tcga_id(x), 1, 12)

clean_na <- function(x) {
  z <- trimws(as.character(x))
  z[toupper(z) %in% c(
    "", "NA", "N/A", "NULL", "[NOT AVAILABLE]", "[NOT APPLICABLE]",
    "[UNKNOWN]", "UNKNOWN", "NOT REPORTED", "[NOT EVALUATED]"
  )] <- NA_character_
  z
}
num_clean <- function(x) suppressWarnings(as.numeric(gsub(",", ".", clean_na(x))))

pick_col <- function(dat, exact = character(0), patterns = character(0)) {
  hit <- exact[exact %in% colnames(dat)]
  if (length(hit) > 0) return(hit[1])

  for (pat in patterns) {
    h <- grep(pat, colnames(dat), value = TRUE, ignore.case = TRUE)
    if (length(h) > 0) return(h[1])
  }
  NA_character_
}

SCORE_MAP <- c(
  "Hallmark_MYC_targets_v1" = "MYC targets v1 (7 genes)",
  "GO_extracellular_exsosome" = "Extracellular exosome (23 genes)"
)
SCORE_COLS <- names(SCORE_MAP)

score_label <- function(x) {
  out <- unname(SCORE_MAP[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

if (!file.exists(INPUT_MERGED)) stop("Missing v1.2 merged table: ", INPUT_MERGED)

dat0 <- read.csv(INPUT_MERGED, check.names = FALSE, stringsAsFactors = FALSE)

cancer_col <- pick_col(
  dat0,
  exact = c("cancer", "cancer_loop", "cancer_from_file", "cancer_type", "CANCER"),
  patterns = c("^cancer")
)
if (is.na(cancer_col)) stop("Cannot detect cancer column in: ", INPUT_MERGED)

if (!("patient_id" %in% colnames(dat0))) {
  patient_source <- pick_col(
    dat0,
    exact = c("PATIENT_ID", "Patient.Identifier", "sample", "sample_raw"),
    patterns = c("PATIENT.*ID")
  )
  if (is.na(patient_source)) stop("Cannot detect patient identifier in merged table.")
  dat0$patient_id <- get_patient12(dat0[[patient_source]])
} else {
  dat0$patient_id <- get_patient12(dat0$patient_id)
}

dat0$.cancer <- toupper(trimws(as.character(dat0[[cancer_col]])))

missing_scores <- setdiff(SCORE_COLS, colnames(dat0))
if (length(missing_scores) > 0) {
  stop("Missing score columns in v1.2 merged table: ",
       paste(missing_scores, collapse = ", "))
}

clean_stage_v13 <- function(x) {
  y <- toupper(clean_na(x))
  y <- gsub("PATHOLOGIC|PATHOLOGICAL|CLINICAL", "", y)
  y <- gsub("STAGE", "", y)
  y <- gsub("[^A-Z0-9]", "", y)

  out <- rep(NA_character_, length(y))
  out[grepl("^0($|A|B|C)", y)] <- "0"
  out[grepl("^I($|A|B|C)$|^1($|A|B|C)", y)] <- "I"
  out[grepl("^II($|A|B|C)$|^2($|A|B|C)", y)] <- "II"
  out[grepl("^III($|A|B|C)$|^3($|A|B|C)", y)] <- "III"
  out[grepl("^IV($|A|B|C)$|^4($|A|B|C)", y)] <- "IV"

  factor(out, levels = c("0", "I", "II", "III", "IV"), ordered = TRUE)
}

clean_grade_v13 <- function(x) {
  y <- toupper(clean_na(x))
  y <- gsub("HISTOLOGIC|HISTOLOGICAL|TUMOR|NEOPLASM|GRADE", "", y)
  y <- gsub("[^A-Z0-9]", "", y)

  out <- rep(NA_character_, length(y))
  out[grepl("(^|G)1$|LOW|WELLDIFFERENTIATED", y)] <- "G1"
  out[grepl("(^|G)2$|MODERATELYDIFFERENTIATED", y)] <- "G2"
  out[grepl("(^|G)3$|POORLYDIFFERENTIATED", y)] <- "G3"
  out[grepl("(^|G)4$|UNDIFFERENTIATED|HIGH", y)] <- "G4"

  out[is.na(out) & grepl("^G1", y)] <- "G1"
  out[is.na(out) & grepl("^G2", y)] <- "G2"
  out[is.na(out) & grepl("^G3", y)] <- "G3"
  out[is.na(out) & grepl("^G4", y)] <- "G4"

  factor(out, levels = c("G1", "G2", "G3", "G4"), ordered = TRUE)
}

clean_sex_v13 <- function(x) {
  y <- toupper(clean_na(x))
  out <- rep(NA_character_, length(y))
  out[grepl("^M$|MALE", y)] <- "Male"
  out[grepl("^F$|FEMALE", y)] <- "Female"
  factor(out, levels = c("Male", "Female"))
}

standardize_clinical_any <- function(clin, source_name = NA_character_) {
  clin <- as.data.frame(clin, stringsAsFactors = FALSE, check.names = FALSE)
  original_names <- colnames(clin)
  colnames(clin) <- make.names(colnames(clin), unique = TRUE)

  patient_col <- pick_col(
    clin,
    exact = c("PATIENT_ID", "Patient.Identifier", "patient_id", "submitter_id",
              "bcr_patient_barcode", "BCR_PATIENT_BARCODE"),
    patterns = c("PATIENT.*ID", "PATIENT.*BARCODE", "SUBMITTER.*ID")
  )
  if (is.na(patient_col)) return(NULL)

  age_col <- pick_col(
    clin,
    exact = c("AGE", "Age", "Diagnosis.Age", "AGE_AT_DIAGNOSIS"),
    patterns = c("^AGE$", "AGE.*DIAG", "DIAG.*AGE")
  )
  sex_col <- pick_col(
    clin,
    exact = c("SEX", "Sex", "Gender", "GENDER"),
    patterns = c("^SEX$", "^GENDER$")
  )
  stage_col <- pick_col(
    clin,
    exact = c(
      "AJCC_PATHOLOGIC_TUMOR_STAGE",
      "AJCC_PATHOLOGIC_STAGE",
      "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code",
      "CLINICAL_STAGE", "PATH_STAGE", "PATHOLOGIC_STAGE", "STAGE",
      "tumor_stage"
    ),
    patterns = c(
      "AJCC.*PATH.*STAGE", "PATHOLOG.*STAGE", "CLINICAL.*STAGE",
      "^STAGE$", "TUMOR.*STAGE"
    )
  )
  grade_col <- pick_col(
    clin,
    exact = c(
      "GRADE", "HISTOLOGIC_GRADE", "NEOPLASM_HISTOLOGIC_GRADE",
      "TUMOR_GRADE", "GRADE_GROUP"
    ),
    patterns = c(
      "HISTOLOG.*GRADE", "NEOPLASM.*GRADE", "TUMOR.*GRADE",
      "^GRADE$", "GRADE.*GROUP"
    )
  )

  out <- data.frame(
    patient_id = get_patient12(clin[[patient_col]]),
    source_name = source_name,
    stringsAsFactors = FALSE
  )

  out$AGE_NEW <- if (!is.na(age_col)) num_clean(clin[[age_col]]) else NA_real_
  
  if (sum(is.finite(out$AGE_NEW)) > 0 &&
      median(out$AGE_NEW, na.rm = TRUE) > 200) {
    out$AGE_NEW <- out$AGE_NEW / 365.25
  }

  out$SEX_NEW <- if (!is.na(sex_col)) as.character(clean_sex_v13(clin[[sex_col]])) else NA_character_
  out$STAGE_RAW_NEW <- if (!is.na(stage_col)) clean_na(clin[[stage_col]]) else NA_character_
  out$STAGE_NEW <- as.character(clean_stage_v13(out$STAGE_RAW_NEW))
  out$GRADE_RAW_NEW <- if (!is.na(grade_col)) clean_na(clin[[grade_col]]) else NA_character_
  out$GRADE_NEW <- as.character(clean_grade_v13(out$GRADE_RAW_NEW))

  attr(out, "detected_columns") <- data.frame(
    source = source_name,
    patient_col = patient_col,
    age_col = age_col,
    sex_col = sex_col,
    stage_col = stage_col,
    grade_col = grade_col,
    original_column_count = length(original_names),
    stringsAsFactors = FALSE
  )
  out
}

clinical_sources <- list()
clinical_source_qc <- list()
cs_i <- 1L

if (!is.na(TCGA_ROOT) && dir.exists(TCGA_ROOT)) {
  all_txt <- list.files(
    TCGA_ROOT,
    pattern = "\\.txt$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  all_txt <- all_txt[
    grepl("clinical", basename(all_txt), ignore.case = TRUE) &
      grepl("patient", basename(all_txt), ignore.case = TRUE)
  ]

  cancers <- sort(unique(dat0$.cancer))
  cancers <- cancers[!is.na(cancers) & nzchar(cancers)]

  for (ca in cancers) {
    ca_lower <- tolower(ca)

    cand <- all_txt[
      grepl(paste0("^", ca_lower, "(_|-).*clinical"), basename(all_txt), ignore.case = TRUE) |
        grepl(paste0("^", ca_lower, "(_|-).*patient"), basename(all_txt), ignore.case = TRUE)
    ]
    if (length(cand) == 0) next

    f <- cand[1]
    raw <- tryCatch(
      read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""),
      error = function(e) NULL
    )

    if (is.null(raw) || !any(grepl("PATIENT", colnames(raw), ignore.case = TRUE))) {
      lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
      if (length(lines) > 0) {
        spl <- strsplit(lines, "\t", fixed = TRUE)
        header_idx <- which(vapply(
          spl,
          function(z) any(grepl("PATIENT_ID|PATIENT.*ID", z, ignore.case = TRUE)),
          logical(1)
        ))
        if (length(header_idx) > 0) {
          raw <- tryCatch(
            read.delim(
              f, skip = header_idx[1] - 1,
              check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""
            ),
            error = function(e) NULL
          )
        }
      }
    }

    if (!is.null(raw) && nrow(raw) > 0) {
      std <- standardize_clinical_any(raw, source_name = paste0("TXT:", basename(f)))
      if (!is.null(std) && nrow(std) > 0) {
        std$.source_cancer <- ca
        clinical_sources[[length(clinical_sources) + 1L]] <- std

        dc <- attr(std, "detected_columns")
        if (!is.null(dc)) {
          dc$cancer <- ca
          dc$file <- f
          clinical_source_qc[[cs_i]] <- dc
          cs_i <- cs_i + 1L
        }
      }
    }
  }

  xlsx_files <- list.files(
    TCGA_ROOT,
    pattern = "TCGA.*Clinical.*\\.xlsx$|Clinical.*TCGA.*\\.xlsx$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(xlsx_files) > 0 && requireNamespace("readxl", quietly = TRUE)) {
    xlsx_file <- xlsx_files[1]
    msg("Optional master clinical XLSX detected: ", xlsx_file)

    sheets <- tryCatch(readxl::excel_sheets(xlsx_file), error = function(e) character(0))
    for (sh in sheets) {
      xr <- tryCatch(
        readxl::read_excel(xlsx_file, sheet = sh, guess_max = 10000),
        error = function(e) NULL
      )
      if (is.null(xr) || nrow(xr) == 0) next

      std <- standardize_clinical_any(
        as.data.frame(xr, stringsAsFactors = FALSE),
        source_name = paste0("XLSX:", basename(xlsx_file), "::", sh)
      )
      if (is.null(std) || nrow(std) == 0) next

      std$.source_cancer <- NA_character_
      clinical_sources[[length(clinical_sources) + 1L]] <- std

      dc <- attr(std, "detected_columns")
      if (!is.null(dc)) {
        dc$cancer <- NA_character_
        dc$file <- paste0(xlsx_file, " [sheet=", sh, "]")
        clinical_source_qc[[cs_i]] <- dc
        cs_i <- cs_i + 1L
      }
    }
  }
}

if (length(clinical_source_qc) > 0) {
  write.csv(
    bind_rows(clinical_source_qc),
    file.path(OUT_CLIN, "clinical_source_detected_columns_v1_3.csv"),
    row.names = FALSE
  )
}

clin_new <- if (length(clinical_sources) > 0) bind_rows(clinical_sources) else data.frame()

if (nrow(clin_new) > 0) {
  useful_new <- c("AGE_NEW", "SEX_NEW", "STAGE_NEW", "GRADE_NEW")
  clin_new$.richness <- rowSums(!is.na(clin_new[, useful_new, drop = FALSE]))
  
  clin_new$.txt_priority <- ifelse(grepl("^TXT:", clin_new$source_name), 1L, 0L)
  clin_new <- clin_new %>%
    arrange(patient_id, desc(.richness), desc(.txt_priority)) %>%
    distinct(patient_id, .keep_all = TRUE)
  clin_new$.richness <- NULL
  clin_new$.txt_priority <- NULL

  dat <- merge(dat0, clin_new, by = "patient_id", all.x = TRUE, sort = FALSE)
} else {
  msg("No reparsable external clinical source was found; using existing v1.2 clinical columns only.")
  dat <- dat0
}

first_nonmissing_chr <- function(old, new) {
  old <- clean_na(old)
  new <- clean_na(new)
  use_new <- !is.na(new)
  out <- old
  out[use_new] <- new[use_new]
  out
}

first_nonmissing_num <- function(old, new) {
  old <- suppressWarnings(as.numeric(old))
  new <- suppressWarnings(as.numeric(new))
  out <- old
  use_new <- is.finite(new)
  out[use_new] <- new[use_new]
  out
}

if (!("AGE" %in% colnames(dat))) dat$AGE <- NA_real_
if (!("SEX" %in% colnames(dat))) dat$SEX <- NA_character_
if (!("STAGE_RAW" %in% colnames(dat))) dat$STAGE_RAW <- NA_character_
if (!("GRADE_RAW" %in% colnames(dat))) dat$GRADE_RAW <- NA_character_

if ("AGE_NEW" %in% colnames(dat)) dat$AGE <- first_nonmissing_num(dat$AGE, dat$AGE_NEW)
if ("SEX_NEW" %in% colnames(dat)) dat$SEX <- first_nonmissing_chr(dat$SEX, dat$SEX_NEW)
if ("STAGE_RAW_NEW" %in% colnames(dat)) dat$STAGE_RAW <- first_nonmissing_chr(dat$STAGE_RAW, dat$STAGE_RAW_NEW)
if ("GRADE_RAW_NEW" %in% colnames(dat)) dat$GRADE_RAW <- first_nonmissing_chr(dat$GRADE_RAW, dat$GRADE_RAW_NEW)

old_stage <- if ("STAGE" %in% colnames(dat)) as.character(dat$STAGE) else rep(NA_character_, nrow(dat))
old_grade <- if ("GRADE" %in% colnames(dat)) as.character(dat$GRADE) else rep(NA_character_, nrow(dat))

stage_reclean <- as.character(clean_stage_v13(dat$STAGE_RAW))
grade_reclean <- as.character(clean_grade_v13(dat$GRADE_RAW))

stage_final <- old_stage
use_stage_new <- !is.na(stage_reclean)
stage_final[use_stage_new] <- stage_reclean[use_stage_new]

grade_final <- old_grade
use_grade_new <- !is.na(grade_reclean)
grade_final[use_grade_new] <- grade_reclean[use_grade_new]

dat$SEX <- as.character(clean_sex_v13(dat$SEX))
dat$SEX <- factor(dat$SEX, levels = c("Male", "Female"))

dat$STAGE <- factor(stage_final, levels = c("0", "I", "II", "III", "IV"), ordered = TRUE)
dat$STAGE_NUM <- as.numeric(dat$STAGE) - 1
dat$STAGE_GROUP <- ifelse(
  dat$STAGE %in% c("0", "I", "II"), "Early",
  ifelse(dat$STAGE %in% c("III", "IV"), "Advanced", NA_character_)
)
dat$STAGE_GROUP <- factor(dat$STAGE_GROUP, levels = c("Early", "Advanced"))

dat$GRADE <- factor(grade_final, levels = c("G1", "G2", "G3", "G4"), ordered = TRUE)

if (!(".cancer" %in% colnames(dat))) {
  dat$.cancer <- toupper(trimws(as.character(dat[[cancer_col]])))
}

REPAIRED_MERGED <- file.path(
  RESULT_DIR, "09_pan_cancer_summary",
  "ALL_cancers_MYC7_Exosome_scores_clinical_merged_v1_3_repaired.csv"
)
write.csv(dat, REPAIRED_MERGED, row.names = FALSE)

clinical_colors <- list(
  SEX = c("Male" = LOW_COLOR, "Female" = HIGH_COLOR),
  STAGE = c(
    "0" = "#D9D9D9", "I" = "#A6CEE3", "II" = "#6BAED6",
    "III" = "#FB6A4A", "IV" = "#CB181D"
  ),
  STAGE_GROUP = c("Early" = LOW_COLOR, "Advanced" = HIGH_COLOR),
  GRADE = c(
    "G1" = "#A6CEE3", "G2" = "#6BAED6",
    "G3" = "#FB6A4A", "G4" = "#CB181D"
  )
)

clinical_label <- c(
  SEX = "Sex",
  STAGE = "Stage",
  STAGE_GROUP = "Stage group",
  GRADE = "Grade"
)

make_score_group <- function(x) {
  med <- median(x, na.rm = TRUE)
  factor(ifelse(x >= med, "High", "Low"), levels = c("Low", "High"))
}

safe_comp_test <- function(tab) {
  if (length(tab) == 0 || any(dim(tab) < 2) || sum(tab) == 0) {
    return(list(test = "NA", p = NA_real_))
  }

  cs <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
  use_fisher <- is.null(cs) || any(cs$expected < 5)

  if (use_fisher) {
    ft <- tryCatch(
      fisher.test(
        tab,
        simulate.p.value = any(dim(tab) > 2),
        B = 10000
      ),
      error = function(e) NULL
    )
    return(list(test = "Fisher", p = if (is.null(ft)) NA_real_ else ft$p.value))
  }

  list(test = "Chi-square", p = cs$p.value)
}

clinical_stats <- list()
clinical_qc <- list()
plot_data_store <- list()
stat_i <- 1L
qc_i <- 1L

cancers <- sort(unique(dat$.cancer))
cancers <- cancers[!is.na(cancers) & nzchar(cancers)]

for (ca in cancers) {
  dc <- dat[dat$.cancer == ca, , drop = FALSE]
  ca_dir <- file.path(OUT_CLIN, ca)
  dir.create(ca_dir, recursive = TRUE, showWarnings = FALSE)

  qc_ca <- list()
  qc_ca_i <- 1L

  for (sc in SCORE_COLS) {
    for (var in c("SEX", "STAGE", "STAGE_GROUP", "GRADE")) {
      score_vec <- suppressWarnings(as.numeric(dc[[sc]]))
      var_vec <- if (var %in% colnames(dc)) dc[[var]] else rep(NA, nrow(dc))

      ok <- is.finite(score_vec) & !is.na(var_vec)
      dd <- dc[ok, , drop = FALSE]
      dd$.score <- score_vec[ok]
      dd$.clinical <- droplevels(factor(as.character(var_vec[ok])))

      n_complete <- nrow(dd)
      n_levels <- if (n_complete > 0) length(unique(dd$.clinical)) else 0L
      level_counts <- if (n_complete > 0) paste(names(table(dd$.clinical)), as.integer(table(dd$.clinical)), sep = ":", collapse = ";") else ""

      eligible <- n_complete >= MIN_CLINICAL_N && n_levels >= 2
      reason <- if (eligible) {
        "eligible"
      } else if (n_complete < MIN_CLINICAL_N) {
        paste0("complete_n_lt_", MIN_CLINICAL_N)
      } else {
        "clinical_variable_has_lt_2_levels"
      }

      qc_row <- data.frame(
        cancer = ca,
        score = sc,
        score_label = score_label(sc),
        variable = var,
        n_cancer_total = nrow(dc),
        n_score_nonmissing = sum(is.finite(score_vec)),
        n_clinical_nonmissing = sum(!is.na(var_vec)),
        n_complete = n_complete,
        n_levels = n_levels,
        level_counts = level_counts,
        eligible = eligible,
        exclusion_reason = reason,
        stringsAsFactors = FALSE
      )
      clinical_qc[[qc_i]] <- qc_row
      qc_ca[[qc_ca_i]] <- qc_row
      qc_i <- qc_i + 1L
      qc_ca_i <- qc_ca_i + 1L

      if (!eligible) next

      dd$.group <- make_score_group(dd$.score)

      if (n_levels == 2) {
        p_score <- tryCatch(
          wilcox.test(.score ~ .clinical, data = dd, exact = FALSE)$p.value,
          error = function(e) NA_real_
        )
        score_test <- "Wilcoxon"
      } else {
        p_score <- tryCatch(
          kruskal.test(.score ~ .clinical, data = dd)$p.value,
          error = function(e) NA_real_
        )
        score_test <- "Kruskal-Wallis"
      }

      ct <- safe_comp_test(table(dd$.group, dd$.clinical))
      p_comp <- ct$p

      clinical_stats[[stat_i]] <- data.frame(
        cancer = ca,
        score = sc,
        score_label = score_label(sc),
        variable = var,
        n = n_complete,
        n_levels = n_levels,
        score_test = score_test,
        score_p = p_score,
        composition_test = ct$test,
        composition_p = p_comp,
        stringsAsFactors = FALSE
      )
      stat_i <- stat_i + 1L

      store_key <- paste(sc, var, ca, sep = "__")
      plot_data_store[[store_key]] <- dd

      pal <- clinical_colors[[var]]

      p1 <- ggplot(dd, aes(x = .clinical, y = .score, fill = .clinical, color = .clinical)) +
        geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82, linewidth = 0.35) +
        geom_jitter(width = 0.14, size = 0.75, alpha = 0.38, show.legend = FALSE) +
        scale_fill_manual(values = pal, drop = FALSE) +
        scale_color_manual(values = pal, drop = FALSE) +
        labs(
          title = paste0(ca, ": ", score_label(sc), " by ", clinical_label[[var]]),
          subtitle = paste0(score_test, " ", format_p(p_score)),
          x = clinical_label[[var]],
          y = "Z-mean score"
        ) +
        theme_step(10) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1)
        )

      save_gg_both(
        p1,
        file.path(ca_dir, paste0(ca, "_", safe_name(sc), "_score_by_", var, "_v1_3")),
        5.8, 5.0
      )

      comp <- as.data.frame(table(dd$.group, dd$.clinical))
      colnames(comp) <- c("group", "clinical", "n")
      comp <- comp %>%
        group_by(group) %>%
        mutate(percent = 100 * n / sum(n)) %>%
        ungroup()

      p2 <- ggplot(comp, aes(x = group, y = percent, fill = clinical)) +
        geom_col(width = 0.70, color = "grey35", linewidth = 0.2) +
        geom_text(
          aes(label = ifelse(percent >= 5, paste0(round(percent, 1), "%"), "")),
          position = position_stack(vjust = 0.5),
          size = 3,
          color = "black",
          check_overlap = TRUE
        ) +
        scale_fill_manual(values = pal, drop = FALSE) +
        scale_y_continuous(
          limits = c(0, 114),
          expand = expansion(mult = c(0, 0))
        ) +
        annotate(
          "text", x = 1.5, y = 108,
          label = paste0(ct$test, " ", format_p(p_comp)),
          size = 3.1
        ) +
        labs(
          title = paste0(ca, ": ", clinical_label[[var]],
                         " composition by ", score_label(sc), " High/Low"),
          x = NULL, y = "Composition (%)", fill = clinical_label[[var]]
        ) +
        theme_step(10) +
        theme(legend.position = "top")

      save_gg_both(
        p2,
        file.path(ca_dir, paste0(ca, "_", safe_name(sc), "_", var,
                                 "_composition_by_high_low_v1_3")),
        5.8, 5.0
      )
    }
  }

  qc_ca_df <- bind_rows(qc_ca)
  write.csv(
    qc_ca_df,
    file.path(ca_dir, paste0(ca, "_clinical_eligibility_QC_v1_3.csv")),
    row.names = FALSE
  )

  if (!any(qc_ca_df$eligible)) {
    writeLines(
      c(
        paste0("No eligible SEX/STAGE/STAGE_GROUP/GRADE plot for ", ca, "."),
        paste0("Minimum complete N = ", MIN_CLINICAL_N, "; at least 2 clinical levels required."),
        "See the CSV in this folder for exact counts and exclusion reasons.",
        "This may reflect unavailable/not-applicable clinical metadata, insufficient matched samples,",
        "or failure of a source column to exist in the available cBioPortal/Excel clinical files."
      ),
      file.path(ca_dir, "NO_ELIGIBLE_CLINICAL_PLOTS_README_v1_3.txt")
    )
  }
}

clin_qc_df <- bind_rows(clinical_qc)
clin_stats_df <- bind_rows(clinical_stats)

write.csv(
  clin_qc_df,
  file.path(OUT_CLIN, "ALL_cancers_clinical_eligibility_QC_v1_3.csv"),
  row.names = FALSE
)

if (nrow(clin_stats_df) > 0) {
  clin_stats_df <- clin_stats_df %>%
    group_by(score, variable) %>%
    mutate(
      score_FDR = p.adjust(score_p, method = "BH"),
      composition_FDR = p.adjust(composition_p, method = "BH")
    ) %>%
    ungroup()

  write.csv(
    clin_stats_df,
    file.path(OUT_CLIN, "ALL_cancers_clinical_statistics_v1_3.csv"),
    row.names = FALSE
  )
  write.csv(
    clin_stats_df,
    file.path(OUT_SUMMARY, "FigureReady_clinical_all_cancers_with_FDR_v1_3.csv"),
    row.names = FALSE
  )
}

make_faceted_score_box <- function(sc, var) {
  keys <- grep(
    paste0("^", sc, "__", var, "__"),
    names(plot_data_store),
    value = TRUE
  )
  if (length(keys) == 0) return(invisible(NULL))

  dd <- bind_rows(lapply(keys, function(k) {
    z <- plot_data_store[[k]]
    z$.cancer_plot <- sub(paste0("^", sc, "__", var, "__"), "", k)
    z
  }))

  st <- clin_stats_df %>%
    filter(score == sc, variable == var) %>%
    transmute(
      .cancer_plot = cancer,
      label = vapply(score_p, format_p, character(1), prefix = "p = ")
    )

  ypos <- dd %>%
    group_by(.cancer_plot) %>%
    summarise(
      ymax = max(.score, na.rm = TRUE),
      ymin = min(.score, na.rm = TRUE),
      nlev = n_distinct(.clinical),
      .groups = "drop"
    ) %>%
    mutate(
      y = ymax + 0.13 * ifelse(ymax > ymin, ymax - ymin, 1),
      x = (nlev + 1) / 2
    )

  ann <- left_join(st, ypos, by = ".cancer_plot")
  pal <- clinical_colors[[var]]

  p <- ggplot(dd, aes(x = .clinical, y = .score, fill = .clinical, color = .clinical)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82, linewidth = 0.25) +
    geom_jitter(width = 0.11, size = 0.32, alpha = 0.25, show.legend = FALSE) +
    geom_text(
      data = ann,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.15,
      color = "black"
    ) +
    facet_wrap(~ .cancer_plot, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = pal, drop = FALSE) +
    scale_color_manual(values = pal, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      title = paste0("All cancers: ", score_label(sc), " by ", clinical_label[[var]]),
      x = NULL,
      y = "Z-mean score"
    ) +
    theme_step(base_size = 8.5, title_size = 10,
                axis_title_size = 8, axis_text_size = 6.8, strip_size = 7.2) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  nca <- n_distinct(dd$.cancer_plot)
  base <- file.path(
    OUT_CLIN_SUM,
    paste0(safe_name(sc), "_ALL_cancers_score_by_", var,
           "_boxplot_faceted_with_p")
  )
  save_gg_both(p, base, width = 13, height = max(6, ceiling(nca / 5) * 2.7))
}

make_faceted_composition <- function(sc, var) {
  keys <- grep(
    paste0("^", sc, "__", var, "__"),
    names(plot_data_store),
    value = TRUE
  )
  if (length(keys) == 0) return(invisible(NULL))

  dd <- bind_rows(lapply(keys, function(k) {
    z <- plot_data_store[[k]]
    z$.cancer_plot <- sub(paste0("^", sc, "__", var, "__"), "", k)
    z
  }))

  comp <- dd %>%
    count(.cancer_plot, .group, .clinical, name = "n") %>%
    group_by(.cancer_plot, .group) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  st <- clin_stats_df %>%
    filter(score == sc, variable == var) %>%
    transmute(
      .cancer_plot = cancer,
      x = 1.5,
      y = 1.08,
      label = vapply(composition_p, format_p, character(1), prefix = "p = ")
    )

  pal <- clinical_colors[[var]]

  p <- ggplot(comp, aes(x = .group, y = prop, fill = .clinical)) +
    geom_col(width = 0.72, color = "grey35", linewidth = 0.15) +
    geom_text(
      aes(label = ifelse(prop >= 0.06, scales::percent(prop, accuracy = 1), "")),
      position = position_stack(vjust = 0.5),
      size = 1.85,
      color = "black",
      check_overlap = TRUE
    ) +
    geom_text(
      data = st,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.2
    ) +
    facet_wrap(~ .cancer_plot, ncol = 5) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.14),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(
      title = paste0(
        "All cancers: ", clinical_label[[var]],
        " composition by ", score_label(sc), " group"
      ),
      x = paste0(score_label(sc), " group"),
      y = "Composition",
      fill = clinical_label[[var]]
    ) +
    theme_step(base_size = 8.5, title_size = 10,
                axis_title_size = 8, axis_text_size = 6.8, strip_size = 7.2) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  nca <- n_distinct(comp$.cancer_plot)
  base <- file.path(
    OUT_CLIN_SUM,
    paste0(safe_name(sc), "_ALL_cancers_", var,
           "_composition_faceted_with_percent_and_p")
  )
  save_gg_both(p, base, width = 13, height = max(6, ceiling(nca / 5) * 2.6))
}

if (nrow(clin_stats_df) > 0) {
  for (sc in SCORE_COLS) {
    for (var in c("SEX", "STAGE", "STAGE_GROUP", "GRADE")) {
      make_faceted_score_box(sc, var)
      make_faceted_composition(sc, var)
    }

    hs <- clin_stats_df %>%
      filter(score == sc) %>%
      mutate(
        analysis = paste0("score_by_", variable),
        neglog10FDR = -log10(pmax(score_FDR, 1e-300)),
        p_label = vapply(score_p, format_p_short, character(1))
      )

    if (nrow(hs) > 0) {
      p_hm <- ggplot(
        hs,
        aes(
          x = factor(cancer, levels = sort(unique(cancer))),
          y = analysis,
          fill = neglog10FDR
        )
      ) +
        geom_tile(color = "white", linewidth = 0.25) +
        geom_text(aes(label = paste0("p=", p_label)), size = 2.0, color = "black") +
        scale_fill_gradient(
          low = "#FFF7EC", high = "#D7301F",
          name = expression(-log[10](FDR))
        ) +
        labs(
          title = paste0(score_label(sc), " clinical associations across cancers"),
          subtitle = "Numbers shown in tiles are raw p values",
          x = "TCGA cancer type",
          y = NULL
        ) +
        theme_step(base_size = 9, title_size = 10,
                    axis_title_size = 8, axis_text_size = 7.2, strip_size = 7) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

      base_hm <- file.path(
        OUT_CLIN_SUM,
        paste0(safe_name(sc), "_ALL_cancers_clinical_association_pvalue_labeled_heatmap")
      )
      save_gg_both(
        p_hm, base_hm,
        width = max(12, length(unique(hs$cancer)) * 0.45),
        height = max(5.5, length(unique(hs$analysis)) * 0.75 + 2)
      )
    }
  }
}

empty_focus <- c(
  "UCS", "UCEC", "THYM", "SARC", "PRAD", "PCPG",
  "LUSC", "LAML", "GBM", "DLBC"
)
focus_qc <- clin_qc_df %>% filter(cancer %in% empty_focus)
write.csv(
  focus_qc,
  file.path(OUT_CLIN, "Previously_empty_cancers_clinical_QC_v1_3.csv"),
  row.names = FALSE
)

RECOMPUTE_COX_FROM_REPAIRED_CLINICAL <- TRUE
REPAIRED_COX_FILE <- file.path(
  RESULT_DIR, "09_pan_cancer_summary",
  "ALL_cancers_all_endpoints_Cox_summary_v1_3_repaired_clinical.csv"
)

normalize_event_v13 <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    y <- suppressWarnings(as.numeric(x))
    y[!(y %in% c(0, 1))] <- NA_real_
    return(y)
  }

  z <- toupper(clean_na(x))
  y <- rep(NA_real_, length(z))
  y[grepl("DECEASED|DEAD|PROGRESSED|RECURRENCE|RECURRED|RELAPSED|^1", z)] <- 1
  y[grepl("LIVING|ALIVE|DISEASE.?FREE|CENSORED|^0", z)] <- 0
  y
}

recompute_one_cox_v13 <- function(d, cancer, score_col, endpoint,
                                  time_col, event_col) {
  if (!all(c(score_col, time_col, event_col) %in% colnames(d))) return(NULL)

  dd <- data.frame(
    time = suppressWarnings(as.numeric(d[[time_col]])),
    event = normalize_event_v13(d[[event_col]]),
    score = suppressWarnings(as.numeric(d[[score_col]])),
    AGE = if ("AGE" %in% colnames(d)) suppressWarnings(as.numeric(d$AGE)) else NA_real_,
    SEX = if ("SEX" %in% colnames(d)) as.character(d$SEX) else NA_character_,
    STAGE_NUM = if ("STAGE_NUM" %in% colnames(d)) suppressWarnings(as.numeric(d$STAGE_NUM)) else NA_real_,
    CPE = if ("CPE" %in% colnames(d)) suppressWarnings(as.numeric(d$CPE)) else NA_real_,
    stringsAsFactors = FALSE
  )

  dd <- dd[
    is.finite(dd$time) & dd$time > 0 &
      dd$event %in% c(0, 1) &
      is.finite(dd$score),
    , drop = FALSE
  ]

  if (nrow(dd) < 20 ||
      length(unique(dd$event)) < 2 ||
      length(unique(dd$score)) < 2) {
    return(NULL)
  }

  cutoff <- median(dd$score, na.rm = TRUE)
  dd$score_scaled <- as.numeric(scale(dd$score))

  rows <- list()
  rr <- 1L

  fit_uni <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ score_scaled, data = dd),
    error = function(e) NULL
  )
  if (!is.null(fit_uni)) {
    su <- summary(fit_uni)
    rows[[rr]] <- data.frame(
      cancer = cancer,
      endpoint = endpoint,
      score = score_col,
      model = "univariate_score",
      covariates = "score_only",
      n = nrow(dd),
      n_events = sum(dd$event == 1, na.rm = TRUE),
      cutoff = cutoff,
      logrank_p = NA_real_,
      HR = su$coefficients["score_scaled", "exp(coef)"],
      p = su$coefficients["score_scaled", "Pr(>|z|)"],
      lower95 = su$conf.int["score_scaled", "lower .95"],
      upper95 = su$conf.int["score_scaled", "upper .95"],
      stringsAsFactors = FALSE
    )
    rr <- rr + 1L
  }

  mv <- dd[
    complete.cases(dd[, c("time", "event", "score_scaled", "AGE", "SEX", "STAGE_NUM")]),
    , drop = FALSE
  ]
  mv$SEX <- factor(mv$SEX, levels = c("Male", "Female"))

  if (nrow(mv) >= 20 &&
      length(unique(mv$event)) >= 2 &&
      length(unique(na.omit(mv$SEX))) >= 2 &&
      length(unique(mv$STAGE_NUM)) >= 2) {

    fit_mv <- tryCatch(
      survival::coxph(
        survival::Surv(time, event) ~ score_scaled + AGE + SEX + STAGE_NUM,
        data = mv
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_mv)) {
      sm <- summary(fit_mv)
      if ("score_scaled" %in% rownames(sm$coefficients)) {
        rows[[rr]] <- data.frame(
          cancer = cancer,
          endpoint = endpoint,
          score = score_col,
          model = "score_scaled",
          covariates = "score + AGE + SEX + STAGE",
          n = nrow(mv),
          n_events = sum(mv$event == 1, na.rm = TRUE),
          cutoff = cutoff,
          logrank_p = NA_real_,
          HR = sm$coefficients["score_scaled", "exp(coef)"],
          p = sm$coefficients["score_scaled", "Pr(>|z|)"],
          lower95 = sm$conf.int["score_scaled", "lower .95"],
          upper95 = sm$conf.int["score_scaled", "upper .95"],
          stringsAsFactors = FALSE
        )
        rr <- rr + 1L
      }
    }
  }

  mv2 <- dd[
    complete.cases(dd[, c("time", "event", "score_scaled", "CPE", "AGE")]),
    , drop = FALSE
  ]

  if (nrow(mv2) >= 20 &&
      length(unique(mv2$event)) >= 2 &&
      length(unique(mv2$CPE)) >= 2) {

    mv2$CPE_scaled <- as.numeric(scale(mv2$CPE))
    fit_mv2 <- tryCatch(
      survival::coxph(
        survival::Surv(time, event) ~ score_scaled + CPE_scaled + AGE,
        data = mv2
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_mv2)) {
      sm <- summary(fit_mv2)
      if ("score_scaled" %in% rownames(sm$coefficients)) {
        rows[[rr]] <- data.frame(
          cancer = cancer,
          endpoint = endpoint,
          score = score_col,
          model = "score_scaled_CPE_AGE",
          covariates = "score + CPE + AGE",
          n = nrow(mv2),
          n_events = sum(mv2$event == 1, na.rm = TRUE),
          cutoff = cutoff,
          logrank_p = NA_real_,
          HR = sm$coefficients["score_scaled", "exp(coef)"],
          p = sm$coefficients["score_scaled", "Pr(>|z|)"],
          lower95 = sm$conf.int["score_scaled", "lower .95"],
          upper95 = sm$conf.int["score_scaled", "upper .95"],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  bind_rows(rows)
}

if (isTRUE(RECOMPUTE_COX_FROM_REPAIRED_CLINICAL)) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    install.packages("survival", dependencies = TRUE)
  }

  endpoint_defs_v13 <- list(
    OS  = c(time = "OS_MONTHS",  event = "OS_EVENT"),
    DFS = c(time = "DFS_MONTHS", event = "DFS_EVENT"),
    PFS = c(time = "PFS_MONTHS", event = "PFS_EVENT"),
    RFS = c(time = "RFS_MONTHS", event = "RFS_EVENT")
  )

  cox_repaired_list <- list()
  cx_i <- 1L

  for (ca in cancers) {
    dc <- dat[dat$.cancer == ca, , drop = FALSE]

    for (sc in SCORE_COLS) {
      for (ep in names(endpoint_defs_v13)) {
        tc <- endpoint_defs_v13[[ep]]["time"]
        ec <- endpoint_defs_v13[[ep]]["event"]

        if (!all(c(tc, ec) %in% colnames(dc))) next

        rr <- tryCatch(
          recompute_one_cox_v13(dc, ca, sc, ep, tc, ec),
          error = function(e) {
            msg("Repaired Cox failed: ", ca, " / ", sc, " / ", ep,
                " | ", conditionMessage(e))
            NULL
          }
        )

        if (!is.null(rr) && nrow(rr) > 0) {
          cox_repaired_list[[cx_i]] <- rr
          cx_i <- cx_i + 1L
        }
      }
    }
  }

  cox_repaired <- bind_rows(cox_repaired_list)

  if (nrow(cox_repaired) > 0) {
    cox_repaired <- cox_repaired %>%
      group_by(score, endpoint, model) %>%
      mutate(FDR = p.adjust(p, method = "BH")) %>%
      ungroup()

    write.csv(cox_repaired, REPAIRED_COX_FILE, row.names = FALSE)
    msg("Recomputed Cox table from repaired clinical data: ", REPAIRED_COX_FILE)
  } else {
    warning("No Cox result could be recomputed from repaired clinical data.")
  }
}

COX_FOR_PLOT <- if (file.exists(REPAIRED_COX_FILE)) REPAIRED_COX_FILE else COX_FILE

if (!file.exists(COX_FOR_PLOT)) {
  warning("Cox summary file not found; survival plotting patch skipped: ", COX_FOR_PLOT)
} else {
  cox <- read.csv(COX_FOR_PLOT, check.names = FALSE, stringsAsFactors = FALSE)
  msg("Survival plots use Cox table: ", COX_FOR_PLOT)

  needed_cox <- c("cancer", "endpoint", "score", "model", "HR", "p")
  miss <- setdiff(needed_cox, colnames(cox))
  if (length(miss) > 0) {
    warning("Cox file missing columns; survival patch skipped: ", paste(miss, collapse = ", "))
  } else {
    cox <- cox %>%
      mutate(
        endpoint = toupper(as.character(endpoint)),
        cancer = toupper(as.character(cancer)),
        HR = suppressWarnings(as.numeric(HR)),
        p = suppressWarnings(as.numeric(p))
      ) %>%
      filter(
        score %in% SCORE_COLS,
        endpoint %in% c("OS", "DFS"),
        is.finite(HR), HR > 0,
        is.finite(p), p > 0
      )

    plot_volcano_p_significance <- function(dat, sc, endpoint_name, model_name, model_label) {
      dd <- dat %>%
        filter(score == sc, endpoint == endpoint_name, model == model_name) %>%
        mutate(
          loge_HR = log(HR),
          neg_log10_p = -log10(p),
          sig = p < 0.05,
          p_text = vapply(p, format_p_short, character(1)),
          label = ifelse(sig, paste0(cancer, "\np=", p_text), "")
        )

      if (nrow(dd) == 0) return(invisible(NULL))

      p <- ggplot(dd, aes(x = loge_HR, y = neg_log10_p)) +
        geom_point(aes(color = sig), size = 2.7, alpha = 0.9) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
        ggrepel::geom_text_repel(
          aes(label = label),
          size = 3.0,
          max.overlaps = Inf,
          box.padding = 0.4,
          show.legend = FALSE
        ) +
        scale_color_manual(
          values = c("FALSE" = "grey75", "TRUE" = "#C74A37"),
          labels = c("P >= 0.05", "P < 0.05")
        ) +
        labs(
          title = paste0(score_label(sc), " - ", endpoint_name, " Cox volcano"),
          subtitle = model_label,
          x = expression(log[e](HR)),
          y = expression(-log[10](P~value)),
          color = NULL
        ) +
        theme_step(base_size = 10, title_size = 9.5,
                    axis_title_size = 9, axis_text_size = 8)

      base <- file.path(
        OUT_SURV,
        paste0(
          safe_name(sc), "_", endpoint_name, "_", model_name,
          "_cox_volcano_RED_is_P_lt_0.05_v1_3"
        )
      )
      save_gg_both(p, base, 5.3, 4.5)
      write.csv(dd, paste0(base, "_data.csv"), row.names = FALSE)
    }

    plot_old_style_bubble <- function(dat, sc, model_name, model_label) {
      dd <- dat %>%
        filter(score == sc, model == model_name, endpoint %in% c("OS", "DFS")) %>%
        mutate(
          loge_HR = log(HR),
          neg_log10_p = -log10(p),
          direction = ifelse(HR > 1, "Risk (HR>1)", "Protective (HR<1)"),
          sig_size = pmin(neg_log10_p, 15),
          p_text = vapply(p, format_p_short, character(1)),
          endpoint = factor(endpoint, levels = c("OS", "DFS")),
          cancer = factor(cancer, levels = sort(unique(cancer)))
        )

      if (nrow(dd) == 0) return(invisible(NULL))

      p <- ggplot(dd, aes(x = cancer, y = endpoint)) +
        geom_point(
          aes(size = sig_size, fill = neg_log10_p, color = direction),
          shape = 21, stroke = 0.7
        ) +
        geom_text(
          aes(label = ifelse(p < 0.05, paste0("p=", p_text), "")),
          size = 2.0,
          vjust = -1.1,
          color = "black",
          check_overlap = TRUE
        ) +
        scale_size_continuous(
          range = c(1.2, 9),
          name = expression(-log[10](P))
        ) +
        scale_fill_gradient(
          low = "#FFFFCC", high = "#D7301F",
          name = expression(-log[10](P))
        ) +
        scale_color_manual(
          values = c(
            "Risk (HR>1)" = "black",
            "Protective (HR<1)" = "#2C7FB8"
          ),
          name = "Direction"
        ) +
        labs(
          title = paste0(score_label(sc), " survival summary"),
          subtitle = paste0(model_label, "; significant bubbles are labeled with raw p values"),
          x = NULL,
          y = NULL
        ) +
        theme_step(base_size = 9.5, title_size = 9.5,
                    axis_title_size = 8.5, axis_text_size = 7.3) +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 9.5),
          plot.subtitle = element_text(size = 8)
        )

      w <- max(10, length(unique(dd$cancer)) * 0.33)
      base <- file.path(
        OUT_SURV,
        paste0(
          safe_name(sc),
          "_survival_endpoint_bubble_",
          ifelse(model_name == "univariate_score", "univariate", "multivariate"),
          "_OS_DFS_v1_3"
        )
      )
      save_gg_both(p, base, w, 3.3)
      write.csv(dd, paste0(base, "_data.csv"), row.names = FALSE)
    }

    for (sc in SCORE_COLS) {
      
      for (ep in c("OS", "DFS")) {
        plot_volcano_p_significance(
          cox, sc, ep,
          "univariate_score",
          "Univariate Cox: score only"
        )
        plot_volcano_p_significance(
          cox, sc, ep,
          "score_scaled",
          "Multivariable Cox: score + AGE + SEX + STAGE when available"
        )
      }

      plot_old_style_bubble(
        cox, sc,
        "univariate_score",
        "Univariate Cox: score only"
      )
      plot_old_style_bubble(
        cox, sc,
        "score_scaled",
        "Multivariable Cox: score term adjusted by AGE + SEX + STAGE where available"
      )

      if (any(cox$score == sc & cox$model == "score_scaled_CPE_AGE")) {
        dd_cpe <- cox %>% filter(score == sc, model == "score_scaled_CPE_AGE")
        dd_cpe <- dd_cpe %>%
          mutate(
            loge_HR = log(HR),
            neg_log10_p = -log10(p),
            direction = ifelse(HR > 1, "Risk (HR>1)", "Protective (HR<1)"),
            sig_size = pmin(neg_log10_p, 15),
            p_text = vapply(p, format_p_short, character(1)),
            endpoint = factor(endpoint, levels = c("OS", "DFS")),
            cancer = factor(cancer, levels = sort(unique(cancer)))
          )

        if (nrow(dd_cpe) > 0) {
          p3 <- ggplot(dd_cpe, aes(x = cancer, y = endpoint)) +
            geom_point(
              aes(size = sig_size, fill = neg_log10_p, color = direction),
              shape = 21, stroke = 0.7
            ) +
            geom_text(
              aes(label = ifelse(p < 0.05, paste0("p=", p_text), "")),
              size = 2.0, vjust = -1.1, color = "black", check_overlap = TRUE
            ) +
            scale_size_continuous(range = c(1.2, 9), name = expression(-log[10](P))) +
            scale_fill_gradient(low = "#FFFFCC", high = "#D7301F",
                                name = expression(-log[10](P))) +
            scale_color_manual(
              values = c("Risk (HR>1)" = "black", "Protective (HR<1)" = "#2C7FB8"),
              name = "Direction"
            ) +
            labs(
              title = paste0(score_label(sc), " survival summary"),
              subtitle = "Multivariable Cox: score + CPE purity + AGE",
              x = NULL, y = NULL
            ) +
            theme_step(base_size = 9.5, title_size = 9.5,
                        axis_title_size = 8.5, axis_text_size = 7.3) +
            theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

          base3 <- file.path(
            OUT_SURV,
            paste0(safe_name(sc),
                   "_survival_endpoint_bubble_multivariate_CPE_AGE_OS_DFS_v1_3")
          )
          save_gg_both(
            p3, base3,
            max(10, length(unique(dd_cpe$cancer)) * 0.33),
            3.3
          )
          write.csv(dd_cpe, paste0(base3, "_data.csv"), row.names = FALSE)
        }
      }
    }
  }
}

manifest <- data.frame(
  item = c(
    "Repaired clinical merged table",
    "Recomputed Cox table from repaired clinical data",
    "Clinical source column QC",
    "All-cancer clinical eligibility QC",
    "Previously empty cancer QC",
    "All-cancer clinical faceted summary folder",
    "Corrected volcano + old-style survival bubble folder"
  ),
  path = c(
    REPAIRED_MERGED,
    REPAIRED_COX_FILE,
    file.path(OUT_CLIN, "clinical_source_detected_columns_v1_3.csv"),
    file.path(OUT_CLIN, "ALL_cancers_clinical_eligibility_QC_v1_3.csv"),
    file.path(OUT_CLIN, "Previously_empty_cancers_clinical_QC_v1_3.csv"),
    OUT_CLIN_SUM,
    OUT_SURV
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(OUT_LOG, "functional_subset_posthoc_fix_v1_3_manifest.csv"),
  row.names = FALSE
)

msg("DONE - functional-subset post-hoc fix v1.3")
msg("Clinical summary folder: ", OUT_CLIN_SUM)
msg("Survival corrected plots folder: ", OUT_SURV)
msg("Repaired merged clinical table: ", REPAIRED_MERGED)
