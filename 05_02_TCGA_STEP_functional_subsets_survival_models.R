# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(123)

RESULT_DIR <- NA_character_  

MIN_MODEL_N <- 20L
MIN_EVENTS  <- 5L
PNG_DPI <- 300

LOW_COLOR  <- "#006FB9"
HIGH_COLOR <- "#E60011"

SIG_RED <- "#D7301F"
NS_GREY <- "grey72"
PROTECTIVE_BLUE <- "#2C7FB8"
RISK_BLACK <- "black"

required_pkgs <- c(
  "survival", "ggplot2", "dplyr", "tidyr",
  "stringr", "scales", "ggrepel", "pheatmap"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
  library(ggrepel)
  library(pheatmap)
})

get_script_dir <- function() {
  p <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)

  if (!is.null(p) && length(p) == 1 && !is.na(p) && nzchar(p)) {
    return(dirname(normalizePath(p, mustWork = FALSE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p2 <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p2)) {
      return(dirname(normalizePath(p2, mustWork = FALSE)))
    }
  }

  normalizePath(getwd(), mustWork = FALSE)
}

has_required_files <- function(x) {
  is.character(x) && length(x) == 1 &&
    !is.na(x) && nzchar(x) &&
    file.exists(file.path(
      x,
      "08_survival_optimal_cutoff_v1_4",
      "00_cutoff_QC",
      "ALL_cancers_OPTIMAL_cutoff_summary_and_eligibility_v1_4.csv"
    )) &&
    file.exists(file.path(
      x,
      "09_pan_cancer_summary",
      "ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv"
    ))
}

resolve_result_dir <- function(manual = RESULT_DIR) {
  if (has_required_files(manual)) {
    return(normalizePath(manual, mustWork = TRUE))
  }

  script_dir <- get_script_dir()
  wd <- normalizePath(getwd(), mustWork = FALSE)

  direct_candidates <- unique(c(
    wd, script_dir, dirname(wd), dirname(script_dir)
  ))

  for (x in direct_candidates) {
    if (has_required_files(x)) {
      message("Detected result folder directly: ", x)
      return(normalizePath(x, mustWork = TRUE))
    }
  }

  roots <- unique(c(
    wd, script_dir, dirname(wd), dirname(script_dir)
  ))
  candidates <- character(0)

  for (rr in roots) {
    if (!dir.exists(rr)) next

    subs <- list.dirs(rr, full.names = TRUE, recursive = FALSE)

    subs <- subs[
      grepl(
        "^TCGA_STEP_functional_subsets_MYC7_Exosome_v1_[0-9]+_",
        basename(subs)
      )
    ]

    if (length(subs) > 0) {
      ok <- vapply(subs, has_required_files, logical(1))
      candidates <- c(candidates, subs[ok])
    }
  }

  candidates <- unique(candidates)

  if (length(candidates) == 0) {
    stop(
      "\nCannot auto-detect a completed functional-subset result folder.\n",
      "Please set RESULT_DIR manually at the top of this script.\n\n",
      "Required files:\n",
      "1) 08_survival_optimal_cutoff_v1_4/00_cutoff_QC/",
      "ALL_cancers_OPTIMAL_cutoff_summary_and_eligibility_v1_4.csv\n",
      "2) 09_pan_cancer_summary/",
      "ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv\n"
    )
  }

  mt <- file.info(candidates)$mtime
  chosen <- candidates[order(mt, decreasing = TRUE)][1]

  message("Auto-detected latest result folder: ", chosen)
  normalizePath(chosen, mustWork = TRUE)
}

RESULT_DIR <- resolve_result_dir(RESULT_DIR)

CUTOFF_FILE <- file.path(
  RESULT_DIR,
  "08_survival_optimal_cutoff_v1_4",
  "00_cutoff_QC",
  "ALL_cancers_OPTIMAL_cutoff_summary_and_eligibility_v1_4.csv"
)

INPUT_REPAIRED <- file.path(
  RESULT_DIR,
  "09_pan_cancer_summary",
  "ALL_cancers_MYC7_Exosome_scores_clinical_merged_v1_3_repaired.csv"
)

INPUT_V12 <- file.path(
  RESULT_DIR,
  "09_pan_cancer_summary",
  "ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv"
)

INPUT_FILE <- if (file.exists(INPUT_REPAIRED)) {
  INPUT_REPAIRED
} else {
  INPUT_V12
}

OUT_ROOT <- file.path(
  RESULT_DIR,
  "08_survival_optimal_cutoff_v1_4",
  "01_downstream_survival_summary_v1_5"
)

OUT_TABLE <- file.path(OUT_ROOT, "00_tables")
OUT_VOLCANO <- file.path(OUT_ROOT, "01_survival_volcano")
OUT_FOREST <- file.path(OUT_ROOT, "02_survival_forest")
OUT_BUBBLE <- file.path(OUT_ROOT, "03_survival_bubble_summary")
OUT_HEATMAP <- file.path(OUT_ROOT, "04_HR_heatmaps")
OUT_FIGREADY <- file.path(OUT_ROOT, "05_figure_ready")
OUT_LOG <- file.path(RESULT_DIR, "00_logs")

for (d in c(
  OUT_ROOT, OUT_TABLE, OUT_VOLCANO, OUT_FOREST,
  OUT_BUBBLE, OUT_HEATMAP, OUT_FIGREADY, OUT_LOG
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(
  OUT_LOG,
  "optimalCutoff_survival_downstream_v1_5.log"
)

msg <- function(...) {
  txt <- paste0(...)
  message(txt)
  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ", txt, "\n",
    file = LOG_FILE, append = TRUE, sep = ""
  )
}

msg("RESULT_DIR = ", RESULT_DIR)
msg("INPUT_FILE = ", INPUT_FILE)
msg("CUTOFF_FILE = ", CUTOFF_FILE)
msg("OUT_ROOT = ", OUT_ROOT)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

format_p <- function(p, prefix = "p = ") {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) {
    return(paste0(prefix, "NA"))
  }

  if (p < 1e-4) {
    return(paste0(prefix, formatC(p, format = "e", digits = 2)))
  }

  paste0(prefix, signif(p, 3))
}

format_p_short <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) {
    return("NA")
  }

  if (p < 1e-4) {
    return(formatC(p, format = "e", digits = 1))
  }

  format(signif(p, 2), trim = TRUE)
}

theme_step <- function(
    base_size = 10,
    title_size = 10,
    axis_title_size = 9,
    axis_text_size = 8,
    strip_size = 8
) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = title_size,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = max(7, title_size - 1),
        hjust = 0.5
      ),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(
        size = axis_text_size,
        color = "black"
      ),
      strip.text = element_text(
        size = strip_size,
        face = "bold"
      ),
      legend.title = element_text(
        size = max(7, axis_title_size - 1)
      ),
      legend.text = element_text(
        size = max(7, axis_text_size - 1)
      ),
      panel.grid.minor = element_blank()
    )
}

save_gg_both <- function(
    p,
    file_base,
    width = 6.5,
    height = 5.2
) {
  print(p)

  ggsave(
    paste0(file_base, ".pdf"),
    p,
    width = width,
    height = height,
    limitsize = FALSE,
    useDingbats = FALSE
  )

  ggsave(
    paste0(file_base, ".png"),
    p,
    width = width,
    height = height,
    dpi = PNG_DPI,
    limitsize = FALSE
  )
}

normalize_event_vector <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    out <- suppressWarnings(as.numeric(x))
    out[!(out %in% c(0, 1))] <- NA_real_
    return(out)
  }

  xx <- toupper(trimws(as.character(x)))
  out <- rep(NA_real_, length(xx))

  out[xx %in% c(
    "1", "TRUE", "DECEASED", "DEAD", "EVENT", "YES",
    "PROGRESSED", "RECURRENCE", "RECURRED", "RELAPSED"
  )] <- 1

  out[xx %in% c(
    "0", "FALSE", "LIVING", "ALIVE", "CENSORED", "NO",
    "DISEASEFREE", "DISEASE_FREE", "NO_RECURRENCE"
  )] <- 0

  out[grepl("^1", xx)] <- 1
  out[grepl("^0", xx)] <- 0

  out
}

find_col <- function(candidates, nms) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) hit[1] else NA_character_
}

SCORE_MAP <- c(
  "Hallmark_MYC_targets_v1" =
    "MYC targets v1 (7 genes)",
  "GO_extracellular_exsosome" =
    "Extracellular exosome (23 genes)"
)

score_label <- function(x) {
  z <- unname(SCORE_MAP[x])
  z[is.na(z)] <- x[is.na(z)]
  z
}

if (!file.exists(INPUT_FILE)) {
  stop("Cannot find merged score/clinical table: ", INPUT_FILE)
}

if (!file.exists(CUTOFF_FILE)) {
  stop("Cannot find v1.4 optimal-cutoff table: ", CUTOFF_FILE)
}

sdat <- read.csv(
  INPUT_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cutoff_df <- read.csv(
  CUTOFF_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cancer_col <- find_col(
  c(
    ".cancer",
    "cancer",
    "cancer_loop",
    "cancer_from_file",
    "cancer_type",
    "CANCER"
  ),
  colnames(sdat)
)

if (is.na(cancer_col)) {
  stop("Cannot detect cancer-type column in merged table.")
}

sdat$.cancer_v15 <- toupper(
  trimws(as.character(sdat[[cancer_col]]))
)

required_cutoff_cols <- c(
  "cancer", "endpoint", "score",
  "cutoff", "eligible"
)

missing_cutoff_cols <- setdiff(
  required_cutoff_cols,
  colnames(cutoff_df)
)

if (length(missing_cutoff_cols) > 0) {
  stop(
    "v1.4 cut-off table is missing columns: ",
    paste(missing_cutoff_cols, collapse = ", ")
  )
}

cutoff_df$cancer <- toupper(
  trimws(as.character(cutoff_df$cancer))
)

cutoff_df$endpoint <- toupper(
  trimws(as.character(cutoff_df$endpoint))
)

cutoff_df$cutoff <- suppressWarnings(
  as.numeric(cutoff_df$cutoff)
)

if (is.logical(cutoff_df$eligible)) {
  cutoff_df$.eligible_v15 <- cutoff_df$eligible
} else {
  cutoff_df$.eligible_v15 <- toupper(
    as.character(cutoff_df$eligible)
  ) %in% c("TRUE", "T", "1", "YES")
}

cutoff_use <- cutoff_df %>%
  filter(
    .eligible_v15,
    is.finite(cutoff),
    score %in% names(SCORE_MAP)
  ) %>%
  distinct(
    cancer,
    endpoint,
    score,
    .keep_all = TRUE
  )

if (nrow(cutoff_use) == 0) {
  stop("No eligible optimal cut-offs were found in v1.4.")
}

write.csv(
  cutoff_use,
  file.path(
    OUT_TABLE,
    "v1_4_optimal_cutoffs_used_by_v1_5.csv"
  ),
  row.names = FALSE
)

msg(
  "Eligible cut-offs loaded: ",
  nrow(cutoff_use)
)

endpoint_defs <- list(
  OS = c(
    time = find_col(
      c(
        "OS_MONTHS", "OS_MONTH",
        "OS.time", "OS_TIME", "OS_DAYS"
      ),
      colnames(sdat)
    ),
    event = find_col(
      c(
        "OS_EVENT", "OS_STATUS",
        "OS", "OS.event"
      ),
      colnames(sdat)
    )
  ),
  DFS = c(
    time = find_col(
      c(
        "DFS_MONTHS", "DFS_MONTH",
        "DFS.time", "DFS_TIME", "DFS_DAYS"
      ),
      colnames(sdat)
    ),
    event = find_col(
      c(
        "DFS_EVENT", "DFS_STATUS",
        "DFS", "DFS.event"
      ),
      colnames(sdat)
    )
  ),
  PFS = c(
    time = find_col(
      c(
        "PFS_MONTHS", "PFS_MONTH",
        "PFS.time", "PFS_TIME", "PFS_DAYS"
      ),
      colnames(sdat)
    ),
    event = find_col(
      c(
        "PFS_EVENT", "PFS_STATUS",
        "PFS", "PFS.event"
      ),
      colnames(sdat)
    )
  ),
  RFS = c(
    time = find_col(
      c(
        "RFS_MONTHS", "RFS_MONTH",
        "RFS.time", "RFS_TIME", "RFS_DAYS"
      ),
      colnames(sdat)
    ),
    event = find_col(
      c(
        "RFS_EVENT", "RFS_STATUS",
        "RFS", "RFS.event"
      ),
      colnames(sdat)
    )
  )
)

endpoint_available <- names(endpoint_defs)[
  vapply(
    endpoint_defs,
    function(x) all(!is.na(x)) && all(nzchar(x)),
    logical(1)
  )
]

msg(
  "Detected survival endpoints: ",
  paste(endpoint_available, collapse = ", ")
)

prepare_covariate_frame <- function(d) {
  out <- data.frame(
    AGE = if ("AGE" %in% colnames(d)) {
      suppressWarnings(as.numeric(d$AGE))
    } else {
      NA_real_
    },

    SEX = if ("SEX" %in% colnames(d)) {
      trimws(as.character(d$SEX))
    } else {
      NA_character_
    },

    STAGE_NUM = if ("STAGE_NUM" %in% colnames(d)) {
      suppressWarnings(as.numeric(d$STAGE_NUM))
    } else if ("STAGE" %in% colnames(d)) {
      
      st <- toupper(trimws(as.character(d$STAGE)))
      stn <- rep(NA_real_, length(st))
      stn[st %in% c("0", "STAGE 0")] <- 0
      stn[grepl("^I($|A|B|C)|^1", st)] <- 1
      stn[grepl("^II($|A|B|C)|^2", st)] <- 2
      stn[grepl("^III($|A|B|C)|^3", st)] <- 3
      stn[grepl("^IV($|A|B|C)|^4", st)] <- 4
      stn
    } else {
      NA_real_
    },

    CPE = if ("CPE" %in% colnames(d)) {
      suppressWarnings(as.numeric(d$CPE))
    } else {
      NA_real_
    },

    stringsAsFactors = FALSE
  )

  out$SEX[toupper(out$SEX) %in% c(
    "", "NA", "N/A", "UNKNOWN",
    "[NOT AVAILABLE]", "NOT REPORTED"
  )] <- NA_character_

  out$SEX[toupper(out$SEX) %in% c(
    "M", "MALE"
  )] <- "Male"

  out$SEX[toupper(out$SEX) %in% c(
    "F", "FEMALE"
  )] <- "Female"

  out
}

covariate_usable <- function(v, type = c("numeric", "factor")) {
  type <- match.arg(type)

  if (type == "numeric") {
    x <- suppressWarnings(as.numeric(v))
    ok <- is.finite(x)

    return(
      sum(ok) >= MIN_MODEL_N &&
        length(unique(x[ok])) >= 2
    )
  }

  x <- trimws(as.character(v))
  x[x == ""] <- NA_character_
  tab <- table(x, useNA = "no")

  length(tab) >= 2 &&
    sum(tab) >= MIN_MODEL_N &&
    all(tab >= 2)
}

choose_best_clinical_covariates <- function(dmod) {
  candidates <- character(0)

  if (covariate_usable(dmod$AGE, "numeric")) {
    candidates <- c(candidates, "AGE")
  }

  if (covariate_usable(dmod$SEX, "factor")) {
    candidates <- c(candidates, "SEX")
  }

  if (covariate_usable(dmod$STAGE_NUM, "numeric")) {
    candidates <- c(candidates, "STAGE_NUM")
  }

  subsets <- list(character(0))

  if (length(candidates) > 0) {
    subsets <- unlist(
      lapply(
        seq_along(candidates),
        function(k) {
          combn(
            candidates,
            k,
            simplify = FALSE
          )
        }
      ),
      recursive = FALSE
    )
  }

  subsets <- subsets[
    order(
      vapply(subsets, length, integer(1)),
      decreasing = TRUE
    )
  ]

  for (covs in subsets) {
    needed <- c(
      "time", "event", "group",
      covs
    )

    dd <- dmod[
      complete.cases(
        dmod[, needed, drop = FALSE]
      ),
      ,
      drop = FALSE
    ]

    if (nrow(dd) < MIN_MODEL_N) next
    if (sum(dd$event == 1, na.rm = TRUE) < MIN_EVENTS) next
    if (length(unique(dd$group)) < 2) next

    if ("SEX" %in% covs) {
      st <- table(dd$SEX)
      if (length(st) < 2 || any(st < 2)) next
    }

    if ("AGE" %in% covs &&
        length(unique(dd$AGE)) < 2) next

    if ("STAGE_NUM" %in% covs &&
        length(unique(dd$STAGE_NUM)) < 2) next

    return(
      list(
        covariates = covs,
        data = dd
      )
    )
  }

  NULL
}

fit_one_optimal_cutoff <- function(
    d,
    cancer,
    endpoint,
    score_col,
    cutoff,
    time_col,
    event_col,
    maxstat_statistic = NA_real_,
    minprop = NA_real_
) {
  if (!all(
    c(score_col, time_col, event_col) %in%
      colnames(d)
  )) {
    return(NULL)
  }

  tt <- suppressWarnings(
    as.numeric(d[[time_col]])
  )

  ee <- normalize_event_vector(
    d[[event_col]]
  )

  ss <- suppressWarnings(
    as.numeric(d[[score_col]])
  )

  ok <- is.finite(tt) &
    tt > 0 &
    ee %in% c(0, 1) &
    is.finite(ss)

  if (sum(ok) < MIN_MODEL_N) {
    return(NULL)
  }

  dcov <- prepare_covariate_frame(
    d[ok, , drop = FALSE]
  )

  dd <- data.frame(
    time = tt[ok],
    event = ee[ok],
    score = ss[ok],
    dcov,
    stringsAsFactors = FALSE
  )

  dd$group <- ifelse(
    dd$score <= cutoff,
    "Low",
    "High"
  )

  dd$group <- factor(
    dd$group,
    levels = c("Low", "High")
  )

  if (length(unique(dd$group)) < 2) {
    return(NULL)
  }

  n_low <- sum(dd$group == "Low")
  n_high <- sum(dd$group == "High")
  n_events <- sum(dd$event == 1)

  if (min(n_low, n_high) < 2 ||
      n_events < MIN_EVENTS) {
    return(NULL)
  }

  lr <- tryCatch(
    survival::survdiff(
      survival::Surv(time, event) ~ group,
      data = dd
    ),
    error = function(e) NULL
  )

  logrank_p <- if (is.null(lr)) {
    NA_real_
  } else {
    1 - pchisq(
      lr$chisq,
      df = length(lr$n) - 1
    )
  }

  rows <- list()
  ri <- 1L

  fit_uni <- tryCatch(
    survival::coxph(
      survival::Surv(time, event) ~ group,
      data = dd
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_uni)) {
    sm <- summary(fit_uni)

    coef_name <- grep(
      "^groupHigh$",
      rownames(sm$coefficients),
      value = TRUE
    )

    if (length(coef_name) == 1) {
      cn <- coef_name[1]

      rows[[ri]] <- data.frame(
        cancer = cancer,
        endpoint = endpoint,
        score = score_col,
        score_label = score_label(score_col),

        cutoff = cutoff,
        cutoff_method =
          "v1.4 surv_cutpoint/maxstat",
        maxstat_statistic =
          maxstat_statistic,
        minprop = minprop,

        comparison = "High_vs_Low",
        model = "univariate_group",
        model_label =
          "Univariate Cox: optimal-cutoff High vs Low",
        covariates = "group only",

        n = nrow(dd),
        n_events = n_events,
        n_low = n_low,
        n_high = n_high,

        logrank_p =
          logrank_p,

        HR =
          sm$coefficients[
            cn, "exp(coef)"
          ],
        p =
          sm$coefficients[
            cn, "Pr(>|z|)"
          ],
        lower95 =
          sm$conf.int[
            cn, "lower .95"
          ],
        upper95 =
          sm$conf.int[
            cn, "upper .95"
          ],

        stringsAsFactors = FALSE
      )

      ri <- ri + 1L
    }
  }

  clin_pick <- choose_best_clinical_covariates(dd)

  if (!is.null(clin_pick)) {
    covs <- clin_pick$covariates
    dm <- clin_pick$data

    if (length(covs) > 0) {
      if ("SEX" %in% covs) {
        dm$SEX <- factor(
          dm$SEX,
          levels = c("Male", "Female")
        )
      }

      rhs <- paste(
        c("group", covs),
        collapse = " + "
      )

      fml <- as.formula(
        paste0(
          "survival::Surv(time, event) ~ ",
          rhs
        )
      )

      fit_mv <- tryCatch(
        survival::coxph(
          fml,
          data = dm
        ),
        error = function(e) NULL
      )

      if (!is.null(fit_mv)) {
        sm <- summary(fit_mv)

        coef_name <- grep(
          "^groupHigh$",
          rownames(sm$coefficients),
          value = TRUE
        )

        if (length(coef_name) == 1) {
          cn <- coef_name[1]

          rows[[ri]] <- data.frame(
            cancer = cancer,
            endpoint = endpoint,
            score = score_col,
            score_label = score_label(score_col),

            cutoff = cutoff,
            cutoff_method =
              "v1.4 surv_cutpoint/maxstat",
            maxstat_statistic =
              maxstat_statistic,
            minprop = minprop,

            comparison = "High_vs_Low",
            model = "multivariable_available_clinical",
            model_label =
              "Multivariable Cox: optimal-cutoff High vs Low + available clinical covariates",
            covariates =
              paste(
                c("group", covs),
                collapse = " + "
              ),

            n = nrow(dm),
            n_events =
              sum(dm$event == 1),
            n_low =
              sum(dm$group == "Low"),
            n_high =
              sum(dm$group == "High"),

            logrank_p =
              logrank_p,

            HR =
              sm$coefficients[
                cn, "exp(coef)"
              ],
            p =
              sm$coefficients[
                cn, "Pr(>|z|)"
              ],
            lower95 =
              sm$conf.int[
                cn, "lower .95"
              ],
            upper95 =
              sm$conf.int[
                cn, "upper .95"
              ],

            stringsAsFactors = FALSE
          )

          ri <- ri + 1L
        }
      }
    }
  }

  cpe_ok <- complete.cases(
    dd[, c(
      "time", "event", "group",
      "CPE", "AGE"
    ), drop = FALSE]
  )

  d2 <- dd[cpe_ok, , drop = FALSE]

  if (
    nrow(d2) >= MIN_MODEL_N &&
    sum(d2$event == 1) >= MIN_EVENTS &&
    length(unique(d2$group)) >= 2 &&
    length(unique(d2$CPE)) >= 2 &&
    length(unique(d2$AGE)) >= 2
  ) {
    d2$CPE_scaled <- as.numeric(
      scale(d2$CPE)
    )

    fit_cpe <- tryCatch(
      survival::coxph(
        survival::Surv(time, event) ~
          group + CPE_scaled + AGE,
        data = d2
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_cpe)) {
      sm <- summary(fit_cpe)

      coef_name <- grep(
        "^groupHigh$",
        rownames(sm$coefficients),
        value = TRUE
      )

      if (length(coef_name) == 1) {
        cn <- coef_name[1]

        rows[[ri]] <- data.frame(
          cancer = cancer,
          endpoint = endpoint,
          score = score_col,
          score_label = score_label(score_col),

          cutoff = cutoff,
          cutoff_method =
            "v1.4 surv_cutpoint/maxstat",
          maxstat_statistic =
            maxstat_statistic,
          minprop = minprop,

          comparison = "High_vs_Low",
          model = "multivariable_CPE_AGE",
          model_label =
            "Multivariable Cox: optimal-cutoff High vs Low + CPE + AGE",
          covariates =
            "group + CPE + AGE",

          n = nrow(d2),
          n_events =
            sum(d2$event == 1),
          n_low =
            sum(d2$group == "Low"),
          n_high =
            sum(d2$group == "High"),

          logrank_p =
            logrank_p,

          HR =
            sm$coefficients[
              cn, "exp(coef)"
            ],
          p =
            sm$coefficients[
              cn, "Pr(>|z|)"
            ],
          lower95 =
            sm$conf.int[
              cn, "lower .95"
            ],
          upper95 =
            sm$conf.int[
              cn, "upper .95"
            ],

          stringsAsFactors = FALSE
        )
      }
    }
  }

  bind_rows(rows)
}

all_models <- list()
mi <- 1L

for (i in seq_len(nrow(cutoff_use))) {
  ca <- cutoff_use$cancer[i]
  ep <- cutoff_use$endpoint[i]
  sc <- cutoff_use$score[i]
  co <- cutoff_use$cutoff[i]

  if (!(ep %in% endpoint_available)) {
    next
  }

  if (!(sc %in% colnames(sdat))) {
    next
  }

  dc <- sdat[
    sdat$.cancer_v15 == ca,
    ,
    drop = FALSE
  ]

  if (nrow(dc) == 0) next

  maxstat_i <- if (
    "maxstat_statistic" %in%
      colnames(cutoff_use)
  ) {
    suppressWarnings(
      as.numeric(
        cutoff_use$maxstat_statistic[i]
      )
    )
  } else {
    NA_real_
  }

  minprop_i <- if (
    "minprop" %in%
      colnames(cutoff_use)
  ) {
    suppressWarnings(
      as.numeric(
        cutoff_use$minprop[i]
      )
    )
  } else {
    NA_real_
  }

  tab <- tryCatch(
    fit_one_optimal_cutoff(
      d = dc,
      cancer = ca,
      endpoint = ep,
      score_col = sc,
      cutoff = co,
      time_col =
        endpoint_defs[[ep]]["time"],
      event_col =
        endpoint_defs[[ep]]["event"],
      maxstat_statistic =
        maxstat_i,
      minprop =
        minprop_i
    ),
    error = function(e) {
      msg(
        "Model failed: ",
        ca, " / ", ep, " / ", sc,
        " | ", conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(tab) &&
      nrow(tab) > 0) {
    all_models[[mi]] <- tab
    mi <- mi + 1L
  }
}

cox_all <- bind_rows(all_models)

if (nrow(cox_all) == 0) {
  stop(
    "No optimal-cutoff Cox models could be fitted. ",
    "Check merged survival data and v1.4 cutoff QC."
  )
}

cox_all <- cox_all %>%
  group_by(
    score,
    endpoint,
    model
  ) %>%
  mutate(
    FDR = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    direction = ifelse(
      HR > 1,
      "Risk (HR>1)",
      "Protective (HR<1)"
    ),
    raw_p_significance = ifelse(
      p < 0.05,
      "P < 0.05",
      "P >= 0.05"
    ),
    FDR_significance = ifelse(
      FDR < 0.05,
      "FDR < 0.05",
      "FDR >= 0.05"
    ),
    log2HR = log2(HR),
    neglog10P = -log10(
      pmax(p, 1e-300)
    ),
    neglog10FDR = -log10(
      pmax(FDR, 1e-300)
    )
  )

write.csv(
  cox_all,
  file.path(
    OUT_TABLE,
    "ALL_optimalCutoff_group_Cox_models_v1_5.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_all %>%
    filter(model == "univariate_group"),
  file.path(
    OUT_TABLE,
    "ALL_optimalCutoff_univariate_Cox_v1_5.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_all %>%
    filter(model == "multivariable_available_clinical"),
  file.path(
    OUT_TABLE,
    "ALL_optimalCutoff_multivariable_availableClinical_Cox_v1_5.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_all %>%
    filter(model == "multivariable_CPE_AGE"),
  file.path(
    OUT_TABLE,
    "ALL_optimalCutoff_multivariable_CPE_AGE_Cox_v1_5.csv"
  ),
  row.names = FALSE
)

plot_optimal_volcano <- function(
    dat,
    sc,
    ep,
    model_name
) {
  dd <- dat %>%
    filter(
      score == sc,
      endpoint == ep,
      model == model_name,
      is.finite(HR),
      HR > 0,
      is.finite(p),
      p > 0
    )

  if (nrow(dd) == 0) {
    return(invisible(NULL))
  }

  dd <- dd %>%
    mutate(
      sig = p < 0.05,
      label = ifelse(
        p < 0.05,
        paste0(
          cancer,
          "\np=",
          vapply(
            p,
            format_p_short,
            character(1)
          )
        ),
        ""
      )
    )

  model_title <- unique(
    dd$model_label
  )[1]

  p <- ggplot(
    dd,
    aes(
      x = log2HR,
      y = neglog10P
    )
  ) +
    geom_point(
      aes(color = sig),
      size = 2.8,
      alpha = 0.90
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 2.9,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.20,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "FALSE" = NS_GREY,
        "TRUE" = SIG_RED
      ),
      labels = c(
        "P >= 0.05",
        "P < 0.05"
      )
    ) +
    labs(
      title = paste0(
        score_label(sc),
        " - ", ep,
        " optimal-cutoff survival volcano"
      ),
      subtitle = paste0(
        model_title,
        "\nHR = High vs Low; RED indicates raw Cox p < 0.05"
      ),
      x = expression(log[2](HR)),
      y = expression(-log[10](P)),
      color = NULL
    ) +
    theme_step(
      base_size = 10,
      title_size = 9.5,
      axis_title_size = 9,
      axis_text_size = 8
    ) +
    theme(
      legend.position = "top"
    )

  base <- file.path(
    OUT_VOLCANO,
    paste0(
      safe_name(sc), "_",
      ep, "_",
      safe_name(model_name),
      "_OPTIMAL_CUTOFF_survival_volcano"
    )
  )

  save_gg_both(
    p,
    base,
    width = 5.7,
    height = 4.8
  )

  write.csv(
    dd,
    paste0(
      base,
      "_data.csv"
    ),
    row.names = FALSE
  )
}

for (sc in unique(cox_all$score)) {
  for (ep in unique(cox_all$endpoint)) {
    for (md in unique(cox_all$model)) {
      plot_optimal_volcano(
        cox_all,
        sc,
        ep,
        md
      )
    }
  }
}

plot_optimal_forest <- function(
    dat,
    sc,
    ep,
    model_name
) {
  dd <- dat %>%
    filter(
      score == sc,
      endpoint == ep,
      model == model_name,
      is.finite(HR),
      HR > 0,
      is.finite(lower95),
      is.finite(upper95),
      lower95 > 0,
      upper95 > 0
    ) %>%
    arrange(HR)

  if (nrow(dd) == 0) {
    return(invisible(NULL))
  }

  dd$cancer <- factor(
    dd$cancer,
    levels = dd$cancer
  )

  dd$sig <- factor(
    ifelse(
      dd$p < 0.05,
      "P < 0.05",
      "P >= 0.05"
    ),
    levels = c(
      "P < 0.05",
      "P >= 0.05"
    )
  )

  model_title <- unique(
    dd$model_label
  )[1]

  p <- ggplot(
    dd,
    aes(
      x = HR,
      y = cancer
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.40
    ) +
    geom_errorbarh(
      aes(
        xmin = lower95,
        xmax = upper95
      ),
      height = 0.16,
      linewidth = 0.38
    ) +
    geom_point(
      aes(fill = sig),
      shape = 21,
      size = 2.5,
      color = "black",
      stroke = 0.35
    ) +
    scale_fill_manual(
      values = c(
        "P < 0.05" = SIG_RED,
        "P >= 0.05" = NS_GREY
      )
    ) +
    scale_x_log10() +
    labs(
      title = paste0(
        score_label(sc),
        " - ", ep,
        " optimal-cutoff Cox forest"
      ),
      subtitle = paste0(
        model_title,
        "\nHR compares High vs Low using each cancer's v1.4 endpoint-specific cut-off"
      ),
      x = "Hazard ratio (High vs Low; log scale)",
      y = NULL,
      fill = NULL
    ) +
    theme_step(
      base_size = 10,
      title_size = 9.5,
      axis_title_size = 9,
      axis_text_size = 7.6
    ) +
    theme(
      legend.position = "top"
    )

  h <- max(
    5.5,
    0.28 * nrow(dd) + 1.9
  )

  base <- file.path(
    OUT_FOREST,
    paste0(
      safe_name(sc), "_",
      ep, "_",
      safe_name(model_name),
      "_OPTIMAL_CUTOFF_forest"
    )
  )

  save_gg_both(
    p,
    base,
    width = 6.5,
    height = h
  )

  write.csv(
    dd,
    paste0(
      base,
      "_data.csv"
    ),
    row.names = FALSE
  )
}

for (sc in unique(cox_all$score)) {
  for (ep in unique(cox_all$endpoint)) {
    for (md in unique(cox_all$model)) {
      plot_optimal_forest(
        cox_all,
        sc,
        ep,
        md
      )
    }
  }
}

plot_optimal_bubble <- function(
    dat,
    sc,
    model_name
) {
  dd <- dat %>%
    filter(
      score == sc,
      model == model_name,
      endpoint %in% c(
        "OS", "DFS"
      ),
      is.finite(HR),
      HR > 0,
      is.finite(p),
      p > 0
    ) %>%
    mutate(
      endpoint = factor(
        endpoint,
        levels = c(
          "OS", "DFS"
        )
      ),
      direction = ifelse(
        HR > 1,
        "Risk (HR>1)",
        "Protective (HR<1)"
      ),
      sig_size = pmin(
        neglog10P,
        15
      ),
      p_text = vapply(
        p,
        format_p_short,
        character(1)
      )
    )

  if (nrow(dd) == 0) {
    return(invisible(NULL))
  }

  cancer_order <- sort(
    unique(dd$cancer)
  )

  dd$cancer <- factor(
    dd$cancer,
    levels = cancer_order
  )

  model_title <- unique(
    dd$model_label
  )[1]

  p <- ggplot(
    dd,
    aes(
      x = cancer,
      y = endpoint
    )
  ) +
    geom_point(
      aes(
        size = sig_size,
        fill = neglog10P,
        color = direction
      ),
      shape = 21,
      stroke = 0.75
    ) +
    geom_text(
      aes(
        label = ifelse(
          p < 0.05,
          paste0(
            "p=",
            p_text
          ),
          ""
        )
      ),
      size = 1.9,
      vjust = -1.15,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_size_continuous(
      range = c(
        1.2, 9
      ),
      name = expression(
        -log[10](P)
      )
    ) +
    scale_fill_gradient(
      low = "#FFFFCC",
      high = SIG_RED,
      name = expression(
        -log[10](P)
      )
    ) +
    scale_color_manual(
      values = c(
        "Risk (HR>1)" =
          RISK_BLACK,
        "Protective (HR<1)" =
          PROTECTIVE_BLUE
      ),
      name = "Direction"
    ) +
    labs(
      title = paste0(
        score_label(sc),
        " optimal-cutoff survival summary"
      ),
      subtitle = paste0(
        model_title,
        "\nEach cancer uses its endpoint-specific v1.4 optimal cut-off; HR = High vs Low"
      ),
      x = NULL,
      y = NULL
    ) +
    theme_step(
      base_size = 9.5,
      title_size = 9.5,
      axis_title_size = 8.5,
      axis_text_size = 7.2
    ) +
    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.minor =
        element_blank(),
      legend.position = "right"
    )

  w <- max(
    10,
    length(cancer_order) * 0.36
  )

  base <- file.path(
    OUT_BUBBLE,
    paste0(
      safe_name(sc),
      "_OPTIMAL_CUTOFF_survival_bubble_",
      safe_name(model_name),
      "_OS_DFS"
    )
  )

  save_gg_both(
    p,
    base,
    width = w,
    height = 3.7
  )

  write.csv(
    dd,
    paste0(
      base,
      "_data.csv"
    ),
    row.names = FALSE
  )
}

for (sc in unique(cox_all$score)) {
  for (md in unique(cox_all$model)) {
    plot_optimal_bubble(
      cox_all,
      sc,
      md
    )
  }
}

fdr_star <- function(x) {
  out <- rep("", length(x))

  out[
    is.finite(x) &
      x < 0.001
  ] <- "***"

  out[
    is.finite(x) &
      x >= 0.001 &
      x < 0.01
  ] <- "**"

  out[
    is.finite(x) &
      x >= 0.01 &
      x < 0.05
  ] <- "*"

  out
}

plot_optimal_hr_heatmap <- function(
    dat,
    sc,
    model_name
) {
  dd <- dat %>%
    filter(
      score == sc,
      model == model_name,
      endpoint %in% c(
        "OS", "DFS"
      ),
      is.finite(log2HR)
    ) %>%
    mutate(
      sig_star = fdr_star(FDR)
    )

  if (nrow(dd) == 0) {
    return(invisible(NULL))
  }

  cancers <- sort(
    unique(dd$cancer)
  )

  endpoints <- intersect(
    c("OS", "DFS"),
    unique(dd$endpoint)
  )

  full_grid <- expand.grid(
    endpoint = endpoints,
    cancer = cancers,
    stringsAsFactors = FALSE
  )

  dd2 <- full_grid %>%
    left_join(
      dd %>%
        select(
          cancer,
          endpoint,
          log2HR,
          FDR,
          sig_star
        ),
      by = c(
        "cancer",
        "endpoint"
      )
    )

  lim <- max(
    abs(dd2$log2HR),
    na.rm = TRUE
  )

  if (!is.finite(lim) ||
      lim == 0) {
    lim <- 1
  }

  model_title <- unique(
    dat$model_label[
      dat$score == sc &
        dat$model == model_name
    ]
  )[1]

  p <- ggplot(
    dd2,
    aes(
      x = cancer,
      y = endpoint,
      fill = log2HR
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.35
    ) +
    geom_text(
      aes(label = sig_star),
      size = 4.0,
      color = "black",
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -lim, lim
      ),
      na.value = "grey92",
      name = expression(
        log[2](HR)
      )
    ) +
    labs(
      title = paste0(
        score_label(sc),
        " optimal-cutoff survival HR heatmap"
      ),
      subtitle = paste0(
        model_title,
        "\nHR = High vs Low; * FDR<0.05, ** FDR<0.01, *** FDR<0.001"
      ),
      x = "TCGA cancer type",
      y = NULL
    ) +
    theme_step(
      base_size = 9,
      title_size = 9.5,
      axis_title_size = 8.5,
      axis_text_size = 7.2
    ) +
    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
      panel.grid = element_blank()
    )

  base <- file.path(
    OUT_HEATMAP,
    paste0(
      safe_name(sc),
      "_OPTIMAL_CUTOFF_HR_heatmap_",
      safe_name(model_name),
      "_OS_DFS"
    )
  )

  save_gg_both(
    p,
    base,
    width = max(
      10,
      length(cancers) * 0.42
    ),
    height = 4.6
  )

  write.csv(
    dd2,
    paste0(
      base,
      "_data.csv"
    ),
    row.names = FALSE
  )
}

for (sc in unique(cox_all$score)) {
  for (md in unique(cox_all$model)) {
    plot_optimal_hr_heatmap(
      cox_all,
      sc,
      md
    )
  }
}

figready_all <- cox_all %>%
  select(
    cancer,
    endpoint,
    score,
    score_label,
    cutoff,
    cutoff_method,
    maxstat_statistic,
    minprop,
    comparison,
    model,
    model_label,
    covariates,
    n,
    n_events,
    n_low,
    n_high,
    logrank_p,
    HR,
    lower95,
    upper95,
    p,
    FDR,
    direction,
    raw_p_significance,
    FDR_significance,
    log2HR,
    neglog10P,
    neglog10FDR
  )

write.csv(
  figready_all,
  file.path(
    OUT_FIGREADY,
    "FigureReady_ALL_optimalCutoff_survival_statistics_v1_5.csv"
  ),
  row.names = FALSE
)

write.csv(
  figready_all %>%
    filter(
      is.finite(p),
      p < 0.05
    ),
  file.path(
    OUT_FIGREADY,
    "FigureReady_optimalCutoff_rawP_lt_0_05_v1_5.csv"
  ),
  row.names = FALSE
)

write.csv(
  figready_all %>%
    filter(
      is.finite(FDR),
      FDR < 0.05
    ),
  file.path(
    OUT_FIGREADY,
    "FigureReady_optimalCutoff_FDR_lt_0_05_v1_5.csv"
  ),
  row.names = FALSE
)

summary_counts <- figready_all %>%
  group_by(
    score,
    score_label,
    endpoint,
    model,
    model_label
  ) %>%
  summarise(
    n_cancers = n_distinct(cancer),
    n_rawP_lt_0_05 =
      sum(p < 0.05, na.rm = TRUE),
    n_FDR_lt_0_05 =
      sum(FDR < 0.05, na.rm = TRUE),
    n_risk_HR_gt_1 =
      sum(HR > 1, na.rm = TRUE),
    n_protective_HR_lt_1 =
      sum(HR < 1, na.rm = TRUE),
    median_HR =
      median(HR, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  summary_counts,
  file.path(
    OUT_FIGREADY,
    "FigureReady_optimalCutoff_survival_summary_counts_v1_5.csv"
  ),
  row.names = FALSE
)

audit <- cutoff_use %>%
  select(
    cancer,
    endpoint,
    score,
    cutoff,
    everything()
  ) %>%
  left_join(
    cox_all %>%
      select(
        cancer,
        endpoint,
        score,
        model,
        n,
        n_events,
        n_low,
        n_high,
        HR,
        p,
        FDR,
        covariates
      ),
    by = c(
      "cancer",
      "endpoint",
      "score"
    )
  )

write.csv(
  audit,
  file.path(
    OUT_TABLE,
    "ALL_optimalCutoff_model_availability_audit_v1_5.csv"
  ),
  row.names = FALSE
)

manifest <- data.frame(
  setting = c(
    "version",
    "result_dir",
    "merged_input",
    "cutoff_input",
    "output_dir",
    "group_definition",
    "univariate_model",
    "clinical_multivariable_model",
    "purity_multivariable_model",
    "volcano_color_logic",
    "bubble_outline_logic",
    "FDR_scope",
    "statistical_note"
  ),
  value = c(
    "v1.5",
    RESULT_DIR,
    INPUT_FILE,
    CUTOFF_FILE,
    OUT_ROOT,
    "Low = score <= v1.4 cutoff; High = score > v1.4 cutoff",
    "Surv ~ High_vs_Low",
    "Surv ~ High_vs_Low + largest viable subset of AGE/SEX/STAGE",
    "Surv ~ High_vs_Low + CPE + AGE, where available",
    "Red = raw Cox p<0.05; grey = p>=0.05; direction is x-axis/HR",
    "Black = Risk HR>1; Blue = Protective HR<1",
    "BH separately within each score x endpoint x model",
    paste(
      "Cut-offs were optimized on the same endpoint/cohort, so",
      "group-based Cox/log-rank p values are exploratory and optimistic."
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(
    OUT_TABLE,
    "optimalCutoff_downstream_survival_manifest_v1_5.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    paste0(
      "Run finished: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "RESULT_DIR: ",
      RESULT_DIR
    ),
    paste0(
      "Merged input: ",
      INPUT_FILE
    ),
    paste0(
      "Cut-off input: ",
      CUTOFF_FILE
    ),
    paste0(
      "Output: ",
      OUT_ROOT
    ),
    "",
    "Core definition:",
    "Low  = score <= endpoint-specific v1.4 optimal cut-off",
    "High = score >  endpoint-specific v1.4 optimal cut-off",
    "",
    "Main outputs:",
    "1) Cox tables with HR / CI / p / FDR",
    "2) survival volcano plots",
    "3) pan-cancer forest plots",
    "4) OS/DFS bubble survival summaries",
    "5) log2(HR) heatmaps",
    "",
    "Interpretation:",
    "Because the same endpoint was used to optimize each cut-off,",
    "the resulting grouped Cox/log-rank p values are selection-biased",
    "and should be treated as exploratory/sensitivity analyses.",
    "Use continuous Cox and median-split results as primary analyses."
  ),
  file.path(
    OUT_TABLE,
    "optimalCutoff_downstream_survival_run_note_v1_5.txt"
  )
)

msg("============================================================")
msg("DONE - optimal-cutoff downstream survival patch v1.5")
msg("Output: ", OUT_ROOT)
msg(
  "Main table: ",
  file.path(
    OUT_TABLE,
    "ALL_optimalCutoff_group_Cox_models_v1_5.csv"
  )
)
msg("============================================================")
