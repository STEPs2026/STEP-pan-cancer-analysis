# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(123)

RESULT_DIR <- NA_character_  

PER_CANCER_ENDPOINTS <- c("OS", "DFS")
OVERVIEW_ENDPOINTS   <- c("OS", "DFS")

MIN_SURV_N <- 20L

MINPROP <- 0.20

OVERVIEW_NCOL <- 4L
PNG_DPI <- 300

LOW_COLOR  <- "#006FB9"
HIGH_COLOR <- "#E60011"
GROUP_COLORS <- c("Low" = LOW_COLOR, "High" = HIGH_COLOR)

required_pkgs <- c(
  "survival", "survminer", "ggplot2", "dplyr",
  "stringr", "ggpubr", "gridExtra"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(stringr)
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

has_v12_merged <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x) &&
    file.exists(file.path(
      x,
      "09_pan_cancer_summary",
      "ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv"
    ))
}

resolve_result_dir <- function(manual = RESULT_DIR) {
  if (has_v12_merged(manual)) {
    return(normalizePath(manual, mustWork = TRUE))
  }

  script_dir <- get_script_dir()
  wd <- normalizePath(getwd(), mustWork = FALSE)

  direct_candidates <- unique(c(
    wd, script_dir, dirname(wd), dirname(script_dir)
  ))

  for (x in direct_candidates) {
    if (has_v12_merged(x)) {
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
      ok <- vapply(subs, has_v12_merged, logical(1))
      candidates <- c(candidates, subs[ok])
    }
  }

  candidates <- unique(candidates)

  if (length(candidates) == 0) {
    stop(
      "\nCannot auto-detect a completed functional-subset result folder.\n",
      "Please set RESULT_DIR manually at the top of this script.\n",
      "The folder must contain:\n",
      "09_pan_cancer_summary/ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv\n"
    )
  }

  mt <- file.info(candidates)$mtime
  chosen <- candidates[order(mt, decreasing = TRUE)][1]

  message("Auto-detected latest result folder: ", chosen)
  normalizePath(chosen, mustWork = TRUE)
}

RESULT_DIR <- resolve_result_dir(RESULT_DIR)

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

INPUT_FILE <- if (file.exists(INPUT_REPAIRED)) INPUT_REPAIRED else INPUT_V12

OUT_DIR <- file.path(RESULT_DIR, "08_survival_optimal_cutoff_v1_4")
OUT_OVERVIEW <- file.path(OUT_DIR, "00_ALL_CANCERS_OVERVIEW")
OUT_QC <- file.path(OUT_DIR, "00_cutoff_QC")
OUT_LOG <- file.path(RESULT_DIR, "00_logs")

for (d in c(OUT_DIR, OUT_OVERVIEW, OUT_QC, OUT_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(OUT_LOG, "KM_optimal_cutoff_patch_v1_4.log")

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
msg("OUT_DIR = ", OUT_DIR)
msg("Optimal cut-off method = survminer::surv_cutpoint / maxstat")
msg("MINPROP = ", MINPROP)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

format_p <- function(p, prefix = "log-rank p = ") {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) {
    return(paste0(prefix, "NA"))
  }
  if (p < 1e-4) {
    return(paste0(prefix, formatC(p, format = "e", digits = 2)))
  }
  paste0(prefix, signif(p, 3))
}

format_cutoff <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  format(signif(x, 4), trim = TRUE, scientific = FALSE)
}

choose_break_time <- function(max_time) {
  if (!is.finite(max_time) || max_time <= 0) return(NULL)
  if (max_time <= 12)  return(3)
  if (max_time <= 24)  return(6)
  if (max_time <= 60)  return(12)
  if (max_time <= 120) return(24)
  if (max_time <= 240) return(48)
  if (max_time <= 360) return(60)
  120
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

save_ggsurvplot_both <- function(g, file_base, width = 6.8, height = 6.6) {
  grDevices::pdf(
    paste0(file_base, ".pdf"),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(g, newpage = FALSE)
  grDevices::dev.off()

  grDevices::png(
    paste0(file_base, ".png"),
    width = round(width * PNG_DPI),
    height = round(height * PNG_DPI),
    res = PNG_DPI
  )
  print(g, newpage = FALSE)
  grDevices::dev.off()

  print(g)
  invisible(TRUE)
}

save_arranged_survplots_both <- function(g, file_base, width, height) {
  grDevices::pdf(
    paste0(file_base, ".pdf"),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(g)
  grDevices::dev.off()

  grDevices::png(
    paste0(file_base, ".png"),
    width = round(width * PNG_DPI),
    height = round(height * PNG_DPI),
    res = PNG_DPI
  )
  print(g)
  grDevices::dev.off()

  invisible(TRUE)
}

if (!file.exists(INPUT_FILE)) {
  stop("Cannot find merged score/clinical table: ", INPUT_FILE)
}

sdat <- read.csv(
  INPUT_FILE,
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

sdat$.cancer_opt <- toupper(trimws(as.character(sdat[[cancer_col]])))

SCORE_MAP <- c(
  "Hallmark_MYC_targets_v1" = "MYC targets v1 (7 genes)",
  "GO_extracellular_exsosome" = "Extracellular exosome (23 genes)"
)

SCORE_COLS <- intersect(names(SCORE_MAP), colnames(sdat))

if (length(SCORE_COLS) == 0) {
  stop(
    "Neither functional-subset score was found in the merged table:\n",
    paste(names(SCORE_MAP), collapse = "\n")
  )
}

score_label <- function(x) {
  lab <- unname(SCORE_MAP[x])
  ifelse(is.na(lab), x, lab)
}

endpoint_defs <- list(
  OS = c(
    time = find_col(
      c("OS_MONTHS", "OS_MONTH", "OS.time", "OS_TIME", "OS_DAYS"),
      colnames(sdat)
    ),
    event = find_col(
      c("OS_EVENT", "OS_STATUS", "OS", "OS.event"),
      colnames(sdat)
    )
  ),
  DFS = c(
    time = find_col(
      c("DFS_MONTHS", "DFS_MONTH", "DFS.time", "DFS_TIME", "DFS_DAYS"),
      colnames(sdat)
    ),
    event = find_col(
      c("DFS_EVENT", "DFS_STATUS", "DFS", "DFS.event"),
      colnames(sdat)
    )
  ),
  PFS = c(
    time = find_col(
      c("PFS_MONTHS", "PFS_MONTH", "PFS.time", "PFS_TIME", "PFS_DAYS"),
      colnames(sdat)
    ),
    event = find_col(
      c("PFS_EVENT", "PFS_STATUS", "PFS", "PFS.event"),
      colnames(sdat)
    )
  ),
  RFS = c(
    time = find_col(
      c("RFS_MONTHS", "RFS_MONTH", "RFS.time", "RFS_TIME", "RFS_DAYS"),
      colnames(sdat)
    ),
    event = find_col(
      c("RFS_EVENT", "RFS_STATUS", "RFS", "RFS.event"),
      colnames(sdat)
    )
  )
)

endpoint_availability <- data.frame(
  endpoint = names(endpoint_defs),
  time_col = vapply(endpoint_defs, function(x) x["time"], character(1)),
  event_col = vapply(endpoint_defs, function(x) x["event"], character(1)),
  available = vapply(
    endpoint_defs,
    function(x) all(!is.na(x)) && all(nzchar(x)),
    logical(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  endpoint_availability,
  file.path(OUT_QC, "detected_survival_endpoints_v1_4.csv"),
  row.names = FALSE
)

endpoint_defs <- endpoint_defs[endpoint_availability$available]

msg("Detected endpoints: ", paste(names(endpoint_defs), collapse = ", "))
msg("Scores: ", paste(SCORE_COLS, collapse = ", "))

find_optimal_cutoff <- function(d, score_col, time_col, event_col,
                                cancer, endpoint,
                                min_n = MIN_SURV_N,
                                minprop = MINPROP) {

  qc <- data.frame(
    cancer = cancer,
    endpoint = endpoint,
    score = score_col,
    score_label = score_label(score_col),
    method = "survminer::surv_cutpoint (maximally selected rank statistic)",
    minprop = minprop,
    n_cancer_total = nrow(d),
    n_time_valid = NA_integer_,
    n_event_valid = NA_integer_,
    n_score_valid = NA_integer_,
    n_complete_valid = NA_integer_,
    n_events = NA_integer_,
    n_censored = NA_integer_,
    n_unique_score = NA_integer_,
    cutoff = NA_real_,
    maxstat_statistic = NA_real_,
    n_low = NA_integer_,
    n_high = NA_integer_,
    prop_low = NA_real_,
    prop_high = NA_real_,
    logrank_p_after_optimization = NA_real_,
    eligible = FALSE,
    exclusion_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!all(c(score_col, time_col, event_col) %in% colnames(d))) {
    qc$exclusion_reason <- "missing_required_column"
    return(list(data = NULL, qc = qc))
  }

  tt <- suppressWarnings(as.numeric(d[[time_col]]))
  ee <- normalize_event_vector(d[[event_col]])
  ss <- suppressWarnings(as.numeric(d[[score_col]]))

  qc$n_time_valid <- sum(is.finite(tt) & tt > 0)
  qc$n_event_valid <- sum(ee %in% c(0, 1))
  qc$n_score_valid <- sum(is.finite(ss))

  ok <- is.finite(tt) & tt > 0 &
    ee %in% c(0, 1) &
    is.finite(ss)

  dd <- data.frame(
    time = tt[ok],
    event = ee[ok],
    score = ss[ok],
    stringsAsFactors = FALSE
  )

  qc$n_complete_valid <- nrow(dd)
  qc$n_events <- sum(dd$event == 1, na.rm = TRUE)
  qc$n_censored <- sum(dd$event == 0, na.rm = TRUE)
  qc$n_unique_score <- length(unique(dd$score))

  if (nrow(dd) < min_n) {
    qc$exclusion_reason <- paste0("n_complete_valid_lt_", min_n)
    return(list(data = NULL, qc = qc))
  }

  if (length(unique(dd$event)) < 2) {
    qc$exclusion_reason <- "only_one_event_status"
    return(list(data = NULL, qc = qc))
  }

  if (length(unique(dd$score)) < 3) {
    qc$exclusion_reason <- "score_has_lt_3_unique_values"
    return(list(data = NULL, qc = qc))
  }

  cp <- tryCatch(
    suppressWarnings(
      survminer::surv_cutpoint(
        data = dd,
        time = "time",
        event = "event",
        variables = "score",
        minprop = minprop,
        progressbar = FALSE
      )
    ),
    error = function(e) {
      attr(e, "v1_4_error") <- conditionMessage(e)
      e
    }
  )

  if (inherits(cp, "error")) {
    qc$exclusion_reason <- paste0(
      "surv_cutpoint_error: ",
      attr(cp, "v1_4_error")
    )
    return(list(data = NULL, qc = qc))
  }

  cutoff <- tryCatch(
    as.numeric(cp$cutpoint["score", "cutpoint"]),
    error = function(e) NA_real_
  )

  statistic <- tryCatch(
    as.numeric(cp$cutpoint["score", "statistic"]),
    error = function(e) NA_real_
  )

  if (!is.finite(cutoff)) {
    qc$exclusion_reason <- "optimal_cutoff_not_finite"
    return(list(data = NULL, qc = qc))
  }

  dd$group <- ifelse(dd$score <= cutoff, "Low", "High")
  dd$group <- factor(dd$group, levels = c("Low", "High"))

  n_low <- sum(dd$group == "Low")
  n_high <- sum(dd$group == "High")

  qc$cutoff <- cutoff
  qc$maxstat_statistic <- statistic
  qc$n_low <- n_low
  qc$n_high <- n_high
  qc$prop_low <- n_low / nrow(dd)
  qc$prop_high <- n_high / nrow(dd)

  if (n_low < 1 || n_high < 1) {
    qc$exclusion_reason <- "optimal_cutoff_created_empty_group"
    return(list(data = NULL, qc = qc))
  }

  min_expected <- floor(nrow(dd) * minprop)
  if (min(n_low, n_high) < max(1L, min_expected)) {
    qc$exclusion_reason <- "group_smaller_than_minprop_guard"
    return(list(data = NULL, qc = qc))
  }

  lr <- tryCatch(
    survival::survdiff(
      survival::Surv(time, event) ~ group,
      data = dd
    ),
    error = function(e) NULL
  )

  if (is.null(lr)) {
    qc$exclusion_reason <- "survdiff_error"
    return(list(data = NULL, qc = qc))
  }

  lr_p <- 1 - pchisq(
    lr$chisq,
    df = length(lr$n) - 1
  )

  qc$logrank_p_after_optimization <- lr_p
  qc$eligible <- TRUE
  qc$exclusion_reason <- "eligible"

  list(data = dd, qc = qc)
}

build_optimal_km <- function(d, score_col, cancer, endpoint,
                             overview_mode = FALSE) {

  if (!(endpoint %in% names(endpoint_defs))) {
    return(list(plot = NULL, qc = NULL))
  }

  time_col <- endpoint_defs[[endpoint]]["time"]
  event_col <- endpoint_defs[[endpoint]]["event"]

  obj <- find_optimal_cutoff(
    d = d,
    score_col = score_col,
    time_col = time_col,
    event_col = event_col,
    cancer = cancer,
    endpoint = endpoint
  )

  if (is.null(obj$qc) || !isTRUE(obj$qc$eligible[1])) {
    return(list(plot = NULL, qc = obj$qc))
  }

  dd <- obj$data
  qc <- obj$qc

  fit <- tryCatch(
    survival::survfit(
      survival::Surv(time, event) ~ group,
      data = dd
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    qc$eligible <- FALSE
    qc$exclusion_reason <- "survfit_error"
    return(list(plot = NULL, qc = qc))
  }

  bt <- choose_break_time(max(dd$time, na.rm = TRUE))

  if (overview_mode) {
    plot_title <- cancer
    legend_setting <- "none"
    title_size <- 8.5
    subtitle_size <- 6.6
    axis_title_size <- 7.2
    axis_text_size <- 6.8
    risk_font <- 2.7
    risk_height <- 0.30
  } else {
    plot_title <- paste0(
      cancer, " ", endpoint, " by ", score_label(score_col)
    )
    legend_setting <- "top"
    title_size <- 10
    subtitle_size <- 8
    axis_title_size <- 9
    axis_text_size <- 8
    risk_font <- 3.2
    risk_height <- 0.27
  }

  subtitle_txt <- paste0(
    "Optimal cut-off = ", format_cutoff(qc$cutoff),
    " | Low n=", qc$n_low,
    ", High n=", qc$n_high,
    " | maxstat, minprop=", MINPROP
  )

  g <- survminer::ggsurvplot(
    fit,
    data = dd,

    risk.table = "absolute",
    risk.table.height = risk_height,
    risk.table.title = "",
    risk.table.fontsize = risk_font,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,

    risk.table.col = "black",

    pval = format_p(
      qc$logrank_p_after_optimization,
      "log-rank p = "
    ),
    pval.size = ifelse(overview_mode, 2.6, 3.5),

    conf.int = FALSE,
    censor = TRUE,
    censor.shape = 124,
    censor.size = ifelse(overview_mode, 2.0, 3.0),

    palette = unname(GROUP_COLORS[c("Low", "High")]),
    legend.title = score_label(score_col),
    legend.labs = c("Low", "High"),
    legend = legend_setting,

    break.time.by = bt,
    xlab = "Time (months)",
    ylab = paste0(endpoint, " probability"),
    title = plot_title,

    ggtheme = ggplot2::theme_classic(
      base_size = ifelse(overview_mode, 8, 10)
    ),

    tables.theme = ggplot2::theme_classic(
      base_size = ifelse(overview_mode, 7, 8)
    )
  )

  g$plot <- g$plot +
    ggplot2::labs(subtitle = subtitle_txt) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_size,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = subtitle_size,
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(
        size = axis_text_size,
        color = "black"
      ),
      legend.position = if (overview_mode) "none" else "top",
      legend.text = ggplot2::element_text(
        size = ifelse(overview_mode, 7, 8)
      ),
      legend.title = ggplot2::element_text(
        size = ifelse(overview_mode, 7, 8)
      )
    )

  g$table <- g$table +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        size = ifelse(overview_mode, 6.8, 8)
      ),
      axis.text.x = ggplot2::element_text(
        size = ifelse(overview_mode, 6.2, 7.2),
        color = "black"
      ),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5)
    )

  list(plot = g, qc = qc)
}

cancers <- sort(unique(sdat$.cancer_opt))
cancers <- cancers[!is.na(cancers) & nzchar(cancers)]

requested_per_cancer <- intersect(
  PER_CANCER_ENDPOINTS,
  names(endpoint_defs)
)

missing_per_cancer <- setdiff(
  PER_CANCER_ENDPOINTS,
  names(endpoint_defs)
)

if (length(missing_per_cancer) > 0) {
  warning(
    "Requested endpoints not detected and therefore skipped: ",
    paste(missing_per_cancer, collapse = ", ")
  )
}

all_qc <- list()
qq <- 1L

for (ep in requested_per_cancer) {
  msg("Per-cancer optimal-cutoff KM: ", ep)

  for (ca in cancers) {
    dc <- sdat[sdat$.cancer_opt == ca, , drop = FALSE]

    ca_dir <- file.path(OUT_DIR, ca)
    dir.create(ca_dir, recursive = TRUE, showWarnings = FALSE)

    for (sc in SCORE_COLS) {
      obj <- tryCatch(
        build_optimal_km(
          d = dc,
          score_col = sc,
          cancer = ca,
          endpoint = ep,
          overview_mode = FALSE
        ),
        error = function(e) {
          data.frame_error <- data.frame(
            cancer = ca,
            endpoint = ep,
            score = sc,
            score_label = score_label(sc),
            method = "survminer::surv_cutpoint (maximally selected rank statistic)",
            minprop = MINPROP,
            n_cancer_total = nrow(dc),
            n_time_valid = NA_integer_,
            n_event_valid = NA_integer_,
            n_score_valid = NA_integer_,
            n_complete_valid = NA_integer_,
            n_events = NA_integer_,
            n_censored = NA_integer_,
            n_unique_score = NA_integer_,
            cutoff = NA_real_,
            maxstat_statistic = NA_real_,
            n_low = NA_integer_,
            n_high = NA_integer_,
            prop_low = NA_real_,
            prop_high = NA_real_,
            logrank_p_after_optimization = NA_real_,
            eligible = FALSE,
            exclusion_reason = paste0(
              "unexpected_error: ",
              conditionMessage(e)
            ),
            stringsAsFactors = FALSE
          )
          list(plot = NULL, qc = data.frame_error)
        }
      )

      if (!is.null(obj$qc)) {
        all_qc[[qq]] <- obj$qc
        qq <- qq + 1L
      }

      if (is.null(obj$plot)) next

      safe_sc <- safe_name(sc)

      file_base <- file.path(
        ca_dir,
        paste0(
          ca, "_", ep, "_", safe_sc,
          "_KM_OPTIMAL_CUTOFF"
        )
      )

      save_ggsurvplot_both(
        obj$plot,
        file_base,
        width = 6.8,
        height = 6.6
      )

      write.csv(
        obj$qc,
        paste0(file_base, "_cutoff_QC.csv"),
        row.names = FALSE
      )
    }
  }
}

qc_df <- if (length(all_qc) > 0) {
  dplyr::bind_rows(all_qc)
} else {
  data.frame()
}

if (nrow(qc_df) > 0) {
  qc_df <- qc_df %>%
    group_by(score, endpoint) %>%
    mutate(
      BH_FDR_across_cancers_after_optimization = ifelse(
        is.finite(logrank_p_after_optimization),
        p.adjust(
          logrank_p_after_optimization,
          method = "BH"
        ),
        NA_real_
      )
    ) %>%
    ungroup()

  write.csv(
    qc_df,
    file.path(
      OUT_QC,
      "ALL_cancers_OPTIMAL_cutoff_summary_and_eligibility_v1_4.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    qc_df %>% filter(eligible),
    file.path(
      OUT_QC,
      "ALL_cancers_OPTIMAL_cutoff_values_ELIGIBLE_only_v1_4.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    qc_df %>% filter(!eligible),
    file.path(
      OUT_QC,
      "ALL_cancers_OPTIMAL_cutoff_exclusions_v1_4.csv"
    ),
    row.names = FALSE
  )
}

requested_overview <- intersect(
  OVERVIEW_ENDPOINTS,
  names(endpoint_defs)
)

overview_qc <- list()
oq <- 1L

for (ep in requested_overview) {
  for (sc in SCORE_COLS) {

    g_list <- list()
    q_list <- list()
    gi <- 1L

    for (ca in cancers) {
      dc <- sdat[sdat$.cancer_opt == ca, , drop = FALSE]

      obj <- tryCatch(
        build_optimal_km(
          d = dc,
          score_col = sc,
          cancer = ca,
          endpoint = ep,
          overview_mode = TRUE
        ),
        error = function(e) NULL
      )

      if (is.null(obj) || is.null(obj$qc)) next

      overview_qc[[oq]] <- obj$qc
      oq <- oq + 1L

      if (is.null(obj$plot)) next

      g_list[[gi]] <- obj$plot
      q_list[[gi]] <- obj$qc
      gi <- gi + 1L
    }

    if (length(g_list) < 3) {
      msg(
        "Skip overview ", ep, " / ", score_label(sc),
        ": fewer than 3 eligible cancers."
      )
      next
    }

    nca <- length(g_list)
    ncol_use <- max(1L, as.integer(OVERVIEW_NCOL))
    nrow_use <- ceiling(nca / ncol_use)

    main_title <- paste0(
      "All cancers: ", ep, " by ", score_label(sc),
      "\nCancer-specific optimal cut-off (maxstat); ",
      "Low = blue (#006FB9), High = red (#E60011), ",
      "minprop = ", MINPROP
    )

    arranged <- survminer::arrange_ggsurvplots(
      x = g_list,
      print = FALSE,
      title = main_title,
      ncol = ncol_use,
      nrow = nrow_use,
      risk.table.height = 0.30,
      surv.plot.height = 0.70
    )

    safe_sc <- safe_name(sc)

    base_out <- file.path(
      OUT_OVERVIEW,
      paste0(
        "ALL_cancers_", ep, "_", safe_sc,
        "_KM_OPTIMAL_CUTOFF_faceted_with_risktable"
      )
    )

    out_width <- max(12, ncol_use * 3.7)
    out_height <- max(6.5, nrow_use * 3.45 + 0.9)

    save_arranged_survplots_both(
      arranged,
      base_out,
      width = out_width,
      height = out_height
    )

    if (length(q_list) > 0) {
      write.csv(
        dplyr::bind_rows(q_list),
        paste0(base_out, "_included_cancers_and_cutoffs.csv"),
        row.names = FALSE
      )
    }

    msg(
      "Finished optimal-cutoff overview: ",
      ep, " / ", score_label(sc),
      " | eligible cancers = ", nca
    )
  }
}

if (length(overview_qc) > 0) {
  write.csv(
    dplyr::bind_rows(overview_qc),
    file.path(
      OUT_OVERVIEW,
      "ALL_cancers_OPTIMAL_cutoff_overview_QC_v1_4.csv"
    ),
    row.names = FALSE
  )
}

if (nrow(qc_df) > 0) {
  figure_ready <- qc_df %>%
    filter(eligible) %>%
    transmute(
      cancer,
      endpoint,
      score,
      score_label,
      cutoff,
      maxstat_statistic,
      n = n_complete_valid,
      n_events,
      n_low,
      n_high,
      prop_low,
      prop_high,
      logrank_p = logrank_p_after_optimization,
      BH_FDR_across_cancers = BH_FDR_across_cancers_after_optimization,
      method,
      minprop
    )

  write.csv(
    figure_ready,
    file.path(
      OUT_QC,
      "FigureReady_OPTIMAL_cutoff_KM_summary_v1_4.csv"
    ),
    row.names = FALSE
  )
}

manifest <- data.frame(
  setting = c(
    "version",
    "result_dir",
    "input_file",
    "output_dir",
    "score_1",
    "score_2",
    "per_cancer_endpoints",
    "overview_endpoints",
    "cutoff_method",
    "minprop",
    "minimum_survival_n",
    "low_color",
    "high_color",
    "cutoff_scope",
    "primary_analysis_note"
  ),
  value = c(
    "v1.4",
    RESULT_DIR,
    INPUT_FILE,
    OUT_DIR,
    "Hallmark_MYC_targets_v1",
    "GO_extracellular_exsosome",
    paste(PER_CANCER_ENDPOINTS, collapse = ";"),
    paste(OVERVIEW_ENDPOINTS, collapse = ";"),
    "survminer::surv_cutpoint / maximally selected rank statistic",
    as.character(MINPROP),
    as.character(MIN_SURV_N),
    LOW_COLOR,
    HIGH_COLOR,
    "cancer x score x endpoint",
    paste(
      "Optimal-cutoff KM is exploratory because cut-off selection and",
      "log-rank testing use the same cohort; retain median-split and",
      "continuous Cox analyses as primary/confirmatory."
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(OUT_QC, "KM_optimal_cutoff_analysis_manifest_v1_4.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste0("Run finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Result directory: ", RESULT_DIR),
    paste0("Input merged table: ", INPUT_FILE),
    paste0("Output: ", OUT_DIR),
    paste0("Method: survminer::surv_cutpoint / maxstat"),
    paste0("MINPROP: ", MINPROP),
    paste0("Per-cancer endpoints: ", paste(PER_CANCER_ENDPOINTS, collapse = ", ")),
    paste0("Low color: ", LOW_COLOR),
    paste0("High color: ", HIGH_COLOR),
    "",
    "Interpretation:",
    "The cut-off is optimized separately for each cancer, score, and endpoint.",
    "Thus OS and DFS may have different cut-offs in the same cancer.",
    "The displayed log-rank p value is calculated after optimization and is",
    "therefore exploratory/optimistic. Median-split KM and continuous Cox",
    "should remain the primary confirmatory analyses."
  ),
  file.path(OUT_QC, "KM_optimal_cutoff_run_note_v1_4.txt")
)

msg("============================================================")
msg("DONE - KM optimal-cutoff patch v1.4")
msg("Output folder: ", OUT_DIR)
msg("Main cutoff table: ",
    file.path(OUT_QC, "ALL_cancers_OPTIMAL_cutoff_summary_and_eligibility_v1_4.csv"))
msg("============================================================")
