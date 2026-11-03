# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
rm(list = ls())
options(stringsAsFactors = FALSE)

ROOT_DIR <- getwd()

INPUT_FILE <- file.path(
  ROOT_DIR,
  "PRECOG_Gide_CTLA4PD1_22rows_to_original_PFS_OS_mapping.csv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "PRECOG_Gide_CTLA4PD1_survival_rescue_KM_22rows_SENSITIVITY"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

COLOR_LOW  <- "#4575B4"
COLOR_HIGH <- "#D73027"
DPI <- 300

pkgs <- c(
  "survival",
  "survminer",
  "ggplot2",
  "dplyr"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

if (!file.exists(INPUT_FILE)) {
  stop(
    "Cannot find input file:\n",
    INPUT_FILE
  )
}

d <- read.csv(
  INPUT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "Expression_ID",
  "Gide_original_patient_no",
  "STEP_UP_resid_CAF3_z",
  "STEP_UP_resid_CAF3_z_group",
  "Progression Free Survival (Days)",
  "PFS_event",
  "Overall Survival (Days)",
  "OS_event"
)

miss <- setdiff(
  required_cols,
  names(d)
)

if (length(miss) > 0) {
  stop(
    "Missing columns: ",
    paste(
      miss,
      collapse = ", "
    )
  )
}

num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

d$Score <- num(
  d$STEP_UP_resid_CAF3_z
)

d$PFS_Days <- num(
  d[[
    "Progression Free Survival (Days)"
  ]]
)

d$OS_Days <- num(
  d[[
    "Overall Survival (Days)"
  ]]
)

d$PFS_event <- num(
  d$PFS_event
)

d$OS_event <- num(
  d$OS_event
)

d$STEP_group <- factor(
  as.character(
    d$STEP_UP_resid_CAF3_z_group
  ),
  levels = c(
    "Low",
    "High"
  )
)

if (
  nrow(d) != 22 ||
  sum(d$STEP_group == "Low") != 11 ||
  sum(d$STEP_group == "High") != 11
) {
  stop(
    "Expected exactly 22 PRECOG rows with 11 Low and 11 High. ",
    "Observed N=", nrow(d),
    ", Low=", sum(d$STEP_group == "Low"),
    ", High=", sum(d$STEP_group == "High")
  )
}

if (
  any(!is.finite(d$PFS_Days)) ||
  any(!is.finite(d$OS_Days)) ||
  any(!d$PFS_event %in% c(0, 1)) ||
  any(!d$OS_event %in% c(0, 1))
) {
  stop(
    "PFS/OS time or event contains invalid values."
  )
}

dup_qc <- d %>%
  group_by(
    Gide_original_patient_no
  ) %>%
  summarise(
    Row_N = n(),
    Groups = paste(
      sort(
        unique(
          as.character(
            STEP_group
          )
        )
      ),
      collapse = ";"
    ),
    Expression_IDs = paste(
      Expression_ID,
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  filter(
    Row_N > 1
  )

write.csv(
  dup_qc,
  file.path(
    OUT_DIR,
    "06_duplicate_patient_QC.csv"
  ),
  row.names = FALSE
)

write.csv(
  d,
  file.path(
    OUT_DIR,
    "05_22rows_survival_input.csv"
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("22-ROW PRECOG SAMPLE-LEVEL SENSITIVITY ANALYSIS\n")
cat("============================================================\n")
cat("Expression rows: ", nrow(d), "\n", sep = "")
cat("Unique Gide patients: ", n_distinct(d$Gide_original_patient_no), "\n", sep = "")
cat("\nPRECOG original group counts:\n")
print(table(d$STEP_group))
cat("\nDuplicated patients:\n")
print(dup_qc)

p_text <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fit_group_cox <- function(
    time,
    event,
    group
) {

  dd <- data.frame(
    .time = time,
    .event = event,
    .group = group
  )

  fit <- tryCatch(
    coxph(
      Surv(
        .time,
        .event
      ) ~
        .group,
      data = dd
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(
      data.frame(
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        P = NA_real_
      )
    )
  }

  s <- summary(fit)

  data.frame(
    HR = unname(
      s$conf.int[
        1,
        "exp(coef)"
      ]
    ),
    CI_low = unname(
      s$conf.int[
        1,
        "lower .95"
      ]
    ),
    CI_high = unname(
      s$conf.int[
        1,
        "upper .95"
      ]
    ),
    P = unname(
      s$coefficients[
        1,
        "Pr(>|z|)"
      ]
    )
  )
}

fit_continuous_cox <- function(
    time,
    event,
    score
) {

  dd <- data.frame(
    .time = time,
    .event = event,
    .score_z = as.numeric(
      scale(
        score
      )
    )
  )

  fit <- tryCatch(
    coxph(
      Surv(
        .time,
        .event
      ) ~
        .score_z,
      data = dd
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(
      data.frame(
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        P = NA_real_
      )
    )
  }

  s <- summary(fit)

  data.frame(
    HR = unname(
      s$conf.int[
        1,
        "exp(coef)"
      ]
    ),
    CI_low = unname(
      s$conf.int[
        1,
        "lower .95"
      ]
    ),
    CI_high = unname(
      s$conf.int[
        1,
        "upper .95"
      ]
    ),
    P = unname(
      s$coefficients[
        1,
        "Pr(>|z|)"
      ]
    )
  )
}

make_km_22 <- function(
    data,
    endpoint = c(
      "PFS",
      "OS"
    )
) {

  endpoint <- match.arg(
    endpoint
  )

  time_col <- if (
    endpoint == "PFS"
  ) {
    "PFS_Days"
  } else {
    "OS_Days"
  }

  event_col <- if (
    endpoint == "PFS"
  ) {
    "PFS_event"
  } else {
    "OS_event"
  }

  km_data <- data.frame(
    .time = num(
      data[[
        time_col
      ]]
    ),
    .event = num(
      data[[
        event_col
      ]]
    ),
    .group = factor(
      data$STEP_group,
      levels = c(
        "Low",
        "High"
      )
    ),
    .score = num(
      data$Score
    ),
    .patient = data$Gide_original_patient_no,
    stringsAsFactors = FALSE
  )

  km_data <- km_data[
    complete.cases(
      km_data
    ) &
      km_data$.time > 0 &
      km_data$.event %in%
        c(
          0,
          1
        ),
    ,
    drop = FALSE
  ]

  fit <- survival::survfit(
    survival::Surv(
      .time,
      .event
    ) ~
      .group,
    data = km_data
  )

  lr <- survival::survdiff(
    survival::Surv(
      .time,
      .event
    ) ~
      .group,
    data = km_data
  )

  logrank_p <- 1 - pchisq(
    lr$chisq,
    df = length(
      lr$n
    ) - 1
  )

  gc <- fit_group_cox(
    km_data$.time,
    km_data$.event,
    km_data$.group
  )

  cc <- fit_continuous_cox(
    km_data$.time,
    km_data$.event,
    km_data$.score
  )

  subtitle <- paste0(
    "NAIVE sample-level analysis: High vs Low HR=",
    ifelse(
      is.finite(gc$HR),
      sprintf(
        "%.2f [%.2f–%.2f]",
        gc$HR,
        gc$CI_low,
        gc$CI_high
      ),
      "NA"
    ),
    "; Cox p=",
    p_text(
      gc$P
    ),
    "\nContinuous HR/1SD=",
    ifelse(
      is.finite(cc$HR),
      sprintf(
        "%.2f [%.2f–%.2f]",
        cc$HR,
        cc$CI_low,
        cc$CI_high
      ),
      "NA"
    ),
    "; p=",
    p_text(
      cc$P
    ),
    "; rows=22; unique patients=16"
  )

  g <- survminer::ggsurvplot(
    fit,
    data = km_data,
    risk.table = TRUE,
    risk.table.height = 0.25,
    pval = paste0(
      "Naive log-rank p = ",
      p_text(
        logrank_p
      )
    ),
    pval.size = 4,
    conf.int = FALSE,
    censor = TRUE,
    palette = c(
      COLOR_LOW,
      COLOR_HIGH
    ),
    legend.title = "PRECOG\nSTEP group",
    legend.labs = c(
      "Low (n=11 rows)",
      "High (n=11 rows)"
    ),
    xlab = paste0(
      endpoint,
      " time (days)"
    ),
    ylab = paste0(
      endpoint,
      " probability"
    ),
    title = paste0(
      "Gide CTLA4 + PD-1: ",
      endpoint,
      " | PRECOG 22-row sensitivity"
    ),
    ggtheme = theme_classic(
      base_size = 12
    ),
    tables.theme = theme_cleantable(
      base_size = 10
    )
  )

  g$plot <- g$plot +
    labs(
      subtitle = subtitle,
      caption = paste0(
        "Descriptive sensitivity only: 22 expression rows represent 16 unique patients; ",
        "patients 33/35/36 contribute samples to both High and Low groups."
      )
    ) +
    theme(
      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = .5
      ),
      plot.subtitle = element_text(
        size = 9.2,
        hjust = .5
      ),
      plot.caption = element_text(
        size = 8,
        hjust = 0
      )
    )

  stats <- data.frame(
    Endpoint = endpoint,
    Expression_rows = nrow(
      km_data
    ),
    Unique_patients = n_distinct(
      km_data$.patient
    ),
    Low_rows = sum(
      km_data$.group == "Low"
    ),
    High_rows = sum(
      km_data$.group == "High"
    ),
    Events = sum(
      km_data$.event == 1
    ),
    Low_events = sum(
      km_data$.event[
        km_data$.group == "Low"
      ] == 1
    ),
    High_events = sum(
      km_data$.event[
        km_data$.group == "High"
      ] == 1
    ),
    Naive_Logrank_P = logrank_p,
    Naive_High_vs_Low_HR = gc$HR,
    Naive_CI_low = gc$CI_low,
    Naive_CI_high = gc$CI_high,
    Naive_Cox_P = gc$P,
    Naive_Continuous_HR_per1SD = cc$HR,
    Naive_Continuous_CI_low = cc$CI_low,
    Naive_Continuous_CI_high = cc$CI_high,
    Naive_Continuous_Cox_P = cc$P,
    Interpretation = "SAMPLE_LEVEL_SENSITIVITY_ONLY_PSEUDOREPLICATION",
    stringsAsFactors = FALSE
  )

  list(
    plot = g,
    stats = stats
  )
}

pfs <- make_km_22(
  d,
  "PFS"
)

os <- make_km_22(
  d,
  "OS"
)

stats <- bind_rows(
  pfs$stats,
  os$stats
)

write.csv(
  stats,
  file.path(
    OUT_DIR,
    "04_22rows_naive_survival_statistics.csv"
  ),
  row.names = FALSE
)

save_surv <- function(
    x,
    filename,
    width = 7.2,
    height = 7.2
) {

  pdf(
    file.path(
      OUT_DIR,
      paste0(
        filename,
        ".pdf"
      )
    ),
    width = width,
    height = height
  )
  print(x)
  dev.off()

  png(
    file.path(
      OUT_DIR,
      paste0(
        filename,
        ".png"
      )
    ),
    width = width,
    height = height,
    units = "in",
    res = DPI
  )
  print(x)
  dev.off()
}

save_surv(
  pfs$plot,
  "01_PFS_KM_22rows"
)

save_surv(
  os$plot,
  "02_OS_KM_22rows"
)

combined <- survminer::arrange_ggsurvplots(
  list(
    pfs$plot,
    os$plot
  ),
  print = FALSE,
  ncol = 2,
  nrow = 1,
  risk.table = TRUE
)

pdf(
  file.path(
    OUT_DIR,
    "03_PFS_OS_KM_22rows_combined.pdf"
  ),
  width = 14.5,
  height = 7.8
)
print(combined)
dev.off()

png(
  file.path(
    OUT_DIR,
    "03_PFS_OS_KM_22rows_combined.png"
  ),
  width = 14.5,
  height = 7.8,
  units = "in",
  res = DPI
)
print(combined)
dev.off()

print(
  pfs$plot
)

print(
  os$plot
)

print(
  combined
)

cat("\n============================================================\n")
cat("22-ROW SENSITIVITY KM COMPLETE\n")
cat("============================================================\n")
print(stats)

cat(
  "\nWARNING:\n",
  "These 22-row p-values/HRs are descriptive only because rows are not independent patients.\n",
  "Use the 16-unique-patient KM as the primary survival analysis.\n",
  sep = ""
)
