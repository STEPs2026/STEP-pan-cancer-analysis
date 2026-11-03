# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(20260713)

ROOT_DIR <- "path/to/ICI_cohorts"

IMVIGOR_DIR <- file.path(
  ROOT_DIR,
  "cBioPortal",
  "blca_iatlas_imvigor210_2017"
)

GSE78220_DIR <- file.path(
  ROOT_DIR,
  "GSE78220"
)

GENESET_DIR <- file.path(
  ROOT_DIR,
  "gene set"
)

GMT_FILE <- file.path(
  GENESET_DIR,
  "STEP_all_human_gene_sets_merged.gmt"
)

IMVIGOR_EXPRESSION_FILE <- file.path(
  IMVIGOR_DIR,
  "data_mrna_seq_tpm.txt"
)

IMVIGOR_SAMPLE_CLINICAL_FILE <- file.path(
  IMVIGOR_DIR,
  "data_clinical_sample.txt"
)

IMVIGOR_PATIENT_CLINICAL_FILE <- file.path(
  IMVIGOR_DIR,
  "data_clinical_patient.txt"
)

GSE78220_EXPRESSION_FILE <- file.path(
  GSE78220_DIR,
  "GSE78220_PatientFPKM.xlsx"
)

GSE78220_METADATA_FILE <- file.path(
  GSE78220_DIR,
  "GSE78220_series_matrix.txt.gz"
)

METHODS_TO_RUN <- c("ZMEAN", "ssGSEA")

MAIN_METHOD <- "ZMEAN"
MAIN_SCORE  <- "STEP_CAF3_adjust"

SCORES_TO_RUN <- c(
  "STEP_UP",
  "CAF3",
  "STEP_CAF3_adjust"
)

GSE78220_BASELINE_ONLY <- TRUE

AVERAGE_REPEATED_BASELINE_BIOPSIES <- TRUE

RUN_CLINICAL_BENEFIT_SENSITIVITY <- TRUE

OUTPUT_FOLDER <- "IMvigor210_GSE78220_V5_6_EXPLICIT_PATHS_results"

cran_packages <- c(
  "data.table",
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "scales",
  "survival",
  "survminer",
  "pROC",
  "ggpubr",
  "cowplot",
  "gridExtra"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {
  stop(
    "Install missing packages first:\ninstall.packages(c(",
    paste(
      sprintf('"%s"', missing_cran),
      collapse = ", "
    ),
    "))"
  )
}

if (
  "ssGSEA" %in% METHODS_TO_RUN &&
  !requireNamespace("GSVA", quietly = TRUE)
) {
  stop(
    "GSVA is required for ssGSEA.\n",
    "Run:\n",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) ",
    "install.packages('BiocManager')\n",
    "BiocManager::install('GSVA')"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(survival)
  library(survminer)
  library(pROC)
  library(ggpubr)
  library(cowplot)
  library(gridExtra)
})

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_sample_id <- function(x) {
  y <- trimws(
    gsub(
      '^"|"$',
      "",
      as.character(x)
    )
  )

  y <- ifelse(
    grepl("^X[0-9]+$", y),
    sub("^X", "", y),
    y
  )

  y <- sub(
    "\\.0$",
    "",
    y
  )

  y
}

first_non_missing <- function(x) {
  y <- x[
    !is.na(x) &
      trimws(as.character(x)) != ""
  ]

  if (length(y) == 0) {
    return(NA)
  }

  y[1]
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p)) {
    return("P = NA")
  }

  if (p < 0.001) {
    return("P < 0.001")
  }

  paste0(
    "P = ",
    format.pval(
      p,
      digits = 3,
      eps = 0.001
    )
  )
}

sanitize_filename <- function(x) {
  gsub(
    "[^A-Za-z0-9_.-]+",
    "_",
    x
  )
}

find_files_recursive <- function(root, pattern) {
  hits <- list.files(
    root,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  unique(
    hits[file.exists(hits)]
  )
}

find_one_file <- function(
  root,
  patterns,
  required = TRUE
) {
  for (pattern in patterns) {
    hits <- find_files_recursive(
      root,
      pattern
    )

    hits <- hits[
      !grepl(
        OUTPUT_FOLDER,
        hits,
        fixed = TRUE
      )
    ]

    if (length(hits) > 0) {
      return(
        normalizePath(
          hits[1],
          mustWork = TRUE
        )
      )
    }
  }

  if (required) {
    stop(
      "Cannot find a file matching:\n",
      paste(
        patterns,
        collapse = "\n"
      )
    )
  }

  NA_character_
}

maybe_log2 <- function(expression_matrix) {
  total_values <- length(expression_matrix)

  if (total_values > 1000000) {
    sampled_index <- sample.int(
      total_values,
      size = 1000000
    )
    values <- as.numeric(
      expression_matrix[sampled_index]
    )
  } else {
    values <- as.numeric(
      expression_matrix
    )
  }

  values <- values[
    is.finite(values) &
      values >= 0
  ]

  q99 <- unname(
    quantile(
      values,
      0.99,
      na.rm = TRUE
    )
  )

  maximum <- max(
    values,
    na.rm = TRUE
  )

  if (
    is.finite(q99) &&
    (
      q99 > 50 ||
        maximum > 100
    )
  ) {
    return(
      list(
        matrix = log2(
          expression_matrix + 1
        ),
        transformation = "log2(x + 1)"
      )
    )
  }

  list(
    matrix = expression_matrix,
    transformation = "No additional log transformation"
  )
}

save_ggplot <- function(
  plot_object,
  directory,
  stem,
  width = 7,
  height = 6
) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  print(plot_object)

  pdf(
    file.path(
      directory,
      paste0(stem, ".pdf")
    ),
    width = width,
    height = height,
    onefile = TRUE
  )

  print(plot_object)
  dev.off()

  if (
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ) {
    ragg::agg_png(
      file.path(
        directory,
        paste0(stem, ".png")
      ),
      width = width,
      height = height,
      units = "in",
      res = 300
    )
  } else {
    png(
      file.path(
        directory,
        paste0(stem, ".png")
      ),
      width = width,
      height = height,
      units = "in",
      res = 300
    )
  }

  print(plot_object)
  dev.off()
}

save_grob <- function(
  grob_object,
  directory,
  stem,
  width = 16,
  height = 12
) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  grid::grid.newpage()
  grid::grid.draw(grob_object)

  pdf(
    file.path(
      directory,
      paste0(stem, ".pdf")
    ),
    width = width,
    height = height,
    onefile = TRUE
  )

  grid::grid.newpage()
  grid::grid.draw(grob_object)
  dev.off()

  if (
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ) {
    ragg::agg_png(
      file.path(
        directory,
        paste0(stem, ".png")
      ),
      width = width,
      height = height,
      units = "in",
      res = 300
    )
  } else {
    png(
      file.path(
        directory,
        paste0(stem, ".png")
      ),
      width = width,
      height = height,
      units = "in",
      res = 300
    )
  }

  grid::grid.newpage()
  grid::grid.draw(grob_object)
  dev.off()
}

read_gmt_simple <- function(path) {
  lines <- readLines(
    path,
    warn = FALSE
  )

  split_lines <- strsplit(
    lines,
    "\t",
    fixed = TRUE
  )

  gene_sets <- lapply(
    split_lines,
    function(x) {
      unique(
        toupper(
          x[-c(1, 2)]
        )
      )
    }
  )

  names(gene_sets) <- vapply(
    split_lines,
    `[`,
    character(1),
    1
  )

  gene_sets
}

select_gene_set <- function(
  gene_sets,
  exact_name,
  regex
) {
  if (
    exact_name %in% names(gene_sets)
  ) {
    return(
      gene_sets[[exact_name]]
    )
  }

  hit <- grep(
    regex,
    names(gene_sets),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(hit) == 0) {
    stop(
      "Cannot find gene set matching: ",
      regex
    )
  }

  message(
    "Using gene set: ",
    hit[1]
  )

  gene_sets[[hit[1]]]
}

score_zmean <- function(
  expression_matrix,
  scoring_sets
) {
  gene_z <- t(
    scale(
      t(expression_matrix)
    )
  )

  gene_z[
    !is.finite(gene_z)
  ] <- NA_real_

  output <- rbind(
    STEP_UP = colMeans(
      gene_z[
        scoring_sets$STEP_UP,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    ),
    CAF3 = colMeans(
      gene_z[
        scoring_sets$CAF3,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  )

  output
}

score_ssgsea <- function(
  expression_matrix,
  scoring_sets
) {
  if (
    exists(
      "ssgseaParam",
      where = asNamespace("GSVA"),
      inherits = FALSE
    )
  ) {
    parameter <- GSVA::ssgseaParam(
      exprData = expression_matrix,
      geneSets = scoring_sets,
      minSize = 1,
      maxSize = Inf,
      normalize = TRUE
    )

    return(
      GSVA::gsva(
        parameter,
        verbose = FALSE
      )
    )
  }

  GSVA::gsva(
    expr = expression_matrix,
    gset.idx.list = scoring_sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE,
    min.sz = 1,
    max.sz = Inf,
    parallel.sz = 1,
    verbose = FALSE
  )
}

calculate_scores <- function(
  expression_matrix,
  method,
  step_genes,
  caf3_genes
) {
  step_present <- intersect(
    step_genes,
    rownames(expression_matrix)
  )

  caf3_present <- intersect(
    caf3_genes,
    rownames(expression_matrix)
  )

  if (length(step_present) < 5) {
    stop(
      "Fewer than five STEP_UP genes were detected."
    )
  }

  if (length(caf3_present) < 2) {
    stop(
      "Fewer than two CAF3 genes were detected."
    )
  }

  scoring_sets <- list(
    STEP_UP = step_present,
    CAF3 = caf3_present
  )

  score_matrix <- switch(
    method,
    ZMEAN = score_zmean(
      expression_matrix,
      scoring_sets
    ),
    ssGSEA = score_ssgsea(
      expression_matrix,
      scoring_sets
    ),
    stop(
      "Unsupported method: ",
      method
    )
  )

  score_table <- data.frame(
    SAMPLE_ID = colnames(
      expression_matrix
    ),
    STEP_UP = as.numeric(
      score_matrix[
        "STEP_UP",
        colnames(expression_matrix)
      ]
    ),
    CAF3 = as.numeric(
      score_matrix[
        "CAF3",
        colnames(expression_matrix)
      ]
    ),
    stringsAsFactors = FALSE
  )

  adjustment_fit <- lm(
    STEP_UP ~ CAF3,
    data = score_table
  )

  score_table$STEP_CAF3_adjust <- residuals(
    adjustment_fit
  )

  coverage <- data.frame(
    GENE_SET = c(
      "STEP_UP",
      "CAF3"
    ),
    N_IN_GMT = c(
      length(step_genes),
      length(caf3_genes)
    ),
    N_PRESENT = c(
      length(step_present),
      length(caf3_present)
    ),
    COVERAGE = c(
      length(step_present) /
        length(step_genes),
      length(caf3_present) /
        length(caf3_genes)
    ),
    PRESENT_GENES = c(
      paste(
        step_present,
        collapse = ";"
      ),
      paste(
        caf3_present,
        collapse = ";"
      )
    ),
    MISSING_GENES = c(
      paste(
        setdiff(
          step_genes,
          step_present
        ),
        collapse = ";"
      ),
      paste(
        setdiff(
          caf3_genes,
          caf3_present
        ),
        collapse = ";"
      )
    ),
    stringsAsFactors = FALSE
  )

  adjustment <- data.frame(
    ALPHA_INTERCEPT = unname(
      coef(adjustment_fit)[1]
    ),
    BETA_CAF3 = unname(
      coef(adjustment_fit)[2]
    ),
    BETA_P = summary(
      adjustment_fit
    )$coefficients[
      "CAF3",
      "Pr(>|t|)"
    ],
    R_SQUARED = summary(
      adjustment_fit
    )$r.squared,
    stringsAsFactors = FALSE
  )

  list(
    scores = score_table,
    coverage = coverage,
    adjustment = adjustment
  )
}

parse_series_matrix <- function(path) {
  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  lines <- readLines(
    con,
    warn = FALSE
  )

  read_values <- function(prefix) {
    index <- which(
      startsWith(
        lines,
        prefix
      )
    )

    if (length(index) == 0) {
      return(character())
    }

    values <- strsplit(
      lines[index[1]],
      "\t",
      fixed = TRUE
    )[[1]][-1]

    gsub(
      '^"|"$',
      "",
      values
    )
  }

  titles <- read_values(
    "!Sample_title"
  )

  gsm <- read_values(
    "!Sample_geo_accession"
  )

  metadata <- data.frame(
    SAMPLE_TITLE = titles,
    GSM = gsm,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  characteristic_lines <- lines[
    startsWith(
      lines,
      "!Sample_characteristics_ch1"
    )
  ]

  used_names <- names(metadata)

  for (line in characteristic_lines) {
    values <- strsplit(
      line,
      "\t",
      fixed = TRUE
    )[[1]][-1]

    values <- gsub(
      '^"|"$',
      "",
      values
    )

    if (
      length(values) !=
        nrow(metadata)
    ) {
      next
    }

    raw_key <- trimws(
      sub(
        ":.*$",
        "",
        values[1]
      )
    )

    key <- toupper(
      gsub(
        "[^A-Za-z0-9]+",
        "_",
        raw_key
      )
    )

    key <- gsub(
      "^_+|_+$",
      "",
      key
    )

    key <- make.unique(
      c(
        used_names,
        key
      )
    )[
      length(used_names) + 1
    ]

    used_names <- c(
      used_names,
      key
    )

    metadata[[key]] <- trimws(
      sub(
        "^[^:]*:\\s*",
        "",
        values
      )
    )
  }

  metadata
}

load_imvigor_expression <- function(path) {
  expression_dt <- data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    data.table = TRUE,
    check.names = FALSE,
    showProgress = TRUE
  )

  if (!"Hugo_Symbol" %in% names(expression_dt)) {
    stop(
      "IMvigor210 expression file must contain a Hugo_Symbol column. ",
      "Detected columns begin with: ",
      paste(
        head(
          names(expression_dt),
          10
        ),
        collapse = ", "
      )
    )
  }

  expression_dt[
    ,
    Hugo_Symbol := toupper(
      trimws(
        as.character(
          Hugo_Symbol
        )
      )
    )
  ]

  sample_columns <- setdiff(
    names(expression_dt),
    "Hugo_Symbol"
  )

  original_sample_ids <- sample_columns
  normalized_sample_ids <- normalize_sample_id(
    original_sample_ids
  )

  if (
    anyDuplicated(
      normalized_sample_ids
    )
  ) {
    stop(
      "Duplicate IMvigor210 sample IDs were created after normalization: ",
      paste(
        unique(
          normalized_sample_ids[
            duplicated(
              normalized_sample_ids
            )
          ]
        ),
        collapse = ", "
      )
    )
  }

  data.table::setnames(
    expression_dt,
    old = sample_columns,
    new = normalized_sample_ids
  )

  sample_columns <- normalized_sample_ids

  expression_dt[
    ,
    (sample_columns) := lapply(
      .SD,
      function(x) {
        suppressWarnings(
          as.numeric(x)
        )
      }
    ),
    .SDcols = sample_columns
  ]

  expression_dt <- expression_dt[
    !is.na(Hugo_Symbol) &
      Hugo_Symbol != ""
  ]

  expression_dt <- expression_dt[
    ,
    lapply(
      .SD,
      mean,
      na.rm = TRUE
    ),
    by = Hugo_Symbol,
    .SDcols = sample_columns
  ]

  gene_names <- expression_dt$Hugo_Symbol

  expression_matrix <- as.matrix(
    expression_dt[
      ,
      sample_columns,
      with = FALSE
    ]
  )

  rownames(
    expression_matrix
  ) <- gene_names

  storage.mode(
    expression_matrix
  ) <- "numeric"

  expression_matrix[
    !is.finite(
      expression_matrix
    )
  ] <- NA_real_

  list(
    matrix = expression_matrix,
    transformation = "No additional log transformation",
    original_sample_ids = original_sample_ids,
    normalized_sample_ids = normalized_sample_ids
  )
}

load_imvigor_clinical <- function(
  sample_file,
  patient_file
) {
  sample_clinical <- data.table::fread(
    sample_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    comment.char = "#",
    check.names = FALSE,
    na.strings = c(
      "",
      "NA"
    )
  )

  patient_clinical <- data.table::fread(
    patient_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    comment.char = "#",
    check.names = FALSE,
    na.strings = c(
      "",
      "NA"
    )
  )

  required_sample_columns <- c(
    "PATIENT_ID",
    "SAMPLE_ID",
    "RESPONSE"
  )

  required_patient_columns <- c(
    "PATIENT_ID",
    "OS_STATUS",
    "OS_MONTHS"
  )

  missing_sample_columns <- setdiff(
    required_sample_columns,
    names(sample_clinical)
  )

  missing_patient_columns <- setdiff(
    required_patient_columns,
    names(patient_clinical)
  )

  if (length(missing_sample_columns) > 0) {
    stop(
      "IMvigor210 sample clinical file is missing: ",
      paste(
        missing_sample_columns,
        collapse = ", "
      )
    )
  }

  if (length(missing_patient_columns) > 0) {
    stop(
      "IMvigor210 patient clinical file is missing: ",
      paste(
        missing_patient_columns,
        collapse = ", "
      )
    )
  }

  sample_clinical <- sample_clinical %>%
    mutate(
      PATIENT_ID = normalize_sample_id(
        PATIENT_ID
      ),
      SAMPLE_ID = normalize_sample_id(
        SAMPLE_ID
      ),
      SAMPLE_KEY = SAMPLE_ID
    )

  patient_clinical <- patient_clinical %>%
    mutate(
      PATIENT_ID = normalize_sample_id(
        PATIENT_ID
      )
    )

  has_tmb <- "TMB_NONSYNONYMOUS" %in%
    names(sample_clinical)

  sample_clinical %>%
    left_join(
      patient_clinical,
      by = "PATIENT_ID",
      suffix = c(
        "_sample",
        "_patient"
      )
    ) %>%
    mutate(
      RESPONSE_SHORT = case_when(
        grepl(
          "Complete",
          RESPONSE,
          ignore.case = TRUE
        ) ~ "CR",
        grepl(
          "Partial",
          RESPONSE,
          ignore.case = TRUE
        ) ~ "PR",
        grepl(
          "Stable",
          RESPONSE,
          ignore.case = TRUE
        ) ~ "SD",
        grepl(
          "Progressive",
          RESPONSE,
          ignore.case = TRUE
        ) ~ "PD",
        TRUE ~ NA_character_
      ),
      RESPONSE_BINARY = case_when(
        RESPONSE_SHORT %in% c(
          "CR",
          "PR"
        ) ~ "Responder",
        RESPONSE_SHORT %in% c(
          "SD",
          "PD"
        ) ~ "Non-responder",
        TRUE ~ NA_character_
      ),
      CLINICAL_BENEFIT_GROUP = case_when(
        RESPONSE_SHORT %in% c(
          "CR",
          "PR",
          "SD"
        ) ~ "Clinical benefit",
        RESPONSE_SHORT == "PD" ~
          "No clinical benefit",
        TRUE ~ NA_character_
      ),
      OS_EVENT = case_when(
        grepl(
          "DECEASED|DEAD",
          OS_STATUS,
          ignore.case = TRUE
        ) ~ 1,
        grepl(
          "LIVING|ALIVE",
          OS_STATUS,
          ignore.case = TRUE
        ) ~ 0,
        TRUE ~ NA_real_
      ),
      OS_MONTHS = safe_numeric(
        OS_MONTHS
      ),
      OS_YEARS = OS_MONTHS / 12,
      TMB_NONSYNONYMOUS =
        if (has_tmb) {
          safe_numeric(
            TMB_NONSYNONYMOUS
          )
        } else {
          NA_real_
        }
    )
}

load_gse78220_expression <- function(path) {
  expression_df <- readxl::read_excel(
    path,
    sheet = 1
  )

  names(
    expression_df
  )[1] <- "Gene"

  expression_df <- expression_df %>%
    mutate(
      Gene = toupper(
        trimws(
          as.character(Gene)
        )
      )
    ) %>%
    filter(
      !is.na(Gene),
      Gene != ""
    ) %>%
    group_by(Gene) %>%
    summarise(
      across(
        where(is.numeric),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      .groups = "drop"
    )

  if (GSE78220_BASELINE_ONLY) {
    keep_columns <- c(
      "Gene",
      names(expression_df)[
        grepl(
          "baseline",
          names(expression_df),
          ignore.case = TRUE
        )
      ]
    )

    expression_df <- expression_df[
      ,
      keep_columns,
      drop = FALSE
    ]
  }

  gene_names <- expression_df$Gene

  expression_matrix <- as.matrix(
    expression_df[
      ,
      -1,
      drop = FALSE
    ]
  )

  rownames(
    expression_matrix
  ) <- gene_names

  storage.mode(
    expression_matrix
  ) <- "numeric"

  transformed <- maybe_log2(
    expression_matrix
  )

  transformed
}

prepare_gse78220_metadata <- function(
  metadata
) {
  find_metadata_column <- function(
    data_names,
    exact_names = character(),
    regex_patterns = character(),
    required = FALSE,
    field_label = "metadata field"
  ) {
    exact_hit <- exact_names[
      exact_names %in% data_names
    ]

    if (length(exact_hit) > 0) {
      return(exact_hit[1])
    }

    for (pattern in regex_patterns) {
      regex_hit <- grep(
        pattern,
        data_names,
        ignore.case = TRUE,
        value = TRUE
      )

      if (length(regex_hit) > 0) {
        return(regex_hit[1])
      }
    }

    if (required) {
      stop(
        "Could not detect ",
        field_label,
        ".\nAvailable GSE78220 metadata columns are:\n",
        paste(
          data_names,
          collapse = ", "
        )
      )
    }

    NA_character_
  }

  metadata_names <- names(
    metadata
  )

  patient_column <- find_metadata_column(
    metadata_names,
    exact_names = c(
      "PATIENT_ID",
      "PATIENT"
    ),
    regex_patterns = c(
      "^PATIENT.*ID",
      "^SUBJECT.*ID"
    ),
    required = TRUE,
    field_label = "patient ID"
  )

  response_column <- find_metadata_column(
    metadata_names,
    exact_names = c(
      "ANTI_PD_1_RESPONSE",
      "ANTI_PD1_RESPONSE",
      "RESPONSE"
    ),
    regex_patterns = c(
      "ANTI.*PD.*1.*RESPONSE",
      "BEST.*RESPONSE",
      "CLINICAL.*RESPONSE"
    ),
    required = TRUE,
    field_label = "anti-PD-1 response"
  )

  os_time_column <- find_metadata_column(
    metadata_names,
    exact_names = c(
      "OVERALL_SURVIVAL_DAYS",
      "OVERALL_SURVIVAL_DAY",
      "OS_DAYS",
      "OS_TIME"
    ),
    regex_patterns = c(
      "OVERALL.*SURVIVAL.*DAY",
      "^OS.*DAY",
      "OVERALL.*SURVIVAL.*TIME",
      "^OS.*TIME"
    ),
    required = FALSE,
    field_label = "overall-survival time"
  )

  vital_column <- find_metadata_column(
    metadata_names,
    exact_names = c(
      "VITAL_STATUS",
      "OS_STATUS",
      "SURVIVAL_STATUS"
    ),
    regex_patterns = c(
      "VITAL.*STATUS",
      "OVERALL.*SURVIVAL.*STATUS",
      "^OS.*STATUS",
      "DEATH.*STATUS"
    ),
    required = FALSE,
    field_label = "vital status"
  )

  age_column <- find_metadata_column(
    metadata_names,
    exact_names = c(
      "AGE_YRS",
      "AGE_YEARS",
      "AGE"
    ),
    regex_patterns = c(
      "^AGE.*YR",
      "^AGE.*YEAR",
      "^AGE$"
    ),
    required = FALSE
  )

  sex_column <- find_metadata_column(
    metadata_names,
    exact_names = c(
      "GENDER",
      "SEX"
    ),
    regex_patterns = c(
      "^GENDER",
      "^SEX"
    ),
    required = FALSE
  )

  biopsy_column <- find_metadata_column(
    metadata_names,
    exact_names = c(
      "BIOPSY_TIME",
      "TREATMENT_TIMEPOINT"
    ),
    regex_patterns = c(
      "BIOPSY.*TIME",
      "TREATMENT.*TIME",
      "TIMEPOINT"
    ),
    required = FALSE
  )

  detected_columns <- data.frame(
    FIELD = c(
      "Patient ID",
      "Response",
      "OS time",
      "Vital status",
      "Age",
      "Sex",
      "Biopsy time"
    ),
    DETECTED_COLUMN = c(
      patient_column,
      response_column,
      os_time_column,
      vital_column,
      age_column,
      sex_column,
      biopsy_column
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    detected_columns,
    file.path(
      output_directory,
      "GSE78220_detected_metadata_columns.csv"
    ),
    row.names = FALSE
  )

  patient_values <- as.character(
    metadata[[patient_column]]
  )

  response_values <- as.character(
    metadata[[response_column]]
  )

  os_days_values <- if (
    !is.na(os_time_column)
  ) {
    safe_numeric(
      metadata[[os_time_column]]
    )
  } else {
    rep(
      NA_real_,
      nrow(metadata)
    )
  }

  vital_values <- if (
    !is.na(vital_column)
  ) {
    as.character(
      metadata[[vital_column]]
    )
  } else {
    rep(
      NA_character_,
      nrow(metadata)
    )
  }

  age_values <- if (
    !is.na(age_column)
  ) {
    safe_numeric(
      metadata[[age_column]]
    )
  } else {
    rep(
      NA_real_,
      nrow(metadata)
    )
  }

  sex_values <- if (
    !is.na(sex_column)
  ) {
    as.character(
      metadata[[sex_column]]
    )
  } else {
    rep(
      NA_character_,
      nrow(metadata)
    )
  }

  biopsy_values <- if (
    !is.na(biopsy_column)
  ) {
    as.character(
      metadata[[biopsy_column]]
    )
  } else {
    rep(
      NA_character_,
      nrow(metadata)
    )
  }

  prepared <- metadata %>%
    mutate(
      PATIENT_ID = patient_values,
      RESPONSE_RAW = response_values,
      RESPONSE_SHORT = case_when(
        grepl(
          "Complete|\\bCR\\b",
          RESPONSE_RAW,
          ignore.case = TRUE
        ) ~ "CR",
        grepl(
          "Partial|\\bPR\\b",
          RESPONSE_RAW,
          ignore.case = TRUE
        ) ~ "PR",
        grepl(
          "Stable|\\bSD\\b",
          RESPONSE_RAW,
          ignore.case = TRUE
        ) ~ "SD",
        grepl(
          "Progressive|\\bPD\\b",
          RESPONSE_RAW,
          ignore.case = TRUE
        ) ~ "PD",
        TRUE ~ NA_character_
      ),
      RESPONSE_BINARY = case_when(
        RESPONSE_SHORT %in% c(
          "CR",
          "PR"
        ) ~ "Responder",
        RESPONSE_SHORT %in% c(
          "SD",
          "PD"
        ) ~ "Non-responder",
        TRUE ~ NA_character_
      ),
      CLINICAL_BENEFIT_GROUP = case_when(
        RESPONSE_SHORT %in% c(
          "CR",
          "PR",
          "SD"
        ) ~ "Clinical benefit",
        RESPONSE_SHORT == "PD" ~
          "No clinical benefit",
        TRUE ~ NA_character_
      ),
      OS_DAYS = os_days_values,
      OS_MONTHS = OS_DAYS /
        30.4375,
      OS_YEARS = OS_DAYS /
        365.25,
      VITAL_STATUS_V5 = vital_values,
      OS_EVENT = case_when(
        grepl(
          "Dead|Deceased|Death",
          VITAL_STATUS_V5,
          ignore.case = TRUE
        ) ~ 1,
        grepl(
          "Alive|Living",
          VITAL_STATUS_V5,
          ignore.case = TRUE
        ) ~ 0,
        TRUE ~ NA_real_
      ),
      AGE = age_values,
      SEX = sex_values,
      BIOPSY_TIME_V5 = biopsy_values
    )

  if (
    sum(
      !is.na(
        prepared$RESPONSE_BINARY
      )
    ) == 0
  ) {
    stop(
      "The response column was detected as '",
      response_column,
      "', but no CR/PR/SD/PD labels could be parsed.\n",
      "Unique response values:\n",
      paste(
        unique(response_values),
        collapse = ", "
      )
    )
  }

  if (
    all(
      !is.finite(
        prepared$OS_DAYS
      )
    )
  ) {
    warning(
      "No usable GSE78220 OS time was detected. ",
      "Response analyses will still run; the GSE78220 KM panel ",
      "will be replaced by an explanatory placeholder."
    )
  }

  prepared
}

aggregate_gse78220_patient <- function(
  sample_table
) {
  if (!AVERAGE_REPEATED_BASELINE_BIOPSIES) {
    sample_table <- sample_table %>%
      arrange(
        PATIENT_ID,
        SAMPLE_ID
      ) %>%
      group_by(PATIENT_ID) %>%
      slice(1) %>%
      ungroup()
  }

  sample_table %>%
    group_by(PATIENT_ID) %>%
    summarise(
      N_BASELINE_BIOPSIES = n(),
      STEP_UP = mean(
        STEP_UP,
        na.rm = TRUE
      ),
      CAF3 = mean(
        CAF3,
        na.rm = TRUE
      ),
      RESPONSE_SHORT = first_non_missing(
        RESPONSE_SHORT
      ),
      RESPONSE_BINARY = first_non_missing(
        RESPONSE_BINARY
      ),
      CLINICAL_BENEFIT_GROUP =
        first_non_missing(
          CLINICAL_BENEFIT_GROUP
        ),
      OS_MONTHS = safe_numeric(
        first_non_missing(
          OS_MONTHS
        )
      ),
      OS_YEARS = safe_numeric(
        first_non_missing(
          OS_YEARS
        )
      ),
      OS_EVENT = safe_numeric(
        first_non_missing(
          OS_EVENT
        )
      ),
      AGE = safe_numeric(
        first_non_missing(
          AGE
        )
      ),
      SEX = first_non_missing(
        SEX
      ),
      .groups = "drop"
    ) %>%
    {
      adjustment_fit <- lm(
        STEP_UP ~ CAF3,
        data = .
      )

      mutate(
        .,
        STEP_CAF3_adjust =
          residuals(
            adjustment_fit
          )
      )
    }
}

add_median_group <- function(
  data,
  score_name
) {
  cutpoint <- median(
    data[[score_name]],
    na.rm = TRUE
  )

  data %>%
    mutate(
      SCORE_GROUP = ifelse(
        .data[[score_name]] >=
          cutpoint,
        "High",
        "Low"
      ),
      SCORE_GROUP = factor(
        SCORE_GROUP,
        levels = c(
          "Low",
          "High"
        )
      )
    )
}

make_km_panel <- function(
  data,
  score_name,
  cohort_title,
  output_directory
) {
  survival_data <- data %>%
    filter(
      is.finite(OS_YEARS),
      OS_YEARS > 0,
      OS_EVENT %in% c(0, 1),
      is.finite(.data[[score_name]])
    )

  if (
    nrow(survival_data) < 4 ||
    length(
      unique(
        survival_data$OS_EVENT
      )
    ) < 2
  ) {
    placeholder <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.58,
        label = cohort_title,
        fontface = "bold",
        size = 5
      ) +
      annotate(
        "text",
        x = 0.5,
        y = 0.42,
        label = "Usable OS time/status not found\nin the current clinical file",
        size = 4
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void()

    save_ggplot(
      placeholder,
      output_directory,
      paste0(
        "Panel_KM_",
        sanitize_filename(cohort_title),
        "_",
        score_name,
        "_OS_unavailable"
      ),
      width = 7.2,
      height = 7
    )

    return(
      list(
        panel = placeholder,
        statistics = data.frame(
          COHORT = cohort_title,
          SCORE = score_name,
          N = nrow(survival_data),
          EVENTS = sum(
            survival_data$OS_EVENT == 1,
            na.rm = TRUE
          ),
          LOGRANK_P = NA_real_,
          HR_HIGH_VS_LOW = NA_real_,
          HR_HIGH_VS_LOW_CI_LOW = NA_real_,
          HR_HIGH_VS_LOW_CI_HIGH = NA_real_,
          HR_HIGH_VS_LOW_P = NA_real_,
          HR_PER_SD = NA_real_,
          HR_PER_SD_CI_LOW = NA_real_,
          HR_PER_SD_CI_HIGH = NA_real_,
          HR_PER_SD_P = NA_real_
        ),
        data = survival_data
      )
    )
  }

  survival_data <- survival_data %>%
    add_median_group(
      score_name
    )

  km_fit <- survfit(
    Surv(
      OS_YEARS,
      OS_EVENT
    ) ~ SCORE_GROUP,
    data = survival_data
  )

  logrank <- survdiff(
    Surv(
      OS_YEARS,
      OS_EVENT
    ) ~ SCORE_GROUP,
    data = survival_data
  )

  logrank_p <- 1 -
    pchisq(
      logrank$chisq,
      df = length(
        logrank$n
      ) - 1
    )

  cox_group <- coxph(
    Surv(
      OS_YEARS,
      OS_EVENT
    ) ~ SCORE_GROUP,
    data = survival_data
  )

  continuous_data <- survival_data %>%
    mutate(
      SCORE_Z = as.numeric(
        scale(
          .data[[score_name]]
        )
      )
    )

  cox_continuous <- coxph(
    Surv(
      OS_YEARS,
      OS_EVENT
    ) ~ SCORE_Z,
    data = continuous_data
  )

  group_summary <- summary(
    cox_group
  )

  continuous_summary <- summary(
    cox_continuous
  )

  statistics <- data.frame(
    COHORT = cohort_title,
    SCORE = score_name,
    N = nrow(survival_data),
    EVENTS = sum(
      survival_data$OS_EVENT == 1
    ),
    LOGRANK_P = logrank_p,
    HR_HIGH_VS_LOW =
      group_summary$coefficients[
        1,
        "exp(coef)"
      ],
    HR_HIGH_VS_LOW_CI_LOW =
      group_summary$conf.int[
        1,
        "lower .95"
      ],
    HR_HIGH_VS_LOW_CI_HIGH =
      group_summary$conf.int[
        1,
        "upper .95"
      ],
    HR_HIGH_VS_LOW_P =
      group_summary$coefficients[
        1,
        "Pr(>|z|)"
      ],
    HR_PER_SD =
      continuous_summary$coefficients[
        1,
        "exp(coef)"
      ],
    HR_PER_SD_CI_LOW =
      continuous_summary$conf.int[
        1,
        "lower .95"
      ],
    HR_PER_SD_CI_HIGH =
      continuous_summary$conf.int[
        1,
        "upper .95"
      ],
    HR_PER_SD_P =
      continuous_summary$coefficients[
        1,
        "Pr(>|z|)"
      ],
    stringsAsFactors = FALSE
  )

  km_plot <- ggsurvplot(
    km_fit,
    data = survival_data,
    pval = paste0(
      "Log-rank ",
      format_p(logrank_p)
    ),
    risk.table = TRUE,
    conf.int = FALSE,
    censor = TRUE,
    xlab = "Time (years)",
    ylab = "Overall survival probability",
    legend.title = score_name,
    legend.labs = c(
      "Low",
      "High"
    ),
    risk.table.height = 0.25,
    ggtheme = theme_classic(
      base_size = 12
    ),
    tables.theme = theme_cleantable(
      base_size = 9
    )
  )

  km_plot$plot <- km_plot$plot +
    ggtitle(
      cohort_title
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      legend.position = c(
        0.68,
        0.86
      )
    )

  km_grob <- gridExtra::arrangeGrob(
    km_plot$plot,
    km_plot$table,
    ncol = 1,
    heights = c(
      3.2,
      1
    )
  )

  panel <- cowplot::ggdraw(
    km_grob
  )

  save_ggplot(
    panel,
    output_directory,
    paste0(
      "Panel_KM_",
      sanitize_filename(
        cohort_title
      ),
      "_",
      score_name
    ),
    width = 7.2,
    height = 7
  )

  list(
    panel = panel,
    statistics = statistics,
    data = survival_data
  )
}

make_response_panel <- function(
  data,
  score_name,
  cohort_title,
  output_directory
) {
  response_data <- data %>%
    filter(
      RESPONSE_BINARY %in% c(
        "Responder",
        "Non-responder"
      ),
      is.finite(
        .data[[score_name]]
      )
    ) %>%
    add_median_group(
      score_name
    ) %>%
    mutate(
      RESPONSE_BINARY = factor(
        RESPONSE_BINARY,
        levels = c(
          "Non-responder",
          "Responder"
        )
      )
    )

  fisher_table <- with(
    response_data,
    table(
      factor(
        SCORE_GROUP,
        levels = c(
          "Low",
          "High"
        )
      ),
      factor(
        RESPONSE_BINARY,
        levels = c(
          "Responder",
          "Non-responder"
        )
      )
    )
  )

  fisher_result <- fisher.test(
    fisher_table
  )

  response_count_raw <- response_data %>%
    transmute(
      SCORE_GROUP = as.character(
        SCORE_GROUP
      ),
      RESPONSE_BINARY = as.character(
        RESPONSE_BINARY
      )
    ) %>%
    count(
      SCORE_GROUP,
      RESPONSE_BINARY,
      name = "n"
    )

  response_grid <- tidyr::expand_grid(
    SCORE_GROUP = c(
      "Low",
      "High"
    ),
    RESPONSE_BINARY = c(
      "Non-responder",
      "Responder"
    )
  )

  response_summary <- response_grid %>%
    left_join(
      response_count_raw,
      by = c(
        "SCORE_GROUP",
        "RESPONSE_BINARY"
      )
    ) %>%
    mutate(
      n = tidyr::replace_na(
        n,
        0L
      ),
      SCORE_GROUP = factor(
        SCORE_GROUP,
        levels = c(
          "Low",
          "High"
        )
      ),
      RESPONSE_BINARY = factor(
        RESPONSE_BINARY,
        levels = c(
          "Non-responder",
          "Responder"
        )
      )
    ) %>%
    group_by(
      SCORE_GROUP
    ) %>%
    mutate(
      TOTAL = sum(n),
      PROPORTION = ifelse(
        TOTAL > 0,
        n / TOTAL,
        0
      ),
      LABEL = paste0(
        percent(
          PROPORTION,
          accuracy = 1
        )
      )
    ) %>%
    ungroup()

  bar_plot <- ggplot(
    response_summary,
    aes(
      SCORE_GROUP,
      PROPORTION,
      fill = RESPONSE_BINARY
    )
  ) +
    geom_col(
      width = 0.68
    ) +
    geom_text(
      aes(label = LABEL),
      position = position_stack(
        vjust = 0.5
      ),
      size = 4
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      limits = c(
        0,
        1
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    labs(
      title = cohort_title,
      x = score_name,
      y = "Percent",
      fill = NULL
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      legend.position = "right"
    )

  count_table_raw <- response_data %>%
    transmute(
      OUTCOME = dplyr::case_when(
        as.character(
          RESPONSE_BINARY
        ) == "Responder" ~ "CR/PR",
        as.character(
          RESPONSE_BINARY
        ) == "Non-responder" ~ "SD/PD",
        TRUE ~ NA_character_
      ),
      SCORE_GROUP = as.character(
        SCORE_GROUP
      )
    ) %>%
    filter(
      !is.na(OUTCOME),
      SCORE_GROUP %in% c(
        "Low",
        "High"
      )
    ) %>%
    count(
      OUTCOME,
      SCORE_GROUP,
      name = "Number"
    )

  count_table_grid <- tidyr::expand_grid(
    OUTCOME = as.character(
      c(
        "CR/PR",
        "SD/PD"
      )
    ),
    SCORE_GROUP = as.character(
      c(
        "Low",
        "High"
      )
    )
  )

  count_table <- count_table_grid %>%
    left_join(
      count_table_raw,
      by = c(
        "OUTCOME",
        "SCORE_GROUP"
      )
    ) %>%
    mutate(
      Number = tidyr::replace_na(
        Number,
        0L
      ),
      OUTCOME = factor(
        OUTCOME,
        levels = c(
          "CR/PR",
          "SD/PD"
        )
      ),
      SCORE_GROUP = factor(
        SCORE_GROUP,
        levels = c(
          "Low",
          "High"
        )
      )
    ) %>%
    arrange(
      OUTCOME,
      SCORE_GROUP
    ) %>%
    pivot_wider(
      names_from = SCORE_GROUP,
      values_from = Number,
      values_fill = 0
    ) %>%
    mutate(
      OUTCOME = as.character(
        OUTCOME
      )
    )

  names(count_table) <- c(
    "Number",
    "Low score",
    "High score"
  )

  table_theme <- gridExtra::ttheme_minimal(
    base_size = 10,
    core = list(
      fg_params = list(
        hjust = 0.5,
        x = 0.5
      )
    ),
    colhead = list(
      fg_params = list(
        fontface = "bold"
      )
    )
  )

  table_grob <- gridExtra::tableGrob(
    count_table,
    rows = NULL,
    theme = table_theme
  )

  p_text <- grid::textGrob(
    paste0(
      "Fisher ",
      format_p(
        fisher_result$p.value
      ),
      "; OR (Low vs High) = ",
      round(
        unname(
          fisher_result$estimate
        ),
        2
      ),
      " [",
      round(
        fisher_result$conf.int[1],
        2
      ),
      ", ",
      round(
        fisher_result$conf.int[2],
        2
      ),
      "]"
    ),
    gp = grid::gpar(
      fontsize = 11,
      fontface = "bold"
    )
  )

  response_grob <- gridExtra::arrangeGrob(
    bar_plot,
    table_grob,
    p_text,
    ncol = 1,
    heights = c(
      3.2,
      1.2,
      0.45
    )
  )

  panel <- cowplot::ggdraw(
    response_grob
  )

  statistics <- data.frame(
    COHORT = cohort_title,
    SCORE = score_name,
    N = nrow(response_data),
    N_RESPONDER = sum(
      response_data$RESPONSE_BINARY ==
        "Responder"
    ),
    N_NON_RESPONDER = sum(
      response_data$RESPONSE_BINARY ==
        "Non-responder"
    ),
    FISHER_OR_LOW_VS_HIGH =
      unname(
        fisher_result$estimate
      ),
    FISHER_CI_LOW =
      fisher_result$conf.int[1],
    FISHER_CI_HIGH =
      fisher_result$conf.int[2],
    FISHER_P =
      fisher_result$p.value,
    stringsAsFactors = FALSE
  )

  save_ggplot(
    panel,
    output_directory,
    paste0(
      "Panel_Response_",
      sanitize_filename(
        cohort_title
      ),
      "_",
      score_name
    ),
    width = 6.5,
    height = 7
  )

  list(
    panel = panel,
    statistics = statistics,
    data = response_data,
    counts = count_table
  )
}

make_pairwise_holm_table <- function(
  data,
  group_column,
  score_column,
  group_levels
) {
  present_levels <- group_levels[
    group_levels %in%
      unique(
        as.character(
          data[[group_column]]
        )
      )
  ]

  if (length(present_levels) < 2) {
    return(
      data.frame()
    )
  }

  comparisons <- combn(
    present_levels,
    2,
    simplify = FALSE
  )

  raw_p <- vapply(
    comparisons,
    function(pair) {
      group_values <- as.character(
        data[[group_column]]
      )

      score_values <- suppressWarnings(
        as.numeric(
          data[[score_column]]
        )
      )

      values_1 <- score_values[
        group_values == pair[1] &
          is.finite(score_values)
      ]

      values_2 <- score_values[
        group_values == pair[2] &
          is.finite(score_values)
      ]

      if (
        length(values_1) == 0 ||
        length(values_2) == 0
      ) {
        return(NA_real_)
      }

      suppressWarnings(
        wilcox.test(
          x = values_1,
          y = values_2,
          exact = FALSE
        )$p.value
      )
    },
    numeric(1)
  )

  adjusted_p <- p.adjust(
    raw_p,
    method = "holm"
  )

  score_range <- range(
    data[[score_column]],
    na.rm = TRUE
  )

  span <- diff(
    score_range
  )

  if (
    !is.finite(span) ||
      span == 0
  ) {
    span <- 1
  }

  data.frame(
    group1 = vapply(
      comparisons,
      `[`,
      character(1),
      1
    ),
    group2 = vapply(
      comparisons,
      `[`,
      character(1),
      2
    ),
    p = raw_p,
    p.adj = adjusted_p,
    label = vapply(
      adjusted_p,
      function(p) {
        if (p < 0.001) {
          "p.adj<0.001"
        } else {
          paste0(
            "p.adj=",
            format.pval(
              p,
              digits = 2
            )
          )
        }
      },
      character(1)
    ),
    y.position = score_range[2] +
      span *
      seq(
        0.12,
        by = 0.12,
        length.out =
          length(comparisons)
      ),
    stringsAsFactors = FALSE
  )
}

make_violin_panel <- function(
  data,
  score_name,
  cohort_title,
  output_directory,
  four_group = FALSE
) {
  plot_data <- data %>%
    filter(
      !is.na(RESPONSE_SHORT),
      is.finite(
        .data[[score_name]]
      )
    )

  if (four_group) {
    response_levels <- c(
      "CR",
      "PR",
      "PD",
      "SD"
    )

    plot_data <- plot_data %>%
      mutate(
        DISPLAY_GROUP = factor(
          RESPONSE_SHORT,
          levels = response_levels
        )
      )
  } else {
    contains_sd <- any(
      plot_data$RESPONSE_SHORT ==
        "SD",
      na.rm = TRUE
    )

    nonresponse_label <- if (
      contains_sd
    ) {
      "SD/PD"
    } else {
      "PD"
    }

    plot_data <- plot_data %>%
      mutate(
        DISPLAY_GROUP = case_when(
          RESPONSE_SHORT %in% c(
            "CR",
            "PR"
          ) ~ "CR/PR",
          RESPONSE_SHORT %in% c(
            "SD",
            "PD"
          ) ~ nonresponse_label,
          TRUE ~ NA_character_
        ),
        DISPLAY_GROUP = factor(
          DISPLAY_GROUP,
          levels = c(
            nonresponse_label,
            "CR/PR"
          )
        )
      ) %>%
      filter(
        !is.na(DISPLAY_GROUP)
      )

    response_levels <- levels(
      plot_data$DISPLAY_GROUP
    )
  }

  overall_test <- if (
    length(
      unique(
        na.omit(
          plot_data$DISPLAY_GROUP
        )
      )
    ) > 2
  ) {
    kruskal.test(
      x = plot_data[[score_name]],
      g = plot_data$DISPLAY_GROUP
    )
  } else {
    group_levels_present <- levels(
      droplevels(
        plot_data$DISPLAY_GROUP
      )
    )

    wilcox.test(
      plot_data[[score_name]][
        plot_data$DISPLAY_GROUP ==
          group_levels_present[1]
      ],
      plot_data[[score_name]][
        plot_data$DISPLAY_GROUP ==
          group_levels_present[2]
      ],
      exact = FALSE
    )
  }

  pairwise_table <- make_pairwise_holm_table(
    plot_data,
    "DISPLAY_GROUP",
    score_name,
    response_levels
  )

  violin_plot <- ggplot(
    plot_data,
    aes(
      DISPLAY_GROUP,
      .data[[score_name]],
      fill = DISPLAY_GROUP
    )
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.92
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = "white"
    ) +
    geom_jitter(
      width = 0.08,
      size = 1.5,
      alpha = 0.45
    ) +
    labs(
      title = cohort_title,
      subtitle = paste0(
        if (
          length(
            unique(
              na.omit(
                plot_data$DISPLAY_GROUP
              )
            )
          ) > 2
        ) {
          "Kruskal-Wallis "
        } else {
          "Wilcoxon "
        },
        format_p(
          overall_test$p.value
        ),
        "; pairwise Holm adjustment"
      ),
      x = NULL,
      y = score_name,
      fill = NULL
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      legend.position = "none",
      axis.text.x = element_text(
        size = 11
      )
    )

  pairwise_table_for_plot <- pairwise_table %>%
    filter(
      is.finite(p.adj),
      !is.na(group1),
      !is.na(group2),
      is.finite(y.position)
    )

  if (nrow(pairwise_table_for_plot) > 0) {
    violin_plot <- violin_plot +
      ggpubr::stat_pvalue_manual(
        pairwise_table_for_plot,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.4,
        size = 3,
        hide.ns = FALSE
      ) +
      scale_y_continuous(
        expand = expansion(
          mult = c(
            0.05,
            0.35
          )
        )
      )
  }

  save_ggplot(
    violin_plot,
    output_directory,
    paste0(
      "Panel_Violin_",
      sanitize_filename(
        cohort_title
      ),
      "_",
      score_name
    ),
    width = 6,
    height = 7
  )

  statistics <- data.frame(
    COHORT = cohort_title,
    SCORE = score_name,
    OVERALL_TEST = if (
      length(
        unique(
          na.omit(
            plot_data$DISPLAY_GROUP
          )
        )
      ) > 2
    ) {
      "Kruskal-Wallis"
    } else {
      "Wilcoxon"
    },
    OVERALL_P =
      overall_test$p.value,
    stringsAsFactors = FALSE
  )

  list(
    panel = violin_plot,
    statistics = statistics,
    pairwise = pairwise_table,
    data = plot_data
  )
}

make_six_panel_figure <- function(
  imvigor_data,
  gse78220_data,
  method,
  score_name,
  output_directory
) {
  combination_dir <- file.path(
    output_directory,
    method,
    score_name
  )

  dir.create(
    combination_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  message("  Building IMvigor210 KM panel: ", method, " / ", score_name)
  panel_a <- make_km_panel(
    imvigor_data,
    score_name,
    "IMvigor210\n(PD-L1 inhibitor cohort)",
    combination_dir
  )

  message("  Building IMvigor210 response panel: ", method, " / ", score_name)
  panel_b <- make_response_panel(
    imvigor_data,
    score_name,
    "IMvigor210",
    combination_dir
  )

  message("  Building IMvigor210 violin panel: ", method, " / ", score_name)
  panel_c <- make_violin_panel(
    imvigor_data,
    score_name,
    "IMvigor210",
    combination_dir,
    four_group = TRUE
  )

  message("  Building GSE78220 KM panel: ", method, " / ", score_name)
  panel_d <- make_km_panel(
    gse78220_data,
    score_name,
    "GSE78220\n(PD-1 inhibitor cohort)",
    combination_dir
  )

  message("  Building GSE78220 response panel: ", method, " / ", score_name)
  panel_e <- make_response_panel(
    gse78220_data,
    score_name,
    "GSE78220",
    combination_dir
  )

  message("  Building GSE78220 violin panel: ", method, " / ", score_name)
  panel_f <- make_violin_panel(
    gse78220_data,
    score_name,
    "GSE78220",
    combination_dir,
    four_group = FALSE
  )

  top_row <- cowplot::plot_grid(
    panel_a$panel,
    panel_b$panel,
    panel_c$panel,
    nrow = 1,
    labels = c(
      "A",
      "B",
      "C"
    ),
    label_size = 18,
    rel_widths = c(
      1.25,
      1,
      1
    ),
    align = "h",
    axis = "tb"
  )

  bottom_row <- cowplot::plot_grid(
    panel_d$panel,
    panel_e$panel,
    panel_f$panel,
    nrow = 1,
    labels = c(
      "D",
      "E",
      "F"
    ),
    label_size = 18,
    rel_widths = c(
      1.25,
      1,
      1
    ),
    align = "h",
    axis = "tb"
  )

  full_figure <- cowplot::plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(
      1,
      1
    )
  )

  full_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "STEP and immunotherapy outcome: ",
        score_name,
        " (",
        method,
        ")"
      ),
      fontface = "bold",
      size = 18,
      x = 0.5,
      hjust = 0.5
    )

  final_figure <- cowplot::plot_grid(
    full_title,
    full_figure,
    ncol = 1,
    rel_heights = c(
      0.055,
      1
    )
  )

  save_ggplot(
    final_figure,
    combination_dir,
    paste0(
      "Figure_6panel_",
      method,
      "_",
      score_name
    ),
    width = 18,
    height = 13
  )

  survival_statistics <- bind_rows(
    panel_a$statistics,
    panel_d$statistics
  )

  response_statistics <- bind_rows(
    panel_b$statistics,
    panel_e$statistics
  )

  violin_statistics <- bind_rows(
    panel_c$statistics,
    panel_f$statistics
  )

  pairwise_statistics <- bind_rows(
    transform(
      panel_c$pairwise,
      COHORT = "IMvigor210",
      SCORE = score_name
    ),
    transform(
      panel_f$pairwise,
      COHORT = "GSE78220",
      SCORE = score_name
    )
  )

  write.csv(
    survival_statistics,
    file.path(
      combination_dir,
      "Statistics_OS_KM_Cox.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    response_statistics,
    file.path(
      combination_dir,
      "Statistics_Response_Fisher.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    violin_statistics,
    file.path(
      combination_dir,
      "Statistics_Response_Distribution.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    pairwise_statistics,
    file.path(
      combination_dir,
      "Statistics_Pairwise_Holm.csv"
    ),
    row.names = FALSE
  )

  list(
    figure = final_figure,
    survival = survival_statistics,
    response = response_statistics,
    violin = violin_statistics,
    pairwise = pairwise_statistics
  )
}

run_three_score_response_summary <- function(
  data,
  cohort_name,
  method,
  output_directory
) {
  output <- list()
  roc_list <- list()

  for (
    score_name in
      SCORES_TO_RUN
  ) {
    d <- data %>%
      filter(
        RESPONSE_BINARY %in% c(
          "Responder",
          "Non-responder"
        ),
        is.finite(
          .data[[score_name]]
        )
      ) %>%
      mutate(
        SCORE_Z = as.numeric(
          scale(
            .data[[score_name]]
          )
        ),
        RESPONSE_Y = ifelse(
          RESPONSE_BINARY ==
            "Responder",
          1,
          0
        )
      )

    logistic_fit <- glm(
      RESPONSE_Y ~ SCORE_Z,
      data = d,
      family = binomial()
    )

    coefficient <- summary(
      logistic_fit
    )$coefficients[
      "SCORE_Z",
      ,
      drop = FALSE
    ]

    roc_object <- pROC::roc(
      response = d$RESPONSE_Y,
      predictor = d[[score_name]],
      levels = c(
        0,
        1
      ),
      direction = "auto",
      quiet = TRUE
    )

    roc_list[[score_name]] <- roc_object

    auc_ci <- pROC::ci.auc(
      roc_object
    )

    output[[score_name]] <- data.frame(
      COHORT = cohort_name,
      METHOD = method,
      SCORE = score_name,
      N = nrow(d),
      OR_PER_SD = exp(
        coefficient[1, "Estimate"]
      ),
      OR_CI_LOW = exp(
        coefficient[1, "Estimate"] -
          1.96 *
          coefficient[1, "Std. Error"]
      ),
      OR_CI_HIGH = exp(
        coefficient[1, "Estimate"] +
          1.96 *
          coefficient[1, "Std. Error"]
      ),
      LOGISTIC_P = coefficient[
        1,
        "Pr(>|z|)"
      ],
      AUC = as.numeric(
        pROC::auc(
          roc_object
        )
      ),
      AUC_CI_LOW = as.numeric(
        auc_ci[1]
      ),
      AUC_CI_HIGH = as.numeric(
        auc_ci[3]
      ),
      stringsAsFactors = FALSE
    )
  }

  summary_table <- bind_rows(
    output
  )

  cohort_dir <- file.path(
    output_directory,
    method,
    "Deep_supplementary",
    cohort_name
  )

  dir.create(
    cohort_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    summary_table,
    file.path(
      cohort_dir,
      "Three_score_response_summary.csv"
    ),
    row.names = FALSE
  )

  combined_roc <- pROC::ggroc(
    roc_list,
    aes = "linetype",
    linewidth = 1.05
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 3
    ) +
    coord_equal() +
    labs(
      title = paste0(
        cohort_name,
        ": three-score ROC comparison"
      ),
      x = "1 - Specificity",
      y = "Sensitivity",
      linetype = "Score"
    ) +
    theme_classic(
      base_size = 13
    )

  save_ggplot(
    combined_roc,
    cohort_dir,
    "Three_score_combined_ROC",
    width = 7,
    height = 6.3
  )

  forest_data <- summary_table %>%
    mutate(
      SCORE = factor(
        SCORE,
        levels = rev(
          SCORES_TO_RUN
        )
      )
    )

  logistic_forest <- ggplot(
    forest_data,
    aes(
      OR_PER_SD,
      SCORE
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = 2
    ) +
    geom_errorbarh(
      aes(
        xmin = OR_CI_LOW,
        xmax = OR_CI_HIGH
      ),
      height = 0.18
    ) +
    geom_point(
      size = 3
    ) +
    scale_x_log10() +
    labs(
      title = paste0(
        cohort_name,
        ": response logistic regression"
      ),
      subtitle = "Odds ratio per 1-SD increase",
      x = "Odds ratio (log scale)",
      y = NULL
    ) +
    theme_classic(
      base_size = 13
    )

  save_ggplot(
    logistic_forest,
    cohort_dir,
    "Three_score_logistic_forest",
    width = 7,
    height = 5.2
  )

  summary_table
}

run_joint_step_caf_analysis <- function(
  data,
  cohort_name,
  method,
  output_directory
) {
  d <- data %>%
    filter(
      RESPONSE_BINARY %in% c(
        "Responder",
        "Non-responder"
      ),
      is.finite(STEP_UP),
      is.finite(CAF3)
    ) %>%
    mutate(
      STEP_GROUP = ifelse(
        STEP_UP >= median(
          STEP_UP,
          na.rm = TRUE
        ),
        "STEP high",
        "STEP low"
      ),
      CAF3_GROUP = ifelse(
        CAF3 >= median(
          CAF3,
          na.rm = TRUE
        ),
        "CAF3 high",
        "CAF3 low"
      ),
      JOINT_GROUP = paste(
        STEP_GROUP,
        CAF3_GROUP,
        sep = " / "
      ),
      JOINT_GROUP = factor(
        JOINT_GROUP,
        levels = c(
          "STEP low / CAF3 low",
          "STEP high / CAF3 low",
          "STEP low / CAF3 high",
          "STEP high / CAF3 high"
        )
      )
    )

  joint_count_raw <- d %>%
    transmute(
      JOINT_GROUP = as.character(
        JOINT_GROUP
      ),
      RESPONSE_BINARY = as.character(
        RESPONSE_BINARY
      )
    ) %>%
    count(
      JOINT_GROUP,
      RESPONSE_BINARY,
      name = "n"
    )

  joint_grid <- tidyr::expand_grid(
    JOINT_GROUP = levels(
      d$JOINT_GROUP
    ),
    RESPONSE_BINARY = c(
      "Non-responder",
      "Responder"
    )
  )

  joint_summary <- joint_grid %>%
    left_join(
      joint_count_raw,
      by = c(
        "JOINT_GROUP",
        "RESPONSE_BINARY"
      )
    ) %>%
    mutate(
      n = tidyr::replace_na(
        n,
        0L
      ),
      JOINT_GROUP = factor(
        JOINT_GROUP,
        levels = levels(
          d$JOINT_GROUP
        )
      ),
      RESPONSE_BINARY = factor(
        RESPONSE_BINARY,
        levels = c(
          "Non-responder",
          "Responder"
        )
      )
    ) %>%
    group_by(
      JOINT_GROUP
    ) %>%
    mutate(
      TOTAL = sum(n),
      PROPORTION = ifelse(
        TOTAL > 0,
        n / TOTAL,
        0
      ),
      LABEL = paste0(
        n,
        "/",
        TOTAL
      )
    ) %>%
    ungroup()

  joint_plot <- ggplot(
    joint_summary,
    aes(
      JOINT_GROUP,
      PROPORTION,
      fill = RESPONSE_BINARY
    )
  ) +
    geom_col() +
    geom_text(
      aes(label = LABEL),
      position = position_stack(
        vjust = 0.5
      ),
      size = 3.5
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      )
    ) +
    labs(
      title = paste0(
        cohort_name,
        ": joint STEP_UP-CAF3 response stratification"
      ),
      x = NULL,
      y = "Proportion",
      fill = NULL
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      axis.text.x = element_text(
        angle = 25,
        hjust = 1
      )
    )

  joint_dir <- file.path(
    output_directory,
    method,
    "Deep_supplementary",
    cohort_name
  )

  save_ggplot(
    joint_plot,
    joint_dir,
    "Joint_STEP_CAF3_response",
    width = 8.5,
    height = 6.2
  )

  write.csv(
    joint_summary,
    file.path(
      joint_dir,
      "Joint_STEP_CAF3_response_counts.csv"
    ),
    row.names = FALSE
  )
}

root_directory <- normalizePath(
  ROOT_DIR,
  mustWork = TRUE
)

output_directory <- file.path(
  root_directory,
  OUTPUT_FOLDER
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

gmt_file <- GMT_FILE

imvigor_expression_file <- IMVIGOR_EXPRESSION_FILE
imvigor_sample_file <- IMVIGOR_SAMPLE_CLINICAL_FILE
imvigor_patient_file <- IMVIGOR_PATIENT_CLINICAL_FILE

gse78220_expression_file <- GSE78220_EXPRESSION_FILE
gse78220_metadata_file <- GSE78220_METADATA_FILE

input_manifest <- data.frame(
  INPUT = c(
    "Root directory",
    "Gene-set GMT",
    "IMvigor210 TPM expression",
    "IMvigor210 sample clinical",
    "IMvigor210 patient clinical",
    "GSE78220 FPKM expression",
    "GSE78220 GEO metadata"
  ),
  PATH = c(
    root_directory,
    gmt_file,
    imvigor_expression_file,
    imvigor_sample_file,
    imvigor_patient_file,
    gse78220_expression_file,
    gse78220_metadata_file
  ),
  EXISTS = c(
    dir.exists(root_directory),
    file.exists(gmt_file),
    file.exists(imvigor_expression_file),
    file.exists(imvigor_sample_file),
    file.exists(imvigor_patient_file),
    file.exists(gse78220_expression_file),
    file.exists(gse78220_metadata_file)
  ),
  stringsAsFactors = FALSE
)

print(input_manifest, row.names = FALSE)

write.csv(
  input_manifest,
  file.path(
    output_directory,
    "00_EXPLICIT_INPUT_FILE_MANIFEST.csv"
  ),
  row.names = FALSE
)

missing_inputs <- input_manifest %>%
  filter(!EXISTS)

if (nrow(missing_inputs) > 0) {
  stop(
    "The following explicitly specified input paths do not exist:\n",
    paste(
      paste0(
        missing_inputs$INPUT,
        ": ",
        missing_inputs$PATH
      ),
      collapse = "\n"
    ),
    "\n\nCorrect ROOT_DIR or the corresponding fixed filename at the top of the script."
  )
}

message(
  "\nV5.4 is reading ONLY these fixed files:\n",
  paste(
    paste0(
      input_manifest$INPUT,
      " -> ",
      input_manifest$PATH
    ),
    collapse = "\n"
  ),
  "\n"
)

gene_sets <- read_gmt_simple(
  gmt_file
)

step_genes <- select_gene_set(
  gene_sets,
  "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN",
  "STEP.*UP.*84|STEP.*UP.*LEADING"
)

caf3_genes <- select_gene_set(
  gene_sets,
  "CAF_3gene_score",
  "CAF.*3"
)

gene_set_input_summary <- data.frame(
  GENE_SET = c(
    "STEP_UP",
    "CAF3"
  ),
  N_GENES_IN_GMT = c(
    length(step_genes),
    length(caf3_genes)
  ),
  GENES = c(
    paste(
      step_genes,
      collapse = ";"
    ),
    paste(
      caf3_genes,
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

print(
  gene_set_input_summary[
    ,
    c(
      "GENE_SET",
      "N_GENES_IN_GMT"
    )
  ],
  row.names = FALSE
)

write.csv(
  gene_set_input_summary,
  file.path(
    output_directory,
    "00c_GENE_SET_INPUT_SUMMARY.csv"
  ),
  row.names = FALSE
)

if (
  length(step_genes) < 5
) {
  stop(
    "STEP_UP was found in the GMT, but fewer than five genes were read."
  )
}

if (
  length(caf3_genes) < 2
) {
  stop(
    "CAF3 was found in the GMT, but fewer than two genes were read."
  )
}

message(
  "\nLoading IMvigor210 expression..."
)

imvigor_expression_info <- load_imvigor_expression(
  imvigor_expression_file
)

message(
  "Loading IMvigor210 clinical data..."
)

imvigor_clinical <- load_imvigor_clinical(
  imvigor_sample_file,
  imvigor_patient_file
)

message(
  "Loading GSE78220 expression..."
)

gse78220_expression_info <- load_gse78220_expression(
  gse78220_expression_file
)

message(
  "Loading GSE78220 metadata..."
)

gse78220_metadata <- parse_series_matrix(
  gse78220_metadata_file
) %>%
  prepare_gse78220_metadata()

input_dimension_summary <- data.frame(
  DATA_OBJECT = c(
    "IMvigor210 expression",
    "IMvigor210 merged clinical",
    "GSE78220 expression",
    "GSE78220 metadata"
  ),
  N_ROWS_OR_GENES = c(
    nrow(imvigor_expression_info$matrix),
    nrow(imvigor_clinical),
    nrow(gse78220_expression_info$matrix),
    nrow(gse78220_metadata)
  ),
  N_COLUMNS_OR_SAMPLES = c(
    ncol(imvigor_expression_info$matrix),
    ncol(imvigor_clinical),
    ncol(gse78220_expression_info$matrix),
    ncol(gse78220_metadata)
  ),
  stringsAsFactors = FALSE
)

print(input_dimension_summary, row.names = FALSE)

write.csv(
  input_dimension_summary,
  file.path(
    output_directory,
    "00b_INPUT_DIMENSION_SUMMARY.csv"
  ),
  row.names = FALSE
)

all_six_panel_results <- list()
all_response_deep <- list()

for (method in METHODS_TO_RUN) {
  message(
    "\n========================================\n",
    "Scoring method: ",
    method,
    "\n========================================"
  )

  method_directory <- file.path(
    output_directory,
    method
  )

  dir.create(
    method_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  imvigor_scoring <- calculate_scores(
    imvigor_expression_info$matrix,
    method,
    step_genes,
    caf3_genes
  )

  imvigor_scores <- imvigor_scoring$scores %>%
    mutate(
      SAMPLE_ID = normalize_sample_id(
        SAMPLE_ID
      ),
      SAMPLE_KEY = SAMPLE_ID
    ) %>%
    left_join(
      imvigor_clinical,
      by = "SAMPLE_KEY",
      suffix = c(
        "_expression",
        "_clinical"
      )
    ) %>%
    mutate(
      SAMPLE_ID = coalesce(
        SAMPLE_ID_expression,
        SAMPLE_ID_clinical,
        SAMPLE_KEY
      )
    )

  imvigor_id_audit <- bind_rows(
    data.frame(
      SOURCE = "Expression",
      RAW_ID = imvigor_expression_info$original_sample_ids,
      NORMALIZED_ID = imvigor_expression_info$normalized_sample_ids,
      stringsAsFactors = FALSE
    ),
    data.frame(
      SOURCE = "Clinical",
      RAW_ID = as.character(
        imvigor_clinical$SAMPLE_ID
      ),
      NORMALIZED_ID = imvigor_clinical$SAMPLE_KEY,
      stringsAsFactors = FALSE
    )
  )

  write.csv(
    imvigor_id_audit,
    file.path(
      method_directory,
      "IMvigor210_SAMPLE_ID_AUDIT.csv"
    ),
    row.names = FALSE
  )

  imvigor_overlap <- intersect(
    imvigor_expression_info$normalized_sample_ids,
    imvigor_clinical$SAMPLE_KEY
  )

  write.csv(
    data.frame(
      MATCHED_SAMPLE_ID = imvigor_overlap,
      stringsAsFactors = FALSE
    ),
    file.path(
      method_directory,
      "IMvigor210_MATCHED_SAMPLE_IDS.csv"
    ),
    row.names = FALSE
  )

  imvigor_match_qc <- data.frame(
    COHORT = "IMvigor210",
    METHOD = method,
    N_EXPRESSION_SAMPLES = nrow(imvigor_scoring$scores),
    N_MATCHED_TO_CLINICAL = sum(
      !is.na(
        imvigor_scores$PATIENT_ID
      )
    ),
    N_WITH_RESPONSE = sum(
      !is.na(imvigor_scores$RESPONSE_BINARY)
    ),
    N_WITH_USABLE_OS = sum(
      is.finite(imvigor_scores$OS_YEARS) &
        imvigor_scores$OS_EVENT %in% c(0, 1)
    ),
    stringsAsFactors = FALSE
  )

  print(imvigor_match_qc, row.names = FALSE)

  write.csv(
    imvigor_match_qc,
    file.path(
      method_directory,
      "IMvigor210_MATCH_QC.csv"
    ),
    row.names = FALSE
  )

  if (length(imvigor_overlap) == 0) {
    stop(
      "IMvigor210 ID overlap is zero even after normalization. ",
      "Inspect IMvigor210_SAMPLE_ID_AUDIT.csv."
    )
  }

  if (
    imvigor_match_qc$N_MATCHED_TO_CLINICAL <
      0.80 * imvigor_match_qc$N_EXPRESSION_SAMPLES
  ) {
    stop(
      "IMvigor210 clinical matching is unexpectedly low: ",
      imvigor_match_qc$N_MATCHED_TO_CLINICAL,
      "/",
      imvigor_match_qc$N_EXPRESSION_SAMPLES,
      ". Check SAMPLE_ID columns before continuing."
    )
  }

  write.csv(
    imvigor_scores,
    file.path(
      method_directory,
      "IMvigor210_scores_clinical.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    transform(
      imvigor_scoring$coverage,
      COHORT = "IMvigor210",
      METHOD = method
    ),
    file.path(
      method_directory,
      "IMvigor210_gene_set_coverage.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    transform(
      imvigor_scoring$adjustment,
      COHORT = "IMvigor210",
      METHOD = method
    ),
    file.path(
      method_directory,
      "IMvigor210_CAF3_adjustment.csv"
    ),
    row.names = FALSE
  )

  gse78220_scoring <- calculate_scores(
    gse78220_expression_info$matrix,
    method,
    step_genes,
    caf3_genes
  )

  gse78220_sample_scores <- gse78220_scoring$scores %>%
    mutate(
      EXPR_SAMPLE = SAMPLE_ID,
      SAMPLE_TITLE = sub(
        "\\..*$",
        "",
        SAMPLE_ID
      )
    ) %>%
    left_join(
      gse78220_metadata,
      by = "SAMPLE_TITLE"
    )

  if (
    GSE78220_BASELINE_ONLY &&
    "BIOPSY_TIME_V5" %in%
      names(gse78220_sample_scores) &&
    any(
      !is.na(
        gse78220_sample_scores$BIOPSY_TIME_V5
      )
    )
  ) {
    gse78220_sample_scores <-
      gse78220_sample_scores %>%
      filter(
        grepl(
          "pre|baseline",
          BIOPSY_TIME_V5,
          ignore.case = TRUE
        )
      )
  }

  gse78220_match_qc <- data.frame(
    COHORT = "GSE78220",
    METHOD = method,
    N_EXPRESSION_SAMPLES = nrow(gse78220_scoring$scores),
    N_MATCHED_TO_METADATA = sum(
      !is.na(gse78220_sample_scores$PATIENT_ID)
    ),
    N_WITH_RESPONSE = sum(
      !is.na(gse78220_sample_scores$RESPONSE_BINARY)
    ),
    N_WITH_USABLE_OS = sum(
      is.finite(gse78220_sample_scores$OS_YEARS) &
        gse78220_sample_scores$OS_EVENT %in% c(0, 1)
    ),
    stringsAsFactors = FALSE
  )

  print(gse78220_match_qc, row.names = FALSE)

  write.csv(
    gse78220_match_qc,
    file.path(
      method_directory,
      "GSE78220_MATCH_QC.csv"
    ),
    row.names = FALSE
  )

  if (
    gse78220_match_qc$N_MATCHED_TO_METADATA == 0
  ) {
    stop(
      "GSE78220 expression samples did not match GEO metadata. ",
      "Inspect GSE78220_sample_level_scores.csv and the sample-name conversion."
    )
  }

  gse78220_patient_scores <-
    aggregate_gse78220_patient(
      gse78220_sample_scores
    )

  if (
    nrow(gse78220_patient_scores) < 5
  ) {
    stop(
      "Fewer than five GSE78220 patients remained after sample matching/aggregation."
    )
  }

  write.csv(
    gse78220_sample_scores,
    file.path(
      method_directory,
      "GSE78220_sample_level_scores.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    gse78220_patient_scores,
    file.path(
      method_directory,
      "GSE78220_patient_level_scores.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    transform(
      gse78220_scoring$coverage,
      COHORT = "GSE78220",
      METHOD = method
    ),
    file.path(
      method_directory,
      "GSE78220_gene_set_coverage.csv"
    ),
    row.names = FALSE
  )

  score_class_audit <- bind_rows(
    data.frame(
      COHORT = "IMvigor210",
      METHOD = method,
      SCORE = SCORES_TO_RUN,
      CLASS = vapply(
        imvigor_scores[
          ,
          SCORES_TO_RUN,
          drop = FALSE
        ],
        function(x) paste(class(x), collapse = ";"),
        character(1)
      ),
      N_FINITE = vapply(
        imvigor_scores[
          ,
          SCORES_TO_RUN,
          drop = FALSE
        ],
        function(x) sum(is.finite(suppressWarnings(as.numeric(x)))),
        integer(1)
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      COHORT = "GSE78220",
      METHOD = method,
      SCORE = SCORES_TO_RUN,
      CLASS = vapply(
        gse78220_patient_scores[
          ,
          SCORES_TO_RUN,
          drop = FALSE
        ],
        function(x) paste(class(x), collapse = ";"),
        character(1)
      ),
      N_FINITE = vapply(
        gse78220_patient_scores[
          ,
          SCORES_TO_RUN,
          drop = FALSE
        ],
        function(x) sum(is.finite(suppressWarnings(as.numeric(x)))),
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
  )

  write.csv(
    score_class_audit,
    file.path(
      method_directory,
      "SCORE_NUMERIC_CLASS_AUDIT.csv"
    ),
    row.names = FALSE
  )

  for (score_name in SCORES_TO_RUN) {
    result_key <- paste(
      method,
      score_name,
      sep = "__"
    )

    all_six_panel_results[[result_key]] <-
      make_six_panel_figure(
        imvigor_scores,
        gse78220_patient_scores,
        method,
        score_name,
        output_directory
      )
  }

  imvigor_response_summary <-
    run_three_score_response_summary(
      imvigor_scores,
      "IMvigor210",
      method,
      output_directory
    )

  gse78220_response_summary <-
    run_three_score_response_summary(
      gse78220_patient_scores,
      "GSE78220",
      method,
      output_directory
    )

  all_response_deep[[
    paste0(
      method,
      "__IMvigor210"
    )
  ]] <- imvigor_response_summary

  all_response_deep[[
    paste0(
      method,
      "__GSE78220"
    )
  ]] <- gse78220_response_summary

  run_joint_step_caf_analysis(
    imvigor_scores,
    "IMvigor210",
    method,
    output_directory
  )

  run_joint_step_caf_analysis(
    gse78220_patient_scores,
    "GSE78220",
    method,
    output_directory
  )

  if (RUN_CLINICAL_BENEFIT_SENSITIVITY) {
    benefit_directory <- file.path(
      output_directory,
      method,
      "Clinical_benefit_sensitivity"
    )

    dir.create(
      benefit_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )

    for (
      cohort_object in
        list(
          IMvigor210 = imvigor_scores,
          GSE78220 = gse78220_patient_scores
        )
    ) {
      invisible(
        cohort_object
      )
    }

    for (
      cohort_name in
        c(
          "IMvigor210",
          "GSE78220"
        )
    ) {
      benefit_data <- if (
        cohort_name ==
          "IMvigor210"
      ) {
        imvigor_scores
      } else {
        gse78220_patient_scores
      }

      benefit_data <- benefit_data %>%
        mutate(
          RESPONSE_BINARY =
            CLINICAL_BENEFIT_GROUP
        )

      valid_classes <- unique(
        na.omit(
          benefit_data$RESPONSE_BINARY
        )
      )

      if (
        length(valid_classes) >= 2
      ) {
        for (
          score_name in
            SCORES_TO_RUN
        ) {
          d <- benefit_data %>%
            filter(
              RESPONSE_BINARY %in% c(
                "Clinical benefit",
                "No clinical benefit"
              ),
              is.finite(
                .data[[score_name]]
              )
            ) %>%
            mutate(
              SCORE_Z = as.numeric(
                scale(
                  .data[[score_name]]
                )
              ),
              Y = ifelse(
                RESPONSE_BINARY ==
                  "Clinical benefit",
                1,
                0
              )
            )

          model <- glm(
            Y ~ SCORE_Z,
            data = d,
            family = binomial()
          )

          coefficient <- summary(
            model
          )$coefficients[
            "SCORE_Z",
            ,
            drop = FALSE
          ]

          sensitivity_output <- data.frame(
            COHORT = cohort_name,
            METHOD = method,
            SCORE = score_name,
            ENDPOINT =
              "Clinical benefit: CR/PR/SD versus PD",
            N = nrow(d),
            OR_PER_SD = exp(
              coefficient[1, "Estimate"]
            ),
            CI_LOW = exp(
              coefficient[1, "Estimate"] -
                1.96 *
                coefficient[
                  1,
                  "Std. Error"
                ]
            ),
            CI_HIGH = exp(
              coefficient[1, "Estimate"] +
                1.96 *
                coefficient[
                  1,
                  "Std. Error"
                ]
            ),
            P_VALUE = coefficient[
              1,
              "Pr(>|z|)"
            ],
            stringsAsFactors = FALSE
          )

          write.table(
            sensitivity_output,
            file.path(
              benefit_directory,
              "Clinical_benefit_logistic_all.csv"
            ),
            sep = ",",
            row.names = FALSE,
            col.names = !file.exists(
              file.path(
                benefit_directory,
                "Clinical_benefit_logistic_all.csv"
              )
            ),
            append = file.exists(
              file.path(
                benefit_directory,
                "Clinical_benefit_logistic_all.csv"
              )
            )
          )
        }
      }
    }
  }

  gc()
}

deep_response_summary <- bind_rows(
  all_response_deep
)

write.csv(
  deep_response_summary,
  file.path(
    output_directory,
    "ALL_DEEP_RESPONSE_SUMMARY.csv"
  ),
  row.names = FALSE
)

main_key <- paste(
  MAIN_METHOD,
  MAIN_SCORE,
  sep = "__"
)

if (
  main_key %in%
    names(
      all_six_panel_results
    )
) {
  main_result <-
    all_six_panel_results[[
      main_key
    ]]

  write.csv(
    main_result$survival,
    file.path(
      output_directory,
      "MAIN_FIGURE_OS_STATISTICS.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    main_result$response,
    file.path(
      output_directory,
      "MAIN_FIGURE_RESPONSE_STATISTICS.csv"
    ),
    row.names = FALSE
  )
}

write.csv(
  data.frame(
    SETTING = c(
      "MAIN_METHOD",
      "MAIN_SCORE",
      "GSE78220_BASELINE_ONLY",
      "AVERAGE_REPEATED_BASELINE_BIOPSIES",
      "IMvigor210_transformation",
      "GSE78220_transformation"
    ),
    VALUE = c(
      MAIN_METHOD,
      MAIN_SCORE,
      GSE78220_BASELINE_ONLY,
      AVERAGE_REPEATED_BASELINE_BIOPSIES,
      imvigor_expression_info$transformation,
      gse78220_expression_info$transformation
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    output_directory,
    "ANALYSIS_SETTINGS.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_directory,
    "sessionInfo.txt"
  )
)

message(
  "\n========================================\n",
  "V5 analysis completed.\n",
  "Main manuscript figure:\n",
  file.path(
    output_directory,
    MAIN_METHOD,
    MAIN_SCORE,
    paste0(
      "Figure_6panel_",
      MAIN_METHOD,
      "_",
      MAIN_SCORE,
      ".pdf"
    )
  ),
  "\n========================================"
)
