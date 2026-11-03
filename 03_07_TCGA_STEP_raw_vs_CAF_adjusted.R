# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

result_dir <- getwd()

pan_dir <- file.path(result_dir, "09_pan_cancer_summary")
cor_file <- file.path(pan_dir, "ALL_cancers_target_vs_all_spearman.csv")

out_dir <- file.path(
  result_dir,
  "10_figure_ready_CAF_adjusted",
  "Supplementary_Dumbbell"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cor_file)) {
  stop(
    "Cannot find:\n",
    cor_file,
    "\nPlease set working directory to the ROOT of the completed TCGA_STEP_full_results_... folder."
  )
}

raw_target <- "STEP_UP_HYBRID_LEADING_EDGE_84_MOUSE_TO_HUMAN"
adjusted_target <- "STEP_UP_residual_CAF3"

variable_info <- data.frame(
  variable = c(
    "CAF_3gene_score",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HUMAN_MESENCHYMAL_SCORE_18"
  ),
  short_name = c(
    "CAF3",
    "HallmarkEMT",
    "MesenchymalScore"
  ),
  plot_title = c(
    "CAF3 score",
    "Hallmark EMT",
    "Mesenchymal score"
  ),
  stringsAsFactors = FALSE
)

cor_all <- read.csv(
  cor_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c("target", "variable", "cancer", "cor")
missing_cols <- setdiff(required_cols, colnames(cor_all))

if (length(missing_cols) > 0) {
  stop(
    "Correlation table is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

build_pair_table <- function(variable_name) {

  raw <- cor_all[
    cor_all$target == raw_target &
      cor_all$variable == variable_name,
    ,
    drop = FALSE
  ]

  adj <- cor_all[
    cor_all$target == adjusted_target &
      cor_all$variable == variable_name,
    ,
    drop = FALSE
  ]

  if (nrow(raw) == 0 || nrow(adj) == 0) return(NULL)

  raw_keep <- raw[, intersect(c("cancer", "cor", "p", "FDR", "n"), colnames(raw)), drop = FALSE]
  adj_keep <- adj[, intersect(c("cancer", "cor", "p", "FDR", "n"), colnames(adj)), drop = FALSE]

  names(raw_keep)[names(raw_keep) == "cor"] <- "raw_rho"
  names(raw_keep)[names(raw_keep) == "p"] <- "raw_p"
  names(raw_keep)[names(raw_keep) == "FDR"] <- "raw_FDR"
  names(raw_keep)[names(raw_keep) == "n"] <- "raw_n"

  names(adj_keep)[names(adj_keep) == "cor"] <- "adjusted_rho"
  names(adj_keep)[names(adj_keep) == "p"] <- "adjusted_p"
  names(adj_keep)[names(adj_keep) == "FDR"] <- "adjusted_FDR"
  names(adj_keep)[names(adj_keep) == "n"] <- "adjusted_n"

  pair <- merge(
    raw_keep,
    adj_keep,
    by = "cancer",
    all = FALSE
  )

  pair <- pair[
    is.finite(pair$raw_rho) &
      is.finite(pair$adjusted_rho),
    ,
    drop = FALSE
  ]

  if (nrow(pair) == 0) return(NULL)

  pair$delta_rho <- pair$adjusted_rho - pair$raw_rho
  pair$abs_change <- abs(pair$delta_rho)

  pair
}

pair_list <- list()

for (i in seq_len(nrow(variable_info))) {
  vv <- variable_info$variable[i]
  pp <- build_pair_table(vv)
  if (!is.null(pp) && nrow(pp) > 0) {
    pair_list[[vv]] <- pp
  }
}

if (length(pair_list) == 0) {
  stop("None of the three requested paired comparisons could be constructed.")
}

stats_list <- list()

for (i in seq_len(nrow(variable_info))) {

  vv <- variable_info$variable[i]
  if (!(vv %in% names(pair_list))) next

  dat <- pair_list[[vv]]

  wt <- tryCatch(
    wilcox.test(
      dat$adjusted_rho,
      dat$raw_rho,
      paired = TRUE,
      exact = FALSE
    ),
    error = function(e) NULL
  )

  p_wilcox <- if (!is.null(wt)) wt$p.value else NA_real_

  stats_list[[vv]] <- data.frame(
    variable = vv,
    short_name = variable_info$short_name[i],
    plot_title = variable_info$plot_title[i],
    n_cancers = nrow(dat),
    median_raw_rho = median(dat$raw_rho, na.rm = TRUE),
    median_adjusted_rho = median(dat$adjusted_rho, na.rm = TRUE),
    median_delta_rho = median(dat$delta_rho, na.rm = TRUE),
    mean_raw_rho = mean(dat$raw_rho, na.rm = TRUE),
    mean_adjusted_rho = mean(dat$adjusted_rho, na.rm = TRUE),
    wilcoxon_p = p_wilcox,
    stringsAsFactors = FALSE
  )
}

stats_all <- do.call(rbind, stats_list)

stats_all$wilcoxon_FDR_BH <- p.adjust(
  stats_all$wilcoxon_p,
  method = "BH"
)

write.csv(
  stats_all,
  file.path(
    out_dir,
    "Supp_RAW_vs_CAF3adjusted_dumbbell_statistics_summary.csv"
  ),
  row.names = FALSE
)

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  format(signif(p, 3), scientific = FALSE)
}

make_dumbbell_plot <- function(
  dat,
  stat_row,
  short_name,
  plot_title
) {

  dat <- dat[
    order(dat$raw_rho, decreasing = FALSE),
    ,
    drop = FALSE
  ]

  dat$cancer <- factor(
    dat$cancer,
    levels = dat$cancer
  )

  subtitle_text <- paste0(
    "Paired Wilcoxon p = ",
    format_p(stat_row$wilcoxon_p),
    "; BH FDR = ",
    format_p(stat_row$wilcoxon_FDR_BH),
    "; n = ",
    stat_row$n_cancers,
    " cancers"
  )

  p <- ggplot(
    dat,
    aes(y = cancer)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      color = "grey55"
    ) +
    geom_segment(
      aes(
        x = raw_rho,
        xend = adjusted_rho,
        yend = cancer
      ),
      linewidth = 0.65,
      color = "grey65"
    ) +
    geom_point(
      aes(
        x = raw_rho,
        color = "Raw STEP_UP"
      ),
      size = 2.4
    ) +
    geom_point(
      aes(
        x = adjusted_rho,
        color = "CAF3-adjusted STEP_UP"
      ),
      size = 2.4
    ) +
    scale_color_manual(
      values = c(
        "Raw STEP_UP" = "#2166AC",
        "CAF3-adjusted STEP_UP" = "#B2182B"
      )
    ) +
    coord_cartesian(
      xlim = c(-1, 1)
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9)
    ) +
    labs(
      title = paste0(
        "Raw vs CAF3-adjusted STEP_UP: ",
        plot_title
      ),
      subtitle = subtitle_text,
      x = "Spearman correlation (rho)",
      y = NULL,
      color = NULL
    )

  write.csv(
    dat,
    file.path(
      out_dir,
      paste0(
        "Supp_",
        short_name,
        "_RAW_vs_adjusted_perCancer_data.csv"
      )
    ),
    row.names = FALSE
  )

  plot_height <- max(
    7.0,
    0.24 * nrow(dat) + 1.8
  )

  file_base <- file.path(
    out_dir,
    paste0(
      "Supp_",
      short_name,
      "_RAW_vs_adjusted_perCancer_dumbbell"
    )
  )

  ggsave(
    paste0(file_base, ".pdf"),
    p,
    width = 7.2,
    height = plot_height,
    limitsize = FALSE
  )

  ggsave(
    paste0(file_base, ".png"),
    p,
    width = 7.2,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )

  print(p)

  invisible(p)
}

plot_list <- list()

for (i in seq_len(nrow(variable_info))) {

  vv <- variable_info$variable[i]
  if (!(vv %in% names(pair_list))) next

  stat_row <- stats_all[
    stats_all$variable == vv,
    ,
    drop = FALSE
  ]

  plot_list[[vv]] <- make_dumbbell_plot(
    dat = pair_list[[vv]],
    stat_row = stat_row,
    short_name = variable_info$short_name[i],
    plot_title = variable_info$plot_title[i]
  )
}

combined_pairs <- do.call(
  rbind,
  lapply(
    names(pair_list),
    function(vv) {
      x <- pair_list[[vv]]
      x$variable <- vv
      x$plot_title <- variable_info$plot_title[
        match(vv, variable_info$variable)
      ]
      x
    }
  )
)

write.csv(
  combined_pairs,
  file.path(
    out_dir,
    "Supp_ALL3_RAW_vs_CAF3adjusted_perCancer_paired_data.csv"
  ),
  row.names = FALSE
)

sink(
  file.path(
    out_dir,
    "sessionInfo_dumbbell_v1_0.txt"
  )
)
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("DONE - Raw vs CAF3-adjusted per-cancer dumbbell plots\n")
cat("============================================================\n")
cat("Output folder:\n", out_dir, "\n", sep = "")
cat("No TCGA expression, scoring, CAF residualization, or correlation analysis was rerun.\n")
