# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

base_dir <- "path/to/GSEA_all gene_EM QC"

old_outdir <- file.path(
  base_dir,
  "02_fgsea_EMT_contrasts_finalGMT_results_v6_full_outputs"
)

combined_file <- file.path(
  old_outdir,
  "ALL_datasets_EMT_contrast_fgsea_results_combined.csv"
)

fig_dir <- file.path(
  old_outdir,
  "figures_PATCH_noCluster_noM2_noDose_noDOWN"
)

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dataset_order <- c(
  "200981_time",
  "213753_time",
  "310220_time",
  "173958_M1",
  "234265_skin",
  "311419_PTC"
)

contrast_order <- c(
  "Hybrid_vs_nonHybrid",
  "Epithelial_vs_Hybrid",
  "Mesenchymal_vs_Hybrid",
  "Epithelial_vs_Mesenchymal"
)

excluded_datasets <- c(
  "173958_M2",
  "213753_dose"
)

if (!file.exists(combined_file)) {
  stop(
    "Cannot find:\n",
    combined_file
  )
}

combined <- suppressMessages(
  readr::read_csv(
    combined_file,
    show_col_types = FALSE
  )
) %>%
  mutate(
    dataset_id = as.character(.data$dataset_id),
    contrast = as.character(.data$contrast),
    pathway = as.character(.data$pathway),
    leadingEdge_count = as.numeric(.data$leadingEdge_count)
  )

up_pattern <- stringr::regex(
  "(^|[_[:space:]-])UP",
  ignore_case = TRUE
)

reference_pattern <- stringr::regex(
  "EPITHELIAL|MESENCHYMAL|KOHN",
  ignore_case = TRUE
)

up_df <- combined %>%
  filter(
    !.data$dataset_id %in% excluded_datasets,
    stringr::str_detect(
      .data$pathway,
      up_pattern
    ),
    !stringr::str_detect(
      .data$pathway,
      reference_pattern
    )
  )

if (nrow(up_df) == 0) {
  stop(
    "No UP gene sets were detected.\n",
    "Please inspect pathway names in ALL_datasets_EMT_contrast_fgsea_results_combined.csv."
  )
}

up_pathway_order <- up_df %>%
  distinct(.data$pathway) %>%
  arrange(.data$pathway) %>%
  pull(.data$pathway)

message(
  "\nDetected UP gene sets (n = ",
  length(up_pathway_order),
  "):"
)

message(
  paste(
    up_pathway_order,
    collapse = "\n"
  )
)

readr::write_csv(
  tibble(
    order = seq_along(up_pathway_order),
    pathway = up_pathway_order
  ),
  file.path(
    fig_dir,
    "Fig8_UP_only_gene_sets_detected.csv"
  )
)

dataset_contrast_order <- as.vector(
  outer(
    dataset_order,
    contrast_order,
    paste,
    sep = "__"
  )
)

plot_up_le <- up_df %>%
  mutate(
    dataset_contrast = paste(
      .data$dataset_id,
      .data$contrast,
      sep = "__"
    ),

    dataset_contrast = factor(
      .data$dataset_contrast,
      levels = dataset_contrast_order
    ),

    pathway = factor(
      .data$pathway,
      levels = up_pathway_order
    )
  )

p_up_le <- ggplot(
  plot_up_le,
  aes(
    x = pathway,
    y = dataset_contrast,
    fill = leadingEdge_count
  )
) +
  geom_tile(
    color = "grey85",
    linewidth = 0.25
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#1B7837",
    na.value = "grey90"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      size = 7
    ),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Leading edge gene count heatmap: STEP UP gene sets",
    x = "STEP UP gene set",
    y = "Dataset__Contrast",
    fill = "Leading edge count"
  )

ggsave(
  file.path(
    fig_dir,
    "Fig8_UP_only_leadingEdge_count_heatmap_EMT_contrasts.pdf"
  ),
  p_up_le,
  width = 12,
  height = 10,
  useDingbats = FALSE
)

ggsave(
  file.path(
    fig_dir,
    "Fig8_UP_only_leadingEdge_count_heatmap_EMT_contrasts.png"
  ),
  p_up_le,
  width = 12,
  height = 10,
  dpi = 300
)

print(p_up_le)

cat("\n============================================================\n")
cat("DONE — UP-only leading-edge count heatmap\n")
cat("============================================================\n")

cat(
  "Detected UP gene-set number:",
  length(up_pathway_order),
  "\n"
)

cat(
  "Excluded datasets: 173958_M2, 213753_dose\n"
)

cat(
  "E/M/KOHN reference sets: excluded\n"
)

cat(
  "Clustering: NO\n"
)

cat(
  "\nOutput folder:\n"
)

cat(
  fig_dir,
  "\n"
)
