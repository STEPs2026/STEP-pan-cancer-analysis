# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

base_dir <- "path/to/GSEA_all gene_EM QC"

outdir <- "02_fgsea_EMT_contrasts_finalGMT_results_v6_full_outputs"
per_dataset_dir <- file.path(outdir, "per_dataset_results")
fig_dir <- file.path(outdir, "figures")
rank_export_dir <- file.path(outdir, "contrast_ranked_lists")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(per_dataset_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rank_export_dir, showWarnings = FALSE, recursive = TRUE)

minSize_use <- 3
maxSize_use <- 5000
set.seed(123)

dataset_order <- c(
  "200981_time", "213753_time", "213753_dose", "310220_time",
  "173958_M1", "173958_M2", "234265_skin", "311419_PTC"
)

species_order <- c("human", "mouse")

contrast_info <- tibble::tribble(
  ~contrast,                         ~group_A,        ~group_B,
  "Hybrid_vs_nonHybrid",             "Hybrid",        "nonHybrid",
  "Epithelial_vs_Hybrid",            "Epithelial",    "Hybrid",
  "Mesenchymal_vs_Hybrid",           "Mesenchymal",   "Hybrid",
  "Epithelial_vs_Mesenchymal",       "Epithelial",    "Mesenchymal"
) %>%
  mutate(contrast = factor(contrast, levels = contrast))

contrast_order <- as.character(contrast_info$contrast)

dataset_info <- tibble::tribble(
  ~dataset_id,      ~species, ~dataset_dir,                                 ~gmt_file,
  "200981_time",    "human", "human/200981_time_gene expr list for GSEA",   "human/final_EMT_human_gene_sets.gmt",
  "213753_time",    "human", "human/213753_time_gene expr list for GSEA",   "human/final_EMT_human_gene_sets.gmt",
  "213753_dose",    "human", "human/213753_dose_gene expr list for GSEA",   "human/final_EMT_human_gene_sets.gmt",
  "310220_time",    "human", "human/310220_time_gene expr list for GSEA",   "human/final_EMT_human_gene_sets.gmt",
  "173958_M1",      "mouse", "mus/173958_M1_gene expr list for GSEA",       "mus/final_EMT_mouse_gene_sets.gmt",
  "173958_M2",      "mouse", "mus/173958_M2_gene expr list for GSEA",       "mus/final_EMT_mouse_gene_sets.gmt",
  "234265_skin",    "mouse", "mus/234265_skin_gene expr list for GSEA",     "mus/final_EMT_mouse_gene_sets.gmt",
  "311419_PTC",     "mouse", "mus/311419_PTC_gene expr list for GSEA",      "mus/final_EMT_mouse_gene_sets.gmt"
) %>%
  mutate(
    dataset_dir = file.path(base_dir, dataset_dir),
    rank_dir = file.path(dataset_dir, "01_group_gene_lists_ranked_lists"),
    gmt_file = file.path(base_dir, gmt_file)
  )

safe_paste <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  paste(x, collapse = ";")
}

find_file <- function(rank_dir, dataset_id, group, suffix) {
  file.path(rank_dir, paste0(dataset_id, "_", group, suffix))
}

read_ranked_csv_df <- function(file) {
  if (!file.exists(file)) stop("Ranked file not found: ", file)
  df <- suppressMessages(readr::read_csv(file, show_col_types = FALSE))

  if (!"gene" %in% colnames(df)) {
    stop("Ranked file lacks 'gene' column: ", file)
  }

  score_col <- NULL
  if ("avg_expr" %in% colnames(df)) score_col <- "avg_expr"
  if (is.null(score_col) && "score" %in% colnames(df)) score_col <- "score"
  if (is.null(score_col) && "rank_score" %in% colnames(df)) score_col <- "rank_score"

  if (is.null(score_col)) {
    stop("Cannot find score column. Expected avg_expr / score / rank_score in: ", file)
  }

  df2 <- df %>%
    dplyr::mutate(
      gene = as.character(.data$gene),
      rank_score_raw = as.numeric(.data[[score_col]])
    ) %>%
    dplyr::filter(!is.na(.data$gene), .data$gene != "", !is.na(.data$rank_score_raw)) %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::summarise(
      avg_expr = max(.data$rank_score_raw, na.rm = TRUE),
      pct_expr = if ("pct_expr" %in% colnames(df)) max(.data$pct_expr, na.rm = TRUE) else NA_real_,
      cell_number = if ("cell_number" %in% colnames(df)) max(.data$cell_number, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$avg_expr))

  return(df2)
}

read_expressed_csv <- function(file) {
  if (!file.exists(file)) return(character(0))
  df <- suppressMessages(readr::read_csv(file, show_col_types = FALSE))
  if (!"gene" %in% colnames(df)) return(character(0))
  unique(as.character(df$gene[!is.na(df$gene) & df$gene != ""]))
}

make_contrast_rank <- function(rank_file_A, rank_file_B, group_A, group_B) {
  df_A <- read_ranked_csv_df(rank_file_A) %>%
    dplyr::rename(
      avg_expr_A = .data$avg_expr,
      pct_expr_A = .data$pct_expr,
      cell_number_A = .data$cell_number
    )

  df_B <- read_ranked_csv_df(rank_file_B) %>%
    dplyr::rename(
      avg_expr_B = .data$avg_expr,
      pct_expr_B = .data$pct_expr,
      cell_number_B = .data$cell_number
    )

  contrast_df <- dplyr::inner_join(df_A, df_B, by = "gene") %>%
    dplyr::mutate(
      rank_score = .data$avg_expr_A - .data$avg_expr_B,
      abs_rank_score = abs(.data$rank_score),
      pct_expr_max = pmax(.data$pct_expr_A, .data$pct_expr_B, na.rm = TRUE),
      cell_number_max = pmax(.data$cell_number_A, .data$cell_number_B, na.rm = TRUE),
      group_A = group_A,
      group_B = group_B
    ) %>%
    dplyr::filter(!is.na(.data$rank_score), is.finite(.data$rank_score)) %>%
    dplyr::arrange(dplyr::desc(.data$rank_score)) %>%
    dplyr::mutate(rank_position = dplyr::row_number())

  ranks <- contrast_df$rank_score
  names(ranks) <- contrast_df$gene
  ranks <- sort(ranks, decreasing = TRUE)

  list(ranks = ranks, rank_meta = contrast_df)
}

run_one_fgsea_contrast <- function(dataset_id, species, contrast, group_A, group_B,
                                   rank_file_A, rank_file_B,
                                   expressed_file_A, expressed_file_B,
                                   gmt_file) {

  cat("\n------------------------------------------------------------\n")
  cat("Running:", dataset_id, species, as.character(contrast), "\n")
  cat("Group A:", group_A, "\n")
  cat("Group B:", group_B, "\n")
  cat("Ranked file A:", rank_file_A, "\n")
  cat("Ranked file B:", rank_file_B, "\n")
  cat("GMT file:", gmt_file, "\n")

  if (!file.exists(rank_file_A)) stop("Group A ranked file not found: ", rank_file_A)
  if (!file.exists(rank_file_B)) stop("Group B ranked file not found: ", rank_file_B)
  if (!file.exists(gmt_file)) stop("GMT file not found: ", gmt_file)

  pathways <- fgsea::gmtPathways(gmt_file)

  ranked <- make_contrast_rank(
    rank_file_A = rank_file_A,
    rank_file_B = rank_file_B,
    group_A = group_A,
    group_B = group_B
  )
  ranks <- ranked$ranks
  rank_meta <- ranked$rank_meta

  expressed_A <- read_expressed_csv(expressed_file_A)
  expressed_B <- read_expressed_csv(expressed_file_B)
  expressed_union <- unique(c(expressed_A, expressed_B))

  pathway_info <- tibble(
    pathway = names(pathways),
    pathway_gene_number_original = lengths(pathways),
    pathway_gene_number_matched = sapply(pathways, function(x) length(intersect(x, names(ranks)))),
    pathway_match_ratio = pathway_gene_number_matched / pmax(pathway_gene_number_original, 1),
    pathway_expressed_A_overlap_count = sapply(pathways, function(x) length(intersect(x, expressed_A))),
    pathway_expressed_B_overlap_count = sapply(pathways, function(x) length(intersect(x, expressed_B))),
    pathway_expressed_union_overlap_count = sapply(pathways, function(x) length(intersect(x, expressed_union))),
    pathway_expressed_A_overlap_genes = sapply(pathways, function(x) safe_paste(intersect(x, expressed_A))),
    pathway_expressed_B_overlap_genes = sapply(pathways, function(x) safe_paste(intersect(x, expressed_B))),
    pathway_expressed_union_overlap_genes = sapply(pathways, function(x) safe_paste(intersect(x, expressed_union)))
  )

  fg <- fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = minSize_use,
    maxSize = maxSize_use,
    eps = 0
  )

  fg_df <- as.data.frame(fg) %>%
    as_tibble() %>%
    dplyr::mutate(
      dataset_id = dataset_id,
      species = species,
      contrast = as.character(contrast),
      group_A = group_A,
      group_B = group_B,
      direction_interpretation = ifelse(NES >= 0,
                                        paste0("enriched_in_", group_A),
                                        paste0("enriched_in_", group_B)),
      FDR = padj,
      leadingEdge_genes = sapply(.data$leadingEdge, safe_paste),
      leadingEdge_count = sapply(.data$leadingEdge, length),
      rank_metric = paste0("avg_expr_", group_A, "_minus_avg_expr_", group_B),
      ranked_gene_number = length(ranks),
      expressed_gene_number_A = length(expressed_A),
      expressed_gene_number_B = length(expressed_B),
      expressed_gene_number_union = length(expressed_union),
      ranked_file_A = rank_file_A,
      ranked_file_B = rank_file_B,
      expressed_file_A = expressed_file_A,
      expressed_file_B = expressed_file_B,
      gmt_file = gmt_file,
      analysis_date = as.character(Sys.Date())
    ) %>%
    dplyr::left_join(pathway_info, by = "pathway") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      leadingEdge_in_expressed_A = safe_paste(intersect(unlist(.data$leadingEdge), expressed_A)),
      leadingEdge_in_expressed_B = safe_paste(intersect(unlist(.data$leadingEdge), expressed_B)),
      leadingEdge_in_expressed_union = safe_paste(intersect(unlist(.data$leadingEdge), expressed_union)),
      leadingEdge_in_expressed_A_count = length(intersect(unlist(.data$leadingEdge), expressed_A)),
      leadingEdge_in_expressed_B_count = length(intersect(unlist(.data$leadingEdge), expressed_B)),
      leadingEdge_in_expressed_union_count = length(intersect(unlist(.data$leadingEdge), expressed_union))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      dataset_id, species, contrast, group_A, group_B, pathway,
      ES, NES, pval, padj, FDR, size,
      direction_interpretation,
      leadingEdge_count, leadingEdge_genes,
      rank_metric, ranked_gene_number,
      pathway_gene_number_original, pathway_gene_number_matched, pathway_match_ratio,
      expressed_gene_number_A, expressed_gene_number_B, expressed_gene_number_union,
      pathway_expressed_A_overlap_count, pathway_expressed_B_overlap_count, pathway_expressed_union_overlap_count,
      pathway_expressed_A_overlap_genes, pathway_expressed_B_overlap_genes, pathway_expressed_union_overlap_genes,
      leadingEdge_in_expressed_A_count, leadingEdge_in_expressed_B_count, leadingEdge_in_expressed_union_count,
      leadingEdge_in_expressed_A, leadingEdge_in_expressed_B, leadingEdge_in_expressed_union,
      ranked_file_A, ranked_file_B, expressed_file_A, expressed_file_B, gmt_file, analysis_date,
      dplyr::everything()
    ) %>%
    dplyr::arrange(.data$FDR, dplyr::desc(abs(.data$NES)))

  out_one <- file.path(per_dataset_dir, paste0(dataset_id, "_", as.character(contrast), "_fgsea_results.csv"))
  readr::write_csv(fg_df, out_one)

  out_rank_csv <- file.path(rank_export_dir, paste0(dataset_id, "_", as.character(contrast), "_ranked_genes.csv"))
  out_rank_rnk <- file.path(rank_export_dir, paste0(dataset_id, "_", as.character(contrast), "_ranked_genes.rnk"))
  readr::write_csv(rank_meta, out_rank_csv)
  readr::write_tsv(
    rank_meta %>% dplyr::select(.data$gene, .data$rank_score),
    out_rank_rnk,
    col_names = FALSE
  )

  le_long <- fg_df %>%
    dplyr::select(dataset_id, species, contrast, group_A, group_B, pathway, ES, NES, pval, padj, FDR, direction_interpretation, leadingEdge) %>%
    tidyr::unnest_longer(.data$leadingEdge, values_to = "gene", keep_empty = FALSE) %>%
    dplyr::mutate(gene = as.character(.data$gene)) %>%
    dplyr::left_join(
      rank_meta %>%
        dplyr::select(
          .data$gene, .data$rank_position, .data$rank_score,
          .data$avg_expr_A, .data$avg_expr_B,
          .data$pct_expr_A, .data$pct_expr_B,
          .data$cell_number_A, .data$cell_number_B
        ),
      by = "gene"
    ) %>%
    dplyr::mutate(
      is_expressed_A = .data$gene %in% expressed_A,
      is_expressed_B = .data$gene %in% expressed_B,
      is_expressed_union = .data$gene %in% expressed_union
    )

  list(fgsea = fg_df, leading_edge_long = le_long)
}

make_matrix <- function(df, value_col) {
  df %>%
    dplyr::mutate(dataset_contrast = paste(.data$dataset_id, .data$contrast, sep = "__")) %>%
    dplyr::select(.data$dataset_contrast, .data$pathway, value = dplyr::all_of(value_col)) %>%
    dplyr::distinct(.data$dataset_contrast, .data$pathway, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = .data$pathway, values_from = .data$value) %>%
    dplyr::arrange(factor(.data$dataset_contrast, levels = as.vector(outer(dataset_order, contrast_order, paste, sep = "__"))))
}

save_plot_pdf_png <- function(p, filename_base, width = 11, height = 7) {
  pdf_file <- file.path(fig_dir, paste0(filename_base, ".pdf"))
  png_file <- file.path(fig_dir, paste0(filename_base, ".png"))
  ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE)
  ggsave(png_file, p, width = width, height = height, dpi = 300)
  print(p)  
}

all_results <- list()
all_leading <- list()
run_log <- list()

for (i in seq_len(nrow(dataset_info))) {
  di <- dataset_info[i, ]

  for (j in seq_len(nrow(contrast_info))) {
    ci <- contrast_info[j, ]

    rank_file_A <- find_file(
      rank_dir = di$rank_dir,
      dataset_id = di$dataset_id,
      group = ci$group_A,
      suffix = "_avg_expression_ranked_genes.csv"
    )

    rank_file_B <- find_file(
      rank_dir = di$rank_dir,
      dataset_id = di$dataset_id,
      group = ci$group_B,
      suffix = "_avg_expression_ranked_genes.csv"
    )

    expressed_file_A <- find_file(
      rank_dir = di$rank_dir,
      dataset_id = di$dataset_id,
      group = ci$group_A,
      suffix = "_expressed_gene_list.csv"
    )

    expressed_file_B <- find_file(
      rank_dir = di$rank_dir,
      dataset_id = di$dataset_id,
      group = ci$group_B,
      suffix = "_expressed_gene_list.csv"
    )

    key <- paste(di$dataset_id, as.character(ci$contrast), sep = "__")

    res <- tryCatch(
      run_one_fgsea_contrast(
        dataset_id = di$dataset_id,
        species = di$species,
        contrast = ci$contrast,
        group_A = ci$group_A,
        group_B = ci$group_B,
        rank_file_A = rank_file_A,
        rank_file_B = rank_file_B,
        expressed_file_A = expressed_file_A,
        expressed_file_B = expressed_file_B,
        gmt_file = di$gmt_file
      ),
      error = function(e) {
        message("FAILED: ", key, " | ", conditionMessage(e))
        run_log[[key]] <<- tibble(
          dataset_id = di$dataset_id,
          species = di$species,
          contrast = as.character(ci$contrast),
          group_A = ci$group_A,
          group_B = ci$group_B,
          status = "failed",
          error_message = conditionMessage(e),
          rank_file_A = rank_file_A,
          rank_file_B = rank_file_B,
          gmt_file = di$gmt_file
        )
        return(NULL)
      }
    )

    if (!is.null(res)) {
      all_results[[key]] <- res$fgsea
      all_leading[[key]] <- res$leading_edge_long
      run_log[[key]] <- tibble(
        dataset_id = di$dataset_id,
        species = di$species,
        contrast = as.character(ci$contrast),
        group_A = ci$group_A,
        group_B = ci$group_B,
        status = "ok",
        error_message = "",
        rank_file_A = rank_file_A,
        rank_file_B = rank_file_B,
        gmt_file = di$gmt_file
      )
    }
  }
}

run_log_df <- bind_rows(run_log)
write_csv(run_log_df, file.path(outdir, "run_log.csv"))

if (length(all_results) == 0) {
  stop("No fgsea result was generated. Check run_log.csv for file path or GMT errors.")
}

combined <- bind_rows(all_results) %>%
  dplyr::mutate(
    dataset_id = factor(.data$dataset_id, levels = dataset_order),
    species = factor(.data$species, levels = species_order),
    contrast = factor(.data$contrast, levels = contrast_order),
    dataset_contrast = paste(as.character(.data$dataset_id), as.character(.data$contrast), sep = "__")
  ) %>%
  dplyr::arrange(.data$species, .data$dataset_id, .data$contrast, .data$FDR, dplyr::desc(abs(.data$NES)))

leading_long <- bind_rows(all_leading) %>%
  dplyr::mutate(
    dataset_id = factor(.data$dataset_id, levels = dataset_order),
    species = factor(.data$species, levels = species_order),
    contrast = factor(.data$contrast, levels = contrast_order),
    dataset_contrast = paste(as.character(.data$dataset_id), as.character(.data$contrast), sep = "__")
  )

write_csv(combined, file.path(outdir, "ALL_datasets_EMT_contrast_fgsea_results_combined.csv"))

sig_all_positive_A <- combined %>%
  dplyr::filter(.data$FDR < 0.05, .data$NES > 0, .data$size >= minSize_use, .data$leadingEdge_count >= 3) %>%
  dplyr::arrange(.data$species, .data$dataset_id, .data$contrast, .data$FDR, dplyr::desc(.data$NES))

sig_all_negative_B <- combined %>%
  dplyr::filter(.data$FDR < 0.05, .data$NES < 0, .data$size >= minSize_use, .data$leadingEdge_count >= 3) %>%
  dplyr::arrange(.data$species, .data$dataset_id, .data$contrast, .data$FDR, .data$NES)

sig_all_abs <- combined %>%
  dplyr::filter(.data$FDR < 0.05, .data$size >= minSize_use, .data$leadingEdge_count >= 3) %>%
  dplyr::arrange(.data$species, .data$dataset_id, .data$contrast, .data$FDR, dplyr::desc(abs(.data$NES)))

write_csv(sig_all_positive_A, file.path(outdir, "ALL_datasets_EMT_contrast_fgsea_significant_FDR_lt_0.05_NES_gt_0_enriched_in_groupA.csv"))
write_csv(sig_all_negative_B, file.path(outdir, "ALL_datasets_EMT_contrast_fgsea_significant_FDR_lt_0.05_NES_lt_0_enriched_in_groupB.csv"))
write_csv(sig_all_abs, file.path(outdir, "ALL_datasets_EMT_contrast_fgsea_significant_FDR_lt_0.05_absNES.csv"))

for (ct in contrast_order) {
  write_csv(
    sig_all_abs %>% dplyr::filter(.data$contrast == ct),
    file.path(outdir, paste0(ct, "_fgsea_significant_FDR_lt_0.05_absNES.csv"))
  )
}

for (ct in contrast_order) {
  write_csv(
    sig_all_positive_A %>% dplyr::filter(.data$contrast == ct),
    file.path(outdir, paste0(ct, "_enriched_in_groupA_FDR_lt_0.05_NES_gt_0.csv"))
  )
  write_csv(
    sig_all_negative_B %>% dplyr::filter(.data$contrast == ct),
    file.path(outdir, paste0(ct, "_enriched_in_groupB_FDR_lt_0.05_NES_lt_0.csv"))
  )
}

write_csv(
  sig_all_positive_A %>% dplyr::filter(.data$contrast == "Hybrid_vs_nonHybrid"),
  file.path(outdir, "Hybrid_vs_nonHybrid_enriched_in_Hybrid_FDR_lt_0.05_NES_gt_0.csv")
)
write_csv(
  sig_all_negative_B %>% dplyr::filter(.data$contrast == "Hybrid_vs_nonHybrid"),
  file.path(outdir, "Hybrid_vs_nonHybrid_enriched_in_nonHybrid_FDR_lt_0.05_NES_lt_0.csv")
)

em_reference_pattern <- regex("EPITHELIAL|MESENCHYMAL|KOHN", ignore_case = TRUE)
write_csv(
  combined %>% dplyr::filter(stringr::str_detect(.data$pathway, em_reference_pattern)),
  file.path(outdir, "ALL_datasets_EM_reference_gene_sets_results.csv")
)
write_csv(
  combined %>% dplyr::filter(!stringr::str_detect(.data$pathway, em_reference_pattern)),
  file.path(outdir, "ALL_datasets_custom_Hybrid_EMT_gene_sets_results.csv")
)

write_csv(leading_long, file.path(outdir, "ALL_leading_edge_genes_long.csv"))

qc_gene_sets <- combined %>%
  dplyr::filter(stringr::str_detect(.data$pathway, regex("EPITHELIAL|MESENCHYMAL|KOHN", ignore_case = TRUE))) %>%
  dplyr::mutate(
    expected_direction_hint = dplyr::case_when(
      stringr::str_detect(.data$pathway, regex("EPITHELIAL", ignore_case = TRUE)) & .data$contrast == "Epithelial_vs_Hybrid" ~ "NES > 0 expected if Epithelial > Hybrid for E score",
      stringr::str_detect(.data$pathway, regex("MESENCHYMAL", ignore_case = TRUE)) & .data$contrast == "Mesenchymal_vs_Hybrid" ~ "NES > 0 expected if Mesenchymal > Hybrid for M score",
      stringr::str_detect(.data$pathway, regex("EPITHELIAL", ignore_case = TRUE)) & .data$contrast == "Epithelial_vs_Mesenchymal" ~ "NES > 0 expected if Epithelial > Mesenchymal for E score",
      stringr::str_detect(.data$pathway, regex("MESENCHYMAL", ignore_case = TRUE)) & .data$contrast == "Epithelial_vs_Mesenchymal" ~ "NES < 0 expected if Mesenchymal > Epithelial for M score",
      TRUE ~ "check context"
    )
  )

write_csv(qc_gene_sets, file.path(outdir, "QC_EM_gene_sets_fgsea_results.csv"))

leading_recurrence <- leading_long %>%
  dplyr::group_by(.data$species, .data$contrast, .data$gene) %>%
  dplyr::summarise(
    occurrence_n = dplyr::n(),
    dataset_n = dplyr::n_distinct(.data$dataset_id),
    pathway_n = dplyr::n_distinct(.data$pathway),
    datasets = safe_paste(unique(as.character(.data$dataset_id))),
    pathways = safe_paste(unique(.data$pathway)),
    mean_NES = mean(.data$NES, na.rm = TRUE),
    median_NES = median(.data$NES, na.rm = TRUE),
    max_abs_NES = max(abs(.data$NES), na.rm = TRUE),
    min_FDR = min(.data$FDR, na.rm = TRUE),
    median_FDR = median(.data$FDR, na.rm = TRUE),
    mean_rank_score = mean(.data$rank_score, na.rm = TRUE),
    median_rank_score = median(.data$rank_score, na.rm = TRUE),
    expressed_A_occurrence_n = sum(.data$is_expressed_A, na.rm = TRUE),
    expressed_B_occurrence_n = sum(.data$is_expressed_B, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$species, .data$contrast, dplyr::desc(.data$dataset_n), dplyr::desc(.data$occurrence_n), .data$min_FDR)

write_csv(leading_recurrence, file.path(outdir, "leading_edge_gene_recurrence_summary.csv"))

tested_n_by_pathway <- combined %>%
  dplyr::group_by(.data$species, .data$contrast, .data$pathway) %>%
  dplyr::summarise(
    tested_dataset_n = dplyr::n_distinct(.data$dataset_id),
    tested_datasets = safe_paste(sort(unique(as.character(.data$dataset_id)))),
    .groups = "drop"
  )

pathway_summary <- combined %>%
  dplyr::group_by(.data$species, .data$contrast, .data$pathway) %>%
  dplyr::summarise(
    tested_dataset_n = dplyr::n_distinct(.data$dataset_id),
    tested_datasets = safe_paste(sort(unique(as.character(.data$dataset_id)))),
    positive_NES_n = sum(.data$NES > 0, na.rm = TRUE),
    negative_NES_n = sum(.data$NES < 0, na.rm = TRUE),
    significant_positive_NES_n = sum(.data$FDR < 0.05 & .data$NES > 0, na.rm = TRUE),
    significant_negative_NES_n = sum(.data$FDR < 0.05 & .data$NES < 0, na.rm = TRUE),
    mean_NES = mean(.data$NES, na.rm = TRUE),
    median_NES = median(.data$NES, na.rm = TRUE),
    mean_abs_NES = mean(abs(.data$NES), na.rm = TRUE),
    min_FDR = min(.data$FDR, na.rm = TRUE),
    median_FDR = median(.data$FDR, na.rm = TRUE),
    mean_leadingEdge_count = mean(.data$leadingEdge_count, na.rm = TRUE),
    median_leadingEdge_count = median(.data$leadingEdge_count, na.rm = TRUE),
    mean_pathway_match_ratio = mean(.data$pathway_match_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    leading_long %>%
      dplyr::group_by(.data$species, .data$contrast, .data$pathway, .data$gene) %>%
      dplyr::summarise(gene_dataset_n = dplyr::n_distinct(.data$dataset_id), .groups = "drop") %>%
      dplyr::filter(.data$gene_dataset_n >= 2) %>%
      dplyr::group_by(.data$species, .data$contrast, .data$pathway) %>%
      dplyr::summarise(
        recurrent_leading_edge_gene_count = dplyr::n_distinct(.data$gene),
        recurrent_leading_edge_genes = safe_paste(sort(unique(.data$gene))),
        .groups = "drop"
      ),
    by = c("species", "contrast", "pathway")
  ) %>%
  dplyr::mutate(
    recurrent_leading_edge_gene_count = ifelse(is.na(.data$recurrent_leading_edge_gene_count), 0, .data$recurrent_leading_edge_gene_count),
    recurrent_leading_edge_genes = ifelse(is.na(.data$recurrent_leading_edge_genes), "", .data$recurrent_leading_edge_genes)
  ) %>%
  dplyr::arrange(.data$species, .data$contrast, dplyr::desc(.data$significant_positive_NES_n + .data$significant_negative_NES_n), dplyr::desc(.data$mean_abs_NES), .data$min_FDR)

write_csv(pathway_summary, file.path(outdir, "leading_edge_by_pathway_summary.csv"))

le_threshold_dir <- file.path(outdir, "leading_edge_threshold_tables")
dir.create(le_threshold_dir, showWarnings = FALSE, recursive = TRUE)

leading_gene_dataset_detail_by_pathway_ALL <- leading_long %>%
  dplyr::group_by(.data$species, .data$contrast, .data$pathway, .data$gene) %>%
  dplyr::summarise(
    gene_dataset_n = dplyr::n_distinct(.data$dataset_id),
    gene_datasets = safe_paste(sort(unique(as.character(.data$dataset_id)))),
    occurrence_n = dplyr::n(),
    mean_NES = mean(.data$NES, na.rm = TRUE),
    median_NES = median(.data$NES, na.rm = TRUE),
    max_abs_NES = max(abs(.data$NES), na.rm = TRUE),
    min_FDR = min(.data$FDR, na.rm = TRUE),
    median_FDR = median(.data$FDR, na.rm = TRUE),
    mean_rank_score = mean(.data$rank_score, na.rm = TRUE),
    median_rank_score = median(.data$rank_score, na.rm = TRUE),
    expressed_A_occurrence_n = sum(.data$is_expressed_A, na.rm = TRUE),
    expressed_B_occurrence_n = sum(.data$is_expressed_B, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(tested_n_by_pathway, by = c("species", "contrast", "pathway")) %>%
  dplyr::mutate(
    dataset_fraction = .data$gene_dataset_n / .data$tested_dataset_n,
    present_in_at_least_1_dataset = .data$gene_dataset_n >= 1,
    present_in_at_least_2_datasets = .data$gene_dataset_n >= 2,
    present_in_at_least_3_datasets = .data$gene_dataset_n >= 3,
    present_in_all_tested_datasets = .data$gene_dataset_n == .data$tested_dataset_n
  ) %>%
  dplyr::arrange(.data$species, .data$contrast, .data$pathway, dplyr::desc(.data$gene_dataset_n), .data$min_FDR, .data$gene)

write_csv(
  leading_gene_dataset_detail_by_pathway_ALL,
  file.path(le_threshold_dir, "leading_edge_gene_detail_by_pathway_ALL_genes_n_ge_1.csv")
)

make_le_threshold_tables <- function(detail_df, threshold_label, filter_fun, out_dir) {

  gene_detail <- detail_df %>%
    filter_fun() %>%
    dplyr::arrange(.data$species, .data$contrast, .data$pathway, dplyr::desc(.data$gene_dataset_n), .data$min_FDR, .data$gene)

  summary_by_pathway <- gene_detail %>%
    dplyr::group_by(.data$species, .data$contrast, .data$pathway) %>%
    dplyr::summarise(
      tested_dataset_n = max(.data$tested_dataset_n, na.rm = TRUE),
      tested_datasets = safe_paste(sort(unique(.data$tested_datasets))),
      leading_edge_gene_count = dplyr::n_distinct(.data$gene),
      leading_edge_genes = safe_paste(sort(unique(.data$gene))),
      mean_gene_dataset_n = mean(.data$gene_dataset_n, na.rm = TRUE),
      median_gene_dataset_n = median(.data$gene_dataset_n, na.rm = TRUE),
      max_gene_dataset_n = max(.data$gene_dataset_n, na.rm = TRUE),
      mean_dataset_fraction = mean(.data$dataset_fraction, na.rm = TRUE),
      median_dataset_fraction = median(.data$dataset_fraction, na.rm = TRUE),
      mean_NES_across_genes = mean(.data$mean_NES, na.rm = TRUE),
      median_NES_across_genes = median(.data$median_NES, na.rm = TRUE),
      min_FDR_across_genes = min(.data$min_FDR, na.rm = TRUE),
      median_FDR_across_genes = median(.data$median_FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$species, .data$contrast, .data$pathway)

  write_csv(
    gene_detail,
    file.path(out_dir, paste0("leading_edge_gene_detail_by_pathway_", threshold_label, ".csv"))
  )

  write_csv(
    summary_by_pathway,
    file.path(out_dir, paste0("leading_edge_by_pathway_summary_", threshold_label, ".csv"))
  )

  return(list(gene_detail = gene_detail, summary_by_pathway = summary_by_pathway))
}

le_n_ge_1 <- make_le_threshold_tables(
  leading_gene_dataset_detail_by_pathway_ALL,
  "n_ge_1",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 1),
  le_threshold_dir
)

le_n_ge_2 <- make_le_threshold_tables(
  leading_gene_dataset_detail_by_pathway_ALL,
  "n_ge_2",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 2),
  le_threshold_dir
)

le_n_ge_3 <- make_le_threshold_tables(
  leading_gene_dataset_detail_by_pathway_ALL,
  "n_ge_3",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 3),
  le_threshold_dir
)

le_n_eq_tested <- make_le_threshold_tables(
  leading_gene_dataset_detail_by_pathway_ALL,
  "n_eq_tested",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n == .data$tested_dataset_n),
  le_threshold_dir
)

message("\nBuilding cross-species merged leading-edge summaries...")

cross_species_dir <- file.path(le_threshold_dir, "cross_species_merged_human_mouse")
dir.create(cross_species_dir, showWarnings = FALSE, recursive = TRUE)

normalize_pathway_for_cross_species <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "^HUMAN_", "")
  x <- stringr::str_replace(x, "^MOUSE_", "")
  x <- stringr::str_replace(x, "_MOUSE_ORTHOLOGS?$", "")
  x <- stringr::str_replace(x, "_MOUSE_ORTHOLOG_STYLE_GENE_SYMBOLS$", "")
  x <- stringr::str_replace(x, "_HUMAN_ORTHOLOGS?$", "")
  x
}

leading_long_cross <- leading_long %>%
  dplyr::mutate(
    dataset_global_id = paste(.data$species, .data$dataset_id, sep = "__"),
    pathway_original = .data$pathway,
    pathway_unified = normalize_pathway_for_cross_species(.data$pathway),
    gene_original = .data$gene,
    gene_unified = dplyr::if_else(.data$species == "human", as.character(.data$gene), NA_character_)
  )

mouse_genes_to_map <- leading_long_cross %>%
  dplyr::filter(.data$species == "mouse") %>%
  dplyr::pull(.data$gene_original) %>%
  unique() %>%
  as.character()

mouse_to_human_le_map <- tibble::tibble(
  mouse_gene = character(),
  human_gene = character()
)

if (length(mouse_genes_to_map) > 0) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    install.packages("babelgene")
  }
  message("Converting mouse leading-edge genes to human ortholog symbols using babelgene...")
  mouse_to_human_raw <- tryCatch(
    babelgene::orthologs(
      genes = mouse_genes_to_map,
      species = "mouse",
      human = FALSE
    ),
    error = function(e) {
      warning("babelgene mouse->human conversion failed. Mouse genes will keep simple upper-case fallback. Error: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(mouse_to_human_raw) && nrow(mouse_to_human_raw) > 0) {
    
    cn <- colnames(mouse_to_human_raw)
    mouse_col <- dplyr::case_when(
      "symbol" %in% cn ~ "symbol",
      "mouse_symbol" %in% cn ~ "mouse_symbol",
      "input" %in% cn ~ "input",
      TRUE ~ cn[1]
    )
    human_col <- dplyr::case_when(
      "human_symbol" %in% cn ~ "human_symbol",
      "human_symbol" %in% cn ~ "human_symbol",
      "human" %in% cn ~ "human",
      TRUE ~ cn[length(cn)]
    )

    mouse_to_human_le_map <- mouse_to_human_raw %>%
      dplyr::transmute(
        mouse_gene = as.character(.data[[mouse_col]]),
        human_gene = as.character(.data[[human_col]])
      ) %>%
      dplyr::filter(!is.na(.data$mouse_gene), !is.na(.data$human_gene), .data$human_gene != "") %>%
      dplyr::distinct(.data$mouse_gene, .data$human_gene)

    mouse_to_human_le_map <- mouse_to_human_le_map %>%
      dplyr::arrange(.data$mouse_gene, .data$human_gene) %>%
      dplyr::group_by(.data$mouse_gene) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }
}

write_csv(mouse_to_human_le_map, file.path(cross_species_dir, "mouse_leading_edge_gene_to_human_symbol_mapping.csv"))

leading_long_cross <- leading_long_cross %>%
  dplyr::left_join(mouse_to_human_le_map, by = c("gene_original" = "mouse_gene")) %>%
  dplyr::mutate(
    gene_unified = dplyr::case_when(
      .data$species == "human" ~ as.character(.data$gene_original),
      .data$species == "mouse" & !is.na(.data$human_gene) & .data$human_gene != "" ~ as.character(.data$human_gene),
      .data$species == "mouse" ~ toupper(as.character(.data$gene_original)),
      TRUE ~ as.character(.data$gene_original)
    ),
    mouse_gene_unmapped_to_human = .data$species == "mouse" & (is.na(.data$human_gene) | .data$human_gene == "")
  ) %>%
  dplyr::select(-.data$human_gene)

write_csv(leading_long_cross, file.path(cross_species_dir, "ALL_leading_edge_genes_long_cross_species_human_symbol.csv"))

tested_n_cross_by_pathway <- combined %>%
  dplyr::mutate(
    dataset_global_id = paste(.data$species, .data$dataset_id, sep = "__"),
    pathway_original = .data$pathway,
    pathway_unified = normalize_pathway_for_cross_species(.data$pathway)
  ) %>%
  dplyr::group_by(.data$contrast, .data$pathway_unified) %>%
  dplyr::summarise(
    tested_dataset_n = dplyr::n_distinct(.data$dataset_global_id),
    tested_human_dataset_n = dplyr::n_distinct(.data$dataset_global_id[.data$species == "human"]),
    tested_mouse_dataset_n = dplyr::n_distinct(.data$dataset_global_id[.data$species == "mouse"]),
    tested_datasets = safe_paste(sort(unique(as.character(.data$dataset_global_id)))),
    original_pathway_names = safe_paste(sort(unique(as.character(.data$pathway_original)))),
    .groups = "drop"
  )

write_csv(tested_n_cross_by_pathway, file.path(cross_species_dir, "cross_species_tested_dataset_n_by_pathway.csv"))

leading_gene_detail_cross_ALL <- leading_long_cross %>%
  dplyr::filter(!is.na(.data$gene_unified), .data$gene_unified != "") %>%
  dplyr::group_by(.data$contrast, .data$pathway_unified, .data$gene_unified) %>%
  dplyr::summarise(
    gene_dataset_n = dplyr::n_distinct(.data$dataset_global_id),
    human_dataset_n = dplyr::n_distinct(.data$dataset_global_id[.data$species == "human"]),
    mouse_dataset_n = dplyr::n_distinct(.data$dataset_global_id[.data$species == "mouse"]),
    gene_datasets = safe_paste(sort(unique(as.character(.data$dataset_global_id)))),
    source_species = safe_paste(sort(unique(as.character(.data$species)))),
    mouse_gene_originals = safe_paste(sort(unique(as.character(.data$gene_original[.data$species == "mouse"])))),
    human_gene_originals = safe_paste(sort(unique(as.character(.data$gene_original[.data$species == "human"])))),
    occurrence_n = dplyr::n(),
    mean_NES = mean(.data$NES, na.rm = TRUE),
    median_NES = median(.data$NES, na.rm = TRUE),
    max_abs_NES = max(abs(.data$NES), na.rm = TRUE),
    min_FDR = min(.data$FDR, na.rm = TRUE),
    median_FDR = median(.data$FDR, na.rm = TRUE),
    mean_rank_score = mean(.data$rank_score, na.rm = TRUE),
    median_rank_score = median(.data$rank_score, na.rm = TRUE),
    expressed_A_occurrence_n = sum(.data$is_expressed_A, na.rm = TRUE),
    expressed_B_occurrence_n = sum(.data$is_expressed_B, na.rm = TRUE),
    unmapped_mouse_occurrence_n = sum(.data$mouse_gene_unmapped_to_human, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(tested_n_cross_by_pathway, by = c("contrast", "pathway_unified")) %>%
  dplyr::mutate(
    dataset_fraction = .data$gene_dataset_n / .data$tested_dataset_n,
    present_in_at_least_1_dataset = .data$gene_dataset_n >= 1,
    present_in_at_least_2_datasets = .data$gene_dataset_n >= 2,
    present_in_at_least_3_datasets = .data$gene_dataset_n >= 3,
    present_in_at_least_4_datasets = .data$gene_dataset_n >= 4,
    present_in_at_least_5_datasets = .data$gene_dataset_n >= 5,
    present_in_at_least_6_datasets = .data$gene_dataset_n >= 6,
    present_in_at_least_7_datasets = .data$gene_dataset_n >= 7,
    present_in_at_least_8_datasets = .data$gene_dataset_n >= 8,
    present_in_all_tested_datasets = .data$gene_dataset_n == .data$tested_dataset_n
  ) %>%
  dplyr::arrange(.data$contrast, .data$pathway_unified, dplyr::desc(.data$gene_dataset_n), .data$min_FDR, .data$gene_unified)

write_csv(leading_gene_detail_cross_ALL, file.path(cross_species_dir, "cross_species_leading_edge_gene_detail_by_pathway_ALL_genes_n_ge_1.csv"))

make_cross_species_le_threshold_tables <- function(detail_df, threshold_label, filter_fun, out_dir) {
  gene_detail <- detail_df %>%
    filter_fun() %>%
    dplyr::arrange(.data$contrast, .data$pathway_unified, dplyr::desc(.data$gene_dataset_n), .data$min_FDR, .data$gene_unified)

  summary_by_pathway <- gene_detail %>%
    dplyr::group_by(.data$contrast, .data$pathway_unified) %>%
    dplyr::summarise(
      tested_dataset_n = max(.data$tested_dataset_n, na.rm = TRUE),
      tested_human_dataset_n = max(.data$tested_human_dataset_n, na.rm = TRUE),
      tested_mouse_dataset_n = max(.data$tested_mouse_dataset_n, na.rm = TRUE),
      tested_datasets = safe_paste(sort(unique(.data$tested_datasets))),
      original_pathway_names = safe_paste(sort(unique(.data$original_pathway_names))),
      leading_edge_gene_count = dplyr::n_distinct(.data$gene_unified),
      leading_edge_genes_human_symbol = safe_paste(sort(unique(.data$gene_unified))),
      mean_gene_dataset_n = mean(.data$gene_dataset_n, na.rm = TRUE),
      median_gene_dataset_n = median(.data$gene_dataset_n, na.rm = TRUE),
      max_gene_dataset_n = max(.data$gene_dataset_n, na.rm = TRUE),
      mean_dataset_fraction = mean(.data$dataset_fraction, na.rm = TRUE),
      median_dataset_fraction = median(.data$dataset_fraction, na.rm = TRUE),
      genes_seen_in_both_species_n = sum(.data$human_dataset_n > 0 & .data$mouse_dataset_n > 0, na.rm = TRUE),
      genes_seen_in_both_species = safe_paste(sort(unique(.data$gene_unified[.data$human_dataset_n > 0 & .data$mouse_dataset_n > 0]))),
      mean_NES_across_genes = mean(.data$mean_NES, na.rm = TRUE),
      median_NES_across_genes = median(.data$median_NES, na.rm = TRUE),
      min_FDR_across_genes = min(.data$min_FDR, na.rm = TRUE),
      median_FDR_across_genes = median(.data$median_FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$contrast, .data$pathway_unified)

  write_csv(
    gene_detail,
    file.path(out_dir, paste0("cross_species_leading_edge_gene_detail_by_pathway_", threshold_label, ".csv"))
  )

  write_csv(
    summary_by_pathway,
    file.path(out_dir, paste0("cross_species_leading_edge_by_pathway_summary_", threshold_label, ".csv"))
  )

  return(list(gene_detail = gene_detail, summary_by_pathway = summary_by_pathway))
}

cross_le_n_ge_1 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_1",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 1),
  cross_species_dir
)

cross_le_n_ge_2 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_2",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 2),
  cross_species_dir
)

cross_le_n_ge_3 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_3",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 3),
  cross_species_dir
)

cross_le_n_ge_4 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_4",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 4),
  cross_species_dir
)

cross_le_n_ge_5 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_5",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 5),
  cross_species_dir
)

cross_le_n_ge_6 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_6",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 6),
  cross_species_dir
)

cross_le_n_ge_7 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_7",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 7),
  cross_species_dir
)

cross_le_n_ge_8 <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_ge_8",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 8),
  cross_species_dir
)

cross_le_n_eq_tested <- make_cross_species_le_threshold_tables(
  leading_gene_detail_cross_ALL,
  "n_eq_tested",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n == .data$tested_dataset_n),
  cross_species_dir
)

message("Cross-species leading-edge summaries written to: ", cross_species_dir)

leading_gene_global_detail_ALL <- leading_long %>%
  dplyr::group_by(.data$species, .data$contrast, .data$gene) %>%
  dplyr::summarise(
    gene_dataset_n = dplyr::n_distinct(.data$dataset_id),
    gene_datasets = safe_paste(sort(unique(as.character(.data$dataset_id)))),
    pathway_n = dplyr::n_distinct(.data$pathway),
    pathways = safe_paste(sort(unique(.data$pathway))),
    occurrence_n = dplyr::n(),
    mean_NES = mean(.data$NES, na.rm = TRUE),
    median_NES = median(.data$NES, na.rm = TRUE),
    max_abs_NES = max(abs(.data$NES), na.rm = TRUE),
    min_FDR = min(.data$FDR, na.rm = TRUE),
    median_FDR = median(.data$FDR, na.rm = TRUE),
    mean_rank_score = mean(.data$rank_score, na.rm = TRUE),
    median_rank_score = median(.data$rank_score, na.rm = TRUE),
    expressed_A_occurrence_n = sum(.data$is_expressed_A, na.rm = TRUE),
    expressed_B_occurrence_n = sum(.data$is_expressed_B, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    combined %>%
      dplyr::group_by(.data$species, .data$contrast) %>%
      dplyr::summarise(
        tested_dataset_n = dplyr::n_distinct(.data$dataset_id),
        tested_datasets = safe_paste(sort(unique(as.character(.data$dataset_id)))),
        .groups = "drop"
      ),
    by = c("species", "contrast")
  ) %>%
  dplyr::mutate(
    dataset_fraction = .data$gene_dataset_n / .data$tested_dataset_n,
    present_in_at_least_1_dataset = .data$gene_dataset_n >= 1,
    present_in_at_least_2_datasets = .data$gene_dataset_n >= 2,
    present_in_at_least_3_datasets = .data$gene_dataset_n >= 3,
    present_in_all_tested_datasets = .data$gene_dataset_n == .data$tested_dataset_n
  ) %>%
  dplyr::arrange(.data$species, .data$contrast, dplyr::desc(.data$gene_dataset_n), dplyr::desc(.data$pathway_n), .data$min_FDR, .data$gene)

write_csv(
  leading_gene_global_detail_ALL,
  file.path(le_threshold_dir, "leading_edge_gene_global_detail_ALL_genes_n_ge_1.csv")
)

make_global_threshold_table <- function(global_df, threshold_label, filter_fun, out_dir) {
  gene_detail <- global_df %>%
    filter_fun() %>%
    dplyr::arrange(.data$species, .data$contrast, dplyr::desc(.data$gene_dataset_n), dplyr::desc(.data$pathway_n), .data$min_FDR, .data$gene)

  group_summary <- gene_detail %>%
    dplyr::group_by(.data$species, .data$contrast) %>%
    dplyr::summarise(
      tested_dataset_n = max(.data$tested_dataset_n, na.rm = TRUE),
      tested_datasets = safe_paste(sort(unique(.data$tested_datasets))),
      leading_edge_gene_count = dplyr::n_distinct(.data$gene),
      leading_edge_genes = safe_paste(sort(unique(.data$gene))),
      mean_gene_dataset_n = mean(.data$gene_dataset_n, na.rm = TRUE),
      median_gene_dataset_n = median(.data$gene_dataset_n, na.rm = TRUE),
      max_gene_dataset_n = max(.data$gene_dataset_n, na.rm = TRUE),
      mean_pathway_n = mean(.data$pathway_n, na.rm = TRUE),
      max_pathway_n = max(.data$pathway_n, na.rm = TRUE),
      mean_dataset_fraction = mean(.data$dataset_fraction, na.rm = TRUE),
      median_dataset_fraction = median(.data$dataset_fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$species, .data$contrast)

  write_csv(
    gene_detail,
    file.path(out_dir, paste0("leading_edge_gene_global_detail_", threshold_label, ".csv"))
  )

  write_csv(
    group_summary,
    file.path(out_dir, paste0("leading_edge_gene_global_summary_", threshold_label, ".csv"))
  )

  return(list(gene_detail = gene_detail, group_summary = group_summary))
}

global_n_ge_1 <- make_global_threshold_table(
  leading_gene_global_detail_ALL,
  "n_ge_1",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 1),
  le_threshold_dir
)

global_n_ge_2 <- make_global_threshold_table(
  leading_gene_global_detail_ALL,
  "n_ge_2",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 2),
  le_threshold_dir
)

global_n_ge_3 <- make_global_threshold_table(
  leading_gene_global_detail_ALL,
  "n_ge_3",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n >= 3),
  le_threshold_dir
)

global_n_eq_tested <- make_global_threshold_table(
  leading_gene_global_detail_ALL,
  "n_eq_tested",
  function(df) df %>% dplyr::filter(.data$gene_dataset_n == .data$tested_dataset_n),
  le_threshold_dir
)

nes_mat <- make_matrix(combined, "NES")
fdr_mat <- make_matrix(combined, "FDR")
le_count_mat <- make_matrix(combined, "leadingEdge_count")

neglog_combined <- combined %>%
  dplyr::mutate(negLog10FDR = -log10(pmax(.data$FDR, 1e-300)))
neglog_mat <- make_matrix(neglog_combined, "negLog10FDR")

write_csv(nes_mat, file.path(outdir, "NES_matrix_dataset_contrast_by_geneSet.csv"))
write_csv(fdr_mat, file.path(outdir, "FDR_matrix_dataset_contrast_by_geneSet.csv"))
write_csv(neglog_mat, file.path(outdir, "negLog10FDR_matrix_dataset_contrast_by_geneSet.csv"))
write_csv(le_count_mat, file.path(outdir, "leadingEdge_count_matrix_dataset_contrast_by_geneSet.csv"))

pathway_order_df <- combined %>%
  dplyr::distinct(.data$pathway) %>%
  dplyr::mutate(
    pathway_group = dplyr::case_when(
      stringr::str_detect(.data$pathway, stringr::regex("EPITHELIAL", ignore_case = TRUE)) ~ "01_Epithelial_reference",
      stringr::str_detect(.data$pathway, stringr::regex("MESENCHYMAL", ignore_case = TRUE)) ~ "02_Mesenchymal_reference",
      stringr::str_detect(.data$pathway, stringr::regex("^UP", ignore_case = TRUE)) ~ "03_Custom_UP",
      stringr::str_detect(.data$pathway, stringr::regex("^DOWN", ignore_case = TRUE)) ~ "04_Custom_DOWN",
      TRUE ~ "99_Other"
    )
  ) %>%
  dplyr::arrange(.data$pathway_group, .data$pathway)

pathway_order <- pathway_order_df$pathway

readr::write_csv(
  pathway_order_df,
  file.path(outdir, "plot_gene_set_x_axis_order_E_M_UP_DOWN.csv")
)

plot_df <- combined %>%
  dplyr::mutate(
    dataset_contrast = paste(as.character(.data$dataset_id), as.character(.data$contrast), sep = "__"),
    dataset_contrast = factor(.data$dataset_contrast, levels = as.vector(outer(dataset_order, contrast_order, paste, sep = "__"))),
    pathway = factor(.data$pathway, levels = pathway_order),
    negLog10FDR = -log10(pmax(.data$FDR, 1e-300)),
    sig_label = ifelse(.data$FDR < 0.05, "*", "")
  )

p_nes <- ggplot(plot_df, aes(x = pathway, y = dataset_contrast, fill = NES)) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  geom_text(aes(label = sig_label), size = 3.2) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, na.value = "grey90") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  ) +
  labs(
    title = "GSEA NES heatmap: EMT pairwise contrast ranked lists",
    subtitle = "NES > 0 enriched in group_A; NES < 0 enriched in group_B; * indicates FDR < 0.05",
    x = "Gene set",
    y = "Dataset__Contrast",
    fill = "NES"
  )

save_plot_pdf_png(p_nes, "Fig1_NES_heatmap_EMT_contrasts", width = 17, height = 10)

p_fdr <- ggplot(plot_df, aes(x = pathway, y = dataset_contrast, fill = negLog10FDR)) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  geom_text(aes(label = sig_label), size = 3.2) +
  scale_fill_gradient(low = "white", high = "#7B3294", na.value = "grey90") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  ) +
  labs(
    title = "GSEA significance heatmap",
    subtitle = "Value = -log10(FDR); * indicates FDR < 0.05",
    x = "Gene set",
    y = "Dataset__Contrast",
    fill = "-log10(FDR)"
  )

save_plot_pdf_png(p_fdr, "Fig2_negLog10FDR_heatmap_EMT_contrasts", width = 17, height = 10)

p_dot <- ggplot(plot_df, aes(x = pathway, y = dataset_contrast)) +
  geom_point(aes(size = negLog10FDR, color = NES)) +
  scale_color_gradient2(low = "#4575B4", mid = "grey90", high = "#D73027", midpoint = 0, na.value = "grey70") +
  scale_size_continuous(range = c(1.2, 7)) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "GSEA NES/FDR dot plot: EMT contrasts",
    x = "Gene set",
    y = "Dataset__Contrast",
    color = "NES",
    size = "-log10(FDR)"
  )

save_plot_pdf_png(p_dot, "Fig3_NES_FDR_dotplot_EMT_contrasts", width = 17, height = 10)

p_le <- ggplot(plot_df, aes(x = pathway, y = dataset_contrast, fill = leadingEdge_count)) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "#1B7837", na.value = "grey90") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Leading edge gene count heatmap: EMT contrasts",
    x = "Gene set",
    y = "Dataset__Contrast",
    fill = "Leading edge count"
  )

save_plot_pdf_png(p_le, "Fig4_leadingEdge_count_heatmap_EMT_contrasts", width = 17, height = 10)

qc_plot_df <- plot_df %>%
  dplyr::filter(stringr::str_detect(.data$pathway, regex("EPITHELIAL|MESENCHYMAL|KOHN", ignore_case = TRUE)))

if (nrow(qc_plot_df) > 0) {
  p_qc <- ggplot(qc_plot_df, aes(x = pathway, y = dataset_contrast, fill = NES)) +
    geom_tile(color = "grey85", linewidth = 0.25) +
    geom_text(aes(label = sig_label), size = 3.2) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, na.value = "grey90") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      title = "QC heatmap: epithelial/mesenchymal reference gene sets",
      subtitle = "Use this to check whether EMT group annotation is consistent with E/M programs",
      x = "E/M reference gene set",
      y = "Dataset__Contrast",
      fill = "NES"
    )

  save_plot_pdf_png(p_qc, "Fig5_QC_EM_reference_gene_sets_NES_heatmap", width = 11, height = 10)
}

hybrid_recur_top <- leading_recurrence %>%
  dplyr::filter(.data$contrast == "Hybrid_vs_nonHybrid") %>%
  dplyr::arrange(dplyr::desc(.data$dataset_n), dplyr::desc(.data$occurrence_n), .data$min_FDR) %>%
  dplyr::group_by(.data$species) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = reorder(.data$gene, .data$occurrence_n))

if (nrow(hybrid_recur_top) > 0) {
  p_recur <- ggplot(hybrid_recur_top, aes(x = gene, y = occurrence_n)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ species, scales = "free_y") +
    theme_bw(base_size = 11) +
    labs(
      title = "Top recurrent leading-edge genes in Hybrid_vs_nonHybrid contrast",
      x = "Gene",
      y = "Leading edge occurrence count"
    )

  save_plot_pdf_png(p_recur, "Fig6_top_recurrent_leading_edge_genes_Hybrid_vs_nonHybrid", width = 10, height = 7)
}

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone!\n")
cat("Main output directory:", outdir, "\n")
cat("Main combined table:\n")
cat(file.path(outdir, "ALL_datasets_EMT_contrast_fgsea_results_combined.csv"), "\n")
cat("QC table for E/M reference gene sets:\n")
cat(file.path(outdir, "QC_EM_gene_sets_fgsea_results.csv"), "\n")
cat("Run log:\n")
cat(file.path(outdir, "run_log.csv"), "\n")
