# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
library(Seurat)
library(dplyr)
library(ggplot2)

base_dir <- "path/to/GSE234265"

input_rds <- file.path(
  base_dir,
  "control_2_reanalysis_V2",
  "GSE234265_control_2_YFP_tumor_final_EMT_annotated.rds"
)

out_dir <- file.path(
  base_dir,
  "control_2_reanalysis_V2",
  "DEG_Hybrid_EMT_states_1_with_GSEA_full_ranked"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

tumor_final <- readRDS(input_rds)

tumor_final

if (!"EMT_state_final" %in% colnames(tumor_final@meta.data)) {
  stop("EMT_state_final not found in metadata.")
}

table(tumor_final$EMT_state_final)

DefaultAssay(tumor_final) <- "RNA"

Idents(tumor_final) <- "EMT_state_final"

levels(Idents(tumor_final))
table(Idents(tumor_final))

standardize_deg_table <- function(deg_df, comparison_name) {
  
  deg_df <- as.data.frame(deg_df)
  deg_df$gene <- rownames(deg_df)
  rownames(deg_df) <- NULL

  if ("avg_log2FC" %in% colnames(deg_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(deg_df)) {
    fc_col <- "avg_logFC"
  } else {
    stop("No logFC column found in DEG table.")
  }
  
  deg_df$comparison <- comparison_name

  deg_df$regulation_in_target_group <- dplyr::case_when(
    deg_df[[fc_col]] > 0 ~ "Up_in_target_group",
    deg_df[[fc_col]] < 0 ~ "Down_in_target_group",
    TRUE ~ "No_change"
  )

  if ("p_val_adj" %in% colnames(deg_df)) {
    deg_df <- deg_df %>%
      arrange(p_val_adj, desc(abs(.data[[fc_col]])))
  } else {
    deg_df <- deg_df %>%
      arrange(p_val, desc(abs(.data[[fc_col]])))
  }
  
  return(deg_df)
}

export_deg_results <- function(deg_df, prefix, padj_cutoff = 0.05, logfc_cutoff = 0.25) {

  if ("avg_log2FC" %in% colnames(deg_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(deg_df)) {
    fc_col <- "avg_logFC"
  } else {
    stop("No logFC column found in DEG table.")
  }

  write.csv(
    deg_df,
    file = paste0(prefix, "_all_genes.csv"),
    row.names = FALSE
  )

  sig_df <- deg_df %>%
    filter(
      p_val_adj < padj_cutoff,
      abs(.data[[fc_col]]) >= logfc_cutoff
    )
  
  write.csv(
    sig_df,
    file = paste0(prefix, "_significant_genes_padj", padj_cutoff, "_abslogFC", logfc_cutoff, ".csv"),
    row.names = FALSE
  )

  up_df <- sig_df %>%
    filter(.data[[fc_col]] > 0)
  
  write.csv(
    up_df,
    file = paste0(prefix, "_UP_in_target_group.csv"),
    row.names = FALSE
  )

  down_df <- sig_df %>%
    filter(.data[[fc_col]] < 0)
  
  write.csv(
    down_df,
    file = paste0(prefix, "_DOWN_in_target_group.csv"),
    row.names = FALSE
  )

  top_up_df <- up_df %>%
    arrange(p_val_adj, desc(.data[[fc_col]])) %>%
    head(50)
  
  write.csv(
    top_up_df,
    file = paste0(prefix, "_TOP50_UP_in_target_group.csv"),
    row.names = FALSE
  )

  top_down_df <- down_df %>%
    arrange(p_val_adj, .data[[fc_col]]) %>%
    head(50)
  
  write.csv(
    top_down_df,
    file = paste0(prefix, "_TOP50_DOWN_in_target_group.csv"),
    row.names = FALSE
  )

  summary_df <- data.frame(
    comparison = unique(deg_df$comparison),
    total_tested_genes = nrow(deg_df),
    significant_genes = nrow(sig_df),
    up_in_target_group = nrow(up_df),
    down_in_target_group = nrow(down_df)
  )
  
  write.csv(
    summary_df,
    file = paste0(prefix, "_summary.csv"),
    row.names = FALSE
  )
  
  return(
    list(
      all = deg_df,
      significant = sig_df,
      up = up_df,
      down = down_df,
      top_up = top_up_df,
      top_down = top_down_df,
      summary = summary_df
    )
  )
}

deg_early_vs_rest_raw <- FindMarkers(
  object = tumor_final,
  ident.1 = "early_hybrid",
  ident.2 = NULL,              
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  only.pos = FALSE
)

deg_early_vs_rest <- standardize_deg_table(
  deg_early_vs_rest_raw,
  comparison_name = "early_hybrid_vs_all_other_EMT_states"
)

res_early_vs_rest <- export_deg_results(
  deg_df = deg_early_vs_rest,
  prefix = "01_early_hybrid_vs_all_other_EMT_states",
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25
)

res_early_vs_rest$summary
head(res_early_vs_rest$top_up, 20)

deg_late_vs_rest_raw <- FindMarkers(
  object = tumor_final,
  ident.1 = "late_hybrid",
  ident.2 = NULL,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  only.pos = FALSE
)

deg_late_vs_rest <- standardize_deg_table(
  deg_late_vs_rest_raw,
  comparison_name = "late_hybrid_vs_all_other_EMT_states"
)

res_late_vs_rest <- export_deg_results(
  deg_df = deg_late_vs_rest,
  prefix = "02_late_hybrid_vs_all_other_EMT_states",
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25
)

res_late_vs_rest$summary
head(res_late_vs_rest$top_up, 20)

tumor_final$EMT_state_hybrid_merged <- as.character(
  tumor_final$EMT_state_final
)

tumor_final$EMT_state_hybrid_merged[
  tumor_final$EMT_state_hybrid_merged %in% c(
    "early_hybrid",
    "late_hybrid"
  )
] <- "Hybrid"

tumor_final$EMT_state_hybrid_merged[
  tumor_final$EMT_state_hybrid_merged == "epithelial"
] <- "Epithelial"

tumor_final$EMT_state_hybrid_merged[
  tumor_final$EMT_state_hybrid_merged == "full_EMT"
] <- "Full_EMT"

tumor_final$EMT_state_hybrid_merged <- factor(
  tumor_final$EMT_state_hybrid_merged,
  levels = c(
    "Epithelial",
    "Hybrid",
    "Full_EMT"
  )
)

table(tumor_final$EMT_state_hybrid_merged)

Idents(tumor_final) <- "EMT_state_hybrid_merged"

deg_hybrid_vs_rest_raw <- FindMarkers(
  object = tumor_final,
  ident.1 = "Hybrid",
  ident.2 = NULL,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  only.pos = FALSE
)

deg_hybrid_vs_rest <- standardize_deg_table(
  deg_hybrid_vs_rest_raw,
  comparison_name = "merged_Hybrid_vs_all_other_cells"
)

res_hybrid_vs_rest <- export_deg_results(
  deg_df = deg_hybrid_vs_rest,
  prefix = "03_merged_Hybrid_vs_all_other_cells",
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25
)

res_hybrid_vs_rest$summary
head(res_hybrid_vs_rest$top_up, 20)

tumor_final$EMT_state_hybrid_pairwise <- as.character(tumor_final$EMT_state_final)

tumor_final$EMT_state_hybrid_pairwise[
  tumor_final$EMT_state_hybrid_pairwise %in% c("early_hybrid", "late_hybrid", "Hybrid")
] <- "Hybrid"

tumor_final$EMT_state_hybrid_pairwise[
  tumor_final$EMT_state_hybrid_pairwise %in% c("epithelial", "Epithelial")
] <- "Epithelial"

tumor_final$EMT_state_hybrid_pairwise[
  tumor_final$EMT_state_hybrid_pairwise %in% c(
    "full_EMT",
    "Full_EMT",
    "full EMT",
    "full-EMT",
    "mesenchymal",
    "Mesenchymal"
  )
] <- "Mesenchymal"

tumor_final$EMT_state_hybrid_pairwise[
  !tumor_final$EMT_state_hybrid_pairwise %in% c(
    "Hybrid",
    "Epithelial",
    "Mesenchymal"
  )
] <- "Other"

tumor_final$EMT_state_hybrid_pairwise <- factor(
  tumor_final$EMT_state_hybrid_pairwise,
  levels = c(
    "Epithelial",
    "Hybrid",
    "Mesenchymal",
    "Other"
  )
)

table(tumor_final$EMT_state_hybrid_pairwise)

export_deg_results_padj_pairwise <- function(
    deg_df,
    prefix,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25
) {

  if ("avg_log2FC" %in% colnames(deg_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(deg_df)) {
    fc_col <- "avg_logFC"
  } else {
    stop("No logFC column found in DEG table.")
  }

  if (!"p_val_adj" %in% colnames(deg_df)) {
    stop("Column p_val_adj not found in DEG table.")
  }

  write.csv(
    deg_df,
    file = paste0(prefix, "_all_genes.csv"),
    row.names = FALSE
  )

  sig_df <- deg_df %>%
    filter(
      p_val_adj < padj_cutoff,
      abs(.data[[fc_col]]) >= logfc_cutoff
    )

  write.csv(
    sig_df,
    file = paste0(prefix, "_significant_genes_padj", padj_cutoff, "_abslog2FC", logfc_cutoff, ".csv"),
    row.names = FALSE
  )

  up_df <- sig_df %>%
    filter(.data[[fc_col]] > 0)

  write.csv(
    up_df,
    file = paste0(prefix, "_UP_in_Hybrid.csv"),
    row.names = FALSE
  )

  down_df <- sig_df %>%
    filter(.data[[fc_col]] < 0)

  write.csv(
    down_df,
    file = paste0(prefix, "_DOWN_in_Hybrid.csv"),
    row.names = FALSE
  )

  top_up_df <- up_df %>%
    arrange(p_val_adj, desc(.data[[fc_col]])) %>%
    head(50)

  write.csv(
    top_up_df,
    file = paste0(prefix, "_TOP50_UP_in_Hybrid.csv"),
    row.names = FALSE
  )

  top_down_df <- down_df %>%
    arrange(p_val_adj, .data[[fc_col]]) %>%
    head(50)

  write.csv(
    top_down_df,
    file = paste0(prefix, "_TOP50_DOWN_in_Hybrid.csv"),
    row.names = FALSE
  )

  summary_df <- data.frame(
    comparison = unique(deg_df$comparison),
    total_tested_genes = nrow(deg_df),
    significant_genes_padj_0.05_abslog2FC_0.25 = nrow(sig_df),
    up_in_Hybrid = nrow(up_df),
    down_in_Hybrid = nrow(down_df),
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  write.csv(
    summary_df,
    file = paste0(prefix, "_summary.csv"),
    row.names = FALSE
  )

  return(
    list(
      all = deg_df,
      significant = sig_df,
      up = up_df,
      down = down_df,
      top_up = top_up_df,
      top_down = top_down_df,
      summary = summary_df
    )
  )
}

run_hybrid_pairwise_deg <- function(
    seurat_obj,
    target_group,
    reference_group,
    prefix,
    comparison_name,
    group_col = "EMT_state_hybrid_pairwise",
    min_cells_per_group = 3,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25
) {

  group_tab <- table(seurat_obj@meta.data[[group_col]])

  if (!target_group %in% names(group_tab) || group_tab[[target_group]] < min_cells_per_group) {
    warning(
      paste0(
        "Skip ", comparison_name,
        ": target group '", target_group,
        "' has fewer than ", min_cells_per_group, " cells."
      )
    )
    return(NULL)
  }

  if (!reference_group %in% names(group_tab) || group_tab[[reference_group]] < min_cells_per_group) {
    warning(
      paste0(
        "Skip ", comparison_name,
        ": reference group '", reference_group,
        "' has fewer than ", min_cells_per_group, " cells or does not exist."
      )
    )
    return(NULL)
  }

  Idents(seurat_obj) <- group_col

  deg_raw <- FindMarkers(
    object = seurat_obj,
    ident.1 = target_group,
    ident.2 = reference_group,
    test.use = "wilcox",
    logfc.threshold = logfc_cutoff,
    min.pct = 0.10,
    only.pos = FALSE
  )

  deg_tbl <- standardize_deg_table(
    deg_raw,
    comparison_name = comparison_name
  )

  res <- export_deg_results_padj_pairwise(
    deg_df = deg_tbl,
    prefix = prefix,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  return(res)
}

Idents(tumor_final) <- "EMT_state_hybrid_pairwise"

res_hybrid_vs_epithelial <- run_hybrid_pairwise_deg(
  seurat_obj = tumor_final,
  target_group = "Hybrid",
  reference_group = "Epithelial",
  prefix = "04_merged_Hybrid_vs_Epithelial_padj0.05_log2FC0.25",
  comparison_name = "merged_Hybrid_vs_Epithelial",
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25
)

res_hybrid_vs_mesenchymal <- run_hybrid_pairwise_deg(
  seurat_obj = tumor_final,
  target_group = "Hybrid",
  reference_group = "Mesenchymal",
  prefix = "05_merged_Hybrid_vs_Mesenchymal_padj0.05_log2FC0.25",
  comparison_name = "merged_Hybrid_vs_Mesenchymal",
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25
)

res_hybrid_vs_other <- run_hybrid_pairwise_deg(
  seurat_obj = tumor_final,
  target_group = "Hybrid",
  reference_group = "Other",
  prefix = "06_merged_Hybrid_vs_Other_padj0.05_log2FC0.25",
  comparison_name = "merged_Hybrid_vs_Other",
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25
)

pairwise_summary_list <- list(
  res_hybrid_vs_epithelial,
  res_hybrid_vs_mesenchymal,
  res_hybrid_vs_other
)

pairwise_summary_all <- bind_rows(
  lapply(
    pairwise_summary_list,
    function(x) {
      if (is.null(x)) {
        return(NULL)
      } else {
        return(x$summary)
      }
    }
  )
)

write.csv(
  pairwise_summary_all,
  file = "07_DEG_summary_merged_Hybrid_pairwise_padj0.05_log2FC0.25.csv",
  row.names = FALSE
)

pairwise_summary_all

gsea_deg_dir <- file.path(out_dir, "GSEA_full_ranked_DEG_tables")
dir.create(gsea_deg_dir, showWarnings = FALSE, recursive = TRUE)

export_gsea_ranked_deg_table <- function(
    deg_df,
    prefix,
    target_group,
    reference_group,
    comparison_name
) {
  deg_df <- as.data.frame(deg_df)

  if ("avg_log2FC" %in% colnames(deg_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(deg_df)) {
    fc_col <- "avg_logFC"
  } else {
    stop("No logFC column found in DEG table.")
  }

  if (!"p_val_adj" %in% colnames(deg_df)) {
    stop("Column p_val_adj not found in DEG table.")
  }

  deg_df <- deg_df %>%
    mutate(
      target_group = target_group,
      reference_group = reference_group,
      comparison = comparison_name,
      p_val_adj_for_rank = case_when(
        is.na(p_val_adj) ~ 1,
        p_val_adj <= 0 ~ 1e-300,
        TRUE ~ p_val_adj
      ),
      rank_metric = sign(.data[[fc_col]]) * (-log10(p_val_adj_for_rank)),
      rank_direction = paste0(
        "rank_metric > 0 means higher in ", target_group,
        "; rank_metric < 0 means higher in ", reference_group
      )
    ) %>%
    arrange(desc(rank_metric))

  write.csv(
    deg_df,
    file = file.path(gsea_deg_dir, paste0(prefix, "_FULL_all_genes_for_ranked_GSEA.csv")),
    row.names = FALSE
  )

  rnk_df <- deg_df %>%
    select(gene, rank_metric) %>%
    filter(!is.na(gene), !is.na(rank_metric)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(desc(rank_metric))

  write.table(
    rnk_df,
    file = file.path(gsea_deg_dir, paste0(prefix, "_rank_metric.rnk")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  summary_df <- data.frame(
    comparison = comparison_name,
    target_group = target_group,
    reference_group = reference_group,
    total_ranked_genes = nrow(deg_df),
    positive_rank_genes_higher_in_target = sum(deg_df$rank_metric > 0, na.rm = TRUE),
    negative_rank_genes_higher_in_reference = sum(deg_df$rank_metric < 0, na.rm = TRUE),
    zero_rank_genes = sum(deg_df$rank_metric == 0, na.rm = TRUE),
    logfc_threshold_for_FindMarkers = 0,
    min_pct_for_FindMarkers = 0.01
  )

  write.csv(
    summary_df,
    file = file.path(gsea_deg_dir, paste0(prefix, "_GSEA_ranked_summary.csv")),
    row.names = FALSE
  )

  return(list(full = deg_df, rnk = rnk_df, summary = summary_df))
}

run_hybrid_pairwise_deg_for_gsea <- function(
    seurat_obj,
    target_group,
    reference_group,
    prefix,
    comparison_name,
    group_col = "EMT_state_hybrid_pairwise",
    min_cells_per_group = 3
) {
  group_tab <- table(seurat_obj@meta.data[[group_col]])
  print(group_tab)

  if (!target_group %in% names(group_tab) || group_tab[[target_group]] < min_cells_per_group) {
    warning(
      paste0(
        "Skip ", comparison_name,
        ": target group '", target_group,
        "' has fewer than ", min_cells_per_group, " cells."
      )
    )
    return(NULL)
  }

  if (!reference_group %in% names(group_tab) || group_tab[[reference_group]] < min_cells_per_group) {
    warning(
      paste0(
        "Skip ", comparison_name,
        ": reference group '", reference_group,
        "' has fewer than ", min_cells_per_group, " cells or does not exist."
      )
    )
    return(NULL)
  }

  Idents(seurat_obj) <- group_col

  deg_raw <- FindMarkers(
    object = seurat_obj,
    ident.1 = target_group,
    ident.2 = reference_group,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.01,
    only.pos = FALSE
  )

  deg_tbl <- standardize_deg_table(
    deg_raw,
    comparison_name = comparison_name
  )

  res <- export_gsea_ranked_deg_table(
    deg_df = deg_tbl,
    prefix = prefix,
    target_group = target_group,
    reference_group = reference_group,
    comparison_name = comparison_name
  )

  return(res)
}

gsea_hybrid_vs_epithelial <- run_hybrid_pairwise_deg_for_gsea(
  seurat_obj = tumor_final,
  target_group = "Hybrid",
  reference_group = "Epithelial",
  prefix = "04_merged_Hybrid_vs_Epithelial",
  comparison_name = "merged_Hybrid_vs_Epithelial"
)

gsea_hybrid_vs_mesenchymal <- run_hybrid_pairwise_deg_for_gsea(
  seurat_obj = tumor_final,
  target_group = "Hybrid",
  reference_group = "Mesenchymal",
  prefix = "05_merged_Hybrid_vs_Mesenchymal",
  comparison_name = "merged_Hybrid_vs_Mesenchymal"
)

gsea_hybrid_vs_other <- run_hybrid_pairwise_deg_for_gsea(
  seurat_obj = tumor_final,
  target_group = "Hybrid",
  reference_group = "Other",
  prefix = "06_merged_Hybrid_vs_Other",
  comparison_name = "merged_Hybrid_vs_Other"
)

gsea_summary_all <- bind_rows(
  lapply(
    list(
      gsea_hybrid_vs_epithelial,
      gsea_hybrid_vs_mesenchymal,
      gsea_hybrid_vs_other
    ),
    function(x) {
      if (is.null(x)) {
        return(NULL)
      } else {
        return(x$summary)
      }
    }
  )
)

write.csv(
  gsea_summary_all,
  file = file.path(gsea_deg_dir, "00_GSEA_full_ranked_DEG_summary_all_pairwise.csv"),
  row.names = FALSE
)

gsea_summary_all

Idents(tumor_final) <- "EMT_state_final"

deg_summary_all <- bind_rows(
  res_early_vs_rest$summary,
  res_late_vs_rest$summary,
  res_hybrid_vs_rest$summary
)

write.csv(
  deg_summary_all,
  file = "08_DEG_summary_original_comparisons_padj0.05_logFC0.25.csv",
  row.names = FALSE
)

deg_summary_all_with_pairwise <- bind_rows(
  deg_summary_all,
  pairwise_summary_all
)

write.csv(
  deg_summary_all_with_pairwise,
  file = "09_DEG_summary_all_comparisons_including_pairwise.csv",
  row.names = FALSE
)

deg_summary_all

saveRDS(
  tumor_final,
  file = "GSE234265_control_2_YFP_tumor_final_with_HybridMerged_annotation.rds"
)

message("DEG analysis for early_hybrid, late_hybrid, and merged Hybrid finished successfully.")

top_n_heatmap <- 50       
cells_per_group <- 120    

set.seed(1234)

get_top_up_genes <- function(up_deg_df, top_n = 50) {
  
  up_deg_df <- as.data.frame(up_deg_df)

  if ("avg_log2FC" %in% colnames(up_deg_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(up_deg_df)) {
    fc_col <- "avg_logFC"
  } else {
    stop("No logFC column found in DEG table.")
  }
  
  top_genes <- up_deg_df %>%
    arrange(p_val_adj, desc(.data[[fc_col]])) %>%
    head(top_n) %>%
    pull(gene)
  
  return(unique(top_genes))
}

downsample_cells_by_group <- function(seurat_obj, group_col, max_cells_per_group = 120) {
  
  meta_df <- seurat_obj@meta.data
  meta_df$cell_barcode <- rownames(meta_df)
  
  sampled_cells <- meta_df %>%
    group_by(.data[[group_col]]) %>%
    slice_sample(n = min(n(), max_cells_per_group)) %>%
    ungroup() %>%
    pull(cell_barcode)
  
  return(sampled_cells)
}

early_top_genes <- get_top_up_genes(
  up_deg_df = res_early_vs_rest$up,
  top_n = top_n_heatmap
)

early_top_genes <- early_top_genes[
  early_top_genes %in% rownames(tumor_final)
]

early_top_genes

write.csv(
  data.frame(gene = early_top_genes),
  file = "05A_early_hybrid_heatmap_top_genes.csv",
  row.names = FALSE
)

tumor_final$EMT_state_final <- factor(
  tumor_final$EMT_state_final,
  levels = c(
    "epithelial",
    "early_hybrid",
    "late_hybrid",
    "full_EMT"
  )
)

Idents(tumor_final) <- "EMT_state_final"

tumor_final <- ScaleData(
  tumor_final,
  features = early_top_genes,
  assay = "RNA"
)

downsample_cells_by_group <- function(
    seurat_obj,
    group_col,
    max_cells_per_group = 120,
    seed = 1234
) {
  
  set.seed(seed)
  
  meta_df <- seurat_obj@meta.data
  meta_df$cell_barcode <- rownames(meta_df)
  
  if (!group_col %in% colnames(meta_df)) {
    stop(paste0("Group column '", group_col, "' not found in Seurat metadata."))
  }
  
  group_list <- split(
    meta_df$cell_barcode,
    meta_df[[group_col]]
  )
  
  sampled_cells <- unlist(
    lapply(
      group_list,
      function(cells_in_group) {
        
        n_take <- min(
          length(cells_in_group),
          max_cells_per_group
        )
        
        sample(
          cells_in_group,
          size = n_take,
          replace = FALSE
        )
      }
    ),
    use.names = FALSE
  )
  
  return(sampled_cells)
}

early_cells_heatmap <- downsample_cells_by_group(
  seurat_obj = tumor_final,
  group_col = "EMT_state_final",
  max_cells_per_group = cells_per_group
)

p_heatmap_early <- DoHeatmap(
  object = tumor_final,
  features = early_top_genes,
  cells = early_cells_heatmap,
  group.by = "EMT_state_final",
  assay = "RNA",
  slot = "scale.data",
  raster = TRUE,
  size = 3
) +
  ggtitle("Top upregulated genes in early_hybrid vs all other EMT states")

print(p_heatmap_early)

ggsave(
  filename = "05B_Heatmap_early_hybrid_UP_genes.pdf",
  plot = p_heatmap_early,
  width = 10,
  height = 10
)

ggsave(
  filename = "05B_Heatmap_early_hybrid_UP_genes.png",
  plot = p_heatmap_early,
  width = 10,
  height = 10,
  dpi = 300
)

late_top_genes <- get_top_up_genes(
  up_deg_df = res_late_vs_rest$up,
  top_n = top_n_heatmap
)

late_top_genes <- late_top_genes[
  late_top_genes %in% rownames(tumor_final)
]

late_top_genes

write.csv(
  data.frame(gene = late_top_genes),
  file = "06A_late_hybrid_heatmap_top_genes.csv",
  row.names = FALSE
)

Idents(tumor_final) <- "EMT_state_final"

tumor_final <- ScaleData(
  tumor_final,
  features = late_top_genes,
  assay = "RNA"
)

late_cells_heatmap <- downsample_cells_by_group(
  seurat_obj = tumor_final,
  group_col = "EMT_state_final",
  max_cells_per_group = cells_per_group
)

p_heatmap_late <- DoHeatmap(
  object = tumor_final,
  features = late_top_genes,
  cells = late_cells_heatmap,
  group.by = "EMT_state_final",
  assay = "RNA",
  slot = "scale.data",
  raster = TRUE,
  size = 3
) +
  ggtitle("Top upregulated genes in late_hybrid vs all other EMT states")

print(p_heatmap_late)

ggsave(
  filename = "06B_Heatmap_late_hybrid_UP_genes.pdf",
  plot = p_heatmap_late,
  width = 10,
  height = 10
)

ggsave(
  filename = "06B_Heatmap_late_hybrid_UP_genes.png",
  plot = p_heatmap_late,
  width = 10,
  height = 10,
  dpi = 300
)

hybrid_top_genes <- get_top_up_genes(
  up_deg_df = res_hybrid_vs_rest$up,
  top_n = top_n_heatmap
)

hybrid_top_genes <- hybrid_top_genes[
  hybrid_top_genes %in% rownames(tumor_final)
]

hybrid_top_genes

write.csv(
  data.frame(gene = hybrid_top_genes),
  file = "07A_merged_Hybrid_heatmap_top_genes.csv",
  row.names = FALSE
)

if (!"EMT_state_hybrid_merged" %in% colnames(tumor_final@meta.data)) {
  stop("EMT_state_hybrid_merged not found. Please run the merged Hybrid DEG section first.")
}

tumor_final$EMT_state_hybrid_merged <- factor(
  tumor_final$EMT_state_hybrid_merged,
  levels = c(
    "Epithelial",
    "Hybrid",
    "Full_EMT"
  )
)

Idents(tumor_final) <- "EMT_state_hybrid_merged"

tumor_final <- ScaleData(
  tumor_final,
  features = hybrid_top_genes,
  assay = "RNA"
)

hybrid_cells_heatmap <- downsample_cells_by_group(
  seurat_obj = tumor_final,
  group_col = "EMT_state_hybrid_merged",
  max_cells_per_group = cells_per_group
)

p_heatmap_hybrid <- DoHeatmap(
  object = tumor_final,
  features = hybrid_top_genes,
  cells = hybrid_cells_heatmap,
  group.by = "EMT_state_hybrid_merged",
  assay = "RNA",
  slot = "scale.data",
  raster = TRUE,
  size = 3
) +
  ggtitle("Top upregulated genes in merged Hybrid vs all other cells")

print(p_heatmap_hybrid)

ggsave(
  filename = "07B_Heatmap_merged_Hybrid_UP_genes.pdf",
  plot = p_heatmap_hybrid,
  width = 10,
  height = 10
)

ggsave(
  filename = "07B_Heatmap_merged_Hybrid_UP_genes.png",
  plot = p_heatmap_hybrid,
  width = 10,
  height = 10,
  dpi = 300
)

Idents(tumor_final) <- "EMT_state_final"

saveRDS(
  tumor_final,
  file = "GSE234265_control_2_YFP_tumor_final_DEG_heatmap_ready.rds"
)

message("All DEG heatmaps finished successfully.")

library(Seurat)
library(dplyr)
library(Matrix)
library(pheatmap)

top_n_avg_heatmap <- 50   

if (!exists("tumor_final")) {
  stop("Object 'tumor_final' not found. Please load your EMT-annotated Seurat object first.")
}

if (!exists("res_early_vs_rest")) {
  stop("Object 'res_early_vs_rest' not found. Please run the DEG analysis section first.")
}

if (!exists("res_late_vs_rest")) {
  stop("Object 'res_late_vs_rest' not found. Please run the DEG analysis section first.")
}

if (!exists("res_hybrid_vs_rest")) {
  stop("Object 'res_hybrid_vs_rest' not found. Please run the DEG analysis section first.")
}

if (!"EMT_state_final" %in% colnames(tumor_final@meta.data)) {
  stop("Metadata column 'EMT_state_final' not found in tumor_final.")
}

if (!"EMT_state_hybrid_merged" %in% colnames(tumor_final@meta.data)) {
  stop("Metadata column 'EMT_state_hybrid_merged' not found. Please run the merged Hybrid DEG section first.")
}

DefaultAssay(tumor_final) <- "RNA"

tumor_final$EMT_state_final <- factor(
  tumor_final$EMT_state_final,
  levels = c(
    "epithelial",
    "early_hybrid",
    "late_hybrid",
    "full_EMT"
  )
)

tumor_final$EMT_state_hybrid_merged <- factor(
  tumor_final$EMT_state_hybrid_merged,
  levels = c(
    "Epithelial",
    "Hybrid",
    "Full_EMT"
  )
)

get_top_up_genes_from_deg <- function(
    up_deg_df,
    top_n = 30
) {
  
  up_deg_df <- as.data.frame(up_deg_df)
  
  if (nrow(up_deg_df) == 0) {
    stop("The UP DEG table is empty.")
  }
  
  if (!"gene" %in% colnames(up_deg_df)) {
    stop("Column 'gene' not found in the DEG table.")
  }
  
  if ("avg_log2FC" %in% colnames(up_deg_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(up_deg_df)) {
    fc_col <- "avg_logFC"
  } else {
    stop("Neither 'avg_log2FC' nor 'avg_logFC' was found in the DEG table.")
  }
  
  if ("p_val_adj" %in% colnames(up_deg_df)) {
    
    top_genes <- up_deg_df %>%
      arrange(p_val_adj, desc(.data[[fc_col]])) %>%
      head(top_n) %>%
      pull(gene)
    
  } else {
    
    top_genes <- up_deg_df %>%
      arrange(desc(.data[[fc_col]])) %>%
      head(top_n) %>%
      pull(gene)
  }
  
  top_genes <- unique(top_genes)
  
  return(top_genes)
}

get_normalized_data_matrix <- function(
    seurat_obj,
    assay = "RNA"
) {
  
  expr_mat <- tryCatch(
    {
      GetAssayData(
        object = seurat_obj,
        assay = assay,
        layer = "data"
      )
    },
    error = function(e) {
      GetAssayData(
        object = seurat_obj,
        assay = assay,
        slot = "data"
      )
    }
  )
  
  return(expr_mat)
}

make_group_average_expression_matrix <- function(
    seurat_obj,
    features,
    group_col,
    group_order,
    assay = "RNA"
) {
  
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Group column '", group_col, "' not found in metadata."))
  }
  
  expr_mat <- get_normalized_data_matrix(
    seurat_obj = seurat_obj,
    assay = assay
  )
  
  features <- unique(features)
  features <- features[features %in% rownames(expr_mat)]
  
  if (length(features) == 0) {
    stop("None of the selected genes are present in the normalized expression matrix.")
  }
  
  expr_sub <- expr_mat[features, , drop = FALSE]
  
  meta_df <- seurat_obj@meta.data[colnames(expr_sub), , drop = FALSE]
  group_vec <- as.character(meta_df[[group_col]])
  
  avg_list <- lapply(
    group_order,
    function(g) {
      
      cells_g <- colnames(expr_sub)[group_vec == g]
      
      if (length(cells_g) == 0) {
        stop(paste0("No cells found for group: ", g))
      }
      
      Matrix::rowMeans(
        expr_sub[, cells_g, drop = FALSE]
      )
    }
  )
  
  avg_mat <- do.call(cbind, avg_list)
  avg_mat <- as.matrix(avg_mat)
  
  rownames(avg_mat) <- features
  colnames(avg_mat) <- group_order
  
  return(avg_mat)
}

prepare_heatmap_matrix <- function(mat) {
  
  mat <- as.matrix(mat)
  
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop(
      paste0(
        "Heatmap matrix is empty: ",
        nrow(mat), " rows x ", ncol(mat), " columns."
      )
    )
  }

  finite_rows <- apply(
    mat,
    1,
    function(x) all(is.finite(x))
  )
  
  mat <- mat[finite_rows, , drop = FALSE]
  
  if (nrow(mat) == 0) {
    stop("All rows were removed because they contained NA/NaN/Inf.")
  }

  row_sd <- apply(
    mat,
    1,
    sd,
    na.rm = TRUE
  )
  
  variable_rows <- is.finite(row_sd) & row_sd > 0
  
  mat <- mat[variable_rows, , drop = FALSE]
  
  if (nrow(mat) == 0) {
    stop(
      "No rows with non-zero variance remain after filtering. ",
      "The selected genes cannot be used for a row-scaled average heatmap."
    )
  }
  
  return(mat)
}

save_average_expression_heatmap <- function(
    mat,
    prefix,
    title_text,
    width = 7,
    height = 10
) {
  
  mat_for_plot <- prepare_heatmap_matrix(mat)

  write.csv(
    mat_for_plot,
    file = paste0(prefix, "_matrix_used_for_heatmap.csv")
  )

  pdf(
    file = paste0(prefix, ".pdf"),
    width = width,
    height = height
  )
  
  pheatmap::pheatmap(
    mat_for_plot,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 10,
    main = title_text
  )
  
  dev.off()

  png(
    filename = paste0(prefix, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  
  pheatmap::pheatmap(
    mat_for_plot,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 10,
    main = title_text
  )
  
  dev.off()
  
  message(
    "Heatmap saved successfully: ",
    prefix,
    " | plotted genes = ",
    nrow(mat_for_plot)
  )
}

early_top_genes_avg <- get_top_up_genes_from_deg(
  up_deg_df = res_early_vs_rest$up,
  top_n = top_n_avg_heatmap
)

early_top_genes_avg <- early_top_genes_avg[
  early_top_genes_avg %in% rownames(tumor_final)
]

write.csv(
  data.frame(gene = early_top_genes_avg),
  file = "13A_early_hybrid_TOP_UP_genes_for_avg_heatmap.csv",
  row.names = FALSE
)

avg_mat_early <- make_group_average_expression_matrix(
  seurat_obj = tumor_final,
  features = early_top_genes_avg,
  group_col = "EMT_state_final",
  group_order = c(
    "epithelial",
    "early_hybrid",
    "late_hybrid",
    "full_EMT"
  ),
  assay = "RNA"
)

write.csv(
  avg_mat_early,
  file = "13B_AvgExpressionMatrix_early_hybrid_UP_genes_raw.csv"
)

save_average_expression_heatmap(
  mat = avg_mat_early,
  prefix = "13C_AvgHeatmap_early_hybrid_UP_genes",
  title_text = "Average expression of early_hybrid-up genes",
  width = 7,
  height = 10
)

late_top_genes_avg <- get_top_up_genes_from_deg(
  up_deg_df = res_late_vs_rest$up,
  top_n = top_n_avg_heatmap
)

late_top_genes_avg <- late_top_genes_avg[
  late_top_genes_avg %in% rownames(tumor_final)
]

write.csv(
  data.frame(gene = late_top_genes_avg),
  file = "14A_late_hybrid_TOP_UP_genes_for_avg_heatmap.csv",
  row.names = FALSE
)

avg_mat_late <- make_group_average_expression_matrix(
  seurat_obj = tumor_final,
  features = late_top_genes_avg,
  group_col = "EMT_state_final",
  group_order = c(
    "epithelial",
    "early_hybrid",
    "late_hybrid",
    "full_EMT"
  ),
  assay = "RNA"
)

write.csv(
  avg_mat_late,
  file = "14B_AvgExpressionMatrix_late_hybrid_UP_genes_raw.csv"
)

save_average_expression_heatmap(
  mat = avg_mat_late,
  prefix = "14C_AvgHeatmap_late_hybrid_UP_genes",
  title_text = "Average expression of late_hybrid-up genes",
  width = 7,
  height = 10
)

hybrid_top_genes_avg <- get_top_up_genes_from_deg(
  up_deg_df = res_hybrid_vs_rest$up,
  top_n = top_n_avg_heatmap
)

hybrid_top_genes_avg <- hybrid_top_genes_avg[
  hybrid_top_genes_avg %in% rownames(tumor_final)
]

write.csv(
  data.frame(gene = hybrid_top_genes_avg),
  file = "15A_merged_Hybrid_TOP_UP_genes_for_avg_heatmap.csv",
  row.names = FALSE
)

avg_mat_hybrid <- make_group_average_expression_matrix(
  seurat_obj = tumor_final,
  features = hybrid_top_genes_avg,
  group_col = "EMT_state_hybrid_merged",
  group_order = c(
    "Epithelial",
    "Hybrid",
    "Full_EMT"
  ),
  assay = "RNA"
)

write.csv(
  avg_mat_hybrid,
  file = "15B_AvgExpressionMatrix_merged_Hybrid_UP_genes_raw.csv"
)

save_average_expression_heatmap(
  mat = avg_mat_hybrid,
  prefix = "15C_AvgHeatmap_merged_Hybrid_UP_genes",
  title_text = "Average expression of merged Hybrid-up genes",
  width = 6,
  height = 10
)

saveRDS(
  tumor_final,
  file = "GSE234265_control_2_YFP_tumor_final_after_average_heatmaps.rds"
)

message("All average-expression heatmaps finished successfully.")
