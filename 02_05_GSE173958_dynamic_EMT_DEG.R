# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
base_dir <- "path/to/GSE173958"
analysis_name <- "GSE173958_M1_removeRiboCluster4_pseudoEMT_v3_4_dynamicHeatmap_forcedOrder"
analysis_name_windows <- paste0(analysis_name, "_Windows_step26plus")

checkpoint_dir <- file.path(base_dir, analysis_name, "RDS")

seurat_checkpoint_file <- file.path(
  checkpoint_dir,
  "GSE173958_M1_checkpoint_step25_Seurat_with_pseudoEMT_EMTstate.rds"
)

cds_checkpoint_file <- file.path(
  checkpoint_dir,
  "GSE173958_M1_checkpoint_step25_Monocle3_cds_with_pseudoEMT.rds"
)

required_packages <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "monocle3",
  "SingleCellExperiment",
  "pheatmap"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "以下 R 包未安装，请先安装后再运行：\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(monocle3)
  library(SingleCellExperiment)
  library(pheatmap)
})

options(future.globals.maxSize = 80 * 1024^3)

out_dir <- file.path(base_dir, analysis_name_windows)

subdirs <- c(
  "00_logs",
  "07_pseudoEMT_dynamic_genes",
  "09_DEG_Hybrid_EMT",
  "09_DEG_Hybrid_EMT/Heatmaps",
  "RDS"
)

for (d in subdirs) {
  dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)
}

sink(file.path(out_dir, "00_logs", "run_log_windows_step26plus.txt"), split = TRUE)

cat("============================================================\n")
cat("GSE173958 M1 Windows step26+ resume analysis\n")
cat("Analysis started at:", format(Sys.time()), "\n")
cat("base_dir:", base_dir, "\n")
cat("input checkpoint_dir:", checkpoint_dir, "\n")
cat("out_dir:", out_dir, "\n")
cat("graph_test_cores:", graph_test_cores, "\n")
cat("============================================================\n\n")

save_plot_pdf <- function(plot_obj, filename, width = 8, height = 6) {
  pdf(filename, width = width, height = height)
  print(plot_obj)
  dev.off()
  print(plot_obj)
}

zscore_rows <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  mat_z[mat_z > 2.5] <- 2.5
  mat_z[mat_z < -2.5] <- -2.5
  return(mat_z)
}

make_pseudotime_binned_matrix <- function(expr_mat, pseudotime_vec, genes, n_bins = 100) {
  common_cells <- intersect(colnames(expr_mat), names(pseudotime_vec))
  genes <- intersect(genes, rownames(expr_mat))

  expr_mat <- expr_mat[genes, common_cells, drop = FALSE]
  pt <- pseudotime_vec[common_cells]

  keep <- is.finite(pt)
  expr_mat <- expr_mat[, keep, drop = FALSE]
  pt <- pt[keep]

  if (length(unique(pt)) < 10) {
    stop("pseudoEMT / pseudotime 的唯一值过少，无法进行动态基因分箱热图。")
  }

  q_breaks <- unique(as.numeric(quantile(
    pt,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  )))

  if (length(q_breaks) < 10) {
    bin_id <- cut(pt, breaks = n_bins, include.lowest = TRUE, labels = FALSE)
  } else {
    bin_id <- cut(pt, breaks = q_breaks, include.lowest = TRUE, labels = FALSE)
  }

  bin_id <- as.integer(bin_id)
  names(bin_id) <- names(pt)
  bin_levels <- sort(unique(bin_id[!is.na(bin_id)]))

  avg_list <- lapply(bin_levels, function(b) {
    cells_b <- names(bin_id)[bin_id == b]
    if (length(cells_b) == 1) {
      as.numeric(expr_mat[, cells_b])
    } else {
      Matrix::rowMeans(expr_mat[, cells_b, drop = FALSE])
    }
  })

  avg_mat <- do.call(cbind, avg_list)
  rownames(avg_mat) <- rownames(expr_mat)
  colnames(avg_mat) <- paste0("Bin_", bin_levels)
  return(avg_mat)
}

top_up_down <- function(deg_df, n_each = 50) {
  deg_df2 <- deg_df %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(p_val_adj))

  up <- deg_df2 %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) %>%
    head(n_each)

  down <- deg_df2 %>%
    dplyr::filter(avg_log2FC < 0) %>%
    dplyr::arrange(p_val_adj, avg_log2FC) %>%
    head(n_each)

  unique(c(up$gene, down$gene))
}

safe_doheatmap <- function(seu_obj, features, group_by, title_text, pdf_file, width = 12, height = 18) {
  features <- intersect(features, rownames(seu_obj))
  if (length(features) < 2) {
    warning("热图基因数少于 2，跳过：", title_text)
    return(invisible(NULL))
  }

  p <- DoHeatmap(
    seu_obj,
    features = features,
    group.by = group_by,
    size = 2.5
  ) + ggtitle(title_text)

  save_plot_pdf(p, pdf_file, width = width, height = height)
}

add_direction_and_filter <- function(deg_df, logfc_cutoff = 0.25, padj_cutoff = 0.05) {
  deg_df %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        avg_log2FC > 0 ~ "up_in_ident1",
        avg_log2FC < 0 ~ "down_in_ident1",
        TRUE ~ "no_change"
      ),
      pass_threshold = !is.na(p_val_adj) &
        p_val_adj < padj_cutoff &
        abs(avg_log2FC) >= logfc_cutoff
    ) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
}

run_deg_pair <- function(
    seu_obj,
    meta_col,
    ident1,
    ident2,
    file_prefix,
    out_deg_dir,
    logfc_cutoff = 0.25,
    padj_cutoff = 0.05,
    min_pct = 0.10,
    max_cells_per_ident = 5000
) {
  message("Running DEG: ", file_prefix)

  meta_vec <- seu_obj[[meta_col]][, 1]
  keep_cells <- colnames(seu_obj)[!is.na(meta_vec)]
  sub_obj <- subset(seu_obj, cells = keep_cells)
  sub_obj[[meta_col]] <- droplevels(factor(sub_obj[[meta_col]][, 1], levels = c(ident1, ident2)))

  tab <- table(sub_obj[[meta_col]][, 1])
  write.csv(
    as.data.frame(tab),
    file.path(out_deg_dir, paste0(file_prefix, "_cell_counts.csv")),
    row.names = FALSE
  )

  if (!all(c(ident1, ident2) %in% names(tab)) || any(tab[c(ident1, ident2)] < 3)) {
    warning("该比较中某组细胞数 < 3，跳过：", file_prefix)
    return(NULL)
  }

  DefaultAssay(sub_obj) <- "RNA"
  Idents(sub_obj) <- meta_col

  deg <- FindMarkers(
    object = sub_obj,
    assay = "RNA",
    ident.1 = ident1,
    ident.2 = ident2,
    logfc.threshold = 0,
    min.pct = min_pct,
    test.use = "wilcox",
    slot = "data",
    densify = FALSE,
    max.cells.per.ident = max_cells_per_ident,
    random.seed = 1234
  )

  deg_all <- add_direction_and_filter(deg, logfc_cutoff, padj_cutoff)

  write.csv(
    deg_all,
    file.path(out_deg_dir, paste0(file_prefix, "_ALL_genes_with_direction.csv")),
    row.names = FALSE
  )

  deg_sig <- deg_all %>% dplyr::filter(pass_threshold)
  write.csv(
    deg_sig,
    file.path(
      out_deg_dir,
      paste0(file_prefix, "_SIG_absLog2FC", logfc_cutoff, "_padj", padj_cutoff, ".csv")
    ),
    row.names = FALSE
  )

  heat_genes <- top_up_down(deg, n_each = 50)
  heat_genes <- intersect(heat_genes, rownames(sub_obj))

  if (length(heat_genes) >= 2) {
    sub_obj <- ScaleData(
      sub_obj,
      assay = "RNA",
      features = heat_genes,
      verbose = FALSE
    )

    safe_doheatmap(
      seu_obj = sub_obj,
      features = heat_genes,
      group_by = meta_col,
      title_text = paste0(file_prefix, ": top 50 up + top 50 down"),
      pdf_file = file.path(out_deg_dir, "Heatmaps", paste0(file_prefix, "_Heatmap_top50UpDown.pdf")),
      width = 12,
      height = 18
    )
  }

  return(list(deg_all = deg_all, deg_sig = deg_sig))
}

if (!file.exists(seurat_checkpoint_file)) {
  stop("找不到 Seurat checkpoint 文件：\n", seurat_checkpoint_file)
}

if (!file.exists(cds_checkpoint_file)) {
  stop("找不到 Monocle3 cds checkpoint 文件：\n", cds_checkpoint_file)
}

m1_seu_clean <- readRDS(seurat_checkpoint_file)
cds_clean <- readRDS(cds_checkpoint_file)

cat("成功读取 checkpoint。\n")
cat("Seurat cells:", ncol(m1_seu_clean), "\n")
cat("Seurat genes:", nrow(m1_seu_clean), "\n")
cat("cds cells:", ncol(cds_clean), "\n")
cat("cds genes:", nrow(cds_clean), "\n\n")

required_meta <- c("pseudoEMT", "EMT_state_initial")
missing_meta <- setdiff(required_meta, colnames(m1_seu_clean@meta.data))
if (length(missing_meta) > 0) {
  stop("Seurat checkpoint 缺少必要 metadata：", paste(missing_meta, collapse = ", "))
}

m1_seu_clean$EMT_state_initial <- factor(
  as.character(m1_seu_clean$EMT_state_initial),
  levels = c("Epithelial", "Early_hybrid", "Late_hybrid", "Mesenchymal")
)

pt_clean <- m1_seu_clean$pseudoEMT
names(pt_clean) <- colnames(m1_seu_clean)
pt_clean <- pt_clean[is.finite(pt_clean)]

cat("EMT_state_initial 细胞数：\n")
print(table(m1_seu_clean$EMT_state_initial, useNA = "ifany"))
cat("\n")

graph_res_clean <- graph_test(
  cds_clean,
  neighbor_graph = "principal_graph",
  cores = graph_test_cores
)

graph_res_clean <- graph_res_clean %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(q_value, dplyr::desc(morans_I))

write.csv(
  graph_res_clean,
  file.path(out_dir, "07_pseudoEMT_dynamic_genes", "PostFilter_graph_test_all_genes_principal_graph.csv"),
  row.names = FALSE
)

dynamic_res_clean <- graph_res_clean %>%
  dplyr::filter(q_value < 0.05, morans_I > 0.1) %>%
  dplyr::arrange(q_value, dplyr::desc(morans_I))

if (nrow(dynamic_res_clean) == 0) {
  stop("没有筛到 q_value < 0.05 且 Moran's I > 0.1 的动态基因，请查看 graph_test 结果。")
}

top_n_dynamic_clean <- min(n_top_dynamic, nrow(dynamic_res_clean))
top_dynamic_clean <- dynamic_res_clean %>% head(top_n_dynamic_clean)
top_dynamic_genes_clean <- top_dynamic_clean$gene

write.csv(
  top_dynamic_clean,
  file.path(
    out_dir,
    "07_pseudoEMT_dynamic_genes",
    paste0("PostFilter_pseudoEMT_dynamic_genes_top", top_n_dynamic_clean, ".csv")
  ),
  row.names = FALSE
)

cat("pseudoEMT 动态基因 q<0.05 & Moran's I>0.1 数量：", nrow(dynamic_res_clean), "\n")
cat("用于动态热图的 top genes 数量：", top_n_dynamic_clean, "\n\n")

norm_counts_cds_clean <- normalized_counts(cds_clean)

common_dynamic_genes_clean <- intersect(
  top_dynamic_genes_clean,
  rownames(norm_counts_cds_clean)
)

if (length(common_dynamic_genes_clean) < 50) {
  stop("可用于热图的动态基因少于 50 个，请检查 normalized_counts(cds_clean)。")
}

bin_mat_clean <- make_pseudotime_binned_matrix(
  expr_mat = norm_counts_cds_clean,
  pseudotime_vec = pt_clean,
  genes = common_dynamic_genes_clean,
  n_bins = n_pseudotime_bins
)

bin_mat_z_clean <- zscore_rows(bin_mat_clean)

hc_rows_clean <- hclust(dist(bin_mat_z_clean), method = "ward.D2")
gene_module_raw_clean <- cutree(hc_rows_clean, k = 6)

gene_module_df_clean <- data.frame(
  gene = names(gene_module_raw_clean),
  raw_module = as.integer(gene_module_raw_clean),
  stringsAsFactors = FALSE
)

module_peak_clean <- lapply(sort(unique(gene_module_raw_clean)), function(mod) {
  genes_mod <- names(gene_module_raw_clean)[gene_module_raw_clean == mod]
  avg_curve <- colMeans(bin_mat_z_clean[genes_mod, , drop = FALSE])
  peak_bin <- which.max(avg_curve)
  data.frame(raw_module = mod, peak_bin = peak_bin)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(peak_bin)

module_order_clean <- c("E", "H1", "H2", "H3", "H4", "M")

module_name_map_clean <- data.frame(
  raw_module = module_peak_clean$raw_module,
  pseudoEMT_gene_module = module_order_clean,
  peak_bin = module_peak_clean$peak_bin,
  stringsAsFactors = FALSE
)

gene_module_df_clean <- gene_module_df_clean %>%
  dplyr::left_join(module_name_map_clean, by = "raw_module") %>%
  dplyr::mutate(
    pseudoEMT_gene_module = factor(pseudoEMT_gene_module, levels = module_order_clean)
  )

write.csv(
  gene_module_df_clean %>%
    dplyr::arrange(pseudoEMT_gene_module, gene),
  file.path(
    out_dir,
    "07_pseudoEMT_dynamic_genes",
    "PostFilter_pseudoEMT_top_dynamic_genes_six_modules_E_H1_H2_H3_H4_M.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_name_map_clean,
  file.path(
    out_dir,
    "07_pseudoEMT_dynamic_genes",
    "PostFilter_pseudoEMT_six_module_peak_order_mapping_E_H1_H2_H3_H4_M.csv"
  ),
  row.names = FALSE
)

hc_gene_order_clean <- rownames(bin_mat_z_clean)[hc_rows_clean$order]

gene_order_for_heatmap_clean <- gene_module_df_clean %>%
  dplyr::mutate(
    module_order = as.integer(pseudoEMT_gene_module),
    hc_order = match(gene, hc_gene_order_clean)
  ) %>%
  dplyr::arrange(module_order, hc_order) %>%
  dplyr::pull(gene)

gene_order_for_heatmap_clean <- intersect(gene_order_for_heatmap_clean, rownames(bin_mat_z_clean))
bin_mat_z_ordered_clean <- bin_mat_z_clean[gene_order_for_heatmap_clean, , drop = FALSE]

anno_row_clean <- gene_module_df_clean %>%
  dplyr::select(gene, Module = pseudoEMT_gene_module) %>%
  dplyr::filter(gene %in% rownames(bin_mat_z_ordered_clean)) %>%
  dplyr::arrange(match(gene, rownames(bin_mat_z_ordered_clean))) %>%
  as.data.frame()

rownames(anno_row_clean) <- anno_row_clean$gene
anno_row_clean$gene <- NULL
anno_row_clean$Module <- factor(anno_row_clean$Module, levels = module_order_clean)

module_sizes_clean <- as.numeric(table(anno_row_clean$Module)[module_order_clean])
module_sizes_clean[is.na(module_sizes_clean)] <- 0
gaps_row_clean <- cumsum(module_sizes_clean)
gaps_row_clean <- gaps_row_clean[-length(gaps_row_clean)]
gaps_row_clean <- gaps_row_clean[gaps_row_clean > 0 & gaps_row_clean < nrow(bin_mat_z_ordered_clean)]

pdf(
  file.path(
    out_dir,
    "07_pseudoEMT_dynamic_genes",
    paste0(
      "PostFilter_Heatmap_top",
      length(common_dynamic_genes_clean),
      "_dynamic_genes_",
      n_pseudotime_bins,
      "bins_E_top_H1_H2_H3_H4_M_bottom.pdf"
    )
  ),
  width = 11,
  height = 16
)

pheatmap::pheatmap(
  bin_mat_z_ordered_clean,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = anno_row_clean,
  gaps_row = gaps_row_clean,
  border_color = NA,
  main = paste0(
    "pseudoEMT dynamic genes top ",
    length(common_dynamic_genes_clean),
    ": E-H1-H2-H3-H4-M"
  )
)

dev.off()

rm(norm_counts_cds_clean, bin_mat_clean, bin_mat_z_clean, bin_mat_z_ordered_clean)
gc()

DefaultAssay(m1_seu_clean) <- "RNA"

m1_seu_clean[["RNA"]] <- JoinLayers(m1_seu_clean[["RNA"]])

m1_seu_clean <- NormalizeData(
  m1_seu_clean,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

m1_seu_clean$Hybrid_all_DEG_group <- dplyr::case_when(
  m1_seu_clean$EMT_state_initial %in% c("Early_hybrid", "Late_hybrid") ~ "Hybrid_all",
  m1_seu_clean$EMT_state_initial == "Epithelial" ~ "Epithelial",
  m1_seu_clean$EMT_state_initial == "Mesenchymal" ~ "Mesenchymal",
  TRUE ~ NA_character_
)

m1_seu_clean$Hybrid_all_vs_Other_group <- dplyr::case_when(
  m1_seu_clean$EMT_state_initial %in% c("Early_hybrid", "Late_hybrid") ~ "Hybrid_all",
  m1_seu_clean$EMT_state_initial %in% c("Epithelial", "Mesenchymal") ~ "Other",
  TRUE ~ NA_character_
)

m1_seu_clean$Early_hybrid_DEG_group <- dplyr::case_when(
  m1_seu_clean$EMT_state_initial == "Early_hybrid" ~ "Early_hybrid",
  m1_seu_clean$EMT_state_initial == "Epithelial" ~ "Epithelial",
  m1_seu_clean$EMT_state_initial == "Mesenchymal" ~ "Mesenchymal",
  TRUE ~ NA_character_
)

m1_seu_clean$Early_hybrid_vs_Other_group <- dplyr::case_when(
  m1_seu_clean$EMT_state_initial == "Early_hybrid" ~ "Early_hybrid",
  m1_seu_clean$EMT_state_initial %in% c("Epithelial", "Late_hybrid", "Mesenchymal") ~ "Other",
  TRUE ~ NA_character_
)

m1_seu_clean$Late_hybrid_DEG_group <- dplyr::case_when(
  m1_seu_clean$EMT_state_initial == "Late_hybrid" ~ "Late_hybrid",
  m1_seu_clean$EMT_state_initial == "Epithelial" ~ "Epithelial",
  m1_seu_clean$EMT_state_initial == "Mesenchymal" ~ "Mesenchymal",
  TRUE ~ NA_character_
)

m1_seu_clean$Late_hybrid_vs_Other_group <- dplyr::case_when(
  m1_seu_clean$EMT_state_initial == "Late_hybrid" ~ "Late_hybrid",
  m1_seu_clean$EMT_state_initial %in% c("Epithelial", "Early_hybrid", "Mesenchymal") ~ "Other",
  TRUE ~ NA_character_
)

deg_plan <- tibble::tribble(
  ~meta_col, ~ident1, ~ident2, ~file_prefix,
  "Hybrid_all_DEG_group", "Hybrid_all", "Epithelial", "PostFilter_RNA_DEG_Hybrid_all_vs_Epithelial",
  "Hybrid_all_DEG_group", "Hybrid_all", "Mesenchymal", "PostFilter_RNA_DEG_Hybrid_all_vs_Mesenchymal",
  "Hybrid_all_vs_Other_group", "Hybrid_all", "Other", "PostFilter_RNA_DEG_Hybrid_all_vs_Other",

  "Early_hybrid_DEG_group", "Early_hybrid", "Epithelial", "PostFilter_RNA_DEG_Early_hybrid_vs_Epithelial",
  "Early_hybrid_DEG_group", "Early_hybrid", "Mesenchymal", "PostFilter_RNA_DEG_Early_hybrid_vs_Mesenchymal",
  "Early_hybrid_vs_Other_group", "Early_hybrid", "Other", "PostFilter_RNA_DEG_Early_hybrid_vs_Other",

  "Late_hybrid_DEG_group", "Late_hybrid", "Epithelial", "PostFilter_RNA_DEG_Late_hybrid_vs_Epithelial",
  "Late_hybrid_DEG_group", "Late_hybrid", "Mesenchymal", "PostFilter_RNA_DEG_Late_hybrid_vs_Mesenchymal",
  "Late_hybrid_vs_Other_group", "Late_hybrid", "Other", "PostFilter_RNA_DEG_Late_hybrid_vs_Other"
)

deg_summary <- list()

for (i in seq_len(nrow(deg_plan))) {
  res_i <- run_deg_pair(
    seu_obj = m1_seu_clean,
    meta_col = deg_plan$meta_col[i],
    ident1 = deg_plan$ident1[i],
    ident2 = deg_plan$ident2[i],
    file_prefix = deg_plan$file_prefix[i],
    out_deg_dir = file.path(out_dir, "09_DEG_Hybrid_EMT"),
    logfc_cutoff = deg_logfc_cutoff,
    padj_cutoff = deg_padj_cutoff,
    min_pct = deg_min_pct,
    max_cells_per_ident = deg_max_cells_per_ident
  )

  if (!is.null(res_i)) {
    deg_summary[[deg_plan$file_prefix[i]]] <- data.frame(
      comparison = deg_plan$file_prefix[i],
      n_all_tested = nrow(res_i$deg_all),
      n_sig = nrow(res_i$deg_sig),
      n_up_in_ident1 = sum(res_i$deg_sig$direction == "up_in_ident1"),
      n_down_in_ident1 = sum(res_i$deg_sig$direction == "down_in_ident1")
    )
  }

  gc()
}

deg_summary_df <- dplyr::bind_rows(deg_summary)

write.csv(
  deg_summary_df,
  file.path(
    out_dir,
    "09_DEG_Hybrid_EMT",
    paste0("PostFilter_RNA_DEG_summary_absLog2FC", deg_logfc_cutoff, "_padj", deg_padj_cutoff, ".csv")
  ),
  row.names = FALSE
)

cat("DEG summary:\n")
print(deg_summary_df)

saveRDS(
  m1_seu_clean,
  file.path(
    out_dir,
    "RDS",
    "GSE173958_M1_Windows_step26plus_Seurat_with_dynamicEMT_DEG.rds"
  )
)

saveRDS(
  cds_clean,
  file.path(
    out_dir,
    "RDS",
    "GSE173958_M1_Windows_step26plus_Monocle3_cds.rds"
  )
)

cat("\n============================================================\n")
cat("Windows step26+ 分析完成。\n")
cat("Analysis finished at:", format(Sys.time()), "\n")
cat("结果输出文件夹：", out_dir, "\n")
cat("============================================================\n")

sink()
