# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
base_dir <- "path/to/GSE173958"

classification_file <- file.path(
  base_dir,
  "pdac_mouse1_final_classifications.txt"
)

analysis_name <- "new_GSE173958_M1_removeRiboCluster4_pseudoEMT_v3_4_dynamicHeatmap_forcedOrder"

initial_monocle_clusters_to_remove <- c("4")
run_dynamic_gene_analysis <- FALSE
set.seed(1234)

required_packages <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "monocle3",
  "SingleCellExperiment",
  "pheatmap",
  "ggrepel"
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
  library(ggrepel)
})

options(future.globals.maxSize = 40 * 1024^3)

out_dir <- file.path(base_dir, analysis_name)

subdirs <- c(
  "00_logs",
  "01_final_singlet_matching",
  "02_QC",
  "03_PreFilter_Seurat_final_singlet",
  "04_FinalTranscriptomeFilter_Monocle_initial",
  "05_PostFilter_Seurat_clean",
  "06_PostFilter_Monocle3_pseudoEMT",
  "07_pseudoEMT_dynamic_genes",
  "08_EMT_state_initial_annotation",
  "09_DEG_Hybrid_EMT",
  "09_DEG_Hybrid_EMT/Heatmaps",
  "RDS"
)

for (d in subdirs) {
  dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)
}

sink(file.path(out_dir, "00_logs", "run_log.txt"), split = TRUE)

cat("============================================================\n")
cat("GSE173958 M1 final_singlet + remove high-ribo cluster 4 pseudoEMT v3.3 RNA-DEG memory-safe\n")
cat("Analysis started at:", format(Sys.time()), "\n")
cat("base_dir:", base_dir, "\n")
cat("classification_file:", classification_file, "\n")
cat("out_dir:", out_dir, "\n")
cat("initial_monocle_clusters_to_remove:", paste(initial_monocle_clusters_to_remove, collapse = ", "), "\n")
cat("run_dynamic_gene_analysis:", run_dynamic_gene_analysis, "\n")
cat("============================================================\n\n")

save_plot_pdf <- function(plot_obj, filename, width = 8, height = 6) {
  pdf(filename, width = width, height = height)
  print(plot_obj)
  dev.off()
  print(plot_obj)
}

find_one_file <- function(folder, pattern) {
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(files) != 1) {
    stop(
      "在文件夹中没有找到唯一匹配文件：\n",
      folder, "\npattern = ", pattern, "\n",
      "当前匹配数量 = ", length(files), "\n",
      "匹配结果：\n",
      paste(files, collapse = "\n")
    )
  }
  return(files)
}

read_one_10x_prefixed <- function(sample_dir) {
  matrix_file   <- find_one_file(sample_dir, "matrix\\.mtx\\.gz$")
  features_file <- find_one_file(sample_dir, "features\\.tsv\\.gz$")
  barcodes_file <- find_one_file(sample_dir, "barcodes\\.tsv\\.gz$")

  mat <- ReadMtx(
    mtx = matrix_file,
    features = features_file,
    cells = barcodes_file,
    feature.column = 2
  )

  return(mat)
}

remove_10x_suffix <- function(x) {
  gsub("-1$", "", x)
}

match_gene_case <- function(target_genes, feature_genes) {
  direct_hits <- intersect(target_genes, feature_genes)
  if (length(direct_hits) >= max(5, 0.3 * length(target_genes))) {
    return(direct_hits)
  }

  to_mouse_case <- function(x) {
    ifelse(
      nchar(x) <= 1,
      toupper(x),
      paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
    )
  }

  converted <- to_mouse_case(target_genes)
  converted_hits <- intersect(converted, feature_genes)

  unique(c(direct_hits, converted_hits))
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

if (!dir.exists(base_dir)) {
  stop("base_dir 不存在：", base_dir)
}

if (!file.exists(classification_file)) {
  stop(
    "未找到 classification_file：\n",
    classification_file, "\n",
    "请把 pdac_mouse1_final_classifications(2).txt 放入 GSE173958 文件夹，",
    "或修改 classification_file 为真实路径。"
  )
}

classification_df <- read.delim(
  classification_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("cellID", "type", "cloneID")
if (!all(required_cols %in% colnames(classification_df))) {
  stop(
    "classification 文件缺少必要列：",
    paste(setdiff(required_cols, colnames(classification_df)), collapse = ", ")
  )
}

classification_df$type <- as.character(classification_df$type)
classification_df$cellID <- as.character(classification_df$cellID)
classification_df$cloneID <- as.character(classification_df$cloneID)

type_summary <- classification_df %>%
  dplyr::count(type, name = "n_cells") %>%
  dplyr::arrange(dplyr::desc(n_cells))

write.csv(
  type_summary,
  file.path(out_dir, "01_final_singlet_matching", "classification_type_summary.csv"),
  row.names = FALSE
)

cat("classification 文件中各 type 数量：\n")
print(type_summary)
cat("\n")

final_singlet_df <- classification_df %>%
  dplyr::filter(type == "final_singlet") %>%
  dplyr::distinct(cellID, .keep_all = TRUE)

cat("final_singlet 唯一 cellID 数量：", nrow(final_singlet_df), "\n\n")

final_singlet_prefix_summary <- final_singlet_df %>%
  dplyr::mutate(
    classification_prefix = sub("_.*$", "", cellID)
  ) %>%
  dplyr::count(classification_prefix, name = "n_final_singlet") %>%
  dplyr::arrange(dplyr::desc(n_final_singlet))

write.csv(
  final_singlet_prefix_summary,
  file.path(out_dir, "01_final_singlet_matching", "final_singlet_prefix_summary.csv"),
  row.names = FALSE
)

cat("final_singlet 按 cellID 前缀统计：\n")
print(final_singlet_prefix_summary)
cat("\n")

all_subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

m1_dirs <- all_subdirs[
  grepl("GSM528348[2-7]_M1-", basename(all_subdirs))
]

if (length(m1_dirs) != 6) {
  cat("当前识别到的 M1 文件夹：\n")
  print(basename(m1_dirs))
  stop("未识别到 6 个 M1 子文件夹，请检查 base_dir 和文件夹命名。")
}

m1_dirs <- sort(m1_dirs)

cat("识别到以下 M1 子文件夹：\n")
print(basename(m1_dirs))
cat("\n")

sample_mapping <- data.frame(
  folder_key = c("M1-PT", "M1-SS", "M1-Met", "M1-Liver", "M1-Lung", "M1-CTC"),
  site_name = c("PT", "SS", "Met", "Liver", "Lung", "CTC"),
  classification_prefix = c("PT", "PTab", "Met", "Liver", "Lung", "Blood"),
  location_label = c("Primary tumor", "Surgical-site metastasis", "Peritoneal metastasis", "Liver", "Lung", "CTCs"),
  stringsAsFactors = FALSE
)

write.csv(
  sample_mapping,
  file.path(out_dir, "01_final_singlet_matching", "sample_prefix_mapping_used.csv"),
  row.names = FALSE
)

cat("本次用于匹配 cellID 的样本前缀映射：\n")
print(sample_mapping)
cat("\n")

seurat_list <- list()
matching_summary_list <- list()
matched_cellids_all <- character(0)

for (sample_dir in m1_dirs) {

  sample_name <- basename(sample_dir)

  map_row <- sample_mapping[
    vapply(
      sample_mapping$folder_key,
      function(k) grepl(k, sample_name, fixed = TRUE),
      logical(1)
    ),
    ,
    drop = FALSE
  ]

  if (nrow(map_row) != 1) {
    stop("无法为文件夹确定唯一前缀映射：", sample_name)
  }

  site_name <- map_row$site_name[1]
  class_prefix <- map_row$classification_prefix[1]
  location_label <- map_row$location_label[1]

  cat("------------------------------------------------------------\n")
  cat("Reading:", sample_name, "\n")
  cat("site_name:", site_name, "\n")
  cat("classification_prefix:", class_prefix, "\n")

  counts <- read_one_10x_prefixed(sample_dir)

  raw_barcode <- colnames(counts)
  barcode_no_suffix <- remove_10x_suffix(raw_barcode)

  unified_cellID <- paste0(class_prefix, "_", barcode_no_suffix)
  colnames(counts) <- unified_cellID

  target_cellids <- final_singlet_df$cellID[
    startsWith(final_singlet_df$cellID, paste0(class_prefix, "_"))
  ]

  matched_cellids <- intersect(colnames(counts), target_cellids)
  unmatched_final_singlet_in_this_sample <- setdiff(target_cellids, colnames(counts))

  cat("该样本 10X 矩阵原始 cell barcode 数：", ncol(counts), "\n")
  cat("该样本 final_singlet 目标 cellID 数：", length(target_cellids), "\n")
  cat("成功匹配到表达矩阵的 final_singlet 数：", length(matched_cellids), "\n")
  cat("未匹配到表达矩阵的 final_singlet 数：", length(unmatched_final_singlet_in_this_sample), "\n")

  if (length(matched_cellids) == 0) {
    stop(
      "样本 ", sample_name, " 未匹配到任何 final_singlet。\n",
      "请检查 barcodes.tsv.gz 是否带 -1 后缀，或 classifications cellID 前缀映射。"
    )
  }

  counts_sub <- counts[, matched_cellids, drop = FALSE]

  seu_i <- CreateSeuratObject(
    counts = counts_sub,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )

  seu_i$sample_id <- sample_name
  seu_i$mouse <- "M1"
  seu_i$site <- site_name
  seu_i$classification_prefix <- class_prefix
  seu_i$location_label <- location_label

  clone_map <- final_singlet_df$cloneID
  names(clone_map) <- final_singlet_df$cellID
  seu_i$cloneID <- unname(clone_map[colnames(seu_i)])

  if (any(is.na(seu_i$cloneID))) {
    warning("样本 ", sample_name, " 中存在 cloneID 未成功映射的细胞。")
  }

  seurat_list[[sample_name]] <- seu_i

  matching_summary_list[[sample_name]] <- data.frame(
    sample_name = sample_name,
    site_name = site_name,
    classification_prefix = class_prefix,
    matrix_cells_before_subset = ncol(counts),
    target_final_singlet_from_classification = length(target_cellids),
    matched_final_singlet_in_matrix = length(matched_cellids),
    unmatched_final_singlet_not_in_matrix = length(unmatched_final_singlet_in_this_sample),
    stringsAsFactors = FALSE
  )

  if (length(unmatched_final_singlet_in_this_sample) > 0) {
    write.csv(
      data.frame(cellID = unmatched_final_singlet_in_this_sample),
      file.path(
        out_dir,
        "01_final_singlet_matching",
        paste0("unmatched_final_singlet_", site_name, ".csv")
      ),
      row.names = FALSE
    )
  }

  matched_cellids_all <- c(matched_cellids_all, matched_cellids)

  rm(counts, counts_sub, seu_i)
  gc()
}

matching_summary <- bind_rows(matching_summary_list)

write.csv(
  matching_summary,
  file.path(out_dir, "01_final_singlet_matching", "final_singlet_matching_summary_by_sample.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("final_singlet 与表达矩阵匹配汇总：\n")
print(matching_summary)
cat("所有样本成功匹配到的 final_singlet 总数：", sum(matching_summary$matched_final_singlet_in_matrix), "\n")
cat("final_classifications 中 final_singlet 总数：", nrow(final_singlet_df), "\n")
cat("============================================================\n\n")

m1_seu <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  project = "GSE173958_M1_final_singlet"
)

rm(seurat_list)
gc()

cat("合并 final_singlet 后细胞数：", ncol(m1_seu), "\n")
cat("合并 final_singlet 后基因数：", nrow(m1_seu), "\n\n")

write.csv(
  data.frame(
    stage = c(
      "classification_final_singlet",
      "matched_to_expression_matrices",
      "merged_final_singlet_seurat_object"
    ),
    n_cells = c(
      nrow(final_singlet_df),
      length(unique(matched_cellids_all)),
      ncol(m1_seu)
    )
  ),
  file.path(out_dir, "01_final_singlet_matching", "cell_number_key_stages_before_QC.csv"),
  row.names = FALSE
)

saveRDS(
  m1_seu,
  file.path(out_dir, "RDS", "GSE173958_M1_final_singlet_merged_before_QC.rds")
)

m1_seu[["percent.mt"]] <- PercentageFeatureSet(m1_seu, pattern = "^mt-")
m1_seu[["percent.ribo"]] <- PercentageFeatureSet(m1_seu, pattern = "^Rpl|^Rps")

qc_summary_before <- m1_seu@meta.data %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    median_percent_ribo = median(percent.ribo),
    .groups = "drop"
  )

write.csv(
  qc_summary_before,
  file.path(out_dir, "02_QC", "QC_summary_final_singlet_before_article_QC.csv"),
  row.names = FALSE
)

cat("final_singlet 合并后、文章 QC 前统计：\n")
print(qc_summary_before)
cat("\n")

p_qc_before <- VlnPlot(
  m1_seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  group.by = "site",
  pt.size = 0.05,
  ncol = 2
)

save_plot_pdf(
  p_qc_before,
  file.path(out_dir, "02_QC", "QC_violin_final_singlet_before_article_QC.pdf"),
  width = 14,
  height = 9
)

cells_before_qc <- ncol(m1_seu)

m1_seu <- subset(
  m1_seu,
  subset = nFeature_RNA > 1000 & percent.mt < 20
)

cells_after_qc <- ncol(m1_seu)

cat("文章参考 QC 前细胞数：", cells_before_qc, "\n")
cat("文章参考 QC 后细胞数：", cells_after_qc, "\n")
cat("QC 剔除细胞数：", cells_before_qc - cells_after_qc, "\n\n")

write.csv(
  data.frame(
    stage = c(
      "final_singlet_matched_before_article_QC",
      "after_article_QC_nFeature_gt1000_percentMT_lt20"
    ),
    n_cells = c(cells_before_qc, cells_after_qc)
  ),
  file.path(out_dir, "02_QC", "cell_number_before_after_article_QC.csv"),
  row.names = FALSE
)

qc_summary_after <- m1_seu@meta.data %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    median_percent_ribo = median(percent.ribo),
    .groups = "drop"
  )

write.csv(
  qc_summary_after,
  file.path(out_dir, "02_QC", "QC_summary_final_singlet_after_article_QC.csv"),
  row.names = FALSE
)

cat("文章参考 QC 后统计：\n")
print(qc_summary_after)
cat("\n")

p_qc_after <- VlnPlot(
  m1_seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  group.by = "site",
  pt.size = 0.05,
  ncol = 2
)

save_plot_pdf(
  p_qc_after,
  file.path(out_dir, "02_QC", "QC_violin_final_singlet_after_article_QC.pdf"),
  width = 14,
  height = 9
)

saveRDS(
  m1_seu,
  file.path(out_dir, "RDS", "GSE173958_M1_final_singlet_after_article_QC.rds")
)

DefaultAssay(m1_seu) <- "RNA"

m1_seu <- NormalizeData(
  m1_seu,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

s_genes <- match_gene_case(cc.genes.updated.2019$s.genes, rownames(m1_seu))
g2m_genes <- match_gene_case(cc.genes.updated.2019$g2m.genes, rownames(m1_seu))

cat("匹配到 S phase genes:", length(s_genes), "\n")
cat("匹配到 G2M phase genes:", length(g2m_genes), "\n\n")

if (length(s_genes) < 10 || length(g2m_genes) < 10) {
  stop("细胞周期基因匹配过少，请检查基因名格式。")
}

m1_seu <- CellCycleScoring(
  m1_seu,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)

m1_seu <- SCTransform(
  m1_seu,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
  vst.flavor = "v2",
  conserve.memory = TRUE,
  ncells = 3000,
  verbose = TRUE
)

m1_seu <- RunPCA(m1_seu, npcs = 50, verbose = FALSE)

p_elbow_prefilter <- ElbowPlot(m1_seu, ndims = 50)
save_plot_pdf(
  p_elbow_prefilter,
  file.path(out_dir, "03_PreFilter_Seurat_final_singlet", "PCA_elbowplot_prefilter.pdf"),
  width = 7,
  height = 5
)

dims_use <- 1:20

m1_seu <- RunUMAP(
  m1_seu,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)

m1_seu <- FindNeighbors(
  m1_seu,
  reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)

m1_seu <- FindClusters(
  m1_seu,
  resolution = 0.5,
  verbose = FALSE
)

p_prefilter_cluster <- DimPlot(
  m1_seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("M1 final_singlet before final transcriptome filters: Seurat clusters")

save_plot_pdf(
  p_prefilter_cluster,
  file.path(out_dir, "03_PreFilter_Seurat_final_singlet", "UMAP_Seurat_clusters_prefilter.pdf"),
  width = 8,
  height = 6.5
)

DefaultAssay(m1_seu) <- "SCT"

epi_genes_raw <- c("Epcam", "Muc1", "Cdh1", "Krt8", "Krt18", "Krt19")
mes_genes_raw <- c("Sparc", "Zeb2", "Col3a1", "Vim", "Fn1", "Col1a1")
emt_tf_genes_raw <- c("Zeb1", "Zeb2", "Snai1", "Snai2", "Prrx1")
hepatic_genes_raw <- c("Alb", "Apoa1", "Apoa2", "Ttr", "Fabp1", "Apob", "Cyp2e1", "Cyp3a11")
marker_genes_raw <- unique(c(epi_genes_raw, mes_genes_raw, emt_tf_genes_raw, hepatic_genes_raw))

epi_genes <- intersect(epi_genes_raw, rownames(m1_seu))
mes_genes <- intersect(mes_genes_raw, rownames(m1_seu))
emt_tf_genes <- intersect(emt_tf_genes_raw, rownames(m1_seu))
hepatic_genes <- intersect(hepatic_genes_raw, rownames(m1_seu))
marker_genes <- intersect(marker_genes_raw, rownames(m1_seu))

cat("Epithelial genes found:", paste(epi_genes, collapse = ", "), "\n")
cat("Mesenchymal genes found:", paste(mes_genes, collapse = ", "), "\n")
cat("EMT TF genes found:", paste(emt_tf_genes, collapse = ", "), "\n")
cat("Hepatic genes found:", paste(hepatic_genes, collapse = ", "), "\n\n")

if (length(epi_genes) >= 2) {
  m1_seu <- AddModuleScore(m1_seu, features = list(epi_genes), name = "EpithelialScore")
} else {
  stop("Epithelial marker 找到少于 2 个，无法继续。")
}

if (length(mes_genes) >= 2) {
  m1_seu <- AddModuleScore(m1_seu, features = list(mes_genes), name = "MesenchymalScore")
} else {
  stop("Mesenchymal marker 找到少于 2 个，无法继续。")
}

if (length(hepatic_genes) >= 2) {
  m1_seu <- AddModuleScore(m1_seu, features = list(hepatic_genes), name = "HepaticScore")
} else {
  m1_seu$HepaticScore1 <- 0
  warning("Hepatic genes 匹配少于 2 个，HepaticScore1 设为 0。")
}

p_prefilter_scores <- FeaturePlot(
  m1_seu,
  features = c("EpithelialScore1", "MesenchymalScore1", "HepaticScore1", "percent.ribo"),
  reduction = "umap",
  ncol = 2,
  order = TRUE
)

save_plot_pdf(
  p_prefilter_scores,
  file.path(out_dir, "03_PreFilter_Seurat_final_singlet", "UMAP_scores_Epi_Mes_Hepatic_Ribo_prefilter.pdf"),
  width = 12,
  height = 10
)

DefaultAssay(m1_seu) <- "RNA"

m1_seu[["RNA"]] <- JoinLayers(m1_seu[["RNA"]])

counts_mat_prefilter <- GetAssayData(
  m1_seu,
  assay = "RNA",
  layer = "counts"
)

cell_metadata_prefilter <- m1_seu@meta.data
gene_metadata_prefilter <- data.frame(
  gene_short_name = rownames(counts_mat_prefilter),
  row.names = rownames(counts_mat_prefilter)
)

cds_initial <- new_cell_data_set(
  expression_data = counts_mat_prefilter,
  cell_metadata = cell_metadata_prefilter,
  gene_metadata = gene_metadata_prefilter
)

rm(counts_mat_prefilter)
gc()

cds_initial <- preprocess_cds(
  cds_initial,
  num_dim = 20
)

cds_initial <- align_cds(
  cds_initial,
  alignment_group = "site",
  residual_model_formula_str = "~ S.Score + G2M.Score + percent.mt"
)

cds_initial <- reduce_dimension(
  cds_initial,
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

cds_initial <- cluster_cells(
  cds_initial,
  reduction_method = "UMAP"
)

initial_monocle_cluster <- as.character(
  clusters(cds_initial, reduction_method = "UMAP")
)

colData(cds_initial)$monocle_initial_cluster <- initial_monocle_cluster
m1_seu$monocle_initial_cluster <- initial_monocle_cluster[colnames(m1_seu)]

p_initial_monocle_cluster <- plot_cells(
  cds_initial,
  color_cells_by = "monocle_initial_cluster",
  label_cell_groups = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = FALSE
) + ggtitle("Initial Monocle3 clusters for final transcriptome filters")

save_plot_pdf(
  p_initial_monocle_cluster,
  file.path(out_dir, "04_FinalTranscriptomeFilter_Monocle_initial", "Initial_Monocle3_clusters_for_filtering.pdf"),
  width = 8,
  height = 6.5
)

p_initial_ribo <- plot_cells(
  cds_initial,
  color_cells_by = "percent.ribo",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = FALSE
) + ggtitle("Initial Monocle3 UMAP colored by percent.ribo")

save_plot_pdf(
  p_initial_ribo,
  file.path(out_dir, "04_FinalTranscriptomeFilter_Monocle_initial", "Initial_Monocle3_percent_ribo.pdf"),
  width = 8,
  height = 6.5
)

p_initial_hepatic <- plot_cells(
  cds_initial,
  color_cells_by = "HepaticScore1",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = FALSE
) + ggtitle("Initial Monocle3 UMAP colored by HepaticScore1")

save_plot_pdf(
  p_initial_hepatic,
  file.path(out_dir, "04_FinalTranscriptomeFilter_Monocle_initial", "Initial_Monocle3_HepaticScore1.pdf"),
  width = 8,
  height = 6.5
)

initial_filter_meta <- as.data.frame(colData(cds_initial)) %>%
  tibble::rownames_to_column("cellID")

monocle_contamination_summary <- initial_filter_meta %>%
  dplyr::group_by(monocle_initial_cluster) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    median_percent_ribo = median(percent.ribo, na.rm = TRUE),
    mean_percent_ribo = mean(percent.ribo, na.rm = TRUE),
    median_HepaticScore1 = median(HepaticScore1, na.rm = TRUE),
    mean_HepaticScore1 = mean(HepaticScore1, na.rm = TRUE),
    fraction_Liver_site = mean(site == "Liver", na.rm = TRUE),
    fraction_CTC_site = mean(site == "CTC", na.rm = TRUE),
    median_EpithelialScore1 = median(EpithelialScore1, na.rm = TRUE),
    median_MesenchymalScore1 = median(MesenchymalScore1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(median_percent_ribo))

write.csv(
  monocle_contamination_summary,
  file.path(out_dir, "04_FinalTranscriptomeFilter_Monocle_initial", "Initial_Monocle_cluster_contamination_summary.csv"),
  row.names = FALSE
)

cat("第一轮 Monocle cluster contamination 汇总：\n")
print(monocle_contamination_summary)
cat("\n")

clusters_to_remove <- as.character(initial_monocle_clusters_to_remove)

clusters_removed_summary <- monocle_contamination_summary %>%
  dplyr::filter(monocle_initial_cluster %in% clusters_to_remove) %>%
  dplyr::mutate(
    removal_reason = "manual_remove_high_ribosomal_candidate_based_on_v3_result"
  )

if (nrow(clusters_removed_summary) == 0) {
  stop(
    "未在初轮 Monocle cluster 中找到设定要删除的 cluster：",
    paste(clusters_to_remove, collapse = ", "),
    "\n请检查 Initial_Monocle_cluster_contamination_summary.csv。"
  )
}

write.csv(
  clusters_removed_summary,
  file.path(out_dir, "04_FinalTranscriptomeFilter_Monocle_initial", "Clusters_selected_for_removal.csv"),
  row.names = FALSE
)

cat("固定计划删除的 Monocle 初轮 cluster：", paste(clusters_to_remove, collapse = ", "), "\n")
cat("删除 cluster 汇总：\n")
print(clusters_removed_summary)
cat("\n")

cells_keep <- initial_filter_meta$cellID[
  !(initial_filter_meta$monocle_initial_cluster %in% clusters_to_remove)
]

cells_remove <- setdiff(colnames(m1_seu), cells_keep)

m1_seu_clean <- subset(m1_seu, cells = cells_keep)

write.csv(
  data.frame(cellID = cells_remove),
  file.path(out_dir, "04_FinalTranscriptomeFilter_Monocle_initial", "Cells_removed_by_final_transcriptome_filters.csv"),
  row.names = FALSE
)

removal_cell_count <- data.frame(
  stage = c(
    "before_final_transcriptome_filters",
    "removed_manual_high_ribo_cluster4",
    "after_final_transcriptome_filters"
  ),
  n_cells = c(
    ncol(m1_seu),
    length(cells_remove),
    ncol(m1_seu_clean)
  )
)

write.csv(
  removal_cell_count,
  file.path(out_dir, "04_FinalTranscriptomeFilter_Monocle_initial", "Cell_number_before_after_final_transcriptome_filters.csv"),
  row.names = FALSE
)

cat("final transcriptome filters 前细胞数：", ncol(m1_seu), "\n")
cat("删除细胞数：", length(cells_remove), "\n")
cat("final transcriptome filters 后细胞数：", ncol(m1_seu_clean), "\n\n")

saveRDS(
  m1_seu_clean,
  file.path(out_dir, "RDS", "GSE173958_M1_after_final_transcriptome_filters_Seurat_prefinal.rds")
)

DefaultAssay(m1_seu_clean) <- "SCT"

m1_seu_clean <- RunPCA(m1_seu_clean, npcs = 50, verbose = FALSE)

p_elbow_clean <- ElbowPlot(m1_seu_clean, ndims = 50)
save_plot_pdf(
  p_elbow_clean,
  file.path(out_dir, "05_PostFilter_Seurat_clean", "PCA_elbowplot_postFilter_clean.pdf"),
  width = 7,
  height = 5
)

m1_seu_clean <- RunUMAP(
  m1_seu_clean,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)

m1_seu_clean <- FindNeighbors(
  m1_seu_clean,
  reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)

m1_seu_clean <- FindClusters(
  m1_seu_clean,
  resolution = 0.5,
  verbose = FALSE
)

p_clean_cluster <- DimPlot(
  m1_seu_clean,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("M1 after final transcriptome filters: Seurat clusters")

save_plot_pdf(
  p_clean_cluster,
  file.path(out_dir, "05_PostFilter_Seurat_clean", "UMAP_Seurat_clusters_postFilter_clean.pdf"),
  width = 8,
  height = 6.5
)

p_clean_site <- DimPlot(
  m1_seu_clean,
  reduction = "umap",
  group.by = "site"
) + ggtitle("M1 after final transcriptome filters: harvest site")

save_plot_pdf(
  p_clean_site,
  file.path(out_dir, "05_PostFilter_Seurat_clean", "UMAP_harvest_site_postFilter_clean.pdf"),
  width = 8,
  height = 6.5
)

p_clean_clone <- DimPlot(
  m1_seu_clean,
  reduction = "umap",
  group.by = "cloneID",
  label = FALSE
) + ggtitle("M1 after final transcriptome filters: cloneID")

save_plot_pdf(
  p_clean_clone,
  file.path(out_dir, "05_PostFilter_Seurat_clean", "UMAP_cloneID_postFilter_clean.pdf"),
  width = 10,
  height = 7.5
)

p_clean_scores <- FeaturePlot(
  m1_seu_clean,
  features = c("EpithelialScore1", "MesenchymalScore1", "HepaticScore1", "percent.ribo"),
  reduction = "umap",
  ncol = 2,
  order = TRUE
)

save_plot_pdf(
  p_clean_scores,
  file.path(out_dir, "05_PostFilter_Seurat_clean", "UMAP_scores_Epi_Mes_Hepatic_Ribo_postFilter_clean.pdf"),
  width = 12,
  height = 10
)

saveRDS(
  m1_seu_clean,
  file.path(out_dir, "RDS", "GSE173958_M1_after_final_transcriptome_filters_Seurat_clean_UMAP.rds")
)

DefaultAssay(m1_seu_clean) <- "RNA"

m1_seu_clean[["RNA"]] <- JoinLayers(m1_seu_clean[["RNA"]])

counts_mat_clean <- GetAssayData(
  m1_seu_clean,
  assay = "RNA",
  layer = "counts"
)

cell_metadata_clean <- m1_seu_clean@meta.data
gene_metadata_clean <- data.frame(
  gene_short_name = rownames(counts_mat_clean),
  row.names = rownames(counts_mat_clean)
)

cds_clean <- new_cell_data_set(
  expression_data = counts_mat_clean,
  cell_metadata = cell_metadata_clean,
  gene_metadata = gene_metadata_clean
)

rm(counts_mat_clean)
gc()

cds_clean <- preprocess_cds(
  cds_clean,
  num_dim = 20
)

cds_clean <- align_cds(
  cds_clean,
  alignment_group = "site",
  residual_model_formula_str = "~ S.Score + G2M.Score + percent.mt"
)

cds_clean <- reduce_dimension(
  cds_clean,
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

cds_clean <- cluster_cells(
  cds_clean,
  reduction_method = "UMAP"
)

cds_clean <- learn_graph(
  cds_clean,
  use_partition = TRUE,
  close_loop = FALSE
)

cds_clean_meta <- as.data.frame(colData(cds_clean))

epi_cut <- quantile(cds_clean_meta$EpithelialScore1, 0.99, na.rm = TRUE)
mes_cut <- quantile(cds_clean_meta$MesenchymalScore1, 0.50, na.rm = TRUE)

root_cells <- rownames(cds_clean_meta)[
  cds_clean_meta$EpithelialScore1 >= epi_cut &
    cds_clean_meta$MesenchymalScore1 <= mes_cut
]

if (length(root_cells) < 10) {
  warning("严格 root 条件得到细胞过少，改用 EpithelialScore top 1% 作为 root。")
  root_cells <- rownames(cds_clean_meta)[cds_clean_meta$EpithelialScore1 >= epi_cut]
}

cat("清理后自动选定 root cells 数量：", length(root_cells), "\n\n")

write.csv(
  data.frame(root_cells = root_cells),
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "pseudoEMT_root_cells_auto_selected_postFilter.csv"),
  row.names = FALSE
)

cds_clean <- order_cells(
  cds_clean,
  reduction_method = "UMAP",
  root_cells = root_cells
)

pt_clean <- pseudotime(cds_clean)

write.csv(
  data.frame(cell = names(pt_clean), pseudoEMT = as.numeric(pt_clean)),
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "pseudoEMT_values_per_cell_postFilter.csv"),
  row.names = FALSE
)

p_clean_cds_cluster <- plot_cells(
  cds_clean,
  color_cells_by = "cluster",
  label_cell_groups = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
) + ggtitle("Post-filter Monocle3 clusters with learned graph")

save_plot_pdf(
  p_clean_cds_cluster,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Monocle3_clusters_with_graph.pdf"),
  width = 8,
  height = 6.5
)

p_clean_cds_site <- plot_cells(
  cds_clean,
  color_cells_by = "site",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
) + ggtitle("Post-filter Monocle3 graph colored by harvest site")

save_plot_pdf(
  p_clean_cds_site,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Monocle3_site_with_graph.pdf"),
  width = 8,
  height = 6.5
)

p_clean_cds_clone <- plot_cells(
  cds_clean,
  color_cells_by = "cloneID",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
) + ggtitle("Post-filter Monocle3 graph colored by cloneID")

save_plot_pdf(
  p_clean_cds_clone,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Monocle3_cloneID_with_graph.pdf"),
  width = 10,
  height = 7.5
)

p_clean_cds_pt <- plot_cells(
  cds_clean,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
) + ggtitle("Post-filter pseudoEMT / pseudotime")

save_plot_pdf(
  p_clean_cds_pt,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Monocle3_pseudoEMT_pseudotime.pdf"),
  width = 8,
  height = 6.5
)

m1_seu_clean$pseudoEMT <- NA_real_
m1_seu_clean$pseudoEMT[names(pt_clean)] <- as.numeric(pt_clean)

m1_seu_clean$pseudoEMT[!is.finite(m1_seu_clean$pseudoEMT)] <- NA_real_

cat("清理后可用 finite pseudoEMT 的细胞数：", sum(is.finite(m1_seu_clean$pseudoEMT)), "\n")
cat("清理后 pseudoEMT 为 NA 的细胞数：", sum(is.na(m1_seu_clean$pseudoEMT)), "\n\n")

DefaultAssay(m1_seu_clean) <- "SCT"

p_clean_seurat_pt <- FeaturePlot(
  m1_seu_clean,
  features = "pseudoEMT",
  reduction = "umap",
  order = TRUE
) + ggtitle("Post-filter Seurat UMAP colored by pseudoEMT")

save_plot_pdf(
  p_clean_seurat_pt,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Seurat_UMAP_pseudoEMT.pdf"),
  width = 8,
  height = 6.5
)

meta_plot_clean <- m1_seu_clean@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::filter(is.finite(pseudoEMT))

p_epi_curve_clean <- ggplot(meta_plot_clean, aes(x = pseudoEMT, y = EpithelialScore1)) +
  geom_point(alpha = 0.12, size = 0.3) +
  geom_smooth(se = FALSE, linewidth = 0.9) +
  labs(
    title = "Post-filter epithelial score along pseudoEMT",
    x = "pseudoEMT",
    y = "EpithelialScore1"
  ) +
  theme_classic()

save_plot_pdf(
  p_epi_curve_clean,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Epithelial_score_along_pseudoEMT.pdf"),
  width = 8,
  height = 5.5
)

p_mes_curve_clean <- ggplot(meta_plot_clean, aes(x = pseudoEMT, y = MesenchymalScore1)) +
  geom_point(alpha = 0.12, size = 0.3) +
  geom_smooth(se = FALSE, linewidth = 0.9) +
  labs(
    title = "Post-filter mesenchymal score along pseudoEMT",
    x = "pseudoEMT",
    y = "MesenchymalScore1"
  ) +
  theme_classic()

save_plot_pdf(
  p_mes_curve_clean,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Mesenchymal_score_along_pseudoEMT.pdf"),
  width = 8,
  height = 5.5
)

p_scores_curve_clean <- p_epi_curve_clean / p_mes_curve_clean

save_plot_pdf(
  p_scores_curve_clean,
  file.path(out_dir, "06_PostFilter_Monocle3_pseudoEMT", "PostFilter_Epi_Mes_scores_along_pseudoEMT_combined.pdf"),
  width = 8,
  height = 10
)

finite_pt_clean <- m1_seu_clean$pseudoEMT[
  is.finite(m1_seu_clean$pseudoEMT)
]

pt_q_clean <- quantile(
  finite_pt_clean,
  probs = c(0.25, 0.50, 0.75),
  na.rm = TRUE
)

m1_seu_clean$EMT_state_initial <- NA_character_

m1_seu_clean$EMT_state_initial[
  is.finite(m1_seu_clean$pseudoEMT) &
    m1_seu_clean$pseudoEMT <= pt_q_clean[1]
] <- "Epithelial"

m1_seu_clean$EMT_state_initial[
  is.finite(m1_seu_clean$pseudoEMT) &
    m1_seu_clean$pseudoEMT > pt_q_clean[1] &
    m1_seu_clean$pseudoEMT <= pt_q_clean[2]
] <- "Early_hybrid"

m1_seu_clean$EMT_state_initial[
  is.finite(m1_seu_clean$pseudoEMT) &
    m1_seu_clean$pseudoEMT > pt_q_clean[2] &
    m1_seu_clean$pseudoEMT <= pt_q_clean[3]
] <- "Late_hybrid"

m1_seu_clean$EMT_state_initial[
  is.finite(m1_seu_clean$pseudoEMT) &
    m1_seu_clean$pseudoEMT > pt_q_clean[3]
] <- "Mesenchymal"

m1_seu_clean$EMT_state_initial <- factor(
  m1_seu_clean$EMT_state_initial,
  levels = c("Epithelial", "Early_hybrid", "Late_hybrid", "Mesenchymal")
)

write.csv(
  data.frame(
    quantile = c("25%", "50%", "75%"),
    pseudoEMT_cutoff = as.numeric(pt_q_clean)
  ),
  file.path(out_dir, "08_EMT_state_initial_annotation", "PostFilter_pseudoEMT_quartile_cutoffs.csv"),
  row.names = FALSE
)

state_counts_clean <- as.data.frame(table(m1_seu_clean$EMT_state_initial))
colnames(state_counts_clean) <- c("EMT_state_initial", "cell_number")

write.csv(
  state_counts_clean,
  file.path(out_dir, "08_EMT_state_initial_annotation", "PostFilter_EMT_state_initial_cell_counts.csv"),
  row.names = FALSE
)

cat("清理后初版 EMT state 细胞数：\n")
print(state_counts_clean)
cat("\n")

p_state_umap_clean <- DimPlot(
  m1_seu_clean,
  reduction = "umap",
  group.by = "EMT_state_initial",
  label = FALSE
) + ggtitle("Post-filter initial EMT state annotation by pseudoEMT quartiles")

save_plot_pdf(
  p_state_umap_clean,
  file.path(out_dir, "08_EMT_state_initial_annotation", "PostFilter_UMAP_EMT_state_initial.pdf"),
  width = 8,
  height = 6.5
)

p_state_pt_clean <- VlnPlot(
  m1_seu_clean,
  features = "pseudoEMT",
  group.by = "EMT_state_initial",
  pt.size = 0.02
) + ggtitle("Post-filter pseudoEMT distribution across initial EMT states")

save_plot_pdf(
  p_state_pt_clean,
  file.path(out_dir, "08_EMT_state_initial_annotation", "PostFilter_Violin_pseudoEMT_by_initial_EMT_state.pdf"),
  width = 8,
  height = 5.5
)

p_state_score_clean <- VlnPlot(
  m1_seu_clean,
  features = c("EpithelialScore1", "MesenchymalScore1"),
  group.by = "EMT_state_initial",
  pt.size = 0.02,
  ncol = 2
)

save_plot_pdf(
  p_state_score_clean,
  file.path(out_dir, "08_EMT_state_initial_annotation", "PostFilter_Violin_Epi_Mes_scores_by_initial_EMT_state.pdf"),
  width = 12,
  height = 5.5
)

saveRDS(
  m1_seu_clean,
  file.path(
    out_dir,
    "RDS",
    "GSE173958_M1_checkpoint_step25_Seurat_with_pseudoEMT_EMTstate.rds"
  )
)

saveRDS(
  cds_clean,
  file.path(
    out_dir,
    "RDS",
    "GSE173958_M1_checkpoint_step25_Monocle3_cds_with_pseudoEMT.rds"
  )
)

cat("\n============================================================\n")
cat("Mac checkpoint 已保存，可以把整个分析文件夹复制到 Windows 后从 step26 继续。\n")
cat("Seurat checkpoint:\n")
cat(file.path(out_dir, "RDS", "GSE173958_M1_checkpoint_step25_Seurat_with_pseudoEMT_EMTstate.rds"), "\n")
cat("Monocle3 checkpoint:\n")
cat(file.path(out_dir, "RDS", "GSE173958_M1_checkpoint_step25_Monocle3_cds_with_pseudoEMT.rds"), "\n")
cat("Analysis stopped before step26 at:", format(Sys.time()), "\n")
cat("============================================================\n")

sink()

stop("已按设置在 step25 后停止。请把输出文件夹复制到 Windows，并运行 windows_resume_step26_dynamicEMT_DEG.R。")

if (run_dynamic_gene_analysis) {

graph_res_clean <- graph_test(
  cds_clean,
  neighbor_graph = "principal_graph",
  cores = 1
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

top_n_dynamic_clean <- min(3000, nrow(dynamic_res_clean))
top_dynamic_clean <- dynamic_res_clean %>% head(top_n_dynamic_clean)
top_dynamic_genes_clean <- top_dynamic_clean$gene

write.csv(
  top_dynamic_clean,
  file.path(out_dir, "07_pseudoEMT_dynamic_genes", paste0("PostFilter_pseudoEMT_dynamic_genes_top", top_n_dynamic_clean, ".csv")),
  row.names = FALSE
)

cat("清理后 pseudoEMT 动态基因满足 q<0.05 & Moran's I>0.1 的数量：", nrow(dynamic_res_clean), "\n")
cat("清理后用于动态热图的 top genes 数量：", top_n_dynamic_clean, "\n\n")

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
  n_bins = 100
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
    "PostFilter_pseudoEMT_top_dynamic_genes_six_modules_E_H1_H2_H3_H4_M_FORCED_ORDER.csv"
  ),
  row.names = FALSE
)

write.csv(
  module_name_map_clean,
  file.path(
    out_dir,
    "07_pseudoEMT_dynamic_genes",
    "PostFilter_pseudoEMT_six_module_peak_order_mapping_FORCED_ORDER.csv"
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
      "_dynamic_genes_100bins_FORCED_E_H1_H2_H3_H4_M_order.pdf"
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
    "Post-filter pseudoEMT dynamic genes top ",
    length(common_dynamic_genes_clean),
    " ordered as E-H1-H2-H3-H4-M"
  )
)

dev.off()

} else {
  cat("\n跳过 pseudoEMT 动态基因分析 graph_test()：run_dynamic_gene_analysis = FALSE。\n")
  cat("主分析仍会继续输出清理后 pseudoEMT、EMT state 初版注释和 hybrid EMT DEG。\n\n")
}

DefaultAssay(m1_seu_clean) <- "RNA"
m1_seu_clean[["RNA"]] <- JoinLayers(m1_seu_clean[["RNA"]])

m1_seu_clean <- NormalizeData(
  m1_seu_clean,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

m1_seu_clean$Hybrid_all_vs_Others <- ifelse(
  m1_seu_clean$EMT_state_initial %in% c("Early_hybrid", "Late_hybrid"),
  "Hybrid_all",
  "Others"
)
m1_seu_clean$Hybrid_all_vs_Others <- factor(
  m1_seu_clean$Hybrid_all_vs_Others,
  levels = c("Hybrid_all", "Others")
)

m1_seu_clean$Late_vs_Epi <- ifelse(
  m1_seu_clean$EMT_state_initial == "Late_hybrid",
  "Late_hybrid",
  ifelse(m1_seu_clean$EMT_state_initial == "Epithelial", "Epithelial", NA)
)
m1_seu_clean$Late_vs_Epi <- factor(
  m1_seu_clean$Late_vs_Epi,
  levels = c("Late_hybrid", "Epithelial")
)

m1_seu_clean$Late_vs_Mes <- ifelse(
  m1_seu_clean$EMT_state_initial == "Late_hybrid",
  "Late_hybrid",
  ifelse(m1_seu_clean$EMT_state_initial == "Mesenchymal", "Mesenchymal", NA)
)
m1_seu_clean$Late_vs_Mes <- factor(
  m1_seu_clean$Late_vs_Mes,
  levels = c("Late_hybrid", "Mesenchymal")
)

Idents(m1_seu_clean) <- "Hybrid_all_vs_Others"

deg_hybrid_vs_others_clean <- FindMarkers(
  object = m1_seu_clean,
  assay = "RNA",
  ident.1 = "Hybrid_all",
  ident.2 = "Others",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  test.use = "wilcox",
  slot = "data",
  densify = FALSE,
  max.cells.per.ident = 2000,
  random.seed = 1234
)

write.csv(
  deg_hybrid_vs_others_clean,
  file.path(out_dir, "09_DEG_Hybrid_EMT", "PostFilter_RNA_DEG_Hybrid_all_vs_Others.csv")
)

m1_late_epi_clean <- subset(m1_seu_clean, subset = !is.na(Late_vs_Epi))
DefaultAssay(m1_late_epi_clean) <- "RNA"
Idents(m1_late_epi_clean) <- "Late_vs_Epi"

deg_late_vs_epi_clean <- FindMarkers(
  object = m1_late_epi_clean,
  assay = "RNA",
  ident.1 = "Late_hybrid",
  ident.2 = "Epithelial",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  test.use = "wilcox",
  slot = "data",
  densify = FALSE,
  max.cells.per.ident = 2000,
  random.seed = 1234
)

write.csv(
  deg_late_vs_epi_clean,
  file.path(out_dir, "09_DEG_Hybrid_EMT", "PostFilter_RNA_DEG_Late_hybrid_vs_Epithelial.csv")
)

m1_late_mes_clean <- subset(m1_seu_clean, subset = !is.na(Late_vs_Mes))
DefaultAssay(m1_late_mes_clean) <- "RNA"
Idents(m1_late_mes_clean) <- "Late_vs_Mes"

deg_late_vs_mes_clean <- FindMarkers(
  object = m1_late_mes_clean,
  assay = "RNA",
  ident.1 = "Late_hybrid",
  ident.2 = "Mesenchymal",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  test.use = "wilcox",
  slot = "data",
  densify = FALSE,
  max.cells.per.ident = 2000,
  random.seed = 1234
)

write.csv(
  deg_late_vs_mes_clean,
  file.path(out_dir, "09_DEG_Hybrid_EMT", "PostFilter_RNA_DEG_Late_hybrid_vs_Mesenchymal.csv")
)

genes_hybrid_heat_clean <- top_up_down(deg_hybrid_vs_others_clean, n_each = 50)
genes_hybrid_heat_clean <- intersect(genes_hybrid_heat_clean, rownames(m1_seu_clean))

if (length(genes_hybrid_heat_clean) >= 2) {
  DefaultAssay(m1_seu_clean) <- "RNA"
  m1_seu_clean <- ScaleData(
    m1_seu_clean,
    assay = "RNA",
    features = genes_hybrid_heat_clean,
    verbose = FALSE
  )
  Idents(m1_seu_clean) <- "Hybrid_all_vs_Others"

  safe_doheatmap(
    seu_obj = m1_seu_clean,
    features = genes_hybrid_heat_clean,
    group_by = "Hybrid_all_vs_Others",
    title_text = "Post-filter RNA DEG: Hybrid_all vs Others, top 50 up + top 50 down",
    pdf_file = file.path(
      out_dir,
      "09_DEG_Hybrid_EMT",
      "Heatmaps",
      "PostFilter_RNA_Heatmap_Hybrid_all_vs_Others_top50UpDown.pdf"
    ),
    width = 12,
    height = 18
  )
}

genes_late_epi_heat_clean <- top_up_down(deg_late_vs_epi_clean, n_each = 50)
genes_late_epi_heat_clean <- intersect(genes_late_epi_heat_clean, rownames(m1_late_epi_clean))

if (length(genes_late_epi_heat_clean) >= 2) {
  DefaultAssay(m1_late_epi_clean) <- "RNA"
  m1_late_epi_clean <- ScaleData(
    m1_late_epi_clean,
    assay = "RNA",
    features = genes_late_epi_heat_clean,
    verbose = FALSE
  )
  Idents(m1_late_epi_clean) <- "Late_vs_Epi"

  safe_doheatmap(
    seu_obj = m1_late_epi_clean,
    features = genes_late_epi_heat_clean,
    group_by = "Late_vs_Epi",
    title_text = "Post-filter RNA DEG: Late_hybrid vs Epithelial, top 50 up + top 50 down",
    pdf_file = file.path(
      out_dir,
      "09_DEG_Hybrid_EMT",
      "Heatmaps",
      "PostFilter_RNA_Heatmap_Late_hybrid_vs_Epithelial_top50UpDown.pdf"
    ),
    width = 12,
    height = 18
  )
}

genes_late_mes_heat_clean <- top_up_down(deg_late_vs_mes_clean, n_each = 50)
genes_late_mes_heat_clean <- intersect(genes_late_mes_heat_clean, rownames(m1_late_mes_clean))

if (length(genes_late_mes_heat_clean) >= 2) {
  DefaultAssay(m1_late_mes_clean) <- "RNA"
  m1_late_mes_clean <- ScaleData(
    m1_late_mes_clean,
    assay = "RNA",
    features = genes_late_mes_heat_clean,
    verbose = FALSE
  )
  Idents(m1_late_mes_clean) <- "Late_vs_Mes"

  safe_doheatmap(
    seu_obj = m1_late_mes_clean,
    features = genes_late_mes_heat_clean,
    group_by = "Late_vs_Mes",
    title_text = "Post-filter RNA DEG: Late_hybrid vs Mesenchymal, top 50 up + top 50 down",
    pdf_file = file.path(
      out_dir,
      "09_DEG_Hybrid_EMT",
      "Heatmaps",
      "PostFilter_RNA_Heatmap_Late_hybrid_vs_Mesenchymal_top50UpDown.pdf"
    ),
    width = 12,
    height = 18
  )
}

saveRDS(
  m1_seu_clean,
  file.path(
    out_dir,
    "RDS",
    "GSE173958_M1_finalTranscriptomeFilter_with_pseudoEMT_initial_EMT_state_DEG_ready.rds"
  )
)

saveRDS(
  cds_initial,
  file.path(
    out_dir,
    "RDS",
    "GSE173958_M1_initial_Monocle3_for_finalTranscriptomeFilter.rds"
  )
)

saveRDS(
  cds_clean,
  file.path(
    out_dir,
    "RDS",
    "GSE173958_M1_postFilter_Monocle3_cds_pseudoEMT.rds"
  )
)

cat("\n============================================================\n")
cat("全部分析完成。\n")
cat("Analysis finished at:", format(Sys.time()), "\n")
cat("结果输出文件夹：", out_dir, "\n")
cat("============================================================\n")

sink()
