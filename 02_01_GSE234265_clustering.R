# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(clustree)

base_dir <- "path/to/GSE234265"

data_dir <- file.path(base_dir, "control_2")

out_dir <- file.path(base_dir, "control_2_reanalysis_V2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

setwd(out_dir)

get_counts <- function(obj, assay = "RNA") {
  counts_mat <- tryCatch(
    {
      GetAssayData(obj, assay = assay, layer = "counts")
    },
    error = function(e) {
      GetAssayData(obj, assay = assay, slot = "counts")
    }
  )
  return(counts_mat)
}

save_plot <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height
  )
}

message("Reading control_2 10X matrix...")

counts <- Read10X(data.dir = data_dir)

sc <- CreateSeuratObject(
  counts = counts,
  project = "GSE234265_control_2",
  min.cells = 3,
  min.features = 200
)

sc$sample_id <- "control_2"
sc$condition <- "control"

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")

sc
table(sc$sample_id)
table(sc$condition)

p_qc_before <- VlnPlot(
  sc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0
) +
  ggtitle("QC before filtering: control_2")

print(p_qc_before)

save_plot(
  p_qc_before,
  "01_QC_violin_before_filtering.pdf",
  width = 12,
  height = 5
)

n_cells_before <- ncol(sc)
message("Cells before QC filtering: ", n_cells_before)

sc <- subset(
  sc,
  subset =
    nFeature_RNA > 600 &
    nFeature_RNA < 8000 &
    percent.mt < 25
)

n_cells_after <- ncol(sc)
message("Cells after QC filtering: ", n_cells_after)

qc_summary <- data.frame(
  dataset = "control_2",
  cells_before_QC = n_cells_before,
  cells_after_QC = n_cells_after
)

write.csv(
  qc_summary,
  file = "02_QC_cell_number_summary.csv",
  row.names = FALSE
)

p_qc_after <- VlnPlot(
  sc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0
) +
  ggtitle("QC after filtering: control_2")

print(p_qc_after)

save_plot(
  p_qc_after,
  "03_QC_violin_after_filtering.pdf",
  width = 12,
  height = 5
)

sc <- NormalizeData(
  sc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

sc <- FindVariableFeatures(
  sc,
  selection.method = "vst",
  nfeatures = 2000
)

sc <- ScaleData(
  sc,
  features = VariableFeatures(sc)
)

sc <- RunPCA(
  sc,
  features = VariableFeatures(sc),
  npcs = 30
)

p_elbow_global <- ElbowPlot(sc, ndims = 30) +
  ggtitle("Global PCA ElbowPlot: control_2")

print(p_elbow_global)

save_plot(
  p_elbow_global,
  "04_Global_ElbowPlot.pdf",
  width = 6,
  height = 5
)

set.seed(1234)

sc <- JackStraw(
  sc,
  num.replicate = 100,
  dims = 30
)

sc <- ScoreJackStraw(
  sc,
  dims = 1:30
)

p_jack_global <- JackStrawPlot(
  sc,
  dims = 1:30
) +
  ggtitle("Global JackStrawPlot: control_2")

print(p_jack_global)

save_plot(
  p_jack_global,
  "05_Global_JackStrawPlot.pdf",
  width = 8,
  height = 6
)

dims_use_global <- 1:30

sc <- RunUMAP(
  sc,
  reduction = "pca",
  dims = dims_use_global,
  seed.use = 1234
)

sc <- FindNeighbors(
  sc,
  reduction = "pca",
  dims = dims_use_global
)

for (res in seq(0.1, 1.0, by = 0.1)) {
  sc <- FindClusters(
    sc,
    resolution = res,
    algorithm = 4,
    random.seed = 1234
  )
}

Idents(sc) <- "RNA_snn_res.0.1"
sc$global_cluster_res0.1 <- Idents(sc)

p_global_cluster <- DimPlot(
  sc,
  reduction = "umap",
  group.by = "global_cluster_res0.1",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Global clustering of control_2, Leiden res = 0.1")

print(p_global_cluster)

save_plot(
  p_global_cluster,
  "06_Global_UMAP_Leiden_res0.1.pdf",
  width = 7,
  height = 6
)

p_global_sample <- DimPlot(
  sc,
  reduction = "umap",
  group.by = "sample_id"
) +
  ggtitle("Sample identity: control_2")

save_plot(
  p_global_sample,
  "07_Global_UMAP_sample.pdf",
  width = 7,
  height = 6
)

p_global_clustree <- clustree(
  sc,
  prefix = "RNA_snn_res."
)

print(p_global_clustree)

ggsave(
  "08_Global_clustree_res0.1_to_1.0.pdf",
  plot = p_global_clustree,
  width = 12,
  height = 10
)

saveRDS(
  sc,
  file = "GSE234265_control_2_global_after_QC_clustering.rds"
)

grep(
  "YFP|yfp|EYFP|eyfp|Venus|venus|GFP|gfp",
  rownames(sc),
  value = TRUE
)

yfp_gene <- "YFP"

if (!yfp_gene %in% rownames(sc)) {
  stop("YFP gene not found. Please check the exact transgene name using grep().")
}

counts_mat <- get_counts(sc, assay = "RNA")

yfp_counts <- counts_mat[yfp_gene, ]

sc$YFP_count <- as.numeric(yfp_counts)
sc$YFP_positive <- sc$YFP_count >= 2

table(sc$YFP_positive)
prop.table(table(sc$YFP_positive))

yfp_summary <- data.frame(
  YFP_status = names(table(sc$YFP_positive)),
  cell_number = as.numeric(table(sc$YFP_positive)),
  proportion = as.numeric(prop.table(table(sc$YFP_positive)))
)

write.csv(
  yfp_summary,
  file = "09_YFP_positive_cell_summary.csv",
  row.names = FALSE
)

yfp_by_global_cluster <- as.data.frame.matrix(
  table(sc$global_cluster_res0.1, sc$YFP_positive)
)

yfp_by_global_cluster$global_cluster_res0.1 <- rownames(yfp_by_global_cluster)

write.csv(
  yfp_by_global_cluster,
  file = "10_YFP_positive_by_global_cluster.csv",
  row.names = FALSE
)

p_yfp_feature <- FeaturePlot(
  sc,
  features = yfp_gene,
  reduction = "umap"
) +
  ggtitle("YFP expression on global UMAP")

p_yfp_status <- DimPlot(
  sc,
  reduction = "umap",
  group.by = "YFP_positive"
) +
  ggtitle("YFP-positive cells on global UMAP")

p_global_cluster_again <- DimPlot(
  sc,
  reduction = "umap",
  group.by = "global_cluster_res0.1",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Global Leiden clusters, res = 0.1")

print(p_yfp_feature)
print(p_yfp_status)
print(p_global_cluster_again)

save_plot(
  p_yfp_feature,
  "11_Global_YFP_FeaturePlot.pdf",
  width = 7,
  height = 6
)

save_plot(
  p_yfp_status,
  "12_Global_YFP_positive_status.pdf",
  width = 7,
  height = 6
)

save_plot(
  p_global_cluster_again,
  "13_Global_clusters_for_YFP_check.pdf",
  width = 7,
  height = 6
)

saveRDS(
  sc,
  file = "GSE234265_control_2_global_with_YFP_info.rds"
)

tumor_initial <- subset(
  sc,
  subset = YFP_positive == TRUE
)

tumor_initial
table(tumor_initial$sample_id)
table(tumor_initial$condition)

message("Initial YFP-positive cell number: ", ncol(tumor_initial))

tumor_initial <- NormalizeData(
  tumor_initial,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

tumor_initial <- FindVariableFeatures(
  tumor_initial,
  selection.method = "vst",
  nfeatures = 2000
)

tumor_initial <- ScaleData(
  tumor_initial,
  features = VariableFeatures(tumor_initial)
)

tumor_initial <- RunPCA(
  tumor_initial,
  features = VariableFeatures(tumor_initial),
  npcs = 30
)

p_elbow_tumor_initial <- ElbowPlot(
  tumor_initial,
  ndims = 30
) +
  ggtitle("Initial YFP+ tumor PCA ElbowPlot")

print(p_elbow_tumor_initial)

save_plot(
  p_elbow_tumor_initial,
  "14_Initial_YFP_tumor_ElbowPlot.pdf",
  width = 6,
  height = 5
)

set.seed(1234)

tumor_initial <- JackStraw(
  tumor_initial,
  num.replicate = 100,
  dims = 30
)

tumor_initial <- ScoreJackStraw(
  tumor_initial,
  dims = 1:30
)

p_jack_tumor_initial <- JackStrawPlot(
  tumor_initial,
  dims = 1:30
) +
  ggtitle("Initial YFP+ tumor JackStrawPlot")

print(p_jack_tumor_initial)

save_plot(
  p_jack_tumor_initial,
  "15_Initial_YFP_tumor_JackStrawPlot.pdf",
  width = 8,
  height = 6
)

dims_use_tumor_initial <- 1:30

tumor_initial <- RunUMAP(
  tumor_initial,
  reduction = "pca",
  dims = dims_use_tumor_initial,
  seed.use = 1234
)

tumor_initial <- FindNeighbors(
  tumor_initial,
  reduction = "pca",
  dims = dims_use_tumor_initial
)

for (res in seq(0.1, 1.0, by = 0.1)) {
  tumor_initial <- FindClusters(
    tumor_initial,
    resolution = res,
    algorithm = 4,
    random.seed = 1234
  )
}

Idents(tumor_initial) <- "RNA_snn_res.0.1"
tumor_initial$tumor_initial_cluster_res0.1 <- Idents(tumor_initial)

p_tumor_initial_res01 <- DimPlot(
  tumor_initial,
  reduction = "umap",
  group.by = "tumor_initial_cluster_res0.1",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Initial YFP+ clustering, Leiden res = 0.1")

print(p_tumor_initial_res01)

save_plot(
  p_tumor_initial_res01,
  "16_Initial_YFP_tumor_UMAP_Leiden_res0.1.pdf",
  width = 7,
  height = 6
)

p_tumor_initial_clustree <- clustree(
  tumor_initial,
  prefix = "RNA_snn_res."
)

print(p_tumor_initial_clustree)

ggsave(
  "17_Initial_YFP_tumor_clustree_res0.1_to_1.0.pdf",
  plot = p_tumor_initial_clustree,
  width = 12,
  height = 10
)

saveRDS(
  tumor_initial,
  file = "GSE234265_control_2_YFP_tumor_initial_clustering.rds"
)

non_tumor_modules <- list(
  Immune = c(
    "Ptprc", "Lyz2", "Cd3d", "Cd3e", "Nkg7", "Ms4a1", "Cd79a"
  ),
  Endothelial = c(
    "Pecam1", "Kdr", "Cdh5", "Emcn", "Klf2", "Robo4"
  ),
  Fibroblast_CAF = c(
    "Col1a1", "Col1a2", "Col3a1", "Col6a3", "Dcn", "Lum", "Pdgfra"
  ),
  Pericyte = c(
    "Rgs5", "Cspg4", "Des", "Notch3", "Pdgfrb", "Acta2"
  ),
  RBC_like = c(
    "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Alas2"
  )
)

non_tumor_modules_filtered <- lapply(
  non_tumor_modules,
  function(gene_vec) {
    gene_vec[gene_vec %in% rownames(tumor_initial)]
  }
)

non_tumor_modules_filtered

non_tumor_modules_filtered <- non_tumor_modules_filtered[
  lengths(non_tumor_modules_filtered) > 0
]

tumor_initial <- AddModuleScore(
  object = tumor_initial,
  features = non_tumor_modules_filtered,
  name = "NonTumorModule"
)

module_score_cols <- paste0(
  "NonTumorModule",
  seq_along(non_tumor_modules_filtered)
)

module_name_map <- data.frame(
  score_column = module_score_cols,
  module_name = names(non_tumor_modules_filtered)
)

write.csv(
  module_name_map,
  file = "18_NonTumor_module_score_column_mapping.csv",
  row.names = FALSE
)

Idents(tumor_initial) <- "tumor_initial_cluster_res0.1"

p_module_feature <- FeaturePlot(
  tumor_initial,
  features = module_score_cols,
  reduction = "umap",
  ncol = 3
)

print(p_module_feature)

ggsave(
  "19_Initial_YFP_NonTumorModule_FeaturePlots.pdf",
  plot = p_module_feature,
  width = 14,
  height = 10
)

p_module_vln <- VlnPlot(
  tumor_initial,
  features = module_score_cols,
  group.by = "tumor_initial_cluster_res0.1",
  pt.size = 0,
  ncol = 2
)

print(p_module_vln)

ggsave(
  "20_Initial_YFP_NonTumorModule_VlnPlots_by_cluster.pdf",
  plot = p_module_vln,
  width = 14,
  height = 10
)

module_score_df <- FetchData(
  tumor_initial,
  vars = c("tumor_initial_cluster_res0.1", module_score_cols)
)

cluster_mean_module_scores <- module_score_df %>%
  group_by(tumor_initial_cluster_res0.1) %>%
  summarise(
    across(
      all_of(module_score_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

write.csv(
  cluster_mean_module_scores,
  file = "21_Initial_YFP_cluster_mean_non_tumor_module_scores.csv",
  row.names = FALSE
)

cluster_mean_module_scores

Idents(tumor_initial) <- "tumor_initial_cluster_res0.1"

markers_tumor_initial <- FindAllMarkers(
  tumor_initial,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers_tumor_initial,
  file = "22_Initial_YFP_cluster_marker_genes_res0.1.csv",
  row.names = FALSE
)

top20_initial_markers <- markers_tumor_initial %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()

write.csv(
  top20_initial_markers,
  file = "23_Initial_YFP_cluster_top20_markers_res0.1.csv",
  row.names = FALSE
)

clusters_to_remove <- c()

if (length(clusters_to_remove) > 0) {
  
  message(
    "Removing suspected contaminating initial YFP+ clusters: ",
    paste(clusters_to_remove, collapse = ", ")
  )
  
  tumor_clean <- subset(
    tumor_initial,
    subset = !(tumor_initial_cluster_res0.1 %in% clusters_to_remove)
  )
  
} else {
  
  message(
    "No contaminating clusters removed. ",
    "tumor_clean will be identical to tumor_initial. ",
    "For article-like analysis, please inspect diagnostic files and update clusters_to_remove if needed."
  )
  
  tumor_clean <- tumor_initial
}

message("YFP+ cell number after contamination removal: ", ncol(tumor_clean))

contamination_removal_summary <- data.frame(
  initial_YFP_positive_cells = ncol(tumor_initial),
  clusters_removed = paste(clusters_to_remove, collapse = ";"),
  remaining_cells_after_removal = ncol(tumor_clean)
)

write.csv(
  contamination_removal_summary,
  file = "24_Contamination_removal_summary.csv",
  row.names = FALSE
)

saveRDS(
  tumor_clean,
  file = "GSE234265_control_2_YFP_tumor_after_contamination_removal_before_rerun.rds"
)

tumor_final <- tumor_clean

tumor_final <- NormalizeData(
  tumor_final,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

tumor_final <- FindVariableFeatures(
  tumor_final,
  selection.method = "vst",
  nfeatures = 2000
)

tumor_final <- ScaleData(
  tumor_final,
  features = VariableFeatures(tumor_final)
)

tumor_final <- RunPCA(
  tumor_final,
  features = VariableFeatures(tumor_final),
  npcs = 30
)

p_elbow_tumor_final <- ElbowPlot(
  tumor_final,
  ndims = 30
) +
  ggtitle("Final YFP+ tumor PCA ElbowPlot after contamination removal")

print(p_elbow_tumor_final)

save_plot(
  p_elbow_tumor_final,
  "25_Final_YFP_tumor_ElbowPlot_after_removal.pdf",
  width = 6,
  height = 5
)

set.seed(1234)

tumor_final <- JackStraw(
  tumor_final,
  num.replicate = 100,
  dims = 30
)

tumor_final <- ScoreJackStraw(
  tumor_final,
  dims = 1:30
)

p_jack_tumor_final <- JackStrawPlot(
  tumor_final,
  dims = 1:30
) +
  ggtitle("Final YFP+ tumor JackStrawPlot after contamination removal")

print(p_jack_tumor_final)

save_plot(
  p_jack_tumor_final,
  "26_Final_YFP_tumor_JackStrawPlot_after_removal.pdf",
  width = 8,
  height = 6
)

dims_use_tumor_final <- 1:30

tumor_final <- RunUMAP(
  tumor_final,
  reduction = "pca",
  dims = dims_use_tumor_final,
  seed.use = 1234
)

tumor_final <- FindNeighbors(
  tumor_final,
  reduction = "pca",
  dims = dims_use_tumor_final
)

for (res in seq(0.1, 1.0, by = 0.1)) {
  tumor_final <- FindClusters(
    tumor_final,
    resolution = res,
    algorithm = 4,
    random.seed = 1234
  )
}

Idents(tumor_final) <- "RNA_snn_res.0.7"
tumor_final$tumor_final_cluster_res0.7 <- Idents(tumor_final)

p_tumor_final_res07 <- DimPlot(
  tumor_final,
  reduction = "umap",
  group.by = "tumor_final_cluster_res0.7",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Final YFP+ tumor clustering after removal, Leiden res = 0.7")

print(p_tumor_final_res07)

save_plot(
  p_tumor_final_res07,
  "27_Final_YFP_tumor_UMAP_Leiden_res0.7.pdf",
  width = 7,
  height = 6
)

p_tumor_final_clustree <- clustree(
  tumor_final,
  prefix = "RNA_snn_res."
)

print(p_tumor_final_clustree)

ggsave(
  "28_Final_YFP_tumor_clustree_res0.1_to_1.0.pdf",
  plot = p_tumor_final_clustree,
  width = 12,
  height = 10
)

Idents(tumor_final) <- "tumor_final_cluster_res0.7"

emt_markers <- c(
  "Epcam",
  "Krt14",
  "Krt8",
  "Vim",
  "Pdgfra"
)

emt_markers_present <- emt_markers[emt_markers %in% rownames(tumor_final)]
emt_markers_missing <- setdiff(emt_markers, rownames(tumor_final))

emt_marker_check <- data.frame(
  EMT_marker = emt_markers,
  present_in_object = emt_markers %in% rownames(tumor_final)
)

write.csv(
  emt_marker_check,
  file = "29_EMT_marker_presence_check.csv",
  row.names = FALSE
)

print(emt_marker_check)

p_emt_feature_final <- FeaturePlot(
  tumor_final,
  features = emt_markers_present,
  reduction = "umap",
  ncol = 3
)

print(p_emt_feature_final)

ggsave(
  "30_Final_YFP_EMT_marker_FeaturePlots.pdf",
  plot = p_emt_feature_final,
  width = 12,
  height = 8
)

p_emt_dot_final <- DotPlot(
  tumor_final,
  features = emt_markers_present,
  group.by = "tumor_final_cluster_res0.7"
) +
  RotatedAxis() +
  ggtitle("Final YFP+ tumor clusters: EMT markers")

print(p_emt_dot_final)

save_plot(
  p_emt_dot_final,
  "31_Final_YFP_EMT_marker_DotPlot_by_cluster.pdf",
  width = 8,
  height = 5
)

p_emt_vln_final <- VlnPlot(
  tumor_final,
  features = emt_markers_present,
  group.by = "tumor_final_cluster_res0.7",
  pt.size = 0,
  ncol = 3
)

print(p_emt_vln_final)

ggsave(
  "32_Final_YFP_EMT_marker_VlnPlots_by_cluster.pdf",
  plot = p_emt_vln_final,
  width = 12,
  height = 8
)

avg_emt_final <- AverageExpression(
  tumor_final,
  features = emt_markers_present,
  group.by = "tumor_final_cluster_res0.7",
  assays = "RNA",
  layer = "data"
)$RNA

avg_emt_final_df <- as.data.frame(avg_emt_final)
avg_emt_final_df$gene <- rownames(avg_emt_final_df)

write.csv(
  avg_emt_final_df,
  file = "33_Final_YFP_EMT_marker_average_expression_by_cluster.csv",
  row.names = FALSE
)

avg_emt_final_df

Idents(tumor_final) <- "tumor_final_cluster_res0.7"

markers_tumor_final <- FindAllMarkers(
  tumor_final,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers_tumor_final,
  file = "34_Final_YFP_cluster_marker_genes_res0.7.csv",
  row.names = FALSE
)

top20_final_markers <- markers_tumor_final %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()

write.csv(
  top20_final_markers,
  file = "35_Final_YFP_cluster_top20_markers_res0.7.csv",
  row.names = FALSE
)

sort(unique(as.character(tumor_final$tumor_final_cluster_res0.7)))

emt_annotation_final <- c(
  "1" = "late_hybrid",
  "2" = "late_hybrid",
  "3" = "full_EMT",
  "4" = "epithelial",
  "5" = "early_hybrid",
  "6" = "full_EMT",
  "7" = "late_hybrid"
)

tumor_final$EMT_state_final <- unname(
  emt_annotation_final[
    as.character(tumor_final$tumor_final_cluster_res0.7)
  ]
)

table(
  tumor_final$tumor_final_cluster_res0.7,
  tumor_final$EMT_state_final,
  useNA = "ifany"
)

table(
  tumor_final$EMT_state_final,
  useNA = "ifany"
)

p_emt_state_final <- DimPlot(
  tumor_final,
  reduction = "umap",
  group.by = "EMT_state_final",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Final EMT-state annotation of control_2 YFP+ tumor cells")

print(p_emt_state_final)

ggsave(
  "36_Final_YFP_EMT_state_annotation_UMAP.pdf",
  p_emt_state_final,
  width = 8,
  height = 6
)

saveRDS(
  tumor_final,
  file = "GSE234265_control_2_YFP_tumor_final_EMT_annotated.rds"
)

message("All analyses finished.")
