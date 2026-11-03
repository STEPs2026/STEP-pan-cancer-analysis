# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(monocle)
library(patchwork)
library(igraph)

packageVersion("igraph")

base_dir <- "path/to/GSE234265"

input_rds <- file.path(
  base_dir,
  "control_2_reanalysis_V2",
  "GSE234265_control_2_YFP_tumor_final_EMT_annotated.rds"
)

out_dir <- file.path(
  base_dir,
  "control_2_reanalysis_V2",
  "pseudotime_monocle2_all_YFP_tumor_cells"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

tumor_final <- readRDS(input_rds)

tumor_final
dim(tumor_final)

colnames(tumor_final@meta.data)

if (!"EMT_state_final" %in% colnames(tumor_final@meta.data)) {
  stop("EMT_state_final not found. Please use the final EMT-annotated Seurat object.")
}

if (!"tumor_final_cluster_res0.7" %in% colnames(tumor_final@meta.data)) {
  stop("tumor_final_cluster_res0.7 not found.")
}

table(tumor_final$EMT_state_final)
table(tumor_final$tumor_final_cluster_res0.7)

DefaultAssay(tumor_final) <- "RNA"

if (!"Phase" %in% colnames(tumor_final@meta.data)) {
  
  s_genes <- Seurat::CaseMatch(
    search = cc.genes.updated.2019$s.genes,
    match = rownames(tumor_final)
  )
  
  g2m_genes <- Seurat::CaseMatch(
    search = cc.genes.updated.2019$g2m.genes,
    match = rownames(tumor_final)
  )
  
  length(s_genes)
  length(g2m_genes)
  
  tumor_final <- CellCycleScoring(
    object = tumor_final,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = FALSE
  )
}

table(tumor_final$Phase)

cell_cycle_summary <- as.data.frame(table(tumor_final$Phase))
colnames(cell_cycle_summary) <- c("Phase", "Cell_number")

write.csv(
  cell_cycle_summary,
  file = "01_CellCycle_phase_summary_all_YFP_tumor_cells.csv",
  row.names = FALSE
)

p_phase_umap <- DimPlot(
  tumor_final,
  reduction = "umap",
  group.by = "Phase"
) +
  ggtitle("Cell-cycle phase of all final YFP+ tumor cells")

print(p_phase_umap)

ggsave(
  "02_CellCycle_phase_on_Seurat_UMAP_all_cells.pdf",
  plot = p_phase_umap,
  width = 7,
  height = 6
)

tumor_all <- tumor_final

message("Number of all final YFP+ tumor cells used for pseudotime: ", ncol(tumor_all))

all_cell_summary <- data.frame(
  all_final_YFP_tumor_cells = ncol(tumor_all)
)

write.csv(
  all_cell_summary,
  file = "03_All_YFP_tumor_cell_number_summary.csv",
  row.names = FALSE
)

ordering_genes <- VariableFeatures(tumor_final)

length(ordering_genes)

ordering_genes <- ordering_genes[
  ordering_genes %in% rownames(tumor_all)
]

length(ordering_genes)

write.csv(
  data.frame(ordering_gene = ordering_genes),
  file = "04_Monocle2_ordering_genes_top2000_HVGs_all_cells.csv",
  row.names = FALSE
)

get_counts <- function(obj, assay = "RNA") {
  tryCatch(
    {
      GetAssayData(obj, assay = assay, layer = "counts")
    },
    error = function(e) {
      GetAssayData(obj, assay = assay, slot = "counts")
    }
  )
}

expression_matrix <- get_counts(tumor_all, assay = "RNA")
expression_matrix <- as(expression_matrix, "sparseMatrix")

pd <- tumor_all@meta.data
pd <- pd[colnames(expression_matrix), , drop = FALSE]

pheno_data <- new(
  "AnnotatedDataFrame",
  data = pd
)

feature_data <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix),
  stringsAsFactors = FALSE
)

feature_data <- new(
  "AnnotatedDataFrame",
  data = feature_data
)

cds_all <- newCellDataSet(
  cellData = expression_matrix,
  phenoData = pheno_data,
  featureData = feature_data,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

cds_all

cds_all <- estimateSizeFactors(cds_all)

ordering_genes <- ordering_genes[
  ordering_genes %in% rownames(fData(cds_all))
]

length(ordering_genes)

cds_all <- setOrderingFilter(
  cds_all,
  ordering_genes
)

set.seed(1234)

cds_all <- reduceDimension(
  cds_all,
  max_components = 2,
  method = "DDRTree"
)

patched_extract_ddrtree_ordering <- function(cds, root_cell, verbose = TRUE) {
  
  dp <- monocle::cellPairwiseDistances(cds)
  dp_mst <- monocle::minSpanningTree(cds)
  
  states <- rep(1, ncol(dp))
  names(states) <- igraph::V(dp_mst)$name
  
  pseudotimes <- rep(0, ncol(dp))
  names(pseudotimes) <- igraph::V(dp_mst)$name
  
  parents <- rep(NA, ncol(dp))
  names(parents) <- igraph::V(dp_mst)$name
  
  mst_traversal <- igraph::dfs(
    graph = dp_mst,
    root = root_cell,
    mode = "all",
    unreachable = FALSE,
    parent = TRUE
  )
  
  parent_vec <- as.numeric(mst_traversal$parent)
  order_vec <- as.numeric(mst_traversal$order)
  
  curr_state <- 1
  
  for (i in seq_along(order_vec)) {
    
    curr_node <- order_vec[i]
    
    if (is.na(curr_node)) {
      next
    }
    
    curr_node_name <- igraph::V(dp_mst)[curr_node]$name
    
    if (!is.na(parent_vec[curr_node])) {
      
      parent_node <- parent_vec[curr_node]
      parent_node_name <- igraph::V(dp_mst)[parent_node]$name
      
      parent_node_pseudotime <- pseudotimes[parent_node_name]
      
      curr_node_pseudotime <- parent_node_pseudotime +
        dp[curr_node_name, parent_node_name]
      
      if (igraph::degree(dp_mst, v = parent_node_name) > 2) {
        curr_state <- curr_state + 1
      }
      
    } else {
      
      parent_node_name <- NA
      curr_node_pseudotime <- 0
    }
    
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_state
    parents[curr_node_name] <- parent_node_name
  }
  
  ordering_df <- data.frame(
    sample_name = names(states),
    cell_state = factor(states),
    pseudo_time = as.vector(pseudotimes),
    parent = parents,
    stringsAsFactors = FALSE
  )
  
  rownames(ordering_df) <- ordering_df$sample_name
  
  return(ordering_df)
}

assignInNamespace(
  x = "extract_ddrtree_ordering",
  value = patched_extract_ddrtree_ordering,
  ns = "monocle"
)

patched_project2MST <- function(cds, Projection_Method) {
  
  dp_mst <- monocle::minSpanningTree(cds)
  Z <- monocle::reducedDimS(cds)
  Y <- monocle::reducedDimK(cds)
  
  cds <- monocle:::findNearestPointOnMST(cds)
  
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex]
  
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
  
  if (!is.function(Projection_Method)) {
    
    P <- Y[, closest_vertex]
    
  } else {
    
    P <- matrix(
      rep(0, length(Z)),
      nrow = nrow(Z)
    )
    
    for (i in 1:length(closest_vertex)) {
      
      neighbor_vs <- igraph::neighbors(
        dp_mst,
        v = closest_vertex_names[i],
        mode = "all"
      )
      
      neighbors <- igraph::as_ids(neighbor_vs)
      
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      
      for (neighbor in neighbors) {
        
        if (closest_vertex_names[i] %in% tip_leaves) {
          
          tmp <- monocle:::projPointOnLine(
            Z_i,
            Y[, c(closest_vertex_names[i], neighbor)]
          )
          
        } else {
          
          tmp <- Projection_Method(
            Z_i,
            Y[, c(closest_vertex_names[i], neighbor)]
          )
        }
        
        projection <- rbind(projection, tmp)
        distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
      }
      
      if (class(projection)[1] != "matrix") {
        projection <- as.matrix(projection)
      }
      
      P[, i] <- projection[
        which(distance == min(distance))[1],
      ]
    }
  }
  
  colnames(P) <- colnames(Z)
  
  dp <- as.matrix(stats::dist(t(P)))
  
  min_dist <- min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0
  
  monocle::cellPairwiseDistances(cds) <- dp
  
  gp <- igraph::graph.adjacency(
    dp,
    mode = "undirected",
    weighted = TRUE
  )
  
  dp_mst <- igraph::minimum.spanning.tree(gp)
  
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  
  return(cds)
}

assignInNamespace(
  x = "project2MST",
  value = patched_project2MST,
  ns = "monocle"
)

cds_all <- orderCells(cds_all)

table(pData(cds_all)$State)
summary(pData(cds_all)$Pseudotime)

state_summary_initial <- as.data.frame(table(pData(cds_all)$State))
colnames(state_summary_initial) <- c("Monocle_State", "Cell_number")

write.csv(
  state_summary_initial,
  file = "05_Monocle2_initial_state_summary_all_cells.csv",
  row.names = FALSE
)

root_state_table <- table(
  pData(cds_all)$State,
  pData(cds_all)$EMT_state_final
)

root_state_table

write.csv(
  as.data.frame.matrix(root_state_table),
  file = "06_Monocle_State_by_EMT_state_before_rooting_all_cells.csv"
)

if ("epithelial" %in% colnames(root_state_table)) {
  
  epithelial_counts_by_state <- root_state_table[, "epithelial"]
  
  root_state <- as.numeric(
    names(epithelial_counts_by_state)[which.max(epithelial_counts_by_state)]
  )
  
} else {
  
  stop("No epithelial cells found. Root state cannot be assigned automatically.")
}

root_state

write.csv(
  data.frame(root_state = root_state),
  file = "07_Selected_root_state_based_on_epithelial_enrichment_all_cells.csv",
  row.names = FALSE
)

cds_all <- orderCells(
  cds_all,
  root_state = root_state
)

summary(pData(cds_all)$Pseudotime)

pseudotime_summary <- data.frame(
  root_state = root_state,
  min_pseudotime = min(pData(cds_all)$Pseudotime, na.rm = TRUE),
  median_pseudotime = median(pData(cds_all)$Pseudotime, na.rm = TRUE),
  mean_pseudotime = mean(pData(cds_all)$Pseudotime, na.rm = TRUE),
  max_pseudotime = max(pData(cds_all)$Pseudotime, na.rm = TRUE)
)

write.csv(
  pseudotime_summary,
  file = "08_Monocle2_pseudotime_summary_all_cells.csv",
  row.names = FALSE
)

plot_cell_trajectory_modern <- function(
    cds,
    color_by = "State",
    show_tree = TRUE,
    show_branch_points = TRUE,
    cell_size = 1.5,
    cell_link_size = 0.75,
    title = NULL
) {
  
  if (!color_by %in% colnames(pData(cds))) {
    stop(paste0("'", color_by, "' not found in pData(cds)."))
  }

  principal_coords <- monocle::reducedDimK(cds)
  
  graph_df <- as.data.frame(t(principal_coords))
  graph_df$sample_name <- rownames(graph_df)
  
  colnames(graph_df)[1:2] <- c(
    "prin_graph_dim_1",
    "prin_graph_dim_2"
  )

  cell_coords <- as.data.frame(t(monocle::reducedDimS(cds)))
  cell_coords$sample_name <- rownames(cell_coords)
  
  colnames(cell_coords)[1:2] <- c(
    "data_dim_1",
    "data_dim_2"
  )
  
  meta_df <- pData(cds)
  meta_df$sample_name <- rownames(meta_df)
  
  data_df <- merge(
    cell_coords,
    meta_df,
    by = "sample_name",
    all.x = TRUE,
    sort = FALSE
  )
  
  data_df <- data_df[
    match(rownames(cell_coords), data_df$sample_name),
  ]

  dp_mst <- monocle::minSpanningTree(cds)
  
  edge_raw <- igraph::as_data_frame(
    dp_mst,
    what = "edges"
  )
  
  edge_df <- data.frame(
    source = edge_raw$from,
    target = edge_raw$to,
    stringsAsFactors = FALSE
  )
  
  source_df <- graph_df[, c(
    "sample_name",
    "prin_graph_dim_1",
    "prin_graph_dim_2"
  )]
  
  colnames(source_df) <- c(
    "source",
    "source_prin_graph_dim_1",
    "source_prin_graph_dim_2"
  )
  
  target_df <- graph_df[, c(
    "sample_name",
    "prin_graph_dim_1",
    "prin_graph_dim_2"
  )]
  
  colnames(target_df) <- c(
    "target",
    "target_prin_graph_dim_1",
    "target_prin_graph_dim_2"
  )
  
  edge_df <- merge(
    edge_df,
    source_df,
    by = "source",
    all.x = TRUE,
    sort = FALSE
  )
  
  edge_df <- merge(
    edge_df,
    target_df,
    by = "target",
    all.x = TRUE,
    sort = FALSE
  )
  
  g <- ggplot(
    data_df,
    aes(
      x = data_dim_1,
      y = data_dim_2
    )
  )
  
  if (show_tree) {
    g <- g +
      geom_segment(
        data = edge_df,
        aes(
          x = source_prin_graph_dim_1,
          y = source_prin_graph_dim_2,
          xend = target_prin_graph_dim_1,
          yend = target_prin_graph_dim_2
        ),
        inherit.aes = FALSE,
        linewidth = cell_link_size
      )
  }
  
  g <- g +
    geom_point(
      aes_string(color = color_by),
      size = cell_size,
      alpha = 0.9
    )

  if (
    show_branch_points &&
    !is.null(cds@auxOrderingData[["DDRTree"]]$branch_points)
  ) {
    
    branch_nodes <- cds@auxOrderingData[["DDRTree"]]$branch_points
    
    branch_df <- graph_df[
      match(branch_nodes, graph_df$sample_name),
      ,
      drop = FALSE
    ]
    
    branch_df <- branch_df[!is.na(branch_df$sample_name), , drop = FALSE]
    
    if (nrow(branch_df) > 0) {
      
      branch_df$branch_point_idx <- seq_len(nrow(branch_df))
      
      g <- g +
        geom_point(
          data = branch_df,
          aes(
            x = prin_graph_dim_1,
            y = prin_graph_dim_2
          ),
          inherit.aes = FALSE,
          size = 5
        ) +
        geom_text(
          data = branch_df,
          aes(
            x = prin_graph_dim_1,
            y = prin_graph_dim_2,
            label = branch_point_idx
          ),
          inherit.aes = FALSE,
          size = 4
        )
    }
  }
  
  g <- g +
    theme_classic() +
    xlab("Component 1") +
    ylab("Component 2") +
    theme(
      legend.position = "top"
    )
  
  if (!is.null(title)) {
    g <- g + ggtitle(title)
  }
  
  return(g)
}

p_traj_pseudotime_all <- plot_cell_trajectory_modern(
  cds_all,
  color_by = "Pseudotime",
  title = "Monocle2 DDRTree trajectory: Pseudotime, all YFP+ tumor cells"
)

p_traj_state_all <- plot_cell_trajectory_modern(
  cds_all,
  color_by = "State",
  title = "Monocle2 DDRTree trajectory: State, all YFP+ tumor cells"
)

p_traj_emt_all <- plot_cell_trajectory_modern(
  cds_all,
  color_by = "EMT_state_final",
  title = "Monocle2 DDRTree trajectory: EMT states, all YFP+ tumor cells"
)

p_traj_cluster_all <- plot_cell_trajectory_modern(
  cds_all,
  color_by = "tumor_final_cluster_res0.7",
  title = "Monocle2 DDRTree trajectory: final Leiden clusters, all cells"
)

p_traj_phase_all <- plot_cell_trajectory_modern(
  cds_all,
  color_by = "Phase",
  title = "Monocle2 DDRTree trajectory: Cell-cycle phase, all cells"
)

print(p_traj_pseudotime_all)
print(p_traj_state_all)
print(p_traj_emt_all)
print(p_traj_cluster_all)
print(p_traj_phase_all)

ggsave(
  "09_Monocle2_DDRTree_Pseudotime_all_cells.pdf",
  plot = p_traj_pseudotime_all,
  width = 7,
  height = 6
)

ggsave(
  "10_Monocle2_DDRTree_State_all_cells.pdf",
  plot = p_traj_state_all,
  width = 7,
  height = 6
)

ggsave(
  "11_Monocle2_DDRTree_EMT_state_all_cells.pdf",
  plot = p_traj_emt_all,
  width = 7,
  height = 6
)

ggsave(
  "12_Monocle2_DDRTree_final_clusters_all_cells.pdf",
  plot = p_traj_cluster_all,
  width = 7,
  height = 6
)

ggsave(
  "13_Monocle2_DDRTree_CellCycle_phase_all_cells.pdf",
  plot = p_traj_phase_all,
  width = 7,
  height = 6
)

pData(cds_all)$EMT_state_article_like <- as.character(
  pData(cds_all)$EMT_state_final
)

pData(cds_all)$EMT_state_article_like[
  pData(cds_all)$EMT_state_article_like %in% c(
    "early_hybrid",
    "late_hybrid"
  )
] <- "Hybrid"

pData(cds_all)$EMT_state_article_like[
  pData(cds_all)$EMT_state_article_like == "epithelial"
] <- "Epithelial"

pData(cds_all)$EMT_state_article_like[
  pData(cds_all)$EMT_state_article_like == "full_EMT"
] <- "Late EMT"

pData(cds_all)$EMT_state_article_like <- factor(
  pData(cds_all)$EMT_state_article_like,
  levels = c(
    "Epithelial",
    "Hybrid",
    "Late EMT"
  )
)

table(pData(cds_all)$EMT_state_article_like)

p_traj_emt_article_like_all <- plot_cell_trajectory_modern(
  cds_all,
  color_by = "EMT_state_article_like",
  title = "Monocle2 DDRTree trajectory: EMT states, article-like grouping"
)

print(p_traj_emt_article_like_all)

ggsave(
  "14_Monocle2_DDRTree_EMT_state_article_like_all_cells.pdf",
  plot = p_traj_emt_article_like_all,
  width = 7,
  height = 6
)

ggsave(
  "14_Monocle2_DDRTree_EMT_state_article_like_all_cells.png",
  plot = p_traj_emt_article_like_all,
  width = 7,
  height = 6,
  dpi = 300
)

pseudotime_meta_all <- pData(cds_all) %>%
  as.data.frame()

pseudotime_meta_all$EMT_state_final <- factor(
  pseudotime_meta_all$EMT_state_final,
  levels = c(
    "epithelial",
    "early_hybrid",
    "late_hybrid",
    "full_EMT"
  )
)

p_pseudotime_by_emt_all <- ggplot(
  pseudotime_meta_all,
  aes(x = EMT_state_final, y = Pseudotime)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  theme_classic() +
  xlab("EMT state") +
  ylab("Monocle2 pseudotime") +
  ggtitle("Pseudotime distribution across EMT states: all cells") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_pseudotime_by_emt_all)

ggsave(
  "15_Pseudotime_distribution_by_EMT_state_all_cells.pdf",
  plot = p_pseudotime_by_emt_all,
  width = 7,
  height = 5
)

pseudotime_by_emt_summary_all <- pseudotime_meta_all %>%
  group_by(EMT_state_final) %>%
  summarise(
    cell_number = n(),
    mean_pseudotime = mean(Pseudotime, na.rm = TRUE),
    median_pseudotime = median(Pseudotime, na.rm = TRUE),
    sd_pseudotime = sd(Pseudotime, na.rm = TRUE),
    .groups = "drop"
  )

pseudotime_by_emt_summary_all

write.csv(
  pseudotime_by_emt_summary_all,
  file = "16_Pseudotime_summary_by_EMT_state_all_cells.csv",
  row.names = FALSE
)

pseudotime_meta_all$EMT_state_article_like <- pData(cds_all)$EMT_state_article_like

p_pseudotime_by_article_group_all <- ggplot(
  pseudotime_meta_all,
  aes(x = EMT_state_article_like, y = Pseudotime)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  theme_classic() +
  xlab("EMT state") +
  ylab("Monocle2 pseudotime") +
  ggtitle("Pseudotime distribution: Epithelial, Hybrid, Late EMT") +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

print(p_pseudotime_by_article_group_all)

ggsave(
  "17_Pseudotime_distribution_article_like_EMT_groups_all_cells.pdf",
  plot = p_pseudotime_by_article_group_all,
  width = 7,
  height = 5
)

emt_markers <- c(
  "Epcam",
  "Krt14",
  "Krt8",
  "Vim",
  "Pdgfra"
)

emt_markers_present <- emt_markers[
  emt_markers %in% rownames(cds_all)
]

emt_markers_present

cds_emt_markers_all <- cds_all[emt_markers_present, ]

p_emt_markers_pseudotime_all <- plot_genes_in_pseudotime(
  cds_emt_markers_all,
  color_by = "EMT_state_final"
) +
  ggtitle("EMT marker expression along Monocle2 pseudotime: all cells")

print(p_emt_markers_pseudotime_all)

ggsave(
  "18_EMT_markers_along_Monocle2_pseudotime_all_cells.pdf",
  plot = p_emt_markers_pseudotime_all,
  width = 10,
  height = 8
)

monocle_pseudotime_all_df <- data.frame(
  cell_barcode = rownames(pData(cds_all)),
  Monocle2_Pseudotime_all_cells = pData(cds_all)$Pseudotime,
  Monocle2_State_all_cells = pData(cds_all)$State,
  EMT_state_final = pData(cds_all)$EMT_state_final,
  EMT_state_article_like = pData(cds_all)$EMT_state_article_like,
  tumor_final_cluster_res0.7 = pData(cds_all)$tumor_final_cluster_res0.7,
  Phase = pData(cds_all)$Phase,
  stringsAsFactors = FALSE
)

write.csv(
  monocle_pseudotime_all_df,
  file = "19_Monocle2_all_cells_pseudotime_per_cell.csv",
  row.names = FALSE
)

tumor_final$Monocle2_Pseudotime_all_cells <- NA
tumor_final$Monocle2_State_all_cells <- NA

common_cells <- intersect(
  colnames(tumor_final),
  monocle_pseudotime_all_df$cell_barcode
)

tumor_final$Monocle2_Pseudotime_all_cells[common_cells] <-
  monocle_pseudotime_all_df$Monocle2_Pseudotime_all_cells[
    match(common_cells, monocle_pseudotime_all_df$cell_barcode)
  ]

tumor_final$Monocle2_State_all_cells[common_cells] <-
  monocle_pseudotime_all_df$Monocle2_State_all_cells[
    match(common_cells, monocle_pseudotime_all_df$cell_barcode)
  ]

p_umap_pseudotime_all <- FeaturePlot(
  tumor_final,
  features = "Monocle2_Pseudotime_all_cells",
  reduction = "umap"
) +
  ggtitle("All-cell Monocle2 pseudotime projected on Seurat UMAP")

print(p_umap_pseudotime_all)

ggsave(
  "20_Monocle2_all_cells_pseudotime_on_Seurat_UMAP.pdf",
  plot = p_umap_pseudotime_all,
  width = 7,
  height = 6
)

saveRDS(
  cds_all,
  file = "21_GSE234265_control_2_YFP_all_cells_Monocle2_DDRTree_final_cds.rds"
)

saveRDS(
  tumor_final,
  file = "22_GSE234265_control_2_YFP_tumor_final_with_all_cell_Monocle2_pseudotime.rds"
)

message("All-cell Monocle2 pseudotime analysis finished successfully.")
