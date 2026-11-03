# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(1234)

required_pkgs <- c(
  "Seurat", "dplyr", "ggplot2", "patchwork",
  "Matrix", "tidyr", "stringr", "scales"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(tidyr)
  library(stringr)
  library(scales)
})

BASE_DIR <- "path/to/GSE234265"

REANALYSIS_DIR <- file.path(
  BASE_DIR,
  "control_2_reanalysis_V2"
)

GLOBAL_RDS <- file.path(
  REANALYSIS_DIR,
  "GSE234265_control_2_global_with_YFP_info.rds"
)

TUMOR_RDS <- file.path(
  REANALYSIS_DIR,
  "GSE234265_control_2_YFP_tumor_final_EMT_annotated.rds"
)

GMT_FILE <- file.path(
  BASE_DIR,
  "gene set",
  "STEP_all_human_gene_sets_merged.gmt"
)

OUT_DIR <- file.path(
  REANALYSIS_DIR,
  "Fig3_scRNA_CAF_STEP_validation_v1_1"
)

DIR_QC <- file.path(OUT_DIR, "00_QC")
DIR_ANNOT <- file.path(OUT_DIR, "01_CAF_annotation")
DIR_SCORE <- file.path(OUT_DIR, "02_scores")
DIR_UMAP <- file.path(OUT_DIR, "03_Fig3_UMAP")
DIR_VLN <- file.path(OUT_DIR, "04_Fig3_violin")
DIR_TABLE <- file.path(OUT_DIR, "05_tables")
DIR_RDS <- file.path(OUT_DIR, "06_RDS")

for (d in c(
  OUT_DIR, DIR_QC, DIR_ANNOT, DIR_SCORE,
  DIR_UMAP, DIR_VLN, DIR_TABLE, DIR_RDS
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

GLOBAL_CLUSTER_COL <- "global_cluster_res0.1"

EMT_COL <- "EMT_state_final"

MIN_COL6A3_PCT <- 0.10

MANUAL_CAF_CLUSTERS <- NULL

MAIN_CAF_SCORE <- "CAF3_Zmean"
MAIN_STEP_SCORE <- "STEP_UP_Zmean"

PNG_DPI <- 300

FIG3_COLORS <- c(
  "Epithelial"    = "#2C7FB8",
  "Hybrid"        = "#F39C34",
  "Mesenchymal"   = "#D73027",
  "Fibroblast/CAF"= "#2A9D55",
  "Other"         = "grey82"
)

FIG3_LEVELS <- c(
  "Epithelial",
  "Hybrid",
  "Mesenchymal",
  "Fibroblast/CAF",
  "Other"
)

for (f in c(GLOBAL_RDS, TUMOR_RDS, GMT_FILE)) {
  if (!file.exists(f)) {
    stop("Required input file not found: ", f)
  }
}

message("Reading global object...")
global <- readRDS(GLOBAL_RDS)

message("Reading final EMT-annotated tumour object...")
tumor <- readRDS(TUMOR_RDS)

DefaultAssay(global) <- "RNA"
DefaultAssay(tumor) <- "RNA"

if (!GLOBAL_CLUSTER_COL %in% colnames(global@meta.data)) {
  stop(
    "Global cluster column not found: ", GLOBAL_CLUSTER_COL,
    "\nAvailable metadata columns:\n",
    paste(colnames(global@meta.data), collapse = ", ")
  )
}

if (!"YFP_positive" %in% colnames(global@meta.data)) {
  stop("YFP_positive not found in global metadata.")
}

if (!EMT_COL %in% colnames(tumor@meta.data)) {
  stop("EMT_state_final not found in tumour metadata.")
}

if (!"umap" %in% names(global@reductions)) {
  stop("UMAP reduction not found in global object.")
}

get_norm_data <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) {
      GetAssayData(obj, assay = assay, slot = "data")
    }
  )
}

get_counts_data <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) {
      GetAssayData(obj, assay = assay, slot = "counts")
    }
  )
}

norm_test <- tryCatch(
  get_norm_data(global),
  error = function(e) NULL
)

if (is.null(norm_test) || nrow(norm_test) == 0) {
  message("Normalized RNA data not found. Running NormalizeData(global)...")
  global <- NormalizeData(
    global,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
}

common_cells <- intersect(
  colnames(tumor),
  colnames(global)
)

barcode_qc <- data.frame(
  global_cells = ncol(global),
  final_tumor_cells = ncol(tumor),
  matched_tumor_cells_in_global = length(common_cells),
  matched_fraction_of_final_tumor =
    length(common_cells) / ncol(tumor),
  stringsAsFactors = FALSE
)

write.csv(
  barcode_qc,
  file.path(DIR_QC, "01_barcode_matching_QC.csv"),
  row.names = FALSE
)

print(barcode_qc)

if (length(common_cells) == 0) {
  stop("No matching barcodes between global and tumour objects.")
}

global$EMT_state_from_final_tumor <- NA_character_

global@meta.data[
  common_cells,
  "EMT_state_from_final_tumor"
] <- as.character(
  tumor@meta.data[
    common_cells,
    EMT_COL
  ]
)

article_lineage_modules <- list(
  CAF = c(
    "Col1a1", "Col1a2", "Col3a1",
    "Col5a1", "Col5a2",
    "Col6a1", "Col6a2", "Col6a3",
    "Dcn", "Lum", "Pdgfra", "Fbln1"
  ),
  Immune = c(
    "Ptprc", "Lyz2", "Tyrobp",
    "Cd3d", "Cd3e", "Nkg7",
    "Ms4a1", "Cd79a", "Itgam"
  ),
  Endothelial = c(
    "Pecam1", "Cdh5", "Kdr",
    "Emcn", "Esam", "Klf2", "Robo4"
  ),
  Pericyte = c(
    "Rgs5", "Cspg4", "Notch3",
    "Pdgfrb", "Des", "Mcam", "Rbp1"
  ),
  RBC_like = c(
    "Hbb-bs", "Hbb-bt",
    "Hba-a1", "Hba-a2", "Alas2"
  )
)

article_lineage_modules_present <- lapply(
  article_lineage_modules,
  function(g) {
    g[g %in% rownames(global)]
  }
)

article_module_gene_qc <- bind_rows(
  lapply(
    names(article_lineage_modules),
    function(nm) {
      data.frame(
        module = nm,
        gene = article_lineage_modules[[nm]],
        present = article_lineage_modules[[nm]] %in% rownames(global),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  article_module_gene_qc,
  file.path(DIR_QC, "02_article_lineage_marker_gene_presence.csv"),
  row.names = FALSE
)

article_lineage_modules_present <-
  article_lineage_modules_present[
    lengths(article_lineage_modules_present) > 0
  ]

if (!"CAF" %in% names(article_lineage_modules_present)) {
  stop("No CAF canonical markers were found in the global object.")
}

global <- AddModuleScore(
  object = global,
  features = article_lineage_modules_present,
  name = "ArticleLineageModule"
)

module_cols <- paste0(
  "ArticleLineageModule",
  seq_along(article_lineage_modules_present)
)

module_map <- data.frame(
  module = names(article_lineage_modules_present),
  score_column = module_cols,
  n_genes_used = lengths(article_lineage_modules_present),
  genes_used = vapply(
    article_lineage_modules_present,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  module_map,
  file.path(DIR_QC, "03_article_lineage_module_score_mapping.csv"),
  row.names = FALSE
)

for (i in seq_len(nrow(module_map))) {
  global@meta.data[
    ,
    paste0(module_map$module[i], "_article_module")
  ] <- global@meta.data[
    ,
    module_map$score_column[i]
  ]
}

canonical_markers <- c(
  "YFP",
  "Col6a3", "Pdgfra", "Dcn", "Lum", "Acta2",
  "Ptprc", "Cd3d",
  "Pecam1",
  "Rgs5", "Pdgfrb"
)

canonical_markers_present <- canonical_markers[
  canonical_markers %in% rownames(global)
]

if (length(canonical_markers_present) > 0) {
  p_marker_features <- FeaturePlot(
    global,
    features = canonical_markers_present,
    reduction = "umap",
    ncol = 3,
    order = TRUE
  )

  print(p_marker_features)

  ggsave(
    file.path(
      DIR_ANNOT,
      "01_global_UMAP_canonical_celltype_markers.pdf"
    ),
    p_marker_features,
    width = 12,
    height = 4 * ceiling(length(canonical_markers_present) / 3),
    limitsize = FALSE
  )

  ggsave(
    file.path(
      DIR_ANNOT,
      "01_global_UMAP_canonical_celltype_markers.png"
    ),
    p_marker_features,
    width = 12,
    height = 4 * ceiling(length(canonical_markers_present) / 3),
    dpi = PNG_DPI,
    limitsize = FALSE
  )
}

dot_genes <- unique(c(
  "Col1a1", "Col1a2", "Col3a1", "Col6a3", "Dcn", "Lum", "Pdgfra",
  "Ptprc", "Lyz2", "Cd3d", "Nkg7",
  "Pecam1", "Cdh5", "Kdr",
  "Rgs5", "Cspg4", "Pdgfrb",
  "YFP"
))

dot_genes <- dot_genes[
  dot_genes %in% rownames(global)
]

if (length(dot_genes) > 0) {
  p_dot <- DotPlot(
    global,
    features = dot_genes,
    group.by = GLOBAL_CLUSTER_COL
  ) +
    RotatedAxis() +
    ggtitle("Global clusters: article-aligned lineage markers")

  print(p_dot)

  ggsave(
    file.path(
      DIR_ANNOT,
      "02_global_cluster_marker_DotPlot.pdf"
    ),
    p_dot,
    width = 12,
    height = 6
  )

  ggsave(
    file.path(
      DIR_ANNOT,
      "02_global_cluster_marker_DotPlot.png"
    ),
    p_dot,
    width = 12,
    height = 6,
    dpi = PNG_DPI
  )
}

meta <- global@meta.data
meta$cell_barcode <- rownames(meta)

yfp_bool <- as.character(meta$YFP_positive)
yfp_bool <- toupper(yfp_bool) %in% c("TRUE", "T", "1")

meta$YFP_positive_bool <- yfp_bool

global$YFP_positive_bool <- yfp_bool

meta_neg <- meta[
  !meta$YFP_positive_bool,
  ,
  drop = FALSE
]

lineage_score_cols <- paste0(
  names(article_lineage_modules_present),
  "_article_module"
)

cluster_module_summary <- meta_neg %>%
  group_by(.data[[GLOBAL_CLUSTER_COL]]) %>%
  summarise(
    n_YFPneg = n(),
    across(
      all_of(lineage_score_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cluster_yfp_summary <- meta %>%
  group_by(.data[[GLOBAL_CLUSTER_COL]]) %>%
  summarise(
    n_total = n(),
    n_YFPpos = sum(YFP_positive_bool, na.rm = TRUE),
    YFPpos_fraction = mean(YFP_positive_bool, na.rm = TRUE),
    .groups = "drop"
  )

norm_mat <- get_norm_data(global)

marker_pct_by_cluster <- function(gene) {
  if (!gene %in% rownames(norm_mat)) {
    return(NULL)
  }

  neg_cells <- rownames(meta_neg)
  values <- as.numeric(
    norm_mat[
      gene,
      neg_cells,
      drop = TRUE
    ]
  )

  tmp <- data.frame(
    cluster = as.character(
      meta_neg[
        neg_cells,
        GLOBAL_CLUSTER_COL
      ]
    ),
    expr = values,
    stringsAsFactors = FALSE
  )

  out <- tmp %>%
    group_by(cluster) %>%
    summarise(
      pct_expr = mean(expr > 0, na.rm = TRUE),
      mean_expr = mean(expr, na.rm = TRUE),
      .groups = "drop"
    )

  names(out)[
    names(out) == "cluster"
  ] <- GLOBAL_CLUSTER_COL

  names(out)[
    names(out) == "pct_expr"
  ] <- paste0("pct_", gene)

  names(out)[
    names(out) == "mean_expr"
  ] <- paste0("mean_", gene)

  out
}

canonical_cluster_tables <- lapply(
  c(
    "Col6a3", "Pdgfra", "Dcn", "Lum",
    "Ptprc", "Pecam1", "Rgs5"
  ),
  marker_pct_by_cluster
)

canonical_cluster_tables <-
  canonical_cluster_tables[
    !vapply(
      canonical_cluster_tables,
      is.null,
      logical(1)
    )
  ]

cluster_qc <- cluster_module_summary %>%
  left_join(
    cluster_yfp_summary,
    by = GLOBAL_CLUSTER_COL
  )

for (tb in canonical_cluster_tables) {
  cluster_qc <- cluster_qc %>%
    left_join(
      tb,
      by = GLOBAL_CLUSTER_COL
    )
}

score_matrix <- as.matrix(
  cluster_qc[
    ,
    lineage_score_cols,
    drop = FALSE
  ]
)

colnames(score_matrix) <- names(
  article_lineage_modules_present
)

cluster_qc$dominant_YFPneg_lineage <- apply(
  score_matrix,
  1,
  function(x) {
    if (all(is.na(x))) {
      return(NA_character_)
    }
    names(which.max(x))[1]
  }
)

if ("pct_Col6a3" %in% colnames(cluster_qc)) {
  cluster_qc$CAF_candidate_auto <-
    cluster_qc$dominant_YFPneg_lineage == "CAF" &
    !is.na(cluster_qc$pct_Col6a3) &
    cluster_qc$pct_Col6a3 >= MIN_COL6A3_PCT
} else {
  warning(
    "Col6a3 not found. Auto CAF call will use CAF-module dominance only."
  )
  cluster_qc$CAF_candidate_auto <-
    cluster_qc$dominant_YFPneg_lineage == "CAF"
}

write.csv(
  cluster_qc,
  file.path(
    DIR_ANNOT,
    "03_YFPnegative_global_cluster_article_aligned_annotation_QC.csv"
  ),
  row.names = FALSE
)

if (is.null(MANUAL_CAF_CLUSTERS)) {

  CAF_CLUSTERS <- as.character(
    cluster_qc[[GLOBAL_CLUSTER_COL]][
      cluster_qc$CAF_candidate_auto
    ]
  )
  CAF_SOURCE <- "automatic_article_aligned"
} else {
  CAF_CLUSTERS <- as.character(
    MANUAL_CAF_CLUSTERS
  )
  CAF_SOURCE <- "manual_override"
}

CAF_CLUSTERS <- unique(CAF_CLUSTERS)

if (length(CAF_CLUSTERS) > 0 &&
    any(grepl("^c\\(", CAF_CLUSTERS))) {
  stop(
    "CAF_CLUSTERS is malformed (collapsed vector). ",
    "Expected values such as c(\"2\", \"3\", \"5\")."
  )
}

write.csv(
  data.frame(
    CAF_cluster = CAF_CLUSTERS,
    source = CAF_SOURCE,
    min_Col6a3_fraction_for_auto =
      MIN_COL6A3_PCT,
    stringsAsFactors = FALSE
  ),
  file.path(
    DIR_ANNOT,
    "04_final_CAF_clusters_used.csv"
  ),
  row.names = FALSE
)

message(
  "CAF clusters used: ",
  ifelse(
    length(CAF_CLUSTERS) == 0,
    "NONE",
    paste(CAF_CLUSTERS, collapse = ", ")
  )
)

if (length(CAF_CLUSTERS) == 0) {
  warning(
    "No CAF cluster was automatically selected. ",
    "Inspect 03_YFPnegative_global_cluster_article_aligned_annotation_QC.csv ",
    "and marker plots. If needed, set MANUAL_CAF_CLUSTERS at the top."
  )
}

global$article_broad_celltype <- "Other"

global$article_broad_celltype[
  global@meta.data$YFP_positive_bool
] <- "Tumour_YFP+"

cluster_to_lineage <- setNames(
  cluster_qc$dominant_YFPneg_lineage,
  as.character(
    cluster_qc[[GLOBAL_CLUSTER_COL]]
  )
)

neg_cells <- rownames(
  global@meta.data
)[
  !global@meta.data$YFP_positive_bool
]

neg_cluster <- as.character(
  global@meta.data[
    neg_cells,
    GLOBAL_CLUSTER_COL
  ]
)

global@meta.data[
  neg_cells,
  "article_broad_celltype"
] <- unname(
  cluster_to_lineage[
    neg_cluster
  ]
)

is_caf <- (
  !global@meta.data$YFP_positive_bool
) & (
  as.character(
    global@meta.data[[GLOBAL_CLUSTER_COL]]
  ) %in% CAF_CLUSTERS
)

global$CAF_article_aligned <- is_caf

global$Fig3_cell_group <- "Other"

emt_raw <- as.character(
  global$EMT_state_from_final_tumor
)

global$Fig3_cell_group[
  emt_raw == "epithelial"
] <- "Epithelial"

global$Fig3_cell_group[
  emt_raw %in% c(
    "early_hybrid",
    "late_hybrid"
  )
] <- "Hybrid"

global$Fig3_cell_group[
  emt_raw %in% c(
    "full_EMT",
    "Full_EMT",
    "mesenchymal",
    "Mesenchymal"
  )
] <- "Mesenchymal"

global$Fig3_cell_group[
  global$CAF_article_aligned
] <- "Fibroblast/CAF"

global$Fig3_cell_group <- factor(
  global$Fig3_cell_group,
  levels = FIG3_LEVELS
)

group_counts <- as.data.frame(
  table(
    global$Fig3_cell_group,
    useNA = "ifany"
  )
)

colnames(group_counts) <- c(
  "Fig3_cell_group",
  "Cell_number"
)

write.csv(
  group_counts,
  file.path(
    DIR_TABLE,
    "01_Fig3_cell_group_counts.csv"
  ),
  row.names = FALSE
)

gmt_lines <- readLines(
  GMT_FILE,
  warn = FALSE
)

caf_line <- gmt_lines[
  grepl(
    "^CAF_3gene_score\\t",
    gmt_lines
  )
]

if (length(caf_line) != 1) {
  stop("CAF_3gene_score was not uniquely found in GMT.")
}

caf_fields <- strsplit(
  caf_line,
  "\t",
  fixed = TRUE
)[[1]]

CAF3_HUMAN <- unique(
  trimws(
    caf_fields[
      seq.int(
        3,
        length(caf_fields)
      )
    ]
  )
)

CAF3_HUMAN <- CAF3_HUMAN[
  nzchar(CAF3_HUMAN)
]

raw_line <- gmt_lines[
  grepl(
    "^RAW_STEP_UP\\b",
    gmt_lines
  )
]

if (length(raw_line) != 1) {
  stop("RAW_STEP_UP was not uniquely found in GMT.")
}

raw_fields <- strsplit(
  raw_line,
  "\t",
  fixed = TRUE
)[[1]]

first_tokens <- strsplit(
  trimws(raw_fields[1]),
  "\\s+"
)[[1]]

RAW_STEP_UP_HUMAN <- unique(c(
  first_tokens[
    first_tokens != "RAW_STEP_UP"
  ],
  raw_fields[
    seq.int(
      2,
      length(raw_fields)
    )
  ]
))

RAW_STEP_UP_HUMAN <- trimws(
  RAW_STEP_UP_HUMAN
)

RAW_STEP_UP_HUMAN <-
  RAW_STEP_UP_HUMAN[
    nzchar(RAW_STEP_UP_HUMAN)
  ]

gene_set_input_qc <- data.frame(
  gene_set = c(
    "CAF_3gene_score",
    "RAW_STEP_UP"
  ),
  human_gene_count = c(
    length(CAF3_HUMAN),
    length(RAW_STEP_UP_HUMAN)
  ),
  genes = c(
    paste(
      CAF3_HUMAN,
      collapse = ";"
    ),
    paste(
      RAW_STEP_UP_HUMAN,
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gene_set_input_qc,
  file.path(
    DIR_QC,
    "04_GMT_gene_sets_parsed.csv"
  ),
  row.names = FALSE
)

map_human_to_mouse_features <- function(
    human_genes,
    mouse_features
) {
  feature_upper <- toupper(
    mouse_features
  )

  feature_map <- setNames(
    mouse_features[
      !duplicated(feature_upper)
    ],
    feature_upper[
      !duplicated(feature_upper)
    ]
  )

  matched <- unname(
    feature_map[
      toupper(human_genes)
    ]
  )

  data.frame(
    human_gene = human_genes,
    mouse_feature = matched,
    matched = !is.na(matched),
    stringsAsFactors = FALSE
  )
}

caf3_map <- map_human_to_mouse_features(
  CAF3_HUMAN,
  rownames(global)
)

step_map <- map_human_to_mouse_features(
  RAW_STEP_UP_HUMAN,
  rownames(global)
)

write.csv(
  caf3_map,
  file.path(
    DIR_QC,
    "05_CAF3_human_to_mouse_feature_matching.csv"
  ),
  row.names = FALSE
)

write.csv(
  step_map,
  file.path(
    DIR_QC,
    "06_RAW_STEP_UP_human_to_mouse_feature_matching.csv"
  ),
  row.names = FALSE
)

CAF3_MOUSE <- unique(
  caf3_map$mouse_feature[
    caf3_map$matched
  ]
)

STEP_UP_MOUSE <- unique(
  step_map$mouse_feature[
    step_map$matched
  ]
)

score_gene_qc <- data.frame(
  gene_set = c(
    "CAF3",
    "RAW_STEP_UP"
  ),
  human_input_n = c(
    length(CAF3_HUMAN),
    length(RAW_STEP_UP_HUMAN)
  ),
  mouse_matched_n = c(
    length(CAF3_MOUSE),
    length(STEP_UP_MOUSE)
  ),
  mouse_genes_used = c(
    paste(CAF3_MOUSE, collapse = ";"),
    paste(STEP_UP_MOUSE, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  score_gene_qc,
  file.path(
    DIR_QC,
    "07_score_gene_matching_summary.csv"
  ),
  row.names = FALSE
)

print(score_gene_qc)

if (length(CAF3_MOUSE) < 2) {
  stop("Fewer than 2 CAF3 genes matched the mouse object.")
}

if (length(STEP_UP_MOUSE) < 10) {
  stop(
    "Too few RAW_STEP_UP genes matched the mouse object. ",
    "Inspect 06_RAW_STEP_UP_human_to_mouse_feature_matching.csv"
  )
}

genes_to_scale <- unique(c(
  CAF3_MOUSE,
  STEP_UP_MOUSE
))

global <- ScaleData(
  global,
  assay = "RNA",
  features = genes_to_scale,
  verbose = FALSE
)

get_scale_data <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(
      obj,
      assay = assay,
      layer = "scale.data"
    ),
    error = function(e) {
      GetAssayData(
        obj,
        assay = assay,
        slot = "scale.data"
      )
    }
  )
}

scale_mat <- get_scale_data(global)

calc_zmean <- function(
    scale_matrix,
    genes
) {
  genes_use <- genes[
    genes %in% rownames(scale_matrix)
  ]

  if (length(genes_use) == 0) {
    return(
      rep(
        NA_real_,
        ncol(scale_matrix)
      )
    )
  }

  Matrix::colMeans(
    scale_matrix[
      genes_use,
      ,
      drop = FALSE
    ]
  )
}

global$CAF3_Zmean <- calc_zmean(
  scale_mat,
  CAF3_MOUSE
)

global$STEP_UP_Zmean <- calc_zmean(
  scale_mat,
  STEP_UP_MOUSE
)

global <- AddModuleScore(
  global,
  features = list(CAF3_MOUSE),
  name = "CAF3_AMS"
)

global <- AddModuleScore(
  global,
  features = list(STEP_UP_MOUSE),
  name = "STEP_UP_AMS"
)

if (!"CAF3_AMS1" %in% colnames(global@meta.data)) {
  stop("CAF3_AMS1 was not created.")
}

if (!"STEP_UP_AMS1" %in% colnames(global@meta.data)) {
  stop("STEP_UP_AMS1 was not created.")
}

p_fig3_umap <- DimPlot(
  global,
  reduction = "umap",
  group.by = "Fig3_cell_group",
  cols = FIG3_COLORS,
  order = c(
    "Other",
    "Fibroblast/CAF",
    "Epithelial",
    "Hybrid",
    "Mesenchymal"
  ),
  pt.size = 0.55
) +
  ggtitle("GSE234265 control_2: Epithelial / Hybrid / Mesenchymal / Fibroblast-CAF") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

print(p_fig3_umap)

ggsave(
  file.path(
    DIR_UMAP,
    "01_Fig3_global_UMAP_E_H_M_Fibroblast_CAF.pdf"
  ),
  p_fig3_umap,
  width = 7.5,
  height = 6
)

ggsave(
  file.path(
    DIR_UMAP,
    "01_Fig3_global_UMAP_E_H_M_Fibroblast_CAF.png"
  ),
  p_fig3_umap,
  width = 7.5,
  height = 6,
  dpi = PNG_DPI
)

p_broad <- DimPlot(
  global,
  reduction = "umap",
  group.by = "article_broad_celltype",
  pt.size = 0.45
) +
  ggtitle("Article-aligned broad cell-type annotation") +
  theme_classic(base_size = 11)

print(p_broad)

ggsave(
  file.path(
    DIR_ANNOT,
    "05_article_aligned_broad_celltype_UMAP.pdf"
  ),
  p_broad,
  width = 7.5,
  height = 6
)

ggsave(
  file.path(
    DIR_ANNOT,
    "05_article_aligned_broad_celltype_UMAP.png"
  ),
  p_broad,
  width = 7.5,
  height = 6,
  dpi = PNG_DPI
)

p_caf_z <- FeaturePlot(
  global,
  features = "CAF3_Zmean",
  reduction = "umap",
  order = TRUE
) +
  ggtitle("CAF3 Z-mean")

p_step_z <- FeaturePlot(
  global,
  features = "STEP_UP_Zmean",
  reduction = "umap",
  order = TRUE
) +
  ggtitle("RAW STEP_UP Z-mean")

p_score_z_pair <- p_caf_z + p_step_z +
  plot_annotation(
    title = "GSE234265 single-cell score distribution"
  )

print(p_score_z_pair)

ggsave(
  file.path(
    DIR_SCORE,
    "01_CAF3_and_RAW_STEP_UP_Zmean_UMAP.pdf"
  ),
  p_score_z_pair,
  width = 12,
  height = 5.5
)

ggsave(
  file.path(
    DIR_SCORE,
    "01_CAF3_and_RAW_STEP_UP_Zmean_UMAP.png"
  ),
  p_score_z_pair,
  width = 12,
  height = 5.5,
  dpi = PNG_DPI
)

p_caf_ams <- FeaturePlot(
  global,
  features = "CAF3_AMS1",
  reduction = "umap",
  order = TRUE
) +
  ggtitle("CAF3 AddModuleScore")

p_step_ams <- FeaturePlot(
  global,
  features = "STEP_UP_AMS1",
  reduction = "umap",
  order = TRUE
) +
  ggtitle("RAW STEP_UP AddModuleScore")

p_score_ams_pair <- p_caf_ams + p_step_ams +
  plot_annotation(
    title = "Sensitivity analysis: Seurat AddModuleScore"
  )

print(p_score_ams_pair)

ggsave(
  file.path(
    DIR_SCORE,
    "02_CAF3_and_RAW_STEP_UP_AddModuleScore_UMAP.pdf"
  ),
  p_score_ams_pair,
  width = 12,
  height = 5.5
)

ggsave(
  file.path(
    DIR_SCORE,
    "02_CAF3_and_RAW_STEP_UP_AddModuleScore_UMAP.png"
  ),
  p_score_ams_pair,
  width = 12,
  height = 5.5,
  dpi = PNG_DPI
)

main_groups <- FIG3_LEVELS[
  FIG3_LEVELS != "Other"
]

plot_df <- global@meta.data %>%
  as.data.frame() %>%
  mutate(
    cell_barcode = rownames(global@meta.data)
  ) %>%
  filter(
    as.character(Fig3_cell_group) %in% main_groups
  ) %>%
  mutate(
    Fig3_cell_group = factor(
      as.character(Fig3_cell_group),
      levels = main_groups
    )
  )

write.csv(
  plot_df[
    ,
    c(
      "cell_barcode",
      GLOBAL_CLUSTER_COL,
      "YFP_positive",
      "EMT_state_from_final_tumor",
      "CAF_article_aligned",
      "Fig3_cell_group",
      "CAF3_Zmean",
      "STEP_UP_Zmean",
      "CAF3_AMS1",
      "STEP_UP_AMS1"
    ),
    drop = FALSE
  ],
  file.path(
    DIR_TABLE,
    "02_Fig3_per_cell_scores_and_groups.csv"
  ),
  row.names = FALSE
)

run_group_stats <- function(
    df,
    score_col,
    group_col = "Fig3_cell_group"
) {
  d <- df[
    is.finite(df[[score_col]]) &
      !is.na(df[[group_col]]),
    ,
    drop = FALSE
  ]

  kw <- kruskal.test(
    d[[score_col]],
    d[[group_col]]
  )

  pw <- pairwise.wilcox.test(
    x = d[[score_col]],
    g = d[[group_col]],
    p.adjust.method = "BH",
    exact = FALSE
  )

  pw_df <- as.data.frame(
    as.table(
      pw$p.value
    )
  )

  colnames(pw_df) <- c(
    "group1",
    "group2",
    "BH_adjusted_p"
  )

  pw_df <- pw_df[
    !is.na(pw_df$BH_adjusted_p),
    ,
    drop = FALSE
  ]

  list(
    kruskal = data.frame(
      score = score_col,
      Kruskal_Wallis_chisq =
        unname(kw$statistic),
      df = unname(kw$parameter),
      p_value = kw$p.value,
      stringsAsFactors = FALSE
    ),
    pairwise = pw_df
  )
}

caf_stats <- run_group_stats(
  plot_df,
  MAIN_CAF_SCORE
)

step_stats <- run_group_stats(
  plot_df,
  MAIN_STEP_SCORE
)

write.csv(
  caf_stats$kruskal,
  file.path(
    DIR_TABLE,
    "03_CAF3_Zmean_Kruskal_Wallis.csv"
  ),
  row.names = FALSE
)

write.csv(
  caf_stats$pairwise,
  file.path(
    DIR_TABLE,
    "04_CAF3_Zmean_pairwise_Wilcoxon_BH.csv"
  ),
  row.names = FALSE
)

write.csv(
  step_stats$kruskal,
  file.path(
    DIR_TABLE,
    "05_STEP_UP_Zmean_Kruskal_Wallis.csv"
  ),
  row.names = FALSE
)

write.csv(
  step_stats$pairwise,
  file.path(
    DIR_TABLE,
    "06_STEP_UP_Zmean_pairwise_Wilcoxon_BH.csv"
  ),
  row.names = FALSE
)

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) {
    return(
      formatC(
        p,
        format = "e",
        digits = 2
      )
    )
  }
  format(
    signif(p, 3),
    trim = TRUE
  )
}

make_violin <- function(
    df,
    score_col,
    ylab_text,
    title_text,
    kw_p
) {
  ggplot(
    df,
    aes(
      x = Fig3_cell_group,
      y = .data[[score_col]],
      fill = Fig3_cell_group
    )
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      alpha = 0.80,
      color = "grey30",
      linewidth = 0.25
    ) +
    geom_boxplot(
      width = 0.14,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_manual(
      values = FIG3_COLORS[
        main_groups
      ],
      drop = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Kruskal-Wallis p = ",
        format_p(kw_p)
      ),
      x = NULL,
      y = ylab_text
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      )
    )
}

p_vln_caf <- make_violin(
  plot_df,
  MAIN_CAF_SCORE,
  "CAF3 Z-mean score",
  "CAF3 across Epithelial / Hybrid / Mesenchymal / Fibroblast-CAF",
  caf_stats$kruskal$p_value
)

p_vln_step <- make_violin(
  plot_df,
  MAIN_STEP_SCORE,
  "RAW STEP_UP Z-mean score",
  "RAW STEP_UP across Epithelial / Hybrid / Mesenchymal / Fibroblast-CAF",
  step_stats$kruskal$p_value
)

p_vln_pair <- p_vln_caf + p_vln_step +
  plot_annotation(
    title = "GSE234265 control_2"
  )

print(p_vln_pair)

ggsave(
  file.path(
    DIR_VLN,
    "01_Fig3_CAF3_STEP_UP_violin_Zmean.pdf"
  ),
  p_vln_pair,
  width = 11,
  height = 5
)

ggsave(
  file.path(
    DIR_VLN,
    "01_Fig3_CAF3_STEP_UP_violin_Zmean.png"
  ),
  p_vln_pair,
  width = 11,
  height = 5,
  dpi = PNG_DPI
)

caf_ams_stats <- run_group_stats(
  plot_df,
  "CAF3_AMS1"
)

step_ams_stats <- run_group_stats(
  plot_df,
  "STEP_UP_AMS1"
)

p_vln_caf_ams <- make_violin(
  plot_df,
  "CAF3_AMS1",
  "CAF3 AddModuleScore",
  "CAF3 AddModuleScore",
  caf_ams_stats$kruskal$p_value
)

p_vln_step_ams <- make_violin(
  plot_df,
  "STEP_UP_AMS1",
  "RAW STEP_UP AddModuleScore",
  "RAW STEP_UP AddModuleScore",
  step_ams_stats$kruskal$p_value
)

p_vln_ams_pair <-
  p_vln_caf_ams +
  p_vln_step_ams +
  plot_annotation(
    title =
      "Sensitivity analysis: AddModuleScore"
  )

print(p_vln_ams_pair)

ggsave(
  file.path(
    DIR_VLN,
    "02_CAF3_STEP_UP_violin_AddModuleScore.pdf"
  ),
  p_vln_ams_pair,
  width = 11,
  height = 5
)

ggsave(
  file.path(
    DIR_VLN,
    "02_CAF3_STEP_UP_violin_AddModuleScore.png"
  ),
  p_vln_ams_pair,
  width = 11,
  height = 5,
  dpi = PNG_DPI
)

group_score_summary <- plot_df %>%
  group_by(Fig3_cell_group) %>%
  summarise(
    n_cells = n(),

    CAF3_Zmean_mean =
      mean(
        CAF3_Zmean,
        na.rm = TRUE
      ),
    CAF3_Zmean_median =
      median(
        CAF3_Zmean,
        na.rm = TRUE
      ),
    CAF3_Zmean_sd =
      sd(
        CAF3_Zmean,
        na.rm = TRUE
      ),

    STEP_UP_Zmean_mean =
      mean(
        STEP_UP_Zmean,
        na.rm = TRUE
      ),
    STEP_UP_Zmean_median =
      median(
        STEP_UP_Zmean,
        na.rm = TRUE
      ),
    STEP_UP_Zmean_sd =
      sd(
        STEP_UP_Zmean,
        na.rm = TRUE
      ),

    CAF3_AMS_mean =
      mean(
        CAF3_AMS1,
        na.rm = TRUE
      ),
    STEP_UP_AMS_mean =
      mean(
        STEP_UP_AMS1,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write.csv(
  group_score_summary,
  file.path(
    DIR_TABLE,
    "07_group_score_summary.csv"
  ),
  row.names = FALSE
)

cor_all <- suppressWarnings(
  cor.test(
    global$CAF3_Zmean,
    global$STEP_UP_Zmean,
    method = "spearman",
    exact = FALSE
  )
)

cor_tumour <- suppressWarnings(
  cor.test(
    plot_df$CAF3_Zmean[
      plot_df$Fig3_cell_group != "CAF"
    ],
    plot_df$STEP_UP_Zmean[
      plot_df$Fig3_cell_group != "CAF"
    ],
    method = "spearman",
    exact = FALSE
  )
)

cor_summary <- data.frame(
  population = c(
    "all_global_cells",
    "tumour_EMT_cells_only"
  ),
  spearman_rho = c(
    unname(cor_all$estimate),
    unname(cor_tumour$estimate)
  ),
  p_value = c(
    cor_all$p.value,
    cor_tumour$p.value
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cor_summary,
  file.path(
    DIR_TABLE,
    "08_CAF3_vs_STEP_UP_singlecell_Spearman.csv"
  ),
  row.names = FALSE
)

umap_df <- as.data.frame(
  Embeddings(
    global,
    reduction = "umap"
  )
)

umap_df$cell_barcode <- rownames(
  umap_df
)

meta_export <- global@meta.data
meta_export$cell_barcode <- rownames(
  meta_export
)

full_export <- meta_export %>%
  left_join(
    umap_df,
    by = "cell_barcode"
  )

write.csv(
  full_export,
  file.path(
    DIR_TABLE,
    "09_ALL_global_cells_Fig3_annotation_scores_UMAP.csv"
  ),
  row.names = FALSE
)

saveRDS(
  global,
  file.path(
    DIR_RDS,
    "GSE234265_control_2_global_Fig3_CAF_STEP_annotated_scored_v1_0.rds"
  )
)

manifest <- data.frame(
  item = c(
    "version",
    "global_RDS",
    "tumour_RDS",
    "GMT",
    "output_dir",
    "CAF_identity",
    "CAF_score",
    "STEP_score",
    "primary_singlecell_scoring",
    "sensitivity_scoring",
    "CAF_auto_rule",
    "manual_CAF_override"
  ),
  value = c(
    "v1.1",
    GLOBAL_RDS,
    TUMOR_RDS,
    GMT_FILE,
    OUT_DIR,
    paste0(
      "YFP-negative + CAF-module-dominant global cluster + canonical Col6a3 support; ",
      "immune/endothelial/pericyte/RBC-like modules used for exclusion"
    ),
    paste(
      CAF3_HUMAN,
      collapse = ";"
    ),
    paste0(
      "RAW_STEP_UP matched genes = ",
      length(STEP_UP_MOUSE),
      "/",
      length(RAW_STEP_UP_HUMAN)
    ),
    "Z-mean",
    "Seurat AddModuleScore",
    paste0(
      "dominant CAF module and pct_Col6a3 >= ",
      MIN_COL6A3_PCT
    ),
    ifelse(
      is.null(MANUAL_CAF_CLUSTERS),
      "NULL / automatic",
      paste(
        MANUAL_CAF_CLUSTERS,
        collapse = ";"
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(
    DIR_QC,
    "08_run_manifest_v1_0.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    "GSE234265 Fig.3 CAF / STEP single-cell validation v1.1",
    "",
    paste0(
      "Finished: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "",
    "CAF identity and CAF3 score are deliberately separated:",
    "  CAF identity = article-aligned YFP-negative cluster annotation.",
    "  CAF3 score = COL1A1/COL1A2/COL3A1, calculated only after identity assignment.",
    "",
    paste0(
      "Final CAF clusters: ",
      ifelse(
        length(CAF_CLUSTERS) == 0,
        "NONE",
        paste(
          CAF_CLUSTERS,
          collapse = ", "
        )
      )
    ),
    paste0(
      "CAF cluster source: ",
      CAF_SOURCE
    ),
    paste0(
      "CAF3 matched mouse genes: ",
      length(CAF3_MOUSE),
      "/",
      length(CAF3_HUMAN)
    ),
    paste0(
      "RAW_STEP_UP matched mouse genes: ",
      length(STEP_UP_MOUSE),
      "/",
      length(RAW_STEP_UP_HUMAN)
    ),
    "",
    "Inspect first:",
    "  01_CAF_annotation/03_YFPnegative_global_cluster_article_aligned_annotation_QC.csv",
    "  01_CAF_annotation/01_global_UMAP_canonical_celltype_markers.pdf",
    "  01_CAF_annotation/02_global_cluster_marker_DotPlot.pdf",
    "  01_CAF_annotation/04_final_CAF_clusters_used.csv",
    "",
    "Main Fig.3 candidate outputs:",
    "  03_Fig3_UMAP/01_Fig3_global_UMAP_E_H_M_Fibroblast_CAF.pdf",
    "  02_scores/01_CAF3_and_RAW_STEP_UP_Zmean_UMAP.pdf",
    "  04_Fig3_violin/01_Fig3_CAF3_STEP_UP_violin_Zmean.pdf"
  ),
  file.path(
    DIR_QC,
    "09_README_results_v1_0.txt"
  )
)

message("============================================================")
message("DONE: GSE234265 Fig.3 single-cell CAF/STEP validation v1.1")
message("Output directory: ", OUT_DIR)
message("Final CAF clusters: ",
        ifelse(
          length(CAF_CLUSTERS) == 0,
          "NONE",
          paste(CAF_CLUSTERS, collapse = ", ")
        ))
message("============================================================")
