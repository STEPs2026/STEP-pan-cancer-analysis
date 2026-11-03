# STEP manuscript analysis.
# Replace path/to placeholders with your local input directories before running.
options(stringsAsFactors = FALSE)
set.seed(1234)

required_pkgs <- c(
  "Seurat", "Matrix", "dplyr", "ggplot2",
  "patchwork", "tidyr", "stringr", "scales", "ggpubr"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(stringr)
  library(scales)
  library(ggpubr)
})

BASE_DIR <- "path/to/GSE173958"

ANALYSIS_DIR <- file.path(
  BASE_DIR,
  "终_GSE173958_M1_removeRiboCluster4_pseudoEMT_v3_4_dynamicHeatmap_forcedOrder"
)

WINDOWS_DIR <- file.path(
  ANALYSIS_DIR,
  "GSE173958_M1_removeRiboCluster4_pseudoEMT_v3_4_dynamicHeatmap_forcedOrder_Windows_step26plus"
)

RDS_CANDIDATES <- c(
  file.path(
    WINDOWS_DIR,
    "RDS",
    "GSE173958_M1_Windows_step26plus_Seurat_with_dynamicEMT_DEG.rds"
  ),
  file.path(
    ANALYSIS_DIR,
    "RDS",
    "GSE173958_M1_checkpoint_step25_Seurat_with_pseudoEMT_EMTstate.rds"
  ),
  file.path(
    ANALYSIS_DIR,
    "RDS",
    "GSE173958_M1_finalTranscriptomeFilter_with_pseudoEMT_initial_EMT_state_DEG_ready.rds"
  )
)

GMT_FILE <- file.path(
  BASE_DIR,
  "gene set",
  "STEP_all_human_gene_sets_merged.gmt"
)

RUN_PSEUDOEMT_TREND <- TRUE

UMAP_WIDTH <- 6.2
UMAP_HEIGHT <- 6.2
UMAP_POINT_SIZE <- 0.38
UMAP_POINT_ALPHA <- 0.88
UMAP_PAD_FRACTION <- 0.04
PNG_DPI <- 300

CAF3_VIOLIN_ZOOM <- TRUE
CAF3_ZOOM_QUANTILES <- c(0.01, 0.99)
CAF3_ZOOM_LOWER_PAD <- 0.08  

EHM_COLORS <- c(
  "Epithelial"  = "#2C7FB8",
  "Hybrid"      = "#F39C34",
  "Mesenchymal" = "#D73027"
)

EHM_LEVELS <- c(
  "Epithelial",
  "Hybrid",
  "Mesenchymal"
)

FOUR_COLORS <- c(
  "Epithelial"   = "#2C7FB8",
  "Early_hybrid" = "#8073AC",
  "Late_hybrid"  = "#F39C34",
  "Mesenchymal"  = "#D73027"
)

FOUR_LEVELS <- c(
  "Epithelial",
  "Early_hybrid",
  "Late_hybrid",
  "Mesenchymal"
)

CONTINUOUS_LOW <- "#F7FBFF"
CONTINUOUS_MID <- "#6BAED6"
CONTINUOUS_HIGH <- "#08306B"

RDS_HITS <- RDS_CANDIDATES[file.exists(RDS_CANDIDATES)]

if (length(RDS_HITS) == 0) {
  stop(
    "Cannot find the expected GSE173958 Seurat RDS.\nChecked:\n",
    paste(RDS_CANDIDATES, collapse = "\n")
  )
}

INPUT_RDS <- RDS_HITS[1]

if (!file.exists(GMT_FILE)) {
  stop("Cannot find GMT file: ", GMT_FILE)
}

PARENT_OUT <- if (dir.exists(WINDOWS_DIR)) WINDOWS_DIR else ANALYSIS_DIR

create_new_output_dir <- function(parent, folder) {
  if (!dir.exists(parent)) stop("Output parent does not exist: ", parent)
  run_id <- 0L
  repeat {
    suffix <- if (run_id == 0L) "" else sprintf("_run%03d", run_id)
    candidate <- file.path(parent, paste0(folder, suffix))
    if (!file.exists(candidate) &&
        dir.create(candidate, recursive = FALSE, showWarnings = FALSE)) {
      return(candidate)
    }
    if (!file.exists(candidate)) {
      stop("Cannot create output directory: ", candidate,
           ". Check permissions and free disk space.")
    }
    run_id <- run_id + 1L
    if (run_id > 9999L) stop("Too many existing v1.3 output folders.")
  }
}

OUT_DIR <- create_new_output_dir(
  PARENT_OUT,
  "Fig3_GSE173958_STEP_CAF_validation_v1_3_UMAP_axes_CAF3_zoom"
)

DIR_QC <- file.path(OUT_DIR, "00_QC")
DIR_UMAP <- file.path(OUT_DIR, "01_UMAP")
DIR_VLN <- file.path(OUT_DIR, "02_violin")
DIR_TREND <- file.path(OUT_DIR, "03_pseudoEMT_trend")
DIR_STATS <- file.path(OUT_DIR, "04_statistics")
DIR_TABLE <- file.path(OUT_DIR, "05_tables")
DIR_RDS <- file.path(OUT_DIR, "06_RDS")

for (d in c(
  OUT_DIR, DIR_QC, DIR_UMAP, DIR_VLN,
  DIR_TREND, DIR_STATS, DIR_TABLE, DIR_RDS
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(
  DIR_QC,
  "run_log_GSE173958_Fig3_validation_v1_3.txt"
)

log_msg <- function(...) {
  z <- paste0(...)
  message(z)
  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ", z, "\n",
    file = LOG_FILE,
    append = TRUE,
    sep = ""
  )
}

log_msg("INPUT_RDS = ", INPUT_RDS)
log_msg("GMT_FILE = ", GMT_FILE)
log_msg("OUT_DIR = ", OUT_DIR)
log_msg("RUN_PSEUDOEMT_TREND = ", RUN_PSEUDOEMT_TREND)

message("Reading existing Seurat object...")
seu <- readRDS(INPUT_RDS)

DefaultAssay(seu) <- "RNA"

required_meta <- c(
  "EMT_state_initial",
  "pseudoEMT"
)

missing_meta <- setdiff(
  required_meta,
  colnames(seu@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Required metadata missing from Seurat object: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (!"umap" %in% Reductions(seu)) {
  stop("Stored UMAP reduction is not present in the Seurat object.")
}

object_qc <- data.frame(
  input_RDS = INPUT_RDS,
  n_cells = ncol(seu),
  n_features = nrow(seu),
  assays = paste(Assays(seu), collapse = ";"),
  reductions = paste(Reductions(seu), collapse = ";"),
  has_EMT_state_initial = "EMT_state_initial" %in% colnames(seu@meta.data),
  has_pseudoEMT = "pseudoEMT" %in% colnames(seu@meta.data),
  has_UMAP = "umap" %in% Reductions(seu),
  stringsAsFactors = FALSE
)

write.csv(
  object_qc,
  file.path(DIR_QC, "01_input_object_QC.csv"),
  row.names = FALSE
)

print(object_qc)

get_norm_data <- function(obj, assay = "RNA") {
  x <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )

  if (is.null(x)) {
    x <- tryCatch(
      GetAssayData(obj, assay = assay, slot = "data"),
      error = function(e) NULL
    )
  }

  x
}

rna_data <- get_norm_data(seu)

if (is.null(rna_data) || nrow(rna_data) == 0 || ncol(rna_data) == 0) {
  log_msg("RNA data layer absent; running a lightweight NormalizeData fallback.")

  if (inherits(seu[["RNA"]], "Assay5")) {
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }

  seu <- NormalizeData(
    seu,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  rna_data <- get_norm_data(seu)
}

if (is.null(rna_data) || nrow(rna_data) == 0) {
  stop("Could not obtain normalized RNA data.")
}

seu$EMT_state_initial <- factor(
  as.character(seu$EMT_state_initial),
  levels = FOUR_LEVELS
)

seu$Fig3_EHM_group <- dplyr::case_when(
  as.character(seu$EMT_state_initial) == "Epithelial" ~ "Epithelial",
  as.character(seu$EMT_state_initial) %in%
    c("Early_hybrid", "Late_hybrid") ~ "Hybrid",
  as.character(seu$EMT_state_initial) == "Mesenchymal" ~ "Mesenchymal",
  TRUE ~ NA_character_
)

seu$Fig3_EHM_group <- factor(
  seu$Fig3_EHM_group,
  levels = EHM_LEVELS
)

counts_4 <- as.data.frame(
  table(
    seu$EMT_state_initial,
    useNA = "ifany"
  )
)
colnames(counts_4) <- c(
  "EMT_state_initial",
  "cell_number"
)

counts_3 <- as.data.frame(
  table(
    seu$Fig3_EHM_group,
    useNA = "ifany"
  )
)
colnames(counts_3) <- c(
  "Fig3_EHM_group",
  "cell_number"
)

write.csv(
  counts_4,
  file.path(DIR_QC, "02_original_four_state_cell_counts.csv"),
  row.names = FALSE
)

write.csv(
  counts_3,
  file.path(DIR_QC, "03_merged_E_H_M_cell_counts.csv"),
  row.names = FALSE
)

gmt_lines <- readLines(
  GMT_FILE,
  warn = FALSE
)

parse_raw_step_up <- function(lines) {
  hit <- lines[
    grepl(
      "^RAW_STEP_UP(\\t|\\s)",
      lines
    )
  ]

  if (length(hit) != 1) {
    stop(
      "RAW_STEP_UP was not uniquely found in GMT. Matches = ",
      length(hit)
    )
  }

  fields <- strsplit(
    hit,
    "\t",
    fixed = TRUE
  )[[1]]

  first_tokens <- unlist(
    strsplit(
      trimws(fields[1]),
      "\\s+"
    )
  )

  genes_first <- first_tokens[
    first_tokens != "RAW_STEP_UP"
  ]

  if (length(genes_first) == 0) {
    genes_rest <- if (length(fields) >= 3) {
      fields[3:length(fields)]
    } else {
      character(0)
    }
  } else {

    genes_rest <- if (length(fields) >= 2) {
      fields[2:length(fields)]
    } else {
      character(0)
    }
  }

  genes <- unique(
    trimws(
      c(
        genes_first,
        genes_rest
      )
    )
  )

  genes[nzchar(genes)]
}

parse_caf3 <- function(lines) {
  hit <- lines[
    grepl(
      "^CAF_3gene_score\\t",
      lines
    )
  ]

  if (length(hit) != 1) {

    warning(
      "CAF_3gene_score row was not uniquely found. ",
      "Using project CAF3 definition COL1A1/COL1A2/COL3A1."
    )
    return(
      c(
        "COL1A1",
        "COL1A2",
        "COL3A1"
      )
    )
  }

  fields <- strsplit(
    hit,
    "\t",
    fixed = TRUE
  )[[1]]

  genes <- if (length(fields) >= 5) {
    fields[3:length(fields)]
  } else {
    c(
      "COL1A1",
      "COL1A2",
      "COL3A1"
    )
  }

  genes <- toupper(
    trimws(
      genes
    )
  )

  expected <- c(
    "COL1A1",
    "COL1A2",
    "COL3A1"
  )

  if (all(expected %in% genes)) {
    expected
  } else {
    warning(
      "CAF3 GMT content did not contain all expected genes. ",
      "Using project CAF3 definition COL1A1/COL1A2/COL3A1."
    )
    expected
  }
}

RAW_STEP_UP_HUMAN <- parse_raw_step_up(
  gmt_lines
)

CAF3_HUMAN <- parse_caf3(
  gmt_lines
)

parsed_gmt <- data.frame(
  gene_set = c(
    "RAW_STEP_UP",
    "CAF_3gene_score"
  ),
  parsed_n = c(
    length(RAW_STEP_UP_HUMAN),
    length(CAF3_HUMAN)
  ),
  expected_n = c(
    84,
    3
  ),
  genes = c(
    paste(
      RAW_STEP_UP_HUMAN,
      collapse = ";"
    ),
    paste(
      CAF3_HUMAN,
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  parsed_gmt,
  file.path(DIR_QC, "04_GMT_target_gene_sets_parsed.csv"),
  row.names = FALSE
)

if (length(RAW_STEP_UP_HUMAN) != 84) {
  warning(
    "RAW_STEP_UP parsed gene count = ",
    length(RAW_STEP_UP_HUMAN),
    " rather than 84. Analysis continues; inspect QC CSV."
  )
}

map_human_to_mouse_features <- function(
    human_genes,
    mouse_features
) {
  feature_upper <- toupper(
    mouse_features
  )

  keep <- !duplicated(
    feature_upper
  )

  fmap <- setNames(
    mouse_features[keep],
    feature_upper[keep]
  )

  matched <- unname(
    fmap[
      toupper(
        human_genes
      )
    ]
  )

  data.frame(
    human_gene = human_genes,
    mouse_feature = matched,
    matched = !is.na(matched),
    stringsAsFactors = FALSE
  )
}

step_map <- map_human_to_mouse_features(
  RAW_STEP_UP_HUMAN,
  rownames(seu)
)

caf3_map <- map_human_to_mouse_features(
  CAF3_HUMAN,
  rownames(seu)
)

write.csv(
  step_map,
  file.path(
    DIR_QC,
    "05_RAW_STEP_UP_human_to_mouse_feature_matching.csv"
  ),
  row.names = FALSE
)

write.csv(
  caf3_map,
  file.path(
    DIR_QC,
    "06_CAF3_human_to_mouse_feature_matching.csv"
  ),
  row.names = FALSE
)

STEP_MOUSE <- unique(
  step_map$mouse_feature[
    step_map$matched
  ]
)

CAF3_MOUSE <- unique(
  caf3_map$mouse_feature[
    caf3_map$matched
  ]
)

matching_summary <- data.frame(
  gene_set = c(
    "RAW_STEP_UP",
    "CAF3"
  ),
  input_human_n = c(
    length(RAW_STEP_UP_HUMAN),
    length(CAF3_HUMAN)
  ),
  matched_mouse_n = c(
    length(STEP_MOUSE),
    length(CAF3_MOUSE)
  ),
  matched_mouse_genes = c(
    paste(STEP_MOUSE, collapse = ";"),
    paste(CAF3_MOUSE, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  matching_summary,
  file.path(
    DIR_QC,
    "07_score_gene_matching_summary.csv"
  ),
  row.names = FALSE
)

print(matching_summary)

if (length(STEP_MOUSE) < 10) {
  stop(
    "Too few RAW_STEP_UP genes matched the mouse object. ",
    "Inspect 05_RAW_STEP_UP_human_to_mouse_feature_matching.csv"
  )
}

if (length(CAF3_MOUSE) < 2) {
  stop(
    "Fewer than 2 CAF3 genes matched the mouse object."
  )
}

calc_zmean_from_norm <- function(
    norm_matrix,
    genes
) {
  genes_use <- intersect(
    genes,
    rownames(norm_matrix)
  )

  if (length(genes_use) == 0) {
    return(
      rep(
        NA_real_,
        ncol(norm_matrix)
      )
    )
  }

  m <- as.matrix(
    norm_matrix[
      genes_use,
      ,
      drop = FALSE
    ]
  )

  mu <- rowMeans(
    m,
    na.rm = TRUE
  )

  sdev <- apply(
    m,
    1,
    sd,
    na.rm = TRUE
  )

  sdev[
    !is.finite(sdev) |
      sdev == 0
  ] <- 1

  mz <- sweep(
    m,
    1,
    mu,
    FUN = "-"
  )

  mz <- sweep(
    mz,
    1,
    sdev,
    FUN = "/"
  )

  colMeans(
    mz,
    na.rm = TRUE
  )
}

seu$RAW_STEP_UP_Zmean <- calc_zmean_from_norm(
  rna_data,
  STEP_MOUSE
)

seu$CAF3_Zmean <- calc_zmean_from_norm(
  rna_data,
  CAF3_MOUSE
)

seu <- AddModuleScore(
  object = seu,
  features = list(
    STEP_MOUSE
  ),
  assay = "RNA",
  name = "RAW_STEP_UP_AMS"
)

seu <- AddModuleScore(
  object = seu,
  features = list(
    CAF3_MOUSE
  ),
  assay = "RNA",
  name = "CAF3_AMS"
)

if (!"RAW_STEP_UP_AMS1" %in%
    colnames(seu@meta.data)) {
  stop("RAW_STEP_UP_AMS1 was not created.")
}

if (!"CAF3_AMS1" %in%
    colnames(seu@meta.data)) {
  stop("CAF3_AMS1 was not created.")
}

umap_mat <- Embeddings(
  seu,
  reduction = "umap"
)

if (ncol(umap_mat) < 2) {
  stop("Stored UMAP does not contain two coordinates.")
}

umap_df <- data.frame(
  cell_barcode = rownames(umap_mat),
  UMAP_1 = as.numeric(umap_mat[, 1]),
  UMAP_2 = as.numeric(umap_mat[, 2]),
  stringsAsFactors = FALSE
)

meta_export <- seu@meta.data
meta_export$cell_barcode <- rownames(
  meta_export
)

plot_df <- umap_df %>%
  left_join(
    meta_export,
    by = "cell_barcode"
  )

xr <- range(
  plot_df$UMAP_1,
  finite = TRUE
)

yr <- range(
  plot_df$UMAP_2,
  finite = TRUE
)

x_center <- mean(xr)
y_center <- mean(yr)

common_span <- max(
  diff(xr),
  diff(yr)
)

common_span <- common_span *
  (1 + 2 * UMAP_PAD_FRACTION)

FIXED_XLIM <- c(
  x_center - common_span / 2,
  x_center + common_span / 2
)

FIXED_YLIM <- c(
  y_center - common_span / 2,
  y_center + common_span / 2
)

umap_geometry_qc <- data.frame(
  raw_x_min = xr[1],
  raw_x_max = xr[2],
  raw_y_min = yr[1],
  raw_y_max = yr[2],
  fixed_x_min = FIXED_XLIM[1],
  fixed_x_max = FIXED_XLIM[2],
  fixed_y_min = FIXED_YLIM[1],
  fixed_y_max = FIXED_YLIM[2],
  fixed_span_x = diff(FIXED_XLIM),
  fixed_span_y = diff(FIXED_YLIM),
  output_width_in = UMAP_WIDTH,
  output_height_in = UMAP_HEIGHT,
  point_size = UMAP_POINT_SIZE,
  stringsAsFactors = FALSE
)

write.csv(
  umap_geometry_qc,
  file.path(
    DIR_QC,
    "08_fixed_UMAP_geometry_settings.csv"
  ),
  row.names = FALSE
)

theme_umap_fixed <- function() {
  theme_classic(
    base_size = 11
  ) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 11, color = "black"),
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 8
        )
      ),
      legend.title = element_text(
        size = 9
      ),
      legend.text = element_text(
        size = 8
      ),
      legend.position = "right",
      plot.margin = margin(
        6, 6, 6, 6
      )
    )
}

add_fixed_umap_coords <- function(p) {
  p +
    coord_fixed(
      ratio = 1,
      xlim = FIXED_XLIM,
      ylim = FIXED_YLIM,
      expand = FALSE,
      clip = "on"
    ) +
    theme_umap_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2")
}

save_umap_both <- function(
    p,
    filename_stem
) {
  print(p)

  ggsave(
    file.path(
      DIR_UMAP,
      paste0(
        filename_stem,
        ".pdf"
      )
    ),
    p,
    width = UMAP_WIDTH,
    height = UMAP_HEIGHT,
    limitsize = FALSE,
    useDingbats = FALSE
  )

  ggsave(
    file.path(
      DIR_UMAP,
      paste0(
        filename_stem,
        ".png"
      )
    ),
    p,
    width = UMAP_WIDTH,
    height = UMAP_HEIGHT,
    dpi = PNG_DPI,
    limitsize = FALSE
  )
}

p_ehm <- ggplot(
  plot_df[
    !is.na(
      plot_df$Fig3_EHM_group
    ),
    ,
    drop = FALSE
  ],
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Fig3_EHM_group
  )
) +
  geom_point(
    size = UMAP_POINT_SIZE,
    alpha = UMAP_POINT_ALPHA
  ) +
  scale_color_manual(
    values = EHM_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE173958: Epithelial / Hybrid / Mesenchymal",
    color = NULL
  )

p_ehm <- add_fixed_umap_coords(
  p_ehm
)

save_umap_both(
  p_ehm,
  "01_UMAP_Epithelial_Hybrid_Mesenchymal_FIXED"
)

p_4state <- ggplot(
  plot_df[
    !is.na(
      plot_df$EMT_state_initial
    ),
    ,
    drop = FALSE
  ],
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = EMT_state_initial
  )
) +
  geom_point(
    size = UMAP_POINT_SIZE,
    alpha = UMAP_POINT_ALPHA
  ) +
  scale_color_manual(
    values = FOUR_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE173958: original four-state EMT annotation",
    color = NULL
  )

p_4state <- add_fixed_umap_coords(
  p_4state
)

save_umap_both(
  p_4state,
  "02_UMAP_Epithelial_EarlyH_LateH_Mesenchymal_FIXED"
)

make_score_umap <- function(
    df,
    score_col,
    title_text
) {
  dd <- df[
    is.finite(
      df[[score_col]]
    ),
    ,
    drop = FALSE
  ]

  p <- ggplot(
    dd,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = .data[[score_col]]
    )
  ) +
    geom_point(
      size = UMAP_POINT_SIZE,
      alpha = UMAP_POINT_ALPHA
    ) +
    scale_color_gradient2(
      low = CONTINUOUS_LOW,
      mid = CONTINUOUS_MID,
      high = CONTINUOUS_HIGH,
      midpoint = median(
        dd[[score_col]],
        na.rm = TRUE
      ),
      name = "Z-mean"
    ) +
    labs(
      title = title_text
    )

  add_fixed_umap_coords(
    p
  )
}

p_step_umap <- make_score_umap(
  plot_df,
  "RAW_STEP_UP_Zmean",
  "RAW_STEP_UP Z-mean"
)

p_caf3_umap <- make_score_umap(
  plot_df,
  "CAF3_Zmean",
  "CAF3 Z-mean"
)

save_umap_both(
  p_step_umap,
  "03_UMAP_RAW_STEP_UP_Zmean_FIXED"
)

save_umap_both(
  p_caf3_umap,
  "04_UMAP_CAF3_Zmean_FIXED"
)

p_score_pair <- p_step_umap +
  p_caf3_umap +
  patchwork::plot_layout(
    ncol = 2
  )

print(p_score_pair)

ggsave(
  file.path(
    DIR_UMAP,
    "05_UMAP_RAW_STEP_UP_CAF3_Zmean_pair_FIXED.pdf"
  ),
  p_score_pair,
  width = UMAP_WIDTH * 2,
  height = UMAP_HEIGHT,
  limitsize = FALSE,
  useDingbats = FALSE
)

ggsave(
  file.path(
    DIR_UMAP,
    "05_UMAP_RAW_STEP_UP_CAF3_Zmean_pair_FIXED.png"
  ),
  p_score_pair,
  width = UMAP_WIDTH * 2,
  height = UMAP_HEIGHT,
  dpi = PNG_DPI,
  limitsize = FALSE
)

p_step_ams_umap <- make_score_umap(
  plot_df,
  "RAW_STEP_UP_AMS1",
  "RAW_STEP_UP AddModuleScore"
)

p_caf3_ams_umap <- make_score_umap(
  plot_df,
  "CAF3_AMS1",
  "CAF3 AddModuleScore"
)

p_ams_pair <- p_step_ams_umap +
  p_caf3_ams_umap +
  patchwork::plot_layout(
    ncol = 2
  )

print(p_ams_pair)

ggsave(
  file.path(
    DIR_UMAP,
    "06_UMAP_RAW_STEP_UP_CAF3_AddModuleScore_pair_FIXED.pdf"
  ),
  p_ams_pair,
  width = UMAP_WIDTH * 2,
  height = UMAP_HEIGHT,
  limitsize = FALSE,
  useDingbats = FALSE
)

ggsave(
  file.path(
    DIR_UMAP,
    "06_UMAP_RAW_STEP_UP_CAF3_AddModuleScore_pair_FIXED.png"
  ),
  p_ams_pair,
  width = UMAP_WIDTH * 2,
  height = UMAP_HEIGHT,
  dpi = PNG_DPI,
  limitsize = FALSE
)

format_p <- function(p) {
  if (length(p) == 0 ||
      is.na(p) ||
      !is.finite(p)) {
    return("NA")
  }

  if (p < 0 || p > 1) return("NA")
  if (p < 1e-300) return("<1e-300")

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

format_p_label <- function(p, prefix = "p", spaced = TRUE) {
  value <- format_p(p)
  gap <- if (spaced) " " else ""
  if (identical(value, "<1e-300")) {
    paste0(prefix, gap, "<", gap, "1e-300")
  } else {
    paste0(prefix, gap, "=", gap, value)
  }
}

run_group_stats <- function(
    df,
    score_col,
    group_col
) {
  d <- df[
    is.finite(
      df[[score_col]]
    ) &
      !is.na(
        df[[group_col]]
      ),
    ,
    drop = FALSE
  ]

  d[[group_col]] <- droplevels(
    factor(
      d[[group_col]]
    )
  )

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
    !is.na(
      pw_df$BH_adjusted_p
    ),
    ,
    drop = FALSE
  ]

  list(
    kw = data.frame(
      score = score_col,
      group = group_col,
      Kruskal_Wallis_chisq =
        unname(
          kw$statistic
        ),
      df =
        unname(
          kw$parameter
        ),
      p_value =
        kw$p.value,
      stringsAsFactors = FALSE
    ),
    pairwise = pw_df
  )
}

scores_to_test <- c(
  "RAW_STEP_UP_Zmean",
  "CAF3_Zmean",
  "RAW_STEP_UP_AMS1",
  "CAF3_AMS1"
)

stats_main <- list()
stats_four <- list()

for (sc in scores_to_test) {
  stats_main[[sc]] <- run_group_stats(
    plot_df,
    sc,
    "Fig3_EHM_group"
  )

  stats_four[[sc]] <- run_group_stats(
    plot_df,
    sc,
    "EMT_state_initial"
  )

  write.csv(
    stats_main[[sc]]$kw,
    file.path(
      DIR_STATS,
      paste0(
        sc,
        "_EHM_Kruskal_Wallis.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    stats_main[[sc]]$pairwise,
    file.path(
      DIR_STATS,
      paste0(
        sc,
        "_EHM_pairwise_Wilcoxon_BH.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    stats_four[[sc]]$kw,
    file.path(
      DIR_STATS,
      paste0(
        sc,
        "_4state_Kruskal_Wallis.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    stats_four[[sc]]$pairwise,
    file.path(
      DIR_STATS,
      paste0(
        sc,
        "_4state_pairwise_Wilcoxon_BH.csv"
      )
    ),
    row.names = FALSE
  )
}

CAF3_ZOOM_QC <- new.env(parent = emptyenv())
CAF3_ZOOM_QC$rows <- list()

get_caf3_zoom <- function(df, dd, score_col, group_col, n_brackets = 0L) {
  if (!isTRUE(CAF3_VIOLIN_ZOOM) ||
      !score_col %in% c("CAF3_Zmean", "CAF3_AMS1")) return(NULL)
  probs <- CAF3_ZOOM_QUANTILES
  if (length(probs) != 2L || any(!is.finite(probs)) ||
      probs[1] < 0 || probs[2] > 1 || probs[1] >= probs[2]) {
    stop("CAF3_ZOOM_QUANTILES must contain two increasing probabilities.")
  }
  x <- dd[[score_col]]
  if (!length(x) || any(!is.finite(x))) {
    stop("No valid CAF3 cells for ", score_col, " / ", group_col,
         "; inspect input scores and state annotations.")
  }
  q <- unname(quantile(x, probs = probs, type = 7))
  fallback <- "none"
  base <- q
  span <- diff(base)
  if (!is.finite(span) || span <= 0) {

    span <- max(abs(q[1]), 1) * 0.02
    base <- q[1] + c(-0.5, 0.5) * span
    fallback <- "coincident_quantiles_symmetric_window"
    warning("Coincident CAF3 quantiles for ", score_col, " / ", group_col,
            "; using recorded symmetric display window.")
  }
  headroom <- if (n_brackets > 0L) (n_brackets + 0.8) * 0.11 * span else 0
  lower_padding <- CAF3_ZOOM_LOWER_PAD * span
  ylim <- c(base[1] - lower_padding, base[2] + headroom)
  positions <- base[2] + seq_len(n_brackets) * 0.11 * span
  known <- !is.na(df[[group_col]])
  valid <- is.finite(df[[score_col]])
  qc <- data.frame(
    score = score_col, grouping = group_col,
    n_input = nrow(df), n_missing_group = sum(!known),
    n_nonfinite_score = sum(!valid), n_ineligible_union = sum(!known | !valid),
    n_used_for_violin = length(x), n_used_for_statistics = length(x),
    n_removed_by_zoom = 0L,
    raw_min = min(x), raw_max = max(x),
    quantile_low_prob = probs[1], quantile_high_prob = probs[2],
    quantile_type = 7L, quantile_low = q[1], quantile_high = q[2],
    base_y_min = base[1], base_y_max = base[2],
    lower_display_padding = lower_padding,
    annotation_headroom = headroom, n_brackets = n_brackets,
    bracket_y_positions = paste(positions, collapse = ";"),
    actual_y_min = ylim[1], actual_y_max = ylim[2],
    n_below_quantile = sum(x < q[1]), n_above_quantile = sum(x > q[2]),
    n_below_display = sum(x < ylim[1]), n_above_display = sum(x > ylim[2]),
    fallback = fallback, zoom_method = "coord_cartesian; no y scale limits",
    stringsAsFactors = FALSE
  )
  CAF3_ZOOM_QC$rows[[paste(score_col, group_col, sep = "__")]] <- qc
  list(ylim = ylim, positions = positions)
}

caf3_zoom_caption <- function(has_brackets = FALSE) {
  paste0(
    "CAF3 display: 1st-99th percentiles + lower padding",
    if (has_brackets) " + bracket space." else ".",
    "\nViolin and statistics use all valid cells; tails are clipped only."
  )
}

make_violin_ehm <- function(
    df,
    score_col,
    ylab_text,
    title_text
) {
  dd <- df %>%
    filter(
      !is.na(
        Fig3_EHM_group
      ),
      is.finite(
        .data[[score_col]]
      )
    ) %>%
    mutate(
      Fig3_EHM_group = factor(
        Fig3_EHM_group,
        levels = EHM_LEVELS
      )
    )

  st <- stats_main[[score_col]]

  y_min <- min(
    dd[[score_col]],
    na.rm = TRUE
  )

  y_max <- max(
    dd[[score_col]],
    na.rm = TRUE
  )

  y_range <- y_max - y_min

  if (!is.finite(y_range) ||
      y_range == 0) {
    y_range <- 1
  }

  pw <- st$pairwise

  pw <- pw[
    pw$group1 %in% EHM_LEVELS &
      pw$group2 %in% EHM_LEVELS,
    ,
    drop = FALSE
  ]

  if (nrow(pw) > 0) {
    pw$label <- vapply(
      pw$BH_adjusted_p, format_p_label, character(1),
      prefix = "BH p", spaced = FALSE
    )

    pw$y.position <- y_max +
      seq_len(
        nrow(pw)
      ) *
      0.11 *
      y_range
  }

  zoom <- get_caf3_zoom(
    df, dd, score_col, "Fig3_EHM_group", n_brackets = nrow(pw)
  )
  if (!is.null(zoom) && nrow(pw) > 0L) pw$y.position <- zoom$positions

  p <- ggplot(
    dd,
    aes(
      x = Fig3_EHM_group,
      y = .data[[score_col]],
      fill = Fig3_EHM_group
    )
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      alpha = 0.82,
      color = "grey30",
      linewidth = 0.28
    ) +
    geom_boxplot(
      width = 0.14,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_manual(
      values = EHM_COLORS,
      drop = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = format_p_label(st$kw$p_value, "Kruskal-Wallis p"),
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
        angle = 0,
        hjust = 0.5
      )
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = if (!is.null(zoom)) c(0, 0) else c(
          0.05,
          ifelse(nrow(pw) > 0, 0.35, 0.10)
        )
      )
    )

  if (nrow(pw) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        pw,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.30,
        size = 2.7
      )
  }

  if (!is.null(zoom)) {
    p <- p +
      coord_cartesian(ylim = zoom$ylim, clip = "on") +
      labs(caption = caf3_zoom_caption(nrow(pw) > 0L)) +
      theme(plot.caption = element_text(size = 7.5, color = "grey35", hjust = 0))
  }

  p
}

p_vln_step <- make_violin_ehm(
  plot_df,
  "RAW_STEP_UP_Zmean",
  "RAW_STEP_UP Z-mean",
  "RAW_STEP_UP across E / H / M"
)

p_vln_caf3 <- make_violin_ehm(
  plot_df,
  "CAF3_Zmean",
  "CAF3 Z-mean",
  "CAF3 across E / H / M"
)

p_vln_pair <- p_vln_step +
  p_vln_caf3 +
  patchwork::plot_layout(
    ncol = 2
  ) +
  patchwork::plot_annotation(
    title = "GSE173958 tumour-cell EMT continuum"
  )

print(p_vln_pair)

ggsave(
  file.path(
    DIR_VLN,
    "01_Fig3_RAW_STEP_UP_CAF3_Zmean_E_H_M_violin.pdf"
  ),
  p_vln_pair,
  width = 11.5,
  height = 5.5,
  limitsize = FALSE,
  useDingbats = FALSE
)

ggsave(
  file.path(
    DIR_VLN,
    "01_Fig3_RAW_STEP_UP_CAF3_Zmean_E_H_M_violin.png"
  ),
  p_vln_pair,
  width = 11.5,
  height = 5.5,
  dpi = PNG_DPI,
  limitsize = FALSE
)

p_vln_step_ams <- make_violin_ehm(
  plot_df,
  "RAW_STEP_UP_AMS1",
  "RAW_STEP_UP AddModuleScore",
  "RAW_STEP_UP AddModuleScore across E / H / M"
)

p_vln_caf3_ams <- make_violin_ehm(
  plot_df,
  "CAF3_AMS1",
  "CAF3 AddModuleScore",
  "CAF3 AddModuleScore across E / H / M"
)

p_vln_ams_pair <- p_vln_step_ams +
  p_vln_caf3_ams +
  patchwork::plot_layout(
    ncol = 2
  ) +
  patchwork::plot_annotation(
    title = "Sensitivity analysis: AddModuleScore"
  )

print(p_vln_ams_pair)

ggsave(
  file.path(
    DIR_VLN,
    "02_RAW_STEP_UP_CAF3_AddModuleScore_E_H_M_violin.pdf"
  ),
  p_vln_ams_pair,
  width = 11.5,
  height = 5.5,
  limitsize = FALSE,
  useDingbats = FALSE
)

ggsave(
  file.path(
    DIR_VLN,
    "02_RAW_STEP_UP_CAF3_AddModuleScore_E_H_M_violin.png"
  ),
  p_vln_ams_pair,
  width = 11.5,
  height = 5.5,
  dpi = PNG_DPI,
  limitsize = FALSE
)

make_violin_four <- function(
    df,
    score_col,
    ylab_text,
    title_text
) {
  dd <- df %>%
    filter(
      !is.na(
        EMT_state_initial
      ),
      is.finite(
        .data[[score_col]]
      )
    ) %>%
    mutate(
      EMT_state_initial = factor(
        EMT_state_initial,
        levels = FOUR_LEVELS
      )
    )

  st <- stats_four[[score_col]]
  zoom <- get_caf3_zoom(df, dd, score_col, "EMT_state_initial")

  p <- ggplot(
    dd,
    aes(
      x = EMT_state_initial,
      y = .data[[score_col]],
      fill = EMT_state_initial
    )
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      alpha = 0.82,
      color = "grey30",
      linewidth = 0.28
    ) +
    geom_boxplot(
      width = 0.13,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_manual(
      values = FOUR_COLORS,
      drop = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = format_p_label(st$kw$p_value, "Kruskal-Wallis p"),
      x = NULL,
      y = ylab_text
    ) +
    theme_classic(
      base_size = 10
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
        angle = 28,
        hjust = 1
      )
    )
  if (!is.null(zoom)) {
    p <- p +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      coord_cartesian(ylim = zoom$ylim, clip = "on") +
      labs(caption = caf3_zoom_caption()) +
      theme(plot.caption = element_text(size = 7.5, color = "grey35", hjust = 0))
  }
  p
}

p4_step <- make_violin_four(
  plot_df,
  "RAW_STEP_UP_Zmean",
  "RAW_STEP_UP Z-mean",
  "RAW_STEP_UP across original four EMT states"
)

p4_caf3 <- make_violin_four(
  plot_df,
  "CAF3_Zmean",
  "CAF3 Z-mean",
  "CAF3 across original four EMT states"
)

p4_pair <- p4_step +
  p4_caf3 +
  patchwork::plot_layout(
    ncol = 2
  )

print(p4_pair)

ggsave(
  file.path(
    DIR_VLN,
    "03_supplementary_RAW_STEP_UP_CAF3_Zmean_4state_violin.pdf"
  ),
  p4_pair,
  width = 11.5,
  height = 5.2,
  limitsize = FALSE,
  useDingbats = FALSE
)

ggsave(
  file.path(
    DIR_VLN,
    "03_supplementary_RAW_STEP_UP_CAF3_Zmean_4state_violin.png"
  ),
  p4_pair,
  width = 11.5,
  height = 5.2,
  dpi = PNG_DPI,
  limitsize = FALSE
)

if (length(CAF3_ZOOM_QC$rows) > 0L) {
  caf3_zoom_qc <- do.call(rbind, unname(CAF3_ZOOM_QC$rows))
  write.csv(
    caf3_zoom_qc,
    file.path(DIR_QC, "11_CAF3_violin_zoom_ranges_v1_3.csv"),
    row.names = FALSE
  )
  print(caf3_zoom_qc)
}

make_pseudoemt_binned_trend <- function(
    df,
    score_col,
    ylab_text,
    title_text,
    n_bins = 40
) {
  dd <- df %>%
    filter(
      is.finite(
        pseudoEMT
      ),
      is.finite(
        .data[[score_col]]
      )
    ) %>%
    arrange(
      pseudoEMT
    )

  if (nrow(dd) < 50) {
    return(NULL)
  }

  dd$pt_bin <- dplyr::ntile(
    dd$pseudoEMT,
    n_bins
  )

  binned <- dd %>%
    group_by(pt_bin) %>%
    summarise(
      pseudoEMT_median =
        median(
          pseudoEMT,
          na.rm = TRUE
        ),
      score_median =
        median(
          .data[[score_col]],
          na.rm = TRUE
        ),
      score_q25 =
        quantile(
          .data[[score_col]],
          0.25,
          na.rm = TRUE
        ),
      score_q75 =
        quantile(
          .data[[score_col]],
          0.75,
          na.rm = TRUE
        ),
      n_cells = n(),
      .groups = "drop"
    )

  ct <- suppressWarnings(
    cor.test(
      dd$pseudoEMT,
      dd[[score_col]],
      method = "spearman",
      exact = FALSE
    )
  )

  p <- ggplot(
    binned,
    aes(
      x = pseudoEMT_median,
      y = score_median
    )
  ) +
    geom_ribbon(
      aes(
        ymin = score_q25,
        ymax = score_q75
      ),
      alpha = 0.18
    ) +
    geom_line(
      linewidth = 0.8
    ) +
    geom_point(
      size = 1.8
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Existing pseudoEMT only; Spearman rho = ",
        signif(
          unname(
            ct$estimate
          ),
          3
        ),
        ", ",
        format_p_label(ct$p.value)
      ),
      x = "Existing pseudoEMT",
      y = ylab_text
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )

  list(
    plot = p,
    binned = binned,
    correlation = data.frame(
      score = score_col,
      n_cells = nrow(dd),
      spearman_rho =
        unname(
          ct$estimate
        ),
      p_value =
        ct$p.value,
      stringsAsFactors = FALSE
    )
  )
}

if (isTRUE(
  RUN_PSEUDOEMT_TREND
)) {
  trend_step <- make_pseudoemt_binned_trend(
    plot_df,
    "RAW_STEP_UP_Zmean",
    "RAW_STEP_UP Z-mean",
    "RAW_STEP_UP along existing pseudoEMT"
  )

  trend_caf3 <- make_pseudoemt_binned_trend(
    plot_df,
    "CAF3_Zmean",
    "CAF3 Z-mean",
    "CAF3 along existing pseudoEMT"
  )

  if (!is.null(trend_step)) {
    print(trend_step$plot)

    ggsave(
      file.path(
        DIR_TREND,
        "01_RAW_STEP_UP_Zmean_existing_pseudoEMT_binned_trend.pdf"
      ),
      trend_step$plot,
      width = 6.2,
      height = 4.8,
      useDingbats = FALSE
    )

    ggsave(
      file.path(
        DIR_TREND,
        "01_RAW_STEP_UP_Zmean_existing_pseudoEMT_binned_trend.png"
      ),
      trend_step$plot,
      width = 6.2,
      height = 4.8,
      dpi = PNG_DPI
    )

    write.csv(
      trend_step$binned,
      file.path(
        DIR_TREND,
        "01_RAW_STEP_UP_Zmean_existing_pseudoEMT_binned_data.csv"
      ),
      row.names = FALSE
    )
  }

  if (!is.null(trend_caf3)) {
    print(trend_caf3$plot)

    ggsave(
      file.path(
        DIR_TREND,
        "02_CAF3_Zmean_existing_pseudoEMT_binned_trend.pdf"
      ),
      trend_caf3$plot,
      width = 6.2,
      height = 4.8,
      useDingbats = FALSE
    )

    ggsave(
      file.path(
        DIR_TREND,
        "02_CAF3_Zmean_existing_pseudoEMT_binned_trend.png"
      ),
      trend_caf3$plot,
      width = 6.2,
      height = 4.8,
      dpi = PNG_DPI
    )

    write.csv(
      trend_caf3$binned,
      file.path(
        DIR_TREND,
        "02_CAF3_Zmean_existing_pseudoEMT_binned_data.csv"
      ),
      row.names = FALSE
    )
  }

  trend_corr <- bind_rows(
    if (!is.null(trend_step)) {
      trend_step$correlation
    } else {
      NULL
    },
    if (!is.null(trend_caf3)) {
      trend_caf3$correlation
    } else {
      NULL
    }
  )

  if (nrow(trend_corr) > 0) {
    write.csv(
      trend_corr,
      file.path(
        DIR_STATS,
        "pseudoEMT_score_Spearman_correlations.csv"
      ),
      row.names = FALSE
    )
  }
}

group_summary_ehm <- plot_df %>%
  filter(
    !is.na(
      Fig3_EHM_group
    )
  ) %>%
  group_by(
    Fig3_EHM_group
  ) %>%
  summarise(
    n_cells = n(),

    RAW_STEP_UP_Zmean_mean =
      mean(
        RAW_STEP_UP_Zmean,
        na.rm = TRUE
      ),
    RAW_STEP_UP_Zmean_median =
      median(
        RAW_STEP_UP_Zmean,
        na.rm = TRUE
      ),
    RAW_STEP_UP_Zmean_sd =
      sd(
        RAW_STEP_UP_Zmean,
        na.rm = TRUE
      ),

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

    RAW_STEP_UP_AMS_mean =
      mean(
        RAW_STEP_UP_AMS1,
        na.rm = TRUE
      ),
    CAF3_AMS_mean =
      mean(
        CAF3_AMS1,
        na.rm = TRUE
      ),

    pseudoEMT_median =
      median(
        pseudoEMT,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write.csv(
  group_summary_ehm,
  file.path(
    DIR_TABLE,
    "01_E_H_M_group_score_summary.csv"
  ),
  row.names = FALSE
)

group_summary_four <- plot_df %>%
  filter(
    !is.na(
      EMT_state_initial
    )
  ) %>%
  group_by(
    EMT_state_initial
  ) %>%
  summarise(
    n_cells = n(),
    RAW_STEP_UP_Zmean_mean =
      mean(
        RAW_STEP_UP_Zmean,
        na.rm = TRUE
      ),
    RAW_STEP_UP_Zmean_median =
      median(
        RAW_STEP_UP_Zmean,
        na.rm = TRUE
      ),
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
    pseudoEMT_median =
      median(
        pseudoEMT,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write.csv(
  group_summary_four,
  file.path(
    DIR_TABLE,
    "02_original_4state_score_summary.csv"
  ),
  row.names = FALSE
)

per_cell_cols <- c(
  "cell_barcode",
  "UMAP_1",
  "UMAP_2",
  "EMT_state_initial",
  "Fig3_EHM_group",
  "pseudoEMT",
  "RAW_STEP_UP_Zmean",
  "CAF3_Zmean",
  "RAW_STEP_UP_AMS1",
  "CAF3_AMS1"
)

per_cell_export <- plot_df[
  ,
  intersect(
    per_cell_cols,
    colnames(plot_df)
  ),
  drop = FALSE
]

write.csv(
  per_cell_export,
  file.path(
    DIR_TABLE,
    "03_ALL_cells_Fig3_annotation_scores_UMAP.csv"
  ),
  row.names = FALSE
)

saveRDS(
  seu,
  file.path(
    DIR_RDS,
    "GSE173958_M1_Fig3_STEP_CAF_scored_v1_3.rds"
  )
)

manifest <- data.frame(
  item = c(
    "version",
    "input_RDS",
    "GMT",
    "output_dir",
    "cell_state_metadata",
    "main_grouping",
    "pseudoEMT",
    "primary_scoring",
    "sensitivity_scoring",
    "RAW_STEP_UP_mouse_genes_matched",
    "CAF3_mouse_genes_matched",
    "UMAP_geometry",
    "UMAP_output_size",
    "UMAP_point_size",
    "UMAP_axes",
    "pseudoEMT_trend_enabled",
    "CAF3_violin_zoom_enabled",
    "CAF3_violin_quantiles",
    "CAF3_violin_zoom_QC",
    "p_display_rule"
  ),
  value = c(
    "v1.3",
    INPUT_RDS,
    GMT_FILE,
    OUT_DIR,
    "EMT_state_initial",
    "Epithelial / Hybrid(Early_hybrid+Late_hybrid) / Mesenchymal",
    "Existing metadata only; NOT recomputed",
    "single-cell Z-mean",
    "Seurat AddModuleScore",
    paste0(
      length(STEP_MOUSE),
      "/",
      length(RAW_STEP_UP_HUMAN)
    ),
    paste0(
      length(CAF3_MOUSE),
      "/",
      length(CAF3_HUMAN)
    ),
    "same stored UMAP + equal x/y span + coord_fixed(ratio=1) for every UMAP",
    paste0(
      UMAP_WIDTH,
      " x ",
      UMAP_HEIGHT,
      " inches"
    ),
    as.character(
      UMAP_POINT_SIZE
    ),
    "Visible axis lines, ticks and numeric labels; UMAP 1 / UMAP 2 titles",
    as.character(
      RUN_PSEUDOEMT_TREND
    ),
    as.character(CAF3_VIOLIN_ZOOM),
    paste(CAF3_ZOOM_QUANTILES, collapse = ";"),
    "00_QC/11_CAF3_violin_zoom_ranges_v1_3.csv",
    "p < 1e-300 for numeric p < 1e-300; numeric tables unchanged"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(
    DIR_QC,
    "09_run_manifest_v1_3.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    "GSE173958 Fig.3 RAW_STEP_UP / CAF3 validation v1.3",
    "",
    paste0(
      "Finished: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Input RDS: ",
      INPUT_RDS
    ),
    paste0(
      "Output: ",
      OUT_DIR
    ),
    "",
    "No QC / clustering / UMAP / Monocle / pseudoEMT / graph_test / DEG was rerun.",
    "",
    "Main Fig.3 group:",
    "  Epithelial",
    "  Hybrid = Early_hybrid + Late_hybrid",
    "  Mesenchymal",
    "",
    "Primary scores:",
    "  RAW_STEP_UP_Zmean",
    "  CAF3_Zmean",
    "",
    "Sensitivity scores:",
    "  RAW_STEP_UP_AMS1",
    "  CAF3_AMS1",
    "",
    "UMAP geometry is locked across ALL panels:",
    paste0(
      "  size = ",
      UMAP_WIDTH,
      " x ",
      UMAP_HEIGHT,
      " inches"
    ),
    "  same stored coordinates",
    "  same x/y limits",
    "  equal x/y span",
    "  coord_fixed(ratio = 1)",
    paste0(
      "  point size = ",
      UMAP_POINT_SIZE
    ),
    "",
    "v1.3 changes from v1.2:",
    "  UMAP axis lines, ticks, numeric labels and UMAP 1 / UMAP 2 titles shown.",
    "  Stored coordinates, fixed limits, ratio = 1 and output dimensions retained.",
    "  Separate v1.3 output folder; earlier outputs are never overwritten.",
    "Retained v1.2 settings:",
    "  CAF3: pooled 1st-99th percentiles + 8% lower padding; E/H/M adds bracket space.",
    "  coord_cartesian only; all valid cells retained for violin and statistics.",
    "  RAW_STEP_UP geometry and all UMAP/scoring/statistical methods unchanged.",
    "  Numeric p < 1e-300 displayed as p < 1e-300; CSV numbers unchanged.",
    "  Zoom QC: 00_QC/11_CAF3_violin_zoom_ranges_v1_3.csv",
    "",
    "Main candidate outputs:",
    "  01_UMAP/01_UMAP_Epithelial_Hybrid_Mesenchymal_FIXED.pdf",
    "  01_UMAP/03_UMAP_RAW_STEP_UP_Zmean_FIXED.pdf",
    "  01_UMAP/04_UMAP_CAF3_Zmean_FIXED.pdf",
    "  01_UMAP/05_UMAP_RAW_STEP_UP_CAF3_Zmean_pair_FIXED.pdf",
    "  02_violin/01_Fig3_RAW_STEP_UP_CAF3_Zmean_E_H_M_violin.pdf",
    "",
    "pseudoEMT trend:",
    paste0(
      "  enabled = ",
      RUN_PSEUDOEMT_TREND
    ),
    "  uses only the existing pseudoEMT metadata; no Monocle rerun."
  ),
  file.path(
    DIR_QC,
    "10_README_results_v1_3.txt"
  )
)

log_msg("============================================================")
log_msg("DONE: GSE173958 Fig.3 STEP/CAF validation v1.3")
log_msg("Output directory: ", OUT_DIR)
log_msg(
  "RAW_STEP_UP mouse matched genes = ",
  length(STEP_MOUSE),
  "/",
  length(RAW_STEP_UP_HUMAN)
)
log_msg(
  "CAF3 mouse matched genes = ",
  length(CAF3_MOUSE),
  "/",
  length(CAF3_HUMAN)
)
log_msg("============================================================")
