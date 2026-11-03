# STEP: Single-cell-defined Transitional Epithelial Plasticity

This repository contains R scripts used in the STEP manuscript. The study derives a hybrid epithelial–mesenchymal transcriptional signature from single-cell RNA-seq data, examines its relationships with tumor composition and cancer-associated fibroblast (CAF) scores, and evaluates CAF-adjusted STEP scores, associated pathways, clinical outcomes, and treatment responses in independent cohorts.

Raw STEP scores, CAF scores, and CAF-adjusted STEP scores represent different quantities. Historical HEMTS terminology has been changed to STEP in this version of the code.

## 1. Scripts and figure correspondence

All 21 R scripts are stored in the repository root. The first two digits identify the manuscript Figure; the last two distinguish scripts within that Figure. Numbering follows figure organization, **not a universal execution order**. Individual scripts may also generate related supplementary figures or additional tables.

|Script|Figure / analysis|
|---|---|
|`02_01_GSE234265_clustering.R`|Figure 2: GSE234265 clustering, YFP-positive tumor selection, and EMT-state annotation|
|`02_02_GSE234265_trajectory.R`|Figure 2: GSE234265 Monocle2 trajectory and pseudotime|
|`02_03_GSE234265_differential_expression.R`|Figure 2: hybrid-state differential expression and heatmaps|
|`02_04_GSE173958_clustering.R`|Figure 2: GSE173958 M1 preprocessing, clustering, Monocle3 pseudoEMT, and checkpoint generation|
|`02_05_GSE173958_dynamic_EMT_DEG.R`|Figure 2: checkpoint-based dynamic-gene analysis and EMT differential expression|
|`03_01_STEP_GSEA.R`|Figure 3: cross-dataset/species candidate-signature enrichment|
|`03_02_STEP_GSEA_leading_edge.R`|Figure 3: leading-edge support summaries and heatmaps|
|`03_03_GSE234265_STEP_CAF3_validation.R`|Figure 3: single-cell STEP and CAF3 scoring in GSE234265|
|`03_04_GSE173958_STEP_CAF3_validation.R`|Figure 3: single-cell STEP and CAF3 scoring in GSE173958|
|`03_05_TCGA_STEP_correlations.R`|Figure 3: TCGA scoring, purity associations, and signature correlations|
|`03_06_TCGA_STEP_CAF_adjustment.R`|Figure 3: CAF-adjusted STEP and associated summary analyses; retained full pipeline contains additional modules|
|`03_07_TCGA_STEP_raw_vs_CAF_adjusted.R`|Figure 3: paired per-cancer correlations before and after CAF adjustment|
|`04_01_TCGA_STEP_Hallmark_Reactome.R`|Figure 4: Hallmark/Reactome enrichment and pan-cancer gene-correlation meta-analysis|
|`05_01_TCGA_STEP_functional_subsets_KM_cutoffs.R`|Figure 5: MYC7/EXO23 optimal survival cutpoints and Kaplan–Meier plots|
|`05_02_TCGA_STEP_functional_subsets_survival_models.R`|Figure 5: cutoff-based Cox models and survival summary figures|
|`05_03_TCGA_STEP_functional_subsets_clinical.R`|Figure 5: clinical-variable repair and functional-subset clinical associations|
|`06_01_TCGA_STEP_TIDE.R`|Figure 6: TCGA/TIDE matching, predicted response, and related summaries; retained full pipeline contains additional modules|
|`06_02_IMvigor210_GSE78220_STEP.R`|Figure 6: STEP scoring, CAF adjustment, survival, and response analyses in IMvigor210/GSE78220|
|`06_03_PRECOG_Kinget_STEP.R`|Figure 6: Kinget response and survival analyses from existing PRECOG score tables|
|`06_04_PRECOG_Gide_CTLA4_PD1_STEP.R`|Figure 6-related supplementary analysis: Gide record-level descriptive survival|
|`06_05_ISPY2_GSE194040_STEP.R`|Figure 6: I-SPY2 pathological complete response analyses, overall and by subtype|

`FILE_MAPPING.csv` records the original filenames and their numbered replacements.

## 2. R environment and dependencies

Earlier module validation used **R 4.6.0**. The terminology-renamed scripts have passed syntax/structure checks; this does not establish a complete rerun of this version.

Major dependencies include:

- **Single-cell analysis:** Seurat, Matrix, sctransform, monocle (Monocle2), monocle3, SingleCellExperiment, DDRTree, igraph, clustree.
- **Scoring/enrichment:** GSVA, fgsea, babelgene, estimate, AnnotationDbi, org.Hs.eg.db.
- **Statistics/survival:** survival, survminer, metafor, pROC.
- **Data handling/plots:** dplyr, tidyr, tibble, readr, stringr, purrr, data.table, readxl, openxlsx, ggplot2, ggpubr, patchwork, pheatmap, scales, RColorBrewer, gridExtra, broom, janitor.

Install CRAN packages with `install.packages()` and Bioconductor packages with `BiocManager::install()`. Monocle3 and the R-Forge `estimate` package require their respective installation procedures and system dependencies. Install only the dependencies required for your selected scripts, using a separate R library when possible. Package requirements are also declared within the scripts. Some scripts automatically install missing packages; review these settings before running.

The package list is a practical setup guide, not a tested environment lockfile. Record `sessionInfo()` for each analysis, and avoid upgrading packages during a reproduction run. Single-cell and dynamic-gene analyses require substantially more memory and runtime than downstream source-table plotting.

## 3. Input data and path configuration

|Analysis|Required inputs / acquisition|
|---|---|
|GSE234265|Obtain the control_2 10X matrices from GEO accession GSE234265; downstream scripts also require the annotated tumor/global Seurat objects generated by preprocessing or supplied as checkpoints.|
|GSE173958|Obtain the six M1 10X inputs and author cell-classification table associated with GEO accession GSE173958. Dynamic analysis requires step25 Seurat and Monocle3 checkpoints; score validation uses the specified processed checkpoint.|
|Candidate GSEA|Provide the original human/mouse candidate GMT files and per-dataset EMT-contrast ranking tables. These processed inputs are not downloaded automatically.|
|TCGA|Provide the original per-cancer GCT expression and clinical tables under `TCGA_201906cBioP/exp_gct/`, the purity table `TCGA_Tumor_purity_TCGAbiolinks.csv`, and the signature GMT files. Use the study's cBioPortal/TCGA snapshot and matching identifiers.|
|Hallmark/Reactome|Provide the Hallmark and Reactome GMT collections used by the study, obtained through MSigDB under its access terms; retain their original release versions. GO GMTs are required only for enabled GO analyses.|
|Functional subsets|Provide the existing `ALL_cancers_MYC7_Exosome_scores_clinical_merged.csv` within `09_pan_cancer_summary/`, together with the clinical source files required for repair. The upstream subset-scoring generator is not included among these 21 scripts.|
|TIDE|Provide the archived TIDE prediction output and the corresponding TCGA score results. Local code analyzes these exports; it does not independently recreate the online prediction service.|
|IMvigor210/GSE78220|Provide the IMvigor210 expression and patient/sample clinical files used by the study, plus the GSE78220 expression workbook and GEO series matrix.|
|PRECOG Kinget|Provide the precomputed `PRECOG_STEP_FORMAL_ANALYSIS_v1/01_MAIN_HIGH_COVERAGE/01_SCORES/` patient-score CSV files. Their upstream PRECOG scoring script is not included here.|
|Gide|Provide `PRECOG_Gide_CTLA4PD1_22rows_to_original_PFS_OS_mapping.csv`. This script analyzes 22 expression records from 16 patients and retains the original record-level grouping.|
|I-SPY2|Obtain GSE194040 gene-level expression and both GPL20078/GPL30493 series matrices from GEO, together with the signature GMT used by the script.|

Replace `path/to/...` placeholders and configure `base_dir`, `BASE_DIR`, `ROOT_DIR`, `RESULT_DIR`, or `PROJECT_DIR` as applicable. Scripts using `getwd()` expect the **data/project working directory**, which need not be the directory containing the R scripts. For I-SPY2, explicitly set `PROJECT_DIR` if its input files are not beside the script. `PATH_PLACEHOLDERS.csv` is not required in this version; inspect the path settings near the top of each numbered script.

**Naming compatibility:** the code now references STEP names, including `STEP_all_human_gene_sets_merged.gmt`, signature names, score columns, and historical result directories. Existing files or checkpoint metadata that still use HEMTS names must be adapted consistently in working copies. Renaming a file alone does not rename GMT entries or data columns. See `EXTERNAL_NAMES_AND_FIELDS.csv` for affected string references. Input data and checkpoints are not bundled with this code folder.

## 4. Execution order

Run analyses in separate R sessions where practical; several scripts clear the workspace or change the working directory. Configure inputs before calling `source()` or `Rscript`.

- **GSE234265:** run `02_01`, then `02_02` and `02_03` from its annotated tumor checkpoint. `03_03` requires both the global and annotated tumor objects.
- **GSE173958:** run `02_04` to generate the step25 checkpoints, then `02_05` if dynamic/DEG recomputation is required. `03_04` requires the processed checkpoint specified in that script. Existing compatible checkpoints can be used directly; dynamic recomputation is optional for checkpoint-based replay. The original `02_04` intentionally stops after saving step25.
- **Signature enrichment:** prepare the ranking tables and GMTs, then run `03_01` → `03_02`.
- **TCGA composition:** run `03_05` or the configured full pipeline `03_06` as required; use the resulting correlation table for `03_07`. These full pipelines overlap and should not be assumed to be interchangeable or automatically sequential.
- **TCGA pathways:** `04_01` starts from TCGA expression and GMT inputs and can run independently of the Figure 3 plots.
- **Functional subsets:** obtain the precomputed merged score table first; run `05_03` for clinical repair, then `05_01` for cutpoints/KM and `05_02` for downstream survival models. Explicitly point all three scripts to the same result directory.
- **Treatment cohorts:** `06_02`, `06_03`, `06_04`, and `06_05` are separate cohort analyses. Kinget/Gide require the processed inputs listed above. Configure `06_01` with its TCGA inputs, archived TIDE export, and optional-module settings; it is a retained full pipeline, not a small TIDE-only plotting script.

For example, after editing its paths, run a selected script with:

```r
source("02_01_GSE234265_clustering.R", encoding = "UTF-8")
```

Do not treat numerical filename order as a one-command complete pipeline.

## 5. Outputs

Outputs are written to script-defined dataset/result directories, often including numbered subdirectories or timestamps. Principal products include:

- Seurat/Monocle checkpoint objects, clustering/EMT UMAPs, trajectory plots, DEG tables, and expression heatmaps.
- STEP/CAF score tables, group comparisons, raw/adjusted correlation tables, and summary plots.
- GSEA NES/P-value/FDR and leading-edge tables; Hallmark/Reactome enrichment plots and gene-level meta-analysis results.
- Clinical association tables; survival cutpoints, Kaplan–Meier curves with risk tables where implemented, Cox hazard ratios/confidence intervals, and summary figures.
- TIDE matched-source and predicted-response tables; treatment-cohort response/survival statistics; I-SPY2 overall/subtype pCR source tables and figures.

Figures are exported in the formats implemented by each script, commonly PDF and PNG, with some TIFF heatmaps. CSV/source-data tables accompany many analyses. Exact filenames and output locations are defined within each script; they are not all redirected to a single repository-wide results folder.
