#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
})

set.seed(20260816)

coculture_file <- "input_data/culture/outputs/019fbe7c_combined_decontX_reclustered_all_cells/Combined_neutrophil_fibroblast_decontX_reclustered.rds"
preculture_file <- "input_data/culture/outputs/019fbe7c_nc1_integration/precoculture_comparison/Neutrophil_only_UMI100_NC1_precoculture_Harmony_integrated.rds"
analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture"
input_dir <- file.path(analysis_dir, "input")
log_dir <- file.path(analysis_dir, "logs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

output_h5ad <- file.path(input_dir, "all_coculture_neutrophil_fibroblast_scvi_input.h5ad")
hvg_file <- file.path(input_dir, "scvi_training_features.csv")
audit_file <- file.path(input_dir, "input_cell_audit.csv")

message("Loading full co-culture neutrophil/fibroblast object")
obj <- readRDS(coculture_file)
stopifnot(ncol(obj) == 128286L)
stopifnot(all(obj$clean_cell_class %in% c("Neutrophil", "Fibroblast_candidate")))
stopifnot(all(obj$ConditionCode %in% c("N", "NNAMPT", "NA5", "CF", "UF", "UN", "UA5", "UD")))

# scVI's negative-binomial likelihood requires the original integer UMI counts.
# DecontX labels and metadata are retained, but fractional DecontX-corrected counts
# are not supplied to the count likelihood.
raw_counts <- LayerData(obj[["RNA"]], layer = "counts")
stopifnot(inherits(raw_counts, "sparseMatrix"))
if (any(raw_counts@x < 0) || any(abs(raw_counts@x - round(raw_counts@x)) > 1e-8)) {
  stop("RNA count layer is not non-negative integer UMI data")
}

message("Loading the independent pre-culture neutrophil reference")
pre_obj <- readRDS(preculture_file)
pre_cells <- rownames(pre_obj@meta.data)[pre_obj$CultureStage == "Pre-coculture"]
stopifnot(length(pre_cells) == 51285L)
stopifnot(length(intersect(colnames(obj), pre_cells)) == 0L)
pre_counts <- LayerData(pre_obj[["RNA"]], layer = "counts")[, pre_cells, drop = FALSE]
if (any(pre_counts@x < 0) || any(abs(pre_counts@x - round(pre_counts@x)) > 1e-8)) {
  stop("Pre-culture RNA count layer is not non-negative integer UMI data")
}

hvg <- VariableFeatures(obj)
required_markers <- c(
  "FCGR3B", "CSF3R", "S100A8", "S100A9", "CXCR4", "OSM", "PADI4", "MX1",
  "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "SPARC", "VIM",
  "FAP", "PDGFRA", "PDGFRB", "THY1", "ACTA2", "TAGLN", "CXCL12", "IL6"
)
hvg <- unique(c(hvg, intersect(required_markers, rownames(raw_counts))))
hvg <- hvg[hvg %in% intersect(rownames(raw_counts), rownames(pre_counts))]
if (length(hvg) < 2000L) stop("Too few training features")
message("Using ", length(hvg), " HVGs/required markers")

md <- obj@meta.data
md$cell <- rownames(md)
md$CultureStage <- "Post-culture"
md$dataset_source <- "Existing_coCulture"
selected_meta_post <- md[, c(
  "cell", "orig.ident", "SampleID", "DonorID", "Condition", "ConditionCode",
  "CultureStage", "dataset_source",
  "clean_cell_class", "clean_state", "clean_annotation", "stromal_confidence",
  "fibroblast_phenotype", "decontX_contamination", "decontX_removed_fraction_observed"
), drop = FALSE]

pre_md <- pre_obj@meta.data[pre_cells, , drop = FALSE]
selected_meta_pre <- data.frame(
  cell = rownames(pre_md),
  orig.ident = as.character(pre_md$orig.ident),
  SampleID = as.character(pre_md$SampleID),
  DonorID = as.character(pre_md$DonorID),
  Condition = "Pre-coculture",
  ConditionCode = "PRE",
  CultureStage = "Pre-coculture",
  dataset_source = as.character(pre_md$dataset_source),
  clean_cell_class = "Neutrophil",
  clean_state = as.character(pre_md$neutrophil_subset_short),
  clean_annotation = paste("Pre-culture", as.character(pre_md$neutrophil_subset_short)),
  stromal_confidence = NA_real_,
  fibroblast_phenotype = NA_character_,
  decontX_contamination = NA_real_,
  decontX_removed_fraction_observed = NA_real_,
  stringsAsFactors = FALSE,
  row.names = rownames(pre_md)
)
selected_meta <- rbind(selected_meta_post, selected_meta_pre)
selected_meta$clean_cell_class[selected_meta$clean_cell_class == "Fibroblast_candidate"] <- "Fibroblast"
rownames(selected_meta) <- selected_meta$cell

combined_counts <- cbind(
  raw_counts[hvg, colnames(obj), drop = FALSE],
  pre_counts[hvg, pre_cells, drop = FALSE]
)
stopifnot(ncol(combined_counts) == 179571L)
stopifnot(identical(colnames(combined_counts), rownames(selected_meta)))

sce <- SingleCellExperiment(
  assays = list(counts = combined_counts),
  colData = S4Vectors::DataFrame(selected_meta, row.names = colnames(combined_counts))
)
rownames(sce) <- hvg
colnames(sce) <- colnames(combined_counts)

write.csv(data.frame(feature = hvg), hvg_file, row.names = FALSE)
write.csv(
  as.data.frame(table(
    ConditionCode = selected_meta$ConditionCode,
    CellClass = selected_meta$clean_cell_class
  )),
  audit_file,
  row.names = FALSE
)

message("Writing sparse H5AD input with native R HDF5 support")
ad <- anndataR::as_AnnData(sce, assay_name = "counts", output_class = "InMemory")
anndataR::write_h5ad(ad, output_h5ad, compression = "gzip", mode = "w")
message("Export complete: ", output_h5ad)
message("Cells: ", ncol(sce), "; features: ", nrow(sce), "; nonzero counts: ", length(assay(sce, "counts")@x))
