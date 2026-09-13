suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(miloR)
  library(data.table)
})

set.seed(20260816)

pkg <- "FAP_ablation_DSS_epithelial_immune_stromal_figures"
out_dir <- file.path(pkg, "miloR_differential_abundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(file.path(pkg, "tables", "included_sample_units.csv"))
manifest <- manifest[sample_uid %chin% fread("config/fap_biological_mouse_units.csv")$sample_uid]
manifest[, group := fcase(
  condition_label == "FAP-TK PBS", "PBS",
  condition_label == "FAP-TK PBS+DSS", "PBS_DSS",
  condition_label == "FAP-TK GCV+DSS", "GCV_DSS"
)]
manifest[, group := factor(group, levels = c("PBS", "PBS_DSS", "GCV_DSS"))]

object_dir <- file.path(
  "MouseColitis_Inhibitors_compartment_subclustering",
  "reassigned_fractions_epithelial_audit_cleaned", "updated_umaps", "level3_reannotated"
)
object_files <- c(
  Epithelial = "MouseColitis_Inhibitors_Epithelial_Level3Reannotated.rds",
  Immune = "MouseColitis_Inhibitors_Immune_Level3DeepResolved.rds",
  Stromal = "MouseColitis_Inhibitors_Stromal_Level3Reannotated.rds"
)

get_level3 <- function(meta) {
  candidates <- c("level_3", "level3", "reassigned_level_3", "level3_resolved", "cell_type_level_3")
  hit <- candidates[candidates %in% colnames(meta)][1]
  if (is.na(hit)) return(rep("Unannotated", nrow(meta)))
  values <- as.character(meta[[hit]])
  values[is.na(values) | values == ""] <- "Unannotated"
  values
}

run_compartment <- function(compartment, object_file) {
  message("Reading ", compartment, " object")
  object <- readRDS(file.path(object_dir, object_file))
  meta <- as.data.table(object@meta.data, keep.rownames = "cell")
  keep_cells <- intersect(meta[sample_uid %chin% manifest$sample_uid, cell], colnames(object))
  if (!length(keep_cells)) stop("No eligible biological cells found for ", compartment)

  latent <- Embeddings(object, "SCVI_Harmony.reconstructed")[keep_cells, , drop = FALSE]
  latent <- latent[, seq_len(min(30L, ncol(latent))), drop = FALSE]
  umap <- Embeddings(object, "umap.updated.compartment")[keep_cells, seq_len(2), drop = FALSE]
  cell_meta <- meta[match(keep_cells, cell)]
  cell_meta[, level3 := get_level3(cell_meta)]
  cell_meta <- cell_meta[, .(cell, sample_uid, level3)]
  if (compartment == "Immune") {
    subtype_file <- file.path(pkg, "tables", "neutrophil_cell_module_subtypes.csv.gz")
    if (file.exists(subtype_file)) {
      subtype_map <- fread(subtype_file)
      cell_meta[subtype_map, on = "cell", level3 := i.neutrophil_subtype]
    }
  }
  cell_meta <- merge(
    cell_meta,
    manifest[, .(sample_uid, condition_label, group)],
    by = "sample_uid", all.x = TRUE, sort = FALSE, suffixes = c("", ".manifest")
  )
  # merge() can reorder rows; restore exact embedding order.
  cell_meta <- cell_meta[match(keep_cells, cell)]

  # The finalized Seurat objects are large; only their coordinates and the
  # three metadata fields above are required downstream.
  rm(object, meta)
  invisible(gc())

  dummy <- Matrix(0, nrow = 1L, ncol = length(keep_cells), sparse = TRUE)
  rownames(dummy) <- "dummy"
  colnames(dummy) <- keep_cells
  sce <- SingleCellExperiment(assays = list(counts = dummy))
  colData(sce) <- S4Vectors::DataFrame(cell_meta, row.names = keep_cells)
  reducedDim(sce, "SCVI_Harmony") <- latent
  reducedDim(sce, "UMAP") <- umap
  milo <- Milo(sce)

  k_use <- min(30L, max(10L, floor(sqrt(ncol(sce)) / 2)))
  message("Building ", compartment, " graph: ", ncol(sce), " cells, k=", k_use)
  milo <- buildGraph(milo, k = k_use, d = ncol(latent), reduced.dim = "SCVI_Harmony")
  milo <- makeNhoods(
    milo, prop = 0.03, k = k_use, d = ncol(latent), refined = TRUE,
    reduced_dims = "SCVI_Harmony"
  )
  message(compartment, " representative neighbourhoods: ", length(nhoodIndex(milo)))
  milo <- countCells(milo, samples = "sample_uid", meta.data = as.data.frame(colData(milo)))

  present_samples <- colnames(nhoodCounts(milo))
  design_df <- as.data.frame(manifest[sample_uid %chin% present_samples, .(sample_uid, group)])
  design_df <- design_df[match(present_samples, design_df$sample_uid), , drop = FALSE]
  rownames(design_df) <- design_df$sample_uid
  design_df$group <- factor(design_df$group, levels = c("PBS", "PBS_DSS", "GCV_DSS"))
  design <- model.matrix(~ 0 + group, data = design_df)
  rownames(design) <- rownames(design_df)

  contrast_defs <- c(
    DSS_induction = "groupPBS_DSS - groupPBS",
    FAP_ablation = "groupGCV_DSS - groupPBS_DSS"
  )
  centers <- rbindlist(lapply(seq_along(nhoodIndex(milo)), function(i) {
    idx <- which(nhoods(milo)[, i] > 0)
    labels <- as.character(colData(milo)$level3[idx])
    dominant <- names(sort(table(labels), decreasing = TRUE))[1]
    data.table(
      Nhood = i,
      UMAP_1 = median(umap[idx, 1], na.rm = TRUE),
      UMAP_2 = median(umap[idx, 2], na.rm = TRUE),
      nhood_cells = length(idx),
      dominant_level3 = dominant
    )
  }))

  results <- list()
  for (contrast_name in names(contrast_defs)) {
    message("Testing ", compartment, ": ", contrast_name)
    result <- as.data.table(testNhoods(
      milo,
      design = design,
      design.df = design_df,
      model.contrasts = contrast_defs[[contrast_name]],
      reduced.dim = "SCVI_Harmony",
      fdr.weighting = "none",
      robust = FALSE
    ))
    result[, Nhood := seq_len(.N)]
    result <- merge(result, centers, by = "Nhood", all.x = TRUE, sort = FALSE)
    result[, `:=`(
      compartment = compartment,
      contrast = contrast_name,
      numerator = ifelse(contrast_name == "DSS_induction", "FAP-TK PBS+DSS", "FAP-TK GCV+DSS"),
      denominator = ifelse(contrast_name == "DSS_induction", "FAP-TK PBS", "FAP-TK PBS+DSS")
    )]
    fwrite(result, file.path(out_dir, paste0(compartment, "_", contrast_name, "_miloR_results.csv")))
    results[[contrast_name]] <- result
  }

  cell_plot <- data.table(
    cell = keep_cells,
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    sample_uid = as.character(colData(milo)$sample_uid),
    level3 = as.character(colData(milo)$level3),
    compartment = compartment
  )
  # Preserve a light, reproducible background layer for the multipanel plot.
  if (nrow(cell_plot) > 12000L) cell_plot <- cell_plot[sample(.N, 12000L)]
  fwrite(cell_plot, file.path(out_dir, paste0(compartment, "_biological_cells_plot_background.csv.gz")))
  saveRDS(milo, file.path(out_dir, paste0(compartment, "_miloR_object.rds")), compress = FALSE)
  rm(sce, milo)
  invisible(gc())
  rbindlist(results, fill = TRUE)
}

all_results <- rbindlist(Map(run_compartment, names(object_files), unname(object_files)), fill = TRUE)
fwrite(all_results, file.path(out_dir, "Figure4C_all_compartments_miloR_results.csv"))

summary_table <- all_results[, .(
  neighborhoods = .N,
  bh_fdr_lt_0_05 = sum(FDR < 0.05, na.rm = TRUE),
  enriched = sum(FDR < 0.05 & logFC > 0, na.rm = TRUE),
  depleted = sum(FDR < 0.05 & logFC < 0, na.rm = TRUE)
), by = .(compartment, contrast, numerator, denominator)]
fwrite(summary_table, file.path(out_dir, "Figure4C_miloR_summary.csv"))

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(out_dir, "sessionInfo.txt"))
message("MiloR analysis complete: ", out_dir)
