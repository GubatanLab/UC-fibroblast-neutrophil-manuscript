suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(miloR)
  library(data.table)
})

set.seed(20260816)

out_dir <- file.path(
  "a5B1_DSS_epithelial_immune_stromal_remodeling",
  "miloR_differential_abundance"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread("MouseColitis_Inhibitors_manifest.csv")
manifest <- manifest[
  (dataset_source == "Merged_a5B1" & source_condition %chin% c("Control", "DSS+ a5B1_inhibitor")) |
    (sample_id %chin% paste0("D", 1:6) & source_condition == "FAP-TK_PBS_DSS"),
  .(sample_uid, sample_id, source_condition)
]
manifest[, group := fcase(
  source_condition == "Control", "Control",
  sample_id %chin% paste0("D", 1:6), "Severe_DSS",
  default = "Blockade"
)]
manifest[, group := factor(group, levels = c("Control", "Severe_DSS", "Blockade"))]

neutrophil_program_file <- file.path(
  "a5B1_DSS_epithelial_immune_stromal_remodeling",
  "neutrophil_trajectory", "program_subsets",
  "neutrophil_program_assignments.csv.gz"
)
neutrophil_programs <- fread(neutrophil_program_file)[
  sample_uid %chin% manifest$sample_uid,
  .(cell, sample_uid, program)
]
expected_neutrophil_programs <- c(
  "Neutrophil OSM", "Neutrophil CXCR4",
  "Neutrophil PADI4", "Neutrophil MX1"
)
stopifnot(
  !anyDuplicated(neutrophil_programs$cell),
  setequal(unique(neutrophil_programs$program), expected_neutrophil_programs)
)

object_dir <- file.path(
  "MouseColitis_Inhibitors_compartment_subclustering",
  "reassigned_fractions_epithelial_audit_cleaned", "updated_umaps", "level3_reannotated"
)
object_files <- c(
  Epithelial = "MouseColitis_Inhibitors_Epithelial_Level3Reannotated.rds",
  Immune = "MouseColitis_Inhibitors_Immune_Level3DeepResolved.rds",
  Stromal = "MouseColitis_Inhibitors_Stromal_Level3Reannotated.rds"
)

run_compartment <- function(compartment, object_file) {
  message("Reading ", compartment)
  object <- readRDS(file.path(object_dir, object_file))
  meta <- as.data.table(object@meta.data, keep.rownames = "cell")
  eligible_meta <- meta[sample_uid %chin% manifest$sample_uid]
  if (compartment != "Stromal") {
    eligible_meta <- eligible_meta[, .SD[sample(.N, min(.N, 5000L))], by = sample_uid]
  }
  keep_cells <- intersect(eligible_meta$cell, colnames(object))
  if (!length(keep_cells)) stop("No eligible cells for ", compartment)

  latent <- Embeddings(object, "SCVI_Harmony.reconstructed")[keep_cells, , drop = FALSE]
  latent <- latent[, seq_len(min(30L, ncol(latent))), drop = FALSE]
  cell_meta <- meta[match(keep_cells, cell), .(
    cell, sample_uid, level3 = as.character(reassigned_level_3)
  )]
  if (compartment == "Immune") {
    original_neutrophils <- grepl("neutrophil", cell_meta$level3, ignore.case = TRUE)
    program_index <- match(cell_meta$cell, neutrophil_programs$cell)
    if (any(original_neutrophils & is.na(program_index))) {
      stop("Every retained immune neutrophil must have a validated OSM/CXCR4/PADI4/MX1 call")
    }
    if (any(!original_neutrophils & !is.na(program_index))) {
      stop("Validated neutrophil program calls matched non-neutrophil Level 3 cells")
    }
    cell_meta[!is.na(program_index), level3 := neutrophil_programs$program[program_index[!is.na(program_index)]]]
  }
  cell_meta <- merge(
    cell_meta, manifest[, .(sample_uid, group)], by = "sample_uid",
    all.x = TRUE, sort = FALSE
  )[match(keep_cells, cell)]
  rm(object, meta)
  invisible(gc())

  dummy <- Matrix(0, nrow = 1L, ncol = length(keep_cells), sparse = TRUE)
  rownames(dummy) <- "dummy"
  colnames(dummy) <- keep_cells
  sce <- SingleCellExperiment(assays = list(counts = dummy))
  colData(sce) <- S4Vectors::DataFrame(cell_meta, row.names = keep_cells)
  reducedDim(sce, "SCVI_Harmony") <- latent
  milo <- Milo(sce)

  k_use <- min(30L, max(10L, floor(sqrt(ncol(sce)) / 2)))
  milo <- buildGraph(milo, k = k_use, d = ncol(latent), reduced.dim = "SCVI_Harmony")
  milo <- makeNhoods(
    milo, prop = 0.03, k = k_use, d = ncol(latent), refined = TRUE,
    reduced_dims = "SCVI_Harmony"
  )
  milo <- countCells(milo, samples = "sample_uid", meta.data = as.data.frame(colData(milo)))

  present_samples <- colnames(nhoodCounts(milo))
  design_df <- as.data.frame(manifest[sample_uid %chin% present_samples, .(sample_uid, group)])
  design_df <- design_df[match(present_samples, design_df$sample_uid), , drop = FALSE]
  rownames(design_df) <- design_df$sample_uid
  design_df$group <- factor(design_df$group, levels = c("Control", "Severe_DSS", "Blockade"))
  design <- model.matrix(~ 0 + group, data = design_df)
  rownames(design) <- rownames(design_df)

  centers <- rbindlist(lapply(seq_along(nhoodIndex(milo)), function(i) {
    idx <- which(nhoods(milo)[, i] > 0)
    labels <- as.character(colData(milo)$level3[idx])
    data.table(
      Nhood = i,
      nhood_cells = length(idx),
      dominant_level3 = names(sort(table(labels), decreasing = TRUE))[1]
    )
  }))

  contrast_defs <- c(
    `Severe DSS colitis - Control` = "groupSevere_DSS - groupControl",
    `Blockade - severe DSS colitis` = "groupBlockade - groupSevere_DSS"
  )
  result_list <- lapply(names(contrast_defs), function(contrast_name) {
    message("Testing ", compartment, ": ", contrast_name)
    result <- as.data.table(testNhoods(
      milo, design = design, design.df = design_df,
      model.contrasts = contrast_defs[[contrast_name]],
      reduced.dim = "SCVI_Harmony", fdr.weighting = "none",
      robust = FALSE
    ))
    result[, Nhood := seq_len(.N)]
    result <- merge(result, centers, by = "Nhood", all.x = TRUE, sort = FALSE)
    result[, `:=`(compartment = compartment, contrast = contrast_name)]
    result
  })
  rm(sce, milo)
  invisible(gc())
  rbindlist(result_list, fill = TRUE)
}

all_results <- rbindlist(
  Map(run_compartment, names(object_files), unname(object_files)), fill = TRUE
)
fwrite(all_results, file.path(out_dir, "Figure5C_all_compartments_miloR_results.csv"))

summary_table <- all_results[, .(
  neighborhoods = .N,
  fdr_lt_0_05 = sum(FDR < 0.05, na.rm = TRUE),
  enriched = sum(FDR < 0.05 & logFC > 0, na.rm = TRUE),
  depleted = sum(FDR < 0.05 & logFC < 0, na.rm = TRUE)
), by = .(compartment, contrast)]
fwrite(summary_table, file.path(out_dir, "Figure5C_miloR_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Figure 5 MiloR analysis complete")
