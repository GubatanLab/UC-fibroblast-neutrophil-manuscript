#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(miloR)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(BiocNeighbors)
  library(BiocParallel)
})

set.seed(20260816)

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture/Neutrophil_Subclustering"
object_dir <- file.path(analysis_dir, "objects")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

contrast_key <- tibble::tribble(
  ~reference, ~comparison, ~contrast, ~paired,
  "PRE", "N", "Neutrophils alone vs pre-culture", FALSE,
  "N", "CF", "Control fibroblast vs neutrophils alone", TRUE,
  "N", "UF", "UC fibroblast vs neutrophils alone", TRUE,
  "CF", "UF", "UC vs control fibroblast", TRUE,
  "N", "NNAMPT", "NAMPTi vs neutrophils alone", TRUE,
  "N", "NA5", "alpha5beta1i vs neutrophils alone", TRUE,
  "UF", "UN", "NAMPTi vs UC", TRUE,
  "UF", "UA5", "alpha5beta1i vs UC", TRUE,
  "UF", "UD", "Dual blockade vs UC", TRUE
)
fwrite(contrast_key, file.path(source_dir, "miloR_neutrophil_DA_contrasts.csv"))

make_minimal_milo <- function(seurat_object, label_col, panel_name) {
  md <- seurat_object@meta.data
  md$DonorID <- as.character(md$DonorID)
  md$ConditionCode <- as.character(md$ConditionCode)
  md$sample_id <- paste(md$DonorID, md$ConditionCode, sep = "__")
  md[[label_col]] <- factor(as.character(md[[label_col]]), levels = unique(as.character(md[[label_col]])))

  placeholder <- Matrix::Matrix(0, nrow = 1, ncol = ncol(seurat_object), sparse = TRUE)
  rownames(placeholder) <- "placeholder"
  colnames(placeholder) <- colnames(seurat_object)
  sce <- SingleCellExperiment(
    assays = list(logcounts = placeholder),
    colData = S4Vectors::DataFrame(md)
  )
  reducedDim(sce, "scvi.harmony") <- Embeddings(seurat_object, "scvi.harmony")[colnames(sce), 1:30, drop = FALSE]
  milo <- Milo(sce)
  message(panel_name, ": building 30-nearest-neighbor graph")
  milo <- buildGraph(
    milo, k = 30, d = 30, reduced.dim = "scvi.harmony",
    BNPARAM = BiocNeighbors::AnnoyParam(ntrees = 50),
    BPPARAM = BiocParallel::SerialParam()
  )
  message(panel_name, ": sampling and refining Milo neighborhoods")
  milo <- makeNhoods(
    milo, prop = 0.10, k = 30, d = 30, refined = TRUE,
    reduced_dims = "scvi.harmony", refinement_scheme = "reduced_dim"
  )
  milo <- countCells(milo, samples = "sample_id", meta.data = as.data.frame(colData(milo)))
  milo <- calcNhoodDistance(milo, d = 30, reduced.dim = "scvi.harmony")
  list(
    milo = milo,
    sample_design = md %>% distinct(sample_id, DonorID, ConditionCode),
    label_col = label_col,
    panel_name = panel_name
  )
}

annotate_neighborhoods <- function(milo, da_res, neighborhood_indices, label_col) {
  labels <- as.character(colData(milo)[[label_col]])
  nh_mat <- nhoods(milo)
  annotation <- lapply(neighborhood_indices, function(idx) {
    members <- which(nh_mat[, idx] != 0)
    tab <- sort(table(labels[members]), decreasing = TRUE)
    data.frame(
      Nhood = idx,
      state = names(tab)[1],
      state_fraction = as.numeric(tab[1]) / sum(tab),
      annotated_nhood_size = length(members),
      stringsAsFactors = FALSE
    )
  })
  annotation <- bind_rows(annotation)
  stopifnot(nrow(annotation) == nrow(da_res))
  bind_cols(as.data.frame(da_res), annotation)
}

run_contrast <- function(milo_bundle, reference, comparison, contrast, paired) {
  milo <- milo_bundle$milo
  sample_design <- milo_bundle$sample_design
  keep_design <- sample_design %>%
    filter(ConditionCode %in% c(reference, comparison)) %>%
    mutate(
      ConditionCode = relevel(factor(ConditionCode), ref = reference),
      DonorID = factor(DonorID)
    )
  keep_samples <- keep_design$sample_id
  keep_samples <- keep_samples[keep_samples %in% colnames(nhoodCounts(milo))]
  keep_design <- keep_design[match(keep_samples, keep_design$sample_id), , drop = FALSE]
  rownames(keep_design) <- keep_design$sample_id

  contrast_milo <- milo
  nhoodCounts(contrast_milo) <- nhoodCounts(contrast_milo)[, keep_samples, drop = FALSE]
  keep_nhoods <- Matrix::rowSums(nhoodCounts(contrast_milo)) > 0
  formula_used <- if (paired) ~ DonorID + ConditionCode else ~ ConditionCode
  message(
    milo_bundle$panel_name, ": ", contrast,
    " (", sum(keep_nhoods), " neighborhoods; ", length(keep_samples), " samples)"
  )
  da_res <- testNhoods(
    contrast_milo,
    design = formula_used,
    design.df = as.data.frame(keep_design),
    fdr.weighting = "k-distance",
    robust = TRUE,
    reduced.dim = "scvi.harmony",
    subset.nhoods = keep_nhoods,
    BPPARAM = BiocParallel::SerialParam()
  )
  neighborhood_indices <- which(keep_nhoods)
  da_res <- annotate_neighborhoods(
    milo_bundle$milo, da_res, neighborhood_indices, milo_bundle$label_col
  )
  da_res$panel <- milo_bundle$panel_name
  da_res$reference <- reference
  da_res$comparison <- comparison
  da_res$contrast <- contrast
  da_res$paired <- paired
  da_res$n_samples <- length(keep_samples)
  da_res$n_pairs <- if (paired) length(unique(keep_design$DonorID)) else NA_integer_
  da_res
}

run_panel <- function(seurat_object, label_col, panel_name) {
  bundle <- make_minimal_milo(seurat_object, label_col, panel_name)
  results <- bind_rows(lapply(seq_len(nrow(contrast_key)), function(i) {
    run_contrast(
      bundle,
      contrast_key$reference[i], contrast_key$comparison[i],
      contrast_key$contrast[i], contrast_key$paired[i]
    )
  }))
  list(bundle = bundle, results = results)
}

figure_a_results_file <- file.path(source_dir, "miloR_FigureA_subcluster_neighborhood_DA_results.csv.gz")
figure_b_results_file <- file.path(source_dir, "miloR_FigureB_canonical_state_neighborhood_DA_results.csv.gz")
figure_a_milo_file <- file.path(object_dir, "miloR_FigureA_subcluster_neighborhoods.rds")
figure_b_milo_file <- file.path(object_dir, "miloR_FigureB_canonical_state_neighborhoods.rds")

resume_from_checkpoints <- all(file.exists(c(
  figure_a_results_file, figure_b_results_file,
  figure_a_milo_file, figure_b_milo_file
)))

if (resume_from_checkpoints) {
  message("Loading completed MiloR neighborhood results from checkpoints")
  figure_a <- list(
    bundle = NULL,
    results = fread(figure_a_results_file)
  )
  figure_b <- list(
    bundle = NULL,
    results = fread(figure_b_results_file)
  )
} else {
  message("Loading Figure A and Figure B neutrophil objects")
  figure_a_object <- readRDS(file.path(object_dir, "neutrophils_scVI_Harmony_subclustered_annotated_clean.rds"))
  figure_b_object <- readRDS(file.path(object_dir, "canonical_neutrophils_scVI_Harmony_reclustered.rds"))

  figure_a <- run_panel(
    figure_a_object,
    label_col = "neutrophil_cluster_annotation",
    panel_name = "Figure A subcluster annotations"
  )
  gc()
  figure_b <- run_panel(
    figure_b_object,
    label_col = "neutrophil_state_annotation",
    panel_name = "Figure B canonical states"
  )
}

fwrite(
  figure_a$results,
  file.path(source_dir, "miloR_FigureA_subcluster_neighborhood_DA_results.csv.gz")
)
fwrite(
  figure_b$results,
  file.path(source_dir, "miloR_FigureB_canonical_state_neighborhood_DA_results.csv.gz")
)

if (!resume_from_checkpoints) {
  saveRDS(
    figure_a$bundle$milo,
    figure_a_milo_file,
    compress = FALSE
  )
  saveRDS(
    figure_b$bundle$milo,
    figure_b_milo_file,
    compress = FALSE
  )
}

summarize_milo <- function(results) {
  results %>%
    filter(state_fraction >= 0.50) %>%
    group_by(panel, contrast, reference, comparison, paired, state) %>%
    summarise(
      neighborhoods = n(),
      median_logFC = median(logFC, na.rm = TRUE),
      mean_logFC = mean(logFC, na.rm = TRUE),
      weighted_mean_logFC = weighted.mean(logFC, w = pmax(annotated_nhood_size, 1), na.rm = TRUE),
      median_state_fraction = median(state_fraction, na.rm = TRUE),
      da_neighborhoods_spatial_FDR_0_05 = sum(SpatialFDR < 0.05, na.rm = TRUE),
      da_up_spatial_FDR_0_05 = sum(SpatialFDR < 0.05 & logFC > 0, na.rm = TRUE),
      da_down_spatial_FDR_0_05 = sum(SpatialFDR < 0.05 & logFC < 0, na.rm = TRUE),
      da_fraction_spatial_FDR_0_05 = mean(SpatialFDR < 0.05, na.rm = TRUE),
      minimum_SpatialFDR = min(SpatialFDR, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_results <- bind_rows(
  summarize_milo(figure_a$results),
  summarize_milo(figure_b$results)
)
fwrite(summary_results, file.path(source_dir, "miloR_neutrophil_state_DA_heatmap_summary.csv"))

figure_a_order <- c(
  "Low-signal neutrophil", "Activated OSM/CXCR4 neutrophil",
  "IL1R2+ inflammatory neutrophil", "DHFR+ proliferative/stress neutrophil",
  "PADI4/translation-high neutrophil", "LTF/BPI immature neutrophil",
  "CCL3/CCL4 inflammatory neutrophil", "RORA+ atypical neutrophil"
)
figure_b_order <- c(
  "PADI4 neutrophil", "OSM neutrophil",
  "CXCR4 neutrophil", "MX1/ISG neutrophil"
)
contrast_order <- contrast_key$contrast

make_heatmap_df <- function(panel_name, state_order) {
  summary_results %>%
    filter(panel == panel_name) %>%
    complete(
      contrast = contrast_order,
      state = state_order,
      fill = list(
        neighborhoods = 0,
        da_neighborhoods_spatial_FDR_0_05 = 0,
        da_up_spatial_FDR_0_05 = 0,
        da_down_spatial_FDR_0_05 = 0
      )
    ) %>%
    mutate(
      contrast = factor(contrast, levels = contrast_order),
      state = factor(state, levels = rev(state_order)),
      tile_label = if_else(
        neighborhoods > 0,
        paste0(
          sprintf("%+.2f", median_logFC), "\n",
          da_neighborhoods_spatial_FDR_0_05, "/", neighborhoods
        ),
        "NA"
      )
    )
}

heat_a <- make_heatmap_df("Figure A subcluster annotations", figure_a_order)
heat_b <- make_heatmap_df("Figure B canonical states", figure_b_order)
all_effects <- c(heat_a$median_logFC, heat_b$median_logFC)
effect_limit <- as.numeric(quantile(abs(all_effects[is.finite(all_effects)]), 0.95, na.rm = TRUE))
effect_limit <- max(effect_limit, 0.25)

plot_heatmap <- function(dat, title, subtitle, show_x = TRUE) {
  p <- ggplot(dat, aes(contrast, state, fill = median_logFC)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = tile_label), size = 2.8, lineheight = 0.9) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
      limits = c(-effect_limit, effect_limit), oob = scales::squish,
      na.value = "#E5E5E5"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL, y = NULL,
      fill = "Median MiloR\nlog2 fold change"
    ) +
    theme_bw(base_family = "Arial", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.5),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank()
    )
  if (show_x) {
    p <- p + theme(axis.text.x = element_text(angle = 42, hjust = 1, size = 8))
  } else {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  p
}

p_a <- plot_heatmap(
  heat_a,
  "Figure A: MiloR differential abundance of neutrophil subcluster annotations",
  "Tile text: median neighborhood log2FC and Spatial FDR <0.05 neighborhoods / tested neighborhoods",
  show_x = FALSE
)
p_b <- plot_heatmap(
  heat_b,
  "Figure B: MiloR differential abundance of canonical neutrophil states",
  "Positive log2FC indicates enrichment in the comparison condition; negative indicates depletion",
  show_x = TRUE
)
p_all <- p_a / p_b +
  plot_layout(heights = c(1.35, 0.95), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption = "MiloR neighborhoods: k = 30, d = 30, 10% refined index sampling. Culture contrasts include donor blocking; pre-culture versus neutrophils alone is unpaired because donor sets differ."
  ) &
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14),
    plot.caption = element_text(size = 8.5, color = "#666666", hjust = 0)
  )

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}

save_plot(
  p_all,
  "Neutrophil_7_MiloR_DA_heatmap_FigureA_subclusters_and_FigureB_canonical_states",
  16, 12
)
message("MiloR Figure A/Figure B differential-abundance analysis complete")
