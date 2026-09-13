suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(SingleCellExperiment)
  library(slingshot)
  library(uwot)
})

set.seed(20260815)

out_dir <- file.path(
  "a5B1_DSS_epithelial_immune_stromal_remodeling",
  "neutrophil_trajectory"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

immune_file <- file.path(
  "MouseColitis_Inhibitors_compartment_subclustering",
  "reassigned_fractions_epithelial_audit_cleaned", "updated_umaps", "level3_reannotated",
  "MouseColitis_Inhibitors_Immune_Level3DeepResolved.rds"
)

condition_levels <- c("Control", "Severe DSS colitis", "DSS + alpha5beta1 blockade")
condition_colors <- c("Control" = "#4C78A8", "Severe DSS colitis" = "#D95F02", "DSS + alpha5beta1 blockade" = "#1B9E77")
state_levels <- c("Cxcr2+ neutrophil", "Inflammatory neutrophil")
state_colors <- c("Cxcr2+ neutrophil" = "#2166AC", "Inflammatory neutrophil" = "#B2182B")

manifest <- fread("MouseColitis_Inhibitors_manifest.csv")
included <- manifest[
  (dataset_source == "Merged_a5B1" & source_condition %chin% c("Control", "DSS+ a5B1_inhibitor")) |
    (sample_id %chin% c("C7", "C8", "C9", "C10") & source_condition == "Control") |
    (sample_id %chin% paste0("D", 1:6) & source_condition == "FAP-TK_PBS_DSS"),
  .(sample_uid, sample_id, source_condition, total_cells = as.numeric(cells))
]
included[, condition := factor(
  fifelse(source_condition == "Control", "Control",
          fifelse(sample_id %chin% paste0("D", 1:6), "Severe DSS colitis", "DSS + alpha5beta1 blockade")),
  levels = condition_levels
)]

message("Loading curated immune object")
immune <- readRDS(immune_file)
meta <- as.data.table(immune@meta.data, keep.rownames = "cell")
neut_meta <- meta[
  sample_uid %chin% included$sample_uid & as.character(reassigned_level_2) == "Neutrophil" &
    as.character(reassigned_level_3) %chin% state_levels,
  .(
    cell, sample_uid, sample_id, source_condition,
    state = factor(as.character(reassigned_level_3), levels = state_levels)
  )
]
neut_meta[, condition := factor(
  fifelse(source_condition == "Control", "Control",
          fifelse(sample_id %chin% paste0("D", 1:6), "Severe DSS colitis", "DSS + alpha5beta1 blockade")),
  levels = condition_levels
)]
if (nrow(neut_meta) < 100L || uniqueN(neut_meta$state) < 2L) {
  stop("Insufficient curated neutrophils or endpoint states for trajectory inference")
}

latent_name <- "SCVI_Harmony.reconstructed"
latent <- Embeddings(immune, reduction = latent_name)[neut_meta$cell, , drop = FALSE]
latent <- latent[, seq_len(min(30L, ncol(latent))), drop = FALSE]

message("Computing neutrophil-specific UMAP and Slingshot continuum")
neut_umap <- uwot::umap(
  latent, n_neighbors = min(30L, nrow(latent) - 1L), min_dist = 0.25,
  metric = "cosine", n_components = 2L, n_threads = 2L, verbose = FALSE,
  ret_model = FALSE, init = "spectral"
)
rownames(neut_umap) <- rownames(latent)

sce <- SingleCellExperiment(assays = list(counts = matrix(0, nrow = 1, ncol = nrow(neut_umap))))
colnames(sce) <- rownames(neut_umap)
reducedDims(sce) <- SimpleList(UMAP = neut_umap)
colData(sce)$state <- neut_meta$state[match(colnames(sce), neut_meta$cell)]
sce <- slingshot(
  sce, clusterLabels = "state", reducedDim = "UMAP",
  start.clus = "Cxcr2+ neutrophil", end.clus = "Inflammatory neutrophil",
  allow.breaks = FALSE, stretch = 0
)
pseudotime <- slingPseudotime(sce)[, 1]
pseudotime <- (pseudotime - min(pseudotime, na.rm = TRUE)) /
  (max(pseudotime, na.rm = TRUE) - min(pseudotime, na.rm = TRUE))
curve <- slingCurves(sce)[[1]]
curve_data <- as.data.table(curve$s[curve$ord, , drop = FALSE])
setnames(curve_data, c("UMAP_1", "UMAP_2"))
fwrite(curve_data, file.path(out_dir, "neutrophil_slingshot_curve.csv"))

cell_data <- copy(neut_meta)
cell_data[, `:=`(
  UMAP_1 = neut_umap[cell, 1], UMAP_2 = neut_umap[cell, 2],
  pseudotime = pseudotime[cell]
)]
cell_data <- cell_data[is.finite(pseudotime)]
fwrite(cell_data, file.path(out_dir, "neutrophil_cell_trajectory.csv.gz"))

# Per-sample-unit summaries preserve the sample unit as the inferential unit.
sample_summary <- cell_data[, .(
  neutrophil_cells = .N,
  median_pseudotime = median(pseudotime),
  mean_pseudotime = mean(pseudotime),
  cxcr2_fraction = mean(state == "Cxcr2+ neutrophil"),
  inflammatory_fraction = mean(state == "Inflammatory neutrophil")
), by = .(sample_uid, sample_id, condition)]
sample_summary <- merge(included[, .(sample_uid, total_cells)], sample_summary, by = "sample_uid", all.y = TRUE)
sample_summary[, total_neutrophil_fraction := neutrophil_cells / total_cells]
fwrite(sample_summary, file.path(out_dir, "neutrophil_trajectory_by_mouse.csv"))
fwrite(sample_summary, file.path(out_dir, "neutrophil_trajectory_by_sample_unit.csv"))

trajectory_contrasts <- data.table(
  contrast = c("Severe DSS colitis - Control", "Blockade - severe DSS colitis"),
  reference = c("Control", "Severe DSS colitis"), treatment = c("Severe DSS colitis", "DSS + alpha5beta1 blockade"),
  x = c(1.5, 2.5)
)
trajectory_test <- rbindlist(lapply(seq_len(nrow(trajectory_contrasts)), function(i) {
  z <- trajectory_contrasts[i]
  ref <- sample_summary[condition == z$reference, median_pseudotime]
  trt <- sample_summary[condition == z$treatment, median_pseudotime]
  wt <- wilcox.test(trt, ref, exact = FALSE)
  data.table(
    contrast = z$contrast, reference = z$reference, treatment = z$treatment,
    reference_sample_units = length(ref), treatment_sample_units = length(trt),
    reference_mean = mean(ref), treatment_mean = mean(trt),
    difference = mean(trt) - mean(ref), p_value = wt$p.value, x = z$x
  )
}))
trajectory_test[, FDR := p.adjust(p_value, method = "BH")]
trajectory_test[, significance := fifelse(FDR < 0.001, "***", fifelse(FDR < 0.01, "**",
                                      fifelse(FDR < 0.05, "*", "ns")))]
fwrite(trajectory_test, file.path(out_dir, "neutrophil_median_pseudotime_test.csv"))

# Six fixed pseudotime bins; occupancy is calculated within each sample unit before group averaging.
cell_data[, pseudotime_bin := cut(
  pseudotime, breaks = seq(0, 1, length.out = 7), include.lowest = TRUE,
  labels = paste0("B", 1:6)
)]
occupancy <- cell_data[, .(cells = .N), by = .(sample_uid, sample_id, condition, pseudotime_bin)]
grid <- CJ(sample_uid = sample_summary$sample_uid, pseudotime_bin = levels(cell_data$pseudotime_bin), unique = TRUE)
grid <- merge(grid, sample_summary[, .(sample_uid, sample_id, condition, neutrophil_cells)], by = "sample_uid")
occupancy <- merge(grid, occupancy, by = c("sample_uid", "sample_id", "condition", "pseudotime_bin"), all.x = TRUE)
occupancy[is.na(cells), cells := 0L]
occupancy[, fraction := cells / neutrophil_cells]
occupancy[, bin_midpoint := (as.integer(factor(pseudotime_bin, levels = paste0("B", 1:6))) - 0.5) / 6]
fwrite(occupancy, file.path(out_dir, "neutrophil_pseudotime_bin_occupancy_by_mouse.csv"))
fwrite(occupancy, file.path(out_dir, "neutrophil_pseudotime_bin_occupancy_by_sample_unit.csv"))

occupancy_summary <- occupancy[, .(
  mean = mean(fraction), sem = sd(fraction) / sqrt(.N), mice = .N
), by = .(condition, pseudotime_bin, bin_midpoint)]

occupancy_tests <- rbindlist(lapply(seq_len(nrow(trajectory_contrasts)), function(i) {
  z <- trajectory_contrasts[i]
  occupancy[, {
    ref <- fraction[condition == z$reference]
    trt <- fraction[condition == z$treatment]
    wt <- wilcox.test(trt, ref, exact = FALSE)
    .(contrast = z$contrast, reference_mean = mean(ref), treatment_mean = mean(trt), p_value = wt$p.value)
  }, by = pseudotime_bin]
}))
occupancy_tests[, FDR := p.adjust(p_value, method = "BH"), by = contrast]
fwrite(occupancy_tests, file.path(out_dir, "neutrophil_pseudotime_bin_tests.csv"))

# Mouse-weighted expression dynamics for interpretable neutrophil markers.
candidate_genes <- c("Cxcr2", "S100a8", "S100a9", "Il1b", "Cxcl2", "Osm")
genes <- intersect(candidate_genes, rownames(immune))
rna_layers <- Layers(immune[["RNA"]])
if ("data" %in% rna_layers) {
  expr <- GetAssayData(immune, assay = "RNA", layer = "data")[genes, cell_data$cell, drop = FALSE]
} else {
  # The merged Seurat v5 object stores one count layer per source. Join only the
  # trajectory cells to avoid materializing a full-object joined matrix.
  neut_expr_object <- subset(immune, cells = cell_data$cell)
  neut_expr_object <- JoinLayers(neut_expr_object, assay = "RNA")
  joined_counts <- GetAssayData(neut_expr_object, assay = "RNA", layer = "counts")
  counts <- joined_counts[genes, cell_data$cell, drop = FALSE]
  totals <- Matrix::colSums(joined_counts[, cell_data$cell, drop = FALSE])
  expr <- log1p(t(t(as.matrix(counts)) / pmax(totals, 1) * 1e4))
  rm(neut_expr_object, joined_counts)
}
expr <- as.matrix(expr)
rownames(expr) <- genes
colnames(expr) <- cell_data$cell
expr_z <- t(scale(t(expr)))
rownames(expr_z) <- genes
colnames(expr_z) <- cell_data$cell
expr_z[!is.finite(expr_z)] <- 0
expr_z[expr_z > 2.5] <- 2.5
expr_z[expr_z < -2.5] <- -2.5
rownames(expr_z) <- genes
colnames(expr_z) <- cell_data$cell
expr_t <- t(expr_z)
expr_long <- as.data.table(expr_t)
expr_long[, cell := rownames(expr_t)]
setcolorder(expr_long, "cell")
expr_long <- melt(expr_long, id.vars = "cell", variable.name = "gene", value.name = "z_expression")
expr_long <- merge(
  expr_long,
  cell_data[, .(cell, sample_uid, condition, pseudotime_bin, bin_midpoint = (as.integer(pseudotime_bin) - 0.5) / 6)],
  by = "cell"
)
expr_mouse_bin <- expr_long[, .(z_expression = mean(z_expression)),
                            by = .(gene, sample_uid, condition, pseudotime_bin, bin_midpoint)]
expr_summary <- expr_mouse_bin[, .(
  mean = mean(z_expression), sem = sd(z_expression) / sqrt(.N), mice = .N
), by = .(gene, condition, pseudotime_bin, bin_midpoint)]
fwrite(expr_summary, file.path(out_dir, "neutrophil_marker_dynamics.csv"))

theme_manuscript <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 7.5, color = "grey25"),
    legend.title = element_blank(), legend.text = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 7.5),
    panel.grid.minor = element_blank()
  )

p_state <- ggplot(cell_data, aes(UMAP_1, UMAP_2, color = state)) +
  geom_point(size = 0.55, alpha = 0.75) +
  geom_path(data = curve_data, aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
            color = "black", linewidth = 0.8,
            arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  scale_color_manual(values = state_colors) +
  labs(title = "Inferred neutrophil state continuum",
       subtitle = "Rooted at Cxcr2+ and terminating in the inflammatory state") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 1))

p_condition <- ggplot(cell_data, aes(UMAP_1, UMAP_2, color = condition)) +
  geom_point(size = 0.55, alpha = 0.75) +
  geom_path(data = curve_data, aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
            color = "black", linewidth = 0.7,
            arrow = arrow(length = unit(0.10, "inches"), type = "closed")) +
  facet_wrap(~condition, nrow = 1) +
  scale_color_manual(values = condition_colors) +
  labs(title = "Condition occupancy along the shared continuum",
       subtitle = "Identical trajectory and UMAP axes across conditions") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5), strip.text = element_text(face = "bold", size = 7),
        legend.position = "none")

p_occupancy <- ggplot(occupancy_summary, aes(bin_midpoint, mean, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = pmax(0, mean - sem), ymax = mean + sem), alpha = 0.16, linewidth = 0) +
  geom_line(linewidth = 0.85) + geom_point(size = 1.8) +
  scale_color_manual(values = condition_colors) + scale_fill_manual(values = condition_colors) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("Cxcr2+", "Intermediate", "Inflammatory"),
                     limits = c(0, 1)) +
  scale_y_continuous(labels = percent) +
  labs(title = "Sample-unit-weighted trajectory occupancy", subtitle = "Mean +/- SEM across sample units",
       x = "Inferred pseudotime", y = "Fraction of neutrophils") +
  theme_manuscript + theme(legend.position = "bottom")

ann_y <- min(1.0, max(sample_summary$median_pseudotime) * 1.03)
trajectory_test[, y := ann_y + c(0, 0.07)]
p_median <- ggplot(sample_summary, aes(condition, median_pseudotime, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.07), size = 1.8, alpha = 0.9) +
  geom_text(data = trajectory_test, aes(x = x, y = y, label = significance),
            inherit.aes = FALSE, fontface = "bold", size = 3.5) +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(limits = c(0, 1.12)) +
  labs(title = "Median pseudotime per sample unit", subtitle = "Control-DSS and DSS-blockade Wilcoxon tests",
       x = NULL, y = "Median inferred pseudotime") +
  theme_manuscript +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1, size = 7))

state_long <- melt(
  sample_summary[, .(sample_uid, condition, `Cxcr2+ neutrophil` = cxcr2_fraction,
                     `Inflammatory neutrophil` = inflammatory_fraction)],
  id.vars = c("sample_uid", "condition"), variable.name = "state", value.name = "fraction"
)
p_states <- ggplot(state_long, aes(condition, fraction, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.45) +
  geom_point(position = position_jitter(width = 0.07), size = 1.5, alpha = 0.9) +
  facet_wrap(~state, nrow = 1) +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(labels = percent) +
  labs(title = "Endpoint-state redistribution", subtitle = "Fractions within total neutrophils",
       x = NULL, y = "Neutrophil-state fraction") +
  theme_manuscript +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_markers <- ggplot(expr_summary, aes(bin_midpoint, mean, color = condition, group = condition)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = condition),
              alpha = 0.13, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  scale_color_manual(values = condition_colors) + scale_fill_manual(values = condition_colors) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("Cxcr2+", "", "Inflammatory"), limits = c(0, 1)) +
  labs(title = "Marker dynamics along inferred pseudotime",
       subtitle = "Gene-wise standardized expression; sample-unit-weighted bin means",
       x = "Inferred pseudotime", y = "Standardized expression") +
  theme_manuscript +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6.5))

layout <- c(
  area(t = 1, l = 1, b = 1, r = 2),
  area(t = 1, l = 3, b = 1, r = 6),
  area(t = 2, l = 1, b = 2, r = 3),
  area(t = 2, l = 4, b = 2, r = 6),
  area(t = 3, l = 1, b = 3, r = 3),
  area(t = 3, l = 4, b = 3, r = 6)
)

trajectory_figure <- p_state + p_condition + p_occupancy + p_median + p_states + p_markers +
  plot_layout(design = layout, heights = c(1.5, 1.15, 1.35)) +
  plot_annotation(
    title = "DSS and alpha5beta1 blockade redistribute neutrophils along an inferred activation continuum",
    subtitle = paste0(
      "Curated neutrophils: ", nrow(cell_data), " cells from ", uniqueN(sample_summary$sample_uid),
      " sample units; arrows indicate orientation, not observed temporal fate"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9, color = "grey25"),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

base <- file.path(out_dir, "Figure_5_supplement_neutrophil_trajectory")
ggsave(paste0(base, ".png"), trajectory_figure, width = 15, height = 13, dpi = 300,
       bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".pdf"), trajectory_figure, width = 15, height = 13,
       bg = "white", device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(base, ".tiff"), trajectory_figure, width = 15, height = 13, dpi = 300,
       bg = "white", compression = "lzw", limitsize = FALSE)

legend <- c(
  "# Figure 5—figure supplement. DSS and alpha5beta1 blockade redistribute neutrophils along an inferred activation continuum.",
  "",
  "(A) Neutrophil-specific UMAP and Slingshot continuum rooted at the curated Cxcr2+ neutrophil state and terminating in the inflammatory neutrophil state. The arrow indicates the imposed biological orientation of the inferred continuum and does not demonstrate observed temporal fate.",
  "(B) Condition-specific occupancy of the shared neutrophil UMAP and trajectory.",
  "(C) Sample-unit-weighted occupancy across six fixed pseudotime bins. Lines and ribbons show the mean and SEM across all included units.",
  "(D) Median inferred pseudotime for each included sample unit.",
  "(E) Sample-unit-level fractions of the two curated endpoint states within total neutrophils.",
  "(F) Sample-unit-weighted standardized expression dynamics for selected neutrophil markers across the inferred continuum.",
  "",
  "Trajectory inference used the first 30 dimensions of the curated SCVI/Harmony representation, a neutrophil-specific UMAP, and Slingshot with Cxcr2+ and inflammatory neutrophils specified as start and terminal clusters. Statistical comparisons use sample units, not individual cells. The trajectory is a cross-sectional state continuum and does not establish lineage direction, transition rates, or cellular fate. Treatment remains aligned with source batch, so all trajectory differences are exploratory."
)
writeLines(legend, file.path(out_dir, "Figure_5_supplement_neutrophil_trajectory_legend.md"))

validation <- data.table(
  cells = nrow(cell_data), sample_units = uniqueN(sample_summary$sample_uid),
  control_sample_units = sample_summary[condition == condition_levels[1], uniqueN(sample_uid)],
  dss_sample_units = sample_summary[condition == condition_levels[2], uniqueN(sample_uid)],
  blockade_sample_units = sample_summary[condition == condition_levels[3], uniqueN(sample_uid)],
  states = uniqueN(cell_data$state), pseudotime_complete = all(is.finite(cell_data$pseudotime)),
  png_exists = file.exists(paste0(base, ".png")), pdf_exists = file.exists(paste0(base, ".pdf")),
  tiff_exists = file.exists(paste0(base, ".tiff"))
)
validation[, passed := states == 2L & pseudotime_complete & png_exists & pdf_exists & tiff_exists]
fwrite(validation, file.path(out_dir, "trajectory_validation.csv"))
message("Neutrophil trajectory figure complete")
