#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
  library(magick)
  library(ggplotify)
  library(grid)
})

set.seed(20260816)

project_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
integrated_dir <- file.path(project_dir, "scVI_Harmony_All_Coculture")
fib_dir <- file.path(integrated_dir, "Fibroblast_Subclustering")
neut_dir <- file.path(integrated_dir, "Neutrophil_Subclustering")
analysis_dir <- file.path(project_dir, "Manuscript_Analyses_1_to_6")
source_dir <- file.path(analysis_dir, "source_data")
output_dir <- file.path(analysis_dir, "Main_And_Supplementary_Figures")
main_dir <- file.path(output_dir, "Main_Figure_3")
supp_dir <- file.path(output_dir, "Supplementary_Figures")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "#444444"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 7.5),
    axis.title = element_text(size = 8.5),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5)
  )

theme_umap <- theme_void(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "#444444"),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 6.8),
    legend.position = "bottom",
    plot.margin = margin(5, 5, 5, 5)
  )

save_figure <- function(plot, stem, width, height, dpi = 400) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         units = "in", dpi = dpi, device = ragg::agg_png, bg = "white")
  ggsave(paste0(stem, ".tiff"), plot, width = width, height = height,
         units = "in", dpi = dpi, device = ragg::agg_tiff,
         compression = "lzw", bg = "white")
}

make_image_panel <- function(path, max_width = 2600) {
  img <- image_read(path) |>
    image_background("white", flatten = TRUE)
  info <- image_info(img)
  if (info$width > max_width) img <- image_resize(img, paste0(max_width, "x"))
  ggplotify::as.ggplot(img) + theme_void() + theme(plot.margin = margin(3, 3, 3, 3))
}

# -----------------------------------------------------------------------------
# Main Figure 3A: coculture study design.
# -----------------------------------------------------------------------------
design_path <- file.path(
  project_dir, "figures", "Figure_3", "source_images",
  "Figure_3A_experimental_coculture_design_revised.png"
)
p_design <- make_image_panel(design_path, max_width = 3600)

# -----------------------------------------------------------------------------
# Main Figure 3B-C: target transcripts in the cleaned fibroblast compartment.
# -----------------------------------------------------------------------------
fib <- readRDS(file.path(fib_dir, "objects", "fibroblasts_subclustered_annotated.rds"))
DefaultAssay(fib) <- "DECONTX"
target_genes <- c("FAP", "ITGA5", "ITGB1")
expr <- as.matrix(GetAssayData(fib, assay = "DECONTX", layer = "data")[target_genes, , drop = FALSE])
expr[expr <= 1e-6] <- 0

condition_order <- c("CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(
  CF = "Control FB", UF = "UC FB", UN = "UC + NAMPTi",
  UA5 = "UC + alpha5beta1i", UD = "UC + dual"
)

target_long <- data.frame(
  cell = rep(colnames(expr), each = length(target_genes)),
  gene = rep(target_genes, times = ncol(expr)),
  expression = as.vector(expr),
  ConditionCode = rep(as.character(fib$ConditionCode), each = length(target_genes)),
  stringsAsFactors = FALSE
) |>
  mutate(
    ConditionCode = factor(ConditionCode, levels = condition_order),
    condition = factor(condition_labels[as.character(ConditionCode)], levels = unname(condition_labels)),
    gene = factor(gene, levels = target_genes)
  ) |>
  group_by(condition, gene) |>
  summarise(
    mean_expression = mean(expression),
    pct_expressed = 100 * mean(expression > 1e-6),
    .groups = "drop"
  ) |>
  group_by(gene) |>
  mutate(scaled_mean = if (sd(mean_expression) > 0) as.numeric(scale(mean_expression)) else 0) |>
  ungroup()

target_scale <- scale_color_gradient2(
  low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
  limits = c(-1.5, 1.5), oob = scales::squish
)
target_size <- scale_size_area(max_size = 8, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

p_a <- ggplot(filter(target_long, gene == "FAP"),
              aes(gene, condition, size = pct_expressed, color = scaled_mean)) +
  geom_point() + target_scale + target_size +
  labs(title = "FAP expression by culture condition", x = NULL, y = NULL,
       color = "Scaled mean", size = "Cells detected (%)") +
  theme_pub + theme(axis.text.x = element_text(face = "italic"))

p_b <- ggplot(filter(target_long, gene %in% c("ITGA5", "ITGB1")),
              aes(gene, condition, size = pct_expressed, color = scaled_mean)) +
  geom_point() + target_scale + target_size +
  labs(title = expression(alpha[5] * beta[1] * " integrin subunits"),
       subtitle = expression(italic(ITGA5) * " and " * italic(ITGB1)),
       x = NULL, y = NULL, color = "Scaled mean", size = "Cells detected (%)") +
  theme_pub + theme(axis.text.x = element_text(face = "italic"))

rm(fib, expr)
invisible(gc())

# -----------------------------------------------------------------------------
# Main Figure 3D-E: original and canonical neutrophil SCVI-Harmony UMAPs.
# -----------------------------------------------------------------------------
original_umap <- fread(
  file.path(neut_dir, "source_data", "neutrophil_cell_metadata_annotations_and_umap.csv.gz"),
  select = c("neutrophil_cluster_annotation", "neutrophil_UMAP_1", "neutrophil_UMAP_2")
)
setnames(original_umap,
         c("neutrophil_cluster_annotation", "neutrophil_UMAP_1", "neutrophil_UMAP_2"),
         c("annotation", "UMAP_1", "UMAP_2"))
original_umap <- original_umap[sample(.N)]

annotation_order <- c(
  "Low-signal neutrophil", "Activated OSM/CXCR4 neutrophil",
  "IL1R2+ inflammatory neutrophil", "DHFR+ proliferative/stress neutrophil",
  "PADI4/translation-high neutrophil", "LTF/BPI immature neutrophil",
  "CCL3/CCL4 inflammatory neutrophil", "RORA+ atypical neutrophil"
)
annotation_colors <- c(
  "Low-signal neutrophil" = "#BDBDBD", "Activated OSM/CXCR4 neutrophil" = "#D55E00",
  "IL1R2+ inflammatory neutrophil" = "#E69F00", "DHFR+ proliferative/stress neutrophil" = "#A6761D",
  "PADI4/translation-high neutrophil" = "#CC79A7", "LTF/BPI immature neutrophil" = "#009E73",
  "CCL3/CCL4 inflammatory neutrophil" = "#F781BF", "RORA+ atypical neutrophil" = "#8C6BB1"
)
original_umap[, annotation := factor(annotation, levels = annotation_order)]

p_c <- ggplot(original_umap, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point_rast(size = 0.15, alpha = 0.55, raster.dpi = 400) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 2)) +
  labs(title = "Neutrophil-only SCVI-Harmony subclusters",
       subtitle = "Original cleaned atlas retained for reference", color = "Subcluster") +
  theme_umap

canonical_umap <- fread(
  file.path(neut_dir, "source_data", "canonical_neutrophil_cell_metadata_and_umap.csv.gz"),
  select = c("cell", "neutrophil_state_annotation", "canonical_UMAP_1", "canonical_UMAP_2")
)
setnames(canonical_umap,
         c("neutrophil_state_annotation", "canonical_UMAP_1", "canonical_UMAP_2"),
         c("state", "UMAP_1", "UMAP_2"))
canonical_umap <- canonical_umap[sample(.N)]

state_order <- c("PADI4 neutrophil", "OSM neutrophil", "CXCR4 neutrophil", "MX1/ISG neutrophil")
state_colors <- c(
  "PADI4 neutrophil" = "#CC79A7", "OSM neutrophil" = "#D55E00",
  "CXCR4 neutrophil" = "#0072B2", "MX1/ISG neutrophil" = "#009E73"
)
canonical_umap[, state := factor(state, levels = state_order)]
state_counts <- table(canonical_umap$state)
state_labels <- setNames(
  paste0(names(state_counts), " (n = ", format(as.integer(state_counts), big.mark = ","), ")"),
  names(state_counts)
)

p_d <- ggplot(canonical_umap, aes(UMAP_1, UMAP_2, color = state)) +
  geom_point_rast(size = 0.20, alpha = 0.65, raster.dpi = 400) +
  scale_color_manual(values = state_colors, labels = state_labels, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1), ncol = 2)) +
  labs(title = "Canonical neutrophil states after reclustering",
       subtitle = "Low-RNA and unresolved neutrophils excluded", color = "Canonical state") +
  theme_umap

# -----------------------------------------------------------------------------
# Main Figure 3F: focused donor-level compositional effects.
# -----------------------------------------------------------------------------
composition <- fread(file.path(source_dir, "analysis2_canonical_state_CLR_differential_abundance.csv"))
main_contrast_ids <- c(
  "N_vs_PRE", "CF_vs_N", "UF_vs_N", "UA5_vs_UF", "UD_vs_UF",
  "alpha5beta1_fibroblast_interaction"
)
main_contrast_labels <- c(
  N_vs_PRE = "N vs PRE", CF_vs_N = "Control FB vs N", UF_vs_N = "UC FB vs N",
  UA5_vs_UF = "alpha5beta1i vs UC", UD_vs_UF = "Dual vs UC",
  alpha5beta1_fibroblast_interaction = "FB-dependent\nalpha5beta1i"
)
composition_main <- composition[contrast_id %in% main_contrast_ids]
composition_main[, contrast_short := factor(main_contrast_labels[contrast_id],
                                             levels = unname(main_contrast_labels))]
composition_main[, state := factor(state, levels = rev(state_order))]

p_e <- ggplot(composition_main, aes(contrast_short, state, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(adj.P.Val < 0.001, "***",
                               ifelse(adj.P.Val < 0.01, "**",
                                      ifelse(adj.P.Val < 0.05, "*", "")))), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-3.1, 3.1), oob = scales::squish) +
  labs(title = "Canonical-state differential abundance",
       subtitle = "Donor-level centered-log-ratio effects; * FDR <0.05",
       x = NULL, y = NULL, fill = "CLR effect") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))

# -----------------------------------------------------------------------------
# Main Figure 3G: focused functional programs for blockade contrasts.
# -----------------------------------------------------------------------------
programs <- fread(file.path(source_dir, "analysis3_functional_program_differential_activity.csv"))
program_ids <- c("CXCR4_aging_retention", "OSM_inflammation", "Type_I_interferon",
                 "Degranulation", "Immature_granulopoiesis", "Proliferation_stress")
program_labels <- c(
  CXCR4_aging_retention = "CXCR4 aging/retention", OSM_inflammation = "OSM inflammation",
  Type_I_interferon = "Type I interferon", Degranulation = "Degranulation",
  Immature_granulopoiesis = "Immature granulopoiesis", Proliferation_stress = "Proliferation/stress"
)
program_contrasts <- c("UA5_vs_UF", "UD_vs_UF", "alpha5beta1_fibroblast_interaction")
program_contrast_labels <- c(
  UA5_vs_UF = "alpha5beta1i vs UC", UD_vs_UF = "Dual vs UC",
  alpha5beta1_fibroblast_interaction = "FB-dependent alpha5beta1i"
)
program_main <- programs[program %in% program_ids & contrast_id %in% program_contrasts]
program_main[, program_label := factor(program_labels[program], levels = rev(unname(program_labels)))]
program_main[, contrast_short := factor(program_contrast_labels[contrast_id],
                                        levels = unname(program_contrast_labels))]
program_main[, state_short := sub(" neutrophil$", "", state)]
program_main[, state_short := factor(state_short, levels = c("PADI4", "OSM", "CXCR4", "MX1/ISG"))]

p_f <- ggplot(program_main, aes(contrast_short, program_label, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(data = program_main[adj.P.Val < 0.05], shape = 8, size = 1.7) +
  facet_wrap(~state_short, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1.25, 1.25), oob = scales::squish,
                       na.value = "#E5E5E5") +
  labs(title = "Functional-program response to alpha5beta1 blockade",
       subtitle = "State-specific donor-level effects; * FDR <0.05; gray = not estimable",
       x = NULL, y = NULL, fill = "Effect") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                    strip.text = element_text(face = "bold", size = 8))

# -----------------------------------------------------------------------------
# Main Figure 3H-J: inferred neutrophil progression.
# -----------------------------------------------------------------------------
trajectory_cells <- fread(file.path(source_dir, "analysis6_cell_pseudotime_and_lineage_weights.csv.gz"),
                          select = c("cell", "state", "weighted_pseudotime"))
trajectory_umap <- merge(
  canonical_umap[, .(cell, UMAP_1, UMAP_2)],
  trajectory_cells,
  by = "cell"
)
centroids <- trajectory_umap[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = state]
centroids[, state := factor(state, levels = state_order)]
setorder(centroids, state)
centroids[, state_label := sub(" neutrophil$", "", as.character(state))]
centroids[, lineage := "Lineage1"]
trajectory_plot_cells <- trajectory_umap[sample(.N, min(.N, 22000))]

p_g <- ggplot(trajectory_plot_cells, aes(UMAP_1, UMAP_2, color = weighted_pseudotime)) +
  geom_point_rast(size = 0.18, alpha = 0.65, raster.dpi = 400) +
  geom_path(data = centroids, aes(UMAP_1, UMAP_2, group = lineage), inherit.aes = FALSE,
            color = "black", linewidth = 0.9,
            arrow = arrow(length = unit(0.10, "inches"), type = "closed")) +
  geom_point(data = centroids, aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 2) +
  geom_text(data = centroids, aes(UMAP_1, UMAP_2, label = state_label), inherit.aes = FALSE,
            nudge_y = 0.22, size = 2.6, fontface = "bold", check_overlap = TRUE) +
  scale_color_viridis_c(option = "C") + coord_equal() +
  labs(title = "SCVI-Harmony neutrophil trajectory",
       subtitle = "Slingshot rooted in the PADI4 state",
       x = "UMAP 1", y = "UMAP 2", color = "Pseudotime") + theme_pub

sample_summary <- fread(file.path(source_dir, "analysis6_trajectory_summary_by_sample.csv"))
trajectory_condition_order <- c("PRE", "N", "NNAMPT", "NA5", "CF", "UF", "UN", "UA5", "UD")
trajectory_condition_labels <- c(
  PRE = "PRE", N = "N", NNAMPT = "N+NAMPTi", NA5 = "N+alpha5beta1i",
  CF = "Control FB", UF = "UC FB", UN = "UC+NAMPTi", UA5 = "UC+alpha5beta1i", UD = "Dual"
)
sample_summary[, condition := factor(trajectory_condition_labels[ConditionCode],
                                     levels = unname(trajectory_condition_labels))]

p_h <- ggplot(sample_summary, aes(condition, median_pseudotime, group = DonorID)) +
  geom_line(color = "#BDBDBD", linewidth = 0.35, alpha = 0.65) +
  geom_point(aes(color = ConditionCode), size = 1.7) +
  geom_boxplot(aes(group = condition), width = 0.5, outlier.shape = NA,
               fill = NA, color = "black", inherit.aes = TRUE) +
  guides(color = "none") +
  labs(title = "Donor-level trajectory position",
       subtitle = "Lines connect paired culture donors; PRE donors are independent",
       x = NULL, y = "Median pseudotime") +
  theme_pub + theme(axis.text.x = element_text(angle = 40, hjust = 1))

trajectory_results <- fread(file.path(source_dir, "analysis6_differential_pseudotime_and_lineage_occupancy.csv"))
trajectory_main <- trajectory_results[contrast_id %in% main_contrast_ids]
trajectory_main[, contrast_short := factor(main_contrast_labels[contrast_id],
                                            levels = unname(main_contrast_labels))]

p_i <- ggplot(trajectory_main, aes(contrast_short, "Median pseudotime", fill = effect)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f%s", effect, ifelse(FDR < 0.05, "*", ""))), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-2.1, 2.1), oob = scales::squish) +
  labs(title = "Differential inferred progression",
       subtitle = "Effect estimates; * FDR <0.05", x = NULL, y = NULL, fill = "Effect") +
  theme_pub + theme(axis.text.x = element_text(angle = 40, hjust = 1),
                    axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Four-row manuscript composition. The design schematic is restored as panel A;
# target-engagement panels B-C share its row, followed by states, phenotypes,
# and trajectory panels D-J.
target_stack <- (p_a / p_b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

figure_3 <-
  ((p_design | target_stack) + plot_layout(widths = c(1.30, 0.70))) /
  (p_c | p_d) /
  (p_e | p_f) /
  (p_g | p_h | p_i) +
  plot_layout(heights = c(1.28, 1.12, 1.05, 0.92)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

main_stem <- file.path(main_dir, "Figure_3_fibroblast_targets_neutrophil_states_and_trajectory")
save_figure(figure_3, main_stem, width = 18, height = 24, dpi = 400)

# -----------------------------------------------------------------------------
# Supplementary Figure S3.1: fibroblast annotation, composition, and targets.
# -----------------------------------------------------------------------------
fib_fig_dir <- file.path(fib_dir, "figures")
s31a <-
  (make_image_panel(file.path(fib_fig_dir, "Fibroblast_UMAP_1_annotations.png")) |
     make_image_panel(file.path(fib_fig_dir, "Fibroblast_UMAP_2_conditions.png"))) /
  (make_image_panel(file.path(fib_fig_dir, "Fibroblast_3_annotation_marker_dotplot.png")) |
     make_image_panel(file.path(fib_fig_dir, "Fibroblast_4_subtype_composition_by_condition.png"))) +
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(s31a, file.path(supp_dir, "Supplementary_Figure_S3_1A_fibroblast_annotation_and_composition"),
            width = 18, height = 14, dpi = 400)

s31b <- make_image_panel(file.path(fib_fig_dir, "Fibroblast_5_FAP_alpha5beta1_expression.png")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(s31b, file.path(supp_dir, "Supplementary_Figure_S3_1B_fibroblast_target_transcripts"),
            width = 14, height = 8, dpi = 400)

# -----------------------------------------------------------------------------
# Supplementary Figure S3.2: integration, neutrophil annotation, and cleanup.
# -----------------------------------------------------------------------------
cleanup <- fread(file.path(neut_dir, "source_data", "neutrophil_cleanup_audit.csv"))
cleanup_summary <- cleanup[, .(cells = sum(cells)), by = .(cleanup_action, cleanup_reason)]
p_cleanup <- ggplot(cleanup_summary,
                    aes(reorder(cleanup_reason, cells), cells, fill = cleanup_action)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_fill_manual(values = c(Retained = "#4DAF4A", Removed = "#999999")) +
  labs(title = "Neutrophil cleanup audit", x = NULL, y = "Cells", fill = "Action") +
  theme_pub

s32a <-
  (make_image_panel(file.path(integrated_dir, "figures", "scVI_Harmony_UMAP_1_cell_compartment.png")) |
     make_image_panel(file.path(integrated_dir, "figures", "scVI_Harmony_UMAP_4_condition_split.png"))) /
  (make_image_panel(file.path(neut_dir, "figures", "Neutrophil_1_scVI_Harmony_subcluster_annotations.png")) |
     make_image_panel(file.path(neut_dir, "figures", "Neutrophil_2_annotation_marker_dotplot.png"))) +
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(s32a, file.path(supp_dir, "Supplementary_Figure_S3_2A_neutrophil_integration_and_annotation"),
            width = 18, height = 14, dpi = 400)

s32b <-
  make_image_panel(file.path(neut_dir, "figures", "Neutrophil_3_OSM_CXCR4_MX1_PADI4_by_condition.png")) |
  p_cleanup +
  plot_layout(widths = c(1.45, 0.75)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(s32b, file.path(supp_dir, "Supplementary_Figure_S3_2B_neutrophil_condition_markers_and_cleanup"),
            width = 18, height = 8, dpi = 400)

# -----------------------------------------------------------------------------
# Supplementary Figure S3.3: complete CLR and MiloR differential abundance.
# -----------------------------------------------------------------------------
s33 <-
  make_image_panel(file.path(analysis_dir, "figures", "Analysis2_compositional_differential_abundance.png"), 3200) /
  make_image_panel(file.path(neut_dir, "figures", "Neutrophil_7_MiloR_DA_heatmap_FigureA_subclusters_and_FigureB_canonical_states.png"), 3200) +
  plot_layout(heights = c(0.9, 1.1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(s33, file.path(supp_dir, "Supplementary_Figure_S3_3_complete_differential_abundance"),
            width = 20, height = 17, dpi = 400)

# -----------------------------------------------------------------------------
# Supplementary Figures S3.4-S3.7: full analysis panels retained without cropping.
# -----------------------------------------------------------------------------
s34 <- make_image_panel(file.path(analysis_dir, "figures", "Analysis1_state_specific_pseudobulk_DE.png"), 3600)
save_figure(s34, file.path(supp_dir, "Supplementary_Figure_S3_4_state_specific_pseudobulk_DE"),
            width = 20, height = 13, dpi = 400)

s35a <- make_image_panel(file.path(analysis_dir, "figures", "Analysis3_functional_program_activity.png"), 3600)
save_figure(s35a, file.path(supp_dir, "Supplementary_Figure_S3_5A_complete_functional_programs"),
            width = 20, height = 13, dpi = 400)
s35b <- make_image_panel(file.path(analysis_dir, "figures", "Analysis4_metacell_gene_modules.png"), 3600)
save_figure(s35b, file.path(supp_dir, "Supplementary_Figure_S3_5B_unbiased_metacell_gene_modules"),
            width = 20, height = 15, dpi = 400)

s36 <- make_image_panel(file.path(analysis_dir, "figures", "Analysis5_fibroblast_neutrophil_signaling.png"), 3600)
save_figure(s36, file.path(supp_dir, "Supplementary_Figure_S3_6_fibroblast_neutrophil_interactome"),
            width = 18, height = 20, dpi = 400)

s37 <- make_image_panel(file.path(analysis_dir, "figures", "Analysis6_neutrophil_trajectory.png"), 3600)
save_figure(s37, file.path(supp_dir, "Supplementary_Figure_S3_7_extended_neutrophil_trajectory"),
            width = 18, height = 20, dpi = 400)

message("Manuscript Figure 3 and supplementary figures completed: ", output_dir)
