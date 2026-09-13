#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
  library(magick)
  library(ggplotify)
})

set.seed(20260816)

project_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
integrated_dir <- file.path(project_dir, "scVI_Harmony_All_Coculture")
fib_dir <- file.path(integrated_dir, "Fibroblast_Subclustering")
neut_dir <- file.path(integrated_dir, "Neutrophil_Subclustering")
base_analysis <- file.path(project_dir, "Manuscript_Analyses_1_to_6")
additional_dir <- file.path(project_dir, "Additional_High_Impact_Analyses_2_to_5")
additional_source <- file.path(additional_dir, "source_data")
main_dir <- file.path(base_analysis, "Main_And_Supplementary_Figures", "Main_Figure_3")
panel_source_dir <- file.path(main_dir, "source_data_revised_A_to_I")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_source_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10.2),
    plot.subtitle = element_text(size = 7.2, color = "#444444"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(face = "bold", size = 7.5),
    legend.text = element_text(size = 6.7),
    strip.text = element_text(face = "bold", size = 7.5),
    plot.margin = margin(4, 4, 4, 4)
  )

theme_umap_compact <- theme_void(base_family = "Arial", base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 7.2, color = "#444444"),
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 5.7),
    legend.key.height = unit(0.12, "in"),
    legend.key.width = unit(0.12, "in"),
    legend.position = "bottom",
    plot.margin = margin(3, 3, 3, 3)
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

make_image_panel <- function(path, max_width = 3600, crop_left = 0) {
  img <- image_read(path) |> image_background("white", flatten = TRUE)
  info <- image_info(img)
  if (crop_left > 0) {
    img <- image_crop(
      img,
      geometry_area(
        width = info$width - crop_left,
        height = info$height,
        x_off = crop_left,
        y_off = 0
      )
    )
    info <- image_info(img)
  }
  if (info$width > max_width) img <- image_resize(img, paste0(max_width, "x"))
  ggplotify::as.ggplot(img) + theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
}

diverging_scale <- function(limits, name = "Effect") {
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = limits, oob = scales::squish, name = name
  )
}

# A — experimental design.
design_path <- file.path(project_dir, "figures", "Figure_3", "source_images",
                         "Figure_3A_experimental_coculture_design_revised.png")
design_img <- image_read(design_path) |> image_background("white", flatten = TRUE)
p_a <- ggplotify::as.ggplot(
  grid::rasterGrob(
    as.raster(design_img),
    x = grid::unit(0.53, "npc"), y = grid::unit(0.50, "npc"),
    width = grid::unit(0.90, "npc"), height = grid::unit(0.90, "npc"),
    interpolate = TRUE
  )
) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

# B — compact combined fibroblast target-transcript dot plot.
fib <- readRDS(file.path(fib_dir, "objects", "fibroblasts_subclustered_annotated.rds"))
DefaultAssay(fib) <- "DECONTX"
target_genes <- c("FAP", "ITGA5", "ITGB1")
expr <- as.matrix(GetAssayData(fib, assay = "DECONTX", layer = "data")[target_genes, , drop = FALSE])
expr[expr <= 1e-6] <- 0
condition_order <- c("CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(CF = "Control FB", UF = "UC FB", UN = "UC + NAMPTi",
                      UA5 = "UC + \u03b15\u03b21i", UD = "UC + dual")
target_long <- data.frame(
  gene = rep(target_genes, times = ncol(expr)), expression = as.vector(expr),
  ConditionCode = rep(as.character(fib$ConditionCode), each = length(target_genes))
) |>
  mutate(ConditionCode = factor(ConditionCode, levels = condition_order),
         condition = factor(condition_labels[as.character(ConditionCode)],
                            levels = unname(condition_labels)),
         gene = factor(gene, levels = target_genes)) |>
  group_by(condition, gene) |>
  summarise(mean_expression = mean(expression),
            pct_expressed = 100 * mean(expression > 1e-6), .groups = "drop") |>
  group_by(gene) |>
  mutate(scaled_mean = if (sd(mean_expression) > 0) as.numeric(scale(mean_expression)) else 0) |>
  ungroup()
fwrite(as.data.table(target_long), file.path(panel_source_dir, "Figure_3B_fibroblast_targets.csv"))

p_b <- ggplot(target_long, aes(gene, condition, size = pct_expressed, color = scaled_mean)) +
  geom_point() +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                        midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish) +
  scale_size_area(max_size = 7, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  labs(title = "Fibroblast target transcripts by culture condition",
       subtitle = expression(italic(FAP) * ", " * italic(ITGA5) * " and " * italic(ITGB1)),
       x = NULL, y = NULL, color = "Scaled mean", size = "Cells detected (%)") +
  theme_pub + theme(axis.text.x = element_text(face = "italic"), legend.position = "right")

rm(fib, expr)
invisible(gc())

# C–D — compact original and canonical neutrophil SCVI-Harmony UMAPs.
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
  geom_point_rast(size = 0.12, alpha = 0.50, raster.dpi = 400) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.0, alpha = 1), ncol = 2)) +
  labs(title = "Neutrophil-only SCVI-Harmony subclusters",
       subtitle = "Original cleaned atlas retained for reference", color = "Subcluster") +
  theme_umap_compact

canonical_umap <- fread(
  file.path(neut_dir, "source_data", "canonical_neutrophil_cell_metadata_and_umap.csv.gz"),
  select = c("neutrophil_state_annotation", "canonical_UMAP_1", "canonical_UMAP_2")
)
setnames(canonical_umap,
         c("neutrophil_state_annotation", "canonical_UMAP_1", "canonical_UMAP_2"),
         c("state", "UMAP_1", "UMAP_2"))
canonical_umap <- canonical_umap[sample(.N)]
state_order <- c("PADI4 neutrophil", "OSM neutrophil", "CXCR4 neutrophil", "MX1/ISG neutrophil")
state_colors <- c("PADI4 neutrophil" = "#CC79A7", "OSM neutrophil" = "#D55E00",
                  "CXCR4 neutrophil" = "#0072B2", "MX1/ISG neutrophil" = "#009E73")
canonical_umap[, state := factor(state, levels = state_order)]
state_legend_labels <- setNames(
  sub(" neutrophil$", "", state_order), state_order
)
p_d <- ggplot(canonical_umap, aes(UMAP_1, UMAP_2, color = state)) +
  geom_point_rast(size = 0.18, alpha = 0.50, raster.dpi = 400) +
  scale_color_manual(values = state_colors, labels = state_legend_labels, drop = FALSE) +
  guides(color = guide_legend(
    override.aes = list(size = 2.3, alpha = 1),
    ncol = 1, byrow = TRUE, title.position = "top", label.position = "right"
  )) +
  coord_equal() +
  labs(title = "Canonical neutrophil states after reclustering",
       subtitle = "Low-RNA and unresolved neutrophils excluded", color = "Canonical state") +
  theme_umap_compact +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = grid::unit(0, "pt"),
    legend.key.height = grid::unit(0.14, "in"),
    legend.key.width = grid::unit(0.09, "in"),
    legend.title = element_text(size = 6.4),
    legend.text = element_text(size = 5.8)
  )

# Keep the legend outside the UMAP data cloud while positioning it tightly
# beside the right edge of the embedding.
p_d_legend <- cowplot::get_legend(p_d)
p_d_compact <- cowplot::ggdraw() +
  cowplot::draw_plot(p_d + theme(legend.position = "none"), 0, 0, 1, 1) +
  cowplot::draw_grob(p_d_legend, x = 0.78, y = 0.38, width = 0.18, height = 0.25)

# D — canonical-state differential abundance restored from the previous Figure 3F.
composition <- fread(file.path(base_analysis, "source_data",
                               "analysis2_canonical_state_CLR_differential_abundance.csv"))
main_contrast_ids <- c("N_vs_PRE", "CF_vs_N", "UF_vs_N", "UA5_vs_UF", "UD_vs_UF",
                       "alpha5beta1_fibroblast_interaction")
main_contrast_labels <- c(
  N_vs_PRE = "N vs PRE", CF_vs_N = "Control FB vs N", UF_vs_N = "UC FB vs N",
  UA5_vs_UF = "\u03b15\u03b21i vs UC", UD_vs_UF = "Dual vs UC",
  alpha5beta1_fibroblast_interaction = "FB-dependent\n\u03b15\u03b21i"
)
panel_d_abundance <- composition[contrast_id %in% main_contrast_ids]
panel_d_abundance[, contrast_label := factor(main_contrast_labels[contrast_id],
                                              levels = unname(main_contrast_labels))]
panel_d_abundance[, state := factor(state, levels = rev(state_order))]
fwrite(panel_d_abundance,
       file.path(panel_source_dir, "Figure_3D_restored_canonical_state_abundance.csv"))
p_previous_f <- ggplot(panel_d_abundance, aes(contrast_label, state, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(adj.P.Val < 0.001, "***",
                               ifelse(adj.P.Val < 0.01, "**",
                                      ifelse(adj.P.Val < 0.05, "*", "")))), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-3.1, 3.1), oob = scales::squish) +
  labs(title = "Canonical-state differential abundance",
       subtitle = "Condition-associated shifts in canonical states",
       x = NULL, y = NULL, fill = "CLR effect") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                    axis.text.y = element_text(size = 6.7),
                    legend.key.height = grid::unit(18, "pt"),
                    legend.key.width = grid::unit(7, "pt"))

# E — functional programs restored from the previous Figure 3G.
programs <- fread(file.path(base_analysis, "source_data",
                            "analysis3_functional_program_differential_activity.csv"))
program_ids <- c("CXCR4_aging_retention", "OSM_inflammation", "Type_I_interferon",
                 "Degranulation", "Immature_granulopoiesis", "Proliferation_stress")
program_labels <- c(
  CXCR4_aging_retention = "CXCR4 aging/retention", OSM_inflammation = "OSM inflammation",
  Type_I_interferon = "Type I interferon", Degranulation = "Degranulation",
  Immature_granulopoiesis = "Immature granulopoiesis", Proliferation_stress = "Proliferation/stress"
)
program_contrasts <- c("UA5_vs_UF", "UD_vs_UF", "alpha5beta1_fibroblast_interaction")
program_contrast_labels <- c(
  UA5_vs_UF = "\u03b15\u03b21i vs UC", UD_vs_UF = "Dual vs UC",
  alpha5beta1_fibroblast_interaction = "FB-dependent \u03b15\u03b21i"
)
panel_e_programs <- programs[program %in% program_ids & contrast_id %in% program_contrasts]
panel_e_programs[, program_label := factor(program_labels[program],
                                            levels = rev(unname(program_labels)))]
panel_e_programs[, contrast_label := factor(program_contrast_labels[contrast_id],
                                             levels = unname(program_contrast_labels))]
panel_e_programs[, state_label := factor(sub(" neutrophil$", "", state),
                                          levels = c("PADI4", "OSM", "CXCR4", "MX1/ISG"))]
fwrite(panel_e_programs,
       file.path(panel_source_dir, "Figure_3E_restored_functional_programs.csv"))
p_previous_g <- ggplot(panel_e_programs, aes(contrast_label, program_label, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(data = panel_e_programs[adj.P.Val < 0.05], shape = 8, size = 1.25) +
  facet_wrap(~state_label, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1.25, 1.25), oob = scales::squish,
                       na.value = "#E5E5E5") +
  labs(title = "Functional-program response to \u03b15\u03b21 blockade",
       subtitle = "State-specific functional responses",
       x = NULL, y = NULL, fill = "Effect") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                    axis.text.y = element_text(size = 6.2),
                    legend.key.height = grid::unit(18, "pt"),
                    legend.key.width = grid::unit(7, "pt"))

# F — formal factorial canonical-state abundance effects.
factorial <- fread(file.path(additional_source, "analysis4_unified_factorial_blockade_effects.csv"))
panel_e <- factorial[family == "Canonical-state CLR abundance"]
contrast_e_labels <- c(
  alpha5beta1_fibroblast_interaction = "FB-dependent \u03b15\u03b21i",
  NAMPT_fibroblast_interaction = "FB-dependent NAMPTi",
  alpha5beta1_NAMPT_synergy_UC = "\u03b15\u03b21i \u00d7 NAMPTi"
)
panel_e[, contrast_label := factor(contrast_e_labels[contrast_id],
                                    levels = unname(contrast_e_labels))]
panel_e[, state_label := factor(sub(" neutrophil$", "", state),
                                levels = rev(c("PADI4", "OSM", "CXCR4", "MX1/ISG")))]
fwrite(panel_e, file.path(panel_source_dir, "Figure_3F_factorial_state_abundance.csv"))
p_e <- ggplot(panel_e, aes(contrast_label, state_label, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%+.2f%s", effect, ifelse(FDR < 0.05, "*", ""))), size = 2.35) +
  diverging_scale(c(-1.35, 1.35), "CLR effect") +
  labs(title = "Factorial effects on canonical-state abundance",
       subtitle = "Fibroblast-dependent and inhibitor interaction effects",
       x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                    legend.key.height = grid::unit(18, "pt"),
                    legend.key.width = grid::unit(7, "pt"))

# G — selected pathway and TF footprints in CXCR4 and OSM neutrophils.
activity <- fread(file.path(additional_source, "analysis2_regulator_pathway_differential_activity.csv"))
selected_regulators <- c("NFkB", "TNFa", "TGFb", "MAPK", "JAK-STAT",
                         "AP1", "RELA", "HIF1A", "CEBPB", "STAT3")
regulator_labels <- c(
  NFkB = "Pathway | NFkB", TNFa = "Pathway | TNFa", TGFb = "Pathway | TGFb",
  MAPK = "Pathway | MAPK", `JAK-STAT` = "Pathway | JAK-STAT",
  AP1 = "TF | JUN/AP-1", RELA = "TF | RELA", HIF1A = "TF | HIF1A",
  CEBPB = "TF | CEBPB", STAT3 = "TF | STAT3"
)
contrast_f_labels <- c(CF_vs_N = "Control FB vs N", UF_vs_N = "UC FB vs N",
                       alpha5beta1_fibroblast_interaction = "FB-dependent \u03b15\u03b21i")
panel_f <- activity[
  state %in% c("CXCR4 neutrophil", "OSM neutrophil") &
    contrast_id %in% names(contrast_f_labels) & regulator %in% selected_regulators
]
panel_f[, regulator_label := factor(regulator_labels[regulator],
                                    levels = rev(unname(regulator_labels)))]
panel_f[, contrast_label := factor(contrast_f_labels[contrast_id],
                                   levels = unname(contrast_f_labels))]
panel_f[, state_label := factor(sub(" neutrophil$", "", state), levels = c("CXCR4", "OSM"))]
fwrite(panel_f, file.path(panel_source_dir, "Figure_3G_selected_regulatory_activity.csv"))
p_f <- ggplot(panel_f, aes(contrast_label, regulator_label, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(data = panel_f[FDR < 0.05], shape = 8, size = 1.25) +
  facet_wrap(~state_label, ncol = 2) +
  diverging_scale(c(-7.6, 7.6), "Activity effect") +
  labs(title = "Selected fibroblast- and \u03b15\u03b21-responsive regulators",
       subtitle = "CXCR4 and OSM neutrophils",
       x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                    axis.text.y = element_text(size = 6.5),
                    legend.key.height = grid::unit(18, "pt"),
                    legend.key.width = grid::unit(7, "pt"))

# H — donor-level pseudotime distribution shifts.
quantile_effects <- fread(file.path(additional_source, "analysis5_donor_pseudotime_quantile_effects.csv"))
contrast_g_labels <- c(CF_vs_N = "Control FB vs N", UF_vs_N = "UC FB vs N",
                       UN_vs_UF = "NAMPTi vs UC", UA5_vs_UF = "\u03b15\u03b21i vs UC",
                       UD_vs_UF = "Dual vs UC")
panel_g <- quantile_effects[contrast_id %in% names(contrast_g_labels)]
panel_g[, contrast_label := factor(contrast_g_labels[contrast_id],
                                   levels = unname(contrast_g_labels))]
panel_g[, quantile := factor(quantile, levels = c("Q90", "Q75", "Q50", "Q25", "Q10"))]
fwrite(panel_g, file.path(panel_source_dir, "Figure_3H_pseudotime_quantile_effects.csv"))
p_g <- ggplot(panel_g, aes(contrast_label, quantile, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.2f%s", effect, ifelse(FDR < 0.05, "*", ""))), size = 2.1) +
  diverging_scale(c(-1.1, 3.0), "Pseudotime effect") +
  labs(title = "Donor-level shifts across the pseudotime distribution",
       subtitle = "Positive values indicate later progression",
       x = NULL, y = "Pseudotime quantile") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                    legend.key.height = grid::unit(18, "pt"),
                    legend.key.width = grid::unit(7, "pt"))

# I — representative alpha5beta1-responsive gene trajectories.
curve_summary <- fread(file.path(additional_source,
                                 "analysis5_condition_specific_pseudotime_gene_curves.csv.gz"))
panel_h <- curve_summary[
  contrast_id == "UA5_vs_UF" & gene %in% c("PRIM1", "GRIP1", "COMMD1") &
    ConditionCode %in% c("UF", "UA5")
]
panel_h[, condition_label := factor(
  fifelse(ConditionCode == "UF", "UC fibroblasts", "UC + \u03b15\u03b21i"),
  levels = c("UC fibroblasts", "UC + \u03b15\u03b21i")
)]
panel_h[, gene := factor(gene, levels = c("PRIM1", "GRIP1", "COMMD1"))]
panel_h_end <- panel_h[order(pt_bin), .SD[.N], by = .(gene, condition_label)]
panel_h_end[, direct_label := fifelse(condition_label == "UC fibroblasts", "UC FB", "UC + \u03b15\u03b21i")]
fwrite(panel_h, file.path(panel_source_dir, "Figure_3I_alpha5beta1_gene_trajectories.csv"))
p_h <- ggplot(panel_h, aes(pt_bin, mean_expression, color = condition_label,
                           fill = condition_label, group = condition_label)) +
  geom_ribbon(aes(ymin = pmax(mean_expression - se_expression, 0),
                  ymax = mean_expression + se_expression),
              alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.95) + geom_point(size = 1.25) +
  geom_text(data = panel_h_end, aes(label = direct_label), hjust = -0.08,
            size = 2.15, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~gene, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("UC fibroblasts" = "#E76F61", "UC + \u03b15\u03b21i" = "#00A9B7")) +
  scale_fill_manual(values = c("UC fibroblasts" = "#E76F61", "UC + \u03b15\u03b21i" = "#00A9B7")) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.20))) +
  labs(title = "\u03b15\u03b21 blockade alters gene-expression trajectories",
       subtitle = "\u03b15\u03b21-responsive expression across pseudotime",
       x = "Pseudotime bin", y = "Mean log-normalized expression",
       color = NULL, fill = NULL) +
  theme_pub + theme(legend.position = "none")

# J — fibroblast-to-neutrophil signaling routes that change in opposite
# directions in UC versus control co-culture and after alpha5beta1 blockade.
# Scores are the consensus of four complementary ligand-receptor metrics.
lr_diff <- fread(file.path(additional_source,
                           "analysis3_differential_consensus_LR_scores.csv.gz"))
lr_reversal <- dcast(
  lr_diff[contrast_id %in% c("UF_vs_CF", "UA5_vs_UF")],
  fibroblast_subtype + neutrophil_state + ligand + receptor + pathway_name ~ contrast_id,
  value.var = "delta_consensus"
)
lr_reversal <- lr_reversal[
  is.finite(UF_vs_CF) & is.finite(UA5_vs_UF) &
    sign(UF_vs_CF) != sign(UA5_vs_UF)
]
lr_reversal[, reversal_strength := pmin(abs(UF_vs_CF), abs(UA5_vs_UF))]
lr_reversal <- lr_reversal[order(-reversal_strength), head(.SD, 1), by = pathway_name]
setorder(lr_reversal, -reversal_strength)
lr_reversal <- head(lr_reversal, 6)
lr_reversal[, sender_short := fifelse(
  grepl("COL1A1", fibroblast_subtype), "ECM-FB", "CCN-FB"
)]
lr_reversal[, state_short := sub(" neutrophil$", "", neutrophil_state)]
lr_reversal[, route_label := paste0(
  "[", sender_short, "] ", ligand, " -> ", receptor, " | ", state_short
)]
route_order <- rev(lr_reversal$route_label)
panel_j <- melt(
  lr_reversal,
  id.vars = c("fibroblast_subtype", "neutrophil_state", "ligand", "receptor",
              "pathway_name", "reversal_strength", "route_label"),
  measure.vars = c("UF_vs_CF", "UA5_vs_UF"),
  variable.name = "contrast_id", value.name = "consensus_change"
)
panel_j[, contrast_label := factor(
  contrast_id,
  levels = c("UF_vs_CF", "UA5_vs_UF"),
  labels = c("UC vs\ncontrol", "\u03b15\u03b21i\nvs UC")
)]
panel_j[, route_label := factor(route_label, levels = route_order)]
fwrite(panel_j, file.path(panel_source_dir,
                          "Figure_3J_fibroblast_neutrophil_signaling_reversals.csv"))
p_j <- ggplot(panel_j, aes(contrast_label, route_label, fill = consensus_change)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%+.2f", consensus_change)), size = 2.05) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-0.85, 0.85), oob = scales::squish, name = "Consensus\nLR change"
  ) +
  labs(
    title = "\u03b15\u03b21 reverses fibroblast-neutrophil signals",
    subtitle = "Signals reversed by \u03b15\u03b21 blockade",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(18, "pt"),
    legend.key.width = grid::unit(7, "pt"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.4),
    axis.text.y = element_text(size = 5.8),
    plot.title = element_text(size = 10.2),
    plot.subtitle = element_text(size = 7.2)
  )

# Revised A-J layout. The enlarged canonical UMAP is paired with differential
# abundance; coord_equal() preserves its native geometry without stretching.
figure_3_revised <-
  ((wrap_elements(full = p_a) | p_b) + plot_layout(widths = c(1.18, 0.82))) /
  ((p_d_compact | p_previous_f) + plot_layout(widths = c(1.16, 0.84))) /
  ((p_previous_g | p_e) + plot_layout(widths = c(1.16, 0.84))) /
  ((p_f | p_g) + plot_layout(widths = c(1.16, 0.84))) /
  ((p_h | p_j) + plot_layout(widths = c(0.64, 0.36))) +
  plot_layout(heights = c(0.96, 1.08, 0.92, 0.92, 0.78)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

main_stem <- file.path(main_dir, "Figure_3_revised_canonical_UMAP_with_restored_abundance_and_program_panels")
save_figure(figure_3_revised, main_stem, width = 18, height = 27.00, dpi = 400)

message("Revised Figure 3 A-J completed: ", main_stem)
