#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
})

set.seed(20260816)

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture/Neutrophil_Subclustering"
object_dir <- file.path(analysis_dir, "objects")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

original <- readRDS(file.path(object_dir, "neutrophils_scVI_Harmony_subclustered_annotated_clean.rds"))
canonical <- readRDS(file.path(object_dir, "canonical_neutrophils_scVI_Harmony_reclustered.rds"))

original_umap <- as.data.frame(Embeddings(original, "umap.neut.scvi.harmony"))
colnames(original_umap) <- c("UMAP_1", "UMAP_2")
original_umap$cell <- rownames(original_umap)
original_umap$annotation <- original$neutrophil_cluster_annotation[original_umap$cell]
original_umap <- original_umap[sample(seq_len(nrow(original_umap))), , drop = FALSE]

canonical_umap <- as.data.frame(Embeddings(canonical, "umap.canonical.scvi.harmony"))
colnames(canonical_umap) <- c("UMAP_1", "UMAP_2")
canonical_umap$cell <- rownames(canonical_umap)
canonical_umap$state <- canonical$neutrophil_state_annotation[canonical_umap$cell]
canonical_umap <- canonical_umap[sample(seq_len(nrow(canonical_umap))), , drop = FALSE]

annotation_order <- c(
  "Low-signal neutrophil", "Activated OSM/CXCR4 neutrophil",
  "IL1R2+ inflammatory neutrophil", "DHFR+ proliferative/stress neutrophil",
  "PADI4/translation-high neutrophil", "LTF/BPI immature neutrophil",
  "CCL3/CCL4 inflammatory neutrophil", "RORA+ atypical neutrophil"
)
annotation_colors <- c(
  "Low-signal neutrophil" = "#BDBDBD",
  "Activated OSM/CXCR4 neutrophil" = "#D55E00",
  "IL1R2+ inflammatory neutrophil" = "#E69F00",
  "DHFR+ proliferative/stress neutrophil" = "#A6761D",
  "PADI4/translation-high neutrophil" = "#CC79A7",
  "LTF/BPI immature neutrophil" = "#009E73",
  "CCL3/CCL4 inflammatory neutrophil" = "#F781BF",
  "RORA+ atypical neutrophil" = "#8C6BB1"
)
original_umap$annotation <- factor(original_umap$annotation, levels = annotation_order)

state_order <- c(
  "PADI4 neutrophil", "OSM neutrophil",
  "CXCR4 neutrophil", "MX1/ISG neutrophil"
)
state_colors <- c(
  "PADI4 neutrophil" = "#CC79A7",
  "OSM neutrophil" = "#D55E00",
  "CXCR4 neutrophil" = "#0072B2",
  "MX1/ISG neutrophil" = "#009E73"
)
canonical_umap$state <- factor(canonical_umap$state, levels = state_order)

theme_umap <- theme_void(base_family = "Arial", base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8.5),
    legend.position = "bottom",
    plot.margin = margin(8, 8, 8, 8)
  )

p_original <- ggplot(original_umap, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point_rast(size = 0.20, alpha = 0.55, raster.dpi = 400) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2)) +
  labs(
    title = "Neutrophil-only scVI-Harmony subclusters",
    subtitle = "Original 92,166-cell cleaned UMAP retained",
    color = "Subcluster annotation"
  ) +
  theme_umap

state_counts <- table(canonical$neutrophil_state_annotation)
state_labels <- setNames(
  paste0(names(state_counts), " (n = ", format(as.integer(state_counts), big.mark = ","), ")"),
  names(state_counts)
)
p_canonical <- ggplot(canonical_umap, aes(UMAP_1, UMAP_2, color = state)) +
  geom_point_rast(size = 0.28, alpha = 0.65, raster.dpi = 400) +
  scale_color_manual(values = state_colors, labels = state_labels, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1), ncol = 2)) +
  labs(
    title = "Canonical neutrophil states after reclustering",
    subtitle = "50,705 cells; Low RNA and unresolved states excluded",
    color = "Canonical neutrophil state"
  ) +
  theme_umap

p_combined <- p_original | p_canonical
p_combined <- p_combined +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    caption = "Left: original cleaned neutrophil UMAP and subclusters. Right: separately recomputed UMAP after retaining only canonical PADI4, OSM, CXCR4, and MX1/ISG states."
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 15),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0)
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
  p_combined,
  "Neutrophil_6_original_subclusters_plus_reclustered_canonical_states_UMAP",
  17, 8.5
)
message("Combined original-subcluster and canonical-state UMAP figure complete")
