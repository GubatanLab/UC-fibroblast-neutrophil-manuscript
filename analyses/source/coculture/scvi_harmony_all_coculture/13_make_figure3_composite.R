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
})

set.seed(20260816)

project_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
integrated_dir <- file.path(project_dir, "scVI_Harmony_All_Coculture")
fib_dir <- file.path(integrated_dir, "Fibroblast_Subclustering")
neut_dir <- file.path(integrated_dir, "Neutrophil_Subclustering")
figure_dir <- file.path(project_dir, "figures", "Figure_3")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

schematic_file <- file.path(
  figure_dir, "source_images", "Figure_3A_experimental_coculture_design.png"
)

# Figure 3A: user-supplied experimental co-culture schematic, flattened to white.
schematic <- image_read(schematic_file) |>
  image_background("white", flatten = TRUE)
p_a <- ggplotify::as.ggplot(schematic) +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))

# Figure 3B: FAP and alpha5beta1-subunit expression by fibroblast culture condition.
fib <- readRDS(file.path(fib_dir, "objects", "fibroblasts_subclustered_annotated.rds"))
DefaultAssay(fib) <- "DECONTX"
target_genes <- c("FAP", "ITGA5", "ITGB1")
stopifnot(all(target_genes %in% rownames(fib)))

condition_order <- c("CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(
  CF = "Control fibroblast\nco-culture",
  UF = "UC fibroblast\nco-culture",
  UN = "UC +\nNAMPTi",
  UA5 = "UC +\nalpha5beta1i",
  UD = "UC +\ndual blockade"
)

expr <- as.matrix(GetAssayData(
  fib, assay = "DECONTX", layer = "data"
)[target_genes, , drop = FALSE])
detection_threshold <- 1e-6
expr[expr <= detection_threshold] <- 0

target_long <- data.frame(
  cell = rep(colnames(expr), each = length(target_genes)),
  gene = rep(target_genes, times = ncol(expr)),
  expression = as.vector(expr),
  ConditionCode = rep(as.character(fib$ConditionCode), each = length(target_genes)),
  stringsAsFactors = FALSE
) |>
  mutate(
    ConditionCode = factor(ConditionCode, levels = condition_order),
    condition_label = factor(
      condition_labels[as.character(ConditionCode)],
      levels = unname(condition_labels)
    ),
    gene = factor(gene, levels = target_genes)
  )

target_dot_condition <- target_long |>
  group_by(ConditionCode, condition_label, gene) |>
  summarise(
    avg_expression = mean(expression),
    pct_expressed = 100 * mean(expression > detection_threshold),
    .groups = "drop"
  ) |>
  group_by(gene) |>
  mutate(scaled_expression = as.numeric(scale(avg_expression))) |>
  ungroup()

theme_pub <- theme_bw(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(7, 7, 7, 7)
  )

p_b <- ggplot(
  target_dot_condition,
  aes(gene, condition_label, size = pct_expressed, color = scaled_expression)
) +
  geom_point() +
  scale_size_area(max_size = 8, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_color_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0
  ) +
  labs(
    title = "Target transcripts by culture condition",
    subtitle = expression(italic(FAP) * " and the " * alpha[5] * beta[1] *
      " integrin subunits " * italic(ITGA5) * " and " * italic(ITGB1)),
    x = NULL, y = NULL,
    color = "Scaled mean", size = "Cells detected (%)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(face = "italic", size = 10.5))

# Figure 3C: original cleaned neutrophil subclusters.
original_umap <- fread(
  file.path(neut_dir, "source_data", "neutrophil_cell_metadata_annotations_and_umap.csv.gz"),
  select = c(
    "neutrophil_cluster_annotation", "neutrophil_UMAP_1", "neutrophil_UMAP_2"
  )
)
setnames(
  original_umap,
  c("neutrophil_cluster_annotation", "neutrophil_UMAP_1", "neutrophil_UMAP_2"),
  c("annotation", "UMAP_1", "UMAP_2")
)
original_umap <- original_umap[sample(.N)]

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

theme_umap <- theme_void(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    plot.margin = margin(7, 7, 7, 7)
  )

p_c <- ggplot(original_umap, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point_rast(size = 0.18, alpha = 0.55, raster.dpi = 400) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.8, alpha = 1), ncol = 2)) +
  labs(
    title = "Neutrophil-only scVI-Harmony subclusters",
    subtitle = "Original 92,166-cell cleaned UMAP",
    color = "Subcluster annotation"
  ) +
  theme_umap

# Figure 3D: canonical neutrophil states after low-RNA/unresolved removal.
canonical_umap <- fread(
  file.path(neut_dir, "source_data", "canonical_neutrophil_cell_metadata_and_umap.csv.gz"),
  select = c("neutrophil_state_annotation", "canonical_UMAP_1", "canonical_UMAP_2")
)
setnames(
  canonical_umap,
  c("neutrophil_state_annotation", "canonical_UMAP_1", "canonical_UMAP_2"),
  c("state", "UMAP_1", "UMAP_2")
)
canonical_umap <- canonical_umap[sample(.N)]

state_order <- c(
  "PADI4 neutrophil", "OSM neutrophil", "CXCR4 neutrophil", "MX1/ISG neutrophil"
)
state_colors <- c(
  "PADI4 neutrophil" = "#CC79A7",
  "OSM neutrophil" = "#D55E00",
  "CXCR4 neutrophil" = "#0072B2",
  "MX1/ISG neutrophil" = "#009E73"
)
canonical_umap$state <- factor(canonical_umap$state, levels = state_order)
state_counts <- table(canonical_umap$state)
state_labels <- setNames(
  paste0(names(state_counts), " (n = ", format(as.integer(state_counts), big.mark = ","), ")"),
  names(state_counts)
)

p_d <- ggplot(canonical_umap, aes(UMAP_1, UMAP_2, color = state)) +
  geom_point_rast(size = 0.25, alpha = 0.65, raster.dpi = 400) +
  scale_color_manual(values = state_colors, labels = state_labels, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1), ncol = 2)) +
  labs(
    title = "Canonical neutrophil states after reclustering",
    subtitle = "50,705 cells; Low RNA and unresolved states excluded",
    color = "Canonical neutrophil state"
  ) +
  theme_umap

figure_3 <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(0.88, 1.12), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 17))

stem <- file.path(
  figure_dir,
  "Figure_3_experimental_design_fibroblast_targets_neutrophil_UMAPs"
)
ggsave(
  paste0(stem, ".pdf"), figure_3, width = 20, height = 14,
  units = "in", device = cairo_pdf, bg = "white"
)
ggsave(
  paste0(stem, ".png"), figure_3, width = 20, height = 14,
  units = "in", dpi = 300, device = ragg::agg_png, bg = "white"
)
ggsave(
  paste0(stem, ".tiff"), figure_3, width = 20, height = 14,
  units = "in", dpi = 600, device = ragg::agg_tiff,
  compression = "lzw", bg = "white"
)

message("Figure 3 composite complete")
