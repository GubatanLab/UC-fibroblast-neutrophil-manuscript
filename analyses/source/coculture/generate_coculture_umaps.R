#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
  library(scales)
})

input_file <- "input_data/culture/outputs/019fbe7c_combined_decontX_reclustered_all_cells/Combined_neutrophil_fibroblast_decontX_reclustered_light.rds"
package_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
figure_dir <- file.path(package_dir, "figures", "UMAPs")
source_dir <- file.path(package_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_file)
reduction_name <- "umap.decontX"
stopifnot(reduction_name %in% names(obj@reductions))

coords <- as.data.frame(Embeddings(obj, reduction_name))
colnames(coords)[1:2] <- c("UMAP_1", "UMAP_2")
coords$cell <- rownames(coords)
md <- obj@meta.data[coords$cell, , drop = FALSE]
coords$ConditionCode <- as.character(md$ConditionCode)
coords$Condition <- as.character(md$Condition)
coords$DonorID <- as.character(md$DonorID)
coords$clean_annotation <- as.character(md$clean_annotation)
coords$clean_cell_class <- as.character(md$clean_cell_class)

condition_order <- c("N", "NNAMPT", "NA5", "CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(
  "N" = "Neutrophils alone",
  "NNAMPT" = "Neutrophils + NAMPTi",
  "NA5" = "Neutrophils + a5b1i",
  "CF" = "Control fibroblast co-culture",
  "UF" = "UC fibroblast co-culture",
  "UN" = "UC co-culture + NAMPTi",
  "UA5" = "UC co-culture + a5b1i",
  "UD" = "UC co-culture + dual blockade"
)
missing_codes <- setdiff(unique(coords$ConditionCode), names(condition_labels))
if (length(missing_codes)) stop("Unmapped condition codes: ", paste(missing_codes, collapse = ", "))
coords$condition_label <- factor(
  unname(condition_labels[coords$ConditionCode]),
  levels = unname(condition_labels[condition_order])
)

condition_colors <- c(
  "Neutrophils alone" = "#4D4D4D",
  "Neutrophils + NAMPTi" = "#56B4E9",
  "Neutrophils + a5b1i" = "#0072B2",
  "Control fibroblast co-culture" = "#009E73",
  "UC fibroblast co-culture" = "#D55E00",
  "UC co-culture + NAMPTi" = "#E69F00",
  "UC co-culture + a5b1i" = "#7B3294",
  "UC co-culture + dual blockade" = "#CC79A7"
)

annotation_colors <- c(
  "CXCR4 neutrophil" = "#0072B2",
  "OSM neutrophil" = "#D55E00",
  "PADI4 neutrophil" = "#E69F00",
  "MX1/ISG-like neutrophil" = "#7B3294",
  "Low-signal/unresolved neutrophil" = "#777777"
)
coords_neut <- coords[coords$clean_cell_class == "Neutrophil", , drop = FALSE]

theme_umap <- theme_void(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0),
    plot.subtitle = element_text(size = 11.5, color = "#333333"),
    strip.text = element_text(face = "bold", size = 10.5, margin = margin(5, 2, 5, 2)),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9.5),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0),
    plot.margin = margin(10, 12, 10, 12)
  )

axis_arrows <- function(d) {
  xr <- range(d$UMAP_1); yr <- range(d$UMAP_2)
  x0 <- xr[1] + 0.04 * diff(xr); y0 <- yr[1] + 0.04 * diff(yr)
  list(
    annotate("segment", x = x0, xend = x0 + 0.13 * diff(xr), y = y0, yend = y0,
             linewidth = 0.55, arrow = arrow(length = unit(0.10, "in"))),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + 0.13 * diff(yr),
             linewidth = 0.55, arrow = arrow(length = unit(0.10, "in"))),
    annotate("text", x = x0 + 0.065 * diff(xr), y = y0 - 0.035 * diff(yr), label = "UMAP 1", size = 3.1),
    annotate("text", x = x0 - 0.025 * diff(xr), y = y0 + 0.065 * diff(yr), label = "UMAP 2", angle = 90, size = 3.1)
  )
}

p_overview <- ggplot(coords, aes(UMAP_1, UMAP_2, color = condition_label)) +
  geom_point_rast(size = 0.24, alpha = 0.42, stroke = 0, raster.dpi = 300) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  labs(
    title = "Shared UMAP of neutrophil-fibroblast co-culture conditions",
    subtitle = "All conditions projected in one DecontX-corrected, Harmony-integrated embedding",
    color = "Condition",
    caption = paste0("n = ", comma(nrow(coords)), " cells; each point is one cell")
  ) + theme_umap + theme(legend.position = "right")
for (layer in axis_arrows(coords)) p_overview <- p_overview + layer

background <- coords[, c("UMAP_1", "UMAP_2")]
p_split_state <- ggplot() +
  geom_point_rast(data = background, aes(UMAP_1, UMAP_2), color = "#E4E4E4", size = 0.08,
                  alpha = 0.14, raster.dpi = 300) +
  geom_point_rast(data = coords_neut, aes(UMAP_1, UMAP_2, color = clean_annotation), size = 0.20,
                  alpha = 0.68, raster.dpi = 300) +
  facet_wrap(~condition_label, ncol = 4) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 1)) +
  labs(
    title = "Neutrophil-state occupancy across co-culture conditions",
    subtitle = "Colored cells indicate the named condition; the shared atlas is shown in light gray",
    color = "Cell state",
    caption = "The same UMAP coordinates and axis limits are used in every panel."
  ) + theme_umap +
  theme(legend.position = "bottom", panel.spacing = unit(7, "pt"))

coords$condition_group <- ifelse(coords$ConditionCode %in% c("N", "CF", "UF"),
                                 "Fibroblast co-culture effect", "Blockade effect")
coords$condition_group <- factor(coords$condition_group,
                                 levels = c("Fibroblast co-culture effect", "Blockade effect"))
p_split_condition <- ggplot() +
  geom_point_rast(data = background, aes(UMAP_1, UMAP_2), color = "#E5E5E5", size = 0.08,
                  alpha = 0.12, raster.dpi = 300) +
  geom_point_rast(data = coords, aes(UMAP_1, UMAP_2, color = condition_label), size = 0.19,
                  alpha = 0.65, raster.dpi = 300) +
  facet_wrap(~condition_label, ncol = 4) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  guides(color = "none") +
  labs(
    title = "Condition-specific distributions in the shared co-culture atlas",
    subtitle = "Untreated, fibroblast-exposed, NAMPT-inhibited, a5b1-inhibited, and dual-blockade cells",
    caption = "Light gray points provide the common embedding reference."
  ) + theme_umap + theme(panel.spacing = unit(7, "pt"))

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}

save_plot(p_overview, "UMAP_1_All_conditions_overlay", 11.5, 7.5)
save_plot(p_split_condition, "UMAP_2_Condition_split", 14, 7.8)
save_plot(p_split_state, "UMAP_3_Condition_split_by_neutrophil_state", 14, 8.4)

write.csv(coords, file.path(source_dir, "Coculture_UMAP_coordinates_and_metadata.csv"), row.names = FALSE)
write.csv(
  as.data.frame(table(Condition = coords$condition_label, Annotation = coords$clean_annotation)),
  file.path(source_dir, "Coculture_UMAP_condition_annotation_counts.csv"), row.names = FALSE
)

cat("Generated three UMAP figure sets in", figure_dir, "\n")
cat("Cells:", nrow(coords), " Conditions:", nlevels(coords$condition_label), "\n")
