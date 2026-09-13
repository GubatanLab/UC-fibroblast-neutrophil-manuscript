#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
  library(scales)
})

input_file <- "input_data/culture/outputs/019fbe7c_nc1_integration/precoculture_comparison/Neutrophil_only_UMI100_NC1_precoculture_Harmony_integrated.rds"
package_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
figure_dir <- file.path(package_dir, "figures", "Preculture_Harmony")
source_dir <- file.path(package_dir, "source_data", "Preculture_Harmony")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_file)
stopifnot("umap.harmony" %in% names(obj@reductions))

d <- as.data.frame(Embeddings(obj, "umap.harmony"))
colnames(d)[1:2] <- c("UMAP_1", "UMAP_2")
d$cell <- rownames(d)
md <- obj@meta.data[d$cell, , drop = FALSE]
d$Condition <- as.character(md$Condition)
d$DonorID <- as.character(md$DonorID)
d$SampleID <- as.character(md$SampleID)
d$CultureStage <- as.character(md$CultureStage)
d$dataset_source <- as.character(md$dataset_source)
d$neutrophil_state <- as.character(md$neutrophil_subset_short)

condition_order <- c(
  "Pre-coculture",
  "Neutrophil_Only",
  "Neutrophil_Only_NAMPTinhibitor",
  "Neutrophil_Only+ a5B1inhibitor",
  "Neutrophil_ControlFibroblasts",
  "Neutrophil_UCFibroblasts",
  "Neutrophil_UCFibroblasts_NAMPTinhibitor",
  "Neutrophil_UCFibroblasts_a5B1inhibitor",
  "Neutrophil_UCFibroblasts_NAMPTa5B1inhibitors"
)
condition_labels <- c(
  "Pre-coculture" = "Pre-culture neutrophils",
  "Neutrophil_Only" = "Neutrophils alone",
  "Neutrophil_Only_NAMPTinhibitor" = "Neutrophils + NAMPTi",
  "Neutrophil_Only+ a5B1inhibitor" = "Neutrophils + a5b1i",
  "Neutrophil_ControlFibroblasts" = "Control fibroblast co-culture",
  "Neutrophil_UCFibroblasts" = "UC fibroblast co-culture",
  "Neutrophil_UCFibroblasts_NAMPTinhibitor" = "UC co-culture + NAMPTi",
  "Neutrophil_UCFibroblasts_a5B1inhibitor" = "UC co-culture + a5b1i",
  "Neutrophil_UCFibroblasts_NAMPTa5B1inhibitors" = "UC co-culture + dual blockade"
)
d$condition_label <- factor(unname(condition_labels[d$Condition]),
                            levels = unname(condition_labels[condition_order]))
d$culture_group <- factor(
  ifelse(d$Condition == "Pre-coculture", "Pre-culture neutrophil", "Cultured neutrophil"),
  levels = c("Pre-culture neutrophil", "Cultured neutrophil")
)

condition_colors <- c(
  "Pre-culture neutrophils" = "#111827",
  "Neutrophils alone" = "#6B7280",
  "Neutrophils + NAMPTi" = "#9CA3AF",
  "Neutrophils + a5b1i" = "#4B5563",
  "Control fibroblast co-culture" = "#2563EB",
  "UC fibroblast co-culture" = "#DC2626",
  "UC co-culture + NAMPTi" = "#F59E0B",
  "UC co-culture + a5b1i" = "#7C3AED",
  "UC co-culture + dual blockade" = "#059669"
)
group_colors <- c("Pre-culture neutrophil" = "#111827", "Cultured neutrophil" = "#4C78A8")
state_colors <- c(
  "LowRNA" = "#9CA3AF", "CXCR4" = "#2563EB", "OSM" = "#DC2626",
  "PADI4" = "#7C3AED", "MX1" = "#059669"
)

theme_umap <- theme_void(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0),
    plot.subtitle = element_text(size = 11.5, color = "#333333"),
    strip.text = element_text(face = "bold", size = 10.2, margin = margin(5, 2, 5, 2)),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9.5),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0),
    plot.margin = margin(10, 12, 10, 12)
  )

pre <- d[d$culture_group == "Pre-culture neutrophil", , drop = FALSE]
post <- d[d$culture_group == "Cultured neutrophil", , drop = FALSE]

p_stage <- ggplot() +
  geom_point_rast(data = post, aes(UMAP_1, UMAP_2, color = culture_group),
                  size = 0.18, alpha = 0.26, raster.dpi = 300) +
  geom_point_rast(data = pre, aes(UMAP_1, UMAP_2, color = culture_group),
                  size = 0.20, alpha = 0.56, raster.dpi = 300) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Pre-culture and cultured neutrophils in joint Harmony space",
    subtitle = "Sparse PCA followed by 30-dimensional Harmony integration and Harmony UMAP",
    color = "Culture stage",
    caption = paste0("n = ", comma(nrow(d)), " neutrophils: ", comma(nrow(pre)),
                     " pre-culture and ", comma(nrow(post)), " cultured")
  ) + theme_umap + theme(legend.position = "right")

background <- d[, c("UMAP_1", "UMAP_2")]
p_condition <- ggplot() +
  geom_point_rast(data = background, aes(UMAP_1, UMAP_2), color = "#E4E4E4",
                  size = 0.06, alpha = 0.10, raster.dpi = 300) +
  geom_point_rast(data = d, aes(UMAP_1, UMAP_2, color = condition_label),
                  size = 0.18, alpha = 0.62, raster.dpi = 300) +
  facet_wrap(~condition_label, ncol = 3) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  guides(color = "none") +
  labs(
    title = "Condition occupancy with a pre-culture reference",
    subtitle = "The same joint Harmony embedding and axis limits are used in every panel",
    caption = "Light gray points show the complete integrated neutrophil atlas."
  ) + theme_umap + theme(panel.spacing = unit(7, "pt"))

p_state <- ggplot() +
  geom_point_rast(data = d, aes(UMAP_1, UMAP_2, color = neutrophil_state),
                  size = 0.18, alpha = 0.58, raster.dpi = 300) +
  facet_wrap(~culture_group, ncol = 2) +
  scale_color_manual(values = state_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 1)) +
  labs(
    title = "Neutrophil states before and after culture",
    subtitle = "Established cultured-cell states were preserved; pre-culture states were transferred in Harmony space",
    color = "Neutrophil state",
    caption = "Pre-culture state transfer used a weighted 50-nearest-neighbor vote."
  ) + theme_umap + theme(legend.position = "bottom", panel.spacing = unit(9, "pt"))

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}

save_plot(p_stage, "Preculture_UMAP_1_Joint_Harmony_stage", 11, 7.6)
save_plot(p_condition, "Preculture_UMAP_2_Joint_Harmony_by_condition", 13.2, 11.2)
save_plot(p_state, "Preculture_UMAP_3_Joint_Harmony_by_state", 12.5, 6.9)

write.csv(d, file.path(source_dir, "Preculture_joint_Harmony_UMAP_coordinates.csv"), row.names = FALSE)
write.csv(
  as.data.frame(table(Condition = d$condition_label, State = d$neutrophil_state)),
  file.path(source_dir, "Preculture_joint_Harmony_condition_state_counts.csv"), row.names = FALSE
)

cat("Generated pre-culture Harmony UMAPs in", figure_dir, "\n")
cat("Cells:", nrow(d), " Pre-culture:", nrow(pre), " Cultured:", nrow(post), "\n")
