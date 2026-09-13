#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
  library(scales)
  library(data.table)
})

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture"
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
object_file <- file.path(analysis_dir, "objects", "all_coculture_scVI_Harmony_clustered_light.rds")
cell_file <- file.path(source_dir, "scVI_Harmony_cell_embeddings_clusters_and_metadata.csv.gz")
latent_file <- file.path(source_dir, "scVI_Harmony_latent_and_metadata.csv.gz")

replace_candidate <- function(z) {
  z <- as.character(z)
  z <- gsub("Fibroblast_candidate", "Fibroblast", z, fixed = TRUE)
  z <- gsub("Fibroblast candidate", "Fibroblast", z, fixed = TRUE)
  gsub("fibroblast candidate", "fibroblast", z, fixed = TRUE)
}

message("Relabeling atlas metadata and source tables")
obj <- readRDS(object_file)
for (nm in colnames(obj@meta.data)) {
  if (is.character(obj@meta.data[[nm]]) || is.factor(obj@meta.data[[nm]])) {
    obj@meta.data[[nm]] <- replace_candidate(obj@meta.data[[nm]])
  }
}
obj$cell_compartment <- ifelse(obj$clean_cell_class == "Neutrophil", "Neutrophil", "Fibroblast")
saveRDS(obj, object_file, compress = FALSE)

d <- fread(cell_file, data.table = FALSE)
for (nm in names(d)) {
  if (is.character(d[[nm]]) || is.factor(d[[nm]])) d[[nm]] <- replace_candidate(d[[nm]])
}
d$clean_cell_class[d$clean_cell_class != "Neutrophil"] <- "Fibroblast"
d$cell_compartment <- factor(
  ifelse(d$clean_cell_class == "Neutrophil", "Neutrophil", "Fibroblast"),
  levels = c("Neutrophil", "Fibroblast")
)
fwrite(d, cell_file)

latent <- fread(latent_file, data.table = FALSE)
for (nm in names(latent)) {
  if (is.character(latent[[nm]]) || is.factor(latent[[nm]])) latent[[nm]] <- replace_candidate(latent[[nm]])
}
latent$clean_cell_class[latent$clean_cell_class != "Neutrophil"] <- "Fibroblast"
fwrite(latent, latent_file)

condition_order <- c("PRE", "N", "NNAMPT", "NA5", "CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(
  PRE = "Pre-culture neutrophils", N = "Neutrophils alone", NNAMPT = "Neutrophils + NAMPTi",
  NA5 = "Neutrophils + a5b1i", CF = "Control fibroblast co-culture",
  UF = "UC fibroblast co-culture", UN = "UC co-culture + NAMPTi",
  UA5 = "UC co-culture + a5b1i", UD = "UC co-culture + dual blockade"
)
d$condition_label <- factor(unname(condition_labels[d$ConditionCode]),
                            levels = unname(condition_labels[condition_order]))
cell_colors <- c("Neutrophil" = "#4C78A8", "Fibroblast" = "#D81B60")
condition_colors <- c(
  "Pre-culture neutrophils" = "#111827", "Neutrophils alone" = "#4D4D4D",
  "Neutrophils + NAMPTi" = "#56B4E9", "Neutrophils + a5b1i" = "#0072B2",
  "Control fibroblast co-culture" = "#009E73", "UC fibroblast co-culture" = "#D55E00",
  "UC co-culture + NAMPTi" = "#E69F00", "UC co-culture + a5b1i" = "#7B3294",
  "UC co-culture + dual blockade" = "#CC79A7"
)
theme_umap <- theme_void(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0),
    plot.subtitle = element_text(size = 11.5, color = "#333333"),
    strip.text = element_text(face = "bold", size = 10.2, margin = margin(5, 2, 5, 2)),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    legend.title = element_text(face = "bold"), legend.text = element_text(size = 9.5),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0),
    plot.margin = margin(10, 12, 10, 12)
  )
neut <- d[d$cell_compartment == "Neutrophil", , drop = FALSE]
fib <- d[d$cell_compartment == "Fibroblast", , drop = FALSE]

p_class <- ggplot() +
  geom_point_rast(data = neut, aes(UMAP_1, UMAP_2, color = cell_compartment),
                  size = 0.18, alpha = 0.40, raster.dpi = 300) +
  geom_point(data = fib, aes(UMAP_1, UMAP_2, color = cell_compartment),
             size = 1.05, alpha = 0.95, shape = 16) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Pre-culture neutrophils and all co-culture cells in scVI-Harmony space",
    subtitle = "scVI latent representation followed by donor-only Harmony, graph clustering, and UMAP",
    color = "Cell compartment",
    caption = paste0("n = ", comma(nrow(d)), " cells; fibroblasts are enlarged for visibility")
  ) + theme_umap + theme(legend.position = "right")

background <- d[, c("UMAP_1", "UMAP_2")]
p_split <- ggplot() +
  geom_point_rast(data = background, aes(UMAP_1, UMAP_2), color = "#E4E4E4",
                  size = 0.06, alpha = 0.10, raster.dpi = 300) +
  geom_point_rast(data = d, aes(UMAP_1, UMAP_2, color = condition_label),
                  size = 0.16, alpha = 0.62, raster.dpi = 300) +
  geom_point(data = fib, aes(UMAP_1, UMAP_2), color = "#D81B60", size = 0.80, alpha = 0.90) +
  facet_wrap(~condition_label, ncol = 3) +
  scale_color_manual(values = condition_colors, drop = FALSE) + guides(color = "none") +
  labs(
    title = "Condition-specific occupancy of the scVI-Harmony atlas",
    subtitle = "Fibroblasts are overlaid in magenta where present",
    caption = "Every panel uses identical UMAP coordinates and axis limits."
  ) + theme_umap + theme(panel.spacing = unit(7, "pt"))

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_plot(p_class, "scVI_Harmony_UMAP_1_cell_compartment", 10.8, 7.5)
save_plot(p_split, "scVI_Harmony_UMAP_4_condition_split", 12.5, 10.5)
write.csv(
  as.data.frame(table(Cluster = d$cluster, Cell_compartment = d$cell_compartment)),
  file.path(source_dir, "primary_cluster_cell_compartment_counts.csv"), row.names = FALSE
)
message("Atlas relabeling and figure refresh complete")

