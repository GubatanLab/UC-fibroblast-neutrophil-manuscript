#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrastr)
})

set.seed(20260816)

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture/Neutrophil_Subclustering"
input_file <- file.path(analysis_dir, "objects", "neutrophils_scVI_Harmony_subclustered_annotated_clean.rds")
object_dir <- file.path(analysis_dir, "objects")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading cleaned neutrophil object and removing Low RNA neutrophils")
neut <- readRDS(input_file)
state_order <- c(
  "PADI4 neutrophil", "OSM neutrophil",
  "CXCR4 neutrophil", "MX1/ISG neutrophil"
)
keep_cells <- rownames(neut@meta.data)[
  as.character(neut$neutrophil_state_annotation) %in% state_order
]
canonical <- subset(neut, cells = keep_cells)
canonical$neutrophil_state_annotation <- factor(
  as.character(canonical$neutrophil_state_annotation), levels = state_order
)
stopifnot(
  ncol(canonical) == 50705L,
  !any(as.character(canonical$neutrophil_state_annotation) == "Low RNA neutrophils"),
  !anyNA(canonical$neutrophil_state_annotation)
)

message("Rebuilding canonical-neutrophil graph and Leiden clusters")
canonical <- FindNeighbors(
  canonical, reduction = "scvi.harmony", dims = 1:30, k.param = 30,
  graph.name = c("canonical_scvi_harmony_nn", "canonical_scvi_harmony_snn"),
  verbose = TRUE
)
for (res in c(0.2, 0.4, 0.6)) {
  canonical <- FindClusters(
    canonical,
    graph.name = "canonical_scvi_harmony_snn",
    resolution = res,
    cluster.name = paste0("canonical_scVI_Harmony_res_", gsub("\\.", "_", res)),
    algorithm = 4,
    random.seed = 20260816,
    verbose = TRUE
  )
}
canonical$canonical_neutrophil_cluster <- factor(canonical$canonical_scVI_Harmony_res_0_4)

message("Recomputing UMAP after Low RNA removal")
canonical <- RunUMAP(
  canonical,
  reduction = "scvi.harmony", dims = 1:30,
  n.neighbors = 30, min.dist = 0.30, spread = 1,
  reduction.name = "umap.canonical.scvi.harmony",
  reduction.key = "CANONICALSCVIHARMONYUMAP_",
  seed.use = 20260816,
  verbose = TRUE
)

state_counts <- as.data.frame(table(
  canonical_state = canonical$neutrophil_state_annotation
))
state_counts$proportion <- state_counts$Freq / sum(state_counts$Freq)
fwrite(state_counts, file.path(source_dir, "canonical_neutrophil_state_counts.csv"))

cluster_counts <- as.data.frame(table(
  resolution_0_2 = canonical$canonical_scVI_Harmony_res_0_2,
  resolution_0_4 = canonical$canonical_scVI_Harmony_res_0_4,
  resolution_0_6 = canonical$canonical_scVI_Harmony_res_0_6
))
fwrite(cluster_counts, file.path(source_dir, "canonical_neutrophil_multiresolution_cluster_counts.csv"))

cluster_state <- as.data.frame(table(
  canonical_cluster = canonical$canonical_neutrophil_cluster,
  canonical_state = canonical$neutrophil_state_annotation
))
fwrite(cluster_state, file.path(source_dir, "canonical_neutrophil_cluster_state_overlap.csv"))

cluster_condition <- as.data.frame(table(
  canonical_cluster = canonical$canonical_neutrophil_cluster,
  condition = canonical$ConditionCode
))
fwrite(cluster_condition, file.path(source_dir, "canonical_neutrophil_cluster_condition_counts.csv"))

um <- as.data.frame(Embeddings(canonical, "umap.canonical.scvi.harmony"))
colnames(um) <- c("canonical_UMAP_1", "canonical_UMAP_2")
um$cell <- rownames(um)
md <- canonical@meta.data
md$cell <- rownames(md)
plot_df <- cbind(um, md[um$cell, , drop = FALSE])
plot_df <- plot_df[, !duplicated(colnames(plot_df)), drop = FALSE]
plot_df <- plot_df[sample(seq_len(nrow(plot_df))), , drop = FALSE]

state_colors <- c(
  "PADI4 neutrophil" = "#CC79A7",
  "OSM neutrophil" = "#D55E00",
  "CXCR4 neutrophil" = "#0072B2",
  "MX1/ISG neutrophil" = "#009E73"
)
legend_labels <- setNames(
  paste0(state_counts$canonical_state, " (n = ", format(state_counts$Freq, big.mark = ","), ")"),
  state_counts$canonical_state
)

p_umap <- ggplot(
  plot_df,
  aes(canonical_UMAP_1, canonical_UMAP_2, color = neutrophil_state_annotation)
) +
  geom_point_rast(size = 0.34, alpha = 0.66, raster.dpi = 400) +
  scale_color_manual(values = state_colors, labels = legend_labels, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(
    title = "Canonical neutrophil states after reclustering",
    subtitle = "50,705 neutrophils; Low RNA and unresolved states excluded before graph and UMAP reconstruction",
    color = "Canonical neutrophil state",
    caption = "Leiden graph and UMAP were recomputed from 30-dimensional scVI-Harmony embeddings; colors show canonical states only."
  ) +
  theme_void(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10.5, color = "#444444"),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_plot(p_umap, "Neutrophil_5_canonical_states_reclustered_scVI_Harmony_UMAP", 11.5, 8)

metadata_output <- canonical@meta.data
metadata_output$cell <- rownames(metadata_output)
coords <- as.data.frame(Embeddings(canonical, "umap.canonical.scvi.harmony"))
colnames(coords) <- c("canonical_UMAP_1", "canonical_UMAP_2")
coords$cell <- rownames(coords)
metadata_output <- merge(metadata_output, coords, by = "cell", sort = FALSE)
fwrite(metadata_output, file.path(source_dir, "canonical_neutrophil_cell_metadata_and_umap.csv.gz"))

saveRDS(
  canonical,
  file.path(object_dir, "canonical_neutrophils_scVI_Harmony_reclustered.rds"),
  compress = FALSE
)
message(
  "Canonical reclustering complete: ", ncol(canonical), " cells; ",
  nlevels(canonical$canonical_neutrophil_cluster), " primary clusters"
)
