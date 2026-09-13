#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(ggrastr)
  library(data.table)
})

set.seed(20260816)

input_file <- "input_data/culture/outputs/019fbe7c_combined_decontX_reclustered_all_cells/Combined_neutrophil_fibroblast_decontX_reclustered.rds"
analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture/Fibroblast_Subclustering"
object_dir <- file.path(analysis_dir, "objects")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading and subsetting fibroblasts")
obj <- readRDS(input_file)
fib <- subset(obj, subset = clean_cell_class == "Fibroblast_candidate")
stopifnot(ncol(fib) == 238L)

# Remove provisional terminology while retaining the original confidence fields.
fib$clean_cell_class <- "Fibroblast"
fib$cell_compartment <- "Fibroblast"
fib$clean_annotation <- gsub("Fibroblast candidate", "Fibroblast", as.character(fib$clean_annotation), fixed = TRUE)
fib$clean_annotation <- gsub("fibroblast candidate", "fibroblast", fib$clean_annotation, fixed = TRUE)
for (nm in colnames(fib@meta.data)) {
  if (is.character(fib@meta.data[[nm]]) || is.factor(fib@meta.data[[nm]])) {
    z <- as.character(fib@meta.data[[nm]])
    z <- gsub("Fibroblast_candidate", "Fibroblast", z, fixed = TRUE)
    z <- gsub("Fibroblast candidate", "Fibroblast", z, fixed = TRUE)
    z <- gsub("fibroblast candidate", "fibroblast", z, fixed = TRUE)
    fib@meta.data[[nm]] <- z
  }
}

# Recompute a fibroblast-only representation from DecontX-corrected expression.
DefaultAssay(fib) <- "DECONTX"
fib@reductions <- list()
fib@graphs <- list()
fib@neighbors <- list()
fib <- NormalizeData(fib, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
fib <- FindVariableFeatures(fib, selection.method = "vst", nfeatures = 2500, verbose = FALSE)
fib <- ScaleData(fib, features = VariableFeatures(fib), verbose = FALSE)
fib <- RunPCA(fib, features = VariableFeatures(fib), npcs = 25, seed.use = 20260816, verbose = FALSE)
fib <- RunHarmony(
  fib, group.by.vars = "DonorID", reduction.use = "pca", dims.use = 1:20,
  reduction.save = "fibro.harmony",
  max_iter = 30, plot_convergence = FALSE, verbose = FALSE
)
fib <- FindNeighbors(
  fib, reduction = "fibro.harmony", dims = 1:20, k.param = 10,
  graph.name = c("fibro_nn", "fibro_snn"), verbose = FALSE
)
for (res in c(0.2, 0.4, 0.6)) {
  fib <- FindClusters(
    fib, graph.name = "fibro_snn", resolution = res,
    cluster.name = paste0("fibro_res_", gsub("\\.", "_", res)),
    algorithm = 4, random.seed = 20260816, verbose = FALSE
  )
}
fib <- RunUMAP(
  fib, reduction = "fibro.harmony", dims = 1:20, n.neighbors = 15,
  min.dist = 0.25, spread = 1, reduction.name = "umap.fibro.harmony",
  reduction.key = "FIBROUMAP_", seed.use = 20260816, verbose = FALSE
)
fib$fibroblast_subcluster <- factor(fib$fibro_res_0_4)

Idents(fib) <- "fibroblast_subcluster"
markers <- FindAllMarkers(
  fib, assay = "DECONTX", slot = "data", only.pos = TRUE,
  min.pct = 0.10, logfc.threshold = 0.20, test.use = "wilcox", verbose = FALSE
)
markers <- markers[order(markers$cluster, markers$p_val_adj, -markers$avg_log2FC), , drop = FALSE]
fwrite(markers, file.path(source_dir, "fibroblast_subcluster_all_markers.csv"))

annotation_map <- data.frame(
  fibroblast_subcluster = as.character(1:4),
  fibroblast_annotation = c(
    "COL1A1+ inflammatory ECM fibroblast",
    "CSF3R+ immune-like fibroblast",
    "CCN1/CCN2+ activated fibroblast",
    "S100A8/A9+ immune-like fibroblast"
  ),
  annotation_confidence = c("Moderate", "Low", "High", "Low"),
  defining_features = c(
    "COL1A1, COL1A2, COL3A1, DCN, LUM, FAP, IL6, CXCL8",
    "CSF3R, C5AR1, TYROBP, FCER1G; weak stromal program",
    "CCN1, CCN2, SPARC, CCL2, FOXF1, TFPI, ACTA2/TAGLN",
    "S100A8, S100A9, CLEC12A, SIGLEC14; weak stromal program"
  ),
  interpretation = c(
    "Matrix-producing fibroblasts with an inflammatory program",
    "Low-confidence immune-like fibroblast profile; possible mixed/ambient signal",
    "Coherent activated and matrix-remodeling stromal program",
    "Low-confidence myeloid-like fibroblast profile; possible mixed/ambient signal"
  ),
  stringsAsFactors = FALSE
)
write.csv(annotation_map, file.path(source_dir, "fibroblast_subcluster_annotation_key.csv"), row.names = FALSE)
fib$fibroblast_annotation <- annotation_map$fibroblast_annotation[
  match(as.character(fib$fibroblast_subcluster), annotation_map$fibroblast_subcluster)
]
fib$fibroblast_annotation <- factor(fib$fibroblast_annotation, levels = annotation_map$fibroblast_annotation)
fib$fibroblast_annotation_confidence <- annotation_map$annotation_confidence[
  match(as.character(fib$fibroblast_subcluster), annotation_map$fibroblast_subcluster)
]
write.csv(
  as.data.frame(table(Subcluster = fib$fibroblast_subcluster, Condition = fib$ConditionCode)),
  file.path(source_dir, "fibroblast_subcluster_condition_counts.csv"), row.names = FALSE
)
write.csv(
  as.data.frame(table(Subcluster = fib$fibroblast_subcluster, Confidence = fib$stromal_confidence)),
  file.path(source_dir, "fibroblast_subcluster_confidence_counts.csv"), row.names = FALSE
)

um <- as.data.frame(Embeddings(fib, "umap.fibro.harmony"))
colnames(um) <- c("UMAP_1", "UMAP_2")
um$cell <- rownames(um)
plot_df <- cbind(um, fib@meta.data[um$cell, , drop = FALSE])
centroids <- aggregate(cbind(UMAP_1, UMAP_2) ~ fibroblast_subcluster, data = plot_df, FUN = median)

annotation_colors <- c("#0072B2", "#999999", "#D55E00", "#CC79A7")
names(annotation_colors) <- annotation_map$fibroblast_annotation
condition_colors <- c(CF = "#009E73", UF = "#D55E00", UN = "#E69F00", UA5 = "#7B3294", UD = "#CC79A7")
theme_umap <- theme_void(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0),
    plot.subtitle = element_text(size = 11.5, color = "#333333"),
    legend.title = element_text(face = "bold"), legend.text = element_text(size = 9.5),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0),
    plot.margin = margin(10, 12, 10, 12)
  )

p_annotation <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = fibroblast_annotation)) +
  geom_point(size = 1.25, alpha = 0.82) +
  geom_label(
    data = centroids, aes(UMAP_1, UMAP_2, label = fibroblast_subcluster), inherit.aes = FALSE,
    size = 4, fontface = "bold", label.size = 0.18, fill = "white"
  ) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Fibroblast subclusters and marker-based annotations",
    subtitle = "DecontX expression, fibroblast-only PCA, donor Harmony, and Leiden resolution 0.4",
    color = "Fibroblast annotation",
    caption = "Immune-like annotations are retained with low confidence because their stromal programs are weak."
  ) + theme_umap + theme(legend.position = "right")

p_condition <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = ConditionCode)) +
  geom_point(size = 1.25, alpha = 0.82) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Fibroblast subclusters by experimental condition",
    subtitle = "All 238 fibroblasts from the five fibroblast-containing conditions",
    color = "Condition"
  ) + theme_umap + theme(legend.position = "right")

dot_features <- c(
  "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "IL6", "CXCL8",
  "CCN1", "CCN2", "SPARC", "CCL2", "FOXF1", "TFPI", "ACTA2", "TAGLN",
  "CSF3R", "C5AR1", "TYROBP", "S100A8", "S100A9", "CLEC12A", "SIGLEC14"
)
Idents(fib) <- "fibroblast_annotation"
p_dot <- DotPlot(
  fib, assay = "DECONTX", features = dot_features, dot.scale = 7,
  cols = c("#2166AC", "#B2182B")
) +
  labs(
    title = "Marker programs supporting fibroblast annotations",
    x = NULL, y = NULL, color = "Scaled expression", size = "Percent expressed"
  ) +
  theme_bw(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8.5),
    axis.text.y = element_text(size = 9.5),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.25),
    panel.grid.minor = element_blank()
  )

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_plot(p_annotation, "Fibroblast_UMAP_1_annotations", 11.5, 7.5)
save_plot(p_condition, "Fibroblast_UMAP_2_conditions", 10.5, 7.5)
save_plot(p_dot, "Fibroblast_3_annotation_marker_dotplot", 14.5, 6.3)

metadata_output <- fib@meta.data
metadata_output$cell <- rownames(metadata_output)
coords_output <- as.data.frame(Embeddings(fib, "umap.fibro.harmony"))
colnames(coords_output) <- c("fibro_UMAP_1", "fibro_UMAP_2")
coords_output$cell <- rownames(coords_output)
metadata_output <- merge(metadata_output, coords_output, by = "cell", sort = FALSE)
fwrite(metadata_output, file.path(source_dir, "fibroblast_cell_metadata_annotations_and_umap.csv.gz"))

saveRDS(fib, file.path(object_dir, "fibroblasts_subclustered_annotated.rds"), compress = FALSE)
message("Fibroblast subclustering and annotation complete: ", ncol(fib), " cells and ", nlevels(fib$fibroblast_subcluster), " primary clusters")
