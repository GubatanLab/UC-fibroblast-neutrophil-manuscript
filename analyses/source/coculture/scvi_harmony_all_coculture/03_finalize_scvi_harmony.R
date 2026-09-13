#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
  library(scales)
  library(data.table)
})

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture"
package_dir <- dirname(analysis_dir)
embedding_file <- file.path(analysis_dir, "source_data", "scVI_Harmony_latent_and_metadata.csv.gz")
figure_dir <- file.path(analysis_dir, "figures")
object_dir <- file.path(analysis_dir, "objects")
source_dir <- file.path(analysis_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading scVI-Harmony embeddings")
d <- fread(embedding_file, data.table = FALSE)
stopifnot(nrow(d) == 179571L)
stopifnot(length(unique(d$cell)) == nrow(d))

message("Creating a lightweight Seurat carrier for all pre-culture and co-culture cells")
d <- d[match(d$cell, d$cell), , drop = FALSE]
dummy <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(1, nrow(d)))
rownames(dummy) <- "SCVI_CARRIER"
colnames(dummy) <- d$cell
obj <- CreateSeuratObject(counts = dummy, meta.data = data.frame(row.names = d$cell))
for (nm in c("orig.ident", "SampleID", "DonorID", "Condition", "ConditionCode",
             "CultureStage", "dataset_source", "clean_cell_class", "clean_state", "clean_annotation")) {
  obj[[nm]] <- d[[nm]]
}

scvi_cols <- paste0("scVI_", 1:30)
harmony_cols <- paste0("scVI_Harmony_", 1:30)
scvi_embedding <- as.matrix(d[, scvi_cols, drop = FALSE])
harmony_embedding <- as.matrix(d[, harmony_cols, drop = FALSE])
rownames(scvi_embedding) <- rownames(harmony_embedding) <- d$cell
colnames(scvi_embedding) <- paste0("SCVI_", 1:30)
colnames(harmony_embedding) <- paste0("SCVIHARMONY_", 1:30)

obj[["scvi.retrained"]] <- CreateDimReducObject(scvi_embedding, assay = DefaultAssay(obj), key = "SCVI_")
obj[["scvi.harmony"]] <- CreateDimReducObject(harmony_embedding, assay = DefaultAssay(obj), key = "SCVIHARMONY_")

message("Building neighbors, clusters, and UMAP from scVI-Harmony embeddings")
obj <- FindNeighbors(
  obj, reduction = "scvi.harmony", dims = 1:30, k.param = 30,
  graph.name = c("scvi_harmony_nn", "scvi_harmony_snn"), verbose = TRUE
)
obj <- FindClusters(
  obj, graph.name = "scvi_harmony_snn", resolution = 0.2,
  cluster.name = "scVI_Harmony_res_0_2", algorithm = 4,
  random.seed = 20260816, verbose = TRUE
)
obj <- FindClusters(
  obj, graph.name = "scvi_harmony_snn", resolution = 0.4,
  cluster.name = "scVI_Harmony_res_0_4", algorithm = 4,
  random.seed = 20260816, verbose = TRUE
)
obj <- FindClusters(
  obj, graph.name = "scvi_harmony_snn", resolution = 0.6,
  cluster.name = "scVI_Harmony_res_0_6", algorithm = 4,
  random.seed = 20260816, verbose = TRUE
)
obj <- RunUMAP(
  obj, reduction = "scvi.harmony", dims = 1:30,
  n.neighbors = 30, min.dist = 0.30, spread = 1,
  reduction.name = "umap.scvi.harmony", reduction.key = "SCVIHARMONYUMAP_",
  seed.use = 20260816, verbose = TRUE
)
obj$scVI_Harmony_primary_cluster <- factor(obj$scVI_Harmony_res_0_4)
saveRDS(obj, file.path(object_dir, "all_coculture_scVI_Harmony_clustered_light.rds"), compress = FALSE)

umap_embedding <- Embeddings(obj, "umap.scvi.harmony")
d$UMAP_1 <- umap_embedding[d$cell, 1]
d$UMAP_2 <- umap_embedding[d$cell, 2]
d$cluster_res_0_2 <- as.character(obj@meta.data[d$cell, "scVI_Harmony_res_0_2"])
d$cluster_res_0_4 <- as.character(obj@meta.data[d$cell, "scVI_Harmony_res_0_4"])
d$cluster_res_0_6 <- as.character(obj@meta.data[d$cell, "scVI_Harmony_res_0_6"])

condition_order <- c("PRE", "N", "NNAMPT", "NA5", "CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(
  PRE = "Pre-culture neutrophils", N = "Neutrophils alone", NNAMPT = "Neutrophils + NAMPTi",
  NA5 = "Neutrophils + a5b1i", CF = "Control fibroblast co-culture",
  UF = "UC fibroblast co-culture", UN = "UC co-culture + NAMPTi",
  UA5 = "UC co-culture + a5b1i", UD = "UC co-culture + dual blockade"
)
d$condition_label <- factor(unname(condition_labels[d$ConditionCode]),
                            levels = unname(condition_labels[condition_order]))
d$cell_compartment <- factor(
  ifelse(d$clean_cell_class == "Neutrophil", "Neutrophil", "Fibroblast"),
  levels = c("Neutrophil", "Fibroblast")
)
d$cluster <- factor(d$cluster_res_0_4, levels = sort(unique(as.integer(d$cluster_res_0_4))))
fwrite(d, file.path(source_dir, "scVI_Harmony_cell_embeddings_clusters_and_metadata.csv.gz"))

cell_colors <- c("Neutrophil" = "#4C78A8", "Fibroblast" = "#D81B60")
condition_colors <- c(
  "Pre-culture neutrophils" = "#111827", "Neutrophils alone" = "#4D4D4D", "Neutrophils + NAMPTi" = "#56B4E9",
  "Neutrophils + a5b1i" = "#0072B2", "Control fibroblast co-culture" = "#009E73",
  "UC fibroblast co-culture" = "#D55E00", "UC co-culture + NAMPTi" = "#E69F00",
  "UC co-culture + a5b1i" = "#7B3294", "UC co-culture + dual blockade" = "#CC79A7"
)
cluster_colors <- setNames(hue_pal(l = 58, c = 95)(nlevels(d$cluster)), levels(d$cluster))

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

p_cluster <- ggplot(d, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point_rast(size = 0.18, alpha = 0.55, raster.dpi = 300) +
  scale_color_manual(values = cluster_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2)) +
  labs(
    title = "Leiden clusters from the scVI-Harmony neighbor graph",
    subtitle = "Primary clustering resolution = 0.4; 30 neighbors in the 30-dimensional scVI-Harmony space",
    color = "Cluster",
    caption = "Experimental condition was not supplied to scVI or Harmony."
  ) + theme_umap + theme(legend.position = "right")

p_condition <- ggplot(d, aes(UMAP_1, UMAP_2, color = condition_label)) +
  geom_point_rast(size = 0.18, alpha = 0.48, raster.dpi = 300) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Condition distribution in scVI-Harmony space",
    subtitle = "Nine experimental conditions share one donor-harmonized embedding",
    color = "Condition",
    caption = "Harmony used DonorID only; condition was preserved."
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
save_plot(p_cluster, "scVI_Harmony_UMAP_2_clusters", 11.2, 7.5)
save_plot(p_condition, "scVI_Harmony_UMAP_3_conditions", 11.5, 7.5)
save_plot(p_split, "scVI_Harmony_UMAP_4_condition_split", 12.5, 10.5)

write.csv(
  as.data.frame(table(Cluster = d$cluster, Cell_compartment = d$cell_compartment)),
  file.path(source_dir, "primary_cluster_cell_compartment_counts.csv"), row.names = FALSE
)
write.csv(
  as.data.frame(table(Condition = d$condition_label, Cluster = d$cluster)),
  file.path(source_dir, "condition_primary_cluster_counts.csv"), row.names = FALSE
)

message("Finalized scVI-Harmony figures and light Seurat object")
