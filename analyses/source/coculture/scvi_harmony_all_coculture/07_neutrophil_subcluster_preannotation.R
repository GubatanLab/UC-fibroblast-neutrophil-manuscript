#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(Seurat)
  library(Matrix)
  library(data.table)
})

set.seed(20260816)

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture"
input_h5ad <- file.path(analysis_dir, "objects", "all_coculture_scVI_Harmony_latent.h5ad")
out_dir <- file.path(analysis_dir, "Neutrophil_Subclustering")
object_dir <- file.path(out_dir, "objects")
source_dir <- file.path(out_dir, "source_data")
figure_dir <- file.path(out_dir, "figures")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading count matrix and scVI-Harmony embeddings")
sce <- read_h5ad(input_h5ad, as = "SingleCellExperiment")
keep <- as.character(colData(sce)$clean_cell_class) == "Neutrophil"
sce <- sce[, keep]
stopifnot(ncol(sce) == 179333L)

counts <- assay(sce, "counts")
md <- as.data.frame(colData(sce))
rownames(md) <- colnames(sce)
neut <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = md, min.cells = 0, min.features = 0)

scvi <- reducedDim(sce, "X_scVI")
harmony <- reducedDim(sce, "X_scVI_Harmony")
rownames(scvi) <- rownames(harmony) <- colnames(neut)
colnames(scvi) <- paste0("SCVI_", seq_len(ncol(scvi)))
colnames(harmony) <- paste0("SCVIHARMONY_", seq_len(ncol(harmony)))
neut[["scvi.retrained"]] <- CreateDimReducObject(scvi, assay = "RNA", key = "SCVI_")
neut[["scvi.harmony"]] <- CreateDimReducObject(harmony, assay = "RNA", key = "SCVIHARMONY_")
rm(sce, counts, scvi, harmony)
gc()

message("Normalizing RNA and constructing neutrophil-only graph")
neut <- NormalizeData(neut, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
neut <- FindNeighbors(
  neut, reduction = "scvi.harmony", dims = 1:30, k.param = 30,
  graph.name = c("neut_scvi_harmony_nn", "neut_scvi_harmony_snn"), verbose = TRUE
)
for (res in c(0.2, 0.4, 0.6)) {
  neut <- FindClusters(
    neut, graph.name = "neut_scvi_harmony_snn", resolution = res,
    cluster.name = paste0("neut_scVI_Harmony_res_", gsub("\\.", "_", res)),
    algorithm = 4, random.seed = 20260816, verbose = TRUE
  )
}
neut <- RunUMAP(
  neut, reduction = "scvi.harmony", dims = 1:30,
  n.neighbors = 30, min.dist = 0.30, spread = 1,
  reduction.name = "umap.neut.scvi.harmony", reduction.key = "NEUTSCVIHARMONYUMAP_",
  seed.use = 20260816, verbose = TRUE
)
neut$neutrophil_subcluster <- factor(neut$neut_scVI_Harmony_res_0_4)

message("Computing cluster diagnostics and downsampled marker discovery")
Idents(neut) <- "neutrophil_subcluster"
marker_genes <- intersect(c(
  "OSM", "CXCR4", "MX1", "PADI4", "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6",
  "FCGR3B", "CSF3R", "S100A8", "S100A9", "FPR1", "CXCL8", "IL1B", "TNF",
  "NAMPT", "MMP9", "OLFM4", "SELL", "ITGAM", "BCL2A1", "NFKBIA", "JUN", "FOS"
), rownames(neut))

avg <- AggregateExpression(neut, assays = "RNA", features = marker_genes, group.by = "neutrophil_subcluster", slot = "data", verbose = FALSE)$RNA
avg_df <- as.data.frame(as.matrix(avg))
avg_df$gene <- rownames(avg_df)
fwrite(avg_df, file.path(source_dir, "neutrophil_cluster_known_marker_average_expression.csv"))

markers <- FindAllMarkers(
  neut, assay = "RNA", slot = "data", only.pos = TRUE,
  min.pct = 0.10, logfc.threshold = 0.20,
  max.cells.per.ident = 2500, random.seed = 20260816,
  test.use = "wilcox", verbose = TRUE
)
markers <- markers[order(markers$cluster, markers$p_val_adj, -markers$avg_log2FC), , drop = FALSE]
fwrite(markers, file.path(source_dir, "neutrophil_subcluster_all_markers_downsampled.csv"))

write.csv(
  as.data.frame(table(Subcluster = neut$neutrophil_subcluster, Existing_state = neut$clean_state)),
  file.path(source_dir, "neutrophil_subcluster_existing_state_overlap.csv"), row.names = FALSE
)
write.csv(
  as.data.frame(table(Subcluster = neut$neutrophil_subcluster, Condition = neut$ConditionCode)),
  file.path(source_dir, "neutrophil_subcluster_condition_counts.csv"), row.names = FALSE
)
write.csv(
  as.data.frame(table(Condition = neut$ConditionCode, Donor = neut$DonorID)),
  file.path(source_dir, "neutrophil_condition_donor_counts.csv"), row.names = FALSE
)

saveRDS(neut, file.path(object_dir, "neutrophils_scVI_Harmony_subclustered_preannotation.rds"), compress = FALSE)
message("Preannotation complete: ", ncol(neut), " neutrophils and ", nlevels(neut$neutrophil_subcluster), " clusters")
