#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(42)

local_library <- normalizePath("R_libs", winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "match"
if (!mode %in% c("extract_reference", "match", "neutrophil_umap")) {
  stop("Mode must be extract_reference, match, or neutrophil_umap")
}

reference_file <- "UCCODEX1_Annotated_withNeutrophil0.8Annotations.rds"
input_file <- "UCCODEX1_Annotated_repaired_withUMAP_neutrophilSubclusters.rds"
output_file <- "UCCODEX1_Annotated_repaired_withUMAP_neutrophil0.8Matched.rds"
output_dir <- "codex_repair_P02_P07/neutrophil_0.8_matched"
reference_meta_file <- file.path(output_dir, "reference_neutrophil_0.8_metadata.rds")
neutrophil_file <- file.path(output_dir, "neutrophil_0.8_matched_umap.rds")
figure_dir <- "giotto_codex_results/figures"
table_dir <- "giotto_codex_results/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

reference_fields <- c(
  "PatientID", "Celltype", "Neutrophil_subtype_0.8", "Neutrophil_annotation_0.8",
  "Neutrophil_cluster_0.8", "Celltype_original", "Celltype_updated"
)

majority_vote <- function(labels) {
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!length(labels)) return(NA_character_)
  names(sort(table(labels), decreasing = TRUE))[1]
}

if (mode == "extract_reference") {
  message("Extracting reference 0.8 annotation metadata")
  reference <- readRDS(reference_file)
  stopifnot(all(reference_fields %in% colnames(reference[[]])))
  reference_meta <- as.data.table(reference[[]], keep.rownames = "cell_ID")[, ..reference_fields]
  reference_meta[, cell_ID := colnames(reference)]
  setcolorder(reference_meta, c("cell_ID", reference_fields))
  saveRDS(reference_meta, reference_meta_file, compress = FALSE)
  message("Saved reference metadata for ", nrow(reference_meta), " cells")
}

if (mode == "match") {
  message("Matching reference annotations and transferring labels to repaired P07")
  obj <- readRDS(input_file)
  reference_meta <- readRDS(reference_meta_file)
  setkey(reference_meta, cell_ID)

  cells <- colnames(obj)
  p07_cells <- cells[obj$PatientID == "P07"]
  direct_cells <- cells[obj$PatientID != "P07" & cells %chin% reference_meta$cell_ID]
  if (length(direct_cells) != sum(obj$PatientID != "P07")) {
    stop("Not all non-P07 cells have direct reference matches")
  }

  direct_meta <- reference_meta[J(direct_cells)]
  harmony_embedding <- Embeddings(obj, "harmony")[, 1:10, drop = FALSE]
  train_labels <- as.character(direct_meta$Celltype_updated)
  neighbors <- RANN::nn2(
    data = harmony_embedding[direct_cells, , drop = FALSE],
    query = harmony_embedding[p07_cells, , drop = FALSE],
    k = 25
  )
  neighbor_labels <- matrix(
    train_labels[neighbors$nn.idx],
    nrow = nrow(neighbors$nn.idx),
    ncol = ncol(neighbors$nn.idx)
  )
  p07_predictions <- apply(neighbor_labels, 1, majority_vote)
  p07_confidence <- apply(neighbor_labels, 1, function(x) {
    tab <- table(x)
    max(tab) / sum(tab)
  })

  celltype_original <- as.character(obj$Celltype)
  celltype_updated <- celltype_original
  names(celltype_updated) <- cells
  celltype_updated[direct_cells] <- as.character(direct_meta$Celltype_updated)
  celltype_updated[p07_cells] <- p07_predictions

  annotation <- setNames(rep(NA_character_, length(cells)), cells)
  subtype <- annotation
  cluster_08 <- annotation
  annotation[direct_cells] <- as.character(direct_meta$Neutrophil_annotation_0.8)
  subtype[direct_cells] <- as.character(direct_meta$Neutrophil_subtype_0.8)
  cluster_08[direct_cells] <- as.character(direct_meta$Neutrophil_cluster_0.8)

  p07_is_subtyped <- grepl("^Neutrophil_", p07_predictions)
  p07_annotation <- rep(NA_character_, length(p07_cells))
  p07_annotation[p07_is_subtyped] <- sub("^Neutrophil_", "", p07_predictions[p07_is_subtyped])
  annotation[p07_cells] <- p07_annotation
  subtype[p07_cells] <- p07_annotation

  train_clusters <- as.character(direct_meta$Neutrophil_cluster_0.8)
  neighbor_clusters <- matrix(
    train_clusters[neighbors$nn.idx],
    nrow = nrow(neighbors$nn.idx),
    ncol = ncol(neighbors$nn.idx)
  )
  p07_clusters <- vapply(seq_along(p07_cells), function(i) {
    if (!p07_is_subtyped[i]) return(NA_character_)
    same_annotation <- neighbor_labels[i, ] == p07_predictions[i]
    majority_vote(neighbor_clusters[i, same_annotation])
  }, character(1))
  cluster_08[p07_cells] <- p07_clusters

  if (!"Neutrophil_de_novo_subcluster" %in% colnames(obj[[]])) {
    obj$Neutrophil_de_novo_subcluster <- obj$Neutrophil_subcluster
  }
  if (!"Neutrophil_de_novo_cluster" %in% colnames(obj[[]])) {
    obj$Neutrophil_de_novo_cluster <- obj$Neutrophil_cluster
  }
  obj$Celltype_original <- celltype_original
  obj$Celltype_updated <- unname(celltype_updated[cells])
  obj$Celltype_refined <- obj$Celltype_updated
  obj$Neutrophil_annotation_0.8 <- unname(annotation[cells])
  obj$Neutrophil_subtype_0.8 <- unname(subtype[cells])
  obj$Neutrophil_cluster_0.8 <- unname(cluster_08[cells])
  obj$Neutrophil_subset_annotation <- ifelse(
    grepl("^Neutrophil", obj$Celltype_updated),
    obj$Celltype_updated,
    NA_character_
  )
  obj$Neutrophil_subcluster <- obj$Neutrophil_subset_annotation
  obj$Neutrophil_cluster <- obj$Neutrophil_cluster_0.8

  direct_check <- reference_meta[J(direct_cells)]
  stopifnot(identical(
    unname(obj$Celltype_updated[match(direct_cells, cells)]),
    unname(as.character(direct_check$Celltype_updated))
  ))

  transfer_qc <- data.table(
    cell_ID = p07_cells,
    predicted_Celltype_updated = p07_predictions,
    predicted_Neutrophil_annotation_0.8 = p07_annotation,
    predicted_Neutrophil_cluster_0.8 = p07_clusters,
    vote_confidence = p07_confidence,
    mean_neighbor_distance = rowMeans(neighbors$nn.dists)
  )
  fwrite(transfer_qc, file.path(table_dir, "P07_neutrophil_0.8_annotation_transfer.csv"))
  fwrite(
    data.table(Celltype_updated = obj$Celltype_updated)[, .N, by = Celltype_updated][order(-N)],
    file.path(table_dir, "celltype_updated_counts_neutrophil_0.8_matched.csv")
  )

  obj@misc$neutrophil_0.8_annotation_match <- list(
    reference_file = normalizePath(reference_file, winslash = "/"),
    direct_match_patients = setdiff(unique(obj$PatientID), "P07"),
    direct_match_cells = length(direct_cells),
    P07_transfer_method = "25-nearest-neighbor majority vote in global Harmony dimensions 1:10",
    P07_cells = length(p07_cells),
    P07_median_vote_confidence = median(p07_confidence),
    date = as.character(Sys.Date())
  )
  saveRDS(obj, output_file)
  message("Saved matched full object: ", output_file)
  print(sort(table(obj$Celltype_updated[grepl("^Neutrophil", obj$Celltype_updated)]), decreasing = TRUE))
}

if (mode == "neutrophil_umap") {
  message("Rebuilding neutrophil-only UMAP with matched 0.8 annotations")
  obj <- readRDS(output_file)
  neutrophil_cells <- colnames(obj)[grepl("^Neutrophil", obj$Celltype_updated)]
  neut <- subset(obj, cells = neutrophil_cells)
  neut <- DietSeurat(
    neut,
    assays = "Akoya",
    layers = "counts",
    dimreducs = NULL,
    graphs = NULL,
    misc = TRUE
  )
  essential <- intersect(
    c(
      "orig.ident", "nCount_Akoya", "nFeature_Akoya", "PatientID", "Diagnosis1", "Diagnosis2",
      "Inflammation", "Medication", "Celltype", "Celltype_original", "Celltype_updated",
      "Neutrophil_subtype_0.8", "Neutrophil_annotation_0.8", "Neutrophil_cluster_0.8",
      "Neutrophil_subset_annotation"
    ),
    colnames(neut[[]])
  )
  neut@meta.data <- neut@meta.data[, essential, drop = FALSE]
  DefaultAssay(neut) <- "Akoya"
  features <- rownames(neut)
  neut <- NormalizeData(neut, normalization.method = "CLR", margin = 2, verbose = FALSE)
  neut <- ScaleData(neut, features = features, verbose = FALSE)
  neut <- RunPCA(neut, features = features, npcs = 30, seed.use = 42, verbose = FALSE)
  neut <- RunHarmony(
    neut,
    group.by.vars = "PatientID",
    reduction.use = "pca",
    dims.use = 1:30,
    project.dim = TRUE,
    verbose = TRUE
  )
  neut <- RunUMAP(
    neut,
    reduction = "harmony",
    dims = 1:30,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    n.epochs = 200,
    seed.use = 42,
    verbose = TRUE
  )
  Idents(neut) <- "Neutrophil_subset_annotation"
  saveRDS(neut, neutrophil_file)

  plot_dt <- as.data.table(Embeddings(neut, "umap"), keep.rownames = "cell_ID")
  setnames(plot_dt, 2:3, c("UMAP_1", "UMAP_2"))
  meta <- as.data.table(neut[[]], keep.rownames = "cell_ID")
  plot_dt <- merge(
    plot_dt,
    meta[, .(cell_ID, Neutrophil_subset_annotation, Diagnosis2, PatientID)],
    by = "cell_ID",
    sort = FALSE
  )
  annotation_levels <- c(
    "Neutrophil_MX1", "Neutrophil_CXCR4", "Neutrophil_PADI4", "Neutrophil_OSM", "Neutrophil"
  )
  annotation_levels <- annotation_levels[annotation_levels %in% plot_dt$Neutrophil_subset_annotation]
  plot_dt[, Neutrophil_subset_annotation := factor(Neutrophil_subset_annotation, levels = annotation_levels)]
  annotation_colors <- c(
    "Neutrophil_MX1" = "#3B82F6",
    "Neutrophil_CXCR4" = "#8B5CF6",
    "Neutrophil_PADI4" = "#EF4444",
    "Neutrophil_OSM" = "#F59E0B",
    "Neutrophil" = "#6B7280"
  )
  group_colors <- c("Control" = "#4C78A8", "UC_Noninflamed" = "#F2CF5B", "UC_Inflamed" = "#E45756")
  centers <- plot_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = Neutrophil_subset_annotation]
  theme_neut <- theme_minimal(base_size = 11) + theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
  p_annotation <- ggplot(plot_dt, aes(UMAP_1, UMAP_2, colour = Neutrophil_subset_annotation)) +
    geom_point(size = 0.62, alpha = 0.78) +
    ggrepel::geom_label_repel(
      data = centers,
      aes(label = Neutrophil_subset_annotation),
      colour = "black", fill = alpha("white", 0.85), linewidth = 0.25,
      fontface = "bold", size = 3.2, seed = 42, min.segment.length = 0,
      max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_colour_manual(values = annotation_colors, drop = FALSE) + coord_equal() +
    labs(
      title = "Neutrophil-only Harmony UMAP",
      subtitle = "Matched to UCCODEX1 neutrophil 0.8 annotations",
      x = "UMAP 1", y = "UMAP 2", colour = "0.8 annotation"
    ) + theme_neut + theme(panel.grid = element_blank())
  p_group <- ggplot(plot_dt, aes(UMAP_1, UMAP_2, colour = Diagnosis2)) +
    geom_point(size = 0.62, alpha = 0.72) +
    scale_colour_manual(values = group_colors, drop = FALSE) + coord_equal() +
    labs(
      title = "Tissue-group distribution",
      subtitle = "Patient-corrected Harmony embedding",
      x = "UMAP 1", y = "UMAP 2", colour = "Tissue group"
    ) + theme_neut + theme(panel.grid = element_blank())
  figure <- p_annotation | p_group
  ggsave(file.path(figure_dir, "Figure_06_neutrophil_subcluster_umap.png"), figure,
         width = 14, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, "Figure_06_neutrophil_subcluster_umap.pdf"), figure,
         width = 14, height = 7, bg = "white")

  counts <- as.data.table(neut[[]], keep.rownames = "cell_ID")[
    , .N,
    by = .(Neutrophil_subset_annotation, Diagnosis2)
  ]
  fwrite(counts, file.path(table_dir, "neutrophil_0.8_annotation_counts_by_tissue_group.csv"))
  fwrite(
    as.data.table(neut[[]], keep.rownames = "cell_ID")[
      , .N,
      by = .(PatientID, Diagnosis2, Neutrophil_subset_annotation)
    ],
    file.path(table_dir, "neutrophil_0.8_annotation_counts_by_patient.csv")
  )
  message("Saved matched neutrophil object and updated Figure 06")
  print(sort(table(neut$Neutrophil_subset_annotation), decreasing = TRUE))
}
