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
mode <- if (length(args)) args[1] else "cluster"
valid_modes <- c("extract", "cluster", "attach", "figures")
if (!mode %in% valid_modes) stop("Mode must be one of: ", paste(valid_modes, collapse = ", "))

input_file <- "UCCODEX1_Annotated_repaired_withUMAP.rds"
output_file <- "UCCODEX1_Annotated_repaired_withUMAP_neutrophilSubclusters.rds"
output_dir <- "codex_repair_P02_P07/neutrophil_subclusters"
input_subset_file <- file.path(output_dir, "neutrophil_input.rds")
neutrophil_file <- file.path(output_dir, "neutrophil_subclustered.rds")
figure_dir <- "giotto_codex_results/figures"
table_dir <- "giotto_codex_results/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (mode == "extract") {
  message("Extracting neutrophils from repaired UMAP-enabled object")
  obj <- readRDS(input_file)
  neutrophil_cells <- colnames(obj)[obj$Celltype == "Neutrophil"]
  if (length(neutrophil_cells) < 100L) stop("Too few neutrophils: ", length(neutrophil_cells))
  neut <- subset(obj, cells = neutrophil_cells)
  neut <- DietSeurat(
    neut,
    assays = "Akoya",
    layers = "counts",
    dimreducs = NULL,
    graphs = NULL,
    misc = TRUE
  )
  essential_metadata <- intersect(
    c(
      "orig.ident", "nCount_Akoya", "nFeature_Akoya", "PatientID", "Diagnosis1",
      "Diagnosis2", "Inflammation", "Medication", "Celltype"
    ),
    colnames(neut[[]])
  )
  neut@meta.data <- neut@meta.data[, essential_metadata, drop = FALSE]
  saveRDS(neut, input_subset_file)
  message("Saved ", ncol(neut), " neutrophils from ", length(unique(neut$PatientID)), " patients")
}

if (mode == "cluster") {
  message("Recomputing neutrophil PCA, Harmony, neighbors, clusters, and UMAP")
  neut <- readRDS(input_subset_file)
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
  neut <- FindNeighbors(neut, reduction = "harmony", dims = 1:15, verbose = FALSE)
  neut <- FindClusters(neut, resolution = 0.4, n.start = 10, random.seed = 0, verbose = FALSE)
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
  neut$Neutrophil_subcluster <- paste0("Neu_", as.character(neut$seurat_clusters))
  neut$Neutrophil_subset_annotation <- neut$Neutrophil_subcluster
  Idents(neut) <- "Neutrophil_subset_annotation"

  cluster_sizes <- as.data.table(neut[[]], keep.rownames = "cell_ID")[
    , .N,
    by = .(Neutrophil_subcluster, Diagnosis2)
  ]
  cluster_totals <- cluster_sizes[, .(n_cells = sum(N)), by = Neutrophil_subcluster]
  patient_counts <- as.data.table(neut[[]], keep.rownames = "cell_ID")[
    , .N,
    by = .(PatientID, Diagnosis2, Neutrophil_subcluster)
  ]

  data_matrix <- LayerData(neut, assay = "Akoya", layer = "data")
  subcluster_levels <- sort(unique(neut$Neutrophil_subcluster))
  average_expression <- sapply(subcluster_levels, function(cluster) {
    rowMeans(data_matrix[, neut$Neutrophil_subcluster == cluster, drop = FALSE])
  })
  rownames(average_expression) <- rownames(neut)
  expression_z <- t(scale(t(average_expression)))
  expression_z[!is.finite(expression_z)] <- 0
  top_markers <- rbindlist(lapply(seq_along(subcluster_levels), function(i) {
    ord <- order(expression_z[, i], decreasing = TRUE)
    data.table(
      Neutrophil_subcluster = subcluster_levels[i],
      rank = seq_len(min(10L, length(ord))),
      marker = rownames(expression_z)[ord[seq_len(min(10L, length(ord)))]],
      average_CLR_expression = average_expression[ord[seq_len(min(10L, length(ord)))], i],
      subcluster_z_score = expression_z[ord[seq_len(min(10L, length(ord)))], i]
    )
  }))

  fwrite(cluster_totals, file.path(table_dir, "neutrophil_subcluster_sizes.csv"))
  fwrite(cluster_sizes, file.path(table_dir, "neutrophil_subcluster_sizes_by_tissue_group.csv"))
  fwrite(patient_counts, file.path(table_dir, "neutrophil_subcluster_counts_by_patient.csv"))
  fwrite(
    as.data.table(average_expression, keep.rownames = "marker"),
    file.path(table_dir, "neutrophil_subcluster_average_marker_expression.csv")
  )
  fwrite(top_markers, file.path(table_dir, "neutrophil_subcluster_top_markers.csv"))

  neut@misc$neutrophil_subclustering <- list(
    source_annotation = "Celltype == Neutrophil",
    n_cells = ncol(neut),
    pca_dimensions = 1:30,
    harmony_group = "PatientID",
    neighbor_dimensions = 1:15,
    resolution = 0.4,
    umap_dimensions = 1:30,
    umap_metric = "cosine",
    umap_neighbors = 30L,
    umap_min_dist = 0.3,
    seed = 42L,
    date = as.character(Sys.Date())
  )
  saveRDS(neut, neutrophil_file)
  message("Saved ", nlevels(Idents(neut)), " neutrophil subclusters")
  print(sort(table(neut$Neutrophil_subcluster), decreasing = TRUE))
}

if (mode == "attach") {
  message("Adding neutrophil subclusters to full Seurat metadata")
  obj <- readRDS(input_file)
  neut <- readRDS(neutrophil_file)
  subclusters <- setNames(as.character(neut$Neutrophil_subcluster), colnames(neut))
  cluster_ids <- setNames(as.character(neut$seurat_clusters), colnames(neut))
  obj$Neutrophil_subcluster <- unname(subclusters[colnames(obj)])
  obj$Neutrophil_cluster <- unname(cluster_ids[colnames(obj)])
  obj$Neutrophil_subset_annotation <- obj$Neutrophil_subcluster
  obj$Celltype_refined <- as.character(obj$Celltype)
  is_neutrophil <- !is.na(obj$Neutrophil_subset_annotation)
  obj$Celltype_refined[is_neutrophil] <- obj$Neutrophil_subset_annotation[is_neutrophil]
  neut$Neutrophil_subset_annotation <- as.character(neut$Neutrophil_subcluster)
  Idents(neut) <- "Neutrophil_subset_annotation"
  stopifnot(
    sum(!is.na(obj$Neutrophil_subcluster)) == ncol(neut),
    all(obj$Celltype[!is.na(obj$Neutrophil_subcluster)] == "Neutrophil"),
    identical(
      obj$Celltype_refined[is_neutrophil],
      obj$Neutrophil_subset_annotation[is_neutrophil]
    )
  )
  obj@misc$neutrophil_subclustering <- neut@misc$neutrophil_subclustering
  obj@misc$neutrophil_subclustering$separate_object <- normalizePath(neutrophil_file, winslash = "/")
  saveRDS(neut, neutrophil_file)
  saveRDS(obj, output_file)
  message("Saved full object with neutrophil subclusters: ", output_file)
}

if (mode == "figures") {
  message("Creating neutrophil-only UMAP figure")
  neut <- readRDS(neutrophil_file)
  plot_dt <- as.data.table(Embeddings(neut, "umap"), keep.rownames = "cell_ID")
  setnames(plot_dt, 2:3, c("UMAP_1", "UMAP_2"))
  meta <- as.data.table(neut[[]], keep.rownames = "cell_ID")
  if (!"Neutrophil_subset_annotation" %in% names(meta)) {
    meta[, Neutrophil_subset_annotation := Neutrophil_subcluster]
  }
  plot_dt <- merge(
    plot_dt,
    meta[, .(cell_ID, Neutrophil_subset_annotation, Diagnosis2, PatientID)],
    by = "cell_ID",
    sort = FALSE
  )
  cluster_levels <- sort(unique(plot_dt$Neutrophil_subset_annotation))
  plot_dt[, Neutrophil_subset_annotation := factor(Neutrophil_subset_annotation, levels = cluster_levels)]
  cluster_colors <- setNames(hue_pal()(length(cluster_levels)), cluster_levels)
  group_colors <- c("Control" = "#4C78A8", "UC_Noninflamed" = "#F2CF5B", "UC_Inflamed" = "#E45756")
  centers <- plot_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = Neutrophil_subset_annotation]
  theme_neut <- theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
      plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold")
    )

  p_cluster <- ggplot(plot_dt, aes(UMAP_1, UMAP_2, colour = Neutrophil_subset_annotation)) +
    geom_point(size = 0.62, alpha = 0.78) +
    ggrepel::geom_label_repel(
      data = centers,
      aes(label = Neutrophil_subset_annotation),
      colour = "black", fill = alpha("white", 0.85), linewidth = 0.25,
      fontface = "bold", size = 3.2, seed = 42, min.segment.length = 0,
      max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_colour_manual(values = cluster_colors) + coord_equal() +
    labs(
      title = "Neutrophil-only Harmony UMAP",
      subtitle = paste0(nrow(plot_dt), " neutrophils; labels from Neutrophil_subset_annotation metadata"),
      x = "UMAP 1", y = "UMAP 2", colour = "Subcluster"
    ) + theme_neut + theme(panel.grid = element_blank())

  p_group <- ggplot(plot_dt, aes(UMAP_1, UMAP_2, colour = Diagnosis2)) +
    geom_point(size = 0.62, alpha = 0.72) +
    scale_colour_manual(values = group_colors, drop = FALSE) + coord_equal() +
    labs(
      title = "Tissue-group distribution",
      subtitle = "Patient-corrected Harmony embedding",
      x = "UMAP 1", y = "UMAP 2", colour = "Tissue group"
    ) + theme_neut + theme(panel.grid = element_blank())

  figure <- p_cluster | p_group
  ggsave(file.path(figure_dir, "Figure_06_neutrophil_subcluster_umap.png"), figure,
         width = 14, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, "Figure_06_neutrophil_subcluster_umap.pdf"), figure,
         width = 14, height = 7, bg = "white")
  message("Saved Figure 06 neutrophil UMAP")
}
