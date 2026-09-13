#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "all"
valid_modes <- c("extract", "compute", "attach", "figures")
if (!mode %in% valid_modes) stop("Mode must be one of: ", paste(valid_modes, collapse = ", "))

input_file <- "UCCODEX1_Annotated_repaired.rds"
output_file <- "UCCODEX1_Annotated_repaired_withUMAP.rds"
work_dir <- "codex_repair_P02_P07"
payload_file <- file.path(work_dir, "harmony_embedding_payload.rds")
umap_file <- file.path(work_dir, "harmony_umap_coordinates.rds")
figure_dir <- file.path("giotto_codex_results", "figures")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (mode == "extract") {
  message("Loading repaired Seurat object and extracting Harmony embeddings")
  obj <- readRDS(input_file)
  harmony <- Embeddings(obj, reduction = "harmony")
  stopifnot(nrow(harmony) == ncol(obj), ncol(harmony) >= 30L)
  payload <- list(
    harmony = harmony[, 1:30, drop = FALSE],
    cell_ID = colnames(obj),
    Celltype = as.character(obj$Celltype),
    seurat_clusters = as.character(obj$seurat_clusters),
    PatientID = as.character(obj$PatientID),
    Diagnosis2 = as.character(obj$Diagnosis2)
  )
  saveRDS(payload, payload_file, compress = FALSE)
  message("Saved ", nrow(payload$harmony), " x ", ncol(payload$harmony), " Harmony payload")
}

if (mode == "compute") {
  message("Computing UMAP from Harmony dimensions 1:30 in an isolated process")
  payload <- readRDS(payload_file)
  umap <- uwot::umap(
    X = payload$harmony,
    n_neighbors = 30,
    n_components = 2,
    metric = "cosine",
    n_epochs = 200,
    learning_rate = 1,
    min_dist = 0.3,
    spread = 1,
    set_op_mix_ratio = 1,
    local_connectivity = 1,
    repulsion_strength = 1,
    negative_sample_rate = 5,
    init = "spectral",
    n_threads = min(4L, parallel::detectCores(logical = TRUE)),
    n_sgd_threads = 1,
    ret_model = FALSE,
    ret_nn = FALSE,
    verbose = TRUE,
    seed = 42
  )
  rownames(umap) <- payload$cell_ID
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  stopifnot(nrow(umap) == 245020L, all(is.finite(umap)))
  saveRDS(umap, umap_file)
  message("Saved UMAP coordinates to ", umap_file)
}

if (mode == "attach") {
  message("Attaching UMAP reduction to repaired Seurat object")
  obj <- readRDS(input_file)
  umap <- readRDS(umap_file)
  umap <- umap[colnames(obj), , drop = FALSE]
  stopifnot(identical(rownames(umap), colnames(obj)), all(is.finite(umap)))
  obj[["umap"]] <- CreateDimReducObject(
    embeddings = umap,
    key = "UMAP_",
    assay = "Akoya",
    global = TRUE,
    misc = list(
      source = "Harmony dimensions 1:30",
      method = "uwot cosine UMAP",
      n_neighbors = 30L,
      min_dist = 0.3,
      n_epochs = 200L,
      seed = 42L
    )
  )
  obj@misc$P02_P07_repair$UMAP <- list(
    source_reduction = "harmony",
    dimensions = 1:30,
    method = "uwot cosine",
    n_neighbors = 30L,
    min_dist = 0.3,
    n_epochs = 200L,
    seed = 42L,
    date = as.character(Sys.Date())
  )
  saveRDS(obj, output_file)
  message("Saved Seurat object with UMAP: ", output_file)
}

if (mode == "figures") {
  message("Regenerating cohort overview with Harmony-derived UMAP")
  payload <- readRDS(payload_file)
  umap <- readRDS(umap_file)
  plot_dt <- data.table(
    cell_ID = payload$cell_ID,
    UMAP_1 = umap[payload$cell_ID, 1],
    UMAP_2 = umap[payload$cell_ID, 2],
    Celltype = payload$Celltype,
    seurat_clusters = factor(payload$seurat_clusters, levels = sort(unique(as.integer(payload$seurat_clusters)))),
    PatientID = payload$PatientID,
    Diagnosis2 = payload$Diagnosis2
  )

  celltype_colors <- c(
    "B Cell" = "#1F78B4", "CD4 T" = "#33A02C", "CD8 T" = "#6A3D9A",
    "Dendritic Cell" = "#3CA3D0", "Endothelial Cell" = "#B15928",
    "Enteroendocrine Cell" = "#E6A01A", "Epithelial Cell" = "#FF7F00",
    "Fibroblast" = "#CAB2D6", "Macrophage" = "#E31A1C", "Neutrophil" = "#F06CA9",
    "Plasma B Cell" = "#8BCB52", "TReg" = "#C8A600"
  )
  group_colors <- c("Control" = "#4C78A8", "UC_Noninflamed" = "#F2CF5B", "UC_Inflamed" = "#E45756")
  theme_codex <- function(base_size = 10) {
    theme_minimal(base_size = base_size) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      )
  }

  sample_counts <- fread(file.path("giotto_codex_results", "tables", "sample_summary_all.csv"))
  sample_counts[, Diagnosis2 := factor(Diagnosis2, levels = names(group_colors))]
  composition <- fread(file.path("giotto_codex_results", "tables", "celltype_composition_by_patient.csv"))
  composition[, Diagnosis2 := factor(Diagnosis2, levels = names(group_colors))]

  p1a <- ggplot(sample_counts[analysis_included == TRUE],
                aes(x = reorder(PatientID, n_cells), y = n_cells, fill = Diagnosis2)) +
    geom_col(width = 0.8) + coord_flip() +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Cells", title = "Cells per patient", fill = "Tissue group") + theme_codex(9)

  p1b <- ggplot(composition, aes(x = PatientID, y = proportion, fill = Celltype)) +
    geom_col(width = 0.9) + scale_fill_manual(values = celltype_colors, drop = FALSE) +
    scale_y_continuous(labels = percent_format()) +
    labs(x = NULL, y = "Cell-type composition", title = "Patient-level composition", fill = "Cell type") +
    theme_codex(9) + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right")

  set.seed(42)
  umap_plot <- if (nrow(plot_dt) > 80000L) plot_dt[sample(.N, 80000L)] else plot_dt
  celltype_centers <- plot_dt[, .(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  ), by = Celltype]
  p1c <- ggplot(umap_plot, aes(x = UMAP_1, y = UMAP_2, colour = Celltype)) +
    geom_point(size = 0.24, alpha = 0.62) + scale_colour_manual(values = celltype_colors, drop = FALSE) +
    ggrepel::geom_label_repel(
      data = celltype_centers,
      aes(label = Celltype),
      colour = "black", fill = alpha("white", 0.82), linewidth = 0.2,
      fontface = "bold", size = 2.5, seed = 42,
      min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
    ) +
    coord_equal() + labs(x = "UMAP 1", y = "UMAP 2", title = "Harmony-derived UMAP", colour = "Cell type") +
    theme_codex(9) + theme(panel.grid = element_blank(), legend.position = "none")

  fig1 <- (p1a | p1c) / p1b + plot_annotation(title = "CODEX cohort overview")
  ggsave(file.path(figure_dir, "Figure_01_cohort_overview.png"), fig1, width = 13, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, "Figure_01_cohort_overview.pdf"), fig1, width = 13, height = 10, bg = "white")

  p_annotation <- ggplot(umap_plot, aes(UMAP_1, UMAP_2, colour = Celltype)) +
    geom_point(size = 0.30, alpha = 0.68) +
    ggrepel::geom_label_repel(
      data = celltype_centers,
      aes(label = Celltype),
      colour = "black", fill = alpha("white", 0.85), linewidth = 0.25,
      fontface = "bold", size = 3.2, seed = 42,
      min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_colour_manual(values = celltype_colors, drop = FALSE) + coord_equal() +
    labs(
      title = "Harmony-derived UMAP annotations",
      subtitle = "Cell-type annotations from Seurat metadata",
      x = "UMAP 1", y = "UMAP 2", colour = "Cell type"
    ) + theme_codex(11) + theme(panel.grid = element_blank(), legend.position = "right")
  ggsave(file.path(figure_dir, "Figure_01b_harmony_umap_annotations.png"), p_annotation, width = 11, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, "Figure_01b_harmony_umap_annotations.pdf"), p_annotation, width = 11, height = 8, bg = "white")

  cluster_centers <- plot_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = seurat_clusters]
  cluster_colors <- setNames(hue_pal()(length(levels(plot_dt$seurat_clusters))), levels(plot_dt$seurat_clusters))
  p_cluster <- ggplot(umap_plot, aes(UMAP_1, UMAP_2, colour = seurat_clusters)) +
    geom_point(size = 0.30, alpha = 0.68) +
    geom_label(
      data = cluster_centers,
      aes(label = seurat_clusters),
      colour = "black", fill = alpha("white", 0.78), linewidth = 0.15,
      fontface = "bold", size = 3, show.legend = FALSE
    ) +
    scale_colour_manual(values = cluster_colors) + coord_equal() +
    labs(
      title = "Harmony-derived UMAP clusters",
      subtitle = "Louvain clustering on Harmony dimensions 1:10; UMAP on Harmony dimensions 1:30",
      x = "UMAP 1", y = "UMAP 2", colour = "Cluster"
    ) + theme_codex(11) + theme(panel.grid = element_blank(), legend.position = "right")
  ggsave(file.path(figure_dir, "Figure_01c_harmony_umap_clusters.png"), p_cluster, width = 10, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, "Figure_01c_harmony_umap_clusters.pdf"), p_cluster, width = 10, height = 8, bg = "white")
  message("Updated Figure 01 and created annotated Figure 01b plus cluster Figure 01c")
}
