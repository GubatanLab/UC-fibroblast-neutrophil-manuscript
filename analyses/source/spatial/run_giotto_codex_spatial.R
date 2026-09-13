#!/usr/bin/env Rscript

# Reproducible Giotto spatial analysis for UCCODEX1_Annotated.rds
# Outputs patient-level Giotto objects, tables, and publication-ready figures.

options(stringsAsFactors = FALSE, giotto.no_python_warn = TRUE)
set.seed(20260710)

suppressPackageStartupMessages({
  library(Seurat)
  library(Giotto)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dbscan)
})

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args
n_sim <- if (quick_mode) 10L else 100L
input_file <- "UCCODEX1_Annotated_repaired_withUMAP_neutrophil0.8Matched.rds"
output_dir <- if (quick_mode) "giotto_codex_results_quick" else "giotto_codex_results"

dirs <- file.path(output_dir, c("tables", "figures", "figures", "giotto_objects", "logs"))
invisible(lapply(unique(dirs), dir.create, recursive = TRUE, showWarnings = FALSE))
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
object_dir <- file.path(output_dir, "giotto_objects")
log_dir <- file.path(output_dir, "logs")
log_file <- file.path(log_dir, "analysis_log.txt")
uv_cache_dir <- file.path(output_dir, ".uv_cache")
dir.create(uv_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(UV_CACHE_DIR = normalizePath(uv_cache_dir, winslash = "/", mustWork = TRUE))

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

write_csv <- function(x, filename) {
  data.table::fwrite(as.data.table(x), file.path(table_dir, filename))
}

save_plot_both <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white", limitsize = FALSE)
}

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

celltype_colors <- c(
  "B Cell" = "#1F78B4",
  "CD4 T" = "#33A02C",
  "CD8 T" = "#6A3D9A",
  "Dendritic Cell" = "#3CA3D0",
  "Endothelial Cell" = "#B15928",
  "Enteroendocrine Cell" = "#E6A01A",
  "Epithelial Cell" = "#FF7F00",
  "Fibroblast" = "#CAB2D6",
  "Macrophage" = "#E31A1C",
  "Neutrophil" = "#F06CA9",
  "Plasma B Cell" = "#8BCB52",
  "TReg" = "#C8A600"
)

group_colors <- c(
  "Control" = "#4C78A8",
  "UC_Noninflamed" = "#F2CF5B",
  "UC_Inflamed" = "#E45756"
)

safe_wilcox <- function(value, group, a, b) {
  xa <- value[group == a & is.finite(value)]
  xb <- value[group == b & is.finite(value)]
  if (length(xa) < 2L || length(xb) < 2L) return(NA_real_)
  tryCatch(wilcox.test(xa, xb, exact = FALSE)$p.value, error = function(e) NA_real_)
}

pairwise_tests <- function(dt, value_col, feature_cols) {
  contrasts <- list(
    c("UC_Noninflamed", "Control"),
    c("UC_Inflamed", "Control"),
    c("UC_Inflamed", "UC_Noninflamed")
  )
  features <- unique(dt[, ..feature_cols])
  out <- vector("list", nrow(features) * length(contrasts))
  idx <- 0L
  for (i in seq_len(nrow(features))) {
    keep <- rep(TRUE, nrow(dt))
    for (fc in feature_cols) keep <- keep & dt[[fc]] == features[[fc]][i]
    sub <- dt[keep]
    for (contrast in contrasts) {
      idx <- idx + 1L
      a <- contrast[1]
      b <- contrast[2]
      va <- sub[Diagnosis2 == a, get(value_col)]
      vb <- sub[Diagnosis2 == b, get(value_col)]
      out[[idx]] <- cbind(
        features[i],
        data.table(
          contrast = paste(a, "vs", b),
          group_1 = a,
          group_2 = b,
          n_group_1 = sum(is.finite(va)),
          n_group_2 = sum(is.finite(vb)),
          median_group_1 = median(va, na.rm = TRUE),
          median_group_2 = median(vb, na.rm = TRUE),
          median_difference = median(va, na.rm = TRUE) - median(vb, na.rm = TRUE),
          p_value = safe_wilcox(sub[[value_col]], sub$Diagnosis2, a, b)
        )
      )
    }
  }
  ans <- rbindlist(out, fill = TRUE)
  ans[, p_adj_BH := p.adjust(p_value, method = "BH"), by = contrast]
  ans[]
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
log_msg("Loading", normalizePath(input_file))
seu <- readRDS(input_file)
log_msg("Loaded", ncol(seu), "cells and", nrow(seu), "features")

if (!"Celltype" %in% colnames(seu[[]])) stop("Required metadata column Celltype is absent")
if (!all(c("PatientID", "Diagnosis1", "Diagnosis2", "Inflammation", "Medication") %in% colnames(seu[[]]))) {
  stop("One or more required patient metadata columns are absent")
}

# Assemble spatial coordinates directly from all Seurat FOV centroid objects.
coord_list <- lapply(Images(seu), function(image_name) {
  x <- as.data.table(GetTissueCoordinates(seu, image = image_name))
  setnames(x, "cell", "cell_ID")
  x[, seurat_image := image_name]
  x
})
coords <- rbindlist(coord_list, use.names = TRUE)
if (nrow(coords) != ncol(seu) || anyDuplicated(coords$cell_ID)) {
  stop("Spatial coordinate assembly failed: coordinates do not map one-to-one to Seurat cells")
}

metadata_cols <- c(
  "PatientID", "Diagnosis1", "Diagnosis2", "Inflammation", "Medication",
  "Celltype", "seurat_clusters", "nCount_Akoya", "nFeature_Akoya"
)
cell_meta <- as.data.table(seu[[]], keep.rownames = "cell_ID")
cell_meta <- cell_meta[, c("cell_ID", intersect(metadata_cols, names(cell_meta))), with = FALSE]
cell_meta <- merge(cell_meta, coords, by = "cell_ID", all.x = TRUE, sort = FALSE)
setkey(cell_meta, cell_ID)
cell_meta <- cell_meta[J(colnames(seu))]
if (anyNA(cell_meta$x) || anyNA(cell_meta$y)) stop("Missing spatial coordinates after metadata join")

# Exact duplicate audit across the formerly duplicated P02/P07 samples.
p02_cells <- cell_meta[PatientID == "P02", cell_ID]
p07_cells <- cell_meta[PatientID == "P07", cell_ID]
same_length <- length(p02_cells) == length(p07_cells) && length(p02_cells) > 0L
coordinates_identical <- same_length && identical(
  unname(as.matrix(cell_meta[p02_cells, on = "cell_ID", .(x, y)])),
  unname(as.matrix(cell_meta[p07_cells, on = "cell_ID", .(x, y)]))
)
expr_all <- LayerData(seu, assay = "Akoya", layer = "data")
expression_identical <- same_length && identical(
  unname(as.matrix(expr_all[, p02_cells, drop = FALSE])),
  unname(as.matrix(expr_all[, p07_cells, drop = FALSE]))
)
duplicate_audit <- data.table(
  patient_kept = "P02",
  patient_excluded = "P07",
  n_cells_each = if (same_length) length(p02_cells) else NA_integer_,
  coordinates_identical = coordinates_identical,
  expression_identical_all_54_markers = expression_identical,
  reason = if (coordinates_identical && expression_identical) {
    "Exact technical duplicate; excluded from inferential and group summaries"
  } else {
    "Not an exact duplicate; retained"
  }
)
write_csv(duplicate_audit, "qc_exact_duplicate_audit.csv")
excluded_patients <- if (coordinates_identical && expression_identical) "P07" else character()

sample_meta <- unique(cell_meta[, .(PatientID, Diagnosis1, Diagnosis2, Inflammation, Medication)])
sample_counts_all <- cell_meta[, .(n_cells = .N), by = .(PatientID, Diagnosis1, Diagnosis2, Inflammation, Medication)]
sample_counts_all[, analysis_included := !PatientID %in% excluded_patients]
setorder(sample_counts_all, PatientID)
write_csv(sample_counts_all, "sample_summary_all.csv")

analysis_cells <- cell_meta[!PatientID %in% excluded_patients, cell_ID]
analysis_meta <- cell_meta[cell_ID %chin% analysis_cells]
analysis_patients <- sort(unique(analysis_meta$PatientID))
if (quick_mode) analysis_patients <- head(analysis_patients, 3L)
analysis_meta <- analysis_meta[PatientID %in% analysis_patients]
log_msg("Analysis includes", length(analysis_patients), "patients and", nrow(analysis_meta), "cells; excluded", paste(excluded_patients, collapse = ","))

# Nearest-neighbor diagnostics support the 50-coordinate-unit Delaunay edge cap.
nn_diagnostics <- rbindlist(lapply(analysis_patients, function(pid) {
  xy <- as.matrix(analysis_meta[PatientID == pid, .(x, y)])
  nn <- dbscan::kNN(xy, k = 1)$dist[, 1]
  data.table(
    PatientID = pid,
    n_cells = nrow(xy),
    nn_min = min(nn),
    nn_q25 = quantile(nn, 0.25),
    nn_median = median(nn),
    nn_q75 = quantile(nn, 0.75),
    nn_q95 = quantile(nn, 0.95),
    nn_q99 = quantile(nn, 0.99),
    nn_max = max(nn)
  )
}))
nn_diagnostics <- merge(nn_diagnostics, sample_meta, by = "PatientID", all.x = TRUE)
write_csv(nn_diagnostics, "qc_nearest_neighbor_distances.csv")

# Patient-level composition analysis.
all_celltypes <- names(celltype_colors)
composition <- analysis_meta[, .(cell_count = .N), by = .(PatientID, Diagnosis1, Diagnosis2, Inflammation, Medication, Celltype)]
composition <- CJ(PatientID = analysis_patients, Celltype = all_celltypes)[composition, on = .(PatientID, Celltype)]
composition <- merge(composition, sample_meta, by = "PatientID", all.x = TRUE, suffixes = c("", ".sample"))
for (nm in c("Diagnosis1", "Diagnosis2", "Inflammation", "Medication")) {
  sample_nm <- paste0(nm, ".sample")
  if (sample_nm %in% names(composition)) {
    composition[is.na(get(nm)), (nm) := get(sample_nm)]
    composition[, (sample_nm) := NULL]
  }
}
composition[is.na(cell_count), cell_count := 0L]
composition[, total_cells := sum(cell_count), by = PatientID]
composition[, proportion := cell_count / total_cells]
composition[, Diagnosis2 := factor(Diagnosis2, levels = names(group_colors))]
write_csv(composition, "celltype_composition_by_patient.csv")

composition_tests <- pairwise_tests(composition, "proportion", "Celltype")
write_csv(composition_tests, "statistics_celltype_composition_pairwise_wilcoxon.csv")

# Figure 1: dataset overview.
sample_counts_plot <- copy(sample_counts_all[analysis_included == TRUE & PatientID %in% analysis_patients])
sample_counts_plot[, Diagnosis2 := factor(Diagnosis2, levels = names(group_colors))]
p1a <- ggplot(sample_counts_plot, aes(x = reorder(PatientID, n_cells), y = n_cells, fill = Diagnosis2)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Cells", title = "Cells per patient", fill = "Tissue group") +
  theme_codex(9)

p1b <- ggplot(composition, aes(x = PatientID, y = proportion, fill = Celltype)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Cell-type composition", title = "Patient-level composition", fill = "Cell type") +
  theme_codex(9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right")

embedding_name <- if ("umap" %in% Reductions(seu)) "umap" else "harmony"
embedding <- as.data.table(Embeddings(seu, embedding_name)[, 1:2, drop = FALSE], keep.rownames = "cell_ID")
setnames(embedding, 2:3, c("embedding_1", "embedding_2"))
embedding <- merge(embedding, analysis_meta[, .(cell_ID, Celltype)], by = "cell_ID")
if (nrow(embedding) > 60000L) embedding <- embedding[sample(.N, 60000L)]
p1c <- ggplot(embedding, aes(x = embedding_1, y = embedding_2, colour = Celltype)) +
  geom_point(size = 0.24, alpha = 0.58) +
  scale_colour_manual(values = celltype_colors, drop = FALSE) +
  coord_equal() +
  labs(
    x = paste(toupper(embedding_name), "1"),
    y = paste(toupper(embedding_name), "2"),
    title = "Phenotypic embedding",
    colour = "Cell type"
  ) +
  theme_codex(9) +
  theme(panel.grid = element_blank(), legend.position = "none")

fig1 <- (p1a | p1c) / p1b + plot_annotation(title = "CODEX cohort overview")
save_plot_both(fig1, "Figure_01_cohort_overview", 13, 10)

# Figure 2: spatial cell maps. Each facet is an independent patient tissue region.
spatial_plot_data <- analysis_meta[, .(cell_ID, PatientID, Diagnosis2, Celltype, x, y)]
spatial_plot_data[, Diagnosis2 := factor(Diagnosis2, levels = names(group_colors))]
p2 <- ggplot(spatial_plot_data, aes(x = x, y = y, colour = Celltype)) +
  geom_point(size = 0.48, alpha = 0.95, stroke = 0) +
  facet_wrap(~ PatientID, scales = "free", ncol = 4) +
  scale_colour_manual(values = celltype_colors, drop = FALSE) +
  scale_y_reverse() +
  labs(x = NULL, y = NULL, title = "Spatial cell-type maps", colour = "Cell type") +
  theme_void(base_size = 8) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    aspect.ratio = 1
  )
save_plot_both(p2, "Figure_02_spatial_celltype_maps", 15, 20)

# Create patient-level Giotto objects and Delaunay networks, then run CPE.
cpe_results <- list()
network_summary <- list()

for (pid in analysis_patients) {
  log_msg("Giotto patient", pid, "started")
  cells <- analysis_meta[PatientID == pid, cell_ID]
  expr <- expr_all[, cells, drop = FALSE]
  meta <- as.data.frame(analysis_meta[cell_ID %chin% cells, .(
    cell_ID, PatientID, Diagnosis1, Diagnosis2, Inflammation, Medication,
    Celltype, seurat_clusters, nCount_Akoya, nFeature_Akoya
  )])
  rownames(meta) <- meta$cell_ID
  xy <- as.data.frame(analysis_meta[cell_ID %chin% cells, .(x, y)])
  rownames(xy) <- cells

  patient_result <- tryCatch({
    g <- createGiottoObject(
      expression = expr,
      expression_feat = "protein",
      spatial_locs = xy,
      cell_metadata = meta,
      verbose = FALSE
    )
    g <- createSpatialNetwork(
      gobject = g,
      name = "Delaunay_network",
      method = "Delaunay",
      maximum_distance_delaunay = 50,
      verbose = FALSE
    )
    cp <- cellProximityEnrichment(
      gobject = g,
      spatial_network_name = "Delaunay_network",
      cluster_column = "Celltype",
      number_of_simulations = n_sim,
      adjust_method = "BH",
      set_seed = TRUE,
      seed_number = 20260710
    )
    saveRDS(g, file.path(object_dir, paste0(pid, "_Giotto.rds")))
    list(g = g, cp = cp, error = NULL)
  }, error = function(e) list(g = NULL, cp = NULL, error = conditionMessage(e)))

  if (!is.null(patient_result$error)) {
    log_msg("Giotto patient", pid, "FAILED:", patient_result$error)
    next
  }

  cp_dt <- as.data.table(patient_result$cp$enrichm_res)
  cp_dt[, PatientID := pid]
  cp_dt <- merge(cp_dt, sample_meta, by = "PatientID", all.x = TRUE)
  cpe_results[[pid]] <- cp_dt
  original_edges <- as.data.table(patient_result$cp$raw_sim_table)[round == "original", sum(V1)]
  network_summary[[pid]] <- data.table(
    PatientID = pid,
    n_cells = length(cells),
    n_undirected_edges = original_edges,
    mean_edges_per_cell = 2 * original_edges / length(cells),
    n_interaction_types = nrow(cp_dt),
    simulations = n_sim
  )
  rm(patient_result, g, cp_dt, expr, meta, xy)
  invisible(gc(FALSE))
  log_msg("Giotto patient", pid, "completed")
}

if (!length(cpe_results)) stop("No patient Giotto analyses completed successfully")
cpe <- rbindlist(cpe_results, fill = TRUE)
network_summary <- rbindlist(network_summary, fill = TRUE)
network_summary <- merge(network_summary, sample_meta, by = "PatientID", all.x = TRUE)
write_csv(network_summary, "giotto_spatial_network_summary.csv")
write_csv(cpe, "giotto_cell_proximity_enrichment_by_patient.csv")

# Parse interaction names and create patient-level statistical comparisons.
pair_parts <- tstrsplit(as.character(cpe$unified_int), "--", fixed = TRUE)
cpe[, celltype_1 := pair_parts[[1]]]
cpe[, celltype_2 := pair_parts[[2]]]
cpe[, Diagnosis2 := factor(Diagnosis2, levels = names(group_colors))]
cpe_tests <- pairwise_tests(cpe, "enrichm", c("celltype_1", "celltype_2"))
write_csv(cpe_tests, "statistics_cell_proximity_pairwise_wilcoxon.csv")

# Group-average proximity enrichment heatmaps.
cpe_group <- cpe[, .(
  mean_log2_enrichment = mean(enrichm, na.rm = TRUE),
  median_log2_enrichment = median(enrichm, na.rm = TRUE),
  n_patients = uniqueN(PatientID)
), by = .(Diagnosis2, celltype_1, celltype_2)]
cpe_group_sym <- rbind(
  cpe_group,
  cpe_group[celltype_1 != celltype_2, .(
    Diagnosis2, celltype_1 = celltype_2, celltype_2 = celltype_1,
    mean_log2_enrichment, median_log2_enrichment, n_patients
  )]
)
write_csv(cpe_group, "giotto_cell_proximity_enrichment_group_summary.csv")

p3 <- ggplot(cpe_group_sym, aes(x = celltype_1, y = celltype_2, fill = mean_log2_enrichment)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~ Diagnosis2, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-2, 2), oob = squish, na.value = "grey90"
  ) +
  coord_equal() +
  labs(
    x = NULL, y = NULL,
    title = "Giotto cell proximity enrichment",
    subtitle = paste0("Patient-level Delaunay networks; mean log2 enrichment; ", n_sim, " label permutations per patient"),
    fill = "Mean log2\nenrichment"
  ) +
  theme_codex(8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1),
    strip.text = element_text(face = "bold", size = 10)
  )
save_plot_both(p3, "Figure_03_cell_proximity_heatmaps", 15, 6)

# Top differential spatial interactions shown as patient-level distributions.
top_interactions <- cpe_tests[is.finite(p_value)][order(p_adj_BH, -abs(median_difference))]
top_interactions[, interaction := paste(celltype_1, celltype_2, sep = " x ")]
top_interactions <- unique(top_interactions[, .(interaction, p_adj_BH, max_abs_difference = abs(median_difference))], by = "interaction")
setorder(top_interactions, p_adj_BH, -max_abs_difference)
top_interactions <- head(top_interactions$interaction, 12L)
cpe[, interaction := paste(celltype_1, celltype_2, sep = " x ")]
plot_cpe <- cpe[interaction %in% top_interactions]
plot_cpe[, interaction := factor(interaction, levels = rev(top_interactions))]

if (nrow(plot_cpe)) {
  p4 <- ggplot(plot_cpe, aes(x = Diagnosis2, y = enrichm, fill = Diagnosis2)) +
    geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.55) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.65) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.9, stroke = 0.3) +
    facet_wrap(~ interaction, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      x = NULL, y = "Giotto log2 proximity enrichment",
      title = "Top group-differential spatial interactions",
      subtitle = "Each point is one independent patient",
      fill = "Tissue group"
    ) +
    theme_codex(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  save_plot_both(p4, "Figure_04_differential_spatial_interactions", 13, 12)
} else {
  log_msg("Differential interaction figure skipped: quick subset does not contain multiple tissue groups")
}

# Spatial marker maps for the largest patient sample from each tissue group.
representatives <- sample_counts_plot[, .SD[which.max(n_cells)], by = Diagnosis2]$PatientID
marker_candidates <- paste0(
  "Cell..",
  c("CD3", "CD4", "CD8", "CD20", "CD68", "CD163", "EpCAM", "aSMA", "CD31", "Ki67", "FoxP3", "HLA-DR"),
  ".mean"
)
markers <- intersect(marker_candidates, rownames(seu))
if (length(markers) >= 3L) {
  marker_dt <- rbindlist(lapply(representatives, function(pid) {
    cells <- analysis_meta[PatientID == pid, cell_ID]
    mat <- as.matrix(expr_all[markers, cells, drop = FALSE])
    long <- as.data.table(as.table(mat))
    setnames(long, c("marker", "cell_ID", "expression"))
    long[, cell_ID := as.character(cell_ID)]
    long <- merge(long, analysis_meta[cell_ID %chin% cells, .(cell_ID, PatientID, Diagnosis2, x, y)], by = "cell_ID")
    long
  }))
  marker_dt[, marker_label := sub("\\.mean$", "", sub("^Cell\\.\\.", "", as.character(marker)))]
  marker_dt[, marker_label := factor(marker_label, levels = sub("\\.mean$", "", sub("^Cell\\.\\.", "", markers)))]
  marker_dt[, display_expression := pmin(expression, quantile(expression, 0.99, na.rm = TRUE)), by = marker]
  p5 <- ggplot(marker_dt, aes(x = x, y = y, colour = display_expression)) +
    geom_point(size = 0.36, alpha = 0.95, stroke = 0) +
    facet_grid(PatientID ~ marker_label, scales = "free") +
    scale_colour_viridis_c(option = "magma", trans = "sqrt") +
    scale_y_reverse() +
    labs(
      x = NULL, y = NULL,
      title = "Representative spatial protein maps",
      subtitle = "Largest cell-count patient in each tissue group; colour capped at marker-specific 99th percentile",
      colour = "Protein\nintensity"
    ) +
    theme_void(base_size = 7) +
    theme(
      strip.text = element_text(face = "bold", size = 7),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9),
      legend.position = "right",
      aspect.ratio = 1
    )
  save_plot_both(p5, "Figure_05_spatial_marker_maps", max(12, length(markers) * 1.45), 7)
}

# Save a lightweight manifest and session record.
manifest <- data.table(
  item = c(
    "input_object", "cells_input", "features", "patients_input",
    "patients_analyzed", "excluded_exact_duplicate", "Giotto_version",
    "Delaunay_maximum_distance", "CPE_permutations", "random_seed"
  ),
  value = c(
    normalizePath(input_file), ncol(seu), nrow(seu), uniqueN(cell_meta$PatientID),
    length(analysis_patients), paste(excluded_patients, collapse = ","),
    as.character(packageVersion("Giotto")), 50, n_sim, 20260710
  )
)
write_csv(manifest, "analysis_manifest.csv")
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

log_msg("Analysis complete. Results:", normalizePath(output_dir))
