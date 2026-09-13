#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, giotto.no_python_warn = TRUE)
set.seed(20260710)

suppressPackageStartupMessages({
  library(Seurat)
  library(Giotto)
  library(data.table)
  library(ggplot2)
})

input_file <- "UCCODEX1_Annotated_repaired_withUMAP_neutrophil0.8Matched.rds"
output_dir <- "giotto_codex_results/neutrophil_fibroblast_subtypes"
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "giotto_objects")
log_dir <- file.path(output_dir, "logs")
figure_dir <- "giotto_codex_results/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
uv_cache_dir <- file.path(output_dir, ".uv_cache")
dir.create(uv_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(UV_CACHE_DIR = normalizePath(uv_cache_dir, winslash = "/", mustWork = TRUE))
log_file <- file.path(log_dir, "analysis_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_wilcox <- function(value, group, a, b) {
  xa <- value[group == a & is.finite(value)]
  xb <- value[group == b & is.finite(value)]
  if (length(xa) < 2L || length(xb) < 2L) return(NA_real_)
  tryCatch(wilcox.test(xa, xb, exact = FALSE)$p.value, error = function(e) NA_real_)
}

pairwise_tests <- function(dt) {
  contrasts <- list(
    c("UC_Noninflamed", "Control"),
    c("UC_Inflamed", "Control"),
    c("UC_Inflamed", "UC_Noninflamed")
  )
  out <- list()
  idx <- 0L
  for (subtype in unique(dt$Neutrophil_subtype)) {
    sub <- dt[Neutrophil_subtype == subtype]
    for (contrast in contrasts) {
      idx <- idx + 1L
      a <- contrast[1]
      b <- contrast[2]
      va <- sub[Diagnosis2 == a, enrichm]
      vb <- sub[Diagnosis2 == b, enrichm]
      out[[idx]] <- data.table(
        Neutrophil_subtype = subtype,
        contrast = paste(a, "vs", b),
        group_1 = a,
        group_2 = b,
        n_group_1 = sum(is.finite(va)),
        n_group_2 = sum(is.finite(vb)),
        median_group_1 = median(va, na.rm = TRUE),
        median_group_2 = median(vb, na.rm = TRUE),
        median_difference = median(va, na.rm = TRUE) - median(vb, na.rm = TRUE),
        p_value = safe_wilcox(sub$enrichm, sub$Diagnosis2, a, b)
      )
    }
  }
  ans <- rbindlist(out)
  ans[, p_adj_BH := p.adjust(p_value, method = "BH"), by = contrast]
  ans[]
}

log_msg("Loading", input_file)
seu <- readRDS(input_file)
stopifnot(all(c("Celltype", "Celltype_updated") %in% colnames(seu[[]])))

coord_list <- lapply(Images(seu), function(image_name) {
  x <- as.data.table(GetTissueCoordinates(seu, image = image_name))
  setnames(x, "cell", "cell_ID")
  x
})
coords <- rbindlist(coord_list)
stopifnot(nrow(coords) == ncol(seu), !anyDuplicated(coords$cell_ID))

meta <- as.data.table(seu[[]], keep.rownames = "cell_ID")
meta <- merge(meta, coords, by = "cell_ID", all.x = TRUE, sort = FALSE)
setkey(meta, cell_ID)
meta <- meta[J(colnames(seu))]
sample_meta <- unique(meta[, .(PatientID, Diagnosis1, Diagnosis2, Inflammation, Medication)])
patients <- sort(unique(meta$PatientID))
expr_one <- LayerData(seu, assay = "Akoya", layer = "data")[1, , drop = FALSE]

subtype_labels <- c(
  "Neutrophil_MX1", "Neutrophil_CXCR4", "Neutrophil_PADI4", "Neutrophil_OSM", "Neutrophil"
)
cpe_results <- list()
network_summary <- list()

for (pid in patients) {
  log_msg("Patient", pid, "started")
  cells <- meta[PatientID == pid, cell_ID]
  patient_meta <- as.data.frame(meta[cell_ID %chin% cells, .(
    cell_ID, PatientID, Diagnosis1, Diagnosis2, Inflammation, Medication,
    Celltype, Celltype_updated
  )])
  rownames(patient_meta) <- patient_meta$cell_ID
  xy <- as.data.frame(meta[cell_ID %chin% cells, .(x, y)])
  rownames(xy) <- cells

  result <- tryCatch({
    g <- createGiottoObject(
      expression = expr_one[, cells, drop = FALSE],
      expression_feat = "protein",
      spatial_locs = xy,
      cell_metadata = patient_meta,
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
      cluster_column = "Celltype_updated",
      number_of_simulations = 100,
      adjust_method = "BH",
      set_seed = TRUE,
      seed_number = 20260710
    )
    saveRDS(g, file.path(object_dir, paste0(pid, "_Giotto_subtypes.rds")))
    list(cp = cp, error = NULL)
  }, error = function(e) list(cp = NULL, error = conditionMessage(e)))

  if (!is.null(result$error)) {
    log_msg("Patient", pid, "FAILED:", result$error)
    next
  }
  cp_dt <- as.data.table(result$cp$enrichm_res)
  cp_dt[, PatientID := pid]
  parts <- tstrsplit(as.character(cp_dt$unified_int), "--", fixed = TRUE)
  cp_dt[, celltype_1 := parts[[1]]]
  cp_dt[, celltype_2 := parts[[2]]]
  keep <- (
    (cp_dt$celltype_1 == "Fibroblast" & cp_dt$celltype_2 %in% subtype_labels) |
      (cp_dt$celltype_2 == "Fibroblast" & cp_dt$celltype_1 %in% subtype_labels)
  )
  nf <- cp_dt[keep]
  if (nrow(nf)) {
    nf[, Neutrophil_subtype := fifelse(celltype_1 == "Fibroblast", celltype_2, celltype_1)]
    nf <- merge(nf, sample_meta, by = "PatientID", all.x = TRUE)
    cpe_results[[pid]] <- nf
  }
  network_summary[[pid]] <- data.table(
    PatientID = pid,
    n_cells = length(cells),
    n_neutrophil_fibroblast_interactions = nrow(nf)
  )
  log_msg("Patient", pid, "completed with", nrow(nf), "subtype-fibroblast interaction types")
  rm(result, patient_meta, xy, cp_dt, nf)
  invisible(gc(FALSE))
}

if (!length(cpe_results)) stop("No subtype-fibroblast results were generated")
cpe <- rbindlist(cpe_results, fill = TRUE)
cpe[, Diagnosis2 := factor(Diagnosis2, levels = c("Control", "UC_Noninflamed", "UC_Inflamed"))]
tests <- pairwise_tests(cpe)
group_summary <- cpe[, .(
  n_patients = uniqueN(PatientID),
  median_log2_enrichment = median(enrichm, na.rm = TRUE),
  mean_log2_enrichment = mean(enrichm, na.rm = TRUE)
), by = .(Diagnosis2, Neutrophil_subtype)]

fwrite(cpe, file.path(table_dir, "neutrophil_subtype_fibroblast_enrichment_by_patient.csv"))
fwrite(tests, file.path(table_dir, "statistics_neutrophil_subtype_fibroblast_pairwise_wilcoxon.csv"))
fwrite(group_summary, file.path(table_dir, "neutrophil_subtype_fibroblast_group_summary.csv"))
fwrite(rbindlist(network_summary), file.path(table_dir, "patient_network_summary.csv"))

display_labels <- c(
  "Neutrophil_MX1" = "MX1",
  "Neutrophil_CXCR4" = "CXCR4",
  "Neutrophil_PADI4" = "PADI4",
  "Neutrophil_OSM" = "OSM",
  "Neutrophil" = "Generic"
)
cpe[, display_subtype := factor(display_labels[Neutrophil_subtype], levels = unname(display_labels))]
group_colors <- c("Control" = "#4C78A8", "UC_Noninflamed" = "#F2CF5B", "UC_Inflamed" = "#E45756")
p <- ggplot(cpe, aes(Diagnosis2, enrichm, fill = Diagnosis2)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.55) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.3, alpha = 0.9, stroke = 0.35) +
  facet_wrap(~ display_subtype, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "Giotto log2 proximity enrichment with fibroblasts",
    title = "Neutrophil subtype-fibroblast spatial interactions",
    subtitle = "Each point is one patient; Delaunay networks with 100 label permutations",
    fill = "Tissue group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(file.path(figure_dir, "Figure_07_neutrophil_subtype_fibroblast_interactions.png"), p,
       width = 15, height = 5, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "Figure_07_neutrophil_subtype_fibroblast_interactions.pdf"), p,
       width = 15, height = 5, bg = "white")

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
log_msg("Analysis complete")
print(tests[order(p_adj_BH, -abs(median_difference))])
