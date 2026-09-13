suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(Matrix)
  library(data.table)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

set.seed(160826)
project_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
analysis_dir <- file.path(project_dir, "Manuscript_Analyses_1_to_6")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
object_dir <- file.path(analysis_dir, "objects")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(project_dir, "scVI_Harmony_All_Coculture/Neutrophil_Subclustering/objects/canonical_neutrophils_scVI_Harmony_reclustered.rds"))
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
meta[, state := as.character(neutrophil_state_annotation)]
state_order <- c("PADI4 neutrophil", "OSM neutrophil", "CXCR4 neutrophil", "MX1/ISG neutrophil")
meta[, state := factor(state, levels = state_order)]

full_embedding <- Embeddings(obj, "scvi.harmony")
meta[, trajectory_stratum := paste(ConditionCode, state, sep = "||")]
landmark_cells <- meta[, .(cell = sample(cell, min(.N, 250L))), by = trajectory_stratum]$cell
landmark_idx <- match(landmark_cells, meta$cell)
sds_landmark <- slingshot(
  full_embedding[landmark_idx, , drop = FALSE],
  clusterLabels = droplevels(meta$state[landmark_idx]),
  start.clus = "PADI4 neutrophil", allow.breaks = FALSE, shrink = TRUE
)
sds_full <- predict(sds_landmark, newdata = full_embedding)
pt <- slingPseudotime(sds_full, na = TRUE)
cw <- slingCurveWeights(sds_full)
lineages <- slingLineages(sds_landmark)
colnames(pt) <- paste0("Lineage", seq_len(ncol(pt)))
colnames(cw) <- colnames(pt)
weighted_pt <- rowSums(pt * cw, na.rm = TRUE) / pmax(rowSums(cw, na.rm = TRUE), 1e-8)
weighted_pt[rowSums(cw, na.rm = TRUE) == 0] <- NA_real_
dominant_lineage <- colnames(cw)[max.col(cw, ties.method = "first")]
trajectory_cells <- cbind(
  meta[, .(cell, DonorID, ConditionCode, state)],
  as.data.table(pt),
  as.data.table(setNames(as.data.frame(cw), paste0(colnames(cw), "_weight")))
)
trajectory_cells[, `:=`(weighted_pseudotime = weighted_pt, dominant_lineage = dominant_lineage)]
fwrite(trajectory_cells, file.path(source_dir, "analysis6_cell_pseudotime_and_lineage_weights.csv.gz"))

lineage_definition <- rbindlist(lapply(seq_along(lineages), function(i) {
  data.table(lineage = paste0("Lineage", i), order = seq_along(lineages[[i]]), state = lineages[[i]])
}))
fwrite(lineage_definition, file.path(source_dir, "analysis6_lineage_definitions.csv"))

sample_summary <- trajectory_cells[, c(
  list(cells = .N, median_pseudotime = median(weighted_pseudotime, na.rm = TRUE), mean_pseudotime = mean(weighted_pseudotime, na.rm = TRUE)),
  lapply(.SD, function(x) mean(x, na.rm = TRUE))
), by = .(DonorID, ConditionCode), .SDcols = patterns("_weight$")]
setnames(sample_summary, grep("_weight$", names(sample_summary), value = TRUE),
         sub("_weight$", "_occupancy", grep("_weight$", names(sample_summary), value = TRUE)))
sample_summary[, sample_id := paste(DonorID, ConditionCode, sep = "||")]
fwrite(sample_summary, file.path(source_dir, "analysis6_trajectory_summary_by_sample.csv"))

contrast_key <- fread(file.path(source_dir, "analysis_contrasts.csv"))
did_key <- fread(file.path(source_dir, "analysis_interaction_contrasts.csv"))
metric_cols <- "median_pseudotime"
if (length(lineages) > 1) metric_cols <- c(metric_cols, grep("_occupancy$", names(sample_summary), value = TRUE))
metric_long <- melt(sample_summary, id.vars = c("sample_id", "DonorID", "ConditionCode", "cells"),
                    measure.vars = metric_cols, variable.name = "metric", value.name = "value")

run_pair <- function(key_row) {
  z <- metric_long[ConditionCode %in% c(key_row$reference, key_row$comparison)]
  if (key_row$paired) {
    donors <- z[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 2, DonorID]
    z <- z[DonorID %in% donors]
  }
  if (uniqueN(z$DonorID) < 3) return(NULL)
  rbindlist(lapply(unique(z$metric), function(metric_name) {
    d <- z[metric == metric_name]
    cond <- relevel(factor(d$ConditionCode), key_row$reference)
    design <- if (key_row$paired) model.matrix(~ factor(d$DonorID) + cond) else model.matrix(~ cond)
    fit <- lmFit(matrix(d$value, nrow = 1), design)
    tab <- topTable(eBayes(fit), coef = ncol(design), number = 1, sort.by = "none")
    data.table(metric = metric_name, effect = tab$logFC, p_value = tab$P.Value,
               contrast_id = key_row$contrast_id, contrast = key_row$contrast,
               n_pairs = if (key_row$paired) uniqueN(d$DonorID) else NA_integer_)
  }))
}

run_did <- function(key_row) {
  conditions <- strsplit(key_row$conditions, ",", fixed = TRUE)[[1]]
  z <- metric_long[ConditionCode %in% conditions]
  donors <- z[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 4, DonorID]
  z <- z[DonorID %in% donors]
  if (uniqueN(z$DonorID) < 3) return(NULL)
  rbindlist(lapply(unique(z$metric), function(metric_name) {
    d <- z[metric == metric_name]
    d$ConditionCode <- factor(d$ConditionCode, levels = conditions)
    design <- model.matrix(~ 0 + ConditionCode + factor(DonorID), data = d)
    colnames(design) <- make.names(colnames(design))
    if (qr(design)$rank < ncol(design)) return(NULL)
    contrast_matrix <- makeContrasts(contrasts = key_row$expression, levels = design)
    fit <- eBayes(contrasts.fit(lmFit(matrix(d$value, nrow = 1), design), contrast_matrix))
    tab <- topTable(fit, coef = 1, number = 1, sort.by = "none")
    data.table(metric = metric_name, effect = tab$logFC, p_value = tab$P.Value,
               contrast_id = key_row$contrast_id, contrast = key_row$contrast,
               n_pairs = uniqueN(d$DonorID))
  }))
}

trajectory_results <- rbindlist(c(
  lapply(seq_len(nrow(contrast_key)), function(i) run_pair(contrast_key[i])),
  lapply(seq_len(nrow(did_key)), function(i) run_did(did_key[i]))
), fill = TRUE)
trajectory_results[, FDR := p.adjust(p_value, method = "BH"), by = metric]
fwrite(trajectory_results, file.path(source_dir, "analysis6_differential_pseudotime_and_lineage_occupancy.csv"))

# tradeSeq on a balanced subset and variable genes to identify dynamic transcripts.
counts <- LayerData(obj, assay = "RNA", layer = "counts")
sampling <- trajectory_cells[is.finite(weighted_pseudotime)]
sampling[, stratum := paste(ConditionCode, state, sep = "||")]
selected_cells <- sampling[, {
  n_take <- min(.N, 350L)
  .(cell = sample(cell, n_take))
}, by = stratum]$cell
selected_idx <- match(selected_cells, trajectory_cells$cell)
sub_counts <- counts[, selected_cells, drop = FALSE]
expressed <- rowSums(sub_counts > 0) >= max(50, round(0.01 * ncol(sub_counts)))
mu <- Matrix::rowMeans(sub_counts[expressed, , drop = FALSE])
sq <- Matrix::rowMeans(sub_counts[expressed, , drop = FALSE]^2)
disp <- (sq - mu^2) / pmax(mu, 0.05)
candidate_genes <- names(sort(disp, decreasing = TRUE))[seq_len(min(220, length(disp)))]
candidate_genes <- unique(c(candidate_genes, intersect(c("OSM", "CXCR4", "MX1", "PADI4", "IL1R2", "CCL3", "CCL4", "LTF", "BPI", "RORA"), rownames(counts))))

association_file <- file.path(source_dir, "analysis6_tradeSeq_trajectory_association.csv")
if (file.exists(association_file)) {
  assoc <- fread(association_file)
  gam <- NULL
} else {
  gam <- fitGAM(
    counts = as.matrix(sub_counts[candidate_genes, , drop = FALSE]),
    pseudotime = pt[selected_idx, , drop = FALSE], cellWeights = cw[selected_idx, , drop = FALSE],
    nknots = 5, verbose = FALSE, parallel = FALSE
  )
  assoc <- as.data.table(associationTest(gam), keep.rownames = "gene")
  assoc[, FDR := p.adjust(pvalue, method = "BH")]
  setorder(assoc, FDR, pvalue)
  fwrite(assoc, association_file)
}

# Smoothed-bin expression for interpretable trajectory-gene heatmap.
top_genes <- head(assoc[FDR < 0.05, gene], 36)
if (length(top_genes) < 20) top_genes <- head(assoc$gene, 36)
log_expr <- log1p(t(t(as.matrix(sub_counts[top_genes, , drop = FALSE])) / pmax(Matrix::colSums(sub_counts), 1) * 10000))
plot_pt <- weighted_pt[selected_idx]
bins <- cut(plot_pt, breaks = quantile(plot_pt, probs = seq(0, 1, length.out = 13), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)
bin_expr <- sapply(sort(unique(bins)), function(b) rowMeans(log_expr[, bins == b, drop = FALSE]))
bin_z <- t(scale(t(bin_expr)))
bin_z[!is.finite(bin_z)] <- 0
colnames(bin_z) <- seq_len(ncol(bin_z))
heat_dt <- as.data.table(as.table(bin_z))
setnames(heat_dt, c("gene", "pseudotime_bin", "z_expression"))
heat_dt[, pseudotime_bin := as.integer(pseudotime_bin)]
fwrite(heat_dt, file.path(source_dir, "analysis6_trajectory_gene_binned_expression.csv"))

theme_pub <- theme_bw(base_family = "Arial", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9.5, color = "#444444"),
        axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
        legend.title = element_text(face = "bold"), plot.margin = margin(7, 7, 7, 7))

umap <- as.data.table(Embeddings(obj, "umap.canonical.scvi.harmony"), keep.rownames = "cell")
setnames(umap, c("cell", "UMAP_1", "UMAP_2"))
umap <- merge(umap, trajectory_cells[, .(cell, state, weighted_pseudotime)], by = "cell")
centroids <- umap[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = state]
centroids[, state_label := sub(" neutrophil$", "", as.character(state))]
curve_dt <- rbindlist(lapply(seq_along(lineages), function(i) {
  z <- centroids[match(lineages[[i]], state)]
  z[, `:=`(lineage = paste0("Lineage", i), order = seq_len(.N))]
  z
}))
plot_cells <- umap[sample(.N, min(.N, 25000))]
p_a <- ggplot(plot_cells, aes(UMAP_1, UMAP_2, color = weighted_pseudotime)) +
  geom_point(size = 0.22, alpha = 0.65) +
  geom_path(data = curve_dt, mapping = aes(UMAP_1, UMAP_2, group = lineage), inherit.aes = FALSE,
            color = "black", linewidth = 1.0,
            arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  geom_point(data = centroids, inherit.aes = FALSE, aes(UMAP_1, UMAP_2), shape = 21, fill = "white", color = "black", size = 2.2) +
  geom_text(data = centroids, inherit.aes = FALSE, aes(UMAP_1, UMAP_2, label = state_label),
            nudge_y = 0.22, size = 3, fontface = "bold", check_overlap = TRUE) +
  scale_color_viridis_c(option = "C") + coord_equal() +
  labs(title = "SCVI-Harmony neutrophil trajectory", subtitle = "Slingshot rooted in the PADI4 state; arrows show inferred lineage direction",
       x = "UMAP 1", y = "UMAP 2", color = "Pseudotime") + theme_pub

condition_order <- c("PRE", "N", "NNAMPT", "NA5", "CF", "UF", "UN", "UA5", "UD")
sample_summary[, ConditionCode := factor(ConditionCode, levels = condition_order)]
p_b <- ggplot(sample_summary, aes(ConditionCode, median_pseudotime, group = DonorID)) +
  geom_line(color = "#BDBDBD", linewidth = 0.4, alpha = 0.6) + geom_point(aes(color = ConditionCode), size = 2) +
  geom_boxplot(data = sample_summary, mapping = aes(ConditionCode, median_pseudotime, group = ConditionCode),
               width = 0.5, outlier.shape = NA, fill = NA, color = "black", inherit.aes = FALSE) +
  guides(color = "none") + labs(title = "Donor-level trajectory position", subtitle = "Lines connect the same culture donor; PRE donors are independent",
                                x = NULL, y = "Median pseudotime") + theme_pub

trajectory_results[, contrast := factor(contrast, levels = c(contrast_key$contrast, did_key$contrast))]
trajectory_results[, metric_label := gsub("_", " ", metric)]
p_c_title <- if (length(lineages) > 1) "Differential progression and lineage occupancy" else "Differential inferred progression"
p_c <- ggplot(trajectory_results, aes(contrast, metric_label, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_point(data = trajectory_results[FDR < 0.05], shape = 8, size = 2, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = p_c_title, subtitle = "Donor-aware effects; asterisks indicate FDR <0.05",
       x = NULL, y = NULL, fill = "Effect") + theme_pub

gene_order <- rev(unique(heat_dt$gene))
heat_dt[, gene := factor(gene, levels = gene_order)]
p_d <- ggplot(heat_dt, aes(pseudotime_bin, gene, fill = z_expression)) +
  geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(title = "Genes varying along inferred pseudotime", subtitle = "Top tradeSeq-associated genes; expression is row-scaled across pseudotime bins",
       x = "Pseudotime quantile bin", y = NULL, fill = "Z score") + theme_pub + theme(axis.text.y = element_text(size = 7))

fig6 <- (p_a | p_b) / p_c / p_d + plot_layout(heights = c(1.05, 0.75, 1.15)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
for (ext in c("pdf", "png", "tiff")) {
  output_file <- file.path(figure_dir, paste0("Analysis6_neutrophil_trajectory.", ext))
  if (ext == "tiff") ggsave(output_file, fig6, width = 14, height = 16, dpi = 400, compression = "lzw")
  else ggsave(output_file, fig6, width = 14, height = 16, dpi = 400)
}

saveRDS(list(slingshot_landmark = sds_landmark, slingshot_projected = sds_full, trajectory_results = trajectory_results, tradeSeq = gam, lineage_definitions = lineage_definition),
        file.path(object_dir, "analysis6_neutrophil_trajectory.rds"), compress = TRUE)
message("Analysis 6 neutrophil trajectory complete")
