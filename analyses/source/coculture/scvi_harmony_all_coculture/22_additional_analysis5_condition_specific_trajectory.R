suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(condiments)
  library(limma)
  library(splines)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260822)
project_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
base_analysis <- file.path(project_dir, "Manuscript_Analyses_1_to_6")
out_dir <- file.path(project_dir, "Additional_High_Impact_Analyses_2_to_5")
source_dir <- file.path(out_dir, "source_data")
figure_dir <- file.path(out_dir, "figures")
object_dir <- file.path(out_dir, "objects")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

save_figure <- function(p, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 400,
         device = ragg::agg_png, bg = "white")
  ggsave(paste0(stem, ".tiff"), p, width = width, height = height, dpi = 400,
         device = ragg::agg_tiff, compression = "lzw", bg = "white")
}

theme_pub <- theme_bw(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, color = "#444444"),
        panel.grid.minor = element_blank(), axis.text = element_text(size = 7.3),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7), strip.text = element_text(face = "bold", size = 8),
        plot.margin = margin(6, 6, 6, 6))

obj <- readRDS(file.path(project_dir,
  "scVI_Harmony_All_Coculture/Neutrophil_Subclustering/objects/canonical_neutrophils_scVI_Harmony_reclustered.rds"))
traj <- fread(file.path(base_analysis, "source_data", "analysis6_cell_pseudotime_and_lineage_weights.csv.gz"))
traj <- traj[is.finite(weighted_pseudotime)]

comparisons <- data.table(
  contrast_id = c("N_vs_PRE", "CF_vs_N", "UF_vs_N", "UN_vs_UF", "UA5_vs_UF", "UD_vs_UF"),
  reference = c("PRE", "N", "N", "UF", "UF", "UF"),
  comparison = c("N", "CF", "UF", "UN", "UA5", "UD"),
  paired = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  contrast = c("N vs PRE", "Control FB vs N", "UC FB vs N", "NAMPTi vs UC",
               "alpha5beta1i vs UC", "Dual vs UC")
)

# Condiments progression sensitivity analysis on donor-balanced cells.
balanced_list <- list()
condiment_rows <- list()
for (i in seq_len(nrow(comparisons))) {
  key <- comparisons[i]
  d <- traj[ConditionCode %in% c(key$reference, key$comparison)]
  if (key$paired) {
    donors <- d[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 2, DonorID]
    d <- d[DonorID %in% donors]
    selected <- d[, {
      n_take <- min(.N, 250L)
      .(cell = sample(cell, n_take))
    }, by = .(DonorID, ConditionCode)]$cell
  } else {
    selected <- d[, .(cell = sample(cell, min(.N, 250L))), by = .(DonorID, ConditionCode)]$cell
  }
  b <- d[cell %in% selected]
  b[, contrast_id := key$contrast_id]
  balanced_list[[i]] <- b
  pt <- matrix(b$weighted_pseudotime, ncol = 1,
               dimnames = list(b$cell, "Lineage1"))
  cw <- matrix(1, nrow = nrow(b), ncol = 1,
               dimnames = list(b$cell, "Lineage1"))
  z <- as.data.table(progressionTest(pt, cellWeights = cw,
                                    conditions = factor(b$ConditionCode,
                                                        levels = c(key$reference, key$comparison)),
                                    global = TRUE, lineages = FALSE, rep = 2000))
  condiment_rows[[i]] <- data.table(contrast_id = key$contrast_id, contrast = key$contrast,
                                    balanced_cells = nrow(b), statistic = z$statistic[1],
                                    p_value = z$p.value[1])
}
balanced <- rbindlist(balanced_list)
condiments_results <- rbindlist(condiment_rows)
condiments_results[, FDR := p.adjust(p_value, method = "BH")]
fwrite(condiments_results, file.path(source_dir, "analysis5_condiments_balanced_progression_tests.csv"))

# Primary donor-level inference across multiple pseudotime quantiles.
quantile_probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
quantiles <- traj[, as.list(setNames(quantile(weighted_pseudotime, probs = quantile_probs, na.rm = TRUE),
                                     paste0("Q", quantile_probs * 100))),
                  by = .(DonorID, ConditionCode)]
quantile_long <- melt(quantiles, id.vars = c("DonorID", "ConditionCode"),
                      variable.name = "quantile", value.name = "pseudotime")

quantile_effects <- rbindlist(lapply(seq_len(nrow(comparisons)), function(i) {
  key <- comparisons[i]
  d <- quantile_long[ConditionCode %in% c(key$reference, key$comparison)]
  if (key$paired) {
    donors <- d[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 2, DonorID]
    d <- d[DonorID %in% donors]
  }
  rbindlist(lapply(unique(d$quantile), function(q) {
    z <- d[quantile == q]
    cond <- relevel(factor(z$ConditionCode), key$reference)
    design <- if (key$paired) model.matrix(~ factor(z$DonorID) + cond) else model.matrix(~ cond)
    fit <- eBayes(lmFit(matrix(z$pseudotime, nrow = 1), design))
    se <- fit$stdev.unscaled[1, ncol(design)] * sqrt(fit$s2.post[1])
    data.table(contrast_id = key$contrast_id, contrast = key$contrast, quantile = q,
               effect = fit$coefficients[1, ncol(design)], SE = se,
               p_value = fit$p.value[1, ncol(design)], n_donors = uniqueN(z$DonorID))
  }))
}), fill = TRUE)
quantile_effects[, FDR := p.adjust(p_value, method = "BH"), by = quantile]
fwrite(quantile_effects, file.path(source_dir, "analysis5_donor_pseudotime_quantile_effects.csv"))

# Donor-condition-pseudotime-bin expression profiles for condition-specific curves.
assoc <- fread(file.path(base_analysis, "source_data", "analysis6_tradeSeq_trajectory_association.csv"))
setorder(assoc, FDR, pvalue)
candidate_genes <- unique(c(head(assoc$gene, 70),
  "PADI4", "OSM", "CXCR4", "MX1", "IL1B", "CXCL8", "CXCL2", "CCL3", "CCL4",
  "LTF", "BPI", "DEFA3", "S100A8", "S100A9", "CD24"))
rna <- LayerData(obj, assay = "RNA", layer = "data")
candidate_genes <- intersect(candidate_genes, rownames(rna))
traj[, pt_bin := cut(weighted_pseudotime,
  breaks = unique(quantile(weighted_pseudotime, probs = seq(0, 1, length.out = 9), na.rm = TRUE)),
  include.lowest = TRUE, labels = FALSE)]
traj <- traj[!is.na(pt_bin)]
traj[, curve_group := paste(DonorID, ConditionCode, pt_bin, sep = "||")]
curve_levels <- unique(traj$curve_group)
indicator <- sparseMatrix(i = seq_len(nrow(traj)), j = match(traj$curve_group, curve_levels), x = 1,
                          dims = c(nrow(traj), length(curve_levels)),
                          dimnames = list(traj$cell, curve_levels))
x <- rna[candidate_genes, traj$cell, drop = FALSE]
curve_mat <- as.matrix(x %*% indicator) /
  rep(as.numeric(table(factor(traj$curve_group, levels = curve_levels))), each = length(candidate_genes))
curve_meta <- data.table(curve_group = curve_levels)
curve_meta[, c("DonorID", "ConditionCode", "pt_bin") := tstrsplit(curve_group, "||", fixed = TRUE)]
curve_meta[, pt_bin := as.numeric(pt_bin)]

curve_results <- list()
curve_summaries <- list()
k <- 0L
for (i in which(comparisons$paired)) {
  key <- comparisons[i]
  idx <- which(curve_meta$ConditionCode %in% c(key$reference, key$comparison))
  sm <- copy(curve_meta[idx])
  donors <- sm[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 2, DonorID]
  keep_idx <- idx[sm$DonorID %in% donors]
  sm <- curve_meta[keep_idx]
  if (uniqueN(sm$DonorID) < 3) next
  cond <- relevel(factor(sm$ConditionCode), key$reference)
  design <- model.matrix(~ factor(sm$DonorID) + cond * ns(sm$pt_bin, df = 3))
  fit <- eBayes(lmFit(curve_mat[, keep_idx, drop = FALSE], design))
  test_cols <- grep("^cond|cond.*:", colnames(design))
  tab <- as.data.table(topTable(fit, coef = test_cols, number = Inf, sort.by = "none"), keep.rownames = "gene")
  mean_curves <- rbindlist(lapply(c(key$reference, key$comparison), function(cc) {
    cols <- keep_idx[curve_meta$ConditionCode[keep_idx] == cc]
    zz <- curve_meta[cols]
    rbindlist(lapply(sort(unique(zz$pt_bin)), function(bb) {
      use <- cols[curve_meta$pt_bin[cols] == bb]
      data.table(gene = candidate_genes, ConditionCode = cc, pt_bin = bb,
                 mean_expression = rowMeans(curve_mat[, use, drop = FALSE]),
                 se_expression = apply(curve_mat[, use, drop = FALSE], 1, sd) / sqrt(length(use)))
    }))
  }))
  auc <- dcast(mean_curves[, .(mean_expression = mean(mean_expression)),
                           by = .(gene, ConditionCode, pt_bin)], gene + pt_bin ~ ConditionCode,
               value.var = "mean_expression")
  auc[, bin_delta := get(key$comparison) - get(key$reference)]
  auc_gene <- auc[, .(curve_AUC_difference = mean(bin_delta, na.rm = TRUE)), by = gene]
  tab <- merge(tab, auc_gene, by = "gene", all.x = TRUE)
  tab[, `:=`(contrast_id = key$contrast_id, contrast = key$contrast,
             reference = key$reference, comparison = key$comparison,
             n_complete_donors = uniqueN(sm$DonorID))]
  k <- k + 1L; curve_results[[k]] <- tab
  mean_curves[, `:=`(contrast_id = key$contrast_id, contrast = key$contrast)]
  curve_summaries[[k]] <- mean_curves
}
curve_tests <- rbindlist(curve_results, fill = TRUE)
curve_tests[, FDR := p.adjust(P.Value, method = "BH"), by = contrast_id]
curve_summary <- rbindlist(curve_summaries, fill = TRUE)
fwrite(curve_tests, file.path(source_dir, "analysis5_condition_specific_pseudotime_gene_curve_tests.csv.gz"))
fwrite(curve_summary, file.path(source_dir, "analysis5_condition_specific_pseudotime_gene_curves.csv.gz"))

contrast_levels <- comparisons$contrast
condition_labels <- c(PRE = "PRE", N = "N", CF = "Control FB", UF = "UC FB",
                      UN = "UC+NAMPTi", UA5 = "UC+alpha5beta1i", UD = "Dual")

density_data <- balanced[contrast_id %in% c("UF_vs_N", "UN_vs_UF", "UA5_vs_UF", "UD_vs_UF")]
density_data <- merge(density_data, comparisons[, .(contrast_id, contrast)], by = "contrast_id")
density_data[, condition := factor(condition_labels[ConditionCode])]
density_data[, contrast := factor(contrast, levels = contrast_levels)]
p_a <- ggplot(density_data, aes(weighted_pseudotime, color = condition, fill = condition)) +
  geom_density(alpha = 0.13, linewidth = 0.65) + facet_wrap(~contrast, ncol = 2, scales = "free_y") +
  labs(title = "Condition-specific progression along the common neutrophil trajectory",
       subtitle = "Donor-balanced cells used for the Condiments progression sensitivity test",
       x = "Pseudotime", y = "Density", color = NULL, fill = NULL) + theme_pub +
  theme(legend.position = "bottom")

quantile_effects[, contrast := factor(contrast, levels = contrast_levels)]
quantile_effects[, quantile := factor(quantile, levels = paste0("Q", quantile_probs * 100))]
p_b <- ggplot(quantile_effects, aes(contrast, quantile, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.2f%s", effect, ifelse(FDR < 0.05, "*", ""))), size = 2.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Donor-level shifts across the pseudotime distribution",
       subtitle = "Positive values indicate later progression; * FDR <0.05",
       x = NULL, y = "Pseudotime quantile", fill = "Effect") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

curve_priority <- curve_tests[, .(priority = max(abs(curve_AUC_difference) *
                                                   (1 + -log10(pmax(FDR, 1e-12))))), by = gene]
top_curve_genes <- head(curve_priority[order(-priority), gene], 24)
curve_plot <- curve_tests[gene %in% top_curve_genes]
curve_plot[, gene := factor(gene, levels = rev(top_curve_genes))]
curve_plot[, contrast := factor(contrast, levels = contrast_levels)]
p_c <- ggplot(curve_plot, aes(contrast, gene, fill = curve_AUC_difference)) +
  geom_tile(color = "white", linewidth = 0.22) +
  geom_point(data = curve_plot[FDR < 0.05], shape = 8, size = 1.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Condition-dependent gene-expression trajectories",
       subtitle = "Donor-blocked spline tests; color is signed curve-AUC difference; * FDR <0.05",
       x = NULL, y = NULL, fill = "Curve AUC") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.5))

alpha_top <- head(curve_tests[contrast_id == "UA5_vs_UF"][order(FDR, -abs(curve_AUC_difference)), gene], 6)
curve_lines <- curve_summary[contrast_id == "UA5_vs_UF" & gene %in% alpha_top]
curve_lines[, gene := factor(gene, levels = alpha_top)]
curve_lines[, condition := factor(condition_labels[ConditionCode],
                                  levels = condition_labels[c("UF", "UA5")])]
p_d <- ggplot(curve_lines, aes(pt_bin, mean_expression, color = condition, group = condition)) +
  geom_ribbon(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression,
                  fill = condition), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(title = "Genes with altered progression after alpha5beta1 blockade",
       subtitle = "Mean donor-bin expression with standard-error bands",
       x = "Pseudotime bin", y = "Mean log-normalized expression", color = NULL, fill = NULL) +
  theme_pub + theme(legend.position = "bottom")

fig <- p_a / p_b / p_c / p_d + plot_layout(heights = c(0.95, 0.70, 1.0, 1.0)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(fig, file.path(figure_dir, "Additional_Analysis_5_condition_specific_trajectory"),
            width = 16, height = 22)

saveRDS(list(condiments = condiments_results, donor_quantiles = quantile_effects,
             gene_curve_tests = curve_tests, gene_curves = curve_summary),
        file.path(object_dir, "additional_analysis5_condition_specific_trajectory.rds"), compress = TRUE)
message("Additional analysis 5 complete")
