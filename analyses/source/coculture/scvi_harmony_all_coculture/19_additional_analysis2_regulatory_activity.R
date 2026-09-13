suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(limma)
  library(decoupleR)
  library(progeny)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
})

set.seed(20260819)
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
        panel.grid = element_blank(), axis.text = element_text(size = 7.5),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        strip.text = element_text(face = "bold", size = 8),
        plot.margin = margin(6, 6, 6, 6))

# Build donor-condition-state pseudobulks from raw neutrophil counts.
obj <- readRDS(file.path(project_dir,
  "scVI_Harmony_All_Coculture/Neutrophil_Subclustering/objects/canonical_neutrophils_scVI_Harmony_reclustered.rds"))
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
meta[, state := as.character(neutrophil_state_annotation)]
meta[, sample_state := paste(DonorID, ConditionCode, state, sep = "||")]
group_levels <- unique(meta$sample_state)
indicator <- sparseMatrix(i = seq_len(nrow(meta)), j = match(meta$sample_state, group_levels), x = 1,
                          dims = c(nrow(meta), length(group_levels)),
                          dimnames = list(meta$cell, group_levels))
counts <- LayerData(obj, assay = "RNA", layer = "counts")[, meta$cell, drop = FALSE]
pb_counts <- counts %*% indicator
pb_meta <- data.table(sample_state = group_levels)
pb_meta[, c("DonorID", "ConditionCode", "state") := tstrsplit(sample_state, "||", fixed = TRUE)]
pb_meta[, sample_id := paste(DonorID, ConditionCode, sep = "||")]

keep <- rowSums(pb_counts >= 10) >= 3
dge <- DGEList(pb_counts[keep, , drop = FALSE])
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# PROGENy footprint network: 500 most informative genes per pathway.
prog_net <- as.data.table(model_human_full)
setorder(prog_net, pathway, p.value)
prog_net <- prog_net[, head(.SD, 500), by = pathway]
prog_net <- prog_net[gene %in% rownames(logcpm), .(source = pathway, target = gene, mor = weight)]

# CollecTRI signed TF-target network from the official rescued OmniPath flat file.
tf_net <- fread(file.path(project_dir, "resources_CollecTRI.csv"))
tf_net <- unique(tf_net[target %in% rownames(logcpm), .(source, target, mor = weight)])
tf_sizes <- tf_net[, .N, by = source]
tf_net <- tf_net[source %in% tf_sizes[N >= 10, source]]

message(sprintf("logCPM: %d genes x %d pseudobulks; PROGENy edges: %d; CollecTRI edges: %d",
                nrow(logcpm), ncol(logcpm), nrow(prog_net), nrow(tf_net)))
stopifnot(nrow(prog_net) > 0, nrow(tf_net) > 0, !anyDuplicated(rownames(logcpm)))
stopifnot(!anyNA(prog_net), !anyNA(tf_net))

# Signed footprint-weighted activity. Gene-wise standardization prevents highly
# expressed granule genes from dominating and the L2 denominator makes scores
# comparable across regulons of different sizes.
score_network <- function(mat, network) {
  gene_z <- t(scale(t(mat)))
  gene_z[!is.finite(gene_z)] <- 0
  sources <- sort(unique(network$source))
  w <- sparseMatrix(i = match(network$target, rownames(mat)),
                    j = match(network$source, sources), x = network$mor,
                    dims = c(nrow(mat), length(sources)),
                    dimnames = list(rownames(mat), sources))
  scores <- as.matrix(crossprod(w, gene_z)) /
    pmax(sqrt(Matrix::colSums(w^2)), 1e-8)
  data.table(source = rep(rownames(scores), times = ncol(scores)),
             condition = rep(colnames(scores), each = nrow(scores)),
             score = as.vector(scores))
}
prog_scores <- score_network(logcpm, prog_net)
message("PROGENy activities calculated")
tf_scores <- score_network(logcpm, tf_net)
activity <- rbindlist(list(
  prog_scores[, .(sample_state = condition, activity_type = "Pathway", regulator = source, activity = score)],
  tf_scores[, .(sample_state = condition, activity_type = "TF", regulator = source, activity = score)]
))
activity <- merge(activity, pb_meta, by = "sample_state", all.x = TRUE)
fwrite(activity, file.path(source_dir, "analysis2_regulator_pathway_activity_by_sample_state.csv.gz"))

contrast_key <- fread(file.path(base_analysis, "source_data", "analysis_contrasts.csv"))
did_key <- fread(file.path(base_analysis, "source_data", "analysis_interaction_contrasts.csv"))

fit_pair <- function(dt, key) {
  d <- dt[ConditionCode %in% c(key$reference, key$comparison)]
  if (isTRUE(key$paired)) {
    donors <- d[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 2, DonorID]
    d <- d[DonorID %in% donors]
  }
  if (uniqueN(d$DonorID) < 3) return(NULL)
  wide <- dcast(d, activity_type + regulator ~ sample_state, value.var = "activity")
  sample_cols <- setdiff(names(wide), c("activity_type", "regulator"))
  sm <- unique(d[, .(sample_state, DonorID, ConditionCode)])
  sm <- sm[match(sample_cols, sample_state)]
  mat <- as.matrix(wide[, ..sample_cols])
  cond <- relevel(factor(sm$ConditionCode), key$reference)
  design <- if (isTRUE(key$paired)) model.matrix(~ factor(sm$DonorID) + cond) else model.matrix(~ cond)
  fit <- eBayes(lmFit(mat, design))
  tab <- as.data.table(topTable(fit, coef = ncol(design), number = Inf, sort.by = "none"))
  cbind(wide[, .(activity_type, regulator)], tab[, .(effect = logFC, t, p_value = P.Value)],
        data.table(contrast_id = key$contrast_id, contrast = key$contrast,
                   n_donors = uniqueN(sm$DonorID)))
}

fit_did <- function(dt, key) {
  conditions <- strsplit(key$conditions, ",", fixed = TRUE)[[1]]
  d <- dt[ConditionCode %in% conditions]
  donors <- d[, .(n = uniqueN(ConditionCode)), by = DonorID][n == length(conditions), DonorID]
  d <- d[DonorID %in% donors]
  if (uniqueN(d$DonorID) < 3) return(NULL)
  wide <- dcast(d, activity_type + regulator ~ sample_state, value.var = "activity")
  sample_cols <- setdiff(names(wide), c("activity_type", "regulator"))
  sm <- unique(d[, .(sample_state, DonorID, ConditionCode)])
  sm <- sm[match(sample_cols, sample_state)]
  mat <- as.matrix(wide[, ..sample_cols])
  sm[, ConditionCode := factor(ConditionCode, levels = conditions)]
  design <- model.matrix(~ 0 + ConditionCode + factor(DonorID), data = sm)
  colnames(design) <- make.names(colnames(design))
  cm <- makeContrasts(contrasts = key$expression, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(mat, design), cm))
  tab <- as.data.table(topTable(fit, coef = 1, number = Inf, sort.by = "none"))
  cbind(wide[, .(activity_type, regulator)], tab[, .(effect = logFC, t, p_value = P.Value)],
        data.table(contrast_id = key$contrast_id, contrast = key$contrast,
                   n_donors = uniqueN(sm$DonorID)))
}

effect_list <- list()
k <- 0L
for (st in unique(activity$state)) {
  d <- activity[state == st]
  for (i in seq_len(nrow(contrast_key))) {
    z <- fit_pair(d, contrast_key[i])
    if (!is.null(z)) { k <- k + 1L; z[, state := st]; effect_list[[k]] <- z }
  }
  for (i in seq_len(nrow(did_key))) {
    z <- fit_did(d, did_key[i])
    if (!is.null(z)) { k <- k + 1L; z[, state := st]; effect_list[[k]] <- z }
  }
}
effects <- rbindlist(effect_list, fill = TRUE)
effects[, FDR := p.adjust(p_value, method = "BH"), by = .(activity_type, state, contrast_id)]
fwrite(effects, file.path(source_dir, "analysis2_regulator_pathway_differential_activity.csv"))

focus_contrasts <- c("CF_vs_N", "UF_vs_N", "UA5_vs_UF", "UD_vs_UF",
                     "alpha5beta1_fibroblast_interaction")
contrast_labels <- c(CF_vs_N = "Control FB vs N", UF_vs_N = "UC FB vs N",
                     UA5_vs_UF = "alpha5beta1i vs UC", UD_vs_UF = "Dual vs UC",
                     alpha5beta1_fibroblast_interaction = "FB-dependent alpha5beta1i")
state_labels <- c("PADI4 neutrophil" = "PADI4", "OSM neutrophil" = "OSM",
                  "CXCR4 neutrophil" = "CXCR4", "MX1/ISG neutrophil" = "MX1/ISG")

path_plot <- effects[activity_type == "Pathway" & contrast_id %in% focus_contrasts]
path_plot[, contrast_short := factor(contrast_labels[contrast_id], levels = unname(contrast_labels))]
path_plot[, state_short := factor(state_labels[state], levels = unname(state_labels))]
p_a <- ggplot(path_plot, aes(contrast_short, regulator, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_point(data = path_plot[FDR < 0.05], shape = 8, size = 1.3) +
  facet_wrap(~state_short, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-3.5, 3.5), oob = scales::squish) +
  labs(title = "Perturbation-responsive pathway activity",
       subtitle = "PROGENy footprint-weighted donor-level effects; * FDR <0.05",
       x = NULL, y = NULL, fill = "Activity effect") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.5))

tf_focus <- effects[activity_type == "TF" & contrast_id %in% focus_contrasts]
tf_priority <- tf_focus[contrast_id %in% c("UA5_vs_UF", "UD_vs_UF", "alpha5beta1_fibroblast_interaction"),
                        .(priority = max(abs(effect) * (1 + -log10(pmax(FDR, 1e-12))))), by = regulator]
top_tfs <- head(tf_priority[order(-priority), regulator], 18)
tf_plot <- tf_focus[regulator %in% top_tfs]
tf_plot[, contrast_short := factor(contrast_labels[contrast_id], levels = unname(contrast_labels))]
tf_plot[, state_short := factor(state_labels[state], levels = unname(state_labels))]
tf_plot[, regulator := factor(regulator, levels = rev(top_tfs))]
p_b <- ggplot(tf_plot, aes(contrast_short, regulator, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_point(data = tf_plot[FDR < 0.05], shape = 8, size = 1.2) +
  facet_wrap(~state_short, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-4, 4), oob = scales::squish) +
  labs(title = "Transcription-factor activity identifies blockade-responsive regulators",
       subtitle = "CollecTRI footprint-weighted donor-level effects; * FDR <0.05",
       x = NULL, y = NULL, fill = "Activity effect") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.5))

interaction_plot <- effects[contrast_id == "alpha5beta1_fibroblast_interaction" &
                              ((activity_type == "Pathway") | regulator %in% top_tfs)]
interaction_plot[, label := paste(state_labels[state], regulator, sep = " | ")]
interaction_plot <- interaction_plot[order(FDR, -abs(effect))][seq_len(min(.N, 24))]
interaction_plot[, label := factor(label, levels = rev(label))]
p_c <- ggplot(interaction_plot, aes(effect, label, color = activity_type, size = -log10(pmax(FDR, 1e-12)))) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = effect, yend = label), linewidth = 0.45, color = "#BDBDBD") +
  geom_point() +
  scale_color_manual(values = c(Pathway = "#0072B2", TF = "#D55E00")) +
  labs(title = "Fibroblast-dependent alpha5beta1 regulatory effect",
       subtitle = "Prioritized effects across canonical neutrophil states",
       x = "Activity effect", y = NULL, color = NULL, size = "-log10 FDR") + theme_pub +
  theme(axis.text.y = element_text(size = 6.5))

fig <- p_a / p_b / p_c + plot_layout(heights = c(1.05, 1.30, 0.90)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(fig, file.path(figure_dir, "Additional_Analysis_2_regulator_and_pathway_activity"),
            width = 16, height = 20)

saveRDS(list(pseudobulk_logCPM = logcpm, pseudobulk_metadata = pb_meta,
             pathway_network = prog_net, TF_network = tf_net,
             activity = activity, effects = effects),
        file.path(object_dir, "additional_analysis2_regulatory_activity.rds"), compress = TRUE)
message("Additional analysis 2 complete")
