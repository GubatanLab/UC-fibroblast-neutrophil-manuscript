suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260821)
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
        panel.grid = element_blank(), axis.text = element_text(size = 7.2),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7), strip.text = element_text(face = "bold", size = 8),
        plot.margin = margin(6, 6, 6, 6))

# Assemble prespecified donor-level outcomes.
comp <- fread(file.path(base_analysis, "source_data", "analysis2_canonical_state_composition_by_sample.csv"))
comp_long <- comp[, .(DonorID, ConditionCode, state, feature = state,
                      family = "Canonical-state CLR abundance", value = clr)]

activity <- fread(file.path(source_dir, "analysis2_regulator_pathway_activity_by_sample_state.csv.gz"))
activity_long <- activity[, .(DonorID, ConditionCode, state, feature = regulator,
                              family = paste(activity_type, "activity"), value = activity)]

modules <- fread(file.path(base_analysis, "source_data", "analysis4_module_eigengenes_by_sample_state.csv"))
module_map <- unique(fread(file.path(base_analysis, "source_data", "analysis4_module_differential_activity.csv"))[
  , .(state, module, module_display)])
modules <- merge(modules, module_map, by = c("state", "module"), all.x = TRUE)
modules[is.na(module_display), module_display := module]
module_long <- modules[, .(DonorID, ConditionCode, state, feature = module_display,
                           family = "Metacell module", value = module_score)]

trajectory <- fread(file.path(base_analysis, "source_data", "analysis6_trajectory_summary_by_sample.csv"))
trajectory_long <- trajectory[, .(DonorID, ConditionCode, state = "All neutrophils",
                                  feature = "Median pseudotime", family = "Trajectory",
                                  value = median_pseudotime)]

outcomes <- rbindlist(list(comp_long, activity_long, module_long, trajectory_long), fill = TRUE)
outcomes <- outcomes[is.finite(value)]
fwrite(outcomes, file.path(source_dir, "analysis4_factorial_model_input_outcomes.csv.gz"))

contrast_defs <- list(
  list(id = "alpha5beta1_fibroblast_interaction",
       label = "FB-dependent alpha5beta1i",
       conditions = c("N", "NA5", "UF", "UA5"),
       weights = c(N = 1, NA5 = -1, UF = -1, UA5 = 1)),
  list(id = "NAMPT_fibroblast_interaction",
       label = "FB-dependent NAMPTi",
       conditions = c("N", "NNAMPT", "UF", "UN"),
       weights = c(N = 1, NNAMPT = -1, UF = -1, UN = 1)),
  list(id = "alpha5beta1_NAMPT_synergy_UC",
       label = "alpha5beta1i x NAMPTi synergy",
       conditions = c("UF", "UN", "UA5", "UD"),
       weights = c(UF = 1, UN = -1, UA5 = -1, UD = 1))
)

fit_feature <- function(d, def) {
  d <- d[ConditionCode %in% def$conditions]
  donors <- d[, .(n = uniqueN(ConditionCode)), by = DonorID][n == length(def$conditions), DonorID]
  d <- d[DonorID %in% donors]
  if (uniqueN(d$DonorID) < 3 || nrow(d) < 8) return(NULL)
  d[, ConditionCode := factor(ConditionCode, levels = def$conditions)]
  design <- model.matrix(~ 0 + ConditionCode + factor(DonorID), data = d)
  if (qr(design)$rank < ncol(design)) return(NULL)
  contrast <- rep(0, ncol(design)); names(contrast) <- colnames(design)
  for (cc in names(def$weights)) contrast[paste0("ConditionCode", cc)] <- def$weights[[cc]]
  fit <- eBayes(contrasts.fit(lmFit(matrix(d$value, nrow = 1), design), contrast))
  se <- fit$stdev.unscaled[1, 1] * sqrt(fit$s2.post[1])
  data.table(effect = fit$coefficients[1, 1], SE = se,
             CI_low = fit$coefficients[1, 1] - qt(0.975, fit$df.total[1]) * se,
             CI_high = fit$coefficients[1, 1] + qt(0.975, fit$df.total[1]) * se,
             t = fit$t[1, 1], p_value = fit$p.value[1, 1],
             n_complete_donors = uniqueN(d$DonorID),
             contrast_id = def$id, contrast = def$label)
}

groups <- split(outcomes, interaction(outcomes$family, outcomes$state, outcomes$feature, drop = TRUE))
results <- rbindlist(lapply(groups, function(d) {
  id <- unique(d[, .(family, state, feature)])
  rbindlist(lapply(contrast_defs, function(def) {
    z <- fit_feature(copy(d), def)
    if (is.null(z)) return(NULL)
    cbind(id, z)
  }), fill = TRUE)
}), fill = TRUE)
results[, FDR := p.adjust(p_value, method = "BH"), by = .(family, contrast_id)]
fwrite(results, file.path(source_dir, "analysis4_unified_factorial_blockade_effects.csv"))

contrast_order <- vapply(contrast_defs, `[[`, character(1), "label")
state_labels <- c("PADI4 neutrophil" = "PADI4", "OSM neutrophil" = "OSM",
                  "CXCR4 neutrophil" = "CXCR4", "MX1/ISG neutrophil" = "MX1/ISG",
                  "All neutrophils" = "All")

comp_plot <- results[family == "Canonical-state CLR abundance"]
comp_plot[, contrast := factor(contrast, levels = contrast_order)]
comp_plot[, state_short := factor(state_labels[state], levels = unname(state_labels[1:4]))]
p_a <- ggplot(comp_plot, aes(contrast, state_short, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%+.2f%s", effect, ifelse(FDR < 0.05, "*", ""))), size = 2.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Factorial effects on canonical-state composition",
       subtitle = "Difference-in-differences and within-UC synergy; * FDR <0.05",
       x = NULL, y = NULL, fill = "CLR effect") + theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

path_plot <- results[family == "Pathway activity"]
path_priority <- path_plot[, .(priority = max(abs(effect) * (1 + -log10(pmax(FDR, 1e-12))))), by = feature]
top_paths <- head(path_priority[order(-priority), feature], 12)
path_plot <- path_plot[feature %in% top_paths]
path_plot[, contrast := factor(contrast, levels = contrast_order)]
path_plot[, state_short := factor(state_labels[state], levels = unname(state_labels[1:4]))]
path_plot[, feature := factor(feature, levels = rev(top_paths))]
p_b <- ggplot(path_plot, aes(contrast, feature, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_point(data = path_plot[FDR < 0.05], shape = 8, size = 1.2) +
  facet_wrap(~state_short, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-3.5, 3.5), oob = scales::squish) +
  labs(title = "Factorial pathway responses",
       subtitle = "PROGENy footprint activities; * FDR <0.05",
       x = NULL, y = NULL, fill = "Effect") + theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6.4))

tf_plot <- results[family == "TF activity"]
tf_priority <- tf_plot[, .(priority = max(abs(effect) * (1 + -log10(pmax(FDR, 1e-12))))), by = feature]
top_tfs <- head(tf_priority[order(-priority), feature], 16)
tf_plot <- tf_plot[feature %in% top_tfs]
tf_plot[, contrast := factor(contrast, levels = contrast_order)]
tf_plot[, state_short := factor(state_labels[state], levels = unname(state_labels[1:4]))]
tf_plot[, feature := factor(feature, levels = rev(top_tfs))]
p_c <- ggplot(tf_plot, aes(contrast, feature, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_point(data = tf_plot[FDR < 0.05], shape = 8, size = 1.1) +
  facet_wrap(~state_short, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-4, 4), oob = scales::squish) +
  labs(title = "Factorial transcription-factor responses",
       subtitle = "Signed CollecTRI footprint activities; * FDR <0.05",
       x = NULL, y = NULL, fill = "Effect") + theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6.2))

module_plot <- results[family %in% c("Metacell module", "Trajectory")]
module_priority <- module_plot[, .(priority = max(abs(effect) * (1 + -log10(pmax(FDR, 1e-12))))),
                               by = .(family, feature)]
top_modules <- head(module_priority[family == "Metacell module"][order(-priority), feature], 16)
module_plot <- module_plot[family == "Trajectory" | feature %in% top_modules]
module_plot[, row_label := ifelse(family == "Trajectory", "All | Median pseudotime",
                                  paste(state_labels[state], feature, sep = " | "))]
row_order <- rev(unique(module_plot[order(family, FDR, -abs(effect)), row_label]))
module_plot[, row_label := factor(row_label, levels = row_order)]
module_plot[, contrast := factor(contrast, levels = contrast_order)]
p_d <- ggplot(module_plot, aes(contrast, row_label, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(data = module_plot[FDR < 0.05], shape = 8, size = 1.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Factorial module and trajectory effects",
       subtitle = "Top unbiased modules plus donor median pseudotime; * FDR <0.05",
       x = NULL, y = NULL, fill = "Effect") + theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6.2))

fig <- p_a / p_b / p_c / p_d + plot_layout(heights = c(0.70, 1.05, 1.20, 0.95)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(fig, file.path(figure_dir, "Additional_Analysis_4_factorial_blockade_models"),
            width = 16, height = 23)

saveRDS(list(outcomes = outcomes, contrast_definitions = contrast_defs, results = results),
        file.path(object_dir, "additional_analysis4_factorial_models.rds"), compress = TRUE)
message("Additional analysis 4 complete")
