#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260816)

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture/Neutrophil_Subclustering"
input_file <- file.path(analysis_dir, "objects", "neutrophils_scVI_Harmony_subclustered_annotated_clean.rds")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

neut <- readRDS(input_file)
md <- neut@meta.data %>%
  transmute(
    DonorID = as.character(DonorID),
    ConditionCode = as.character(ConditionCode),
    neutrophil_state_annotation = as.character(neutrophil_state_annotation)
  )

condition_order <- c("PRE", "N", "NNAMPT", "NA5", "CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(
  PRE = "Pre-culture", N = "Neutrophils alone", NNAMPT = "Neutrophils + NAMPTi",
  NA5 = "Neutrophils + alpha5beta1i", CF = "Control fibroblast co-culture",
  UF = "UC fibroblast co-culture", UN = "UC + NAMPTi",
  UA5 = "UC + alpha5beta1i", UD = "UC + dual blockade"
)
condition_colors <- c(
  PRE = "#111827", N = "#4D4D4D", NNAMPT = "#56B4E9", NA5 = "#0072B2",
  CF = "#009E73", UF = "#D55E00", UN = "#E69F00", UA5 = "#7B3294", UD = "#CC79A7"
)
state_order <- c(
  "Low RNA neutrophils", "PADI4 neutrophil", "OSM neutrophil",
  "CXCR4 neutrophil", "MX1/ISG neutrophil"
)
state_colors <- c(
  "Low RNA neutrophils" = "#999999", "PADI4 neutrophil" = "#CC79A7",
  "OSM neutrophil" = "#D55E00", "CXCR4 neutrophil" = "#0072B2",
  "MX1/ISG neutrophil" = "#009E73"
)

sample_info <- md %>% distinct(DonorID, ConditionCode)
sample_totals <- md %>% count(DonorID, ConditionCode, name = "total_neutrophils")
observed_counts <- md %>%
  count(DonorID, ConditionCode, neutrophil_state_annotation, name = "cells")

state_abundance <- sample_info %>%
  crossing(neutrophil_state_annotation = state_order) %>%
  left_join(observed_counts, by = c("DonorID", "ConditionCode", "neutrophil_state_annotation")) %>%
  left_join(sample_totals, by = c("DonorID", "ConditionCode")) %>%
  mutate(
    cells = replace_na(cells, 0L),
    proportion = cells / total_neutrophils,
    condition_label = unname(condition_labels[ConditionCode]),
    ConditionCode = factor(ConditionCode, levels = condition_order),
    condition_label = factor(condition_label, levels = unname(condition_labels)),
    neutrophil_state_annotation = factor(neutrophil_state_annotation, levels = state_order)
  ) %>%
  arrange(ConditionCode, DonorID, neutrophil_state_annotation)
fwrite(state_abundance, file.path(source_dir, "neutrophil_state_differential_abundance_donor_condition.csv"))

contrast_key <- tibble::tribble(
  ~reference, ~comparison, ~contrast, ~paired,
  "PRE", "N", "Neutrophils alone vs pre-culture", FALSE,
  "N", "CF", "Control fibroblast vs neutrophils alone", TRUE,
  "N", "UF", "UC fibroblast vs neutrophils alone", TRUE,
  "CF", "UF", "UC vs control fibroblast", TRUE,
  "N", "NNAMPT", "NAMPTi vs neutrophils alone", TRUE,
  "N", "NA5", "alpha5beta1i vs neutrophils alone", TRUE,
  "UF", "UN", "NAMPTi vs UC", TRUE,
  "UF", "UA5", "alpha5beta1i vs UC", TRUE,
  "UF", "UD", "Dual blockade vs UC", TRUE
)
fwrite(contrast_key, file.path(source_dir, "neutrophil_state_differential_abundance_contrasts.csv"))

safe_wilcox <- function(x, y, paired) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  if (paired && all(abs(x - y) < .Machine$double.eps^0.5)) return(1)
  if (!paired && length(unique(c(x, y))) == 1) return(1)
  tryCatch(
    suppressWarnings(wilcox.test(x, y, paired = paired, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

test_one <- function(state, reference, comparison, contrast, paired) {
  dat <- state_abundance %>%
    filter(
      as.character(neutrophil_state_annotation) == state,
      as.character(ConditionCode) %in% c(reference, comparison)
    )
  if (paired) {
    wide <- dat %>%
      select(DonorID, ConditionCode, proportion) %>%
      pivot_wider(names_from = ConditionCode, values_from = proportion) %>%
      filter(!is.na(.data[[reference]]), !is.na(.data[[comparison]]))
    ref_values <- wide[[reference]]
    comp_values <- wide[[comparison]]
    n_reference <- nrow(wide)
    n_comparison <- nrow(wide)
    n_pairs <- nrow(wide)
    effect <- median(comp_values - ref_values) * 100
    mean_effect <- mean(comp_values - ref_values) * 100
    p_value <- safe_wilcox(comp_values, ref_values, paired = TRUE)
  } else {
    ref_values <- dat$proportion[as.character(dat$ConditionCode) == reference]
    comp_values <- dat$proportion[as.character(dat$ConditionCode) == comparison]
    n_reference <- length(ref_values)
    n_comparison <- length(comp_values)
    n_pairs <- NA_integer_
    effect <- (median(comp_values) - median(ref_values)) * 100
    mean_effect <- (mean(comp_values) - mean(ref_values)) * 100
    p_value <- safe_wilcox(comp_values, ref_values, paired = FALSE)
  }
  tibble(
    neutrophil_state_annotation = state,
    reference = reference,
    comparison = comparison,
    contrast = contrast,
    paired = paired,
    n_reference = n_reference,
    n_comparison = n_comparison,
    n_pairs = n_pairs,
    reference_median_pct = median(ref_values) * 100,
    comparison_median_pct = median(comp_values) * 100,
    median_change_percentage_points = effect,
    mean_change_percentage_points = mean_effect,
    p_value = p_value
  )
}

da_results <- bind_rows(lapply(seq_len(nrow(contrast_key)), function(i) {
  bind_rows(lapply(state_order, function(state) {
    test_one(
      state,
      contrast_key$reference[i], contrast_key$comparison[i],
      contrast_key$contrast[i], contrast_key$paired[i]
    )
  }))
})) %>%
  group_by(contrast) %>%
  mutate(FDR_within_contrast = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    significance = case_when(
      FDR_global < 0.001 ~ "***",
      FDR_global < 0.01 ~ "**",
      FDR_global < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
fwrite(da_results, file.path(source_dir, "neutrophil_state_differential_abundance_tests.csv"))

pooled_composition <- md %>%
  count(ConditionCode, neutrophil_state_annotation, name = "cells") %>%
  group_by(ConditionCode) %>%
  mutate(total_neutrophils = sum(cells), proportion = cells / total_neutrophils) %>%
  ungroup() %>%
  mutate(
    condition_label = factor(unname(condition_labels[ConditionCode]), levels = unname(condition_labels)),
    neutrophil_state_annotation = factor(neutrophil_state_annotation, levels = state_order)
  )
fwrite(pooled_composition, file.path(source_dir, "neutrophil_state_differential_abundance_pooled_composition.csv"))

p_composition <- ggplot(
  pooled_composition,
  aes(condition_label, proportion, fill = neutrophil_state_annotation)
) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Neutrophil state composition across culture conditions",
    subtitle = "Pre-culture is included; bars show pooled cell fractions after all quality exclusions",
    x = NULL, y = "Neutrophils (%)", fill = "Neutrophil state"
  ) +
  theme_bw(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text.x = element_text(angle = 38, hjust = 1, size = 8.5),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2))

p_donor <- ggplot(
  state_abundance,
  aes(condition_label, proportion, color = ConditionCode)
) +
  geom_boxplot(aes(group = condition_label), width = 0.62, outlier.shape = NA, color = "#777777", fill = NA) +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 1.9, alpha = 0.85) +
  facet_wrap(~neutrophil_state_annotation, ncol = 3, scales = "free_y") +
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Donor-level state abundance",
    subtitle = "Each point is one donor-condition; facet-specific y-scales expose both common and rare states",
    x = NULL, y = "State fraction within donor-condition"
  ) +
  theme_bw(base_family = "Arial", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 42, hjust = 1, size = 7.5),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()
  )

contrast_order <- contrast_key$contrast
heatmap_df <- da_results %>%
  mutate(
    contrast = factor(contrast, levels = contrast_order),
    neutrophil_state_annotation = factor(neutrophil_state_annotation, levels = rev(state_order)),
    effect_label = paste0(sprintf("%+.1f", median_change_percentage_points), significance)
  )
effect_limit <- max(abs(heatmap_df$median_change_percentage_points), na.rm = TRUE)
p_effect <- ggplot(
  heatmap_df,
  aes(contrast, neutrophil_state_annotation, fill = median_change_percentage_points)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = effect_label), size = 3.1, fontface = ifelse(heatmap_df$significance == "", "plain", "bold")) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(-effect_limit, effect_limit)
  ) +
  labs(
    title = "Differential abundance of neutrophil states",
    subtitle = "Values are median percentage-point changes (comparison minus reference); global Wilcoxon FDR symbols: * <0.05, ** <0.01, *** <0.001",
    x = NULL, y = NULL, fill = "Median change\n(percentage points)"
  ) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 42, hjust = 1, size = 8.2),
    axis.text.y = element_text(size = 9.5),
    panel.grid = element_blank()
  )

p_all <- p_composition / p_donor / p_effect +
  plot_layout(heights = c(0.95, 1.45, 1.05)) +
  plot_annotation(tag_levels = "A")

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_plot(p_all, "Neutrophil_4_state_differential_abundance_by_condition", 15.5, 18)
save_plot(p_effect, "Neutrophil_4B_state_differential_abundance_heatmap", 15, 6.5)

message("Neutrophil state differential-abundance analysis complete")
