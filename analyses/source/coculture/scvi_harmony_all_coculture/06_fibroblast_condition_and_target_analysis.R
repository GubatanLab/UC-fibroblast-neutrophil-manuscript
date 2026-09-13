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

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture/Fibroblast_Subclustering"
object_file <- file.path(analysis_dir, "objects", "fibroblasts_subclustered_annotated.rds")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

fib <- readRDS(object_file)
DefaultAssay(fib) <- "DECONTX"

condition_order <- c("CF", "UF", "UN", "UA5", "UD")
condition_labels <- c(
  CF = "Control fibroblast\nco-culture",
  UF = "UC fibroblast\nco-culture",
  UN = "UC +\nNAMPTi",
  UA5 = "UC +\nalpha5beta1i",
  UD = "UC +\ndual blockade"
)
condition_colors <- c(CF = "#009E73", UF = "#D55E00", UN = "#E69F00", UA5 = "#7B3294", UD = "#CC79A7")

annotation_order <- levels(fib$fibroblast_annotation)
annotation_short <- c(
  "COL1A1+ inflammatory ECM fibroblast" = "COL1A1+ inflammatory ECM",
  "CSF3R+ immune-like fibroblast" = "CSF3R+ immune-like",
  "CCN1/CCN2+ activated fibroblast" = "CCN1/CCN2+ activated",
  "S100A8/A9+ immune-like fibroblast" = "S100A8/A9+ immune-like"
)
annotation_colors <- c("#0072B2", "#999999", "#D55E00", "#CC79A7")
names(annotation_colors) <- annotation_order

target_genes <- c("FAP", "ITGA5", "ITGB1")
stopifnot(all(target_genes %in% rownames(fib)))
expr <- as.matrix(GetAssayData(fib, assay = "DECONTX", layer = "data")[target_genes, , drop = FALSE])
detection_threshold <- 1e-6
expr[expr <= detection_threshold] <- 0
alpha_double_positive <- expr["ITGA5", ] > detection_threshold & expr["ITGB1", ] > detection_threshold
alpha_score <- ifelse(alpha_double_positive, sqrt(expr["ITGA5", ] * expr["ITGB1", ]), 0)
fib$alpha5beta1_transcript_score <- alpha_score[colnames(fib)]
fib$alpha5beta1_double_positive <- alpha_double_positive[colnames(fib)]

md <- fib@meta.data
md$cell <- rownames(md)
md$ConditionCode <- factor(md$ConditionCode, levels = condition_order)
md$condition_label <- factor(condition_labels[as.character(md$ConditionCode)], levels = unname(condition_labels))
md$annotation_short <- factor(annotation_short[as.character(md$fibroblast_annotation)], levels = unname(annotation_short))
for (g in target_genes) md[[g]] <- expr[g, md$cell]
md$alpha5beta1_score <- alpha_score[md$cell]
md$alpha5beta1_double_positive <- as.numeric(alpha_double_positive[md$cell])

# Sample-level subtype composition, retaining explicit zeros.
sample_info <- md %>% distinct(SampleID, DonorID, ConditionCode, condition_label) %>%
  mutate(total_fibroblasts = as.integer(table(md$SampleID)[SampleID]))
sample_comp <- tidyr::expand_grid(
  SampleID = unique(md$SampleID),
  fibroblast_annotation = annotation_order
) %>%
  left_join(sample_info, by = "SampleID") %>%
  left_join(
    md %>% count(SampleID, fibroblast_annotation, name = "cells"),
    by = c("SampleID", "fibroblast_annotation")
  ) %>%
  mutate(
    cells = replace_na(cells, 0L),
    proportion = cells / total_fibroblasts,
    annotation_short = factor(annotation_short[as.character(fibroblast_annotation)], levels = unname(annotation_short)),
    ConditionCode = factor(ConditionCode, levels = condition_order),
    condition_label = factor(condition_label, levels = unname(condition_labels))
  )

condition_comp <- sample_comp %>%
  group_by(ConditionCode, condition_label, fibroblast_annotation, annotation_short) %>%
  summarise(
    cells = sum(cells), total_fibroblasts = sum(total_fibroblasts),
    pooled_proportion = cells / total_fibroblasts,
    sample_median_proportion = median(proportion),
    samples = n(), .groups = "drop"
  )

contrasts <- tibble::tribble(
  ~reference, ~comparison, ~contrast,
  "CF", "UF", "UF vs CF",
  "UF", "UN", "NAMPTi vs UF",
  "UF", "UA5", "alpha5beta1i vs UF",
  "UF", "UD", "Dual vs UF"
)

paired_test <- function(dat, value_col, reference, comparison) {
  wide <- dat %>%
    filter(as.character(ConditionCode) %in% c(reference, comparison)) %>%
    select(DonorID, ConditionCode, value = all_of(value_col)) %>%
    group_by(DonorID, ConditionCode) %>% summarise(value = mean(value), .groups = "drop") %>%
    pivot_wider(names_from = ConditionCode, values_from = value)
  wide <- wide[complete.cases(wide[, c(reference, comparison)]), , drop = FALSE]
  if (nrow(wide) < 3) return(tibble(n_pairs = nrow(wide), median_paired_change = NA_real_, p_value = NA_real_))
  tst <- suppressWarnings(wilcox.test(wide[[comparison]], wide[[reference]], paired = TRUE, exact = FALSE))
  tibble(
    n_pairs = nrow(wide),
    median_paired_change = median(wide[[comparison]] - wide[[reference]]),
    p_value = unname(tst$p.value)
  )
}

composition_tests <- bind_rows(lapply(annotation_order, function(a) {
  dat <- sample_comp %>% filter(fibroblast_annotation == a)
  bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
    bind_cols(
      tibble(fibroblast_annotation = a), contrasts[i, ],
      paired_test(dat, "proportion", contrasts$reference[i], contrasts$comparison[i])
    )
  }))
})) %>% mutate(FDR = p.adjust(p_value, method = "BH"))

# Sample-level target transcript summaries. Means include zero-expression cells.
sample_expression <- md %>%
  group_by(SampleID, DonorID, ConditionCode, condition_label) %>%
  summarise(
    fibroblast_cells = n(),
    FAP = mean(FAP), ITGA5 = mean(ITGA5), ITGB1 = mean(ITGB1),
    alpha5beta1_score = mean(alpha5beta1_score),
    alpha5beta1_double_positive_pct = 100 * mean(alpha5beta1_double_positive),
    .groups = "drop"
  )

expression_metrics <- c("FAP", "ITGA5", "ITGB1", "alpha5beta1_score", "alpha5beta1_double_positive_pct")
expression_tests <- bind_rows(lapply(expression_metrics, function(metric) {
  bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
    bind_cols(
      tibble(metric = metric), contrasts[i, ],
      paired_test(sample_expression, metric, contrasts$reference[i], contrasts$comparison[i])
    )
  }))
})) %>% mutate(FDR = p.adjust(p_value, method = "BH"))

expression_condition_summary <- sample_expression %>%
  pivot_longer(all_of(expression_metrics), names_to = "metric", values_to = "sample_value") %>%
  group_by(ConditionCode, condition_label, metric) %>%
  summarise(
    samples = n(), median = median(sample_value), mean = mean(sample_value),
    q1 = quantile(sample_value, 0.25), q3 = quantile(sample_value, 0.75),
    .groups = "drop"
  )

subtype_condition_expression <- md %>%
  group_by(ConditionCode, condition_label, fibroblast_annotation, annotation_short) %>%
  summarise(
    cells = n(),
    FAP_mean = mean(FAP), FAP_pct = 100 * mean(FAP > detection_threshold),
    ITGA5_mean = mean(ITGA5), ITGA5_pct = 100 * mean(ITGA5 > detection_threshold),
    ITGB1_mean = mean(ITGB1), ITGB1_pct = 100 * mean(ITGB1 > detection_threshold),
    alpha5beta1_score_mean = mean(alpha5beta1_score),
    alpha5beta1_double_positive_pct = 100 * mean(alpha5beta1_double_positive),
    .groups = "drop"
  )

fwrite(sample_comp, file.path(source_dir, "fibroblast_subtype_composition_by_sample.csv"))
fwrite(condition_comp, file.path(source_dir, "fibroblast_subtype_composition_by_condition.csv"))
fwrite(composition_tests, file.path(source_dir, "fibroblast_subtype_condition_paired_tests.csv"))
fwrite(sample_expression, file.path(source_dir, "fibroblast_FAP_alpha5beta1_expression_by_sample.csv"))
fwrite(expression_condition_summary, file.path(source_dir, "fibroblast_FAP_alpha5beta1_expression_condition_summary.csv"))
fwrite(expression_tests, file.path(source_dir, "fibroblast_FAP_alpha5beta1_condition_paired_tests.csv"))
fwrite(subtype_condition_expression, file.path(source_dir, "fibroblast_FAP_alpha5beta1_by_subtype_and_condition.csv"))

theme_pub <- theme_bw(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "#444444"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_comp_pooled <- ggplot(condition_comp, aes(condition_label, pooled_proportion, fill = fibroblast_annotation)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = annotation_colors, labels = annotation_short, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Pooled fibroblast subtype composition",
    subtitle = "Bars show the fraction of fibroblast calls within each culture condition",
    x = NULL, y = "Fibroblast subtype (%)", fill = "Fibroblast subtype"
  ) + theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))

p_comp_sample <- ggplot(sample_comp, aes(condition_label, proportion, color = ConditionCode)) +
  geom_boxplot(aes(group = condition_label), width = 0.56, outlier.shape = NA, color = "#777777", fill = NA) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.9, alpha = 0.82) +
  facet_wrap(~annotation_short, ncol = 2, scales = "free_y") +
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Sample-level fibroblast subtype proportions",
    subtitle = "Each point is one donor-condition sample; inference is exploratory because only 238 fibroblasts were recovered",
    x = NULL, y = "Within-sample fibroblast proportion"
  ) + theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5))

p_comp <- p_comp_pooled / p_comp_sample + plot_layout(heights = c(0.9, 1.35)) +
  plot_annotation(tag_levels = "A")

target_long <- md %>%
  select(cell, ConditionCode, condition_label, fibroblast_annotation, annotation_short, all_of(target_genes)) %>%
  pivot_longer(all_of(target_genes), names_to = "gene", values_to = "expression")
target_dot_condition <- target_long %>%
  group_by(ConditionCode, condition_label, gene) %>%
  summarise(avg_expression = mean(expression), pct_expressed = 100 * mean(expression > detection_threshold), .groups = "drop") %>%
  group_by(gene) %>% mutate(scaled_expression = as.numeric(scale(avg_expression))) %>% ungroup()
target_dot_subtype <- target_long %>%
  group_by(fibroblast_annotation, annotation_short, gene) %>%
  summarise(avg_expression = mean(expression), pct_expressed = 100 * mean(expression > detection_threshold), .groups = "drop") %>%
  group_by(gene) %>% mutate(scaled_expression = as.numeric(scale(avg_expression))) %>% ungroup()

p_target_condition <- ggplot(target_dot_condition, aes(gene, condition_label, size = pct_expressed, color = scaled_expression)) +
  geom_point() +
  scale_size_area(max_size = 9, limits = c(0, 100)) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Target transcripts by culture condition",
    subtitle = "FAP and the alpha5beta1 integrin subunits ITGA5 and ITGB1",
    x = NULL, y = NULL, color = "Scaled mean", size = "Cells detected (%)"
  ) + theme_pub + theme(axis.text.x = element_text(face = "italic", size = 11))

p_target_subtype <- ggplot(target_dot_subtype, aes(gene, annotation_short, size = pct_expressed, color = scaled_expression)) +
  geom_point() +
  scale_size_area(max_size = 9, limits = c(0, 100)) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Target transcripts by fibroblast subtype",
    subtitle = "Expression is based on DecontX-corrected log-normalized RNA",
    x = NULL, y = NULL, color = "Scaled mean", size = "Cells detected (%)"
  ) + theme_pub + theme(axis.text.x = element_text(face = "italic", size = 11))

sample_expression_long <- sample_expression %>%
  select(SampleID, DonorID, ConditionCode, condition_label, FAP, ITGA5, ITGB1, alpha5beta1_score) %>%
  pivot_longer(c(FAP, ITGA5, ITGB1, alpha5beta1_score), names_to = "metric", values_to = "mean_expression") %>%
  mutate(metric = factor(metric,
    levels = c("FAP", "ITGA5", "ITGB1", "alpha5beta1_score"),
    labels = c("FAP", "ITGA5", "ITGB1", "alpha5beta1 joint score")
  ))

p_target_samples <- ggplot(sample_expression_long, aes(condition_label, mean_expression, color = ConditionCode)) +
  geom_boxplot(aes(group = condition_label), width = 0.56, outlier.shape = NA, color = "#777777", fill = NA) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.9, alpha = 0.82) +
  facet_wrap(~metric, ncol = 2, scales = "free_y") +
  scale_color_manual(values = condition_colors, guide = "none") +
  labs(
    title = "Sample-level FAP and alpha5beta1 transcript expression",
    subtitle = "Joint score is the geometric mean of ITGA5 and ITGB1 log-expression and is zero unless both are detected",
    x = NULL, y = "Mean log-normalized expression"
  ) + theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5))

coexpr_heat <- subtype_condition_expression %>%
  mutate(label = ifelse(cells < 3, paste0(round(alpha5beta1_double_positive_pct), "%\n(n=", cells, ")"), paste0(round(alpha5beta1_double_positive_pct), "%")))
p_coexpr <- ggplot(coexpr_heat, aes(condition_label, annotation_short, fill = alpha5beta1_double_positive_pct)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 3.1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 100), na.value = "#EEEEEE") +
  labs(
    title = "ITGA5/ITGB1 double-positive fibroblasts",
    subtitle = "Percent co-expressing both alpha5beta1 subunit transcripts; n is shown where fewer than three cells were available",
    x = NULL, y = NULL, fill = "Double-positive (%)"
  ) + theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_targets <- (p_target_condition | p_target_subtype) / p_target_samples / p_coexpr +
  plot_layout(heights = c(0.8, 1.2, 0.9)) + plot_annotation(tag_levels = "A")

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_plot(p_comp, "Fibroblast_4_subtype_composition_by_condition", 12.5, 12.5)
save_plot(p_targets, "Fibroblast_5_FAP_alpha5beta1_expression", 13.5, 16.0)

saveRDS(fib, object_file, compress = FALSE)
message("Condition composition and FAP/alpha5beta1 analysis complete")
