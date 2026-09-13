#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
})

set.seed(20260816)

analysis_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture/Neutrophil_Subclustering"
input_file <- file.path(analysis_dir, "objects", "neutrophils_scVI_Harmony_subclustered_preannotation.rds")
object_dir <- file.path(analysis_dir, "objects")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

neut <- readRDS(input_file)
DefaultAssay(neut) <- "RNA"

cluster_key <- data.frame(
  neutrophil_subcluster = as.character(1:31),
  neutrophil_cluster_annotation = "Singleton/outlier",
  annotation_confidence = "Low",
  defining_features = "Very small graph community",
  stringsAsFactors = FALSE
)
set_cluster <- function(ids, label, confidence, features) {
  idx <- cluster_key$neutrophil_subcluster %in% as.character(ids)
  cluster_key$neutrophil_cluster_annotation[idx] <<- label
  cluster_key$annotation_confidence[idx] <<- confidence
  cluster_key$defining_features[idx] <<- features
}
set_cluster(1, "Low-signal neutrophil", "Moderate", "Low RNA signal; majority prior low-signal/LowRNA calls")
set_cluster(2, "Activated OSM/CXCR4 neutrophil", "High", "FCGR3B, FPR1, FOS, OSM, CXCR4; activated mature program")
set_cluster(3, "IL1R2+ inflammatory neutrophil", "Moderate", "IL1R2 with inflammatory/stress program")
set_cluster(4, "DHFR+ proliferative/stress neutrophil", "Low", "DHFR, DNAJA1 and stress-associated program")
set_cluster(5, "PADI4/translation-high neutrophil", "Moderate", "PADI4 enrichment with EEF1A1, PTMA, RACK1")
set_cluster(6, "HLA-II antigen-presenting/monocyte-like", "Low", "HLA-DRA, HLA-DRB1, CD74, LYZ, VCAN")
set_cluster(7, "CLC+ granulocyte-like", "Low", "CLC, PDE4D; possible eosinophil-like contamination")
set_cluster(8, "SPARC/CLU vascular-like outlier", "Low", "SPARC, CLU, ITGB3, GNG11")
set_cluster(9, "LTF/BPI immature neutrophil", "High", "LTF, BPI, CD24, S100A8")
set_cluster(10, "RORA+ atypical neutrophil", "Low", "RORA, SETBP1, ZBTB20")
set_cluster(11, "B-cell-like doublet", "Low", "CD79A, MS4A1, HLA-DQA1, HLA-DPA1")
set_cluster(28, "CCL3/CCL4 inflammatory neutrophil", "Moderate", "CCL3, CCL3L3, CCL4L2")
write.csv(cluster_key, file.path(source_dir, "neutrophil_subcluster_annotation_key.csv"), row.names = FALSE)

neut$neutrophil_cluster_annotation <- cluster_key$neutrophil_cluster_annotation[
  match(as.character(neut$neutrophil_subcluster), cluster_key$neutrophil_subcluster)
]
annotation_order <- c(
  "Low-signal neutrophil", "Activated OSM/CXCR4 neutrophil", "IL1R2+ inflammatory neutrophil",
  "DHFR+ proliferative/stress neutrophil", "PADI4/translation-high neutrophil",
  "LTF/BPI immature neutrophil", "CCL3/CCL4 inflammatory neutrophil",
  "HLA-II antigen-presenting/monocyte-like", "CLC+ granulocyte-like",
  "RORA+ atypical neutrophil", "SPARC/CLU vascular-like outlier", "B-cell-like doublet",
  "Singleton/outlier"
)
neut$neutrophil_cluster_annotation <- factor(neut$neutrophil_cluster_annotation, levels = annotation_order)
neut$neutrophil_cluster_annotation_confidence <- cluster_key$annotation_confidence[
  match(as.character(neut$neutrophil_subcluster), cluster_key$neutrophil_subcluster)
]

state_map <- c(
  "OSM" = "OSM neutrophil",
  "CXCR4" = "CXCR4 neutrophil",
  "MX1" = "MX1/ISG neutrophil",
  "MX1/ISG-like" = "MX1/ISG neutrophil",
  "PADI4" = "PADI4 neutrophil",
  "LowRNA" = "Low RNA neutrophils",
  "Low-signal/unresolved" = "Low-signal/unresolved neutrophil"
)
state_order <- c(
  "Low RNA neutrophils", "PADI4 neutrophil", "OSM neutrophil",
  "CXCR4 neutrophil", "MX1/ISG neutrophil", "Low-signal/unresolved neutrophil"
)
neut$neutrophil_state_annotation <- factor(unname(state_map[as.character(neut$clean_state)]), levels = state_order)

# Exclude graph communities whose marker programs support doublets, non-neutrophil
# contaminants, or technical/singleton outliers. The preannotation object remains
# unchanged and provides the complete unfiltered record.
remove_annotations <- c(
  "HLA-II antigen-presenting/monocyte-like",
  "CLC+ granulocyte-like",
  "SPARC/CLU vascular-like outlier",
  "B-cell-like doublet",
  "Singleton/outlier"
)
remove_states <- "Low-signal/unresolved neutrophil"
cleanup_audit <- neut@meta.data %>%
  mutate(
    neutrophil_cluster_annotation = as.character(neutrophil_cluster_annotation),
    neutrophil_state_annotation = as.character(neutrophil_state_annotation),
    cleanup_reason = case_when(
      neutrophil_cluster_annotation %in% remove_annotations ~ "Doublet, contaminant-like, or singleton/outlier community",
      neutrophil_state_annotation %in% remove_states ~ "Low-signal/unresolved neutrophil state",
      TRUE ~ "Retained neutrophil"
    ),
    cleanup_action = if_else(cleanup_reason == "Retained neutrophil", "Retained", "Removed")
  ) %>%
  count(
    neutrophil_subcluster, neutrophil_cluster_annotation,
    neutrophil_state_annotation, cleanup_action, cleanup_reason,
    name = "cells"
  ) %>%
  arrange(desc(cleanup_action), cleanup_reason, as.numeric(as.character(neutrophil_subcluster)))
fwrite(cleanup_audit, file.path(source_dir, "neutrophil_cleanup_audit.csv"))

keep_cells <- rownames(neut@meta.data)[
  !as.character(neut$neutrophil_cluster_annotation) %in% remove_annotations &
    !as.character(neut$neutrophil_state_annotation) %in% remove_states
]
neut <- subset(neut, cells = keep_cells)
annotation_order <- setdiff(annotation_order, remove_annotations)
state_order <- setdiff(state_order, remove_states)
neut$neutrophil_cluster_annotation <- factor(
  as.character(neut$neutrophil_cluster_annotation), levels = annotation_order
)
neut$neutrophil_state_annotation <- factor(
  as.character(neut$neutrophil_state_annotation), levels = state_order
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
neut$ConditionCode <- factor(neut$ConditionCode, levels = condition_order)
neut$condition_label <- factor(unname(condition_labels[as.character(neut$ConditionCode)]), levels = unname(condition_labels))

requested_genes <- c("OSM", "CXCR4", "MX1", "PADI4")
stopifnot(all(requested_genes %in% rownames(neut)))
expr <- LayerData(neut[["RNA"]], layer = "data")[requested_genes, , drop = FALSE]
detection_threshold <- 1e-6
md <- neut@meta.data
md$cell <- rownames(md)
for (g in requested_genes) md[[g]] <- as.numeric(expr[g, md$cell])

# Donor-condition summaries are the inferential units; multiple sample tags for a donor-condition are combined.
sample_expression <- md %>%
  group_by(DonorID, ConditionCode, condition_label) %>%
  summarise(
    neutrophil_cells = n(),
    OSM_mean = mean(OSM), OSM_pct = 100 * mean(OSM > detection_threshold),
    CXCR4_mean = mean(CXCR4), CXCR4_pct = 100 * mean(CXCR4 > detection_threshold),
    MX1_mean = mean(MX1), MX1_pct = 100 * mean(MX1 > detection_threshold),
    PADI4_mean = mean(PADI4), PADI4_pct = 100 * mean(PADI4 > detection_threshold),
    .groups = "drop"
  )
fwrite(sample_expression, file.path(source_dir, "neutrophil_OSM_CXCR4_MX1_PADI4_by_donor_condition.csv"))

condition_expression <- md %>%
  group_by(ConditionCode, condition_label) %>%
  summarise(
    cells = n(), donors = n_distinct(DonorID),
    OSM_mean = mean(OSM), OSM_pct = 100 * mean(OSM > detection_threshold),
    CXCR4_mean = mean(CXCR4), CXCR4_pct = 100 * mean(CXCR4 > detection_threshold),
    MX1_mean = mean(MX1), MX1_pct = 100 * mean(MX1 > detection_threshold),
    PADI4_mean = mean(PADI4), PADI4_pct = 100 * mean(PADI4 > detection_threshold),
    .groups = "drop"
  )
fwrite(condition_expression, file.path(source_dir, "neutrophil_OSM_CXCR4_MX1_PADI4_condition_summary.csv"))

contrasts <- tibble::tribble(
  ~reference, ~comparison, ~contrast,
  "PRE", "N", "N vs pre-culture",
  "N", "NNAMPT", "NAMPTi vs neutrophils alone",
  "N", "NA5", "alpha5beta1i vs neutrophils alone",
  "N", "CF", "Control fibroblast vs neutrophils alone",
  "CF", "UF", "UC vs control fibroblast",
  "UF", "UN", "NAMPTi vs UC",
  "UF", "UA5", "alpha5beta1i vs UC",
  "UF", "UD", "Dual blockade vs UC"
)

paired_test <- function(dat, metric, reference, comparison) {
  wide <- dat %>%
    filter(as.character(ConditionCode) %in% c(reference, comparison)) %>%
    select(DonorID, ConditionCode, value = all_of(metric)) %>%
    pivot_wider(names_from = ConditionCode, values_from = value)
  wide <- wide[complete.cases(wide[, c(reference, comparison)]), , drop = FALSE]
  if (nrow(wide) < 3) return(tibble(n_pairs = nrow(wide), median_paired_change = NA_real_, p_value = NA_real_))
  tst <- suppressWarnings(wilcox.test(wide[[comparison]], wide[[reference]], paired = TRUE, exact = FALSE))
  tibble(n_pairs = nrow(wide), median_paired_change = median(wide[[comparison]] - wide[[reference]]), p_value = unname(tst$p.value))
}

metrics <- unlist(lapply(requested_genes, function(g) c(paste0(g, "_mean"), paste0(g, "_pct"))))
paired_tests <- bind_rows(lapply(metrics, function(metric) {
  bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
    bind_cols(tibble(metric = metric), contrasts[i, ], paired_test(
      sample_expression, metric, contrasts$reference[i], contrasts$comparison[i]
    ))
  }))
})) %>% group_by(grepl("_pct$", metric)) %>% mutate(FDR = p.adjust(p_value, method = "BH")) %>% ungroup()
fwrite(paired_tests, file.path(source_dir, "neutrophil_OSM_CXCR4_MX1_PADI4_paired_condition_tests.csv"))

state_sample <- md %>% count(DonorID, ConditionCode, condition_label, neutrophil_state_annotation, name = "cells") %>%
  group_by(DonorID, ConditionCode, condition_label) %>%
  mutate(total_neutrophils = sum(cells), proportion = cells / total_neutrophils) %>% ungroup()
state_condition <- md %>% count(ConditionCode, condition_label, neutrophil_state_annotation, name = "cells") %>%
  group_by(ConditionCode, condition_label) %>% mutate(total_neutrophils = sum(cells), proportion = cells / total_neutrophils) %>% ungroup()
fwrite(state_sample, file.path(source_dir, "neutrophil_state_composition_by_donor_condition.csv"))
fwrite(state_condition, file.path(source_dir, "neutrophil_state_composition_by_condition.csv"))

um <- as.data.frame(Embeddings(neut, "umap.neut.scvi.harmony"))
colnames(um) <- c("UMAP_1", "UMAP_2")
um$cell <- rownames(um)
plot_df <- cbind(um, md[um$cell, , drop = FALSE])
plot_df <- plot_df[, !duplicated(colnames(plot_df)), drop = FALSE]

annotation_colors <- c(
  "Low-signal neutrophil" = "#BDBDBD",
  "Activated OSM/CXCR4 neutrophil" = "#D55E00",
  "IL1R2+ inflammatory neutrophil" = "#E69F00",
  "DHFR+ proliferative/stress neutrophil" = "#A6761D",
  "PADI4/translation-high neutrophil" = "#CC79A7",
  "LTF/BPI immature neutrophil" = "#009E73",
  "CCL3/CCL4 inflammatory neutrophil" = "#F781BF",
  "RORA+ atypical neutrophil" = "#8C6BB1"
)
state_colors <- c(
  "Low RNA neutrophils" = "#999999", "PADI4 neutrophil" = "#CC79A7",
  "OSM neutrophil" = "#D55E00", "CXCR4 neutrophil" = "#0072B2",
  "MX1/ISG neutrophil" = "#009E73", "Low-signal/unresolved neutrophil" = "#D9D9D9"
)
theme_umap <- theme_void(base_family = "Arial", base_size = 12) + theme(
  plot.title = element_text(face = "bold", size = 16),
  plot.subtitle = element_text(size = 10.5, color = "#444444"),
  legend.title = element_text(face = "bold"), legend.text = element_text(size = 8.5),
  plot.caption = element_text(size = 8.5, color = "#666666", hjust = 0)
)

p_cluster_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = neutrophil_cluster_annotation)) +
  geom_point_rast(size = 0.16, alpha = 0.52, raster.dpi = 300) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  labs(
    title = "Neutrophil-only scVI-Harmony subclusters",
    subtitle = "Marker-supported communities after removal of doublets, contaminants, and unresolved cells",
    color = "Cluster annotation",
    caption = "Doublet, contaminant-like, singleton/outlier, and low-signal/unresolved cells were excluded before downstream analysis."
  ) + theme_umap + theme(legend.position = "right")

p_state_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = neutrophil_state_annotation)) +
  geom_point_rast(size = 0.16, alpha = 0.52, raster.dpi = 300) +
  scale_color_manual(values = state_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Canonical neutrophil state annotations",
    subtitle = "OSM, CXCR4, MX1/ISG, PADI4, and Low RNA neutrophil programs",
    color = "Neutrophil state"
  ) + theme_umap + theme(legend.position = "right")

p_umap <- p_cluster_umap | p_state_umap

marker_panel <- intersect(c(
  "FCGR3B", "FPR1", "OSM", "CXCR4", "MX1", "ISG15", "IFIT2", "PADI4",
  "IL1R2", "LTF", "BPI", "CCL3", "CCL4", "HLA-DRA", "CD74", "CLC", "SPARC", "MS4A1"
), rownames(neut))
marker_expr <- LayerData(neut[["RNA"]], layer = "data")[marker_panel, , drop = FALSE]
ann <- as.character(neut$neutrophil_cluster_annotation)
marker_dot <- rbindlist(lapply(annotation_order, function(a) {
  idx <- ann == a
  if (!any(idx)) return(NULL)
  data.table(
    annotation = a, gene = marker_panel,
    mean_expression = Matrix::rowMeans(marker_expr[, idx, drop = FALSE]),
    pct_expressed = 100 * Matrix::rowMeans(marker_expr[, idx, drop = FALSE] > detection_threshold),
    cells = sum(idx)
  )
}))
marker_dot[, scaled_expression := as.numeric(scale(mean_expression)), by = gene]
marker_dot$annotation <- factor(marker_dot$annotation, levels = rev(annotation_order))
marker_dot$gene <- factor(marker_dot$gene, levels = marker_panel)
fwrite(marker_dot, file.path(source_dir, "neutrophil_annotation_marker_dotplot_source.csv"))

p_markers <- ggplot(marker_dot, aes(gene, annotation, size = pct_expressed, color = scaled_expression)) +
  geom_point() +
  scale_size_area(max_size = 8, limits = c(0, 100)) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Marker programs supporting neutrophil cluster annotations",
    subtitle = "Dot size is detection frequency; color is scaled mean log-normalized expression",
    x = NULL, y = NULL, size = "Cells detected (%)", color = "Scaled mean"
  ) +
  theme_bw(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text.x = element_text(angle = 50, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 9), panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.25)
  )

condition_dot <- rbindlist(lapply(condition_order, function(cc) {
  idx <- as.character(md$ConditionCode) == cc
  data.table(
    ConditionCode = cc, condition_label = condition_labels[cc], gene = requested_genes,
    mean_expression = Matrix::rowMeans(expr[, idx, drop = FALSE]),
    pct_expressed = 100 * Matrix::rowMeans(expr[, idx, drop = FALSE] > detection_threshold)
  )
}))
condition_dot[, scaled_expression := as.numeric(scale(mean_expression)), by = gene]
condition_dot$condition_label <- factor(condition_dot$condition_label, levels = unname(condition_labels))
condition_dot$gene <- factor(condition_dot$gene, levels = requested_genes)

p_condition_dot <- ggplot(condition_dot, aes(gene, condition_label, size = pct_expressed, color = scaled_expression)) +
  geom_point() +
  scale_size_area(max_size = 10, limits = c(0, 100)) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(
    title = "OSM, CXCR4, MX1, and PADI4 across culture conditions",
    subtitle = paste0(format(ncol(neut), big.mark = ","), " cleaned neutrophils; RNA values are log-normalized raw UMI counts"),
    x = NULL, y = NULL, size = "Cells detected (%)", color = "Scaled mean"
  ) +
  theme_bw(base_family = "Arial", base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 15), axis.text.x = element_text(face = "italic", size = 11), panel.grid.minor = element_blank())

sample_long <- sample_expression %>%
  select(DonorID, ConditionCode, condition_label, ends_with("_mean")) %>%
  pivot_longer(ends_with("_mean"), names_to = "gene", values_to = "mean_expression") %>%
  mutate(gene = factor(sub("_mean$", "", gene), levels = requested_genes))
p_sample <- ggplot(sample_long, aes(condition_label, mean_expression, color = ConditionCode)) +
  geom_boxplot(aes(group = condition_label), width = 0.58, outlier.shape = NA, color = "#777777", fill = NA) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.7, alpha = 0.80) +
  facet_wrap(~gene, ncol = 2, scales = "free_y") +
  scale_color_manual(values = condition_colors, guide = "none") +
  labs(
    title = "Donor-level expression across conditions",
    subtitle = "Each point is one donor-condition mean; paired donor tests are provided in the source tables",
    x = NULL, y = "Mean log-normalized expression"
  ) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14), strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 8), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()
  )

p_state_comp <- ggplot(state_condition, aes(condition_label, proportion, fill = neutrophil_state_annotation)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Neutrophil state composition across conditions",
    subtitle = "Pooled cell fractions shown descriptively; donor-level composition is exported separately",
    x = NULL, y = "Neutrophils (%)", fill = "Neutrophil state"
  ) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14), axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom"
  ) + guides(fill = guide_legend(nrow = 2))

p_expression <- p_condition_dot / p_sample / p_state_comp + plot_layout(heights = c(0.8, 1.25, 0.9)) + plot_annotation(tag_levels = "A")

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_plot(p_umap, "Neutrophil_1_scVI_Harmony_subcluster_annotations", 16, 7.5)
save_plot(p_markers, "Neutrophil_2_annotation_marker_dotplot", 14.5, 8.5)
save_plot(p_expression, "Neutrophil_3_OSM_CXCR4_MX1_PADI4_by_condition", 13.5, 17)

metadata_output <- neut@meta.data
metadata_output$cell <- rownames(metadata_output)
coords <- as.data.frame(Embeddings(neut, "umap.neut.scvi.harmony"))
colnames(coords) <- c("neutrophil_UMAP_1", "neutrophil_UMAP_2")
coords$cell <- rownames(coords)
metadata_output <- merge(metadata_output, coords, by = "cell", sort = FALSE)
fwrite(metadata_output, file.path(source_dir, "neutrophil_cell_metadata_annotations_and_umap.csv.gz"))

saveRDS(neut, file.path(object_dir, "neutrophils_scVI_Harmony_subclustered_annotated.rds"), compress = FALSE)
saveRDS(neut, file.path(object_dir, "neutrophils_scVI_Harmony_subclustered_annotated_clean.rds"), compress = FALSE)
message("Neutrophil annotation and condition marker analysis complete")
