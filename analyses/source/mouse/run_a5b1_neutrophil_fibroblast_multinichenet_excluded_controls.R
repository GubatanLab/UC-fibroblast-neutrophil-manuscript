suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(nichenetr)
  library(multinichenetr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(20260816)

out_dir <- "a5B1_DSS_neutrophil_fibroblast_MultiNicheNet_excluded_controls"
dir.create(file.path(out_dir, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

input_root <- file.path(
  "MouseColitis_Inhibitors_compartment_subclustering",
  "reassigned_fractions_epithelial_audit_cleaned", "updated_umaps", "level3_reannotated"
)
immune_file <- file.path(input_root, "MouseColitis_Inhibitors_Immune_Level3DeepResolved.rds")
stromal_file <- file.path(input_root, "MouseColitis_Inhibitors_Stromal_Level3Reannotated.rds")
manifest_file <- "MouseColitis_Inhibitors_manifest.csv"
resource_dir <- "FAP_ablation_DSS_fibroblast_myeloid_MultiNicheNet/resources"
lr_file <- file.path(resource_dir, "lr_network_mouse_allInfo_30112033.rds")
lt_file <- file.path(resource_dir, "ligand_target_matrix_nsga2r_final_mouse.rds")
stopifnot(all(file.exists(c(immune_file, stromal_file, manifest_file, lr_file, lt_file))))

manifest <- fread(manifest_file)[
  (dataset_source == "Merged_a5B1" & source_condition %chin% c("Control", "DSS+ a5B1_inhibitor")) |
    (sample_id %chin% paste0("D", 1:6) & source_condition == "FAP-TK_PBS_DSS")
]
manifest[, group := fcase(
  source_condition == "Control", "Control",
  sample_id %chin% paste0("D", 1:6), "DSS",
  source_condition == "DSS+ a5B1_inhibitor", "Blockade",
  default = NA_character_
)]
manifest <- manifest[!is.na(group)]
manifest[, condition_label := factor(group, levels = c("Control", "DSS", "Blockade"),
                                     labels = c("Control", "DSS colitis", "DSS + alpha5beta1 blockade"))]
manifest[, sample_uid := make.names(sample_uid)]

extract_compartment <- function(file, celltype_name) {
  object <- readRDS(file)
  metadata <- as.data.table(object@meta.data, keep.rownames = "cell")
  metadata[, sample_uid_match := make.names(as.character(sample_uid))]
  metadata <- metadata[sample_uid_match %in% manifest$sample_uid]
  if (celltype_name == "Fibroblast") {
    metadata <- metadata[reassigned_level_2 == "Fibroblast"]
  } else {
    metadata <- metadata[reassigned_level_2 == "Neutrophil"]
  }
  metadata <- merge(
    metadata,
    manifest[, .(sample_uid_match = sample_uid, group, condition_label)],
    by = "sample_uid_match", all.x = TRUE, sort = FALSE
  )
  metadata[, `:=`(
    celltype = celltype_name,
    level3 = make.names(as.character(reassigned_level_3)),
    sample = sample_uid_match
  )]
  object <- JoinLayers(object, assay = "RNA")
  counts <- LayerData(object[["RNA"]], layer = "counts")
  counts <- counts[, metadata$cell, drop = FALSE]
  metadata <- metadata[match(colnames(counts), cell)]
  stopifnot(identical(metadata$cell, colnames(counts)))
  list(counts = counts, metadata = metadata)
}

message("Extracting fibroblasts and neutrophils from finalized annotations")
fibro <- extract_compartment(stromal_file, "Fibroblast")
neut <- extract_compartment(immune_file, "Neutrophil")
common_genes <- intersect(rownames(fibro$counts), rownames(neut$counts))
counts <- cbind(fibro$counts[common_genes, , drop = FALSE], neut$counts[common_genes, , drop = FALSE])
counts@x <- round(counts@x)
counts <- drop0(counts)
metadata <- rbindlist(list(fibro$metadata, neut$metadata), fill = TRUE)
metadata <- metadata[match(colnames(counts), cell)]
stopifnot(identical(metadata$cell, colnames(counts)))
nonzero_cells <- Matrix::colSums(counts) > 0
counts <- counts[, nonzero_cells, drop = FALSE]
metadata <- metadata[nonzero_cells]
stopifnot(identical(metadata$cell, colnames(counts)))

sce_all <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(as.data.frame(metadata[, .(
    cell, sample, sample_id, group, celltype, level3, condition_label
  )]), row.names = metadata$cell)
)
rownames(sce_all) <- make.names(rownames(sce_all), unique = TRUE)
saveRDS(sce_all, file.path(out_dir, "objects", "fibroblast_neutrophil_input_sce.rds"), compress = FALSE)

cell_counts <- as.data.table(as.data.frame(colData(sce_all)))[, .N, by = .(
  sample, sample_id, group, condition_label, celltype, level3
)]
fwrite(cell_counts, file.path(out_dir, "tables", "input_cell_counts_by_sample_and_level3.csv"))
fwrite(cell_counts[, .(cells = sum(N)), by = .(
  sample, sample_id, group, condition_label, celltype
)], file.path(out_dir, "tables", "input_cell_counts_by_sample_and_celltype.csv"))

message("Preparing official mouse NicheNet priors")
lr_network_all <- readRDS(lr_file) |>
  mutate(
    ligand = make.names(convert_alias_to_symbols(ligand, organism = "mouse")),
    receptor = make.names(convert_alias_to_symbols(receptor, organism = "mouse"))
  )
lr_network <- lr_network_all |> distinct(ligand, receptor)
ligand_target_matrix <- readRDS(lt_file)
colnames(ligand_target_matrix) <- make.names(convert_alias_to_symbols(colnames(ligand_target_matrix), organism = "mouse"))
rownames(ligand_target_matrix) <- make.names(convert_alias_to_symbols(rownames(ligand_target_matrix), organism = "mouse"))
lr_network <- lr_network |> filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix <- ligand_target_matrix[, unique(lr_network$ligand), drop = FALSE]

run_pair <- function(comparator, comparison_id) {
  groups <- c("DSS", comparator)
  contrast_forward <- paste0("DSS-", comparator)
  contrast_reverse <- paste0(comparator, "-DSS")
  # Blockade-associated fibroblasts are sparse. Requiring >1 cell retains three
  # blockade units with both lineages after zero-library cells are removed.
  min_cells_pair <- 1
  pair_meta <- as.data.table(as.data.frame(colData(sce_all)))[group %in% groups]
  pair_counts <- pair_meta[, .N, by = .(sample, group, celltype)]
  eligible_pair_samples <- pair_counts[N > min_cells_pair, .N, by = .(sample, group)][N == 2, sample]
  sce <- sce_all[, colData(sce_all)$group %in% groups & colData(sce_all)$sample %in% eligible_pair_samples]
  fwrite(
    pair_counts[, retained_for_pair := sample %in% eligible_pair_samples],
    file.path(out_dir, "tables", paste0(comparison_id, "_sample_celltype_eligibility.csv"))
  )
  message("Running ", comparison_id, ": ", contrast_forward)
  result <- multi_nichenet_analysis(
    sce = sce,
    celltype_id = "celltype", sample_id = "sample", group_id = "group",
    batches = NA, covariates = NA,
    lr_network = lr_network, ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c(paste0("'", contrast_forward, "','", contrast_reverse, "'")),
    contrast_tbl = tibble(
      contrast = c(contrast_forward, contrast_reverse),
      group = c("DSS", comparator)
    ),
    senders_oi = c("Fibroblast", "Neutrophil"),
    receivers_oi = c("Fibroblast", "Neutrophil"),
    fraction_cutoff = 0.05, min_sample_prop = 0.50,
    scenario = "regular", ligand_activity_down = FALSE,
    min_cells = min_cells_pair, logFC_threshold = 0.50, p_val_threshold = 0.05,
    p_val_adj = FALSE, empirical_pval = FALSE,
    top_n_target = 250, verbose = TRUE, n.cores = 2,
    return_lr_prod_matrix = TRUE, top_n_LR = 2500
  )
  saveRDS(result, file.path(out_dir, "objects", paste0("multinichenet_", comparison_id, ".rds")), compress = FALSE)
  group_tbl <- as.data.table(result$prioritization_tables$group_prioritization_tbl)
  sample_tbl <- as.data.table(result$prioritization_tables$sample_prioritization_tbl)
  fwrite(group_tbl, file.path(out_dir, "tables", paste0(comparison_id, "_group_prioritization_all.csv")))
  fwrite(sample_tbl, file.path(out_dir, "tables", paste0(comparison_id, "_sample_prioritization_all.csv")))
  fwrite(as.data.table(result$celltype_de), file.path(out_dir, "tables", paste0(comparison_id, "_celltype_DE.csv")))

  cross <- group_tbl[sender != receiver]
  cross[, direction := paste(sender, "->", receiver)]
  wide <- dcast(
    cross,
    sender + receiver + direction + ligand + receptor + lr_interaction + id ~ group,
    value.var = c("prioritization_score", "fraction_expressing_ligand_receptor")
  )
  score_a <- paste0("prioritization_score_", "DSS")
  score_b <- paste0("prioritization_score_", comparator)
  frac_a <- paste0("fraction_expressing_ligand_receptor_", "DSS")
  frac_b <- paste0("fraction_expressing_ligand_receptor_", comparator)
  wide[, `:=`(
    comparison = comparison_id,
    comparator = comparator,
    interaction = paste0(ligand, " - ", receptor),
    score_DSS = get(score_a),
    score_comparator = get(score_b),
    fraction_DSS = get(frac_a),
    fraction_comparator = get(frac_b)
  )]
  wide[, score_delta := score_DSS - score_comparator]
  wide[, effect := fifelse(score_delta > 0, "Higher in DSS", "Lower in DSS")]
  wide[, abs_delta := abs(score_delta)]
  setorder(wide, -abs_delta)
  wide[, rank_by_direction_effect := frank(-abs_delta, ties.method = "first"), by = .(direction, effect)]
  fwrite(wide, file.path(out_dir, "tables", paste0(comparison_id, "_cross_lineage_score_changes.csv")))
  wide
}

reuse_core <- identical(Sys.getenv("A5B1_MN_REUSE_CORE"), "1")
disease_table_file <- file.path(out_dir, "tables", "DSS_vs_Control_cross_lineage_score_changes.csv")
blockade_table_file <- file.path(out_dir, "tables", "DSS_vs_Blockade_cross_lineage_score_changes.csv")
if (file.exists(disease_table_file)) {
  message("Reusing completed DSS-versus-Control MultiNicheNet result")
  disease <- fread(disease_table_file)
} else {
  disease <- run_pair("Control", "DSS_vs_Control")
}
if (reuse_core && file.exists(blockade_table_file)) {
  message("Reusing completed DSS-versus-blockade MultiNicheNet result")
  blockade <- fread(blockade_table_file)
} else {
  blockade <- run_pair("Blockade", "DSS_vs_Blockade")
}
combined <- rbindlist(list(disease, blockade), fill = TRUE)
combined[, interaction := paste0(ligand, " - ", receptor)]

ambient_epithelial <- c(
  "Epcam", "Cldn7", "Cldn2", "Cftr", "Guca2a", "Dsc2", "Dsg2", "Ocln",
  "Krt8", "Krt18", "Krt19", "Muc1", "Muc13", "Alpi", "Vil1"
)
combined[, ambient_epithelial_flag := ligand %in% ambient_epithelial | receptor %in% ambient_epithelial]
focused <- combined[ambient_epithelial_flag == FALSE]
fwrite(combined, file.path(out_dir, "tables", "both_comparisons_all_cross_lineage_interactions.csv"))
fwrite(focused, file.path(out_dir, "tables", "both_comparisons_nonambient_interactions.csv"))
fwrite(focused[rank_by_direction_effect <= 25], file.path(out_dir, "tables", "top25_by_comparison_direction_effect.csv"))

comparison_labels <- c(
  DSS_vs_Control = "DSS colitis vs Control\n(DSS effect)",
  DSS_vs_Blockade = "DSS colitis vs blockade\n(alpha5beta1-blockade effect)"
)
direction_colors <- c(
  "Fibroblast -> Neutrophil" = "#0072B2",
  "Neutrophil -> Fibroblast" = "#CC79A7"
)

top_plot_data <- focused[rank_by_direction_effect <= 7]
top_plot_data[, display_group := paste(comparison, direction, effect, sep = "__")]
top_plot_data[, interaction_ordered := reorder(interaction, score_delta)]

p_rank <- ggplot(top_plot_data, aes(score_delta, interaction_ordered, color = direction)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = score_delta, yend = interaction_ordered), linewidth = 0.7) +
  geom_point(size = 2.3) +
  facet_grid(effect + direction ~ comparison, scales = "free", space = "free_y",
             labeller = labeller(comparison = comparison_labels)) +
  scale_color_manual(values = direction_colors, guide = "none") +
  labs(
    title = "Largest changes in reciprocal neutrophil-fibroblast signaling",
    subtitle = "Positive values indicate higher MultiNicheNet priority in DSS colitis",
    x = "DSS colitis minus comparator priority score", y = NULL, tag = "A"
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold", size = 7.2),
    axis.text.y = element_text(size = 6.5),
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.25),
    plot.tag = element_text(face = "bold", size = 13)
  )

selected <- focused[rank_by_direction_effect <= 5, .(
  comparison, comparator, sender, receiver, direction, ligand, receptor, interaction, effect
)]
long_scores <- rbindlist(list(
  merge(
    combined[, .(comparison, comparator, sender, receiver, direction, ligand, receptor,
                 score = score_DSS, fraction = fraction_DSS)],
    selected, by = c("comparison", "comparator", "sender", "receiver", "direction", "ligand", "receptor")
  )[, condition := "DSS colitis"],
  merge(
    combined[, .(comparison, comparator, sender, receiver, direction, ligand, receptor,
                 score = score_comparator, fraction = fraction_comparator)],
    selected, by = c("comparison", "comparator", "sender", "receiver", "direction", "ligand", "receptor")
  )[, condition := fifelse(comparator == "Control", "Control", "DSS + alpha5beta1 blockade")]
), fill = TRUE)
long_scores[, condition := factor(
  condition,
  levels = c("Control", "DSS colitis", "DSS + alpha5beta1 blockade")
)]
long_scores[, interaction_ordered := reorder(interaction, score)]

p_support <- ggplot(long_scores, aes(condition, interaction_ordered)) +
  geom_point(aes(size = fraction, color = score), alpha = 0.95) +
  facet_grid(direction + effect ~ comparison, scales = "free", space = "free_y",
             labeller = labeller(comparison = comparison_labels)) +
  scale_color_viridis_c(option = "C", limits = c(0, 1), name = "Priority score") +
  scale_size_continuous(range = c(1.0, 4.8), labels = percent_format(accuracy = 1),
                        name = "Fraction expressing\nligand and receptor") +
  labs(
    title = "Condition-specific support for leading interactions",
    subtitle = "Color integrates differential expression, specificity, prevalence, and ligand activity",
    x = NULL, y = NULL, tag = "B"
  ) +
  theme_classic(base_size = 8.2) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold", size = 7),
    axis.text.y = element_text(size = 6.2), axis.text.x = element_text(size = 7.5),
    panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.22),
    legend.position = "bottom", plot.tag = element_text(face = "bold", size = 13)
  )

counts_plot <- cell_counts[, .(cells = sum(N)), by = .(sample, group, celltype)]
counts_plot[, group_label := factor(
  group, levels = c("Control", "DSS", "Blockade"),
  labels = c("Control", "DSS colitis", "DSS + alpha5beta1 blockade")
)]
p_counts <- ggplot(counts_plot, aes(group_label, cells, color = group_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.4, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.06, seed = 20260816), size = 1.8, alpha = 0.9) +
  facet_wrap(~celltype, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c(
    "Control" = "#4C78A8", "DSS colitis" = "#D95F02",
    "DSS + alpha5beta1 blockade" = "#1B9E77"
  ), guide = "none") +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), labels = comma) +
  labs(title = "Cell recovery by condition", subtitle = "All retained units with recovered lineage cells", x = NULL,
       y = "Cells retained", tag = "C") +
  theme_classic(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "#F2F2F2", color = NA),
        strip.text = element_text(face = "bold"), plot.tag = element_text(face = "bold", size = 13))

figure <- p_rank / (p_support | p_counts) +
  plot_layout(heights = c(1.25, 1), widths = c(3.2, 1)) +
  plot_annotation(
    title = "DSS colitis and alpha5beta1 blockade reshape reciprocal neutrophil-fibroblast communication",
    subtitle = "MultiNicheNet 2.1.0; pairwise analyses require recovered fibroblasts and neutrophils",
    caption = paste(
      "Scores prioritize expression-supported ligand-receptor hypotheses and do not measure signaling flux.",
      "Three blockade units pass the >1-cell threshold for both lineages after zero-library filtering; blockade-associated fibroblast recovery is sparse and condition is aligned with source batch."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 10, color = "#4D4D4D"),
      plot.caption = element_text(size = 7.5, hjust = 0, color = "#4D4D4D"),
      plot.margin = margin(10, 12, 8, 12)
    )
  )

png_file <- file.path(out_dir, "figures", "a5B1_DSS_neutrophil_fibroblast_MultiNicheNet_two_contrasts.png")
pdf_file <- file.path(out_dir, "figures", "a5B1_DSS_neutrophil_fibroblast_MultiNicheNet_two_contrasts.pdf")
ggsave(png_file, figure, width = 18, height = 15, dpi = 320, bg = "white", limitsize = FALSE)
ggsave(pdf_file, figure, width = 18, height = 15, device = cairo_pdf, bg = "white", limitsize = FALSE)

disease_eligibility <- fread(file.path(out_dir, "tables", "DSS_vs_Control_sample_celltype_eligibility.csv"))
blockade_eligibility <- fread(file.path(out_dir, "tables", "DSS_vs_Blockade_sample_celltype_eligibility.csv"))
disease_pair_units <- unique(disease_eligibility[retained_for_pair == TRUE, .(sample, group)])
blockade_pair_units <- unique(blockade_eligibility[retained_for_pair == TRUE, .(sample, group)])

validation <- data.table(
  package_version = as.character(packageVersion("multinichenetr")),
  cells = ncol(sce_all), genes = nrow(sce_all),
  fibroblast_cells = sum(colData(sce_all)$celltype == "Fibroblast"),
  neutrophil_cells = sum(colData(sce_all)$celltype == "Neutrophil"),
  eligible_sample_units = uniqueN(colData(sce_all)$sample),
  retained_control_units = uniqueN(colData(sce_all)$sample[colData(sce_all)$group == "Control"]),
  retained_dss_units = uniqueN(colData(sce_all)$sample[colData(sce_all)$group == "DSS"]),
  retained_blockade_units = uniqueN(colData(sce_all)$sample[colData(sce_all)$group == "Blockade"]),
  disease_pair_control_units = disease_pair_units[group == "Control", uniqueN(sample)],
  disease_pair_dss_units = disease_pair_units[group == "DSS", uniqueN(sample)],
  blockade_pair_dss_units = blockade_pair_units[group == "DSS", uniqueN(sample)],
  blockade_pair_blockade_units = blockade_pair_units[group == "Blockade", uniqueN(sample)],
  disease_interactions = nrow(disease), blockade_interactions = nrow(blockade),
  png_exists = file.exists(png_file), pdf_exists = file.exists(pdf_file),
  passed = nrow(disease) > 0 && nrow(blockade) > 0 &&
    file.exists(png_file) && file.exists(pdf_file) &&
    !any(manifest$sample_id %chin% c("C7", "C8", "C9", "C10"))
)
fwrite(validation, file.path(out_dir, "analysis_validation.csv"))
message("Two-contrast a5B1 neutrophil-fibroblast MultiNicheNet analysis complete")
