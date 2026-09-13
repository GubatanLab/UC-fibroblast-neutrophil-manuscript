suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

set.seed(20260815)

trajectory_dir <- file.path(
  "a5B1_DSS_epithelial_immune_stromal_remodeling", "neutrophil_trajectory"
)
out_dir <- file.path(trajectory_dir, "program_subsets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

immune_file <- file.path(
  "MouseColitis_Inhibitors_compartment_subclustering",
  "reassigned_fractions_epithelial_audit_cleaned", "updated_umaps", "level3_reannotated",
  "MouseColitis_Inhibitors_Immune_Level3DeepResolved.rds"
)

condition_levels <- c("Control", "Severe DSS colitis", "DSS + alpha5beta1 blockade")
condition_colors <- c("Control" = "#4C78A8", "Severe DSS colitis" = "#D95F02", "DSS + alpha5beta1 blockade" = "#1B9E77")
program_levels <- c("Neutrophil OSM", "Neutrophil CXCR4", "Neutrophil PADI4", "Neutrophil MX1")
program_colors <- c(
  "Neutrophil OSM" = "#B2182B", "Neutrophil CXCR4" = "#2166AC",
  "Neutrophil PADI4" = "#E08214", "Neutrophil MX1" = "#762A83"
)

# Cross-species, reference-guided programs based on the existing human UC neutrophil
# subclustering and restricted to one-to-one mouse symbols present in this object.
program_definitions <- list(
  `Neutrophil OSM` = c("Osm", "Il1b", "Socs3", "Pde4b", "Hcar2", "Prok2", "Pfkfb3", "Nfkbia", "Cxcl2"),
  `Neutrophil CXCR4` = c("Cxcr4", "Hilpda", "Plin2", "Vegfa", "Hmox1", "Bhlhe40", "Sqstm1", "Ccrl2", "Cd83"),
  `Neutrophil PADI4` = c("Padi4", "Mmp9", "Cebpd", "Cebpb", "Cxcr2", "Sorl1", "Rgs2", "Camp", "Ngp"),
  `Neutrophil MX1` = c("Mx1", "Isg15", "Ifit1", "Ifit2", "Ifit3", "Rsad2", "Xaf1", "Gbp4", "Gbp5", "Ifi6")
)

cell_data <- fread(file.path(trajectory_dir, "neutrophil_cell_trajectory.csv.gz"))
cell_data[, condition := factor(condition, levels = condition_levels)]
sample_summary <- fread(file.path(trajectory_dir, "neutrophil_trajectory_by_mouse.csv"))
sample_summary[, condition := factor(condition, levels = condition_levels)]

message("Loading curated immune object and scoring neutrophil programs")
immune <- readRDS(immune_file)
feature_lookup <- setNames(rownames(immune), toupper(rownames(immune)))
mapped_programs <- lapply(program_definitions, function(x) {
  unique(na.omit(unname(feature_lookup[toupper(x)])))
})
gene_map <- rbindlist(lapply(names(mapped_programs), function(program) {
  data.table(program, requested_gene = program_definitions[[program]],
             mapped_gene = unname(feature_lookup[toupper(program_definitions[[program]])]))
}))
fwrite(gene_map, file.path(out_dir, "neutrophil_program_gene_map.csv"))
if (any(lengths(mapped_programs) < 4L)) {
  stop("Fewer than four mapped genes in at least one neutrophil program")
}

neut_object <- subset(immune, cells = cell_data$cell)
neut_object <- JoinLayers(neut_object, assay = "RNA")
counts <- GetAssayData(neut_object, assay = "RNA", layer = "counts")
totals <- Matrix::colSums(counts)
all_genes <- unique(unlist(mapped_programs))
expr <- log1p(t(t(as.matrix(counts[all_genes, cell_data$cell, drop = FALSE])) / pmax(totals[cell_data$cell], 1) * 1e4))
rownames(expr) <- all_genes
colnames(expr) <- cell_data$cell
expr_z <- t(scale(t(expr)))
expr_z[!is.finite(expr_z)] <- 0

program_scores <- sapply(mapped_programs, function(genes) {
  colMeans(expr_z[genes, , drop = FALSE])
})
program_scores <- as.matrix(program_scores)
rownames(program_scores) <- cell_data$cell
score_dt <- as.data.table(program_scores, keep.rownames = "cell")

score_matrix <- as.matrix(score_dt[, names(mapped_programs), with = FALSE])
top_index <- max.col(score_matrix, ties.method = "first")
sorted_scores <- t(apply(score_matrix, 1, sort, decreasing = TRUE))
score_dt[, `:=`(
  top_program = names(mapped_programs)[top_index],
  top_score = sorted_scores[, 1],
  score_margin = sorted_scores[, 1] - sorted_scores[, 2]
)]
# Assign each cell to the dominant one of the four requested neutrophil programs.
score_dt[, program := top_program]
score_dt[, program := factor(program, levels = program_levels)]

cell_program <- merge(cell_data, score_dt, by = "cell")
fwrite(cell_program, file.path(out_dir, "neutrophil_program_assignments.csv.gz"))

assignment_summary <- cell_program[, .(
  cells = .N, fraction = .N / nrow(cell_program),
  median_margin = median(score_margin), median_top_score = median(top_score)
), by = program][order(program)]
fwrite(assignment_summary, file.path(out_dir, "neutrophil_program_assignment_summary.csv"))

# Sample-unit-level fractions with explicit zeros for absent programs.
requested_programs <- program_levels
program_counts <- cell_program[, .(cells = .N), by = .(sample_uid, sample_id, condition, program)]
program_grid <- CJ(sample_uid = sample_summary$sample_uid, program = program_levels, unique = TRUE)
program_grid <- merge(
  program_grid,
  sample_summary[, .(sample_uid, sample_id, condition, neutrophil_cells)],
  by = "sample_uid"
)
program_by_mouse <- merge(
  program_grid, program_counts,
  by = c("sample_uid", "sample_id", "condition", "program"), all.x = TRUE
)
program_by_mouse[is.na(cells), cells := 0L]
program_by_mouse[, fraction := cells / neutrophil_cells]
fwrite(program_by_mouse, file.path(out_dir, "neutrophil_program_fractions_by_mouse.csv"))

program_tests <- program_by_mouse[program %chin% requested_programs, {
  ref <- fraction[condition == condition_levels[2]]
  trt <- fraction[condition == condition_levels[3]]
  wt <- wilcox.test(trt, ref, exact = FALSE)
  .(
    dss_sample_units = length(ref), blockade_sample_units = length(trt),
    dss_mean = mean(ref), blockade_mean = mean(trt),
    difference_percentage_points = 100 * (mean(trt) - mean(ref)),
    p_value = wt$p.value
  )
}, by = program]
program_tests[, FDR := p.adjust(p_value, method = "BH")]
program_tests[, significance := fifelse(FDR < 0.001, "***", fifelse(FDR < 0.01, "**",
                                      fifelse(FDR < 0.05, "*", "ns")))]
fwrite(program_tests, file.path(out_dir, "neutrophil_program_fraction_tests.csv"))

# Sample-unit median pseudotime within each program.
program_pseudotime <- cell_program[program %chin% requested_programs, .(
  cells = .N, median_pseudotime = median(pseudotime), mean_pseudotime = mean(pseudotime)
), by = .(sample_uid, sample_id, condition, program)]
fwrite(program_pseudotime, file.path(out_dir, "neutrophil_program_pseudotime_by_mouse.csv"))
fwrite(program_pseudotime, file.path(out_dir, "neutrophil_program_pseudotime_by_sample_unit.csv"))

program_pseudotime_tests <- program_pseudotime[, {
  ref <- median_pseudotime[condition == condition_levels[2]]
  trt <- median_pseudotime[condition == condition_levels[3]]
  wt <- wilcox.test(trt, ref, exact = FALSE)
  .(
    dss_sample_units = length(ref), blockade_sample_units = length(trt),
    dss_mean = mean(ref), blockade_mean = mean(trt),
    difference = mean(trt) - mean(ref), p_value = wt$p.value
  )
}, by = program]
program_pseudotime_tests[, FDR := p.adjust(p_value, method = "BH")]
program_pseudotime_tests[, significance := fifelse(FDR < 0.001, "***", fifelse(FDR < 0.01, "**",
                                                 fifelse(FDR < 0.05, "*", "ns")))]
program_pseudotime_tests <- merge(
  program_pseudotime_tests,
  program_pseudotime[, .(ymax = max(median_pseudotime)), by = program],
  by = "program"
)
program_pseudotime_tests[, `:=`(x = 2.5, y = pmin(1.04, ymax * 1.07 + 0.01))]
fwrite(program_pseudotime_tests, file.path(out_dir, "neutrophil_program_pseudotime_tests.csv"))

# Program-score dynamics use the existing six pseudotime bins and sample-unit-weighted means.
cell_program[, pseudotime_bin := cut(
  pseudotime, breaks = seq(0, 1, length.out = 7), include.lowest = TRUE,
  labels = paste0("B", 1:6)
)]
score_long <- melt(
  cell_program,
  id.vars = c("cell", "sample_uid", "condition", "pseudotime", "pseudotime_bin"),
  measure.vars = names(mapped_programs), variable.name = "score_program", value.name = "program_score"
)
score_long[, bin_midpoint := (as.integer(pseudotime_bin) - 0.5) / 6]
score_mouse_bin <- score_long[, .(program_score = mean(program_score)),
                              by = .(sample_uid, condition, score_program, pseudotime_bin, bin_midpoint)]
score_dynamics <- score_mouse_bin[, .(
  mean = mean(program_score), sem = sd(program_score) / sqrt(.N), mice = .N
), by = .(condition, score_program, pseudotime_bin, bin_midpoint)]
fwrite(score_dynamics, file.path(out_dir, "neutrophil_program_scores_across_pseudotime.csv"))

# Program composition by pseudotime bin, calculated within each mouse first.
bin_program_counts <- cell_program[, .(cells = .N),
  by = .(sample_uid, condition, pseudotime_bin, program)]
bin_totals <- cell_program[, .(bin_cells = .N), by = .(sample_uid, condition, pseudotime_bin)]
bin_program_counts <- merge(bin_program_counts, bin_totals,
                            by = c("sample_uid", "condition", "pseudotime_bin"))
bin_program_counts[, fraction := cells / bin_cells]
bin_program_summary <- bin_program_counts[, .(
  mean_fraction = mean(fraction), sem = sd(fraction) / sqrt(.N), mice = .N
), by = .(condition, pseudotime_bin, program)]
fwrite(bin_program_summary, file.path(out_dir, "neutrophil_program_composition_across_pseudotime.csv"))

theme_manuscript <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 7.5, color = "grey25"),
    legend.title = element_blank(), legend.text = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 7.3),
    panel.grid.minor = element_blank()
  )

p_umap <- ggplot(cell_program, aes(UMAP_1, UMAP_2, color = program)) +
  geom_point(size = 0.65, alpha = 0.78) +
  scale_color_manual(values = program_colors, drop = FALSE) +
  labs(title = "Reference-guided neutrophil programs",
       subtitle = "Dominant expression program; ambiguous cells retained as mixed") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 2, byrow = TRUE))

p_condition <- ggplot(cell_program, aes(UMAP_1, UMAP_2, color = program)) +
  geom_point(size = 0.6, alpha = 0.78) +
  facet_wrap(~condition, nrow = 1) +
  scale_color_manual(values = program_colors, drop = FALSE) +
  labs(title = "Program occupancy by treatment",
       subtitle = "Shared neutrophil UMAP and program definitions") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5), strip.text = element_text(face = "bold", size = 7),
        legend.position = "none")

plot_fraction_data <- program_by_mouse[program %chin% requested_programs]
ann <- merge(
  program_tests,
  plot_fraction_data[, .(ymax = max(fraction)), by = program],
  by = "program"
)
ann[, `:=`(x = 2.5, y = pmin(1.03, ymax * 1.12 + 0.01))]
p_fractions <- ggplot(plot_fraction_data, aes(condition, fraction, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.54, linewidth = 0.45) +
  geom_point(position = position_jitter(width = 0.07), size = 1.45, alpha = 0.88) +
  geom_text(data = ann, aes(x = x, y = y, label = significance), inherit.aes = FALSE,
            size = 3, fontface = "bold") +
  facet_wrap(~program, scales = "free_y", ncol = 2) +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0.03, 0.22))) +
  labs(title = "Neutrophil-program redistribution",
       subtitle = "Fractions within total neutrophils; each point is one sample unit",
       x = NULL, y = "Program fraction") +
  theme_manuscript +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_program_pt <- ggplot(program_pseudotime, aes(condition, median_pseudotime, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.54, linewidth = 0.45) +
  geom_point(aes(size = cells), position = position_jitter(width = 0.07), alpha = 0.85) +
  geom_text(data = program_pseudotime_tests,
            aes(x = x, y = y, label = significance), inherit.aes = FALSE,
            size = 3, fontface = "bold") +
  facet_wrap(~program, scales = "free_y", ncol = 2) +
  scale_color_manual(values = condition_colors) +
  scale_size_continuous(range = c(1, 2.7)) +
  labs(title = "Program position along inferred pseudotime",
       subtitle = "Sample-unit medians; point size reflects program-assigned cells",
       x = NULL, y = "Median inferred pseudotime", size = "Cells") +
  theme_manuscript +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1, size = 6.5))

p_dynamics <- ggplot(score_dynamics, aes(bin_midpoint, mean, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.14, color = NA, na.rm = TRUE) +
  geom_line(linewidth = 0.75) + geom_point(size = 1.25) +
  facet_wrap(~score_program, scales = "free_y", nrow = 1) +
  scale_color_manual(values = condition_colors) + scale_fill_manual(values = condition_colors) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("Cxcr2+", "Intermediate", "Inflammatory"),
                     limits = c(0, 1)) +
  labs(title = "Program activity across the inferred neutrophil continuum",
       subtitle = "Sample-unit-weighted mean +/- SEM",
       x = "Inferred pseudotime", y = "Standardized program score") +
  theme_manuscript + theme(legend.position = "bottom", axis.text.x = element_text(size = 6.5))

p_effect <- ggplot(program_tests, aes(difference_percentage_points, reorder(program, difference_percentage_points),
                                      color = program)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = difference_percentage_points,
                   yend = reorder(program, difference_percentage_points)), linewidth = 0.75) +
  geom_point(size = 2.5) +
  geom_text(aes(label = significance), color = "black", nudge_y = 0.15,
            size = 3, fontface = "bold") +
  scale_color_manual(values = program_colors) +
  labs(title = "Blockade-associated program effects",
       subtitle = "Positive values indicate enrichment with blockade",
  x = "Blockade - severe DSS colitis (percentage points)", y = NULL) +
  theme_manuscript + theme(legend.position = "none", panel.grid.major.y = element_blank())

layout <- c(
  area(t = 1, l = 1, b = 1, r = 2),
  area(t = 1, l = 3, b = 1, r = 6),
  area(t = 2, l = 1, b = 2, r = 3),
  area(t = 2, l = 4, b = 2, r = 6),
  area(t = 3, l = 1, b = 3, r = 2),
  area(t = 3, l = 3, b = 3, r = 6)
)

figure <- p_umap + p_condition + p_fractions + p_program_pt + p_effect + p_dynamics +
  plot_layout(design = layout, heights = c(1.45, 1.45, 1.2)) +
  plot_annotation(
    title = "DSS and alpha5beta1 blockade reshape neutrophil programs along pseudotime",
    subtitle = paste0(
      "Reference-guided program scoring of ", nrow(cell_program),
      " curated neutrophils positioned on the existing inferred pseudotime continuum"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9, color = "grey25"),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

base <- file.path(out_dir, "Figure_5_supplement_neutrophil_programs_pseudotime")
ggsave(paste0(base, ".png"), figure, width = 15, height = 13.5, dpi = 300,
       bg = "white", limitsize = FALSE)
ggsave(paste0(base, ".pdf"), figure, width = 15, height = 13.5,
       bg = "white", device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(base, ".tiff"), figure, width = 15, height = 13.5, dpi = 300,
       bg = "white", compression = "lzw", limitsize = FALSE)

legend <- c(
  "# Figure 5—figure supplement. Severe DSS colitis and alpha5beta1 blockade reshape Neutrophil OSM, CXCR4, PADI4, and MX1 subsets along inferred pseudotime.",
  "",
  "(A) Neutrophil-specific UMAP colored by the dominant reference-guided Neutrophil OSM, CXCR4, PADI4, or MX1 subset assignment.",
  "(B) Program occupancy in Control, DSS, and alpha5beta1-blockade conditions on the shared neutrophil UMAP.",
  "(C) Sample-unit-level fractions of OSM-, CXCR4-, PADI4-, and MX1-associated programs within total neutrophils. Stars denote Benjamini-Hochberg-adjusted FDR: * <0.05, ** <0.01, *** <0.001.",
  "(D) Median inferred pseudotime of each program in each included sample unit. Point size represents the number of program-assigned cells contributing to the median; stars denote Benjamini-Hochberg-adjusted FDR.",
  "(E) Treatment-associated differences in sample-unit-level program fractions.",
  "(F) Sample-unit-weighted mean program activity across six fixed bins of the previously inferred Cxcr2-to-inflammatory neutrophil continuum.",
  "",
  "Subsets were assigned from the dominant score among four programs transferred from the existing human UC neutrophil reference using mapped mouse orthologous markers. They should not be interpreted as independently validated discrete mouse lineages. Pseudotime is a cross-sectional state continuum, not observed temporal fate. Treatment remains aligned with source batch, so all effects are exploratory."
)
writeLines(legend, file.path(out_dir, "Figure_5_supplement_neutrophil_programs_pseudotime_legend.md"))

validation <- data.table(
  cells = nrow(cell_program), sample_units = uniqueN(cell_program$sample_uid),
  requested_programs_detected = sum(assignment_summary$program %chin% requested_programs & assignment_summary$cells > 0),
  png_exists = file.exists(paste0(base, ".png")), pdf_exists = file.exists(paste0(base, ".pdf")),
  tiff_exists = file.exists(paste0(base, ".tiff"))
)
validation[, passed := requested_programs_detected == 4L & png_exists & pdf_exists & tiff_exists]
fwrite(validation, file.path(out_dir, "program_pseudotime_validation.csv"))
message("Neutrophil program pseudotime analysis complete")
