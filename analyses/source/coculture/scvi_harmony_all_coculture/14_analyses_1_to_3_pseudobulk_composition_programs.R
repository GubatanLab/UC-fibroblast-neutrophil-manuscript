#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260816)

project_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
neut_dir <- file.path(project_dir, "scVI_Harmony_All_Coculture", "Neutrophil_Subclustering")
analysis_dir <- file.path(project_dir, "Manuscript_Analyses_1_to_6")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
object_dir <- file.path(analysis_dir, "objects")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

contrast_key <- data.frame(
  contrast_id = c("N_vs_PRE", "CF_vs_N", "UF_vs_N", "UF_vs_CF", "NNAMPT_vs_N",
                  "NA5_vs_N", "UN_vs_UF", "UA5_vs_UF", "UD_vs_UF"),
  reference = c("PRE", "N", "N", "CF", "N", "N", "UF", "UF", "UF"),
  comparison = c("N", "CF", "UF", "UF", "NNAMPT", "NA5", "UN", "UA5", "UD"),
  contrast = c(
    "Neutrophils alone vs pre-culture",
    "Control fibroblast vs neutrophils alone",
    "UC fibroblast vs neutrophils alone",
    "UC vs control fibroblast",
    "NAMPTi vs neutrophils alone",
    "alpha5beta1i vs neutrophils alone",
    "NAMPTi vs UC", "alpha5beta1i vs UC", "Dual blockade vs UC"
  ),
  paired = c(FALSE, rep(TRUE, 8)),
  stringsAsFactors = FALSE
)
did_key <- data.frame(
  contrast_id = c("NAMPT_fibroblast_interaction", "alpha5beta1_fibroblast_interaction"),
  contrast = c(
    "Fibroblast-dependent NAMPTi effect",
    "Fibroblast-dependent alpha5beta1i effect"
  ),
  expression = c(
    "ConditionCodeUN-ConditionCodeUF-ConditionCodeNNAMPT+ConditionCodeN",
    "ConditionCodeUA5-ConditionCodeUF-ConditionCodeNA5+ConditionCodeN"
  ),
  conditions = c("N,NNAMPT,UF,UN", "N,NA5,UF,UA5"),
  stringsAsFactors = FALSE
)
fwrite(contrast_key, file.path(source_dir, "analysis_contrasts.csv"))
fwrite(did_key, file.path(source_dir, "analysis_interaction_contrasts.csv"))

aggregate_sparse_counts <- function(obj, state_col, panel, min_cells = 10L) {
  md <- obj@meta.data
  stopifnot(all(c(state_col, "DonorID", "ConditionCode") %in% colnames(md)))
  state <- as.character(md[[state_col]])
  donor <- as.character(md$DonorID)
  condition <- as.character(md$ConditionCode)
  valid <- !is.na(state) & !is.na(donor) & !is.na(condition)
  obj <- obj[, valid]
  state <- state[valid]
  donor <- donor[valid]
  condition <- condition[valid]
  group <- factor(paste(state, donor, condition, sep = "||"))
  indicator <- sparseMatrix(
    i = seq_along(group), j = as.integer(group), x = 1,
    dims = c(length(group), nlevels(group)),
    dimnames = list(NULL, levels(group))
  )
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  aggregated <- counts %*% indicator
  parts <- tstrsplit(colnames(aggregated), "\\|\\|", fixed = FALSE)
  group_meta <- data.table(
    group = colnames(aggregated), state = parts[[1]], DonorID = parts[[2]],
    ConditionCode = parts[[3]], cells = as.integer(colSums(indicator)), panel = panel
  )
  keep <- group_meta$cells >= min_cells
  list(counts = aggregated[, keep, drop = FALSE], meta = group_meta[keep])
}

run_edger_pair <- function(counts, meta, state_name, key_row) {
  sm <- meta[state == state_name & ConditionCode %in% c(key_row$reference, key_row$comparison)]
  if (key_row$paired) {
    complete_donors <- sm[, .(n_conditions = uniqueN(ConditionCode)), by = DonorID][n_conditions == 2, DonorID]
    sm <- sm[DonorID %in% complete_donors]
  }
  group_sizes <- sm[, .N, by = ConditionCode]
  if (nrow(group_sizes) < 2 || min(group_sizes$N) < 3) return(NULL)
  y <- DGEList(as.matrix(counts[, sm$group, drop = FALSE]))
  cond <- relevel(factor(sm$ConditionCode), ref = key_row$reference)
  if (key_row$paired) {
    design <- model.matrix(~ factor(sm$DonorID) + cond)
  } else {
    design <- model.matrix(~ cond)
  }
  keep_gene <- filterByExpr(y, design = design, min.count = 3)
  if (sum(keep_gene) < 20) return(NULL)
  y <- calcNormFactors(y[keep_gene, , keep.lib.sizes = FALSE])
  fit <- glmQLFit(y, design, robust = TRUE)
  tab <- topTags(glmQLFTest(fit, coef = ncol(design)), n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$state <- state_name
  tab$panel <- unique(sm$panel)
  tab$contrast_id <- key_row$contrast_id
  tab$contrast <- key_row$contrast
  tab$reference <- key_row$reference
  tab$comparison <- key_row$comparison
  tab$paired <- key_row$paired
  tab$n_reference <- group_sizes[ConditionCode == key_row$reference, N]
  tab$n_comparison <- group_sizes[ConditionCode == key_row$comparison, N]
  if (key_row$paired) tab$n_pairs <- uniqueN(sm$DonorID) else tab$n_pairs <- NA_integer_
  as.data.table(tab)
}

run_edger_did <- function(counts, meta, state_name, key_row) {
  conditions <- strsplit(key_row$conditions, ",", fixed = TRUE)[[1]]
  sm <- meta[state == state_name & ConditionCode %in% conditions]
  complete_donors <- sm[, .(n_conditions = uniqueN(ConditionCode)), by = DonorID][n_conditions == 4, DonorID]
  sm <- sm[DonorID %in% complete_donors]
  if (uniqueN(sm$DonorID) < 3) return(NULL)
  sm$ConditionCode <- factor(sm$ConditionCode, levels = conditions)
  design <- model.matrix(~ 0 + ConditionCode + factor(DonorID), data = sm)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design)) return(NULL)
  y <- DGEList(as.matrix(counts[, sm$group, drop = FALSE]))
  keep_gene <- filterByExpr(y, design = design, min.count = 3)
  if (sum(keep_gene) < 20) return(NULL)
  y <- calcNormFactors(y[keep_gene, , keep.lib.sizes = FALSE])
  fit <- glmQLFit(y, design, robust = TRUE)
  contrast_vec <- makeContrasts(contrasts = key_row$expression, levels = design)
  tab <- topTags(glmQLFTest(fit, contrast = contrast_vec), n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$state <- state_name
  tab$panel <- unique(sm$panel)
  tab$contrast_id <- key_row$contrast_id
  tab$contrast <- key_row$contrast
  tab$reference <- "inhibitor-alone and untreated backgrounds"
  tab$comparison <- "fibroblast-dependent inhibitor interaction"
  tab$paired <- TRUE
  tab$n_reference <- NA_integer_
  tab$n_comparison <- NA_integer_
  tab$n_pairs <- uniqueN(sm$DonorID)
  as.data.table(tab)
}

run_pseudobulk_panel <- function(bundle) {
  states <- unique(bundle$meta$state)
  ordinary <- rbindlist(lapply(states, function(st) {
    rbindlist(lapply(seq_len(nrow(contrast_key)), function(i) {
      run_edger_pair(bundle$counts, bundle$meta, st, contrast_key[i, ])
    }), fill = TRUE)
  }), fill = TRUE)
  interactions <- rbindlist(lapply(states, function(st) {
    rbindlist(lapply(seq_len(nrow(did_key)), function(i) {
      run_edger_did(bundle$counts, bundle$meta, st, did_key[i, ])
    }), fill = TRUE)
  }), fill = TRUE)
  rbindlist(list(ordinary, interactions), fill = TRUE)
}

message("Loading neutrophil objects for donor-aware pseudobulk analysis")
original <- readRDS(file.path(neut_dir, "objects", "neutrophils_scVI_Harmony_subclustered_annotated_clean.rds"))
canonical <- readRDS(file.path(neut_dir, "objects", "canonical_neutrophils_scVI_Harmony_reclustered.rds"))
original <- subset(
  original,
  cells = rownames(original@meta.data)[original@meta.data$neutrophil_cluster_annotation != "Low-signal neutrophil"]
)

original_bundle <- aggregate_sparse_counts(
  original, "neutrophil_cluster_annotation", "Figure A detailed subclusters", min_cells = 10L
)
canonical_bundle <- aggregate_sparse_counts(
  canonical, "neutrophil_state_annotation", "Figure B canonical states", min_cells = 10L
)

message("Running state-specific edgeR pseudobulk models")
de_results <- rbindlist(list(
  run_pseudobulk_panel(original_bundle),
  run_pseudobulk_panel(canonical_bundle)
), fill = TRUE)
de_results[, FDR_global := p.adjust(PValue, method = "BH"), by = panel]
fwrite(de_results, file.path(source_dir, "analysis1_state_specific_pseudobulk_DE.csv.gz"))

de_summary <- de_results[, .(
  tested_genes = .N,
  DE_FDR_0_05 = sum(FDR < 0.05),
  DE_FDR_0_05_logFC_0_5 = sum(FDR < 0.05 & abs(logFC) >= 0.5),
  up_FDR_0_05_logFC_0_5 = sum(FDR < 0.05 & logFC >= 0.5),
  down_FDR_0_05_logFC_0_5 = sum(FDR < 0.05 & logFC <= -0.5),
  median_significant_logFC = if (any(FDR < 0.05)) median(logFC[FDR < 0.05]) else NA_real_,
  minimum_FDR = min(FDR, na.rm = TRUE)
), by = .(panel, state, contrast_id, contrast, n_pairs, n_reference, n_comparison)]
fwrite(de_summary, file.path(source_dir, "analysis1_state_specific_pseudobulk_DE_summary.csv"))

# Generic limma testing for sample-level feature matrices.
run_limma_pair <- function(value_matrix, sample_meta, key_row, analysis_label) {
  sm <- sample_meta[ConditionCode %in% c(key_row$reference, key_row$comparison)]
  if (key_row$paired) {
    complete_donors <- sm[, .(n_conditions = uniqueN(ConditionCode)), by = DonorID][n_conditions == 2, DonorID]
    sm <- sm[DonorID %in% complete_donors]
  }
  gs <- sm[, .N, by = ConditionCode]
  if (nrow(gs) < 2 || min(gs$N) < 3) return(NULL)
  cond <- relevel(factor(sm$ConditionCode), ref = key_row$reference)
  if (key_row$paired) design <- model.matrix(~ factor(sm$DonorID) + cond) else design <- model.matrix(~ cond)
  fit <- eBayes(lmFit(value_matrix[, sm$sample_id, drop = FALSE], design), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
  tab$feature <- rownames(tab)
  tab$analysis <- analysis_label
  tab$contrast_id <- key_row$contrast_id
  tab$contrast <- key_row$contrast
  tab$n_pairs <- if (key_row$paired) uniqueN(sm$DonorID) else NA_integer_
  as.data.table(tab)
}

run_limma_did <- function(value_matrix, sample_meta, key_row, analysis_label) {
  conditions <- strsplit(key_row$conditions, ",", fixed = TRUE)[[1]]
  sm <- sample_meta[ConditionCode %in% conditions]
  complete_donors <- sm[, .(n_conditions = uniqueN(ConditionCode)), by = DonorID][n_conditions == 4, DonorID]
  sm <- sm[DonorID %in% complete_donors]
  if (uniqueN(sm$DonorID) < 3) return(NULL)
  sm$ConditionCode <- factor(sm$ConditionCode, levels = conditions)
  design <- model.matrix(~ 0 + ConditionCode + factor(DonorID), data = sm)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design)) return(NULL)
  fit <- lmFit(value_matrix[, sm$sample_id, drop = FALSE], design)
  fit <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = key_row$expression, levels = design)), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  tab$feature <- rownames(tab)
  tab$analysis <- analysis_label
  tab$contrast_id <- key_row$contrast_id
  tab$contrast <- key_row$contrast
  tab$n_pairs <- uniqueN(sm$DonorID)
  as.data.table(tab)
}

run_limma_all <- function(value_matrix, sample_meta, analysis_label) {
  a <- rbindlist(lapply(seq_len(nrow(contrast_key)), function(i) {
    run_limma_pair(value_matrix, sample_meta, contrast_key[i, ], analysis_label)
  }), fill = TRUE)
  b <- rbindlist(lapply(seq_len(nrow(did_key)), function(i) {
    run_limma_did(value_matrix, sample_meta, did_key[i, ], analysis_label)
  }), fill = TRUE)
  rbindlist(list(a, b), fill = TRUE)
}

# Analysis 2: donor-condition centered-log-ratio composition of canonical states.
canonical_md <- as.data.table(canonical@meta.data, keep.rownames = "cell")
composition_counts <- canonical_md[, .N, by = .(
  DonorID, ConditionCode, state = neutrophil_state_annotation
)]
sample_totals <- canonical_md[, .(total_cells = .N), by = .(DonorID, ConditionCode)]
all_states <- sort(unique(composition_counts$state))
composition_complete <- CJ(
  DonorID = unique(sample_totals$DonorID),
  ConditionCode = unique(sample_totals$ConditionCode),
  state = all_states,
  unique = TRUE
)[sample_totals, on = .(DonorID, ConditionCode), nomatch = 0]
composition_complete <- composition_counts[
  composition_complete, on = .(DonorID, ConditionCode, state)
]
composition_complete[is.na(N), N := 0L]
composition_complete[, proportion := N / total_cells]
composition_complete[, sample_id := paste(DonorID, ConditionCode, sep = "||")]
composition_complete[, clr := log(N + 0.5) - mean(log(N + 0.5)), by = sample_id]
fwrite(composition_complete, file.path(source_dir, "analysis2_canonical_state_composition_by_sample.csv"))

composition_wide <- dcast(composition_complete, state ~ sample_id, value.var = "clr")
composition_matrix <- as.matrix(composition_wide[, -1])
rownames(composition_matrix) <- composition_wide$state
composition_sample_meta <- unique(composition_complete[, .(sample_id, DonorID, ConditionCode, total_cells)])
composition_results <- run_limma_all(
  composition_matrix, composition_sample_meta, "Canonical-state CLR composition"
)
composition_results[, state := feature]
fwrite(composition_results, file.path(source_dir, "analysis2_canonical_state_CLR_differential_abundance.csv"))

# Analysis 3: curated neutrophil functional programs, scored from log-normalized RNA.
programs <- list(
  NETosis = c("PADI4", "MPO", "ELANE", "PRTN3", "CTSG", "AZU1", "CYBB", "NCF1", "NCF2", "S100A8", "S100A9"),
  Degranulation = c("LTF", "BPI", "CAMP", "LCN2", "MMP8", "MMP9", "OLFM4", "CEACAM8", "FCGR3B"),
  Oxidative_burst = c("CYBB", "CYBA", "NCF1", "NCF2", "NCF4", "RAC2", "MPO"),
  Phagocytosis = c("FCGR3B", "FCGR2A", "ITGAM", "ITGB2", "FPR1", "FPR2", "C5AR1", "TREM1", "SYK", "LYN"),
  Chemotaxis = c("CXCR1", "CXCR2", "CXCR4", "FPR1", "FPR2", "C5AR1", "CCR1", "CCRL2"),
  OSM_inflammation = c("OSM", "IL1B", "TNF", "NFKBIA", "CXCL8", "CCL3", "CCL4", "IL1R2", "SOCS3"),
  Type_I_interferon = c("MX1", "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "OAS1", "OAS2", "OASL", "RSAD2", "IRF7", "STAT1"),
  Immature_granulopoiesis = c("LTF", "BPI", "CAMP", "LCN2", "OLFM4", "CEACAM8", "DEFA4", "MS4A3"),
  Survival_antiapoptosis = c("BCL2A1", "MCL1", "BCL2L1", "CFLAR", "FAS", "BAX", "BBC3", "CASP3"),
  Hypoxia_glycolysis = c("HIF1A", "SLC2A1", "HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "LDHA", "ENO1"),
  ER_stress = c("XBP1", "ATF4", "DDIT3", "HSPA5", "HERPUD1", "DNAJB9", "ERO1A"),
  Proliferation_stress = c("MKI67", "PCNA", "TYMS", "DHFR", "TUBA1B", "STMN1"),
  Integrin_adhesion = c("ITGA5", "ITGB1", "ITGAM", "ITGB2", "SELPLG", "ICAM1"),
  CXCR4_aging_retention = c("CXCR4", "SELL", "FCGR3B", "PECAM1", "C5AR1", "BCL2A1")
)

data_matrix <- GetAssayData(canonical, assay = "RNA", layer = "data")
available_programs <- lapply(programs, intersect, y = rownames(data_matrix))
available_programs <- available_programs[lengths(available_programs) >= 3]
program_gene_table <- rbindlist(lapply(names(available_programs), function(nm) {
  data.table(program = nm, gene = available_programs[[nm]])
}))
fwrite(program_gene_table, file.path(source_dir, "analysis3_curated_functional_program_genes.csv"))

union_genes <- unique(program_gene_table$gene)
expr_union <- as.matrix(data_matrix[union_genes, , drop = FALSE])
gene_means <- rowMeans(expr_union)
gene_sds <- apply(expr_union, 1, sd)
gene_sds[gene_sds == 0 | is.na(gene_sds)] <- 1
expr_z <- sweep(sweep(expr_union, 1, gene_means, "-"), 1, gene_sds, "/")
program_scores <- sapply(available_programs, function(gs) colMeans(expr_z[gs, , drop = FALSE]))
program_scores <- as.data.table(program_scores)
program_scores[, cell := colnames(canonical)]
program_scores <- cbind(
  canonical_md[, .(cell, DonorID, ConditionCode, state = neutrophil_state_annotation)],
  program_scores[, -"cell"]
)
program_long <- melt(
  program_scores,
  id.vars = c("cell", "DonorID", "ConditionCode", "state"),
  variable.name = "program", value.name = "score"
)
program_sample <- program_long[, .(
  cells = .N, mean_score = mean(score), median_score = median(score)
), by = .(DonorID, ConditionCode, state, program)]
program_sample <- program_sample[cells >= 10]
program_sample[, sample_id := paste(DonorID, ConditionCode, sep = "||")]
fwrite(program_sample, file.path(source_dir, "analysis3_functional_program_scores_by_sample_state.csv"))

program_feature <- paste(program_sample$state, program_sample$program, sep = "||")
program_sample[, feature := paste(state, program, sep = "||")]
program_wide <- dcast(program_sample, feature ~ sample_id, value.var = "mean_score")
program_matrix <- as.matrix(program_wide[, -1])
rownames(program_matrix) <- program_wide$feature
program_sample_meta <- unique(program_sample[, .(sample_id, DonorID, ConditionCode)])
program_results <- run_limma_all(program_matrix, program_sample_meta, "State-specific functional programs")
program_results[, c("state", "program") := tstrsplit(feature, "\\|\\|", fixed = FALSE)]
fwrite(program_results, file.path(source_dir, "analysis3_functional_program_differential_activity.csv"))

# Manuscript figures for analyses 1-3.
contrast_order <- c(contrast_key$contrast, did_key$contrast)
canonical_order <- c("PADI4 neutrophil", "OSM neutrophil", "CXCR4 neutrophil", "MX1/ISG neutrophil")
detail_order <- c(
  "Low-signal neutrophil", "Activated OSM/CXCR4 neutrophil",
  "IL1R2+ inflammatory neutrophil", "DHFR+ proliferative/stress neutrophil",
  "PADI4/translation-high neutrophil", "LTF/BPI immature neutrophil",
  "CCL3/CCL4 inflammatory neutrophil", "RORA+ atypical neutrophil"
)
theme_pub <- theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(7, 7, 7, 7)
  )

plot_de_summary <- function(panel_name, state_order, title) {
  dd <- de_summary[panel == panel_name]
  dd[, state := factor(state, levels = rev(state_order))]
  dd[, contrast := factor(contrast, levels = contrast_order)]
  dd[, signed_score := sign(up_FDR_0_05_logFC_0_5 - down_FDR_0_05_logFC_0_5) *
       log10(1 + up_FDR_0_05_logFC_0_5 + down_FDR_0_05_logFC_0_5)]
  dd[, label := paste0(up_FDR_0_05_logFC_0_5, "/", down_FDR_0_05_logFC_0_5)]
  ggplot(dd, aes(contrast, state, fill = signed_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(
      title = title,
      subtitle = "Tile text: up/down genes with FDR <0.05 and |log2FC| >=0.5",
      x = NULL, y = NULL, fill = "Signed\nDE burden"
    ) + theme_pub
}

p_de_a <- plot_de_summary(
  "Figure A detailed subclusters", detail_order,
  "Differential state: detailed neutrophil subclusters"
)
p_de_b <- plot_de_summary(
  "Figure B canonical states", canonical_order,
  "Differential state: canonical neutrophil states"
)
fig1 <- p_de_a / p_de_b + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

composition_results[, state := factor(state, levels = rev(canonical_order))]
composition_results[, contrast := factor(contrast, levels = contrast_order)]
composition_results[, star := fifelse(adj.P.Val < 0.001, "***", fifelse(adj.P.Val < 0.01, "**", fifelse(adj.P.Val < 0.05, "*", "")))]
p_comp <- ggplot(composition_results, aes(contrast, state, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = star), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Donor-level compositional differential abundance",
    subtitle = "Centered-log-ratio state abundance; * adjusted P <0.05",
    x = NULL, y = NULL, fill = "CLR effect"
  ) + theme_pub

program_results[, state := factor(state, levels = rev(canonical_order))]
program_results[, contrast := factor(contrast, levels = contrast_order)]
program_results[, program := gsub("_", " ", program)]
program_results[, sig := adj.P.Val < 0.05]
p_program <- ggplot(program_results, aes(contrast, program, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_point(data = program_results[sig == TRUE], shape = 21, size = 1.8, fill = "black", color = "black") +
  facet_wrap(~state, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Neutrophil functional-program responses",
    subtitle = "State-specific donor-level effects; dots mark adjusted P <0.05",
    x = NULL, y = NULL, fill = "Score effect"
  ) + theme_pub + theme(axis.text.y = element_text(size = 8))

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}

save_plot(fig1, "Analysis1_state_specific_pseudobulk_DE", 16, 12)
save_plot(p_comp, "Analysis2_compositional_differential_abundance", 14, 5.8)
save_plot(p_program, "Analysis3_functional_program_activity", 16, 10)

saveRDS(
  list(
    contrast_key = contrast_key, did_key = did_key,
    programs = available_programs,
    composition_matrix = composition_matrix,
    program_matrix = program_matrix
  ),
  file.path(object_dir, "analyses_1_to_3_supporting_objects.rds"),
  compress = TRUE
)

message("Analyses 1-3 complete")
