#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(edgeR)
  library(limma)
  library(igraph)
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

contrast_key <- fread(file.path(source_dir, "analysis_contrasts.csv"))
did_key <- fread(file.path(source_dir, "analysis_interaction_contrasts.csv"))
support <- readRDS(file.path(object_dir, "analyses_1_to_3_supporting_objects.rds"))

canonical <- readRDS(file.path(neut_dir, "objects", "canonical_neutrophils_scVI_Harmony_reclustered.rds"))
md <- as.data.table(canonical@meta.data, keep.rownames = "cell")
md[, state := as.character(neutrophil_state_annotation)]

# Construct approximately 50-cell metacells within donor-condition-state strata.
md[, random_order := sample(.N), by = .(DonorID, ConditionCode, state)]
setorder(md, DonorID, ConditionCode, state, random_order)
md[, metacell_index := ceiling(seq_len(.N) / 50), by = .(DonorID, ConditionCode, state)]
md[, metacell_id := paste(DonorID, ConditionCode, state, metacell_index, sep = "||")]
metacell_sizes <- md[, .N, by = metacell_id]
valid_metacells <- metacell_sizes[N >= 20, metacell_id]
md <- md[metacell_id %in% valid_metacells]

cell_index <- match(md$cell, colnames(canonical))
group <- factor(md$metacell_id)
indicator <- sparseMatrix(
  i = cell_index, j = as.integer(group), x = 1,
  dims = c(ncol(canonical), nlevels(group)),
  dimnames = list(colnames(canonical), levels(group))
)
counts <- GetAssayData(canonical, assay = "RNA", layer = "counts")
metacell_counts <- counts %*% indicator
metacell_meta <- unique(md[, .(
  metacell_id, DonorID, ConditionCode, state, metacell_index
)])
metacell_meta[, cells := metacell_sizes[.SD, on = "metacell_id", x.N]]
setkey(metacell_meta, metacell_id)
metacell_meta <- metacell_meta[colnames(metacell_counts)]

y <- DGEList(as.matrix(metacell_counts))
y <- calcNormFactors(y)
log_cpm <- cpm(y, log = TRUE, prior.count = 2)

detected_fraction <- rowMeans(metacell_counts >= 3)
gene_variance <- apply(log_cpm, 1, var)
lineage_contaminant_genes <- c(
  "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "IL7R", "LTB",
  "GNLY", "NKG7", "KLRD1", "CCL5", "GZMB", "PRF1",
  "CD79A", "CD79B", "MS4A1", "CD37", "MZB1", "JCHAIN",
  "PPBP", "PF4", "NRGN", "HBB", "HBA1", "HBA2"
)
eligible <- detected_fraction >= 0.10 &
  !grepl("^(MT-|RPL|RPS|ENSG)", rownames(log_cpm)) &
  !rownames(log_cpm) %in% lineage_contaminant_genes
eligible_genes <- names(sort(gene_variance[eligible], decreasing = TRUE))
network_genes <- head(eligible_genes, 900)
network_expression <- log_cpm[network_genes, , drop = FALSE]
network_expression <- network_expression[apply(network_expression, 1, sd) > 0, , drop = FALSE]

message("Calculating metacell gene-gene correlations")
gene_cor <- cor(t(network_expression), method = "pearson", use = "pairwise.complete.obs")
diag(gene_cor) <- -Inf
gene_names <- rownames(gene_cor)
edges <- rbindlist(lapply(seq_along(gene_names), function(i) {
  ix <- head(order(gene_cor[i, ], decreasing = TRUE), 8)
  data.table(from = gene_names[i], to = gene_names[ix], weight = gene_cor[i, ix])
}))
edges <- edges[is.finite(weight) & weight >= 0.25 & from != to]
edges[, pair := ifelse(from < to, paste(from, to, sep = "||"), paste(to, from, sep = "||"))]
edges <- edges[, .(from = first(from), to = first(to), weight = max(weight)), by = pair]

graph <- graph_from_data_frame(edges[, .(from, to, weight)], directed = FALSE, vertices = gene_names)
communities <- cluster_louvain(graph, weights = E(graph)$weight)
membership_dt <- data.table(gene = names(membership(communities)), raw_module = as.integer(membership(communities)))
module_sizes <- membership_dt[, .N, by = raw_module][order(-N)]
valid_raw_modules <- module_sizes[N >= 10, raw_module]
module_map <- setNames(paste0("M", seq_along(valid_raw_modules)), valid_raw_modules)
membership_dt[, module := ifelse(raw_module %in% valid_raw_modules, module_map[as.character(raw_module)], "grey")]
membership_dt[, module := factor(module, levels = c(paste0("M", seq_along(valid_raw_modules)), "grey"))]

# Module eigengenes and intramodular hub correlations.
valid_modules <- setdiff(levels(membership_dt$module), "grey")
module_eigengenes <- matrix(
  NA_real_, nrow = length(valid_modules), ncol = ncol(network_expression),
  dimnames = list(valid_modules, colnames(network_expression))
)
hub_tables <- vector("list", length(valid_modules))
for (i in seq_along(valid_modules)) {
  mod <- valid_modules[i]
  genes <- membership_dt[module == mod, gene]
  z <- t(scale(t(network_expression[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  pc <- prcomp(t(z), center = FALSE, scale. = FALSE)$x[, 1]
  if (cor(pc, colMeans(z)) < 0) pc <- -pc
  module_eigengenes[mod, ] <- pc
  hubs <- apply(z, 1, function(v) cor(v, pc, use = "pairwise.complete.obs"))
  hub_tables[[i]] <- data.table(module = mod, gene = names(hubs), hub_correlation = as.numeric(hubs))
}
hub_table <- rbindlist(hub_tables)
membership_dt <- merge(membership_dt, hub_table, by = c("module", "gene"), all.x = TRUE)

# Annotate modules by overlap with the curated functional programs.
programs <- support$programs
background <- unique(membership_dt[module != "grey", gene])
enrichment <- rbindlist(lapply(valid_modules, function(mod) {
  mod_genes <- membership_dt[module == mod, gene]
  rbindlist(lapply(names(programs), function(program) {
    gs <- intersect(programs[[program]], background)
    overlap <- intersect(mod_genes, gs)
    p <- phyper(length(overlap) - 1, length(gs), length(background) - length(gs), length(mod_genes), lower.tail = FALSE)
    data.table(
      module = mod, program = program, module_genes = length(mod_genes),
      program_genes = length(gs), overlap = length(overlap),
      overlap_genes = paste(overlap, collapse = ";"), p_value = p
    )
  }))
}))
enrichment[, FDR := p.adjust(p_value, method = "BH")]
module_labels <- enrichment[order(module, FDR, -overlap), .SD[1], by = module]
module_labels[, display_label := paste0(module, ": ", gsub("_", " ", program))]

membership_dt <- merge(
  membership_dt,
  module_labels[, .(module, module_annotation = program, module_display = display_label)],
  by = "module", all.x = TRUE
)
fwrite(membership_dt, file.path(source_dir, "analysis4_metacell_gene_module_membership_and_hubs.csv"))
fwrite(enrichment, file.path(source_dir, "analysis4_module_functional_enrichment.csv"))
fwrite(edges[, .(from, to, weight)], file.path(source_dir, "analysis4_gene_coexpression_network_edges.csv.gz"))

# Aggregate module eigengenes to donor-condition-state before differential testing.
me_long <- as.data.table(t(module_eigengenes), keep.rownames = "metacell_id")
me_long <- melt(me_long, id.vars = "metacell_id", variable.name = "module", value.name = "eigengene")
me_long <- merge(me_long, metacell_meta, by = "metacell_id")
me_sample <- me_long[, .(
  metacells = .N, cells = sum(cells), module_score = mean(eigengene)
), by = .(DonorID, ConditionCode, state, module)]
me_sample[, sample_id := paste(DonorID, ConditionCode, sep = "||")]
me_sample[, feature := paste(state, module, sep = "||")]
fwrite(me_sample, file.path(source_dir, "analysis4_module_eigengenes_by_sample_state.csv"))

me_wide <- dcast(me_sample, feature ~ sample_id, value.var = "module_score")
me_matrix <- as.matrix(me_wide[, -1])
rownames(me_matrix) <- me_wide$feature
sample_meta <- unique(me_sample[, .(sample_id, DonorID, ConditionCode)])

run_pair <- function(key_row) {
  sm <- sample_meta[ConditionCode %in% c(key_row$reference, key_row$comparison)]
  if (key_row$paired) {
    donors <- sm[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 2, DonorID]
    sm <- sm[DonorID %in% donors]
  }
  gs <- sm[, .N, by = ConditionCode]
  if (nrow(gs) < 2 || min(gs$N) < 3) return(NULL)
  cond <- relevel(factor(sm$ConditionCode), key_row$reference)
  if (key_row$paired) design <- model.matrix(~ factor(sm$DonorID) + cond) else design <- model.matrix(~ cond)
  test_matrix <- me_matrix[, sm$sample_id, drop = FALSE]
  keep <- rowSums(!is.finite(test_matrix)) == 0 & apply(test_matrix, 1, sd) > 1e-08
  if (!any(keep) || nrow(design) <= ncol(design)) return(NULL)
  fit <- eBayes(lmFit(test_matrix[keep, , drop = FALSE], design), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
  tab$feature <- rownames(tab)
  tab$contrast_id <- key_row$contrast_id
  tab$contrast <- key_row$contrast
  tab$n_pairs <- if (key_row$paired) uniqueN(sm$DonorID) else NA_integer_
  as.data.table(tab)
}

run_did <- function(key_row) {
  conditions <- strsplit(key_row$conditions, ",", fixed = TRUE)[[1]]
  sm <- sample_meta[ConditionCode %in% conditions]
  donors <- sm[, .(n = uniqueN(ConditionCode)), by = DonorID][n == 4, DonorID]
  sm <- sm[DonorID %in% donors]
  if (uniqueN(sm$DonorID) < 3) return(NULL)
  sm$ConditionCode <- factor(sm$ConditionCode, levels = conditions)
  design <- model.matrix(~ 0 + ConditionCode + factor(DonorID), data = sm)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design)) return(NULL)
  test_matrix <- me_matrix[, sm$sample_id, drop = FALSE]
  keep <- rowSums(!is.finite(test_matrix)) == 0 & apply(test_matrix, 1, sd) > 1e-08
  if (!any(keep) || nrow(design) <= ncol(design)) return(NULL)
  fit <- lmFit(test_matrix[keep, , drop = FALSE], design)
  contrast_matrix <- makeContrasts(contrasts = key_row$expression, levels = design)
  fit <- eBayes(contrasts.fit(fit, contrast_matrix), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  tab$feature <- rownames(tab)
  tab$contrast_id <- key_row$contrast_id
  tab$contrast <- key_row$contrast
  tab$n_pairs <- uniqueN(sm$DonorID)
  as.data.table(tab)
}

module_results <- rbindlist(list(
  rbindlist(lapply(seq_len(nrow(contrast_key)), function(i) run_pair(contrast_key[i]))),
  rbindlist(lapply(seq_len(nrow(did_key)), function(i) run_did(did_key[i])))
), fill = TRUE)
module_results[, c("state", "module") := tstrsplit(feature, "\\|\\|", fixed = FALSE)]
module_results <- merge(module_results, module_labels[, .(module, module_display = display_label)], by = "module", all.x = TRUE)
fwrite(module_results, file.path(source_dir, "analysis4_module_differential_activity.csv"))

# Manuscript figure: functional annotation, condition responses, and hub genes.
theme_pub <- theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(7, 7, 7, 7)
  )

module_order <- module_labels$display_label
enrichment_plot <- merge(enrichment, module_labels[, .(module, module_display = display_label)], by = "module")
enrichment_plot <- enrichment_plot[overlap > 0]
enrichment_plot[, program_label := gsub("_", " ", program)]
enrichment_plot[, module_display := factor(module_display, levels = rev(module_order))]
p_a <- ggplot(enrichment_plot, aes(program_label, module_display, size = overlap, color = -log10(pmax(FDR, 1e-12)))) +
  geom_point() +
  scale_color_gradient(low = "#FDD49E", high = "#B30000") +
  labs(
    title = "Unbiased metacell modules and functional annotation",
    subtitle = "Point size is overlapping genes; color is enrichment strength",
    x = NULL, y = NULL, size = "Genes", color = "-log10 FDR"
  ) + theme_pub

canonical_order <- c("PADI4 neutrophil", "OSM neutrophil", "CXCR4 neutrophil", "MX1/ISG neutrophil")
contrast_order <- c(contrast_key$contrast, did_key$contrast)
module_results[, state := factor(state, levels = canonical_order)]
module_results[, contrast := factor(contrast, levels = contrast_order)]
module_results[, module_display := factor(module_display, levels = rev(module_order))]
module_results[, significant := adj.P.Val < 0.05]
p_b <- ggplot(module_results, aes(contrast, module_display, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_point(data = module_results[significant == TRUE], size = 1.4, color = "black") +
  facet_wrap(~state, ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Condition-dependent module activity",
    subtitle = "Donor-level module eigengenes; dots mark adjusted P <0.05",
    x = NULL, y = NULL, fill = "Effect"
  ) + theme_pub + theme(axis.text.y = element_text(size = 7.5))

hub_display <- membership_dt[module != "grey"][order(module, -hub_correlation), head(.SD, 5), by = module]
hub_display[, rank := seq_len(.N), by = module]
hub_display[, module_display := factor(module_display, levels = rev(module_order))]
p_c <- ggplot(hub_display, aes(rank, module_display)) +
  geom_text(aes(label = gene, color = hub_correlation), fontface = "italic", size = 3.1) +
  scale_x_continuous(breaks = 1:5, labels = paste0("Hub ", 1:5)) +
  scale_color_gradient(low = "#969696", high = "#7A0177") +
  labs(
    title = "Top intramodular hub genes",
    subtitle = "Ranked by correlation with the module eigengene",
    x = NULL, y = NULL, color = "Hub correlation"
  ) + theme_pub + theme(panel.grid.major.y = element_line(color = "#EEEEEE"))

fig4 <- p_a / p_b / p_c + plot_layout(heights = c(0.8, 1.45, 0.8)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_plot(fig4, "Analysis4_metacell_gene_modules", 16, 18)

saveRDS(
  list(
    module_membership = membership_dt,
    module_eigengenes = module_eigengenes,
    metacell_metadata = metacell_meta,
    module_graph = graph,
    module_labels = module_labels
  ),
  file.path(object_dir, "analysis4_metacell_gene_modules.rds"),
  compress = TRUE
)

message("Analysis 4 metacell gene modules complete")
