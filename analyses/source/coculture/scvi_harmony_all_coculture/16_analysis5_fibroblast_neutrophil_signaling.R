suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
})

set.seed(160825)
project_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
analysis_dir <- file.path(project_dir, "Manuscript_Analyses_1_to_6")
source_dir <- file.path(analysis_dir, "source_data")
figure_dir <- file.path(analysis_dir, "figures")
object_dir <- file.path(analysis_dir, "objects")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

fib <- readRDS(file.path(project_dir, "scVI_Harmony_All_Coculture/Fibroblast_Subclustering/objects/fibroblasts_subclustered_annotated.rds"))
neu <- readRDS(file.path(project_dir, "scVI_Harmony_All_Coculture/Neutrophil_Subclustering/objects/canonical_neutrophils_scVI_Harmony_reclustered.rds"))
fib <- subset(
  fib,
  cells = rownames(fib@meta.data)[!grepl("immune-like", fib@meta.data$fibroblast_annotation, ignore.case = TRUE)]
)

data(CellChatDB.human)
lr <- as.data.table(CellChatDB.human$interaction)
lr <- unique(lr[, .(
  interaction_name, pathway_name, annotation,
  ligand_components = ligand.symbol,
  receptor_components = receptor.symbol
)])
split_components <- function(x) trimws(unlist(strsplit(x, ",", fixed = TRUE)))
ligand_genes <- unique(unlist(lapply(lr$ligand_components, split_components)))
receptor_genes <- unique(unlist(lapply(lr$receptor_components, split_components)))

fib_mat <- LayerData(fib, assay = "DECONTX", layer = "data")
neu_mat <- LayerData(neu, assay = "RNA", layer = "data")
ligand_genes <- intersect(ligand_genes, rownames(fib_mat))
receptor_genes <- intersect(receptor_genes, rownames(neu_mat))
lr <- lr[vapply(ligand_components, function(x) all(split_components(x) %in% ligand_genes), logical(1)) &
           vapply(receptor_components, function(x) all(split_components(x) %in% receptor_genes), logical(1))]

summarize_expression <- function(mat, metadata, genes, group_cols, cell_name = "cell") {
  metadata <- as.data.table(metadata, keep.rownames = cell_name)
  groups <- interaction(metadata[, ..group_cols], drop = TRUE, lex.order = TRUE, sep = "||")
  group_levels <- levels(groups)
  indicator <- sparseMatrix(
    i = seq_along(groups), j = as.integer(groups), x = 1,
    dims = c(length(groups), length(group_levels)),
    dimnames = list(rownames(metadata), group_levels)
  )
  x <- mat[genes, metadata[[cell_name]], drop = FALSE]
  n <- as.numeric(table(groups)[group_levels])
  means <- as.matrix(x %*% indicator) / rep(n, each = length(genes))
  pct <- as.matrix((x > 0) %*% indicator) / rep(n, each = length(genes))
  out <- as.data.table(as.table(means))
  setnames(out, c("gene", "group", "mean_expression"))
  pct_dt <- as.data.table(as.table(pct))
  setnames(pct_dt, c("gene", "group", "pct_expressed"))
  out <- merge(out, pct_dt, by = c("gene", "group"))
  group_meta <- unique(metadata[, c(group_cols), with = FALSE])
  group_meta[, group := interaction(.SD, drop = TRUE, lex.order = TRUE, sep = "||"), .SDcols = group_cols]
  group_meta[, cells := n]
  merge(out, group_meta, by = "group")
}

fib_meta <- fib@meta.data
fib_meta$fibroblast_annotation <- droplevels(factor(fib_meta$fibroblast_annotation))
neu_meta <- neu@meta.data
neu_meta$neutrophil_state_annotation <- droplevels(factor(neu_meta$neutrophil_state_annotation))

lig_expr <- summarize_expression(
  fib_mat, fib_meta, ligand_genes,
  c("ConditionCode", "fibroblast_annotation")
)
rec_expr <- summarize_expression(
  neu_mat, neu_meta, receptor_genes,
  c("ConditionCode", "neutrophil_state_annotation")
)
fwrite(lig_expr, file.path(source_dir, "analysis5_fibroblast_ligand_expression.csv.gz"))
fwrite(rec_expr, file.path(source_dir, "analysis5_neutrophil_receptor_expression.csv.gz"))

complex_signal <- function(expr_dt, component_string, condition, celltype_col, celltype) {
  components <- split_components(component_string)
  z <- expr_dt[ConditionCode == condition & get(celltype_col) == celltype & gene %in% components]
  if (nrow(z) != length(components)) return(c(mean = 0, pct = 0, cells = 0))
  c(mean = min(z$mean_expression), pct = min(z$pct_expressed), cells = min(z$cells))
}

conditions <- intersect(unique(lig_expr$ConditionCode), unique(rec_expr$ConditionCode))
fib_types <- unique(lig_expr$fibroblast_annotation)
neu_states <- unique(rec_expr$neutrophil_state_annotation)
score_rows <- vector("list", length(conditions) * length(fib_types) * length(neu_states) * nrow(lr))
k <- 0L
for (condition in conditions) for (ft in fib_types) for (ns in neu_states) {
  for (i in seq_len(nrow(lr))) {
    ls <- complex_signal(lig_expr, lr$ligand_components[i], condition, "fibroblast_annotation", ft)
    rs <- complex_signal(rec_expr, lr$receptor_components[i], condition, "neutrophil_state_annotation", ns)
    if (ls["cells"] < 3 || rs["cells"] < 20 || ls["pct"] < 0.05 || rs["pct"] < 0.05) next
    k <- k + 1L
    score_rows[[k]] <- data.table(
      ConditionCode = condition, fibroblast_subtype = ft, neutrophil_state = ns,
      interaction_name = lr$interaction_name[i], pathway_name = lr$pathway_name[i],
      annotation = lr$annotation[i], ligand = lr$ligand_components[i], receptor = lr$receptor_components[i],
      ligand_mean = ls["mean"], ligand_pct = ls["pct"], receptor_mean = rs["mean"], receptor_pct = rs["pct"],
      fibroblast_cells = ls["cells"], neutrophil_cells = rs["cells"],
      interaction_score = sqrt((ls["mean"] * ls["pct"]) * (rs["mean"] * rs["pct"]))
    )
  }
}
scores <- rbindlist(score_rows[seq_len(k)], fill = TRUE)
fwrite(scores, file.path(source_dir, "analysis5_condition_subtype_state_LR_scores.csv.gz"))

contrast_key <- data.table(
  contrast_id = c("UF_vs_CF", "UN_vs_UF", "UA5_vs_UF", "UD_vs_UF"),
  reference = c("CF", "UF", "UF", "UF"), comparison = c("UF", "UN", "UA5", "UD"),
  contrast = c("UC vs control fibroblast", "NAMPTi vs UC", "alpha5beta1i vs UC", "Dual blockade vs UC")
)
id_cols <- c("fibroblast_subtype", "neutrophil_state", "interaction_name", "pathway_name", "annotation", "ligand", "receptor")
diff_list <- lapply(seq_len(nrow(contrast_key)), function(i) {
  ck <- contrast_key[i]
  ref <- scores[ConditionCode == ck$reference, c(id_cols, "interaction_score"), with = FALSE]
  cmp <- scores[ConditionCode == ck$comparison, c(id_cols, "interaction_score"), with = FALSE]
  setnames(ref, "interaction_score", "score_reference")
  setnames(cmp, "interaction_score", "score_comparison")
  z <- merge(ref, cmp, by = id_cols, all = TRUE)
  z[is.na(score_reference), score_reference := 0]
  z[is.na(score_comparison), score_comparison := 0]
  z[, `:=`(
    delta_score = score_comparison - score_reference,
    log2_fold_change = log2((score_comparison + 0.01) / (score_reference + 0.01)),
    contrast_id = ck$contrast_id, contrast = ck$contrast
  )]
  z
})
diff_scores <- rbindlist(diff_list, fill = TRUE)
fwrite(diff_scores, file.path(source_dir, "analysis5_differential_LR_scores.csv.gz"))

# Couple changing fibroblast ligands to unbiased neutrophil response modules.
module_effects <- fread(file.path(source_dir, "analysis4_module_differential_activity.csv"))
module_effects <- module_effects[contrast_id %in% contrast_key$contrast_id]
ligand_effects <- diff_scores[, .(
  ligand_delta = weighted.mean(delta_score, w = pmax(score_reference + score_comparison, 0.01)),
  max_abs_LR_delta = max(abs(delta_score)), pathways = paste(unique(pathway_name), collapse = ";")
), by = .(contrast_id, contrast, neutrophil_state, ligand)]
coupling <- merge(
  ligand_effects,
  module_effects[, .(contrast_id, neutrophil_state = state, module, module_display, module_logFC = logFC, module_FDR = adj.P.Val)],
  by = c("contrast_id", "neutrophil_state"), allow.cartesian = TRUE
)
coupling[, response_coupling := ligand_delta * module_logFC]
fwrite(coupling, file.path(source_dir, "analysis5_ligand_to_neutrophil_module_response_coupling.csv"))

theme_pub <- theme_bw(base_family = "Arial", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9.5, color = "#444444"),
        axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
        legend.title = element_text(face = "bold"), plot.margin = margin(7, 7, 7, 7))

top_diff <- diff_scores[, head(.SD[order(-abs(delta_score))], 8), by = contrast_id]
top_diff[, interaction_label := paste0(ligand, " -> ", receptor)]
top_diff[, route := paste(fibroblast_subtype, neutrophil_state, sep = " -> ")]
top_diff[, contrast := factor(contrast, levels = contrast_key$contrast)]
p_a <- ggplot(top_diff, aes(contrast, interaction_label, size = abs(delta_score), color = delta_score)) +
  geom_point(alpha = 0.9) + facet_wrap(~neutrophil_state, scales = "free_y", ncol = 2) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(title = "Fibroblast-to-neutrophil signaling changes", subtitle = "CellChatDB-supported interactions; positive values increase in the comparison condition",
       x = NULL, y = NULL, size = "|Score change|", color = "Score change") + theme_pub + theme(axis.text.y = element_text(size = 7))

pathway_summary <- diff_scores[, .(pathway_change = sum(delta_score), interaction_count = .N),
                               by = .(contrast, fibroblast_subtype, neutrophil_state, pathway_name)]
pathway_summary[, route_pathway := paste(fibroblast_subtype, "->", neutrophil_state, "|", pathway_name)]
pathway_top <- pathway_summary[, head(.SD[order(-abs(pathway_change))], 8), by = contrast]
p_b <- ggplot(pathway_top, aes(contrast, route_pathway, size = interaction_count, color = pathway_change)) +
  geom_point() +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(title = "Sender subtype and receiver-state pathway map", subtitle = "Aggregated ligand-receptor score changes",
       x = NULL, y = NULL, size = "Interactions", color = "Pathway change") + theme_pub +
  theme(axis.text.y = element_text(size = 6.5), strip.text.y = element_text(size = 7))

top_coupling <- coupling[is.finite(response_coupling), head(.SD[order(-abs(response_coupling))], 12), by = contrast_id]
top_coupling[, pair := paste(ligand, module_display, sep = " -> ")]
top_coupling[, contrast := factor(contrast, levels = contrast_key$contrast)]
p_c <- ggplot(top_coupling, aes(contrast, pair, fill = response_coupling)) +
  geom_tile(color = "white", linewidth = 0.25) + facet_wrap(~neutrophil_state, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Ligand-to-neutrophil module response coupling", subtitle = "Descriptive coupling of ligand-receptor change with donor-level module response",
       x = NULL, y = NULL, fill = "Coupling") + theme_pub + theme(axis.text.y = element_text(size = 6.5))

fig5 <- p_a / p_b / p_c + plot_layout(heights = c(1.15, 1.35, 1.05)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
for (ext in c("pdf", "png", "tiff")) {
  output_file <- file.path(figure_dir, paste0("Analysis5_fibroblast_neutrophil_signaling.", ext))
  if (ext == "tiff") {
    ggsave(output_file, fig5, width = 14, height = 18, dpi = 400, compression = "lzw")
  } else {
    ggsave(output_file, fig5, width = 14, height = 18, dpi = 400)
  }
}
saveRDS(list(scores = scores, differential_scores = diff_scores, response_coupling = coupling),
        file.path(object_dir, "analysis5_fibroblast_neutrophil_signaling.rds"), compress = TRUE)
message("Analysis 5 fibroblast-neutrophil signaling complete")
