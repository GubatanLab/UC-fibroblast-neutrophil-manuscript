suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggraph)
})

set.seed(20260820)
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

sig_obj <- readRDS(file.path(base_analysis, "objects", "analysis5_fibroblast_neutrophil_signaling.rds"))
scores <- as.data.table(sig_obj[["scores"]])
tf_effects <- fread(file.path(source_dir, "analysis2_regulator_pathway_differential_activity.csv"))

# Four complementary expression and detection metrics are combined as percentile
# ranks within each sender-receiver context. This prevents any single scoring
# convention from dominating and is robust to scale differences between metrics.
scores[, `:=`(
  metric_geometric = interaction_score,
  metric_expression_product = ligand_mean * receptor_mean,
  metric_detection_product = ligand_pct * receptor_pct,
  metric_balanced_min = pmin(ligand_mean * ligand_pct, receptor_mean * receptor_pct)
)]
metric_cols <- c("metric_geometric", "metric_expression_product",
                 "metric_detection_product", "metric_balanced_min")
for (m in metric_cols) {
  rank_col <- paste0(m, "_percentile")
  scores[, (rank_col) := frank(get(m), ties.method = "average") / .N,
         by = .(ConditionCode, fibroblast_subtype, neutrophil_state)]
}
rank_cols <- paste0(metric_cols, "_percentile")
scores[, consensus_score := rowMeans(.SD), .SDcols = rank_cols]
scores[, route := paste(fibroblast_subtype, neutrophil_state, sep = " -> ")]
fwrite(scores, file.path(source_dir, "analysis3_multimetric_consensus_LR_scores.csv.gz"))

contrast_key <- data.table(
  contrast_id = c("UF_vs_CF", "UN_vs_UF", "UA5_vs_UF", "UD_vs_UF"),
  reference = c("CF", "UF", "UF", "UF"), comparison = c("UF", "UN", "UA5", "UD"),
  contrast = c("UC vs control FB", "NAMPTi vs UC", "alpha5beta1i vs UC", "Dual vs UC")
)
id_cols <- c("fibroblast_subtype", "neutrophil_state", "interaction_name",
             "pathway_name", "annotation", "ligand", "receptor")
diff_list <- lapply(seq_len(nrow(contrast_key)), function(i) {
  key <- contrast_key[i]
  ref <- scores[ConditionCode == key$reference,
                c(id_cols, "consensus_score", "interaction_score"), with = FALSE]
  cmp <- scores[ConditionCode == key$comparison,
                c(id_cols, "consensus_score", "interaction_score"), with = FALSE]
  setnames(ref, c("consensus_score", "interaction_score"),
           c("consensus_reference", "expression_reference"))
  setnames(cmp, c("consensus_score", "interaction_score"),
           c("consensus_comparison", "expression_comparison"))
  z <- merge(ref, cmp, by = id_cols, all = TRUE)
  fill_cols <- c("consensus_reference", "expression_reference",
                 "consensus_comparison", "expression_comparison")
  for (cc in fill_cols) set(z, which(is.na(z[[cc]])), cc, 0)
  z[, `:=`(
    delta_consensus = consensus_comparison - consensus_reference,
    delta_expression = expression_comparison - expression_reference,
    mean_consensus = (consensus_comparison + consensus_reference) / 2,
    contrast_id = key$contrast_id, contrast = key$contrast
  )]
  z
})
diff_scores <- rbindlist(diff_list, fill = TRUE)
fwrite(diff_scores, file.path(source_dir, "analysis3_differential_consensus_LR_scores.csv.gz"))

# Couple consensus LR changes to the strongest donor-level TF responses in the
# matching neutrophil state. This is response coupling, not a causal path claim.
comm_ids <- contrast_key$contrast_id
tf_top <- tf_effects[activity_type == "TF" & contrast_id %in% comm_ids]
tf_top[, abs_effect := abs(effect)]
setorder(tf_top, contrast_id, state, FDR, -abs_effect)
tf_top <- tf_top[, head(.SD, 5), by = .(contrast_id, state)]
coupling <- merge(
  diff_scores,
  tf_top[, .(contrast_id, neutrophil_state = state, TF = regulator,
             TF_effect = effect, TF_FDR = FDR)],
  by = c("contrast_id", "neutrophil_state"), allow.cartesian = TRUE
)
coupling[, `:=`(
  directional_coupling = delta_consensus * TF_effect,
  route_confidence = mean_consensus * abs(TF_effect) * (1 + pmin(-log10(pmax(TF_FDR, 1e-12)), 6)),
  interaction_label = paste(ligand, receptor, sep = " -> ")
)]
fwrite(coupling, file.path(source_dir, "analysis3_consensus_LR_to_TF_response_coupling.csv.gz"))

contrast_labels <- setNames(contrast_key$contrast, contrast_key$contrast_id)
state_labels <- c("PADI4 neutrophil" = "PADI4", "OSM neutrophil" = "OSM",
                  "CXCR4 neutrophil" = "CXCR4", "MX1/ISG neutrophil" = "MX1/ISG")

top_lr <- diff_scores[contrast_id %in% comm_ids]
top_lr[, priority := abs(delta_consensus) * (0.5 + mean_consensus)]
top_lr <- top_lr[order(-priority), head(.SD, 7), by = .(contrast_id, neutrophil_state)]
top_lr[, interaction_route := paste(ligand, "->", receptor)]
top_lr[, contrast_short := factor(contrast_labels[contrast_id], levels = contrast_key$contrast)]
top_lr[, state_short := factor(state_labels[neutrophil_state], levels = unname(state_labels))]
p_a <- ggplot(top_lr, aes(contrast_short, interaction_route,
                          color = delta_consensus, size = mean_consensus)) +
  geom_point(alpha = 0.9) + facet_wrap(~state_short, scales = "free_y", ncol = 2) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  labs(title = "Consensus fibroblast-to-neutrophil interaction changes",
       subtitle = "Agreement across expression, detection, balance, and geometric metrics",
       x = NULL, y = NULL, color = "Consensus change", size = "Mean rank") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.4))

coupling_top <- coupling[contrast_id %in% c("UA5_vs_UF", "UD_vs_UF")]
coupling_top[, priority := route_confidence * abs(directional_coupling)]
coupling_top <- coupling_top[order(-priority), head(.SD, 12),
                             by = .(contrast_id, neutrophil_state)]
coupling_top[, pair := paste(interaction_label, TF, sep = " -> ")]
coupling_top[, contrast_short := factor(contrast_labels[contrast_id], levels = contrast_key$contrast)]
coupling_top[, state_short := factor(state_labels[neutrophil_state], levels = unname(state_labels))]
p_b <- ggplot(coupling_top, aes(contrast_short, pair, fill = directional_coupling)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~state_short, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Blockade-sensitive ligand-receptor-to-TF response coupling",
       subtitle = "Positive values indicate concordant LR and TF activity changes",
       x = NULL, y = NULL, fill = "Coupling") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.0))

# Compact three-layer network for the highest-confidence alpha5beta1/dual routes.
network_routes <- coupling[contrast_id %in% c("UA5_vs_UF", "UD_vs_UF")]
network_routes[, priority := route_confidence * abs(directional_coupling)]
network_routes <- unique(network_routes[order(-priority),
  .SD[!duplicated(paste(fibroblast_subtype, interaction_label, TF, sep = "||"))][1:min(.N, 12)]])
network_routes[, sender_node := sub(" fibroblast$", "", fibroblast_subtype)]
network_routes[, lr_node := paste(interaction_label, state_labels[neutrophil_state], sep = " | ")]
network_routes[, tf_node := paste0(TF, " activity")]

node_names <- unique(c(network_routes$sender_node, network_routes$lr_node, network_routes$tf_node))
node_type <- c(setNames(rep("Fibroblast", uniqueN(network_routes$sender_node)), unique(network_routes$sender_node)),
               setNames(rep("Ligand-receptor", uniqueN(network_routes$lr_node)), unique(network_routes$lr_node)),
               setNames(rep("Neutrophil TF", uniqueN(network_routes$tf_node)), unique(network_routes$tf_node)))
nodes <- data.table(name = node_names, type = unname(node_type[node_names]))
nodes[, x := match(type, c("Fibroblast", "Ligand-receptor", "Neutrophil TF"))]
nodes[, y := seq(0, 1, length.out = .N), by = type]
edges1 <- network_routes[, .(from = sender_node, to = lr_node,
                             coupling = directional_coupling, confidence = route_confidence)]
edges2 <- network_routes[, .(from = lr_node, to = tf_node,
                             coupling = directional_coupling, confidence = route_confidence)]
edges <- unique(rbind(edges1, edges2))
graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
p_c <- ggraph(graph, layout = "manual", x = V(graph)$x, y = V(graph)$y) +
  geom_edge_link(aes(color = coupling, width = confidence), alpha = 0.48,
                 arrow = arrow(length = unit(0.055, "inches")), end_cap = circle(2, "mm")) +
  geom_node_point(aes(shape = type), size = 2.3, color = "#222222") +
  geom_node_text(aes(label = name), size = 2.4, repel = TRUE) +
  scale_edge_color_gradient2(low = "#2166AC", mid = "#BDBDBD", high = "#B2182B", midpoint = 0) +
  scale_edge_width(range = c(0.25, 1.3)) +
  scale_shape_manual(values = c(Fibroblast = 15, `Ligand-receptor` = 16, `Neutrophil TF` = 17)) +
  scale_x_continuous(breaks = 1:3, labels = c("Fibroblast sender", "LR interaction / receiver", "TF response"),
                     expand = expansion(mult = c(0.12, 0.12))) +
  labs(title = "Prioritized alpha5beta1-sensitive signaling-to-regulator model",
       subtitle = "Exploratory consensus routes; fibroblast sampling is limited",
       edge_color = "Directional coupling", edge_width = "Confidence", shape = NULL) +
  theme_void(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, color = "#444444"),
        axis.text.x = element_text(size = 8, face = "bold"),
        legend.position = "bottom", plot.margin = margin(8, 12, 8, 12))

fig <- p_a / p_b / p_c + plot_layout(heights = c(1.0, 1.15, 0.95)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
save_figure(fig, file.path(figure_dir, "Additional_Analysis_3_consensus_interactome_and_TF_coupling"),
            width = 16, height = 20)

saveRDS(list(consensus_scores = scores, differential_scores = diff_scores,
             TF_response_coupling = coupling, prioritized_routes = network_routes),
        file.path(object_dir, "additional_analysis3_consensus_interactome_TF.rds"), compress = TRUE)
message("Additional analysis 3 complete")
