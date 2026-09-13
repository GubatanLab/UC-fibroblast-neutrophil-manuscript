suppressPackageStartupMessages(library(data.table))
d <- "giotto_codex_results/tables"

samples <- fread(file.path(d, "sample_summary_all.csv"))
net <- fread(file.path(d, "giotto_spatial_network_summary.csv"))
comp <- fread(file.path(d, "statistics_celltype_composition_pairwise_wilcoxon.csv"))
cpe <- fread(file.path(d, "statistics_cell_proximity_pairwise_wilcoxon.csv"))
cpe_group <- fread(file.path(d, "giotto_cell_proximity_enrichment_group_summary.csv"))

cat("Analyzed samples by group:\n")
print(samples[analysis_included == TRUE, .(patients = .N, cells = sum(n_cells)), by = Diagnosis2])
cat("\nNetwork totals/range:\n")
print(net[, .(
  patients = .N,
  cells = sum(n_cells),
  edges = sum(n_undirected_edges),
  mean_degree_min = min(mean_edges_per_cell),
  mean_degree_median = median(mean_edges_per_cell),
  mean_degree_max = max(mean_edges_per_cell)
)])

cat("\nComposition tests BH < 0.05:\n")
print(comp[p_adj_BH < 0.05][order(p_adj_BH, -abs(median_difference))])

cat("\nCPE tests BH < 0.05 (top 30):\n")
print(head(cpe[p_adj_BH < 0.05][order(p_adj_BH, -abs(median_difference))], 30))
cat("Significant CPE tests by contrast:\n")
print(cpe[p_adj_BH < 0.05, .N, by = contrast])

cat("\nStrongest group-average enrichments/depletions:\n")
for (grp in unique(cpe_group$Diagnosis2)) {
  cat("\n", grp, "enriched:\n", sep = "")
  print(head(cpe_group[Diagnosis2 == grp][order(-mean_log2_enrichment)], 8))
  cat(grp, "depleted:\n")
  print(head(cpe_group[Diagnosis2 == grp][order(mean_log2_enrichment)], 8))
}
