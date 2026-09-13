suppressPackageStartupMessages(library(data.table))

out_dir <- "a5B1_DSS_neutrophil_fibroblast_MultiNicheNet_excluded_controls"
table_dir <- file.path(out_dir, "tables")
all_results <- fread(file.path(table_dir, "both_comparisons_nonambient_interactions.csv"))
key <- c("sender", "receiver", "direction", "ligand", "receptor", "interaction")

disease <- all_results[comparison == "DSS_vs_Control", c(
  key, "score_DSS", "score_comparator", "fraction_DSS", "fraction_comparator", "score_delta"
), with = FALSE]
setnames(disease, c("score_DSS", "score_comparator", "fraction_DSS", "fraction_comparator", "score_delta"),
         c("score_DSS_disease_model", "score_Control", "fraction_DSS_disease_model", "fraction_Control", "delta_DSS_vs_Control"))

blockade <- all_results[comparison == "DSS_vs_Blockade", c(
  key, "score_DSS", "score_comparator", "fraction_DSS", "fraction_comparator", "score_delta"
), with = FALSE]
setnames(blockade, c("score_DSS", "score_comparator", "fraction_DSS", "fraction_comparator", "score_delta"),
         c("score_DSS_blockade_model", "score_Blockade", "fraction_DSS_blockade_model", "fraction_Blockade", "delta_DSS_vs_Blockade"))

paired <- merge(disease, blockade, by = key)
paired[, `:=`(
  concordance_score = pmin(abs(delta_DSS_vs_Control), abs(delta_DSS_vs_Blockade)),
  combined_delta = delta_DSS_vs_Control + delta_DSS_vs_Blockade
)]

attenuated <- paired[delta_DSS_vs_Control > 0 & delta_DSS_vs_Blockade > 0]
setorder(attenuated, -concordance_score, -combined_delta)
attenuated[, rank := .I]
fwrite(attenuated, file.path(table_dir, "DSS_induced_blockade_attenuated_ranked.csv"))

restored <- paired[delta_DSS_vs_Control < 0 & delta_DSS_vs_Blockade < 0]
setorder(restored, -concordance_score, combined_delta)
restored[, rank := .I]
fwrite(restored, file.path(table_dir, "DSS_suppressed_blockade_restored_ranked.csv"))

discordant <- paired[sign(delta_DSS_vs_Control) != sign(delta_DSS_vs_Blockade)]
setorder(discordant, -concordance_score)
discordant[, rank := .I]
fwrite(discordant, file.path(table_dir, "discordant_condition_effects_ranked.csv"))

key_genes <- c(
  "Osm", "Osmr", "Il6st", "Cxcl1", "Cxcl2", "Cxcl5", "Cxcr2", "Nampt",
  "Itga5", "Itgb1", "Il1b", "Il1r1", "Il1rap", "Il1r2", "S100a8", "S100a9",
  "Tlr4", "Alcam", "Vegfa", "Nrp1", "Thbs1", "Itga4", "Col1a1", "Col1a2",
  "Col3a1", "Fn1", "Icam1", "Pecam1", "Tgfb1", "Csf1", "Csf1r"
)
key_tbl <- all_results[ligand %chin% key_genes | receptor %chin% key_genes]
setorder(key_tbl, comparison, direction, -abs_delta)
fwrite(key_tbl, file.path(table_dir, "biologically_prioritized_pathway_interactions.csv"))

summary_tbl <- paired[, .(
  tested_in_both = .N,
  DSS_induced_blockade_attenuated = sum(delta_DSS_vs_Control > 0 & delta_DSS_vs_Blockade > 0),
  DSS_suppressed_blockade_restored = sum(delta_DSS_vs_Control < 0 & delta_DSS_vs_Blockade < 0),
  discordant = sum(sign(delta_DSS_vs_Control) != sign(delta_DSS_vs_Blockade))
), by = direction]
fwrite(summary_tbl, file.path(table_dir, "paired_effect_direction_summary.csv"))

message("Wrote cross-contrast interpretation tables")
