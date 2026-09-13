suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
})

mn_dir <- "a5B1_DSS_neutrophil_fibroblast_MultiNicheNet_excluded_controls"
out_dir <- file.path(
  "a5B1_DSS_epithelial_immune_stromal_remodeling",
  "manuscript_figure5"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Prespecified signaling families represented in Figure 5E. Ligand-target
# evidence is aggregated across the DSS-vs-control and DSS-vs-blockade models.
lr_families <- data.table(
  direction = c(
    rep("Fibroblast -> Neutrophil", 8),
    rep("Neutrophil -> Fibroblast", 8)
  ),
  receiver = c(rep("Neutrophil", 8), rep("Fibroblast", 8)),
  ligand = c(
    "Il1b", "Cxcl1", "Cxcl2", "Cxcl5", "Csf1", "Fn1", "Nampt", "Icam1",
    "Osm", "Osm", "Il1b", "Vegfa", "Nampt", "Thbs1", "S100a9", "Col1a1"
  ),
  lr_pathway = c(
    "IL1", rep("CXCL-CXCR2", 3), "CSF1", "alpha5-integrin/ECM",
    "alpha5-integrin/ECM", "ICAM adhesion",
    "OSM", "OSM", "IL1", "VEGF", "alpha5-integrin/ECM",
    "Integrin/ECM", "S100A9-ALCAM", "Integrin/ECM"
  )
)

models <- file.path(
  mn_dir, "objects",
  c("multinichenet_DSS_vs_Control.rds", "multinichenet_DSS_vs_Blockade.rds")
)

target_parts <- vector("list", length(models))
background_parts <- vector("list", length(models))
for (i in seq_along(models)) {
  model <- readRDS(models[i])
  activity <- as.data.table(
    model$ligand_activities_targets_DEgenes$ligand_activities
  )
  activity <- activity[
    as.character(direction_regulation) == "up" &
      !is.na(ligand_target_weight)
  ]
  background_parts[[i]] <- unique(activity[, .(receiver, target)])
  selected <- merge(activity, lr_families, by = c("receiver", "ligand"), all = FALSE)
  target_parts[[i]] <- selected[, .(
    target_weight = max(ligand_target_weight, na.rm = TRUE)
  ), by = .(direction, receiver, lr_pathway, target)]
  rm(model, activity, selected)
  invisible(gc())
}

target_evidence <- rbindlist(target_parts)[, .(
  target_weight = max(target_weight, na.rm = TRUE)
), by = .(direction, receiver, lr_pathway, target)]
background <- unique(rbindlist(background_parts))

hallmark <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
hallmark[, downstream_pathway := gsub("^HALLMARK_", "", gs_name)]
hallmark[, downstream_pathway := gsub("_", " ", downstream_pathway)]
hallmark[, downstream_pathway := tools::toTitleCase(tolower(downstream_pathway))]
hallmark_sets <- split(hallmark$gene_symbol, hallmark$downstream_pathway)

ora_one <- function(direction_i, receiver_i, lr_pathway_i) {
  selected_genes <- unique(target_evidence[
    direction == direction_i & lr_pathway == lr_pathway_i, target
  ])
  universe <- unique(background[receiver == receiver_i, target])
  selected_genes <- intersect(selected_genes, universe)

  rbindlist(lapply(names(hallmark_sets), function(pathway_i) {
    pathway_genes <- intersect(unique(hallmark_sets[[pathway_i]]), universe)
    overlap_genes <- intersect(selected_genes, pathway_genes)
    a <- length(overlap_genes)
    b <- length(selected_genes) - a
    c <- length(pathway_genes) - a
    d <- length(universe) - a - b - c
    p_value <- if (length(selected_genes) >= 3L && length(pathway_genes) >= 3L && d >= 0L) {
      fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
    } else 1
    data.table(
      direction = direction_i,
      receiver = receiver_i,
      lr_pathway = lr_pathway_i,
      downstream_pathway = pathway_i,
      overlap = a,
      target_genes = length(selected_genes),
      pathway_genes = length(pathway_genes),
      gene_ratio = if (length(selected_genes)) a / length(selected_genes) else 0,
      p_value = p_value,
      leading_targets = paste(head(overlap_genes, 5), collapse = ", ")
    )
  }))
}

families <- unique(lr_families[, .(direction, receiver, lr_pathway)])
enrichment <- rbindlist(lapply(seq_len(nrow(families)), function(i) {
  ora_one(families$direction[i], families$receiver[i], families$lr_pathway[i])
}))
enrichment[, fdr := p.adjust(p_value, method = "BH"), by = .(direction, lr_pathway)]
enrichment[, enrichment_strength := pmin(-log10(pmax(fdr, 1e-6)), 6)]

selected_rows <- enrichment[overlap >= 2][
  order(direction, lr_pathway, fdr, -overlap, -gene_ratio),
  head(.SD, 2),
  by = .(direction, lr_pathway)
]
selected_pathways <- selected_rows[, .(
  best_fdr = min(fdr, na.rm = TRUE),
  total_overlap = sum(overlap, na.rm = TRUE)
), by = .(direction, downstream_pathway)][
  order(direction, best_fdr, -total_overlap),
  head(.SD, 6),
  by = direction
][, .(direction, downstream_pathway)]
plot_data <- merge(
  enrichment, selected_pathways,
  by = c("direction", "downstream_pathway"), all = FALSE
)[overlap > 0]

fwrite(plot_data, file.path(out_dir, "Figure_5F_downstream_target_pathways.csv"))
fwrite(target_evidence, file.path(out_dir, "Figure_5F_ligand_target_evidence.csv"))
message("Wrote ", nrow(plot_data), " Figure 5F pathway-family cells.")
