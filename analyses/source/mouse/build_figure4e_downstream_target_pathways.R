suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
})

out_dir <- file.path(
  "FAP_ablation_DSS_epithelial_immune_stromal_figures",
  "multipanel_figure"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lr_families <- data.table(
  direction = c(
    rep("Fibroblast -> Neutrophil", 8),
    rep("Neutrophil -> Fibroblast", 8)
  ),
  receiver = c(rep("Neutrophil", 8), rep("Fibroblast", 8)),
  ligand = c(
    "Il1b", "Il33", "Cxcl1", "Cxcl2", "Cxcl5", "Tgfbi", "Vcam1", "H2.DMb1",
    "Thbs1", "Tgm2", "Nampt", "S100a9", "Osm", "Tnf", "Il1a", "Il10"
  ),
  lr_pathway = c(
    "IL1/IL33", "IL1/IL33", rep("CXCL-CXCR2", 3),
    "TGF-beta/ECM", "VCAM adhesion", "MHC-II/CD74",
    rep("Integrin/ECM", 3), "S100A9-ALCAM", "OSM", "TNF", "IL1", "IL10"
  )
)

models <- data.table(
  model_file = file.path(
    "FAPTK_neutrophil_fibroblast_MultiNicheNet", "objects",
    c(
      "multinichenet_PBS_DSS_vs_PBS.rds",
      "multinichenet_PBS_DSS_vs_GCV_DSS.rds"
    )
  ),
  contrast = c("PBS_DSS-PBS", "PBS_DSS-GCV_DSS")
)

target_parts <- vector("list", nrow(models))
background_parts <- vector("list", nrow(models))
for (i in seq_len(nrow(models))) {
  model <- readRDS(models$model_file[i])
  activity <- as.data.table(
    model$ligand_activities_targets_DEgenes$ligand_activities
  )
  activity <- activity[
    contrast == models$contrast[i] &
      as.character(direction_regulation) == "up" &
      !is.na(ligand_target_weight)
  ]

  background_parts[[i]] <- unique(activity[, .(receiver, target)])
  selected <- merge(
    activity,
    lr_families,
    by = c("receiver", "ligand"),
    all = FALSE
  )
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
    direction == direction_i & lr_pathway == lr_pathway_i,
    target
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
    p_value <- if (
      length(selected_genes) >= 3L && length(pathway_genes) >= 3L && d >= 0L
    ) {
      fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
    } else {
      1
    }
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

# Retain the strongest two receiver programs per ligand-receptor family, then
# show their union within each signaling direction for a compact shared matrix.
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
  enrichment,
  selected_pathways,
  by = c("direction", "downstream_pathway"),
  all = FALSE
)
plot_data <- plot_data[overlap > 0]

direction_order <- c("Fibroblast -> Neutrophil", "Neutrophil -> Fibroblast")
family_order <- c(
  "IL1/IL33", "CXCL-CXCR2", "TGF-beta/ECM", "VCAM adhesion", "MHC-II/CD74",
  "Integrin/ECM", "S100A9-ALCAM", "OSM", "TNF", "IL1", "IL10"
)
plot_data[, direction := factor(direction, levels = direction_order)]
plot_data[, lr_pathway := factor(lr_pathway, levels = family_order)]
plot_data[, downstream_pathway := factor(
  downstream_pathway,
  levels = rev(unique(selected_rows[
    order(factor(direction, levels = direction_order), fdr, -overlap),
    downstream_pathway
  ]))
)]

fwrite(
  plot_data,
  file.path(out_dir, "Figure_4F_downstream_target_pathways.csv")
)
fwrite(
  target_evidence,
  file.path(out_dir, "Figure_4F_ligand_target_evidence.csv")
)

message(
  "Wrote ", nrow(plot_data), " pathway-family cells across ",
  uniqueN(plot_data$downstream_pathway), " downstream pathways."
)
