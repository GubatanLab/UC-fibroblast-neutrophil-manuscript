options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Matrix); library(SeuratObject)})
root <- normalizePath(getwd(), winslash = "/")
out <- file.path(root, "output/pdf/Figure_1_Additional_Panels_2026-09-12/source_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
src <- "input_data/human"
score_path <- file.path(src, "UC_fibroblast_neutrophil_analysis/objects/neutrophil_cell_scores.rds")
scores <- readRDS(score_path)
write.csv(scores, file.path(out, "neutrophil_original_cell_scores.csv"), row.names = FALSE)
exprs <- parse(file.path(src, "run_uc_fibro_neutrophil_analysis.R"))
for (ex in exprs) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.symbol(ex[[2]]) &&
      as.character(ex[[2]]) %in% c("gene_sets_fib", "gene_sets_neut")) eval(ex)
}
defs <- rbind(do.call(rbind, lapply(names(gene_sets_fib), function(k) data.frame(compartment="Fibroblast", program=k, gene=gene_sets_fib[[k]]))),
              do.call(rbind, lapply(names(gene_sets_neut), function(k) data.frame(compartment="Neutrophil", program=k, gene=gene_sets_neut[[k]]))))
write.csv(defs, file.path(out, "original_program_gene_definitions.csv"), row.names=FALSE)
raw_path <- "input_data/mouse/UCGNE_fibroblast_neutrophil_subclustering/objects/UCGNE_colon_neutrophils_subclustered.rds"
message("Loading existing neutrophil expression object")
obj <- readRDS(raw_path)
message("Object loaded: ", nrow(obj), " genes, ", ncol(obj), " cells")
assay <- obj[["RNA"]]
dat <- if (inherits(assay, "Assay5")) LayerData(assay, layer="data") else slot(assay, "data")
overlap <- intersect(scores$cell, colnames(dat))
message("Original-score cells matched: ",length(overlap)," / ",nrow(scores))
missing <- scores[!scores$cell %in% overlap,]
write.csv(missing,file.path(out,"cells_absent_from_neutrophil_expression_object.csv"),row.names=FALSE)
stopifnot(all(missing$neut_state %in% c("Erythrocyte","Monocyte-derived macrophage")))
scores <- scores[scores$cell %in% overlap,]
dat <- dat[,scores$cell,drop=FALSE]
cover <- do.call(rbind,lapply(names(gene_sets_neut),function(k) {
  gg <- gene_sets_neut[[k]]
  data.frame(program=k,gene=gg,present=gg %in% rownames(dat),state_marker=gg %in% c("OSM","CXCR4","PADI4","MX1"))
}))
write.csv(cover,file.path(out,"neutrophil_program_gene_coverage.csv"),row.names=FALSE)
errors <- sapply(names(gene_sets_neut),function(k) {
  gg <- intersect(gene_sets_neut[[k]],rownames(dat))
  max(abs(Matrix::colMeans(dat[gg,,drop=FALSE])-scores[[k]]))
})
write.csv(data.frame(program=names(errors),max_abs_difference=errors),file.path(out,"expression_reproduction_audit.csv"),row.names=FALSE)
print(errors)
stopifnot(all(errors < 1e-7))
clean <- scores[,c("PatientID","condition","biopsy_id","neut_state","cell")]
for (k in names(gene_sets_neut)) {
  gg <- setdiff(intersect(gene_sets_neut[[k]],rownames(dat)),c("OSM","CXCR4","PADI4","MX1"))
  clean[[k]] <- Matrix::colMeans(dat[gg,,drop=FALSE])
}
write.csv(clean,file.path(out,"neutrophil_marker_excluded_cell_scores.csv"),row.names=FALSE)
genes <- unique(c("OSM","CXCR4","PADI4","MX1","IL1B","CXCL8","NFKBIA","CXCR2","FPR1","ITGAM","ICAM1","CD44","BCL2A1","MMP9","ELANE","LTF","CYBB","NCF2","RAC2","ISG15","IFIT3","IFI6"))
genes <- intersect(genes,rownames(dat))
cellgenes <- cbind(scores[,c("PatientID","condition","biopsy_id","neut_state","cell")],as.data.frame(t(as.matrix(dat[genes,,drop=FALSE]))))
write.csv(cellgenes,file.path(out,"neutrophil_selected_gene_expression.csv"),row.names=FALSE)
writeLines(c(paste("Original scores:",score_path),paste("Expression:",raw_path),paste("Matched cells:",length(overlap)),paste("Maximum score reproduction error:",max(errors))),file.path(out,"extraction_provenance.txt"))
message("Extraction and marker-exclusion sensitivity complete")
