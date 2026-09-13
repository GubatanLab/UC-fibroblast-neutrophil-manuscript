options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table); library(dplyr)})
out <- "UC_fibroblast_neutrophil_analysis/objects"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- readRDS("UCGNE.RDS")
m <- x@meta.data
m$cell <- rownames(m)
m$condition <- dplyr::case_when(
  m$tissue.inf == "Colon Inf UC" ~ "Inflamed UC",
  m$tissue.inf == "Colon Uninf UC" ~ "Uninflamed UC",
  m$tissue.inf == "Colon HD" ~ "Healthy",
  TRUE ~ NA_character_)
m$fib_state <- ifelse(m$Annotations_Level_1 == "Stromal", m$Annotations_Level_3, NA_character_)
m$neut_state <- ifelse(m$Annotations_Level_1 %in% c("Neutrophil", "Immature Neutrophil"), m$Annotations_Level_3, NA_character_)
colon <- !is.na(m$condition)
clinical_cols <- intersect(c("PatientID","condition","MAYO_ES_SCORE","HISTOLOGY_SCORE","FECAL_CAL","CRP",
  "ENDOSCOPIC_SCORE","AGE","SEX","BIOPSY_LOCATION","AntiTNF","Vedolizumab","MEDICATION_AT_BIOPSY","RESPONSE"), names(m))
clinical <- m[colon, clinical_cols, drop = FALSE] %>% distinct(PatientID, condition, .keep_all = TRUE)
data.table::fwrite(clinical, file.path(out, "clinical_by_biopsy.tsv"), sep = "\t", na = "NA")

reds <- names(x@reductions)
cand <- reds[grepl("umap", reds, ignore.case = TRUE)]
red <- cand[which.max(vapply(cand, function(z) nrow(x@reductions[[z]]@cell.embeddings), numeric(1)))]
emb <- x@reductions[[red]]@cell.embeddings
sel <- which(colon & (m$Annotations_Level_1 == "Stromal" | m$Annotations_Level_1 %in% c("Neutrophil","Immature Neutrophil")))
e <- emb[m$cell[sel], 1:2, drop = FALSE]
emb_out <- data.frame(cell = rownames(e), UMAP1 = e[,1], UMAP2 = e[,2],
  PatientID = m$PatientID[sel], condition = m$condition[sel], fib_state = m$fib_state[sel], neut_state = m$neut_state[sel])
saveRDS(emb_out, file.path(out, "fibro_neutrophil_embeddings.rds"), compress = FALSE)

lr_pairs <- data.frame(
  ligand = c("CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CSF3","IL6","CCL2","CCL7",
             "OSM","IL1B","TNF","TGFB1","SPP1","AREG","HBEGF","VEGFA","FN1","COL1A1","ICAM1","VCAM1"),
  receptor = c("CXCR2","CXCR2","CXCR2","CXCR2","CXCR1","CXCR1","CSF3R","IL6R","CCR2","CCR2",
               "OSMR","IL1R1","TNFRSF1A","TGFBR2","CD44","EGFR","EGFR","FLT1","ITGA5","ITGB1","ITGAL","ITGA4"),
  stringsAsFactors = FALSE)
data.table::fwrite(lr_pairs, file.path(out, "curated_lr_pairs.tsv"), sep = "\t")
genes <- intersect(unique(c(lr_pairs$ligand, lr_pairs$receptor, "ITGB1","LIFR","IL6ST","TNFRSF1B","TGFBR1")), rownames(x))
cts <- x[["RNA"]]@counts; dat <- x[["RNA"]]@data
groups <- list(
  Fibroblast_Inflamed = which(m$Annotations_Level_1 == "Stromal" & m$condition == "Inflamed UC"),
  Fibroblast_Uninflamed = which(m$Annotations_Level_1 == "Stromal" & m$condition == "Uninflamed UC"),
  Neutrophil_Inflamed = which(m$Annotations_Level_1 == "Neutrophil" & m$condition == "Inflamed UC"),
  Neutrophil_Uninflamed = which(m$Annotations_Level_1 == "Neutrophil" & m$condition == "Uninflamed UC"))
expr <- bind_rows(lapply(names(groups), function(g) {
  ii <- groups[[g]]
  data.frame(group = g, gene = genes, n_cells = length(ii),
    pct = Matrix::rowMeans(cts[genes, ii, drop = FALSE] > 0),
    avg = Matrix::rowMeans(dat[genes, ii, drop = FALSE]))
}))
data.table::fwrite(expr, file.path(out, "lr_expression_summary.tsv"), sep = "\t")
writeLines(paste("UMAP reduction:", red), file.path(out, "extraction_note.txt"))
