suppressPackageStartupMessages({library(Seurat); library(data.table)})
obj <- readRDS("UCCODEX1_Annotated_repaired_withUMAP_neutrophil0.8Matched.rds")
out <- file.path("giotto_codex_results", "additional_analyses")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
coords <- rbindlist(lapply(Images(obj), function(im) {
  z <- as.data.table(GetTissueCoordinates(obj, image=im)); setnames(z,"cell","cell_ID"); z[,region:=im]; z
}), fill=TRUE)
metadata <- obj[[]]
cell_ids <- rownames(metadata)
m <- as.data.table(metadata, keep.rownames="cell_ID")
a5_counts <- as.numeric(LayerData(obj, assay="Akoya", layer="counts")["Cell..a5B1.mean", cell_ids])
a5_clr <- as.numeric(LayerData(obj, assay="Akoya", layer="data")["Cell..a5B1.mean", cell_ids])
fap_counts <- as.numeric(LayerData(obj, assay="Akoya", layer="counts")["Cell..FAPa.mean", cell_ids])
fap_clr <- as.numeric(LayerData(obj, assay="Akoya", layer="data")["Cell..FAPa.mean", cell_ids])
m[, `:=`(a5B1_cell_raw=a5_counts, a5B1_cell_clr=a5_clr,
         FAPa_cell_raw=fap_counts, FAPa_cell_clr=fap_clr)]
base <- c("cell_ID","PatientID","Diagnosis2","Medication","Celltype_updated",
          "Neutrophil_subtype_0.8","Cell..Area","Nucleus..Area",
          "a5B1_cell_raw","a5B1_cell_clr","FAPa_cell_raw","FAPa_cell_clr")
markers <- c("PADI4","CXCR4","MX.1","OSM","CD66b","CD15","CD16","FAPa","CD140a",
             "aSMA","Vimentin","Collagen.4","Podoplanin","CD34","CD68","CD163","CD14",
             "HLA.DR","CD11b","CD4","CD8","FoxP3","CD25","CD279","TIGIT","LAG.3",
             "EpCAM","CK7","P53","Ki67","HLA.ABC")
marker_cols <- paste0("Nucleus..", markers, ".mean")
keep <- intersect(c(base, marker_cols), names(m))
m <- m[, ..keep]
setnames(m, "Celltype_updated", "cell_type")
setnames(m, intersect(marker_cols,names(m)), sub("^Nucleus\\.\\.","",sub("\\.mean$","",intersect(marker_cols,names(m)))))
m <- merge(m, coords[,.(cell_ID,x,y,region)], by="cell_ID", all.x=TRUE)
stopifnot(!anyNA(m$x), !anyNA(m$y), !anyDuplicated(m$cell_ID))
fwrite(m, file.path(out,"codex_extended_cells.csv"))
cat("Exported",nrow(m),"cells and",ncol(m),"columns\n")
