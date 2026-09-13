suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

input_file <- "UCCODEX1_Annotated_repaired_withUMAP_neutrophil0.8Matched.rds"
output_dir <- file.path("giotto_codex_results", "nolan_neighborhoods")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_file)
coords <- rbindlist(lapply(Images(obj), function(image_name) {
  z <- as.data.table(GetTissueCoordinates(obj, image = image_name))
  setnames(z, "cell", "cell_ID")
  z[, region := image_name]
  z
}), use.names = TRUE, fill = TRUE)

meta <- as.data.table(obj[[]], keep.rownames = "cell_ID")
keep <- c("cell_ID", "PatientID", "Diagnosis2", "Celltype_updated")
stopifnot(all(keep %in% names(meta)))
meta <- merge(meta[, ..keep], coords[, .(cell_ID, x, y, region)], by = "cell_ID", all.x = TRUE)
if (anyNA(meta$x) || anyNA(meta$y) || anyDuplicated(meta$cell_ID)) stop("Invalid spatial coordinate mapping")
meta[, cell_type := as.character(Celltype_updated)]
meta[is.na(cell_type) | cell_type == "", cell_type := "Unclassified"]
meta[, Celltype_updated := NULL]
setorder(meta, PatientID, region, cell_ID)
fwrite(meta, file.path(output_dir, "nolan_input_cells.csv"))
cat("Exported", nrow(meta), "cells from", uniqueN(meta$PatientID), "patients\n")
