suppressPackageStartupMessages({library(Seurat);library(Matrix);library(edgeR)})
root <- '.'
out <- file.path(root,'output/Priority_Figure_Additions_2026-09-07')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source <- 'input_data/mouse/MouseColitis_Inhibitors_compartment_subclustering/reassigned_fractions_epithelial_audit_cleaned/updated_umaps/level3_reannotated/MouseColitis_Inhibitors_Epithelial_Level3Reannotated.rds'
cat('Loading epithelial object\n');flush.console()
obj <- readRDS(source)
md <- obj[[]];md$cell <- rownames(md)
cat('Loaded',nrow(md),'cells; assays:',names(obj@assays),'\n');print(colnames(md));flush.console()
units <- read.csv(file.path(root,'config/biological_mouse_units.csv'))
md <- md[md$sample_uid %in% units$sample_uid,,drop=FALSE]
write.csv(md,file.path(out,'epithelial_cell_metadata.csv'),row.names=FALSE)
write.csv(units,file.path(out,'mouse_units.csv'),row.names=FALSE)
assay <- if('RNA' %in% names(obj@assays)) 'RNA' else stop('No RNA assay')
layers <- Layers(obj[[assay]],search='^counts')
cat('Count layers:',layers,'\n');flush.console()
genes <- rownames(obj[[assay]])
keys <- rbind(data.frame(cell=md$cell,sample_uid=md$sample_uid,scope='All epithelial',state='All epithelial'),data.frame(cell=md$cell,sample_uid=md$sample_uid,scope='Level 2',state=md$reassigned_level_2),data.frame(cell=md$cell,sample_uid=md$sample_uid,scope='Level 3',state=md$reassigned_level_3))
keys$key <- paste(keys$scope,keys$state,keys$sample_uid,sep='___')
info <- unique(keys[,c('key','scope','state','sample_uid')]);info$n_cells <- as.integer(table(factor(keys$key,levels=info$key)))
pb <- matrix(0,nrow=length(genes),ncol=nrow(info),dimnames=list(genes,info$key))
seen <- character()
for(lay in layers){
 cat('Aggregating',lay,'\n');flush.console()
 mat <- LayerData(obj,assay=assay,layer=lay)
 cells <- intersect(colnames(mat),md$cell)
 stopifnot(length(intersect(cells,seen))==0);seen <- c(seen,cells)
 if(!length(cells))next
 k <- keys[keys$cell %in% cells,]
 design <- sparseMatrix(i=match(k$cell,cells),j=match(k$key,info$key),x=1,dims=c(length(cells),nrow(info)))
 sums <- mat[,cells,drop=FALSE] %*% design
 pb[match(rownames(mat),genes),] <- pb[match(rownames(mat),genes),]+as.matrix(sums)
 rm(mat,sums);gc()
}
stopifnot(setequal(seen,md$cell))
write.csv(info,file.path(out,'epithelial_pseudobulk_inventory.csv'),row.names=FALSE)
con <- gzfile(file.path(out,'epithelial_pseudobulk_counts.csv.gz'),'wt');write.csv(pb,con);close(con)
norm <- matrix(NA_real_,nrow=nrow(pb),ncol=ncol(pb),dimnames=dimnames(pb))
for(scope_state in unique(paste(info$scope,info$state,sep='___'))){
 ix <- which(paste(info$scope,info$state,sep='___')==scope_state & info$n_cells>=20 & colSums(pb)>0)
 if(length(ix)<2)next
 d <- calcNormFactors(DGEList(counts=pb[,ix,drop=FALSE]),method='TMM')
 norm[,ix] <- cpm(d,log=TRUE,prior.count=1)
}
con <- gzfile(file.path(out,'epithelial_logcpm.csv.gz'),'wt');write.csv(norm,con);close(con)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
cat('DONE epithelial extraction\n')
