p <- 'input_data/culture/outputs/019fbe7c_donor_paired_mechanism/nichenet_prior/ligand_target_matrix_nsga2r_final.rds'
x <- readRDS(p)
out <- 'input_data/flow/outputs/15-Figure3-revised'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ligands <- c('COL1A2','FN1','IL1B','CEACAM1','CXCL12','C5')
print(dim(x)); print(ligands %in% colnames(x))
rows <- list()
de <- read.csv('./output/Nature_Extended_Data_Consolidated_2026-09-03/source_data/updated_analyses/full_gene_pseudobulk_differential_expression.csv')
eligible <- unique(de$gene[de$population=='Resolved'])
for(l in ligands){
 if(l %in% colnames(x)){
  z <- x[,l]; z <- z[is.finite(z) & z>0 & names(z) %in% eligible]; z <- sort(z,decreasing=TRUE)
  rows[[l]] <- data.frame(ligand=l,target=names(z)[seq_len(min(3,length(z)))],prior_weight=head(z,3),rank=seq_len(min(3,length(z))))
 }else rows[[l]] <- data.frame(ligand=l,target=NA,prior_weight=NA,rank=NA)
}
write.csv(do.call(rbind,rows),file.path(out,'Figure3g_prior_targets.csv'),row.names=FALSE)
print(do.call(rbind,rows))
