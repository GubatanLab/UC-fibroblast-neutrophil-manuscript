p <- 'input_data/culture/outputs/019fbe7c_donor_paired_mechanism/nichenet_prior/ligand_target_matrix_nsga2r_final.rds'
x <- readRDS(p)
root <- 'input_data/flow'
out <- file.path(root,'outputs/17-Figure3g-module-preview');dir.create(out,recursive=TRUE,showWarnings=FALSE)
de <- read.csv('./output/Nature_Extended_Data_Consolidated_2026-09-03/source_data/updated_analyses/full_gene_pseudobulk_differential_expression.csv')
genes <- intersect(rownames(x),unique(de$gene[de$population=='Resolved']))
defs <- read.csv(file.path(root,'outputs/12-RNA-module-heatmap/RNA_module_gene_membership.csv'))
ligands <- c('COL1A2','FN1','IL1B','CEACAM1','CXCL12','C5')
mods <- c('OSM_inflammation','Type_I_interferon','Degranulation','Oxidative_burst','Chemotaxis','Hypoxia_glycolysis')
rows <- list(); details <- list(); n <- length(genes)
for(l in ligands){
 z <- x[genes,l];stopifnot(all(is.finite(z)))
 percentile <- (rank(z,ties.method='average')-.5)/n
 for(m in mods){
  requested <- unique(defs$gene[defs$program==m & tolower(defs$in_full)=="true"]);present <- intersect(requested,genes)
  stopifnot(length(present)>0)
  ix <- match(present,genes);r <- percentile[ix];weights <- z[ix]
  ranked <- order(-weights,present)
  rows[[paste(l,m)]] <- data.frame(ligand=l,module=m,score=2*mean(r)-1,n_module_genes=length(present),n_requested=length(requested),top_decile_targets=sum(r>=.90),examples=paste(present[head(ranked,2)],collapse=';'),background_genes=n)
  details[[paste(l,m)]] <- data.frame(ligand=l,module=m,gene=present,prior_weight=weights,percentile=r,top_decile=r>=.90)
 }
}
write.csv(do.call(rbind,rows),file.path(out,'Ligand_module_scores.csv'),row.names=FALSE)
write.csv(do.call(rbind,details),file.path(out,'Ligand_module_gene_support.csv'),row.names=FALSE)
writeLines(c(p,paste('Background genes:',n),'Score: 2 x mean full-background percentile rank - 1; average ranks for ties.','Top decile: average percentile rank >=0.90. Descriptive prior support; no statistical significance test.'),file.path(out,'Analysis_provenance.txt'))
print(do.call(rbind,rows)[,c('ligand','module','score','top_decile_targets')])

