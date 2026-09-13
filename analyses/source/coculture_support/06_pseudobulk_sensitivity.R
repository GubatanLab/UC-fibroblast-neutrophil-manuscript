source('output/additional_analyses_2026-09-02/code/common.R')
suppressPackageStartupMessages({library(edgeR);library(limma)})
bun<-readRDS(file.path(OUT,'objects/full_gene_pseudobulk.rds'))
pr<-readRDS(file.path(OUT,'objects/rna_analysis_provenance.rds'))$programs
res<-list()
for(pop in names(bun))for(cn in c('alpha_interaction','alpha_in_UC')){
 b<-bun[[pop]];w<-contrasts_list[[cn]];sm<-b$meta[ConditionCode%in%names(w)&cells>=20]
 donors<-sm[,.(n=.N),by=DonorID][n==length(w),DonorID];sm<-sm[DonorID%in%donors];sm[,ConditionCode:=factor(ConditionCode,levels=names(w))]
 design<-model.matrix(~0+ConditionCode+factor(DonorID),sm);cv<-c(unname(w),rep(0,ncol(design)-length(w)))
 y<-DGEList(as.matrix(b$counts[,sm$sample,drop=FALSE]));keep<-filterByExpr(y,design,min.count=3);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);vo<-voom(y,design,plot=FALSE);inds<-ids2indices(pr,rownames(y),remove.empty=TRUE)
 for(cor in c(.01,NA_real_)){
  ca<-as.data.table(camera(vo,inds,design,contrast=cv,inter.gene.cor=cor),keep.rownames='program');ca[,`:=`(population=pop,contrast=cn,n_donors=length(donors),correlation_method=if(is.na(cor))'Estimated within gene set' else 'Fixed 0.01')];res[[paste(pop,cn,cor)]]<-ca
 }
}
write_tab(rbindlist(res,fill=TRUE),'pseudobulk_program_correlation_sensitivity')
# Eligibility audit of the existing external spatial files; no downloads or new data.
spbase<-'input_data/mouse/Figure 6 Spatial External Validation/data'
rocha<-file.path(spbase,'Rocha2024/longitudinal_exprMat.csv')
meta<-file.path(spbase,'Rocha2024/longitudinal_metadata.csv')
rocha_genes<-names(fread(rocha,nrows=0))
scfiles<-list.files(file.path(spbase,'SCP3818'),pattern='^UC1_inflamed_.*\\.tsv$',full.names=FALSE)
scgenes<-sub('^UC1_inflamed_|\\.tsv$','',scfiles);scgenes<-gsub('\\.tsv$','',scgenes)
coverage<-rbindlist(lapply(names(pr),function(nm)rbindlist(list(data.table(dataset='Existing Rocha/Mennillo CosMx',program=nm,defined=length(pr[[nm]]),available=sum(pr[[nm]]%in%rocha_genes),genes=paste(intersect(pr[[nm]],rocha_genes),collapse=';')),data.table(dataset='Existing SCP3818 one-donor exports',program=nm,defined=length(pr[[nm]]),available=sum(pr[[nm]]%in%scgenes),genes=paste(intersect(pr[[nm]],scgenes),collapse=';'))))))
coverage[,fraction:=available/defined];write_tab(coverage,'external_program_measurement_coverage')
mt<-fread(meta,select=c('Fine_annotation_3','HS','Condition'))
write_tab(mt[,.(cells=.N),by=.(Fine_annotation_3,HS,Condition)],'external_spatial_annotation_inventory')
message('Pseudobulk sensitivity and external eligibility audit complete')
