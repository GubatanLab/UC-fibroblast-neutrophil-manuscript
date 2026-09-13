source('output/additional_analyses_2026-09-02/code/common.R')
suppressPackageStartupMessages({library(SeuratObject);library(edgeR);library(limma)})
message('Reproduce composition, then evaluate independent donor differences')
comp<-fread(file.path(CO,'Manuscript_Analyses_1_to_6/source_data/analysis2_canonical_state_composition_by_sample.csv'))
comp<-comp[ConditionCode!='PRE'];comp[,cells:=total_cells]
stopifnot(max(abs(comp[,sum(proportion),by=.(DonorID,ConditionCode)]$V1-1))<1e-10)
all_comp<-list();all_ind<-list()
for(th in c(0,20,50)){
 z<-calc_contrast(comp,keys='state',value='clr',mincells=th);all_comp[[as.character(th)]]<-z$stats;all_ind[[as.character(th)]]<-z$individual
}
cs<-rbindlist(all_comp);ci<-rbindlist(all_ind)
cs[,q:=p.adjust(p,'BH'),by=.(contrast,min_cells)];cs[,q_signflip:=p.adjust(p_signflip,'BH'),by=.(contrast,min_cells)]
write_tab(cs,'composition_donor_contrasts');write_tab(ci,'composition_individual_differences');write_tab(unique(comp[,.(DonorID,ConditionCode,cells)]),'canonical_cell_counts')
lod<-rbindlist(lapply(unique(ci$DonorID),function(d){x<-ci[DonorID!=d & min_cells==0 & contrast=='alpha_interaction'];y<-x[,mean_test(delta),by=state];y[,omitted_donor:=d];y}))
write_tab(lod,'composition_leave_one_donor_out')
orig<-fread(file.path(CO,'Manuscript_Analyses_1_to_6/Main_And_Supplementary_Figures/Main_Figure_3/source_data_revised_A_to_I/Figure_3F_factorial_state_abundance.csv'))
audit<-merge(orig[contrast_id=='alpha5beta1_fibroblast_interaction',.(state,original_effect=effect,original_q=FDR)],cs[contrast=='alpha_interaction' & min_cells==0],by='state')
stopifnot(max(abs(audit$original_effect-audit$effect))<1e-10)
write_tab(audit,'original_vs_donor_variance_comparison')
message('Loading full-gene original coculture object')
rawfile<-'input_data/culture/outputs/019fbe7c_combined_decontX_reclustered_all_cells/Combined_neutrophil_fibroblast_decontX_reclustered.rds'
raw<-readRDS(rawfile)
print(dim(raw));print(Assays(raw));print(Layers(raw[['RNA']]))
md<-as.data.table(raw@meta.data,keep.rownames='barcode');cts<-LayerData(raw[['RNA']],layer='counts')
stopifnot(all(colnames(cts)==md$barcode),all(cts@x>=0),all(abs(cts@x-round(cts@x))<1e-8))
canmd<-fread(file.path(ND,'source_data/canonical_neutrophil_cell_metadata_and_umap.csv.gz'))[ConditionCode!='PRE']
md[,state:=canmd$neutrophil_state_annotation[match(barcode,canmd$cell)]]
md[,resolved:=!is.na(state)]
lineage<-fread(file.path(OUT,'tables/lineage_eligibility_before_unresolved_removal.csv'))
md[,lineage_eligible:=clean_cell_class=='Neutrophil' & barcode%in%lineage[lineage_eligible==TRUE,barcode]]
stopifnot(all(md[resolved==TRUE,lineage_eligible]))
stopifnot(sum(md$resolved)==nrow(canmd));stopifnot(all(md[resolved==TRUE,DonorID]==canmd$DonorID[match(md[resolved==TRUE,barcode],canmd$cell)]))
md[,full_UMI:=as.numeric(Matrix::colSums(cts))];md[,full_genes:=as.numeric(Matrix::colSums(cts>0))]
qca<-md[,.(cells=.N,resolved=sum(resolved),median_UMI=median(full_UMI),median_genes=median(full_genes)),by=.(DonorID,ConditionCode,clean_cell_class,clean_state)]
write_tab(qca,'full_gene_cell_quality_by_sample_state')
write_tab(md[,.(barcode,DonorID,SampleID,ConditionCode,clean_cell_class,clean_state,state,resolved,lineage_eligible,full_UMI,full_genes)],'rna_cell_annotation_and_quality')
neut<-md[lineage_eligible==TRUE];qc<-neut[,.(all_neutrophils=.N,resolved_cells=sum(resolved),unresolved_fraction=mean(!resolved),median_UMI=median(full_UMI),median_genes=median(full_genes)),by=.(DonorID,ConditionCode)]
write_tab(qc,'all_neutrophil_sample_quality')
allcomp<-merge(comp,qc[,.(DonorID,ConditionCode,all_neutrophils)],by=c('DonorID','ConditionCode'))
allcomp[,fraction_all:=N/all_neutrophils]
den<-calc_contrast(allcomp,'state',value='fraction_all',mincells=0);den$stats[,q:=p.adjust(p,'BH'),by=contrast]
write_tab(allcomp,'composition_both_denominators');write_tab(den$stats,'all_neutrophil_denominator_contrasts')
programs<-list(
 NETosis=c('PADI4','MPO','ELANE','PRTN3','CTSG','AZU1','CYBB','NCF1','NCF2','S100A8','S100A9'),
 Degranulation=c('LTF','BPI','CAMP','LCN2','MMP8','MMP9','OLFM4','CEACAM8','FCGR3B'),
 Oxidative_burst=c('CYBB','CYBA','NCF1','NCF2','NCF4','RAC2','MPO'),
 Phagocytosis=c('FCGR3B','FCGR2A','ITGAM','ITGB2','FPR1','FPR2','C5AR1','TREM1','SYK','LYN'),
 Chemotaxis=c('CXCR1','CXCR2','CXCR4','FPR1','FPR2','C5AR1','CCR1','CCRL2'),
 OSM_inflammation=c('OSM','IL1B','TNF','NFKBIA','CXCL8','CCL3','CCL4','IL1R2','SOCS3'),
 Type_I_interferon=c('MX1','ISG15','IFIT1','IFIT2','IFIT3','IFI6','OAS1','OAS2','OASL','RSAD2','IRF7','STAT1'),
 Immature_granulopoiesis=c('LTF','BPI','CAMP','LCN2','OLFM4','CEACAM8','DEFA4','MS4A3'),
 Survival_antiapoptosis=c('BCL2A1','MCL1','BCL2L1','CFLAR','FAS','BAX','BBC3','CASP3'),
 Hypoxia_glycolysis=c('HIF1A','SLC2A1','HK2','PFKP','ALDOA','GAPDH','PGK1','LDHA','ENO1'),
 ER_stress=c('XBP1','ATF4','DDIT3','HSPA5','HERPUD1','DNAJB9','ERO1A'),
 Proliferation_stress=c('MKI67','PCNA','TYMS','DHFR','TUBA1B','STMN1'),
 Integrin_adhesion=c('ITGA5','ITGB1','ITGAM','ITGB2','SELPLG','ICAM1'),
 CXCR4_aging_retention=c('CXCR4','SELL','FCGR3B','PECAM1','C5AR1','BCL2A1'))
pg<-rbindlist(lapply(names(programs),function(p)data.table(program=p,gene=programs[[p]])))
reduced<-readRDS(file.path(OUT,'tables/input_audit.rds'))$genes
pg[,`:=`(in_full=gene%in%rownames(cts),in_reduced=gene%in%reduced,identity_marker=gene%in%c('OSM','PADI4','CXCR4','MX1'))]
write_tab(pg,'program_gene_coverage')
genes<-intersect(unique(pg$gene),rownames(cts));use<-which(md$lineage_eligible & md$full_UMI>0)
e<-as.matrix(cts[genes,use,drop=FALSE]);e<-log1p(sweep(e,2,md$full_UMI[use],'/')*1e4)
nm<-copy(md[use]);nm[,sample:=paste(DonorID,ConditionCode,sep='|')]
# Equal-weight untreated donor-condition groups define a fixed reference scale.
ref<-which(nm$resolved & nm$ConditionCode%in%c('N','CF','UF'))
groups<-split(ref,nm$sample[ref]);mu<-rowMeans(sapply(groups,function(ii)rowMeans(e[,ii,drop=FALSE])))
second<-rowMeans(sapply(groups,function(ii)rowMeans(e[,ii,drop=FALSE]^2)))
sdref<-sqrt(pmax(second-mu^2,0));sdref[sdref<1e-8]<-1
zz<-sweep(sweep(e,1,mu,'-'),1,sdref,'/')
write_tab(data.table(gene=genes,reference_mean=mu,reference_sd=sdref),'program_reference_scaling')
scorelists<-list();pops<-list();sts<-list()
for(ver in c('full_gene','identity_markers_removed')){
 avail<-lapply(programs,intersect,y=genes)
 if(ver=='identity_markers_removed')avail<-lapply(avail,setdiff,y=c('OSM','PADI4','CXCR4','MX1'))
 avail<-avail[lengths(avail)>=3]
 sc<-sapply(avail,function(g)colMeans(zz[g,,drop=FALSE]));colnames(sc)<-names(avail)
 sd<-cbind(nm[,.(barcode,DonorID,ConditionCode,state,resolved,full_UMI,full_genes)],as.data.table(sc))
 lo<-melt(sd,id.vars=c('barcode','DonorID','ConditionCode','state','resolved','full_UMI','full_genes'),variable.name='program',value.name='value')
 st<-lo[resolved==TRUE,.(value=mean(value),cells=.N),by=.(DonorID,ConditionCode,state,program)];st[,version:=ver];sts[[ver]]<-st
 pp<-rbindlist(list(lo[resolved==TRUE,.(value=mean(value),cells=.N,population='Resolved'),by=.(DonorID,ConditionCode,program)],lo[,.(value=mean(value),cells=.N,population='All annotated'),by=.(DonorID,ConditionCode,program)],lo[resolved & full_genes>=100,.(value=mean(value),cells=.N,population='Resolved; >=100 genes'),by=.(DonorID,ConditionCode,program)]))
 pp[,version:=ver];pops[[ver]]<-pp
}
ps<-rbindlist(pops);ss<-rbindlist(sts)
write_tab(ps,'population_program_scores');write_tab(ss,'state_program_scores')
stats<-list();inds<-list()
for(th in c(0,20,50)){
 z<-calc_contrast(ps,c('program','population','version'),mincells=th);stats[[as.character(th)]]<-z$stats;inds[[as.character(th)]]<-z$individual
}
pt<-rbindlist(stats);pi<-rbindlist(inds)
pt[,q:=p.adjust(p,'BH'),by=.(contrast,min_cells,population,version)];pt[,q_signflip:=p.adjust(p_signflip,'BH'),by=.(contrast,min_cells,population,version)]
write_tab(pt,'population_program_contrasts');write_tab(pi,'population_program_donor_differences')
zs<-calc_contrast(ss,c('state','program','version'),mincells=10);zs$stats[,q:=p.adjust(p,'BH'),by=.(contrast,version)]
write_tab(zs$stats,'within_state_program_contrasts')
# Symmetric composition/activity decomposition on states with adequate support in both arms.
dec<-list();k<-0
for(ver in unique(ss$version))for(d in unique(ss$DonorID))for(th in c(1,10)){
 l<-ss[DonorID==d & version==ver & ConditionCode=='UF'];r<-ss[DonorID==d & version==ver & ConditionCode=='UA5']
 a<-merge(l,r,by=c('DonorID','state','program','version'),suffixes=c('_0','_1'))
 a<-a[cells_0>=th & cells_1>=th];if(!nrow(a))next
 n0<-comp[DonorID==d & ConditionCode=='UF',unique(total_cells)];n1<-comp[DonorID==d & ConditionCode=='UA5',unique(total_cells)]
 a[,`:=`(p0=cells_0/sum(cells_0),p1=cells_1/sum(cells_1),retained_0=sum(cells_0)/n0,retained_1=sum(cells_1)/n1),by=program]
 a[,`:=`(composition=(p1-p0)*(value_1+value_0)/2,activity=(p1+p0)*(value_1-value_0)/2,total=p1*value_1-p0*value_0)]
 y<-a[,.(composition=sum(composition),activity=sum(activity),total=sum(total),retained_0=retained_0[1],retained_1=retained_1[1],n_states=.N,n0=n0,n1=n1),by=.(DonorID,program,version)]
 y[,min_state_cells:=th];k<-k+1;dec[[k]]<-y
}
de<-rbindlist(dec);stopifnot(max(abs(de$total-de$composition-de$activity))<1e-10)
write_tab(de,'composition_activity_decomposition_by_donor')
dl<-melt(de[n0>=20 & n1>=20 & retained_0>=.8 & retained_1>=.8],id.vars=c('DonorID','program','version','min_state_cells'),measure.vars=c('total','composition','activity'),variable.name='component',value.name='delta')
ds<-dl[,mean_test(delta),by=.(program,version,min_state_cells,component)]
ds[,q:=p.adjust(p,'BH'),by=.(version,min_state_cells,component)];write_tab(ds,'composition_activity_decomposition_statistics')
message('Full-gene donor-level pseudobulk models')
agg<-function(ids){
 m<-md[ids];g<-factor(paste(m$DonorID,m$ConditionCode,sep='|'))
 ind<-sparseMatrix(i=seq_along(g),j=as.integer(g),x=1,dims=c(length(g),nlevels(g)))
 cc<-cts[,ids,drop=FALSE]%*%ind;colnames(cc)<-levels(g)
 sm<-m[,.(cells=.N),by=.(DonorID,ConditionCode)];sm[,sample:=paste(DonorID,ConditionCode,sep='|')]
 list(counts=cc,meta=sm)
}
bundles<-list(Resolved=agg(which(md$resolved)),All_annotated=agg(which(md$lineage_eligible)))
saveRDS(bundles,file.path(OUT,'objects/full_gene_pseudobulk.rds'))
rm(raw,cts,e,zz,lo,sd);gc()
destats<-list();enrich<-list()
for(pop in names(bundles))for(cn in c('alpha_interaction','alpha_in_UC','UC_vs_control_fibroblasts')){
 b<-bundles[[pop]];w<-contrasts_list[[cn]];sm<-b$meta[ConditionCode%in%names(w) & cells>=20]
 donors<-sm[,.(n=.N),by=DonorID][n==length(w),DonorID];sm<-sm[DonorID%in%donors]
 if(length(donors)<3)next
 sm[,ConditionCode:=factor(ConditionCode,levels=names(w))]
 design<-model.matrix(~0+ConditionCode+factor(DonorID),sm)
 yy<-DGEList(as.matrix(b$counts[,sm$sample,drop=FALSE]));keep<-filterByExpr(yy,design,min.count=3)
 yy<-calcNormFactors(yy[keep,,keep.lib.sizes=FALSE]);yy<-estimateDisp(yy,design,robust=TRUE)
 fit<-glmQLFit(yy,design,robust=TRUE);cv<-rep(0,ncol(design));cv[seq_along(w)]<-w
 fitq<-glmQLFTest(fit,contrast=cv);de1<-as.data.table(topTags(fitq,n=Inf,sort.by='none')$table,keep.rownames='gene')
 de1[,`:=`(population=pop,contrast=cn,n_donors=length(donors),tested_genes=sum(keep))];destats[[paste(pop,cn)]]<-de1
 # Competitive gene-set testing uses all filtered genes as the measured background.
 inds1<-ids2indices(programs,rownames(yy),remove.empty=TRUE)
 if(length(inds1)){
  vo<-voom(yy,design,plot=FALSE);ca<-camera(vo,inds1,design,contrast=cv,inter.gene.cor=.01)
  ca<-as.data.table(ca,keep.rownames='program');ca[,`:=`(population=pop,contrast=cn,n_donors=length(donors))];enrich[[paste(pop,cn)]]<-ca
 }
}
deall<-rbindlist(destats);write_tab(deall,'full_gene_pseudobulk_differential_expression');write_tab(rbindlist(enrich),'full_gene_competitive_program_tests')
write_tab(deall[,.(tested_genes=.N,significant_FDR_05=sum(FDR<.05),significant_FDR_05_absLFC05=sum(FDR<.05 & abs(logFC)>=.5)),by=.(population,contrast,n_donors)],'pseudobulk_summary')
saveRDS(list(programs=programs,source=rawfile,full_gene_count=length(readRDS(file.path(OUT,'objects/full_gene_pseudobulk.rds'))[[1]]$counts[,1]),session=sessionInfo()),file.path(OUT,'objects/rna_analysis_provenance.rds'))
message('RNA analysis complete')
