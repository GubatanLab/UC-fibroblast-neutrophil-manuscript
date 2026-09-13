source('output/additional_analyses_2026-09-02/code/common.R')
suppressPackageStartupMessages(library(SeuratObject))
o<-readRDS(file.path(ND,'objects/neutrophils_scVI_Harmony_subclustered_before_unresolved_removal.rds'))
md<-as.data.table(o@meta.data,keep.rownames='barcode');print(dim(o));print(names(md));print(table(md$neutrophil_cluster_annotation,useNA='ifany'))
bad<-c('HLA-II antigen-presenting/monocyte-like','CLC+ granulocyte-like','SPARC/CLU vascular-like outlier','B-cell-like doublet','Singleton/outlier')
md[,lineage_eligible:=!neutrophil_cluster_annotation%in%bad]
write_tab(md[ConditionCode!='PRE',.(barcode,DonorID,ConditionCode,neutrophil_cluster_annotation,neutrophil_state_annotation,lineage_eligible)],'lineage_eligibility_before_unresolved_removal')
write_tab(md[ConditionCode!='PRE',.(cells=.N),by=.(ConditionCode,neutrophil_cluster_annotation,lineage_eligible)],'lineage_audit_counts')
message('Lineage audit complete')
