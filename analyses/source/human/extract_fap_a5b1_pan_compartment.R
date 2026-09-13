options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(dplyr); library(tidyr); library(edgeR); library(CellChat)
})
set.seed(20260815)
root <- "UC_fibroblast_neutrophil_analysis/FAP_A5B1_pan_compartment"
obj_dir <- file.path(root, "objects"); tab_dir <- file.path(root, "tables")
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
wr <- function(x,n) data.table::fwrite(as.data.frame(x),file.path(tab_dir,n),sep="\t",na="NA")

message("Loading source object")
x <- readRDS("UCGNE.RDS"); m <- x@meta.data; m$cell <- rownames(m)
cts <- x[["RNA"]]@counts; dat <- x[["RNA"]]@data
m$condition <- dplyr::case_when(m$tissue.inf=="Colon Inf UC"~"Inflamed UC",m$tissue.inf=="Colon Uninf UC"~"Uninflamed UC",m$tissue.inf=="Colon HD"~"Healthy",TRUE~NA_character_)
colon <- which(!is.na(m$condition))

fib <- m$Annotations_Level_1=="Stromal" & grepl("fibroblast",m$Annotations_Level_3,ignore.case=TRUE) & !is.na(m$condition)
fap <- as.vector(cts["FAP",]>0); a5 <- as.vector(cts["ITGA5",]>0 & cts["ITGB1",]>0)
m$fib_class <- NA_character_
m$fib_class[fib & fap & a5] <- "FAP+A5B1 double-positive"
m$fib_class[fib & fap & !a5] <- "FAP-only"
m$fib_class[fib & !fap & a5] <- "A5B1-only"
m$fib_class[fib & !fap & !a5] <- "Double-negative"
m$FAP_status <- ifelse(fib,ifelse(fap,"FAP+","FAP-"),NA_character_)
m$A5B1_status <- ifelse(fib,ifelse(a5,"A5B1+","A5B1-"),NA_character_)

epithelial <- c("Enterocytes","Enteroendocrine","Goblet","Stem Cycling","Tuft")
immune <- c("B Cell","Eosinophil","Immature Neutrophil","Macrophage and DC","Mast Cell","Monocyte and DC","Neutrophil","Plasma","T Cell and NK")
m$broad_compartment <- case_when(m$Annotations_Level_1 %in% epithelial~"Epithelial",m$Annotations_Level_1 %in% immune~"Immune",
  m$Annotations_Level_1=="Endothelial"~"Vascular",m$Annotations_Level_1 %in% c("Stromal","Pericyte")~"Stromal",TRUE~"Other")
m$receiver_type <- case_when(
  m$Annotations_Level_1 %in% c("Neutrophil","Immature Neutrophil")~"Neutrophil",
  m$Annotations_Level_1 %in% c("Macrophage and DC","Monocyte and DC")~"Myeloid APC",
  m$Annotations_Level_1=="T Cell and NK"~"T/NK",
  m$Annotations_Level_1=="B Cell"~"B cell",
  m$Annotations_Level_1=="Plasma"~"Plasma",
  m$Annotations_Level_1=="Mast Cell"~"Mast",
  m$Annotations_Level_1=="Eosinophil"~"Eosinophil",
  m$Annotations_Level_1 %in% epithelial~m$Annotations_Level_1,
  m$Annotations_Level_1=="Endothelial"~"Endothelial",
  fib~"Fibroblast",
  m$Annotations_Level_1 %in% c("Pericyte","Stromal")~"Other stromal",
  TRUE~NA_character_)

fib_ab <- m[fib,] %>% count(PatientID,condition,fib_class,FAP_status,A5B1_status,name="n") %>%
  group_by(PatientID,condition) %>% mutate(total_fib=sum(n),fraction=n/total_fib) %>% ungroup()
wr(fib_ab,"01_fibroblast_class_abundance.tsv")
overlap_ab <- bind_rows(
  m[fib,] %>% count(PatientID,condition,FAP_status,name="n") %>% group_by(PatientID,condition)%>%mutate(total=sum(n),fraction=n/total,definition="FAP")%>%ungroup()%>%rename(status=FAP_status),
  m[fib,] %>% count(PatientID,condition,A5B1_status,name="n") %>% group_by(PatientID,condition)%>%mutate(total=sum(n),fraction=n/total,definition="A5B1")%>%ungroup()%>%rename(status=A5B1_status))
wr(overlap_ab,"02_FAP_A5B1_overlapping_abundance.tsv")

data(CellChatDB.human)
lr <- CellChatDB.human$interaction %>% transmute(interaction_name,pathway_name,annotation,ligand_symbols=ligand.symbol,receptor_symbols=receptor.symbol) %>% distinct()
split_symbols <- function(z) unique(trimws(unlist(strsplit(z,",",fixed=TRUE))))
sig_genes <- unique(c(split_symbols(lr$ligand_symbols),split_symbols(lr$receptor_symbols)))
sig_genes <- intersect(sig_genes,rownames(cts)); wr(lr,"03_CellChatDB_interactions.tsv")

expr_one <- function(idx,label,condition,kind) {
  if(length(idx)<5)return(NULL)
  data.frame(group=label,condition=condition,kind=kind,gene=sig_genes,n_cells=length(idx),
    pct=Matrix::rowMeans(cts[sig_genes,idx,drop=FALSE]>0),avg=Matrix::rowMeans(dat[sig_genes,idx,drop=FALSE]))
}
message("Aggregating signaling expression")
expr <- list(); k <- 0
for(co in c("Healthy","Uninflamed UC","Inflamed UC")) {
  fi <- which(fib & m$condition==co)
  groups <- list("FAP+ fibroblast"=fi[fap[fi]],"FAP- fibroblast"=fi[!fap[fi]],
    "A5B1+ fibroblast"=fi[a5[fi]],"A5B1- fibroblast"=fi[!a5[fi]],
    "FAP+A5B1 double-positive"=fi[fap[fi]&a5[fi]],"FAP-only"=fi[fap[fi]&!a5[fi]],
    "A5B1-only"=fi[!fap[fi]&a5[fi]],"Double-negative"=fi[!fap[fi]&!a5[fi]])
  for(g in names(groups)){k<-k+1;expr[[k]]<-expr_one(groups[[g]],g,co,"Fibroblast class")}
  for(g in sort(unique(na.omit(m$receiver_type[colon])))){
    ii<-which(m$condition==co & m$receiver_type==g);k<-k+1;expr[[k]]<-expr_one(ii,g,co,"Receiver cell type")
  }
}
expr <- bind_rows(expr); wr(expr,"04_signaling_expression_summary.tsv")

programs <- list(
  Epithelial_barrier=c("EPCAM","TJP1","OCLN","CLDN1","CLDN4","CLDN7","MUC2","TFF3","KRT8","KRT18"),
  Epithelial_regeneration=c("AREG","HBEGF","SOX9","YAP1","CTGF","KRT17","KRT19","OLFM4","MKI67","PCNA"),
  Epithelial_inflammatory=c("CXCL1","CXCL2","CXCL8","CCL20","IL18","NFKBIA","REG1A","REG3A"),
  Epithelial_stress=c("XBP1","DDIT3","HSPA5","ATF4","HIF1A","ERN1"),
  Immune_inflammatory=c("IL1B","TNF","OSM","IFNG","IL6","NFKBIA","CXCL8"),
  Myeloid_remodeling=c("SPP1","MMP9","MMP12","CTSB","CTSD","VEGFA","TGFB1"),
  Cytotoxicity=c("NKG7","GNLY","GZMB","PRF1","IFNG"),
  Stromal_ECM=c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","POSTN","MMP2","MMP3","MMP14","TIMP1"),
  Stromal_inflammatory=c("FAP","PDPN","CXCL1","CXCL2","CXCL5","CXCL6","CXCL8","CSF3","IL6","ICAM1"),
  TGFb_response=c("TGFBR1","TGFBR2","SMAD2","SMAD3","SERPINE1","CTGF","COL1A1"),
  OSM_response=c("OSMR","LIFR","IL6ST","STAT3","SOCS3","CEBPD","CXCL1","CXCL8"),
  Angiogenesis=c("VEGFA","KDR","FLT1","ANGPT1","ANGPT2","TEK","VWF","EMCN")
)
score_mat <- do.call(rbind,lapply(programs,function(gs){gs<-intersect(gs,rownames(dat));Matrix::colMeans(dat[gs,colon,drop=FALSE])}))
rownames(score_mat)<-names(programs)
group_id <- paste(m$PatientID[colon],m$condition[colon],m$broad_compartment[colon],sep="__")
fac <- factor(group_id); des <- Matrix::sparseMatrix(i=seq_along(fac),j=as.integer(fac),x=1,dims=c(length(fac),nlevels(fac)),dimnames=list(NULL,levels(fac)))
agg <- score_mat %*% des; nn<-Matrix::colSums(des); agg<-sweep(agg,2,nn,"/")
prog_biopsy <- as.data.frame(t(agg)); prog_biopsy$group_id<-rownames(prog_biopsy);prog_biopsy$n_cells<-nn
prog_biopsy <- prog_biopsy %>% separate(group_id,c("PatientID","condition","compartment"),sep="__",remove=FALSE)
wr(prog_biopsy,"05_biopsy_compartment_programs.tsv")

fib_program_idx <- which(fib)
score_fib <- do.call(rbind,lapply(programs,function(gs){gs<-intersect(gs,rownames(dat));Matrix::colMeans(dat[gs,fib_program_idx,drop=FALSE])}))
fg <- paste(m$condition[fib_program_idx],m$fib_class[fib_program_idx],sep="__"); ff<-factor(fg)
fd<-Matrix::sparseMatrix(i=seq_along(ff),j=as.integer(ff),x=1,dims=c(length(ff),nlevels(ff)),dimnames=list(NULL,levels(ff)))
fa<-score_fib%*%fd; fn<-Matrix::colSums(fd);fa<-sweep(fa,2,fn,"/")
fib_prog <- as.data.frame(t(fa));fib_prog$group<-rownames(fib_prog);fib_prog$n_cells<-fn
fib_prog<-fib_prog%>%separate(group,c("condition","fib_class"),sep="__",remove=FALSE)
wr(fib_prog,"06_fibroblast_class_programs.tsv")

aggregate_pb <- function(idx,status,min_cells=20){
  gp<-paste(m$PatientID[idx],status,sep="__");tt<-table(gp);valid<-names(tt)[tt>=min_cells];use<-gp%in%valid
  if(length(valid)<4)return(NULL);f<-factor(gp[use],levels=valid);m0<-Matrix::sparseMatrix(i=seq_along(f),j=as.integer(f),x=1,dims=c(length(f),nlevels(f)),dimnames=list(NULL,levels(f)))
  pb<-cts[,idx[use],drop=FALSE]%*%m0;list(pb=pb,celln=tt[valid])
}
run_de <- function(co,status_vec,positive,label){
  idx<-which(fib&m$condition==co);st<-status_vec[idx];ag<-aggregate_pb(idx,st)
  if(is.null(ag))return(NULL);sm<-data.frame(sample=colnames(ag$pb))%>%separate(sample,c("PatientID","status"),sep="__",remove=FALSE)
  paired<-intersect(sm$PatientID[sm$status==positive],sm$PatientID[sm$status!=positive]);sm<-sm%>%filter(PatientID%in%paired)
  if(length(paired)<3)return(NULL);y<-DGEList(ag$pb[,sm$sample,drop=FALSE]);keep<-filterByExpr(y,group=sm$status,min.count=5);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
  sm$PatientID<-factor(sm$PatientID);sm$status<-factor(sm$status);design<-model.matrix(~PatientID+status,sm);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
  coefn<-grep("^status",colnames(design),value=TRUE)[1];q<-glmQLFTest(fit,coef=coefn);o<-topTags(q,n=Inf)$table;o$gene<-rownames(o);o$condition<-co;o$contrast<-label;o$n_pairs<-length(paired)
  # Ensure positive minus negative orientation.
  lev<-levels(sm$status);if(lev[2]!=positive)o$logFC<--o$logFC
  o[,c("condition","contrast","n_pairs","gene","logFC","logCPM","F","PValue","FDR")]
}
message("Running donor-paired fibroblast-state pseudobulk contrasts")
de <- bind_rows(run_de("Inflamed UC",m$FAP_status,"FAP+","FAP+ vs FAP-"),run_de("Uninflamed UC",m$FAP_status,"FAP+","FAP+ vs FAP-"),
  run_de("Inflamed UC",m$A5B1_status,"A5B1+","A5B1+ vs A5B1-"),run_de("Uninflamed UC",m$A5B1_status,"A5B1+","A5B1+ vs A5B1-"))
wr(de,"07_fibroblast_positive_state_pseudobulk_DE.tsv")

clinical_cols<-intersect(c("PatientID","condition","MAYO_ES_SCORE","HISTOLOGY_SCORE","FECAL_CAL","CRP","ENDOSCOPIC_SCORE","MEDICATION_AT_BIOPSY"),names(m))
wr(m[colon,clinical_cols,drop=FALSE]%>%distinct(PatientID,condition,.keep_all=TRUE),"08_clinical_by_biopsy.tsv")
writeLines(c(paste("Colon cells:",length(colon)),paste("Fibroblasts:",sum(fib)),paste("Signaling genes:",length(sig_genes))),file.path(root,"EXTRACTION_SUMMARY.txt"))
