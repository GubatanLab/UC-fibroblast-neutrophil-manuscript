options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(ggrepel);library(pheatmap)})
set.seed(20260815);theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank()))
root<-"UC_fibroblast_neutrophil_analysis/FAP_A5B1_pan_compartment";td<-file.path(root,"tables");fd<-file.path(root,"figures")
dir.create(fd,recursive=TRUE,showWarnings=FALSE)
rd<-function(n)data.table::fread(file.path(td,n),data.table=FALSE,na.strings=c("NA",""));wr<-function(x,n)data.table::fwrite(as.data.frame(x),file.path(td,n),sep="\t",na="NA")
savep<-function(p,n,w=9,h=6){ggsave(file.path(fd,paste0(n,".png")),p,width=w,height=h,dpi=300,bg="white");ggsave(file.path(fd,paste0(n,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white")}

ab<-rd("01_fibroblast_class_abundance.tsv");ov<-rd("02_FAP_A5B1_overlapping_abundance.tsv");lr<-rd("03_CellChatDB_interactions.tsv")
ex<-rd("04_signaling_expression_summary.tsv");bp<-rd("05_biopsy_compartment_programs.tsv");fp<-rd("06_fibroblast_class_programs.tsv");de<-rd("07_fibroblast_positive_state_pseudobulk_DE.tsv")
splitg<-function(z)trimws(strsplit(z,",",fixed=TRUE)[[1]])
lig_list<-lapply(lr$ligand_symbols,splitg);rec_list<-lapply(lr$receptor_symbols,splitg)

# Precompute multimer-aware ligand and receptor expression for each cell profile.
profiles<-split(ex,interaction(ex$condition,ex$kind,ex$group,drop=TRUE,sep="|||"))
profile_scores<-lapply(profiles,function(z){
  av<-setNames(z$avg,z$gene);pc<-setNames(z$pct,z$gene)
  calc<-function(lst,vals,fun){vapply(lst,function(gs){v<-vals[gs];if(any(is.na(v)))return(0);fun(v)},numeric(1))}
  data.frame(lig_avg=calc(lig_list,av,function(v)exp(mean(log(pmax(v,1e-8))))),lig_pct=calc(lig_list,pc,min),
    rec_avg=calc(rec_list,av,function(v)exp(mean(log(pmax(v,1e-8))))),rec_pct=calc(rec_list,pc,min))
})
key<-function(co,kind,g)paste(co,kind,g,sep="|||")
fib_groups<-c("FAP+ fibroblast","A5B1+ fibroblast","FAP+A5B1 double-positive","FAP-only","A5B1-only","Double-negative")
receiver_groups<-sort(unique(ex$group[ex$kind=="Receiver cell type"]))
ctype_comp<-function(g)case_when(g%in%c("Enterocytes","Enteroendocrine","Goblet","Stem Cycling","Tuft")~"Epithelial",
  g%in%c("Neutrophil","Myeloid APC","T/NK","B cell","Plasma","Mast","Eosinophil")~"Immune",g=="Endothelial"~"Vascular",TRUE~"Stromal")

score_pair<-function(co,sender,sender_kind,receiver,receiver_kind,direction){
  sk<-key(co,sender_kind,sender);rk<-key(co,receiver_kind,receiver);if(!sk%in%names(profile_scores)||!rk%in%names(profile_scores))return(NULL)
  s<-profile_scores[[sk]];r<-profile_scores[[rk]]
  data.frame(condition=co,direction=direction,sender=sender,receiver=receiver,target_compartment=if(direction=="Outgoing")ctype_comp(receiver)else"Fibroblast",
    interaction_name=lr$interaction_name,pathway_name=lr$pathway_name,annotation=lr$annotation,ligand=lr$ligand_symbols,receptor=lr$receptor_symbols,
    ligand_avg=s$lig_avg,ligand_pct=s$lig_pct,receptor_avg=r$rec_avg,receptor_pct=r$rec_pct,
    score=sqrt(s$lig_avg*r$rec_avg)*sqrt(s$lig_pct*r$rec_pct)) %>% filter(ligand_pct>=.05,receptor_pct>=.05,score>0)
}
message("Scoring pan-compartment ligand-receptor interactions")
res<-list();k<-0
for(co in c("Uninflamed UC","Inflamed UC"))for(fg in fib_groups)for(rg in receiver_groups){
  k<-k+1;res[[k]]<-score_pair(co,fg,"Fibroblast class",rg,"Receiver cell type","Outgoing")
  k<-k+1;res[[k]]<-score_pair(co,rg,"Receiver cell type",fg,"Fibroblast class","Incoming")
}
ints<-bind_rows(res);wr(ints,"09_pan_compartment_interaction_scores.tsv")
delta<-ints%>%select(-ligand_avg,-ligand_pct,-receptor_avg,-receptor_pct)%>%pivot_wider(names_from=condition,values_from=score,values_fill=0)%>%
  mutate(log2FC_inflamed=log2((`Inflamed UC`+1e-5)/(`Uninflamed UC`+1e-5)),score_change=`Inflamed UC`-`Uninflamed UC`)%>%arrange(desc(abs(score_change)))
wr(delta,"10_interaction_inflammation_changes.tsv")

pathway<-ints%>%group_by(condition,direction,sender,receiver,target_compartment,pathway_name,annotation)%>%
  summarise(pathway_score=sum(score),n_interactions=n(),.groups="drop")
wr(pathway,"11_pathway_scores.tsv")

# Correct paired abundance tests, separating absent state from absent biopsy.
paired_ab_test<-function(dat,state_col,value_col="fraction"){
  samples<-dat%>%filter(condition%in%c("Inflamed UC","Uninflamed UC"))%>%distinct(PatientID,condition)
  paired<-samples%>%count(PatientID)%>%filter(n==2)%>%pull(PatientID);states<-unique(dat[[state_col]])
  bind_rows(lapply(states,function(st){
    z<-samples%>%filter(PatientID%in%paired)%>%cross_join(data.frame(state=st));names(z)[ncol(z)]<-state_col
    z<-z%>%left_join(dat%>%select(PatientID,condition,all_of(state_col),value=all_of(value_col)),by=c("PatientID","condition",state_col))%>%mutate(value=replace_na(value,0))%>%
      select(PatientID,condition,value)%>%pivot_wider(names_from=condition,values_from=value)
    wt<-wilcox.test(z$`Inflamed UC`,z$`Uninflamed UC`,paired=TRUE,exact=FALSE)
    data.frame(state=st,n_pairs=nrow(z),mean_uninflamed=mean(z$`Uninflamed UC`),mean_inflamed=mean(z$`Inflamed UC`),median_change=median(z$`Inflamed UC`-z$`Uninflamed UC`),p_value=wt$p.value)
  }))%>%mutate(FDR=p.adjust(p_value,"BH"))
}
class_test<-paired_ab_test(ab,"fib_class");overlap_test<-ov%>%group_by(definition)%>%group_modify(~paired_ab_test(.x,"status"))%>%ungroup()
wr(class_test,"12_fibroblast_class_paired_tests.tsv");wr(overlap_test,"13_FAP_A5B1_paired_tests.tsv")

# Paired epithelial, immune, stromal, and vascular remodeling programs.
program_cols<-setdiff(names(bp),c("group_id","PatientID","condition","compartment","n_cells"))
applicable<-list(Epithelial=grep("^Epithelial",program_cols,value=TRUE),Immune=c("Immune_inflammatory","Myeloid_remodeling","Cytotoxicity"),
  Stromal=c("Stromal_ECM","Stromal_inflammatory","TGFb_response","OSM_response"),Vascular=c("Angiogenesis","TGFb_response","OSM_response"))
paired_program<-bind_rows(lapply(names(applicable),function(cp)bind_rows(lapply(applicable[[cp]],function(pr){
  z<-bp%>%filter(compartment==cp,condition%in%c("Inflamed UC","Uninflamed UC"),n_cells>=20)%>%select(PatientID,condition,value=all_of(pr))%>%pivot_wider(names_from=condition,values_from=value)%>%drop_na()
  if(nrow(z)<3)return(NULL);wt<-wilcox.test(z$`Inflamed UC`,z$`Uninflamed UC`,paired=TRUE,exact=FALSE)
  data.frame(compartment=cp,program=pr,n_pairs=nrow(z),mean_uninflamed=mean(z$`Uninflamed UC`),mean_inflamed=mean(z$`Inflamed UC`),median_change=median(z$`Inflamed UC`-z$`Uninflamed UC`),p_value=wt$p.value)
}))))%>%mutate(FDR=p.adjust(p_value,"BH"))%>%arrange(FDR)
wr(paired_program,"14_paired_compartment_remodeling_tests.tsv")

# Condition-adjusted, biopsy-level association of FAP/A5B1 fractions with remodeling.
pred<-ov%>%filter(status%in%c("FAP+","A5B1+"))%>%select(PatientID,condition,predictor=definition,fraction)%>%pivot_wider(names_from=predictor,values_from=fraction,names_prefix="frac_")
assoc_dat<-bp%>%filter(condition%in%c("Inflamed UC","Uninflamed UC"),n_cells>=20)%>%inner_join(pred,by=c("PatientID","condition"))
assoc<-bind_rows(lapply(c("frac_FAP","frac_A5B1"),function(px)bind_rows(lapply(names(applicable),function(cp)bind_rows(lapply(applicable[[cp]],function(pr){
  z<-assoc_dat%>%filter(compartment==cp)%>%select(condition,x=all_of(px),y=all_of(pr))%>%drop_na();if(nrow(z)<8||sd(z$x)==0||sd(z$y)==0)return(NULL)
  raw<-cor.test(z$x,z$y,method="spearman",exact=FALSE);infl<-as.integer(z$condition=="Inflamed UC")
  rx<-resid(lm(rank(x)~infl,z));ry<-resid(lm(rank(y)~infl,z));adj<-cor.test(rx,ry,method="pearson")
  data.frame(predictor=px,compartment=cp,program=pr,n=nrow(z),rho_raw=unname(raw$estimate),p_raw=raw$p.value,rho_condition_adjusted=unname(adj$estimate),p_adjusted=adj$p.value)
}))))))%>%mutate(FDR_adjusted=p.adjust(p_adjusted,"BH"))%>%arrange(FDR_adjusted)
wr(assoc,"15_condition_adjusted_remodeling_associations.tsv")

# Figures
cols<-c("FAP+A5B1 double-positive"="#B2182B","FAP-only"="#EF8A62","A5B1-only"="#2166AC","Double-negative"="#BDBDBD")
ovp<-ov%>%filter(status%in%c("FAP+","A5B1+"),condition%in%c("Uninflamed UC","Inflamed UC"))
p1<-ggplot(ovp,aes(condition,fraction,group=PatientID,color=status))+geom_line(alpha=.35)+geom_point(size=1.8)+facet_wrap(~status,scales="free_y")+
  scale_color_manual(values=c("FAP+"="#D73027","A5B1+"="#2166AC"),guide="none")+labs(x=NULL,y="Fraction of fibroblasts",title="Paired FAP+ and alpha5-beta1+ fibroblast abundance")+theme(axis.text.x=element_text(angle=25,hjust=1))
savep(p1,"01_paired_FAP_A5B1_abundance",9,5)
p2<-ggplot(ab%>%filter(condition%in%c("Uninflamed UC","Inflamed UC")),aes(condition,fraction,fill=fib_class))+geom_col()+facet_wrap(~PatientID,ncol=6)+
  scale_fill_manual(values=cols)+labs(x=NULL,y="Fibroblast composition",fill="Exclusive class",title="Exclusive FAP/alpha5-beta1 fibroblast classes by biopsy")+theme(axis.text.x=element_text(angle=55,hjust=1),legend.position="bottom")
savep(p2,"02_exclusive_fibroblast_classes",13,8)

hm<-as.matrix(fp%>%filter(condition=="Inflamed UC")%>%select(all_of(program_cols)));rownames(hm)<-fp$fib_class[fp$condition=="Inflamed UC"];hm<-t(scale(t(hm)));hm[!is.finite(hm)]<-0
png(file.path(fd,"03_fibroblast_class_program_heatmap.png"),2600,1800,res=300);pheatmap(hm,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),border_color=NA,main="Inflamed UC fibroblast programs by FAP/alpha5-beta1 class");dev.off()
pdf(file.path(fd,"03_fibroblast_class_program_heatmap.pdf"),9,6);pheatmap(hm,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),border_color=NA,main="Inflamed UC fibroblast programs by FAP/alpha5-beta1 class");dev.off()

plot_path_heat<-function(direction,name,title){
  z<-pathway%>%filter(condition=="Inflamed UC",.data$direction==direction,sender%in%c("FAP+ fibroblast","A5B1+ fibroblast","FAP-only","A5B1-only"))
  if(direction=="Incoming")z<-pathway%>%filter(condition=="Inflamed UC",.data$direction==direction,receiver%in%c("FAP+ fibroblast","A5B1+ fibroblast","FAP-only","A5B1-only"))%>%rename(fib_group=receiver,source=sender)
  else z<-z%>%rename(fib_group=sender,source=receiver)
  top<-z%>%group_by(pathway_name)%>%summarise(s=sum(pathway_score),.groups="drop")%>%slice_max(s,n=20)%>%pull(pathway_name)
  zz<-z%>%filter(pathway_name%in%top)%>%group_by(fib_group,source,pathway_name)%>%summarise(score=sum(pathway_score),.groups="drop")
  p<-ggplot(zz,aes(source,pathway_name,fill=log1p(score)))+geom_tile(color="white")+facet_wrap(~fib_group,scales="free_x")+scale_fill_viridis_c()+labs(x=NULL,y=NULL,fill="log1p score",title=title)+theme(axis.text.x=element_text(angle=60,hjust=1))
  savep(p,name,14,9)
}
plot_path_heat("Outgoing","04_outgoing_pathway_heatmap","Fibroblast outgoing signaling across epithelial, immune, vascular and stromal targets")
plot_path_heat("Incoming","05_incoming_pathway_heatmap","Incoming signals to FAP/alpha5-beta1 fibroblasts")

topd<-delta%>%filter(direction=="Outgoing",sender%in%c("FAP+ fibroblast","A5B1+ fibroblast"),`Inflamed UC`>.02)%>%group_by(sender,target_compartment)%>%slice_max(abs(score_change),n=8)%>%ungroup()
p6<-ggplot(topd,aes(score_change,reorder(paste(ligand,receptor,receiver,sep=" -> "),score_change),color=target_compartment,size=`Inflamed UC`))+geom_point()+facet_wrap(~sender,scales="free_y")+
  geom_vline(xintercept=0,color="grey60")+labs(x="Interaction score change: inflamed - uninflamed",y=NULL,size="Inflamed score",color="Target",title="Inflammation-enhanced outgoing fibroblast interactions")
savep(p6,"06_top_outgoing_interaction_changes",13,10)

longbp<-bp%>%filter(condition%in%c("Uninflamed UC","Inflamed UC"),n_cells>=20)%>%pivot_longer(all_of(program_cols),names_to="program",values_to="score")%>%
  filter((compartment=="Epithelial"&grepl("^Epithelial",program))|(compartment=="Immune"&program%in%applicable$Immune)|(compartment=="Stromal"&program%in%applicable$Stromal)|(compartment=="Vascular"&program%in%applicable$Vascular))
p7<-ggplot(longbp%>%filter(compartment=="Epithelial"),aes(condition,score,group=PatientID))+geom_line(alpha=.35,color="grey50")+geom_point(aes(color=condition),size=1.5)+facet_wrap(~program,scales="free_y")+
  scale_color_manual(values=c("Uninflamed UC"="#2166AC","Inflamed UC"="#B2182B"),guide="none")+labs(x=NULL,y="Mean program expression",title="Paired epithelial remodeling")+theme(axis.text.x=element_text(angle=25,hjust=1))
savep(p7,"07_paired_epithelial_remodeling",10,7)
p8<-ggplot(longbp%>%filter(compartment%in%c("Immune","Stromal","Vascular")),aes(condition,score,group=PatientID))+geom_line(alpha=.3,color="grey55")+geom_point(aes(color=condition),size=1.3)+facet_grid(compartment~program,scales="free_y")+
  scale_color_manual(values=c("Uninflamed UC"="#2166AC","Inflamed UC"="#B2182B"),guide="none")+labs(x=NULL,y="Mean program expression",title="Paired immune, stromal and vascular remodeling")+theme(axis.text.x=element_text(angle=25,hjust=1))
savep(p8,"08_paired_immune_stromal_vascular_remodeling",14,9)

am<-assoc%>%select(predictor,compartment,program,rho_condition_adjusted)%>%unite(feature,compartment,program,sep=": ")%>%pivot_wider(names_from=predictor,values_from=rho_condition_adjusted)%>%as.data.frame();rownames(am)<-am$feature;am$feature<-NULL;am<-as.matrix(am)
png(file.path(fd,"09_condition_adjusted_association_heatmap.png"),2200,2500,res=300);pheatmap(am,cluster_rows=TRUE,cluster_cols=FALSE,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),breaks=seq(-1,1,length.out=102),display_numbers=TRUE,number_format="%.2f",main="Condition-adjusted fibroblast-remodeling associations");dev.off()
pdf(file.path(fd,"09_condition_adjusted_association_heatmap.pdf"),7,9);pheatmap(am,cluster_rows=TRUE,cluster_cols=FALSE,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),breaks=seq(-1,1,length.out=102),display_numbers=TRUE,number_format="%.2f",main="Condition-adjusted fibroblast-remodeling associations");dev.off()

vol<-function(z,title){z<-z%>%mutate(sig=FDR<.05&abs(logFC)>.5,label=ifelse(FDR<.005&abs(logFC)>1,gene,NA));ggplot(z,aes(logFC,-log10(pmax(FDR,1e-300)),color=sig))+geom_point(alpha=.5,size=.8)+geom_text_repel(aes(label=label),max.overlaps=16,size=2.4)+scale_color_manual(values=c("FALSE"="grey75","TRUE"="#B2182B"),guide="none")+labs(x="log2 fold-change: positive / negative",y="-log10 FDR",title=title)}
vd<-de%>%filter(condition=="Inflamed UC");p10<-vol(vd%>%filter(contrast=="FAP+ vs FAP-"),"FAP+ vs FAP- fibroblasts")|vol(vd%>%filter(contrast=="A5B1+ vs A5B1-"),"Alpha5-beta1+ vs negative fibroblasts")
savep(p10,"10_positive_state_pseudobulk_volcanoes",13,5.5)

keyres<-bind_rows(overlap_test%>%transmute(category="Abundance",feature=paste(definition,state),effect=median_change,p_value,FDR),
  paired_program%>%transmute(category="Remodeling",feature=paste(compartment,program),effect=median_change,p_value,FDR),
  assoc%>%slice_head(n=20)%>%transmute(category="Adjusted association",feature=paste(predictor,compartment,program),effect=rho_condition_adjusted,p_value=p_adjusted,FDR=FDR_adjusted))
wr(keyres,"00_key_results.tsv")
writeLines(c("Pan-compartment FAP/alpha5-beta1 fibroblast analysis complete.",paste("Interaction rows:",nrow(ints)),paste("Paired remodeling tests:",nrow(paired_program)),"CellChatDB-compatible scores are expression-supported and descriptive; biopsies are inferential units for paired and association tests."),file.path(root,"RUN_SUMMARY.txt"))
