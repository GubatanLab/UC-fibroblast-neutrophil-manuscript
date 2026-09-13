suppressPackageStartupMessages({library(data.table);library(Matrix);library(ggplot2);library(patchwork)})
options(stringsAsFactors=FALSE)
set.seed(20260902)
ROOT <- '.'
OUT <- file.path(ROOT,'output/additional_analyses_2026-09-02')
for(z in c('tables','figures','objects','logs'))dir.create(file.path(OUT,z),recursive=TRUE,showWarnings=FALSE)
CO <- 'input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures'
ND <- file.path(CO,'scVI_Harmony_All_Coculture/Neutrophil_Subclustering')
write_tab <- function(x,name)fwrite(x,file.path(OUT,'tables',paste0(name,'.csv')))
states <- c('PADI4 neutrophil','OSM neutrophil','CXCR4 neutrophil','MX1/ISG neutrophil')
short_state <- setNames(c('PADI4','OSM','CXCR4','MX1/ISG'),states)
conditions <- c('N','CF','UF','NA5','UA5','NNAMPT','UN','UD')
condition_labels <- c(N='Neutrophils alone',CF='Control fibroblasts',UF='UC fibroblasts',NA5='Alone + blockade',UA5='UC + blockade',NNAMPT='Alone + NAMPTi',UN='UC + NAMPTi',UD='UC + dual blockade')
contrasts_list <- list(alpha_interaction=c(N=1,NA5=-1,UF=-1,UA5=1),alpha_in_UC=c(UF=-1,UA5=1),alpha_alone=c(N=-1,NA5=1),UC_vs_control_fibroblasts=c(CF=-1,UF=1),UC_vs_alone=c(N=-1,UF=1),control_fibroblasts_vs_alone=c(N=-1,CF=1))
mean_test <- function(x){
 x<-x[is.finite(x)];n<-length(x)
 if(n<3)return(data.table(n=n,effect=if(n)mean(x) else NA_real_,se=NA_real_,low=NA_real_,high=NA_real_,p=NA_real_,p_signflip=NA_real_,same_sign=NA_integer_))
 m<-mean(x);se<-sd(x)/sqrt(n);crit<-qt(.975,n-1)
 p<-if(se>0)2*pt(-abs(m/se),n-1) else if(m==0)1 else 0
 ps<-if(n<=16){signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));mean(abs(as.vector(signs%*%x)/n)>=abs(m)-1e-12)} else NA_real_
 data.table(n=n,effect=m,se=se,low=m-crit*se,high=m+crit*se,p=p,p_signflip=ps,same_sign=sum(sign(x)==sign(m)))
}
calc_contrast <- function(dt,keys,value='value',mincells=0){
 value_column<-value
 ans<-list();ind<-list();a<-0
 for(cn in names(contrasts_list)){
  w<-contrasts_list[[cn]];d<-dt[ConditionCode%in%names(w) & cells>=mincells]
  complete<-d[,.(nc=uniqueN(ConditionCode)),by=c(keys,'DonorID')][nc==length(w)]
  d<-merge(d,complete[,!'nc'],by=c(keys,'DonorID'))
  if(!nrow(d))next
  d[,v:=get(value_column)*unname(w[ConditionCode])]
  dd<-d[,.(delta=sum(v)),by=c(keys,'DonorID')]
  dd[,contrast:=cn];dd[,min_cells:=mincells];ind[[cn]]<-dd
  ss<-dd[,mean_test(delta),by=keys];ss[,`:=`(contrast=cn,min_cells=mincells)];ans[[cn]]<-ss
 }
 list(stats=rbindlist(ans,fill=TRUE),individual=rbindlist(ind,fill=TRUE))
}
theme_pub<-theme_classic(base_size=10,base_family='Arial')+theme(plot.title=element_text(size=12,face='bold'),plot.subtitle=element_text(size=9,color='#4C5665'),axis.title=element_text(size=10),axis.text=element_text(size=9,color='#27313D'),strip.background=element_rect(fill='#F0F3F6',color=NA),strip.text=element_text(size=10,face='bold'),legend.position='bottom',plot.margin=margin(9,13,9,9),plot.tag=element_text(size=16,face='bold'))
cols<-c('PADI4'='#8056A2','OSM'='#D65A4A','CXCR4'='#2E7FAD','MX1/ISG'='#248E84')
fmtq<-function(x)ifelse(is.na(x),'NA',ifelse(x<.001,formatC(x,format='e',digits=1),sprintf('%.3f',x)))
save_plot<-function(p,name,w=11,h=8.5){
 ggsave(file.path(OUT,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=320,device=ragg::agg_png,bg='white',limitsize=FALSE)
 ggsave(file.path(OUT,'figures',paste0(name,'.svg')),p,width=w,height=h,device=svglite::svglite,bg='white',limitsize=FALSE)
 saveRDS(list(plot=p,width=w,height=h),file.path(OUT,'objects',paste0(name,'_plot.rds')))
}
