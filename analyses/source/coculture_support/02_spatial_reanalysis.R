source('output/additional_analyses_2026-09-02/code/common.R')
suppressPackageStartupMessages({library(RANN);library(dbscan)})
file<-'input_data/codex/giotto_codex_results/additional_analyses/codex_extended_cells.csv'
d<-fread(file);markers<-c('CXCR4','OSM','CD16','CD11b');B<-499L
fibuc<-d[cell_type=='Fibroblast' & Diagnosis2!='Control']
cuts<-c(FAP=unname(quantile(fibuc$FAPa_cell_clr,.75)),a5=unname(quantile(fibuc$a5B1_cell_clr,.75)))
write_tab(data.table(marker=names(cuts),cutoff=cuts,definition='Pooled UC fibroblasts; 75th percentile CLR intensity'),'spatial_thresholds')
scale1<-function(x){s<-sd(x);if(!is.finite(s)||s<1e-10)return(rep(0,length(x)));(x-mean(x))/s}
nn<-function(ref,query){if(!nrow(ref))return(rep(NA_real_,nrow(query)));as.vector(nn2(ref,query,k=1)$nn.dists)}
rows<-list();counts<-list();cells<-list();nulls<-list();l<-0
for(id in unique(d$PatientID)){
 a<-d[PatientID==id];a[,isneut:=grepl('^Neutrophil',cell_type)];nf<-a[isneut==TRUE];fb<-a[cell_type=='Fibroblast']
 counts[[id]]<-data.table(PatientID=id,Diagnosis2=a$Diagnosis2[1],total_cells=nrow(a),neutrophils=nrow(nf),fibroblasts=nrow(fb),included=nrow(nf)>=30)
 if(nrow(nf)<30||nrow(fb)<20)next
 xy<-as.matrix(a[,.(x,y)]);qxy<-as.matrix(nf[,.(x,y)]);fxy<-as.matrix(fb[,.(x,y)])
 exy<-as.matrix(a[cell_type%in%c('Epithelial Cell','Enteroendocrine Cell'),.(x,y)]);vxy<-as.matrix(a[cell_type=='Endothelial Cell',.(x,y)])
 target<-fb$FAPa_cell_clr>=cuts['FAP'] & fb$a5B1_cell_clr<cuts['a5'];if(sum(target)<5)next
 neigh<-frNN(xy,eps=50,query=qxy,sort=FALSE)$id
 nFib<-vapply(neigh,function(ii)sum(a$cell_type[ii]=='Fibroblast'),numeric(1));nNeu<-vapply(neigh,function(ii)max(0,sum(a$isneut[ii])-1),numeric(1));nOther<-lengths(neigh)-nFib-nNeu-1
 nearest<-nn2(fxy,qxy,k=1)
 dd<-data.table(PatientID=id,Diagnosis2=a$Diagnosis2[1],cell_ID=nf$cell_ID,x=nf$x,y=nf$y,
  target_distance=nn(fxy[target,,drop=FALSE],qxy),fib_distance=as.vector(nearest$nn.dists),epi_distance=nn(exy,qxy),endo_distance=nn(vxy,qxy),fib_count50=nFib,neut_count50=nNeu,other_count50=nOther,
  nearest_FAP=fb$FAPa_cell_clr[as.vector(nearest$nn.idx)],nearest_a5=fb$a5B1_cell_clr[as.vector(nearest$nn.idx)])
 if(anyNA(dd))next
 covar<-cbind(fib=scale1(log1p(nFib)),neut=scale1(log1p(nNeu)),other=scale1(log1p(pmax(nOther,0))),epi=scale1(log1p(dd$epi_distance)),endo=scale1(log1p(dd$endo_distance)),dist_fib=scale1(log1p(dd$fib_distance)))
 covar<-covar[,apply(covar,2,sd)>1e-10,drop=FALSE];X<-cbind(Intercept=1,covar);qx<-qr(X)
 y<-sapply(markers,function(m)scale1(log1p(pmax(0,nf[[m]]))));colnames(y)<-markers
 yr<-qr.resid(qx,y);prox<-scale1(-log1p(dd$target_distance));xr<-qr.resid(qx,prox)
 observed<-as.vector(crossprod(xr,yr)/sum(xr^2));unadjusted<-as.vector(crossprod(prox,y)/sum(prox^2))
 for(j in seq_along(markers)){
  rows[[length(rows)+1]]<-data.table(PatientID=id,Diagnosis2=a$Diagnosis2[1],marker=markers[j],model='Unadjusted',effect=unadjusted[j],n_neut=nrow(nf))
  rows[[length(rows)+1]]<-data.table(PatientID=id,Diagnosis2=a$Diagnosis2[1],marker=markers[j],model='Density + anatomy adjusted',effect=observed[j],n_neut=nrow(nf))
 }
 cf<-scale1(dd$nearest_FAP);ca<-scale1(dd$nearest_a5)
 Z<-cbind(Intercept=1,FAP=cf,a5=ca,FAP_a5=cf*ca,covar)
 fit<-lm.fit(Z,y)
 for(term in c('FAP','a5','FAP_a5'))for(j in seq_along(markers))rows[[length(rows)+1]]<-data.table(PatientID=id,Diagnosis2=a$Diagnosis2[1],marker=markers[j],model=paste('Continuous nearest fibroblast',term),effect=fit$coefficients[term,j],n_neut=nrow(nf))
 dd[,proximity_z:=prox];for(m in markers)dd[,(m):=nf[[m]]];cells[[id]]<-dd
 # Preserve fibroblast coordinates and target counts within 200-um tiles and nearest anatomical compartment.
 fdE<-nn(exy,fxy);fdV<-nn(vxy,fxy)
 strata<-interaction(floor((fb$x-min(a$x))/200),floor((fb$y-min(a$y))/200),fdE<fdV,drop=TRUE)
 chunks<-split(seq_along(target),strata);movable<-sum(vapply(chunks,function(ii)if(length(unique(target[ii]))>1)length(ii) else 0L,integer(1)))
 nv<-matrix(NA_real_,B,length(markers))
 for(b in seq_len(B)){
  tp<-target
  for(ii in chunks)if(length(ii)>1)tp[ii]<-target[sample(ii,length(ii),replace=FALSE)]
  pp<-scale1(-log1p(nn(fxy[tp,,drop=FALSE],qxy)));rr<-qr.resid(qx,pp)
  if(sum(rr^2)>1e-10)nv[b,]<-as.vector(crossprod(rr,yr)/sum(rr^2))
 }
 no<-as.data.table(nv);setnames(no,markers);no[,`:=`(iteration=seq_len(B),PatientID=id,Diagnosis2=a$Diagnosis2[1],movable_fibroblasts=movable,total_fibroblasts=nrow(fb))];nulls[[id]]<-no
 message('Spatial ROI ',id,'; neutrophils=',nrow(nf),'; permutable fibroblasts=',movable)
}
rr<-rbindlist(rows);cc<-rbindlist(counts);nl<-rbindlist(nulls)
write_tab(cc,'spatial_sample_coverage');write_tab(rr,'spatial_ROI_effects');write_tab(rbindlist(cells),'spatial_neutrophil_context_features');write_tab(nl,'spatial_restricted_label_null')
pooled<-rr[Diagnosis2!='Control',mean_test(effect),by=.(marker,model)];pooled[,q:=p.adjust(p,'BH'),by=model]
bygroup<-rr[,mean_test(effect),by=.(Diagnosis2,marker,model)];bygroup[,q:=p.adjust(p,'BH'),by=.(Diagnosis2,model)]
write_tab(pooled,'spatial_UC_pooled_effects');write_tab(bygroup,'spatial_group_effects')
lo<-rbindlist(lapply(unique(rr[Diagnosis2!='Control',PatientID]),function(id){z<-rr[Diagnosis2!='Control' & PatientID!=id & model=='Density + anatomy adjusted',mean_test(effect),by=marker];z[,omitted_id:=id];z}))
write_tab(lo,'spatial_leave_one_ROI_out')
nullLong<-melt(nl[Diagnosis2!='Control'],id.vars=c('iteration','PatientID','Diagnosis2','movable_fibroblasts','total_fibroblasts'),measure.vars=markers,variable.name='marker',value.name='effect')
nullMean<-nullLong[,.(effect=mean(effect)),by=.(iteration,marker)]
test<-rbindlist(lapply(markers,function(m){obs<-pooled[marker==m & model=='Density + anatomy adjusted',effect];v<-nullMean[marker==m,effect];data.table(marker=m,observed=obs,null_mean=mean(v),null_low=quantile(v,.025),null_high=quantile(v,.975),p_restricted=min(1,2*min((1+sum(v>=obs))/(length(v)+1),(1+sum(v<=obs))/(length(v)+1))),B=length(v))}))
test[,q_restricted:=p.adjust(p_restricted,'BH')];write_tab(test,'spatial_restricted_label_tests');write_tab(nullMean,'spatial_pooled_null_distribution')
message('Spatial analysis complete')
