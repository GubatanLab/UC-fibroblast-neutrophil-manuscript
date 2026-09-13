source('output/additional_analyses_2026-09-02/code/common.R')
TA<-'input_data/mouse/Figure 6 TAURUS External Validation/tables'
fl<-'input_data/flow/outputs/02-gating-reanalysis'
welch<-function(a,b){
 a<-a[is.finite(a)];b<-b[is.finite(b)];na<-length(a);nb<-length(b);eff<-mean(a)-mean(b)
 va<-var(a)/na;vb<-var(b)/nb;se<-sqrt(va+vb);df<-(va+vb)^2/(va^2/(na-1)+vb^2/(nb-1))
 data.table(n_a=na,n_b=nb,effect=eff,se=se,low=eff-qt(.975,df)*se,high=eff+qt(.975,df)*se,p=2*pt(-abs(eff/se),df))
}
exact_unpaired<-function(a,b){
 z<-c(a,b);k<-length(a);n<-length(z);ix<-combn(n,k);sa<-colSums(matrix(z[ix],nrow=k));del<-sa/k-(sum(z)-sa)/(n-k)
 mean(abs(del)>=abs(mean(a)-mean(b))-1e-12)
}
hc3<-function(fit,term){X<-model.matrix(fit);ii<-solve(crossprod(X));u<-resid(fit)/(1-hatvalues(fit));vc<-ii%*%crossprod(X,X*u^2)%*%ii;s<-sqrt(diag(vc));b<-coef(fit)[term];se<-s[term];df<-df.residual(fit);data.table(effect=unname(b),se=unname(se),low=unname(b-qt(.975,df)*se),high=unname(b+qt(.975,df)*se),p=unname(2*pt(-abs(b/se),df)))}
pat<-fread(file.path(TA,'taurus_uc_patient_paired_deltas.tsv'));site<-fread(file.path(TA,'taurus_uc_site_matched_deltas.tsv'))
features<-c('activated FAP fibroblast fraction','FAP inflammatory','alpha5beta1 adhesion','neutrophil recruitment','OSM response')
pat<-pat[feature%in%features & Remission_status%in%c('Remission','Non_Remission')]
recalc<-site[feature%in%features,.(pre=mean(Pre),post=mean(Post),delta=mean(delta),n_sites=.N),by=.(Patient,Remission_status,feature)]
audit<-merge(pat,recalc,by=c('Patient','Remission_status','feature'),suffixes=c('_source','_recomputed'))
stopifnot(max(abs(audit$delta_source-audit$delta_recomputed),na.rm=TRUE)<1e-10)
write_tab(pat,'clinical_paired_patient_values');write_tab(site[feature%in%features],'clinical_site_values');write_tab(audit,'clinical_site_aggregation_audit')
res<-list();loo<-list()
for(f in features){
 d<-pat[feature==f];a<-d[Remission_status=='Remission',delta];b<-d[Remission_status=='Non_Remission',delta]
 if(length(a)<2||length(b)<2)next
 w<-welch(a,b);w[,`:=`(feature=f,model='Change difference',p_exact=exact_unpaired(a,b))];res[[paste(f,'raw')]]<-w
 d[,remission:=as.integer(Remission_status=='Remission')];fit<-lm(delta~remission+pre,data=d)
 adj<-hc3(fit,'remission');adj[,`:=`(feature=f,model='Baseline-adjusted change',n_a=length(a),n_b=length(b),p_exact=NA_real_)];res[[paste(f,'adjusted')]]<-adj
 for(id in d$Patient){x<-d[Patient!=id];aa<-x[Remission_status=='Remission',delta];bb<-x[Remission_status=='Non_Remission',delta];if(min(length(aa),length(bb))>=2){z<-welch(aa,bb);z[,`:=`(feature=f,omitted_patient=id)];loo[[paste(f,id)]]<-z}}
}
cr<-rbindlist(res,fill=TRUE);cr[,q:=p.adjust(p,'BH'),by=model];cr[,q_exact:=p.adjust(p_exact,'BH'),by=model]
write_tab(cr,'clinical_change_comparisons');write_tab(rbindlist(loo),'clinical_leave_one_patient_out')
message('Clinical analysis complete')
# Use the corrected existing flow audit, with its original biological-sample identities and caveats.
fm<-fread(file.path(fl,'revised_marker_measurements.csv'))[quantile==.99]
fq<-fread(file.path(fl,'sample_qc.csv'));fs<-fread(file.path(fl,'revised_statistical_comparisons.csv'))
write_tab(fm,'flow_existing_review_measurements');write_tab(fq,'flow_existing_sample_QC')
write_tab(fs[quantile==.99 & metric=='median_fi'],'flow_existing_complete_intensity_tests')
fr<-list();find<-list();j<-0
for(m in c('CXCR4','OSM')){
 a<-fm[experiment=='co' & marker==m]
 cc<-list(UC_fibroblasts_vs_alone=c('N+IF','N'),UC_vs_noninflamed_fibroblasts=c('N+IF','N+NF'),ATN_in_UC=c('N+IF+ATN','N+IF'),ATN_alone=c('N+ATN','N'))
 for(cn in names(cc)){
  g<-cc[[cn]];x<-a[group==g[1],median_fi];y<-a[group==g[2],median_fi];if(min(length(x),length(y))<2)next
  z<-welch(x,y);z[,`:=`(marker=m,contrast=cn,group_a=g[1],group_b=g[2],p_exact=exact_unpaired(x,y),inference='Unpaired biological samples; exploratory review measurements')]
  old<-fs[experiment=='co' & marker==m & quantile==.99 & metric=='median_fi' & group_a==g[1] & group_b==g[2]]
  z[,original_Holm_36:=if(nrow(old))old$p_holm[1] else NA_real_];j<-j+1;fr[[j]]<-z
 }
 # Four-group difference-in-differences, without assuming unconfirmed donor pairing.
 groups<-c('N','N+ATN','N+IF','N+IF+ATN');w<-c(1,-1,-1,1);dd<-a[group%in%groups]
 g<-dd[,.(n=.N,mean=mean(median_fi),var=var(median_fi)),by=group];g<-g[match(groups,group)]
 if(nrow(g)==4 && all(g$n>=2)){
  v<-g$var/g$n;eff<-sum(w*g$mean);se<-sqrt(sum(v));df<-sum(v)^2/sum(v^2/(g$n-1));j<-j+1
  fr[[j]]<-data.table(marker=m,contrast='ATN_context_interaction',effect=eff,se=se,low=eff-qt(.975,df)*se,high=eff+qt(.975,df)*se,p=2*pt(-abs(eff/se),df),n_a=sum(g$n[3:4]),n_b=sum(g$n[1:2]),inference='Unpaired four-group contrast; heteroscedastic Satterthwaite interval; exploratory')
 }
}
ff<-rbindlist(fr,fill=TRUE);ff[,targeted_q:=p.adjust(p,'BH')]
write_tab(ff,'flow_continuous_effects');write_tab(fm[experiment=='co' & marker%in%c('CXCR4','OSM')],'flow_continuous_figure_values')
message('Flow summary and new unpaired interaction analysis complete')
