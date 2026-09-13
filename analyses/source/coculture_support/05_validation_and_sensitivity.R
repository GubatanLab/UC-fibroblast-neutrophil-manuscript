source('output/additional_analyses_2026-09-02/code/common.R')
rd<-function(n)fread(file.path(OUT,'tables',paste0(n,'.csv')))
# Remove inferential statistics where fewer than three biological units are available.
for(f in c('composition_donor_contrasts','population_program_contrasts','within_state_program_contrasts','composition_activity_decomposition_statistics')){
 d<-rd(f);bad<-which(d$n<3);for(k in intersect(c('se','low','high','p','p_signflip','q','q_signflip'),names(d)))if(length(bad))set(d,i=bad,j=k,value=NA_real_)
 family_cols<-switch(f,composition_donor_contrasts=c('contrast','min_cells'),population_program_contrasts=c('contrast','min_cells','population','version'),within_state_program_contrasts=c('contrast','version'),composition_activity_decomposition_statistics=c('version','min_state_cells','component'))
 d[,q:=p.adjust(p,'BH'),by=family_cols]
 if('q_signflip'%in%names(d))d[,q_signflip:=p.adjust(p_signflip,'BH'),by=family_cols]
 stopifnot(all(is.na(d[n<3,p])))
 write_tab(d,f)
}
# Pairing sensitivity: use exactly the same recorded complete donor sets, but treat arms as independent.
ps<-rd('population_program_scores')[population=='Resolved' & version=='full_gene']
comp<-fread(file.path(CO,'Manuscript_Analyses_1_to_6/source_data/analysis2_canonical_state_composition_by_sample.csv'))[ConditionCode!='PRE'];comp[,`:=`(cells=total_cells,value=clr,feature=state)]
pa<-ps[,.(DonorID,ConditionCode,cells,value,feature=program)]
ind<-list();kk<-0
for(fam in c('Canonical composition','Full-gene program'))for(th in c(0,20,50))for(cn in c('alpha_interaction','alpha_in_UC')){
 dat<-if(fam=='Canonical composition')copy(comp) else copy(pa);w<-contrasts_list[[cn]]
 dat<-dat[ConditionCode%in%names(w) & cells>=th]
 complete<-dat[,.(nc=uniqueN(ConditionCode)),by=.(feature,DonorID)][nc==length(w)]
 dat<-merge(dat,complete[,.(feature,DonorID)],by=c('feature','DonorID'))
 for(ft in unique(dat$feature)){
  z<-dat[feature==ft,.(n=.N,m=mean(value),v=var(value)/.N),by=ConditionCode];z<-z[match(names(w),ConditionCode)]
  if(nrow(z)!=length(w)||any(z$n<3))next
  eff<-sum(w*z$m);vv<-sum(w^2*z$v);df<-vv^2/sum((w^2*z$v)^2/(z$n-1));se<-sqrt(vv)
  kk<-kk+1;ind[[kk]]<-data.table(family=fam,feature=ft,min_cells=th,contrast=cn,n_per_arm=z$n[1],effect=eff,low=eff-qt(.975,df)*se,high=eff+qt(.975,df)*se,p=2*pt(-abs(eff/se),df))
 }
}
ip<-rbindlist(ind);ip[,q:=p.adjust(p,'BH'),by=.(family,min_cells,contrast)];write_tab(ip,'recorded_pairing_sensitivity_independent_arms')
de<-rd('composition_activity_decomposition_by_donor');stopifnot(max(abs(de$total-de$composition-de$activity))<1e-10)
cs<-rd('original_vs_donor_variance_comparison');stopifnot(max(abs(cs$original_effect-cs$effect))<1e-10)
ca<-rd('clinical_site_aggregation_audit');stopifnot(max(abs(ca$delta_source-ca$delta_recomputed))<1e-10)
cov<-rd('spatial_sample_coverage');stopifnot(sum(cov$included & cov$Diagnosis2!='Control')==17)
qc<-rd('all_neutrophil_sample_quality');stopifnot(all(qc$resolved_cells<=qc$all_neutrophils),max(abs(qc$unresolved_fraction-(1-qc$resolved_cells/qc$all_neutrophils)))<1e-10)
cnt<-rd('canonical_cell_counts');stopifnot(all(qc$resolved_cells[match(paste(cnt$DonorID,cnt$ConditionCode),paste(qc$DonorID,qc$ConditionCode))]==cnt$cells))
pg<-rd('program_gene_coverage');coverage<-pg[,.(defined_genes=.N,full_genes=sum(in_full),mapping_genes=sum(in_reduced)),by=program];write_tab(coverage,'program_coverage_summary')
sp<-rd('spatial_UC_pooled_effects');sl<-rd('spatial_leave_one_ROI_out')
lr<-sl[,.(minimum_effect=min(effect),maximum_effect=max(effect),sign_stable=all(effect>0)|all(effect<0)),by=marker];write_tab(lr,'spatial_leave_one_out_summary')
checks<-data.table(check=c('Original composition effect reproduced','Decomposition sums exactly','Clinical site means reproduced','Canonical counts match full source','Unresolved fractions sum correctly','17 UC spatial IDs included','At least 3 units for inferential rows'),status='PASS')
write_tab(checks,'validation_checks')
capture.output(sessionInfo(),file=file.path(OUT,'R_session_information.txt'))
message('Validation and pairing sensitivity complete')
