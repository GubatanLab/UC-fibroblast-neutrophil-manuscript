options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(ggrepel); library(pheatmap)
})
set.seed(20260815)
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

root <- "UC_fibroblast_neutrophil_analysis"
tab_dir <- file.path(root, "tables"); fig_dir <- file.path(root, "figures"); obj_dir <- file.path(root, "objects")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(n) data.table::fread(file.path(tab_dir, n), data.table = FALSE, na.strings = c("NA", ""))
wr <- function(x, n) data.table::fwrite(as.data.frame(x), file.path(tab_dir, n), sep = "\t", na = "NA")
savep <- function(p, n, w = 8, h = 6) {
  ggsave(file.path(fig_dir, paste0(n, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(n, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
}
safe_num <- function(x) suppressWarnings(as.numeric(sub(".*?([0-9]+(?:\\.[0-9]+)?).*", "\\1", as.character(x))))

marker <- rd("03_fibroblast_marker_summary.tsv")
fib_ab <- rd("05_fibroblast_state_abundance.tsv")
neut_ab <- rd("06_neutrophil_state_abundance.tsv")
fib_sample <- rd("09_fibroblast_biopsy_programs.tsv")
neut_sample <- rd("10_neutrophil_biopsy_programs.tsv")
prog_tests <- rd("11_paired_program_tests.tsv")
fib_de <- rd("12_fibroblast_paired_pseudobulk_DE.tsv")
neut_de <- rd("13_neutrophil_paired_pseudobulk_DE.tsv")
cross <- rd("14_cross_compartment_biopsy_matrix.tsv")
cross_cor <- rd("15_cross_compartment_correlations.tsv")
clinical <- data.table::fread(file.path(obj_dir, "clinical_by_biopsy.tsv"), data.table = FALSE, na.strings = c("NA", "N/A", "None"))
emb <- readRDS(file.path(obj_dir, "fibro_neutrophil_embeddings.rds"))
fib_scores <- readRDS(file.path(obj_dir, "fibroblast_cell_scores.rds"))
neut_scores <- readRDS(file.path(obj_dir, "neutrophil_cell_scores.rds"))

fib_programs <- c("FAP_inflammatory","Alpha5Beta1_adhesion","Neutrophil_recruitment","OSM_response",
                  "ECM_remodeling","TGFb_response","NFkB_AP1","YAP_mechanotransduction")
neut_programs <- c("Recruitment_migration","Retention_aging","Degranulation","NET_associated","Oxidative_burst",
                   "Inflammatory","Interferon","Survival_immaturity","Tissue_injury")

# Complete genuinely observed biopsies with zeroes for states not seen in that biopsy,
# while never treating a missing biopsy as a zero-abundance biopsy.
complete_abundance <- function(dat, state_col) {
  samples <- dat %>% distinct(PatientID, condition)
  states <- sort(unique(dat[[state_col]]))
  base <- tidyr::crossing(samples, state = states)
  names(base)[names(base) == "state"] <- state_col
  base %>% left_join(dat, by = c("PatientID", "condition", state_col)) %>%
    group_by(PatientID, condition) %>%
    mutate(n = replace_na(n, 0), total = max(replace_na(if ("total_fib" %in% names(.)) total_fib else total_neut, 0)),
           fraction = replace_na(fraction, 0)) %>% ungroup()
}
paired_composition <- function(dat, state_col) {
  paired <- dat %>% filter(condition %in% c("Inflamed UC", "Uninflamed UC")) %>%
    distinct(PatientID, condition) %>% count(PatientID) %>% filter(n == 2) %>% pull(PatientID)
  bind_rows(lapply(unique(dat[[state_col]]), function(st) {
    z <- dat %>% filter(PatientID %in% paired, .data[[state_col]] == st,
                        condition %in% c("Inflamed UC", "Uninflamed UC")) %>%
      select(PatientID, condition, fraction) %>% pivot_wider(names_from = condition, values_from = fraction)
    z <- z %>% drop_na()
    if (nrow(z) < 3) return(NULL)
    wt <- wilcox.test(z$`Inflamed UC`, z$`Uninflamed UC`, paired = TRUE, exact = FALSE)
    data.frame(state = st, n_pairs = nrow(z), mean_uninflamed = mean(z$`Uninflamed UC`),
      mean_inflamed = mean(z$`Inflamed UC`), median_paired_change = median(z$`Inflamed UC` - z$`Uninflamed UC`),
      p_value = wt$p.value)
  })) %>% mutate(FDR = p.adjust(p_value, "BH"))
}
fib_ab_complete <- complete_abundance(fib_ab, "fib_state")
neut_ab_complete <- complete_abundance(neut_ab, "neut_state")
fib_comp <- paired_composition(fib_ab_complete, "fib_state")
neut_comp <- paired_composition(neut_ab_complete, "neut_state")
wr(fib_comp, "07_fibroblast_paired_composition_tests.tsv")
wr(neut_comp, "08_neutrophil_paired_composition_tests.tsv")

# Clinical associations
clinical <- clinical %>% mutate(
  Mayo = safe_num(MAYO_ES_SCORE), Histology = safe_num(HISTOLOGY_SCORE),
  Fecal_calprotectin = safe_num(FECAL_CAL), CRP_num = safe_num(CRP), Endoscopic = safe_num(ENDOSCOPIC_SCORE))
cdat <- cross %>% left_join(clinical, by = c("PatientID", "condition"))
pv <- intersect(c(fib_programs, neut_programs, "frac_Inflammatory fibroblast", "frac_Neutrophil OSM",
                  "frac_Neutrophil PADI4", "frac_Neutrophil CXCR4"), names(cdat))
cv <- c("Mayo","Histology","Fecal_calprotectin","CRP_num","Endoscopic")
clinical_res <- bind_rows(lapply(pv, function(p) bind_rows(lapply(cv, function(cn) {
  z <- cdat[, c(p, cn)] %>% drop_na()
  if (nrow(z) < 5 || sd(z[[1]]) == 0 || sd(z[[2]]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(z[[1]], z[[2]], method = "spearman", exact = FALSE))
  data.frame(program = p, clinical_variable = cn, n = nrow(z), rho = unname(ct$estimate), p_value = ct$p.value)
}))))
if (nrow(clinical_res)) clinical_res <- clinical_res %>% mutate(FDR = p.adjust(p_value, "BH")) %>% arrange(FDR)
wr(clinical_res, "17_clinical_correlations.tsv")

# Expression-supported, curated bidirectional ligand-receptor analysis
lr <- data.table::fread(file.path(obj_dir, "curated_lr_pairs.tsv"), data.table = FALSE)
ex <- data.table::fread(file.path(obj_dir, "lr_expression_summary.tsv"), data.table = FALSE)
score_lr <- function(sender_group, receiver_group, direction) {
  se <- ex %>% filter(group == sender_group) %>% select(gene, ligand_pct = pct, ligand_avg = avg)
  re <- ex %>% filter(group == receiver_group) %>% select(gene, receptor_pct = pct, receptor_avg = avg)
  lr %>% inner_join(se, by = c("ligand" = "gene")) %>% inner_join(re, by = c("receptor" = "gene")) %>%
    mutate(direction = direction,
           expression_score = sqrt(pmax(ligand_avg,0) * pmax(receptor_avg,0)) * sqrt(ligand_pct * receptor_pct),
           supported = ligand_pct >= .05 & receptor_pct >= .05) %>% arrange(desc(expression_score))
}
lr_screen <- bind_rows(
  score_lr("Fibroblast_Inflamed", "Neutrophil_Inflamed", "Fibroblast -> Neutrophil"),
  score_lr("Neutrophil_Inflamed", "Fibroblast_Inflamed", "Neutrophil -> Fibroblast"))
wr(lr_screen, "18_curated_ligand_receptor_screen.tsv")

# NicheNet ligand activity using paired-pseudobulk response genes when packaged matrices are available.
nn <- data.frame()
if (requireNamespace("nichenetr", quietly = TRUE)) {
  try({
    env <- new.env(); data("ligand_target_matrix", package = "nichenetr", envir = env)
    if (exists("ligand_target_matrix", envir = env)) {
      ltm <- get("ligand_target_matrix", envir = env)
      run_nn <- function(de, analysis, potential_ligands, direction) {
        d <- de %>% filter(.data$analysis == analysis)
        geneset <- intersect(d$gene[d$FDR < .10 & d$logFC > 0], colnames(ltm))
        background <- intersect(d$gene, colnames(ltm))
        lig <- intersect(potential_ligands, rownames(ltm))
        if (length(geneset) >= 10 && length(background) >= 100 && length(lig) >= 2)
          nichenetr::predict_ligand_activities(geneset, background, ltm, lig) %>% mutate(direction = direction)
        else NULL
      }
      f2n_lig <- lr_screen$ligand[lr_screen$direction == "Fibroblast -> Neutrophil" & lr_screen$supported]
      n2f_lig <- lr_screen$ligand[lr_screen$direction == "Neutrophil -> Fibroblast" & lr_screen$supported]
      nn <- bind_rows(run_nn(neut_de, "All mature neutrophils", f2n_lig, "Fibroblast -> Neutrophil"),
                      run_nn(fib_de, "All fibroblasts", n2f_lig, "Neutrophil -> Fibroblast"))
    }
  }, silent = TRUE)
}
if (nrow(nn)) wr(nn, "19_nichenet_ligand_activities.tsv") else
  writeLines("No NicheNet activity table was estimable: too few FDR<0.10 receiver response genes and/or packaged ligand-target matrix unavailable.",
             file.path(tab_dir, "19_nichenet_not_estimable.txt"))

# Cross-sectional exploratory mediation: condition -> fibroblast program -> neutrophil outcome.
# Small n and observational data make this hypothesis-generating only.
cross$inflamed <- as.integer(cross$condition == "Inflamed UC")
med_pairs <- expand.grid(mediator = intersect(c("FAP_inflammatory","Alpha5Beta1_adhesion","Neutrophil_recruitment","OSM_response"), names(cross)),
                         outcome = intersect(c("Inflammatory","NET_associated","Retention_aging","Survival_immaturity",
                           "frac_Neutrophil OSM","frac_Neutrophil PADI4"), names(cross)), stringsAsFactors = FALSE)
boot_med <- function(m, y, B = 2000) {
  z <- cross[, c("inflamed",m,y)] %>% drop_na(); names(z) <- c("x","m","y")
  if (nrow(z) < 10 || length(unique(z$x)) < 2) return(NULL)
  a <- coef(lm(m ~ x, z))[2]; b <- coef(lm(y ~ m + x, z))[2]; est <- a*b
  bs <- replicate(B, {ii <- sample(seq_len(nrow(z)), replace = TRUE); zz <- z[ii,];
    if (length(unique(zz$x)) < 2) return(NA_real_)
    coef(lm(m ~ x, zz))[2] * coef(lm(y ~ m + x, zz))[2]})
  data.frame(mediator=m, outcome=y, n=nrow(z), indirect_product=est,
             ci_low=quantile(bs,.025,na.rm=TRUE), ci_high=quantile(bs,.975,na.rm=TRUE))
}
med <- bind_rows(lapply(seq_len(nrow(med_pairs)), function(i) boot_med(med_pairs$mediator[i], med_pairs$outcome[i])))
wr(med, "20_exploratory_mediation.tsv")

# Descriptive trajectories on integrated UMAP, restricted to colon.
run_sling <- function(dat, state_col, start, end, label) {
  if (!requireNamespace("slingshot", quietly = TRUE)) return(NULL)
  d <- dat %>% filter(!is.na(.data[[state_col]]), is.finite(UMAP1), is.finite(UMAP2))
  state_n <- table(d[[state_col]])
  d <- d[d[[state_col]] %in% names(state_n)[state_n >= 10], , drop = FALSE]
  cl <- factor(d[[state_col]])
  args <- list(data = as.matrix(d[,c("UMAP1","UMAP2")]), clusterLabels = cl)
  if (start %in% levels(cl)) args$start.clus <- start
  if (end %in% levels(cl)) args$end.clus <- end
  s <- do.call(slingshot::slingshot, args)
  pt <- slingshot::slingPseudotime(s, na = TRUE)
  d$pseudotime <- apply(pt, 1, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm=TRUE))
  wr(d, paste0("21_",label,"_trajectory.tsv")); d
}
fib_path_states <- c("Resting crypt-top fibroblast", "Activated crypt-top fibroblast", "Inflammatory fibroblast")
tf <- tryCatch(run_sling(emb %>% filter(fib_state %in% fib_path_states), "fib_state",
                         "Resting crypt-top fibroblast", "Inflammatory fibroblast", "fibroblast"), error=function(e) NULL)
tn <- tryCatch(run_sling(emb, "neut_state", "Neutrophil lowRNA", "Neutrophil CXCR4", "neutrophil"), error=function(e) NULL)

# Figures
fib_cols <- c("Inflammatory fibroblast"="#D73027","Activated crypt-top fibroblast"="#FC8D59",
  "Resting crypt-top fibroblast"="#91BFDB","Crypt-bottom fibroblast"="#4575B4","LP fibroblast"="#74ADD1",
  "Adventitial fibroblast"="#984EA3","SMC-1"="#666666","SMC-2"="#999999")
neut_cols <- c("Neutrophil OSM"="#D73027","Neutrophil PADI4"="#FC8D59","Neutrophil CXCR4"="#984EA3",
  "Neutrophil MX1"="#4575B4","Neutrophil lowRNA"="#999999","Myelocyte"="#66C2A5","Promyelocyte"="#1B9E77")

p1 <- ggplot(fib_ab, aes(condition,n,fill=fib_state)) + geom_col() + facet_wrap(~PatientID,scales="free_y",ncol=7) +
  scale_fill_manual(values=fib_cols,na.value="grey80") + labs(x=NULL,y="Fibroblast cells",fill="State",title="Fibroblast recovery by biopsy") +
  theme(axis.text.x=element_text(angle=60,hjust=1),legend.position="bottom")
savep(p1,"01_fibroblast_sample_inventory",14,9)
p2 <- ggplot(marker,aes(gene,fib_state,size=pct_detected,color=mean_logexpr)) + geom_point() +
  scale_color_viridis_c(option="magma") + scale_size(range=c(.5,8)) +
  labs(x=NULL,y=NULL,size="% detected",color="Mean log expression",title="FAP, alpha5-beta1 and inflammatory fibroblast features") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
savep(p2,"02_fibroblast_marker_dotplot",11,6)
p3 <- ggplot(fib_ab_complete %>% filter(condition %in% c("Uninflamed UC","Inflamed UC")), aes(condition,fraction,group=PatientID,color=fib_state)) +
  geom_line(alpha=.45)+geom_point(size=1.5)+facet_wrap(~fib_state,scales="free_y")+scale_color_manual(values=fib_cols,guide="none")+
  labs(x=NULL,y="Fraction of fibroblasts",title="Paired fibroblast-state changes")+theme(axis.text.x=element_text(angle=30,hjust=1))
savep(p3,"03_paired_fibroblast_composition",12,8)
fl <- fib_sample %>% filter(condition %in% c("Uninflamed UC","Inflamed UC")) %>%
  select(PatientID,condition,all_of(intersect(fib_programs,names(.)))) %>% pivot_longer(-c(PatientID,condition),names_to="program",values_to="score")
p4 <- ggplot(fl,aes(condition,score,group=PatientID))+geom_line(alpha=.4,color="grey45")+geom_point(aes(color=condition),size=1.6)+
  facet_wrap(~program,scales="free_y",ncol=4)+scale_color_manual(values=c("Uninflamed UC"="#4575B4","Inflamed UC"="#D73027"),guide="none")+
  labs(x=NULL,y="Mean module expression",title="Paired fibroblast program changes")+theme(axis.text.x=element_text(angle=30,hjust=1))
savep(p4,"04_paired_fibroblast_programs",13,7)
p5 <- ggplot(neut_ab_complete %>% filter(condition %in% c("Uninflamed UC","Inflamed UC")),aes(condition,fraction,group=PatientID,color=neut_state))+
  geom_line(alpha=.4)+geom_point(size=1.4)+facet_wrap(~neut_state,scales="free_y")+scale_color_manual(values=neut_cols,guide="none")+
  labs(x=NULL,y="Fraction of neutrophils",title="Paired neutrophil-state changes")+theme(axis.text.x=element_text(angle=30,hjust=1))
savep(p5,"05_paired_neutrophil_composition",12,7)

nsp <- neut_scores %>% group_by(neut_state) %>% summarise(n=n(),across(all_of(neut_programs),mean),.groups="drop")
wr(nsp,"22_neutrophil_state_program_summary.tsv")
hm <- as.matrix(nsp[,neut_programs]); rownames(hm)<-nsp$neut_state; hm<-t(scale(t(hm))); hm[!is.finite(hm)]<-0
png(file.path(fig_dir,"06_neutrophil_state_program_heatmap.png"),2700,1800,res=300); pheatmap(hm,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),border_color=NA,main="Neutrophil functional programs by state"); dev.off()
pdf(file.path(fig_dir,"06_neutrophil_state_program_heatmap.pdf"),9,6); pheatmap(hm,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),border_color=NA,main="Neutrophil functional programs by state"); dev.off()

cm <- cross_cor %>% select(fibroblast_predictor,neutrophil_outcome,rho) %>% pivot_wider(names_from=neutrophil_outcome,values_from=rho) %>% as.data.frame()
rownames(cm)<-cm$fibroblast_predictor; cm$fibroblast_predictor<-NULL; cm<-as.matrix(cm)
png(file.path(fig_dir,"07_cross_compartment_correlation_heatmap.png"),3000,2400,res=300); pheatmap(cm,cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),breaks=seq(-1,1,length.out=102),display_numbers=TRUE,number_format="%.2f",main="Biopsy-level fibroblast-neutrophil associations"); dev.off()
pdf(file.path(fig_dir,"07_cross_compartment_correlation_heatmap.pdf"),10,8); pheatmap(cm,cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),breaks=seq(-1,1,length.out=102),display_numbers=TRUE,number_format="%.2f",main="Biopsy-level fibroblast-neutrophil associations"); dev.off()

spairs <- list(c("FAP_inflammatory","Inflammatory"),c("Neutrophil_recruitment","frac_Neutrophil OSM"),c("Alpha5Beta1_adhesion","NET_associated"),c("OSM_response","frac_Neutrophil PADI4"))
ps <- lapply(spairs,function(z) ggplot(cross,aes(x=.data[[z[1]]],y=.data[[z[2]]],color=condition,label=PatientID))+geom_point(size=2)+geom_smooth(method="lm",se=TRUE,linewidth=.6)+scale_color_manual(values=c("Uninflamed UC"="#4575B4","Inflamed UC"="#D73027"))+labs(x=z[1],y=z[2])+theme(legend.position="none"))
savep(wrap_plots(ps,ncol=2)+plot_annotation(title="Fibroblast programs versus neutrophil phenotypes"),"08_key_cross_compartment_scatterplots",11,9)

lrt <- lr_screen %>% filter(supported) %>% group_by(direction) %>% slice_max(expression_score,n=18) %>% ungroup()
p9 <- ggplot(lrt,aes(expression_score,reorder(paste(ligand,receptor,sep=" -> "),expression_score),color=direction,size=100*sqrt(ligand_pct*receptor_pct)))+geom_point()+facet_wrap(~direction,scales="free_y")+labs(x="Expression-supported interaction score",y=NULL,size="Geometric mean\n% detected",title="Prioritized bidirectional ligand-receptor interactions")+theme(legend.position="bottom")
savep(p9,"09_ligand_receptor_screen",13,8)

vol <- function(de,an,title){z<-de%>%filter(analysis==an)%>%mutate(sig=FDR<.05&abs(logFC)>=.5,label=ifelse(FDR<.01&abs(logFC)>=1,gene,NA));ggplot(z,aes(logFC,-log10(pmax(FDR,1e-300)),color=sig))+geom_point(alpha=.55,size=1)+geom_text_repel(aes(label=label),max.overlaps=20,size=2.5)+scale_color_manual(values=c("FALSE"="grey70","TRUE"="#D73027"),guide="none")+geom_vline(xintercept=c(-.5,.5),linetype=2,color="grey60")+labs(x="log2 fold-change: inflamed / uninflamed",y="-log10 FDR",title=title)}
savep(vol(fib_de,"All fibroblasts","Fibroblast paired pseudobulk DE")|vol(neut_de,"All mature neutrophils","Neutrophil paired pseudobulk DE"),"10_paired_pseudobulk_volcanoes",13,5.5)

if(nrow(clinical_res)){ch<-clinical_res%>%select(program,clinical_variable,rho)%>%pivot_wider(names_from=clinical_variable,values_from=rho)%>%as.data.frame();rownames(ch)<-ch$program;ch$program<-NULL;ch<-as.matrix(ch);png(file.path(fig_dir,"11_clinical_correlation_heatmap.png"),2500,2600,res=300);pheatmap(ch,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),breaks=seq(-1,1,length.out=102),display_numbers=TRUE,number_format="%.2f",main="Exploratory clinical associations");dev.off();pdf(file.path(fig_dir,"11_clinical_correlation_heatmap.pdf"),8,9);pheatmap(ch,color=colorRampPalette(c("#2166AC","white","#B2182B"))(101),breaks=seq(-1,1,length.out=102),display_numbers=TRUE,number_format="%.2f",main="Exploratory clinical associations");dev.off()}
if(!is.null(tf)){p12<-ggplot(tf,aes(UMAP1,UMAP2,color=pseudotime))+geom_point(size=.35,alpha=.75)+scale_color_viridis_c(na.value="grey85")+coord_equal()+labs(title="Descriptive fibroblast state trajectory",color="Pseudotime")+theme_void();savep(p12,"12_fibroblast_trajectory",7,6)}
if(!is.null(tn)){p13<-ggplot(tn,aes(UMAP1,UMAP2,color=pseudotime))+geom_point(size=.2,alpha=.65)+scale_color_viridis_c(na.value="grey85")+coord_equal()+labs(title="Descriptive colon neutrophil state trajectory",color="Pseudotime")+theme_void();savep(p13,"13_neutrophil_trajectory",7,6)}

# Compact key-results table
key <- bind_rows(
  prog_tests %>% transmute(category="Paired program",feature=paste(compartment,program,sep=": "),effect=median_change,p_value,FDR),
  cross_cor %>% slice_head(n=20) %>% transmute(category="Cross-compartment",feature=paste(fibroblast_predictor,neutrophil_outcome,sep=" vs "),effect=rho,p_value,FDR),
  clinical_res %>% slice_head(n=10) %>% transmute(category="Clinical",feature=paste(program,clinical_variable,sep=" vs "),effect=rho,p_value,FDR))
wr(key,"00_key_results.tsv")
writeLines(c("Completed donor-aware paired pseudobulk, composition and program analyses.",
  "Completed biopsy-level fibroblast-neutrophil associations, clinical correlations, curated ligand-receptor screening, descriptive trajectories and exploratory mediation.",
  "Primary inference: paired inflamed versus uninflamed UC biopsies. Cross-sectional interactions, clinical associations, trajectories and mediation are hypothesis-generating.",
  "ITGA5+ITGB1 RNA co-detection does not establish alpha5-beta1 surface protein."),file.path(root,"RUN_SUMMARY.txt"))
writeLines(capture.output(sessionInfo()),file.path(root,"sessionInfo_finalize.txt"))
