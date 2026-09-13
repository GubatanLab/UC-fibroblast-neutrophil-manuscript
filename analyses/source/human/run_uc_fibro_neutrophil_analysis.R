options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(edgeR)
  library(limma)
  library(pheatmap)
})

set.seed(20260815)
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

input_file <- "UCGNE.RDS"
out_dir <- "UC_fibroblast_neutrophil_analysis"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
obj_dir <- file.path(out_dir, "objects")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

write_tab <- function(x, name) {
  data.table::fwrite(as.data.frame(x), file.path(tab_dir, name), sep = "\t", na = "NA")
}

save_plot <- function(p, name, width = 8, height = 6) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height,
         device = cairo_pdf, bg = "white")
}

safe_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(x))))
clamp <- function(x, lo = -4, hi = 4) pmax(lo, pmin(hi, x))

message("Loading ", input_file)
obj <- readRDS(input_file)
meta <- obj@meta.data
meta$cell <- rownames(meta)
counts <- obj[["RNA"]]@counts
logdata <- obj[["RNA"]]@data
stopifnot(identical(colnames(counts), rownames(meta)))

meta$condition <- dplyr::case_when(
  meta$tissue.inf == "Colon Inf UC" ~ "Inflamed UC",
  meta$tissue.inf == "Colon Uninf UC" ~ "Uninflamed UC",
  meta$tissue.inf == "Colon HD" ~ "Healthy",
  TRUE ~ NA_character_
)
meta$condition <- factor(meta$condition, levels = c("Healthy", "Uninflamed UC", "Inflamed UC"))
meta$biopsy_id <- ifelse(is.na(meta$condition), NA_character_, paste(meta$PatientID, meta$condition, sep = "__"))
meta$fib_state <- ifelse(meta$Annotations_Level_1 == "Stromal", meta$Annotations_Level_3, NA_character_)
meta$neut_state <- ifelse(meta$Annotations_Level_1 %in% c("Neutrophil", "Immature Neutrophil"),
                          meta$Annotations_Level_3, NA_character_)

fib_idx <- which(meta$Annotations_Level_1 == "Stromal" & !is.na(meta$condition))
neut_idx <- which(meta$Annotations_Level_1 == "Neutrophil" & !is.na(meta$condition))
all_neut_idx <- which(meta$Annotations_Level_1 %in% c("Neutrophil", "Immature Neutrophil") &
                        !is.na(meta$condition))
colon_idx <- which(!is.na(meta$condition))

sample_inventory <- meta[colon_idx, ] %>%
  count(PatientID, condition, Annotations_Level_1, name = "cells") %>%
  arrange(condition, PatientID, desc(cells))
write_tab(sample_inventory, "01_sample_cell_inventory.tsv")

paired_status <- meta[colon_idx, ] %>%
  distinct(PatientID, condition) %>%
  count(PatientID, name = "n_conditions") %>%
  mutate(paired_uc = PatientID %in% intersect(
    unique(meta$PatientID[meta$condition == "Inflamed UC"]),
    unique(meta$PatientID[meta$condition == "Uninflamed UC"])
  ))
write_tab(paired_status, "02_patient_pairing.tsv")

gene_sets_fib <- list(
  FAP_inflammatory = c("FAP", "PDPN", "CXCL1", "CXCL5", "CXCL6", "CXCL8", "CSF3", "IL6", "ICAM1"),
  Alpha5Beta1_adhesion = c("ITGA5", "ITGB1", "FN1", "PTK2", "SRC", "PXN", "VCL", "TLN1", "ACTN1", "RHOA", "ROCK1", "MYL9"),
  Neutrophil_recruitment = c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CSF3"),
  OSM_response = c("OSMR", "LIFR", "IL6ST", "STAT3", "SOCS3", "CEBPD", "JUNB", "FOSL2", "CXCL1", "CXCL8"),
  ECM_remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "FN1", "POSTN", "MMP2", "MMP3", "MMP14", "TIMP1"),
  TGFb_response = c("TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF", "COL1A1", "COL3A1"),
  NFkB_AP1 = c("NFKBIA", "NFKB1", "RELA", "TNFAIP3", "JUN", "JUNB", "FOS", "FOSB", "ATF3"),
  YAP_mechanotransduction = c("YAP1", "WWTR1", "TEAD1", "CTGF", "CYR61", "ANKRD1", "AMOTL2")
)

gene_sets_neut <- list(
  Recruitment_migration = c("CXCR1", "CXCR2", "FPR1", "FPR2", "SELL", "ITGAM", "ITGB2"),
  Retention_aging = c("CXCR4", "ICAM1", "CD44", "SELL", "BCL2A1"),
  Degranulation = c("MMP8", "MMP9", "ELANE", "CTSG", "AZU1", "LTF", "DEFA3", "OLFM4"),
  NET_associated = c("PADI4", "MPO", "ELANE", "CTSG", "HIST1H2AC", "HIST1H4C"),
  Oxidative_burst = c("CYBA", "CYBB", "NCF1", "NCF2", "NCF4", "RAC2"),
  Inflammatory = c("OSM", "IL1B", "TNF", "PTGS2", "NFKBIA", "CXCL8"),
  Interferon = c("MX1", "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "OAS1"),
  Survival_immaturity = c("CSF3R", "BCL2A1", "OLFM4", "LCN2", "LTF", "MPO", "ELANE"),
  Tissue_injury = c("MMP8", "MMP9", "S100A8", "S100A9", "LCN2", "CTSG", "ELANE")
)

score_cells <- function(idx, gene_sets) {
  ans <- data.frame(cell = colnames(logdata)[idx], stringsAsFactors = FALSE)
  for (nm in names(gene_sets)) {
    gs <- intersect(gene_sets[[nm]], rownames(logdata))
    ans[[nm]] <- if (length(gs)) Matrix::colMeans(logdata[gs, idx, drop = FALSE]) else NA_real_
  }
  ans
}

message("Scoring fibroblast and neutrophil programs")
fib_scores <- score_cells(fib_idx, gene_sets_fib)
fib_scores <- cbind(meta[fib_idx, c("PatientID", "condition", "biopsy_id", "fib_state")], fib_scores)
neut_scores <- score_cells(all_neut_idx, gene_sets_neut)
neut_scores <- cbind(meta[all_neut_idx, c("PatientID", "condition", "biopsy_id", "neut_state")], neut_scores)
saveRDS(fib_scores, file.path(obj_dir, "fibroblast_cell_scores.rds"), compress = FALSE)
saveRDS(neut_scores, file.path(obj_dir, "neutrophil_cell_scores.rds"), compress = FALSE)

marker_genes <- c("FAP", "ITGA5", "ITGB1", "FN1", "PDGFRA", "PDPN", "CXCL1", "CXCL2",
                  "CXCL5", "CXCL6", "CXCL8", "CSF3", "IL6", "OSMR", "ICAM1", "VCAM1")
marker_genes <- intersect(marker_genes, rownames(counts))
marker_summ <- lapply(split(seq_along(fib_idx), meta$fib_state[fib_idx]), function(j) {
  ii <- fib_idx[j]
  data.frame(
    fib_state = meta$fib_state[ii[1]], gene = marker_genes,
    pct_detected = 100 * Matrix::rowMeans(counts[marker_genes, ii, drop = FALSE] > 0),
    mean_logexpr = Matrix::rowMeans(logdata[marker_genes, ii, drop = FALSE]),
    n_cells = length(ii)
  )
}) %>% bind_rows()
write_tab(marker_summ, "03_fibroblast_marker_summary.tsv")

alpha_joint_genes <- intersect(c("FAP", "ITGA5", "ITGB1"), rownames(counts))
alpha_joint <- as.data.frame(t(as.matrix(counts[alpha_joint_genes, fib_idx, drop = FALSE] > 0)))
alpha_joint$cell <- colnames(counts)[fib_idx]
alpha_joint <- cbind(meta[fib_idx, c("PatientID", "condition", "fib_state")], alpha_joint)
alpha_joint$FAP_pos <- if ("FAP" %in% names(alpha_joint)) alpha_joint$FAP else FALSE
alpha_joint$A5B1_RNA_pos <- alpha_joint$ITGA5 & alpha_joint$ITGB1
joint_summary <- alpha_joint %>%
  group_by(condition, fib_state) %>%
  summarise(n = n(), FAP_pct = 100 * mean(FAP_pos), A5B1_RNA_pct = 100 * mean(A5B1_RNA_pos),
            joint_pct = 100 * mean(FAP_pos & A5B1_RNA_pos), .groups = "drop")
write_tab(joint_summary, "04_FAP_ITGA5_ITGB1_coexpression.tsv")

fib_abund <- meta[fib_idx, ] %>%
  count(PatientID, condition, biopsy_id, fib_state, name = "n") %>%
  group_by(PatientID, condition, biopsy_id) %>%
  mutate(total_fib = sum(n), fraction = n / total_fib) %>% ungroup()
neut_abund <- meta[all_neut_idx, ] %>%
  count(PatientID, condition, biopsy_id, neut_state, name = "n") %>%
  group_by(PatientID, condition, biopsy_id) %>%
  mutate(total_neut = sum(n), fraction = n / total_neut) %>% ungroup()
write_tab(fib_abund, "05_fibroblast_state_abundance.tsv")
write_tab(neut_abund, "06_neutrophil_state_abundance.tsv")

paired_test <- function(dat, state_col) {
  state_sym <- rlang::sym(state_col)
  states <- unique(dat[[state_col]])
  bind_rows(lapply(states, function(st) {
    z <- dat %>% filter(!!state_sym == st, condition %in% c("Inflamed UC", "Uninflamed UC")) %>%
      select(PatientID, condition, fraction) %>%
      pivot_wider(names_from = condition, values_from = fraction, values_fill = 0)
    if (nrow(z) < 3) return(NULL)
    wt <- suppressWarnings(wilcox.test(z$`Inflamed UC`, z$`Uninflamed UC`, paired = TRUE, exact = FALSE))
    data.frame(state = st, n_pairs = nrow(z), mean_uninflamed = mean(z$`Uninflamed UC`),
               mean_inflamed = mean(z$`Inflamed UC`),
               median_paired_change = median(z$`Inflamed UC` - z$`Uninflamed UC`), p_value = wt$p.value)
  })) %>% mutate(FDR = p.adjust(p_value, "BH"))
}
fib_comp_test <- paired_test(fib_abund, "fib_state")
neut_comp_test <- paired_test(neut_abund, "neut_state")
write_tab(fib_comp_test, "07_fibroblast_paired_composition_tests.tsv")
write_tab(neut_comp_test, "08_neutrophil_paired_composition_tests.tsv")

summarize_scores <- function(score_df, program_names, state_col) {
  score_df %>%
    group_by(PatientID, condition, biopsy_id) %>%
    summarise(n_cells = n(), across(all_of(program_names), mean, na.rm = TRUE), .groups = "drop")
}
fib_sample <- summarize_scores(fib_scores, names(gene_sets_fib), "fib_state")
neut_sample <- summarize_scores(neut_scores, names(gene_sets_neut), "neut_state")

add_state_fraction <- function(sample_df, abund, state_col, prefix) {
  wide <- abund %>% select(PatientID, condition, !!rlang::sym(state_col), fraction) %>%
    pivot_wider(names_from = !!rlang::sym(state_col), values_from = fraction, values_fill = 0,
                names_prefix = prefix)
  left_join(sample_df, wide, by = c("PatientID", "condition"))
}
fib_sample <- add_state_fraction(fib_sample, fib_abund, "fib_state", "frac_")
neut_sample <- add_state_fraction(neut_sample, neut_abund, "neut_state", "frac_")
write_tab(fib_sample, "09_fibroblast_biopsy_programs.tsv")
write_tab(neut_sample, "10_neutrophil_biopsy_programs.tsv")

program_paired_tests <- function(dat, programs, compartment) {
  bind_rows(lapply(programs, function(pr) {
    z <- dat %>% filter(condition %in% c("Inflamed UC", "Uninflamed UC")) %>%
      select(PatientID, condition, value = all_of(pr)) %>%
      pivot_wider(names_from = condition, values_from = value) %>% drop_na()
    if (nrow(z) < 3) return(NULL)
    wt <- wilcox.test(z$`Inflamed UC`, z$`Uninflamed UC`, paired = TRUE, exact = FALSE)
    data.frame(compartment = compartment, program = pr, n_pairs = nrow(z),
               mean_uninflamed = mean(z$`Uninflamed UC`), mean_inflamed = mean(z$`Inflamed UC`),
               median_change = median(z$`Inflamed UC` - z$`Uninflamed UC`), p_value = wt$p.value)
  })) %>% mutate(FDR = p.adjust(p_value, "BH"))
}
program_tests <- bind_rows(
  program_paired_tests(fib_sample, names(gene_sets_fib), "Fibroblast"),
  program_paired_tests(neut_sample, names(gene_sets_neut), "Neutrophil")
)
write_tab(program_tests, "11_paired_program_tests.tsv")

aggregate_counts <- function(cell_idx, group, min_cells = 20) {
  tt <- table(group)
  valid <- names(tt)[tt >= min_cells]
  if (!length(valid)) return(list(counts = Matrix::Matrix(0, nrow = nrow(counts), ncol = 0, sparse = TRUE), cell_counts = tt[FALSE]))
  use <- group %in% valid
  group <- factor(group[use], levels = valid)
  design <- Matrix::sparseMatrix(i = seq_along(group), j = as.integer(group), x = 1,
                                 dims = c(length(group), nlevels(group)),
                                 dimnames = list(NULL, levels(group)))
  pb <- counts[, cell_idx[use], drop = FALSE] %*% design
  list(counts = pb, cell_counts = tt[valid])
}

run_paired_edger <- function(cell_idx, state = NULL, state_vec = NULL, min_cells = 20, label = "all") {
  if (!is.null(state)) cell_idx <- cell_idx[state_vec[cell_idx] == state]
  mm <- meta[cell_idx, ]
  keep_cond <- mm$condition %in% c("Inflamed UC", "Uninflamed UC")
  cell_idx <- cell_idx[keep_cond]
  mm <- mm[keep_cond, ]
  group <- paste(mm$PatientID, mm$condition, sep = "__")
  ag <- aggregate_counts(cell_idx, group, min_cells)
  if (ncol(ag$counts) == 0) return(NULL)
  sm <- data.frame(sample = colnames(ag$counts)) %>%
    separate(sample, into = c("PatientID", "condition"), sep = "__", extra = "merge", remove = FALSE)
  paired <- intersect(sm$PatientID[sm$condition == "Inflamed UC"], sm$PatientID[sm$condition == "Uninflamed UC"])
  sm <- sm %>% filter(PatientID %in% paired)
  if (length(paired) < 3) return(NULL)
  y <- DGEList(ag$counts[, sm$sample, drop = FALSE])
  keep_gene <- filterByExpr(y, group = sm$condition, min.count = 5)
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  sm$PatientID <- factor(sm$PatientID)
  sm$condition <- factor(sm$condition, levels = c("Uninflamed UC", "Inflamed UC"))
  design <- model.matrix(~ PatientID + condition, sm)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "conditionInflamed UC")
  out <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  out$gene <- rownames(out)
  out$analysis <- label
  out$n_pairs <- length(paired)
  out[, c("analysis", "n_pairs", "gene", "logFC", "logCPM", "F", "PValue", "FDR")]
}

message("Running paired donor-pseudobulk differential expression")
fib_de <- run_paired_edger(fib_idx, label = "All fibroblasts")
for (st in unique(na.omit(meta$fib_state[fib_idx]))) {
  z <- run_paired_edger(fib_idx, state = st, state_vec = meta$fib_state, label = st)
  if (!is.null(z)) fib_de <- bind_rows(fib_de, z)
}
neut_de <- run_paired_edger(neut_idx, label = "All mature neutrophils")
for (st in unique(na.omit(meta$neut_state[neut_idx]))) {
  z <- run_paired_edger(neut_idx, state = st, state_vec = meta$neut_state, label = st)
  if (!is.null(z)) neut_de <- bind_rows(neut_de, z)
}
write_tab(fib_de, "12_fibroblast_paired_pseudobulk_DE.tsv")
write_tab(neut_de, "13_neutrophil_paired_pseudobulk_DE.tsv")

message("Cross-compartment biopsy-level associations")
cross <- inner_join(fib_sample, neut_sample, by = c("PatientID", "condition"), suffix = c("_fib", "_neut")) %>%
  filter(condition %in% c("Inflamed UC", "Uninflamed UC"), n_cells_fib >= 20, n_cells_neut >= 20)
fib_vars <- intersect(c(names(gene_sets_fib), "frac_Inflammatory fibroblast", "frac_Activated crypt-top fibroblast"), names(cross))
neut_vars <- intersect(c(names(gene_sets_neut), "frac_Neutrophil OSM", "frac_Neutrophil PADI4",
                         "frac_Neutrophil CXCR4", "frac_Neutrophil MX1"), names(cross))
cross_cor <- bind_rows(lapply(fib_vars, function(fv) bind_rows(lapply(neut_vars, function(nv) {
  z <- cross[, c(fv, nv)] %>% drop_na()
  if (nrow(z) < 5 || sd(z[[1]]) == 0 || sd(z[[2]]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(z[[1]], z[[2]], method = "spearman", exact = FALSE))
  data.frame(fibroblast_predictor = fv, neutrophil_outcome = nv, n_biopsies = nrow(z),
             rho = unname(ct$estimate), p_value = ct$p.value)
}))))
if (nrow(cross_cor)) {
  cross_cor <- cross_cor %>% mutate(FDR = p.adjust(p_value, "BH")) %>% arrange(FDR)
} else {
  cross_cor <- data.frame(fibroblast_predictor = character(), neutrophil_outcome = character(),
                          n_biopsies = integer(), rho = numeric(), p_value = numeric(), FDR = numeric())
}
write_tab(cross, "14_cross_compartment_biopsy_matrix.tsv")
write_tab(cross_cor, "15_cross_compartment_correlations.tsv")

paired_delta <- function(dat, vars) {
  dat %>% select(PatientID, condition, all_of(vars)) %>%
    pivot_wider(names_from = condition, values_from = all_of(vars), names_sep = "__") %>%
    mutate(across(ends_with("__Inflamed UC"), identity))
}
delta_rows <- list()
for (fv in fib_vars) for (nv in neut_vars) {
  z <- cross %>% select(PatientID, condition, all_of(c(fv, nv))) %>%
    pivot_wider(names_from = condition, values_from = all_of(c(fv, nv)), names_sep = "__")
  needed <- c(paste0(fv, "__Inflamed UC"), paste0(fv, "__Uninflamed UC"),
              paste0(nv, "__Inflamed UC"), paste0(nv, "__Uninflamed UC"))
  if (!all(needed %in% names(z))) next
  dx <- z[[needed[1]]] - z[[needed[2]]]
  dy <- z[[needed[3]]] - z[[needed[4]]]
  ok <- is.finite(dx) & is.finite(dy)
  if (sum(ok) < 5 || sd(dx[ok]) == 0 || sd(dy[ok]) == 0) next
  ct <- suppressWarnings(cor.test(dx[ok], dy[ok], method = "spearman", exact = FALSE))
  delta_rows[[length(delta_rows) + 1]] <- data.frame(fibroblast_predictor = fv,
    neutrophil_outcome = nv, n_pairs = sum(ok), rho_delta = unname(ct$estimate), p_value = ct$p.value)
}
delta_cor <- bind_rows(delta_rows)
if (nrow(delta_cor)) {
  delta_cor <- delta_cor %>% mutate(FDR = p.adjust(p_value, "BH")) %>% arrange(FDR)
} else {
  delta_cor <- data.frame(fibroblast_predictor = character(), neutrophil_outcome = character(),
                          n_pairs = integer(), rho_delta = numeric(), p_value = numeric(), FDR = numeric())
}
write_tab(delta_cor, "16_within_patient_delta_correlations.tsv")

message("Clinical associations")
clinical <- meta[colon_idx, ] %>%
  distinct(PatientID, condition, .keep_all = TRUE) %>%
  transmute(PatientID, condition,
            Mayo = safe_num(MAYO_ES_SCORE), Histology = safe_num(HISTOLOGY_SCORE),
            Fecal_calprotectin = safe_num(FECAL_CAL), CRP = safe_num(CRP),
            Endoscopic = safe_num(ENDOSCOPIC_SCORE), Age = safe_num(AGE),
            AntiTNF = as.character(AntiTNF), Vedolizumab = as.character(Vedolizumab))
clinical_dat <- cross %>% left_join(clinical, by = c("PatientID", "condition"))
prog_vars <- c(fib_vars, neut_vars)
clinical_vars <- c("Mayo", "Histology", "Fecal_calprotectin", "CRP", "Endoscopic")
clinical_cor <- bind_rows(lapply(prog_vars, function(pv) bind_rows(lapply(clinical_vars, function(cv) {
  z <- clinical_dat[, c(pv, cv)] %>% drop_na()
  if (nrow(z) < 5 || sd(z[[1]]) == 0 || sd(z[[2]]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(z[[1]], z[[2]], method = "spearman", exact = FALSE))
  data.frame(program = pv, clinical_variable = cv, n = nrow(z), rho = unname(ct$estimate), p_value = ct$p.value)
}))))
if (nrow(clinical_cor)) {
  clinical_cor <- clinical_cor %>% mutate(FDR = p.adjust(p_value, "BH")) %>% arrange(FDR)
} else {
  clinical_cor <- data.frame(program = character(), clinical_variable = character(), n = integer(),
                             rho = numeric(), p_value = numeric(), FDR = numeric())
}
write_tab(clinical_cor, "17_clinical_correlations.tsv")

message("Ligand-receptor and NicheNet-compatible screening")
lr <- NULL
try({
  data("lr_network", package = "nichenetr", envir = environment())
  if (exists("lr_network")) lr <- get("lr_network")
}, silent = TRUE)
if (is.null(lr)) {
  lr <- data.frame(from = c("CXCL1","CXCL2","CXCL5","CXCL6","CXCL8","CSF3","IL6","OSM","IL1B","TNF","TGFB1","FN1"),
                   to = c("CXCR2","CXCR2","CXCR2","CXCR1","CXCR1","CSF3R","IL6R","OSMR","IL1R1","TNFRSF1A","TGFBR2","ITGA5"))
}
names(lr)[1:2] <- c("ligand", "receptor")
lr <- distinct(lr, ligand, receptor)

expr_summary <- function(idx) {
  genes <- intersect(unique(c(lr$ligand, lr$receptor)), rownames(counts))
  data.frame(gene = genes,
             pct = Matrix::rowMeans(counts[genes, idx, drop = FALSE] > 0),
             avg = Matrix::rowMeans(logdata[genes, idx, drop = FALSE]))
}
fib_inf_idx <- fib_idx[meta$condition[fib_idx] == "Inflamed UC"]
neut_inf_idx <- neut_idx[meta$condition[neut_idx] == "Inflamed UC"]
fib_expr <- expr_summary(fib_inf_idx)
neut_expr <- expr_summary(neut_inf_idx)

score_lr_direction <- function(sender, receiver, direction) {
  lr %>% inner_join(sender, by = c("ligand" = "gene")) %>%
    rename(ligand_pct = pct, ligand_avg = avg) %>%
    inner_join(receiver, by = c("receptor" = "gene")) %>%
    rename(receptor_pct = pct, receptor_avg = avg) %>%
    filter(ligand_pct >= 0.05, receptor_pct >= 0.05) %>%
    mutate(direction = direction,
           expression_score = sqrt(pmax(ligand_avg, 0) * pmax(receptor_avg, 0)) *
             sqrt(ligand_pct * receptor_pct)) %>%
    arrange(desc(expression_score))
}
lr_screen <- bind_rows(score_lr_direction(fib_expr, neut_expr, "Fibroblast -> Neutrophil"),
                       score_lr_direction(neut_expr, fib_expr, "Neutrophil -> Fibroblast"))
write_tab(lr_screen, "18_ligand_receptor_expression_screen.tsv")

nichenet_activity <- list()
try({
  data("ligand_target_matrix", package = "nichenetr", envir = environment())
  if (exists("ligand_target_matrix") && requireNamespace("nichenetr", quietly = TRUE)) {
    ltm <- get("ligand_target_matrix")
    run_nn <- function(de_tab, receiver_idx, sender_expr, direction) {
      de0 <- de_tab %>% filter(analysis == unique(analysis)[1])
      geneset <- de0$gene[de0$FDR < 0.1 & de0$logFC > 0]
      bg <- rownames(counts)[Matrix::rowSums(counts[, receiver_idx, drop = FALSE] > 0) >= max(10, 0.05 * length(receiver_idx))]
      lig <- intersect(sender_expr$gene[sender_expr$pct >= 0.05], rownames(ltm))
      geneset <- intersect(geneset, colnames(ltm)); bg <- intersect(bg, colnames(ltm))
      if (length(geneset) >= 10 && length(lig) >= 2)
        nichenetr::predict_ligand_activities(geneset, bg, ltm, lig) %>% mutate(direction = direction)
      else NULL
    }
    nichenet_activity[[1]] <- run_nn(neut_de, neut_inf_idx, fib_expr, "Fibroblast -> Neutrophil")
    nichenet_activity[[2]] <- run_nn(fib_de, fib_inf_idx, neut_expr, "Neutrophil -> Fibroblast")
  }
}, silent = TRUE)
nichenet_activity <- bind_rows(nichenet_activity)
if (nrow(nichenet_activity)) write_tab(nichenet_activity, "19_nichenet_ligand_activities.tsv")

message("Exploratory mediation")
mediation_pairs <- expand.grid(
  mediator = intersect(c("FAP_inflammatory", "Alpha5Beta1_adhesion", "Neutrophil_recruitment", "OSM_response"), fib_vars),
  outcome = intersect(c("Inflammatory", "NET_associated", "Retention_aging", "frac_Neutrophil OSM", "frac_Neutrophil PADI4"), neut_vars),
  stringsAsFactors = FALSE
)
boot_product <- function(dat, mediator, outcome, B = 1000) {
  z <- dat %>% select(PatientID, condition, all_of(c(mediator, outcome))) %>%
    pivot_wider(names_from = condition, values_from = all_of(c(mediator, outcome)), names_sep = "__")
  req <- c(paste0(mediator, "__Inflamed UC"), paste0(mediator, "__Uninflamed UC"),
           paste0(outcome, "__Inflamed UC"), paste0(outcome, "__Uninflamed UC"))
  if (!all(req %in% names(z))) return(NULL)
  dm <- z[[req[1]]] - z[[req[2]]]; dy <- z[[req[3]]] - z[[req[4]]]
  ok <- is.finite(dm) & is.finite(dy); dm <- dm[ok]; dy <- dy[ok]
  if (length(dm) < 6 || sd(dm) == 0 || sd(dy) == 0) return(NULL)
  # In a paired change model, a is the condition-to-mediator mean change and b is
  # the slope relating mediator change to outcome change. This remains exploratory.
  a <- mean(dm); b <- coef(lm(dy ~ dm))[2]; est <- a * b
  boots <- replicate(B, {i <- sample(seq_along(dm), replace = TRUE); mean(dm[i]) * coef(lm(dy[i] ~ dm[i]))[2]})
  data.frame(mediator = mediator, outcome = outcome, n_pairs = length(dm), indirect_product = est,
             ci_low = quantile(boots, 0.025, na.rm = TRUE), ci_high = quantile(boots, 0.975, na.rm = TRUE))
}
mediation_res <- bind_rows(lapply(seq_len(nrow(mediation_pairs)), function(i)
  boot_product(cross, mediation_pairs$mediator[i], mediation_pairs$outcome[i])))
if (nrow(mediation_res)) write_tab(mediation_res, "20_exploratory_mediation.tsv")

message("Trajectory analyses")
trajectory_results <- list()
run_sling <- function(idx, cluster, start, end, name) {
  if (!requireNamespace("slingshot", quietly = TRUE)) return(NULL)
  reds <- names(obj@reductions)
  cand <- reds[grepl("umap", reds, ignore.case = TRUE)]
  if (!length(cand)) return(NULL)
  red <- cand[which.max(vapply(cand, function(z) nrow(obj@reductions[[z]]@cell.embeddings), numeric(1)))]
  emb <- obj@reductions[[red]]@cell.embeddings[colnames(obj)[idx], 1:2, drop = FALSE]
  keep <- complete.cases(emb) & !is.na(cluster)
  emb <- emb[keep, , drop = FALSE]; cl <- factor(cluster[keep])
  args <- list(data = emb, clusterLabels = cl)
  if (start %in% levels(cl)) args$start.clus <- start
  if (end %in% levels(cl)) args$end.clus <- end
  sds <- do.call(slingshot::slingshot, args)
  pt <- slingshot::slingPseudotime(sds, na = TRUE)
  out <- data.frame(cell = rownames(emb), cluster = cl, UMAP1 = emb[,1], UMAP2 = emb[,2],
                    pseudotime = apply(pt, 1, function(z) if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)))
  write_tab(out, paste0("21_", name, "_trajectory.tsv"))
  out
}
trajectory_results$fib <- tryCatch(run_sling(fib_idx, meta$fib_state[fib_idx],
  "Resting crypt-top fibroblast", "Inflammatory fibroblast", "fibroblast"), error = function(e) NULL)
trajectory_results$neut <- tryCatch(run_sling(neut_idx, meta$neut_state[neut_idx],
  "Neutrophil lowRNA", "Neutrophil CXCR4", "neutrophil"), error = function(e) NULL)

message("Generating figures")
state_cols_fib <- c("Inflammatory fibroblast"="#D73027", "Activated crypt-top fibroblast"="#FC8D59",
                    "Resting crypt-top fibroblast"="#91BFDB", "Crypt-bottom fibroblast"="#4575B4",
                    "LP fibroblast"="#74ADD1", "Adventitial fibroblast"="#984EA3",
                    "SMC-1"="#666666", "SMC-2"="#999999")
state_cols_neut <- c("Neutrophil OSM"="#D73027", "Neutrophil PADI4"="#FC8D59",
                     "Neutrophil CXCR4"="#984EA3", "Neutrophil MX1"="#4575B4",
                     "Neutrophil lowRNA"="#999999", "Myelocyte"="#66C2A5", "Promyelocyte"="#1B9E77")

p1 <- ggplot(fib_abund, aes(condition, n, fill = fib_state)) +
  geom_col(position = "stack") + facet_wrap(~PatientID, scales = "free_y", ncol = 7) +
  scale_fill_manual(values = state_cols_fib, na.value = "grey80") +
  labs(x = NULL, y = "Fibroblast cells", fill = "State", title = "Fibroblast recovery by biopsy") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")
save_plot(p1, "01_fibroblast_sample_inventory", 14, 9)

p2 <- ggplot(marker_summ, aes(gene, fib_state, size = pct_detected, color = mean_logexpr)) +
  geom_point() + scale_color_viridis_c(option = "magma") + scale_size(range = c(0.5, 8)) +
  labs(x = NULL, y = NULL, size = "% detected", color = "Mean log expression",
       title = "FAP, alpha5-beta1 and inflammatory fibroblast features") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p2, "02_fibroblast_marker_dotplot", 11, 6)

fib_plot_dat <- fib_abund %>% filter(condition %in% c("Uninflamed UC", "Inflamed UC"))
p3 <- ggplot(fib_plot_dat, aes(condition, fraction, group = PatientID, color = fib_state)) +
  geom_line(alpha = 0.45) + geom_point(size = 1.5) + facet_wrap(~fib_state, scales = "free_y") +
  scale_color_manual(values = state_cols_fib, guide = "none") +
  labs(x = NULL, y = "Fraction of fibroblasts", title = "Paired fibroblast-state changes") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p3, "03_paired_fibroblast_composition", 12, 8)

fib_long <- fib_sample %>% filter(condition %in% c("Uninflamed UC", "Inflamed UC")) %>%
  select(PatientID, condition, all_of(names(gene_sets_fib))) %>%
  pivot_longer(all_of(names(gene_sets_fib)), names_to = "program", values_to = "score")
p4 <- ggplot(fib_long, aes(condition, score, group = PatientID)) +
  geom_line(alpha = 0.4, color = "grey45") + geom_point(aes(color = condition), size = 1.6) +
  facet_wrap(~program, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Uninflamed UC"="#4575B4", "Inflamed UC"="#D73027"), guide = "none") +
  labs(x = NULL, y = "Mean module expression", title = "Paired fibroblast program changes") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p4, "04_paired_fibroblast_programs", 13, 7)

neut_plot_dat <- neut_abund %>% filter(condition %in% c("Uninflamed UC", "Inflamed UC"))
p5 <- ggplot(neut_plot_dat, aes(condition, fraction, group = PatientID, color = neut_state)) +
  geom_line(alpha = 0.4) + geom_point(size = 1.4) + facet_wrap(~neut_state, scales = "free_y") +
  scale_color_manual(values = state_cols_neut, guide = "none") +
  labs(x = NULL, y = "Fraction of neutrophils", title = "Paired neutrophil-state changes") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p5, "05_paired_neutrophil_composition", 12, 7)

neut_state_program <- neut_scores %>% group_by(neut_state) %>%
  summarise(n = n(), across(all_of(names(gene_sets_neut)), mean), .groups = "drop")
write_tab(neut_state_program, "22_neutrophil_state_program_summary.tsv")
hm <- as.matrix(neut_state_program[, names(gene_sets_neut)])
rownames(hm) <- neut_state_program$neut_state
hm <- t(scale(t(hm))); hm[!is.finite(hm)] <- 0
png(file.path(fig_dir, "06_neutrophil_state_program_heatmap.png"), width = 2700, height = 1800, res = 300)
pheatmap(hm, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("#2166AC","white","#B2182B"))(101),
         main = "Neutrophil functional programs by state", border_color = NA)
dev.off()
pdf(file.path(fig_dir, "06_neutrophil_state_program_heatmap.pdf"), width = 9, height = 6)
pheatmap(hm, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("#2166AC","white","#B2182B"))(101),
         main = "Neutrophil functional programs by state", border_color = NA)
dev.off()

if (nrow(cross_cor)) {
  cm <- cross_cor %>% select(fibroblast_predictor, neutrophil_outcome, rho) %>%
    pivot_wider(names_from = neutrophil_outcome, values_from = rho) %>% as.data.frame()
  rownames(cm) <- cm$fibroblast_predictor; cm$fibroblast_predictor <- NULL; cm <- as.matrix(cm)
  png(file.path(fig_dir, "07_cross_compartment_correlation_heatmap.png"), width = 2800, height = 2200, res = 300)
  pheatmap(cm, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("#2166AC","white","#B2182B"))(101), breaks = seq(-1,1,length.out=102),
           display_numbers = TRUE, number_format = "%.2f", border_color = "grey90",
           main = "Biopsy-level fibroblast-neutrophil associations")
  dev.off()
  pdf(file.path(fig_dir, "07_cross_compartment_correlation_heatmap.pdf"), width = 10, height = 8)
  pheatmap(cm, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("#2166AC","white","#B2182B"))(101), breaks = seq(-1,1,length.out=102),
           display_numbers = TRUE, number_format = "%.2f", border_color = "grey90",
           main = "Biopsy-level fibroblast-neutrophil associations")
  dev.off()
}

scatter_pairs <- list(c("FAP_inflammatory", "Inflammatory"), c("Neutrophil_recruitment", "frac_Neutrophil OSM"),
                      c("Alpha5Beta1_adhesion", "NET_associated"), c("OSM_response", "frac_Neutrophil PADI4"))
scatters <- lapply(scatter_pairs, function(z) {
  if (!all(z %in% names(cross))) return(ggplot() + theme_void())
  ggplot(cross, aes(x = .data[[z[1]]], y = .data[[z[2]]], color = condition, label = PatientID)) +
    geom_point(size = 2) + geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
    scale_color_manual(values = c("Uninflamed UC"="#4575B4", "Inflamed UC"="#D73027")) +
    labs(x = z[1], y = z[2]) + theme(legend.position = "none")
})
p8 <- wrap_plots(scatters, ncol = 2) + plot_annotation(title = "Fibroblast programs versus neutrophil phenotypes")
save_plot(p8, "08_key_cross_compartment_scatterplots", 11, 9)

lr_top <- lr_screen %>% group_by(direction) %>% slice_max(expression_score, n = 25) %>% ungroup()
p9 <- ggplot(lr_top, aes(expression_score, reorder(paste(ligand, receptor, sep = " -> "), expression_score),
                         color = direction, size = 100 * sqrt(ligand_pct * receptor_pct))) +
  geom_point() + facet_wrap(~direction, scales = "free_y") +
  labs(x = "Expression-supported interaction score", y = NULL, size = "Geometric mean\n% detected",
       title = "Prioritized bidirectional ligand-receptor interactions") +
  theme(legend.position = "bottom")
save_plot(p9, "09_ligand_receptor_screen", 13, 10)

volcano_one <- function(de, analysis_name, title) {
  z <- de %>% filter(analysis == analysis_name) %>%
    mutate(sig = FDR < 0.05 & abs(logFC) >= 0.5, label = ifelse(FDR < 0.01 & abs(logFC) >= 1, gene, NA))
  ggplot(z, aes(logFC, -log10(pmax(FDR, 1e-300)), color = sig)) + geom_point(alpha = 0.55, size = 1) +
    geom_text_repel(aes(label = label), max.overlaps = 20, size = 2.5) +
    scale_color_manual(values = c("FALSE"="grey70", "TRUE"="#D73027"), guide = "none") +
    geom_vline(xintercept = c(-0.5,0.5), linetype = 2, color = "grey60") +
    labs(x = "log2 fold-change: inflamed / uninflamed", y = "-log10 FDR", title = title)
}
p10 <- volcano_one(fib_de, "All fibroblasts", "Fibroblast paired pseudobulk DE") |
  volcano_one(neut_de, "All mature neutrophils", "Neutrophil paired pseudobulk DE")
save_plot(p10, "10_paired_pseudobulk_volcanoes", 13, 5.5)

if (nrow(clinical_cor)) {
  clinm <- clinical_cor %>% select(program, clinical_variable, rho) %>%
    pivot_wider(names_from = clinical_variable, values_from = rho) %>% as.data.frame()
  rownames(clinm) <- clinm$program; clinm$program <- NULL; clinm <- as.matrix(clinm)
  png(file.path(fig_dir, "11_clinical_correlation_heatmap.png"), width = 2500, height = 2600, res = 300)
  pheatmap(clinm, cluster_rows = TRUE, cluster_cols = FALSE,
           color = colorRampPalette(c("#2166AC","white","#B2182B"))(101), breaks = seq(-1,1,length.out=102),
           display_numbers = TRUE, number_format = "%.2f", border_color = "grey90",
           main = "Exploratory clinical associations")
  dev.off()
  pdf(file.path(fig_dir, "11_clinical_correlation_heatmap.pdf"), width = 8, height = 9)
  pheatmap(clinm, cluster_rows = TRUE, cluster_cols = FALSE,
           color = colorRampPalette(c("#2166AC","white","#B2182B"))(101), breaks = seq(-1,1,length.out=102),
           display_numbers = TRUE, number_format = "%.2f", border_color = "grey90",
           main = "Exploratory clinical associations")
  dev.off()
}

if (!is.null(trajectory_results$fib)) {
  tf <- trajectory_results$fib
  p12 <- ggplot(tf, aes(UMAP1, UMAP2, color = pseudotime)) + geom_point(size = 0.35, alpha = 0.75) +
    scale_color_viridis_c(na.value = "grey85") + coord_equal() +
    labs(title = "Descriptive fibroblast state trajectory", color = "Pseudotime") + theme_void()
  save_plot(p12, "12_fibroblast_trajectory", 7, 6)
}
if (!is.null(trajectory_results$neut)) {
  tn <- trajectory_results$neut
  p13 <- ggplot(tn, aes(UMAP1, UMAP2, color = pseudotime)) + geom_point(size = 0.2, alpha = 0.65) +
    scale_color_viridis_c(na.value = "grey85") + coord_equal() +
    labs(title = "Descriptive colon neutrophil state trajectory", color = "Pseudotime") + theme_void()
  save_plot(p13, "13_neutrophil_trajectory", 7, 6)
}

summary_lines <- c(
  paste("Object:", ncol(obj), "cells and", nrow(obj), "genes"),
  paste("Colon fibroblasts:", length(fib_idx)),
  paste("Colon mature neutrophils:", length(neut_idx)),
  paste("Paired UC patients:", sum(paired_status$paired_uc)),
  paste("Biopsies in cross-compartment analysis after cell-count filters:", nrow(cross)),
  paste("Fibroblast DE tests:", paste(unique(fib_de$analysis), collapse = "; ")),
  paste("Neutrophil DE tests:", paste(unique(neut_de$analysis), collapse = "; ")),
  "Inference uses biopsy/patient-level aggregation; trajectories and mediation are exploratory.",
  "ITGA5+ITGB1 RNA co-detection is not proof of surface alpha5-beta1 protein."
)
writeLines(summary_lines, file.path(out_dir, "RUN_SUMMARY.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
saveRDS(list(fib_sample = fib_sample, neut_sample = neut_sample, cross = cross,
             fib_comp_test = fib_comp_test, neut_comp_test = neut_comp_test,
             program_tests = program_tests, cross_cor = cross_cor, delta_cor = delta_cor,
             clinical_cor = clinical_cor, mediation = mediation_res),
        file.path(obj_dir, "analysis_summaries.rds"))
message("Analysis complete: ", normalizePath(out_dir))
