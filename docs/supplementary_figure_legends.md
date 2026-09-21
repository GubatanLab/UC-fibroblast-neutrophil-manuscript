Supplementary Figure 1. Stromal programs, neutrophil annotation controls, and paired coexpression

(A) Eight fibroblast-associated programs in the broader stromal compartment, including smooth-muscle cells and pericytes; 14 paired UC patients. Lines connect samples, bars indicate medians, and open points mark <10 cells. Paired Wilcoxon tests are BH-adjusted across eight programs. (B) Pooled original program expression, z-scored across seven annotations within each program. Counts: erythrocytes, 2; macrophages, 7; CXCR4, 9,614; MX1, 5,222; OSM, 8,432; PADI4, 5,911; low RNA, 9,505. The non-neutrophil annotations are controls. This cell-weighted display retains naming genes and differs from Figure 1E. (C) Coexpression in nine patient pairs; paired Wilcoxon tests retain BH correction across five original routes. Displayed q values are 0.055, 0.055, and 0.076; none is significant at q<0.05. The NAMPT receptor summary uses the lower ITGA5/ITGB1 expression and does not establish receptor assembly or binding.

Supplementary Figure 2. CODEX annotations and patient composition

(A) UMAP of 245,020 cells across 12 broad annotations. (B) Eight-marker median fluorescence across five neutrophil phenotypes, z-scored within marker; 14,906 neutrophils. (C) Cell-type fractions across 24 independent patients: 6 control, 9 noninflamed UC, and 9 inflamed UC. These descriptive summaries support Figure 2; cells are not independent replicates.

Supplementary Figure 3. Spatial specificity and marker-threshold sensitivity

(A) Patient-median log₂ observed/expected enrichment of six immune lineages within 50 μm of four fibroblast categories; nulls retain coordinates and target abundance. (B) FAP-positive/α5β1-low proximity across 49 threshold combinations; outlines mark the 75th/75th percentiles. A,B include 6 control, 9 noninflamed, and 9 inflamed UC patients. (C) Adjusted coefficients and pointwise 95% CIs from 17 UC patients with ≥30 neutrophils; BH q values for covariate-adjusted models and 499 restricted permutations. Only CD16 passes both corrections. These reference definitions differ from Figure 2E.

Supplementary Figure 4. Fibroblast annotation and coculture transcript controls

(A) Fibroblast embedding after DecontX, principal-component analysis, Harmony, and Leiden clustering; 238 audit cells, including low-confidence immune-like controls. Communication analyses use 170 bona fide fibroblasts. (B) FAP, ITGA5, and ITGB1 across five fibroblast-containing conditions; size, detection; color, scaled expression. (C) Marker-program support for four annotations. Values use corrected log-normalized RNA. These descriptive controls have no hypothesis tests; sparse recovery limits condition-specific inference, and transcript codetection does not establish surface integrins.

Supplementary Figure 5. Neutrophil annotation and quality control

(A) Neutrophil-only scVI–Harmony embedding, including the low-signal audit group. (B) Eight subcluster annotations; size, marker detection; color, scaled expression (−2 to +2). (C) Cleanup inventory: 92,166 resolved neutrophils, 80,111 low-signal/unresolved cells excluded from state analyses, and 7,056 contaminant-, doublet-, or outlier-like cells removed. Counts describe different annotation stages and differ from protein-assay denominators.

Supplementary Figure 6. Recorded-donor neutrophil state contrasts

(A) Four-state centered-log-ratio effects across six contrasts; means and 95% parametric CIs using donor-specific variance. Six recorded donors contribute per comparison; no positive minimum-cell filter was applied. Filled/open symbols indicate BH q<0.05/≥0.05. The α5β1-by-context interaction is the inhibitor effect with UC fibroblasts minus its effect without fibroblasts; none of four state interactions passes correction. Biological matching remains unverified. These composition analyses are separate from Figure 3 protein and RNA-program tests.

Supplementary Figure 7. Neutrophil recovery and protein markers after FAP-lineage ablation

(A) Shared manifold occupancy: DSS, 1,074 cells from six mice; ablation, two cells from five mice. Gray denotes the reference; curves do not measure time. (B) State fractions among all recovered cells per mouse; mean±SEM; Wilcoxon rank-sum tests with BH correction across four states. (C–F) Osm-, Padi4-, Mx1-, and Cxcr4-positive neutrophil measurements; black, PBS; dark red, PBS+DSS; purple, ganciclovir+DSS. Source observations, summaries, and significance annotations are retained. Sparse residual cells preclude within-state inference. Treatment–batch confounding and assay-specific statistical provenance are described in Supplementary Methods.

Supplementary Figure 8. Neutrophil occupancy and transferred programs after ATN-161

(A) Condition occupancy and Slingshot curve rooted in the Cxcr2-associated state. (B) Fractions assigned to two endpoints. (C) Fractions assigned to transferred human neutrophil programs. Points indicate source units; boxes, median/interquartile range; whiskers, 1.5×interquartile range. DSS-versus-blockade Wilcoxon tests use BH correction across four programs. The inventory contains 2,147 cells and 28 units (10 control, 6 DSS, 12 blockade); Figure 5C excludes four control units upstream. (D) Program scores across pseudotime; 60 source-unit-by-bin means±SE. Treatment and batch are confounded; cross-sectional curves do not establish temporal transitions.

Supplementary Figure 9. Marker-positive neutrophils in human fibroblast coculture

(A–D) PADI4-, OSM-, MX1-, and CXCR4-positive neutrophil percentages after culture alone, with control fibroblasts, or with UC fibroblasts. Source observations, mean±SEM designation, and significance annotations are retained. Samples and gates are not assumed to match Figure 3B,C. Marker-positive percentages do not establish mutually exclusive RNA states or effector activity.

Supplementary Figure 10. Stromal detection and complementary neutrophil protein measurements after ATN-161

(A,B) Detectable α5β1-positive fractions among PDPN-positive and PDPN-positive/FAP-positive fibroblasts. (C–F) OSM-, PADI4-, CXCR4- and MX1-positive percentages within neutrophils. (G) MPO fluorescence within CXCR4-positive neutrophils; OSM- and PADI4-subset MPO appear in Figure 5B. Black, control; dark red, DSS+PBS; purple, DSS+ATN-161. Source observations, error bars and annotations are retained; MX1 has no added error bars. Control-versus-DSS annotations in A,B,G retain reconstructed-point Welch/BH tests, with the three-endpoint MPO correction family shared with Figure 5B. Other annotations retain source tests. Lower integrin detection may reflect receptor occupancy or epitope interference. Marker-positive percentages do not establish absolute subset expansion or effector function.

Supplementary Figure 11. Complete epithelial program gene expression in mature absorptive colonocytes

Stored pseudobulk expression of all 39 genes from the six epithelial programs in mature absorptive colonocytes. Columns represent six DSS, four FAP-ablation and twelve α5β1-blockade biological mice, each with ≥20 cells in this state. Genes are grouped by program; mice are ordered by condition and identifier. Colors show gene-wise z-scores across the 22 equally weighted mice, calculated from existing log₂ counts per million using the population SD and saturated at ±2.5. White separators distinguish programs and conditions. All genes were detected in at least three mice; †Tacstd2 was detected in only 4/22. The original normalization includes a prior count. This descriptive display does not imply gene-level significance, protein abundance, epithelial lineage conversion or functional repair. Treatment and acquisition batch are confounded.

Supplementary Figure 12. Cell-state concordance and predicted communication responses to stromal interventions

(A) Within-compartment proportion differences for 47 shared cell states; dashed line, equal differences; Spearman correlation describes concordance. These comparisons use six DSS, five ablation and twelve blockade biological mice. (B) Matched communication-priority contrasts for 691 routes; orange, fibroblast-to-neutrophil; blue, neutrophil-to-fibroblast. Of 502 jointly DSS-induced routes, 490 have lower inferred priority after both interventions. (C) Selected reciprocal-route differences; purple circles, ablation; teal triangles, blockade. Positive values indicate lower inferred priorities relative to DSS. Communication analyses retain their separate source eligibility and sparse-cell limitations, detailed in Supplementary Methods. Comparisons derive from separate experiments sharing a DSS reference; treatment and acquisition batch are confounded. Inferred priorities do not establish signaling activity or direct mechanisms.
