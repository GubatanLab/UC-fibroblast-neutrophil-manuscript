# Methods and interpretation of deposited data

## Study components and units

Human discovery displays contain 5,187 fibroblasts and 38,684 neutrophils. Paired state analyses use 13 fibroblast and 17 neutrophil biopsy pairs, and receptor-expression analyses use 10 biopsies. CODEX comprises 24 different patients: 6 controls, 9 noninflamed UC and 9 inflamed UC. Eligibility differs by endpoint; the analytical-populations table records the manuscript denominators. Patient/sample summaries are the inferential units where specified. Cells, events and overlapping neighborhoods are not independent biological replicates.

## Ethics approvals

Collection of human blood and biopsies for coculture, single-cell RNA sequencing and CODEX analyses was approved by Stanford University (IRB 28437, 60958 and 52317). Mouse procedures were approved by the Stanford University IACUC (22000 and 27715) and Mayo Clinic IACUC (A00008640-26). These identifiers reproduce the author-supplied final Methods. No consent or waiver terms are inferred from an approval number.

## Human transcriptomic programs

The discovery atlas is available at SCP3755 and described by Eshghi, Gubatan, Mazrooei and colleagues (Front Immunol 2026;17:1705328; https://doi.org/10.3389/fimmu.2026.1705328). Figure 1D uses paired stromal-program summaries with Wilcoxon tests and BH correction. Figure 1E uses matched-patient neutrophil programs excluding naming markers, with Friedman tests and BH correction across nine programs. Gene membership, patient-state summaries and threshold sensitivity are retained in source_data/figure_01/. Original cell-weighted summaries and intermediate tables are separately identifiable by filename; they are not interchangeable with current patient-weighted panels.

## CODEX

The staining panel comprised 53 antibodies plus DAPI. The current reagent table records clone, catalog, manufacturer, barcode, reporter, acquisition cycle, exposure and dilution. Substance P was not used in the staining panel. The α5β1-specific antibody was clone M200 (Volociximab), rabbit IgG chimeric, Abcam ab275977. CD4 is ATTO550 cycle 4. The conflicting CD274 acquisition cycle is unspecified rather than inferred.

The recorded Seurat object contains 245,020 cells. Cell-level fluorescence exported from QuPath was analyzed using Seurat 5.5.0 and SeuratObject 5.4.0 in R 4.5.3, CLR normalization with margin=2, scaling, 30-component PCA and Harmony 2.0.5 adjustment by PatientID. SNN/Louvain clustering used Harmony dimensions 1–10 at resolution 1.0; the cohort UMAP used dimensions 1–30. Lineage-protein combinations and selected marker relationships from the UC atlas guided biological annotation. Neutrophil identity used CD66b/CD15/CD16/CD11b; OSM, CXCR4, PADI4 and MX1 informed phenotypes. Existing reference labels were matched by cell identifier for 226,659 cells; 18,361 P07 cells were assigned by 25-nearest-neighbor majority vote in Harmony dimensions 1–10. This is CODEX-to-CODEX label transfer, not direct RNA-to-protein transfer.

The saved normalization/scaling/PCA retained 54 feature columns, including a legacy Substance P column. Its presence does not establish measured Substance P staining. The object was not reprocessed during this deposit preparation. CD31 is present in the clustering object but absent from the downstream proximity export. Neutrophil marker summaries used exported nuclear mean-fluorescence fields; these are operational intensity readouts, not proof of subcellular localization.

Current Figure 2A–C uses 50-µm cellular neighborhoods and patient-level summaries. Figure 2D uses upper-quartile FAP/α5β1 categories and availability-adjusted proximity. Figure 2E uses a separate posterior gate, adjusted marker-proximity coefficients and eight eligible inflamed UC patients. Files are separated into neighborhoods/, quartile_proximity/ and posterior_gate/. Older codex_high_yield and earlier_figure_sources tables include distinct sensitivity definitions; use the figure map to avoid substituting them for current panels. Spatial associations do not establish direct contact or signaling.

## Coculture measurements

Control-blood neutrophils were cocultured with control or UC fibroblasts for 16 hours. ATN-161 denotes the α5β1-directed treatment. Sample-median compensated fluorescence, marker-positive fractions, NET-associated elastase activity and RNA-program scores are distinct endpoints. Figure 3C uses exact unpaired permutation tests with BH correction across 45 comparisons. Figure 3D includes six independent experiments per arm, with a common scale and normalization confirmed by the investigator, and exact unpaired comparisons with Holm correction. The activity measurement follows removal of free elastase and nuclease digestion; it does not by itself establish NET structure or mechanism. Some source values originated from supplied plots, with provenance retained; author verification of plotted values does not convert them into original instrument exports. RNA-program t tests and exact sign-flip sensitivity checks are separate analyses.

## Mouse and epithelial readouts

Biological-mouse source tables include a core DSS/ablation/blockade comparison of 6/5/12 mice; whole-epithelium RNA scoring requires at least 20 cells and uses 6/4/12 eligible mice. Retained inventories define the relevant population for each comparison. FAP-lineage ablation and ATN-161 experiments share a DSS reference but have separate experimental origins. Treatment and acquisition batch are aligned in the transcriptomic datasets; observed RNA patterns are intervention-associated. Existing neighborhood and communication estimates are supplied as stored model summaries. Their original graphs, fitting populations and model resources are not distributed here, so these summaries cannot be refit from the deposited biological-sample tables alone.

Figure 6A,B retains existing epithelial RNA-program estimates; counts were aggregated by biological mouse, TMM-normalized and summarized as log2 CPM. Figure 6C,D describes recovered-cell proportions, not absolute abundance. Figure 6E and Supplementary Figure 11 provide selected or full 39-gene mature absorptive-colonocyte displays; gene-wise z scores derive from the stored log2 CPM data. Figure 6F provides individual-mouse values for existing chemokine, absorptive and mucus/secretory programs. RNA programs and marker transcripts do not establish epithelial function or repair. Figure 5 display-reconstruction tables are explicitly approximate source-image measurements; they are not original flow-cytometry events.

## Independent human evidence and statistical interpretation

TAURUS and SCP3818 provide independent human expression/spatial evidence. Transcript detection does not establish assembled receptors; two SCP3818 sections represent one donor. Ligand–receptor and ligand–target scores identify candidate networks, not measured signaling flux. Tests use endpoint-specific multiple-testing families and assay-specific denominators. Source units, contrast directions and missing-value codes are retained. File names containing pseudotime refer to computational ordering of cells and do not denote longitudinal animal sampling.

The included verification script checks selected point estimates and P/q values. It does not create new manuscript analyses or validate all mechanistic interpretations. Exact source-data transformations predating this package are recorded in metadata/source_provenance.csv where available.

## Verified bibliography and final figure labels

The final manuscript reference audit covers 36 main and 18 supplementary entries, including 39 distinct verified cited DOIs. The two Broad datasets remain accession citations because dataset-specific DOIs were not identified. The manuscript archive DOI was checked separately. Full citation details and verification sources are in `metadata/reference_verification_2026-09-21.csv`. MultiNicheNet is cited as a preprint (https://doi.org/10.1101/2023.06.13.544751), alongside software version 2.1.0.

Visible panel letters use uppercase throughout the corrected figures. Figure 3D is labeled NET-associated elastase activity; Figure 7C describes receptor components rather than measured receptor availability. Typography and label changes do not alter plotted measurements or statistical results. Histology scale calibration and field identities and the MX1/PADI4 flow-channel assignments still require their original acquisition records; provisional qualifications are retained in the legends.
