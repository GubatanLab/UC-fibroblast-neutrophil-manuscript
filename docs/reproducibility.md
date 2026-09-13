# Reproducibility and interpretation

## Executable result checks

`analyses/reproduce_summary_analyses.py` checks the following against deposited results:

1. Figure 1d: paired Wilcoxon fibroblast-program tests with continuity correction and BH adjustment across eight programs.
2. Figure 1e: within-patient Friedman tests across four neutrophil states and BH adjustment across nine programs.
3. Figure 3c: exact unpaired permutation tests of sample-median compensated fluorescence, with BH adjustment across 45 comparisons.
4. Figure 3d: exact unpaired elastase comparisons and Holm adjustment across three comparisons.
5. Figure 6d–e: epithelial program mean differences, Mann–Whitney tests and scope-specific BH adjustment.
6. Figure 7c: patient-average receptor-component contrasts with exact sign-flip tests and BH adjustment across six features.

Reproduction checks compare point estimates and relevant P/q values. They do not refit expression integration, recalibrate images, rerun flow gating, or independently reconstruct bootstrap confidence intervals. Source workflows contain the corresponding upstream procedures and uncertainty calculations. Validation results are written to `results/reproduction/verification.json`.

## Analytical units

Human tissue comparisons use patients or matched patient/biopsy summaries as specified by the analysis. Coculture source tables retain recorded donor and condition labels. Compensated fluorescence is a sample-median measurement; individual events are not treated as independent biological replicates. Biological-mouse comparisons use the deposited inventory and assay-specific eligibility thresholds. The core comparison contains six shared DSS mice, five ablation mice and twelve blockade mice; epithelial whole-compartment scoring includes four eligible ablation mice.

Overlapping Milo neighborhoods are not independent animals. Existing neighborhood-model and communication summaries are supplied as fixed source estimates. Their graph construction, eligibility and input population are properties of the original fit. The deposited biological-mouse refit procedures should be treated as new fits when their inputs differ; agreement with a stored model is not assumed. Very low residual neutrophil recovery constrains the ablation signaling analysis.

## Measurement and design limits

Mouse transcriptomic treatment groups are aligned with acquisition batch, so comparisons describe intervention-associated patterns. State fractions are relative representation, not absolute abundance. RNA state labels do not define protein-positive gates. Spatial proximity and ligand–receptor predictions are associative. Receptor-component RNA does not establish an assembled receptor or the cellular target of ATN-161.

RNA programs, compensated fluorescence and elastase activity are separate readouts. Some flow-channel identities and source replicate/normalization metadata remain provisional. Several elastase and historical flow values were transcribed from supplied displays; their provenance columns are retained. Elastase activity alone does not establish NET structure or formation mechanism. Original image calibration and assay records are needed for acquisition-level validation.

## Figure integrity

All deposited canonical figure PDFs are byte-identical to the manuscript workspace copies. The repository validation script checks stored hashes, all 17 individual figures, the 7-page main bundle and the 10-page extended data bundle. Rebuilding bundles writes separate outputs under `results/figure_bundles/` and does not alter the canonical files.

The originating source scripts have been parsed for syntax where a local interpreter is available. They have not all been executed end to end in a fresh environment; raw inputs and some model resources are external. Package requirements for the executable deposited-table workflow are separate from the larger source-workflow dependency list.
