# Fibroblast regulation of neutrophils in ulcerative colitis

Analysis code, derived source tables and canonical figures accompanying **Single-cell profiling reveals stromal niche and α5β1 blockade regulation of neutrophil phenotypes in ulcerative colitis**. This release updates the repository to the September 21, 2026 manuscript and figure set.

The study integrates human colonic single-cell transcriptomics, CODEX spatial protein profiling, fibroblast–neutrophil coculture, DSS colitis experiments and independent human datasets. The epithelial analyses describe shared reductions in chemokine programs and distinct barrier-associated and regenerative expression patterns after FAP-lineage ablation and α5β1-directed treatment. Protein, NET-associated elastase and RNA endpoints are distinguished throughout.

## Figures and tables

- [Main Figures 1–7](figures/Main_Figures_1_to_7.pdf)
- [Supplementary Figures 1–12](figures/Supplementary_Figures_1_to_12.pdf)
- [Individual main figures](figures/main) and [supplementary figures](figures/supplementary)
- [Supplementary Table 1: CODEX markers](tables/Supplementary_Table_1_CODEX_Markers.csv)
- [Supplementary Table 2: Analytical populations](tables/Supplementary_Table_2_Analytical_Populations.csv)
- [Figure and analysis guide](docs/figure_guide.md), [panel-to-data map](docs/FIGURE_DATA_MAP.csv) and [main figure legends](docs/main_figure_legends.md)

Figure 3 includes the updated coculture schematic. Figure 5 includes OSM/PADI4-subset MPO in panel B. Figure 6 presents six epithelial panels, with the full gene heatmap and supporting concordance/communication displays in Supplementary Figures 11–12. Supplementary Figures 1–10 also retain matching Extended Data aliases for existing links.

## Repository contents

| Location | Contents |
|---|---|
| `source_data/`, `tables/` | 347 derived data tables, including current spatial, epithelial and CODEX records |
| `figures/` | Canonical figure PDFs, PNG previews and combined PDFs |
| `analyses/source/` | Originating Python and R workflows grouped by study component |
| `analyses/reproduce_summary_analyses.py` | Executable checks of selected reported results |
| `config/`, `environment/` | Analysis inventories, gene programs and dependencies |
| `docs/` | Methods, interpretation, data access and figure/data mapping |
| `metadata/` | Column dictionary, provenance, figure/analysis indices and checksums |
| `tools/` | Repository verification and PDF assembly |

## Quick start

Python 3.11 or later is recommended for the deposited-table workflow.

```bash
python -m venv .venv
# Activate the environment for your operating system.
python -m pip install -r environment/requirements-summary.txt
python tools/validate_repository.py
python analyses/reproduce_summary_analyses.py
```

The reproduction workflow recomputes selected fibroblast-program, neutrophil-state, coculture fluorescence, elastase, epithelial-program and human receptor-expression statistics from deposited tables. It writes verification results under `results/reproduction/`; canonical source data and figures are not overwritten. The release passed 17 selected statistical checks. To assemble additional figure bundles from the individual pages, run `python tools/assemble_figures.py`.

## Reproducibility scope

The deposit contains selected derived data and reproducible summary-statistic checks. Upstream source workflows require their specified external inputs and configuration and do not constitute a single end-to-end pipeline. The current epithelial figure assembly is in `analyses/source/figure_assembly/figure6_current.py` (Arial and the documented plotting dependencies are required). Historical assembly scripts may use earlier panel labels; the current PDFs and panel-to-data map define the release. Raw reads, full Seurat/AnnData objects, FCS files, whole-slide images and original model-fitting resources are held separately. See [Methods and interpretation](docs/METHODS_AND_INTERPRETATION.md), [Data access](docs/DATA_ACCESS.md) and [Reproducibility notes](docs/reproducibility.md).

Biological-sample inventories and assay-specific denominators define the analysis units. Stored neighborhood and communication summaries cannot be refit from the deposited biological-sample tables alone. Mouse transcriptomic treatment groups are aligned with acquisition batch. Spatial proximity and predicted communication identify candidate mechanisms; RNA programs do not establish epithelial function or repair. No new experiments or inferential analyses were generated for this release.

## Citation and reuse

See [CITATION.cff](CITATION.cff). Release DOI: [10.5281/zenodo.22881923](https://doi.org/10.5281/zenodo.22881923). The manuscript itself remains unpublished. The original code is licensed under [MIT](LICENSE-CODE), and original data tables and figures under [CC BY 4.0](LICENSE-DATA.md). Third-party materials retain their original terms. See [Licensing](LICENSE.md).
