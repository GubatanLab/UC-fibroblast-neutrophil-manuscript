# Fibroblast regulation of neutrophils in ulcerative colitis

Analysis code, source tables and canonical figures for the September 12, 2026 Gastroenterology manuscript examining fibroblast niches, neutrophil phenotypes and α5β1-directed treatment in ulcerative colitis.

The study integrates human colonic single-cell transcriptomics, CODEX spatial protein profiling, fibroblast–neutrophil coculture, mouse colitis interventions and independent human validation datasets.

## Figures

- [Main Figures 1–7](figures/Main_Figures_1_to_7.pdf)
- [Extended Data Figures 1–10](figures/Extended_Data_Figures_1_to_10.pdf)
- [Individual main figures](figures/main)
- [Individual extended data figures](figures/extended_data)
- [Figure and analysis guide](docs/figure_guide.md)

The deposited PDFs are unchanged copies of the canonical manuscript artwork. Extended Data Figures 1–10 correspond to Supplementary Figures 1–10 in the Gastroenterology manuscript. Figure 1 contains panels a–f; the paired coexpression display is in Extended Data Figure 1c.

## Repository contents

| Location | Contents |
|---|---|
| `figures/` | Canonical PDFs, PNG previews and the two combined PDF sets |
| `source_data/` | Figure-specific tables, analysis summaries and selected biological-sample measurements |
| `analyses/source/` | Originating Python and R workflows, grouped by study component |
| `analyses/reproduce_summary_analyses.py` | Executable checks of key results using deposited tables |
| `config/` | Biological-mouse analysis inventories and epithelial gene programs |
| `environment/` | Package requirements and source-workflow dependencies |
| `metadata/` | Figure index, analysis index, source provenance and file checksums |
| `tools/` | Repository verification and PDF assembly |

## Quick start

Python 3.11 or later is recommended for the deposited-table workflow.

```bash
python -m venv .venv
# Activate the environment using the command appropriate for your operating system.
python -m pip install -r environment/requirements-summary.txt
python tools/validate_repository.py
python analyses/reproduce_summary_analyses.py
```

The executable workflow recomputes fibroblast paired-program tests, neutrophil state-program tests, coculture fluorescence comparisons, elastase comparisons, epithelial program effects and human receptor-expression contrasts. It compares its results with the deposited estimates and writes a report under `results/reproduction/`. The canonical figures and source tables are not overwritten.

To recreate the combined PDFs from the deposited individual figures:

```bash
python tools/assemble_figures.py
```

## Reproducibility scope

The deposited-table workflow is runnable from the repository. The source workflows document upstream processing and figure development, but require their specified input objects and configuration; they are not a single end-to-end pipeline. Some source figure scripts represent an earlier assembly stage, while the PDFs under `figures/` define the manuscript figure version. The figure guide identifies the relevant components.

Raw sequencing reads, full Seurat/AnnData objects, FCS files and whole-slide images are not distributed here. Public accession links and input requirements are in [Data access](docs/data_access.md). Biological-mouse comparisons use the deposited sample inventories. Existing neighborhood and communication estimates are supplied as fitted summaries; the source object, graph and full modeling inputs are required to refit them. Rerunning a model with a different input population does not reproduce a stored model estimate.

Treatment and acquisition batch are aligned in the mouse transcriptomic comparisons. Spatial proximity and predicted communication support candidate mechanisms; they do not establish direct contact, receptor specificity or clinical efficacy. RNA programs, protein fluorescence and elastase activity are distinct measurements. The [reproducibility notes](docs/reproducibility.md) describe these interpretation limits and the validation performed for this repository.

## Citation and reuse

See [CITATION.cff](CITATION.cff). The manuscript is unpublished and has no manuscript DOI assigned in this deposit. Related public datasets retain their original citations and access conditions. No software or data redistribution license is granted by this private repository; contact the study authors about reuse.
