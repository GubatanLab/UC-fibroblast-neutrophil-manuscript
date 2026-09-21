# Figure and analysis guide

## Main figures

| Figure | Analysis | Source workflows | Principal deposited tables |
|---|---|---|---|
| 1 | Human fibroblast and neutrophil states, paired fibroblast programs, matched-patient neutrophil programs and predicted signals | `analyses/source/human`, `analyses/source/figure_01` | `source_data/figure_01`, `source_data/earlier_figure_sources/human_discovery` |
| 2 | CODEX neighborhoods and neutrophil proximity to FAP-high/α5β1-high fibroblasts | `analyses/source/spatial`, `analyses/source/figure_02` | `source_data/figure_02`, `source_data/codex_high_yield`, `source_data/earlier_figure_sources/figure2` |
| 3 | Coculture fluorescence, elastase activity, RNA programs and candidate ligand–module associations | `analyses/source/coculture`, `analyses/source/flow`, `analyses/source/coculture_support` | `source_data/figure_03`, `source_data/coculture_fluorescence`, `source_data/coculture_support` |
| 4 | FAP-lineage ablation, tissue inflammation, neighborhood estimates and candidate signaling | `analyses/source/mouse`, `analyses/source/figure_assembly` | `source_data/figure_04`, `source_data/earlier_figure_sources/figure4`, `source_data/mouse_comparison` |
| 5 | ATN-161, tissue inflammation, subset MPO, neutrophil neighborhoods and predicted signaling | `analyses/source/mouse`, `analyses/source/figure_assembly` | `source_data/figure_05`, `source_data/earlier_figure_sources/figure5`, `source_data/mouse_comparison` |
| 6 | Epithelial programs, cell composition, selected genes and individual-mouse profiles | `analyses/source/mouse`, `analyses/source/epithelial_and_validation` | `source_data/figure_06`, `source_data/epithelial_programs`, `source_data/mouse_comparison` |
| 7 | Independent human fibroblast receptor components, longitudinal observations and spatial validation | `analyses/source/external_validation` | `source_data/figure_07`, `source_data/taurus`, `source_data/taurus_receptors`, `source_data/spatial_validation` |

## Supplementary figures

| Figure | Content | Relevant study component |
|---|---|---|
| 1 | Stromal programs, neutrophil annotation controls and paired coexpression | Human discovery; `source_data/extended_data/ED1_revision_2026-09-12` |
| 2 | CODEX annotations and patient composition | Spatial profiling |
| 3 | Spatial specificity and marker-threshold sensitivity | Spatial profiling |
| 4 | Fibroblast annotation and coculture transcript controls | Coculture and annotation tables |
| 5 | Neutrophil annotation and quality control | Coculture and annotation tables |
| 6 | Recorded-donor neutrophil state contrasts | Coculture RNA analyses |
| 7 | Neutrophil recovery and protein markers after FAP-lineage ablation | Mouse and flow analyses |
| 8 | Neutrophil occupancy and transferred programs after ATN-161 | Mouse trajectory and neighborhood analyses |
| 9 | Marker-positive neutrophils in human fibroblast coculture | Human flow analyses |
| 10 | Stromal detection and remaining neutrophil proteins after ATN-161; OSM/PADI4-subset MPO is in main 5B | Mouse flow analyses |
| 11 | Full 39-gene epithelial heatmap | Epithelial expression data |
| 12 | State-level concordance and predicted reciprocal communication | Mouse comparative analyses |

Assembly workflows are in `analyses/source/extended_data`. Earlier figure-source directories preserve intermediate table names, which may use the panel numbering from that analysis stage. Use this guide and the current figure PDFs to establish the manuscript numbering; do not infer it from an older table filename alone.

`metadata/figure_index.csv` lists every deposited figure. `metadata/analysis_index.csv` lists each source script, and `metadata/source_provenance.csv` records the originating analysis location and checksum. These indices separate deposited evidence from the executable result checks described in the README.

See `docs/FIGURE_DATA_MAP.csv` for the current panel-level mapping. Supplementary Figures 1–10 retain matching Extended Data aliases.
