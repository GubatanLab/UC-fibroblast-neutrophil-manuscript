# Data access

## Public human datasets

| Dataset | Access | Use in this manuscript |
|---|---|---|
| Published UC single-cell atlas | [Broad Single Cell Portal SCP3755](https://singlecell.broadinstitute.org/single_cell/study/SCP3755) and [Frontiers in Immunology publication](https://doi.org/10.3389/fimmu.2026.1705328) | Discovery fibroblast and neutrophil states, paired programs and candidate communication |
| TAURUS IBD atlas | [Version 3 dataset](https://doi.org/10.5281/zenodo.14007626) and [Nature Immunology publication](https://doi.org/10.1038/s41590-024-01994-8) | Independent UC fibroblast receptor components and paired tissue/longitudinal comparisons |
| UC Xenium spatial dataset | [Broad Single Cell Portal SCP3818](https://singlecell.broadinstitute.org/single_cell/study/SCP3818) | Two spatial sections from one human UC donor |

Follow each dataset's original access terms. Authentication or an access agreement may be required at the data provider. No credentials or access tokens are included here.

## Study inputs supplied separately

- Human discovery: annotated expression objects, sample/condition mapping, and normalized expression used for module scoring.
- CODEX: registered multichannel images, cell coordinates, protein intensities, annotations and patient mapping.
- Coculture: compensated FCS measurements, gating definitions, sample/condition mapping, annotated RNA objects and elastase assay records.
- Mouse studies: annotated compartment objects, biological-mouse inventories, original histology images, flow measurements and model resources.
- Communication models: the exact species-specific ligand–receptor network and ligand–target matrix used by the source workflow, together with compatible MultiNicheNet/NicheNet versions.

Selected derived tables are deposited under `source_data/`; the figure guide and source provenance table map them to analyses. These tables do not replace the original acquisition records.

## Local configuration

Source scripts use relative placeholders such as `input_data/human`, `input_data/codex`, `input_data/flow`, and `input_data/mouse`. Configure the input and output variables near the beginning of a source script before running it. Paths inferred from the original working-directory structure also require adjustment. The source archive preserves the analysis structure rather than providing a relocated raw-data mirror.

Use a separate working directory for upstream reruns. The executable deposited-table workflow and PDF assembly commands are already configured to read repository paths and write to `results/`.
