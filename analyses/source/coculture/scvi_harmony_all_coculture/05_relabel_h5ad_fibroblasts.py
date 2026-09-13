#!/usr/bin/env python

import json
from pathlib import Path

import scvi  # noqa: F401; import first in this Windows environment
import anndata as ad

analysis_dir = Path(r"input_data/mouse\Neutrophil_Fibroblast_Coculture_Blockade_Figures\scVI_Harmony_All_Coculture")
files = [
    analysis_dir / "input" / "all_coculture_neutrophil_fibroblast_scvi_input.h5ad",
    analysis_dir / "objects" / "all_coculture_scVI_Harmony_latent.h5ad",
]

for path in files:
    x = ad.read_h5ad(path)
    labels = x.obs["clean_cell_class"].astype(str).replace({"Fibroblast_candidate": "Fibroblast"})
    x.obs["clean_cell_class"] = labels.astype("category")
    x.write_h5ad(path, compression="gzip")

metadata_path = analysis_dir / "source_data" / "scVI_Harmony_run_metadata.json"
metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
metadata["fibroblasts"] = metadata.pop("fibroblast_candidates", metadata.get("fibroblasts", 238))
metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
print("H5AD and run metadata relabeling complete")

