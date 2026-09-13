#!/usr/bin/env python

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Import scvi first in this Windows environment.
import scvi
import anndata as ad
import harmonypy as hm
import torch


SEED = 20260816
ANALYSIS_DIR = Path(
    r"input_data/mouse\Neutrophil_Fibroblast_Coculture_Blockade_Figures\scVI_Harmony_All_Coculture"
)
INPUT_FILE = ANALYSIS_DIR / "input" / "all_coculture_neutrophil_fibroblast_scvi_input.h5ad"
MODEL_DIR = ANALYSIS_DIR / "model" / "scvi_30d"
OBJECT_DIR = ANALYSIS_DIR / "objects"
SOURCE_DIR = ANALYSIS_DIR / "source_data"
LOG_DIR = ANALYSIS_DIR / "logs"
for directory in (MODEL_DIR.parent, OBJECT_DIR, SOURCE_DIR, LOG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

scvi.settings.seed = SEED
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.set_num_threads(max(1, min(12, (torch.get_num_threads() or 1))))


def log(message: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"{stamp} | {message}"
    print(line, flush=True)
    with (LOG_DIR / "scvi_harmony_training.log").open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


log(f"Python {sys.version.split()[0]}; scvi-tools {scvi.__version__}; torch {torch.__version__}")
log(f"Loading {INPUT_FILE}")
adata = ad.read_h5ad(INPUT_FILE)
adata.var_names_make_unique()
if "counts" not in adata.layers:
    raise RuntimeError("H5AD input does not contain the counts layer")
if not sp.issparse(adata.layers["counts"]):
    adata.layers["counts"] = sp.csr_matrix(adata.layers["counts"])
else:
    adata.layers["counts"] = adata.layers["counts"].tocsr()
count_matrix = adata.layers["counts"]

if adata.n_obs != 179571:
    raise RuntimeError(f"Unexpected cell count: {adata.n_obs}")
if set(adata.obs["clean_cell_class"].astype(str).unique()) != {"Neutrophil", "Fibroblast"}:
    raise RuntimeError("Unexpected cell classes")
if np.min(count_matrix.data) < 0 or not np.allclose(count_matrix.data, np.round(count_matrix.data)):
    raise RuntimeError("scVI input is not non-negative integer count data")

for key in ["orig.ident", "SampleID", "DonorID", "Condition", "ConditionCode", "CultureStage", "dataset_source", "clean_cell_class"]:
    adata.obs[key] = adata.obs[key].astype("category")

log(f"Input verified: {adata.n_obs} cells, {adata.n_vars} features, {count_matrix.nnz} nonzero counts")
scvi.model.SCVI.setup_anndata(adata, layer="counts")
model = scvi.model.SCVI(
    adata,
    n_hidden=128,
    n_latent=30,
    n_layers=2,
    dropout_rate=0.10,
    dispersion="gene",
    gene_likelihood="nb",
    latent_distribution="normal",
)

log("Training scVI on CPU: 30 latent dimensions, maximum 80 epochs, early stopping patience 10")
train_start = time.time()
model.train(
    max_epochs=80,
    accelerator="cpu",
    devices=1,
    batch_size=512,
    train_size=0.90,
    validation_size=0.10,
    early_stopping=True,
    early_stopping_patience=10,
    check_val_every_n_epoch=1,
    enable_progress_bar=True,
)
training_seconds = time.time() - train_start
log(f"scVI training completed in {training_seconds / 60:.1f} minutes")
model.save(MODEL_DIR, overwrite=True, save_anndata=False)

history = model.history
if isinstance(history, dict):
    history_frames = []
    for metric, values in history.items():
        frame = pd.DataFrame(values).reset_index().rename(columns={"index": "epoch"})
        value_cols = [c for c in frame.columns if c != "epoch"]
        if value_cols:
            frame = frame.rename(columns={value_cols[0]: "value"})
        frame["metric"] = metric
        history_frames.append(frame[["epoch", "metric", "value"]])
    if history_frames:
        pd.concat(history_frames, ignore_index=True).to_csv(SOURCE_DIR / "scvi_training_history.csv", index=False)
else:
    pd.DataFrame(history).to_csv(SOURCE_DIR / "scvi_training_history.csv")

log("Extracting 30-dimensional scVI latent representation")
latent = np.asarray(model.get_latent_representation(), dtype=np.float32)
if latent.shape != (adata.n_obs, 30):
    raise RuntimeError(f"Unexpected scVI latent shape: {latent.shape}")
adata.obsm["X_scVI"] = latent

log("Running Harmony on scVI latent space using DonorID only; condition is preserved")
harmony_start = time.time()
harmony = hm.run_harmony(
    latent,
    adata.obs,
    vars_use=["DonorID"],
    theta=2.0,
    lamb=1.0,
    max_iter_harmony=30,
    max_iter_kmeans=20,
    epsilon_harmony=1e-4,
    epsilon_cluster=1e-5,
    random_state=SEED,
    verbose=True,
)
harmony_latent = np.asarray(harmony.Z_corr, dtype=np.float32)
if harmony_latent.shape == (30, adata.n_obs):
    harmony_latent = harmony_latent.T
if harmony_latent.shape != latent.shape:
    raise RuntimeError(f"Unexpected Harmony shape: {harmony_latent.shape}")
adata.obsm["X_scVI_Harmony"] = harmony_latent
log(f"Harmony completed in {(time.time() - harmony_start) / 60:.1f} minutes")

log("Writing scVI and scVI-Harmony embeddings for graph clustering and UMAP in Seurat")
cell_output = pd.DataFrame(
    {
        "cell": adata.obs_names.astype(str),
        "orig.ident": adata.obs["orig.ident"].astype(str).to_numpy(),
        "SampleID": adata.obs["SampleID"].astype(str).to_numpy(),
        "DonorID": adata.obs["DonorID"].astype(str).to_numpy(),
        "Condition": adata.obs["Condition"].astype(str).to_numpy(),
        "ConditionCode": adata.obs["ConditionCode"].astype(str).to_numpy(),
        "CultureStage": adata.obs["CultureStage"].astype(str).to_numpy(),
        "dataset_source": adata.obs["dataset_source"].astype(str).to_numpy(),
        "clean_cell_class": adata.obs["clean_cell_class"].astype(str).to_numpy(),
        "clean_state": adata.obs["clean_state"].astype(str).to_numpy(),
        "clean_annotation": adata.obs["clean_annotation"].astype(str).to_numpy(),
    }
)
for idx in range(30):
    cell_output[f"scVI_{idx + 1}"] = latent[:, idx]
for idx in range(30):
    cell_output[f"scVI_Harmony_{idx + 1}"] = harmony_latent[:, idx]
cell_output.to_csv(SOURCE_DIR / "scVI_Harmony_latent_and_metadata.csv.gz", index=False, compression="gzip")

run_metadata = {
    "cells": int(adata.n_obs),
    "features": int(adata.n_vars),
    "neutrophils": int((adata.obs["clean_cell_class"].astype(str) == "Neutrophil").sum()),
    "fibroblasts": int((adata.obs["clean_cell_class"].astype(str) == "Fibroblast").sum()),
    "preculture_neutrophils": int((adata.obs["ConditionCode"].astype(str) == "PRE").sum()),
    "scvi_latent_dimensions": 30,
    "scvi_max_epochs": 80,
    "scvi_training_minutes": training_seconds / 60,
    "scvi_batch_covariate": None,
    "harmony_covariate": "DonorID",
    "harmony_condition_corrected": False,
    "downstream_graph_runtime": "Seurat FindNeighbors/FindClusters and RunUMAP",
    "planned_neighbors": 30,
    "planned_primary_resolution": 0.4,
    "random_seed": SEED,
    "python": sys.version,
    "platform": platform.platform(),
    "scvi_tools": scvi.__version__,
    "torch": torch.__version__,
    "cuda_available": bool(torch.cuda.is_available()),
}
with (SOURCE_DIR / "scVI_Harmony_run_metadata.json").open("w", encoding="utf-8") as handle:
    json.dump(run_metadata, handle, indent=2)

adata.write_h5ad(OBJECT_DIR / "all_coculture_scVI_Harmony_latent.h5ad", compression="gzip")
log("scVI and Harmony stages completed successfully")
