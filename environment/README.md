# Source workflow environments

The smaller executable deposited-table workflow uses `requirements-summary.txt`. Its tested package versions are recorded in `requirements-summary-tested.txt`.

Upstream scripts span several analysis environments. Install dependencies only for the workflow you intend to rerun; there is no validated single environment for every source script. R source syntax was checked with R 4.5.3. MultiNicheNet outputs in the manuscript use version 2.1.0. The communication workflows also require compatible NicheNet resources; Seurat, scVI/Harmony and trajectory workflows require their respective input formats and model environments.

R package names detected in source workflows:

BiocNeighbors, BiocParallel, CellChat, Giotto, Matrix, RANN, S4Vectors, Seurat, SeuratObject, SingleCellExperiment, SummarizedExperiment, anndataR, condiments, data.table, dbscan, decoupleR, dplyr, edgeR, ggplot2, ggplotify, ggraph, ggrastr, ggrepel, grid, harmony, igraph, limma, magick, miloR, msigdbr, multinichenetr, nichenetr, patchwork, pheatmap, progeny, scales, slingshot, splines, tibble, tidyr, tradeSeq, uwot.

Python import names detected in source workflows (these include standard-library and local modules, so this is an inventory, not a pip installation command):

PIL, __future__, anndata, collections, concurrent, copy, csv, flowkit, h5py, harmonypy, hashlib, html, itertools, json, lxml, math, matplotlib, numpy, pandas, pathlib, platform, pymupdf, pypdf, re, scipy, scvi, seaborn, shutil, sklearn, statsmodels, sys, tempfile, textwrap, time, torch, urllib, warnings, xml, zipfile.

Arial is used in the original figure layouts. Install it under an appropriate license or select an available font when rerendering. Fonts are not redistributed. Canonical PDF files already contain their exported text appearance.
