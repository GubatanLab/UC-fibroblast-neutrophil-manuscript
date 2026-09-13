from pathlib import Path
import anndata as ad
import h5py
import numpy as np

ROOT = Path(r"input_data/mouse\Figure 6 TAURUS External Validation")
for stem in ("fibperi_final", "myeloid_final"):
    path = ROOT / "data" / f"{stem}.h5ad"
    obj = ad.read_h5ad(path, backed="r")
    obs = obj.obs
    uc = obs.loc[obs["Disease"].eq("UC")]
    sample = uc.drop_duplicates("sample_id")
    print(f"\n### {stem}")
    print("cells", len(uc), "samples", sample['sample_id'].nunique(), "patients", sample['Patient'].nunique())
    print(sample.groupby(["Treatment", "Remission_status"], observed=True).size())
    print("inflammation by treatment/outcome")
    print(sample.groupby(["Treatment", "Remission_status", "Inflammation"], observed=True).size())
    print("sample IDs repeated across site?", int(sample.duplicated(["Patient", "Treatment", "Site"]).sum()))
    print("genes present:")
    genes = ["FAP","PDPN","THY1","ITGA5","ITGB1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CSF3","IL6","ICAM1","OSMR","NAMPT","OSM","IL1B","TNF","S100A9"]
    print({g: (g in obj.var_names) for g in genes})
    print("first var names", list(obj.var_names[:8]))
    print("uns keys", list(obj.uns.keys()))
    obj.file.close()
    with h5py.File(path, "r") as h:
        x = h["X"]
        print("X type", type(x).__name__, "attrs", dict(x.attrs))
        if hasattr(x, "keys"):
            print("X keys", list(x.keys()), "data range", float(np.min(x['data'][:1000000])), float(np.max(x['data'][:1000000])))
