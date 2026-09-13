from pathlib import Path
import anndata as ad
import pandas as pd

ROOT = Path(r"input_data/mouse\Figure 6 TAURUS External Validation")

for stem in ("fibperi_final", "myeloid_final"):
    path = ROOT / "data" / f"{stem}.h5ad"
    obj = ad.read_h5ad(path, backed="r")
    print(f"\n### {stem}: {obj.shape}")
    print("obs columns:", list(obj.obs.columns))
    print("var columns:", list(obj.var.columns))
    print("layers:", list(obj.layers.keys()))
    print("raw:", None if obj.raw is None else obj.raw.shape)
    for col in obj.obs.columns:
        s = obj.obs[col]
        n = s.nunique(dropna=False)
        if n <= 50:
            vc = s.value_counts(dropna=False).head(50)
            print(f"\n[{col}] n={n}\n{vc.to_string()}")
    obj.file.close()

pair = pd.read_csv(ROOT / "data" / "paired_sample_list.csv")
print("\n### paired_sample_list")
print(pair.columns.tolist())
print(pair.to_string(index=False))
