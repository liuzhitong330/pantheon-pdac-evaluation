from pathlib import Path
import pandas as pd
import scanpy as sc
import anndata as ad

project = Path(__file__).resolve().parents[2]
data = project / "data"

ann = pd.read_csv(data / "GSE242230_annotations.txt.gz", sep="\t")
ann["analysis_group"] = ann["cell_type_specific"].map({
    "Malignant - Basal": "Malignant",
    "Malignant - Classical": "Malignant",
    "Normal Epithelial": "Normal_epithelial",
})
selected = ann.loc[ann["analysis_group"].notna()].set_index("cell_id")

parts = []
for sample, meta in selected.groupby("sample_id", sort=True):
    print(f"Reading {sample}...", flush=True)

    files = list(data.glob(f"*_{sample}_matrix.mtx.gz"))
    if len(files) != 1:
        raise ValueError(f"{sample}: expected one matrix, found {len(files)}")

    prefix = files[0].name.removesuffix("matrix.mtx.gz")
    a = sc.read_10x_mtx(
        data, prefix=prefix, var_names="gene_symbols", make_unique=True
    )
    a.obs_names = sample + "_" + a.obs_names

    missing = meta.index.difference(a.obs_names)
    if len(missing):
        raise ValueError(f"{sample}: {len(missing)} cell IDs did not match")

    a = a[meta.index].copy()
    a.obs = meta.copy()

    if parts and not a.var["gene_ids"].equals(parts[0].var["gene_ids"]):
        raise ValueError(f"{sample}: gene information differs; check before merging")

    parts.append(a)

combined = ad.concat(parts, join="inner", merge="same")
assert combined.n_obs == len(selected)
assert combined.obs_names.is_unique

combined.uns["source_accession"] = "GSE242230"
output = project / "pdac_epithelial_counts.h5ad"
combined.write_h5ad(output)

print("\nCombined shape (cells, genes):", combined.shape)
print(combined.obs["analysis_group"].value_counts().to_string())
print("Samples:", combined.obs["sample_id"].nunique())
print("Saved:", output)
