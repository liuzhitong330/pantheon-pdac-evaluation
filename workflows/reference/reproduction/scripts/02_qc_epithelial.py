from pathlib import Path
import pandas as pd
import scanpy as sc

project = Path(__file__).resolve().parents[2]
a = sc.read_h5ad(project / "pdac_epithelial_counts.h5ad")

a.var["mt"] = a.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    a, qc_vars=["mt"], percent_top=None,
    log1p=False, inplace=True
)

keep = (
    (a.obs["n_genes_by_counts"] >= 200)
    & (a.obs["pct_counts_mt"] <= 25)
)

print("Cells before basic QC:", a.n_obs)
print("Cells removed:", int((~keep).sum()))
a = a[keep].copy()
print("Cells after basic QC:", a.n_obs)

counts = pd.crosstab(
    a.obs["sample_id"], a.obs["analysis_group"]
).reindex(columns=["Malignant", "Normal_epithelial"], fill_value=0)

min_cells = 10
paired_samples = counts.index[
    counts.ge(min_cells).all(axis=1)
].tolist()

print("\nSamples with at least 10 cells in each group:")
print(counts.loc[paired_samples].to_string())
print("\nNumber of paired samples:", len(paired_samples))

a.uns["paired_samples_min10"] = paired_samples
a.uns["basic_qc_rules"] = {
    "min_genes": 200,
    "max_pct_counts_mt": 25,
}

output = project / "pdac_epithelial_qc.h5ad"
a.write_h5ad(output)
print("Saved:", output)
