from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

folder = Path(__file__).resolve().parents[2]
a = ad.read_h5ad(folder / "pdac_epithelial_qc.h5ad")

# 1. 读取原始计数
if "counts" in a.layers:
    X = a.layers["counts"]
    source = "layers['counts']"
else:
    if "log1p" in a.uns:
        raise ValueError("X可能已经做过log转换，请先确认原始counts的位置。")
    X = a.X
    source = "X"

print("Count source:", source)

X = X.tocsr() if sparse.issparse(X) else np.asarray(X)
values = X.data if sparse.issparse(X) else X

if not np.isfinite(values).all() or (values < 0).any():
    raise ValueError("计数中存在缺失值、无穷值或负数。")
if not np.equal(values, np.floor(values)).all():
    raise ValueError("发现非整数：请使用原始counts，而不是归一化后的数据。")
if not a.var_names.is_unique:
    raise ValueError("基因标识存在重复，需要先检查a.var。")

# 2. 找出两组都至少有10个细胞的样本
groups = ["Normal_epithelial", "Malignant"]
obs = a.obs[["sample_id", "analysis_group"]].copy()

if obs.isna().any().any():
    raise ValueError("sample_id或analysis_group存在缺失值。")

obs = obs.astype(str)
sizes = pd.crosstab(obs["sample_id"], obs["analysis_group"])
sizes = sizes.reindex(columns=groups, fill_value=0)
paired = sizes.index[sizes.ge(10).all(axis=1)].tolist()

if not paired:
    raise ValueError("没有找到两组均至少有10个细胞的样本。")

# 3. 对每个样本、每个组的原始计数求和
rows = []
records = []

for sample in paired:
    for group in groups:
        mask = (
            (obs["sample_id"] == sample)
            & (obs["analysis_group"] == group)
        ).to_numpy()

        summed = np.asarray(
            X[mask].sum(axis=0, dtype=np.int64)
        ).ravel()

        total = int(summed.sum())
        if total == 0:
            raise ValueError(f"{sample}的{group}组总计数为0。")

        rows.append(summed)
        records.append({
            "profile_id": f"{sample}__{group}",
            "sample_id": sample,
            "analysis_group": group,
            "n_cells": int(mask.sum()),
            "total_counts": total,
        })

# 4. 生成表达矩阵和对应的样本信息
metadata = pd.DataFrame(records).set_index("profile_id")
counts = pd.DataFrame(
    np.vstack(rows),
    index=metadata.index,
    columns=a.var_names,
)

# 5. 保存结果
outputs = {
    "pdac_pseudobulk_counts.csv": counts,
    "pdac_pseudobulk_metadata.csv": metadata,
    "pdac_pseudobulk_genes.csv": a.var,
}

for filename, table in outputs.items():
    table.to_csv(folder / filename)

print("\nPaired samples:", len(paired))
print("Cells used:", int(metadata["n_cells"].sum()))
print("Pseudobulk shape (profiles, genes):", counts.shape)
print()
print(metadata.to_string())
print()

for filename in outputs:
    print("Saved:", folder / filename)
