from pathlib import Path
from importlib.metadata import version

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference


def main():
    folder = Path(__file__).resolve().parents[2]
    if version("pydeseq2") != "0.4.12":
        raise ValueError("This script requires pydeseq2==0.4.12.")

    counts = pd.read_csv(
        folder / "pdac_pseudobulk_counts.csv", index_col=0
    )
    metadata = pd.read_csv(
        folder / "pdac_pseudobulk_metadata.csv", index_col=0
    )

    # 检查并对齐表达矩阵与样本信息
    if not counts.index.is_unique or not metadata.index.is_unique:
        raise ValueError("Duplicate profile identifiers.")
    if set(counts.index) != set(metadata.index):
        raise ValueError("Counts and metadata contain different profiles.")
    metadata = metadata.loc[counts.index].copy()

    # 使用简短标签，避免旧版本解析下划线时出错
    metadata["pair"] = metadata["sample_id"].astype(str)
    metadata["condition"] = metadata["analysis_group"].map({
        "Normal_epithelial": "Normal",
        "Malignant": "Malignant",
    })

    if metadata[["sample_id", "condition"]].isna().any().any():
        raise ValueError("Missing sample ID or unexpected group label.")

    pairs = pd.crosstab(metadata["pair"], metadata["condition"])
    pairs = pairs.reindex(
        columns=["Normal", "Malignant"], fill_value=0
    )
    if len(pairs) != 8 or not pairs.eq(1).all().all():
        raise ValueError(
            "Expected 8 pairs, each with one Normal and one Malignant profile."
        )

    # 确认输入为原始整数计数
    values = counts.to_numpy()
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.equal(values, np.floor(values)).all()
    ):
        raise ValueError("Input must contain nonnegative integer raw counts.")

    # 低表达过滤：至少8份数据中counts达到10
    n_before = counts.shape[1]
    keep = counts.ge(10).sum(axis=0).ge(8)
    pd.DataFrame({"kept": keep}).to_csv(
        folder / "pdac_de_gene_filter.csv"
    )
    counts = counts.loc[:, keep].astype("int64")

    if counts.shape[1] == 0 or counts.sum(axis=1).eq(0).any():
        raise ValueError("No usable genes or an empty profile after filtering.")

    print("Paired samples:", len(pairs), flush=True)
    print("Genes before filtering:", n_before, flush=True)
    print("Genes after filtering:", counts.shape[1], flush=True)
    print(
        "Model: pair + condition; contrast: Malignant vs Normal",
        flush=True,
    )

    # 配对模型：控制样本差异，比较恶性与正常上皮
    inference = DefaultInference(n_cpus=2)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata[["pair", "condition"]].copy(),
        design_factors=["pair", "condition"],
        ref_level=["condition", "Normal"],
        refit_cooks=False,
        inference=inference,
    )

    # 内部完成归一化、离散度估计和表达差异估计
    dds.deseq2()

    stats = DeseqStats(
        dds,
        contrast=["condition", "Malignant", "Normal"],
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
        inference=inference,
    )
    stats.summary()

    # 保存全部经过检验的基因及原始注释
    results = stats.results_df.copy()
    results.index.name = "gene"
    results = results.sort_values("padj", na_position="last")

    gene_info = pd.read_csv(
        folder / "pdac_pseudobulk_genes.csv", index_col=0
    )
    results = results.join(gene_info.add_prefix("annotation_"))
    results.to_csv(folder / "pdac_de_all.csv")

    # 初步候选：校正后p值<0.05，估计变化至少2倍
    significant = results[results["padj"].lt(0.05)]
    candidates = significant[
        significant["log2FoldChange"].abs().ge(1)
    ]
    candidates.to_csv(folder / "pdac_de_candidates.csv")

    # 保存归一化计数和样本信息，供后续检查和画图使用
    normalized = pd.DataFrame(
        dds.layers["normed_counts"],
        index=dds.obs_names,
        columns=dds.var_names,
    )
    normalized.to_csv(folder / "pdac_de_normalized_counts.csv")

    diagnostics = metadata.loc[dds.obs_names].copy()
    diagnostics["size_factor"] = np.asarray(
        dds.obsm["size_factors"]
    ).ravel()
    diagnostics.to_csv(folder / "pdac_de_sample_qc.csv")

    print("\nGenes with padj < 0.05:", len(significant))
    print(
        "Candidates with estimated |log2FC| >= 1:",
        len(candidates),
    )
    print(
        "Up in Malignant:",
        int(candidates["log2FoldChange"].ge(1).sum()),
    )
    print(
        "Down in Malignant:",
        int(candidates["log2FoldChange"].le(-1).sum()),
    )
    print("Genes with missing padj:", int(results["padj"].isna().sum()))

    print("\nTop 20 candidates, ordered by padj:")
    print(
        candidates[["baseMean", "log2FoldChange", "padj"]]
        .head(20)
        .to_string()
    )

    print("\nSaved results to:", folder)
    print("pdac_de_all.csv, pdac_de_candidates.csv")
    print("pdac_de_normalized_counts.csv, pdac_de_sample_qc.csv")
    print("pdac_de_gene_filter.csv")


if __name__ == "__main__":
    main()
