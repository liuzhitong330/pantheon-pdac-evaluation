# Evaluation design and interpretation

The evaluation was completed on 22 October 2026; this public curation was prepared on 23 October 2026. It uses the saved results, not a new analysis.

## Input-blinded round

PantheonOS 0.6.4, official source commit `8cc6ae86f624322377b799173424e6bada5e3321`, was run through its Python Agent API with `codex/gpt-6.1-sol` and high reasoning. A restricted shell exposed the dataset and an isolated analysis workspace. This is a particular recorded configuration, not a benchmark of every Pantheon interface or model.

The [initial prompt](../protocol/blind_initial_prompt.txt) supplied the biological question and asked for a plan before inference. Inputs were the original count matrices, cell annotations and public sample metadata. Reference code, results, expected genes and preferred statistical settings were withheld from the blind-round session.

The [sealed initial plan](../protocol/blind_initial_plan.txt) correctly distinguished specific epithelial labels from a misleading broad cell-type field, matched 43,430 annotated barcodes and mapped 25 samples to distinct patients. It proposed pseudobulk aggregation and paired inference. Because moderated count-model packages were unavailable initially, it proposed an approximate paired log-CPM test. The evaluator offered dependency provisioning without prescribing a method; Pantheon then selected robust edgeR quasi-likelihood inference. The initial approximate test was not the reported primary analysis.

Before testing, [Pantheon's QC report](../protocol/pantheon_pre_DE_QC.md) identified the hemoglobin-dominated wu16 malignant profile. The evaluator accepted Pantheon's recommendation to exclude the entire pair. The saved [amendment](../workflows/pantheon/outputs/analysis/PLAN_AMENDMENT_2026-10-22.md) uses the historical phrase “User authorized”; in the evaluation record this was the evaluator's acceptance within the authorized task, not an additional direct human QC decision. Its original wording is preserved as provenance.

## Statistical comparison

The reference cohort consists of wu03, wu15, wu16, wu26, wu27, wu28, wu29 and wu42; the Pantheon primary omits wu16. Positive log2 fold changes indicate **higher expression in malignant than annotated normal cells**. Counts were summed within patient and epithelial group; cells were not treated as independent patients.

The reference uses ≥200 detected genes, ≤25% mitochondrial counts, ≥10 cells per group, and genes with ≥10 counts in at least 8 of 16 pseudobulks. Its paired PyDESeq2 0.4.12 workflow and exact settings are in the [four saved scripts](../workflows/reference/reproduction/scripts). The saved parametric dispersion-trend failure and fallback to a mean fit reproduce; it was not silently retuned.

Pantheon uses supplied `filtered=False` specific epithelial labels and nonzero-count cells, ≥10 cells per group, and pre-DE exclusion of wu16. The primary uses raw integer counts, TMM, `~patient+status`, Normal as reference, `filterByExpr(group=status,min.count=10,min.total.count=15,large.n=10,min.prop=0.7)`, robust dispersion estimation and robust `glmQLFit(legacy=FALSE)` followed by `glmQLFTest`. [Saved code](../workflows/pantheon/outputs/analysis/code/run_edger.R) is authoritative. Each discovery cohort has its own BH FDR family. Leave-one-patient-out fits provide influence diagnostics without new discovery p-values.

Codex compared the completed saved outputs using unique Ensembl IDs, each workflow's original FDR values, and FDR <0.05 without a fold-change cutoff. It separated filtered genes from tested nonsignificant genes and checked pseudobulk equality. Correlations are descriptive comparisons on shared data, not independent validation statistics. The [saved full report](SAVED_COMPARISON_REPORT.txt) documents method, filtering and sample-inclusion effects separately.

## Assistance and reproduction evidence

Codex provisioned dependencies, resolved packaging/model-access and sandbox issues, reviewed the pre-DE plan/QC and performed separate audits. Pantheon repaired its own edgeR export/validation-field error. No barcode, pairing or comparison-direction correction was needed. Full private conversations and execution logs remain local; this repository supplies the scientific code, safe plans and derived results rather than those logs.

A fresh second Pantheon session received the reference scripts and a current package snapshot, but no result tables. Its path-only reproduction matched ten required reference artifacts. A separate evaluator replay of the first round matched 372 scientific tables after column alignment, with maximum log2FC/FDR differences of zero. These are **previously completed local checks**, not reruns performed while publishing. Host paths, timestamps and binary metadata were not claimed identical.

Known limits include the small predominantly classical cohort, annotation dependence, wu16's unresolved cause, an extreme HOXB6 estimate with near-zero adjusted residual degrees of freedom, and the preparation checksum scan's old inspection-directory reference. The scientific replay reads the fresh inspection output correctly; the separate prior audit verified its inputs. A model-generated gene nomination, including EFS, is outside this expression comparison's validation scope.
