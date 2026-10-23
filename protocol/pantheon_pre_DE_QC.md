# Raw-count QC and pseudobulk preparation: completed, inference on hold

## Scope and provenance
Executed raw-count validation, descriptive cell QC, patient/status and patient/subtype count summation, eligibility manifests, and design/contrast validation. No DE inference, abundance gene filtering, TMM normalization, dispersion fitting, ranking, or DE-selected gene list was produced. All 33,538 genes remain in each count matrix. Earlier plans and inspection outputs are unchanged; checksums of all supplied inputs and prior artifacts are in input_and_prior_artifact_checksums.tsv.

## Integrity and cohort results
- All 25 integer matrices parsed with declared entry totals and row/column bounds checked; no negative values or explicit zero entries found. Feature order is identical. No duplicate nonzero gene-cell coordinates occurred among annotated cells.
- All 43,430 annotated IDs match unique sample barcodes. Cells outside the supplied annotation are not called or added to the analysis.
- All 11,133 selected epithelial cells have positive counts; none were excluded. Preserve 10,783 malignant (5,123 basal + 5,660 classical) and 350 Normal Epithelial cells.
- Created 46 observed patient/status pseudobulks (25 malignant and 21 normal) and 58 observed patient/specific-subtype pseudobulks. Missing normal groups are not encoded as zero-count comparator libraries.
- Saved integer sums exactly reconcile with retained-cell UMI sums; combined malignant counts reconcile with basal-plus-classical sums. Sample/patient mapping is one-to-one and supplied clinical fields agree after whitespace trimming.
- Primary eligibility remains 8 patients / 16 libraries, 1,983 malignant cells (1,964 classical + 19 basal), 318 normal. Minimum-5 sensitivity: 10 pairs; minimum-20: 5 pairs; classical minimum-10: 8 pairs.
- Saved metadata has a deterministic Normal-then-Malignant ordering within sample. The primary design has rank 9 with 16 rows and 7 nominal residual df. Each within-patient Malignant-minus-Normal design-row difference is exactly the saved statusMalignant contrast. Reference status: Normal; positive log2FC: higher in Malignant. No clinical covariates constant within patient are added to the paired design.

## Severe QC issue requiring review BEFORE testing: wu16
Observed count composition, not DE:
- 127 malignant cells are all labeled Malignant - Classical. Median UMI=5,657 but median detected genes=234, versus 11,764 UMI / 3,476 genes in 69 normal epithelial cells from the same patient.
- HBB is the top-count gene in 105/127 malignant cells. HBB alone exceeds 50% of all counts in 103/127 cells; median HBB fraction across malignant cells is 84.88%. Median combined HBA1/HBA2/HBB fraction is 89.53%.
- The malignant pseudobulk contains 656,628 UMIs; HBB accounts for 80.30%, and HBA1/HBA2/HBB together 84.49%. The normal pseudobulk hemoglobin fraction is 2.16% (median individual normal-cell hemoglobin fraction 0.058%).
- 111/127 malignant cells have fewer than 500 detected genes. Other primary malignant groups have median gene counts 2,346.5-4,121.5 and median single-gene fractions approximately 2.9-3.6%.

Interpretation: this is a severe, group-specific low-complexity/hemoglobin-dominance anomaly, consistent with erythroid contamination and/or annotation error; ambient RNA, mixed cells, and other explanations cannot be distinguished from supplied counts alone. No relabeling or doublet classification is justified here. TMM may address some composition shifts but cannot establish that these are valid malignant epithelial cells. Simply removing hemoglobin genes would conceal the QC concern without restoring epithelial information.

Recommended decision for review: exclude the entire wu16 pair from primary inference, retain all original-count files, and treat the originally planned eight-pair cohort only as a clearly labeled QC sensitivity. This would leave seven primary pairs, 1,856 malignant and 249 normal cells. This is a RECOMMENDATION ONLY and has NOT been applied. Alternatively, retaining the pair requires explicit acceptance of this quality concern; targeted cell reclassification/filtering would need a separately approved, count-QC-defined rule. Testing remains on hold until this decision is made.

## Other descriptive QC flags (not exclusions)
Among all retained epithelial cells, 969 have mitochondrial fraction >20%, none >50%; 440 exceed 50,000 UMIs and 566 exceed 7,000 genes. In the original primary cohort, these counts are 263, 0, 141, and 120, respectively. No retained epithelial cell has fewer than 200 genes. These are transparent review flags only, not new filtering cutoffs or evidence sufficient to call doublets. High-count cells and moderate mitochondrial fractions vary across patients; no additional automatic exclusion is proposed. Supplied filtering is retained. Per-cell and group quantiles are saved, including all filtered/non-epithelial annotated cells for context.

## Software / execution limitation
Python 3.12.14, NumPy 2.5.4, pandas 3.0.6, SciPy 1.18.1 executed preparation successfully. Full installed-Python manifest is in python_versions.json.
The newly provided .r-env executable cannot initialize its library under current sandbox permissions. First startup failed because R_TempDir contains a space; relative TMPDIR avoids that issue, but R then reports Operation not permitted when resolving its provisioned library and cannot load base packages. Full errors are saved in logs/r_initial_probe.log, logs/r_relative_tmp_probe.log and logs/check_r_inputs.log. No R package version can yet be verified and no R sessionInfo was obtainable. Do not infer successful installation from the requested versions. An input-only R verification script is saved for when this runtime access is repaired; Python checks already establish the intended design algebra and ordering. No attempt was made to access other project/user data.

## Key artifacts
- retained_epithelial_cells.tsv.gz: exact IDs, sample/patient, specific and broad labels, grouping, barcode column, QC and both pseudobulk memberships.
- all_annotated_cell_qc.tsv.gz: same provenance and selection/exclusion flags for every annotated cell.
- feature_mapping.tsv: 1-based matrix row, unique Ensembl ID, nonunique symbol, feature type, mitochondrial/ribosomal symbol flags.
- raw_pseudobulk_counts_all.tsv.gz and raw_subtype_pseudobulk_counts_all.tsv.gz: all observed group sums, genes as rows, explicit pseudobulk IDs as columns, raw integers.
- primary_min10_raw_counts.tsv.gz and primary_min10_metadata.tsv: original-plan primary inputs, provisional pending wu16 review. classical_min10 equivalents also saved.
- sample_patient_metadata.tsv, pseudobulk_metadata_all.tsv, subtype_pseudobulk_metadata_all.tsv: clinical sources, patient mapping, cell numbers, library sizes.
- sample_eligibility.tsv: decisions/reasons for all 25 samples at each threshold.
- cohort_pseudobulk_manifest.tsv and cohort_definitions.json: deterministic input order and membership for planned cohorts and leave-one-patient-out checks.
- design_*.tsv, contrast_*.tsv, design_validation.tsv: explicit full-rank paired designs and contrast definitions; no model fitting.
- cell_qc_by_sample_specific_label.tsv, selected_group_composition_qc.tsv, qc_review_flag_counts.tsv, raw_matrix_validation.tsv: complete descriptive summaries and validation results.
- wu16_cell_composition_review.tsv.gz, wu16_composition_summary.tsv, hemoglobin_composition_all_pseudobulks.tsv: traceable composition-review evidence.
- code/, reproduce.sh, execution_commands.txt, python_versions.json, logs/: executable code, parameters, commands, versions and full run/error logs. checksums.tsv covers final artifacts except itself.
