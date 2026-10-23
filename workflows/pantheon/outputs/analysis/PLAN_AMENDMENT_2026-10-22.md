# Approved pre-inference amendment: exclude wu16 from primary

Recorded UTC: 2026-10-22T05:40:54.799659+00:00
User authorized this amendment before any dataset DE test. No dataset DE inference has been performed at the time of writing.

## Rationale and preservation
Independent review confirmed mapping, raw aggregation and the hemoglobin-dominance anomaly. The entire wu16 pair is excluded from primary inference for count-QC reasons, not in response to DE results. Median malignant-cell complexity is 234 genes; median HBB fraction is 84.88%; 103/127 malignant cells have HBB >50% of counts. Preserve every original cell label and count. Do not remove hemoglobin genes. All earlier plans, failure logs and eight-patient preparation files remain unchanged.

## Exact cohorts, roles and FDR families
The machine-readable cohort_specifications.tsv below is authoritative; each listed sample is a distinct mapped patient with a Normal and Malignant pseudobulk. Normal reference; positive statusMalignant log2FC means Malignant minus Normal.

                  analysis_id                     role                                           samples  contains_wu16    source filter_mode                       parent                    fdr_family
       primary7_min10_no_wu16                  primary                wu03;wu15;wu26;wu27;wu28;wu29;wu42          False  combined independent                                     primary7_min10_no_wu16
 qc_original8_min10_with_wu16           qc_sensitivity           wu03;wu15;wu16;wu26;wu27;wu28;wu29;wu42           True  combined independent                               qc_original8_min10_with_wu16
             sensitivity_min5 prespecified_sensitivity wu03;wu07;wu15;wu16;wu26;wu27;wu28;wu29;wu39;wu42           True  combined independent                                           sensitivity_min5
            sensitivity_min20 prespecified_sensitivity                          wu03;wu16;wu26;wu27;wu42           True  combined independent                                          sensitivity_min20
              classical_min10 prespecified_sensitivity           wu03;wu15;wu16;wu26;wu27;wu28;wu29;wu42           True classical independent                                            classical_min10
sensitivity_min5_without_wu39 prespecified_sensitivity      wu03;wu07;wu15;wu16;wu26;wu27;wu28;wu29;wu42           True  combined independent                              sensitivity_min5_without_wu39
sensitivity_min5_without_wu07 prespecified_sensitivity      wu03;wu15;wu16;wu26;wu27;wu28;wu29;wu39;wu42           True  combined independent                              sensitivity_min5_without_wu07
    loo_primary7_without_wu03     influence_diagnostic                     wu15;wu26;wu27;wu28;wu29;wu42          False  combined parent_mask       primary7_min10_no_wu16    NONE_effect_stability_only
    loo_primary7_without_wu15     influence_diagnostic                     wu03;wu26;wu27;wu28;wu29;wu42          False  combined parent_mask       primary7_min10_no_wu16    NONE_effect_stability_only
    loo_primary7_without_wu26     influence_diagnostic                     wu03;wu15;wu27;wu28;wu29;wu42          False  combined parent_mask       primary7_min10_no_wu16    NONE_effect_stability_only
    loo_primary7_without_wu27     influence_diagnostic                     wu03;wu15;wu26;wu28;wu29;wu42          False  combined parent_mask       primary7_min10_no_wu16    NONE_effect_stability_only
    loo_primary7_without_wu28     influence_diagnostic                     wu03;wu15;wu26;wu27;wu29;wu42          False  combined parent_mask       primary7_min10_no_wu16    NONE_effect_stability_only
    loo_primary7_without_wu29     influence_diagnostic                     wu03;wu15;wu26;wu27;wu28;wu42          False  combined parent_mask       primary7_min10_no_wu16    NONE_effect_stability_only
    loo_primary7_without_wu42     influence_diagnostic                     wu03;wu15;wu26;wu27;wu28;wu29          False  combined parent_mask       primary7_min10_no_wu16    NONE_effect_stability_only
   loo_original8_without_wu03     influence_diagnostic                wu15;wu16;wu26;wu27;wu28;wu29;wu42           True  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only
   loo_original8_without_wu15     influence_diagnostic                wu03;wu16;wu26;wu27;wu28;wu29;wu42           True  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only
   loo_original8_without_wu16     influence_diagnostic                wu03;wu15;wu26;wu27;wu28;wu29;wu42          False  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only
   loo_original8_without_wu26     influence_diagnostic                wu03;wu15;wu16;wu27;wu28;wu29;wu42           True  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only
   loo_original8_without_wu27     influence_diagnostic                wu03;wu15;wu16;wu26;wu28;wu29;wu42           True  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only
   loo_original8_without_wu28     influence_diagnostic                wu03;wu15;wu16;wu26;wu27;wu29;wu42           True  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only
   loo_original8_without_wu29     influence_diagnostic                wu03;wu15;wu16;wu26;wu27;wu28;wu42           True  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only
   loo_original8_without_wu42     influence_diagnostic                wu03;wu15;wu16;wu26;wu27;wu28;wu29           True  combined parent_mask qc_original8_min10_with_wu16    NONE_effect_stability_only

Primary7 comprises seven pairs, 1,856 malignant (1,837 classical + 19 basal) and 249 normal cells. Original8 is only a QC sensitivity, NOT a competing primary answer. Threshold and subtype sensitivities retain their original membership, explicitly including wu16. Classical sensitivity uses Malignant - Classical only versus the same Normal Epithelial cells; it is also the original eight-patient cohort. The min5-without-wu39 and min5-without-wu07 analyses remove only the named patient and still contain wu16. No sensitivity is chosen based on discoveries. Every discovery analysis has its own BH family, with no pooled FDR assertion across families.

## Unchanged primary/test settings
Raw integer pseudobulk count input, unique Ensembl rows, no imputation/integration, downsampling or further cell relabeling. DGEList; full-rank ~patient+status; Normal reference. Stage/location/pathology/survival main effects not added because of patient blocking. filterByExpr(group=status,min.count=10,min.total.count=15,large.n=10,min.prop=0.7). In Primary7 this uses a seven-library minimum and a CPM cutoff determined from median raw library size; save the mask and derived cutoff. Subset keep.lib.sizes=FALSE; calcNormFactors(method="TMM") with edgeR 4.4.0 defaults (logratioTrim=.3,sumTrim=.05,doWeighting=TRUE,Acutoff=-1e10,automatic reference); record factors and effective library offsets.

estimateDisp(design,robust=TRUE); glmQLFit(design,robust=TRUE,abundance.trend=TRUE,legacy=FALSE); glmQLFTest(coef="statusMalignant",poisson.bound=TRUE). Other package defaults from the pinned edgeR version are retained and saved via function-formals manifest. BH FDR <.05 across all genes tested within each discovery family. abs(log2FC)>=1 is a descriptive effect-size prioritization flag, not a minimum-effect hypothesis. No ad hoc confidence intervals. Preserve all 33,538 gene rows in full-results files with explicit not-tested statuses. Preserve fits, dispersions, QL df, offsets, metadata, masks, warnings and commands.

## Influence diagnostics
Leave each of the seven Primary7 patients out in turn: six-patient fits using the Primary7 gene mask. Separately leave each Original8 patient out: seven-patient fits using the Original8 mask. Recompute TMM and dispersions per reduced cohort; retain status direction and report effects/estimability. These are influence diagnostics, not new discovery lists: no new p-values/FDR are computed for them. Removing wu16 from Original8 matches Primary7 patient membership but may use a different mask; it is not a new primary. Any all-zero or not-estimable genes in a reduced cohort are explicitly flagged, not assigned spurious significance.

## Software and readiness
Requested runtime now accessible. Run saved check_r_inputs.R directly with tool-supplied temporary environment, never old TMPDIR override. Verify R 4.4.2, Bioconductor release 3.20 via BiocVersion 3.20.0, edgeR 4.4.0, limma 3.62.0 and statmod 1.5.0. Full software manifest/sessionInfo, function defaults, RDS models and logs will accompany results. Existing original R failure logs remain untouched. Tests may proceed under this approved amendment without another pause, with invalid fits reported rather than silently substituted.
