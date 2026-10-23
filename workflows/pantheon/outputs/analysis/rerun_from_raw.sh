#!/bin/zsh
# Run from workspace root; NEVER overwrite original preparation/plans.
set -eu
ROOT=${1:-outputs/full_rerun_$(date -u +%Y%m%dT%H%M%SZ)}
if [[ -e "$ROOT" ]]; then print -u2 "Refusing existing output directory: $ROOT"; exit 1; fi
mkdir -p "$ROOT/inspection" "$ROOT/preparation/logs" "$ROOT/analysis/logs"
export INSPECTION_OUTPUT="$ROOT/inspection"
export PREP_OUTPUT="$ROOT/preparation"
export ANALYSIS_OUTPUT="$ROOT/analysis"
# Runtime tool provides accessible TMPDIR: do not override it.
.analysis-venv/bin/python outputs/analysis/code/rerun_inspection.py > "$ROOT/inspection/inspection_report.txt" 2>&1
.analysis-venv/bin/python outputs/analysis/code/rerun_prepare_counts.py > "$PREP_OUTPUT/logs/preparation.log" 2>&1
.analysis-venv/bin/python outputs/analysis/code/rerun_validate_and_review.py > "$PREP_OUTPUT/logs/validation.log" 2>&1
cp outputs/analysis/cohort_specifications.tsv "$ANALYSIS_OUTPUT/"
cp outputs/analysis/PLAN_AMENDMENT_2026-10-22.md "$ANALYSIS_OUTPUT/"
.r-env/bin/Rscript --vanilla outputs/analysis/code/run_edger.R > "$ANALYSIS_OUTPUT/logs/run_edger.log" 2>&1
.r-env/bin/Rscript --vanilla outputs/analysis/code/audit_fits.R > "$ANALYSIS_OUTPUT/logs/audit_fits.log" 2>&1
.analysis-venv/bin/python outputs/analysis/code/summarize_results.py > "$ANALYSIS_OUTPUT/logs/summarize_results.log" 2>&1
.r-env/bin/Rscript --vanilla outputs/analysis/code/add_presentation_plots.R > "$ANALYSIS_OUTPUT/logs/presentation.log" 2>&1
print "Raw-input-to-results rerun complete: $ROOT"
