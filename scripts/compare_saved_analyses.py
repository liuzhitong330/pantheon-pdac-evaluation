#!/usr/bin/env python3
"""Independent evaluator comparison of frozen results. Does not fit any model.

Run with round2/.analysis-venv/bin/python. Default paths derive from this file.
All gene joins use unique Ensembl IDs; FDR threshold is strictly <0.05, with no
fold-change cutoff. Missing statistics for filtered genes remain missing.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--evaluation', type=Path, default=Path(__file__).resolve().parents[1] / 'work/saved/evaluation')
parser.add_argument('--original', type=Path, default=Path(__file__).resolve().parents[1] / 'work/saved/reference')
parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'work/comparison')
args = parser.parse_args()
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=True)
analysis = args.evaluation / 'round1/outputs/analysis'
prep = args.evaluation / 'round1/outputs/preparation'
sources = {}

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read(path, **kwargs):
    path = path.resolve()
    st = path.stat()
    sources[str(path)] = {'sha256': sha(path), 'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns}
    return pd.read_csv(path, **kwargs)

def save(df, name, index=False):
    df.to_csv(out / name, sep='\t', index=index, na_rep='NA', float_format='%.17g')

def write_json(value, name):
    (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def state(tested, fdr):
    return pd.Series(np.select([~tested, fdr.isna(), fdr < .05],
                     ['filtered_out', 'tested_no_adjusted_p', 'significant'],
                     default='tested_not_significant'), index=fdr.index)

def effect_stats(x, y):
    keep = x.notna() & y.notna()
    x, y = x[keep], y[keep]
    d = (y - x).abs()
    return {'n': len(x), 'pearson': float(x.corr(y)) if len(x) > 1 else None,
            'spearman': float(spearmanr(x, y).statistic) if len(x) > 1 else None,
            'same_sign': int((np.sign(x) == np.sign(y)).sum()),
            'opposite_sign': int((x * y < 0).sum()),
            'median_absolute_log2FC_difference': float(d.median()) if len(x) else None,
            'mean_absolute_log2FC_difference': float(d.mean()) if len(x) else None,
            'p90_absolute_log2FC_difference': float(d.quantile(.9)) if len(x) else None,
            'p95_absolute_log2FC_difference': float(d.quantile(.95)) if len(x) else None,
            'max_absolute_log2FC_difference': float(d.max()) if len(x) else None}

genes = read(args.original / 'pdac_pseudobulk_genes.csv', index_col=0)
filt = read(args.original / 'pdac_de_gene_filter.csv', index_col=0)
old = read(args.original / 'pdac_de_all.csv')
assert genes.index.is_unique and genes.gene_ids.is_unique
assert filt.index.equals(genes.index) and filt.kept.dtype == bool
assert old.annotation_gene_ids.is_unique and old.gene.is_unique
assert set(old.gene) == set(filt.index[filt.kept])
assert (old.set_index('gene').annotation_gene_ids == genes.loc[old.gene, 'gene_ids']).all()
universe = genes.gene_ids.tolist()
old = old.set_index('annotation_gene_ids').reindex(universe)
base = pd.DataFrame(index=pd.Index(universe, name='ensembl_id'))
base['original_gene_label'] = genes.index.to_numpy()
base['original_tested'] = filt.kept.to_numpy()
for src, dst in [('log2FoldChange', 'log2FC'), ('pvalue', 'pvalue'), ('padj', 'FDR'), ('lfcSE', 'log2FC_SE')]:
    base['original_' + dst] = old[src]
base['original_status'] = state(base.original_tested, base.original_FDR)

specs = read(analysis / 'cohort_specifications.tsv', sep='\t')
discovery = specs[specs.fdr_family != 'NONE_effect_stability_only'].copy()
cohorts = {}
for row in discovery.itertuples():
    d = read(analysis / row.analysis_id / 'all_genes_results.tsv', sep='\t').set_index('ensembl_id')
    assert d.index.is_unique and set(d.index) == set(universe)
    d = d.reindex(universe)
    d['status'] = state(d.test_status.eq('tested'), d.FDR)
    assert d.loc[d.test_status.eq('tested'), ['logFC', 'PValue', 'FDR']].notna().all().all()
    assert d.loc[~d.test_status.eq('tested'), ['logFC', 'PValue', 'FDR']].isna().all().all()
    cohorts[row.analysis_id] = d
primary_id = 'primary7_min10_no_wu16'
eight_id = 'qc_original8_min10_with_wu16'
primary, eight = cohorts[primary_id], cohorts[eight_id]
base.insert(0, 'gene_symbol', primary.gene_symbol)
for prefix, d in [('pantheon', primary), ('edgeR8_sensitivity', eight)]:
    base[prefix + '_tested'] = d.test_status.eq('tested')
    base[prefix + '_status'] = d.status
    for src, dst in [('logFC', 'log2FC'), ('PValue', 'pvalue'), ('FDR', 'FDR')]:
        base[prefix + '_' + dst] = d[src]
osig = base.original_status.eq('significant')
psig = base.pantheon_status.eq('significant')
base['discovery_category'] = np.select([osig & psig, osig, psig],
                 ['shared_significant', 'original_only', 'pantheon_only'], default='neither_significant')
base['other_analysis_status_for_unique_gene'] = np.select([osig & ~psig, psig & ~osig],
                 [base.pantheon_status, base.original_status], default='not_applicable')
base['log2FC_difference_Pantheon_minus_original'] = base.pantheon_log2FC - base.original_log2FC
base['absolute_log2FC_difference'] = base.log2FC_difference_Pantheon_minus_original.abs()
common = base.original_tested & base.pantheon_tested
base['direction_agrees'] = pd.Series(pd.NA, index=base.index, dtype='boolean')
base.loc[common, 'direction_agrees'] = np.sign(base.loc[common, 'original_log2FC']) == np.sign(base.loc[common, 'pantheon_log2FC'])
annot = read(analysis / 'primary_all_genes_annotated.tsv', sep='\t').set_index('ensembl_id')
for col in ['n_patients_direction_agree', 'all_loo_same_direction', 'loo_logFC_min', 'loo_logFC_max',
            'most_influential_omission', 'max_abs_loo_effect_change', 'flag_near_zero_adjusted_residual_df']:
    base['pantheon_' + col] = annot[col]
save(base, 'gene_comparison_all_33538.tsv', True)
save(base[common].sort_values('absolute_log2FC_difference', ascending=False), 'common_tested_effects.tsv', True)
save(base[common & ~base.direction_agrees.fillna(False)].sort_values('absolute_log2FC_difference', ascending=False),
     'opposite_direction_genes.tsv', True)
save(base[osig & psig].sort_values('pantheon_FDR'), 'shared_significant_genes.tsv', True)
save(base[osig & ~psig].sort_values('original_FDR'), 'original_only_significant_genes.tsv', True)
save(base[psig & ~osig].sort_values('pantheon_FDR'), 'pantheon_only_significant_genes.tsv', True)
save(pd.crosstab(base.original_status, base.pantheon_status), 'test_status_cross_table.tsv', True)

summary = {'comparison_type': 'Independent evaluator; saved results only, no model refits',
           'gene_key': 'unique Ensembl ID', 'FDR_threshold': '<0.05', 'fold_change_cutoff': None,
           'direction': 'positive log2FC means higher in malignant than annotated normal epithelial cells',
           'n_genes_in_input': len(base), 'original_tested': int(base.original_tested.sum()),
           'pantheon_tested': int(base.pantheon_tested.sum()), 'common_tested': int(common.sum()),
           'original_only_tested': int((base.original_tested & ~base.pantheon_tested).sum()),
           'pantheon_only_tested': int((~base.original_tested & base.pantheon_tested).sum()),
           'neither_tested': int((~base.original_tested & ~base.pantheon_tested).sum()),
           'original_significant': int(osig.sum()), 'pantheon_significant': int(psig.sum()),
           'shared_significant': int((osig & psig).sum()),
           'original_only_significant': int((osig & ~psig).sum()),
           'pantheon_only_significant': int((psig & ~osig).sum()),
           'original_only_split_by_pantheon_status': base.loc[osig & ~psig, 'pantheon_status'].value_counts().to_dict(),
           'pantheon_only_split_by_original_status': base.loc[psig & ~osig, 'original_status'].value_counts().to_dict(),
           'effect_agreement': {}}
for name, mask in [('common_tested', common), ('shared_significant', osig & psig),
                   ('significant_in_either_and_tested_in_both', common & (osig | psig)),
                   ('both_absolute_log2FC_at_least_0.5', common & base.original_log2FC.abs().ge(.5) & base.pantheon_log2FC.abs().ge(.5))]:
    summary['effect_agreement'][name] = effect_stats(base.loc[mask, 'original_log2FC'], base.loc[mask, 'pantheon_log2FC'])
flips = common & base.original_log2FC.mul(base.pantheon_log2FC).lt(0)
summary['direction_flip_diagnostics'] = {
    'n': int(flips.sum()), 'max_abs_original_log2FC': float(base.loc[flips, 'original_log2FC'].abs().max()),
    'max_abs_pantheon_log2FC': float(base.loc[flips, 'pantheon_log2FC'].abs().max()),
    'minimum_original_FDR': float(base.loc[flips, 'original_FDR'].min()),
    'minimum_pantheon_FDR': float(base.loc[flips, 'pantheon_FDR'].min())}
summary_rows, transition_rows = [], []
for row in discovery.itertuples():
    d = cohorts[row.analysis_id]
    sig = d.status.eq('significant')
    tested = d.test_status.eq('tested')
    summary_rows.append({'analysis_id': row.analysis_id, 'role': row.role, 'patients': row.samples,
        'n_tested': int(tested.sum()), 'n_significant': int(sig.sum()),
        'n_overlap_original_significant': int((sig & osig).sum()),
        'n_overlap_primary_significant': int((sig & psig).sum()),
        'n_original_only_filtered_here': int((osig & ~tested).sum()),
        'n_original_only_tested_nonsignificant_here': int((osig & tested & ~sig).sum()),
        **{'vs_original_' + k: v for k, v in effect_stats(base.original_log2FC, d.logFC).items()}})
    for label, mask in [('original_only', osig & ~psig), ('pantheon_only', psig & ~osig), ('shared', osig & psig)]:
        for status, n in d.loc[mask, 'status'].value_counts().items():
            transition_rows.append({'analysis_id': row.analysis_id, 'primary_comparison_category': label, 'status': status, 'n': int(n)})
save(pd.DataFrame(summary_rows), 'saved_cohort_comparisons.tsv')
save(pd.DataFrame(transition_rows), 'unique_gene_status_across_saved_cohorts.tsv')
pairwise = []
for left, right, change in [
    (eight_id, primary_id, 'exclude wu16 pair for pre-DE composition QC'),
    (eight_id, 'sensitivity_min5', 'lower minimum cells from 10 to 5; add wu07 and wu39'),
    (eight_id, 'sensitivity_min20', 'raise minimum cells from 10 to 20; exclude wu15, wu28 and wu29'),
    ('sensitivity_min5', 'sensitivity_min5_without_wu39', 'exclude basal-rich wu39 from min5 cohort'),
    ('sensitivity_min5', 'sensitivity_min5_without_wu07', 'exclude adenosquamous wu07 from min5 cohort'),
    (eight_id, 'classical_min10', 'restrict malignant group to classical-labelled cells')]:
    a, b = cohorts[left], cohorts[right]
    sa, sb = a.status.eq('significant'), b.status.eq('significant')
    pairwise.append({'analysis_A': left, 'analysis_B': right, 'change': change,
        'A_significant': int(sa.sum()), 'B_significant': int(sb.sum()),
        'shared_significant': int((sa & sb).sum()), 'A_only_significant': int((sa & ~sb).sum()),
        'B_only_significant': int((sb & ~sa).sum()),
        'A_only_filtered_by_B': int((sa & b.status.eq('filtered_out')).sum()),
        'B_only_filtered_by_A': int((sb & a.status.eq('filtered_out')).sum()),
        'opposite_direction_significant_in_both': int((sa & sb & (a.logFC * b.logFC < 0)).sum()),
        **effect_stats(a.logFC, b.logFC)})
save(pd.DataFrame(pairwise), 'sample_and_subtype_sensitivity_comparisons.tsv')
save(pd.crosstab(eight.status, primary.status, rownames=['edgeR8_status'], colnames=['primary7_status']),
     'wu16_exclusion_status_transitions.tsv', True)
effect = base[['gene_symbol', 'original_status', 'original_log2FC', 'original_FDR',
               'pantheon_status', 'pantheon_log2FC', 'pantheon_FDR',
               'edgeR8_sensitivity_status', 'edgeR8_sensitivity_log2FC', 'edgeR8_sensitivity_FDR']].copy()
loo = read(analysis / 'loo_original8_without_wu16/all_genes_results.tsv', sep='\t').set_index('ensembl_id').reindex(universe)
assert loo.FDR.isna().all() and loo.PValue.isna().all()
assert loo.test_status.eq('effect_diagnostic_only_no_test').equals(eight.test_status.eq('tested'))
assert loo.logFC.notna().equals(eight.test_status.eq('tested'))
effect['fixed_8patient_gene_mask_without_wu16_log2FC'] = loo.logFC
effect['fixed_mask_exclusion_effect_change'] = loo.logFC - eight.logFC
effect['primary7_minus_edgeR8_log2FC'] = primary.logFC - eight.logFC
effect['primary7_minus_fixedmask7_log2FC'] = primary.logFC - loo.logFC
effect['abs_primary7_minus_edgeR8_log2FC'] = effect.primary7_minus_edgeR8_log2FC.abs()
save(effect.sort_values('abs_primary7_minus_edgeR8_log2FC', ascending=False), 'wu16_exclusion_gene_effects.tsv', True)
summary['wu16_exclusion'] = {
    'edgeR8_significant': int(eight.status.eq('significant').sum()),
    'shared_edgeR8_primary7_significant': int((eight.status.eq('significant') & psig).sum()),
    'primary7_only_vs_edgeR8': int((psig & ~eight.status.eq('significant')).sum()),
    'edgeR8_only_vs_primary7': int((~psig & eight.status.eq('significant')).sum()),
    'edgeR8_vs_primary7_effects': effect_stats(eight.logFC, primary.logFC),
    'edgeR8_vs_fixed_mask_without_wu16_effects': effect_stats(eight.logFC, loo.logFC),
    'primary7_vs_fixed_mask_without_wu16_effects': effect_stats(loo.logFC, primary.logFC),
    'interpretation': 'Saved sensitivity only. Exclusion refits normalization/dispersion and changes independent abundance filtering. Fixed-mask leave-one-out supplies effects only, no discovery p-values.'}

# Verify actual pseudobulk inputs to distinguish pipeline effects from cell/count changes.
old_counts = read(args.original / 'pdac_pseudobulk_counts.csv', index_col=0).T
old_counts.index = genes.loc[old_counts.index, 'gene_ids'].to_numpy()
old_counts.columns = old_counts.columns.str.replace('Normal_epithelial', 'Normal', regex=False)
count_checks = {}
raw = {}
for label in [eight_id, primary_id]:
    c = read(analysis / label / 'raw_counts.tsv.gz', sep='\t', index_col=0)
    aligned = old_counts.loc[c.index, c.columns]
    assert np.array_equal(aligned.to_numpy(), c.to_numpy())
    count_checks[label] = {'shape_genes_by_profiles': list(c.shape), 'exact_original_subset': True,
                           'total_counts': int(c.to_numpy().sum()), 'different_entries': 0}
    raw[label] = c.reindex(universe)
fp = read(analysis / primary_id / 'filter_parameters.tsv', sep='\t').iloc[0]
raw7 = raw[primary_id]
filter_audit = base[['gene_symbol', 'original_status', 'pantheon_status', 'discovery_category']].copy()
filter_audit['original_profiles_counts_ge10'] = (old_counts.reindex(universe) >= 10).sum(axis=1)
filter_audit['original_required_profiles'] = 8
filter_audit['primary7_profiles_CPM_ge_cutoff'] = (raw7.div(raw7.sum(axis=0), axis=1) * 1e6 >= fp.derived_CPM_cutoff).sum(axis=1)
filter_audit['primary7_CPM_cutoff'] = fp.derived_CPM_cutoff
filter_audit['primary7_required_profiles'] = fp.independent_min_libraries
filter_audit['primary7_total_counts'] = raw7.sum(axis=1)
filter_audit['primary7_required_total_counts'] = fp.min_total_count
assert filter_audit.original_profiles_counts_ge10.ge(8).equals(base.original_tested)
recomputed_primary = filter_audit.primary7_profiles_CPM_ge_cutoff.ge(fp.independent_min_libraries) & filter_audit.primary7_total_counts.ge(fp.min_total_count)
assert recomputed_primary.equals(base.pantheon_tested)
save(filter_audit, 'gene_filter_audit.tsv', True)
examples = ['S100P', 'MIA', 'MUC5AC', 'FOLR1', 'NR5A2', 'CRISP3', 'CYP2S1', 'PLAC8',
            'IL2RG', 'TNS4', 'SERPINA6', 'HABP2', 'ANPEP', 'TPM2', 'GABRP', 'OLFM4',
            'GABARAP', 'HBB', 'HOXB6', 'RGS10', 'CRP']
example_table = base[base.gene_symbol.isin(examples)].copy()
for col in filter_audit.columns:
    if col not in example_table:
        example_table[col] = filter_audit[col]
example_table['display_order'] = example_table.gene_symbol.map({g: i for i, g in enumerate(examples)})
save(example_table.sort_values('display_order').drop(columns='display_order'), 'representative_gene_examples.tsv', True)
old_meta = read(args.original / 'pdac_pseudobulk_metadata.csv')
elig = read(prep / 'sample_eligibility.tsv', sep='\t')
for row in discovery.itertuples():
    elig['included_' + row.analysis_id] = elig.sample_id.isin(row.samples.split(';'))
elig['included_original'] = elig.sample_id.isin(old_meta.sample_id)
elig['primary_exclusion_reason'] = np.where(elig['included_' + primary_id], 'included',
    np.where(elig.sample_id.eq('wu16'), 'pre-DE hemoglobin/composition QC; entire pair excluded', elig.reason_min10))
save(elig, 'patient_inclusion_comparison.tsv')
summary['pseudobulk_input_validation'] = count_checks
summary['original_patients'] = sorted(old_meta.sample_id.unique().tolist())
summary['original_cells'] = {str(k): int(v) for k, v in old_meta.groupby('analysis_group').n_cells.sum().items()}
primary_meta = read(analysis / primary_id / 'metadata.tsv', sep='\t')
summary['pantheon_cells'] = {str(k): int(v) for k, v in primary_meta.groupby('group').n_cells.sum().items()}
assert summary['original_significant'] == summary['shared_significant'] + summary['original_only_significant']
assert summary['pantheon_significant'] == summary['shared_significant'] + summary['pantheon_only_significant']
write_json(summary, 'comparison_summary.json')
save(pd.DataFrame([{'analysis': 'Original PyDESeq2, 8 patients', 'tested': int(base.original_tested.sum()),
                   'FDR_lt_0.05': int(osig.sum()), 'higher_in_malignant': int((osig & base.original_log2FC.gt(0)).sum()),
                   'lower_in_malignant': int((osig & base.original_log2FC.lt(0)).sum())},
                  {'analysis': 'Pantheon edgeR primary, 7 patients', 'tested': int(base.pantheon_tested.sum()),
                   'FDR_lt_0.05': int(psig.sum()), 'higher_in_malignant': int((psig & base.pantheon_log2FC.gt(0)).sum()),
                   'lower_in_malignant': int((psig & base.pantheon_log2FC.lt(0)).sum())}]), 'headline_counts.tsv')

# A static scientific figure; all plotted values come from the saved result tables.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.5), gridspec_kw={'width_ratios': [1.2, 1]})
colors = {'neither_significant': '#b7bfc8', 'shared_significant': '#147d92',
          'original_only': '#d77d20', 'pantheon_only': '#8254a2'}
names = {'neither_significant': 'Neither significant', 'shared_significant': 'Shared significant',
         'original_only': 'Original only', 'pantheon_only': 'Pantheon only'}
for cat in colors:
    d = base[common & base.discovery_category.eq(cat)]
    axes[0].scatter(d.original_log2FC, d.pantheon_log2FC, s=8, alpha=.45 if cat == 'neither_significant' else .75,
                    c=colors[cat], label=f'{names[cat]} ({len(d):,})', rasterized=True, linewidths=0)
lo = min(base.loc[common, 'original_log2FC'].min(), base.loc[common, 'pantheon_log2FC'].min()) - .3
hi = max(base.loc[common, 'original_log2FC'].max(), base.loc[common, 'pantheon_log2FC'].max()) + .3
axes[0].plot([lo, hi], [lo, hi], color='#394653', lw=.8, ls='--')
axes[0].axhline(0, color='#dedede', lw=.5); axes[0].axvline(0, color='#dedede', lw=.5)
axes[0].set(xlabel='Original log2 fold change (8 patients)', ylabel='Pantheon log2 fold change (7 patients)',
            title='6,601 genes tested in both: r = 0.988', xlim=(lo, hi), ylim=(lo, hi))
axes[0].legend(frameon=False, fontsize=8, loc='upper left')
groups = ['Original\n1,373 significant', 'Pantheon\n1,007 significant']
bottom = np.array([0, 0])
segments = [('Shared significant', [int((osig & psig).sum())] * 2, '#147d92'),
            ('Unique: other tested, not significant', [int((osig & ~psig & common).sum()), int((psig & ~osig & common).sum())], '#d77d20'),
            ('Unique: filtered out by other', [int((osig & ~base.pantheon_tested).sum()), int((psig & ~base.original_tested).sum())], '#8f9baa')]
for label, vals, color in segments:
    axes[1].bar(groups, vals, bottom=bottom, color=color, width=.57, label=label)
    for i, value in enumerate(vals):
        center = bottom[i] + value/2
        if value < 100:
            axes[1].annotate(f'{value:,}', xy=(i+.285, center),
                             xytext=(i+.43, center+(55 if value == 43 else -20)),
                             fontsize=9, va='center', color='#17212b',
                             arrowprops={'arrowstyle': '-', 'color': '#657383', 'lw': .7})
        else:
            axes[1].text(i, center, f'{value:,}', ha='center', va='center', fontsize=10,
                         color='white' if label == 'Shared significant' else '#17212b')
    bottom += vals
axes[1].set(ylabel='Number of genes', title='FDR < 0.05; no fold-change cutoff', ylim=(0, 1510), xlim=(-.65, 1.8))
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, fontsize=8, loc='lower left', bbox_to_anchor=(.56, .015))
fig.suptitle('GSE242230: saved blind Pantheon analysis vs original workflow', fontsize=13, y=.99)
fig.tight_layout(rect=(0, .15, 1, .94))
fig.savefig(out / 'comparison_overview.png', dpi=200)
fig.savefig(out / 'comparison_overview.pdf')
plt.close(fig)

versions = {'python': sys.version, 'platform': platform.platform(), 'numpy': np.__version__,
            'pandas': pd.__version__, 'scipy': scipy.__version__, 'matplotlib': matplotlib.__version__}
write_json(versions, 'comparison_software_versions.json')
for path, recorded in sources.items():
    p = Path(path)
    assert p.stat().st_mtime_ns == recorded['mtime_ns'] and sha(p) == recorded['sha256'], path
write_json({'inputs': sources, 'all_inputs_unchanged_after_comparison': True,
            'comparison_script_sha256': sha(Path(__file__).resolve())}, 'comparison_input_manifest.json')
print(json.dumps({k: summary[k] for k in ['original_significant', 'pantheon_significant', 'shared_significant',
     'original_only_significant', 'pantheon_only_significant', 'original_only_split_by_pantheon_status',
     'pantheon_only_split_by_original_status', 'effect_agreement', 'wu16_exclusion']}, indent=2))
