# Provenance and responsibility

This repository curates a completed evaluation; it does not claim that publication reran the biological analysis. Existing source projects remain unchanged. Local absolute paths are replaced by relative source labels in public manifests.

| Component | Origin | Public location |
|---|---|---|
| Blind plan, count preparation, paired edgeR analysis, sensitivities and presentation | PantheonOS, with evaluator setup/QC assistance | [protocol](../protocol), [Pantheon workflow](../workflows/pantheon/outputs/analysis) |
| Existing four-script PyDESeq2 reference | Existing PDAC project; saved path-only copies used by Pantheon in a fresh reproduction session | [reference workflow](../workflows/reference/reproduction/scripts) |
| Gene overlap, effect comparison, filtering audit and comparison figure | Codex evaluator, after blind results were saved | [comparison script](../scripts/compare_saved_analyses.py), [results](../results), [figures](../figures) |
| Public packaging, unpacker and replay convenience wrapper | Codex curation for this repository | [scripts](../scripts), [environment](../environment), this documentation |

Scientific files copied here are identified in [copied_files.json](../provenance/copied_files.json) by original evaluation-relative location, source SHA-256 and public SHA-256. The comparison script has only three path-default changes, recorded in [comparison_script_provenance.json](../results/comparison_script_provenance.json); its statistical/comparison logic is preserved. The twenty bundled saved inputs retain their original bytes, documented in [input_manifest.json](../results/input_manifest.json). Fresh wrappers select runtimes and output directories without changing scientific settings.

The raw-data hashes are from the completed input audit. The repository's [SHA256SUMS.txt](../SHA256SUMS.txt) checks the curated payload. The [environment records](../environment) distinguish the observed successful local environments from cross-platform reproducibility claims.

Full agent conversations, OAuth/authentication material, private execution logs, installed runtimes, raw count datasets and unrelated EFS/network projects are deliberately outside this public repository. Safe saved plans and scientific result tables are included. The [methods record](METHODS.md) summarizes assistance and limitations without publishing private logs.

Previously saved verification found exact scientific agreement on replay: 372 first-round tables and ten required reference artifacts. Those checks establish computational consistency in the recorded environments. The new public wrappers have only been statically checked, and this curation did not rerun fits. The source-data checksum gap noted in the original preparation script is documented in [METHODS.md](METHODS.md), rather than silently altering the original workflow.

This is a personal evaluation using [PantheonOS](https://github.com/aristoteleo/PantheonOS), not an official Pantheon benchmark or a Qiu lab endorsement. AI-generated scientific code and evaluator work are explicitly attributed; this repository is not evidence of unassisted implementation or experimental validation.
