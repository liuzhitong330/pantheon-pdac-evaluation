# Replaying the saved scientific workflows

These instructions replay the saved code, including its fixed cohort decisions and statistical settings. They do not start Pantheon, call a language model, repeat the blind method-selection exercise, or require model credentials. The public convenience wrapper was written by Codex; the scientific scripts under `workflows/` were preserved from the evaluation.

The saved evaluation included successful local replays. This public package has received static checks and a `--help` check only. **No fresh environment installation or scientific replay was run during public curation.** The recorded platform was macOS 14.2.1 on Apple Silicon. Numerical identity on other operating systems, package builds, or BLAS implementations has not been tested.

Run commands below from the repository root. Use Python 3.9 or newer for the convenience wrapper. It refuses an existing destination and copies the scientific code unchanged into a new scratch directory. Source hashes, commands, exit status, and scientific execution logs are written under that directory's `replay/` folder. The wrapper links data files; it does not modify them or copy the original project.

## Reproduce the comparison from saved results

This route only reads the twenty bundled derived inputs; it fits no DE model.
The [actual comparison script](../scripts/compare_saved_analyses.py) is preserved
except for three public path defaults. Its saved runtime was Python 3.9.6.

```sh
python3.9 -m venv .venvs/comparison
.venvs/comparison/bin/python -m pip install -r environment/comparison-requirements.txt
python3 scripts/unpack_saved_results.py --destination work/saved
.venvs/comparison/bin/python scripts/compare_saved_analyses.py \
  --evaluation work/saved/evaluation --original work/saved/reference \
  --output work/comparison
```

The unpacker verifies the archive and every member against
[input_manifest.json](../results/input_manifest.json) and refuses an existing
destination. It materializes about 86 MB of saved results and pseudobulks from
the 13 MB archive. The comparison writes complete gene tables, filtering and
sample-inclusion checks, sensitivity comparisons and the figure to
`work/comparison`. It asserts exact saved count equality. Compare its summary
with [comparison_summary.json](../results/comparison_summary.json); plots can
have binary/rendering differences across platforms. Figure labels are specific
to this frozen evaluation; this is not a generic DE comparison package.

Neither this comparison nor either scientific workflow was rerun for public
curation. The unpacker was checked with temporary extraction; the scientific
scripts were inspected without executing them.

## Inputs for a full raw-data replay

Obtain the original public [GSE242230 data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242230). Place the decompressed archive's **compressed member files** together in one directory:

- 25 `GSM…_wu…_matrix.mtx.gz` files;
- each matrix's matching `_features.tsv.gz` and `_barcodes.tsv.gz` files;
- `GSE242230_annotations.txt.gz`.

Do not rename sample prefixes or decompress the individual `.gz` files. There should be 76 files needed by the replay: 75 matrix triplet members and one annotation file. Verify their recorded input checksums separately. The wrapper checks names, file count, unique sample IDs, and nonempty files; it does not establish byte identity by those checks alone.

The Pantheon workflow also reads the two supplied public metadata files in this repository:

- `data/metadata/GSE242230_sample_records.soft`
- `data/metadata/Table_S1_sample_metadata.tsv`

The wrapper exposes these at `data/sample_metadata/` in its scratch directory, matching the saved code. The reference workflow does not consume these metadata files. Both workflows consume the same count matrices and annotations.

## Recorded environments

`environment/observed-versions.json` distinguishes recorded versions from testing performed for this public package. The dependency specifications contain no installed runtime, authentication, or private machine paths.

| Workflow | Recorded scientific runtime |
| --- | --- |
| Pantheon saved code | Python 3.12.14; NumPy 2.5.4, pandas 3.0.6, SciPy 1.18.1; R 4.4.2, Bioconductor 3.20, edgeR 4.4.0, limma 3.62.0, statmod 1.5.0 |
| Reference saved code | Python 3.9.6; 43 package pins, including PyDESeq2 0.4.12, Scanpy 1.10.3, AnnData 0.10.9, NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1 |

The reference package snapshot came from the original project's environment at inspection time. It is not a complete contemporaneous lock from the historical analysis date. Preserve that limitation when describing reproduction.

Using Python executables of the recorded versions, create separate environments:

```sh
# Verify these report 3.12.14 and 3.9.6 respectively before proceeding.
python3.12 --version
python3.9 --version
python3.12 -m venv .venvs/pantheon
.venvs/pantheon/bin/python -m pip install -r environment/pantheon-python-requirements.txt
python3.9 -m venv .venvs/reference
.venvs/reference/bin/python -m pip install -r environment/reference-python-requirements.txt
```

These commands document an intended clean setup; they were not executed for this release. A package version alone does not fix the interpreter build, compiler, or numerical libraries. Retain installation output and the resulting environment inventory when performing a new replay. Do not silently substitute versions to obtain agreement.

## Exact R package sources

`environment/r-osx-arm64-explicit.txt` is the observed explicit conda package/build lock. It is **specific to macOS arm64** and does not define a Linux environment. The observed micromamba version was 2.9.0; see the [official micromamba installation instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

Use a new R prefix without spaces. The observed conda R installation could not compile reliably in a prefix containing spaces, and required a no-space temporary directory. A compatible compiler toolchain is needed for the source packages.

```sh
# Replace this with your own new absolute path, containing no spaces.
export R_PREFIX="/absolute/no-space/path/to/r-runtime"
micromamba create --yes --prefix "$R_PREFIX" --file environment/r-osx-arm64-explicit.txt
```

The explicit environment supplies R and its resolved dependencies. Install the exact additional source packages in this order: **BiocVersion 3.20.0, statmod 1.5.0, limma 3.62.0, edgeR 4.4.0**. The official package archive/mirror URLs and observed SHA-256 digests are in `environment/r-source-archives.json`. This avoids a current release resolver choosing a newer patch version. The observed transitive R packages were locfit 1.5-9.12, Rcpp 1.1.2, and lattice 0.23-1; the complete inventory is `environment/r-package-manifest.tsv`.

The following setup command downloads only those four recorded archives, verifies each digest, and installs them. It does not read dataset files or run an analysis:

```sh
python3 - "$R_PREFIX" <<'PY'
from pathlib import Path
import hashlib, json, os, subprocess, sys, tempfile, urllib.request

prefix = Path(sys.argv[1]).expanduser().absolute()
if " " in str(prefix):
    raise SystemExit("Choose an R prefix without spaces.")
sources = json.loads(Path("environment/r-source-archives.json").read_text())
downloads = Path(tempfile.mkdtemp(prefix="pdac-r-sources-", dir="/tmp"))
env = dict(os.environ)
env.update(TMPDIR=str(downloads), R_ENVIRON_USER="/dev/null",
           R_PROFILE_USER="/dev/null", R_LIBS_USER=str(prefix / "lib/R/library"))
for source in sources:
    archive = downloads / source["filename"]
    blob = urllib.request.urlopen(source["url"], timeout=60).read()
    if hashlib.sha256(blob).hexdigest() != source["sha256"]:
        raise SystemExit("Source checksum mismatch: " + source["filename"])
    archive.write_bytes(blob)
    subprocess.run(["micromamba", "run", "--prefix", str(prefix),
                    "R", "CMD", "INSTALL", "--no-multiarch", str(archive)],
                   env=env, check=True)
print("Verified source archives retained in", downloads)
PY
```

The saved R analysis itself asserts R 4.4.2, edgeR 4.4.0, limma 3.62.0, statmod 1.5.0, and BiocVersion 3.20.0 before fitting. Do not remove those checks to force a run with different versions.

## Replay commands

Inspect the wrapper safely:

```sh
python3 scripts/replay_workflow.py --help
```

For the Pantheon workflow, supply the matching Python and Rscript executables:

```sh
python3 scripts/replay_workflow.py pantheon runs/pantheon-replay \
  --data-dir /absolute/path/to/GSE242230-files \
  --python "$PWD/.venvs/pantheon/bin/python" \
  --rscript "$R_PREFIX/bin/Rscript"
```

This preserves the saved entrypoint `outputs/analysis/rerun_from_raw.sh`, run with `zsh`. Use `--shell /path/to/zsh` if needed. The wrapper supplies executable shims at the two runtime locations expected by that entrypoint; it does not rewrite its scientific code. The wrapper creates a temporary no-space directory and sets numerical thread limits to two.

The entrypoint runs inspection, raw-count preparation, validation, edgeR fitting, fitted-result auditing, summaries, and plots. Its output is `runs/pantheon-replay/outputs/replay/`. It uses the frozen `cohort_specifications.tsv` and plan amendment, including the seven-patient primary analysis, eight-patient QC sensitivity, and other recorded sensitivity/influence analyses. Leave-one-patient-out fits retain their diagnostic role rather than creating new significance families.

For the reference workflow:

```sh
python3 scripts/replay_workflow.py reference runs/reference-replay \
  --data-dir /absolute/path/to/GSE242230-files \
  --python "$PWD/.venvs/reference/bin/python"
```

This runs the four saved `01`–`04` scripts in order. Their preserved `reproduction/scripts/` depth determines the scratch output root. H5AD files, pseudobulk tables, and the `pdac_de_*.csv` results are written directly under `runs/reference-replay/`. The public wrapper invokes those scripts directly and does not include the evaluation's separate fitted-object capture harness. Its commands, source hashes, and stdout/stderr are still recorded under `replay/`.

Add `--prepare-only` to either command to validate required paths and create the copies/links without executing scientific code. A prepared destination cannot subsequently be reused by this wrapper: choose another new destination for execution. Neither mode invokes a model or changes methodological settings.

## Comparing a new replay

Compare cell identifiers and membership, sample/patient labels, pseudobulk counts, gene identifiers and masks, design and contrast, then aligned fold changes and adjusted p-values. Excluded genes must remain distinguishable from tested nonsignificant genes. Compare tables by identifiers and column names; timestamps, host paths, compressed-container bytes, and diagnostic-column order need not match even when scientific values do. A changed platform or dependency may cause numerical differences; investigate before claiming exact reproduction.

The original Pantheon preparation script has a recorded provenance gap: its checksum scan names `outputs/inspection` and `outputs/planning_revision`, although scientific reads correctly use the fresh `INSPECTION_OUTPUT` location. This is preserved in the scientific source. The wrapper records the copied code hashes and linked input locations; it does not silently repair that checksum scope. Verify input byte checksums separately.

The saved `validate_reproduction.py` depends on the original evaluation's particular output/checksum layout and is not this public wrapper's validation command. Full source/environment provenance and external replay validation remain separate from Pantheon's own scientific outputs.
