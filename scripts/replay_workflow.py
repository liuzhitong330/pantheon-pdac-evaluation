#!/usr/bin/env python3
"""Replay saved scientific code in a new directory; never call a model.

This Codex convenience wrapper is separate from the saved Pantheon code. It
copies that code unchanged, links locally supplied public data, and records
commands and source hashes. --help and --prepare-only do not run the analyses.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


PANTHEON_FILES = [
    "outputs/analysis/rerun_from_raw.sh",
    "outputs/analysis/cohort_specifications.tsv",
    "outputs/analysis/PLAN_AMENDMENT_2026-10-22.md",
] + ["outputs/analysis/code/" + name for name in (
    "rerun_inspection.py", "rerun_prepare_counts.py", "rerun_validate_and_review.py",
    "run_edger.R", "audit_fits.R", "summarize_results.py", "add_presentation_plots.R",
)]
REFERENCE_STAGES = [
    "01_load_epithelial.py", "02_qc_epithelial.py", "03_pseudobulk.py", "04_paired_de.py",
]
METADATA_FILES = ["GSE242230_sample_records.soft", "Table_S1_sample_metadata.tsv"]


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def executable(value, label):
    expanded = os.path.expanduser(value)
    found = shutil.which(expanded)
    if found is None:
        raise ValueError(f"{label} is not an executable file: {value}")
    # Do not resolve the final symlink: doing so can bypass a Python venv.
    return Path(os.path.abspath(found))


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Required nonempty file is missing: {path}")


def interpreter_shim(path, target):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nexec " + shlex.quote(str(target)) + ' "$@"\n')
    path.chmod(0o755)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", choices=["pantheon", "reference"])
    parser.add_argument("destination", type=Path, help="New scratch directory; existing paths are refused")
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Directory containing the 25 original matrix triplets and annotation file")
    parser.add_argument("--python", required=True,
                        help="Python executable in the matching scientific environment")
    parser.add_argument("--rscript", help="Rscript executable in the matching R environment; required for pantheon")
    parser.add_argument("--shell", default="zsh", help="Shell for Pantheon's saved entrypoint (default: zsh)")
    parser.add_argument("--prepare-only", action="store_true",
                        help="Validate paths and prepare copies/links; do not execute scientific code")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    destination = args.destination.expanduser().absolute()
    data = args.data_dir.expanduser().resolve()
    try:
        if destination.exists() or destination.is_symlink():
            raise ValueError("Destination must not already exist; saved results are never overwritten.")
        python = executable(args.python, "Python")
        rscript = shell = None
        if args.workflow == "pantheon":
            if not args.rscript:
                raise ValueError("--rscript is required for the pantheon workflow.")
            rscript = executable(args.rscript, "Rscript")
            shell = executable(args.shell, "Shell")
        elif args.rscript:
            raise ValueError("--rscript applies only to the pantheon workflow.")
        annotation = data / "GSE242230_annotations.txt.gz"
        require_file(annotation)
        matrices = sorted(data.glob("*_matrix.mtx.gz"))
        if len(matrices) != 25:
            raise ValueError(f"Expected 25 raw matrix files, found {len(matrices)} in {data}")
        data_files = [annotation]
        sample_ids = []
        for matrix in matrices:
            match = re.fullmatch(r"GSM\d+_(wu\d+)_matrix\.mtx\.gz", matrix.name)
            if not match:
                raise ValueError(f"Unexpected matrix filename: {matrix.name}")
            sample_ids.append(match.group(1))
            prefix = matrix.name.removesuffix("matrix.mtx.gz")
            for suffix in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
                path = data / (prefix + suffix)
                require_file(path)
                data_files.append(path)
        if len(set(sample_ids)) != 25:
            raise ValueError("Matrix sample identifiers are not unique.")
        relative_files = PANTHEON_FILES if args.workflow == "pantheon" else [
            "reproduction/scripts/" + name for name in REFERENCE_STAGES]
        source_root = repo / "workflows" / args.workflow
        for rel in relative_files:
            require_file(source_root / rel)
        if args.workflow == "pantheon":
            for name in METADATA_FILES:
                require_file(repo / "data/metadata" / name)
    except ValueError as exc:
        parser.error(str(exc))

    destination.mkdir(parents=True)
    (destination / "data").mkdir()
    (destination / "replay/logs").mkdir(parents=True)
    for path in data_files:
        (destination / "data" / path.name).symlink_to(path)
    if args.workflow == "pantheon":
        meta = destination / "data/sample_metadata"
        meta.mkdir()
        for name in METADATA_FILES:
            (meta / name).symlink_to(repo / "data/metadata" / name)
    hashes = {}
    for rel in relative_files:
        source = source_root / rel
        target = destination / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[rel] = sha256(source)
        if sha256(target) != hashes[rel]:
            raise RuntimeError(f"Copied source hash mismatch: {rel}")
    interpreter_shim(destination / ".analysis-venv/bin/python", python)
    if rscript:
        interpreter_shim(destination / ".r-env/bin/Rscript", rscript)
    manifest = {
        "wrapper_role": "Codex convenience wrapper; saved scientific code copied unchanged",
        "workflow": args.workflow, "created_utc": utc_now(), "status": "prepared",
        "source_root": str(source_root), "destination": str(destination),
        "python": str(python), "rscript": str(rscript) if rscript else None,
        "source_sha256": hashes,
        "data_links": [{"name": path.name, "source": str(path), "bytes": path.stat().st_size}
                       for path in data_files],
        "input_validation": "Names, count, uniqueness and nonempty files only; verify release checksums separately.",
        "commands": [],
    }
    manifest_path = destination / "replay/manifest.json"

    def save_manifest():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    save_manifest()
    if args.prepare_only:
        print(f"Prepared without executing analysis: {destination}")
        print("To execute, call this wrapper again with a different new destination and omit --prepare-only.")
        return

    env = dict(os.environ)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", TZ="UTC",
               OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2",
               NUMBA_NUM_THREADS="2", R_ENVIRON_USER="/dev/null", R_PROFILE_USER="/dev/null")
    for key, rel in [("MPLCONFIGDIR", ".cache/matplotlib"), ("NUMBA_CACHE_DIR", ".cache/numba"),
                     ("PIP_CACHE_DIR", ".cache/pip"), ("JOBLIB_TEMP_FOLDER", ".tmp/joblib")]:
        path = destination / rel
        path.mkdir(parents=True, exist_ok=True)
        env[key] = str(path)
    commands = ([[str(shell), "outputs/analysis/rerun_from_raw.sh", "outputs/replay"]]
                if args.workflow == "pantheon" else
                [[str(python), "reproduction/scripts/" + name] for name in REFERENCE_STAGES])
    # The recorded R installation requires a temporary directory without spaces.
    with tempfile.TemporaryDirectory(prefix="pdac-replay-", dir="/tmp") as tmp:
        env.update(TMPDIR=tmp, TMP=tmp, TEMP=tmp, TMPPREFIX=tmp + "/zsh")
        manifest.update(status="running", started_utc=utc_now(), temporary_directory=tmp)
        save_manifest()
        try:
            for index, command in enumerate(commands, 1):
                record = {"argv": command, "started_utc": utc_now()}
                manifest["commands"].append(record)
                save_manifest()
                stem = f"{index:02d}_{args.workflow}"
                with (destination / "replay/logs" / (stem + ".stdout.log")).open("w") as stdout, \
                     (destination / "replay/logs" / (stem + ".stderr.log")).open("w") as stderr:
                    result = subprocess.run(command, cwd=destination, env=env, stdout=stdout, stderr=stderr)
                record.update(exit_code=result.returncode, ended_utc=utc_now())
                save_manifest()
                if result.returncode:
                    raise subprocess.CalledProcessError(result.returncode, command)
            manifest["status"] = "completed"
        except BaseException as exc:
            manifest.update(status="failed", error=repr(exc))
            raise
        finally:
            manifest["ended_utc"] = utc_now()
            save_manifest()
    print(f"Saved-code replay completed: {destination}")


if __name__ == "__main__":
    main()
