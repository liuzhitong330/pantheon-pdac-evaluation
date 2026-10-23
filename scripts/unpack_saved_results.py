#!/usr/bin/env python3
"""Verify and unpack the frozen comparison inputs into a new directory.

This copies saved tables only. It does not run Pantheon or fit any model.
The destination must not already exist, including as a symbolic link.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import re
import shutil
import sys
import tarfile


def sha256_file(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def validate_member_name(name):
    if not isinstance(name, str) or "\\" in name or "\x00" in name:
        raise ValueError("Invalid archive member name")
    path = PurePosixPath(name)
    if (path.is_absolute() or not path.parts or
            any(part in (".", "..", "") or ":" in part for part in path.parts) or
            path.as_posix() != name or
            path.parts[0] not in ("evaluation", "reference")):
        raise ValueError("Unsafe archive member name: " + repr(name))
    return path.parts


def unpack(archive_path, manifest_path, destination):
    # Do not resolve the destination: resolving could conceal a dangling symlink.
    destination = destination.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Destination already exists; choose a new directory")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries = manifest["members"]
    if not isinstance(entries, list) or not entries:
        raise ValueError("Manifest must contain a nonempty members list")
    expected = {}
    for entry in entries:
        name = entry["member"]
        validate_member_name(name)
        if name in expected:
            raise ValueError("Duplicate manifest member: " + name)
        if (type(entry["bytes"]) is not int or entry["bytes"] < 0 or
                not isinstance(entry["sha256"], str) or
                not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"])):
            raise ValueError("Invalid size or checksum for " + name)
        expected[name] = entry
    archive_info = manifest["archive"]
    if archive_path.stat().st_size != archive_info["bytes"]:
        raise ValueError("Archive byte size does not match the manifest")
    if sha256_file(archive_path) != archive_info["sha256"]:
        raise ValueError("Archive checksum does not match the manifest")

    with tarfile.open(archive_path, mode="r:gz") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate archive members are not allowed")
        if set(names) != set(expected):
            raise ValueError("Archive members differ from the manifest")
        for member in members:
            validate_member_name(member.name)
            if not member.isfile() or member.issym() or member.islnk():
                raise ValueError("Only regular files are allowed: " + member.name)
            if member.size != expected[member.name]["bytes"]:
                raise ValueError("Member byte size differs: " + member.name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.mkdir()  # Exclusive creation: never reuse an existing directory.
        try:
            for member in members:
                target = destination.joinpath(*validate_member_name(member.name))
                target.parent.mkdir(parents=True, exist_ok=True)
                checksum = hashlib.sha256()
                written = 0
                with archive.extractfile(member) as source, target.open("xb") as output:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        written += len(chunk)
                        if written > member.size:
                            raise ValueError("Member exceeds declared byte size")
                        output.write(chunk)
                        checksum.update(chunk)
                if (written != member.size or
                        checksum.hexdigest() != expected[member.name]["sha256"]):
                    raise ValueError("Member checksum differs: " + member.name)
        except BaseException:
            shutil.rmtree(destination)
            raise
    return len(expected), sum(entry["bytes"] for entry in entries)


def main():
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path,
                        default=repository / "results/comparison_inputs.tar.gz")
    parser.add_argument("--manifest", type=Path,
                        default=repository / "results/input_manifest.json")
    parser.add_argument("--destination", type=Path,
                        default=repository / "work/saved")
    args = parser.parse_args()
    try:
        count, size = unpack(args.archive, args.manifest, args.destination)
    except (OSError, ValueError, KeyError, TypeError, tarfile.TarError) as error:
        parser.exit(1, "Unpack failed: " + str(error) + "\n")
    print("Verified and unpacked {} saved files ({} bytes).".format(count, size))
    print("No comparison or analysis was run.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
