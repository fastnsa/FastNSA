"""Validate and assemble FastNSA release artifacts; this tool never uploads."""

from __future__ import annotations

import argparse
from email.parser import BytesParser
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import sys
import tarfile
import tomllib
import urllib.error
import urllib.request
import zipfile


TARGETS = {"linux-x86_64", "windows-amd64", "macos-x86_64", "macos-arm64"}
REQUIRED_RUNTIME = {
    "fastnsa/__init__.py", "fastnsa/runtime.py", "fastnsa/enums.py",
    "fastnsa/exceptions.py", "fastnsa/bindings/capi.py",
    "fastnsa/bindings/_lib.py", "fastnsa/stats/network_kde.py",
}
BUNDLE_FILES = {"manifest.json", "SHA256SUMS", "RELEASE_NOTES.md"}
STABLE_VERSION = re.compile(r"(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)")


class ReleaseError(ValueError):
    """A version, distribution, or release bundle is unsafe to publish."""


def require(condition, message):
    if not condition:
        raise ReleaseError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_version(root, ref_type="", ref_name=""):
    root = Path(root)
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    require(project.get("name") == "fastnsa", "Wrong project name")
    version = project["version"]
    require(STABLE_VERSION.fullmatch(version), "Only canonical X.Y.Z stable versions are supported")
    cmake = (root / "cpp/CMakeLists.txt").read_text(encoding="utf-8")
    match = re.search(r"project\(fastnsa\s+VERSION\s+([^\s)]+)", cmake)
    require(match and match.group(1) == version, "CMake and Python versions differ")
    if ref_type == "tag":
        require(ref_name == "v" + version, "Release tag must exactly match v" + version)
    notes = root / "python/docs/releases" / (version + ".md")
    require(notes.is_file(), "Missing release notes: " + str(notes))
    return project


def safe_member(name):
    path = PurePosixPath(name)
    require(name and not path.is_absolute() and ".." not in path.parts
            and "\\" not in name and ":" not in name, "Unsafe archive path: " + name)


def check_metadata(data, project):
    from packaging.specifiers import SpecifierSet

    metadata = BytesParser().parsebytes(data)
    require(metadata["Name"] == "fastnsa", "Wrong distribution project name")
    require(metadata["Version"] == project["version"], "Distribution version differs from project")
    require(SpecifierSet(metadata.get("Requires-Python", "")) ==
            SpecifierSet(project["requires-python"]), "Requires-Python differs from project")


def inspect_wheel(path, project):
    from packaging.utils import parse_wheel_filename

    name, version, _, tags = parse_wheel_filename(path.name)
    require(name == "fastnsa" and str(version) == project["version"], "Wrong wheel name/version")
    targets = set()
    for tag in tags:
        require(tag.interpreter == "py3" and tag.abi == "none", "Unexpected wheel Python/ABI tag")
        platform = tag.platform
        if platform in {"manylinux2014_x86_64", "manylinux_2_17_x86_64"}:
            targets.add("linux-x86_64")
        elif platform == "win_amd64":
            targets.add("windows-amd64")
        elif re.fullmatch(r"macosx_\d+_\d+_x86_64", platform):
            targets.add("macos-x86_64")
        elif re.fullmatch(r"macosx_\d+_\d+_arm64", platform):
            targets.add("macos-arm64")
        else:
            raise ReleaseError("Unsupported release wheel platform: " + platform)
    require(len(targets) == 1, "A wheel must identify one release target")
    target = targets.pop()
    library = ("fastnsa.dll" if target == "windows-amd64" else
               "libfastnsa.dylib" if target.startswith("macos") else "libfastnsa.so")
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)), "Duplicate wheel members")
        for entry in archive.infolist():
            safe_member(entry.filename)
            require(not stat.S_ISLNK(entry.external_attr >> 16), "Wheel symlinks are not allowed")
        require(REQUIRED_RUNTIME | {"fastnsa/bindings/" + library} <= set(names),
                "Wheel is missing the current native/Python runtime")
        metadata = [n for n in names if n.endswith(".dist-info/METADATA")]
        require(len(metadata) == 1, "Expected exactly one wheel METADATA file")
        check_metadata(archive.read(metadata[0]), project)
        require(archive.testzip() is None, "Corrupt wheel ZIP member")
    return target


def inspect_sdist(path, project):
    from packaging.utils import parse_sdist_filename

    name, version = parse_sdist_filename(path.name)
    require(name == "fastnsa" and str(version) == project["version"], "Wrong source package name/version")
    prefix = "fastnsa-" + project["version"] + "/"
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        names = [m.name for m in members]
        require(len(names) == len(set(names)), "Duplicate source package members")
        for member in members:
            safe_member(member.name)
            require(member.name.startswith(prefix) or member.name == prefix[:-1],
                    "Unexpected source package root")
            require(member.isfile() or member.isdir(), "Source package links/devices are not allowed")
            require(not member.name.endswith((".so", ".dll", ".dylib")) and ".so." not in member.name,
                    "Source package contains a native binary")
        required = {"PKG-INFO", "pyproject.toml", "cpp/CMakeLists.txt",
                    "cpp/bindings/fastnsa_capi.cpp", "cpp/bindings/fastnsa_capi.h",
                    "python/LICENSE", "python/README.md", "python/fastnsa/runtime.py"}
        require({prefix + n for n in required} <= set(names), "Incomplete source package")
        check_metadata(archive.extractfile(prefix + "PKG-INFO").read(), project)


def assemble(root, artifacts, output, commit, ref):
    project = check_version(root)
    require(re.fullmatch(r"[0-9a-f]{40}", commit), "Expected a full Git commit SHA")
    require(ref.startswith(("refs/heads/", "refs/tags/")), "Expected a Git branch or tag ref")
    if ref.startswith("refs/tags/"):
        check_version(root, "tag", ref.removeprefix("refs/tags/"))
    artifacts, output = Path(artifacts), Path(output)
    require(artifacts.is_dir() and not artifacts.is_symlink(), "Missing artifact directory")
    files = []
    for path in artifacts.rglob("*"):
        require(not path.is_symlink(), "Artifact symlinks are not allowed")
        if path.is_file():
            require(path.name.endswith((".whl", ".tar.gz")), "Unexpected artifact file: " + path.name)
            files.append(path)
    require(len(files) == 5, "Expected exactly four platform wheels and one source package")
    records, targets, names = [], set(), set()
    for path in sorted(files):
        require(path.name not in names, "Duplicate artifact filename")
        names.add(path.name)
        if path.suffix == ".whl":
            target = inspect_wheel(path, project)
            require(target not in targets, "Duplicate wheel target: " + target)
            targets.add(target)
        else:
            inspect_sdist(path, project)
            target = "sdist"
        records.append({"name": path.name, "target": target, "sha256": sha256(path),
                        "size": path.stat().st_size})
    require(targets == TARGETS and sum(r["target"] == "sdist" for r in records) == 1,
            "Missing release platform/source package")
    require(not output.exists() and not output.is_symlink(), "Refusing to overwrite a release bundle")
    output.mkdir(parents=True)
    (output / "dist").mkdir()
    for path in files:
        shutil.copyfile(path, output / "dist" / path.name)
    manifest = {"schema_version": 1, "package": "fastnsa", "version": project["version"],
                "source_commit": commit, "source_ref": ref, "files": records,
                "release_notes_sha256": sha256(Path(root) / "python/docs/releases" / (project["version"] + ".md"))}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output / "SHA256SUMS").write_text(
        "".join(r["sha256"] + "  " + r["name"] + "\n" for r in records), encoding="ascii"
    )
    shutil.copyfile(Path(root) / "python/docs/releases" / (project["version"] + ".md"),
                    output / "RELEASE_NOTES.md")
    verify(output, project["version"], commit, ref)
    return manifest


def verify(bundle, version, commit, ref):
    """Verify handoff integrity using only the standard library, without building code."""
    root = Path(bundle)
    require(re.fullmatch(r"[0-9a-f]{40}", commit), "Expected a full Git commit SHA")
    require(ref.startswith(("refs/heads/", "refs/tags/")), "Expected a Git branch or tag ref")
    require(root.is_dir() and not root.is_symlink(), "Missing release bundle")
    for path in root.rglob("*"):
        require(not path.is_symlink(), "Release bundle symlinks are not allowed")
        require(path.is_file() or path.is_dir(), "Unexpected release bundle entry")
    require({p.name for p in root.iterdir()} == BUNDLE_FILES | {"dist"},
            "Unexpected release bundle root files")
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    require(isinstance(manifest, dict), "Malformed release manifest")
    require(manifest.get("schema_version") == 1 and manifest.get("package") == "fastnsa",
            "Unsupported release manifest")
    require(manifest.get("version") == version and STABLE_VERSION.fullmatch(version),
            "Release manifest version mismatch")
    require(manifest.get("source_commit") == commit and manifest.get("source_ref") == ref,
            "Release bundle commit/ref mismatch")
    if ref.startswith("refs/tags/"):
        require(ref == "refs/tags/v" + version, "Release bundle tag mismatch")
    records = manifest.get("files")
    require(isinstance(records, list) and len(records) == 5,
            "Expected five release distributions")
    for record in records:
        require(isinstance(record, dict) and
                all(isinstance(record.get(k), str) for k in ("name", "target", "sha256")) and
                type(record.get("size")) is int and record["size"] > 0 and
                re.fullmatch(r"[0-9a-f]{64}", record["sha256"]), "Malformed distribution record")
    require({r["target"] for r in records} == TARGETS | {"sdist"}, "Wrong platform coverage")
    names = {r["name"] for r in records}
    require(len(names) == 5 and {p.name for p in (root / "dist").iterdir()} == names,
            "Unexpected, missing, or duplicate distribution files")
    for record in records:
        name = record["name"]
        require(PurePosixPath(name).name == name and "\\" not in name, "Unsafe distribution filename")
        require(name.startswith("fastnsa-" + version + "-") or name == "fastnsa-" + version + ".tar.gz",
                "Distribution filename version mismatch")
        path = root / "dist" / name
        require(path.is_file() and path.stat().st_size == record["size"] and
                sha256(path) == record["sha256"], "Distribution checksum/size mismatch: " + name)
    expected = "".join(r["sha256"] + "  " + r["name"] + "\n" for r in records)
    require((root / "SHA256SUMS").read_text(encoding="ascii") == expected, "Checksum list mismatch")
    require((root / "RELEASE_NOTES.md").is_file() and
            sha256(root / "RELEASE_NOTES.md") == manifest.get("release_notes_sha256"),
            "Missing or modified release notes")
    return manifest


def check_index(version):
    """Fail closed on existing releases or network errors; never upload or delete anything."""
    require(STABLE_VERSION.fullmatch(version), "Invalid release version")
    request = urllib.request.Request("https://pypi.org/pypi/fastnsa/" + version + "/json",
                                     headers={"User-Agent": "FastNSA-release-preflight"})
    try:
        with urllib.request.urlopen(request, timeout=30):
            pass
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return
        raise ReleaseError("PyPI availability check failed: HTTP " + str(error.code)) from error
    raise ReleaseError("PyPI version already exists; inspect partial/previous uploads before retrying")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    version = sub.add_parser("version")
    version.add_argument("--root", type=Path, default=Path("."))
    version.add_argument("--github-output", type=Path)
    assembly = sub.add_parser("assemble")
    assembly.add_argument("--root", type=Path, default=Path("."))
    assembly.add_argument("--artifacts", type=Path, required=True)
    assembly.add_argument("--output", type=Path, required=True)
    verification = sub.add_parser("verify")
    verification.add_argument("--bundle", type=Path, required=True)
    verification.add_argument("--version", required=True)
    for command in (assembly, verification):
        command.add_argument("--commit", required=True)
        command.add_argument("--ref", required=True)
    index = sub.add_parser("check-index")
    index.add_argument("--version", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "version":
            project = check_version(args.root, os.environ.get("GITHUB_REF_TYPE", ""),
                                    os.environ.get("GITHUB_REF_NAME", ""))
            values = {"version": project["version"], "tag": "v" + project["version"]}
            if args.github_output:
                with args.github_output.open("a", encoding="utf-8") as output:
                    output.write("".join(k + "=" + v + "\n" for k, v in values.items()))
            print(json.dumps(values))
        elif args.command == "assemble":
            print(json.dumps(assemble(args.root, args.artifacts, args.output, args.commit, args.ref), indent=2))
        elif args.command == "verify":
            verify(args.bundle, args.version, args.commit, args.ref)
            print("Release bundle verified.")
        else:
            check_index(args.version)
            print("PyPI version is available.")
        return 0
    except (ReleaseError, OSError, ValueError, KeyError, zipfile.BadZipFile, tarfile.TarError) as error:
        print("Release check failed: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
