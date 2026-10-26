"""Release validation tests: no upload, credentials, or native library needed."""

from contextlib import redirect_stdout, redirect_stderr
import importlib.util
import io
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import urllib.error
import zipfile


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("release_tool", ROOT / "tools/ci/release.py")
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)
COMMIT = "a" * 40
REF = "refs/tags/v0.3.0"
PLATFORMS = {
    "linux-x86_64": "manylinux2014_x86_64.manylinux_2_17_x86_64",
    "windows-amd64": "win_amd64",
    "macos-x86_64": "macosx_10_13_x86_64",
    "macos-arm64": "macosx_11_0_arm64",
}


class ReleaseToolTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="fastnsa-release-test-")
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.root = self.base / "project"
        self.root.mkdir()
        (self.root / "cpp").mkdir()
        (self.root / "python/docs/releases").mkdir(parents=True)
        self.set_version("0.3.0")
        (self.root / "python/docs/releases/0.3.0.md").write_text("# FastNSA 0.3.0\n")
        self.artifacts = self.base / "artifacts"
        self.artifacts.mkdir()
        self.bundle = self.base / "bundle"
        for platform in PLATFORMS.values():
            self.make_wheel(platform)
        self.make_sdist()

    def set_version(self, version):
        (self.root / "pyproject.toml").write_text(
            '[project]\nname = "fastnsa"\nversion = "' + version + '"\nrequires-python = ">=3.11"\n'
        )
        (self.root / "cpp/CMakeLists.txt").write_text("project(fastnsa VERSION " + version + ")\n")

    @staticmethod
    def metadata(version="0.3.0", python=">=3.11"):
        return ("Metadata-Version: 2.4\nName: fastnsa\nVersion: " + version +
                "\nRequires-Python: " + python + "\n").encode()

    def make_wheel(self, platform, version="0.3.0", metadata_version=None, python=">=3.11",
                   omit=None, extra=None):
        path = self.artifacts / ("fastnsa-" + version + "-py3-none-" + platform + ".whl")
        library = ("fastnsa.dll" if platform == "win_amd64" else
                   "libfastnsa.dylib" if platform.startswith("macos") else "libfastnsa.so")
        with zipfile.ZipFile(path, "w") as archive:
            for name in sorted(release.REQUIRED_RUNTIME | {"fastnsa/bindings/" + library}):
                if name != omit:
                    archive.writestr(name, b"test-only placeholder")
            archive.writestr("fastnsa-" + version + ".dist-info/METADATA",
                             self.metadata(metadata_version or version, python))
            if extra is not None:
                archive.writestr(extra, b"unexpected")
        return path

    def make_sdist(self, extra=None, omit=None):
        path = self.artifacts / "fastnsa-0.3.0.tar.gz"
        names = {"PKG-INFO", "pyproject.toml", "cpp/CMakeLists.txt",
                 "cpp/bindings/fastnsa_capi.cpp", "cpp/bindings/fastnsa_capi.h",
                 "python/LICENSE", "python/README.md", "python/fastnsa/runtime.py"}
        with tarfile.open(path, "w:gz") as archive:
            for name in sorted(names | ({extra} if extra else set())):
                if name == omit:
                    continue
                data = self.metadata() if name == "PKG-INFO" else b"source fixture"
                info = tarfile.TarInfo("fastnsa-0.3.0/" + name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
        return path

    def assemble(self, ref=REF):
        return release.assemble(self.root, self.artifacts, self.bundle, COMMIT, ref)

    def verify(self, **changes):
        args = {"version": "0.3.0", "commit": COMMIT, "ref": REF}
        args.update(changes)
        return release.verify(self.bundle, **args)

    def replace_linux(self, platform, **kwargs):
        (self.artifacts / ("fastnsa-0.3.0-py3-none-" + PLATFORMS["linux-x86_64"] + ".whl")).unlink()
        return self.make_wheel(platform, **kwargs)

    def test_matching_tag_and_branch_dry_run(self):
        self.assertEqual(release.check_version(self.root, "tag", "v0.3.0")["version"], "0.3.0")
        self.assertEqual(release.check_version(self.root, "branch", "main")["version"], "0.3.0")

    def test_mismatched_tag_rejected(self):
        for tag in ("v0.2.0", "0.3.0", "v0.3.0-rc1"):
            with self.subTest(tag=tag), self.assertRaisesRegex(release.ReleaseError, "exactly match"):
                release.check_version(self.root, "tag", tag)

    def test_noncanonical_and_prerelease_versions_rejected(self):
        for version in ("0.3.0rc1", "0.3.0.dev0", "0.3.0+local", "00.3.0", "0.3"):
            self.set_version(version)
            with self.subTest(version=version), self.assertRaisesRegex(release.ReleaseError, "canonical"):
                release.check_version(self.root)

    def test_cmake_mismatch_rejected(self):
        (self.root / "cpp/CMakeLists.txt").write_text("project(fastnsa VERSION 0.2.0)")
        with self.assertRaisesRegex(release.ReleaseError, "CMake"):
            release.check_version(self.root)

    def test_missing_notes_rejected(self):
        (self.root / "python/docs/releases/0.3.0.md").unlink()
        with self.assertRaisesRegex(release.ReleaseError, "notes"):
            release.check_version(self.root)

    def test_version_cli_writes_only_validated_outputs(self):
        output = self.base / "outputs"
        with patch.dict(os.environ, {"GITHUB_REF_TYPE": "branch", "GITHUB_REF_NAME": "main"}), \
                redirect_stdout(io.StringIO()):
            self.assertEqual(release.main(["version", "--root", str(self.root),
                                          "--github-output", str(output)]), 0)
        self.assertEqual(output.read_text(), "version=0.3.0\ntag=v0.3.0\n")

    def test_complete_bundle_and_hashes(self):
        manifest = self.assemble()
        self.assertEqual({r["target"] for r in manifest["files"]}, release.TARGETS | {"sdist"})
        self.assertEqual(self.verify(), manifest)

    def test_dry_run_bundle_is_not_a_tag_bundle(self):
        self.assemble(ref="refs/heads/main")
        release.verify(self.bundle, "0.3.0", COMMIT, "refs/heads/main")
        with self.assertRaisesRegex(release.ReleaseError, "commit/ref"):
            self.verify()

    def test_missing_platform_rejected(self):
        next(self.artifacts.glob("*win_amd64.whl")).unlink()
        with self.assertRaisesRegex(release.ReleaseError, "exactly four"):
            self.assemble()

    def test_duplicate_platform_rejected(self):
        next(self.artifacts.glob("*arm64.whl")).unlink()
        self.make_wheel("macosx_11_0_x86_64")
        with self.assertRaisesRegex(release.ReleaseError, "Duplicate wheel target"):
            self.assemble()

    def test_universal_and_unrepaired_linux_wheels_rejected(self):
        for platform in ("any", "linux_x86_64"):
            with self.subTest(platform=platform):
                path = self.replace_linux(platform)
                with self.assertRaisesRegex(release.ReleaseError, "Unsupported"):
                    self.assemble()
                path.unlink()
                self.make_wheel(PLATFORMS["linux-x86_64"])

    def test_old_version_rejected(self):
        self.replace_linux(PLATFORMS["linux-x86_64"], version="0.2.0")
        with self.assertRaisesRegex(release.ReleaseError, "name/version"):
            self.assemble()

    def test_wheel_metadata_mismatch_rejected(self):
        for values in ({"python": ">=3.9"}, {"metadata_version": "0.3.1"}):
            with self.subTest(values=values):
                self.make_wheel(PLATFORMS["linux-x86_64"], **values)
                with self.assertRaises(release.ReleaseError):
                    self.assemble()

    def test_incomplete_runtime_rejected(self):
        self.make_wheel(PLATFORMS["linux-x86_64"], omit="fastnsa/runtime.py")
        with self.assertRaisesRegex(release.ReleaseError, "runtime"):
            self.assemble()

    def test_wheel_traversal_and_symlinks_rejected(self):
        for extra in ("../escape", "/absolute"):
            self.make_wheel(PLATFORMS["linux-x86_64"], extra=extra)
            with self.subTest(extra=extra), self.assertRaisesRegex(release.ReleaseError, "Unsafe"):
                self.assemble()
        entry = zipfile.ZipInfo("fastnsa/linked")
        entry.create_system = 3
        entry.external_attr = (stat.S_IFLNK | 0o777) << 16
        self.make_wheel(PLATFORMS["linux-x86_64"], extra=entry)
        with self.assertRaisesRegex(release.ReleaseError, "symlinks"):
            self.assemble()

    def test_source_traversal_binary_and_missing_header_rejected(self):
        for kwargs in ({"extra": "../escape"}, {"extra": "libfastnsa.so"},
                       {"omit": "cpp/bindings/fastnsa_capi.h"}):
            self.make_sdist(**kwargs)
            with self.subTest(kwargs=kwargs), self.assertRaises(release.ReleaseError):
                self.assemble()

    def test_unexpected_artifact_rejected(self):
        (self.artifacts / "untrusted.exe").write_bytes(b"not a release file")
        with self.assertRaisesRegex(release.ReleaseError, "Unexpected artifact"):
            self.assemble()

    def test_existing_bundle_not_overwritten(self):
        self.bundle.mkdir()
        sentinel = self.bundle / "keep"
        sentinel.write_text("keep")
        with self.assertRaisesRegex(release.ReleaseError, "overwrite"):
            self.assemble()
        self.assertEqual(sentinel.read_text(), "keep")

    def test_tampered_distribution_rejected(self):
        manifest = self.assemble()
        path = self.bundle / "dist" / manifest["files"][0]["name"]
        path.write_bytes(path.read_bytes() + b"tampered")
        with self.assertRaisesRegex(release.ReleaseError, "checksum/size"):
            self.verify()

    def test_malformed_manifest_and_distribution_records_rejected(self):
        manifest = self.assemble()
        path = self.bundle / "manifest.json"
        for broken in ([], {**manifest, "files": None}, {**manifest, "files": [None] * 5}):
            path.write_text(json.dumps(broken))
            with self.subTest(broken=broken), self.assertRaises(release.ReleaseError):
                self.verify()

    def test_duplicate_artifact_names_in_separate_downloads_rejected(self):
        original = next(self.artifacts.glob("*arm64.whl"))
        duplicate = self.artifacts / "second-download" / original.name
        duplicate.parent.mkdir()
        original.rename(duplicate)
        next(self.artifacts.glob("*win_amd64.whl")).unlink()
        misplaced = self.artifacts / "third-download"
        misplaced.mkdir()
        linux = next(self.artifacts.glob("*manylinux*.whl"))
        (misplaced / linux.name).write_bytes(linux.read_bytes())
        with self.assertRaisesRegex(release.ReleaseError, "Duplicate artifact filename"):
            self.assemble()

    def test_commit_ref_and_version_mismatch_rejected(self):
        self.assemble()
        for values in ({"commit": "b" * 40}, {"ref": "refs/heads/main"}, {"version": "0.4.0"}):
            with self.subTest(values=values), self.assertRaises(release.ReleaseError):
                self.verify(**values)

    def test_extra_bundle_file_and_modified_notes_rejected(self):
        self.assemble()
        extra = self.bundle / "dist/untrusted.whl"
        extra.write_bytes(b"extra")
        with self.assertRaises(release.ReleaseError):
            self.verify()
        extra.unlink()
        (self.bundle / "RELEASE_NOTES.md").write_text("altered")
        with self.assertRaisesRegex(release.ReleaseError, "release notes"):
            self.verify()

    def test_verifier_needs_no_third_party_packages(self):
        self.assemble()
        result = subprocess.run(
            [sys.executable, "-B", "-I", "-S", str(ROOT / "tools/ci/release.py"), "verify",
             "--bundle", str(self.bundle), "--version", "0.3.0", "--commit", COMMIT, "--ref", REF],
            capture_output=True, text=True, timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_bundle_symlink_rejected(self):
        self.assemble()
        path = self.bundle / "link"
        try:
            path.symlink_to(self.root, target_is_directory=True)
        except OSError as error:
            self.skipTest("Symlinks unavailable: " + str(error))
        with self.assertRaisesRegex(release.ReleaseError, "symlinks"):
            self.verify()

    def test_existing_index_version_rejected(self):
        with patch.object(release.urllib.request, "urlopen"):
            with self.assertRaisesRegex(release.ReleaseError, "already exists"):
                release.check_index("0.3.0")

    def test_only_index_404_is_accepted(self):
        for code in (404, 403, 500):
            error = urllib.error.HTTPError("https://pypi.org/example", code, "test", {}, None)
            with self.subTest(code=code), patch.object(release.urllib.request, "urlopen", side_effect=error):
                if code == 404:
                    release.check_index("0.3.0")
                else:
                    with self.assertRaises(release.ReleaseError):
                        release.check_index("0.3.0")

    def test_index_network_error_fails_closed(self):
        with patch.object(release.urllib.request, "urlopen", side_effect=urllib.error.URLError("offline")), \
                redirect_stderr(io.StringIO()):
            self.assertEqual(release.main(["check-index", "--version", "0.3.0"]), 1)


class ReleaseWorkflowTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import yaml

        cls.workflow = yaml.load((ROOT / ".github/workflows/publish.yml").read_text(), Loader=yaml.BaseLoader)
        cls.matrix = yaml.load((ROOT / ".github/workflows/multiplatform.yml").read_text(), Loader=yaml.BaseLoader)

    def test_only_tag_push_can_publish(self):
        self.assertEqual(set(self.workflow["on"]), {"push", "workflow_dispatch"})
        self.assertEqual(self.workflow["on"]["push"]["tags"], ["v*"])
        for name in ("publish", "github-release"):
            condition = self.workflow["jobs"][name]["if"]
            self.assertIn("github.event_name == 'push'", condition)
            self.assertIn("github.repository == 'fastnsa/FastNSA'", condition)
            self.assertIn("startsWith(github.ref, 'refs/tags/v')", condition)

    def test_only_publisher_has_oidc_permission(self):
        self.assertEqual(self.workflow["permissions"], {"contents": "read"})
        jobs = self.workflow["jobs"]
        self.assertEqual(jobs["publish"]["environment"]["name"], "pypi")
        for name, job in jobs.items():
            self.assertEqual(job.get("permissions", {}).get("id-token"),
                             "write" if name == "publish" else None)
            if name != "github-release":
                self.assertNotEqual(job.get("permissions", {}).get("contents"), "write")
        action = next(s["uses"] for s in jobs["publish"]["steps"]
                      if s.get("uses", "").startswith("pypa/"))
        self.assertRegex(action, r"^pypa/gh-action-pypi-publish@[0-9a-f]{40}$")
        self.assertNotIn("secrets", json.dumps(jobs["publish"]))

    def test_all_platform_checks_gate_assembly_and_publish(self):
        jobs = self.workflow["jobs"]
        self.assertIn("workflow_call", self.matrix["on"])
        self.assertEqual(jobs["checks"]["uses"], "./.github/workflows/multiplatform.yml")
        self.assertTrue({"native", "wheel", "installed-wheel", "linux-geo", "release-checks"}
                        <= set(self.matrix["jobs"]))
        self.assertEqual(set(jobs["assemble"]["needs"]), {"preflight", "checks", "sdist"})
        self.assertIn("assemble", jobs["publish"]["needs"])
        self.assertIn("publish", jobs["github-release"]["needs"])

    def test_github_release_stays_draft_and_source_smoke_wheel_not_uploaded(self):
        jobs = self.workflow["jobs"]
        command = next(s["run"] for s in jobs["github-release"]["steps"]
                       if s.get("name") == "Create a draft GitHub Release")
        self.assertIn("--draft", command)
        self.assertIn("--verify-tag", command)
        upload = next(s for s in jobs["sdist"]["steps"]
                      if s.get("uses", "").startswith("actions/upload-artifact"))
        self.assertEqual(upload["with"]["path"], "source-dist/*.tar.gz")
        self.assertEqual(self.workflow["concurrency"]["cancel-in-progress"], "false")


if __name__ == "__main__":
    unittest.main()
