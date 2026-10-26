"""Check that package, CMake, and source fallback versions stay consistent."""

import importlib.util
from importlib.metadata import PackageNotFoundError
from pathlib import Path
import re
import tomllib
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]


class ReleaseVersionTest(unittest.TestCase):
    def test_cmake_version_matches_package_metadata(self):
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        cmake = (ROOT / "cpp/CMakeLists.txt").read_text(encoding="utf-8")
        match = re.search(r"project\(fastnsa\s+VERSION\s+([^\s)]+)", cmake)
        self.assertIsNotNone(match)
        self.assertEqual(match.group(1), project["project"]["version"])

    def test_source_fallback_matches_package_version(self):
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        spec = importlib.util.spec_from_file_location(
            "_fastnsa_version_check", ROOT / "python/fastnsa/__init__.py"
        )
        module = importlib.util.module_from_spec(spec)
        with patch("importlib.metadata.version", side_effect=PackageNotFoundError("fastnsa")):
            spec.loader.exec_module(module)
        self.assertEqual(module.__version__, project["project"]["version"] + ".dev0")


if __name__ == "__main__":
    unittest.main()
