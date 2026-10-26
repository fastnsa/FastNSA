"""Audit the native runtime inside a Linux release wheel using readelf."""

import argparse
from pathlib import Path
import re
import subprocess
import tempfile
import zipfile


def check_wheel(path):
    with zipfile.ZipFile(path) as wheel:
        if not any(
            tag in path.name
            for tag in ("manylinux_2_17_x86_64", "manylinux2014_x86_64")
        ):
            raise ValueError(f"Expected a GLIBC 2.17 manylinux x86_64 wheel: {path.name}")
        libraries = [
            name for name in wheel.namelist()
            if name.endswith("/bindings/libfastnsa.so")
        ]
        if len(libraries) != 1:
            raise ValueError(f"Expected exactly one bundled libfastnsa.so: {libraries}")
        with tempfile.TemporaryDirectory(prefix="fastnsa-wheel-audit-") as temporary:
            library = Path(temporary) / "libfastnsa.so"
            library.write_bytes(wheel.read(libraries[0]))
            report = subprocess.run(
                ["readelf", "-d", str(library)],
                check=True, capture_output=True, text=True,
            ).stdout
        dependencies = re.findall(r"\(NEEDED\).*?\[(.*?)\]", report)
        if not dependencies:
            raise ValueError("Could not identify native dynamic dependencies.")
        print(f"{path.name}: {libraries[0]}")
        for dependency in dependencies:
            print(f"  NEEDED {dependency}")
            if re.match(r"lib(?:geos|gomp|omp|iomp)", dependency, re.IGNORECASE):
                raise ValueError(f"Unexpected GEOS/OpenMP dependency: {dependency}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheels", type=Path, nargs="+")
    arguments = parser.parse_args()
    for path in arguments.wheels:
        check_wheel(path)


if __name__ == "__main__":
    main()
