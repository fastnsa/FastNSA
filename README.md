# FastNSA

Native-accelerated network spatial analysis for Python. This repository owns
the Python package, C++ core, C ABI, and platform wheels. QGIS plugin sources
are maintained in a separate repository.

- [Python package and quick start](python/README.md)
- [User documentation](python/docs/index.md)
- [Algorithm references](python/docs/references.md)
- [License](LICENSE)

## Build and Test

Python 3.11 or newer, CMake 3.18 or newer, and a C++17 compiler are required
for a source build. Compatible wheels include the native library.

```bash
python -m pip install .
python -m fastnsa._selftest
python -m unittest discover -s tests/python -t . -v
cmake -S cpp -B build/native -DCMAKE_BUILD_TYPE=Release -DFASTNSA_BUILD_TESTS=ON
cmake --build build/native --config Release --parallel 2
ctest --test-dir build/native -C Release --output-on-failure
```

Network-dependent integration tests are opt-in; see their module docstrings.
The native/wheel matrix covers Linux x86_64, Windows AMD64, and macOS Intel
and Apple Silicon. CI success does not certify every external host environment.

`EXPORT_PROVENANCE.json` records the upstream snapshot and exported file hashes.
It is a provenance record, not a public release approval or a secret audit.
