# Installation and Deployment

## Requirements

- Python 3.11 or newer; release CI targets Python 3.11 and 3.12.
- A 64-bit operating system and interpreter for published wheels.
- NumPy 1.23.5 or newer.
- A compatible FastNSA native library bundled in the wheel or supplied by the
  developer.

## Install from a Wheel

~~~bash
python -m pip install --upgrade pip
python -m pip install fastnsa
python -m fastnsa._selftest
~~~

The self-test loads the native library and runs a small analysis. It is a more
useful deployment check than importing the package alone.

Binary wheel availability varies by release. Check the files attached to the
release or package index before assuming that a particular Python and platform
combination is available. If a compatible wheel does not exist, pip may fall
back to a source distribution.

## Optional Dependencies

The base package intentionally depends only on NumPy. Optional extras keep
QGIS and embedded deployments small.

| Extra | Packages | Use case |
| --- | --- | --- |
| dataframe | pandas | Result DataFrame export |
| xarray | xarray | Spatiotemporal surface export |
| geo | pandas, GeoPandas, Shapely | Point files and geometry-based workflows |
| osm | geo packages, NetworkX, OSMnx | OSM download and graph construction |
| viz | Matplotlib | User-provided plotting workflows |
| all | All of the above | Research and notebook environment |

Install an extra with:

~~~bash
python -m pip install "fastnsa[geo,dataframe]"
~~~

## Build from Source

A source build requires CMake, a C++17 compiler, Python development files, and
the build tools declared in pyproject.toml.

~~~bash
git clone https://github.com/fastnsa/FastNSA.git
cd FastNSA
python -m pip install --upgrade pip
python -m pip install .
python -m fastnsa._selftest
~~~

For editable development:

~~~bash
python -m pip install -e ".[all]"
python -m unittest discover -s tests/python -t . -v
~~~

Release builds enable the CMake Release configuration. The project CMake
configuration also applies supported link-time optimization, section garbage
collection, and hidden symbol visibility to reduce the native binary surface.

## Native Library Discovery

The loader expects the platform-specific library beside the Python bindings:

| Platform | Filename |
| --- | --- |
| Linux | libfastnsa.so |
| macOS | libfastnsa.dylib |
| Windows | fastnsa.dll |

Development builds can override that location:

~~~bash
export FASTNSA_LIBRARY_PATH=/absolute/path/to/libfastnsa.so
python -m fastnsa._selftest
~~~

On Windows, set the same environment variable to the DLL path. The loader adds
its parent directory to the DLL search path for Python 3.8 and newer.

Do not point FASTNSA_LIBRARY_PATH at a library built for a different CPU,
operating system, C runtime, or pointer width.

## QGIS and Embedded Python

QGIS uses its own Python interpreter and package search path. Installing
FastNSA into a system Python does not make it available inside QGIS. Use the
Python executable belonging to the QGIS environment, and verify from the QGIS
Python console:

~~~python
import fastnsa
print(fastnsa.__version__)
print(fastnsa.__file__)
~~~

For a plugin, prefer the base runtime and NativeNetwork/NativePointSet path
when the host already supplies geometries and attributes. This avoids pulling
GeoPandas, OSMnx, and NetworkX into QGIS solely for data conversion.

## Release Verification

Test each produced wheel in a clean environment:

~~~bash
python -m venv wheel-check
wheel-check/bin/python -m pip install dist/fastnsa-*.whl
wheel-check/bin/python -m fastnsa._selftest
~~~

Use the appropriate Scripts path on Windows. Also test one optional geospatial
workflow separately because the base self-test deliberately avoids optional
packages.

