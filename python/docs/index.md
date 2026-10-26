# FastNSA Python Documentation

FastNSA combines a NumPy-facing Python API with a C++ backend for
network-constrained spatial point analysis. This documentation covers the
standalone Python package. QGIS-specific workflows are documented separately
in the plugin documentation.

## Choose a Starting Point

- New users: follow the [quick start](quickstart.md).
- Package maintainers: read [installation and deployment](installation.md).
- GIS users: read [data model and preprocessing](data-model.md).
- Analysts: compare methods in [analysis methods](analyses.md).
- Application developers: use the [API reference](api-reference.md).
- Chinese readers: start with the [中文使用指南](README_ZH.md).

## Documentation Map

| Guide | Contents |
| --- | --- |
| [Installation](installation.md) | Wheels, optional dependencies, source builds, native library loading |
| [Quick start](quickstart.md) | A complete network, event, analysis, and export workflow |
| [Data model](data-model.md) | Networks, mapped events, CRS, time, preprocessing, direct arrays |
| [Analyses](analyses.md) | NKDE, network K, ST K, Monte Carlo inference, performance |
| [Results](results.md) | Result fields, DataFrame, CSV, NPZ, xarray, reproducibility |
| [API reference](api-reference.md) | Public functions, classes, enums, exceptions, callbacks |
| [Algorithm references](references.md) | Original papers, DOIs, and research artifacts |
| [Troubleshooting](troubleshooting.md) | Loading failures, validation errors, crashes, performance |

## Supported Analysis Pipeline

The normal workflow has four stages:

1. Construct a projected road network.
2. Map event points to edge IDs and along-edge offsets.
3. Run NKDE, network K-function, or spatiotemporal K-function analysis.
4. inspect, export, or persist the typed result object.

Mapping can be reused. If a dataset already contains valid edge IDs and
offsets for the same network, use preprocess_events or
PointEvents.from_mapped_arrays instead of repeating nearest-edge matching.

## Units and Statistical Interpretation

FastNSA treats all network distances as numeric map units. For usual road
analysis, use a projected CRS in metres. A geographic CRS such as EPSG:4326
stores coordinates in degrees and is not suitable for metre-valued bandwidths
or thresholds.

K-function theoretical curves and envelopes are simulation-derived. When
simulations is zero, no Monte Carlo theoretical curve or envelope is returned.
The result metadata records the method, random seed, envelope type, confidence
level, network length, and number of points needed to interpret an analysis.

## Version Scope

These documents describe the public API of FastNSA 0.2.x. Native result files
include a format marker and should be loaded with fastnsa.load_result rather
than by reading NPZ members directly.

