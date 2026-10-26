# FastNSA

FastNSA is a Python toolkit for network-constrained spatial point analysis.
It combines a C++ computational backend with Python APIs for studying events
along road networks and other linear networks.

## Features

- Point-to-network matching and reusable event preprocessing.
- Network kernel density estimation with ADA and LION, including event weights.
- Network K-functions with exact and approximate methods.
- Spatiotemporal network K-functions for events with time information.
- Monte Carlo envelopes and reproducible simulations.
- CSV and NPZ result storage, with optional pandas and xarray integration.

The lightweight array API requires NumPy. Optional integrations support
GeoPandas, Shapely, NetworkX, OSMnx, and Matplotlib.

## Installation

Python 3.11 or newer is required. To install the current source version:

```bash
git clone https://github.com/fastnsa/FastNSA.git
cd FastNSA
python -m pip install .
```

Building from source also requires CMake 3.18 or newer and a C++17 compiler.
See the [installation guide](python/docs/installation.md) for wheel installation
and optional dependencies.

## Quick Start

This example uses a small built-in network and requires no GIS data files:

```python
import fastnsa as nsa

with nsa.make_line_network(length=1000.0) as network:
    with nsa.make_line_events(network, offsets=[100.0, 350.0, 700.0, 900.0]) as events:
        result = nsa.network_k_function(
            network=network,
            points=events,
            r_values=[50.0, 100.0, 250.0],
            simulations=0,
        )
        print(result.to_records())
```

For real-world data, use a projected coordinate reference system with map
units appropriate for distance analysis.

## Documentation

- [Package overview](python/README.md)
- [User guide](python/docs/index.md)
- [Quick start](python/docs/quickstart.md)
- [Analysis methods](python/docs/analyses.md)
- [API reference](python/docs/api-reference.md)
- [Algorithm references](python/docs/references.md)
- [Chinese guide](python/docs/README_ZH.md)

## Support and License

Report bugs and feature requests through the
[issue tracker](https://github.com/fastnsa/FastNSA/issues).
FastNSA is licensed under [LGPL-2.1-only](LICENSE).
