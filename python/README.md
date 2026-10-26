# FastNSA Python

FastNSA is a native-accelerated Python toolkit for network-constrained spatial
point analysis. It provides point-to-network mapping, network kernel density
estimation (NKDE), network K-functions, and spatiotemporal network K-functions.

The package has two integration levels:

- A lightweight array API based on NumPy and the native FastNSA library.
- A geospatial API for GeoPandas, Shapely, NetworkX, and OSMnx workflows.

## Installation

Python 3.11 or newer is required. Install the base NumPy/native runtime with:

~~~bash
python -m pip install fastnsa
~~~

Published binary wheels determine the supported operating systems, CPU
architectures, and Python versions for each release. If no compatible wheel is
available, pip may attempt a source build, which requires CMake and a C++17
compiler.

Install only the optional integrations you need:

~~~bash
python -m pip install "fastnsa[dataframe]"  # pandas exports
python -m pip install "fastnsa[xarray]"     # ST result DataArray exports
python -m pip install "fastnsa[geo]"        # GeoPandas and Shapely workflows
python -m pip install "fastnsa[osm]"        # NetworkX and OSMnx workflows
python -m pip install "fastnsa[viz]"        # Matplotlib visualizations
python -m pip install "fastnsa[all]"        # all optional integrations
~~~

## Minimal Example

The bundled line dataset is useful for installation checks and examples. It
does not require GeoPandas, NetworkX, or OSMnx.

~~~python
import fastnsa as nsa

network = nsa.make_line_network(length=1000.0)
events = nsa.make_line_events(network, offsets=[100.0, 350.0, 700.0, 900.0])

result = nsa.network_k_function(
    network=network,
    points=events,
    r_values=[50.0, 100.0, 250.0],
    method=nsa.KFunctionMethod.STARS,
    simulations=5,
    random_state=2026,
)

print(result.to_records())
result.save_npz("network_k_result.npz")
restored = nsa.load_result("network_k_result.npz")

events.close()
network.close()
~~~

Native resources are context managers, so larger applications can use a
deterministic lifecycle:

~~~python
with nsa.make_line_network() as network:
    with nsa.make_line_events(network) as events:
        result = nsa.network_k_function(
            network=network,
            points=events,
            r_max=500.0,
            r_step=50.0,
        )
~~~

## Geospatial Workflow

With the osm extra installed, a network can be downloaded through OSMnx and
event files can be loaded through GeoPandas:

~~~python
import fastnsa as nsa

network = nsa.Network.from_osm(place="Central Park, New York", network_type="drive")
points = nsa.load_points("events.gpkg", layer="events", time_col="event_time")
events = nsa.PointEvents(points, target_network=network, time_col="event_time")

k_result = nsa.network_k_function(
    network=network,
    points=events,
    r_min=0.0,
    r_max=2000.0,
    r_step=50.0,
    method=nsa.KFunctionMethod.STARS,
    simulations=99,
    envelope="global",
)
~~~

Use a projected coordinate reference system whose map units are metres for
distance-based analyses. FastNSA does not interpret degree coordinates as
geodesic distances.

## Main Features

- Point mapping with R-tree, exact comparison, and approximate comparison.
- NKDE methods including ADA and LION with optional event weights.
- Network K-function methods SPS, NS, and STARS.
- Spatiotemporal methods RQS, SPS, EAR, MTS, and EARTH.
- Pointwise or global Monte Carlo envelopes, reproducible seeds, and p-values.
- CSV, pandas, NPZ, and optional xarray result exports.
- Progress callbacks and cooperative cancellation for host applications.
- Direct NumPy/native inputs for low-dependency QGIS and server integrations.

## Documentation

- [Documentation index](docs/index.md)
- [Installation and deployment](docs/installation.md)
- [Quick start](docs/quickstart.md)
- [Data model and preprocessing](docs/data-model.md)
- [Analysis methods](docs/analyses.md)
- [Results and persistence](docs/results.md)
- [API reference](docs/api-reference.md)
- [Algorithm references](docs/references.md)
- [Troubleshooting](docs/troubleshooting.md)
- [中文使用指南](docs/README_ZH.md)

## License

FastNSA is licensed under GNU LGPL v2.1 only. Applications may dynamically
link to the library under the LGPL terms. Modifications to the FastNSA library
itself must remain available under the applicable LGPL requirements. See the
packaged LICENSE file for the complete license text.
