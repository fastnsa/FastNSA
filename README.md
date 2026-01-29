# FastNSA

FastNSA is a high-performance toolkit for **network-constrained spatial analysis**,
designed for large-scale road networks and massive point datasets.
It provides a Python interface backed by a native implementation for efficient
spatial indexing and network-based statistical computation.



---

## Installation

FastNSA is distributed via PyPI with prebuilt native backends.

```bash
pip install fastnsa
```
Supported platforms:

* Linux (x86_64)

* Python ≥ 3.9

## Quick Example

```python
import fastnsa as nsa

# 1. Build a road network from OpenStreetMap
# (internally constructs spatial and network indexes)
net = nsa.Network.from_osm(place="Central Park, NY")

# 2. Load event points and project them onto the network
gdf = nsa.load_points("accidents.csv")
events = nsa.PointEvents(gdf, target_network=net)

# 3. Compute a network-based K-function
results = nsa.network_k_function(
    network=net,
    points=events,
    r_values=[10, 50, 100]
)

print(results)
```

##  Design Overview

The Python layer provides:

* Network construction and I/O

* Event projection to networks

* tatistical analysis interfaces

Performance-critical components such as spatial indexing,
point-to-network mapping, and distance evaluation are implemented
in a native backend and exposed through Python bindings.


## License

FastNSA is licensed under the GNU Lesser General Public License v3.0 (LGPL-2.1).

You are free to use this library in proprietary and open-source software,
but any modifications to the FastNSA library itself must be released under LGPL.
