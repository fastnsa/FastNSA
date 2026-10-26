# Data Model and Preprocessing

## Core Objects

FastNSA analyses consume a network and a mapped point set. Each has a
lightweight class and a compatibility class.

| Purpose | Lightweight | Geospatial compatibility |
| --- | --- | --- |
| Network | NativeNetwork | Network |
| Events | NativePointSet | PointEvents |

NativeNetwork and NativePointSet require only NumPy and the bundled native
library. Network and PointEvents add graph, geometry, CRS, and file-oriented
conveniences. Both expose the native pointers expected by analysis functions.

## Native Network Arrays

NativeNetwork accepts:

- edges: float array with shape (n_edges, 3), containing start node index, end
  node index, and edge length.
- offsets and indices: int32 CSR adjacency arrays.
- coordinates: float array with shape (n_vertices, 2), containing packed edge
  vertices. This is the preferred dependency-free input path.
- geometry_offsets: int32 positions with n_edges + 1 entries delimiting each
  edge in coordinates.
- wkb_blob and wkb_offsets: optional compatibility input for concatenated
  LineString WKB geometries.
- crs: optional application metadata.

The final CSR offset must equal the number of adjacency indices. The final
geometry offset must equal the number of coordinate rows. Node and edge
identifiers passed to the native layer must fit signed 32-bit integers.

~~~python
network = nsa.NativeNetwork(
    edges=edges,
    offsets=csr_offsets,
    indices=csr_indices,
    coordinates=coordinates,
    geometry_offsets=geometry_offsets,
    crs="EPSG:32618",
)
~~~

Edge row order defines the native edge ID. Preserve that order when storing or
reusing mapped event columns. FastNSA builds and retains a native BVH over the
packed polylines, so repeated event mapping does not rebuild the spatial index.

## Graph-backed Networks

Network accepts a NetworkX-style multigraph. It copies the graph, aligns edge
geometry direction with graph direction, assigns contiguous native node and
edge indices, and packs edge geometry.

Required graph properties are:

- graph.graph["crs"] for CRS-aware workflows;
- a finite edge length attribute;
- LineString geometry, or finite x/y attributes on both endpoint nodes.

Network.from_osm can construct and project a graph from a place or bounding
box. It requires the osm optional extra. Network.from_prepared_arrays wraps
arrays already prepared by QGIS or another host without constructing a full
NetworkX graph.

## Point Mapping

NativePointSet.map_xy maps coordinate arrays to the nearest network edge:

~~~python
events = nsa.NativePointSet.map_xy(
    network,
    x,
    y,
    mapping_method=nsa.MappingMethod.RTREE,
    n_threads=4,
    times=times,
)
~~~

Set `n_threads=0` to use the bounded automatic policy. Automatic mapping uses
at most eight workers and scales down for small point sets. Explicit requests
are capped by detected CPU availability and a default hard limit of 32; set
`FASTNSA_MAX_THREADS` before loading FastNSA to override that hard limit.

Available mapping methods are:

| Enum | String | Guidance |
| --- | --- | --- |
| MappingMethod.RTREE | R-tree | General default for spatial lookup |
| MappingMethod.COMP_EXACT | comp_e | Exhaustive comparison reference |
| MappingMethod.COMP_APPROX | comp_a | Comparison with comp_eps tolerance |

Use exact comparison mainly for validation or small datasets. Benchmark on the
target network before changing the R-tree default.

The matching algorithms are based on Ye et al., *An Efficient and Accurate
Grid Compression Solution for Point-to-Line Nearest Neighbor Search* (KDD
2026). See the [matching reference and artifact](references.md#point-to-line-nearest-edge-matching).

PointEvents performs the same operation from a point GeoDataFrame. If point and
network CRS values differ, it reprojects the points to the network CRS before
extracting x/y coordinates.

## Reusing Mapped Events

Mapping is often a noticeable part of repeated analyses. Reuse valid edge IDs
and offsets when the network has not changed:

~~~python
with nsa.preprocess_events(
    points,
    network,
    edge_id_col="edge_id",
    offset_col="offset",
    time_col="event_time",
) as prepared:
    result = nsa.network_k_function(
        network=network,
        points=prepared.events,
        r_max=5000.0,
        r_step=100.0,
    )
    print(prepared.reused_mapping)
    print(prepared.timings)
~~~

edge_id_col and offset_col must be supplied together. These values are only
valid for the exact network topology and edge order used during mapping. An
application should store a network signature beside them and remap if edge
geometry, order, or topology changes.

The lower-level equivalent is:

~~~python
events = nsa.NativePointSet.from_mapped(
    network,
    edge_ids=edge_ids,
    offsets=offsets,
    times=times,
)
~~~

## CRS and Distance Units

FastNSA computes planar distances from network edge lengths and geometry. It
does not perform geodesic conversion for longitude/latitude coordinates.

- Project the network before constructing the FastNSA object.
- Use a local projected CRS appropriate for the study region.
- Express bandwidth, lixel length, and spatial thresholds in the same unit as
  edge lengths.
- Recompute stored mappings after changing the network CRS or geometry.

## Time Values

NativePointSet accepts finite float64 times. The backend does not attach a unit
or timezone to those numbers. Spatial-temporal thresholds must use the same
unit.

PointEvents supports either:

- times: an explicit numeric array; or
- time_col: a numeric or datetime-compatible GeoDataFrame column.

Datetime columns are parsed as UTC and converted to Unix seconds. Numeric
columns are passed through without unit conversion. Applications should set
time_unit, time_origin, and time_mode metadata when those concepts are known.

Spatiotemporal analysis requires at least two events and a positive observed
time span. MTS and EARTH additionally require evenly spaced threshold axes.

## Loading Point Files

load_points supports:

- CSV with configurable x_col and y_col, defaulting to lon and lat;
- GeoPackage, optionally selecting a layer;
- GeoJSON and JSON;
- Parquet and GeoParquet.

~~~python
points = nsa.load_points(
    "events.parquet",
    time_col="timestamp",
)
~~~

The result must be a non-empty point GeoDataFrame with a CRS. CSV coordinates
use the crs argument, which defaults to EPSG:4326. The legacy
return_bounding_box option returns max_y, min_y, max_x, min_x; new code should
prefer points.total_bounds to avoid relying on that historical ordering.

## Resource Ownership

Network objects own native memory. Point sets refer to their network, so close
point sets before closing the network. Context managers are the safest pattern.
Calling an analysis after either object is closed raises RuntimeError.

