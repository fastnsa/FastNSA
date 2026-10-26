# Public API Reference

This is a compact reference for the supported top-level API. Import public
objects from fastnsa rather than private modules.

~~~python
import fastnsa as nsa
~~~

## Data Construction

### NativeNetwork

~~~text
NativeNetwork(
    *, edges, offsets, indices, coordinates=None, geometry_offsets=None,
    wkb_blob=b"", wkb_offsets=None, crs=None
)
~~~

Creates an owning native network from contiguous-compatible arrays. Packed
coordinates are preferred; LineString WKB remains available for compatibility.
Properties: num_edges, num_nodes, and total_length. Supports close and context
management.

### NativePointSet

~~~text
NativePointSet.map_xy(
    network, x, y, *, mapping_method="R-tree", comp_eps=10.0,
    n_threads=0, times=None
)
NativePointSet.from_mapped(network, edge_ids, offsets, *, times=None)
~~~

Methods: set_times, edge_ids_array, offsets_array, times_array, and close.

### Network

~~~text
Network(graph)
Network.from_prepared_arrays(
    *, crs, edges, offsets, indices, coordinates=None, geometry_offsets=None,
    wkb_blob=b"", wkb_offsets=None, graph=None
)
Network.from_osm(
    bbox=None, place=None, network_type="drive", crs=None,
    consolidate_tolerance=20
)
~~~

Graph-compatible wrapper around NativeNetwork. save_graph_gpkg requires a
graph-backed object and the osm extra.

### PointEvents

~~~text
PointEvents(
    points_gdf, target_network, mapping_method="R-tree", comp_eps=10.0,
    n_threads=0, time_col=None, times=None, time_unit=None,
    time_origin=None, time_mode=None
)
PointEvents.from_mapped_arrays(
    points_gdf, target_network, edge_ids, offsets, *, time_col=None,
    times=None, time_unit=None, time_origin=None, time_mode=None
)
~~~

Methods: edge_ids_array, offsets_array, times_array, to_mapping_records,
dump_to_file, dump_to_file_st, and close. timings records conversion phases.

### load_points

~~~text
load_points(
    points_path, return_bounding_box=False, time_col=None, *,
    x_col="lon", y_col="lat", crs="EPSG:4326", layer=None
)
~~~

Loads CSV, GeoPackage, GeoJSON/JSON, or GeoParquet through the geo extra.

### preprocess_events

~~~text
preprocess_events(
    points, network, *, edge_id_col=None, offset_col=None, time_col=None,
    times=None, time_unit=None, time_origin=None, time_mode=None,
    mapping_method="R-tree", comp_eps=10.0, n_threads=0
)
~~~

Returns PreprocessResult with events, reused_mapping, n_points, timings, and
mapping_method. Supports context management.

### Example Datasets

~~~text
make_line_network(length=1000.0, *, crs="EPSG:3857")
make_line_events(network, offsets=None, *, times=None)
~~~

Return deterministic NativeNetwork and NativePointSet objects for examples and
smoke tests.

## Analysis Functions

### network_kde

~~~text
network_kde(
    *, network, points, bandwidth, lixel_length, method="LION",
    kernel="EPANECHNIKOV", normalize=False, progress_callback=None,
    weights=None
)
~~~

Returns NetworkKDEPointResult.

### network_k_function

~~~text
network_k_function(
    *, network, points, r_min=0.0, r_max=2000.0, r_step=50,
    r_values=None, method="NS", simulations=0, alpha=0.05,
    return_simulations=False, epsilon=0.1, progress_callback=None,
    random_state=5489, envelope="pointwise"
)
~~~

Returns KFunctionResult.

### spatiotemporal_network_k_function

~~~text
spatiotemporal_network_k_function(
    *, network, points, s_min=0.0, s_max=2000.0, s_step=50,
    t_min=0.0, t_max=2000.0, t_step=50, s_values=None, t_values=None,
    method="EARTH", simulations=0, alpha=0.05,
    return_simulations=False, progress_callback=None, random_state=5489,
    envelope="pointwise", return_xarray=False
)
~~~

Returns SpatiotemporalKFunctionResult, or the observed xarray DataArray when
return_xarray is true.

## Result Classes

- KFunctionResult
- SpatiotemporalKFunctionResult
- NetworkKDEPointResult

All provide to_records, to_dataframe, save_csv/to_csv, save_npz, load_npz, and
summary. K-function results additionally provide envelope_status and
envelope_description. ST results provide to_xarray.

### load_result

~~~text
load_result(filename)
~~~

Loads a versioned NPZ result and returns the appropriate result class.

## Stable Enums

| Enum | Values |
| --- | --- |
| MappingMethod | RTREE, COMP_APPROX, COMP_EXACT |
| KFunctionMethod | SPS, NS, STARS |
| SpatiotemporalKFunctionMethod | RQS, SPS, EAR, MTS, EARTH |
| NKDVMethod | BASIC, SHARED_SPD, ADA, IA, ADA_IA, LION |
| NKDVKernel | EPANECHNIKOV |

Functions also accept the corresponding exported string values.

## Exceptions

| Exception | Meaning |
| --- | --- |
| FastNSAError | Base FastNSA exception |
| FastNSAValidationError | Invalid input or parameter |
| NativeBackendError | Native backend reported an error |
| FastNSACancelledError | Cooperative cancellation completed |

Standard ValueError, RuntimeError, ImportError, FileNotFoundError, and OSError
may still be raised for conventional Python failures described in messages.

## Progress Callback

~~~text
callback(stage: str, completed: int, total: int) -> bool
~~~

Return true to request cancellation. Exceptions raised by a callback cannot
cross the C boundary and are converted into a cancellation request.

## Logging

FastNSA uses the standard logging hierarchy under fastnsa. Library code does
not install handlers. Structured extra fields may include fastnsa_stage,
fastnsa_method, fastnsa_points, fastnsa_simulations, and fastnsa_duration_s.

~~~python
import logging

logging.basicConfig(level=logging.INFO)
logging.getLogger("fastnsa").setLevel(logging.INFO)
~~~

## Algorithm Citations

See [Algorithm references](references.md) for the original matching, NKDE,
Network K-function, STARS, and EARTH publications and research artifacts.

