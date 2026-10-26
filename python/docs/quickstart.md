# Quick Start

This guide starts with deterministic in-memory data, then shows the equivalent
GeoPandas workflow.

## 1. Verify the Runtime

~~~python
import fastnsa as nsa

print(nsa.__version__)
~~~

For a complete installation check, run python -m fastnsa._selftest from a
terminal.

## 2. Create a Small Network and Events

~~~python
import fastnsa as nsa

network = nsa.make_line_network(length=1000.0, crs="EPSG:3857")
events = nsa.make_line_events(
    network,
    offsets=[100.0, 250.0, 550.0, 850.0],
    times=[0.0, 2.0, 5.0, 9.0],
)

print(network.num_edges)
print(network.total_length)
print(events.edge_ids_array())
print(events.offsets_array())
~~~

An event location is represented by an edge ID and an offset measured from the
edge start. This mapped representation is reused by every analysis.

## 3. Compute Network KDE

~~~python
kde = nsa.network_kde(
    network=network,
    points=events,
    bandwidth=250.0,
    lixel_length=50.0,
    method=nsa.NKDVMethod.LION,
    kernel=nsa.NKDVKernel.EPANECHNIKOV,
    normalize=True,
)

print(kde.summary())
kde.to_csv("line_kde.csv")
~~~

The output contains one row per lixel segment, including edge ID, start/end
offsets, and density.

## 4. Compute a Network K-function

~~~python
k_result = nsa.network_k_function(
    network=network,
    points=events,
    r_values=[50.0, 100.0, 200.0, 400.0],
    method=nsa.KFunctionMethod.STARS,
    simulations=19,
    alpha=0.05,
    envelope="pointwise",
    random_state=2026,
)

print(k_result.summary())
print(k_result.envelope_description)
k_result.to_csv("line_k.csv")
k_result.save_npz("line_k.npz")
~~~

Use a larger simulation count for final inference. Small counts such as 5 or
19 are useful for functional tests but provide coarse Monte Carlo estimates.

## 5. Compute a Spatiotemporal K-function

~~~python
st_result = nsa.spatiotemporal_network_k_function(
    network=network,
    points=events,
    s_values=[100.0, 300.0, 500.0],
    t_values=[1.0, 3.0, 5.0],
    method=nsa.SpatiotemporalKFunctionMethod.EARTH,
    simulations=5,
    random_state=2026,
)

print(st_result.k_obs.shape)
st_result.to_csv("line_st_k.csv")
~~~

Time values are numeric and their unit is determined by the application. For
example, if input times are days since an origin, temporal thresholds are also
days. Store the chosen unit and origin in surrounding workflow metadata.

## 6. Restore Results

~~~python
restored = nsa.load_result("line_k.npz")
assert isinstance(restored, nsa.KFunctionResult)
print(restored.method, restored.random_seed)
~~~

NPZ preserves arrays and analysis metadata. CSV is intended for interchange
and does not reconstruct a typed result.

For scientific reporting, cite the publication associated with the selected
method from the [algorithm references](references.md), and record the FastNSA
version and parameters alongside the result.

## 7. Clean Up Native Resources

~~~python
events.close()
network.close()
~~~

Prefer context managers in applications that repeatedly create analyses:

~~~python
with nsa.make_line_network() as network:
    with nsa.make_line_events(network) as events:
        result = nsa.network_k_function(network=network, points=events)
~~~

## Working with Real Files

Install the geo extra, load a point file, and construct or download a network:

~~~python
import fastnsa as nsa

points = nsa.load_points(
    "events.gpkg",
    layer="events",
    time_col="occurred_at",
)
network = nsa.Network.from_osm(
    place="Gainesville, Florida, USA",
    network_type="drive",
)
events = nsa.PointEvents(
    points,
    target_network=network,
    time_col="occurred_at",
    mapping_method=nsa.MappingMethod.RTREE,
)
~~~

The point layer is reprojected to the network CRS by PointEvents. The network
must itself use a projected CRS suitable for the chosen distance values.

