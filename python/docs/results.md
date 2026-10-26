# Results and Persistence

FastNSA returns typed dataclasses instead of unstructured tuples. The classes
retain analysis metadata, provide tabular exports, and support versioned NPZ
persistence.

## Common Operations

~~~python
print(result.summary())
rows = result.to_records()
frame = result.to_dataframe()
result.to_csv("result.csv")
path = result.save_npz("result.npz")
restored = nsa.load_result(path)
~~~

to_dataframe requires the dataframe optional extra. CSV export uses the Python
standard library and does not require pandas. save_csv and to_csv are aliases.

## KFunctionResult

| Field | Meaning |
| --- | --- |
| r | Network distance thresholds |
| k_obs | Observed K(r) values |
| k_theo | Mean Monte Carlo reference, if simulations were run |
| k_lo, k_hi | Envelope bounds |
| k_sim | Simulation matrix when return_simulations is true |
| method | Backend algorithm name |
| n_points | Event count |
| intensity | Event count divided by network length |
| network_length | Total network length |
| random_seed | Effective Monte Carlo seed |
| envelope_method | pointwise or global |
| confidence_level | One minus alpha |
| p_value | Global maximum-deviation Monte Carlo p-value |

envelope_status is an int8 array with -1 below the lower envelope, 0 inside,
and 1 above the upper envelope. It is None when no envelope exists.
envelope_description is a human-readable property suitable for report notes.

## SpatiotemporalKFunctionResult

The observed, theoretical, and envelope surfaces have shape (n_s, n_t).
Simulation output has shape (n_simulations, n_s, n_t).

In addition to K-function metadata, the object stores time_unit, time_origin,
and time_mode copied from the point set. These values are descriptive metadata;
the backend always operates on numeric time values.

Export a surface to xarray:

~~~python
observed = result.to_xarray("k_obs")
upper = result.to_xarray("k_hi")
~~~

The resulting DataArray dimensions are spatial_threshold and
temporal_threshold. The xarray extra is required.

## NetworkKDEPointResult

| Field | Meaning |
| --- | --- |
| edge_ids | Native edge ID for each lixel |
| offsets_n1 | Distance from edge start to evaluation position |
| offsets_n2 | Distance from evaluation position to edge end |
| density | Estimated density |
| method, kernel | Native analysis choices |
| bandwidth, lixel_length | Distance parameters |
| normalized | Whether network-integral normalization was applied |
| weighted, total_weight | Event-weight metadata |

Derived segment_start_n1, segment_end_n1, and segment_lengths properties turn
evaluation offsets into lixel intervals. Tabular exports include these derived
columns and can be joined to host network geometry by edge_id.

## CSV versus NPZ

Use CSV when people or other software need a flat table. Use NPZ when the
result will be reopened by FastNSA.

CSV does not preserve all metadata, simulation arrays, exact dtypes, or result
class identity. NPZ preserves these values and uses allow_pickle=False while
loading. load_result reads its result_type marker and dispatches to the correct
class.

Class-specific loading is also available:

~~~python
restored = nsa.KFunctionResult.load_npz("network_k.npz")
~~~

Do not manually modify NPZ members or rely on undocumented member names. The
format is versioned so future releases can validate compatibility.

## Reproducible Reports

Archive at least:

- the NPZ result;
- network identity, CRS, and edge ordering or a stable network signature;
- point source and mapping parameters;
- method and all distance/time parameters;
- simulation count, random seed, alpha, and envelope type;
- FastNSA version and native platform;
- time unit and origin for spatiotemporal work.

summary provides a compact machine-readable starting point but is not a full
provenance record.

