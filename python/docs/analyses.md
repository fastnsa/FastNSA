# Analysis Methods

## Network Kernel Density Estimation

The ADA/IA augmentation family and LION method come from the PVLDB 2021 and
PVLDB 2024 papers listed in the
[NKDE references](references.md#network-kernel-density-visualization).

network_kde estimates density along lixel segments of a network.

~~~python
result = nsa.network_kde(
    network=network,
    points=events,
    bandwidth=500.0,
    lixel_length=50.0,
    method=nsa.NKDVMethod.LION,
    kernel=nsa.NKDVKernel.EPANECHNIKOV,
    normalize=False,
)
~~~

Parameters:

- bandwidth controls the network-distance support of the kernel.
- lixel_length controls output segmentation. Smaller values produce more rows
  and increase memory and computation.
- method selects the native evaluation strategy.
- normalize divides by the density integral over actual lixel lengths.
- weights supplies one finite non-negative value per event.

The public API currently enables only the Epanechnikov kernel. Although the
native ABI defines identifiers for additional kernels, they remain disabled
until all accelerated methods have numerical regression coverage.

Supported method names are BASIC, SHARED_SPD, ADA, IA, ADA_IA, and LION. ADA
and LION are the primary optimized choices exposed by the QGIS workflow.
Different methods should agree within numerical tolerance; benchmark both
speed and agreement on representative network structures before selecting a
production default.

Weighted example:

~~~python
result = nsa.network_kde(
    network=network,
    points=events,
    bandwidth=500.0,
    lixel_length=50.0,
    method="ADA",
    weights=severity,
)
~~~

Weights are applied during native kernel aggregation. At least one weight must
be positive. Tiny negative values caused by floating-point subtractive updates
are clamped to zero before results are returned.

## Network K-function

NS/ASPS and STARS are documented in the PVLDB 2022 and KDD 2026 papers listed
in the [Network K references](references.md#network-k-function).

network_k_function evaluates K(r) for ordered event pairs whose network
distance is within each threshold.

~~~python
result = nsa.network_k_function(
    network=network,
    points=events,
    r_min=0.0,
    r_max=5000.0,
    r_step=100.0,
    method=nsa.KFunctionMethod.STARS,
    simulations=99,
    alpha=0.05,
    envelope="global",
    random_state=2026,
)
~~~

Thresholds can instead be provided as a finite, non-negative, ordered
r_values array. Supported methods are SPS, NS, and STARS. STARS is generally
the preferred optimized path for large threshold sets, while the other
methods remain useful for comparison and validation. Performance depends on
network topology, threshold range, event distribution, and backend build, so
benchmark with the intended workload.

The estimator is:

~~~text
network_length * ordered_pair_count / (n_points * (n_points - 1))
~~~

At least two events are required.

## Spatiotemporal Network K-function

EAR, MTS, and EARTH are introduced in the KDD 2026 paper listed in the
[spatiotemporal reference](references.md#spatiotemporal-network-k-function).

spatiotemporal_network_k_function evaluates K(s, t) on spatial and temporal
threshold axes. Output arrays have shape (number of s values, number of t
values).

~~~python
result = nsa.spatiotemporal_network_k_function(
    network=network,
    points=events,
    s_values=[250.0, 500.0, 750.0],
    t_values=[1.0, 4.0, 7.0],
    method=nsa.SpatiotemporalKFunctionMethod.EARTH,
    simulations=19,
    envelope="pointwise",
)
~~~

Supported methods are RQS, SPS, EAR, MTS, and EARTH. EARTH is the default
optimized method. MTS and EARTH require evenly spaced s and t values; use the
min/max/step arguments when that requirement applies. Threshold units must
match the network and event time units respectively.

The estimator is:

~~~text
network_length * observed_time_span * ordered_pair_count
--------------------------------------------------------
              n_points * (n_points - 1)
~~~

Set return_xarray=True to return an xarray DataArray of the observed surface
directly. Otherwise the function returns SpatiotemporalKFunctionResult, whose
to_xarray method can export any available surface later.

## Monte Carlo Inference

K-function simulations sample events uniformly with respect to network length.
Spatiotemporal simulations also sample over the observed temporal range.

- simulations=0 disables theoretical curves and envelopes.
- k_theo is the mean Monte Carlo result, not an analytic closed-form curve.
- alpha controls the envelope significance level.
- random_state controls reproducibility; None uses the stable default seed.
- return_simulations=True retains every simulation and uses more memory.

Pointwise envelopes take quantiles independently at every threshold. They show
local departures but do not control error across the entire curve or surface.
Global envelopes use maximum absolute deviation from the simulation mean and
provide a single global Monte Carlo p-value.

For quick smoke tests, 5 simulations can confirm that a path works. It is not
adequate for fine significance levels. Use a simulation count justified by the
required p-value resolution and report the count, seed, alpha, and envelope
type with results.

## Progress and Cancellation

All analysis functions accept progress_callback. The callback receives:

~~~python
def progress(stage, completed, total):
    print(stage, completed, total)
    return cancel_requested
~~~

Return a truthy value to request cooperative cancellation. Cancellation is
checked at backend-defined safe points, so response is not necessarily
instantaneous inside a long native phase. A cancelled operation raises
FastNSACancelledError when the backend reports cancellation.

Callbacks may execute from native worker context. Keep them short, thread-safe,
and free of direct GUI updates; dispatch UI work to the host event loop.

## Performance Planning

The main growth factors are:

| Analysis | Important cost drivers |
| --- | --- |
| Mapping | Number of edges, point count, geometry complexity, method |
| NKDE | Lixel count, bandwidth, event count, topology, method |
| Network K | Event count, threshold count/range, simulations, method |
| ST K | Event count, s by t grid size, simulations, method |

Each simulation repeats a substantial part of K-function computation, so
runtime is approximately proportional to simulation count after fixed setup
costs. Avoid unnecessarily fine threshold grids during exploration. Reuse
network construction and mapped point sets across compatible analyses.

