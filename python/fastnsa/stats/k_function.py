"""Network and spatiotemporal K-functions with Monte Carlo inference."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Optional
import logging
import numpy as np

from fastnsa.core.types import KFunctionResult, SpatiotemporalKFunctionResult
from fastnsa.bindings import capi
from fastnsa.enums import enum_value
from fastnsa.exceptions import FastNSAValidationError, backend_error


logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from fastnsa.runtime import NativeNetwork, NativePointSet


def network_k_function(
    *,
    network: NativeNetwork,
    points: NativePointSet,
    r_min: float = 0.0,
    r_max: float = 2000.0,
    r_step: int = 50,
    r_values: Optional[np.ndarray] = None,
    method: Literal["SPS", "NS", "STARS"] = "NS",
    simulations: int = 0, 
    alpha: float = 0.05,
    return_simulations: bool = False,
    epsilon: float = 0.1,
    progress_callback=None,
    random_state: Optional[int] = 5489,
    envelope: str = "pointwise",
) -> KFunctionResult:
    """
    Compute a network K-function over ordered event pairs.

    Parameters
    ----------
    network : Network or NativeNetwork
        Open network on which shortest-path distances are measured.
    points : PointEvents or NativePointSet
        At least two event locations mapped to ``network``. A point set from a
        different network is rejected.
    r_min, r_max, r_step : float, optional
        Range and step size for network distances.
        Ignored if `r_values` is provided.
        Values use the same distance unit as network edge lengths. The generated
        range includes the upper endpoint when it lies on the step grid.
    r_values : array-like, optional
        Explicit finite, non-negative, non-decreasing thresholds. When given,
        ``r_min``, ``r_max``, and ``r_step`` are ignored.
    method : {"SPS", "NS", "STARS"}, default "NS"
        Native computation algorithm. Stable :class:`KFunctionMethod` values
        are also accepted.
    simulations : int, default 0
        Number of complete-spatial-randomness simulations. Zero disables the
        theoretical curve, envelope, and p-value.
    alpha : float, default 0.05
        Envelope significance level in the open interval (0, 1).
    return_simulations : bool, default False
        Retain the full ``(simulations, n_thresholds)`` matrix. Simulations are
        still computed for inference when this is false, but not returned.
    epsilon : float, default 0.1
        STARS source-event sampling control: at most ``ceil(1 / epsilon)``
        source events are sampled per occupied edge. Smaller values increase
        sampling work and usually improve accuracy; this is not a guaranteed
        percentage error bound. NS and SPS do not use epsilon for counting.
    progress_callback : callable, optional
        Called as ``callback(stage, completed, total)``. Return a truthy value
        to request cooperative cancellation.
    random_state : int or None, optional
        Monte Carlo seed. ``None`` uses the stable default seed 5489.
    envelope : {"pointwise", "global"}, default "pointwise"
        Pointwise quantiles or a simultaneous maximum-deviation envelope.

    Returns
    -------
    KFunctionResult
        Observed values, analysis metadata, and optional Monte Carlo outputs.

    Raises
    ------
    ValueError
        If objects are incompatible or thresholds and simulation parameters are
        invalid.
    RuntimeError
        If native objects are closed or a requested backend feature is absent.
    NativeBackendError
        If native computation fails.
    FastNSACancelledError
        If cooperative cancellation completes.

    Notes
    -----
    Values use the estimator ``network_length * ordered_pair_count /
    (n_points * (n_points - 1))``. Monte Carlo points are sampled uniformly
    with respect to network length.
    ``k_theo`` is the simulation mean under CSR on the network, not an analytic
    curve. A global p-value is calculated from maximum absolute deviation for
    either envelope choice. If ``simulations == 0``, Monte Carlo fields are
    ``None``.

    References
    ----------
    Tsz Nam Chan et al. "Fast Network K-function-based Spatial Analysis."
    PVLDB 15(11), 2022. https://doi.org/10.14778/3551793.3551836

    Hongwei Ye et al. "STARS: A Sampling and Threshold Sharing Solution for
    Network K-function Analytics." KDD 2026.
    https://doi.org/10.1145/3770855.3817709
    """
    if not network or not points:
        raise ValueError("Both `network` and `points` must be provided and non-empty.")
    if not getattr(network, "_ptr", None) or not getattr(points, "_ptr", None):
        raise RuntimeError("The network and point set must remain open during analysis.")
    if getattr(getattr(points, "network", None), "_ptr", None) != network._ptr:
        raise ValueError("The point set belongs to a different network.")
    if points.n_points < 2 or network.total_length <= 0.0:
        raise ValueError("Network K-function requires a non-empty network and at least two points.")

    if simulations < 0:
        raise ValueError("`simulations` must be non-negative.")

    if not (0.0 < alpha < 1.0):
        raise ValueError("`alpha` must be in (0, 1).")
    
    if r_values is None:
        if r_step <= 0:
            raise ValueError("`r_step` must be positive.")
        if r_max <= r_min:
            raise ValueError("`r_max` must be larger than `r_min`.")

        r = np.arange(r_min, r_max + r_step, r_step, dtype=np.float64)
    else:
        r = np.ascontiguousarray(r_values, dtype=np.float64)
        if r.ndim != 1 or r.size == 0:
            raise ValueError("`r_values` must be a non-empty 1D array.")
    if not np.isfinite(r).all() or np.any(r < 0.0) or np.any(np.diff(r) < 0.0):
        raise ValueError("K-function thresholds must be finite, non-negative, and ordered.")

    n_r = r.size
    n_points = points.n_points

    method_map = {
        "NS": capi.FASTNSA_K_NS,
        "SPS": capi.FASTNSA_K_SPS,
        "STARS": capi.FASTNSA_K_STARS,
    }
    method_key = str(enum_value(method)).upper()
    if method_key not in method_map:
        raise ValueError(f"Unsupported network K-function method: {method!r}.")
    method_flag = method_map[method_key]


    k_obs = np.zeros(n_r, dtype=np.float64)

    k_sim_all = np.zeros((simulations, n_r), dtype=np.float64) if simulations > 0 else None
    k_sim = k_sim_all if return_simulations else None

    seed = _normalize_random_state(random_state)
    native_callback = capi.make_progress_callback(progress_callback)
    callback_pointer = native_callback or capi.PROGRESS_CALLBACK()
    use_seeded_api = capi.HAS_SEEDED_K_FUNCTION
    if not use_seeded_api and seed != 5489:
        raise RuntimeError(
            "The loaded native backend does not support custom random seeds."
        )
    use_progress_api = native_callback is not None and capi.HAS_PROGRESS_K_FUNCTION
    if use_seeded_api:
        function = capi.lib.fastnsa_compute_k_function_seeded
    elif use_progress_api:
        function = capi.lib.fastnsa_compute_k_function_with_progress
    else:
        function = (
            capi.lib.fastnsa_compute_k_function_checked
            if capi.HAS_CHECKED_K_FUNCTION
            else capi.lib.fastnsa_compute_k_function
        )
    arguments = [
        network._ptr,
        points._ptr,
        r.ctypes.data_as(capi.P_DOUBLE),
        n_r,
        method_flag,
        simulations,
        k_obs.ctypes.data_as(capi.P_DOUBLE),
        k_sim_all.ctypes.data_as(capi.P_DOUBLE) if k_sim_all is not None else None,
        epsilon,
    ]
    if use_seeded_api:
        arguments.extend([seed, callback_pointer, None])
    elif use_progress_api:
        arguments.extend([native_callback, None])
    status = function(*arguments)
    if (capi.HAS_CHECKED_K_FUNCTION or native_callback is not None) and status != 0:
        raise backend_error(capi.last_error())

    k_theo, k_lo, k_hi, p_value, envelope_method = _monte_carlo_inference(
        k_obs, k_sim_all, alpha, envelope
    )

    result = KFunctionResult(
            r=r,
            k_obs=k_obs,
            k_theo=k_theo,
            k_lo=k_lo,
            k_hi=k_hi,
            k_sim=k_sim if return_simulations else None,
            method=method_key,
            n_points=n_points,
            intensity=n_points / network.total_length,
            network_length=network.total_length,
            random_seed=seed,
            envelope_method=envelope_method,
            confidence_level=(1.0 - alpha) if simulations > 0 else None,
            p_value=p_value,
        )
    logger.info(
        "Computed network K-function",
        extra={
            "fastnsa_stage": "network_k_function",
            "fastnsa_method": method_key,
            "fastnsa_points": n_points,
            "fastnsa_simulations": simulations,
        },
    )
    return result


def spatiotemporal_network_k_function(
    *,
    network: NativeNetwork,
    points: NativePointSet,
    s_min: float = 0.0,
    s_max: float = 2000.0,
    s_step: int = 50,
    t_min: float = 0.0,
    t_max: float = 2000.0,
    t_step: int = 50,
    s_values: Optional[np.ndarray] = None,
    t_values: Optional[np.ndarray] = None,
    method: Literal["RQS", "SPS", "EAR", "MTS", "EARTH"] = "EARTH",
    simulations: int = 0,
    alpha: float = 0.05,
    return_simulations: bool = False,
    progress_callback=None,
    random_state: Optional[int] = 5489,
    envelope: str = "pointwise",
    return_xarray: bool = False,
):
    """
    Compute a spatiotemporal network K-function surface.

    This function evaluates K(s, t), where:
    - `s` is the network distance threshold.
    - `t` is the temporal threshold.

    Parameters
    ----------
    network : Network or NativeNetwork
        Open target network.
    points : PointEvents or NativePointSet
        At least two events mapped to ``network`` with finite numeric times and
        a positive observed time span.
    s_min, s_max, s_step : float, optional
        Spatial threshold range. Ignored if `s_values` is provided.
    t_min, t_max, t_step : float, optional
        Temporal threshold range. Ignored if `t_values` is provided.
        Spatial values use network distance units; temporal values use the same
        unit as event times.
    s_values : array-like, optional
        Explicit finite, non-negative, non-decreasing spatial thresholds.
    t_values : array-like, optional
        Explicit finite, non-negative, non-decreasing temporal thresholds.
    method : {"RQS", "SPS", "EAR", "MTS", "EARTH"}, optional
        Backend algorithm. Default is `"EARTH"`.
    simulations : int, default 0
        Number of network/time CSR simulations. Zero disables Monte Carlo
        reference and envelope surfaces.
    alpha : float, default 0.05
        Envelope significance level in the open interval (0, 1).
    return_simulations : bool, default False
        Retain the full ``(simulations, n_s, n_t)`` simulation array.
    progress_callback : callable, optional
        Called as ``callback(stage, completed, total)``. Return true to request
        cooperative cancellation.
    random_state : int or None, optional
        Monte Carlo seed. ``None`` uses the stable default seed 5489.
    envelope : {"pointwise", "global"}, default "pointwise"
        Pointwise quantiles or a simultaneous maximum-deviation envelope.
    return_xarray : bool, default False
        Return the observed surface directly as ``xarray.DataArray``. This
        requires ``fastnsa[xarray]`` and omits access to the full result object.

    Returns
    -------
    SpatiotemporalKFunctionResult or xarray.DataArray
        Contains observed `k_obs` with shape `(n_s, n_t)` and optional
        envelope surfaces (`k_lo`, `k_hi`) and simulations (`k_sim`), unless
        ``return_xarray`` requests the observed DataArray.

    Raises
    ------
    ValueError
        If events, time span, thresholds, or simulation parameters are invalid.
        MTS and EARTH require at least two evenly spaced thresholds per axis.
    ImportError
        If ``return_xarray`` is true and xarray is unavailable.
    NativeBackendError
        If native computation fails or is cancelled.

    Notes
    -----
    Output arrays are indexed as ``[s_index, t_index]``.
    Values use the estimator ``network_length * time_span * ordered_pair_count /
    (n_points * (n_points - 1))``. Monte Carlo events are sampled uniformly
    with respect to network length and over the observed temporal range.
    All five methods use the input time unit without implicit conversion.
    Multiplying both event times and temporal thresholds by a positive factor
    multiplies observed, reference, and envelope K values by that factor.
    Shifting the time origin leaves K unchanged. Compare results only after
    aligning time units, thresholds, and Monte Carlo seeds.

    References
    ----------
    Tsz Nam Chan et al. "EARTH: Accelerating Spatiotemporal Network K-function-
    based Analytics." KDD 2026. https://doi.org/10.1145/3770855.3817706

    Examples
    --------
    >>> result = spatiotemporal_network_k_function(
    ...     network=net,
    ...     points=events,
    ...     s_values=np.array([50.0, 100.0, 150.0]),
    ...     t_values=np.array([1.0, 4.0, 7.0]),
    ...     method="EARTH",
    ...     simulations=99,
    ... )
    >>> result.k_obs.shape
    (3, 3)
    """
    if not network or not points:
        raise ValueError("Both `network` and `points` must be provided and non-empty.")
    if not getattr(network, "_ptr", None) or not getattr(points, "_ptr", None):
        raise RuntimeError("The network and point set must remain open during analysis.")
    if getattr(getattr(points, "network", None), "_ptr", None) != network._ptr:
        raise ValueError("The point set belongs to a different network.")
    if points.n_points < 2 or network.total_length <= 0.0:
        raise ValueError(
            "Spatiotemporal K-function requires a non-empty network and at least two points."
        )
    if simulations < 0:
        raise ValueError("`simulations` must be non-negative.")
    if not (0.0 < alpha < 1.0):
        raise ValueError("`alpha` must be in (0, 1).")

    if s_values is None:
        if s_step <= 0 or s_max <= s_min:
            raise ValueError("Invalid s-range: require s_step>0 and s_max>s_min.")
        s = np.arange(s_min, s_max + s_step, s_step, dtype=np.float64)
    else:
        s = np.ascontiguousarray(s_values, dtype=np.float64)
        if s.ndim != 1 or s.size == 0:
            raise ValueError("`s_values` must be a non-empty 1D array.")

    if t_values is None:
        if t_step <= 0 or t_max <= t_min:
            raise ValueError("Invalid t-range: require t_step>0 and t_max>t_min.")
        t = np.arange(t_min, t_max + t_step, t_step, dtype=np.float64)
    else:
        t = np.ascontiguousarray(t_values, dtype=np.float64)
        if t.ndim != 1 or t.size == 0:
            raise ValueError("`t_values` must be a non-empty 1D array.")
    if not np.isfinite(s).all() or np.any(s < 0.0) or np.any(np.diff(s) < 0.0):
        raise ValueError("Spatial thresholds must be finite, non-negative, and ordered.")
    if not np.isfinite(t).all() or np.any(t < 0.0) or np.any(np.diff(t) < 0.0):
        raise ValueError("Temporal thresholds must be finite, non-negative, and ordered.")

    n_s = s.size
    n_t = t.size
    n_points = points.n_points
    point_times = points.times_array()
    time_span = float(np.ptp(point_times))
    if not np.isfinite(time_span) or time_span <= 0.0:
        raise ValueError(
            "Spatiotemporal K-function requires event times with a positive time span."
        )

    method_map = {
        "RQS": capi.FASTNSA_K_RQS,
        "SPS": capi.FASTNSA_K_SPS,
        "EAR": capi.FASTNSA_K_EAR,
        "MTS": capi.FASTNSA_K_MTS,
        "EARTH": capi.FASTNSA_K_EARTH,
    }
    method_key = str(enum_value(method)).upper()
    if method_key not in method_map:
        raise ValueError(f"Unsupported spatiotemporal K-function method: {method!r}.")
    method_flag = method_map[method_key]
    if method_key in {"MTS", "EARTH"}:
        if len(s) < 2 or len(t) < 2:
            raise ValueError(f"{method_key} requires at least two thresholds per axis.")
        if not np.allclose(np.diff(s), np.diff(s)[0]) or not np.allclose(
            np.diff(t),
            np.diff(t)[0],
        ):
            raise ValueError(f"{method_key} requires evenly spaced thresholds.")

    k_obs = np.zeros((n_s, n_t), dtype=np.float64)
    k_sim_all = np.zeros((simulations, n_s, n_t), dtype=np.float64) if simulations > 0 else None
    k_sim = k_sim_all if return_simulations else None

    seed = _normalize_random_state(random_state)
    native_callback = capi.make_progress_callback(progress_callback)
    callback_pointer = native_callback or capi.PROGRESS_CALLBACK()
    use_seeded_api = capi.HAS_SEEDED_ST_K_FUNCTION
    if not use_seeded_api and seed != 5489:
        raise RuntimeError(
            "The loaded native backend does not support custom random seeds."
        )
    use_progress_api = native_callback is not None and capi.HAS_PROGRESS_ST_K_FUNCTION
    if use_seeded_api:
        function = capi.lib.fastnsa_compute_spatiotemporal_k_function_seeded
    elif use_progress_api:
        function = capi.lib.fastnsa_compute_spatiotemporal_k_function_with_progress
    else:
        function = (
            capi.lib.fastnsa_compute_spatiotemporal_k_function_checked
            if capi.HAS_CHECKED_ST_K_FUNCTION
            else capi.lib.fastnsa_compute_spatiotemporal_k_function
        )
    arguments = [
        network._ptr,
        points._ptr,
        s.ctypes.data_as(capi.P_DOUBLE),
        n_s,
        t.ctypes.data_as(capi.P_DOUBLE),
        n_t,
        method_flag,
        simulations,
        k_obs.ctypes.data_as(capi.P_DOUBLE),
        k_sim_all.ctypes.data_as(capi.P_DOUBLE) if k_sim_all is not None else None,
    ]
    if use_seeded_api:
        arguments.extend([seed, callback_pointer, None])
    elif use_progress_api:
        arguments.extend([native_callback, None])
    status = function(*arguments)
    if (capi.HAS_CHECKED_ST_K_FUNCTION or native_callback is not None) and status != 0:
        raise backend_error(capi.last_error())

    k_theo, k_lo, k_hi, p_value, envelope_method = _monte_carlo_inference(
        k_obs, k_sim_all, alpha, envelope
    )

    result = SpatiotemporalKFunctionResult(
        s=s,
        t=t,
        k_obs=k_obs,
        k_theo=k_theo,
        k_lo=k_lo,
        k_hi=k_hi,
        k_sim=k_sim if return_simulations else None,
        method=method_key,
        n_points=n_points,
        intensity=n_points / network.total_length,
        network_length=network.total_length,
        time_unit=getattr(points, "time_unit", None),
        time_origin=getattr(points, "time_origin", None),
        time_mode=getattr(points, "time_mode", None),
        random_seed=seed,
        envelope_method=envelope_method,
        confidence_level=(1.0 - alpha) if simulations > 0 else None,
        p_value=p_value,
    )
    logger.info(
        "Computed spatiotemporal network K-function",
        extra={
            "fastnsa_stage": "spatiotemporal_network_k_function",
            "fastnsa_method": method_key,
            "fastnsa_points": n_points,
            "fastnsa_simulations": simulations,
        },
    )
    return result.to_xarray() if return_xarray else result


def _monte_carlo_inference(observed, simulations, alpha: float, envelope: str):
    envelope_method = str(envelope).casefold()
    if envelope_method not in {"pointwise", "global"}:
        raise FastNSAValidationError("`envelope` must be 'pointwise' or 'global'.")
    if simulations is None or len(simulations) == 0:
        return None, None, None, None, None

    theoretical = np.mean(simulations, axis=0)
    simulation_deviation = np.max(
        np.abs(simulations - theoretical), axis=tuple(range(1, simulations.ndim))
    )
    observed_deviation = float(np.max(np.abs(observed - theoretical)))
    p_value = float(
        (1 + np.count_nonzero(simulation_deviation >= observed_deviation)) /
        (len(simulations) + 1)
    )

    if envelope_method == "pointwise":
        lower = np.quantile(simulations, alpha / 2.0, axis=0)
        upper = np.quantile(simulations, 1.0 - alpha / 2.0, axis=0)
    else:
        critical = float(np.quantile(simulation_deviation, 1.0 - alpha, method="higher"))
        lower = theoretical - critical
        upper = theoretical + critical
    return theoretical, lower, upper, p_value, envelope_method


def _normalize_random_state(random_state: Optional[int]) -> int:
    if random_state is None:
        return 5489
    if isinstance(random_state, bool) or not isinstance(random_state, (int, np.integer)):
        raise FastNSAValidationError("`random_state` must be an integer or None.")
    if random_state < 0 or random_state > np.iinfo(np.uint32).max:
        raise FastNSAValidationError("`random_state` must be in the uint32 range.")
    return int(random_state)
