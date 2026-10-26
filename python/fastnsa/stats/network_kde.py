"""Native network kernel-density estimation and weighted lixel output."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal
import logging

import numpy as np

from fastnsa.bindings import capi
from fastnsa.core.types import NetworkKDEPointResult
from fastnsa.enums import enum_value
from fastnsa.exceptions import FastNSAValidationError, backend_error


logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from fastnsa.runtime import NativeNetwork, NativePointSet


def network_kde(
    *,
    network: NativeNetwork,
    points: NativePointSet,
    bandwidth: float,
    lixel_length: float,
    method: Literal["BASIC", "SHARED_SPD", "ADA", "IA", "ADA_IA", "LION"] = "LION",
    kernel: Literal["EPANECHNIKOV"] = "EPANECHNIKOV",
    normalize: bool = False,
    progress_callback=None,
    weights=None,
) -> NetworkKDEPointResult:
    """
    Estimate network kernel density on lixel segments.

    Parameters
    ----------
    network : Network or NativeNetwork
        Open target network. Edge lengths define all distance units.
    points : PointEvents or NativePointSet
        Non-empty event locations mapped to ``network``.
    bandwidth : float
        Kernel bandwidth in network distance units.
    lixel_length : float
        Length of lixel segmentation used for evaluation.
    method : {"BASIC", "SHARED_SPD", "ADA", "IA", "ADA_IA", "LION"}, default "LION"
        Native NKDV algorithm. Stable :class:`NKDVMethod` values are accepted.
    kernel : {"EPANECHNIKOV"}, default "EPANECHNIKOV"
        Kernel family. Other kernel families remain disabled until their
        accelerated backend implementations have numerical regression tests.
    normalize : bool, default False
        Divide values by their network integral, using each lixel's actual
        segment length. This is not a simple array-sum normalization.
    progress_callback : callable, optional
        Called as ``callback(stage, completed, total)``. Return true to request
        cooperative cancellation.
    weights : array-like, optional
        One finite non-negative event weight per point. At least one weight must
        be positive. Weights are applied inside native kernel aggregation.

    Returns
    -------
    NetworkKDEPointResult
        Lixel edge IDs, offsets, density values, and analysis metadata.

    Raises
    ------
    ValueError
        If bandwidth, lixel length, kernel, or object compatibility is invalid.
    FastNSAValidationError
        If the method or weights are invalid.
    RuntimeError
        If objects are closed or weighted NKDV is unavailable in the backend.
    NativeBackendError
        If native computation fails or is cancelled.

    Notes
    -----
    Tiny negative values introduced by floating-point subtractive updates are
    clamped to zero. Smaller lixels create finer output but increase the number
    of returned segments and typically increase computation time.

    References
    ----------
    Tsz Nam Chan et al. "Fast Augmentation Algorithms for Network Kernel
    Density Visualization." PVLDB 14(9), 2021.
    https://doi.org/10.14778/3461535.3461540

    Tsz Nam Chan et al. "LION: Fast and High-Resolution Network Kernel Density
    Visualization." PVLDB 17(6), 2024.
    https://doi.org/10.14778/3648160.3648168
    """
    if network is None or points is None:
        raise ValueError("Both `network` and `points` must be provided.")
    if not getattr(network, "_ptr", None) or not getattr(points, "_ptr", None):
        raise RuntimeError("The network and point set must remain open during analysis.")
    if getattr(getattr(points, "network", None), "_ptr", None) != network._ptr:
        raise ValueError("The point set belongs to a different network.")
    if points.n_points <= 0 or network.total_length <= 0.0:
        raise ValueError("The network and point set must be non-empty.")
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError("`bandwidth` must be positive.")
    if not np.isfinite(lixel_length) or lixel_length <= 0:
        raise ValueError("`lixel_length` must be positive.")

    method_key = str(enum_value(method)).upper()
    kernel_key = str(enum_value(kernel)).upper()
    method_map = {
        "BASIC": capi.FASTNSA_NKDV_BASIC,
        "SHARED_SPD": capi.FASTNSA_NKDV_SHARED_SPD,
        "ADA": capi.FASTNSA_NKDV_ADA,
        "IA": capi.FASTNSA_NKDV_IA,
        "ADA_IA": capi.FASTNSA_NKDV_ADA_IA,
        "LION": capi.FASTNSA_NKDV_LION,
    }
    kernel_map = {
        "GAUSSIAN": capi.FASTNSA_NKDV_KERNEL_GAUSSIAN,
        "TRIANGULAR": capi.FASTNSA_NKDV_KERNEL_TRIANGULAR,
        "EPANECHNIKOV": capi.FASTNSA_NKDV_KERNEL_EPANECHNIKOV,
        "QUARTIC": capi.FASTNSA_NKDV_KERNEL_QUARTIC,
    }

    if method_key not in method_map:
        raise FastNSAValidationError(f"Unsupported NKDV method: {method!r}.")
    if kernel_key not in kernel_map:
        raise ValueError(f"Unsupported NKDV kernel: {kernel!r}.")
    if kernel_key != "EPANECHNIKOV":
        raise ValueError("The current NKDV backend only supports the Epanechnikov kernel.")

    n_lixels = capi.lib.fastnsa_nkdv_count_lixels(network._ptr, float(lixel_length))
    if n_lixels <= 0:
        raise backend_error(capi.last_error())

    edge_ids = np.empty(n_lixels, dtype=np.int32)
    offsets_n1 = np.empty(n_lixels, dtype=np.float64)
    offsets_n2 = np.empty(n_lixels, dtype=np.float64)
    density = np.empty(n_lixels, dtype=np.float64)

    native_callback = capi.make_progress_callback(progress_callback)
    callback_pointer = native_callback or capi.PROGRESS_CALLBACK()
    weight_values = None
    if weights is not None:
        weight_values = np.ascontiguousarray(weights, dtype=np.float64)
        if weight_values.ndim != 1 or len(weight_values) != points.n_points:
            raise FastNSAValidationError(
                "`weights` must be a 1D array matching the point count."
            )
        if not np.isfinite(weight_values).all() or np.any(weight_values < 0.0):
            raise FastNSAValidationError("`weights` must be finite and non-negative.")
        if not np.any(weight_values > 0.0):
            raise FastNSAValidationError("At least one event weight must be positive.")
        if not capi.HAS_WEIGHTED_NKDV:
            raise RuntimeError("The loaded native backend does not support event weights.")

    use_progress_api = native_callback is not None and capi.HAS_PROGRESS_NKDV
    if weight_values is not None:
        function = capi.lib.fastnsa_compute_nkdv_weighted
    elif use_progress_api:
        function = capi.lib.fastnsa_compute_nkdv_with_progress
    else:
        function = capi.lib.fastnsa_compute_nkdv
    arguments = [
        network._ptr,
        points._ptr,
        method_map[method_key],
        float(lixel_length),
        kernel_map[kernel_key],
        float(bandwidth),
        int(n_lixels),
        edge_ids.ctypes.data_as(capi.P_INT),
        offsets_n1.ctypes.data_as(capi.P_DOUBLE),
        offsets_n2.ctypes.data_as(capi.P_DOUBLE),
        density.ctypes.data_as(capi.P_DOUBLE),
    ]
    if weight_values is not None:
        arguments = [
            network._ptr,
            points._ptr,
            weight_values.ctypes.data_as(capi.P_DOUBLE),
            len(weight_values),
        ] + arguments[2:] + [callback_pointer, None]
    elif use_progress_api:
        arguments.extend([native_callback, None])
    ret = function(*arguments)
    if ret != 0:
        raise backend_error(capi.last_error())

    # Kernel density estimates are non-negative. Native implementations may
    # leave tiny negative round-off values after subtractive updates.
    np.maximum(density, 0.0, out=density)

    if normalize:
        edge_lengths = offsets_n1 + offsets_n2
        segment_starts = np.floor(offsets_n1 / float(lixel_length)) * float(lixel_length)
        segment_lengths = np.minimum(float(lixel_length), edge_lengths - segment_starts)
        segment_lengths = np.maximum(segment_lengths, 0.0)
        density_integral = float(np.sum(density * segment_lengths))
        if density_integral > 0:
            density = density / density_integral

    result = NetworkKDEPointResult(
        edge_ids=edge_ids,
        offsets_n1=offsets_n1,
        offsets_n2=offsets_n2,
        density=density,
        method=method_key,
        kernel=kernel_key,
        bandwidth=float(bandwidth),
        lixel_length=float(lixel_length),
        n_points=points.n_points,
        network_length=network.total_length,
        normalized=normalize,
        normalization="network_integral" if normalize else None,
        weighted=weight_values is not None,
        total_weight=(
            float(np.sum(weight_values)) if weight_values is not None else float(points.n_points)
        ),
    )
    logger.info(
        "Computed network KDE",
        extra={
            "fastnsa_stage": "network_kde",
            "fastnsa_method": method_key,
            "fastnsa_points": points.n_points,
            "fastnsa_lixels": n_lixels,
            "fastnsa_weighted": weight_values is not None,
        },
    )
    return result
