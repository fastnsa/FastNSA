"""Small deterministic datasets for examples and smoke tests."""

import struct

import numpy as np

from fastnsa.runtime import NativeNetwork, NativePointSet


def make_line_network(length: float = 1000.0, *, crs: str = "EPSG:3857") -> NativeNetwork:
    """Create a deterministic one-edge network for examples and smoke tests.

    Parameters
    ----------
    length : float, default 1000.0
        Positive edge length and x-coordinate of the ending node.
    crs : str, default "EPSG:3857"
        Metadata attached to the network; no coordinate conversion is applied.

    Returns
    -------
    NativeNetwork
        Open two-node network with one horizontal LineString edge.

    Notes
    -----
    This synthetic network is not intended as a scientific benchmark dataset.
    """
    length = float(length)
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("`length` must be finite and positive.")
    line_wkb = struct.pack("<BII4d", 1, 2, 2, 0.0, 0.0, length, 0.0)
    return NativeNetwork(
        edges=np.array([[0, 1, length]], dtype=np.float64),
        offsets=np.array([0, 1, 2], dtype=np.int32),
        indices=np.array([0, 0], dtype=np.int32),
        wkb_blob=line_wkb,
        wkb_offsets=np.array([0, len(line_wkb)], dtype=np.int32),
        crs=crs,
    )


def make_line_events(network: NativeNetwork, offsets=None, *, times=None) -> NativePointSet:
    """Create mapped events on a one-edge sample network.

    Parameters
    ----------
    network : NativeNetwork
        Open network containing exactly one edge.
    offsets : array-like, optional
        Along-edge event offsets. Defaults to 10%, 40%, and 80% of total length.
    times : array-like, optional
        Finite numeric times matching ``offsets``.

    Returns
    -------
    NativePointSet
        Open point set whose edge IDs are all zero.
    """
    if network.num_edges != 1:
        raise ValueError("`make_line_events` requires a single-edge network.")
    if offsets is None:
        offsets = np.asarray(network.total_length) * np.array([0.1, 0.4, 0.8])
    offset_values = np.ascontiguousarray(offsets, dtype=np.float64)
    return NativePointSet.from_mapped(
        network,
        np.zeros(len(offset_values), dtype=np.int32),
        offset_values,
        times=times,
    )
