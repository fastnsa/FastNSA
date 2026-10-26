"""Minimal installed-wheel diagnostic used by release automation."""

import struct

import numpy as np

from fastnsa.runtime import NativeNetwork, NativePointSet
from fastnsa.datasets import make_line_events, make_line_network
from fastnsa.stats.k_function import network_k_function
from fastnsa.stats.network_kde import network_kde


def main() -> None:
    """Run mapping, K-function, NKDE, and sample-data installation checks.

    The function is intentionally silent on success and raises ``RuntimeError``
    when a native result is missing or numerically invalid. It is invoked by
    ``python -m fastnsa._selftest`` in wheel release verification.
    """
    wkb = struct.pack("<BII4d", 1, 2, 2, 0.0, 0.0, 100.0, 0.0)
    with NativeNetwork(
        edges=np.array([[0, 1, 100.0]], dtype=np.float64),
        offsets=np.array([0, 1, 2], dtype=np.int32),
        indices=np.array([0, 0], dtype=np.int32),
        wkb_blob=wkb,
        wkb_offsets=np.array([0, len(wkb)], dtype=np.int32),
        crs="EPSG:3857",
    ) as network:
        with NativePointSet.map_xy(
            network, [10.0, 40.0, 80.0], [1.0, -1.0, 2.0]
        ) as points:
            mapped = points.offsets_array()
            if not np.allclose(mapped, [10.0, 40.0, 80.0]):
                raise RuntimeError(f"Unexpected nearest-edge offsets: {mapped!r}")
            k_result = network_k_function(
                network=network,
                points=points,
                r_values=np.array([10.0, 50.0]),
                method="NS",
                simulations=1,
                random_state=1234,
            )
            if not np.isfinite(k_result.k_obs).all():
                raise RuntimeError("K-function self-test returned non-finite values.")
            kde_result = network_kde(
                network=network,
                points=points,
                bandwidth=40.0,
                lixel_length=20.0,
                method="ADA",
            )
            if not np.isfinite(kde_result.density).all():
                raise RuntimeError("NKDV self-test returned non-finite values.")

    with make_line_network(100.0) as network:
        with make_line_events(network) as points:
            if points.n_points != 3:
                raise RuntimeError("Built-in sample events are invalid.")


if __name__ == "__main__":
    main()
