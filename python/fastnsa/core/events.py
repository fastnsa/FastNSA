"""GeoDataFrame event mapping and reusable mapped-event accessors."""

from __future__ import annotations

import ctypes
import os
from time import perf_counter
from typing import Literal, Optional

import numpy as np

from fastnsa.bindings import capi
from fastnsa.runtime import NativePointSet


class PointEvents:
    """Map a point GeoDataFrame to a FastNSA network.

    Point geometries are reprojected to the network CRS when necessary, then
    mapped to native edge IDs and along-edge offsets. The instance records
    per-stage timings and owns the resulting :class:`NativePointSet`.

    Parameters
    ----------
    points_gdf : geopandas.GeoDataFrame
        Non-empty point geometries with a defined CRS.
    target_network : Network or NativeNetwork
        Open target network. Its CRS and coordinates should use a suitable
        projected distance unit.
    mapping_method : {"R-tree", "comp_a", "comp_e"}, default "R-tree"
        Native nearest-edge mapping strategy.
    comp_eps : float, default 10.0
        Approximation tolerance for ``comp_a``.
    n_threads : int, default 0
        Requested native mapping worker count. Zero uses the bounded automatic
        policy.
    time_col : str, optional
        Numeric or datetime-compatible event-time column.
    times : array-like, optional
        Explicit numeric times. Takes precedence over ``time_col``.
    time_unit, time_origin, time_mode : str, optional
        Descriptive metadata copied to spatiotemporal result objects. FastNSA
        does not convert numeric values based on these labels.

    Notes
    -----
    Non-numeric datetime columns are parsed as UTC Unix seconds. Numeric time
    columns are preserved. Close events before closing their network.

    See Also
    --------
    fastnsa.NativePointSet.map_xy : Native matching API and source publication.
    """

    def __init__(
        self,
        points_gdf,
        target_network,
        mapping_method: Literal["R-tree", "comp_a", "comp_e"] = "R-tree",
        comp_eps: float = 10.0,
        n_threads: int = 0,
        time_col: Optional[str] = None,
        times: Optional[np.ndarray] = None,
        time_unit: Optional[str] = None,
        time_origin: Optional[str] = None,
        time_mode: Optional[str] = None,
    ):
        total_started = perf_counter()
        self.network = target_network
        self.n_points = len(points_gdf)
        self.timings = {}
        self.time_unit = time_unit
        self.time_origin = time_origin
        self.time_mode = time_mode

        crs_started = perf_counter()
        if points_gdf.crs != self.network.crs:
            points_gdf = points_gdf.to_crs(self.network.crs)
        self.timings["crs_align_s"] = perf_counter() - crs_started

        xy_started = perf_counter()
        x_values = np.ascontiguousarray(points_gdf.geometry.x.values, dtype=np.float64)
        y_values = np.ascontiguousarray(points_gdf.geometry.y.values, dtype=np.float64)
        self.timings["xy_extract_s"] = perf_counter() - xy_started

        mapping_started = perf_counter()
        self._native = NativePointSet.map_xy(
            self.network,
            x_values,
            y_values,
            mapping_method=mapping_method,
            comp_eps=comp_eps,
            n_threads=n_threads,
        )
        self.timings["cpp_mapping_s"] = perf_counter() - mapping_started

        times_started = perf_counter()
        try:
            time_values = self._prepare_times(time_col, times, points_gdf)
            if time_values is not None:
                self._native.set_times(time_values)
        except Exception:
            self.close()
            raise
        self.timings["time_update_s"] = perf_counter() - times_started
        self.timings["total_s"] = perf_counter() - total_started

    @classmethod
    def from_mapped_arrays(
        cls,
        points_gdf,
        target_network,
        edge_ids,
        offsets,
        *,
        time_col: Optional[str] = None,
        times: Optional[np.ndarray] = None,
        time_unit: Optional[str] = None,
        time_origin: Optional[str] = None,
        time_mode: Optional[str] = None,
    ) -> "PointEvents":
        """Create events from mapping columns without nearest-edge matching.

        ``edge_ids`` and ``offsets`` must correspond to ``points_gdf`` row order
        and to the exact native edge order of ``target_network``. CRS alignment
        is still performed for consistency, but point geometry is not remapped.

        Parameters are otherwise equivalent to :class:`PointEvents`.
        """
        total_started = perf_counter()
        instance = cls.__new__(cls)
        instance.network = target_network
        instance.n_points = len(points_gdf)
        instance.timings = {}
        instance.time_unit = time_unit
        instance.time_origin = time_origin
        instance.time_mode = time_mode

        crs_started = perf_counter()
        if points_gdf.crs != instance.network.crs:
            points_gdf = points_gdf.to_crs(instance.network.crs)
        instance.timings["crs_align_s"] = perf_counter() - crs_started

        create_started = perf_counter()
        instance._native = NativePointSet.from_mapped(
            target_network, edge_ids, offsets
        )
        instance.timings["mapped_set_s"] = perf_counter() - create_started

        times_started = perf_counter()
        try:
            time_values = instance._prepare_times(time_col, times, points_gdf)
            if time_values is not None:
                instance._native.set_times(time_values)
        except Exception:
            instance.close()
            raise
        instance.timings["time_update_s"] = perf_counter() - times_started
        instance.timings["total_s"] = perf_counter() - total_started
        return instance

    @property
    def _ptr(self):
        return self._native._ptr

    def close(self):
        """Release the underlying native point set."""
        native = getattr(self, "_native", None)
        if native is not None:
            native.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def _prepare_times(self, time_col, times, points_gdf):
        if times is not None:
            time_values = np.ascontiguousarray(times, dtype=np.float64)
            self.time_mode = self.time_mode or "numeric"
        elif time_col is not None:
            if time_col not in points_gdf.columns:
                raise ValueError(f"`time_col` {time_col!r} not found in points data.")
            try:
                import pandas as pd
            except ImportError as exc:
                raise ImportError(
                    'Datetime conversion requires the "geo" extra: '
                    'pip install "fastnsa[geo]".'
                ) from exc
            column = points_gdf[time_col]
            if pd.api.types.is_numeric_dtype(column):
                time_values = np.ascontiguousarray(column.to_numpy(), dtype=np.float64)
                self.time_mode = self.time_mode or "numeric"
            else:
                timestamps = pd.to_datetime(column, errors="coerce", utc=True)
                if timestamps.isna().any():
                    raise ValueError(
                        f"Cannot parse some values in time column {time_col!r}."
                    )
                time_values = np.ascontiguousarray(
                    (timestamps.astype("int64") / 1e9).to_numpy(),
                    dtype=np.float64,
                )
                self.time_mode = self.time_mode or "datetime"
                self.time_unit = self.time_unit or "seconds"
                self.time_origin = self.time_origin or "Unix epoch (1970-01-01T00:00:00Z)"
        else:
            return None

        if time_values.ndim != 1 or len(time_values) != len(points_gdf):
            raise ValueError(
                "Time input must be 1D and have length equal to number of points."
            )
        if not np.isfinite(time_values).all():
            raise ValueError("Time input contains non-finite values (NaN or inf).")
        return time_values

    def dump_to_file(self, filename):
        """Write mapped edge IDs and offsets using the native text format."""
        return self._dump(filename, capi.lib.fastnsa_pointset_dump_text)

    def dump_to_file_st(self, filename):
        """Write mapped locations and elapsed times using the native text format.

        Times are relative to the earliest event and retain the input unit.
        Convert seconds to days explicitly before exporting if days are needed;
        older native libraries unconditionally divided exported times by 86400.
        """
        return self._dump(filename, capi.lib.fastnsa_pointset_dump_text_st)

    def _dump(self, filename, function):
        path = os.fspath(filename)
        status = function(self._ptr, ctypes.c_char_p(os.fsencode(path)))
        if status != 0:
            raise OSError(capi.last_error())
        return path

    def edge_ids_array(self) -> np.ndarray:
        """Return a copy of mapped native edge IDs."""
        return self._native.edge_ids_array()

    def offsets_array(self) -> np.ndarray:
        """Return a copy of along-edge offsets in network distance units."""
        return self._native.offsets_array()

    def times_array(self) -> np.ndarray:
        """Return a copy of numeric event times stored by the native backend."""
        return self._native.times_array()

    def to_mapping_records(self) -> list[dict[str, float | int]]:
        """Return mapping data as one plain dictionary per input point.

        Each record contains ``point_index``, ``edge_id``, ``offset``, and
        ``time``. The dictionaries are convenient for host-layer attribute
        creation but are less memory-efficient than the array accessors.
        """
        edge_ids = self.edge_ids_array()
        offsets = self.offsets_array()
        times = self.times_array()
        return [
            {
                "point_index": int(index),
                "edge_id": int(edge_ids[index]),
                "offset": float(offsets[index]),
                "time": float(times[index]),
            }
            for index in range(self.n_points)
        ]
