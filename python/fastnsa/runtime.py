"""Lightweight NumPy/ctypes runtime for the FastNSA native library."""

from __future__ import annotations

from typing import Any

import numpy as np

from fastnsa.bindings import capi
from fastnsa.enums import enum_value
from fastnsa.exceptions import FastNSAValidationError, backend_error


class NativeNetwork:
    """Own a native network created directly from topology and geometry arrays.

    This is the lowest-dependency network interface. Input arrays are converted
    to contiguous native dtypes and retained for the lifetime of the object.
    Edge row order defines the edge IDs used by :class:`NativePointSet`.

    Parameters
    ----------
    edges : array-like, shape (n_edges, 3)
        Rows contain ``(start_node, end_node, length)``. Node indices and edge
        lengths must be finite; lengths and geometries must use the same planar
        coordinate unit expected by the analyses.
    offsets, indices : array-like
        CSR adjacency arrays. ``offsets`` has ``n_nodes + 1`` entries and its
        final value must equal ``len(indices)``. Values are converted to int32.
    wkb_blob : bytes, optional
        Concatenated WKB edge geometries in the same order as ``edges``.
    wkb_offsets : array-like, shape (n_edges + 1,), optional
        Byte offsets delimiting each geometry in ``wkb_blob``. The final entry
        must equal ``len(wkb_blob)``. Used when packed coordinates are absent.
    coordinates : array-like, shape (n_vertices, 2), optional
        Packed two-dimensional edge vertices. This is the preferred path for
        host integrations because it avoids serializing and parsing WKB.
    geometry_offsets : array-like, shape (n_edges + 1,), optional
        Vertex offsets delimiting each edge in ``coordinates``.
    crs : object, optional
        CRS metadata retained for host integrations. The native backend does
        not reproject coordinates.

    Notes
    -----
    The object owns a native allocation. Use it as a context manager or call
    :meth:`close` when finished. Point sets that reference this network must be
    closed before the network.
    """

    def __init__(
        self,
        *,
        edges: np.ndarray,
        offsets: np.ndarray,
        indices: np.ndarray,
        wkb_blob: bytes = b"",
        wkb_offsets: np.ndarray | None = None,
        coordinates: np.ndarray | None = None,
        geometry_offsets: np.ndarray | None = None,
        crs: Any = None,
    ):
        self._ptr = None
        edges_array = np.ascontiguousarray(edges, dtype=np.float64)
        offsets_array = np.ascontiguousarray(offsets, dtype=np.int32)
        indices_array = np.ascontiguousarray(indices, dtype=np.int32)
        use_coordinates = coordinates is not None or geometry_offsets is not None
        if use_coordinates and (coordinates is None or geometry_offsets is None):
            raise ValueError("coordinates and geometry_offsets must be provided together.")

        if edges_array.ndim != 2 or edges_array.shape[1] != 3:
            raise FastNSAValidationError("edges must have shape (n_edges, 3).")
        if len(edges_array) == 0 or not np.isfinite(edges_array).all():
            raise ValueError("edges must be non-empty and contain only finite values.")
        if offsets_array.ndim != 1 or len(offsets_array) < 2:
            raise ValueError("offsets must be a 1D CSR offset array.")
        if indices_array.ndim != 1:
            raise ValueError("indices must be a 1D CSR index array.")
        if int(offsets_array[-1]) != len(indices_array):
            raise ValueError("The final CSR offset must equal the indices size.")

        self.crs = crs
        self._edges = edges_array
        self._offsets = offsets_array
        self._indices = indices_array
        self._wkb_blob = bytes(wkb_blob)
        self._wkb_offsets = None
        self._coordinates = None
        self._geometry_offsets = None
        self._ptr = capi.lib.fastnsa_network_create(
            self._edges.ctypes.data_as(capi.P_DOUBLE),
            len(self._edges),
            self._offsets.ctypes.data_as(capi.P_INT),
            self._indices.ctypes.data_as(capi.P_INT),
            len(self._offsets) - 1,
        )
        if not self._ptr:
            raise backend_error(capi.last_error())

        if use_coordinates:
            coordinates_array = np.ascontiguousarray(coordinates, dtype=np.float64)
            geometry_offsets_array = np.ascontiguousarray(
                geometry_offsets, dtype=np.int32
            )
            if coordinates_array.ndim != 2 or coordinates_array.shape[1] != 2:
                self.close()
                raise ValueError("coordinates must have shape (n_vertices, 2).")
            if geometry_offsets_array.ndim != 1:
                self.close()
                raise ValueError("geometry_offsets must be a 1D array.")
            if len(geometry_offsets_array) != len(edges_array) + 1:
                self.close()
                raise ValueError("geometry_offsets must contain n_edges + 1 entries.")
            if (
                int(geometry_offsets_array[0]) != 0
                or int(geometry_offsets_array[-1]) != len(coordinates_array)
            ):
                self.close()
                raise ValueError("geometry_offsets must span the packed coordinate array.")
            self._coordinates = coordinates_array
            self._geometry_offsets = geometry_offsets_array
            status = capi.lib.fastnsa_network_set_polylines(
                self._ptr,
                self._coordinates.ctypes.data_as(capi.P_DOUBLE),
                self._geometry_offsets.ctypes.data_as(capi.P_INT),
                len(self._edges),
            )
        else:
            if wkb_offsets is None:
                self.close()
                raise ValueError("Provide packed coordinates or WKB edge geometries.")
            wkb_offsets_array = np.ascontiguousarray(wkb_offsets, dtype=np.int32)
            if wkb_offsets_array.ndim != 1:
                self.close()
                raise ValueError("wkb_offsets must be a 1D array.")
            if len(wkb_offsets_array) != len(edges_array) + 1:
                self.close()
                raise ValueError("wkb_offsets must contain n_edges + 1 entries.")
            if int(wkb_offsets_array[-1]) != len(wkb_blob):
                self.close()
                raise ValueError("The final WKB offset must equal the WKB buffer size.")
            self._wkb_offsets = wkb_offsets_array
            wkb_array = np.frombuffer(self._wkb_blob, dtype=np.uint8)
            status = capi.lib.fastnsa_network_set_wkbs_checked(
                self._ptr,
                wkb_array.ctypes.data_as(capi.P_UINT8),
                self._wkb_offsets.ctypes.data_as(capi.P_INT),
                len(self._edges),
            )
        if status != 0:
            self.close()
            raise backend_error(capi.last_error())

    @property
    def num_edges(self) -> int:
        """Number of edges in native edge-ID order."""
        return len(self._edges)

    @property
    def num_nodes(self) -> int:
        """Number of nodes represented by the CSR topology."""
        return len(self._offsets) - 1

    @property
    def total_length(self) -> float:
        """Sum of edge lengths in the input network distance unit."""
        return float(np.sum(self._edges[:, 2]))

    def close(self):
        """Release the native network allocation.

        Calling this method more than once is safe. Any point set or analysis
        still using the network becomes invalid after it is closed.
        """
        if getattr(self, "_ptr", None):
            capi.lib.fastnsa_network_destroy(self._ptr)
            self._ptr = None

    def _require_open(self):
        if not getattr(self, "_ptr", None):
            raise RuntimeError("The native network has already been closed.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


class NativePointSet:
    """Own event locations mapped to a :class:`NativeNetwork`.

    Events are represented by a native edge ID and an along-edge offset, with
    an optional numeric time value. Construct instances with :meth:`map_xy` or
    :meth:`from_mapped`; the direct constructor is intended for internal use.

    The point set keeps a Python reference to its network but does not own the
    network allocation. Close the point set before closing its network.
    """

    _METHODS = {
        "R-tree": capi.FASTNSA_PLNNS_RTREE,
        "comp_a": capi.FASTNSA_PLNNS_COMP_APPROX,
        "comp_e": capi.FASTNSA_PLNNS_COMP_EXACT,
    }

    def __init__(self, network: NativeNetwork, pointer, n_points: int):
        if not pointer:
            raise RuntimeError(capi.last_error())
        self.network = network
        self._ptr = pointer
        self.n_points = int(n_points)

    @classmethod
    def map_xy(
        cls,
        network: NativeNetwork,
        x,
        y,
        *,
        mapping_method: str = "R-tree",
        comp_eps: float = 10.0,
        n_threads: int = 0,
        times=None,
    ) -> "NativePointSet":
        """Project planar x/y coordinates onto the nearest network edges.

        Parameters
        ----------
        network : NativeNetwork
            Open network whose geometry and topology define the mapping.
        x, y : array-like
            Equal-length one-dimensional coordinate arrays in the network CRS.
        mapping_method : {"R-tree", "comp_a", "comp_e"}, default "R-tree"
            ``R-tree`` is the general indexed method. ``comp_e`` performs exact
            comparison and ``comp_a`` uses the ``comp_eps`` approximation.
        comp_eps : float, default 10.0
            Approximation tolerance used by ``comp_a``.
        n_threads : int, default 0
            Requested mapping worker count. Zero uses FastNSA's bounded
            automatic policy. Negative values are invalid.
        times : array-like, optional
            Finite numeric values with one entry per event. FastNSA does not
            infer a unit or origin for this low-level input.

        Returns
        -------
        NativePointSet
            Owning mapped point set associated with ``network``.

        Raises
        ------
        ValueError
            If coordinate or time arrays are empty, mismatched, or non-finite.
        FastNSAValidationError
            If ``mapping_method`` is unsupported.
        RuntimeError
            If the network is closed or the native mapper fails.

        References
        ----------
        Hongwei Ye et al. "An Efficient and Accurate Grid Compression Solution
        for Point-to-Line Nearest Neighbor Search." KDD 2026.
        https://doi.org/10.1145/3770855.3817711
        """
        if not getattr(network, "_ptr", None):
            raise RuntimeError("Cannot map points using a closed network.")
        x_array = np.ascontiguousarray(x, dtype=np.float64)
        y_array = np.ascontiguousarray(y, dtype=np.float64)
        if x_array.ndim != 1 or y_array.ndim != 1 or len(x_array) != len(y_array):
            raise ValueError("x and y must be 1D arrays with equal length.")
        if len(x_array) == 0:
            raise ValueError("At least one point coordinate is required.")
        if not np.isfinite(x_array).all() or not np.isfinite(y_array).all():
            raise ValueError("Point coordinates contain NaN or infinite values.")
        mapping_method = enum_value(mapping_method)
        if mapping_method not in cls._METHODS:
            raise FastNSAValidationError(f"Unsupported mapping method: {mapping_method!r}.")
        thread_count = int(n_threads)
        if thread_count < 0:
            raise ValueError("n_threads must be zero (automatic) or a positive integer.")

        pointer = capi.lib.fastnsa_map_points(
            network._ptr,
            x_array.ctypes.data_as(capi.P_DOUBLE),
            y_array.ctypes.data_as(capi.P_DOUBLE),
            len(x_array),
            cls._METHODS[mapping_method],
            float(comp_eps),
            thread_count,
        )
        pointset = cls(network, pointer, len(x_array))
        if times is not None:
            try:
                pointset.set_times(times)
            except Exception:
                pointset.close()
                raise
        return pointset

    @classmethod
    def from_mapped(
        cls,
        network: NativeNetwork,
        edge_ids,
        offsets,
        *,
        times=None,
    ) -> "NativePointSet":
        """Create events from precomputed edge IDs and offsets.

        Use this method to avoid repeating nearest-edge mapping. The values are
        meaningful only for the exact network edge order used to produce them.

        Parameters
        ----------
        network : NativeNetwork
            Open network associated with the precomputed mapping.
        edge_ids : array-like
            Finite integer edge IDs in signed int32 range.
        offsets : array-like
            Finite along-edge offsets with one value per edge ID, measured from
            the stored start of each edge in network distance units.
        times : array-like, optional
            Finite numeric event times with matching length.

        Returns
        -------
        NativePointSet
            Owning native point set associated with ``network``.

        Notes
        -----
        Applications should compare a stable network signature before reusing
        mappings. Reprojection, filtering, sorting, or editing edges can make
        stored IDs and offsets invalid.
        Endpoint round-off within ``max(1e-9, 1e-12 * edge_length)`` network
        units is clamped to the endpoint. Larger errors are rejected.
        """
        if not getattr(network, "_ptr", None):
            raise RuntimeError("Cannot create points on a closed network.")
        raw_edge_ids = np.asarray(edge_ids)
        if not np.issubdtype(raw_edge_ids.dtype, np.number):
            raise ValueError("edge_ids must contain numeric integer values.")
        if not np.isfinite(raw_edge_ids).all() or not np.equal(
            raw_edge_ids,
            np.floor(raw_edge_ids),
        ).all():
            raise ValueError("edge_ids must contain finite integer values.")
        int32_info = np.iinfo(np.int32)
        if np.any(raw_edge_ids < int32_info.min) or np.any(
            raw_edge_ids > int32_info.max
        ):
            raise ValueError("edge_ids contains a value outside the int32 range.")
        edge_ids_array = np.ascontiguousarray(raw_edge_ids, dtype=np.int32)
        offsets_array = np.ascontiguousarray(offsets, dtype=np.float64)
        if edge_ids_array.ndim != 1 or offsets_array.ndim != 1:
            raise ValueError("edge_ids and offsets must be 1D arrays.")
        if len(edge_ids_array) != len(offsets_array):
            raise ValueError("edge_ids and offsets must have equal length.")
        if len(edge_ids_array) == 0:
            raise ValueError("At least one mapped point is required.")
        if not np.isfinite(offsets_array).all():
            raise ValueError("offsets contains NaN or infinite values.")

        pointer = capi.lib.fastnsa_pointset_create(network._ptr)
        pointset = cls(network, pointer, len(edge_ids_array))
        status = capi.lib.fastnsa_pointset_set(
            pointset._ptr,
            edge_ids_array.ctypes.data_as(capi.P_INT),
            offsets_array.ctypes.data_as(capi.P_DOUBLE),
            len(edge_ids_array),
        )
        if status != 0:
            pointset.close()
            raise backend_error(capi.last_error())
        if times is not None:
            try:
                pointset.set_times(times)
            except Exception:
                pointset.close()
                raise
        return pointset

    def set_times(self, times):
        """Replace event times with a finite one-dimensional numeric array.

        The array length must equal :attr:`n_points`. Time units are not
        converted or interpreted by the native runtime.
        """
        self._require_open()
        times_array = np.ascontiguousarray(times, dtype=np.float64)
        if times_array.ndim != 1 or len(times_array) != self.n_points:
            raise ValueError("times must be a 1D array matching the point count.")
        if not np.isfinite(times_array).all():
            raise ValueError("times contains NaN or infinite values.")
        status = capi.lib.fastnsa_pointset_set_times(
            self._ptr,
            times_array.ctypes.data_as(capi.P_DOUBLE),
            self.n_points,
        )
        if status != 0:
            raise backend_error(capi.last_error())

    def _copy_array(self, function, dtype, label: str) -> np.ndarray:
        self._require_open()
        normalized_dtype = np.dtype(dtype)
        values = np.empty(self.n_points, dtype=normalized_dtype)
        pointer_type = (
            capi.P_INT if normalized_dtype == np.dtype(np.int32) else capi.P_DOUBLE
        )
        status = function(self._ptr, values.ctypes.data_as(pointer_type), self.n_points)
        if status != 0:
            raise RuntimeError(f"Could not copy {label}; backend returned {status}.")
        return values

    def edge_ids_array(self) -> np.ndarray:
        """Return a new int32 array containing mapped edge IDs."""
        return self._copy_array(capi.lib.fastnsa_pointset_copy_edge_ids, np.int32, "edge IDs")

    def offsets_array(self) -> np.ndarray:
        """Return a new float64 array of offsets from stored edge starts."""
        return self._copy_array(capi.lib.fastnsa_pointset_copy_offsets, np.float64, "offsets")

    def times_array(self) -> np.ndarray:
        """Return a new float64 array containing native event times."""
        return self._copy_array(capi.lib.fastnsa_pointset_copy_times, np.float64, "times")

    def close(self):
        """Release the native point-set allocation; repeated calls are safe."""
        if getattr(self, "_ptr", None):
            capi.lib.fastnsa_pointset_destroy(self._ptr)
            self._ptr = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def _require_open(self):
        if not getattr(self, "_ptr", None):
            raise RuntimeError("The native point set has already been closed.")
