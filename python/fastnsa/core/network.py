"""Graph-compatible network construction backed by the native runtime."""

from __future__ import annotations

from typing import Any

import numpy as np

from fastnsa.runtime import NativeNetwork


class Network:
    """Build a native FastNSA network from a NetworkX-style multigraph.

    The input graph is copied. Edge geometries are aligned with graph direction,
    graph node IDs are mapped to contiguous native indices, and graph edge order
    defines the edge IDs exposed to mapped events.

    Parameters
    ----------
    graph : networkx.MultiGraph-like
        Graph with CRS metadata in ``graph.graph["crs"]``. Every edge must have
        a LineString geometry, or both endpoint nodes must provide finite ``x``
        and ``y`` coordinates from which a straight geometry can be created.

    Notes
    -----
    Network distances use geometry lengths in the graph CRS. Use a suitable
    projected CRS rather than longitude/latitude degrees. This object owns a
    native network and should be closed or used as a context manager.
    """

    def __init__(self, graph):
        if graph is None:
            raise ValueError("`graph` must be a NetworkX-style graph.")
        graph = graph.copy()
        self.crs = graph.graph.get("crs")
        self.graph = graph
        self.fix_direction(graph)

        self.nodes_list = list(graph.nodes())
        self._node_to_idx = {
            node_id: index for index, node_id in enumerate(self.nodes_list)
        }
        ir_data = self._prepare_network_ir(graph)
        geometries = [
            data["geometry"] for _, _, _, data in graph.edges(keys=True, data=True)
        ]
        coordinates, geometry_offsets = self._pack_coordinates(geometries)
        self._initialize_backend(
            ir_data,
            coordinates=coordinates,
            geometry_offsets=geometry_offsets,
        )

    @classmethod
    def from_prepared_arrays(
        cls,
        *,
        crs: Any,
        edges: np.ndarray,
        offsets: np.ndarray,
        indices: np.ndarray,
        wkb_blob: bytes = b"",
        wkb_offsets: np.ndarray | None = None,
        coordinates: np.ndarray | None = None,
        geometry_offsets: np.ndarray | None = None,
        graph: Any = None,
    ) -> "Network":
        """Create a compatibility wrapper from host-prepared native arrays.

        Parameters follow :class:`fastnsa.NativeNetwork`. ``graph`` is optional;
        methods requiring original graph attributes, such as
        :meth:`save_graph_gpkg`, are unavailable when it is omitted. This path
        is intended for QGIS and other hosts that already own network features.
        """
        instance = cls.__new__(cls)
        instance.crs = crs
        instance.graph = graph
        instance.nodes_list = list(range(int(len(offsets) - 1)))
        instance._node_to_idx = {
            node_id: node_id for node_id in instance.nodes_list
        }
        instance.edge_id_map = {}
        ir_data = {
            "edges": np.ascontiguousarray(edges, dtype=np.float64),
            "offsets": np.ascontiguousarray(offsets, dtype=np.int32),
            "indices": np.ascontiguousarray(indices, dtype=np.int32),
            "n_edges": int(len(edges)),
            "n_nodes": int(len(offsets) - 1),
        }
        instance._initialize_backend(
            ir_data,
            wkb_blob=wkb_blob,
            wkb_offsets=wkb_offsets,
            coordinates=coordinates,
            geometry_offsets=geometry_offsets,
        )
        return instance

    def _initialize_backend(
        self,
        ir_data,
        wkb_blob=b"",
        wkb_offsets=None,
        *,
        coordinates=None,
        geometry_offsets=None,
    ):
        self._data_cache = ir_data
        self._wkb_blob = bytes(wkb_blob)
        self._wkb_offsets = (
            None if wkb_offsets is None
            else np.ascontiguousarray(wkb_offsets, dtype=np.int32)
        )
        self._coordinates = (
            None if coordinates is None
            else np.ascontiguousarray(coordinates, dtype=np.float64)
        )
        self._geometry_offsets = (
            None if geometry_offsets is None
            else np.ascontiguousarray(geometry_offsets, dtype=np.int32)
        )
        self._native = NativeNetwork(
            edges=ir_data["edges"],
            offsets=ir_data["offsets"],
            indices=ir_data["indices"],
            wkb_blob=self._wkb_blob,
            wkb_offsets=self._wkb_offsets,
            coordinates=self._coordinates,
            geometry_offsets=self._geometry_offsets,
            crs=self.crs,
        )

    @property
    def _ptr(self):
        return self._native._ptr

    def close(self):
        """Release the underlying native network allocation."""
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

    @staticmethod
    def _pack_wkbs(wkb_list):
        wkb_bytes = [value if value is not None else b"" for value in wkb_list]
        lengths = np.fromiter((len(value) for value in wkb_bytes), dtype=np.int64)
        offsets = np.empty(len(wkb_bytes) + 1, dtype=np.int64)
        offsets[0] = 0
        np.cumsum(lengths, out=offsets[1:])
        if offsets[-1] > np.iinfo(np.int32).max:
            raise ValueError("Packed network WKB exceeds the native int32 offset limit.")
        return b"".join(wkb_bytes), offsets.astype(np.int32)

    @staticmethod
    def _pack_coordinates(geometries):
        """Pack LineString coordinates for the dependency-free native index."""
        coordinate_rows = []
        offsets = [0]
        for geometry in geometries:
            if geometry is None or not hasattr(geometry, "coords"):
                raise ValueError("Every network edge must provide a LineString geometry.")
            rows = [(float(row[0]), float(row[1])) for row in geometry.coords]
            if len(rows) < 2:
                raise ValueError("Every network edge must contain at least two vertices.")
            coordinate_rows.extend(rows)
            offsets.append(len(coordinate_rows))
        return (
            np.ascontiguousarray(coordinate_rows, dtype=np.float64),
            np.ascontiguousarray(offsets, dtype=np.int32),
        )

    def save_graph_gpkg(self, filename):
        """Write graph edge geometries to an ``edges`` GeoPackage layer.

        Requires a graph-backed instance and the ``fastnsa[osm]`` extra.

        Parameters
        ----------
        filename : path-like
            Destination GeoPackage path.
        """
        if self.graph is None:
            raise RuntimeError("save_graph_gpkg requires a graph-backed Network.")
        try:
            import geopandas as gpd
            import osmnx as ox
        except ImportError as exc:
            raise ImportError(
                'Saving graph data requires the "osm" extra: pip install "fastnsa[osm]".'
            ) from exc

        edges_gdf = ox.graph_to_gdfs(self.graph, nodes=False)
        edges = gpd.GeoDataFrame(geometry=edges_gdf.geometry, crs=edges_gdf.crs)
        edges.to_file(filename, driver="GPKG", layer="edges")

    @staticmethod
    def fix_direction(graph):
        """Align each edge geometry start with its graph start node in place.

        Missing geometries are constructed as straight LineStrings from finite
        endpoint ``x``/``y`` attributes. This helper is primarily used during
        :class:`Network` construction.
        """
        coordinates = {
            node_id: (data.get("x"), data.get("y"))
            for node_id, data in graph.nodes(data=True)
        }
        line_string = None
        for u, v, data in graph.edges(data=True):
            geometry = data.get("geometry")
            if geometry is None:
                start = coordinates.get(u)
                end = coordinates.get(v)
                coordinate_values = () if start is None or end is None else (*start, *end)
                if (
                    not coordinate_values
                    or any(value is None for value in coordinate_values)
                    or not np.isfinite(np.asarray(coordinate_values, dtype=np.float64)).all()
                ):
                    raise ValueError(
                        "Edges without geometry require finite x/y coordinates on both nodes."
                    )
                if line_string is None:
                    try:
                        from shapely.geometry import LineString
                    except ImportError as exc:
                        raise ImportError(
                            'Graph conversion requires the "geo" extra: '
                            'pip install "fastnsa[geo]".'
                        ) from exc
                    line_string = LineString
                geometry = line_string([start, end])
                data["geometry"] = geometry
            x_coords = geometry.xy[0]
            node_x = coordinates.get(u, (None, None))[0]
            if node_x is not None and abs(x_coords[0] - node_x) > 1e-5:
                if line_string is None:
                    try:
                        from shapely.geometry import LineString
                    except ImportError as exc:
                        raise ImportError(
                            'Graph conversion requires the "geo" extra: '
                            'pip install "fastnsa[geo]".'
                        ) from exc
                    line_string = LineString
                data["geometry"] = line_string(list(geometry.coords)[::-1])

    # Compatibility alias for the misspelled pre-0.2 private helper.
    fix_direciton = fix_direction

    @classmethod
    def from_osm(
        cls,
        bbox=None,
        place=None,
        network_type="drive",
        crs=None,
        consolidate_tolerance=20,
    ):
        """Download, project, and prepare an OpenStreetMap network.

        Parameters
        ----------
        bbox : sequence, optional
            Bounding-box arguments accepted by the installed OSMnx
            ``graph_from_bbox`` version. Ignored when ``place`` is provided.
        place : str or mapping, optional
            OSMnx geocoding query. Either ``place`` or ``bbox`` is required.
        network_type : str, default "drive"
            OSMnx network type, such as ``drive``, ``walk``, or ``bike``.
        crs : CRS-like, optional
            Target projected CRS. If omitted, OSMnx chooses a projected CRS.
        consolidate_tolerance : float, default 20
            Intersection consolidation tolerance in projected network units.

        Returns
        -------
        Network
            Undirected, projected FastNSA network.

        Raises
        ------
        ImportError
            If the ``fastnsa[osm]`` dependencies are unavailable.
        ValueError
            If neither ``place`` nor ``bbox`` is supplied.
        """
        try:
            import osmnx as ox
        except ImportError as exc:
            raise ImportError(
                'OSM network construction requires: pip install "fastnsa[osm]".'
            ) from exc

        if place:
            graph = ox.graph_from_place(place, simplify=True, network_type=network_type)
        else:
            if bbox is None:
                raise ValueError("Either place or bbox must be provided.")
            graph = ox.graph_from_bbox(*bbox, simplify=True, network_type=network_type)
        graph = ox.project_graph(graph, to_crs=crs)
        graph = ox.consolidate_intersections(
            graph, tolerance=consolidate_tolerance, rebuild_graph=True
        )
        return cls(graph.to_undirected())

    @property
    def num_edges(self) -> int:
        """Number of edges in native edge-ID order."""
        return self._native.num_edges

    @property
    def num_nodes(self) -> int:
        """Number of contiguous native nodes."""
        return self._native.num_nodes

    @property
    def total_length(self) -> float:
        """Total edge length in projected network units."""
        return self._native.total_length

    def _prepare_network_ir(self, graph):
        edge_list = []
        self.edge_id_map = {}
        for edge_id, (u, v, key, data) in enumerate(
            graph.edges(keys=True, data=True)
        ):
            geometry = data.get("geometry")
            if geometry is None:
                raise ValueError("Every network edge must contain a geometry.")
            u_index = self._node_to_idx[u]
            v_index = self._node_to_idx[v]
            edge_list.append((u_index, v_index, geometry.length))
            self.edge_id_map[(u, v, key)] = edge_id

        edges = np.asarray(edge_list, dtype=np.float64).reshape((-1, 3))
        adjacency = [[] for _ in self.nodes_list]
        for edge_id, (u_index, v_index, _) in enumerate(edge_list):
            adjacency[int(u_index)].append(edge_id)
            adjacency[int(v_index)].append(edge_id)

        offsets = np.zeros(len(adjacency) + 1, dtype=np.int32)
        indices = np.empty(sum(map(len, adjacency)), dtype=np.int32)
        cursor = 0
        for node_index, edge_ids in enumerate(adjacency):
            offsets[node_index] = cursor
            indices[cursor : cursor + len(edge_ids)] = edge_ids
            cursor += len(edge_ids)
        offsets[-1] = cursor
        return {
            "edges": edges,
            "offsets": offsets,
            "indices": indices,
            "n_edges": len(edges),
            "n_nodes": len(adjacency),
        }
