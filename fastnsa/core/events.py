import ctypes
from shapely.geometry import LineString
import numpy as np
from fastnsa.bindings import capi
from typing import Literal, Optional
import os

class PointEvents:
    def __init__(self, points_gdf, target_network, 
                 mapping_method:Literal["R-tree", "comp_a", "comp_e"] = "R-tree",
                 comp_eps: float = 10.0,
                 n_threads: int = 1,
                 time_col: Optional[str] = None,
                 times: Optional[np.ndarray] = None
                 ):
        """
        points_gdf: GeoDataFrame with point geometries to be mapped
        target_network: Network object representing the target network
        mapping_method: Method to use for mapping points onto the network
            "R-tree": Use R-tree based nearest edge search
            "comp_a": Use Grid Compression approximate mapping method
            "comp_e": Use Grid Compression exact mapping method 
        comp_eps: Tolerance for approximate methods
        n_threads: Number of threads to use for mapping
        """
        self.network = target_network  # keep a reference to avoid premature GC
        self.n_points = len(points_gdf)
                
        # make sure CRS matches
        if points_gdf.crs != self.network.crs:
            points_gdf = points_gdf.to_crs(self.network.crs)
            
        x_ptr = np.ascontiguousarray(points_gdf.geometry.x.values, dtype=np.float64)
        y_ptr = np.ascontiguousarray(points_gdf.geometry.y.values, dtype=np.float64)
        
        method_map = {
            "R-tree": capi.FASTNSA_PLNNS_RTREE,
            "comp_a": capi.FASTNSA_PLNNS_COMP_APPROX,
            "comp_e": capi.FASTNSA_PLNNS_COMP_EXACT,
        }
        method_flag = method_map.get(mapping_method, capi.FASTNSA_PLNNS_RTREE)

        self._perform_cpp_mapping(x_ptr, y_ptr, method_flag, comp_eps,n_threads)
        self._update_times(time_col=time_col, times=times, points_gdf=points_gdf)
        
        
    def __del__(self):
        if hasattr(self, '_ptr') and self._ptr:
            capi.lib.fastnsa_pointset_destroy(self._ptr)
            
    def dump_to_file(self, filename):
        path = os.fspath(filename)
        c_path = ctypes.c_char_p(filename.encode("utf-8"))        
        capi.lib.fastnsa_pointset_dump_text(self._ptr, c_path)


    def _perform_cpp_mapping(self, x_ptr, y_ptr, method_flag, comp_eps=10, n_threads=1):
        self._ptr = capi.lib.fastnsa_map_points(
            self.network._ptr, 
            x_ptr.ctypes.data_as(capi.P_DOUBLE),
            y_ptr.ctypes.data_as(capi.P_DOUBLE),
            self.n_points,
            method_flag,
            comp_eps,
            n_threads
        )


    def _perform_python_mapping_fallback(self, points_gdf):
        import osmnx as ox
        
        self.edge_ids = np.zeros(self.n_points, dtype=np.int32)
        self.offsets = np.zeros(self.n_points, dtype=np.float64)

        nearest = ox.nearest_edges(
            self.network.graph,
            points_gdf.geometry.x.values,
            points_gdf.geometry.y.values,
        )

        for i, (u, v, k) in enumerate(nearest):
            eid = self.network.edge_id_map.get((u, v, k))
            if eid is None:
                eid = self.network.edge_id_map.get((v, u, k))
            
            self.edge_ids[i] = eid
            
            geom:LineString= self.network.graph[u][v][k]['geometry']
            self.offsets[i] = geom.project(points_gdf.geometry.iloc[i])
        
        self._ptr = capi.lib.fastnsa_pointset_create(self.network._ptr)
        self._update_data(self.edge_ids, self.offsets)
    
    def _update_data(self, edge_ids, offsets):
            capi.lib.fastnsa_pointset_set(
                self._ptr,
                edge_ids.ctypes.data_as(capi.P_INT),
                offsets.ctypes.data_as(capi.P_DOUBLE),
                len(edge_ids)
            )

    def _update_times(self, time_col: Optional[str], times: Optional[np.ndarray], points_gdf):
        if times is not None:
            time_arr = np.ascontiguousarray(times, dtype=np.float64)
        elif time_col is not None:
            if time_col not in points_gdf.columns:
                raise ValueError(f"`time_col` {time_col!r} not found in points data.")
            time_arr = np.ascontiguousarray(points_gdf[time_col].to_numpy(), dtype=np.float64)
        else:
            return

        if time_arr.ndim != 1 or len(time_arr) != self.n_points:
            raise ValueError("`times` must be a 1D array with length equal to number of points.")

        ret = capi.lib.fastnsa_pointset_set_times(
            self._ptr,
            time_arr.ctypes.data_as(capi.P_DOUBLE),
            self.n_points
        )
        if ret != 0:
            raise RuntimeError(f"fastnsa_pointset_set_times failed with code {ret}.")
        

