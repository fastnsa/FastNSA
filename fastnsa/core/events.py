from shapely.geometry import LineString
import numpy as np
from fastnsa.bindings import capi

class PointEvents:
    def __init__(self, points_gdf, target_network):
        """
        points_gdf: 包含几何点的 GeoDataFrame
        target_network: 已经构建好的 Network 对象
        """
        self.network = target_network  # 保持引用，防止 Network 先被回收
        self.n_points = len(points_gdf)
                
        # 确保点数据与路网 CRS 一致
        if points_gdf.crs != self.network.crs:
            points_gdf = points_gdf.to_crs(self.network.crs)
            
        x_ptr = np.ascontiguousarray(points_gdf.geometry.x.values, dtype=np.float64)
        y_ptr = np.ascontiguousarray(points_gdf.geometry.y.values, dtype=np.float64)

        self._perform_cpp_mapping(x_ptr, y_ptr)
        
        

    def __del__(self):
        if hasattr(self, '_ptr') and self._ptr:
            capi.lib.fastnsa_pointset_destroy(self._ptr)


    def _perform_cpp_mapping(self, x_ptr, y_ptr):
        self._ptr = capi.lib.fastnsa_map_points(
            self.network._ptr, 
            x_ptr.ctypes.data_as(capi.P_DOUBLE),
            y_ptr.ctypes.data_as(capi.P_DOUBLE),
            self.n_points
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
        


