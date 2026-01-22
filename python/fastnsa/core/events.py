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
            
        # coords = np.ascontiguousarray(
        #     np.vstack([points_gdf.geometry.x, points_gdf.geometry.y]).T, 
        #     dtype=np.float64
        # )

        # TODO： 考虑大数据量时，分批处理以避免内存问题
        x_ptr = np.ascontiguousarray(points_gdf.geometry.x.values, dtype=np.float64)
        y_ptr = np.ascontiguousarray(points_gdf.geometry.y.values, dtype=np.float64)

        # 3. 调用 C++ 执行 Mapping (假设 C++ 侧实现了空间索引加速)
        # self._perform_cpp_mapping(x_ptr, y_ptr)
        self._perform_python_mapping_fallback(points_gdf)
        
        
        

    def __del__(self):
        if hasattr(self, '_ptr') and self._ptr:
            capi.lib.fastnsa_pointset_destroy(self._ptr)


    def _perform_cpp_mapping(self, x_ptr, y_ptr):
        # 注意：这里直接传入 network 的 C++ 句柄
        self._ptr = capi.lib.fastnsa_map_points( # 直接返回PointSet指针
            self.network._ptr, # todo graph geometry对接
            x_ptr.ctypes.data_as(capi.P_DOUBLE),
            y_ptr.ctypes.data_as(capi.P_DOUBLE),
            self.n_points
        )


    def _perform_python_mapping_fallback(self, points_gdf):
        import osmnx as ox

        # 1. 预分配 C++ 映射结果的内存空间
        self.edge_ids = np.zeros(self.n_points, dtype=np.int32)
        self.offsets = np.zeros(self.n_points, dtype=np.float64)

        nearest = ox.nearest_edges(
            self.network.graph, # todo 对接
            points_gdf.geometry.x.values,
            points_gdf.geometry.y.values,
        )

        # for i in range(len(points_gdf)):
        #     u, v, _ = nearest[i]
        #     if u > v:
        #         u, v = v, u

        #     edge_geom = self.network.graph[u][v][0]['geometry']
        #     proj_dist = edge_geom.project(points_gdf.geometry.iloc[i])
        #     self.edge_ids[i] = self.network._node_to_idx[(u, v)]
        #     self.offsets[i] = proj_dist

        for i, (u, v, k) in enumerate(nearest):
            # 这里的关键：必须使用 Network 类在构建时生成的 edge_id_map
            eid = self.network.edge_id_map.get((u, v, k))
            if eid is None:
                # 针对无向图的兜底判断
                eid = self.network.edge_id_map.get((v, u, k))
            
            self.edge_ids[i] = eid
            
            # 计算偏移 (这里还是调用了 shapely，在大数据量下会很慢)
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
        


    @property
    def ir_data(self):
        """返回符合你之前设计的 points_ir 格式: [edge_id, offset]"""
        return self.edge_ids, self.offsets
