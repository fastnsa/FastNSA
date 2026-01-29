import numpy as np
import networkx as nx
from fastnsa.bindings import capi
from shapely.geometry import LineString
import osmnx as ox


class Network: 
    def __init__(self, graph:nx.MultiGraph):
        self.crs = graph.graph.get("crs")
        self.graph = graph  # 保持引用，防止被回收
        
        self.fix_direciton(graph)
    
        self.nodes_list = list(graph.nodes())
        self._node_to_idx = {node_id: i for i, node_id in enumerate(self.nodes_list)}

        ir_data = self._prepare_network_ir(graph)
        self._data_cache = ir_data

        self._ptr = capi.lib.fastnsa_network_create(
            ir_data['edges'].ctypes.data_as(capi.P_DOUBLE),
            ir_data['n_edges'],
            ir_data['offsets'].ctypes.data_as(capi.P_INT),
            ir_data['indices'].ctypes.data_as(capi.P_INT),
            ir_data['n_nodes']     
        )
        self._set_wkbs(graph)
        
    
    def __del__(self):
        if hasattr(self, '_ptr') and self._ptr:
            capi.lib.fastnsa_network_destroy(self._ptr)
            self._ptr = None
            
        
    def _set_wkbs(self, graph):    
        edges_gdf = ox.graph_to_gdfs(graph, nodes=False)
        wkb_list = edges_gdf.geometry.apply(
            lambda g: g.wkb if g is not None else None
        )
        wkb_bytes = []
        offsets = [0]
        for b in wkb_list:
            if b is None:
                b = b''

            wkb_bytes.append(b)
            offsets.append(offsets[-1] + len(b))

        wkb_blob = b''.join(wkb_bytes)

        wkb_array = np.frombuffer(wkb_blob, dtype=np.uint8)
        offsets = np.array(offsets, dtype=np.int32)
        
        capi.lib.fastnsa_network_set_wkbs(
            self._ptr,
            wkb_array.ctypes.data_as(capi.P_UINT8),
            offsets.ctypes.data_as(capi.P_INT),
            len(offsets) - 1
        )
    
    @staticmethod
    def fix_direciton(g:nx.MultiGraph):
        x_dic = {node_id: data['x'] for node_id, data in g.nodes(data=True)}

        for u, v, data in g.edges(data=True):
            geom = data.get('geometry')

            x_coords = geom.xy[0]
            if abs(x_coords[0] - x_dic.get(u, x_coords[0])) > 1e-5:
                reversed_coords = list(geom.coords)[::-1]
                data['geometry'] = LineString(reversed_coords)


    @classmethod
    def from_osm(cls, bbox=None, place=None, network_type="drive", 
                 crs=None, consolidate_tolerance=20):
        import osmnx as ox
        
        if place:
            g = ox.graph_from_place(place, simplify=True, network_type=network_type)
        else:
            g = ox.graph_from_bbox(*bbox, simplify=True, network_type=network_type)

        g = ox.project_graph(g, to_crs=crs)
        g = ox.consolidate_intersections(
            g, tolerance=consolidate_tolerance, rebuild_graph=True
        )

        g = g.to_undirected()

        return cls(g)
    

    @property
    def num_edges(self) -> int:
        return self._data_cache['n_edges']
    
    @property
    def num_nodes(self) -> int:
        return self._data_cache['n_nodes']

    @property
    def total_length(self) -> float:
        return np.sum(self._data_cache['edges'][:, 2])

    
    def _prepare_network_ir(self, graph):
            # ---- edge indexing (critical) ----
        edge_list = []
        self.edge_id_map = {}   # (u, v, k) -> eid

        for eid, (u, v, k, data) in enumerate(graph.edges(keys=True, data=True)):
            # geom: LineString = data["geometry"]
            length = data["geometry"].length

            u_idx = self._node_to_idx[u]
            v_idx = self._node_to_idx[v]
            edge_list.append((u_idx, v_idx, length))
            self.edge_id_map[(u, v, k)] = eid

        # edges = sorted(edges, key=lambda x: (x[0], x[1]))  
        n_edges = len(edge_list)
        edges_arr = np.array(edge_list, dtype=np.float64) # [u_idx, v_idx, len]

        n_nodes = len(self.nodes_list)
        adj = [[] for _ in range(n_nodes)]
        for eid, (u_idx, v_idx, _) in enumerate(edge_list):
            adj[int(u_idx)].append(eid)
            adj[int(v_idx)].append(eid)

        offsets = np.zeros(n_nodes + 1, dtype=np.int32)
        indices = np.zeros(n_edges * 2, dtype=np.int32)

        curr = 0
        for i in range(n_nodes):
            offsets[i] = curr
            for eid in adj[i]:
                indices[curr] = eid
                curr += 1
        offsets[n_nodes] = curr

        return {
            'edges': edges_arr,
            'offsets': offsets,
            'indices': indices,
            'n_edges': n_edges,
            'n_nodes': n_nodes
        }

