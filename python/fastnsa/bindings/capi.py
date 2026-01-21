from ._lib import lib
import ctypes



FASTNSA_K_NS = 0
FASTNSA_K_SPS = 1
FASTNSA_K_STARS = 2

P_DOUBLE = ctypes.POINTER(ctypes.c_double)
P_INT = ctypes.POINTER(ctypes.c_int)
VOID_P = ctypes.c_void_p
P_UINT8 = ctypes.POINTER(ctypes.c_uint8)


# --- Network 相关绑定 ---
# 参数：edges(double*), n_edges(int), offsets(int*), indices(int*), n_nodes(int)
lib.fastnsa_network_create.argtypes = [
    P_DOUBLE,
    ctypes.c_int,
    P_INT,
    P_INT,
    ctypes.c_int
]
lib.fastnsa_network_create.restype = VOID_P

lib.fastnsa_network_destroy.argtypes = [VOID_P]
lib.fastnsa_network_destroy.restype = None

lib.fastnsa_network_set_wkbs.argtypes = [VOID_P, P_DOUBLE, P_UINT8,ctypes.c_int]
lib.fastnsa_network_set_wkbs.restype = None


# --- PointSet 相关绑定 ---
lib.fastnsa_pointset_create.argtypes = [VOID_P] # 传入 network_ptr
lib.fastnsa_pointset_restype = VOID_P

lib.fastnsa_pointset_destroy.argtypes = [VOID_P]
lib.fastnsa_pointset_destroy.restype = None

# 将数据从 NumPy 填充进 C++ PointSet 
lib.fastnsa_pointset_set.argtypes = [VOID_P, P_INT, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_set.restype = ctypes.c_int # 返回状态码



# 参数：net_ptr(void*), points(double*), n_points(int), edge_ids(int

    # capi.lib.fastnsa_map_points(
    #     self.network._ptr,
    #     coords.ctypes.data_as(capi.P_DOUBLE),
    #     self.n_points,
    #     self.edge_ids.ctypes.data_as(capi.P_INT),
    #     self.offsets.data_as(capi.P_DOUBLE)
    # )

lib.fastnsa_map_points.argtypes = [
    ctypes.c_void_p,  # net_ptr - pointer to the network object
    P_DOUBLE,
    P_DOUBLE,  
    ctypes.c_int,  # n_points(int) - number of points to map
]
lib.fastnsa_map_points.restype = VOID_P # 返回一个 PointSet 对象的指针



# 参数：net_ptr, , n_points, radii, n_radii, out_results
lib.fastnsa_compute_k_function.argtypes = [
    VOID_P,  
    VOID_P,
    ctypes.c_int,
    P_DOUBLE,      # radii
    ctypes.c_int, # n_radii
    ctypes.c_int, # method
    P_DOUBLE
]
lib.fastnsa_compute_k_function.restype = None



