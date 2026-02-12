from ._lib import lib
import ctypes


FASTNSA_K_RQS = 1
FASTNSA_K_SPS = 2
FASTNSA_K_EAR = 3
FASTNSA_K_MTS = 4
FASTNSA_K_EARTH = 5
FASTNSA_K_NS = 6
FASTNSA_K_STARS = 7

FASTNSA_PLNNS_RTREE = 0
FASTNSA_PLNNS_COMP_APPROX = 1
FASTNSA_PLNNS_COMP_EXACT = 3

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

# 参数：network_ptr(void*), wkb_data(uint8*), offsets(int*), n_edges(int)
lib.fastnsa_network_set_wkbs.argtypes = [VOID_P, P_UINT8, P_INT, ctypes.c_int]
lib.fastnsa_network_set_wkbs.restype = None


# --- PointSet 相关绑定 ---
lib.fastnsa_pointset_create.argtypes = [VOID_P] # 传入 network_ptr
lib.fastnsa_pointset_restype = VOID_P

lib.fastnsa_pointset_destroy.argtypes = [VOID_P]
lib.fastnsa_pointset_destroy.restype = None

# 将数据从 NumPy 填充进 C++ PointSet 
# 参数：points_ptr(void*), edge_ids(int*), offsets(double*), n_points(int)
lib.fastnsa_pointset_set.argtypes = [VOID_P, P_INT, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_set.restype = ctypes.c_int # 返回状态码
lib.fastnsa_pointset_set_st.argtypes = [VOID_P, P_INT, P_DOUBLE, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_set_st.restype = ctypes.c_int


lib.fastnsa_map_points.argtypes = [
    ctypes.c_void_p,  # net_ptr - pointer to the network object
    P_DOUBLE,
    P_DOUBLE,  
    ctypes.c_int,  # n_points(int) - number of points to map
    ctypes.c_int,  # method(int) - mapping method
    ctypes.c_double, # comp_eps(double) - computation epsilon for approximate methods
    ctypes.c_int   # n_threads(int) - number of threads to use
]
lib.fastnsa_map_points.restype = VOID_P # 返回一个 PointSet 对象的指针

lib.fastnsa_pointset_dump_text.restype = ctypes.c_int # 状态码：0成功， -1 空指针， -2文件打开失败
lib.fastnsa_pointset_dump_text.argtypes = [VOID_P, ctypes.c_char_p] # void point_set_ptr， const char* filenme



# 参数：net_ptr, event_ptr, r_values, n_r, k_obs,k_sim
lib.fastnsa_compute_k_function.argtypes = [
    VOID_P,   # network_ptr
    VOID_P,   # point_set_ptr
    P_DOUBLE,      # r_values
    ctypes.c_int, # n_r
    ctypes.c_int, # method
    ctypes.c_int, # n_simulations
    P_DOUBLE, # out k_obs
    P_DOUBLE, # out k_sim
    ctypes.c_double # tolerance for approximate methods
]
lib.fastnsa_compute_k_function.restype = None

lib.fastnsa_compute_k_function_earth.argtypes = [
    VOID_P,   # network_ptr
    VOID_P,   # point_set_ptr
    P_DOUBLE, # s_values
    ctypes.c_int, # n_s
    P_DOUBLE, # t_values
    ctypes.c_int, # n_t
    ctypes.c_int, # method
    ctypes.c_int, # n_simulations
    P_DOUBLE, # out k_obs
    P_DOUBLE, # out k_sim
]
lib.fastnsa_compute_k_function_earth.restype = None


