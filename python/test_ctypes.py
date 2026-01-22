import ctypes
import numpy as np
from ctypes import c_void_p, c_int, c_double

from fastnsa.bindings.capi import lib



# ---- build network ----
edges = np.array([
    [0, 1, 10.0]
], dtype=np.float64)

csr_offsets = np.array([0, 1, 2], dtype=np.int32)
csr_indices = np.array([0, 0], dtype=np.int32)

net = lib.fastnsa_network_create(
    edges.ctypes.data_as(ctypes.POINTER(c_double)),
    1,
    csr_offsets.ctypes.data_as(ctypes.POINTER(c_int)),
    csr_indices.ctypes.data_as(ctypes.POINTER(c_int)),
    2
)

print("Network created:", net)

# ---- pointset ----
pts = lib.fastnsa_pointset_create(net)

edge_ids = np.array([0], dtype=np.int32)
offsets = np.array([5.0], dtype=np.float64)

lib.fastnsa_pointset_set(
    pts,
    edge_ids.ctypes.data_as(ctypes.POINTER(c_int)),
    offsets.ctypes.data_as(ctypes.POINTER(c_double)),
    1
)
print("PointSet created and set:", pts)
# ---- K-function ----
r_vals = np.array([1.0, 2.0, 5.0], dtype=np.float64)
out_k = np.zeros(3, dtype=np.float64)

lib.fastnsa_compute_k_function(
    net, pts,
    r_vals.ctypes.data_as(ctypes.POINTER(c_double)),
    len(r_vals),
    6,  # method
    0,  # simulations
    out_k.ctypes.data_as(ctypes.POINTER(c_double)),
    None
)
print("K-function computed.")

print("K:", out_k)
print("Python ctypes test OK")

