import numpy as np
import ctypes
from fastnsa.bindings.capi import lib

def test_map_points(graph_network):
    xs = np.ascontiguousarray([5.0], dtype=np.float64)
    ys = np.ascontiguousarray([0.0], dtype=np.float64)

    # pts = lib.fastnsa_map_points(
    #     graph_network,
    #     xs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    #     ys.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    #     1
    # )

    # assert pts is not None and pts != 0

    # lib.fastnsa_pointset_destroy(pts)
