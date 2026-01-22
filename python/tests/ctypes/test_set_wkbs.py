import numpy as np
import ctypes
from fastnsa.bindings.capi import lib

def test_network_set_wkbs(simple_network):
    # 假的 WKB 数据（哪怕是空，也要合法）
    wkb_data = np.ascontiguousarray([1, 2, 3, 4], dtype=np.uint8)
    offsets = np.ascontiguousarray([0, 4], dtype=np.int32)

    lib.fastnsa_network_set_wkbs(
        simple_network,
        wkb_data.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        1
    )
