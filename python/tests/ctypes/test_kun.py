import numpy as np
import ctypes

from fastnsa.bindings.capi import (
    lib,
    FASTNSA_K_NS,
)

def test_compute_k_function(simple_network, simple_pointset):
    r_values = np.ascontiguousarray([1.0, 2.0, 5.0], dtype=np.float64)
    n_r = len(r_values)

    k_obs = np.zeros(n_r, dtype=np.float64)
    k_sim = np.zeros(n_r, dtype=np.float64)

    lib.fastnsa_compute_k_function(
        simple_network,
        simple_pointset,
        r_values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        n_r,
        6, #ns method
        0,  # n_simulations
        k_obs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        k_sim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )

    # 不验证数值，只验证没崩
    print("K-function observed:", k_obs)
    assert np.all(np.isfinite(k_obs))
