import numpy as np
import pytest
import ctypes
import sys
sys.path.insert(0, '/home/yeswin/project/fast-nsa/python')

from fastnsa.bindings.capi import lib


@pytest.fixture(scope="function")
def simple_network():
    """
    构造一个最小可用 network:
    2 nodes, 1 edge
    """
    # edge: u, v, length
    edges = np.ascontiguousarray(
        [0.0, 1.0, 10.0],
        dtype=np.float64
    )

    # CSR: 0 -> 1
    offsets = np.ascontiguousarray([0, 1, 1], dtype=np.int32)
    indices = np.ascontiguousarray([1], dtype=np.int32)

    net = lib.fastnsa_network_create(
        edges.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        1,
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        indices.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        2
    )

    assert net is not None and net != 0

    yield net

    lib.fastnsa_network_destroy(net)


@pytest.fixture(scope="function")
def graph_network():
    """
    构造一个最小可用 network：
    2 nodes, 1 edge
    """
    # edge: u, v, length
    edges = np.ascontiguousarray(
        [0.0, 1.0, 10.0],
        dtype=np.float64
    )

    # CSR: 0 -> 1
    offsets = np.ascontiguousarray([0, 1, 1], dtype=np.int32)
    indices = np.ascontiguousarray([1], dtype=np.int32)

    net = lib.fastnsa_network_create(
        edges.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        1,
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        indices.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        2
    )
    
    
    wkb_data = np.ascontiguousarray([1, 2, 3, 4], dtype=np.uint8)
    offsets = np.ascontiguousarray([0, 4], dtype=np.int32)

    lib.fastnsa_network_set_wkbs(
        net,
        wkb_data.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        1
    )

    assert net is not None and net != 0

    yield net

    lib.fastnsa_network_destroy(net)


@pytest.fixture(scope="function")
def simple_pointset(simple_network):
    """
    手动构造一个 PointSet（edge_id + offset）
    """
    pts = lib.fastnsa_pointset_create(simple_network)
    assert pts is not None and pts != 0

    edge_ids = np.ascontiguousarray([0], dtype=np.int32)
    offsets = np.ascontiguousarray([5.0], dtype=np.float64)

    ret = lib.fastnsa_pointset_set(
        pts,
        edge_ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        1
    )

    assert ret == 0

    yield pts

    lib.fastnsa_pointset_destroy(pts)
