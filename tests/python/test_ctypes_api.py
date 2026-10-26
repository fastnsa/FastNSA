"""Focused checks for the low-level ctypes binding contract."""

import unittest

import numpy as np

from fastnsa.bindings import capi


class CtypesSmokeTest(unittest.TestCase):
    def test_checked_k_function_call(self):
        edges = np.array([[0, 1, 10.0]], dtype=np.float64)
        csr_offsets = np.array([0, 1, 2], dtype=np.int32)
        csr_indices = np.array([0, 0], dtype=np.int32)

        network = capi.lib.fastnsa_network_create(
            edges.ctypes.data_as(capi.P_DOUBLE),
            1,
            csr_offsets.ctypes.data_as(capi.P_INT),
            csr_indices.ctypes.data_as(capi.P_INT),
            2,
        )
        self.assertTrue(network, capi.last_error())

        points = capi.lib.fastnsa_pointset_create(network)
        self.assertTrue(points, capi.last_error())
        try:
            edge_ids = np.array([0, 0], dtype=np.int32)
            offsets = np.array([2.0, 8.0], dtype=np.float64)
            status = capi.lib.fastnsa_pointset_set(
                points,
                edge_ids.ctypes.data_as(capi.P_INT),
                offsets.ctypes.data_as(capi.P_DOUBLE),
                2,
            )
            self.assertEqual(status, 0, capi.last_error())

            r_values = np.array([1.0, 2.0, 5.0], dtype=np.float64)
            observed = np.zeros(len(r_values), dtype=np.float64)
            function = (
                capi.lib.fastnsa_compute_k_function_checked
                if capi.HAS_CHECKED_K_FUNCTION
                else capi.lib.fastnsa_compute_k_function
            )
            status = function(
                network,
                points,
                r_values.ctypes.data_as(capi.P_DOUBLE),
                len(r_values),
                capi.FASTNSA_K_NS,
                0,
                observed.ctypes.data_as(capi.P_DOUBLE),
                None,
                0.1,
            )
            if capi.HAS_CHECKED_K_FUNCTION:
                self.assertEqual(status, 0, capi.last_error())
            self.assertTrue(np.isfinite(observed).all())
        finally:
            capi.lib.fastnsa_pointset_destroy(points)
            capi.lib.fastnsa_network_destroy(network)


if __name__ == "__main__":
    unittest.main()

