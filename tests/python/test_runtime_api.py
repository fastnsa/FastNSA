"""Tests for the public Python wrappers around the native runtime."""

import struct
import unittest

import numpy as np

from fastnsa import Network, NativeNetwork, NativePointSet, PointEvents, preprocess_events
from fastnsa.bindings import capi
from fastnsa.stats.k_function import network_k_function


def _network_arrays(length=100.0):
    wkb = struct.pack("<BII4d", 1, 2, 2, 0.0, 0.0, length, 0.0)
    return {
        "edges": np.array([[0, 1, length]], dtype=np.float64),
        "offsets": np.array([0, 1, 2], dtype=np.int32),
        "indices": np.array([0, 0], dtype=np.int32),
        "wkb_blob": wkb,
        "wkb_offsets": np.array([0, len(wkb)], dtype=np.int32),
    }


class _Points:
    crs = "EPSG:3857"
    columns = ()

    def __len__(self):
        return 3


class _MappedPoints(_Points):
    columns = ("edge_id", "offset")

    def __getitem__(self, name):
        return {
            "edge_id": np.array([0, 0, 0]),
            "offset": np.array([10.0, 50.0, 90.0]),
        }[name]


class PythonRuntimeApiTest(unittest.TestCase):
    def test_preprocess_events_reuses_mapping_columns(self):
        with Network.from_prepared_arrays(
            crs="EPSG:3857", **_network_arrays()
        ) as network:
            with preprocess_events(
                _MappedPoints(),
                network,
                edge_id_col="edge_id",
                offset_col="offset",
            ) as prepared:
                self.assertTrue(prepared.reused_mapping)
                self.assertEqual(prepared.mapping_method, "precomputed")
                self.assertEqual(prepared.events.edge_ids_array().tolist(), [0, 0, 0])
    def test_native_abi_matches_python_bindings(self):
        self.assertEqual(capi.NATIVE_ABI_VERSION, capi.REQUIRED_ABI_VERSION)

    def test_legacy_wrappers_delegate_to_native_runtime(self):
        with Network.from_prepared_arrays(
            crs="EPSG:3857", **_network_arrays()
        ) as network:
            self.assertIsInstance(network._native, NativeNetwork)
            with PointEvents.from_mapped_arrays(
                _Points(),
                network,
                [0, 0, 0],
                [10.0, 50.0, 90.0],
                times=[0.0, 1.0, 2.0],
            ) as points:
                self.assertIsInstance(points._native, NativePointSet)
                self.assertEqual(points.edge_ids_array().tolist(), [0, 0, 0])
                result = network_k_function(
                    network=network,
                    points=points,
                    r_values=[0.0, 50.0],
                )
                self.assertEqual(result.k_obs.shape, (2,))

    def test_closed_native_objects_are_rejected(self):
        network = NativeNetwork(**_network_arrays())
        points = NativePointSet.from_mapped(network, [0, 0], [10.0, 20.0])
        points.close()
        with self.assertRaisesRegex(RuntimeError, "remain open"):
            network_k_function(
                network=network,
                points=points,
                r_values=[0.0, 10.0],
            )
        network.close()
        with self.assertRaisesRegex(RuntimeError, "closed network"):
            NativePointSet.from_mapped(network, [0], [1.0])

    def test_point_set_cannot_be_used_with_another_network(self):
        with NativeNetwork(**_network_arrays()) as first:
            with NativeNetwork(**_network_arrays()) as second:
                with NativePointSet.from_mapped(first, [0, 0], [10.0, 20.0]) as points:
                    with self.assertRaisesRegex(ValueError, "different network"):
                        network_k_function(
                            network=second,
                            points=points,
                            r_values=[0.0, 10.0],
                        )


if __name__ == "__main__":
    unittest.main()
