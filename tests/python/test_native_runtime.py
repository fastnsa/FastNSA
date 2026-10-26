"""Regression tests for the installed FastNSA native runtime."""

import struct
import heapq
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from fastnsa.bindings import capi
from fastnsa.runtime import NativeNetwork, NativePointSet
from fastnsa.stats.k_function import (
    network_k_function,
    spatiotemporal_network_k_function,
)
from fastnsa.stats.network_kde import network_kde
from fastnsa.core.types import (
    KFunctionResult,
    NetworkKDEPointResult,
    SpatiotemporalKFunctionResult,
)
from fastnsa.io.results import load_result


def _line_wkb(length: float) -> bytes:
    return struct.pack("<BII4d", 1, 2, 2, 0.0, 0.0, length, 0.0)


def _network(length: float = 100.0) -> NativeNetwork:
    wkb = _line_wkb(length)
    return NativeNetwork(
        edges=np.array([[0, 1, length]], dtype=np.float64),
        offsets=np.array([0, 1, 2], dtype=np.int32),
        indices=np.array([0, 0], dtype=np.int32),
        wkb_blob=wkb,
        wkb_offsets=np.array([0, len(wkb)], dtype=np.int32),
        crs="EPSG:3857",
    )


class NativeRuntimeTest(unittest.TestCase):
    def test_endpoint_roundoff_mapping_and_kde_reuse(self):
        length = np.nextafter(100.0, -np.inf)
        with NativeNetwork(
            edges=[[0, 1, length]], offsets=[0, 1, 2], indices=[0, 0],
            coordinates=[[0.0, 0.0], [100.0, 0.0]], geometry_offsets=[0, 2],
        ) as network:
            with NativePointSet.map_xy(network, [100.0, 0.0], [0.0, 0.0]) as mapped:
                np.testing.assert_array_equal(mapped.offsets_array(), [length, 0.0])
                with NativePointSet.from_mapped(
                    network, mapped.edge_ids_array(), mapped.offsets_array()
                ) as reused:
                    for method in ("ADA", "LION"):
                        with self.subTest(method=method):
                            direct = network_kde(network=network, points=mapped,
                                bandwidth=30, lixel_length=20, method=method)
                            restored = network_kde(network=network, points=reused,
                                bandwidth=30, lixel_length=20, method=method)
                            self.assertTrue(np.isfinite(direct.density).all())
                            np.testing.assert_allclose(restored.density, direct.density)
            with NativePointSet.from_mapped(
                network, [0, 0], [-1e-11, length + 1e-11], times=[0.0, 1.0]
            ) as legacy:
                np.testing.assert_array_equal(legacy.offsets_array(), [0.0, length])
            for offsets in ([-1e-5, length], [0.0, length + 1e-5]):
                with self.subTest(offsets=offsets):
                    with self.assertRaisesRegex(RuntimeError, "invalid edge offset"):
                        NativePointSet.from_mapped(network, [0, 0], offsets)

    def test_mapping_thread_policy_preserves_results_and_rejects_negative_values(self):
        query_x = np.linspace(0.0, 100.0, 1025)
        query_y = np.linspace(-5.0, 5.0, 1025)
        with _network() as network:
            outputs = []
            for thread_count in (0, 1, 4, 32, 512):
                with NativePointSet.map_xy(
                    network, query_x, query_y, n_threads=thread_count
                ) as points:
                    outputs.append(
                        (points.edge_ids_array(), points.offsets_array())
                    )
            for edge_ids, offsets in outputs[1:]:
                np.testing.assert_array_equal(edge_ids, outputs[0][0])
                np.testing.assert_allclose(offsets, outputs[0][1], atol=1e-12)
            with self.assertRaisesRegex(ValueError, "n_threads"):
                NativePointSet.map_xy(network, [1.0], [1.0], n_threads=-1)

    def test_packed_polylines_match_wkb_and_preserve_bent_edge_offsets(self):
        edges = np.array([[0, 1, 20.0], [1, 2, 10.0]], dtype=np.float64)
        adjacency_offsets = np.array([0, 1, 3, 4], dtype=np.int32)
        adjacency_indices = np.array([0, 0, 1, 1], dtype=np.int32)
        coordinates = np.array(
            [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [10.0, 10.0], [20.0, 10.0]],
            dtype=np.float64,
        )
        geometry_offsets = np.array([0, 3, 5], dtype=np.int32)
        bent_wkb = struct.pack(
            "<BII6d", 1, 2, 3, 0.0, 0.0, 10.0, 0.0, 10.0, 10.0
        )
        straight_wkb = struct.pack(
            "<BII4d", 1, 2, 2, 10.0, 10.0, 20.0, 10.0
        )
        packed_wkb = bent_wkb + straight_wkb
        wkb_offsets = np.array([0, len(bent_wkb), len(packed_wkb)], dtype=np.int32)

        common = {
            "edges": edges,
            "offsets": adjacency_offsets,
            "indices": adjacency_indices,
        }
        with NativeNetwork(
            **common,
            coordinates=coordinates,
            geometry_offsets=geometry_offsets,
        ) as packed_network, NativeNetwork(
            **common,
            wkb_blob=packed_wkb,
            wkb_offsets=wkb_offsets,
        ) as wkb_network:
            query_x = np.array([9.0, 11.0, 16.0])
            query_y = np.array([7.0, 9.0, 12.0])
            with NativePointSet.map_xy(packed_network, query_x, query_y) as packed_points, \
                    NativePointSet.map_xy(wkb_network, query_x, query_y) as wkb_points:
                np.testing.assert_array_equal(
                    packed_points.edge_ids_array(), wkb_points.edge_ids_array()
                )
                np.testing.assert_allclose(
                    packed_points.offsets_array(), wkb_points.offsets_array(), atol=1e-12
                )
                np.testing.assert_array_equal(packed_points.edge_ids_array(), [0, 0, 1])
                np.testing.assert_allclose(packed_points.offsets_array(), [17.0, 19.0, 6.0])

    def test_wkb_compatibility_accepts_big_endian_and_z_linestrings(self):
        edges = np.array([[0, 1, 10.0]], dtype=np.float64)
        common = {
            "edges": edges,
            "offsets": np.array([0, 1, 2], dtype=np.int32),
            "indices": np.array([0, 0], dtype=np.int32),
        }
        wkbs = (
            struct.pack(">BII4d", 0, 2, 2, 0.0, 0.0, 10.0, 0.0),
            struct.pack("<BII6d", 1, 1002, 2, 0.0, 0.0, 7.0, 10.0, 0.0, 9.0),
        )
        for wkb in wkbs:
            with self.subTest(first_bytes=wkb[:5]):
                with NativeNetwork(
                    **common,
                    wkb_blob=wkb,
                    wkb_offsets=np.array([0, len(wkb)], dtype=np.int32),
                ) as network:
                    with NativePointSet.map_xy(network, [4.0], [2.0]) as points:
                        self.assertEqual(points.edge_ids_array().tolist(), [0])
                        np.testing.assert_allclose(points.offsets_array(), [4.0])

    def test_comp_exact_matches_exact_bvh(self):
        rng = np.random.default_rng(410)
        starts = rng.uniform(-200.0, 200.0, size=(120, 2))
        ends = starts + rng.uniform(-20.0, 20.0, size=(120, 2))
        lengths = np.linalg.norm(ends - starts, axis=1)
        edges = np.column_stack((np.arange(120), np.arange(1, 121), lengths))
        coordinates = np.empty((240, 2), dtype=np.float64)
        coordinates[0::2] = starts
        coordinates[1::2] = ends
        queries = rng.uniform(-210.0, 210.0, size=(200, 2))
        with NativeNetwork(
            edges=edges,
            offsets=np.arange(122, dtype=np.int32),
            indices=np.concatenate((np.array([0]), np.arange(120))).astype(np.int32),
            coordinates=coordinates,
            geometry_offsets=np.arange(0, 241, 2, dtype=np.int32),
        ) as network:
            with NativePointSet.map_xy(
                network, queries[:, 0], queries[:, 1], mapping_method="R-tree"
            ) as exact, NativePointSet.map_xy(
                network,
                queries[:, 0],
                queries[:, 1],
                mapping_method="comp_e",
                comp_eps=8.0,
            ) as comparison:
                np.testing.assert_array_equal(
                    comparison.edge_ids_array(), exact.edge_ids_array()
                )
                np.testing.assert_allclose(
                    comparison.offsets_array(), exact.offsets_array(), atol=1e-10
                )

    def test_exact_bvh_matches_brute_force_for_random_segments(self):
        rng = np.random.default_rng(20261020)
        starts = rng.uniform(-1000.0, 1000.0, size=(500, 2))
        ends = starts + rng.uniform(-40.0, 40.0, size=(500, 2))
        lengths = np.linalg.norm(ends - starts, axis=1)
        edges = np.column_stack((np.arange(500), np.arange(1, 501), lengths))
        adjacency_offsets = np.arange(502, dtype=np.int32)
        adjacency_indices = np.concatenate((np.array([0]), np.arange(500))).astype(np.int32)
        coordinates = np.empty((1000, 2), dtype=np.float64)
        coordinates[0::2] = starts
        coordinates[1::2] = ends
        geometry_offsets = np.arange(0, 1001, 2, dtype=np.int32)
        queries = rng.uniform(-1100.0, 1100.0, size=(300, 2))

        segment_vectors = ends - starts
        segment_lengths_squared = np.sum(segment_vectors * segment_vectors, axis=1)
        relative = queries[:, None, :] - starts[None, :, :]
        fractions = np.sum(relative * segment_vectors[None, :, :], axis=2)
        fractions /= segment_lengths_squared[None, :]
        fractions = np.clip(fractions, 0.0, 1.0)
        projections = starts[None, :, :] + fractions[:, :, None] * segment_vectors[None, :, :]
        distances_squared = np.sum((queries[:, None, :] - projections) ** 2, axis=2)
        expected_edges = np.argmin(distances_squared, axis=1).astype(np.int32)
        expected_offsets = fractions[np.arange(len(queries)), expected_edges] * lengths[expected_edges]

        with NativeNetwork(
            edges=edges,
            offsets=adjacency_offsets,
            indices=adjacency_indices,
            coordinates=coordinates,
            geometry_offsets=geometry_offsets,
        ) as network:
            with NativePointSet.map_xy(network, queries[:, 0], queries[:, 1]) as points:
                np.testing.assert_array_equal(points.edge_ids_array(), expected_edges)
                np.testing.assert_allclose(points.offsets_array(), expected_offsets, atol=1e-10)

    def test_result_archives_round_trip_metadata_and_optional_arrays(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            k_result = KFunctionResult(
                r=np.array([0.0, 10.0]),
                k_obs=np.array([0.0, 12.0]),
                k_lo=np.array([0.0, 4.0]),
                k_hi=np.array([0.0, 8.0]),
                method="STARS",
                n_points=3,
                intensity=0.03,
                network_length=100.0,
                random_seed=42,
            )
            restored_k = KFunctionResult.load_npz(k_result.save_npz(root / "k_result"))
            self.assertIsNone(restored_k.k_theo)
            self.assertEqual(restored_k.method, "STARS")
            self.assertEqual(restored_k.random_seed, 42)
            self.assertEqual(restored_k.envelope_status.tolist(), [0, 1])
            self.assertEqual(restored_k.summary()["above_envelope"], 1)
            self.assertIsInstance(load_result(root / "k_result.npz"), KFunctionResult)
            csv_path = Path(restored_k.save_csv(root / "k_result.csv"))
            self.assertEqual(csv_path.read_text(encoding="utf-8").splitlines()[0], "r,k_obs,k_lo,k_hi")

            st_result = SpatiotemporalKFunctionResult(
                s=np.array([0.0, 10.0]),
                t=np.array([0.0, 1.0]),
                k_obs=np.array([[0.0, 2.0], [8.0, 4.0]]),
                k_lo=np.ones((2, 2)),
                k_hi=np.full((2, 2), 6.0),
                method="EARTH",
                n_points=4,
                time_unit="days",
                time_origin="2026-01-01",
                time_mode="datetime",
                random_seed=7,
            )
            restored_st = SpatiotemporalKFunctionResult.load_npz(
                st_result.save_npz(root / "st_result.npz")
            )
            np.testing.assert_array_equal(
                restored_st.envelope_status,
                np.array([[-1, 0], [1, 0]], dtype=np.int8),
            )
            self.assertEqual(len(restored_st.to_records()), 4)
            self.assertEqual(restored_st.time_unit, "days")

            kde_result = NetworkKDEPointResult(
                edge_ids=np.array([0, 0], dtype=np.int32),
                offsets_n1=np.array([5.0, 12.5]),
                offsets_n2=np.array([10.0, 2.5]),
                density=np.array([0.05, 0.1]),
                method="ADA",
                kernel="EPANECHNIKOV",
                bandwidth=20.0,
                lixel_length=10.0,
                n_points=2,
                network_length=15.0,
                normalized=True,
                normalization="network_integral",
            )
            restored_kde = NetworkKDEPointResult.load_npz(
                kde_result.save_npz(root / "kde_result")
            )
            self.assertEqual(restored_kde.n_points, 2)
            self.assertEqual(restored_kde.normalization, "network_integral")
            self.assertAlmostEqual(restored_kde.summary()["density_integral"], 1.0)

            with self.assertRaisesRegex(ValueError, "Expected a 'network_kde'"):
                NetworkKDEPointResult.load_npz(root / "k_result.npz")

    def test_nkdv_unsorted_points_self_loops_and_shortcuts(self):
        # Edge 0 is a self-loop, edges 1 and 2 are parallel with a shorter
        # alternative route. Both break a direct-distance-only same-edge path.
        edges = np.array([[0, 0, 100.0], [0, 1, 100.0], [0, 1, 10.0], [1, 2, 30.0]])
        adjacency = [[0, 0, 1, 2], [1, 2, 3], [3]]
        offsets = np.array([0, 4, 7, 8], dtype=np.int32)
        indices = np.array([e for row in adjacency for e in row], dtype=np.int32)
        wkbs = [_line_wkb(length) for _, _, length in edges]
        wkb_offsets = np.array([0] + list(np.cumsum([len(w) for w in wkbs])), dtype=np.int32)
        point_edges = np.array([1, 0, 2, 0, 3, 1, 0, 1, 0], dtype=np.int32)
        point_offsets = np.array([95, 95, 5, 5, 10, 5, 95, 50, 60], dtype=float)
        bandwidth = 40.0

        def exact(edge_id, x):
            u, v, length = edges[edge_id]
            distances = [float('inf')] * 3
            queue = []
            for node, distance in ((int(u), x), (int(v), length - x)):
                if distance < distances[node]:
                    distances[node] = distance
                    heapq.heappush(queue, (distance, node))
            while queue:
                distance, node = heapq.heappop(queue)
                if distance != distances[node]:
                    continue
                for other_edge in adjacency[node]:
                    a, b, weight = edges[other_edge]
                    other = int(b) if int(a) == node else int(a)
                    candidate = distance + weight
                    if candidate < distances[other]:
                        distances[other] = candidate
                        heapq.heappush(queue, (candidate, other))
            total = 0.0
            for host, point in zip(point_edges, point_offsets):
                a, b, weight = edges[host]
                distance = min(distances[int(a)] + point, distances[int(b)] + weight - point)
                if host == edge_id:
                    distance = min(distance, abs(x - point))
                if distance < bandwidth:
                    total += 1.0 - (distance / bandwidth) ** 2
            return total

        with NativeNetwork(edges=edges, offsets=offsets, indices=indices,
                           wkb_blob=b''.join(wkbs), wkb_offsets=wkb_offsets) as network:
            expected = None
            for permutation in (np.arange(len(point_edges)), np.arange(len(point_edges))[::-1]):
                with NativePointSet.from_mapped(network, point_edges[permutation], point_offsets[permutation]) as points:
                    for method in ('BASIC', 'ADA', 'LION'):
                        result = network_kde(network=network, points=points, bandwidth=bandwidth,
                                             lixel_length=20, kernel='EPANECHNIKOV', method=method)
                        reference = np.array([exact(int(e), x) for e, x in zip(result.edge_ids, result.offsets_n1)])
                        np.testing.assert_allclose(result.density, reference, rtol=1e-12, atol=1e-9,
                                                   err_msg=f'{method} with permutation {permutation}')
                        if expected is None:
                            expected = result.density.copy()
                        else:
                            np.testing.assert_allclose(result.density, expected, rtol=1e-12, atol=1e-9)

    def test_maps_points_and_runs_statistics(self):
        with _network() as network:
            with NativePointSet.map_xy(
                network,
                [10.0, 40.0, 80.0],
                [1.0, -1.0, 2.0],
            ) as points:
                self.assertEqual(points.edge_ids_array().tolist(), [0, 0, 0])
                self.assertTrue(
                    np.allclose(points.offsets_array(), [10.0, 40.0, 80.0])
                )

                kde = network_kde(
                    network=network,
                    points=points,
                    bandwidth=30.0,
                    lixel_length=10.0,
                    method="ADA",
                    kernel="EPANECHNIKOV",
                )
                self.assertEqual(len(kde.edge_ids), 10)
                self.assertTrue(np.all(kde.density >= 0.0))

                result = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([10.0, 20.0]),
                    simulations=2,
                    return_simulations=True,
                )
                self.assertEqual(result.k_obs.shape, (2,))
                self.assertEqual(result.k_sim.shape, (2, 2))
                self.assertEqual(result.k_lo.shape, (2,))
                self.assertEqual(result.k_hi.shape, (2,))

    def test_k_function_is_normalized_by_network_length(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                result = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([10.0, 30.0, 40.0, 70.0]),
                    simulations=0,
                )
                self.assertTrue(
                    np.allclose(result.k_obs, [0.0, 100.0 / 3.0, 200.0 / 3.0, 100.0])
                )
                sps_result = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([10.0, 30.0, 40.0, 70.0]),
                    method="SPS",
                    simulations=0,
                )
                self.assertTrue(np.allclose(result.k_obs, sps_result.k_obs))
                stars_result = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([10.0, 30.0, 40.0, 70.0]),
                    method="STARS",
                    simulations=0,
                    epsilon=0.1,
                )
                self.assertTrue(np.allclose(result.k_obs, stars_result.k_obs))

    def test_nkdv_rejects_unvalidated_kernel(self):
        with _network() as network:
            with NativePointSet.from_mapped(network, [0, 0], [20.0, 80.0]) as points:
                with self.assertRaisesRegex(ValueError, "only supports the Epanechnikov"):
                    network_kde(
                        network=network,
                        points=points,
                        bandwidth=30.0,
                        lixel_length=10.0,
                        method="ADA",
                        kernel="QUARTIC",
                    )

    def test_native_nkdv_rejects_unvalidated_kernel(self):
        with _network() as network:
            with NativePointSet.from_mapped(network, [0, 0], [20.0, 80.0]) as points:
                n_lixels = capi.lib.fastnsa_nkdv_count_lixels(network._ptr, 10.0)
                edge_ids = np.empty(n_lixels, dtype=np.int32)
                offsets_n1 = np.empty(n_lixels, dtype=np.float64)
                offsets_n2 = np.empty(n_lixels, dtype=np.float64)
                density = np.empty(n_lixels, dtype=np.float64)
                status = capi.lib.fastnsa_compute_nkdv(
                    network._ptr,
                    points._ptr,
                    capi.FASTNSA_NKDV_ADA,
                    10.0,
                    capi.FASTNSA_NKDV_KERNEL_QUARTIC,
                    30.0,
                    n_lixels,
                    edge_ids.ctypes.data_as(capi.P_INT),
                    offsets_n1.ctypes.data_as(capi.P_DOUBLE),
                    offsets_n2.ctypes.data_as(capi.P_DOUBLE),
                    density.ctypes.data_as(capi.P_DOUBLE),
                )
                self.assertNotEqual(status, 0)
                self.assertIn("only supports the Epanechnikov", capi.last_error())

    def test_nkdv_rejects_excessive_lixel_output(self):
        with _network() as network:
            with NativePointSet.from_mapped(network, [0, 0], [20.0, 80.0]) as points:
                with self.assertRaisesRegex(RuntimeError, "excessive output size"):
                    network_kde(
                        network=network,
                        points=points,
                        bandwidth=30.0,
                        lixel_length=1e-300,
                        method="ADA",
                    )

    def test_lion_matches_ada_with_a_short_final_lixel(self):
        first_wkb = _line_wkb(95.0)
        second_wkb = _line_wkb(67.0)
        with NativeNetwork(
            edges=np.array([[0, 1, 95.0], [1, 2, 67.0]], dtype=np.float64),
            offsets=np.array([0, 1, 3, 4], dtype=np.int32),
            indices=np.array([0, 0, 1, 1], dtype=np.int32),
            wkb_blob=first_wkb + second_wkb,
            wkb_offsets=np.array(
                [0, len(first_wkb), len(first_wkb) + len(second_wkb)],
                dtype=np.int32,
            ),
        ) as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 1, 1],
                offsets=[2.0, 90.0, 3.0, 65.0],
            ) as points:
                kwargs = {
                    "network": network,
                    "points": points,
                    "bandwidth": 35.0,
                    "lixel_length": 20.0,
                    "kernel": "EPANECHNIKOV",
                }
                ada = network_kde(method="ADA", **kwargs)
                lion = network_kde(method="LION", **kwargs)
                self.assertTrue(np.all(lion.density >= 0.0))
                self.assertTrue(np.allclose(ada.density, lion.density, atol=1e-10))

    def test_weighted_nkdv_is_linear_and_ada_matches_lion(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                kwargs = {
                    "network": network,
                    "points": points,
                    "bandwidth": 35.0,
                    "lixel_length": 10.0,
                    "kernel": "EPANECHNIKOV",
                }
                unweighted = network_kde(method="ADA", **kwargs)
                ones = network_kde(method="ADA", weights=[1.0, 1.0, 1.0], **kwargs)
                doubled = network_kde(method="ADA", weights=[2.0, 2.0, 2.0], **kwargs)
                lion = network_kde(method="LION", weights=[1.0, 3.0, 0.5], **kwargs)
                ada = network_kde(method="ADA", weights=[1.0, 3.0, 0.5], **kwargs)
                basic = network_kde(method="BASIC", weights=[1.0, 3.0, 0.5], **kwargs)
                np.testing.assert_allclose(ones.density, unweighted.density, atol=1e-10)
                np.testing.assert_allclose(doubled.density, 2.0 * unweighted.density, atol=1e-10)
                np.testing.assert_allclose(lion.density, ada.density, atol=1e-10)
                np.testing.assert_allclose(basic.density, ada.density, atol=1e-10)
                self.assertTrue(ada.weighted)
                self.assertEqual(ada.total_weight, 4.5)

    def test_monte_carlo_samples_proportionally_to_edge_length(self):
        short_wkb = _line_wkb(1.0)
        long_wkb = _line_wkb(99.0)
        with NativeNetwork(
            edges=np.array([[0, 1, 1.0], [2, 3, 99.0]], dtype=np.float64),
            offsets=np.array([0, 1, 2, 3, 4], dtype=np.int32),
            indices=np.array([0, 0, 1, 1], dtype=np.int32),
            wkb_blob=short_wkb + long_wkb,
            wkb_offsets=np.array([0, len(short_wkb), len(short_wkb) + len(long_wkb)], dtype=np.int32),
        ) as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[1] * 200,
                offsets=np.linspace(0.0, 99.0, 200),
            ) as points:
                result = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([1.0]),
                    simulations=1,
                    return_simulations=True,
                )
                # Uniform edge selection would place roughly half the events on
                # the one-unit edge and inflate this statistic above 20.
                self.assertLess(float(result.k_sim[0, 0]), 10.0)

    def test_rejects_invalid_network_arrays(self):
        wkb = _line_wkb(100.0)
        with self.assertRaises((RuntimeError, ValueError)):
            NativeNetwork(
                edges=np.array([[0.5, 1, 100.0]], dtype=np.float64),
                offsets=np.array([0, 1, 2], dtype=np.int32),
                indices=np.array([0, 0], dtype=np.int32),
                wkb_blob=wkb,
                wkb_offsets=np.array([0, len(wkb)], dtype=np.int32),
            )

    def test_rejects_non_finite_point_coordinates(self):
        with _network() as network:
            with self.assertRaisesRegex(ValueError, "NaN or infinite"):
                NativePointSet.map_xy(network, [float("nan")], [0.0])

    def test_earth_rejects_single_threshold_axes(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0],
                offsets=[10.0, 20.0],
                times=[0.0, 1.0],
            ) as points:
                with self.assertRaisesRegex(ValueError, "at least two thresholds"):
                    spatiotemporal_network_k_function(
                        network=network,
                        points=points,
                        s_values=np.array([10.0]),
                        t_values=np.array([1.0]),
                        method="EARTH",
                    )

    def test_spatiotemporal_k_function_requires_two_points(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0],
                offsets=[10.0],
                times=[0.0],
            ) as points:
                with self.assertRaisesRegex(ValueError, "at least two points"):
                    spatiotemporal_network_k_function(
                        network=network,
                        points=points,
                        s_values=np.array([0.0, 10.0]),
                        t_values=np.array([0.0, 1.0]),
                    )

    def test_spatiotemporal_k_function_is_normalized_by_space_and_time(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0],
                offsets=[10.0, 90.0],
                times=[0.0, 10.0],
            ) as points:
                result = spatiotemporal_network_k_function(
                    network=network,
                    points=points,
                    s_values=np.array([0.0, 100.0]),
                    t_values=np.array([0.0, 10.0]),
                    method="RQS",
                )
                self.assertAlmostEqual(float(result.k_obs[-1, -1]), 1000.0)

    def test_spatiotemporal_k_function_rejects_constant_times(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0],
                offsets=[10.0, 90.0],
                times=[5.0, 5.0],
            ) as points:
                with self.assertRaisesRegex(ValueError, "positive time span"):
                    spatiotemporal_network_k_function(
                        network=network,
                        points=points,
                        s_values=np.array([0.0, 100.0]),
                        t_values=np.array([0.0, 10.0]),
                        method="RQS",
                    )

    def test_spatiotemporal_methods_agree_on_small_network(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 90.0],
                times=[0.0, 4.0, 10.0],
            ) as points:
                kwargs = {
                    "network": network,
                    "points": points,
                    "s_values": np.array([0.0, 50.0, 100.0]),
                    "t_values": np.array([0.0, 5.0, 10.0]),
                }
                expected = spatiotemporal_network_k_function(method="RQS", **kwargs)
                for method in ("SPS", "EAR", "MTS", "EARTH"):
                    actual = spatiotemporal_network_k_function(method=method, **kwargs)
                    self.assertTrue(
                        np.allclose(expected.k_obs, actual.k_obs),
                        msg=f"{method} differs from RQS",
                    )

    def test_spatiotemporal_text_export_preserves_native_time_unit(self):
        with _network() as network, tempfile.TemporaryDirectory() as directory:
            for scale in (1.0, 86400.0):
                with NativePointSet.from_mapped(
                    network, edge_ids=[0, 0], offsets=[90.0, 10.0],
                    times=np.array([0.125, 0.375]) * scale,
                ) as points:
                    path = Path(directory) / "events.txt"
                    status = capi.lib.fastnsa_pointset_dump_text_st(
                        points._ptr, os.fsencode(path)
                    )
                    self.assertEqual(status, 0)
                    row = path.read_text().splitlines()[1].split()
                    self.assertEqual(float(row[4]), 10.0)
                    self.assertEqual(float(row[5]), 0.25 * scale)
                    self.assertEqual(float(row[6]), 90.0)
                    self.assertEqual(float(row[7]), 0.0)

    def test_spatiotemporal_methods_match_pair_counts_across_time_units(self):
        # Multiple populated edges also exercise EARTH's MTS counting path.
        lengths = np.array([100.0, 80.0, 120.0])
        edge_ids = np.repeat(np.arange(3), 12)
        offsets = np.concatenate([
            np.linspace(2.0, length - 2.0, 12) for length in lengths
        ])
        positions = np.concatenate(([0.0], np.cumsum(lengths)))
        locations = positions[edge_ids] + offsets
        times = ((np.arange(36) * 7) % 23).astype(float) / 2.0
        spatial = np.array([15.0, 75.0, 135.0, 195.0])
        pair_distances = np.abs(locations[:, None] - locations)
        pair_times = np.abs(times[:, None] - times)
        off_diagonal = ~np.eye(len(times), dtype=bool)
        with NativeNetwork(
            edges=np.array([[0, 1, 100.0], [1, 2, 80.0], [2, 3, 120.0]]),
            offsets=np.array([0, 1, 3, 5, 6], dtype=np.int32),
            indices=np.array([0, 0, 1, 1, 2, 2], dtype=np.int32),
            coordinates=np.column_stack((
                positions[[0, 1, 1, 2, 2, 3]], np.zeros(6)
            )),
            geometry_offsets=np.array([0, 2, 4, 6], dtype=np.int32),
        ) as network:
            for temporal in (np.arange(4.0), np.arange(1.0, 5.0)):
                expected = np.array([
                    [np.count_nonzero(
                        off_diagonal & (pair_distances <= s) & (pair_times <= t)
                    ) for t in temporal] for s in spatial
                ]) * lengths.sum() * np.ptp(times) / (len(times) * (len(times) - 1))
                for scale, origin in ((1.0, 0.0), (86400.0, 0.0),
                                      (86400000.0, 1700000000000.0)):
                    with NativePointSet.from_mapped(
                        network, edge_ids=edge_ids, offsets=offsets,
                        times=times * scale + origin,
                    ) as points:
                        for method in ("RQS", "SPS", "EAR", "MTS", "EARTH"):
                            with self.subTest(method=method, scale=scale,
                                              origin=origin, start=temporal[0]):
                                result = spatiotemporal_network_k_function(
                                    network=network, points=points,
                                    s_values=spatial, t_values=temporal * scale,
                                    method=method,
                                )
                                np.testing.assert_allclose(result.k_obs / scale, expected)

    def test_spatiotemporal_monte_carlo_methods_and_time_units_agree(self):
        with _network() as network:
            baseline = None
            for scale in (1.0, 86400.0):
                with NativePointSet.from_mapped(
                    network, edge_ids=[0] * 12,
                    offsets=np.linspace(1.0, 99.0, 12),
                    times=np.arange(12.0) * scale,
                ) as points:
                    for method in ("RQS", "SPS", "EAR", "MTS", "EARTH"):
                        with self.subTest(method=method, scale=scale):
                            result = spatiotemporal_network_k_function(
                                network=network, points=points,
                                s_values=np.array([10.0, 40.0, 70.0, 100.0]),
                                t_values=np.array([1.0, 4.0, 7.0, 10.0]) * scale,
                                method=method, simulations=3, random_state=42,
                                return_simulations=True,
                            )
                            surfaces = [result.k_obs, result.k_sim,
                                        result.k_theo, result.k_lo, result.k_hi]
                            if baseline is None:
                                baseline = [surface.copy() for surface in surfaces]
                            for actual, expected in zip(surfaces, baseline):
                                np.testing.assert_allclose(actual / scale, expected)

    def test_rqs_and_sps_progress_is_global_and_monotonic(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 90.0],
                times=[0.0, 4.0, 10.0],
            ) as points:
                for method in ("RQS", "SPS"):
                    updates = []

                    def progress(_stage, completed, total):
                        updates.append((completed, total))
                        return False

                    spatiotemporal_network_k_function(
                        network=network,
                        points=points,
                        s_values=np.array([0.0, 50.0, 100.0]),
                        t_values=np.array([0.0, 5.0, 10.0]),
                        method=method,
                        simulations=1,
                        progress_callback=progress,
                    )

                    self.assertTrue(updates)
                    totals = {total for _completed, total in updates}
                    self.assertEqual(len(totals), 1)
                    self.assertGreater(next(iter(totals)), 0)
                    completed = [value for value, _total in updates]
                    self.assertEqual(completed, sorted(completed))

    def test_checked_c_api_is_available(self):
        self.assertGreaterEqual(capi.lib.fastnsa_abi_version(), 1)
        self.assertTrue(capi.HAS_CHECKED_WKB_SETTER)
        self.assertTrue(capi.HAS_CHECKED_K_FUNCTION)
        self.assertTrue(capi.HAS_CHECKED_ST_K_FUNCTION)
        self.assertTrue(capi.HAS_PROGRESS_K_FUNCTION)
        self.assertTrue(capi.HAS_PROGRESS_ST_K_FUNCTION)
        self.assertTrue(capi.HAS_PROGRESS_NKDV)

    def test_monte_carlo_is_deterministic_for_the_default_seed(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                first = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([10.0, 30.0, 70.0]),
                    simulations=3,
                    return_simulations=True,
                )
                second = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([10.0, 30.0, 70.0]),
                    simulations=3,
                    return_simulations=True,
                )
                self.assertTrue(np.array_equal(first.k_sim, second.k_sim))

    def test_monte_carlo_accepts_reproducible_custom_seed(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                kwargs = {
                    "network": network,
                    "points": points,
                    "r_values": np.array([10.0, 30.0, 70.0]),
                    "simulations": 4,
                    "return_simulations": True,
                    "random_state": 12345,
                }
                first = network_k_function(**kwargs)
                second = network_k_function(**kwargs)
                self.assertEqual(first.random_seed, 12345)
                self.assertTrue(np.array_equal(first.k_sim, second.k_sim))

    def test_global_envelope_reports_confidence_and_simulation_p_value(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                result = network_k_function(
                    network=network,
                    points=points,
                    r_values=np.array([10.0, 30.0, 70.0]),
                    simulations=5,
                    envelope="global",
                    alpha=0.1,
                )
                self.assertEqual(result.envelope_method, "global")
                self.assertAlmostEqual(result.confidence_level, 0.9)
                self.assertGreaterEqual(result.p_value, 1.0 / 6.0)
                self.assertLessEqual(result.p_value, 1.0)
                self.assertIn("90.0% global", result.envelope_description)

    def test_st_xarray_export_has_an_actionable_optional_dependency(self):
        result = SpatiotemporalKFunctionResult(
            s=np.array([0.0, 1.0]),
            t=np.array([0.0, 1.0]),
            k_obs=np.zeros((2, 2)),
        )
        try:
            data = result.to_xarray()
        except ImportError as exc:
            self.assertIn("fastnsa[xarray]", str(exc))
        else:
            self.assertEqual(data.dims, ("spatial_threshold", "temporal_threshold"))

    def test_normalized_nkdv_integrates_using_short_final_lixels(self):
        with _network(length=95.0) as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                result = network_kde(
                    network=network,
                    points=points,
                    bandwidth=30.0,
                    lixel_length=20.0,
                    method="ADA",
                    normalize=True,
                )
                self.assertEqual(result.normalization, "network_integral")
                self.assertAlmostEqual(
                    float(np.sum(result.density * result.segment_lengths)),
                    1.0,
                    places=12,
                )

    def test_kfunction_search_cancellation_raises_typed_error_and_can_rerun(self):
        from fastnsa.exceptions import FastNSACancelledError

        with _network() as network:
            with NativePointSet.from_mapped(
                network, edge_ids=[0, 0, 0], offsets=[10.0, 40.0, 80.0],
                times=[0.0, 1.0, 3.0],
            ) as points:
                for temporal, methods in (
                    (False, ("NS", "STARS", "SPS")),
                    (True, ("RQS", "SPS", "EAR", "MTS", "EARTH")),
                ):
                    function = spatiotemporal_network_k_function if temporal else network_k_function
                    thresholds = {"s_values": [100.0, 200.0], "t_values": [0.0, 3.0]} \
                        if temporal else {"r_values": [100.0, 200.0]}
                    for method in methods:
                        with self.subTest(temporal=temporal, method=method):
                            kwargs = dict(network=network, points=points, method=method,
                                          simulations=2, return_simulations=True, **thresholds)
                            baseline = function(**kwargs)
                            polls = {}

                            def cancel(stage, completed, total):
                                if stage != "Monte Carlo batch" or completed <= 0:
                                    return False
                                key = (completed, total)
                                polls[key] = polls.get(key, 0) + 1
                                # Repeated progress occurs inside searches, not
                                # at the advancing source-edge boundary.
                                return polls[key] >= 2

                            with self.assertRaises(FastNSACancelledError):
                                function(progress_callback=cancel, **kwargs)
                            self.assertTrue(any(count >= 2 for count in polls.values()))
                            recovered = function(**kwargs)
                            np.testing.assert_array_equal(recovered.k_obs, baseline.k_obs)
                            np.testing.assert_array_equal(recovered.k_sim, baseline.k_sim)

    def test_native_progress_callback_can_cancel(self):
        stages = []

        def cancel(stage, completed, total):
            stages.append((stage, completed, total))
            return True

        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                with self.assertRaisesRegex(RuntimeError, "canceled"):
                    network_k_function(
                        network=network,
                        points=points,
                        r_values=np.array([10.0, 30.0, 70.0]),
                        simulations=3,
                        progress_callback=cancel,
                    )
        self.assertTrue(stages)

    def test_native_nkdv_progress_callback_can_cancel(self):
        with _network() as network:
            with NativePointSet.from_mapped(
                network,
                edge_ids=[0, 0, 0],
                offsets=[10.0, 40.0, 80.0],
            ) as points:
                with self.assertRaisesRegex(RuntimeError, "canceled"):
                    network_kde(
                        network=network,
                        points=points,
                        bandwidth=30.0,
                        lixel_length=10.0,
                        method="LION",
                        progress_callback=lambda *_args: True,
                    )


if __name__ == "__main__":
    unittest.main()
