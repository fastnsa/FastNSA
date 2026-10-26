"""Tests requiring the optional NetworkX and Shapely dependencies."""

import unittest

try:
    import networkx as nx
    from shapely.geometry import LineString
except ImportError:  # pragma: no cover - exercised in minimal wheel environments
    nx = None
    LineString = None

from fastnsa import Network


@unittest.skipUnless(nx is not None and LineString is not None, "geo extras are not installed")
class OptionalNetworkTest(unittest.TestCase):
    def test_multigraph_builds_native_topology_and_geometry(self):
        graph = nx.MultiGraph(crs="EPSG:3857")
        graph.add_node("west", x=0.0, y=0.0)
        graph.add_node("center", x=10.0, y=0.0)
        graph.add_node("north", x=10.0, y=10.0)
        graph.add_edge(
            "west",
            "center",
            length=10.0,
            geometry=LineString([(0.0, 0.0), (10.0, 0.0)]),
        )
        graph.add_edge(
            "center",
            "north",
            length=10.0,
            geometry=LineString([(10.0, 0.0), (10.0, 10.0)]),
        )

        with Network(graph) as network:
            self.assertEqual(network.num_nodes, 3)
            self.assertEqual(network.num_edges, 2)
            self.assertEqual(network.total_length, 20.0)
            self.assertEqual(network._data_cache["offsets"].shape, (4,))


if __name__ == "__main__":
    unittest.main()
