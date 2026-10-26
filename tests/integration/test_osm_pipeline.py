"""Opt-in OSM download and analysis pipeline test."""

import os
import unittest

import numpy as np

try:
    import geopandas as gpd
    from shapely.geometry import Point
except ImportError:  # pragma: no cover - depends on integration environment
    gpd = None
    Point = None

import fastnsa as nsa


RUN_NETWORK_TESTS = os.environ.get("FASTNSA_RUN_NETWORK_TESTS") == "1"


@unittest.skipUnless(RUN_NETWORK_TESTS, "set FASTNSA_RUN_NETWORK_TESTS=1")
@unittest.skipUnless(gpd is not None and Point is not None, "OSM/geo extras are not installed")
class OsmPipelineIntegrationTest(unittest.TestCase):
    def test_osm_mapping_and_k_function(self):
        with nsa.Network.from_osm(place="Central Park, New York, USA") as network:
            points = gpd.GeoDataFrame(
                {"event_id": [1, 2]},
                geometry=[Point(-73.9712, 40.7831), Point(-73.9680, 40.7815)],
                crs="EPSG:4326",
            )
            with nsa.PointEvents(points, target_network=network) as events:
                result = nsa.network_k_function(
                    network=network,
                    points=events,
                    r_values=[10.0, 50.0, 100.0],
                    method="STARS",
                )
                self.assertEqual(result.r.shape, (3,))
                self.assertEqual(result.k_obs.shape, (3,))
                self.assertTrue(np.isfinite(result.k_obs).all())


if __name__ == "__main__":
    unittest.main()
