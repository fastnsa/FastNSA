"""Optional geospatial file I/O tests."""

import tempfile
import unittest
from pathlib import Path

try:
    import geopandas as gpd
except ImportError:  # pragma: no cover - exercised in minimal wheel environments
    gpd = None

from fastnsa.io.points import load_points


@unittest.skipUnless(gpd is not None, "GeoPandas is not installed")
class PointIoTest(unittest.TestCase):
    def test_load_csv_with_custom_coordinates_and_datetime(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.csv"
            path.write_text(
                "east,north,occurred\n1,2,2026-01-01T00:00:00Z\n"
                "3,4,2026-01-02T00:00:00Z\n",
                encoding="utf-8",
            )
            points, bounds = load_points(
                path,
                return_bounding_box=True,
                time_col="occurred",
                x_col="east",
                y_col="north",
                crs="EPSG:3857",
            )
            self.assertEqual(str(points.crs), "EPSG:3857")
            self.assertEqual(bounds, (4.0, 2.0, 3.0, 1.0))
            self.assertEqual(
                float(points["occurred"].iloc[1] - points["occurred"].iloc[0]),
                86400.0,
            )

    def test_load_geojson_preserves_geometry_and_crs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.geojson"
            source = gpd.GeoDataFrame(
                {"value": [1, 2]},
                geometry=gpd.points_from_xy([1.0, 3.0], [2.0, 4.0]),
                crs="EPSG:4326",
            )
            source.to_file(path, driver="GeoJSON")
            loaded = load_points(path)
            self.assertEqual(len(loaded), 2)
            self.assertEqual(str(loaded.crs), "EPSG:4326")


if __name__ == "__main__":
    unittest.main()
