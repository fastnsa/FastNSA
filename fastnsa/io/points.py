import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point


def load_points(points_path, return_bounding_box=False):
    """
    Load input points and infer bounding box.

    Parameters
    ----------
    points_path : str
        Path to CSV file containing lon/lat columns.

    Returns
    -------
    points : GeoSeries
        Point geometries (EPSG:4326).
    bounding_box : tuple
        (lat_max, lat_min, lon_max, lon_min)
    """

    ext = os.path.splitext(points_path)[-1].lower()

    if ext == ".csv":
        df = pd.read_csv(points_path)
        if "lon" not in df.columns or "lat" not in df.columns:
            raise ValueError("CSV must contain 'lon' and 'lat' columns")
        
        if return_bounding_box:

            lat_max = df["lat"].max()
            lat_min = df["lat"].min()
            lon_max = df["lon"].max()
            lon_min = df["lon"].min()
            bounding_box = (lat_max, lat_min, lon_max, lon_min)

        points = gpd.GeoSeries(
            [Point(xy) for xy in zip(df["lon"], df["lat"])],
            crs="EPSG:4326",
        )
    else:
        raise ValueError(f"Unsupported file format: {ext}")
    
    if return_bounding_box:
        return points, bounding_box
    
    return points
