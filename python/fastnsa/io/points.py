"""Point dataset loading with optional GeoPandas integration."""

from pathlib import Path
import logging

import numpy as np


logger = logging.getLogger(__name__)


def _geo_dependencies():
    try:
        import geopandas as gpd
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            'Point file loading requires the "geo" extra: pip install "fastnsa[geo]".'
        ) from exc
    return gpd, pd


def _prepare_time_column(points, time_col, pd):
    if time_col is None:
        return points
    if time_col not in points.columns:
        raise ValueError(f"Point data must contain time column {time_col!r}.")
    points = points.copy()
    if pd.api.types.is_numeric_dtype(points[time_col]):
        points[time_col] = points[time_col].astype("float64")
    else:
        timestamps = pd.to_datetime(points[time_col], errors="coerce", utc=True)
        if timestamps.isna().any():
            raise ValueError(f"Cannot parse some values in time column {time_col!r}.")
        epoch = pd.Timestamp("1970-01-01", tz="UTC")
        points[time_col] = (timestamps - epoch).dt.total_seconds().astype("float64")
    return points


def _validate_points(points, source) -> None:
    if points.empty:
        raise ValueError(f"Point dataset is empty: {source!s}.")
    if points.crs is None:
        raise ValueError("Point dataset has no CRS; pass crs for CSV input or define it first.")
    if points.geometry.isna().any() or points.geometry.is_empty.any():
        raise ValueError("Point dataset contains null or empty geometries.")
    geometry_types = set(points.geometry.geom_type)
    if not geometry_types.issubset({"Point"}):
        raise ValueError(
            "Point dataset must contain only Point geometries; found "
            f"{', '.join(sorted(geometry_types))}."
        )


def load_points(
    points_path,
    return_bounding_box=False,
    time_col=None,
    *,
    x_col="lon",
    y_col="lat",
    crs="EPSG:4326",
    layer=None,
):
    """Load and validate point events through GeoPandas.

    Parameters
    ----------
    points_path : path-like
        CSV, GeoPackage, GeoJSON/JSON, Parquet, or GeoParquet source.
    return_bounding_box : bool, default False
        If true, return ``(points, max_y, min_y, max_x, min_x)`` for backward
        compatibility. New code should prefer ``points.total_bounds`` because
        the legacy tuple order is unusual.
    time_col : str, optional
        Time column to validate. Numeric values are converted to float64;
        other values are parsed as UTC datetimes and converted to Unix seconds.
    x_col, y_col : str, default "lon", "lat"
        CSV coordinate column names. Ignored for geospatial formats.
    crs : CRS-like, default "EPSG:4326"
        CRS assigned to CSV coordinates. Vector and Parquet inputs preserve
        their stored CRS.
    layer : str, optional
        GeoPackage layer name. Ignored for formats without layers.

    Returns
    -------
    geopandas.GeoDataFrame or tuple
        Validated, non-empty Point data, optionally followed by the legacy
        bounding-box values.

    Raises
    ------
    ImportError
        If the ``fastnsa[geo]`` dependencies are unavailable.
    ValueError
        If the format, coordinate columns, CRS, geometry, or time data is
        invalid.

    Notes
    -----
    This function does not project points. :class:`fastnsa.PointEvents` aligns
    them to the target network during mapping.
    """
    gpd, pd = _geo_dependencies()
    path = Path(points_path)
    extension = path.suffix.casefold()

    if extension == ".csv":
        frame = pd.read_csv(path)
        missing = [name for name in (x_col, y_col) if name not in frame.columns]
        if missing:
            raise ValueError(f"CSV is missing coordinate columns: {', '.join(missing)}.")
        coordinates = frame[[x_col, y_col]].apply(pd.to_numeric, errors="coerce")
        if coordinates.isna().any().any() or not np.isfinite(coordinates.to_numpy()).all():
            raise ValueError("CSV coordinate columns contain invalid or non-finite values.")
        points = gpd.GeoDataFrame(
            frame.copy(),
            geometry=gpd.points_from_xy(coordinates[x_col], coordinates[y_col]),
            crs=crs,
        )
    elif extension in {".gpkg", ".geojson", ".json"}:
        read_options = {"layer": layer} if layer is not None else {}
        points = gpd.read_file(path, **read_options)
    elif extension in {".parquet", ".geoparquet"}:
        points = gpd.read_parquet(path)
    else:
        raise ValueError(
            f"Unsupported point format {extension!r}; expected CSV, GPKG, GeoJSON, or Parquet."
        )

    points = _prepare_time_column(points, time_col, pd)
    _validate_points(points, path)
    logger.info(
        "Loaded point dataset",
        extra={
            "fastnsa_stage": "point_load",
            "fastnsa_source": str(path),
            "fastnsa_format": extension,
            "fastnsa_points": len(points),
            "fastnsa_crs": str(points.crs),
        },
    )

    if not return_bounding_box:
        return points
    min_x, min_y, max_x, max_y = (float(value) for value in points.total_bounds)
    return points, (max_y, min_y, max_x, min_x)
