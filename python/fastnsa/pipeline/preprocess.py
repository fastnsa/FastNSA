"""Reusable event preprocessing workflows."""

from dataclasses import dataclass
import logging
from time import perf_counter
from typing import Optional

import numpy as np

from fastnsa.core.events import PointEvents


logger = logging.getLogger(__name__)


@dataclass
class PreprocessResult:
    """Prepared events plus mapping provenance and timing information.

    Attributes
    ----------
    events : PointEvents
        Open mapped events consumed by analysis functions.
    reused_mapping : bool
        Whether edge-ID and offset columns were reused.
    n_points : int
        Number of prepared events.
    timings : dict[str, float]
        Stage durations in seconds, including ``preprocess_total_s``.
    mapping_method : str
        Selected mapping method, or ``precomputed`` when columns were reused.
    """

    events: PointEvents
    reused_mapping: bool
    n_points: int
    timings: dict[str, float]
    mapping_method: str

    def close(self) -> None:
        """Close the owned :class:`PointEvents` object."""
        self.events.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def preprocess_events(
    points,
    network,
    *,
    edge_id_col: Optional[str] = None,
    offset_col: Optional[str] = None,
    time_col: Optional[str] = None,
    times=None,
    time_unit: Optional[str] = None,
    time_origin: Optional[str] = None,
    time_mode: Optional[str] = None,
    mapping_method: str = "R-tree",
    comp_eps: float = 10.0,
    n_threads: int = 0,
) -> PreprocessResult:
    """Map event geometries or reuse validated mapping columns.

    ``edge_id_col`` and ``offset_col`` must be provided together. When set,
    their values are passed to the native runtime without repeating nearest-
    edge matching; native validation still checks edge IDs and offsets.

    Parameters
    ----------
    points : geopandas.GeoDataFrame
        Input point rows. Geometry and CRS are required when mapping is not
        reused; row count and optional time fields are always used.
    network : Network or NativeNetwork
        Open target network.
    edge_id_col, offset_col : str, optional
        Existing mapping columns. Supply both or neither. Stored mappings are
        valid only for the exact same network and edge order.
    time_col : str, optional
        Numeric or datetime-compatible point attribute containing event time.
    times : array-like, optional
        Explicit numeric times, taking precedence over ``time_col``.
    time_unit, time_origin, time_mode : str, optional
        Descriptive metadata for downstream spatiotemporal results.
    mapping_method : {"R-tree", "comp_a", "comp_e"}, default "R-tree"
        Method used only when mapping columns are absent.
    comp_eps : float, default 10.0
        Approximation tolerance for ``comp_a``.
    n_threads : int, default 0
        Requested point-mapping worker count. Zero uses the bounded automatic
        policy.

    Returns
    -------
    PreprocessResult
        Context-manageable result containing the open events and timings.

    Raises
    ------
    ValueError
        If only one mapping column is specified, a named column is missing, or
        event/time values are invalid.

    See Also
    --------
    fastnsa.NativePointSet.map_xy : Matching algorithm reference.
    fastnsa.PointEvents.from_mapped_arrays : Direct mapping reuse.
    """
    started = perf_counter()
    if (edge_id_col is None) != (offset_col is None):
        raise ValueError("`edge_id_col` and `offset_col` must be provided together.")

    reused_mapping = edge_id_col is not None
    if reused_mapping:
        missing = [name for name in (edge_id_col, offset_col) if name not in points.columns]
        if missing:
            raise ValueError(f"Point data is missing mapping columns: {', '.join(missing)}.")
        edge_ids = np.asarray(points[edge_id_col])
        offsets = np.asarray(points[offset_col])
        events = PointEvents.from_mapped_arrays(
            points,
            network,
            edge_ids,
            offsets,
            time_col=time_col,
            times=times,
            time_unit=time_unit,
            time_origin=time_origin,
            time_mode=time_mode,
        )
    else:
        events = PointEvents(
            points,
            network,
            mapping_method=mapping_method,
            comp_eps=comp_eps,
            n_threads=n_threads,
            time_col=time_col,
            times=times,
            time_unit=time_unit,
            time_origin=time_origin,
            time_mode=time_mode,
        )

    timings = dict(events.timings)
    timings["preprocess_total_s"] = perf_counter() - started
    logger.info(
        "Prepared network events",
        extra={
            "fastnsa_stage": "event_preprocess",
            "fastnsa_points": len(points),
            "fastnsa_reused_mapping": reused_mapping,
            "fastnsa_duration_s": timings["preprocess_total_s"],
        },
    )
    return PreprocessResult(
        events=events,
        reused_mapping=reused_mapping,
        n_points=len(points),
        timings=timings,
        mapping_method="precomputed" if reused_mapping else mapping_method,
    )
