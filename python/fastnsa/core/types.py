"""Typed result containers, exports, summaries, and NPZ persistence."""

from dataclasses import dataclass
import csv
from pathlib import Path
from typing import Optional
import numpy as np


_RESULT_FORMAT_VERSION = 1


def _array_or_empty(value, *, dimensions: int = 1) -> np.ndarray:
    if value is not None:
        return np.asarray(value)
    return np.empty((0,) * dimensions, dtype=np.float64)


def _npz_path(filename) -> Path:
    path = Path(filename)
    return path if path.suffix.casefold() == ".npz" else Path(f"{path}.npz")


def _optional_array(archive, name: str) -> Optional[np.ndarray]:
    if name not in archive:
        return None
    value = np.asarray(archive[name])
    return value if value.size else None


def _scalar(archive, name: str, default=None):
    if name not in archive:
        return default
    value = np.asarray(archive[name]).item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return value


def _optional_text(archive, name: str) -> Optional[str]:
    value = _scalar(archive, name, "")
    return str(value) if value not in (None, "") else None


def _optional_number(archive, name: str, *, integer: bool = False):
    value = _scalar(archive, name, -1 if integer else np.nan)
    if integer:
        return None if int(value) < 0 else int(value)
    return None if not np.isfinite(value) else float(value)


def _validate_result_type(archive, expected: str) -> None:
    result_type = _optional_text(archive, "result_type")
    if result_type is not None and result_type != expected:
        raise ValueError(
            f"Expected a {expected!r} result archive, found {result_type!r}."
        )


def _to_dataframe(records):
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            'DataFrame export requires pandas: pip install "fastnsa[dataframe]".'
        ) from exc
    return pd.DataFrame.from_records(records)


def _save_records_csv(records, filename) -> str:
    path = Path(filename)
    rows = list(records)
    if not rows:
        raise ValueError("Cannot export an empty result to CSV.")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def _envelope_status(observed, lower, upper) -> Optional[np.ndarray]:
    if lower is None or upper is None:
        return None
    status = np.zeros(np.asarray(observed).shape, dtype=np.int8)
    status[np.asarray(observed) < np.asarray(lower)] = -1
    status[np.asarray(observed) > np.asarray(upper)] = 1
    return status


@dataclass
class KFunctionResult:
    """Network K-function curve and Monte Carlo inference metadata.

    ``r`` and ``k_obs`` always have shape ``(n_r,)``. Theoretical and envelope
    arrays are ``None`` when no simulations were requested. ``k_sim`` is kept
    only when ``return_simulations=True`` was passed to the analysis function.

    Instances support dictionary-style field access, tabular export, summaries,
    and versioned NPZ persistence. They do not own native resources.
    """

    # --- distance axis ---
    r: np.ndarray
    """
    1D array of network distances at which K(r) is evaluated.
    Shape: (n_r,)
    """

    # --- observed statistics ---
    k_obs: np.ndarray
    """
    Observed network K-function values.
    Shape: (n_r,)
    """

    # --- theoretical reference (CSR on network) ---
    k_theo: Optional[np.ndarray] = None
    """
    Theoretical K-function under complete spatial randomness (CSR).
    Shape: (n_r,)
    """

    # --- Monte Carlo envelopes ---
    k_lo: Optional[np.ndarray] = None
    """
    Lower envelope from Monte Carlo simulations.
    Shape: (n_r,)
    """

    k_hi: Optional[np.ndarray] = None
    """
    Upper envelope from Monte Carlo simulations.
    Shape: (n_r,)
    """

    # --- optional simulation details ---
    k_sim: Optional[np.ndarray] = None
    """
    Simulated K-functions.
    Shape: (n_sim, n_r)
    Only present if return_simulations=True.
    """

    # --- metadata ---
    method: Optional[str] = None
    """
    Algorithm used to compute the K-function
    (e.g., 'SPS', 'NS', 'STARS').
    """

    n_points: Optional[int] = None
    """
    Number of events used in the analysis.
    """

    intensity: Optional[float] = None
    """
    Estimated intensity (lambda) on the network.
    """

    network_length: Optional[float] = None
    """
    Total length of the network used for analysis.
    """
    
    random_seed: Optional[int] = None
    envelope_method: Optional[str] = None
    confidence_level: Optional[float] = None
    p_value: Optional[float] = None

    def __getitem__(self, key):
        """Return a named dataclass field, raising ``KeyError`` if unknown."""
        if key not in self.__dataclass_fields__:
            raise KeyError(
                f"{key!r} is not a valid KFunctionResult field"
            )
        return getattr(self, key)

    def to_records(self) -> list[dict[str, float]]:
        """Convert the curve to one plain dictionary per distance threshold."""
        records = []
        for index, distance in enumerate(self.r):
            record = {"r": float(distance), "k_obs": float(self.k_obs[index])}
            for name in ("k_theo", "k_lo", "k_hi"):
                values = getattr(self, name)
                if values is not None:
                    record[name] = float(values[index])
            records.append(record)
        return records

    def to_dataframe(self):
        """Return curve records as a pandas DataFrame.

        Requires the ``fastnsa[dataframe]`` optional dependency.
        """
        return _to_dataframe(self.to_records())

    def save_csv(self, filename) -> str:
        """Write threshold and available K-value columns to CSV.

        Returns the output path as a string. CSV does not preserve all result
        metadata or retained simulation arrays; use :meth:`save_npz` for a
        restorable FastNSA result.
        """
        return _save_records_csv(self.to_records(), filename)

    to_csv = save_csv

    @property
    def envelope_status(self) -> Optional[np.ndarray]:
        """Classify observed values relative to the Monte Carlo envelope.

        Returns int8 values ``-1`` below, ``0`` within, and ``1`` above the
        envelope, or ``None`` when envelope arrays are unavailable.
        """
        return _envelope_status(self.k_obs, self.k_lo, self.k_hi)

    def summary(self) -> dict:
        """Return compact analysis and envelope counts for reports or logs."""
        status = self.envelope_status
        return {
            "result_type": "network_k_function",
            "method": self.method,
            "n_points": self.n_points,
            "n_thresholds": int(len(self.r)),
            "below_envelope": int(np.count_nonzero(status == -1)) if status is not None else 0,
            "above_envelope": int(np.count_nonzero(status == 1)) if status is not None else 0,
            "has_envelope": status is not None,
            "random_seed": self.random_seed,
            "envelope_method": self.envelope_method,
            "confidence_level": self.confidence_level,
            "p_value": self.p_value,
        }

    @property
    def envelope_description(self) -> Optional[str]:
        """Human-readable confidence, envelope, and p-value description."""
        if self.envelope_method is None or self.confidence_level is None:
            return None
        return (
            f"{self.confidence_level:.1%} {self.envelope_method} Monte Carlo envelope; "
            f"simulation p-value={self.p_value:.4g}."
            if self.p_value is not None
            else f"{self.confidence_level:.1%} {self.envelope_method} Monte Carlo envelope."
        )

    def save_npz(self, filename) -> str:
        """Save arrays and metadata to a compressed, versioned NPZ archive.

        The ``.npz`` suffix is added when absent. Returns the resulting path.
        """
        path = _npz_path(filename)
        np.savez_compressed(
            path, format_version=_RESULT_FORMAT_VERSION,
            result_type="network_k_function",
            r=self.r, k_obs=self.k_obs, k_theo=_array_or_empty(self.k_theo),
            k_lo=_array_or_empty(self.k_lo), k_hi=_array_or_empty(self.k_hi),
            k_sim=_array_or_empty(self.k_sim, dimensions=2),
            method=self.method or "", random_seed=-1 if self.random_seed is None else self.random_seed,
            n_points=-1 if self.n_points is None else self.n_points,
            intensity=np.nan if self.intensity is None else self.intensity,
            network_length=np.nan if self.network_length is None else self.network_length,
            envelope_method=self.envelope_method or "",
            confidence_level=np.nan if self.confidence_level is None else self.confidence_level,
            p_value=np.nan if self.p_value is None else self.p_value,
        )
        return str(path)

    @classmethod
    def load_npz(cls, filename) -> "KFunctionResult":
        """Restore a network K-function result from its NPZ archive."""
        with np.load(filename, allow_pickle=False) as archive:
            _validate_result_type(archive, "network_k_function")
            return cls(
                r=np.asarray(archive["r"]),
                k_obs=np.asarray(archive["k_obs"]),
                k_theo=_optional_array(archive, "k_theo"),
                k_lo=_optional_array(archive, "k_lo"),
                k_hi=_optional_array(archive, "k_hi"),
                k_sim=_optional_array(archive, "k_sim"),
                method=_optional_text(archive, "method"),
                n_points=_optional_number(archive, "n_points", integer=True),
                intensity=_optional_number(archive, "intensity"),
                network_length=_optional_number(archive, "network_length"),
                random_seed=_optional_number(archive, "random_seed", integer=True),
                envelope_method=_optional_text(archive, "envelope_method"),
                confidence_level=_optional_number(archive, "confidence_level"),
                p_value=_optional_number(archive, "p_value"),
            )


@dataclass
class SpatiotemporalKFunctionResult:
    """Spatiotemporal K-function surfaces and inference metadata.

    Surface arrays use shape ``(n_s, n_t)`` and are indexed by spatial then
    temporal threshold. Retained simulations use shape ``(n_sim, n_s, n_t)``.
    ``time_unit``, ``time_origin``, and ``time_mode`` describe input times but
    do not alter numeric values.
    """
    s: np.ndarray
    t: np.ndarray
    k_obs: np.ndarray
    k_theo: Optional[np.ndarray] = None
    k_lo: Optional[np.ndarray] = None
    k_hi: Optional[np.ndarray] = None
    k_sim: Optional[np.ndarray] = None
    method: Optional[str] = None
    n_points: Optional[int] = None
    intensity: Optional[float] = None
    network_length: Optional[float] = None
    time_unit: Optional[str] = None
    time_origin: Optional[str] = None
    time_mode: Optional[str] = None
    random_seed: Optional[int] = None
    envelope_method: Optional[str] = None
    confidence_level: Optional[float] = None
    p_value: Optional[float] = None

    def __getitem__(self, key):
        """Return a named dataclass field, raising ``KeyError`` if unknown."""
        if key not in self.__dataclass_fields__:
            raise KeyError(
                f"{key!r} is not a valid SpatiotemporalKFunctionResult field"
            )
        return getattr(self, key)

    def to_records(self) -> list[dict[str, float]]:
        """Flatten surfaces to one record per spatial/temporal threshold pair."""
        records = []
        for s_index, spatial_threshold in enumerate(self.s):
            for t_index, temporal_threshold in enumerate(self.t):
                record = {
                    "s": float(spatial_threshold),
                    "t": float(temporal_threshold),
                    "k_obs": float(self.k_obs[s_index, t_index]),
                }
                for name in ("k_theo", "k_lo", "k_hi"):
                    values = getattr(self, name)
                    if values is not None:
                        record[name] = float(values[s_index, t_index])
                records.append(record)
        return records

    def to_dataframe(self):
        """Return flattened surface records as a pandas DataFrame."""
        return _to_dataframe(self.to_records())

    def save_csv(self, filename) -> str:
        """Write flattened thresholds and available surfaces to CSV."""
        return _save_records_csv(self.to_records(), filename)

    to_csv = save_csv

    def to_xarray(self, field: str = "k_obs"):
        """Export one available surface as an xarray DataArray.

        Parameters
        ----------
        field : {"k_obs", "k_theo", "k_lo", "k_hi"}, default "k_obs"
            Surface to export. Theoretical and envelope fields require Monte
            Carlo simulations.

        Returns
        -------
        xarray.DataArray
            Array with ``spatial_threshold`` and ``temporal_threshold`` axes
            and result metadata in ``attrs``.
        """
        if field not in {"k_obs", "k_theo", "k_lo", "k_hi"}:
            raise ValueError(f"Unsupported surface field: {field!r}.")
        values = getattr(self, field)
        if values is None:
            raise ValueError(f"Surface {field!r} is not available in this result.")
        try:
            import xarray as xr
        except ImportError as exc:
            raise ImportError(
                'xarray export requires: pip install "fastnsa[xarray]".'
            ) from exc
        data = xr.DataArray(
            values,
            coords={"spatial_threshold": self.s, "temporal_threshold": self.t},
            dims=("spatial_threshold", "temporal_threshold"),
            name=field,
        )
        data.attrs.update(
            method=self.method or "",
            time_unit=self.time_unit or "",
            time_origin=self.time_origin or "",
            time_mode=self.time_mode or "",
            envelope_method=self.envelope_method or "",
            confidence_level=self.confidence_level,
            p_value=self.p_value,
        )
        return data

    @property
    def envelope_status(self) -> Optional[np.ndarray]:
        """Return -1/0/1 classifications for every surface cell, if available."""
        return _envelope_status(self.k_obs, self.k_lo, self.k_hi)

    def summary(self) -> dict:
        """Return compact shape, time, and envelope metadata."""
        status = self.envelope_status
        return {
            "result_type": "spatiotemporal_network_k_function",
            "method": self.method,
            "n_points": self.n_points,
            "surface_shape": tuple(int(value) for value in self.k_obs.shape),
            "below_envelope": int(np.count_nonzero(status == -1)) if status is not None else 0,
            "above_envelope": int(np.count_nonzero(status == 1)) if status is not None else 0,
            "has_envelope": status is not None,
            "time_unit": self.time_unit,
            "time_origin": self.time_origin,
            "time_mode": self.time_mode,
            "random_seed": self.random_seed,
            "envelope_method": self.envelope_method,
            "confidence_level": self.confidence_level,
            "p_value": self.p_value,
        }

    @property
    def envelope_description(self) -> Optional[str]:
        """Human-readable confidence, envelope, and p-value description."""
        if self.envelope_method is None or self.confidence_level is None:
            return None
        return (
            f"{self.confidence_level:.1%} {self.envelope_method} Monte Carlo envelope; "
            f"simulation p-value={self.p_value:.4g}."
            if self.p_value is not None
            else f"{self.confidence_level:.1%} {self.envelope_method} Monte Carlo envelope."
        )

    def save_npz(self, filename) -> str:
        """Save all surfaces and metadata to a compressed versioned archive."""
        path = _npz_path(filename)
        np.savez_compressed(
            path, format_version=_RESULT_FORMAT_VERSION,
            result_type="spatiotemporal_network_k_function",
            s=self.s, t=self.t, k_obs=self.k_obs,
            k_theo=_array_or_empty(self.k_theo, dimensions=2),
            k_lo=_array_or_empty(self.k_lo, dimensions=2),
            k_hi=_array_or_empty(self.k_hi, dimensions=2),
            k_sim=_array_or_empty(self.k_sim, dimensions=3),
            method=self.method or "", time_unit=self.time_unit or "",
            time_origin=self.time_origin or "", time_mode=self.time_mode or "",
            random_seed=-1 if self.random_seed is None else self.random_seed,
            n_points=-1 if self.n_points is None else self.n_points,
            intensity=np.nan if self.intensity is None else self.intensity,
            network_length=np.nan if self.network_length is None else self.network_length,
            envelope_method=self.envelope_method or "",
            confidence_level=np.nan if self.confidence_level is None else self.confidence_level,
            p_value=np.nan if self.p_value is None else self.p_value,
        )
        return str(path)

    @classmethod
    def load_npz(cls, filename) -> "SpatiotemporalKFunctionResult":
        """Restore a spatiotemporal K-function result from NPZ."""
        with np.load(filename, allow_pickle=False) as archive:
            _validate_result_type(archive, "spatiotemporal_network_k_function")
            return cls(
                s=np.asarray(archive["s"]),
                t=np.asarray(archive["t"]),
                k_obs=np.asarray(archive["k_obs"]),
                k_theo=_optional_array(archive, "k_theo"),
                k_lo=_optional_array(archive, "k_lo"),
                k_hi=_optional_array(archive, "k_hi"),
                k_sim=_optional_array(archive, "k_sim"),
                method=_optional_text(archive, "method"),
                n_points=_optional_number(archive, "n_points", integer=True),
                intensity=_optional_number(archive, "intensity"),
                network_length=_optional_number(archive, "network_length"),
                time_unit=_optional_text(archive, "time_unit"),
                time_origin=_optional_text(archive, "time_origin"),
                time_mode=_optional_text(archive, "time_mode"),
                random_seed=_optional_number(archive, "random_seed", integer=True),
                envelope_method=_optional_text(archive, "envelope_method"),
                confidence_level=_optional_number(archive, "confidence_level"),
                p_value=_optional_number(archive, "p_value"),
            )


@dataclass
class NetworkKDEPointResult:
    """Network KDE values and lixel location metadata.

    Each array entry describes one lixel evaluation on ``edge_ids``. The
    ``offsets_n1`` and ``offsets_n2`` values are distances from the evaluation
    position to the stored start and end of its edge. Derived segment properties
    expose lixel start, end, and actual length in start-node coordinates.
    """
    edge_ids: np.ndarray
    offsets_n1: np.ndarray
    offsets_n2: np.ndarray
    density: np.ndarray
    method: Optional[str] = None
    kernel: Optional[str] = None
    bandwidth: Optional[float] = None
    lixel_length: Optional[float] = None
    n_points: Optional[int] = None
    network_length: Optional[float] = None
    normalized: bool = False
    normalization: Optional[str] = None
    weighted: bool = False
    total_weight: Optional[float] = None

    def __getitem__(self, key):
        """Return a named dataclass field, raising ``KeyError`` if unknown."""
        if key not in self.__dataclass_fields__:
            raise KeyError(
                f"{key!r} is not a valid NetworkKDEPointResult field"
            )
        return getattr(self, key)

    @property
    def segment_start_n1(self) -> np.ndarray:
        """Lixel start offsets measured from each stored edge start."""
        if not self.lixel_length or self.lixel_length <= 0:
            return self.offsets_n1.copy()
        return np.floor(self.offsets_n1 / self.lixel_length) * self.lixel_length

    @property
    def segment_end_n1(self) -> np.ndarray:
        """Lixel end offsets clipped to each edge's actual length."""
        edge_lengths = self.offsets_n1 + self.offsets_n2
        return np.minimum(edge_lengths, self.segment_start_n1 + float(self.lixel_length or 0.0))

    @property
    def segment_lengths(self) -> np.ndarray:
        """Actual non-negative lixel lengths, including short final segments."""
        return np.maximum(self.segment_end_n1 - self.segment_start_n1, 0.0)

    def save_npz(self, filename) -> str:
        """Save lixel arrays and analysis metadata to versioned NPZ."""
        path = _npz_path(filename)
        np.savez_compressed(
            path, format_version=_RESULT_FORMAT_VERSION,
            result_type="network_kde",
            edge_ids=self.edge_ids, offsets_n1=self.offsets_n1,
            offsets_n2=self.offsets_n2, density=self.density,
            segment_start_n1=self.segment_start_n1,
            segment_end_n1=self.segment_end_n1,
            method=self.method or "", kernel=self.kernel or "",
            bandwidth=self.bandwidth, lixel_length=self.lixel_length,
            normalized=self.normalized, normalization=self.normalization or "",
            n_points=-1 if self.n_points is None else self.n_points,
            network_length=np.nan if self.network_length is None else self.network_length,
            weighted=self.weighted,
            total_weight=np.nan if self.total_weight is None else self.total_weight,
        )
        return str(path)

    @classmethod
    def load_npz(cls, filename) -> "NetworkKDEPointResult":
        """Restore a Network KDE result from NPZ."""
        with np.load(filename, allow_pickle=False) as archive:
            _validate_result_type(archive, "network_kde")
            return cls(
                edge_ids=np.asarray(archive["edge_ids"]),
                offsets_n1=np.asarray(archive["offsets_n1"]),
                offsets_n2=np.asarray(archive["offsets_n2"]),
                density=np.asarray(archive["density"]),
                method=_optional_text(archive, "method"),
                kernel=_optional_text(archive, "kernel"),
                bandwidth=_optional_number(archive, "bandwidth"),
                lixel_length=_optional_number(archive, "lixel_length"),
                n_points=_optional_number(archive, "n_points", integer=True),
                network_length=_optional_number(archive, "network_length"),
                normalized=bool(_scalar(archive, "normalized", False)),
                normalization=_optional_text(archive, "normalization"),
                weighted=bool(_scalar(archive, "weighted", False)),
                total_weight=_optional_number(archive, "total_weight"),
            )

    def to_dataframe(self):
        """Return lixel records as a pandas DataFrame."""
        return _to_dataframe(self.to_records())

    def save_csv(self, filename) -> str:
        """Write lixel locations, segment lengths, and density values to CSV."""
        return _save_records_csv(self.to_records(), filename)

    to_csv = save_csv

    def summary(self) -> dict:
        """Return density range, integral, weighting, and lixel-count metadata."""
        finite_density = self.density[np.isfinite(self.density)]
        return {
            "result_type": "network_kde",
            "method": self.method,
            "kernel": self.kernel,
            "n_points": self.n_points,
            "n_lixels": int(len(self.edge_ids)),
            "density_min": float(np.min(finite_density)) if finite_density.size else None,
            "density_max": float(np.max(finite_density)) if finite_density.size else None,
            "density_integral": float(np.sum(self.density * self.segment_lengths)),
            "normalized": self.normalized,
            "normalization": self.normalization,
            "weighted": self.weighted,
            "total_weight": self.total_weight,
        }

    def to_records(self) -> list[dict[str, float | int]]:
        """Return one plain dictionary per output lixel segment."""
        return [
            {
                "edge_id": int(self.edge_ids[i]),
                "offset_n1": float(self.offsets_n1[i]),
                "offset_n2": float(self.offsets_n2[i]),
                "segment_start_n1": float(self.segment_start_n1[i]),
                "segment_end_n1": float(self.segment_end_n1[i]),
                "segment_length": float(self.segment_lengths[i]),
                "density": float(self.density[i]),
            }
            for i in range(len(self.edge_ids))
        ]
    
