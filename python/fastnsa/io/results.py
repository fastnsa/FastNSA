"""Persistence helpers for FastNSA analysis results."""

from pathlib import Path

import numpy as np

from fastnsa.core.types import (
    KFunctionResult,
    NetworkKDEPointResult,
    SpatiotemporalKFunctionResult,
)


_RESULT_TYPES = {
    "network_k_function": KFunctionResult,
    "spatiotemporal_network_k_function": SpatiotemporalKFunctionResult,
    "network_kde": NetworkKDEPointResult,
}


def load_result(filename):
    """Load a typed FastNSA result from a versioned NPZ archive.

    Parameters
    ----------
    filename : path-like
        Archive created by a result object's ``save_npz`` method.

    Returns
    -------
    KFunctionResult, SpatiotemporalKFunctionResult, or NetworkKDEPointResult
        Result class selected from the archive's ``result_type`` metadata.

    Raises
    ------
    ValueError
        If the file cannot be read safely, lacks type metadata, or contains an
        unsupported result type.

    Notes
    -----
    Archives are loaded with ``allow_pickle=False``. CSV exports cannot be
    restored with this function because they do not retain full metadata.
    """
    path = Path(filename)
    try:
        with np.load(path, allow_pickle=False) as archive:
            if "result_type" not in archive:
                raise ValueError(
                    "The archive has no result_type metadata; use the specific "
                    "result class to load legacy archives."
                )
            result_type = str(np.asarray(archive["result_type"]).item())
    except (OSError, ValueError) as exc:
        if isinstance(exc, ValueError) and "result_type" in str(exc):
            raise
        raise ValueError(f"Could not read FastNSA result archive {path!s}: {exc}") from exc

    result_class = _RESULT_TYPES.get(result_type)
    if result_class is None:
        raise ValueError(f"Unsupported FastNSA result type: {result_type!r}.")
    return result_class.load_npz(path)
