from typing import Literal, Optional
import numpy as np

from fastnsa.core.types import KFunctionResult
from fastnsa.core.network import Network
from fastnsa.core.events import PointEvents
from fastnsa.bindings import capi


# todo
def network_k_function(
    *,
    network: Network,
    points: PointEvents,
    r_min: float = 0.0,
    r_max: float = 2000.0,
    r_step: int = 50,
    r_values: Optional[np.ndarray] = None,
    method: Literal["SPS", "NS", "STARS"] = "NS",
    simulations: int = 0, 
    alpha: float = 0.05,
    return_simulations: bool = False,
) -> KFunctionResult:
    """
    Compute the network K-function for a set of point events on a network.

    Parameters
    ----------
    points : PointEvents
        Event locations projected onto the network.
    network : Network
        Network on which distances are measured.
    r_min, r_max, r_step : float, optional
        Range and step size for network distances.
        Ignored if `r_values` is provided.
    r_values : np.ndarray, optional
        Explicit distance values at which K(r) is evaluated.
    method : {"SPS", "NS", "STARS"}
        Algorithm used to compute the K-function.
    simulations : int, optional
        Number of Monte Carlo simulations for envelope construction.
        Set to 0 to disable simulations.
    alpha : float, optional
        Significance level for Monte Carlo envelopes.
    return_simulations : bool, optional
        Whether to return all simulated K-functions.

    Returns
    -------
    KFunctionResult
        Result object containing observed K-function values and
        optional theoretical reference and simulation envelopes.
    """
    if not network or not points:
        raise ValueError("Both `network` and `points` must be provided and non-empty.")

    if simulations < 0:
        raise ValueError("`simulations` must be non-negative.")

    if not (0.0 < alpha < 1.0):
        raise ValueError("`alpha` must be in (0, 1).")
    
    if r_values is None:
        if r_step <= 0:
            raise ValueError("`r_step` must be positive.")
        if r_max <= r_min:
            raise ValueError("`r_max` must be larger than `r_min`.")

        r = np.arange(r_min, r_max + r_step, r_step, dtype=np.float64)
    else:
        r = np.ascontiguousarray(r_values, dtype=np.float64)
        if r.ndim != 1:
            raise ValueError("`r_values` must be a 1D array.")

    # todo: call the C++ backend to compute K-function values
    n_r = r.size
    n_points = points.n_points

    method_map = {
        "NS": 6,
        "SPS": capi.FASTNSA_K_SPS,
        "STARS": capi.FASTNSA_K_STARS,
    }
    method_flag = method_map.get(method, 0)


    k_obs = np.zeros(n_r, dtype=np.float64)

    if return_simulations and simulations > 0:
        k_sim = np.zeros((simulations, n_r), dtype=np.float64)
    else:
        k_sim = None

    # k_sim = np.zeros((simulations, n_r), dtype=np.float64) if simulations > 0 else None

    capi.lib.fastnsa_compute_k_function(
        network._ptr,
        points._ptr,
        r.ctypes.data_as(capi.P_DOUBLE),
        n_r,
        method_flag,
        simulations,
        k_obs.ctypes.data_as(capi.P_DOUBLE),
        k_sim.ctypes.data_as(capi.P_DOUBLE) if k_sim is not None else None
    )

    k_lo = k_hi = None

    if simulations > 0:
        # Compute envelopes
        lo_q = alpha / 2.0
        hi_q = 1.0 - alpha / 2.0
        k_lo = np.quantile(k_sim, lo_q, axis=0)
        k_hi = np.quantile(k_sim, hi_q, axis=0)

        if not return_simulations:
            k_sim = None
    
    k_theo = r.copy()

    return KFunctionResult(
            r=r,
            k_obs=k_obs,
            k_theo=k_theo,
            k_lo=k_lo,
            k_hi=k_hi,
            k_sim=k_sim if return_simulations else None,
            method=method,
            n_points=n_points,
            intensity=n_points / network.total_length,
            network_length=network.total_length
        )
