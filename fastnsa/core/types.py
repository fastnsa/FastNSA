from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class KFunctionResult:
    """
    Result of a network K-function computation.
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
    
    def __getitem__(self, key):
        if key not in self.__dataclass_fields__:
            raise KeyError(
                f"{key!r} is not a valid KFunctionResult field"
            )
        return getattr(self, key)
    
