"""Stable string enums accepted by the public API."""

from enum import Enum


class _StringEnum(str, Enum):
    def __str__(self) -> str:
        return self.value


class MappingMethod(_StringEnum):
    """Nearest-edge mapping algorithms accepted by point constructors."""
    RTREE = "R-tree"
    COMP_APPROX = "comp_a"
    COMP_EXACT = "comp_e"


class KFunctionMethod(_StringEnum):
    """Native algorithms available for the network K-function."""
    SPS = "SPS"
    NS = "NS"
    STARS = "STARS"


class SpatiotemporalKFunctionMethod(_StringEnum):
    """Native algorithms available for spatiotemporal K-functions."""
    RQS = "RQS"
    SPS = "SPS"
    EAR = "EAR"
    MTS = "MTS"
    EARTH = "EARTH"


class NKDVMethod(_StringEnum):
    """Native network KDE/NKDV evaluation algorithms."""
    BASIC = "BASIC"
    SHARED_SPD = "SHARED_SPD"
    ADA = "ADA"
    IA = "IA"
    ADA_IA = "ADA_IA"
    LION = "LION"


class NKDVKernel(_StringEnum):
    """Validated kernel families exposed by the public NKDV API."""
    EPANECHNIKOV = "EPANECHNIKOV"


def enum_value(value):
    """Return an enum's string value while preserving regular strings."""
    return value.value if isinstance(value, Enum) else value
