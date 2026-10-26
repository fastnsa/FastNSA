"""FastNSA public API with optional geospatial dependencies loaded lazily."""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("fastnsa")
except PackageNotFoundError:
    __version__ = "0.2.0.dev0"


_EXPORTS = {
    "NativeNetwork": ("fastnsa.runtime", "NativeNetwork"),
    "NativePointSet": ("fastnsa.runtime", "NativePointSet"),
    "Network": ("fastnsa.core.network", "Network"),
    "PointEvents": ("fastnsa.core.events", "PointEvents"),
    "network_k_function": ("fastnsa.stats.k_function", "network_k_function"),
    "spatiotemporal_network_k_function": (
        "fastnsa.stats.k_function",
        "spatiotemporal_network_k_function",
    ),
    "network_kde": ("fastnsa.stats.network_kde", "network_kde"),
    "KFunctionResult": ("fastnsa.core.types", "KFunctionResult"),
    "NetworkKDEPointResult": ("fastnsa.core.types", "NetworkKDEPointResult"),
    "SpatiotemporalKFunctionResult": (
        "fastnsa.core.types",
        "SpatiotemporalKFunctionResult",
    ),
    "load_points": ("fastnsa.io.points", "load_points"),
    "load_result": ("fastnsa.io.results", "load_result"),
    "preprocess_events": ("fastnsa.pipeline.preprocess", "preprocess_events"),
    "PreprocessResult": ("fastnsa.pipeline.preprocess", "PreprocessResult"),
    "make_line_network": ("fastnsa.datasets", "make_line_network"),
    "make_line_events": ("fastnsa.datasets", "make_line_events"),
    "MappingMethod": ("fastnsa.enums", "MappingMethod"),
    "KFunctionMethod": ("fastnsa.enums", "KFunctionMethod"),
    "SpatiotemporalKFunctionMethod": ("fastnsa.enums", "SpatiotemporalKFunctionMethod"),
    "NKDVMethod": ("fastnsa.enums", "NKDVMethod"),
    "NKDVKernel": ("fastnsa.enums", "NKDVKernel"),
    "FastNSAError": ("fastnsa.exceptions", "FastNSAError"),
    "FastNSAValidationError": ("fastnsa.exceptions", "FastNSAValidationError"),
    "NativeBackendError": ("fastnsa.exceptions", "NativeBackendError"),
    "FastNSACancelledError": ("fastnsa.exceptions", "FastNSACancelledError"),
}


def __getattr__(name):
    target = _EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attribute_name = target
    value = getattr(import_module(module_name), attribute_name)
    globals()[name] = value
    return value

__all__ = [
    "Network",
    "PointEvents",
    "NativeNetwork",
    "NativePointSet",
    "network_k_function",
    "spatiotemporal_network_k_function",
    "network_kde",
    "KFunctionResult",
    "NetworkKDEPointResult",
    "SpatiotemporalKFunctionResult",
    "load_points",
    "load_result",
    "preprocess_events",
    "PreprocessResult",
    "make_line_network",
    "make_line_events",
    "MappingMethod",
    "KFunctionMethod",
    "SpatiotemporalKFunctionMethod",
    "NKDVMethod",
    "NKDVKernel",
    "FastNSAError",
    "FastNSAValidationError",
    "NativeBackendError",
    "FastNSACancelledError",
    "__version__",
]
