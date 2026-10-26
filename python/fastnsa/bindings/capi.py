"""ctypes declarations for the stable FastNSA C ABI.

This module is internal. Public code should use :mod:`fastnsa.runtime` and the
analysis functions instead of calling these symbols directly.
"""

from ._lib import lib
import ctypes

REQUIRED_ABI_VERSION = 2

FASTNSA_K_RQS = 1
FASTNSA_K_SPS = 2
FASTNSA_K_EAR = 3
FASTNSA_K_MTS = 4
FASTNSA_K_EARTH = 5
FASTNSA_K_NS = 6
FASTNSA_K_STARS = 7

FASTNSA_PLNNS_RTREE = 0
FASTNSA_PLNNS_COMP_APPROX = 1
FASTNSA_PLNNS_COMP_EXACT = 3

FASTNSA_NKDV_BASIC = 1
FASTNSA_NKDV_SHARED_SPD = 2
FASTNSA_NKDV_ADA = 3
FASTNSA_NKDV_IA = 4
FASTNSA_NKDV_ADA_IA = 5
FASTNSA_NKDV_LION = 6

FASTNSA_NKDV_KERNEL_GAUSSIAN = 0
FASTNSA_NKDV_KERNEL_TRIANGULAR = 1
FASTNSA_NKDV_KERNEL_EPANECHNIKOV = 2
FASTNSA_NKDV_KERNEL_QUARTIC = 3

P_DOUBLE = ctypes.POINTER(ctypes.c_double)
P_INT = ctypes.POINTER(ctypes.c_int)
VOID_P = ctypes.c_void_p
P_UINT8 = ctypes.POINTER(ctypes.c_uint8)
PROGRESS_CALLBACK = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_void_p,
)


try:
    _abi_version = lib.fastnsa_abi_version
except AttributeError as exc:
    raise ImportError(
        "The FastNSA native library does not expose fastnsa_abi_version; "
        "reinstall the Python package and its matching native wheel."
    ) from exc

_abi_version.argtypes = []
_abi_version.restype = ctypes.c_int
NATIVE_ABI_VERSION = int(_abi_version())
if NATIVE_ABI_VERSION != REQUIRED_ABI_VERSION:
    raise ImportError(
        "Incompatible FastNSA native library ABI: "
        f"expected {REQUIRED_ABI_VERSION}, found {NATIVE_ABI_VERSION}. "
        "Reinstall the Python package and remove stale libfastnsa binaries."
    )

_REQUIRED_SYMBOLS = (
    "fastnsa_last_error",
    "fastnsa_network_create",
    "fastnsa_network_destroy",
    "fastnsa_network_set_wkbs_checked",
    "fastnsa_network_set_polylines",
    "fastnsa_pointset_create",
    "fastnsa_pointset_destroy",
    "fastnsa_pointset_set",
    "fastnsa_pointset_set_times",
    "fastnsa_pointset_size",
    "fastnsa_pointset_copy_edge_ids",
    "fastnsa_pointset_copy_offsets",
    "fastnsa_pointset_copy_times",
    "fastnsa_map_points",
    "fastnsa_pointset_dump_text",
    "fastnsa_pointset_dump_text_st",
    "fastnsa_compute_k_function_checked",
    "fastnsa_compute_k_function_with_progress",
    "fastnsa_compute_spatiotemporal_k_function_checked",
    "fastnsa_compute_spatiotemporal_k_function_with_progress",
    "fastnsa_nkdv_count_lixels",
    "fastnsa_compute_nkdv",
    "fastnsa_compute_nkdv_with_progress",
)
_missing_symbols = [name for name in _REQUIRED_SYMBOLS if not hasattr(lib, name)]
if _missing_symbols:
    raise ImportError(
        "The FastNSA native library is incomplete for ABI "
        f"{REQUIRED_ABI_VERSION}; missing: {', '.join(_missing_symbols)}."
    )

if hasattr(lib, "fastnsa_last_error"):
    lib.fastnsa_last_error.argtypes = []
    lib.fastnsa_last_error.restype = ctypes.c_char_p


def last_error() -> str:
    """Return the latest native error for the current thread."""
    if not hasattr(lib, "fastnsa_last_error"):
        return "Native FastNSA operation failed."
    message = lib.fastnsa_last_error()
    return message.decode("utf-8", errors="replace") if message else "Native FastNSA operation failed."


# Network lifetime and packed-geometry bindings. Keep these declarations in
# lockstep with cpp/bindings/fastnsa_capi.h; ctypes cannot validate ABI types.
lib.fastnsa_network_create.argtypes = [
    P_DOUBLE,
    ctypes.c_int,
    P_INT,
    P_INT,
    ctypes.c_int
]
lib.fastnsa_network_create.restype = VOID_P

lib.fastnsa_network_destroy.argtypes = [VOID_P]
lib.fastnsa_network_destroy.restype = None

# WKB offsets contain n_edges + 1 entries and delimit a packed byte buffer.
lib.fastnsa_network_set_wkbs.argtypes = [VOID_P, P_UINT8, P_INT, ctypes.c_int]
lib.fastnsa_network_set_wkbs.restype = None

lib.fastnsa_network_set_polylines.argtypes = [VOID_P, P_DOUBLE, P_INT, ctypes.c_int]
lib.fastnsa_network_set_polylines.restype = ctypes.c_int

HAS_CHECKED_WKB_SETTER = hasattr(lib, "fastnsa_network_set_wkbs_checked")
if HAS_CHECKED_WKB_SETTER:
    lib.fastnsa_network_set_wkbs_checked.argtypes = list(
        lib.fastnsa_network_set_wkbs.argtypes
    )
    lib.fastnsa_network_set_wkbs_checked.restype = ctypes.c_int


# Point sets borrow their network pointer; Python wrappers keep the network
# alive and require point sets to be destroyed first.
lib.fastnsa_pointset_create.argtypes = [VOID_P]
lib.fastnsa_pointset_create.restype = VOID_P

lib.fastnsa_pointset_destroy.argtypes = [VOID_P]
lib.fastnsa_pointset_destroy.restype = None

# Populate native events from contiguous edge-ID, offset, and optional time
# arrays. Functions return zero on success unless their declaration says void.
lib.fastnsa_pointset_set.argtypes = [VOID_P, P_INT, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_set.restype = ctypes.c_int
lib.fastnsa_pointset_set_st.argtypes = [VOID_P, P_INT, P_DOUBLE, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_set_st.restype = ctypes.c_int
lib.fastnsa_pointset_set_times.argtypes = [VOID_P, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_set_times.restype = ctypes.c_int
lib.fastnsa_pointset_size.argtypes = [VOID_P]
lib.fastnsa_pointset_size.restype = ctypes.c_int
lib.fastnsa_pointset_copy_edge_ids.argtypes = [VOID_P, P_INT, ctypes.c_int]
lib.fastnsa_pointset_copy_edge_ids.restype = ctypes.c_int
lib.fastnsa_pointset_copy_offsets.argtypes = [VOID_P, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_copy_offsets.restype = ctypes.c_int
lib.fastnsa_pointset_copy_times.argtypes = [VOID_P, P_DOUBLE, ctypes.c_int]
lib.fastnsa_pointset_copy_times.restype = ctypes.c_int


lib.fastnsa_map_points.argtypes = [
    ctypes.c_void_p,  # network pointer
    P_DOUBLE,
    P_DOUBLE,
    ctypes.c_int,  # point count
    ctypes.c_int,  # mapping method enum
    ctypes.c_double,  # approximation epsilon
    ctypes.c_int,  # requested worker count
]
lib.fastnsa_map_points.restype = VOID_P

lib.fastnsa_pointset_dump_text.restype = ctypes.c_int
lib.fastnsa_pointset_dump_text.argtypes = [VOID_P, ctypes.c_char_p]
lib.fastnsa_pointset_dump_text_st.restype = ctypes.c_int
lib.fastnsa_pointset_dump_text_st.argtypes = [VOID_P, ctypes.c_char_p]



# Base K-function entry points predate checked/progress/seeded variants and
# return void. Public wrappers prefer the checked variants when available.
lib.fastnsa_compute_k_function.argtypes = [
    VOID_P,  # network pointer
    VOID_P,  # point-set pointer
    P_DOUBLE,  # distance thresholds
    ctypes.c_int,  # threshold count
    ctypes.c_int,  # method enum
    ctypes.c_int,  # simulation count
    P_DOUBLE,  # observed output
    P_DOUBLE,  # optional simulation output
    ctypes.c_double,  # approximate-method tolerance
]
lib.fastnsa_compute_k_function.restype = None

HAS_CHECKED_K_FUNCTION = hasattr(lib, "fastnsa_compute_k_function_checked")
if HAS_CHECKED_K_FUNCTION:
    lib.fastnsa_compute_k_function_checked.argtypes = list(lib.fastnsa_compute_k_function.argtypes)
    lib.fastnsa_compute_k_function_checked.restype = ctypes.c_int

HAS_PROGRESS_K_FUNCTION = hasattr(lib, "fastnsa_compute_k_function_with_progress")
if HAS_PROGRESS_K_FUNCTION:
    lib.fastnsa_compute_k_function_with_progress.argtypes = list(
        lib.fastnsa_compute_k_function.argtypes
    ) + [PROGRESS_CALLBACK, VOID_P]
    lib.fastnsa_compute_k_function_with_progress.restype = ctypes.c_int

HAS_SEEDED_K_FUNCTION = hasattr(lib, "fastnsa_compute_k_function_seeded")
if HAS_SEEDED_K_FUNCTION:
    lib.fastnsa_compute_k_function_seeded.argtypes = list(
        lib.fastnsa_compute_k_function.argtypes
    ) + [ctypes.c_uint32, PROGRESS_CALLBACK, VOID_P]
    lib.fastnsa_compute_k_function_seeded.restype = ctypes.c_int

lib.fastnsa_compute_spatiotemporal_k_function.argtypes = [
    VOID_P,   # network_ptr
    VOID_P,   # point_set_ptr
    P_DOUBLE, # s_values
    ctypes.c_int, # n_s
    P_DOUBLE, # t_values
    ctypes.c_int, # n_t
    ctypes.c_int, # method
    ctypes.c_int, # n_simulations
    P_DOUBLE, # out k_obs
    P_DOUBLE, # out k_sim
]
lib.fastnsa_compute_spatiotemporal_k_function.restype = None

HAS_CHECKED_ST_K_FUNCTION = hasattr(
    lib,
    "fastnsa_compute_spatiotemporal_k_function_checked",
)
if HAS_CHECKED_ST_K_FUNCTION:
    lib.fastnsa_compute_spatiotemporal_k_function_checked.argtypes = list(
        lib.fastnsa_compute_spatiotemporal_k_function.argtypes
    )
    lib.fastnsa_compute_spatiotemporal_k_function_checked.restype = ctypes.c_int

HAS_PROGRESS_ST_K_FUNCTION = hasattr(
    lib, "fastnsa_compute_spatiotemporal_k_function_with_progress"
)
if HAS_PROGRESS_ST_K_FUNCTION:
    lib.fastnsa_compute_spatiotemporal_k_function_with_progress.argtypes = list(
        lib.fastnsa_compute_spatiotemporal_k_function.argtypes
    ) + [PROGRESS_CALLBACK, VOID_P]
    lib.fastnsa_compute_spatiotemporal_k_function_with_progress.restype = ctypes.c_int

HAS_SEEDED_ST_K_FUNCTION = hasattr(
    lib, "fastnsa_compute_spatiotemporal_k_function_seeded"
)
if HAS_SEEDED_ST_K_FUNCTION:
    lib.fastnsa_compute_spatiotemporal_k_function_seeded.argtypes = list(
        lib.fastnsa_compute_spatiotemporal_k_function.argtypes
    ) + [ctypes.c_uint32, PROGRESS_CALLBACK, VOID_P]
    lib.fastnsa_compute_spatiotemporal_k_function_seeded.restype = ctypes.c_int

lib.fastnsa_nkdv_count_lixels.argtypes = [
    VOID_P,
    ctypes.c_double,
]
lib.fastnsa_nkdv_count_lixels.restype = ctypes.c_int

lib.fastnsa_compute_nkdv.argtypes = [
    VOID_P,
    VOID_P,
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_int,
    P_INT,
    P_DOUBLE,
    P_DOUBLE,
    P_DOUBLE,
]
lib.fastnsa_compute_nkdv.restype = ctypes.c_int

HAS_PROGRESS_NKDV = hasattr(lib, "fastnsa_compute_nkdv_with_progress")
if HAS_PROGRESS_NKDV:
    lib.fastnsa_compute_nkdv_with_progress.argtypes = list(
        lib.fastnsa_compute_nkdv.argtypes
    ) + [PROGRESS_CALLBACK, VOID_P]
    lib.fastnsa_compute_nkdv_with_progress.restype = ctypes.c_int

HAS_WEIGHTED_NKDV = hasattr(lib, "fastnsa_compute_nkdv_weighted")
if HAS_WEIGHTED_NKDV:
    lib.fastnsa_compute_nkdv_weighted.argtypes = [
        VOID_P, VOID_P, P_DOUBLE, ctypes.c_int,
        ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double,
        ctypes.c_int, P_INT, P_DOUBLE, P_DOUBLE, P_DOUBLE,
        PROGRESS_CALLBACK, VOID_P,
    ]
    lib.fastnsa_compute_nkdv_weighted.restype = ctypes.c_int


def make_progress_callback(callback):
    """Wrap a Python progress/cancel callable for the native ABI.

    The Python callable receives ``(stage, completed, total)`` and should return
    a truthy value to request cancellation. Exceptions cannot cross the C ABI;
    an exception raised by the callable is therefore converted to a cancellation
    request. The returned ctypes callback must remain alive for the complete
    native call.
    """
    if callback is None:
        return None

    def bridge(stage, completed, total, _user_data):
        try:
            stage_text = stage.decode("utf-8", errors="replace") if stage else "backend compute"
            return int(bool(callback(stage_text, int(completed), int(total))))
        except Exception:
            # Never allow Python exceptions to cross the C callback boundary.
            return 1

    return PROGRESS_CALLBACK(bridge)
