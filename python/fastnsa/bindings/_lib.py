"""Locate and load the platform-specific FastNSA shared library."""

import ctypes
import os
import platform
from pathlib import Path

__all__ = ['lib']

_dll_directory_handle = None


def _get_library_filename() -> str:
    """
    Return platform-specific shared library filename.
    """
    system = platform.system()

    if system == "Linux":
        return "libfastnsa.so"
    elif system == "Darwin":
        return "libfastnsa.dylib"
    elif system == "Windows":
        return "fastnsa.dll"
    else:
        raise RuntimeError(f"Unsupported platform: {system}")

def _load_library() -> ctypes.CDLL:
    """
    Load the fastnsa shared library from the same directory as this file.
    """
    global _dll_directory_handle

    lib_name = _get_library_filename()
    this_dir = Path(__file__).resolve().parent
    override = os.environ.get("FASTNSA_LIBRARY_PATH")
    lib_path = Path(override).expanduser().resolve() if override else this_dir / lib_name

    if not lib_path.exists():
        raise FileNotFoundError(
            f"Cannot find {lib_name} in {this_dir}\n"
            f"Expected path: {lib_path}\n"
            "Install a wheel matching this operating system and CPU architecture, "
            "or set FASTNSA_LIBRARY_PATH for a development build."
        )

    # Windows needs explicit DLL search path handling (Python >= 3.8)
    if os.name == "nt":
        _dll_directory_handle = os.add_dll_directory(str(lib_path.parent))

    try:
        return ctypes.CDLL(str(lib_path))
    except OSError as e:
        raise OSError(
            f"Failed to load shared library: {lib_path}\n"
            f"Platform: {platform.system()} {platform.machine()}\n"
            f"Python architecture: {ctypes.sizeof(ctypes.c_void_p) * 8}-bit\n"
            f"Original error: {e}"
        ) from e


lib = _load_library()

