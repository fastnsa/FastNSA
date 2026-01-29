import ctypes
import os
import platform
from pathlib import Path

__all__ = ['lib']


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
    lib_name = _get_library_filename()
    this_dir = Path(__file__).resolve().parent
    lib_path = this_dir / lib_name

    if not lib_path.exists():
        raise FileNotFoundError(
            f"Cannot find {lib_name} in {this_dir}\n"
            f"Expected path: {lib_path}\n"
            f"Make sure the C++ core library is built and installed correctly."
        )

    # Windows needs explicit DLL search path handling (Python >= 3.8)
    if os.name == "nt":
        os.add_dll_directory(str(this_dir))

    try:
        return ctypes.CDLL(str(lib_path))
    except OSError as e:
        raise OSError(
            f"Failed to load shared library: {lib_path}\n"
            f"Original error: {e}"
        )


lib = _load_library()

