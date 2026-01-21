import ctypes
import os

def test_load_lib():
    here = os.path.dirname(__file__)
    lib_path = os.path.join(
        here, "..","fastnsa", "bindings", "libfastnsa.so"
    )

    lib = ctypes.CDLL(lib_path)

    assert hasattr(lib, "fastnsa_network_create")
    assert hasattr(lib, "fastnsa_network_destroy")
    print("Library loaded and functions found successfully.")
    
if __name__ == "__main__":
    test_load_lib()
