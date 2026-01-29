import numpy as np
import pandas as pd
import  sys
sys.path.insert(0, '/home/yeswin/project/fast-nsa/python')

import fastnsa as nsa


def test_python_mapping_k_function_pipeline():
    net = nsa.Network.from_osm(
        place="Central Park, NY, USA"
    )

    assert net is not None
    assert hasattr(net, "_ptr")
    assert net._ptr is not None
    
    
