import numpy as np
import pytest
import  sys
sys.path.insert(0, '/home/yeswin/project/fast-nsa/python')

import fastnsa as nsa


@pytest.mark.integration
def test_python_mapping_k_function_pipeline():
    """
    End-to-end test:
    Network.from_osm
      -> Python-based point-line mapping
      -> network K-function (C++ core)
    """

    # --------------------------------------------------
    # 1. 构建 Network（Python -> C++）
    # --------------------------------------------------
    net = nsa.Network.from_osm(
        place="Central Park, NY, USA"
    )

    assert net is not None
    assert hasattr(net, "_ptr")
    assert net._ptr is not None

    # --------------------------------------------------
    # 2. 构建事件点（Python mapping）
    # --------------------------------------------------
    # 用最小 mock 数据，避免真实 accidents.csv 带来的不确定性
    import geopandas as gpd
    from shapely.geometry import Point

    gdf = gpd.GeoDataFrame(
        {
            "id": [1, 2],
            "geometry": [
                Point(-73.9712, 40.7831),
                Point(-73.9680, 40.7815),
            ]
        },
        crs="EPSG:4326"
    )

    events = nsa.PointEvents(
        gdf,
        target_network=net
    )

    assert events is not None
    assert hasattr(events, "_ptr")
    assert events._ptr is not None

    # 基本 sanity check（Python mapping 已经完成）
    # assert len(events.edge_ids) == len(events.offsets)
    # assert len(events.edge_ids) > 0

    # --------------------------------------------------
    # 3. 调用 Network K-function
    # --------------------------------------------------
    r_values = [10.0, 50.0, 100.0]

    results = nsa.stats.k_function.network_k_function(
        network=net,
        points=events,
        r_values=r_values,
        method="NS",          # 或 FASTNSA_K_NS
    )

    # --------------------------------------------------
    # 4. 结果检查（不验证数值，只验证结构）
    # --------------------------------------------------
    assert results is not None

    # network_k_function returns a dataclass, access attributes directly
    assert hasattr(results, "r")
    assert hasattr(results, "k_obs")

    r = results.r
    k_obs = results.k_obs

    assert len(r) == len(r_values)
    assert len(k_obs) == len(r_values)

    # 不允许 NaN / inf（这是非常关键的稳定性检查）
    assert np.all(np.isfinite(k_obs))
