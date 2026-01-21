import geopandas as gpd
from shapely.geometry import Point
from fastnsa import PointEvents

def test_point_mapping(simple_grid_network):
    net = simple_grid_network
    # 创建一个刚好在边上的点
    df = gpd.GeoDataFrame({
        'geometry': [Point(10, 50)] # 在节点 0 和 1 之间
    }, crs="EPSG:3857")
    
    events = PointEvents(df, net)
    
    # 验证映射结果
    # 10, 50 应该落在横向那条边上，offset 为 10
    assert events.edge_ids[0] in [0, 1] # 取决于你 edge_id_map 的生成顺序
    assert events.offsets[0] == 10.0