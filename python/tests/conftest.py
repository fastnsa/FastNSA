import pytest
import networkx as nx
from shapely.geometry import LineString
from fastnsa import Network

@pytest.fixture
def simple_grid_network():
    """创建一个 100x100 的简单十字路口路网"""
    g = nx.MultiGraph(crs="EPSG:3857")
    # 节点坐标
    g.add_node(0, x=0, y=50)
    g.add_node(1, x=100, y=50)
    g.add_node(2, x=50, y=0)
    g.add_node(3, x=50, y=100)
    # 边
    g.add_edge(0, 1, weight=100, geometry=LineString([(0, 50), (100, 50)]))
    g.add_edge(2, 3, weight=100, geometry=LineString([(50, 0), (50, 100)]))
    return Network(g)