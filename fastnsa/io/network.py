
import osmnx as ox
from shapely.geometry import LineString


def build_osm_network(
    bbox=None,
    place=None,
    network_type="drive",
    crs=None,
    consolidate_tolerance=20
):
    if place:
        g = ox.graph_from_place(place, simplify=True, network_type=network_type)
    else:
        g = ox.graph_from_bbox(*bbox, simplify=True, network_type=network_type)

    g = ox.project_graph(g, to_crs=crs)
    g = ox.consolidate_intersections(
        g, tolerance=consolidate_tolerance, rebuild_graph=True
    )
    g = g.to_undirected()

    # Fix edge direction based on node coordinates 
    x_dic = {node_id: data['x'] for node_id, data in g.nodes(data=True)}

    for u, v, data in g.edges(data=True):
        geom = data.get('geometry')

        x_coords = geom.xy[0]
        if abs(x_coords[0] - x_dic.get(u, x_coords[0])) > 1e-5:
            reversed_coords = list(geom.coords)[::-1]
            data['geometry'] = LineString(reversed_coords)

    return g
