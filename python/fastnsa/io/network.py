"""Legacy OSMnx graph construction helpers."""


def build_osm_network(
    bbox=None,
    place=None,
    network_type="drive",
    crs=None,
    consolidate_tolerance=20
):
    """Build an undirected projected OSMnx graph for FastNSA conversion.

    This compatibility helper returns a graph rather than a FastNSA
    :class:`Network`. New code can usually call ``Network.from_osm`` directly.

    Parameters
    ----------
    bbox, place, network_type, crs, consolidate_tolerance
        Passed to the equivalent OSMnx download, projection, and intersection
        consolidation operations. ``place`` takes precedence over ``bbox``.

    Returns
    -------
    networkx.MultiGraph
        Undirected projected graph with edge geometry aligned to node direction.

    Raises
    ------
    ImportError
        If the ``fastnsa[osm]`` dependencies are unavailable.
    """
    try:
        import osmnx as ox
        from shapely.geometry import LineString
    except ImportError as exc:
        raise ImportError(
            'OSM network construction requires: pip install "fastnsa[osm]".'
        ) from exc

    if place:
        g = ox.graph_from_place(place, simplify=True, network_type=network_type)
    else:
        g = ox.graph_from_bbox(*bbox, simplify=True, network_type=network_type)

    g = ox.project_graph(g, to_crs=crs)
    g = ox.consolidate_intersections(
        g, tolerance=consolidate_tolerance, rebuild_graph=True
    )
    g = g.to_undirected()

    # Native offsets are measured from the graph's start node, so keep edge
    # geometry direction consistent with topology before conversion.
    x_dic = {node_id: data['x'] for node_id, data in g.nodes(data=True)}

    for u, v, data in g.edges(data=True):
        geom = data.get('geometry')

        x_coords = geom.xy[0]
        if abs(x_coords[0] - x_dic.get(u, x_coords[0])) > 1e-5:
            reversed_coords = list(geom.coords)[::-1]
            data['geometry'] = LineString(reversed_coords)

    return g
