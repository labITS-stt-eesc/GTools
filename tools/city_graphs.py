# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 15:05:17 2020

@author: ANDRE BORGATO MORELLI
"""

import pickle
import igraph as ig
import networkx as nx
import osmnx as ox
import os
from shapely.geometry import Point
from  shapely.ops import cascaded_union
import pyproj
from functools import partial
from shapely.ops import transform, unary_union
import geopandas as gpd
from tools.graph_operations import *

def read_pickle(path):
    return pickle.load(open(path,'rb'))

def graph_from_circle(query, radius=1000, network_type='all_private', dual = False, return_igraph = False,
                      save_pickle=False, fname='graphs\\city_graph', osmnx_query_kws = {}):
    """
    Like fetching a graph with osmnx but with a circle instead of a square.

    Parameters
    ----------
    query : string or dict 
        Location to query Nominatim.
    radius: float
        Radius of the circle.
    network_type: string 
        see osmnx network types
    dual: bool
        if true converts graph to its dual form
    return_igraph: bool
        if true retruns the graph as iGraph
    save_pickle: bool
        if True saves file as a pickle in fname directory
    fname: string
        Directory to save.
    osmnx_query_kws: dict
        options for osmnx query. See osmnx properties at 
        https://osmnx.readthedocs.io/en/stable/osmnx.html.

    Returns
    -------
    iGraph or NetworkX graph
    """
    pt = ox.geocode(query)
    poly = circle_from_lat_lon(*pt, radius)
    G = ox.graph_from_polygon(poly, network_type=network_type, **osmnx_query_kws)
    G.graph['kind'] = 'primal'
    if dual:
        G = get_dual(G)
    if return_igraph:
        G=get_full_igraph(G)
    if save_pickle and return_igraph:
        _save_pickle_file(G,fname, extention='ig')
    elif save_pickle and not return_igraph:
        _save_pickle_file(G,fname, extention='nx')
    return G

def graph_from_address(query, distance=1000, network_type='all_private', dual = False, return_igraph = False,
                       save_pickle=False, fname='graphs\\city_graph', osmnx_query_kws = {}):
    """
    Like fetching a graph with osmnx but with additional functionality.

    Parameters
    ----------
    query : string or dict 
        Location to query Nominatim.
    distance: float
        distance of the sides of square from the center.
    network_type: string 
        see osmnx network types
    dual: bool
        if true converts graph to its dual form
    return_igraph: bool
        if true retruns the graph as iGraph
    save_pickle: bool
        if True saves file as a pickle in fname directory
    fname: string
        Directory to save.
    osmnx_query_kws: dict
        options for osmnx query. See osmnx properties at 
        https://osmnx.readthedocs.io/en/stable/osmnx.html.

    Returns
    -------
    iGraph or NetworkX graph
    """
    G = ox.graph_from_address(query, distance=distance, network_type=network_type, **osmnx_query_kws)
    G.graph['kind'] = 'primal'
    if dual:
        G = get_dual(G)
    if return_igraph:
        G=get_full_igraph(G)
    if save_pickle and return_igraph:
        _save_pickle_file(G,fname, extention='ig')
    elif save_pickle and not return_igraph:
        _save_pickle_file(G,fname, extention='nx')
    return G

def graph_from_place(query, network_type='all_private', dual = False, return_igraph = False,
                     save_pickle=False, fname='graphs\\city_graph', osmnx_query_kws = {}):
    """
    Like fetching a graph with osmnx but with additional functionality.

    Parameters
    ----------
    query : string or dict 
        Location to query Nominatim.
    network_type: string 
        see osmnx network types
    dual: bool
        if true converts graph to its dual form
    return_igraph: bool
        if true retruns the graph as iGraph
    save_pickle: bool
        if True saves file as a pickle in fname directory
    fname: string
        Directory to save.
    osmnx_query_kws: dict
        options for osmnx query. See osmnx properties at 
        https://osmnx.readthedocs.io/en/stable/osmnx.html.

    Returns
    -------
    iGraph or NetworkX graph
    """
    G = ox.graph_from_place(query, network_type=network_type, **osmnx_query_kws)
    G.graph['kind'] = 'primal'
    if dual:
        G = get_dual(G)
    if return_igraph:
        G=get_full_igraph(G)
    if save_pickle and return_igraph:
        _save_pickle_file(G,fname, extention='ig')
    elif save_pickle and not return_igraph:
        _save_pickle_file(G,fname, extention='nx')
    return G

def graph_from_point(point, buffer=0, buffer_type='circle', network_type='all_private', dual = False, return_igraph = False,
                     save_pickle=False, fname='graphs\\city_graph', osmnx_query_kws = {}):
    """
    Like fetching a graph with osmnx but with option of a circle or square.

    Parameters
    ----------
    point : float 
        Central point as (lat, lon).
    buffer: float
        Radius of the circle or half the side of the square.
    network_type: string 
        see osmnx network types
    dual: bool
        if true converts graph to its dual form
    return_igraph: bool
        if true retruns the graph as iGraph
    save_pickle: bool
        if True saves file as a pickle in fname directory
    fname: string
        Directory to save.
    osmnx_query_kws: dict
        options for osmnx query. See osmnx properties at 
        https://osmnx.readthedocs.io/en/stable/osmnx.html.

    Returns
    -------
    iGraph or NetworkX graph
    """
    if buffer_type=='circle':
        poly = circle_from_lat_lon(*point, buffer)
        G = ox.graph_from_polygon(poly, network_type=network_type, **osmnx_query_kws)
    elif buffer_type=='square':
        G = ox.graph_from_point(point, network_type=network_type, **osmnx_query_kws)
    G.graph['kind'] = 'primal'
    if dual:
        G = get_dual(G)
    if return_igraph:
        G=get_full_igraph(G)
    if save_pickle and return_igraph:
        _save_pickle_file(G,fname, extention='ig')
    elif save_pickle and not return_igraph:
        _save_pickle_file(G,fname, extention='nx')
    return G

def graph_from_traffic_zones(shp_directory, network_type='all_private', mark_traffic_zones_to_nodes = False, 
                             zone_column=None, convex_hull=False, dual = False, return_igraph = False,
                             save_pickle=False, fname='graphs\\city_graph', osmnx_query_kws={}):
    """
    Get graph from a polygon shapefile of traffic zones.

    Parameters
    ----------
    shp_directory : string
        Shapefile directory.
    network_type: string 
        See osmnx network types
    mark_traffic_zones_to_nodes: bool
        If True, add an attribute to nodes marking what zone they belong to.
    zone_column: bool
        The column of GeoDataFrame with the names of the zones. required to
        mark_traffic_zones_to_nodes.
    convex_hull: bool
        If True, don't use only traffic zones, but the convex hull of the shapes.
        This may correct imperfections in traffic zones, but may add nodes not
        contained in the area analysed.
    dual: bool
        If true converts graph to its dual form
    return_igraph: bool
        If true retruns the graph as iGraph
    save_pickle: bool
        If True saves file as a pickle in fname directory
    fname: string
        Directory to save.
    osmnx_query_kws: dict
        options for osmnx query. See osmnx properties at 
        https://osmnx.readthedocs.io/en/stable/osmnx.html.

    Returns
    -------
    iGraph or NetworkX graph
    """
    gdf = gpd.read_file('ZTs')
    gdf = gdf.to_crs(epsg=4326)
    polygons = [poly for poly in gdf.geometry]
    if convex_hull:
        boundary = gpd.GeoSeries(unary_union(polygons)).convex_hull[0]
    else:
        boundary = gpd.GeoSeries(unary_union(polygons))[0]
    G = ox.graph_from_polygon(boundary, network_type=network_type, **osmnx_query_kws)
    G.graph['kind'] = 'primal'
    if mark_traffic_zones_to_nodes:
        if zone_column == None:
            raise ValueError('Missing zone column name. If mark_traffic_zones_to_nodes is True, must pass zone_column attribute.')
        G = add_traffic_zones_to_nodes(G, gdf, zone_column)

    if dual:
        G = get_dual(G)
    if return_igraph:
        G=get_full_igraph(G)
    if save_pickle and return_igraph:
        _save_pickle_file(G,fname, extention='ig')
    elif save_pickle and not return_igraph:
        _save_pickle_file(G,fname, extention='nx')
    return G
    
def add_traffic_zones_to_nodes(G, gdf, zone_column):
    """
    Adds traffic zones to nodes of a graph.

    Parameters
    ----------
    G : NetworkX Graph or DiGraph
        Graph for info to be added.
    gdf: GeoDataFrame 
        GeoDataFrame of the traffic zones
    zone_column: string
        Name of the column.
    zone_column: bool
        The column of GeoDataFrame with the names of the zones.

    Returns
    -------
    NetworkX graph
    """
    gdf_g = ox.graph_to_gdfs(G)[0]
    gdf_g = gdf_g.to_crs(epsg=4326)
    join = gpd.sjoin(gdf,gdf_g, op='contains')
    
    for zone, osmid in zip(join[zone_column], join.osmid):
        G.nodes[osmid]['zone'] = zone
    for node in G.nodes:
        try:
            G.nodes[node]['zone']
        except:
            G.nodes[node]['zone']=None
    return G


def _save_pickle_file(G,fname, extention):
    """
    Saves graph into a pickle.

    Parameters
    ----------
    G : NetworkX Graph or DiGraph or iGraph Graph
        Graph to save.
    fname: string 
        directory + file name
    extention: string
        may be anything, if you think about it.

    Returns
    -------
    None
    """
    if len(fname.split('\\'))>len(fname.split('/')):
        if not os.path.exists('\\'.join(fname.split('\\')[:-1])):
            try:
                os.mkdir('\\'.join(fname.split('\\')[:-1]))
            except OSError:
                raise OSError("Creation of the directory %s failed" % path)
    elif len(fname.split('\\'))<len(fname.split('/')):
        if not os.path.exists('\\'.join(fname.split('/')[:-1])):
            try:
                os.mkdir('\\'.join(fname.split('/')[:-1]))
            except OSError:
                raise OSError("Creation of the directory %s failed" % path)
    pickle.dump(G,open(fname+'.'+extention, 'wb'))
    return None

def circle_from_lat_lon(lat, lon, radius):
    """
    Get a projected circle in lat, lon coordinates.

    Parameters
    ----------
    lat : float
        Central point latitude.
    lon : float
        Central point longitude.
    radius: float
        Radius of circle.

    Returns
    -------
    shapely Polygon
    """
    local_azimuthal_projection = f"+proj=aeqd +R=6371000 +units=m +lat_0={lat} +lon_0={lon}"
    point = Point(lon, lat)
    wgs84_to_aeqd = partial(
        pyproj.transform,
        pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'),
        pyproj.Proj(local_azimuthal_projection),
    )

    aeqd_to_wgs84 = partial(
        pyproj.transform,
        pyproj.Proj(local_azimuthal_projection),
        pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'),
    )

    point_transformed = transform(wgs84_to_aeqd, point)

    buffer = point_transformed.buffer(radius)
    buffer_wgs84 = transform(aeqd_to_wgs84, buffer)
    return buffer_wgs84