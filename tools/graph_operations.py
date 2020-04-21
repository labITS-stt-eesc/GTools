# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 10:48:09 2020

@author: ANDRE BORGATO MORELLI
"""

import igraph as ig
import networkx as nx
import numpy as np
import operator
import osmnx as ox
import math

def get_igraph(G, edge_weights=None, node_weights=None):
    """
    Transforms a NetworkX graph into an iGraph graph.

    Parameters
    ----------
    G : NetworkX DiGraph or Graph 
        The graph to be converted.
    edge_weights: list or string
        weights stored in edges in the original graph to be kept in new graph. 
        If None, no weight will be carried. See get_full_igraph to get all 
        weights and attributes into the graph.
    node_weights: list or string
        weights stored in nodes in the original graph to be kept in new graph. 
        If None, no weight will be carried. See get_full_igraph to get all 
        weights and attributes into the graph.

    Returns
    -------
    iGraph graph
    """
    if type(edge_weights) == str:
        edge_weights = [edge_weights]
    if type(node_weights) == str:
        node_weights = [node_weights]
    G = G.copy()
    G = nx.relabel.convert_node_labels_to_integers(G)
    Gig = ig.Graph(directed=True)
    Gig.add_vertices(list(G.nodes()))
    Gig.add_edges(list(G.edges()))
    if 'kind' not in G.graph.keys():
        G.graph['kind']=primal # if not specified, assume graph id primal
    if G.graph['kind']=='primal':
        Gig.vs['osmid'] = list(nx.get_node_attributes(G, 'osmid').values())
    elif G.graph['kind']=='dual':
        Gig.vs['osmid'] = list(G.edges)
    if edge_weights != None:
        for weight in edge_weights:
            Gig.es[weight] = [n for _,_,n in G.edges(data=weight)]
            
    if node_weights != None:
        for weight in node_weights:
            Gig.vs[weight] = [n for _,n in G.nodes(data=weight)]
    for v in Gig.vs:
        v['name'] = v['osmid']
    return Gig

def get_full_igraph(G):
    """
    Transforms a NetworkX graph into an iGraph graph keeping all possible info.

    Parameters
    ----------
    G : NetworkX DiGraph or Graph 
        The graph to be converted.
    
    Returns
    -------
    iGraph graph
    """
    all_edge_attrs = []
    all_node_attrs = []
    for edge in G.edges:
        for attr in G.edges[edge].keys():
            if attr in all_edge_attrs:
                continue
            all_edge_attrs.append(attr)
                
    for node in G.nodes:
        G.nodes[node]['osmid'] = str(G.nodes[node]['osmid'])
        for attr in G.nodes[node].keys():
            if attr in all_node_attrs:
                continue
            all_node_attrs.append(attr)       
    
    return get_igraph(G, all_edge_attrs, all_node_attrs)

def get_dual(G, node_to_edge='first'):
    """
    Transforms a NetworkX primal graph into a NetworkX dual graph.
    This function keeps edge properties in the nodes of dual graph
    and adds angular weights, keeping the x and y central positions 
    of edges so this graph can be plotted with OSMnx. 
    
    obs: This function is partially based on previous work by 
    Alceu Dal Bosco Jr.

    Parameters
    ----------
    G : NetworkX DiGraph or Graph 
        The graph to be converted.
    node_to_edge: 'first','last' or None
        gets attributes from nodes of the original graph to add to dual
        graph. 'first' gets the attribute from the first node of the edge, 
        'last' gets from the last. If None, ignore orginal graphs node
        attributes
    Returns
    -------
    NetworkX graph
    """
    Gd = nx.line_graph(G)
    Gd.graph = G.graph
    all_edge_attrs = []
    all_node_attrs = []
    for edge in G.edges:
        for attr in G.edges[edge].keys():
            if attr in all_edge_attrs:
                continue
            all_edge_attrs.append(attr)
    for node in G.nodes:
        G.nodes[node]['osmid'] = str(G.nodes[node]['osmid'])
        for attr in G.nodes[node].keys():
            if attr in all_node_attrs:
                continue
            all_node_attrs.append(attr)

    for node in Gd.nodes:
        if node_to_edge=='first':
            for attr in all_node_attrs:
                if attr in ['x', 'y']:
                    continue
                try:
                    Gd.nodes[node][attr] = G.nodes[node[0]][attr]
                except:
                    Gd.nodes[node][attr] = None
        for attr in all_edge_attrs:
            try:
                Gd.nodes[node][attr] = G.edges[node][attr]
            except:
                    Gd.nodes[node][attr] = None

    G_ = ox.project_graph(G)
    for e in Gd.nodes:
        Gd.nodes[e]['x'] = (G.nodes[e[0]]['x']+G.nodes[e[1]]['x'])/2
        Gd.nodes[e]['y'] = (G.nodes[e[0]]['y']+G.nodes[e[1]]['y'])/2
        Gd.nodes[e]['osmid'] = (G.nodes[e[0]]['osmid'], G.nodes[e[1]]['osmid'])
    for e1,e2,l in Gd.edges:
        ang_dif = _dif_angle((G_.nodes[e1[0]]['x'],G_.nodes[e1[0]]['y']), 
                            (G_.nodes[e1[1]]['x'],G_.nodes[e1[1]]['y']), 
                            (G_.nodes[e2[1]]['x'],G_.nodes[e2[1]]['y']))
        if ang_dif == 0:
            ang_dif += .0001 #avoid null values
        Gd.edges[(e1,e2,l)]['ang_dif'] = ang_dif
    Gd.graph['kind'] = 'dual'
    return Gd

def fast_betweenness(G, weight=None, kind = 'edge', norm=True):
    """
    Gets betweenness centrality. For relativelly large graphs, this func is 
    faster than networkx

    Parameters
    ----------
    G : NetworkX DiGraph or Graph 
        The graph to be considered.
    weight: string
        edge weights for shortest paths.
    kind: 'edge' or 'node'
        Betweenness for edges or nodes.
    norm: bool
        If True, returns norm betweenness (bet/((N-1)*(N-2))).

    Returns
    -------
    dict
    """
    if weight != None:
        Gig = get_igraph(G, edge_weights = weight)
    else:
        Gig = get_igraph(G)
    norm_val = len(G.nodes)*(len(G.nodes)-1)
    if kind=='edge':
        bet = Gig.edge_betweenness(weights=weight)
        if norm==True:
            return {e:b/norm_val for e,b in zip(G.edges,bet)}
        else:
            return {e:b for e,b in zip(G.edges,bet)}
    elif kind=='node':
        bet = Gig.betweenness(weights=weight)
        if norm==True:
            return {e:b/norm_val for e,b in zip(G.nodes,bet)}
        else:
            return {e:b for e,b in zip(G.nodes,bet)}
    
    
    
def fast_closeness(G, kind = 'edge', weight=None, norm=True):
    """
    Gets closeness centrality. For relativelly large graphs, this func is 
    faster than networkx

    Parameters
    ----------
    G : NetworkX DiGraph or Graph 
        The graph to be considered.
    weight: string
        edge weights for shortest paths.
    kind: 'edge' or 'node'
        closeness for edges or nodes. for edges, the closeness id the average
        of extreme nodes. 
    norm: bool
        If True, returns norm betweenness (clo*N).

    Returns
    -------
    dict
    """
    if weight != None:
        Gig = get_igraph(G, edge_weights = weight)
    else:
        Gig = get_igraph(G)
    
    clo = Gig.closeness(weights=weight)
    if kind == 'node':
        if norm:
            return {n:c for n,c in zip(G.nodes,clo)}
        else:
            return {n:c/len(G.nodes) for n,c in zip(G.nodes,clo)}
    elif kind == 'edge':
        
        n_clo = {n:c for n,c in zip(G.nodes,clo)}
        if norm:
            e_clo = {edge:(n_clo[edge[0]]+n_clo[edge[1]]) for edge in G.edges}
        else:
            e_clo = {edge:(n_clo[edge[0]]+n_clo[edge[1]])/len(G.nodes)/2 for edge in G.edges}
        return e_clo

def fast_shortest_paths_dijkstra(G, od_pairs, weights=None):
    pass

def gini(values, weights=None):
    """Calculate the Gini coefficient of a numpy array."""
    # addapted from https://github.com/oliviaguest/gini
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    if weights == None:
        array = np.array(values) #all values are treated equally
    else:
        if type(weights) != list:
            raise ValueError('weights must be a list')
        elif len(weights)!=len(values):
            raise ValueError('values and weights must be lists of same size')
        array = []
        for val, weight in zip(values, weights):
            array = array + [val]*weight
        array = np.array(array) 
    if np.amin(array) < 0:
        array -= np.amin(array) #values cannot be negative
    array += 0.0000001 #values cannot be 0
    array = np.sort(array) #values must be sorted
    index = np.arange(1,array.shape[0]+1) #index per array element
    n = array.shape[0]#number of array elements
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) #Gini coefficient

def concentration(values, upper_strata=1, weights=None):
    # upper_strata in percentage
    if weights != None:
        if type(weights) != list:
            raise ValueError('weights must be a list')
        elif len(weights)!=len(values):
            raise ValueError('values and weights must be lists of same size')
        array = []
        for val, weight in zip(values, weights):
            array = array + [val]*weight
        array = np.array(array)
        array
    else:
        array = np.array(values)
    array = np.sort(array)
    upper_concentration = array[-int(upper_strata/100*len(values)):].sum()
    lower_strata = ((array.cumsum() <= upper_concentration)*1).sum()
    lower_strata = (lower_strata+1)/len(array)*100
    return (lower_strata/upper_strata)
    

def get_attr_gini_coef(G, attribute, weight=None, kind = 'edge'):
    if kind == 'edge':
        attrs = [attr for n1, n2, attr in G.edges(data=attribute)]
        if weight != None:
            weight = [int(attr) for n1, n2, attr in G.edges(data=weight)]
    elif kind == 'node':
        attrs = [attr for n, attr in G.nodes(data=attribute)]
        if weight != None:
            weight = [int(attr) for n, attr in G.nodes(data=weight)]
    return gini(attrs, weights = weight)

def get_attr_concentration_coef(G, attribute, upper_strata=1, weight=None, kind = 'edge'):
    if kind == 'edge':
        attrs = [attr for n1, n2, attr in G.edges(data=attribute)]
        if weight != None:
            weight = [int(attr) for n1, n2, attr in G.edges(data=weight)]
    elif kind == 'node':
        attrs = [attr for n, attr in G.nodes(data=attribute)]
        if weight != None:
            weight = [int(attr) for n, attr in G.nodes(data=weight)]
    return concentration(attrs, weights = weight)
    

def _dif_angle(a, b, c):
    """
    Gets angle difference

    Parameters
    ----------
    a, b, c : tuples
        sequential node coords on projected graph.
    
    Returns
    -------
    float
    """
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return abs(ang) - 180 if abs(ang) > 180 else 180 - abs(ang)
