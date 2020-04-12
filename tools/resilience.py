# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 12:29:18 2020

@author: ANDRE BORGATO MORELLI
"""

import igraph as ig
import networkx as nx
import numpy as np
import operator
from tools.graph_operations import *

def get_efficiency(Gig, weight = None, track_progress=False):
    if weight != None:
        Gig = get_igraph(G, edge_weights = [weight])
    else:
        Gig = get_igraph(G)
    efficiency = {}
    for node in Gig.vs:
        shrt = Gig.shortest_paths_dijkstra(node,weights=weight)
        total_proximity=0
        for l in shrt[0]:
            if l == 0:
                continue
            total_proximity += 1/l
        total_proximity = total_proximity/len(Gig.vs)
        efficiency[int(node['osmid'])] = total_proximity
    return efficiency

def get_connectivity(Gig,track_progress=True):
    valid_paths = {}
    for node in Gig.vs:
        shrt = Gig.shortest_paths(node)
        valid = 0
        for l in shrt[0]:
            if l == 0:
                continue
            if l!=np.inf:
                valid+=1
        valid_paths[int(node['osmid'])] = valid
    return valid_paths

def get_targeted_efficiency_and_connectivity(Gig, target, target_true_label=0, weight = None, track_progress=True):
    efficiency = {}
    valid_paths = {}
    nodes = Gig.vs
    target_list = [node for node in nodes if node[target]==target_true_label]
    for node in nodes:
        shrt = Gig.shortest_paths_dijkstra(node,target=target_list,weights=weight)
        l = min(shrt)
        efficiency[int(node['osmid'])] = 1/l
        if 1/l==0:
            valid_paths[int(node['osmid'])] = 0
        else:
            valid_paths[int(node['osmid'])] = 1
    return efficiency, valid_paths

def get_efficiency_and_connectivity(Gig, weight = None, track_progress=True):
    efficiency = {}
    valid_paths = {}
    for node in Gig.vs:
        shrt = Gig.shortest_paths_dijkstra(node,weights=weight)
        total_proximity = 0
        valid = 0
        for l in shrt[0]:
            if l == 0:
                continue
            if l!=np.inf:
                valid+=1
            total_proximity += 1/l
        total_proximity = total_proximity/(len(Gig.vs)-1)
        efficiency[int(node['osmid'])] = total_proximity
        valid_paths[int(node['osmid'])] = valid
    return efficiency, valid_paths

def remove_nodes_by_attr(G, attr, remove, ascending=False):
    G_new = G.copy()
    lst = [(G.nodes[n][attr], n) for n in G.nodes]
    if ascending:
        lst= sorted(lst)
    else:
        lst= sorted(lst, reverse=True)
    remove = [n for m,n in lst[:int(remove*len(lst))]]
    for node in remove:
        G_new.remove_node(node)
    return G_new

def remove_edges_by_attr(G, attr, remove, ascending=False):
    G_new = G.copy()
    lst = [(G.edges[e][attr], e) for e in G.edges]
    if ascending:
        lst= sorted(lst)
    else:
        lst= sorted(lst, reverse=True)
    remove = [n for m,n in lst[:int(remove*len(lst))]]
    for edge in remove:
        G_new.remove_edge(*edge)
    return G_new

def remove_edges_random(G, remove, random_seed=None):
    random.seed(random_seed)
    G_new = G.copy()
    remove = random.sample(list(G.edges), int(remove*len(G.edges)))
    for edge in remove:
        G_new.remove_edge(*edge)
    return G_new

def remove_edge_random_stochastic(G, p, random_seed=None, copy=True):
    """
    Recieves a Graph and p (uniform probability of failure).
    
    Returns a degraded Graph
    """
    
    if copy:
        G_ = G.copy()
    else:
        G_ = G
    
    if random_seed is not None:
        random.seed(random_seed)
    for edge in list(G.edges):
        if random.random()<=p:
            G_.remove_edge(*edge)
    
    return(G_)


def remove_edge_stochastic_function(G, parameter, prob_func, prob_func_kws={}, random_seed=None, copy=True):
    """
    Recieves a Graph and p.
    p is function of a defined parameter

    Returns a degraded Graph
    """
    if random_seed is not None:
        random.seed(random_seed)
        
    if copy:
        G_ = G.copy()
    else:
        G_ = G
    
    
    lst = [G.edges[n][parameter] for n in G.edges]
    vmax, vmin = max(lst), min(lst)
    prob_func_kws['vmax'] = vmax
    prob_func_kws['vmin'] = vmin
    lst=None
    for edge in list(G.edges):
        p = prob_func(G.edges[edge][parameter], **prob_func_kws)
        if random.random()<=p:
            G_.remove_edge(*edge)
    return(G_)


# From here: set of stochastic functions for degradation
def _betweenness_degradation(bet, alpha, vmax=1, vmin=0):
    vmin, vmax = [n**.3 for n in [vmin, vmax]]
    norm = (bet**.3 - vmin)/(vmax-vmin)
    return alpha*norm