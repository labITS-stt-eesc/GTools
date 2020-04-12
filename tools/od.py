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

import pandas as pd
import numpy as np
import random

def get_edges_in_route(G, z_origin, z_destination, weight): #função para gerar as rotas de uma zona para outra
    origins = [n for n, zone in G.nodes(data="zone") if zone==z_origin]
    dests = [n for n, zone in G.nodes(data="zone") if zone==z_destination]
    if len(origins)==0 or len(dests)==0:
        return[]
    while True: # encontra uma rota válida
        o = random.choice(origins)
        d = random.choice(dests)
        
        if o==d:
            continue
        elif not nx.has_path(G, o, d):
            continue
        break
    
    path = nx.dijkstra_path(G, o, d, weight) #Nota p/ futuro - trocando para o iGraph isso fica mais rápido
    edges = [(path[i], path[i+1], 0) for i in range(len(path)-1)]
    return edges


def allocate_od_trips_from_bulk_data(G, od_df, origin_column, destination_column, trip_reason_column = None, 
                                     trip_mode_column = None, weight=None, reasons='all', modes='all', k=1):
    edge_counts = {}.fromkeys(G.edges,0)
    for row in tqdm_notebook(df.index):
        
        if reasons=='all': # Pula a iteração se o motivo não pertence aos fornecidos
            pass
        elif (trip_reason_column is not None) and (od_df[trip_reason_column][row] not in reasons):
            continue
        if modes=='all': # Pula a iteração se o modo não pertence aos fornecidos
            pass
        elif (trip_reason_column is not None) and (od_df[trip_mode_column][row] not in modes):
            continue
        for i in range(k):
            try:
                origin = int(od_df['Zorigem'][row]) 
                destination = int(od_df['Zdestino'][row])
            except ValueError: # alguns valores não podem ser convertidos p/ int como strings com "#N/A"
                break
            edges = get_edges_in_route(G, origin, destination, weight = weight)
            for edge in edges: #adiciona as arestas da rota na contagem total
                edge_counts[edge] += 1
    return edge_counts