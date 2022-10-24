import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constant import *
import random
import networkx as nx
import HogswartStat as stat
import HogswartDistance as distance
import HogswartPDist as pdist
import HogswartTarget as target


# add random edges to a network with a given percentage
def add_random_edge(G, percentage):
    no_of_edge = G.number_of_edges()
    no_to_add = no_of_edge * percentage
    nodes = list(G.nodes()).copy()
    count = 0
    while count < no_to_add:
        u = random.choice(nodes)
        v = random.choice(nodes)
        if not G.has_edge(u, v):
            G.add_edge(u, v)
            count += 1

    return G


# add the edges and then compute and plot statistics for the new network
def run_add(network,percentage,iter_no):
    added_edge = add_random_edge(network, percentage)
    nx.write_gexf(added_edge,f'{random_path}\\{iter_no}\\added_edge_{percentage}\\added_edge_{percentage}.gexf')
    random_stat_path = f'{random_path}\\{iter_no}\\added_edge_{percentage}\\statistics'
    random_ppr_path = f'{random_path}\\{iter_no}\\added_edge_{percentage}\\ppr'
    random_dppr_path = f'{random_path}\\{iter_no}\\added_edge_{percentage}\\dppr'
    random_pdist_path = f'{random_path}\\{iter_no}\\added_edge_{percentage}\\pdist'
    random_target_path = f'{random_path}\\{iter_no}\\added_edge_{percentage}\\target'
    random_stat_chart_path = f'{random_path}\\{iter_no}\\added_edge_{percentage}\\statistics\\charts'
    random_distance_path = f'{random_path}\\{iter_no}\\added_edge_{percentage}\\distance'
    for p in [random_stat_chart_path, random_stat_path, random_ppr_path, random_dppr_path,
              random_pdist_path, random_distance_path, random_target_path]:
        if not os.path.exists(p):
            os.makedirs(p)
    # compute statistics of the new network
    stat.compute_network_stats(added_edge, f'added_edge_{percentage}', random_stat_path)
    stat.plot_network_stats(f'added_edge_{percentage}', random_stat_path,
                            f'{random_stat_path}\\charts')
    # execute the target-related functions in the new network
    target.write_target_info(f'added_edge_{percentage}', random_stat_path, random_target_path)
    target.target_chart(random_target_path, f'added_edge_{percentage}')
    # compute the ppr, dppr, pdist of the new network
    pdist.pdist_alpha(added_edge, random_ppr_path, random_stat_path, f'added_edge_{percentage}', random_dppr_path,
                      random_pdist_path, 2)
    # compute the shortest distance for each pair of nodes in the new network

    distance.compute_shortest_distance(added_edge, random_distance_path)
    target.plot_pdist(random_pdist_path,random_target_path,2)
    target.plot_distance(random_distance_path,random_target_path)
    target.st_diff(random_target_path, f'added_edge_{percentage}',2)


# remove random edges from a network with a given percentage
def remove_random_edge(G, percentage):
    no_of_edge = G.number_of_edges()
    no_to_remove = no_of_edge * percentage
    nodes = list(G.nodes()).copy()
    count = 0
    while count < no_to_remove:
        u = random.choice(nodes)
        v = random.choice(nodes)
        if G.has_edge(u, v):
            G.remove_edge(u, v)
            count += 1

    # after removing edges, some nodes will have no neighbors, we remove such nodes also
    nodes_to_remove = []
    for u in G.nodes():
        if G.out_degree(u) == 0 and G.in_degree(u) == 0:
            nodes_to_remove.append(u)
    G.remove_nodes_from(nodes_to_remove)

    return G


# remove random edges and compute/plot the statistics for the new network
def run_remove(network,percentage,iter_no):
    removed_edge = remove_random_edge(network, percentage)
    nx.write_gexf(removed_edge, f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\removed_edge_{percentage}.gexf')
    random_stat_path = f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\statistics'
    random_ppr_path = f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\ppr'
    random_dppr_path = f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\dppr'
    random_pdist_path = f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\pdist'
    random_target_path = f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\target'
    random_stat_chart_path = f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\statistics\\charts'
    random_distance_path = f'{random_path}\\{iter_no}\\removed_edge_{percentage}\\distance'
    for p in [random_stat_chart_path, random_stat_path, random_ppr_path, random_dppr_path,
              random_pdist_path, random_distance_path, random_target_path]:
        if not os.path.exists(p):
            os.makedirs(p)
    stat.compute_network_stats(removed_edge, f'removed_edge_{percentage}', random_stat_path)
    stat.plot_network_stats(f'removed_edge_{percentage}', random_stat_path,
                            f'{random_stat_path}\\charts')
    target.write_target_info(f'removed_edge_{percentage}', random_stat_path, random_target_path)
    target.target_chart(random_target_path, f'removed_edge_{percentage}')
    pdist.pdist_alpha(removed_edge, random_ppr_path, random_stat_path, f'removed_edge_{percentage}', random_dppr_path,
                      random_pdist_path, 2)
    distance.compute_shortest_distance(removed_edge, random_distance_path)
    target.plot_pdist(random_pdist_path,random_target_path,2)
    target.plot_distance(random_distance_path,random_target_path)
    target.st_diff(random_target_path, f'removed_edge_{percentage}',2)