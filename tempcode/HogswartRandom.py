import os.path
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
def run_add(network):
    added_edge = add_random_edge(network, 0.2)
    random_stat_path = f'{random_path}\\added_edge\\statistics'
    random_ppr_path = f'{random_path}\\added_edge\\ppr'
    random_dppr_path = f'{random_path}\\added_edge\\dppr'
    random_pdist_path = f'{random_path}\\added_edge\\pdist'
    random_target_path = f'{random_path}\\added_edge\\target'
    random_stat_chart_path = f'{random_path}\\added_edge\\statistics\\charts'
    random_distance_path = f'{random_path}\\added_edge\\distance'
    for p in [random_stat_chart_path, random_stat_path, random_ppr_path, random_dppr_path,
              random_pdist_path, random_distance_path, random_target_path]:
        if not os.path.exists(p):
            os.makedirs(p)
    # compute statistics of the new network
    stat.compute_network_stats(added_edge, 'added_edge', f'{random_path}\\added_edge\\statistics')
    stat.plot_network_stats('added_edge', f'{random_path}\\added_edge\\statistics',
                            f'{random_path}\\added_edge\\statistics\\charts')
    # execute the target-related functions in the new network
    target.write_target_info('added_edge', random_stat_path, random_target_path)
    target.target_chart(random_target_path, 'added_edge')
    # compute the ppr, dppr, pdist of the new network
    pdist.pdist_alpha(added_edge, random_ppr_path, random_stat_path, 'added_edge', random_dppr_path,
                      random_pdist_path, 9)
    # compute the shortest distance for each pair of nodes in the new network
    for source in sorted(added_edge.nodes()):
        distance.compute_shortest_distance(added_edge, source, random_distance_path)
    target.target_pdist(random_pdist_path, random_target_path, added_edge)
    target.target_distance(random_distance_path, random_target_path, added_edge)
    target.plot_pdist(random_target_path)
    target.plot_distance(random_target_path)
    target.st_diff(random_target_path, 'added_edge')


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
def run_remove(network):
    removed_edge1 = remove_random_edge(network, 0.2)
    random_stat_path = f'{random_path}\\removed_edge1\\statistics'
    random_ppr_path = f'{random_path}\\removed_edge1\\ppr'
    random_dppr_path = f'{random_path}\\removed_edge1\\dppr'
    random_pdist_path = f'{random_path}\\removed_edge1\\pdist'
    random_target_path = f'{random_path}\\removed_edge1\\target'
    random_stat_chart_path = f'{random_path}\\removed_edge1\\statistics\\charts'
    random_distance_path = f'{random_path}\\removed_edge1\\distance'
    for p in [random_stat_chart_path, random_stat_path, random_ppr_path, random_dppr_path,
              random_pdist_path, random_distance_path, random_target_path]:
        if not os.path.exists(p):
            os.makedirs(p)
    stat.compute_network_stats(removed_edge1, 'removed_edge1', f'{random_path}\\removed_edge1\\statistics')
    stat.plot_network_stats('removed_edge1', f'{random_path}\\removed_edge1\\statistics',
                            f'{random_path}\\removed_edge1\\statistics\\charts')
    target.write_target_info('removed_edge1', random_stat_path, random_target_path)
    target.target_chart(random_target_path, 'removed_edge1')
    pdist.pdist_alpha(removed_edge1, random_ppr_path, random_stat_path, 'removed_edge1', random_dppr_path,
                      random_pdist_path, 9)
    for source in sorted(removed_edge1.nodes()):
        distance.compute_shortest_distance(removed_edge1, source, random_distance_path)
    target.target_pdist(random_pdist_path, random_target_path, removed_edge1)
    target.target_distance(random_distance_path, random_target_path, removed_edge1)
    target.plot_pdist(random_target_path)
    target.plot_distance(random_target_path)
    target.st_diff(random_target_path, 'removed_edge1')
