import sys
import os
sys.path.insert(1, 'D://Code//Hogwarts//code//signaling network')
import pandas as pd
import HogwartsHallmarkAnalysis as hallmark
from constant import *
import networkx as nx
import numpy as np
import itertools
from statistics import mean
from HogwartsCanceType import *


# check if any pair of nodes in a given cancer_genes are connected in a graph
# def check_connected(cancer_genes, g):
#     flag = True
#     for u in cancer_genes:
#         for v in cancer_genes:
#             if not u == v:
#                 if not nx.has_path(g, u, v):
#                     flag = False
#     return flag


# construct a smallest connected network containing all the input cancer genes
def constructNetwork(cancer_genes,nodetype, targetset, g, cancer_name):  # construct cancer network
    temp1 = cancer_genes.copy()  # cancer_genes is the set of cancer_genes of this kind of cancer type
    cancer_genes1 = []
    for u in temp1:
        if u in g.nodes():
            cancer_genes1.append(u)  # remove the cancer_genes which are not in the network
    temp2 = targetset.copy()  # target set is the set of drug targets for the input cancer type
    targetset1 = []
    for v in temp2:
        if v in g.nodes():
            targetset1.append(v)  # remove the targets which are not in the network
    print('There are {} cancer genes and {} targets'.format(len(cancer_genes1), len(targetset1)))
    new_nodes = set()
    count = 0
    # if there is at least one path from a target u to cancer_gene v
    # if the shortest path length <= 5, add the paths with the length <=5 to the network
    # else, add all the paths from u to v to the network.
    for u in targetset1:
        for v in cancer_genes1:
            if nx.has_path(g, u, v):
                lp = nx.shortest_path_length(g, u, v)
                # print('length of the shortest path between {} and {} is {}'.format(u,v,lp))
                if lp <= 5:
                    temp = nx.all_simple_paths(g, u, v, 5)
                else:
                    temp = nx.all_simple_paths(g, u, v, lp)
                for p in temp:
                    new_nodes.update(p)
            # if there is no path from target u to cancer_gene v, add u and v to network only.
            else:
                new_nodes.update({u, v})

                # construct cancer network
    cancer_network = nx.subgraph(g, new_nodes)
    cancer_network_1 = cancer_network.copy()

    # compute statistics and pdist within the cancer network
    for i in [f'{output_path}/prune/{cancer_name}_{nodetype}/ppr',
              f'{output_path}/prune/{cancer_name}_{nodetype}/dppr',
              f'{output_path}/prune/{cancer_name}_{nodetype}/pdist',
              f'{output_path}/prune/{cancer_name}_{nodetype}/stats',
              f'{output_path}/prune/{cancer_name}_{nodetype}/distance']:
        if not os.path.exists(i):
            os.makedirs(i)
    # write the network to gexf file
    nx.write_gexf(cancer_network_1, f'{output_path}/prune/{cancer_name}_{nodetype}/{cancer_name}_{nodetype}.gexf')
    import HogwartsStat as stats
    # compute the topological features of the network
    stats.compute_network_stats(cancer_network_1, cancer_name, f'{output_path}/prune/{cancer_name}_{nodetype}/stats')
    import HogwartsPDist as pdist
    # compute the pdist of the network
    pdist.pdist_alpha(cancer_network_1,
                      f'{output_path}/prune/{cancer_name}_{nodetype}/ppr',
                      f'{output_path}/prune/{cancer_name}_{nodetype}/stats',
                      cancer_name, f'{output_path}/prune/{cancer_name}_{nodetype}/dppr',
                      f'{output_path}/prune/{cancer_name}_{nodetype}/pdist', iter=2)
    import HogwartsDistance as dist
    dist.compute_shortest_distance(cancer_network_1, f'{output_path}/prune/{cancer_name}_{nodetype}/distance')

    return cancer_network_1