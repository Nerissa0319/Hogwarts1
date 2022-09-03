from ConstantNeu import *
import pandas as pd
import os.path
import networkx as nx
from sklearn import preprocessing
import numpy as np
import json


# the algorithm to compute ppr, dppr and pdist
def compute_pdist(graph, alpha, ppr_path, stat_path, filename, dppr_path, pdist_path):
    # create ppr path for different alphas
    ppr_alpha_path = os.path.join(ppr_path, f'alpha = {alpha}')
    for paths in [ppr_alpha_path, dppr_path, pdist_path]:
        if not os.path.exists(paths):
            os.makedirs(paths)

    # for each node in the graph, compute ppr, dppr and pdist
    for src in range(graph.number_of_nodes()):
        print("src: " + str(src))
        print(" ")
        nodename = sorted(graph.nodes())[src]

        # the algorithm for computing ppr
        ppr_dict = {}
        numIter = 100
        A = nx.adjacency_matrix(graph, nodelist=sorted(graph.nodes()))
        P = preprocessing.normalize(A, norm='l1', axis=1).T
        e = np.zeros(graph.number_of_nodes())
        e[src] = 1.0
        ppr = e.copy()
        for i in range(numIter):
            ppr = (1 - alpha) * P.dot(ppr) + e
        ppr_ls = alpha * ppr
        count = 0
        for j in sorted(graph.nodes()):
            ppr_dict[str((nodename, j))] = ppr_ls[count]
            count += 1

        # sort ppr
        sorted_ppr_dict = {k: v for k, v in sorted(ppr_dict.items(),
                                                   key=lambda item: item[1],
                                                   reverse=True)}
        # write ppr in files
        with open(os.path.join(ppr_alpha_path, '_'.join((nodename, 'ppr.txt'))), 'w') as f:
            json_str = json.dumps(sorted_ppr_dict, indent=0)
            f.write(json_str)
            f.write('\n')
        f.close()
        # compute dppr
        dppr_dict = {}
        out_degree_df = pd.read_csv(os.path.join(
            stat_path, f'{filename}_sorted out_degree.txt'),
            sep='\t', header=None)
        node = list(out_degree_df.iloc[:, 0])
        value = list(out_degree_df.iloc[:, 1])
        out_degree_dict = {}
        for i in range(len(node)):
            out_degree_dict[node[i]] = value[i]
        # dppr = out_degree * ppr
        for key, value in sorted_ppr_dict.items():
            dppr_dict[key] = out_degree_dict[str(nodename)] * value
        # write dppr into files
        with open(os.path.join(dppr_path, f'{nodename}_dppr.txt'), 'w') as f:
            json_str = json.dumps(dppr_dict, indent=0)
            f.write(json_str)
            f.write('\n')
        f.close()
        # compute pdist
        pdist_dict = {}
        # pdist = 1 - log(dppr + 0.00001) and is rounded to 4 decimal points
        for key, value in dppr_dict.items():
            pdist_dict[key] = round((1 - np.log(value + (1e-05))), 4)
        # write pdist to files
        with open(os.path.join(pdist_path, f'{nodename}_PDist.txt'), 'w') as f:
            json_str = json.dumps(pdist_dict, indent=0)
            f.write(json_str)
            f.write('\n')
        f.close()
        src += 1


# compute ppr, dppr and pdist for a specific network
# iter is the number of different alphas
# for example, when iter = 3, the function will compute pdist for alpha = 0.1, 0.2 and 0.3
# when iter = 9, the function will compute pdist for alpha = 0.1 to 0.9
def pdist_alpha(network,network_ppr_path,stat_network_path,network_name,network_dppr_path,network_pdist_path,iter=9):
    alpha = 0
    for i in range(iter):
        alpha = (i + 1) / 10
        dppr_path = os.path.join(network_dppr_path, f'alpha = {alpha}')
        pdist_path = os.path.join(network_pdist_path, f'alpha = {alpha}')
        compute_pdist(network, alpha, network_ppr_path, stat_network_path, network_name, dppr_path,
                      pdist_path)
