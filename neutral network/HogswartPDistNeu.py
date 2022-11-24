from ConstantNeu import *
import pandas as pd
import os.path
import networkx as nx
from sklearn import preprocessing
import numpy as np
import json


# the algorithm to compute ppr, dppr and pdist
def compute_pdist(network, alpha, ppr_path, stat_path, filename, dppr_path, pdist_path):
    # create ppr path for different alphas
    ppr_alpha_path = os.path.join(ppr_path, f'alpha = {alpha}')
    for paths in [ppr_alpha_path, dppr_path, pdist_path]:
        if not os.path.exists(paths):
            os.makedirs(paths)
    out_degree_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted out_degree.txt'),
        sep='\t', header=None)
    node = list(out_degree_df.iloc[:, 0])
    value = list(out_degree_df.iloc[:, 1])
    out_degree_dict = {}
    for i in range(len(node)):
        out_degree_dict[node[i]] = value[i]
    ppr_df = pd.DataFrame(index=pd.Index(sorted(network.nodes())), columns=pd.Index(sorted(network.nodes())))
    dppr_df = pd.DataFrame(index=pd.Index(sorted(network.nodes())), columns=pd.Index(sorted(network.nodes())))
    pdist_df = pd.DataFrame(index=pd.Index(sorted(network.nodes())), columns=pd.Index(sorted(network.nodes())))

    # for each node in the network, compute ppr, dppr and pdist
    for src in range(network.number_of_nodes()):
        nodename = sorted(network.nodes())[src]

        # the algorithm for computing ppr
        ppr_dict = {}
        numIter = 100
        A = nx.adjacency_matrix(network, nodelist=sorted(network.nodes()))
        P = preprocessing.normalize(A, norm='l1', axis=1).T
        e = np.zeros(network.number_of_nodes())
        e[src] = 1.0
        ppr = e.copy()
        for i in range(numIter):
            ppr = (1 - alpha) * P.dot(ppr) + e
        ppr_ls = alpha * ppr
        count = 0
        for j in sorted(network.nodes()):
            # ppr_dict[str((nodename, j))] = ppr_ls[count]
            ppr_df.loc[nodename, j] = ppr_ls[count]
            dppr_df.loc[nodename, j] = ppr_ls[count] * out_degree_dict[str(nodename)]
            pdist_df.loc[nodename, j] = round((1 - np.log(dppr_df.loc[nodename, j] + (1e-05))), 4)
            count += 1
        print(src)
        src += 1
    ppr_df.to_csv(f'{ppr_alpha_path}/ppr.txt', header=True, index=True, sep='\t')
    dppr_df.to_csv(f'{dppr_path}/dppr.txt', header=True, index=True, sep='\t')
    pdist_df.to_csv(f'{pdist_path}/pdist.txt', header=True, index=True, sep='\t')


# compute ppr, dppr and pdist for a specific network
# iter is the number of different alphas
# for example, when iter = 3, the function will compute pdist for alpha = 0.1, 0.2 and 0.3
# when iter = 9, the function will compute pdist for alpha = 0.1 to 0.9
def pdist_alpha(network, network_ppr_path, stat_network_path, network_name, network_dppr_path, network_pdist_path,
                iter=9):
    alpha = 0
    for i in range(iter):
        alpha = (i + 1) / 10
        dppr_path = os.path.join(network_dppr_path, f'alpha = {alpha}')
        pdist_path = os.path.join(network_pdist_path, f'alpha = {alpha}')
        compute_pdist(network, alpha, network_ppr_path, stat_network_path, network_name, dppr_path,
                      pdist_path)

