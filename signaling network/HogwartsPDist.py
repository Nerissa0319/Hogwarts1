from constant import *
import pandas as pd
import os.path
import networkx as nx
from sklearn import preprocessing
import numpy as np
import json


# the algorithm to compute ppr, dppr and pdist
def compute_pdist(network, alpha, ppr_path, stat_path, filename, dppr_path, pdist_path):
    out_degree_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted out_degree.txt'),
        sep='\t', header=None) # read out-degree of the network
    node = list(out_degree_df.iloc[:, 0]) # nodes in the network
    value = list(out_degree_df.iloc[:, 1]) # out-degree values of the nodes
    out_degree_dict = {}
    for i in range(len(node)):
        out_degree_dict[node[i]] = value[i] # create an out-degree dictionary
    # create empty dataframes to store ppr, dppr and pdist
    ppr_df = pd.DataFrame(index=pd.Index(sorted(network.nodes())), columns=pd.Index(sorted(network.nodes())))
    dppr_df = pd.DataFrame(index=pd.Index(sorted(network.nodes())), columns=pd.Index(sorted(network.nodes())))
    pdist_df = pd.DataFrame(index=pd.Index(sorted(network.nodes())), columns=pd.Index(sorted(network.nodes())))

    # for each node in the network, compute ppr, dppr and pdist
    for src in range(network.number_of_nodes()): # src is the position of the current node in network
        nodename = sorted(network.nodes())[src]
        # the algorithm for computing ppr
        ppr_dict = {}
        numIter = 100 # iteration number
        A = nx.adjacency_matrix(network, nodelist=sorted(network.nodes())) # adjacency matrix of the network
        P = preprocessing.normalize(A, norm='l1', axis=1).T # Scale input vectors individually to unit norm (vector length)
        e = np.zeros(network.number_of_nodes())
        e[src] = 1.0
        ppr = e.copy() # initialize the ppr list
        for i in range(numIter): # compute ppr
            ppr = (1 - alpha) * P.dot(ppr) + e
        ppr_ls = alpha * ppr
        count = 0
        for j in sorted(network.nodes()):
            ppr_df.loc[nodename, j] = ppr_ls[count] # store ppr values
            dppr_df.loc[nodename, j] = ppr_ls[count] * out_degree_dict[str(nodename)] # store dppr values which is ppr * outdegree
            pdist_df.loc[nodename, j] = round((1 - np.log(dppr_df.loc[nodename, j] + (1e-05))), 4) # compute pdist
            count += 1
        print(src)
        src += 1
    # write ppr, dppr, pdist to txt files
    ppr_df.to_csv(f'{ppr_path}/ppr.txt', header=True, index=True, sep='\t')
    dppr_df.to_csv(f'{dppr_path}/dppr.txt', header=True, index=True, sep='\t')
    pdist_df.to_csv(f'{pdist_path}/pdist.txt', header=True, index=True, sep='\t')


# compute ppr, dppr and pdist for a specific network
# iter is the number of different alphas
# for example, when iter = 3, the function will compute pdist for alpha = 0.1, 0.2 and 0.3
# when iter = 9, the function will compute pdist for alpha = 0.1 to 0.9
def pdist_alpha(network, network_ppr_path, stat_network_path, network_name, network_dppr_path, network_pdist_path,
                iter=2):
    alpha = 0
    for i in range(iter):
        alpha = (i + 1) / 10 # when iter = 1, alpha = 0.1
        dppr_path = os.path.join(network_dppr_path, f'alpha = {alpha}')
        pdist_path = os.path.join(network_pdist_path, f'alpha = {alpha}')
        ppr_path = os.path.join(network_ppr_path, f'alpha = {alpha}')
        for path in [dppr_path,ppr_path,pdist_path]:
            if not os.path.exists(path):
                os.makedirs(path)
        compute_pdist(network, alpha, ppr_path, stat_network_path, network_name, dppr_path,
                      pdist_path)

