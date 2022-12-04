import json
import operator
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats as st
from ConstantNeu import *
import networkx as nx
import HogwartsOtherFuncNeu as others


def add_legend():
    # add legend
    handles, labels = plt.gca().get_legend_handles_labels()
    i = 1
    while i < len(labels):
        if labels[i] in labels[:i]:
            del (labels[i])
            del (handles[i])
        else:
            i += 1
    plt.legend(handles, labels)
    return


def compute_network_stats(graph, filename, save_to):
    f = open(os.path.join(save_to, '.'.join((filename, 'txt'))), 'w',
             encoding='utf-8')
    f.write(str(nx.info(graph)))
    f.write('\n')

    # density
    density = nx.density(graph)
    f.write('Density of the network: ')
    f.write(str(round(density, 4)))
    f.write('\n')
    print('density computed')
    # average clustering coefficient
    a_clust = nx.average_clustering(graph)
    f.write("Average clustering coefficient of the network: ")
    f.write(str(round(a_clust, 4)))
    f.write("\n")
    print('average clustering coefficient computed')
    # number of strongly connected components and giant comoponent
    components = nx.strongly_connected_components(graph)
    components_no = nx.number_strongly_connected_components(graph)
    f.write("Number of strongly connected components: ")
    f.write(str(components_no))
    f.write("\n")
    largest_component = max(components, key=len)  # largest connected components
    f.write("Size of the largest strongly connected component: ")
    f.write(str(len(largest_component)))
    f.write("\n")
    # number of weakly connected components
    w_components = nx.weakly_connected_components(graph)
    w_components_no = nx.number_weakly_connected_components(graph)
    f.write("Number of weakly connected components: ")
    f.write(str(w_components_no))
    f.write("\n")
    w_largest_component = max(w_components, key=len)  # largest connected components
    f.write("Size of the largest weakly connected component: ")
    f.write(str(len(w_largest_component)))
    f.write("\n")
    print('component computed')
    # compute degree centrality
    degree_centrality_dict = dict(nx.degree_centrality(graph))
    sorted_degree_centrality = others.rank(degree_centrality_dict)
    f.write("\n")
    f.write("Top 10 nodes by degree centrality: ")
    f.write("\n")
    for index, row in sorted_degree_centrality[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(row['Value']))
        f.write("\n")
    f.write("\n")
    with open(os.path.join(save_to, '_'.join((filename, 'sorted degree centrality.txt'))), 'w') as f1:
        for index, row in sorted_degree_centrality.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()
    print('degree centrality computed')

    closeness_centrality_dict = dict(nx.closeness_centrality(graph))
    sorted_closeness_centrality = others.rank(closeness_centrality_dict)
    f.write("\n")
    f.write("Top 10 nodes by closeness centrality: ")
    f.write("\n")
    for index, row in sorted_closeness_centrality[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(row['Value']))
        f.write("\n")
    f.write("\n")
    with open(os.path.join(save_to, '_'.join((filename, 'sorted closeness centrality.txt'))), 'w') as f1:
        for index, row in sorted_closeness_centrality.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()
    print('closeness centrality computed')

    # compute betweenness and eigenvector centralities
    betweenness_dict = nx.betweenness_centrality(graph, weight='weight')
    eigenvector_dict = nx.eigenvector_centrality(graph, max_iter=600, weight='weight')

    # add betweenness and eigenvector to node attributes
    sorted_betweenness = others.rank(betweenness_dict)
    with open(os.path.join(save_to, '_'.join((filename, 'sorted betweenness centrality.txt'))), 'w') as f1:
        for index, row in sorted_betweenness.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()

    f.write("Top 10 nodes by betweenness centrality: ")
    f.write("\n")
    for index, row in sorted_betweenness[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(round(row['Value'], 4)))
        f.write("\n")
    f.write("\n")
    print('betweenness centrality computed')
    sorted_eigenvector = others.rank(eigenvector_dict)
    with open(os.path.join(save_to, '_'.join((filename, 'sorted eigenvector centrality.txt'))), 'w') as f1:
        for index, row in sorted_eigenvector.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()
    f.write("Top 10 nodes by eigenvector centrality: ")
    f.write("\n")
    for index, row in sorted_eigenvector[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(round(row['Value'], 4)))
        f.write("\n")
    print('eigenvector centrality computed')
    # degree, in-degree and out-degree
    degree_dict = dict(graph.degree(graph.nodes()))
    sorted_degree = others.rank(degree_dict)
    f.write("\n")
    f.write("Top 10 nodes by degree: ")
    f.write("\n")
    for index, row in sorted_degree[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(row['Value']))
        f.write("\n")
    f.write("\n")

    with open(os.path.join(save_to, '_'.join((filename, 'sorted degree.txt'))), 'w') as f1:
        for index, row in sorted_degree.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()
    print('degree computed')
    in_degree_dict = dict(graph.in_degree(graph.nodes()))
    sorted_in_degree = others.rank(in_degree_dict)
    f.write("Top 10 nodes by in_degree: ")
    f.write("\n")
    for index, row in sorted_in_degree[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(row['Value']))
        f.write("\n")
    f.write("\n")

    with open(os.path.join(save_to, '_'.join((filename, 'sorted in_degree.txt'))), 'w') as f1:
        for index, row in sorted_in_degree.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()
    print('in-degree computed')
    out_degree_dict = dict(graph.out_degree(graph.nodes()))
    sorted_out_degree = others.rank(out_degree_dict)
    f.write("Top 10 nodes by out_degree: ")
    f.write("\n")
    for index, row in sorted_out_degree[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(row['Value']))
        f.write("\n")
    f.write("\n")

    with open(os.path.join(save_to, '_'.join((filename, 'sorted out_degree.txt'))), 'w') as f1:
        for index, row in sorted_out_degree.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()
    print('out-degree computed')
    # computer pagerank
    pr = nx.pagerank(graph, weight='weight')
    sorted_pr = others.rank(pr)
    with open(os.path.join(save_to, '_'.join((filename, 'sorted pagerank.txt'))),
              'w') as f1:
        for index, row in sorted_pr.iterrows():
            f1.write(str(str(row['Gene'])))
            f1.write('\t')
            f1.write(str(row['Value']))
            f1.write('\t')
            f1.write(str(row['Rank']))
            f1.write("\n")
    f1.close()

    f.write("\nTop 10 nodes by pagerank: ")
    f.write("\n")
    for index, row in sorted_pr[:10].iterrows():
        f.write("\t")
        f.write(str(row['Gene']))
        f.write(" ")
        f.write(str(round(row['Value'], 4)))
        f.write("\n")
    print("pagerank computed")

    # close the file
    f.close()


def plot_histogram(data, dataname, filename, save_to, mark_genes=True):
    node = list(data.iloc[:, 0])
    value = list(data.iloc[:, 1])
    dict = {}
    for i in range(len(node)):
        dict[node[i]] = value[i]
    # q25, q75 = np.percentile(value, [25, 75])
    # bin_width = 2 * (q75 - q25) * len(value) ** (-1 / 3)
    # bins = round((max(value) - min(value)) / bin_width)
    plt.hist(value, bins=50, edgecolor='black')
    plt.yscale('log')
    mn, mx = plt.xlim()
    plt.xlim(mn, mx)
    plt.ylabel("Count")
    plt.xlabel(dataname)
    plt.title(f'{dataname} Distribution')
    plt.savefig(os.path.join(save_to, f'{filename} {dataname} Distribution.png'))
    plt.close('all')


def read_stats(filename, stat_path):
    # betweenness centrality
    betweenness_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted betweenness centrality.txt'),
        sep='\t', header=None)

    # degree centrality
    degree_centrality_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted degree centrality.txt'),
        sep='\t', header=None)

    # closeness centrality
    closeness_centrality_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted closeness centrality.txt'),
        sep='\t', header=None)

    # eigenvector
    eigen_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted eigenvector centrality.txt'),
        sep='\t', header=None)
    # degree
    degree_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted degree.txt'),
        sep='\t', header=None)

    # in-degree
    in_degree_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted in_degree.txt'),
        sep='\t', header=None)

    # out-degree
    out_degree_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted out_degree.txt'),
        sep='\t', header=None)

    # pagerank
    # this is the pagerank computed by the built-in function, not the ppr
    pagerank_df = pd.read_csv(os.path.join(
        stat_path, f'{filename}_sorted pagerank.txt'),
        sep='\t', header=None)
    return betweenness_df, degree_centrality_df, closeness_centrality_df, eigen_df, degree_df, in_degree_df, out_degree_df, pagerank_df


def plot_network_stats(filename, stat_path, save_to):
    betweenness_df, degree_centrality_df, closeness_centrality_df, eigen_df, degree_df, in_degree_df, out_degree_df, pagerank_df = read_stats(
        filename, stat_path)
    plot_histogram(betweenness_df, 'Betweenness Centrality', filename, save_to)
    plot_histogram(degree_centrality_df, 'Degree Centrality', filename, save_to)
    plot_histogram(closeness_centrality_df, 'Closeness Centrality', filename, save_to)
    plot_histogram(degree_df, 'Degree', filename, save_to)
    plot_histogram(in_degree_df, 'in-Degree', filename, save_to)
    plot_histogram(out_degree_df, 'out-Degree', filename, save_to)
    plot_histogram(pagerank_df, 'pagerank', filename, save_to)
    plot_histogram(eigen_df, 'Eigenvector Centrality', filename, save_to)
