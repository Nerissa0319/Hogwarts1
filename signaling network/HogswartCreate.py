import os

import networkx as nx
import pandas as pd
from constant import *


# remove nerual links from a network
def remove_phy_edge(graph):
    removed_edge = []
    for u, v in graph.edges():
        if graph.edges[u, v]['edge_types'] == 'Phy':
            removed_edge.append([u, v])

    for e in removed_edge:
        u = e[0]
        v = e[1]
        graph.remove_edge(u, v)

    new_G = nx.DiGraph()
    new_G.add_edges_from(graph.edges())
    node_ls = new_G.nodes()
    result_graph = graph.subgraph(node_ls)
    # nx.write_gexf(result, os.path.join(output_path, '_'.join((filename, 'signaling.gexf'))))
    return result_graph


# read expression level file and add the data to node attributes
def ER_positive_breast_cancer(G):
    mutation_frequency = pd.read_csv(os.path.join(ERBC_input_path, 'COSMIC_ERBC_GeneMutationFrequency_14July2016.csv'),
                                     encoding='unicode_escape')
    gene_expression_fold_change = pd.read_csv(os.path.join(ERBC_input_path, 'GEO2R_ERPos_GSE20437_NormalVsER.csv'),
                                              encoding='unicode_escape')

    for u in list(mutation_frequency['Gene name']):
        if u in G.nodes:
            G.nodes[u]['ER Positive Breast Cancer_Mutation_Frequency'] = {
                "Mutated Samples": list(
                    mutation_frequency.loc[mutation_frequency['Gene name'] == u, "Mutated samples"])[0],
                "Samples Tested": list(
                    mutation_frequency.loc[mutation_frequency['Gene name'] == u, "Samples tested"])[0],
                "Percentage": list(
                    mutation_frequency.loc[mutation_frequency['Gene name'] == u, "Percentage"])[0]
            }
            if len(list(gene_expression_fold_change.loc[gene_expression_fold_change['Gene'] == u, "Fold change"])) != 0:
                G.nodes[u]['ER Positive Breast Cancer_Fold_Change'] = {
                    "Fold Change":
                        list(gene_expression_fold_change.loc[gene_expression_fold_change['Gene'] == u, "Fold change"])[
                            0]
                }


def Triple_negative_cancer(G):
    mutation_frequency = pd.read_csv(os.path.join(TNBC_input_path, 'COSMIC_TNBC_GeneMutationFrequency_12Nov2015.csv'),
                                     encoding='unicode_escape')
    gene_expression_fold_change = pd.read_csv(os.path.join(TNBC_input_path, 'GEO2R_TNBC_GSE38959_TNBCvsNormal.csv'),
                                              encoding='unicode_escape')

    for u in list(mutation_frequency['Gene name']):
        if u in G.nodes:
            G.nodes[u]['Triple Negative Cancer_Mutation_Frequency'] = {
                "Mutated Samples": list(
                    mutation_frequency.loc[mutation_frequency['Gene name'] == u, "Mutated samples"])[0],
                "Samples Tested": list(
                    mutation_frequency.loc[mutation_frequency['Gene name'] == u, "Samples tested"])[0],
                "Percentage": list(
                    mutation_frequency.loc[mutation_frequency['Gene name'] == u, "Percentage"])[0]
            }
            if len(list(gene_expression_fold_change.loc[gene_expression_fold_change['Gene'] == u, "Fold change"])) != 0:
                G.nodes[u]['Triple Negative Cancer_Fold_Change'] = {
                    "Fold Change":
                        list(gene_expression_fold_change.loc[gene_expression_fold_change['Gene'] == u, "Fold change"])[
                            0]
                }


# add the data from protease output file to node attribute
def add_attri(protease_file_name, G):
    protease_df = pd.read_csv(os.path.join(protease_output_path, '.'.join((protease_file_name, 'txt'))),
                              sep="\t", header=0)
    col_names = protease_df.columns

    # add the data of the protease output file to node attribute
    for u in list(protease_df['Gene.symbol']):
        if u in G.nodes:
            if ((len(list(protease_df.loc[protease_df['Gene.symbol'] == u, "ID"])) != 0)
                    # (len(list(protease_df.loc[protease_df['Gene'] == u, "ID"])) != 0) &
                    # (len(list(protease_df.loc[protease_df['Gene'] == u, "ID"])) != 0) &
                    # (len(list(protease_df.loc[protease_df['Gene'] == u, "ID"])) != 0) &
                    # (len(list(protease_df.loc[protease_df['Gene'] == u, "ID"])) != 0) &
            ):

                G.nodes[u][protease_file_name] = {
                }
                for name in col_names:
                    G.nodes[u][protease_file_name][name] = list(
                        protease_df.loc[protease_df['Gene.symbol'] == u, name])[0]

    print(protease_file_name + " expression levels added")
    print(" ")


# add pathways to each gene's node attribute
# find out how many genes do not belong to any pathway
def add_pathway(G):
    # read the pathway gene list from the file downloaded from reactome.org
    df = pd.read_csv(os.path.join(input_path, 'PathwayGeneList.gmt'), sep="|", names=["A"])
    cols = ["pathways", "ID", "Gene"]
    df[cols] = df.A.str.split("\t", n=2, expand=True)
    df['Gene'] = df['Gene'].str.split('\t')
    df = df.drop(columns='A', axis=1)
    all_pathway = pd.read_csv(os.path.join(input_path, 'PathwayRelationship.txt'), sep="\t", header=0)
    parent_pathway_ls = list(all_pathway['Parent'])
    number_no_pathway = 0
    number_no_pathway_list = []

    # stored the genes which do not belong to any pathway
    g1 = open(os.path.join(pathway_txt_path, 'number_no_pathway.txt'), 'w', encoding='utf-8')
    for u in sorted(G.nodes()):
        pathway_dict = {}
        g = open(os.path.join(pathway_txt_path, '.'.join((str(u), 'txt'))), 'w', encoding="utf-8")
        g.write('Pathways\t')
        g.write('ID\n')
        a = 0
        for i in range(len(df)):
            for v in df['Gene'][i]:
                if str(u) == str(v):

                    if df['ID'][i] in parent_pathway_ls:
                        a = 1
                    else:
                        pathway_dict[i] = str(df['ID'][i])
                        g.write(str(df['pathways'][i]))
                        g.write('\t')
                        g.write(str(df['ID'][i]))
                        g.write('\n')

        g.close()
        count = 0

        file = open(os.path.join(pathway_txt_path, '.'.join((str(u), 'txt'))), 'r', encoding="utf-8")
        for line in file:
            if line != "\n":
                count += 1
        if count == 1:
            number_no_pathway_list.append(str(u))
            g1.write(str(u))
            g1.write('\n')
            number_no_pathway += 1

        # add pathways to node attribute
        G.nodes[u]['pathway'] = pathway_dict
    print('Pathways of the genes have been added')

    g1.write('There are ' + str(number_no_pathway) + ' nodes that do not belong to any pathway')
    g1.write('\n')
    percent_str = '{:.1%}'.format(number_no_pathway / G.number_of_nodes())
    g1.write('The fraction is ' + str(percent_str))
    g1.close()


def create_whole_signaling():
    # create a digraph
    whole_graph = nx.DiGraph()
    # read the edges from csv file
    whole_network = pd.read_csv(
        os.path.join(common_path, 'HumanSignalingNet_v6.csv'),
        header=None, skiprows=[0, 1]
    )
    columns = ['Source node Entrez gene ID', 'Source node Entrez gene name',
               'Target node Entrez gene ID', 'Target node Entrez gene name',
               'Edge type']
    whole_network.columns = columns
    source_name = list(whole_network['Source node Entrez gene name'])
    target_name = list(whole_network['Target node Entrez gene name'])
    edge_types = whole_network['Edge type']

    # add edge attributes: edgetype, edgecolor
    # Positive edge = red, negative edge = blue, physical binding edge = green
    edge_colors = []
    types = ['Pos', 'Neg', 'Phy']
    for i in range(len(edge_types)):

        if edge_types[i] == types[0]:
            edge_colors.append('red')
        if edge_types[i] == types[1]:
            edge_colors.append('blue')
        if edge_types[i] == types[2]:
            edge_colors.append('green')

    for i in range(len(source_name)):
        whole_graph.add_edge(source_name[i], target_name[i], edge_types=edge_types[i], edge_colors=edge_colors[i])

    whole_signaling = remove_phy_edge(whole_graph)  # the signaling network which removes all neutral edges

    # add attributes to the network
    # read hallmarks from AllHallmarks.csv
    hallMarks_df = pd.read_csv(os.path.join(hallmark_path, 'AllHallmarks.csv'))
    hallMarks_name = hallMarks_df.columns
    all_nodes = list(hallMarks_df[hallMarks_name[0]])

    for u in all_nodes:
        if u in whole_signaling.nodes():
            for j in range(len(hallMarks_name)):
                if j > 0:
                    whole_signaling.nodes[u][hallMarks_name[j]] = \
                        list(hallMarks_df.loc[hallMarks_df[hallMarks_name[0]] == u,
                                              hallMarks_name[j]])[0]  # add all hallmarks to node attributes

    # whether a node is an essential gene
    # read the human essential genes from humanEssentialGenes.csv
    essential_df = pd.read_csv(os.path.join(common_path, 'humanEssentialGenes.csv'))
    essential_name = essential_df.columns
    essential_nodes = list(essential_df[essential_name[0]])
    # add the attribute to identify whether the node is human essential gene
    # if a node is human essential gene, the value of this attribute will be 1, otherwise, 0
    for u in whole_signaling.nodes:
        if u in essential_nodes:
            whole_signaling.nodes[u]['essential'] = 1
        else:
            whole_signaling.nodes[u]['essential'] = 0

    # add the mutation frequency and fold change to node attributes
    ER_positive_breast_cancer(whole_signaling)
    Triple_negative_cancer(whole_signaling)

    # add attributes from the PROTEASE OUTPUT file
    for filename in protease_filename:
        add_attri(filename, whole_signaling)

    # add pathways of nodes to the attributes
    add_pathway(whole_signaling)

    # save the graph and attributes of the graph to a gexf file
    nx.write_gexf(whole_signaling, os.path.join(output_path, 'whole_signaling.gexf'))

    return whole_signaling


# def create_disease_signaling(graph):
#     ERBC_list = []
#     ERBC_Fold_Change = pd.read_csv(os.path.join(ERBC_input_path, 'GEO2R_ERPos_GSE20437_NormalVsER.csv'),
#                                    encoding='unicode_escape')
#     ERBC_list.extend(list(ERBC_Fold_Change['Gene']))
#     ERBC_list = list(set(ERBC_list))
#     ERBC_nodes = []
#     for u in ERBC_list:
#         if u in list(graph.nodes()):
#             ERBC_nodes.append(u)
#     ERBC_signaling = graph.subgraph(ERBC_nodes)
#     nx.write_gexf(ERBC_signaling, os.path.join(output_path, 'ERBC_signaling.gexf'))
#
#     TNBC_Fold_Change = pd.read_csv(os.path.join(TNBC_input_path, 'GEO2R_TNBC_GSE38959_TNBCvsNormal.csv'),
#                                    encoding="unicode_escape")
#     TNBC_list = []
#     TNBC_list.extend(list(TNBC_Fold_Change['Gene']))
#     TNBC_list = list(set(TNBC_list))
#     TNBC_nodes = []
#     for u in TNBC_list:
#         if u in list(graph.nodes()):
#             TNBC_nodes.append(u)
#     # create the subgraph and compute the statistics of the subgraph
#     TNBC_signaling = graph.subgraph(TNBC_nodes)
#     # write the subgraph to gexf file
#     nx.write_gexf(TNBC_signaling, os.path.join(output_path, 'TNBC_signaling.gexf'))
#
#     return ERBC_signaling, TNBC_signaling
