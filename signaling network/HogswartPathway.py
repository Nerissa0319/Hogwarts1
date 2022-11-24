from constant import *
import networkx as nx
import HogswartStat as stat
import matplotlib.pyplot as plt

'''create and visualize pathway graph'''


def pathway_graph(G, filename):  # G is the network, filename is the network name
    pathway_list = []  # create an empty list
    for u in sorted(G.nodes()):
        # read pathways of each node from file
        pathway_df = pd.read_csv(os.path.join(pathway_txt_path, '.'.join((str(u), 'txt'))), sep="\t", encoding='utf-8',
                                 header=0)  # each node has a txt file to store all the pathways it belongs to
        pathway_id_df = list(pathway_df['ID'])
        pathway_list.extend(pathway_id_df)  # pathway_list.extend(pathway_id_df)

    pathway_list = list(set(pathway_list))  # remove redundant elements
    df = pd.read_csv(os.path.join(input_path, 'PathwayGeneList.gmt'), sep="|", names=["A"])  # read pathway - gene list
    cols = ["pathways", "ID", "Gene"]
    df[cols] = df.A.str.split("\t", n=2, expand=True)
    df['Gene'] = df['Gene'].str.split('\t')
    df = df.drop(columns='A', axis=1)

    G_pathway = nx.DiGraph()

    # add node to pathway graph
    G_pathway.add_nodes_from(pathway_list)  # add pathway ids as a node to a graph
    for u in G_pathway.nodes():  # add attributes
        # 'gene' is the attraibute of a pathway which store all the genes belong to the pathway
        G_pathway.nodes[u]['Gene'] = {'Gene List': list(df.loc[df['ID'] == u, 'Gene'])[0]}
        # 'name' is the name of the pathway
        G_pathway.nodes[u]['Name'] = {'Name': list(df.loc[df['ID'] == u, 'pathways'])[0]}

    count = 0
    # add edges in pathway graph
    for p in pathway_list:
        p_gene = list(df['Gene'][df['ID'] == p])[0]
        p_gene_set = set(p_gene)
        # if there is at least one gene that belong to pathway p and q at the same time,
        # then we add an edge between p and q
        for q in pathway_list:
            if p == q:
                continue
            else:

                if G_pathway.has_edge(p, q):
                    continue
                else:
                    q_gene = list(df['Gene'][df['ID'] == q])[0]
                    q_gene_set = set(q_gene)
                    if p_gene_set & q_gene_set:
                        G_pathway.add_edge(p, q)
        count += 1
    print('Edges of the pathway graph have been added')
    # compute statistics of the pathway network
    stat.compute_network_stats(G_pathway, '_'.join((filename, 'pathway', 'stats')), pathway_path)

    # write the graph to gexf file

    nx.write_gexf(G_pathway, os.path.join(pathway_path, f'{filename}_pathway.gexf'))

    return


# visualize the pathway graph
def visualize_pathway(filename):
    G_pathway = nx.read_gexf(f'{pathway_path}//{filename}_pathway.gexf')
    node_colors = []
    edge_colors = []
    for u in range(G_pathway.number_of_nodes()):
        node_colors.append('red')
    for e in range(G_pathway.number_of_edges()):
        edge_colors.append('blue')

    plt.figure(figsize=(200, 200))
    pos = nx.spring_layout(G_pathway, k=0.15, iterations=20)
    nx.draw(G_pathway,
            pos=pos,
            node_color=node_colors,
            edge_color=edge_colors,
            node_size=3000,
            with_labels=True,
            alpha=0.75,
            font_color='grey',
            font_size=150)
    plt.axis('equal')
    plt.savefig(f'{pathway_path}//{filename}_pathway.png')
    plt.close()
    return
