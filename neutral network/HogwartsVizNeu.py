from ConstantNeu import *
import os
import matplotlib.pyplot as plt
import networkx as nx


# this function visualizes a networkx graph and save the visualization to file
def viz_network(graph, filename, save_to):
    # add edge colors according to the edge types
    colors = []
    for e in graph.edges():
        colors.append(graph.edges[e]['edge_colors'])

    plt.figure(figsize=(100, 100))
    pos = nx.spring_layout(graph, k=0.15, iterations=20)
    nx.draw_networkx_nodes(graph, pos, node_size=1000, alpha=0.75)
    nx.draw_networkx_edges(graph, pos, width=0.7, edge_color=colors, arrows=True, arrowsize=10)
    nx.draw_networkx_labels(graph, pos, font_size=15, font_color='grey')
    plt.axis('equal')
    print(f'{filename} visualized')

    # write to png file
    plt.savefig(os.path.join(save_to, '.'.join((filename, 'png'))))
    print(f'{filename} visualization saved to file')
    plt.close('all')
    return




