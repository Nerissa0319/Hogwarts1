from constant import *
import networkx as nx
import json
import pandas as pd


# compute the shortest distance between any pairs of nodes in a network
def compute_shortest_distance(G, save_to):
    distance_df = pd.DataFrame(index=pd.Index(sorted(G.nodes())), columns=pd.Index(sorted(G.nodes())))
    for u in sorted(G.nodes()):
        for v in sorted(G.nodes()):
            try:
                distance_df.loc[u, v] = nx.shortest_path_length(G, u, v)
            except nx.NetworkXNoPath:
                distance_df.loc[u, v] = 0

    distance_df.to_csv(f'{save_to}\distance.txt', index=True, header=True, sep='\t')
