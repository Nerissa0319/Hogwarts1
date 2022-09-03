from constant import *
import networkx as nx
import json

# compute the shortest distance between any pairs of nodes in a network
def compute_shortest_distance(G, source, save_to):
    if source in G.nodes():
        shortest_distance_dict = {}
        sorted_st_dict = {}
        for t in G.nodes():
            if not t == source:
                try:
                    # for any u, v in disease nodes
                    # if there is a path between u and v
                    # we add all the nodes on the shortest path from u to v, to the graph
                    shortest_distance_dict[str((source, t))] = nx.shortest_path_length(G, source, t)
                except nx.NetworkXNoPath:
                    # if there is no path between u and v
                    # print on screen
                    shortest_distance_dict[str((source, t))] = 0
                    # sort the shortest distance
                sorted_st_dict = {k: v for k, v in
                                  sorted(shortest_distance_dict.items(), key=lambda item: item[1], reverse=True)}

        with open(os.path.join(save_to, f'{source}_shortest distance.txt'), 'w') as f:
            json_str = json.dumps(sorted_st_dict, indent=0)
            f.write(json_str)
            f.write('\n')
        f.close()

