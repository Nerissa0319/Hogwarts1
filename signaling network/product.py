import networkx as nx
g = nx.DiGraph()
g.add_edge(1,2)
g.add_edge(1,3)
g.add_edge(2,4)
g.add_edge(2,1)
g.add_edge(1,6)
g.add_edge(5,6)

print(g.in_degree[4])