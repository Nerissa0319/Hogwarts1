from constant import *
import networkx as nx
from pygraphviz import *
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
whole_signaling = nx.read_gexf(os.path.join(output_path, 'whole_signaling.gexf'))
target = ['CDK6','BRCA1','PIK3C3']
cancer= ['ERBB2','ESR1']
G=nx.DiGraph()
for u in target:
    G.add_node(u,color='green')
for u in cancer:
    G.add_node(u,color='red')
G.add_edge('CDK6','ERBB2',minlen =3.0,pdist=5.9667)
G.add_edge('CDK6','ESR1',minlen =3.0,pdist=4.8956)
G.add_edge('BRCA1','ERBB2',minlen =3.0,pdist=6.0348)
# G.add_edge('BRCA1','ESR1',minlen =1.0,pdist=2.6792)
G.add_edge('PIK3C3','ERBB2',minlen =3.0,pdist=8.7261)
G.add_edge('PIK3C3','ESR1',minlen =3.0,pdist=7.9536)
# G.add_edge('PGR','ERBB2',minlen =2.0,pdist=4.4482)
# G.add_edge('PGR','ESR1',minlen =1.0,pdist=2.6518)

# pos = nx.spring_layout(G)
pos = graphviz_layout(G,prog='circo')
labels = nx.get_edge_attributes(G,'pdist')
colors = nx.get_node_attributes(G,'color').values()
nx.draw(G,pos, with_labels=True,node_color=colors)
nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
plt.show()