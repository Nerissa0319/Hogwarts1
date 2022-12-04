import powerlaw
from scipy import stats
import pandas as pd
import numpy as np
# read data
neutral_path = r'D:\Hogwarts\Output\Neutral\statistics\whole_neutral'
signaling_path = r'D:\Hogwarts\Output\Signaling Network\statistics\whole_signaling'
signal_betweenness = pd.read_csv(f'{signaling_path}/whole_signaling_sorted betweenness centrality.txt',
                                 sep='\t', header=None)
signal_closeness = pd.read_csv(f'{signaling_path}/whole_signaling_sorted closeness centrality.txt',
                                 sep='\t', header=None)
signal_degree_cen = pd.read_csv(f'{signaling_path}/whole_signaling_sorted degree centrality.txt',
                                 sep='\t', header=None)
signal_eigen = pd.read_csv(f'{signaling_path}/whole_signaling_sorted eigenvector centrality.txt',
                                 sep='\t', header=None)
signal_pagerank = pd.read_csv(f'{signaling_path}/whole_signaling_sorted pagerank.txt',
                                 sep='\t', header=None)
signal_d = pd.read_csv(f'{signaling_path}/whole_signaling_sorted degree.txt',
                                 sep='\t', header=None)
neutral_betweenness = pd.read_csv(f'{neutral_path}/whole_neutral_sorted betweenness centrality.txt',
                                 sep='\t', header=None)
neutral_closeness = pd.read_csv(f'{neutral_path}/whole_neutral_sorted closeness centrality.txt',
                                 sep='\t', header=None)
neutral_degree_cen = pd.read_csv(f'{neutral_path}/whole_neutral_sorted degree centrality.txt',
                                 sep='\t', header=None)
neutral_eigen = pd.read_csv(f'{neutral_path}/whole_neutral_sorted eigenvector centrality.txt',
                                 sep='\t', header=None)
neutral_pagerank = pd.read_csv(f'{neutral_path}/whole_neutral_sorted pagerank.txt',
                                 sep='\t', header=None)
neutral_d = pd.read_csv(f'{neutral_path}/whole_neutral_sorted degree.txt',
                                 sep='\t', header=None)

# degree centrality
print(stats.ks_2samp(signal_degree_cen.iloc[:,1],neutral_degree_cen.iloc[:,1]))
print(stats.ks_2samp(signal_betweenness.iloc[:,1],neutral_betweenness.iloc[:,1]))
print(stats.ks_2samp(signal_closeness.iloc[:,1],neutral_closeness.iloc[:,1]))
print(stats.ks_2samp(signal_eigen.iloc[:,1],neutral_eigen.iloc[:,1]))
print(stats.ks_2samp(signal_pagerank.iloc[:,1],neutral_pagerank.iloc[:,1]))

# normal test
print(stats.normaltest(neutral_closeness.iloc[0:4925,1]))
print(stats.normaltest(signal_closeness.iloc[0:4276,1]))
#neutral eigenvector
large = neutral_eigen.iloc[0:75,0]
import networkx as nx
g = nx.read_gexf(r'D:\Hogwarts\Output\Neutral\whole_neutral.gexf')
# for u in large:
#     temp = list(g.out_edges(u,data=True))
#     for neigh in temp:
#
#         print(neigh[2]['edge_types'])
#     print('___________________________________________________________________')
phy = []
for u,v,data in g.edges(data = True):
    if data['edge_types'] == 'Phy':
        phy.append(u)
        phy.append(v)
phy = list(set(phy))
node_neu=[]
gsig = nx.read_gexf(r'D:\Hogwarts\Output\Signaling Network\whole_signaling.gexf')
# for u in g.nodes():
#     if u not in gsig.nodes():
#         node_neu.append(u)
# print(node_neu)
# for u in large:
#     print(g.degree()[u])
# for u in large:
#     temp = g.neighbors(u)
#     neigh = []
#     for v in temp:
#         neigh.append(v)
#     for v in neigh:
#         print(g.degree(v))
#     print('_____________________________________________________________________')
# count = 0
# for u in large:
#     if u in gsig.nodes():
#         count += 1
#         print(gsig.degree(u))

# for u in ['UBB','RPS13','RPS3','RPS11','RPSA','RPS6','RPS20','UBC','RPS10','HBA2','HBA1']:
#     temp = list(signal_betweenness.loc[neutral_d.iloc[:,0]==u,:][1])
#     for n in temp:
#         print(n)

# targets for alpha
import matplotlib.pyplot as plt
plt.figure()
for i in range(9):
    alpha = (i+1)/10
    with open(f'D:\\Hogwarts\\Output\\Signaling Network\\target\\whole_signaling\\pdist\\alpha '
              f'= {alpha}\\target\\cancer_genes_pdist_mean.txt','r') as f:
        temp_dict = f.read()
    f.close()
    import ast
    temp_dict = ast.literal_eval(temp_dict)
    name = list(temp_dict.keys())
    value = list(temp_dict.values())
    plt.plot(value,label=f'alpha={alpha}')

plt.title('Pdist for Different Alpha')
plt.ylabel('PDistance')
plt.legend()
plt.show()
