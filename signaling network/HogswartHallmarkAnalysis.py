from constant import *
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

'''This code analyzed the hallmarks of drug targets in order to help prune candidate targets'''


def hallmark_pruning(cancer_genes, candidate_targets, k):
    hallMarks_df = pd.read_csv(os.path.join(hallmark_path, 'AllHallmarks.csv'), index_col=0, header=0)
    hallMarks_name = list(hallMarks_df.columns[1:])
    cancer_index = []
    for u in cancer_genes:
        if u in hallMarks_df.index:
            cancer_index.append(u)
    cancer_hallmark = hallMarks_df.loc[cancer_index, :]
    cancer_hallmarks_frequency = pd.Series(index=hallMarks_name)
    for h in hallMarks_name:
        if len(cancer_hallmark[h].value_counts()) == 1:
            cancer_hallmarks_frequency[h] = 0 / len(cancer_index)
        else:
            cancer_hallmarks_frequency[h] = cancer_hallmark[h].value_counts()[1] / len(cancer_index)
    cancer_hallmarks_frequency.sort_values(axis=0, inplace=True, ascending=False)
    k_hallmarks = list(cancer_hallmarks_frequency.index[0:k])
    target_index = []
    for u in candidate_targets:
        if u in hallMarks_df.index:
            target_index.append(u)
    target_hallmark = hallMarks_df.loc[target_index, :]
    k_sets = []
    for i in range(k):
        k_sets.append([])

    for u in target_hallmark.index:
        for i in range(k):
            if target_hallmark.loc[u, k_hallmarks[i]] == 1:
                k_sets[i].append(u)
    # form combination
    for i in range(k):
        print(len(k_sets[i]))
    # print(cancer_hallmarks_frequency)
    return k_sets


from HogswartCancerType import *

breast_onco, prostate_onco, breast_target, prostate_target, breast_nontarget, prostate_nontarget = find_subgene(
    whole_signaling)


# hallmark_pruning(prostate_onco,prostate_target,3)

def paths_pruning(cancer_network, cancer_genes, gene_sets):
    candidates = []
    for u in gene_sets:
        flag = True
        if u not in cancer_network.nodes():
            continue
        else:
            for v in cancer_genes:
                if v not in cancer_network.nodes():
                    continue
                if cancer_network.in_degree(v) == 0:
                    continue
                if not nx.has_path(cancer_network, u, v):
                    flag = False
                    break
            if flag:
                candidates.append(u)

    return candidates


# a = paths_pruning(whole_signaling, prostate_onco, whole_signaling.nodes())
# b = hallmark_pruning(prostate_onco,a,3)
# print(len(prostate_target))
# for i in range(len(b)):
#     print(len(b[i]))
#     count = 0
#     for u in b[i]:
#         if u in prostate_target:
#             count += 1
#     print(count)
#
# seta = []
# for u in b[0]:
#     if u in prostate_target:
#         seta.append(u)
# for u in b[1]:
#     if u in prostate_target:
#         seta.append(u)
# for u in b[2]:
#     if u in prostate_target:
#         seta.append(u)
# seta = list(set(seta))
# print(len(seta))
# for u in seta:
#     print(u)
# print('----------------------------------')
# for u in prostate_target:
#     print(u)


ERBC_disease_nodes = ["ESR1", "PGR", "ERBB2", "TP53", "MUC1", "CEACAM5", "BRCA1", "BRCA2"]  # list of disease nodes
TNBC_disease_nodes = ["EGFR", "KIT", "KRT5", "TP53", "TOP2A", "PARP1",
                      "HSP90AA1", "HSP90AB1", "MTOR"]  # list of disease nodes
whole_disease_nodes = list(set(ERBC_disease_nodes + TNBC_disease_nodes))  # list of disease nodes
a = paths_pruning(whole_signaling, whole_disease_nodes, whole_signaling.nodes())
b = hallmark_pruning(whole_disease_nodes,a,3)
print(len(breast_target))
for i in range(len(b)):
    print(len(b[i]))
    count = 0
    for u in b[i]:
        if u in breast_target:
            count += 1
    print(count)

seta = []
for u in b[0]:
    if u in breast_target:
        seta.append(u)
for u in b[1]:
    if u in breast_target:
        seta.append(u)
for u in b[2]:
    if u in breast_target:
        seta.append(u)
seta = list(set(seta))
print(len(seta))
for u in seta:
    print(u)
print('----------------------------------')
for u in breast_target:
    print(u)