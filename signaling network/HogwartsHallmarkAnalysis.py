from constant import *
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

'''This code analyzed the hallmarks of drug targets in order to help prune candidate targets'''


def hallmark_pruning(cancer_genes, candidate_targets, k, network):
    # summarize the hallmarks of the input cancer genes
    hallMarks_df = pd.read_csv(os.path.join(hallmark_path, 'AllHallmarks.csv'), index_col=0, header=0)
    hallMarks_name = list(hallMarks_df.columns[1:])
    cancer_index = []
    for u in cancer_genes:
        if u in hallMarks_df.index:
            cancer_index.append(u)
    cancer_hallmark = hallMarks_df.loc[cancer_index, :]
    # summarize the frequency of each hallmarks in the input cancer genes
    cancer_hallmarks_frequency = pd.Series(index=hallMarks_name)
    for h in hallMarks_name:
        if len(cancer_hallmark[h].value_counts()) == 1:
            cancer_hallmarks_frequency[h] = 0 / len(cancer_index)
        else:
            cancer_hallmarks_frequency[h] = cancer_hallmark[h].value_counts()[1] / len(cancer_index)
    cancer_hallmarks_frequency.sort_values(axis=0, inplace=True, ascending=False)
    # extract the first k hallmarks that occur most frequently in cancer genes
    k_hallmarks = list(cancer_hallmarks_frequency.index[0:k])
    # clean the input targets data
    target_index = []
    for u in candidate_targets:
        if u in hallMarks_df.index:
            target_index.append(u)
    target_hallmark = hallMarks_df.loc[target_index, :]  # summarize the hallmarks of input candidate targets
    # prune the nodes which have none of the first k hallmarks
    # then assign the rest of genes to set n if it has the nth frequent hallmark, where n <= k
    k_sets = []  # k is a list of list, where each element is a set of genes with nth hallmark
    for i in range(k):
        k_sets.append([])

    for u in target_hallmark.index:
        for i in range(k):
            if target_hallmark.loc[u, k_hallmarks[i]] == 1:
                k_sets[i].append(u)  # the first element in k_sets is a list of genes with the most frequent hallmarks
    print('k_sets completed')
    # form combination
    # a k-combination is a set of k genes selected from candidate targets
    # genes in one combination are from k different sets based on their hallmarks
    from itertools import product  # to generate combinations
    combination = list(product(*k_sets))  # combine genes from k_sets
    print(len(combination))
    print('first combination completed')
    # some genes have more than one hallmarks, so they might belong to different hallmark groups
    # thus the combination formed might have duplicated elements such as (TP53, TP53, AR)
    # prune the combination with duplicated genes to make sure each candidate combo
    # has k different frequent hallmarks
    combo_k = []
    for c in combination:
        temp = set(c)
        if len(temp) == k:
            combo_k.append(temp)
    print('combo_k completed')
    # then prune the duplicated combinations such as (TP53, ERBB2, AR) and (AR, ERBB2, TP53)
    combo_unique = []
    for c in combo_k:
        if c not in combo_unique:
            combo_unique.append(c)
    print('unique combination completed')
    return k_sets,combo_unique  # combo_unique is the final candidate combo set based on hallmarks


# prune the candidates based on their pdistance to cancer genes
def pdist_pruning(cancer_genes, network,candidate_combo, k, pdist_df, pdist_threshold=12.5129):

    # for a candidate combination, the genes should have path to each of the connected components or
    # the average pdistance should be larger than a threshold
    # the default pdist threshold is 12.5129 which means no connections
    import statistics  # compute mean pdist
    final_candidate = candidate_combo.copy()
    cancer_genes_new = cancer_genes.copy()
    for c in cancer_genes:
        if network.in_degree[c] == 0:
            cancer_genes_new.remove(c)
    for combo in candidate_combo:
        for v in cancer_genes_new:
            combo_pdist = []
            for u in combo:
                combo_pdist.append(pdist_df.loc[u,v])
            if min(combo_pdist) > pdist_threshold:
                final_candidate.remove(combo)
                break

    return final_candidate


from datetime import datetime

# track the time of running the project
now = datetime.now()
start_time = now.strftime('%d/%m/%Y %H:%M:%S')
print(f'The start time is {start_time} \n')

from HogwartsCanceType import *

breast_onco, prostate_onco, breast_target, prostate_target, breast_nontarget, prostate_nontarget = find_subgene(
    whole_signaling)
import random
random_candidate = random.sample(list(whole_signaling.nodes()),70)
sample_targets = random.sample(breast_target,70)
random_candidate = list(set(random_candidate).union(set(sample_targets)))
k_sets, combination1 = hallmark_pruning(breast_onco, random_candidate, 3, whole_signaling)
print(len(combination1))

count = 0
for r in random_candidate:
    if r in breast_target:
        count+=1
print(f'there are {count} targets in {len(random_candidate)} random nodes')
all_target = set()
for ls in k_sets:
    all_target = all_target.union(set(ls))
count = 0
for u in all_target:
    if u in breast_target:
        count += 1
print(f'there are {count} targets in {len(all_target)} candidates')
# record the time of running the code
now = datetime.now()
end_time = now.strftime('%d/%m/%Y %H:%M:%S')
print(f'time completing first phase: {end_time} \n')

combo2 = pdist_pruning(breast_onco,whole_signaling,combination1,3,whole_pdist_df,11)
print(len(combo2))

now = datetime.now()
end_time = now.strftime('%d/%m/%Y %H:%M:%S')
print(f'end time is: {end_time} \n')

# ERBC_disease_nodes = ["ESR1", "PGR", "ERBB2", "TP53", "MUC1", "CEACAM5", "BRCA1", "BRCA2"]  # list of disease nodes
# TNBC_disease_nodes = ["EGFR", "KIT", "KRT5", "TP53", "TOP2A", "PARP1",
#                       "HSP90AA1", "HSP90AB1", "MTOR"]  # list of disease nodes
# whole_disease_nodes = list(set(ERBC_disease_nodes + TNBC_disease_nodes))  # list of disease nodes
#
