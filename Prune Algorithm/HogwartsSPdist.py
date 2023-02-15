import sys
import os
sys.path.insert(1, 'D://Code//Hogwarts//code//signaling network')
import pandas as pd
import HogwartsHallmarkAnalysis as hallmark
from constant import *
import networkx as nx
import numpy as np
import itertools
from statistics import mean
from HogwartsCanceType import *


def pruneByPdist_s(nodeset, targetset, candidate, network, cancer_name):  # prune target candidates by pdistance
    temp = nodeset.copy()
    for u in temp:
        if u not in network.nodes():
            nodeset.remove(u)
    temp = targetset.copy()
    for u in temp:
        if u not in network.nodes():
            targetset.remove(u)
    candidate_by_pdist = candidate.copy()
    pdist = pd.read_csv(f'{prune_path}/{cancer_name}/pdist/alpha = 0.2/pdist.txt', sep='\t',
                        header=0, index_col=0)  # read pdist.txt file to dataframe
    target_pdist = pdist.loc[targetset, nodeset]  # extract target-oncogene
    target_check_bytarget = (target_pdist == 12.5129).all(
        axis=1)  # check if no connection from a target to all oncogenes
    # prune out the targets which have no connection to any of the oncogenes
    for t in targetset:
        if target_check_bytarget[t]:
            target_pdist.drop(index=t, inplace=True)
    max_pdist = pd.Series(index=nodeset, dtype='float64')
    min_pdist = pd.Series(index=nodeset, dtype='float64')
    abs_max = pd.Series(index=nodeset, dtype='float64')
    # if the pdist between the candidate and oncogenes out of the threshold range, prune it out
    for u in nodeset:
        q1 = target_pdist.loc[:, u].quantile(0.25)  # 1st quantile
        q3 = target_pdist.loc[:, u].quantile(0.75)  # 3rd quantile
        iqr = q3 - q1  # inter quantile range
        max_temp = target_pdist.loc[:, u].max()  # maximum
        min_temp = target_pdist.loc[:, u].min()  # minimum
        upperbound = q3 + 1.5 * iqr
        lowerbound = q1 - 1.5 * iqr
        max_pdist[u] = min(max_temp, upperbound)  # upper threshold
        min_pdist[u] = max(min_temp, lowerbound)  # lower threshold
        abs_max[u] = 10
    for v in candidate:
        for p in nodeset:
            if (pdist.loc[v, p] == 12.5129) and (max_pdist[p] != 12.5129):
                candidate_by_pdist.remove(v)
                break

    return candidate_by_pdist