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


def pruneByHallMarks_s(cancer_network, targetset):  # prune single candidate by hallmarks
    candidate = []
    temp = targetset.copy()  # target set for the input cancer type
    for u in temp:
        if u not in cancer_network.nodes():
            targetset.remove(u)  # remove oncogene which are not in cancer network.
    cancer_target_hallmark = hallmark.target_hallmarks.loc[targetset, :]
    hallmark_binary = list(
        hallmark.targets_analysis(target_hallmarks_df=cancer_target_hallmark).index)  # hallmarks for known targets
    for u in cancer_network.nodes():
        temp_str = ''
        for h in hallmark.hallMarks_df.columns[1:]:
            temp_str = ''.join([temp_str, str(hallmark.hallMarks_df.loc[u, h])])
        if temp_str in hallmark_binary:  # if nodes in cancer network have the same hallmarks as known targets
            candidate.append(u)  # add to candidates list
    candidate = list(set(candidate))
    return candidate