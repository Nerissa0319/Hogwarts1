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
from HogwartsSortCombo import *
from HogwartsSPdist import *
from HogwartsSHallmark import *
from HogwartsConstruct import *


def prune(tumour_type,cancer_name,k):
    onco1, target1, nontarget1 = find_subgene(whole_signaling, tumour_type, cancer_name)
    # create new path if path does not exist
    if not os.path.exists(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf'):
        onco_network = constructNetwork(onco1, target1, whole_signaling, cancer_name)
    else:
        # construct the cancer network
        onco_network = nx.read_gexf(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf')

    target_in_network = target1.copy()
    t_in_network = 0
    for u in target1:
        if u not in onco_network.nodes():
            target_in_network.remove(u)
        else:
            t_in_network += 1

    candidate1 = pruneByHallMarks_s(onco_network, targetset=target1)  # prune by hallmarks

    t1_in_network = 0
    for u in candidate1:
        if u in target1:
            t1_in_network += 1

    # prune by pdistance
    candidate3 = pruneByPdist_s(onco1, target1, candidate1, onco_network, cancer_name)
    sortByPdist(candidate3, k, onco1, target1, onco_network, cancer_name, 10)
    t3_in_network = 0
    for u in candidate3:
        if u in target1:
            t3_in_network += 1
    # print basic information about the cancer
    with open(f'{prune_path}//{cancer_name}//basic_info.csv', 'w') as f:
        print('Cancer Subtype:{}'.format(cancer_name))
        print('k-set:{}'.format(k))
        print('There are {} known targets for this cancer subtype'.format(len(target1)))
        print('There are {} known targets in the cancer network'.format(t_in_network))
        print('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
        print('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(onco_network)))
        print('Number of candidates after pruning by single target hallmarks:{}'.format(len(candidate1)))
        print(f'Number of known targets left after pruning by hallmarks: {t1_in_network}, '
              f'the percentage is {t1_in_network / len(target1):.2%}')
        print('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
        print(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
              f'the percentage is {t3_in_network / t1_in_network:.2%}')
        f.write('Cancer Subtype:{}'.format(cancer_name))
        f.write('\n')
        f.write('k-set:{}'.format(k))
        f.write('\n')
        f.write('There are {} known targets for this cancer subtype'.format(len(target1)))
        f.write('\n')
        f.write('There are {} known targets in the cancer network'.format(t_in_network))
        f.write('\n')
        f.write('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
        f.write('\n')
        f.write('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(onco_network)))
        f.write('\n')
        f.write('Number of candidates after pruning by single target hallmarks:{}'.format(len(candidate1)))
        f.write('\n')
        f.write(f'Number of known targets left after pruning by hallmarks: {t1_in_network}, '
                f'the percentage is {t1_in_network / len(target1):.2%}')
        f.write('\n')
        f.write('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
        f.write('\n')
        f.write(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
                f'the percentage is {t3_in_network / t1_in_network:.2%}')
        f.write('\n')
    f.close()


if __name__ == "__main__":
    prune('prostate','Prostate Cancer',2)
    prune('breast','Breast Cancer',2)
