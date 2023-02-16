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
from math import floor

'''
This file is for analyzing the performance of prune algorithm.
'''

def analyze_prune(cancer_name,known_targets_path,combo_candidates_path,k):
    known_targets = pd.read_csv(known_targets_path,sep = ',',header = 0, index_col = 0)
    combo_candidates = pd.read_csv(combo_candidates_path,sep=',',header=0,index_col=[0,1])
    with open(f'{prune_path}//{cancer_name}//known targets analysis.txt', 'w') as f:
        combo_len = len(combo_candidates)
        known_len = len(known_targets)
        # similarity
        na_count0 = combo_candidates['similarity'].isna().sum()
        na_count1 = known_targets['similarity'].isna().sum()
        na_percentage0 = na_count0 / combo_len
        na_percentage1 = na_count1 / known_len
        # drop na values
        dropna0 = combo_candidates['similarity'].sort_values(ascending = False).dropna()
        dropna1 = known_targets['similarity'].sort_values(ascending = False).dropna()
        f.write('There are {} target combinations'.format(combo_len))
        f.write('\n')
        f.write('There are {} known target combinations'.format(known_len))
        f.write('\n')
        f.write('Similarity:\n')
        f.write('Similarity scores are sorted in descending order\n')
        f.write('There are {} NA values in similarity column in combination candidates\n'.format(na_percentage0))
        f.write('There are {} NA values in similarity column in known target combination\n'.format(na_percentage1))
        f.write('After dropping NA values, \n')
        similarity_percentage = {}
        for pc in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]:
            count = (dropna1 > dropna0[floor(len(dropna0) * pc)]).sum()
            percent = count / len(dropna1)
            similarity_percentage[pc] = percent
            f.write(f'{count} ({percent:.2%}) known target combinations are in the top {pc:.2%} combination candidates\n')
        f.write('\n\n')
        # pdist
        na_count0 = combo_candidates['pdist'].isna().sum()
        na_count1 = known_targets['pdist'].isna().sum()
        na_percentage0 = na_count0 / combo_len
        na_percentage1 = na_count1 / known_len
        dropna0 = combo_candidates['pdist'].sort_values(ascending=False).dropna()
        dropna1 = known_targets['pdist'].sort_values(ascending=False).dropna()
        f.write('Pdist scores:\n')
        f.write('Pdist scores are sorted in descending order\n')
        f.write('There are {} NA values in pdist column in combination candidates\n'.format(na_percentage0))
        f.write('There are {} NA values in pdist column in known target combination\n'.format(na_percentage1))
        f.write('After dropping NA values, \n')
        pdist_percentage = {}
        for pc in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]:
            count = (dropna1 > dropna0[floor(len(dropna0) * pc)]).sum()
            percent = count / len(dropna1)
            pdist_percentage[pc] = percent
            f.write(f'{count} ({percent:.2%}) known target combinations are in the top {pc:.2%} combination candidates\n')
        f.write('\n\n')
        # pdist to oncogenes
        na_count0 = combo_candidates['pdist_oncogenes'].isna().sum()
        na_count1 = known_targets['pdist_oncogenes'].isna().sum()
        na_percentage0 = na_count0 / combo_len
        na_percentage1 = na_count1 / known_len
        dropna0 = combo_candidates['pdist_oncogenes'].sort_values(ascending=True).dropna()
        dropna1 = known_targets['pdist_oncogenes'].sort_values(ascending=True).dropna()
        f.write('Pdist to oncogenes:\n')
        f.write('Pdist to oncogenes are sorted in ascending order\n')
        f.write('There are {} NA values in pdist_to_oncogenes column in combination candidates\n'.format(na_percentage0))
        f.write('There are {} NA values in pdist_to_oncogenes column in known target combination\n'.format(na_percentage1))
        f.write('After dropping NA values, \n')
        pdistonco_percentage = {}
        for pc in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]:
            count = (dropna1 < dropna0[floor(len(dropna0) * pc)]).sum()
            percent = count / len(dropna1)
            pdistonco_percentage[pc] = percent
            f.write(f'{count} ({percent:.2%}) known target combinations are in the top {pc:.2%} combination candidates\n')
        f.write('\n\n')
        # pdist to other genes
        na_count0 = combo_candidates['pdist_othergenes'].isna().sum()
        na_count1 = known_targets['pdist_othergenes'].isna().sum()
        na_percentage0 = na_count0 / combo_len
        na_percentage1 = na_count1 / known_len
        dropna0 = combo_candidates['pdist_othergenes'].sort_values(ascending=True).dropna()
        dropna1 = known_targets['pdist_othergenes'].sort_values(ascending=True).dropna()
        f.write('Pdist to other genes:\n')
        f.write('Pdist to other genes are sorted in ascending order\n')
        f.write(
            'There are {} NA values in pdist_to_other_genes column in combination candidates\n'.format(na_percentage0))
        f.write(
            'There are {} NA values in pdist_to_other_genes column in known target combination\n'.format(na_percentage1))
        f.write('After dropping NA values, \n')
        pdistother_percentage = {}
        for pc in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]:
            count = (dropna1 < dropna0[floor(len(dropna0) * pc)]).sum()
            percent = count / len(dropna1)
            pdistother_percentage[pc] = percent
            f.write(
                f'{count} ({percent:.2%}) known target combinations are in the top {pc:.2%} combination candidates\n')
        f.write('\n\n')
        # distance within combination
        na_count0 = combo_candidates['distance'].isna().sum()
        na_count1 = known_targets['distance'].isna().sum()
        na_percentage0 = na_count0 / combo_len
        na_percentage1 = na_count1 / known_len
        dropna0 = combo_candidates['distance'].sort_values(ascending=True).dropna()
        dropna1 = known_targets['distance'].sort_values(ascending=True).dropna()
        f.write('Distance:\n')
        f.write('Distance are sorted in ascending order\n')
        f.write(
            'There are {} NA values in Distance column in combination candidates\n'.format(na_percentage0))
        f.write(
            'There are {} NA values in Distance column in known target combination\n'.format(
                na_percentage1))
        f.write('After dropping NA values, \n')
        distance_percentage = {}
        for pc in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]:
            count = (dropna1 < dropna0[floor(len(dropna0) * pc)]).sum()
            percent = count / len(dropna1)
            distance_percentage[pc] = percent
            f.write(f'{count} ({percent:.2%}) known target combinations are in the top {pc:.2%} combination candidates\n')
        f.write('\n\n')
        # overlapping
        na_count0 = combo_candidates['overlapping'].isna().sum()
        na_count1 = known_targets['overlapping'].isna().sum()
        na_percentage0 = na_count0 / combo_len
        na_percentage1 = na_count1 / known_len
        dropna0 = combo_candidates['overlapping'].sort_values(ascending=True).dropna()
        dropna1 = known_targets['overlapping'].sort_values(ascending=True).dropna()
        f.write('Overlapping:\n')
        f.write(
            'There are {} NA values in Overlapping column in combination candidates\n'.format(na_percentage0))
        f.write(
            'There are {} NA values in Overlapping column in known target combination\n'.format(
                na_percentage1))
        f.write('After dropping NA values, \n')
        overlapping_percentage = pd.DataFrame(index = [i for i in range(k+1)], columns = ['known targets','combo candidates'])
        for i in range(k+1):
            count0 = (dropna0 == i).sum()
            count1 = (dropna1 == i).sum()
            overlapping_percentage.loc[i,'combo candidates'] = count0/len(dropna0)
            overlapping_percentage.loc[i,'known targets'] = count1/len(dropna1)
            f.write(f'{count0} ({count0/len(dropna0):.2%}) combination candidates have {i} targets overlapped with oncogenes\n')
            f.write(f'{count1} ({count1/len(dropna1):.2%}) known target combinations have {i} targets overlapped with oncogenes\n')
        f.close()

cancer_name = 'Prostate Cancer'
k = 2
analyze_prune(cancer_name,f'{prune_path}//{cancer_name}//known_targets.csv',f'{prune_path}//{cancer_name}//{k}set_combo.csv',k)