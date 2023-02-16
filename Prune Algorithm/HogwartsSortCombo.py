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


def jaccard_binary(x, y):
    """A function for finding the similarity between two binary vectors"""
    intersection = np.logical_and(x, y)
    union = np.logical_or(x, y)
    similarity = intersection.sum() / float(union.sum())
    return similarity


def sortByPdist(candidates, k, nodeset, targetset, cancer_network, cancer_name, m):  # sort the candidate combination
    import HogwartsHallmarkAnalysis as hallmark
    # remove the nodes not in the cancer network
    temp = nodeset.copy()
    for u in temp:
        if u not in cancer_network.nodes():
            nodeset.remove(u)
    targets = targetset.copy()
    for u in targets:
        if u not in cancer_network.nodes():
            targetset.remove(u)
    # import pdist dataset
    cancer_pdist_path = f'{prune_path}/{cancer_name}/pdist/alpha = 0.2'
    cancer_dist_path = f'{prune_path}/{cancer_name}/distance'
    pdist_df = pd.read_csv(f'{cancer_pdist_path}/pdist.txt', sep='\t', index_col=0, header=0)
    dist_df = pd.read_csv(f'{cancer_dist_path}/distance.txt', sep='\t', index_col=0, header=0)
    oncogene_pdist = pdist_df.loc[:, nodeset] # pdistance to oncogenes
    othernodes = []
    for u in cancer_network.nodes():
        if u not in nodeset:
            othernodes.append(u)
    other_pdist = pdist_df.loc[:, othernodes] # pdistance to other nodes except oncogenes
    # analyze the hallmarks of known drug targets
    # for each drug with k targets for the input cancer type
    # compute jaccard similarity index for each pair of targets of this drug
    # then compute the mean similarity of the targets in this drug
    # hallmark_k = hallmark.drug_hallmarks[hallmark.drug_hallmarks['No. of Targets'] == k]  # drugs with k targets
    hallmark_k = hallmark.drug_hallmarks.copy()
    drug_hallmark_k = hallmark.drug_hallmarks.copy() # hallmarks for each drug's targets

    # select only rows with the cancer type we want
    for i in hallmark_k.index:
        if type(hallmark_k.loc[i, 'Indications']) == str:
            if cancer_name.lower() not in hallmark_k.loc[i, 'Indications'].lower():
                drug_hallmark_k.drop(index=i, inplace=True)  # drugs with k targets for the input cancer type
        else:
            drug_hallmark_k.drop(index=i, inplace=True)

    # summary of drug hallmarks similarity
    drug_summary = pd.DataFrame(index=drug_hallmark_k.index)  # summary of drug's hallmark similarity and pdist
    for drug in drug_hallmark_k.index:
        temp_similarity = []
        temp_target1 = hallmark.drug_targets.loc[
            drug, 'Targets']  # store the similarity for each pair of targets in drug

        temp_target = temp_target1.copy() # all the targets for the current drug
        for t in temp_target1:
            if t not in cancer_network.nodes():
                temp_target.remove(t) # remove drug targets not in cancer network
        po1 = [] # pdistance to oncogenes
        pr1 = [] # pdistance to other genes
        if len(temp_target) > 0: # for drugs with at least 1 target
            # compute the similarity of each pair of targets in the current drug
            for tc in itertools.combinations(temp_target, 2):  # for each pair of targets in a drug

                h1 = list(hallmark.target_hallmarks.loc[tc[0]][0:-1])  # hallmarks for one target
                h2 = list(hallmark.target_hallmarks.loc[tc[1]][0:-1])  # hallmarks for another target
                h1_h2 = jaccard_binary(h1, h2)  # compute similarity
                temp_similarity.append(h1_h2)  # add the similarity to the temp dataframe
        if not len(temp_similarity) == 0:
            # the hallmark score is the mean of hallmark similarity of all target-pairs in this drug
            # for example:
            # drug A has target t1, t2, t3
            # compute s1 = similarity(t1, t2), s2 = similarity(t1,t3), s3 = similarity(t2, t3)
            # hallmark score for drug A is mean(s1,s2,s3)
            drug_summary.loc[drug, 'Mean Similarity'] = statistics.mean(temp_similarity)  # compute the mean

        # distance between targets
        temp_distance = []
        if len(temp_target) > 0:
            for t in temp_target:
                po1.append(mean(oncogene_pdist.loc[t, :]))
                pr1.append(mean(other_pdist.loc[t, :]))
            if (len(po1) > 0) and (len(pr1) > 0):
                mean_po = mean(po1)
                mean_pr = mean(pr1)
                diff = mean_pr - mean_po  # difference between pdistance to oncogenes and other nodes
                drug_summary.loc[drug, 'pdist_to_oncogenes'] = mean_po
                drug_summary.loc[drug, 'pdist_to_othergenes'] = mean_pr
                drug_summary.loc[drug, 'pdist_score'] = (diff ** 2) / mean_po

            for tar1 in temp_target:
                for tar2 in temp_target:
                    temp_distance.append(dist_df.loc[tar1, tar2])
            drug_summary.loc[drug, 'mean distance'] = mean(temp_distance)

        # overlapping
        overlapping_count = 0
        for t in temp_target:
            if t in nodeset:
                overlapping_count += 1
        drug_summary.loc[drug, 'overlapping count'] = overlapping_count
    drug_summary.to_csv(f'{prune_path}//{cancer_name}//drug_summary.csv', header=True, index=True, sep=',')
    # total_score = pd.DataFrame(index=drug_hallmark_k.index)
    # for drug in drug_hallmark_k.index:
    #     total_score.loc[drug, 'similarity_score'] = target_similarity.loc[drug, 'Mean Similarity']
    #     total_score.loc[drug, 'pdist_score'] = (target_pdist.loc[drug, 'pdist_score'] - target_pdist.loc[:,
    #                                                                                     'pdist_score'].min()) / (
    #                                                    target_pdist.loc[:, 'pdist_score'].max() - target_pdist.loc[
    #                                                                                               :,
    #                                                                                               'pdist_score'].min())
    #     total_score.loc[drug, 'total_score'] = (total_score.loc[drug, 'similarity_score'] + total_score.loc[
    #         drug, 'pdist_score']) / 2
    #
    # total_mean = total_score.loc[:, 'total_score'].mean()

    pdist_oncogenes = {}
    pdist_othergenes = {}
    pdist_scores = {}
    similarity_scores = {}
    distance_scores = {}
    overlapping_scores = {}
    count = 0
    for subset in itertools.combinations(candidates, k):
        # hallmark score:
        temp_similarity1 = []
        for tc in itertools.combinations(subset, 2):
            h1 = list(hallmark.hallMarks_df.loc[tc[0]][1:])  # hallmarks for one target
            h2 = list(hallmark.hallMarks_df.loc[tc[1]][1:])  # hallmarks for another target
            h1_h2 = jaccard_binary(h1, h2)  # compute similarity
            temp_similarity1.append(h1_h2)  # add the similarity to the temp dataframe
        similarity_scores[subset] = statistics.mean(temp_similarity1)  # compute the mean
        po = []  # pdistance to oncogenes
        pr = []  # pdistance to rest genes
        distance_1 = []
        overlapping_count1 = 0
        for u in subset:
            po.append(mean(oncogene_pdist.loc[u, :]))
            pr.append(mean(other_pdist.loc[u, :]))
            if u in nodeset:
                overlapping_count1 += 1
            for v in subset:
                distance_1.append(dist_df.loc[u, v])
        mean_po = mean(po)
        mean_pr = mean(pr)
        diff = mean_pr - mean_po  # difference between pdistance to oncogenes and other nodes
        mean_distance = mean(distance_1)
        pdist_oncogenes[subset] = mean_po
        pdist_othergenes[subset] = mean_pr
        pdist_scores[subset] = (diff ** 2) / mean_po
        distance_scores[subset] = mean_distance
        overlapping_scores[subset] = overlapping_count1

        count += 1
        print(f'{count},{subset}')

    cs_df = pd.concat([pd.Series(d) for d in
                       [similarity_scores, pdist_scores, pdist_oncogenes, pdist_othergenes, distance_scores,
                        overlapping_scores]], axis=1)
    cs_df.columns = ['similarity', 'pdist', 'pdist_oncogenes', 'pdist_othergenes', 'distance', 'overlapping']
    cs_df.to_csv(f'{prune_path}//{cancer_name}//{k}set_combo.csv', header=True, index=True, sep=',')
    # cs_df = pd.read_csv(f'{prune_path}//{cancer_name}//{k}set_combo.csv', header=0, index_col=[0,1], sep=',')
    known_targets = pd.DataFrame()
    for drug in drug_hallmark_k.index:
        temp_target1 = hallmark.drug_targets.loc[
            drug, 'Targets']  # store the similarity for each pair of targets in drug

        temp_target = temp_target1.copy()
        for t in temp_target1:
            if t not in cancer_network.nodes():
                temp_target.remove(t)
        for subset in itertools.combinations(temp_target, k):
            if subset in cs_df.index:
                known_targets[subset] = cs_df.T[subset]
    known_targets = known_targets.T
    known_targets.to_csv(f'{prune_path}//{cancer_name}//known_targets.csv', header=True, index=True, sep=',')

