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


def sortCombo(candidates, k, cancer_genes, nodetype, targetset, cancer_network, cancer_name,
              m):  # sort the candidate combination
    import HogwartsHallmarkAnalysis as hallmark
    # remove the nodes not in the cancer network
    temp = cancer_genes.copy()
    for u in temp:
        if u not in cancer_network.nodes():
            cancer_genes.remove(u)
    targets = targetset.copy()
    for u in targets:
        if u not in cancer_network.nodes():
            targetset.remove(u)
    # import pdist dataset
    cancer_pdist_path = f'{prune_path}/{cancer_name}_{nodetype}/pdist/alpha = 0.2'
    cancer_dist_path = f'{prune_path}/{cancer_name}_{nodetype}/distance'
    pdist_df = pd.read_csv(f'{cancer_pdist_path}/pdist.txt', sep='\t', index_col=0, header=0)
    dist_df = pd.read_csv(f'{cancer_dist_path}/distance.txt', sep='\t', index_col=0, header=0)
    cancergene_pdist = pdist_df.loc[:, cancer_genes]  # pdistance to cancergenes
    othernodes = []
    for u in cancer_network.nodes():
        if u not in cancer_genes:
            othernodes.append(u)
    other_pdist = pdist_df.loc[:, othernodes]  # pdistance to other nodes except cancergenes
    # analyze the hallmarks of known drug targets
    # for each drug with k targets for the input cancer type
    # compute jaccard similarity index for each pair of targets of this drug
    # then compute the mean similarity of the targets in this drug
    # hallmark_k = hallmark.drug_hallmarks[hallmark.drug_hallmarks['No. of Targets'] == k]  # drugs with k targets
    hallmark_k = hallmark.drug_hallmarks.copy()
    drug_hallmark_k = hallmark.drug_hallmarks.copy()  # hallmarks for each drug's targets

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

        temp_target = temp_target1.copy()  # all the targets for the current drug
        for t in temp_target1:
            if t not in cancer_network.nodes():
                temp_target.remove(t)  # remove drug targets not in cancer network
        po1 = []  # pdistance to cancergenes
        pr1 = []  # pdistance to other genes
        if len(temp_target) > 0:  # for drugs with at least 1 target
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
                po1.append(mean(cancergene_pdist.loc[t, :]))
                pr1.append(mean(other_pdist.loc[t, :]))
            if (len(po1) > 0) and (len(pr1) > 0):
                mean_po = mean(po1)
                mean_pr = mean(pr1)
                diff = mean_pr - mean_po  # difference between pdistance to cancergenes and other nodes
                drug_summary.loc[drug, 'pdist_to_cancergenes'] = mean_po
                drug_summary.loc[drug, 'pdist_to_othergenes'] = mean_pr
                drug_summary.loc[drug, 'pdist_score'] = (diff ** 2) / mean_po
                drug_summary.loc[drug, 'pdist_diff'] = diff
            for tar1 in temp_target:
                for tar2 in temp_target:
                    temp_distance.append(dist_df.loc[tar1, tar2])
            drug_summary.loc[drug, 'mean distance'] = mean(temp_distance)

        # overlapping
        overlapping_count = 0
        for t in temp_target:
            if t in cancer_genes:
                overlapping_count += 1
        drug_summary.loc[drug, 'overlapping count'] = overlapping_count
    drug_summary.to_csv(f'{prune_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_drug_summary.csv',
                        header=True, index=True, sep=',')
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
    pdist_cancergenes = {}
    pdist_othergenes = {}
    pdist_scores = {}
    pdist_diff = {}
    similarity_scores = {}
    distance_scores = {}
    overlapping_scores = {}
    for drug in drug_hallmark_k.index:
        temp_target1 = hallmark.drug_targets.loc[
            drug, 'Targets']  # store the similarity for each pair of targets in drug
        temp_target = temp_target1.copy()
        for t in temp_target1:
            if t not in cancer_network.nodes():
                temp_target.remove(t)

        for subset in itertools.combinations(temp_target, k):
            # hallmark score:
            temp_similarity1 = []
            for tc in itertools.combinations(subset, 2):
                h1 = list(hallmark.hallMarks_df.loc[tc[0]][1:])  # hallmarks for one target
                h2 = list(hallmark.hallMarks_df.loc[tc[1]][1:])  # hallmarks for another target
                h1_h2 = jaccard_binary(h1, h2)  # compute similarity
                temp_similarity1.append(h1_h2)  # add the similarity to the temp dataframe
            similarity_scores[subset] = statistics.mean(temp_similarity1)  # compute the mean
            po = []  # pdistance to cancergenes
            pr = []  # pdistance to rest genes
            distance_1 = []
            overlapping_count1 = 0
            for u in subset:
                po.append(mean(cancergene_pdist.loc[u, :]))
                pr.append(mean(other_pdist.loc[u, :]))
                if u in cancer_genes:
                    overlapping_count1 += 1
                for v in subset:
                    distance_1.append(dist_df.loc[u, v])
            mean_po = mean(po)
            mean_pr = mean(pr)
            diff = mean_pr - mean_po  # difference between pdistance to cancergenes and other nodes
            mean_distance = mean(distance_1)
            pdist_cancergenes[subset] = mean_po
            pdist_othergenes[subset] = mean_pr
            pdist_diff[subset] = diff
            pdist_scores[subset] = (diff ** 2) / mean_po
            distance_scores[subset] = mean_distance
            overlapping_scores[subset] = overlapping_count1
    known_targets = pd.concat([pd.Series(d) for d in
                               [similarity_scores, pdist_scores, pdist_cancergenes, pdist_othergenes, pdist_diff,
                                distance_scores,
                                overlapping_scores]], axis=1)
    known_targets.columns = ['similarity', 'pdist', 'pdist_cancergenes', 'pdist_othergenes', 'pdist_diff', 'distance',
                             'overlapping']
    known_targets.to_csv(f'{prune_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_known_targets.csv',
                         header=True, index=True, sep=',')
    if not os.path.exists(f'{prune_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}set_combo.txt'):
        pdist_cancergenes = {}
        pdist_othergenes = {}
        pdist_scores = {}
        pdist_diff = {}
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
            po = []  # pdistance to cancergenes
            pr = []  # pdistance to rest genes
            distance_1 = []
            overlapping_count1 = 0
            for u in subset:
                po.append(mean(cancergene_pdist.loc[u, :]))
                pr.append(mean(other_pdist.loc[u, :]))
                if u in cancer_genes:
                    overlapping_count1 += 1
                for v in subset:
                    distance_1.append(dist_df.loc[u, v])
            mean_po = mean(po)
            mean_pr = mean(pr)
            diff = mean_pr - mean_po  # difference between pdistance to cancergenes and other nodes
            mean_distance = mean(distance_1)
            pdist_cancergenes[subset] = mean_po
            pdist_othergenes[subset] = mean_pr
            pdist_diff[subset] = diff
            pdist_scores[subset] = (diff ** 2) / mean_po
            distance_scores[subset] = mean_distance
            overlapping_scores[subset] = overlapping_count1
            # if the pdistance_difference or distance within targets are beyond the range of known targets, then prune it
            # out


            count += 1
            print(f'{count},{subset}')

        cs_df = pd.concat([pd.Series(d) for d in
                           [similarity_scores, pdist_scores, pdist_cancergenes, pdist_othergenes, pdist_diff,
                            distance_scores,
                            overlapping_scores]], axis=1)

        cs_df.columns = ['similarity', 'pdist', 'pdist_cancergenes', 'pdist_othergenes', 'pdist_diff', 'distance',
                         'overlapping']

        cs_df.to_csv(f'{prune_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}set_combo.txt', header=True,
                     index=True, sep='\t')
    else:
        cs_df = pd.read_csv(f'{prune_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}set_combo.txt', sep='\t',header=0,index_col=[0,1])
    cs_df1 = cs_df[(cs_df['pdist_diff'] <= known_targets['pdist_diff'].max()) & (
                cs_df['pdist_diff'] >= known_targets['pdist_diff'].min()) & (
                               cs_df['distance'] <= known_targets['distance'].max()) & (
                               cs_df['distance'] >= known_targets['distance'].min())]
    cs_df1.to_csv(f'{prune_path}//{cancer_name}_{nodetype}//{cancer_name}_{nodetype}_{k}_pruned.txt',
                     header=True,
                     index=True, sep='\t')
