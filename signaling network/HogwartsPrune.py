import pandas as pd
import HogwartsHallmarkAnalysis as hallmark
from constant import *
import networkx as nx
from typing import List
from math import comb
import itertools
from statistics import mean
from sklearn import preprocessing


# construct a smallest connected network containing all the input cancer genes
def constructNetwork(nodeset, targetset, g, cancer_name):  # construct cancer network
    temp1 = nodeset.copy()  # nodeset is the set of oncogenes of this kind of cancer type
    nodeset1 = []
    for u in temp1:
        if u in g.nodes():
            nodeset1.append(u)  # remove the oncogenes which are not in the network
    temp2 = targetset.copy()  # target set is the set of drug targets for the input cancer type
    targetset1 = []
    for v in temp2:
        if v in g.nodes():
            targetset1.append(v)  # remove the targets which are not in the network
    print('There are {} cancer genes and {} targets'.format(len(nodeset1), len(targetset1)))
    new_nodes = set()
    count = 0
    # if there is at least one path from a target u to oncogene v
    # if the shortest path length <= 5, add the paths with the length <=5 to the network
    # else, add all the paths from u to v to the network.
    for u in targetset1:
        for v in nodeset1:
            if nx.has_path(g, u, v):
                lp = nx.shortest_path_length(g, u, v)
                # print('length of the shortest path between {} and {} is {}'.format(u,v,lp))
                if lp <= 5:
                    temp = nx.all_simple_paths(g, u, v, 5)
                else:
                    temp = nx.all_simple_paths(g, u, v, lp)
                for p in temp:
                    new_nodes.update(p)
            # if there is no path from target u to oncogene v, add u and v to network only.
            else:
                new_nodes.update({u, v})

                # construct cancer network
    cancer_network = nx.subgraph(g, new_nodes)
    cancer_network_1 = cancer_network.copy()
    # cancer_network_1.remove_nodes_from(list(nx.isolates(cancer_network)))
    # scc_set = nx.strongly_connected_components(g)
    # sorted_scc = sorted(scc_set, key=len, reverse=True)
    # gscc = nx.subgraph(g, sorted_scc[0])  # find the largest strongly connected components
    # node_in_gscc = []
    # node_notin_gscc = []
    # for u in nodeset:
    #     if u in gscc.nodes():
    #         node_in_gscc.append(u)  # cancer nodes that are also in the greatest connected component
    #     else:
    #         node_notin_gscc.append(u)  # cancer nodes that are not in the greatest connected component
    # gscc_nodes = list(gscc.nodes())  # all the nodes in the greatest connected component
    # # for any node in the connected component,
    # # if the removal of it will not affect the connection between cancer nodes,
    # # then we will remove it from the subgraph
    # for u in gscc_nodes:
    #     if not u in node_in_gscc:
    #         if not u in targetset:
    #             temp_graph = gscc.copy()
    #             temp_graph.remove_node(u)
    #             if nx.is_strongly_connected(temp_graph):
    #                 gscc = temp_graph.copy()
    # nodes_to_be_added = []  # add path between nodes_in_gscc and nodes_not_in_gscc
    # # for node not in gscc, add the shortest path between these nodes and the nodes in gscc
    # for v in node_notin_gscc:
    #     for u in node_in_gscc:
    #         if nx.has_path(g, u, v):
    #             nodes_to_be_added.extend(nx.shortest_path(g, u, v))
    # nodes_to_be_added = list(set(nodes_to_be_added))
    # nodes_of_network = list(gscc.nodes()).copy()
    # nodes_of_network.extend(nodes_to_be_added)
    # cancer_network = nx.subgraph(g, nodes_of_network)

    # compute statistics and pdist within the cancer network
    for i in [f'{output_path}/prune/{cancer_name}/ppr',
              f'{output_path}/prune/{cancer_name}/dppr',
              f'{output_path}/prune/{cancer_name}/pdist',
              f'{output_path}/prune/{cancer_name}/stats',
              f'{output_path}/prune/{cancer_name}/distance']:
        if not os.path.exists(i):
            os.makedirs(i)
    # write the network to gexf file
    nx.write_gexf(cancer_network_1, f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf')
    import HogwartsStat as stats
    # compute the topological features of the network
    stats.compute_network_stats(cancer_network_1, cancer_name, f'{output_path}/prune/{cancer_name}/stats')
    import HogwartsPDist as pdist
    # compute the pdist of the network
    pdist.pdist_alpha(cancer_network_1,
                      f'{output_path}/prune/{cancer_name}/ppr',
                      f'{output_path}/prune/{cancer_name}/stats',
                      cancer_name, f'{output_path}/prune/{cancer_name}/dppr',
                      f'{output_path}/prune/{cancer_name}/pdist', iter=2)
    import HogwartsDistance as dist
    dist.compute_shortest_distance(cancer_network_1,f'{output_path}/prune/{cancer_name}/distance')

    return cancer_network_1


# check if any pair of nodes in a given nodeset are connected in a graph
# def check_connected(nodeset, g):
#     flag = True
#     for u in nodeset:
#         for v in nodeset:
#             if not u == v:
#                 if not nx.has_path(g, u, v):
#                     flag = False
#     return flag


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


# def pruneByTopo_s(nodeset, targetset, candidate, network, cancer_name):
#     topological_names = ['degree centrality', 'betweenness centrality',
#                          'eigenvector centrality', 'pagerank', 'closeness centrality']
#     temp1 = nodeset.copy()
#     for u in temp1:
#         if u not in network.nodes():
#             nodeset.remove(u)
#     temp2 = targetset.copy()
#     for u in temp2:
#         if u not in network.nodes():
#             targetset.remove(u)
#     candidate_by_tp = candidate.copy()
#     for tp in topological_names:
#         tp_df = pd.read_csv(f'{prune_path}//{cancer_name}//stats//{cancer_name}_sorted {tp}.txt',
#                             sep='\t', header=None, index_col=0)
#         tp_df.columns = ['Value', 'Rank']
#         tp_df_target = tp_df.loc[targetset, :]
#         q1 = tp_df_target.loc[:, 'Value'].quantile(0.25)
#         q3 = tp_df_target.loc[:, 'Value'].quantile(0.75)
#         max_tp = tp_df_target.loc[:, 'Value'].max()
#         min_tp = tp_df_target.loc[:, 'Value'].min()
#         iqr = q3 - q1
#         lower_bound = max(min_tp, q1 - 1.5 * iqr)
#         upper_bound = min(max_tp, q3 + 1.5 * iqr)
#         for c in candidate:
#             if c in candidate_by_tp:
#                 if tp_df.loc[c, 'Value'] < lower_bound or tp_df.loc[c, 'Value'] > upper_bound:
#                     candidate_by_tp.remove(c)
#
#     return candidate_by_tp


# def pruneByHallMarks_c(candidate_targets, k, cancer_name):
#     import itertools
#     candidate_combination = []
#     hallmark_k = hallmark.drug_hallmarks[hallmark.drug_hallmarks['No. of Targets'] == k]
#     drug_hallmark_k = hallmark_k.copy()
#     for i in hallmark_k.index:
#         if type(hallmark_k.loc[i, 'Indications']) == list:
#             if cancer_name not in hallmark_k.loc[i, 'Indications']:
#                 drug_hallmark_k.drop(index=i, inplace=True)
#     # find hallmarks never targeted by the drug targets of this cancer type
#     hallmark_never_targeted = []
#     if_targeted = (drug_hallmark_k == 0).all()
#     for h in hallmark.hallmark_names:
#         if if_targeted[h]:
#             hallmark_never_targeted.append(h)
#
#     if len(hallmark_never_targeted) == 0:
#         print('No hallmarks that are never targeted by the drugs of this cancer')
#     else:
#         print('Hallmarks never targeted:')
#         for h in hallmark_never_targeted:
#             print('----------{}-----------'.format(h))
#
#     # prune out candidates which target the hallamarks in hallmark_never_targeted
#     candidate_targets_1 = candidate_targets.copy()
#     if len(hallmark_never_targeted) > 0:
#         for c in candidate_targets_1:
#             for h in hallmark_never_targeted:
#                 if hallmark.hallMarks_df[c, h] == 1:
#                     candidate_targets.remove(c)
#                     break
#     max_h = pd.Series(index=hallmark.hallmark_names, dtype='float64')
#     min_h = pd.Series(index=hallmark.hallmark_names, dtype='float64')
#     score_combo = []
#     # for h in hallmark.hallmark_names:
#     #     max_h[h] = drug_hallmark_k[h].max()
#     #     min_h[h] = drug_hallmark_k[h].min()
#     for i in drug_hallmark_k.index:
#         temp_str = ''
#         for h in hallmark.hallmark_names:
#             temp_str = ''.join((temp_str, str(drug_hallmark_k.loc[i, h])))
#         score_combo.append(temp_str)
#     with open(f'{prune_path}//{cancer_name}//{k}set_target_hallmarks_combo.txt', 'w') as f:
#         for i in score_combo:
#             f.write(i)
#             f.write('\n')
#     f.close()
#     c = 0
#     with open(f'{prune_path}//{cancer_name}//{k}set_combination.txt', 'w') as f:
#         for subset in itertools.permutations(candidate_targets, k):
#             if set(subset) not in candidate_combination:  # skip if duplicated
#                 score = pd.Series(0, index=hallmark.hallmark_names)
#                 for t in subset:
#                     for h in hallmark.hallmark_names:
#                         if hallmark.hallMarks_df.loc[t][h] == 1:
#                             score[h] += 1
#
#                 # score = score / k
#                 flag = True
#                 temp_str1 = ''
#                 for h in hallmark.hallmark_names:
#                     temp_str1 = ''.join((temp_str1, str(score[h])))
#
#                 if temp_str1 not in score_combo:
#                     flag = False
#                 # for h in hallmark.hallmark_names:
#                 #     if score[h] < min_h[h] or score[h] > max_h[h]:
#                 #         flag = False
#                 #         break
#                 if flag:
#                     candidate_combination.append(set(subset))
#                     f.write(str(subset))
#                     f.write('\n')
#             print(c)
#             c += 1
#             # if c > 50000:
#             #     break
#         f.close()
#     return candidate_targets, candidate_combination

def jaccard_binary(x, y):
    """A function for finding the similarity between two binary vectors"""
    intersection = np.logical_and(x, y)
    union = np.logical_or(x, y)
    similarity = intersection.sum() / float(union.sum())
    return similarity


def sortByPdist(candidates, k, nodeset, targetset,cancer_network, cancer_name, m):  # sort the candidate combination
    import HogwartsHallmarkAnalysis as hallmark
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
    dist_df = pd.read_csv(f'{cancer_dist_path}/distance.txt',sep='\t', index_col = 0, header = 0)
    oncogene_pdist = pdist_df.loc[:, nodeset]
    othernodes = []
    for u in cancer_network.nodes():
        if u not in nodeset:
            othernodes.append(u)
    other_pdist = pdist_df.loc[:, othernodes]
    # analyze the hallmarks of known drug targets
    # for each drug with k targets for the input cancer type
    # compute jaccard similarity index for each pair of targets of this drug
    # then compute the mean similarity of the targets in this drug
    # hallmark_k = hallmark.drug_hallmarks[hallmark.drug_hallmarks['No. of Targets'] == k]  # drugs with k targets
    hallmark_k = hallmark.drug_hallmarks.copy()
    drug_hallmark_k = hallmark.drug_hallmarks.copy()

    for i in hallmark_k.index:
        if type(hallmark_k.loc[i, 'Indications']) == str:
            if cancer_name.lower() not in hallmark_k.loc[i, 'Indications'].lower():
                drug_hallmark_k.drop(index=i, inplace=True)  # drugs with k targets for the input cancer type
        else:
            drug_hallmark_k.drop(index=i, inplace=True)

    drug_summary = pd.DataFrame(index=drug_hallmark_k.index) # summary of drug's hallmark similarity and pdist
    for drug in drug_hallmark_k.index:
        temp_similarity = []
        temp_target1 = hallmark.drug_targets.loc[
            drug, 'Targets']  # store the similarity for each pair of targets in drug

        temp_target = temp_target1.copy()
        for t in temp_target1:
            if t not in cancer_network.nodes():
                temp_target.remove(t)
        po1 = []
        pr1 = []
        if len(temp_target) > 0:
            for tc in itertools.combinations(temp_target, 2):  # for each pair of targets in a drug

                h1 = list(hallmark.target_hallmarks.loc[tc[0]][0:-1])  # hallmarks for one target
                h2 = list(hallmark.target_hallmarks.loc[tc[1]][0:-1])  # hallmarks for another target
                h1_h2 = jaccard_binary(h1, h2)  # compute similarity
                temp_similarity.append(h1_h2)  # add the similarity to the temp dataframe
        if not len(temp_similarity) == 0:
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
                drug_summary.loc[drug,'pdist_to_oncogenes'] = mean_po
                drug_summary.loc[drug,'pdist_to_othergenes'] = mean_pr
                drug_summary.loc[drug, 'pdist_score'] = (diff ** 2) / mean_po

            for tar1 in temp_target:
                for tar2 in temp_target:
                    temp_distance.append(dist_df.loc[tar1,tar2])
            drug_summary.loc[drug,'mean distance'] = mean(temp_distance)

        # overlapping
        overlapping_count = 0
        for t in temp_target:
            if t in nodeset:
                overlapping_count += 1
        drug_summary.loc[drug,'overlapping count'] = overlapping_count
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

    combination = []
    visited = []
    candidate_summary = pd.DataFrame()
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
                distance_1.append(dist_df.loc[u,v])
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

    cs_df = pd.concat([pd.Series(d) for d in [similarity_scores, pdist_scores, pdist_oncogenes,pdist_othergenes,distance_scores,overlapping_scores]], axis=1)
    cs_df.columns = ['similarity', 'pdist', 'pdist_oncogenes', 'pdist_othergenes','distance','overlapping']
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
        for subset in itertools.combinations(temp_target,k):
            if subset in cs_df.index:
                known_targets[subset] = cs_df.T[subset]
    known_targets = known_targets.T
    known_targets.to_csv(f'{prune_path}//{cancer_name}//known_targets.csv', header=True, index=True, sep=',')

if __name__ == "__main__":
    from HogwartsCanceType import *

    breast_onco, prostate_onco, breast_target, prostate_target, breast_nontarget, prostate_nontarget = find_subgene(
        whole_signaling)
    # arguments
    cancer_name = 'Prostate Cancer'
    k = 2

    # create new path if path does not exist
    if not os.path.exists(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf'):
        onco_network = constructNetwork(prostate_onco, prostate_target, whole_signaling, cancer_name)
    else:
        # construct the cancer network
        onco_network = nx.read_gexf(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf')

    target_in_network = prostate_target.copy()
    t_in_network = 0
    for u in prostate_target:
        if u not in onco_network.nodes():
            target_in_network.remove(u)
        else:
            t_in_network += 1

    candidate1 = pruneByHallMarks_s(onco_network, targetset=prostate_target)

    t1_in_network = 0
    for u in candidate1:
        if u in prostate_target:
            t1_in_network += 1

    # candidate2 = pruneByTopo_s(nodeset=prostate_onco, targetset=prostate_target,
    #                            candidate=candidate1, network=onco_network, cancer_name=cancer_name)
    # print('number of targets after pruning by single topological features:{}'.format(len(candidate2)))
    #
    # count = 0
    # for u in candidate2:
    #     if u in prostate_target:
    #         count += 1
    # print(f'  where {count} are real targets for {cancer_name} drugs, the percentage is {count / len(prostate_target):.2%}')
    # prune by pdistance
    candidate3 = pruneByPdist_s(prostate_onco, prostate_target, candidate1, onco_network, cancer_name)
    sortByPdist(candidate3, k, prostate_onco, prostate_target, onco_network, 'Prostate Cancer', 10)
    t3_in_network = 0
    for u in candidate3:
        if u in prostate_target:
            t3_in_network += 1
    # print basic information about the cancer
    with open(f'{prune_path}//{cancer_name}//basic_info.csv', 'w') as f:
        print('Cancer Subtype:{}'.format(cancer_name))
        print('k-set:{}'.format(k))
        print('There are {} known targets for this cancer subtype'.format(len(prostate_target)))
        print('There are {} known targets in the cancer network'.format(t_in_network))
        print('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
        print('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(onco_network)))
        print('Number of candidates after pruning by single target hallmarks:{}'.format(len(candidate1)))
        print(f'Number of known targets left after pruning by hallmarks: {t1_in_network}, '
              f'the percentage is {t1_in_network / len(prostate_target):.2%}')
        print('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
        print(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
              f'the percentage is {t3_in_network / t1_in_network:.2%}')
        f.write('Cancer Subtype:{}'.format(cancer_name))
        f.write('\n')
        f.write('k-set:{}'.format(k))
        f.write('\n')
        f.write('There are {} known targets for this cancer subtype'.format(len(prostate_target)))
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
                f'the percentage is {t1_in_network / len(prostate_target):.2%}')
        f.write('\n')
        f.write('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
        f.write('\n')
        f.write(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
                f'the percentage is {t3_in_network / t1_in_network:.2%}')
        f.write('\n')
    f.close()
    # candidate4, can_combo = pruneByHallMarks_c(candidate3, k, cancer_name=cancer_name)
    # count = 0
    # for u in candidate4:
    #     if u in prostate_target:
    #         count += 1
    # print('number of candidate targets in the end:{},'.format(len(candidate4)))
    # print(f'  where {count} are real targets for {cancer_name} drugs, the percentage is {count / len(prostate_target):.2%}')
    # print('number of possible combinations:{}'.format(comb(len(candidate4), k)))
    # print('number of target combinations:{}'.format(len(can_combo)))
    # print('the percentage is:{0:.2%}'.format(len(can_combo) / comb(len(candidate4), k)))
    # print('-----------------------------------------------------------------')

# # prune on Prostate Cancer
# # import the targets and oncogenes of prostate cancer
# from HogwartsCanceType import *
#
# breast_onco, prostate_onco, breast_target, prostate_target, breast_nontarget, prostate_nontarget = find_subgene(
#     whole_signaling)
# # arguments
# cancer_name = 'Prostate Cancer'
# k = 2
#
# # create new path if path does not exist
# if not os.path.exists(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf'):
#     onco_network = constructNetwork(prostate_onco, prostate_target, whole_signaling, cancer_name)
# else:
#     # construct the cancer network
#     onco_network = nx.read_gexf(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf')
#
# target_in_network = prostate_target.copy()
# t_in_network = 0
# for u in prostate_target:
#     if u not in onco_network.nodes():
#         target_in_network.remove(u)
#     else:
#         t_in_network += 1
#
# candidate1 = pruneByHallMarks_s(onco_network, targetset=prostate_target)
#
# t1_in_network = 0
# for u in candidate1:
#     if u in prostate_target:
#         t1_in_network += 1
#
# # candidate2 = pruneByTopo_s(nodeset=prostate_onco, targetset=prostate_target,
# #                            candidate=candidate1, network=onco_network, cancer_name=cancer_name)
# # print('number of targets after pruning by single topological features:{}'.format(len(candidate2)))
# #
# # count = 0
# # for u in candidate2:
# #     if u in prostate_target:
# #         count += 1
# # print(f'  where {count} are real targets for {cancer_name} drugs, the percentage is {count / len(prostate_target):.2%}')
# # prune by pdistance
# candidate3 = pruneByPdist_s(prostate_onco, prostate_target, candidate1, onco_network, cancer_name)
# sortByPdist(candidate3, k, prostate_onco, prostate_target,onco_network, 'Prostate Cancer', 10)
# t3_in_network = 0
# for u in candidate3:
#     if u in prostate_target:
#         t3_in_network += 1
# # print basic information about the cancer
# with open(f'{prune_path}//{cancer_name}//basic_info.csv','w') as f:
#     print('Cancer Subtype:{}'.format(cancer_name))
#     print('k-set:{}'.format(k))
#     print('There are {} known targets for this cancer subtype'.format(len(prostate_target)))
#     print('There are {} known targets in the cancer network'.format(t_in_network))
#     print('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
#     print('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(onco_network)))
#     print('Number of candidates after pruning by single target hallmarks:{}'.format(len(candidate1)))
#     print(f'Number of known targets left after pruning by hallmarks: {t1_in_network}, '
#           f'the percentage is {t1_in_network / len(prostate_target):.2%}')
#     print('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
#     print(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
#           f'the percentage is {t3_in_network / t1_in_network:.2%}')
#     f.write('Cancer Subtype:{}'.format(cancer_name))
#     f.write('\n')
#     f.write('k-set:{}'.format(k))
#     f.write('\n')
#     f.write('There are {} known targets for this cancer subtype'.format(len(prostate_target)))
#     f.write('\n')
#     f.write('There are {} known targets in the cancer network'.format(t_in_network))
#     f.write('\n')
#     f.write('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
#     f.write('\n')
#     f.write('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(onco_network)))
#     f.write('\n')
#     f.write('Number of candidates after pruning by single target hallmarks:{}'.format(len(candidate1)))
#     f.write('\n')
#     f.write(f'Number of known targets left after pruning by hallmarks: {t1_in_network}, '
#           f'the percentage is {t1_in_network / len(prostate_target):.2%}')
#     f.write('\n')
#     f.write('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
#     f.write('\n')
#     f.write(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
#           f'the percentage is {t3_in_network / t1_in_network:.2%}')
#     f.write('\n')
# f.close()
# # candidate4, can_combo = pruneByHallMarks_c(candidate3, k, cancer_name=cancer_name)
# # count = 0
# # for u in candidate4:
# #     if u in prostate_target:
# #         count += 1
# # print('number of candidate targets in the end:{},'.format(len(candidate4)))
# # print(f'  where {count} are real targets for {cancer_name} drugs, the percentage is {count / len(prostate_target):.2%}')
# # print('number of possible combinations:{}'.format(comb(len(candidate4), k)))
# # print('number of target combinations:{}'.format(len(can_combo)))
# # print('the percentage is:{0:.2%}'.format(len(can_combo) / comb(len(candidate4), k)))
# # print('-----------------------------------------------------------------')
#
#
# cancer_name = 'Breast Cancer'
# k = 2
#
# # create new path if path does not exist
# if not os.path.exists(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf'):
#     onco_network = constructNetwork(breast_onco, breast_target, whole_signaling, cancer_name)
# else:
#     # construct the cancer network
#     onco_network = nx.read_gexf(f'{output_path}/prune/{cancer_name}/{cancer_name}.gexf')
#
# target_in_network = breast_target.copy()
# t_in_network = 0
# for u in breast_target:
#     if u not in onco_network.nodes():
#         target_in_network.remove(u)
#     else:
#         t_in_network += 1
#
# candidate1 = pruneByHallMarks_s(onco_network, targetset=breast_target)
#
# t1_in_network = 0
# for u in candidate1:
#     if u in breast_target:
#         t1_in_network += 1
#
# # candidate2 = pruneByTopo_s(nodeset=breast_onco, targetset=breast_target,
# #                            candidate=candidate1, network=onco_network, cancer_name=cancer_name)
# # print('number of targets after pruning by single topological features:{}'.format(len(candidate2)))
# #
# # count = 0
# # for u in candidate2:
# #     if u in breast_target:
# #         count += 1
# # print(f'  where {count} are real targets for {cancer_name} drugs, the percentage is {count / len(breast_target):.2%}')
# # prune by pdistance
# candidate3 = pruneByPdist_s(breast_onco, breast_target, candidate1, onco_network, cancer_name)
# t3_in_network = 0
# for u in candidate3:
#     if u in breast_target:
#         t3_in_network += 1
# sortByPdist(candidate3, k, breast_onco,breast_target, onco_network, 'Breast Cancer', 10)
# t3_in_network = 0
# for u in candidate3:
#     if u in breast_target:
#         t3_in_network += 1
# # print basic information about the cancer
# with open(f'{prune_path}//{cancer_name}//basic_info.csv','w') as f:
#     print('Cancer Subtype:{}'.format(cancer_name))
#     print('k-set:{}'.format(k))
#     print('There are {} known targets for this cancer subtype'.format(len(breast_target)))
#     print('There are {} known targets in the cancer network'.format(t_in_network))
#     print('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
#     print('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(onco_network)))
#     print('Number of candidates after pruning by single target hallmarks:{}'.format(len(candidate1)))
#     print(f'Number of known targets left after pruning by hallmarks: {t1_in_network}, '
#           f'the percentage is {t1_in_network / len(breast_target):.2%}')
#     print('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
#     print(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
#           f'the percentage is {t3_in_network / t1_in_network:.2%}')
#     f.write('Cancer Subtype:{}'.format(cancer_name))
#     f.write('\n')
#     f.write('k-set:{}'.format(k))
#     f.write('\n')
#     f.write('There are {} known targets for this cancer subtype'.format(len(breast_target)))
#     f.write('\n')
#     f.write('There are {} known targets in the cancer network'.format(t_in_network))
#     f.write('\n')
#     f.write('Number of nodes in signaling network:{}'.format(nx.number_of_nodes(whole_signaling)))
#     f.write('\n')
#     f.write('Number of nodes in the cancer network:{}'.format(nx.number_of_nodes(onco_network)))
#     f.write('\n')
#     f.write('Number of candidates after pruning by single target hallmarks:{}'.format(len(candidate1)))
#     f.write('\n')
#     f.write(f'Number of known targets left after pruning by hallmarks: {t1_in_network}, '
#           f'the percentage is {t1_in_network / len(breast_target):.2%}')
#     f.write('\n')
#     f.write('Number of targets after pruning by single pdistance:{}'.format(len(candidate3)))
#     f.write('\n')
#     f.write(f'Number of known targets left after pruning by pdistance: {t3_in_network}, '
#           f'the percentage is {t3_in_network / t1_in_network:.2%}')
#     f.write('\n')
# f.close()
# # candidate4, can_combo = pruneByHallMarks_c(candidate3, k, cancer_name=cancer_name)
# # count = 0
# # for u in candidate4:
# #     if u in breast_target:
# #         count += 1
# # print('number of candidate targets in the end:{},'.format(len(candidate4)))
# # print(f'  where {count} are real targets for {cancer_name} drugs, the percentage is {count / len(breast_target):.2%}')
# # print('number of possible combinations:{}'.format(comb(len(candidate4), k)))
# # print('number of target combinations:{}'.format(len(can_combo)))
# # print('the percentage is:{0:.2%}'.format(len(can_combo) / comb(len(candidate4), k)))
# # print('-----------------------------------------------------------------')