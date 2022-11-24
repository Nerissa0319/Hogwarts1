import os
import networkx as nx
import pandas as pd
from datetime import datetime
from constant import *


def main():
    # for some time-consuming functions like ppr/distance computing, ALWAYS KEEP FALSE
    # after creating the network for the first time, switch to False
    # remember NOT TO DELETE any .gexf file in output folder
    is_create = True  # create the signaling networks
    is_visualize = True  # visualize the signaling network and save to png file
    is_statistics = True  # compute the statistics of the network
    is_plot_statistics = True  # plot the distribution of the statistics for each network
    is_pathway = True
    # DO NOT DELETE ppr folder under Output, it will take quite a long time to compute the ppr for all networkx
    # KEEP FALSE
    is_compute_PDist = True

    # KEEP FALSE
    is_distance = True  # whether to compute shortest_distance of each network

    is_target_info = True  # whether to extract the information of targeted genes
    is_target_chart = True  # whether to plot the degree, eigen, etc. of targeted genes

    is_target_pdist_plot = True  # plot mean value of dppr between targets/non-targets and cancer/disease genes
    is_target_distance_plot = True  # plot mean value of distance between targets/non-targets and cancer/disease genes
    is_target_test = True  # perform K-S test to check whether targets and non targets are stastitically different

    is_add_random = True  # add random edges
    is_remove_random = True  # remove random edges
    is_cancer_comparison = True  # compare all topological features in different cancer subtypes
    if is_create:
        import HogswartCreate as cr
        cr.create_whole_signaling()  # create signaling networks without neural links and add attributes to nodes
        print('Signaling network creation finished')

    whole_signaling = nx.read_gexf(os.path.join(output_path, 'whole_signaling.gexf'))  # read graph
    # ERBC_signaling = nx.read_gexf(os.path.join(output_path, 'ERBC_signaling.gexf'))
    # TNBC_signaling = nx.read_gexf(os.path.join(output_path, 'TNBC_signaling.gexf'))
    if is_visualize:
        import HogswartViz as viz
        viz.viz_network(whole_signaling, 'whole_signaling',
                        visualization_path)  # visualization network and save to png file
        print('\nvisualization finished\n')

    import HogswartStat as stats
    if is_statistics:
        # summarize the statistics of the network (degree, betweennss centrality...)
        # results are written to txt files under stat_whole_path
        stats.compute_network_stats(whole_signaling, 'whole_signaling', stat_whole_path)
    if is_plot_statistics:
        # plot scatter, histogram and line charts of the statistics of the network
        stats.plot_network_stats('whole_signaling', stat_whole_path, stat_chart_whole_path)

    if is_pathway:
        import HogswartPathway as pathway
        # create and visualize pathway graph
        pathway.pathway_graph(whole_signaling, 'whole_signaling')
        pathway.visualize_pathway('whole_signaling')

    if is_compute_PDist:
        import HogswartPDist as pdist
        # compute pdistance and write the resutls to txt files
        pdist.pdist_alpha(whole_signaling, whole_ppr_path, stat_whole_path, 'whole_signaling', whole_dppr_path,
                          whole_pdist_path, iter=9)
        # the argument 'iter' determines the maximum alpha value
        # for example,
        # if iter = 9, the function will compute the ppr for three networks from alpha = 0.1 to 0.9
        # if iter = 3, the function will comopute the ppr for three networks from alpha = 0.1 to 0.3
        print('ppr computed')

    # compute distance and write the results to txt files
    if is_distance:
        import HogswartDistance as distance
        distance.compute_shortest_distance(whole_signaling, whole_st_path)
        print('shortest distance of whole signaling network computed')

    import HogswartTarget as target
    if is_target_info:
        target.write_target_info('whole_signaling', stat_whole_path, whole_target_path)
        print('information of targets and non-targets written into files')
    if is_target_chart:
        # plot charts of target/non-target to oncogenes
        target.target_chart(whole_target_path, 'whole_signaling')
    # plot charts to show the pdistance between targets/nontargets to oncogenes
    if is_target_pdist_plot:
        alpha = 0
        for i in range(9):
            alpha = (i + 1) / 10
            target.plot_pdist(whole_pdist_path, whole_target_path, alpha)
    # plot charts to show shortest distance between targets/nontargets to oncogenes
    if is_target_distance_plot:
        target.plot_distance(whole_st_path, whole_target_path)
    # perform K-S test to check whether targets and non targets are statistically different
    if is_target_test:
        target.st_diff(whole_st_path, whole_pdist_path, whole_target_path, 'whole_signaling', 9)
    # perform analysis on edge_remvoed_graph and edge_added_graph
    if is_remove_random:
        import HogswartRandom as hogRan
        for i in range(4):
            percentage = 5 * (i + 1) / 100 # remove 5%, 10%, 15%, 20% edges
            hogRan.run_remove(whole_signaling, percentage)
    if is_add_random:
        import HogswartRandom as hogRan
        for i in range(4):
            percentage = 5 * (i + 1) / 100 # add 5%, 10%, 15%, 20% edges
            hogRan.run_add(whole_signaling, percentage)

    if is_cancer_comparison:  # compare data in cancer subtype
        import HogswartCanceType as cancer
        breast_onco, prostate_onco, breast_target, prostate_target, breast_nontarget, prostate_nontarget = cancer.find_subgene(
            whole_signaling)  # prepare oncogenes and targets for specific cancer type
        # if the txt files exist, read the files directly, else, run the function 'find_feature'
        if not os.path.exists(f'{cancer_comparison_path}//prostate_cancer_target_distance.txt'):
            cancer.find_feature(breast_onco, prostate_onco, breast_target, prostate_target, breast_nontarget,
                                prostate_nontarget)  # find pdistance and distance for genes in cancer subtype
        # read the distance and pdistance from files
        prostate_target_distance = pd.read_csv(f'{cancer_comparison_path}//prostate_cancer_target_distance.txt',
                                               sep='\t', header=0, index_col=0)
        prostate_non_distance = pd.read_csv(f'{cancer_comparison_path}//prostate_cancer_nontarget_distance.txt',
                                            sep='\t', header=0, index_col=0)
        prostate_target_pdist = pd.read_csv(f'{cancer_comparison_path}//prostate_cancer_target_pdistance.txt',
                                            sep='\t', header=0, index_col=0)
        prostate_non_pdist = pd.read_csv(f'{cancer_comparison_path}//prostate_cancer_nontarget_pdistance.txt',
                                         sep='\t', header=0, index_col=0)
        breast_target_distance = pd.read_csv(f'{cancer_comparison_path}//breast_cancer_target_distance.txt',
                                             sep='\t', header=0, index_col=0)
        breast_non_distance = pd.read_csv(f'{cancer_comparison_path}//breast_cancer_nontarget_distance.txt',
                                          sep='\t', header=0, index_col=0)
        breast_target_pdist = pd.read_csv(f'{cancer_comparison_path}//breast_cancer_target_pdistance.txt',
                                          sep='\t', header=0, index_col=0)
        breast_non_pdist = pd.read_csv(f'{cancer_comparison_path}//breast_cancer_nontarget_pdistance.txt',
                                       sep='\t', header=0, index_col=0)
        # use boxplot of distance and pdistance to compare targets and non-targets in cancer subtypes
        # target/nontarget_df_ is the distance/pdist dataframe for targets/nontargets,
        # cancertype is 'Prostate Cancer'/'Breast Cancer', featurename is 'Distance'/'Pdistance'
        cancer.box_plot(prostate_target_distance, prostate_non_distance, 'Prostate Cancer', 'Distance', prostate_target)
        cancer.box_plot(prostate_target_pdist, prostate_non_pdist, 'Prostate Cancer', 'PDistance', prostate_target)
        cancer.box_plot(breast_target_distance, breast_non_distance, 'Breast Cancer', 'Distance', breast_target)
        cancer.box_plot(breast_target_pdist, breast_non_pdist, 'Breast Cancer', 'PDistance', breast_target)
        # scatterplot of distance/pdistance mean vs other topolotical features(degree,eigenvector centrality...)
        cancer.scatter_plot(prostate_target_pdist, prostate_non_pdist, 'PDistance', 'Prostate Cancer')
        cancer.scatter_plot(prostate_target_distance, prostate_non_distance, 'Distance', 'Prostate Cancer')
        cancer.scatter_plot(breast_target_pdist, breast_non_pdist, 'PDistance', 'Breast Cancer')
        cancer.scatter_plot(breast_target_distance, breast_non_distance, 'Distance', 'Breast Cancer')
        # plot scatter plot of pdistance vs distance for targets and non targets in cancer subtype
        cancer.pdist_vs_distance(prostate_target_pdist, prostate_non_pdist, prostate_target_distance,
                                 prostate_non_distance,
                                 'Prostate Cancer')
        cancer.pdist_vs_distance(breast_target_pdist, breast_non_pdist, breast_target_distance, breast_non_distance,
                                 'Breast Cancer')


if __name__ == '__main__':
    # track the time of running the project
    now = datetime.now()
    start_time = now.strftime('%d/%m/%Y %H:%M:%S')
    print(f'The start time is {start_time} \n')
    main()
    # record the time of running the code
    print(f'start date and time is: {start_time} \n')
    now = datetime.now()
    end_time = now.strftime('%d/%m/%Y %H:%M:%S')
    print(f'end date and time is: {end_time} \n')
