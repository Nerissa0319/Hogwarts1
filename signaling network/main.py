import os
import networkx as nx
import pandas as pd
from datetime import datetime
from constant import *


def main():
    # for some time-consuming functions like ppr/distance computing, ALWAYS KEEP FALSE
    # after creating the network for the first time, switch to False
    # remember NOT TO DELETE any .gexf file in output folder
    is_create = False  # create the signaling networks
    is_visualize = False  # visualize the signaling network and save to png file
    is_statistics = False  # compute the statistics of the network
    is_plot_statistics = False  # plot the distribution of the statistics for each network
    is_pathway = False
    # DO NOT DELETE ppr folder under Output, it will take quite a long time to compute the ppr for all networkx
    # KEEP FALSE
    is_compute_PDist = False

    # KEEP FALSE
    is_distance = False  # whether to compute shortest_distance of each network

    is_target_info = False # whether to extract the information of targeted genes
    is_target_chart = False  # whether to plot the degree, eigen, etc. of targeted genes

    is_target_pdist_plot = False # plot mean value of dppr between targets/non-targets and cancer/disease genes
    is_target_distance_plot = False  # plot mean value of distance between targets/non-targets and cancer/disease genes
    is_target_test = True

    is_add_random = False
    is_remove_random = False
    is_cancer_comparison = False
    if is_create:
        import HogswartCreate as cr
        cr.create_whole_signaling()
        print('Signaling network creation finished')

    whole_signaling = nx.read_gexf(os.path.join(output_path, 'whole_signaling.gexf'))
    # ERBC_signaling = nx.read_gexf(os.path.join(output_path, 'ERBC_signaling.gexf'))
    # TNBC_signaling = nx.read_gexf(os.path.join(output_path, 'TNBC_signaling.gexf'))
    if is_visualize:
        import HogswartViz as viz
        viz.viz_network(whole_signaling, 'whole_signaling', visualization_path)
        print('\nvisualization finished\n')

    import HogswartStat as stats
    if is_statistics:
        stats.compute_network_stats(whole_signaling, 'whole_signaling', stat_whole_path)
    if is_plot_statistics:
        stats.plot_network_stats('whole_signaling', stat_whole_path, stat_chart_whole_path)

    if is_pathway:
        import HogswartPathway as pathway
        pathway.pathway_graph(whole_signaling, 'whole_signaling')
        pathway.visualize_pathway('whole_signaling')

    if is_compute_PDist:
        import HogswartPDist as pdist
        pdist.pdist_alpha(whole_signaling, whole_ppr_path, stat_whole_path, 'whole_signaling', whole_dppr_path,
                          whole_pdist_path, iter=9)
        # the argument 'iter' determines the maximum alpha value
        # for example,
        # if iter = 9, the function will compute the ppr for three networks from alpha = 0.1 to 0.9
        # if iter = 3, the function will comopute the ppr for three networks from alpha = 0.1 to 0.3
        print('ppr computed')

    if is_distance:
        import HogswartDistance as distance
        distance.compute_shortest_distance(whole_signaling, whole_st_path)
        print('shortest distance of whole signaling network computed')

    import HogswartTarget as target
    if is_target_info:
        target.write_target_info('whole_signaling', stat_whole_path, whole_target_path)
        print('information of targets and non-targets written into files')
    if is_target_chart:
        target.target_chart(whole_target_path, 'whole_signaling')
    if is_target_pdist_plot:
        alpha = 0
        for i in range(9):
            alpha = (i+1)/10
            target.plot_pdist(whole_pdist_path,whole_target_path,alpha)
    if is_target_distance_plot:
        target.plot_distance(whole_st_path,whole_target_path)
    if is_target_test:
        target.st_diff(whole_st_path,whole_pdist_path,whole_target_path, 'whole_signaling', 9)

    if is_remove_random:
        import HogswartRandom as hogRan
        for i in range(4):
            percentage = 5 * (i+1)/100
            for j in range(5):
                hogRan.run_remove(whole_signaling, percentage,j+1)
    if is_add_random:
        import HogswartRandom as hogRan
        for i in range(4):
            percentage = 5 * (i+1)/100
            for j in range(5):
                hogRan.run_add(whole_signaling, percentage,j+1)

    if is_cancer_comparison:
        import HogswartCanceType as cancer
        cancer.pdist_vs_distance()
        cancer.pdist_vs_other()
        cancer.distance_vs_other()

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
