import os
import networkx as nx
import pandas as pd
from datetime import datetime
from ConstantNeu import *


def main():
    # for some time-consuming functions like ppr/distance computing, ALWAYS KEEP FALSE

    # after creating the network for the first time, switch to False
    # remember NOT TO DELETE any .gexf file in output folder
    is_create = False  # create the signaling networks
    is_visualize = False  # visualize the signaling network and save to png file
    is_statistics = False  # compute the statistics of the network
    is_plot_statistics = False  # plot the distribution of the statistics for each network

    # DO NOT DELETE ppr folder under Output, it will take quite a long time to compute the ppr for all networkx
    # KEEP FALSE
    is_compute_PDist = False

    # KEEP FALSE
    is_distance = False  # whether to compute shortest_distance of each network

    is_target_info = False  # whether to extract the information of targeted genes
    is_target_chart = False  # whether to plot the degree, eigen, etc. of targeted genes

    # KEEP FALSE
    is_target_pdist = False  # whether to extract dppr between target/non-target and cancer/disease genes

    is_target_pdist_plot = False  # plot mean value of dppr between targets/non-targets and cancer/disease genes
    # KEEP FALSE
    is_target_distance = False

    is_target_distance_plot = False  # plot mean value of distance between targets/non-targets and cancer/disease genes
    is_target_test = True
    is_add_random = True
    is_remove_random = True
    if is_create:
        import HogswartCreateNeu as cr
        cr.create_whole_neutral()
        print('Signaling network creation finished')

    whole_neutral = nx.read_gexf(os.path.join(output_path, 'whole_neutral.gexf'))
    # ERBC_neutral = nx.read_gexf(os.path.join(output_path, 'ERBC_neutral.gexf'))
    # TNBC_neutral = nx.read_gexf(os.path.join(output_path, 'TNBC_neutral.gexf'))
    if is_visualize:
        import HogswartVizNeu as viz
        viz.viz_network(whole_neutral, 'whole_neutral', visualization_path)
        print('\nvisualization finished\n')

    import HogswartStatNeu as stats
    if is_statistics:
        stats.compute_network_stats(whole_neutral, 'whole_neutral', stat_whole_path)
    if is_plot_statistics:
        stats.plot_network_stats('whole_neutral', stat_whole_path, stat_chart_whole_path)

    if is_compute_PDist:
        import HogswartPDistNeu as pdist
        pdist.pdist_alpha(whole_neutral, whole_ppr_path, stat_whole_path, 'whole_neutral',whole_dppr_path,whole_pdist_path,
                          iter=9)  # the argument 'iter' determines the maximum alpha value
        # for example,
        # if iter = 9, the function will compute the ppr for three networks from alpha = 0.1 to 0.9
        # if iter = 3, the function will comopute the ppr for three networks from alpha = 0.1 to 0.3
        print('ppr computed')

    if is_distance:
        import HogswartDistanceNeu as distance
        for source in sorted(whole_neutral.nodes()):
            distance.compute_shortest_distance(whole_neutral, source, whole_st_path)
        print('shortest distance of whole signaling network computed')

    import HogswartTargetNeu as target
    if is_target_info:
        target.write_target_info('whole_neutral', stat_whole_path, whole_target_path)
        print('information of targets and non-targets written into files')
    if is_target_chart:
        target.target_chart(whole_target_path, 'whole_neutral')
    if is_target_pdist:
        target.target_pdist(whole_pdist_path, whole_target_path, whole_neutral,9)
    if is_target_distance:
        target.target_distance(whole_st_path, whole_target_path, whole_neutral)
    if is_target_pdist_plot:
        target.plot_pdist(whole_target_path,9)
    if is_target_distance_plot:
        target.plot_distance(whole_target_path)
    if is_target_test:
        target.st_diff(whole_target_path, 'whole_neutral',9)
    if is_add_random:
        import HogswartRandomNeu as hogRan
        for i in range(6):
            percentage = 0.05 * (i+1)
            hogRan.run_add(whole_neutral,percentage)
    if is_remove_random:
        import HogswartRandomNeu as hogRan
        for i in range(6):
            percentage = 0.05 * (i+1)
            hogRan.run_remove(whole_neutral,percentage)


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
