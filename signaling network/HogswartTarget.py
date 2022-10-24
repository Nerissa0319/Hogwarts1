import os.path
import statistics
from HogswartStat import *
import pandas as pd
import ast
import csv
import numpy as np


# write all statistical values for targets, non-targets and cancer genes to files
def write_target_info(filename, stat_path, save_to):
    # read the data from files
    betweenness_df, degree_centrality_df, closeness_centrality_df, eigen_df, degree_df, in_degree_df, out_degree_df, pagerank_df = \
        read_stats(filename, stat_path)
    df_list = [betweenness_df, degree_centrality_df, closeness_centrality_df, eigen_df, degree_df, in_degree_df,
               out_degree_df, pagerank_df]
    # write the statistical data to files for targets
    with open(os.path.join(save_to, f'Targeted Genes Info_{filename}.txt'), 'w') as g:
        # write the column names
        data = ['Gene Name', 'Betweenness Centrality', 'Betweenness Rank', 'Degree Centrality',
                'Degree Centrality Rank', 'Closeness Centrality', 'Closeness Centrality Rank', 'Eigenvector Centrality',
                'Eigenvector Rank', 'Degree', 'Degree Rank', 'In-Degree', 'In-Degree Rank', 'Out-Degree',
                'Out-Degree Rank', 'Pagerank', 'Pagerank Rank']
        for col in data:
            g.write(str(col))
            g.write('\t')
        g.write('\n')
        # write the data values
        for u in target_ls:
            if u in betweenness_df.iloc[:, 0].values:
                g.write(str(u))

                for df in df_list:
                    g.write('\t')
                    g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 1].values[0]))
                    g.write('\t')
                    g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 2].values[0]))

                g.write('\n')
    g.close()
    # write the statistical data to files for non-targets
    with open(os.path.join(save_to, f'Non-Targeted Genes Info_{filename}.txt'), 'w') as g:
        data = ['Gene Name', 'Betweenness Centrality', 'Betweenness Rank', 'Degree Centrality',
                'Degree Centrality Rank',
                'Closeness Centrality',
                'Closeness Centrality Rank', 'Eigenvector Centrality',
                'Eigenvector Rank', 'Degree', 'Degree Rank', 'In-Degree', 'In-Degree Rank', 'Out-Degree',
                'Out-Degree Rank', 'Pagerank', 'Pagerank Rank']
        for col in data:
            g.write(str(col))
            g.write('\t')
        g.write('\n')
        for u in betweenness_df.iloc[:, 0].values:
            if u not in target_ls:
                g.write(str(u))
                for df in df_list:
                    g.write('\t')
                    g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 1].values[0]))
                    g.write('\t')
                    g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 2].values[0]))

                g.write('\n')
    g.close()
    # write the statistical data to files for all genes
    with open(os.path.join(save_to, f'All Genes Info_{filename}.txt'), 'w') as g:
        data = ['Gene Name', 'Betweenness Centrality', 'Betweenness Rank', 'Degree Centrality',
                'Degree Centrality Rank', 'Closeness Centrality', 'Closeness Centrality Rank', 'Eigenvector Centrality',
                'Eigenvector Rank', 'Degree', 'Degree Rank', 'In-Degree', 'In-Degree Rank', 'Out-Degree',
                'Out-Degree Rank', 'Pagerank', 'Pagerank Rank']
        for col in data:
            g.write(str(col))
            g.write('\t')
        g.write('\n')
        for u in betweenness_df.iloc[:, 0].values:

            g.write(str(u))
            for df in df_list:
                g.write('\t')
                g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 1].values[0]))
                g.write('\t')
                g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 2].values[0]))

            g.write('\n')
    g.close()
    # write the statistical data to files for cancer genes
    with open(os.path.join(save_to, f'Cancer Genes Info_{filename}.txt'), 'w') as g:
        data = ['Gene Name', 'Betweenness Centrality', 'Betweenness Rank', 'Degree Centrality',
                'Degree Centrality Rank', 'Closeness Centrality', 'Closeness Centrality Rank', 'Eigenvector Centrality',
                'Eigenvector Rank', 'Degree', 'Degree Rank', 'In-Degree', 'In-Degree Rank', 'Out-Degree',
                'Out-Degree Rank', 'Pagerank', 'Pagerank Rank']
        for col in data:
            g.write(str(col))
            g.write('\t')
        g.write('\n')
        for u in cancer_ls:
            if u in betweenness_df.iloc[:, 0].values:
                g.write(str(u))

                for df in df_list:
                    g.write('\t')
                    g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 1].values[0]))
                    g.write('\t')
                    g.write(str(df.loc[df.iloc[:, 0] == u].iloc[:, 2].values[0]))

                g.write('\n')
    g.close()


# plot histogram of all statistical values for targets and non-targets
def plot_target_histogram(target_path, filename):
    # read data from files
    tar_df = pd.read_csv(os.path.join(target_path, f'Targeted Genes Info_{filename}.txt'), sep='\t', header=0)
    non_target_df = pd.read_csv(os.path.join(target_path, f'Non-Targeted Genes Info_{filename}.txt'), sep='\t',
                                header=0)
    cancer_info_df = pd.read_csv(os.path.join(target_path, f'Cancer Genes Info_{filename}.txt'), sep='\t',
                                 header=0)

    # function to plot histogram
    # data argument takes the dataframe of the gene information with all statistical values
    # topology is the string of statistical values (betweenness centrality, degree...)
    # if_target takes 'targets','non-targets','all_genes'
    # file_name is the name of the network ('whole_signaling' in this case)
    # save_to is the directory where the output will be saved
    def plot_target_hist(data, topology, if_target, file_name, save_to, log=True):
        node = list(data.iloc[:, 0])
        value = list(data.loc[:, topology])
        stat_dict = {}
        for i in range(len(node)):
            stat_dict[node[i]] = value[i]
        if log:
            plt.hist(value, bins=50, edgecolor='black')
            plt.yscale('log')
        mn, mx = plt.xlim()
        plt.xlim(mn, mx)
        # set x,y labels
        plt.ylabel("Count")
        plt.xlabel(topology)
        # set title
        plt.title(f'{if_target} {topology} Histogram')
        fig = plt.gcf()
        fig.set_size_inches(12, 9)
        # save to png
        save_to = os.path.join(save_to, topology, 'histogram')
        if not os.path.exists(save_to):
            os.makedirs(save_to)
        plt.savefig(os.path.join(save_to, f'{file_name}_{if_target} {topology} Histogram.png'))
        plt.close('all')

    dataname = ['Betweenness Centrality',
                'Degree Centrality',
                'Closeness Centrality',
                'Eigenvector Centrality',
                'Degree',
                'In-Degree',
                'Out-Degree',
                'Pagerank']
    # execute the plot histogram function
    for stat in dataname:
        plot_target_hist(tar_df, stat, 'Target', filename, target_path, True)
        plot_target_hist(cancer_info_df, stat, 'Cancer Genes', filename, target_path, True)
        # plot_target_hist(erbc_info_df, stat, 'ERBC Disease Nodes', filename, target_path, False)
        # plot_target_hist(tnbc_info_df, stat, 'TNBC Disease Nodes', filename, target_path, False)
        plot_target_hist(non_target_df, stat, 'Non-Target', filename,
                         target_path, True)


# plot scatter of all statistical values for targets and non-targets
def plot_target_scatter(target_path, filename):
    tar_df = pd.read_csv(os.path.join(target_path, f'Targeted Genes Info_{filename}.txt'), sep='\t', header=0)

    non_target_df = pd.read_csv(os.path.join(target_path, f'Non-Targeted Genes Info_{filename}.txt'), sep='\t',
                                header=0)
    cancer_info_df = pd.read_csv(os.path.join(target_path, f'Cancer Genes Info_{filename}.txt'), sep='\t',
                                 header=0)
    all_df = pd.read_csv(os.path.join(target_path, f'All Genes Info_{filename}.txt'), sep='\t', header=0)

    # function to plot scatters
    def plot_target_scat(data, topology, if_target, file_name, save_to):
        node = list(data.iloc[:, 0])
        value = list(data.loc[:, topology])
        stat_dict = {}
        for i in range(len(node)):
            stat_dict[node[i]] = value[i]

        stat_dict = {k: v for k, v in sorted(stat_dict.items(), key=lambda item: item[1], reverse=True)}
        node = list(stat_dict.keys())
        value = list(stat_dict.values())
        # if the data is for all genes (including targets and non-targets)
        # set different markers for targets and non-targets, and show in the legend
        if if_target == 'All Genes':
            tarx = []
            tary = []
            nonx = []
            nony = []
            x0 = 0
            for u in node:
                x0 += 1
                if u in target_ls:
                    tarx.append(x0)
                    tary.append(stat_dict[u])
                else:
                    nonx.append(x0)
                    nony.append(stat_dict[u])

            plt.scatter(tarx, tary, s=5, alpha=0.7, c='red', marker='v', label='Targets')
            plt.scatter(nonx, nony, s=1, alpha=0.4, c='blue', marker='x', label='Non-Targets')
            plt.yscale('log')
            plt.legend()
            # set the xticks to percentage
            start = 1
            end = len(node)
            diff = end - start
            plt.xticks(
                [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
                ['0', '20%', '40%', '60%', '80%', '100%'])
        # if the data is not for all genes, just plot scatter
        else:
            plt.scatter(node, value, s=5, alpha=0.6)
            plt.yscale('log')
            # set the xticks to percentage
            locs, labels = plt.xticks()
            start = locs[0]
            end = locs[-1]
            diff = end - start
            plt.xticks(
                [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
                ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.ylabel(f'log({topology})')
        plt.xlabel('Nodes')
        plt.title(f'{if_target} {topology} Scatter')
        save_to = os.path.join(save_to, topology, 'scatter')
        if not os.path.exists(save_to):
            os.makedirs(save_to)
        plt.savefig(os.path.join(save_to, f'{file_name}_{if_target} {topology} Scatter.png'))
        plt.close('all')

    stat_ls = ['Betweenness Centrality',
               'Degree Centrality',
               'Closeness Centrality',
               'Eigenvector Centrality',
               'Degree',
               'In-Degree',
               'Out-Degree',
               'Pagerank']
    for stat in stat_ls:
        plot_target_scat(tar_df, stat, 'Target', filename, target_path)
        plot_target_scat(cancer_info_df, stat, 'Cancer Genes', filename, target_path)
        plot_target_scat(non_target_df, stat, 'Non-Target', filename, target_path)
        plot_target_scat(all_df, stat, 'All Genes', filename, target_path)


# plot line charts of all statistical values for targets and non-targets
def plot_target_line(target_path, filename):
    tar_df = pd.read_csv(os.path.join(target_path, f'Targeted Genes Info_{filename}.txt'), sep='\t', header=0)

    non_target_df = pd.read_csv(os.path.join(target_path, f'Non-Targeted Genes Info_{filename}.txt'), sep='\t',
                                header=0)
    cancer_info_df = pd.read_csv(os.path.join(target_path, f'Cancer Genes Info_{filename}.txt'), sep='\t',
                                 header=0)
    all_info_df = pd.read_csv(os.path.join(target_path, f'All Genes Info_{filename}.txt'), sep='\t',
                              header=0)

    def plot_target_stats_line(data, dfname, topology='', file_name='', save_to=''):
        plt.figure(figsize=(12, 9))
        for i in range(len(data)):
            temp_data = data[i]
            temp_name = dfname[i]
            node = list(temp_data.iloc[:, 0])
            value = list(temp_data.loc[:, topology])
            stat_dict = {}
            for j in range(len(node)):
                stat_dict[node[j]] = value[j]

            stat_dict = {k: v for k, v in sorted(stat_dict.items(), key=lambda item: item[1], reverse=True)}
            # node = list(stat_dict.keys())
            value = list(stat_dict.values())
            x1 = np.linspace(0, 1, len(value))
            plt.plot(x1, value, label=temp_name)
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.yscale('log')
        plt.legend()
        plt.ylabel(f'log({topology})')
        plt.title(f'{topology}')

        save_to = os.path.join(save_to, topology, 'line')
        if not os.path.exists(save_to):
            os.makedirs(save_to)
        plt.savefig(os.path.join(save_to, f'{file_name}_{topology} Comparison.png'))
        plt.close('all')

    stat_ls = ['Betweenness Centrality',
               'Degree Centrality',
               'Closeness Centrality',
               'Eigenvector Centrality',
               'Degree',
               'In-Degree',
               'Out-Degree',
               'Pagerank']
    # datals = [tar_df, cancer_info_df, erbc_info_df, tnbc_info_df, non_target_df,all_info_df]
    datals = [tar_df, cancer_info_df, non_target_df, all_info_df]
    dataname_ls = ['Target', 'Cancer Genes', 'Non-Target', 'All-Genes']
    # dataname_ls = ['Target', 'Cancer Genes', 'ERBC Disease Nodes', 'TNBC Disease Nodes', 'Non-Target','All-Genes']
    for stat in stat_ls:
        plot_target_stats_line(datals, dataname_ls, stat, filename, target_path)


# execute the plot function of histogram, scatter, line plot
def target_chart(network_target_path, network_name):
    plot_target_histogram(network_target_path, network_name)
    plot_target_scatter(network_target_path, network_name)
    plot_target_line(network_target_path, network_name)
    print('statistics charts of targets and non-targets finished')


# plot histogram and scatter for mean values of pdist for both targets and non-targets
def plot_pdist(network_pdist_path,network_target_path,alpha):
    pdist = pd.read_csv(f'{network_pdist_path}/alpha = {alpha}/pdist.txt',sep='\t',
                        header=0,index_col=0)  # read pdist.txt file to dataframe
    target_pdist = pdist.loc[pdist.index.isin(target_ls),pdist.columns.isin(cancer_ls)] # extract target-oncogene
    nontarget_pdist = pdist.loc[~pdist.index.isin(target_ls),pdist.columns.isin(cancer_ls)] # extract nontarget - oncogene
    all_pdist = pdist.loc[:,pdist.columns.isin(cancer_ls)] # extract all genes - oncogene
    # target_pdist = target_pdist.replace(to_replace = 12.5129,value = np.NaN)
    # nontarget_pdist = nontarget_pdist.replace(to_replace = 12.5129,value = np.NaN)
    # all_pdist = all_pdist.replace(to_replace = 12.5129,value = np.NaN)
    target_mean = target_pdist.mean(axis=1).sort_values(ascending=False) # compute mean value for targets
    nontarget_mean = nontarget_pdist.mean(axis=1).sort_values(ascending=False) # compute mean value for non-targets
    all_mean = all_pdist.mean(axis=1).sort_values(ascending=False) # compute mean value for all genes

    num = 0
    for series in [target_mean,nontarget_mean]:
        plt.figure(figsize=(12,9))
        plt.scatter(series.index,series.values,s=5,alpha=0.6)
        locs,labels=plt.xticks()
        start = locs[0]
        end = locs[-1]
        diff = end - start
        plt.xticks(
            [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
            ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.ylabel('Mean pdist')
        plt.xlabel('Nodes')
        saveto = f'{network_target_path}/pdist/alpha = {alpha}'
        if not os.path.exists(saveto):
            os.makedirs(saveto)
        if num == 0:
            plt.title(r'Target - Oncogenes Mean Pdist Scatter')
            plt.savefig(f'{network_target_path}/pdist/alpha = {alpha}/Target - Oncogenes Mean Pdist Scatter.png')
        else:
            plt.title(r'NonTarget - Oncogenes Mean Pdist Scatter')
            plt.savefig(f'{network_target_path}/pdist/alpha = {alpha}/NonTarget - Oncogenes Mean Pdist Scatter.png')
        plt.close('all')

        plt.figure(figsize=(12, 9))
        plt.hist(series.values, bins=50, edgecolor='black')
        plt.yscale('log')
        plt.ylabel('Count')
        plt.xlabel('Mean pdist')
        if num == 0:
            plt.title(r'Target - Oncogenes Mean Pdist Histogram')
            plt.savefig(f'{network_target_path}/pdist/alpha = {alpha}/Target - Oncogenes Mean Pdist Histogram.png')
        else:
            plt.title(r'NonTarget - Oncogenes Mean Pdist Histogram')
            plt.savefig(f'{network_target_path}/pdist/alpha = {alpha}/NonTarget - Oncogenes Mean Pdist Histogramo.png')

        plt.close('all')
        num += 1

    plt.figure(figsize=(12,9))
    x1=np.linspace(0,1,len(target_mean))
    x2=np.linspace(0,1,len(nontarget_mean))
    x3=np.linspace(0,1,len(all_mean))
    plt.plot(x1,target_mean.values,label='Targets')
    plt.plot(x2,nontarget_mean.values,label='Non Targets')
    plt.plot(x3,all_mean.values,label='All Genes')
    plt.legend()
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
    plt.xlabel('Node')
    plt.title(f'Target vs Non-Target Mean Pdist Line')
    plt.savefig(f'{network_target_path}/pdist/alpha = {alpha}/NonTarget - Oncogenes Mean Pdist Comparison.png')
    plt.close('all')


def plot_distance(network_distance_path,network_target_path):
    distance = pd.read_csv(f'{network_distance_path}/distance.txt',sep='\t',
                        header=0,index_col=0)  # read distance.txt file to dataframe
    target_distance = distance.loc[distance.index.isin(target_ls),distance.columns.isin(cancer_ls)] # extract target-oncogene
    nontarget_distance = distance.loc[~distance.index.isin(target_ls),distance.columns.isin(cancer_ls)] # extract nontarget - oncogene
    all_distance = distance.loc[:,distance.columns.isin(cancer_ls)] # extract all genes - oncogene
    # target_distance = target_distance.replace(to_replace = 0,value = np.NaN)
    # nontarget_distance = nontarget_distance.replace(to_replace = 0,value = np.NaN)
    # all_distance = all_distance.replace(to_replace = 0,value = np.NaN)
    target_mean = target_distance.mean(axis=1).sort_values(ascending=False) # compute mean value for targets
    nontarget_mean = nontarget_distance.mean(axis=1).sort_values(ascending=False) # compute mean value for non-targets
    all_mean = all_distance.mean(axis=1).sort_values(ascending=False) # compute mean value for all genes
    saveto = f'{network_target_path}/distance'
    if not os.path.exists(saveto):
        os.makedirs(saveto)
    num = 0
    for series in [target_mean,nontarget_mean]:
        plt.figure(figsize=(12,9))
        plt.scatter(series.index,series.values,s=5,alpha=0.6)
        locs,labels=plt.xticks()
        start = locs[0]
        end = locs[-1]
        diff = end - start
        plt.xticks(
            [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
            ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.ylabel('Mean distance')
        plt.xlabel('Nodes')

        if num == 0:
            plt.title(r'Target - Oncogenes Mean Distance Scatter')
            plt.savefig(f'{network_target_path}/distance/Target - Oncogenes Mean Distance Scatter.png')
        else:
            plt.title(r'NonTarget - Oncogenes Mean Distance Scatter')
            plt.savefig(f'{network_target_path}/distance/NonTarget - Oncogenes Mean Distance Scatter.png')
        plt.close('all')

        plt.figure(figsize=(12, 9))
        plt.hist(series.values, bins=50, edgecolor='black')
        plt.yscale('log')
        plt.ylabel('Count')
        plt.xlabel('Mean distance')
        if num == 0:
            plt.title(r'Target - Oncogenes Mean Distance Histogram')
            plt.savefig(f'{network_target_path}/distance/Target - Oncogenes Mean Distance Histogram.png')
        else:
            plt.title(r'NonTarget - Oncogenes Mean Distance Histogram')
            plt.savefig(f'{network_target_path}/distance/NonTarget - Oncogenes Mean Distance Histogramo.png')

        plt.close('all')
        num += 1

    plt.figure(figsize=(12,9))
    x1=np.linspace(0,1,len(target_mean))
    x2=np.linspace(0,1,len(nontarget_mean))
    x3=np.linspace(0,1,len(all_mean))
    plt.plot(x1,target_mean.values,label='Targets')
    plt.plot(x2,nontarget_mean.values,label='Non Targets')
    plt.plot(x3,all_mean.values,label='All Genes')
    plt.legend()
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
    plt.xlabel('Node')
    plt.title(f'Target vs Non-Target Mean Distance Line')
    plt.savefig(f'{network_target_path}/distance/NonTarget - Oncogenes Mean Distance Comparison.png')
    plt.close('all')


# conduct the two-sample ks test for targets and non-targets
def st_diff(network_distance_path,network_pdist_path,target_path, filename, iter):
    def ks_test_target(target_data, non_target_data, statistics_name):
        target_value = target_data.loc[:, statistics_name]
        non_target_value = non_target_data.loc[:, statistics_name]
        from scipy import stats
        s, pvalue = stats.ks_2samp(target_value, non_target_value)
        if pvalue < 0.05:
            print(f'pvalue of {pvalue} is lower than the threshold 0.05, reject the null hypothesis.')
            print(
                f'Thus {statistics_name} for targets and non-targets in signaling network is statistically different\n')
        else:
            print(f'pvalue of {pvalue} is not lower than the threshold 0.05, we do not reject the null hypothesis.')
            print(f'Thus {statistics_name} for targets and non-targets in signaling network has no statistical '
                  f'difference\n')

    tar_df = pd.read_csv(os.path.join(target_path, f'Targeted Genes Info_{filename}.txt'), sep='\t', header=0)
    non_target_df = pd.read_csv(os.path.join(target_path, f'Non-Targeted Genes Info_{filename}.txt'), sep='\t',
                                header=0)
    dataname = ['Betweenness Centrality',
                'Degree Centrality',
                'Closeness Centrality',
                'Eigenvector Centrality',
                'Degree',
                'In-Degree',
                'Out-Degree',
                'Pagerank']
    for stat in dataname:
        ks_test_target(tar_df, non_target_df, stat)

    # test for pdist and distance
    import numpy as np
    distance = pd.read_csv(f'{network_distance_path}/distance.txt', sep='\t',
                           header=0, index_col=0)  # read distance.txt file to dataframe
    target_distance = distance.loc[
        distance.index.isin(target_ls), distance.columns.isin(cancer_ls)]  # extract target-oncogene
    nontarget_distance = distance.loc[
        ~distance.index.isin(target_ls), distance.columns.isin(cancer_ls)]  # extract nontarget - oncogene
    all_distance = distance.loc[:, distance.columns.isin(cancer_ls)]  # extract all genes - oncogene
    target_mean = target_distance.mean(axis=1).sort_values(ascending=False)  # compute mean value for targets
    nontarget_mean = nontarget_distance.mean(axis=1).sort_values(ascending=False)  # compute mean value for non-targets
    from scipy import stats
    s, pvalue = stats.ks_2samp(target_mean, nontarget_mean)
    if pvalue < 0.05:
        print(f'pvalue of {pvalue} is lower than the threshold 0.05, reject the null hypothesis.')
        print(
            f'Thus shortest distance for targets and non-targets in signaling network is statistically different\n')
    else:
        print(f'pvalue of {pvalue} is not lower than the threshold 0.05, we do not reject the null hypothesis.')
        print(f'Thus shortest distance for targets and non-targets in signaling network has no statistical '
              f'difference\n')
    for i in range(iter):
        alpha = (i+1)/10
        pdist = pd.read_csv(f'{network_pdist_path}/alpha = {alpha}/pdist.txt', sep='\t',
                            header=0, index_col=0)  # read pdist.txt file to dataframe
        target_pdist = pdist.loc[pdist.index.isin(target_ls), pdist.columns.isin(cancer_ls)]  # extract target-oncogene
        nontarget_pdist = pdist.loc[
            ~pdist.index.isin(target_ls), pdist.columns.isin(cancer_ls)]  # extract nontarget - oncogene
        all_pdist = pdist.loc[:, pdist.columns.isin(cancer_ls)]  # extract all genes - oncogene
        # target_pdist = target_pdist.replace(to_replace = 12.5129,value = np.NaN)
        # nontarget_pdist = nontarget_pdist.replace(to_replace = 12.5129,value = np.NaN)
        # all_pdist = all_pdist.replace(to_replace = 12.5129,value = np.NaN)
        target_mean = target_pdist.mean(axis=1).sort_values(ascending=False)  # compute mean value for targets
        nontarget_mean = nontarget_pdist.mean(axis=1).sort_values(ascending=False)  # compute mean value for non-targets

        from scipy import stats
        s, pvalue = stats.ks_2samp(target_mean, nontarget_mean)
        if pvalue < 0.05:
            print(f'pvalue of {pvalue} is lower than the threshold 0.05, reject the null hypothesis.')
            print(
                f'Thus pdist when alpha = {alpha} for targets and non-targets in signaling network is '
                f'statistically different\n')
        else:
            print(f'pvalue of {pvalue} is not lower than the threshold 0.05, we do not reject the null hypothesis.')
            print(
                f'Thus pdist when alpha = {alpha} for targets and non-targets in signaling network has no '
                f'statistical '
                f'difference\n')
