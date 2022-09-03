import os.path
import statistics
from HogswartStatNeu import *
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


def target_pdist(network_pdist_path, network_target_path, network,iter):
    # find pdist between all target/non-target nodes and cancer nodes
    def find_pdist(pdist_path, save_to, alpha, graph, target, gene_list, filename, genetype):
        if target not in graph.nodes():
            return
        pdist_alpha_path = os.path.join(pdist_path, f'alpha = {alpha}')
        save_to = os.path.join(save_to, 'pdist', f'alpha = {alpha}', genetype, filename)
        if not os.path.exists(save_to):
            os.makedirs(save_to)

        pairs = []
        for gene in gene_list:
            if gene in graph.nodes():
                temp = f"('{target}', '{gene}')"
                pairs.append(temp)
        pdist_txt = os.path.join(pdist_alpha_path, f'{target}_pdist.txt')
        with open(pdist_txt, 'r') as file:
            dict = file.read()
        file.close()
        dict = ast.literal_eval(dict)
        selected_ppr_dict = {}
        for key in dict.keys():
            if key in pairs:
                selected_ppr_dict[key] = dict[key]
        with open(os.path.join(save_to, f'{target}_{filename}_pdist.txt'), 'w') as f:
            json_str = json.dumps(selected_ppr_dict, indent=0)
            f.write(json_str)
            f.write('\n')
        f.close()
    # execute computing pdist for different alpha values
    for i in range(iter):
        alpha = (i + 1) / 10
        for target in target_ls:
            find_pdist(network_pdist_path, network_target_path, alpha, network, target, cancer_ls, 'cancer_genes',
                       'target')

        for nontar in network:
            if nontar not in target_ls:
                find_pdist(network_pdist_path, network_target_path, alpha, network, nontar, cancer_ls,
                           'cancer_genes',
                           'non-target')
        for all in network:
            find_pdist(network_pdist_path, network_target_path, alpha, network,
                       all, cancer_ls, 'cancer_genes', 'all_genes')
    print('PDist or targets and non-targets computed')


    # function to compute the mean value of pdist
    import statistics
    def compute_pdist_mean(pdist_path, save_to, filename, genetype, alpha):
        pdist_path = os.path.join(pdist_path, 'pdist', f'alpha = {alpha}',
                                  genetype, filename)
        save_to = os.path.join(save_to, 'pdist', f'alpha = {alpha}',
                               genetype)
        mean_dict = {}
        count = 0
        file_ls = list(os.listdir(pdist_path))
        for file in file_ls:
            idx = str(file).find('_')
            name = str(file)[:idx]
            with open(os.path.join(pdist_path, file), 'r') as f:
                temp_dict = f.read()
            f.close()
            pdist_dict = ast.literal_eval(temp_dict)

            pdist_value = list(pdist_dict.values())
            pdist_mean = statistics.mean(pdist_value)
            mean_dict[name] = pdist_mean
            count += 1
        sorted_mean_dict = {k: v for k, v in sorted(mean_dict.items(), key=lambda item: item[1], reverse=True)}
        with open(os.path.join(save_to, f'{filename}_pdist_mean.txt'), 'w') as f:
            json_str = json.dumps(sorted_mean_dict, indent=0)
            f.write(json_str)
            f.write('\n')
        f.close()

    for i in range(iter):
        alpha = (i + 1) / 10
        # for network_path in [network_target_path, erbc_target_path, tnbc_target_path]:
        for network_path in [network_target_path]:
            # for gene_ls in ['cancer_genes', 'ERBC_disease_nodes', 'TNBC_disease_nodes']:
            for gene_ls in ['cancer_genes']:
                for genetype in ['target', 'non-target', 'all_genes']:
                    compute_pdist_mean(network_path, network_path,
                                       gene_ls, genetype, alpha)


def target_distance(network_st_path, network_target_path, network):
    # find distance between all target/non-target nodes and cancer/disease nodes
    def find_distance(distance_path, save_to, graph, tar, gene_list, filename, genetype):
        if tar not in graph.nodes():
            return
        else:
            distance_path = os.path.join(distance_path)
            save_to = os.path.join(save_to, 'distance', genetype, filename)
            if not os.path.exists(save_to):
                os.makedirs(save_to)
            pairs = []
            for gene in gene_list:
                if gene in graph.nodes():
                    temp = f"('{tar}', '{gene}')"
                    pairs.append(temp)
            distance_txt = os.path.join(distance_path, f'{tar}_shortest distance.txt')
            with open(distance_txt, 'r') as file:
                distance_dict = file.read()
            file.close()
            distance_dict = ast.literal_eval(distance_dict)
            selected_distance_dict = {}
            for key in distance_dict.keys():
                if key in pairs:
                    selected_distance_dict[key] = distance_dict[key]
            with open(os.path.join(save_to, f'{tar}_{filename}_distance.txt'), 'w') as f:
                json_str = json.dumps(selected_distance_dict, indent=0)
                f.write(json_str)
                f.write('\n')
            f.close()

    for target in target_ls:
        find_distance(network_st_path, network_target_path, network, target, cancer_ls, 'cancer_genes', 'target')
    for non_tar in network:
        find_distance(network_st_path, network_target_path, network, non_tar, cancer_ls, 'cancer_genes',
                      'all_genes')
        if non_tar not in target_ls:
            find_distance(network_st_path, network_target_path, network, non_tar, cancer_ls, 'cancer_genes',
                          'non-target')

    def compute_distance_mean(distance_path, save_to, filename, genetype):
        distance_path = os.path.join(distance_path, 'distance', genetype, filename)
        save_to = os.path.join(save_to, 'distance', genetype)
        mean_dict = {}
        count = 0
        file_ls = list(os.listdir(distance_path))
        for file in file_ls:
            idx = str(file).find('_')
            name = str(file)[:idx]
            with open(os.path.join(distance_path, file), 'r') as f:
                temp_dict = f.read()
            f.close()
            distance_dict = ast.literal_eval(temp_dict)

            distance_value = list(distance_dict.values())
            distance_mean = statistics.mean(distance_value)
            mean_dict[name] = distance_mean
            count += 1
        sorted_mean_dict = {k: v for k, v in sorted(mean_dict.items(), key=lambda item: item[1], reverse=True)}
        with open(os.path.join(save_to, f'{filename}_distance_mean.txt'), 'w') as f:
            json_str = json.dumps(sorted_mean_dict, indent=0)
            f.write(json_str)
            f.write('\n')
        f.close()

    # for network_path in [network_target_path, erbc_target_path, tnbc_target_path]:
    for network_path in [network_target_path]:
        # for gene_ls in ['cancer_genes', 'ERBC_disease_nodes', 'TNBC_disease_nodes']:
        for gene_ls in ['cancer_genes']:
            for genetype in ['target', 'non-target', 'all_genes']:
                compute_distance_mean(network_path, network_path,
                                      gene_ls, genetype)


# plot histogram and scatter for mean values of pdist for both targets and non-targets
def plot_pdist(network_target_path,iter):
    def plot_pdist_mean(mean_dict, save_to, filename, genetypes, alpha):
        save_to = os.path.join(save_to, 'pdist', f'alpha = {alpha}',
                               genetype)
        node = list(mean_dict.keys())
        value = list(mean_dict.values())
        plt.figure(figsize=(12, 9))
        if genetypes == 'all_genes':
            tarx = []
            tary = []
            nonx = []
            nony = []
            x0 = 0
            for u in node:
                x0 += 1
                if u in target_ls:
                    tarx.append(x0)
                    tary.append(mean_dict[u])
                else:
                    nonx.append(x0)
                    nony.append(mean_dict[u])

            plt.scatter(tarx, tary, s=5, alpha=0.7, c='red', marker='v', label='Targets')
            plt.scatter(nonx, nony, s=1, alpha=0.4, c='blue', marker='x', label='Non-Targets')
            plt.legend()
            start = 1
            end = len(node)
            diff = end - start
            plt.xticks(
                [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
                ['0', '20%', '40%', '60%', '80%', '100%'])
        else:
            plt.scatter(node, value, s=5, alpha=0.6)
            locs, labels = plt.xticks()
            start = locs[0]
            end = locs[-1]
            diff = end - start
            plt.xticks(
                [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
                ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.ylabel('Mean pdist')
        plt.xlabel('Nodes')
        plt.title(f'{genetype} - {filename} Mean pdist Scatter')
        save_to_scatter = os.path.join(save_to, 'mean pdist scatter')
        if not os.path.exists(save_to_scatter):
            os.makedirs(save_to_scatter)
        plt.savefig(os.path.join(save_to_scatter, f'{genetype}_{filename} Mean pdist Scatter.png'))
        plt.close('all')
        plt.figure(figsize=(12, 9))
        plt.hist(value, bins=50, edgecolor='black')
        plt.yscale('log')
        plt.ylabel('Count')
        plt.xlabel('Mean pdist')
        plt.title(f'{genetype} - {filename} Mean pdist Histogram')
        save_to_hist = os.path.join(save_to, 'mean pdist histogram')
        if not os.path.exists(save_to_hist):
            os.makedirs(save_to_hist)
        plt.savefig(os.path.join(save_to_hist, f'{genetype}_{filename} Mean pdist Histogram.png'))
        plt.close('all')

    def plot_pdist_mean_line(data, dictname, alpha, filename='', save_to=''):
        save_to = os.path.join(save_to, 'pdist', f'alpha = {alpha}')
        plt.figure(figsize=(12, 9))
        for i in range(len(data)):
            temp_data = data[i]
            temp_name = dictname[i]
            node = list(temp_data.keys())
            value = list(temp_data.values())
            x1 = np.linspace(0, 1, len(value))
            plt.plot(x1, value, label=temp_name)
        plt.legend()
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.xlabel('Node')
        plt.title(f'Target vs Non-Target Mean pdist Line ({filename})')
        save_to_line = os.path.join(save_to, 'mean pdist line')
        if not os.path.exists(save_to_line):
            os.makedirs(save_to_line)
        plt.savefig(os.path.join(save_to_line, f'{filename}_Mean pdist Comparison.png'))
        plt.close('all')

    for i in range(iter):
        alpha = (i + 1) / 10
        # for network_path in [network_target_path, erbc_target_path, tnbc_target_path]:
        for network_path in [network_target_path]:
            # for gene_ls in ['cancer_genes', 'ERBC_disease_nodes', 'TNBC_disease_nodes']:
            for gene_ls in ['cancer_genes']:
                for genetype in ['target', 'non-target', 'all_genes']:
                    with open(os.path.join(network_path, 'pdist', f'alpha = {alpha}', genetype,
                                           f'{gene_ls}_pdist_mean.txt'), 'r') as f:
                        temp_dict = f.read()
                    f.close()

                    temp_dict = ast.literal_eval(temp_dict)

                    plot_pdist_mean(temp_dict, network_path, gene_ls, genetype, alpha)

                with open(
                        os.path.join(network_path, 'pdist', f'alpha = {alpha}', 'target', f'{gene_ls}_pdist_mean.txt'),
                        'r') as f:
                    temp_dict = f.read()
                f.close()
                temp_target_dict = ast.literal_eval(temp_dict)
                with open(os.path.join(network_path, 'pdist', f'alpha = {alpha}', 'non-target',
                                       f'{gene_ls}_pdist_mean.txt'), 'r') as f:
                    temp_dict = f.read()
                f.close()
                temp_nontarget_dict = ast.literal_eval(temp_dict)
                with open(os.path.join(network_path, 'pdist', f'alpha = {alpha}', 'all_genes',
                                       f'{gene_ls}_pdist_mean.txt'), 'r') as f:
                    temp_dict = f.read()
                f.close()
                temp_all_dict = ast.literal_eval(temp_dict)
                data = [temp_target_dict, temp_nontarget_dict, temp_all_dict]
                plot_pdist_mean_line(data, ['target', 'non-target', 'all_genes'], alpha, gene_ls, network_path)
        print(f'plot of mean pdist for alpha = {alpha} finished')


# plot histogram and scatter for mean value of distance for both targets and non-targets
def plot_distance(network_target_path):
    def plot_distance_mean(mean_dict, save_to, filename, genetypes):
        save_to = os.path.join(save_to, 'distance', genetype)
        node = list(mean_dict.keys())
        value = list(mean_dict.values())
        plt.figure(figsize=(12, 9))
        if genetypes == 'all_genes':
            tarx = []
            tary=[]
            nonx = []
            nony=[]
            x0 = 0
            for u in node:
                x0 += 1
                if u in target_ls:
                    tarx.append(x0)
                    tary.append(mean_dict[u])
                else:
                    nonx.append(x0)
                    nony.append(mean_dict[u])

            plt.scatter(tarx,tary,s=5,alpha=0.7,c='red',marker='v',label='Targets')
            plt.scatter(nonx,nony,s=1,alpha=0.4,c='blue',marker='x',label='Non-Targets')
            plt.legend()
            start = 1
            end = len(node)
            diff = end - start
            plt.xticks(
                [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
                ['0', '20%', '40%', '60%', '80%', '100%'])
        else:
            plt.scatter(node, value, s=5, alpha=0.6)
            locs, labels = plt.xticks()
            start = locs[0]
            end = locs[-1]
            diff = end - start
            plt.xticks(
                [start + diff * 0, start + diff * 0.2, start + diff * 0.4, start + diff * 0.6, start + diff * 0.8, end],
                ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.ylabel('Mean distance')
        plt.xlabel('Nodes')
        plt.title(f'{genetype} - {filename} Mean Distance Scatter')
        save_to_scatter = os.path.join(save_to, 'mean distance scatter')
        if not os.path.exists(save_to_scatter):
            os.makedirs(save_to_scatter)
        plt.savefig(os.path.join(save_to_scatter, f'{genetype}_{filename} Mean Distance Scatter.png'))
        plt.close('all')
        plt.figure(figsize=(12, 9))
        plt.hist(value, bins=50, edgecolor='black')
        plt.yscale('log')

        plt.ylabel('Count')
        plt.xlabel('Mean distance')
        plt.title(f'{genetype} - {filename} Mean Distance Histogram')
        save_to_hist = os.path.join(save_to, 'mean distance histogram')
        if not os.path.exists(save_to_hist):
            os.makedirs(save_to_hist)
        plt.savefig(os.path.join(save_to_hist, f'{genetype}_{filename} Mean Distance Histogram.png'))
        plt.close('all')

    def plot_distance_mean_line(data_ls, dictname, filename='', save_to=''):
        save_to = os.path.join(save_to, 'distance')
        plt.figure(figsize=(12, 9))
        for i in range(len(data_ls)):
            temp_data = data_ls[i]
            temp_name = dictname[i]
            node = list(temp_data.keys())
            value = list(temp_data.values())
            x1 = np.linspace(0, 1, len(value))
            plt.plot(x1, value, label=temp_name)
        plt.legend()
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
        plt.xlabel('Node')
        plt.title(f'Target vs Non-Target Mean Distance Line ({filename})')
        save_to_line = os.path.join(save_to, 'mean distance line')
        if not os.path.exists(save_to_line):
            os.makedirs(save_to_line)
        plt.savefig(os.path.join(save_to_line, f'{filename}_Mean Distance Comparison.png'))
        plt.close('all')

    # for network_path in [network_target_path, erbc_target_path, tnbc_target_path]:
    for network_path in [network_target_path]:
        # for gene_ls in ['cancer_genes', 'ERBC_disease_nodes', 'TNBC_disease_nodes']:
        for gene_ls in ['cancer_genes']:
            for genetype in ['target', 'non-target', 'all_genes']:
                with open(os.path.join(network_path, 'distance', genetype,
                                       f'{gene_ls}_distance_mean.txt'), 'r') as f:
                    temp_dict = f.read()
                f.close()
                temp_target_dict = ast.literal_eval(temp_dict)

                plot_distance_mean(temp_target_dict, network_path, gene_ls, genetype)

            with open(os.path.join(network_path, 'distance', 'target', f'{gene_ls}_distance_mean.txt'),
                      'r') as f:
                temp_dict = f.read()
            f.close()
            temp_target_dict = ast.literal_eval(temp_dict)
            with open(os.path.join(network_path, 'distance', 'non-target',
                                   f'{gene_ls}_distance_mean.txt'), 'r') as f:
                temp_dict = f.read()
            f.close()
            temp_nontarget_dict = ast.literal_eval(temp_dict)
            with open(os.path.join(network_path, 'distance', 'all_genes',
                                   f'{gene_ls}_distance_mean.txt'), 'r') as f:
                temp_dict = f.read()
            f.close()
            temp_all_dict = ast.literal_eval(temp_dict)
            data = [temp_target_dict, temp_nontarget_dict, temp_all_dict]
            plot_distance_mean_line(data, ['target', 'non-target', 'all_genes'], gene_ls, network_path)
    print(f'plot of mean distance finished')


# conduct the two-sample ks test for targets and non-targets
def st_diff(target_path, filename,iter):
    def ks_test_target(target_data, non_target_data, statistics_name):
        target_value = target_data.loc[:, statistics_name]
        non_target_value = non_target_data.loc[:, statistics_name]
        from scipy import stats
        s, pvalue = stats.ks_2samp(target_value, non_target_value)
        if pvalue < 0.05:
            print(f'pvalue of {pvalue} is lower than the threshold 0.05, reject the null hypothesis.')
            print(f'Thus {statistics_name} for targets and non-targets in signaling network is statistically different\n')
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
    for gene_ls in ['cancer_genes']:
        with open(os.path.join(target_path, 'distance', 'target',
                               f'{gene_ls}_distance_mean.txt'), 'r') as f:
            temp_dict = f.read()
        f.close()
        temp_target_dict = ast.literal_eval(temp_dict)
        target_distance_arr = np.array(list(temp_target_dict.values()))
        with open(os.path.join(target_path, 'distance', 'non-target',
                               f'{gene_ls}_distance_mean.txt'), 'r') as f:
            temp_dict = f.read()
        f.close()
        temp_nontarget_dict = ast.literal_eval(temp_dict)
        non_distance_arr = np.array(list(temp_nontarget_dict.values()))
        from scipy import stats
        s, pvalue = stats.ks_2samp(target_distance_arr, non_distance_arr)
        if pvalue < 0.05:
            print(f'pvalue of {pvalue} is lower than the threshold 0.05, reject the null hypothesis.')
            print(f'Thus shortest distance for targets and non-targets in signaling network is statistically different\n')
        else:
            print(f'pvalue of {pvalue} is not lower than the threshold 0.05, we do not reject the null hypothesis.')
            print(f'Thus shortest distance for targets and non-targets in signaling network has no statistical '
                  f'difference\n')

    for i in range(iter):
        alpha = (i + 1) / 10
        for gene_ls in ['cancer_genes']:
            with open(os.path.join(target_path, 'pdist', f'alpha = {alpha}', 'target',
                                   f'{gene_ls}_pdist_mean.txt'), 'r') as f:
                temp_dict = f.read()
            f.close()
            temp_target_dict = ast.literal_eval(temp_dict)
            target_pdist_arr = np.array(list(temp_target_dict.values()))
            with open(os.path.join(target_path, 'pdist', f'alpha = {alpha}', 'non-target',
                                   f'{gene_ls}_pdist_mean.txt'), 'r') as f:
                temp_dict = f.read()
            f.close()
            temp_nontarget_dict = ast.literal_eval(temp_dict)
            non_pdist_arr = np.array(list(temp_nontarget_dict.values()))
            from scipy import stats
            s, pvalue = stats.ks_2samp(target_pdist_arr, non_pdist_arr)
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
