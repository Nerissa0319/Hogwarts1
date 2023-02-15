from constant import *
import pandas as pd
import ast
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random
import math
import statistics
from scipy.stats import spearmanr

'''Comparing topological features in different cancer subtypes'''

# read graph
whole_signaling = nx.read_gexf(os.path.join(output_path, 'whole_signaling.gexf'))
# read data from pdistance(alpha=0.2) and distance txt files
whole_pdist_df = pd.read_csv(f'{whole_pdist_path}\\alpha = 0.2\\pdist.txt', sep='\t', index_col=0, header=0)
whole_distance_df = pd.read_csv(f'{whole_st_path}\\distance.txt', sep='\t', index_col=0, header=0)


# find genes for cancer subtypes:
def find_subgene(network, tumour_type,cancer_name): # tumour_type: 'prostate', cancer_name:'Prostate Cancer'
    onco_ = []
    target_ = []
    # find oncogenes for each cancer subtype from cancer gene census dataframe ()
    for index, row in cancer_df.iterrows():
        if tumour_type.lower() in str(row['Tumour Types(Somatic)']).lower() or tumour_type.lower() in str(row['Tumour Types(Germline)']).lower():
            if row['Gene Symbol'] in network.nodes():  # remove those not in signaling network
                onco_.append(row['Gene Symbol'])

    # find target genes for each cancer subtype
    for index, row in target_df.iterrows():
        if cancer_name.lower() in str(row['Indications']).lower():
            temp = str(row['Targets']).split('; ')
            for gene in temp:
                if gene in network.nodes():
                    target_.append(gene)

    target_ = list(set(target_))  # remove redundant elements from the list
    nontarget_ = []  # genes in signaling network which are not prostate cancer targets
    for u in network.nodes():
        if not u in target_:
            nontarget_.append(u)
    return onco_,target_,nontarget_


# extract topological features for different cancer subtypes
def find_feature(breast_onco_, prostate_onco_, breast_target_, prostate_target_, breast_nontarget_,
                 prostate_nontarget_):
    # create a dataframe to store pdist and distance for cancer subtype
    prostate_target_pdist_ = pd.DataFrame(index=pd.Index(prostate_target_), columns=pd.Index(prostate_onco_))
    prostate_non_pdist_ = pd.DataFrame(index=pd.Index(prostate_nontarget_), columns=pd.Index(prostate_onco_))
    breast_target_pdist_ = pd.DataFrame(index=pd.Index(breast_target_), columns=pd.Index(breast_onco_))
    breast_non_pdist_ = pd.DataFrame(index=pd.Index(breast_nontarget_), columns=pd.Index(breast_onco_))

    prostate_target_distance_ = pd.DataFrame(index=pd.Index(prostate_target_), columns=pd.Index(prostate_onco_))
    prostate_non_distance_ = pd.DataFrame(index=pd.Index(prostate_nontarget_), columns=pd.Index(prostate_onco_))
    breast_target_distance_ = pd.DataFrame(index=pd.Index(breast_target_), columns=pd.Index(breast_onco_))
    breast_non_distance_ = pd.DataFrame(index=pd.Index(breast_nontarget_), columns=pd.Index(breast_onco_))
    # assign values to dataframes
    for columns in prostate_onco_:
        for index in prostate_target_:
            prostate_target_pdist_.loc[index, columns] = whole_pdist_df.loc[index, columns]
            prostate_target_distance_.loc[index, columns] = whole_distance_df.loc[index, columns]
        for index in prostate_nontarget_:
            prostate_non_pdist_.loc[index, columns] = whole_pdist_df.loc[index, columns]
            prostate_non_distance_.loc[index, columns] = whole_distance_df.loc[index, columns]
    for columns in breast_onco_:
        for index in breast_target_:
            breast_target_pdist_.loc[index, columns] = whole_pdist_df.loc[index, columns]
            breast_target_distance_.loc[index, columns] = whole_distance_df.loc[index, columns]
        for index in breast_nontarget_:
            breast_non_pdist_.loc[index, columns] = whole_pdist_df.loc[index, columns]
            breast_non_distance_.loc[index, columns] = whole_distance_df.loc[index, columns]
    # write the dataframes to txt files
    prostate_target_distance_.to_csv(f'{cancer_comparison_path}//prostate_cancer_target_distance.txt', sep='\t',
                                     header=True,
                                     index=True)
    prostate_target_pdist_.to_csv(f'{cancer_comparison_path}//prostate_cancer_target_pdistance.txt', sep='\t',
                                  header=True,
                                  index=True)
    prostate_non_distance_.to_csv(f'{cancer_comparison_path}//prostate_cancer_nontarget_distance.txt', sep='\t',
                                  header=True,
                                  index=True)
    prostate_non_pdist_.to_csv(f'{cancer_comparison_path}//prostate_cancer_nontarget_pdistance.txt', sep='\t',
                               header=True,
                               index=True)
    breast_target_distance_.to_csv(f'{cancer_comparison_path}//breast_cancer_target_distance.txt', sep='\t',
                                   header=True,
                                   index=True)
    breast_target_pdist_.to_csv(f'{cancer_comparison_path}//breast_cancer_target_pdistance.txt', sep='\t', header=True,
                                index=True)
    breast_non_distance_.to_csv(f'{cancer_comparison_path}//breast_cancer_nontarget_distance.txt', sep='\t',
                                header=True,
                                index=True)
    breast_non_pdist_.to_csv(f'{cancer_comparison_path}//breast_cancer_nontarget_pdistance.txt', sep='\t', header=True,
                             index=True)

    return


# use boxplot of distance and pdistance to compare targets and non-targets in cancer subtypes
# target/nontarget_df_ is the distance/pdist dataframe for targets/nontargets,
# cancertype is 'Prostate Cancer'/'Breast Cancer', featurename is 'Distance'/'Pdistance'
def box_plot(target_df_, nontarget_df_, cancertype_, featurename_, targets):
    # we remove data of all the pairs of nodes which have no connections when plotting boxplots
    # that is, remove all the pairs of nodes whose distance = 0 or pdistance = 12.5129
    if featurename_ == 'PDistance':
        target_df_ = target_df_.replace(12.5129,
                                        np.NaN)  # replace 12.5129 with NaN so that the data will not be plotted
        nontarget_df_ = nontarget_df_.replace(12.5129, np.NaN)
    # we remove data of all the (drug target, oncogene) pairs where the target and the oncogene are the same
    for df in [target_df_, nontarget_df_]:
        for i in df.index:
            for c in df.columns:
                if i == c:  # if target name == oncogene name, replace the value with NaN
                    df.replace(df.loc[i, c], np.NaN, inplace=True)
    # record the percentage of NaN values (the percentage of targets/nontargets with no connections to oncogenes)
    # and the percentage of outliers (values > Q3 + 1.5IQR)
    target_stat = pd.DataFrame(columns=['NaN', 'NaN %', 'outliers', 'outlier %'], index=target_df_.columns)
    non_stat = pd.DataFrame(columns=['NaN', 'NaN %', 'outliers', 'outlier %'], index=nontarget_df_.columns)
    # in target genes
    for index, row in target_df_.transpose().iterrows():
        nan_no = row.isna().sum()  # count the number of NaN values
        target_stat.loc[index, 'NaN'] = nan_no
        iqr = row.quantile(0.75) - row.quantile(0.25)  # IQR(Interquartile range) = third Quartile - first Quartile
        upper = min(row.quantile(0.75) + 1.5 * iqr, row.max())  # upper limit of the boxplot = 3rd quartile + 1.5iqr
        lower = max(row.quantile(0.25) - 1.5 * iqr, row.min())  # lower limit of the boxplot = 1st quartile - 1.5iqr
        # datas beyond the upper and lower limit are called outliers
        outlier_no = row[(row < lower) | (row > upper)].count()  # count the number of outliers
        target_stat.loc[index, 'outliers'] = outlier_no
        target_stat.loc[index, 'NaN %'] = '{:.2%}'.format(nan_no / len(row))
        target_stat.loc[index, 'outlier %'] = '{:.2%}'.format(outlier_no / (len(row) - nan_no))
    # in nontarget genes
    for index, row in nontarget_df_.transpose().iterrows():
        nan_no = row.isna().sum()
        non_stat.loc[index, 'NaN'] = nan_no
        iqr = row.quantile(0.75) - row.quantile(0.25)  # IQR(Interquartile range) = third Quartile - first Quartile
        upper = min(row.quantile(0.75) + 1.5 * iqr, row.max())  # upper limit of the boxplot = 3rd quartile + 1.5iqr
        lower = max(row.quantile(0.25) - 1.5 * iqr, row.min())  # lower limit of the boxplot = 1st quartile - 1.5iqr
        # datas beyond the upper and lower limit are called outliers
        outlier_no = row[(row < lower) | (row > upper)].count()  # count the number of outliers
        non_stat.loc[index, 'outliers'] = outlier_no
        non_stat.loc[index, 'NaN %'] = '{:.2%}'.format(nan_no / len(row))
        non_stat.loc[index, 'outlier %'] = '{:.2%}'.format(outlier_no / (len(row) - nan_no))
    # write the summary of NaN values and outliers to txt files
    target_stat.to_csv(f'{cancer_comparison_path}//{cancertype_} targets {featurename_} statistics.txt',
                       sep='\t', header=True, index=True)
    non_stat.to_csv(f'{cancer_comparison_path}//{cancertype_} non targets {featurename_} statistics.txt',
                    sep='\t', header=True, index=True)

    whole_df = pd.concat([target_df_, nontarget_df_])  # concatenate the both targets and nontargets dataframes into one
    # the argument of cancertype is 'Prostate Cancer'/'Breast Cancer'
    # choose the data according to the cancertype argument

    for u in whole_df.index:  # in whole_df, add a columns called 'type' to record whether the gene is a target
        if u in targets:
            whole_df.loc[u, 'type'] = 'Targets'
        else:
            whole_df.loc[u, 'type'] = 'Non-Targets'
    # start plotting
    plt.figure(figsize=(16, 9))
    sns.stripplot(x='variable', y='value', hue='type', data=pd.melt(whole_df, id_vars='type'),
                  jitter=0.1, alpha=0.5, size=2, dodge=True)  # the scatter plot of all the values
    ax = sns.boxplot(x='variable', y='value', hue='type', data=pd.melt(whole_df, id_vars='type'),
                     medianprops=dict(color='red'),
                     showmeans=True, meanline=True, meanprops={'color': 'red', 'ls': ':'})  # boxplot

    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:2], labels[0:2], loc='upper center', ncol=2, frameon=False)  # legend
    plt.xlabel('oncogenes')
    plt.ylabel(featurename_)
    plt.xticks(rotation=30)  # rotate the xtick labels
    plt.title(f'{featurename_} for {cancertype_}')
    plt.tight_layout()
    plt.savefig(f'{cancer_comparison_path}//{featurename_} for {cancertype_}')
    plt.close()


# scatterplot of distance/pdistance mean vs other topolotical features(degree,eigenvector centrality...)
def scatter_plot(target_df_, nontarget_df_, featurename_, cancertype_):
    # we remove data of all the pairs of nodes which have no connections when plotting boxplots
    # that is, remove all the pairs of nodes whose distance = 0 or pdistance = 12.5129
    if featurename_ == 'PDistance':
        target_df_ = target_df_.replace(12.5129,
                                        np.NaN)  # replace 12.5129 with NaN so that the data will not be plotted
        nontarget_df_ = nontarget_df_.replace(12.5129, np.NaN)
    # we remove data of all the (drug target, oncogene) pairs where the target and the oncogene are the same
    for df in [target_df_, nontarget_df_]:
        for i in df.index:
            for c in df.columns:
                if i == c:  # if target name == oncogene name, replace the value with NaN
                    df.replace(df.loc[i, c], np.NaN, inplace=True)
    # compute the average pdistance or distance
    target_mean = target_df_.mean(axis=1)  # compute mean value for distance and pdistance
    non_mean = nontarget_df_.mean(axis=1)
    # read topological features (degree centrality ...) from txt files
    all_info = pd.read_csv(f'{whole_target_path}\\All Genes Info_whole_signaling.txt', sep='\t', header=0, index_col=0)
    # for each feature, output a scatter plot
    for centrality in ['Betweenness Centrality', 'Degree Centrality',
                       'Pagerank', 'Eigenvector Centrality', 'Closeness Centrality']:
        target_centrality = all_info.loc[target_mean.index, centrality]
        non_centrality = all_info.loc[non_mean.index, centrality]
        fig = plt.figure()
        plt.scatter(target_mean, target_centrality)  # plot scatter for targets
        plt.title(f'{featurename_} vs {centrality} for {cancertype_} targets')
        plt.ylabel(f'{centrality}')
        plt.xlabel(f'{featurename_}')
        plt.savefig(f'{cancer_comparison_path}/{featurename_} vs {centrality} for {cancertype_} targets')
        plt.close()
        fig2 = plt.figure()
        plt.scatter(non_mean, non_centrality)  # plot scatter for nontargets
        plt.title(f'{featurename_} vs {centrality} for {cancertype_} non-targets')
        plt.ylabel(f'{centrality}')
        plt.xlabel(f'{featurename_}')
        plt.savefig(f'{cancer_comparison_path}/{featurename_} vs {centrality} for {cancertype_} non-targets')
        plt.close()


# plot scatter plot of pdistance vs distance for targets and non targets in cancer subtype
def pdist_vs_distance(target_pdist_df, non_pdist_df, target_distance_df, non_distance_df, cancertype):
    target_pdist_df = target_pdist_df.replace(12.5129, np.NaN)
    non_pdist_df = non_pdist_df.replace(12.5129, np.NaN)
    for df in [target_pdist_df, non_pdist_df, target_distance_df, non_distance_df]:
        for i in df.index:
            for c in df.columns:
                if i == c:
                    df.replace(df.loc[i, c], np.NaN, inplace=True)
    target_pdist_mean = target_pdist_df.mean(axis=1)  # compute mean values
    non_pdist_mean = non_pdist_df.mean(axis=1)  # compute mean values for pdist
    target_distance_mean = target_distance_df.mean(axis=1)  # compute mean values
    non_distance_mean = non_distance_df.mean(axis=1)  # compute mean values for distance
    plt.figure()
    plt.scatter(target_distance_mean, target_pdist_mean)  # plot distance vs pdistance for targets
    plt.ylabel('pdistance')
    plt.xlabel('distance')
    plt.title(f'pdistance vs distance for {cancertype} targets')
    plt.savefig(f'{cancer_comparison_path}/pdistance vs distance for {cancertype} targets')  # save plot
    plt.close()
    plt.figure()
    plt.scatter(non_distance_mean, non_pdist_mean)  # plot distance vs pdistance for nontargets
    plt.ylabel('pdistance')
    plt.xlabel('distance')
    plt.title(f'pdistance vs distance for {cancertype} non-targets')
    plt.savefig(f'{cancer_comparison_path}/pdistance vs distance for {cancertype} non-targets')  # save plot
    plt.close()
