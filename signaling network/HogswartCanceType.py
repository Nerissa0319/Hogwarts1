from constant import *
import pandas as pd
import ast
import matplotlib.pyplot as plt
import numpy as np
import random
import math
whole_signaling = nx.read_gexf(os.path.join(output_path, 'whole_signaling.gexf'))
target_in_network = []
nontar_in_network = []
for u in whole_signaling.nodes():
    if u in target_ls:
        target_in_network.append(u)
    else:
        nontar_in_network.append(u)


# find genes for cancer subtypes:
def find_subgene(network):
    prostate_onco = []
    breast_onco = []
    prostate_target = []
    breast_target = []
    # extract oncogenes for specific cancer types
    for index, row in cancer_df.iterrows():
        if 'prostate' in str(row['Tumour Types(Somatic)']) or 'prostate' in str(row['Tumour Types(Germline)']):
            if row['Gene Symbol'] in network.nodes():  # remove those not in signaling network
                prostate_onco.append(row['Gene Symbol'])
        if 'breast' in str(row['Tumour Types(Somatic)']) or 'prostate' in str(row['Tumour Types(Germline)']):
            if row['Gene Symbol'] in network.nodes():  # remove those not in signaling network
                breast_onco.append(row['Gene Symbol'])
    for index, row in target_df.iterrows():
        if 'Prostate Cancer' in str(row['Indications']) or 'prostate cancer' in str(row['Indications']):
            temp = str(row['Targets']).split('; ')
            for gene in temp:
                if gene in network.nodes():
                    prostate_target.append(gene)
        if 'Breast Cancer' in str(row['Indications']) or 'breast cancer' in str(row['Indications']):
            temp = str(row['Targets']).split('; ')
            for gene in temp:
                if gene in network.nodes():
                    breast_target.append(gene)
    prostate_target = list(set(prostate_target))
    breast_target = list(set(breast_target))

    return breast_onco, prostate_onco, breast_target, prostate_target


bc_onco, pro_onco, bc_tar, pro_tar = find_subgene(whole_signaling)


def find_feature():
    # create a dataframe to store pdist value between targets/non-targets and prostate oncogenes
    pro_target_pdist = pd.DataFrame(index=pd.Index(target_in_network), columns=pd.Index(pro_onco))
    pro_non_pdist = pd.DataFrame(index=pd.Index(nontar_in_network), columns=pd.Index(pro_onco))
    breast_target_pdist = pd.DataFrame(index=pd.Index(target_in_network), columns=pd.Index(bc_onco))
    breast_non_pdist = pd.DataFrame(index=pd.Index(nontar_in_network), columns=pd.Index(bc_onco))

    pro_target_distance = pd.DataFrame(index=pd.Index(target_in_network), columns=pd.Index(pro_onco))
    pro_non_distance = pd.DataFrame(index=pd.Index(nontar_in_network), columns=pd.Index(pro_onco))
    breast_target_distance = pd.DataFrame(index=pd.Index(target_in_network), columns=pd.Index(bc_onco))
    breast_non_distance = pd.DataFrame(index=pd.Index(nontar_in_network), columns=pd.Index(bc_onco))
    whole_pdist_df = pd.read_csv(f'{whole_pdist_path}\\alpha = 0.2\\pdist.txt', sep='\t', index_col=0, header=0)
    whole_distance_df = pd.read_csv(f'{whole_st_path}\\distance.txt', sep='\t', index_col=0, header=0)

    for columns in pro_onco:
        for index in target_in_network:
            pro_target_pdist.loc[index, columns] = whole_pdist_df.loc[index, columns]
            pro_target_distance.loc[index, columns] = whole_distance_df.loc[index, columns]
        for index in nontar_in_network:
            pro_non_pdist.loc[index, columns] = whole_pdist_df.loc[index, columns]
            pro_non_distance.loc[index, columns] = whole_distance_df.loc[index, columns]
    for columns in bc_onco:
        for index in target_in_network:
            breast_target_pdist.loc[index, columns] = whole_pdist_df.loc[index, columns]
            breast_target_distance.loc[index, columns] = whole_distance_df.loc[index, columns]
        for index in nontar_in_network:
            breast_non_pdist.loc[index, columns] = whole_pdist_df.loc[index, columns]
            breast_non_distance.loc[index, columns] = whole_distance_df.loc[index, columns]

    # breast_target_distance.replace(0,np.NaN,inplace=True)
    # breast_non_distance.replace(0,np.NaN,inplace=True)
    # pro_target_distance.replace(0,np.NaN,inplace=True)
    # pro_non_distance.replace(0,np.NaN,inplace=True)
    # breast_target_pdist.replace(12.5129,np.NaN,inplace=True)
    # breast_non_pdist.replace(12.5129,np.NaN,inplace=True)
    # pro_target_pdist.replace(12.5129,np.NaN,inplace=True)
    # pro_non_pdist.replace(12.5129,np.NaN,inplace=True)

    return breast_target_pdist, breast_non_pdist, pro_target_pdist, pro_non_pdist, breast_target_distance, breast_non_distance, pro_target_distance, pro_non_distance


bc_target_pdist, bc_non_pdist, prostate_target_pdist, prostate_non_pdist, bc_target_distance, bc_non_distance, prostate_target_distance, prostate_non_distance = find_feature()


def pdist_vs_distance():
    # prostate cancer
    tar_pdist_mean = []
    non_pdist_mean = []
    tar_distance_mean = []
    non_distance_mean = []
    for index, row in prostate_target_pdist.iterrows():
        tar_pdist_mean.append(row.mean())
    sorted_tar_pdist_mean = sorted(tar_pdist_mean,key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    for index, row in prostate_non_pdist.iterrows():
        non_pdist_mean.append(row.mean())
    sorted_non_pdist_mean = sorted(non_pdist_mean, key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    plt.figure(figsize=(12, 9))
    x1 = np.linspace(0, 1, len(sorted_tar_pdist_mean))
    plt.plot(x1, sorted_tar_pdist_mean, label='targets')
    x2 = np.linspace(0, 1, len(sorted_non_pdist_mean))
    plt.plot(x2, sorted_non_pdist_mean, label='non-targets')
    plt.legend()
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
    plt.xlabel('Node')
    plt.title(f'PDist to Prostate Cancer Oncogenes')
    plt.savefig(f'{cancer_comparison_path}/PDist to Prostate Cancer Oncogenes.png')
    plt.close('all')
    for index, row in prostate_target_distance.iterrows():
        tar_distance_mean.append(row.mean())
    sorted_tar_distance_mean = sorted(tar_distance_mean, key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    for index, row in prostate_non_distance.iterrows():
        non_distance_mean.append(row.mean())
    sorted_non_distance_mean = sorted(non_distance_mean, key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    plt.figure(figsize=(12, 9))
    x1 = np.linspace(0, 1, len(sorted_tar_distance_mean))
    plt.plot(x1, sorted_tar_distance_mean, label='targets')
    x2 = np.linspace(0, 1, len(sorted_non_distance_mean))
    plt.plot(x2, sorted_non_distance_mean, label='non-targets')
    plt.legend()
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
    plt.xlabel('Node')
    plt.title(f'Distance to Prostate Cancer Oncogenes')
    plt.savefig(f'{cancer_comparison_path}/Distance to Prostate Cancer Oncogenes.png')
    plt.close('all')
    plt.figure(figsize=(12, 9))
    plt.scatter(tar_distance_mean, tar_pdist_mean, label='target')
    n = range(len(nontar_in_network))
    n_list = [*n]
    random_pos = random.sample(n_list,len(target_in_network))
    distance_sample = []
    pdist_sample = []
    for i in random_pos:
        distance_sample.append(non_distance_mean[i])
        pdist_sample.append(non_pdist_mean[i])
    plt.scatter(distance_sample, pdist_sample, label='non-target')
    plt.ylabel('Pdist')
    plt.xlabel('Distance')
    plt.title('Prostate Cancer')
    plt.legend()
    plt.savefig(f'{cancer_comparison_path}/Pdist vs Distance (Prostate Cancer).png')

    # breast cancer
    tar_pdist_mean = []
    non_pdist_mean = []
    tar_distance_mean = []
    non_distance_mean = []
    for index, row in bc_target_pdist.iterrows():
        tar_pdist_mean.append(row.mean())
    sorted_tar_pdist_mean = sorted(tar_pdist_mean, key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    for index, row in bc_non_pdist.iterrows():
        non_pdist_mean.append(row.mean())
    sorted_non_pdist_mean = sorted(non_pdist_mean, key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    plt.figure(figsize=(12, 9))
    x1 = np.linspace(0, 1, len(sorted_tar_pdist_mean))
    plt.plot(x1, sorted_tar_pdist_mean, label='targets')
    x2 = np.linspace(0, 1, len(sorted_non_pdist_mean))
    plt.plot(x2, sorted_non_pdist_mean, label='non-targets')
    plt.legend()
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
    plt.xlabel('Node')
    plt.title(f'PDist to Breast Cancer Oncogenes')
    plt.savefig(f'{cancer_comparison_path}/PDist to Breast Cancer Oncogenes.png')
    plt.close('all')
    for index, row in bc_target_distance.iterrows():
        tar_distance_mean.append(row.mean())
    sorted_tar_distance_mean = sorted(tar_distance_mean, key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    for index, row in bc_non_distance.iterrows():
        non_distance_mean.append(row.mean())
    sorted_non_distance_mean = sorted(non_distance_mean, key=lambda f: float('-inf') if math.isnan(f) else f,reverse=True)
    plt.figure(figsize=(12, 9))
    x1 = np.linspace(0, 1, len(sorted_tar_distance_mean))
    plt.plot(x1, sorted_tar_distance_mean, label='targets')
    x2 = np.linspace(0, 1, len(sorted_non_distance_mean))
    plt.plot(x2, sorted_non_distance_mean, label='non-targets')
    plt.legend()
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '20%', '40%', '60%', '80%', '100%'])
    plt.xlabel('Node')
    plt.title(f'Distance to Breast Cancer Oncogenes')
    plt.savefig(f'{cancer_comparison_path}/Distance to Breast Cancer Oncogenes.png')
    plt.close('all')
    plt.figure(figsize=(12, 9))
    plt.scatter(tar_distance_mean, tar_pdist_mean, label='target')
    n = range(len(nontar_in_network))
    n_list = [*n]
    random_pos = random.sample(n_list, len(target_in_network))
    distance_sample = []
    pdist_sample = []
    for i in random_pos:
        distance_sample.append(non_distance_mean[i])
        pdist_sample.append(non_pdist_mean[i])
    plt.scatter(distance_sample, pdist_sample, label='non-target')
    plt.ylabel('Pdist')
    plt.xlabel('Distance')
    plt.title('Breast Cancer')
    plt.legend()
    plt.savefig(f'{cancer_comparison_path}/Pdist vs Distance (Breast Cancer).png')

target_info = pd.read_csv(r'D:\Hogswart\Output\Signaling Network\target\whole_signaling\Targeted Genes '
                          r'Info_whole_signaling.txt',sep='\t',header=0,index_col=0)
nontar_info = pd.read_csv(r'D:\Hogswart\Output\Signaling Network\target\whole_signaling\Non-Targeted Genes '
                          r'Info_whole_signaling.txt',sep='\t',header=0,index_col=0)
def pdist_vs_other():


    # prostate cancer
    tar_pdist_mean = []
    non_pdist_mean = []
    for index, row in prostate_target_pdist.iterrows():
        tar_pdist_mean.append(row.mean())
    for index, row in prostate_non_pdist.iterrows():
        non_pdist_mean.append(row.mean())
    for feature in ['Degree Centrality','Eigenvector Centrality','Closeness Centrality','Betweenness Centrality','Pagerank']:
        plt.figure()
        plt.scatter(tar_pdist_mean,target_info.loc[prostate_target_pdist.index.tolist(),'Degree Centrality'],label='Targets')
        plt.scatter(non_pdist_mean,nontar_info.loc[prostate_non_pdist.index.tolist(),'Degree Centrality'],label='Non Targets',alpha=0.2)
        plt.xlabel('Pdist')
        plt.ylabel(f'{feature}')
        plt.title(f'pdist vs {feature} (Prostate Cancer)')
        plt.legend()
        plt.savefig(f'D:\\Hogswart\\Output\\Signaling Network\\cancer_subtype\\pdist vs {feature} (Prostate Cancer).png')
        plt.close()

    # breast cancer
    tar_pdist_mean = []
    non_pdist_mean = []
    for index, row in bc_target_pdist.iterrows():
        tar_pdist_mean.append(row.mean())
    for index, row in bc_non_pdist.iterrows():
        non_pdist_mean.append(row.mean())
    for feature in ['Degree Centrality', 'Eigenvector Centrality', 'Closeness Centrality', 'Betweenness Centrality',
                    'Pagerank']:
        plt.figure()
        plt.scatter(tar_pdist_mean, target_info.loc[bc_target_pdist.index.tolist(), 'Degree Centrality'],
                    label='Targets')
        plt.scatter(non_pdist_mean, nontar_info.loc[bc_non_pdist.index.tolist(), 'Degree Centrality'],
                    label='Non Targets', alpha=0.2)
        plt.xlabel('Pdist')
        plt.ylabel(f'{feature}')
        plt.title(f'pdist vs {feature} (Breast Cancer)')
        plt.legend()
        plt.savefig(f'D:\\Hogswart\\Output\\Signaling Network\\cancer_subtype\\pdist vs {feature} (Breast Cancer).png')
        plt.close()
def distance_vs_other():
    # prostate cancer
    tar_distance_mean = []
    non_distance_mean = []
    for index, row in prostate_target_distance.iterrows():
        tar_distance_mean.append(row.mean())
    for index, row in prostate_non_distance.iterrows():
        non_distance_mean.append(row.mean())

    for feature in ['Degree Centrality', 'Eigenvector Centrality', 'Closeness Centrality', 'Betweenness Centrality',
                    'Pagerank']:
        plt.figure()
        plt.scatter(tar_distance_mean, target_info.loc[prostate_target_distance.index.tolist(), 'Degree Centrality'],
                    label='Targets')
        plt.scatter(non_distance_mean, nontar_info.loc[prostate_non_distance.index.tolist(), 'Degree Centrality'],
                    label='Non Targets', alpha=0.2)
        plt.xlabel('Distance')
        plt.ylabel(f'{feature}')
        plt.title(f'Distance vs {feature} (Prostate Cancer)')
        plt.legend()
        plt.savefig(f'D:\\Hogswart\\Output\\Signaling Network\\cancer_subtype\\Distance vs {feature} (Prostate Cancer).png')

    # Breast Cancer
    tar_distance_mean = []
    non_distance_mean = []
    for index, row in bc_target_distance.iterrows():
        tar_distance_mean.append(row.mean())
    for index, row in bc_non_distance.iterrows():
        non_distance_mean.append(row.mean())

    for feature in ['Degree Centrality', 'Eigenvector Centrality', 'Closeness Centrality', 'Betweenness Centrality',
                    'Pagerank']:
        plt.figure()
        plt.scatter(tar_distance_mean, target_info.loc[bc_target_distance.index.tolist(), 'Degree Centrality'],
                    label='Targets')
        plt.scatter(non_distance_mean, nontar_info.loc[bc_non_distance.index.tolist(), 'Degree Centrality'],
                    label='Non Targets', alpha=0.2)
        plt.xlabel('Distance')
        plt.ylabel(f'{feature}')
        plt.title(f'Distance vs {feature} (Breast Cancer)')
        plt.legend()
        plt.savefig(
            f'D:\\Hogswart\\Output\\Signaling Network\\cancer_subtype\\Distance vs {feature} (Breast Cancer).png')
