from constant import *
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

'''This code analyzed the hallmarks of drug targets in order to help prune candidate targets'''

# clean the drug targets data
drug_targets = target_df.copy()
target_raw = list(target_df['Targets'])
targets_in_list = []
indications = []
for tg in target_raw:
    if type(tg) == str:
        temp = tg.split('; ')
        targets_in_list.append(temp)
    else:
        targets_in_list.append(tg)
for ind in list(target_df['Indications']):
    if type(ind) == str:
        temp = ind.split('; ')
        indications.append(temp)
    else:
        indications.append(ind)
del drug_targets['Targets'] # delete the original Targets column
del drug_targets['Indications'] # delete the original Indications column
drug_targets['Indications'] = indications # assign the new indications column
drug_targets['Targets'] = targets_in_list # assign the new targets column
drug_targets.dropna(subset=['Targets'],inplace=True) # delete the rows with no targets
drug_targets.set_index(['Product'],inplace=True)

hallmark_names = ['promote proliferation',
                  'evade growth repressor',
                  'resist cell death',
                  'enable replicative immortality',
                  'Induce angiogenesis',
                  'Activate invasion & metastasis',
                  'Reprogramme energy metabolism',
                  'Evade immune destruction',
                  'Genome instability',
                  'Tumor promoting inflammation']

# drug hallmarks summary
drug_columns = hallmark_names
drug_hallmarks = pd.DataFrame(0,index=list(drug_targets.index),columns = drug_columns)
hallMarks_df = pd.read_csv(os.path.join(hallmark_path, 'AllHallmarks.csv'),index_col=0,header=0)
no_of_targets = []
if not os.path.exists(f'{whole_target_path}/drug_targets_hallmarks.csv'):
    eff_targets = [0]*len(drug_targets.index)
    ind = 0
    indi = []
    # drug's hallmark score = no. of targets with this hallmark/no. of targets of the drug
    for drug in drug_targets.index:
        for target in drug_targets.loc[drug]['Targets']:
            if target in hallMarks_df.index:
                eff_targets[ind] += 1 # eff_targets is the genes that can be found in our dataset
                for h in drug_columns:
                    if hallMarks_df.loc[target][h] == 1:
                        drug_hallmarks.loc[drug][h] += 1 # compute the frequency of each hallmark
        no_of_targets.append(eff_targets[ind])
        indi.append(drug_targets.loc[drug,'Indications'])
        ind += 1

    # drug_hallmarks = drug_hallmarks.div(eff_targets,axis='index') # divide frequency by no. of targets
    drug_hallmarks['No. of Targets'] = no_of_targets
    sum_score = []
    for row in drug_hallmarks.iterrows():
        sum_score.append(sum(row[1][0:9]))
    drug_hallmarks['Sum of Score'] = sum_score
    drug_hallmarks['Indications'] = indi
    pd.DataFrame.to_csv(drug_hallmarks,f'{whole_target_path}/drug_targets_hallmarks.csv') # save as a file

else:
    # drug_hallmarks is a dataframe with all the drugs for cancer and their hallmark score
    drug_hallmarks = pd.read_csv(f'{whole_target_path}/drug_targets_hallmarks.csv',sep=',',
                                 index_col=0,header=0)

# target hallmarks summary
# target_hallmarks is the dataframe of each drug_targets and their hallmarks
if not os.path.exists(f'{whole_target_path}/target_hallmarks.csv'):
    target_hallmarks = pd.DataFrame(index=target_ls,columns=hallmark_names)
    for u in target_ls:
        if u in hallMarks_df.index:
            target_hallmarks.loc[u] = hallMarks_df.loc[u]
    target_hallmarks.dropna(inplace=True)
    drug_of_targets = []
    for u in target_hallmarks.index:
        temp = []
        for d in drug_targets.index:
            if u in drug_targets.loc[d]['Targets']:
                temp.append(d)
        drug_of_targets.append(temp)
    target_hallmarks['Drugs'] = drug_of_targets
    target_hallmarks.to_csv(f'{whole_target_path}/target_hallmarks.csv')
else:
    target_hallmarks = pd.read_csv(f'{whole_target_path}/target_hallmarks.csv',sep=',',
                                   index_col=0,header=0)


# analyzes the patterns of hallmarks for drug targets
def targets_analysis(target_hallmarks_df):
    hallmark_combo = {}
    for target in target_hallmarks_df.index:
        temp_str = ''
        for h in target_hallmarks_df.columns[:-1]:
            temp_str = temp_str + str(target_hallmarks_df.loc[target,h])
        if temp_str not in hallmark_combo.keys():
            hallmark_combo[temp_str] = 1
        else:
            hallmark_combo[temp_str] += 1

    hallmark_combo = {k:v for k,v in sorted(hallmark_combo.items(),reverse=True,key=lambda item:item[1])}
    # target_hallmark_combo is a dataframe with each drug_target and its 10-digit binary code of hallmarks
    target_hallmark_combo = pd.DataFrame.from_dict(hallmark_combo,orient='index')

    return target_hallmark_combo


