# this file list all the constant, dataset, paths we will use in the whole project
# so that we will  not have to write the whole string every time when we want to use it

import os
import pandas as pd
import networkx as nx

# list the path we will use in the project
project_path = os.path.dirname(os.path.dirname(os.getcwd()))
output_path = os.path.join(project_path, 'Output', 'Neutral')  # output path
input_path = os.path.join(project_path, 'Input')  # input path
common_path = os.path.join(input_path, 'common files')  # the dir of common files
hallmark_path = os.path.join(common_path, 'node', 'hallmark annotation')
ERBC_input_path = os.path.join(input_path, 'ER positive breast cancer')
TNBC_input_path = os.path.join(input_path, 'Triple negative breast cancer')
protease_output_path = os.path.join(input_path, 'PROTEASE output')
pathway_path = os.path.join(output_path, 'pathway')
pathway_txt_path = os.path.join(pathway_path,'pathways')
visualization_path = os.path.join(output_path, 'visualization')
stat_whole_path = os.path.join(output_path, 'statistics', 'whole_neutral')
stat_chart_whole_path = os.path.join(stat_whole_path, 'charts')
# stat_erbc_path = os.path.join(output_path,'statistics','erbc_neutral')
# stat_chart_erbc_path = os.path.join(stat_erbc_path,'charts')
# stat_tnbc_path = os.path.join(output_path,'statistics','tnbc_neutral')
# stat_chart_tnbc_path = os.path.join(stat_tnbc_path,'charts')
whole_ppr_path = os.path.join(output_path, 'ppr', 'whole_ppr')
# erbc_ppr_path = os.path.join(output_path,'ppr','erbc_ppr')
# tnbc_ppr_path = os.path.join(output_path,'ppr','tnbc_ppr')
whole_st_path = os.path.join(output_path, 'distance', 'whole_distance')
# erbc_st_path = os.path.join(output_path,'distance','erbc_distance')
# tnbc_st_path = os.path.join(output_path,'distance','tnbc_distance')
whole_target_path = os.path.join(output_path, 'target', 'whole_neutral')
# erbc_target_path = os.path.join(output_path,'target','erbc_neutral')
# tnbc_target_path = os.path.join(output_path,'target','tnbc_neutral')
whole_dppr_path = os.path.join(output_path, 'dppr', 'whole_neutral')
# erbc_dppr_path = os.path.join(output_path,'dppr','erbc_neutral')
# tnbc_dppr_path = os.path.join(output_path,'dppr','tnbc_neutral')
whole_pdist_path = os.path.join(output_path, 'pdist', 'whole_neutral')
random_path = os.path.join(output_path,'random')
for paths in [output_path, pathway_path, visualization_path,
              # stat_tnbc_path,
              # stat_erbc_path,
              stat_whole_path,
              whole_ppr_path,
              # erbc_ppr_path,tnbc_ppr_path,
              stat_chart_whole_path,
              # stat_chart_erbc_path,
              # stat_chart_tnbc_path,
              whole_st_path,
              # erbc_st_path,
              # tnbc_st_path,
              whole_target_path,
              # erbc_target_path,tnbc_target_path,
              whole_dppr_path,
              # erbc_dppr_path,tnbc_dppr_path,
              whole_pdist_path,random_path,pathway_txt_path
              ]:
    os.makedirs(paths, exist_ok=True)

# define the constant and dataset we will use in the project
target_df = pd.read_csv(os.path.join(input_path, 'Drug targets.csv'), encoding="unicode_escape")
target_raw = list(target_df['Targets'])
target_ls = []
for tg in target_raw:
    if type(tg) == str:
        temp = tg.split('; ')
        target_ls.extend(temp)
target_ls = list(set(target_ls))
cancer_df = pd.read_csv(os.path.join(input_path, 'oncogene.csv'), encoding='unicode_escape')
cancer_ls = list(cancer_df['Gene Symbol'])  # list of cancer genes
# ERBC_disease_nodes = ["ESR1", "PGR", "ERBB2", "TP53", "MUC1", "CEACAM5", "BRCA1", "BRCA2"]  # list of disease nodes
# TNBC_disease_nodes = ["EGFR", "KIT", "KRT5", "TP53", "TOP2A", "PARP1",
#                       "HSP90AA1", "HSP90AB1", "MTOR"]  # list of disease nodes
# whole_disease_nodes = list(set(ERBC_disease_nodes + TNBC_disease_nodes))  # list of disease nodes
protease_filename = ['Homo sapiens_GPL96_Lung Carcinoma_cell line_gemcitabine resistant derivative',
                     'Homo sapiens_GPL96_Lung Carcinoma_gender_female[1]',
                     'Homo sapiens_GPL96_Lung Carcinoma_gender_male[2]',
                     'Homo sapiens_GPL570_Lung Carcinoma_cell line_Beas2B',
                     'Homo sapiens_GPL570_Lung Carcinoma_time',
                     'Homo sapiens_GPL570_Prostate Carcinoma_cell line',
                     'Homo sapiens_GPL2986_Prostate Carcinoma_cell line',
                     'Homo sapiens_GPL6244_Lung Carcinoma_cell line',
                     'Homo sapiens_GPL6480_Lung Carcinoma_time',
                     'Homo sapiens_Prostate Carcinoma_time']
