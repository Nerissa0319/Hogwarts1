import pandas as pd
import os.path
from ConstantNeu import *
import networkx as nx
from sklearn import preprocessing
import numpy as np
import json


def rank(dict):
    sorted_dict = {k: v for k, v in sorted(dict.items(), key=lambda item: item[1], reverse=True)}
    key = list(sorted_dict.keys())
    value = list(sorted_dict.values())
    ranking_ls = []
    current_value = 10000
    current_index = 0
    for index, p in enumerate(value, start=1):
        if p < current_value:
            current_index = index
        ranking_ls.append(current_index)
        current_value = p

    ranking_df = pd.DataFrame()
    for i in range(len(key)):
        data = pd.DataFrame({'Gene': [key[i]], 'Value': [value[i]], 'Rank': [ranking_ls[i]]})
        ranking_df = pd.concat([ranking_df, data])

    return ranking_df



