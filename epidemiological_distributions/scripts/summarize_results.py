# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 2022

@author: dsquevedo
@author: ntorres
"""
import yaml
import pandas as pd
import numpy as np

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'epidemiological_distributions')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'epidemiological_distributions')

import sys
sys.path.append(UTILS_PATH)
import utilities as ut

df_best_models = pd.read_csv(OUT_PATH + 'best_models.csv')
df_best_models = df_best_models.set_index(df_best_models.columns[0]).transpose()

# 1. Prepare the data
all_dfs, columns = ut.prepare_confirmed_cases_data()

stat = ['mean', ut.q975, ut.q025]
##############################################################################
# print n samples and range of data
for df in all_dfs:
    col = str(df.columns[4])
    print(col, len(df[col].index), df[col].min(), '-', df[col].max())
##############################################################################
# load the samples (models fits for every epidemiological distribution)

dist_posteriors  = {'icu_stay':{},
                    'hosp_stay':{},
                    'onset_icu':{},
                    'onset_hosp':{},
                    'onset_death':{}
                   }

for col in columns:
    dist_posteriors[col].update({'Gamma': pd.read_csv(OUT_PATH + col +'-samples-gamma.csv').agg(stat)})
    dist_posteriors[col].update({'Lognormal': pd.read_csv(OUT_PATH + col +'-samples-log_normal.csv').agg(stat)})
    dist_posteriors[col].update({'Weibull': pd.read_csv(OUT_PATH + col +'-samples-weibull.csv').agg(stat)})
    dist_posteriors[col].update({'Exponential': pd.read_csv(OUT_PATH + col +'-samples-exponential.csv').agg(stat)})
    dist_posteriors[col].update({'Gen Lognormal': pd.read_csv(OUT_PATH + col +'-samples-gln.csv').agg(stat)})
    
df_res = pd.DataFrame({})
for dist in dist_posteriors:
    best = df_best_models[df_best_models[dist] == 0].index[0]
    df_temp = dist_posteriors[dist][best]
    df_temp = df_temp.reset_index()
    df_temp = df_temp.rename(columns = {'index':'stat'})
    cols = ['t', 'best', 'stat'] + [col for col in df_temp.columns.tolist() if 'mu' in col]
    df_temp['t'] = dist
    df_temp['best'] = best
    df_temp = df_temp[cols]
    df_res = pd.concat([df_res, df_temp])
    
for col in df_res.columns:
    if 'mu' in col:
        df_res[col] = np.exp(df_res[col]).round(1)
        
df_res.to_csv(OUT_PATH + 'best_fit_summary.csv', index = False)