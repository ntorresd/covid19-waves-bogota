# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 2022

@author: dsquevedo
@author: cwhittaker
@author: ntorres
"""
import yaml
import pandas as pd

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'genomics')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'genomics')

# Load raw fitting data
df_fit_raw = pd.read_csv(OUT_PATH + 'fit_raw.csv')
df_variants = pd.read_csv(DATA_PATH + 'variants_pivot.csv')

# 2.5th Percentile
def q025(x):
    return x.quantile(0.025)
# 97.5th Percentile
def q975(x):
    return x.quantile(0.975)

# Prepare results
#theta
variants_dict = {'1' : 'Alpha',
                 '2' : 'Delta',
                 '3' : 'Gamma',
                 '4' : 'Mu',
                 '5' : 'Omicron'
            }
theta_cols = [col for col in df_fit_raw.columns if 'theta' in col]
df_theta = df_fit_raw[theta_cols]
df_theta_mean =  df_theta.agg(['mean', q025, q975])
df_theta_mean.reset_index(inplace = True)
df_theta_mean.rename(columns = {'index' : 'stat'}, inplace = True)
df_theta_mean = pd.melt(df_theta_mean, id_vars=['stat'])
df_theta_mean = df_theta_mean.rename(columns = {'value' : 'theta'})
df_theta_mean[['trash1','week_var']] = df_theta_mean['variable'].str.split('[', 1, expand = True)
df_theta_mean[['week','variant']] = df_theta_mean['week_var'].str.split(',', 1, expand = True)
df_theta_mean[['variant', 'thash2']] = df_theta_mean['variant'].str.split(']', 1, expand = True)
df_theta_mean = df_theta_mean[['week','variant','stat','theta']]
df_theta_mean['week'] = df_theta_mean['week'].astype(int) - 1
df_theta_mean.replace({"variant": variants_dict}, inplace = True)

# Variant counts - week
df_variants['weekly_count_variants'] = df_variants['Alpha'] + df_variants['Delta'] + df_variants['Gamma'] + df_variants['Mu'] + df_variants['Omicron']
df_variants['t'] = df_variants['t'].astype(int)

df_theta_mean = df_theta_mean.merge(df_variants, how = 'left', right_on = 't', left_on = 'week')
del(df_theta_mean['t'])
df_theta_mean.sort_values(by = ['variant', 'week'], inplace = True)
df_theta_mean.to_csv(OUT_PATH + 'theta.csv', index = False)