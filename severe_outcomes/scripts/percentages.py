# -*- coding: utf-8 -*-
"""
Created on Thr Jul 31 2022

@author: dsquevedo
@author: ntorresd
"""

import sys, os
import yaml
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

# Load configuration
ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

# Paths
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'severe_outcomes')

# Plot style
plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

df_hosp = pd.read_csv(DATA_PATH+'hosp_waves_bog.csv')
df_icu = pd.read_csv(DATA_PATH+'icu_waves_bog.csv')
df_death = pd.read_csv(DATA_PATH+'death_waves_bog.csv')

# Percentages calculation
def calculate_percentage(df, strat='wave', group='age_group', print_sum=False):
    df_percentage = df.groupby(by=[strat, group]).size().reset_index(name='counts')
    df_percentage['percentage'] = 0
    strat_list = df_percentage[strat].unique()
    for strat_ in strat_list:
        mask = (df_percentage[strat]==strat_)
        df_percentage.loc[mask, 'percentage'] = 100 * df_percentage.loc[mask, 'counts'] / df_percentage.loc[mask, 'counts'].sum()
        if print_sum:
            print('The sum over the percentages gives: ', df_percentage.loc[mask ,'percentage'].sum())
    return df_percentage

df_hosp_perc = calculate_percentage(df_hosp)
df_icu_perc = calculate_percentage(df_icu)
df_death_perc = calculate_percentage(df_death)  

df_hosp_perc.to_csv(OUT_PATH+'hosp_percentages.csv', index=False)
df_icu_perc.to_csv(OUT_PATH+'icu_percentages.csv', index=False)
df_death_perc.to_csv(OUT_PATH+'deaths_percentages.csv', index=False)