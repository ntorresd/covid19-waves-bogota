# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 2022

@author: dsquevedo
@author: ntorresd
""" 

import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plot_scripts')

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# Read confirmed cases information
df_confirmed_bogota = pd.read_csv(DATA_PATH + 'confirmed_cases.csv')
df_confirmed_bogota = df_confirmed_bogota.astype({'age':int})

# functions
def age_group_dec(df,var,var_unit,dic): 
    conditions=[(df[var]<10)|(df[var_unit].isin([2,3]))] + [(df[var]>= dic[ge][0]) & (df[var] <= dic[ge][1]) for ge in list(dic.keys())[1:]]
    choices=[ge for ge in list(dic.keys())]
    df['age_group']=np.select(conditions,choices,default=np.nan)
    return df  

# build age groups
age_group_dic={
    '0-9':[0,9],
    '10-19':[10,19],
    '20-29':[20,29], 
    '30-39':[30,39], 
    '40-49':[40,49], 
    '50-59':[50,59],
    '60-69':[60,69],
    '70-79':[70,79], 
    '80+':[80,9999]
}

df_confirmed_bogota=age_group_dec(df_confirmed_bogota,'age','age_unit', age_group_dic)

def plot_pyramid():
    df_counts = df_confirmed_bogota.groupby(by=['age_group', 'sex']).size().reset_index(name='counts')

    df_female = df_counts[df_counts['sex']=='F']
    df_male = df_counts[df_counts['sex']=='M']
    df_prop = pd.crosstab(index=df_confirmed_bogota['age_group'],
                                columns=df_confirmed_bogota['sex'],
                                normalize="index")*100

    fig, ax1 = plt.subplots(figsize = (15,5))

    hist_f = ax1.barh(df_female['age_group'], -df_female['counts'], align='center', label='Female')
    hist_m = ax1.barh(df_male['age_group'], df_male['counts'], align='center', label='Male')

    ax1.bar_label(hist_m, labels=['{:.1f}%'.format(prop) for prop in df_prop['M'].round(1)], label_type='edge')
    ax1.bar_label(hist_f, labels=['{:.1f}%'.format(prop) for prop in df_prop['F'].round(1)], label_type='edge')

    ticks =  ax1.get_xticks()
    ax1.set_xticklabels([int(abs(tick)) for tick in ticks])

    ax1.set_xlabel('Cases')
    ax1.set_ylabel('Age group')

    fig.legend()
    plt.savefig(FIG_PATH + 'population_pyramid.png')
    return fig, ax1

plot_pyramid()
