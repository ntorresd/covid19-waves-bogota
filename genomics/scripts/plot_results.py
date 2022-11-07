# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 2022

@author: dsquevedo
@author: cwhittaker
@author: ntorres
"""     
import json
import pandas as pd
import matplotlib.pyplot as plt

config_path = open('../../config.json')
config = json.load(config_path)
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH']

plt.style.use(config['PLOTS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

df_results = pd.read_csv(OUT_PATH + 'theta.csv')

stat = "mean"

def plot1():
    fig, ax = plt.subplots()
    
    n = 0 
    variant = 'Alpha'
    mask = (df_results.stat == stat) & (df_results.variant == variant)
    mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
    mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
    ax.plot(df_results[mask].week, df_results[mask].theta, color  = colors[n], label = variant)
    ax.plot(df_results[mask].week, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
            color  = colors[n], marker = '*', linestyle = '')
    ax.fill_between(df_results[mask].week, df_results[mask1].theta, df_results[mask2].theta,
                    color  = colors[n], alpha = 0.2)
    
    n += 1
    variant = 'Delta'
    mask = (df_results.stat == stat) & (df_results.variant == variant)
    mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
    mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
    ax.plot(df_results[mask].week, df_results[mask].theta, color  = colors[n], label = variant)
    ax.plot(df_results[mask].week, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
            color  = colors[n], marker = '*', linestyle = '')
    ax.fill_between(df_results[mask].week, df_results[mask1].theta, df_results[mask2].theta,
                    color  = colors[n], alpha = 0.2)
    
    n += 1
    variant = 'Gamma'
    mask = (df_results.stat == stat) & (df_results.variant == variant)
    mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
    mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
    ax.plot(df_results[mask].week, df_results[mask].theta, color  = colors[n], label = variant)
    ax.plot(df_results[mask].week, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
            color  = colors[n], marker = '*', linestyle = '')
    ax.fill_between(df_results[mask].week, df_results[mask1].theta, df_results[mask2].theta,
                    color  = colors[n], alpha = 0.2)
    
    n += 1
    variant = 'Mu'
    mask = (df_results.stat == stat) & (df_results.variant == variant)
    mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
    mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
    ax.plot(df_results[mask].week, df_results[mask].theta, color  = colors[n], label = variant)
    ax.plot(df_results[mask].week, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
            color  = colors[n], marker = '*', linestyle = '')
    ax.fill_between(df_results[mask].week, df_results[mask1].theta, df_results[mask2].theta,
                    color  = colors[n], alpha = 0.2)
    
    n += 1
    variant = 'Omicron'
    mask = (df_results.stat == stat) & (df_results.variant == variant)
    mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
    mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
    ax.plot(df_results[mask].week, df_results[mask].theta, color  = colors[n], label = variant)
    ax.plot(df_results[mask].week, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
            color  = colors[n], marker = '*', linestyle = '')
    ax.fill_between(df_results[mask].week, df_results[mask1].theta, df_results[mask2].theta,
                    color  = colors[n], alpha = 0.2)
    ax.set_xlabel('Week')
    ax.set_ylabel('Prevalence per variant')
    fig.legend(loc = 'center right')
    
    fig.show()
    fig.savefig(FIG_PATH + 'Variants_multinomial.png', dpi = 100)
    return fig, ax

plot1()