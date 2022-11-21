# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 2022

@author: dsquevedo
@author: ntorres
"""     
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'waves')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'waves')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'waves')
# The roots were selected by visual inspection of /figures/roots_confirmed_cases.png 
# and set by default in config.yml
b = config['WAVES']['GAUSSIAN_KERNEL'] 

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

sys.path.append(UTILS_PATH)
import utilities as ut

# Confirmed cases
df_confirmed_bogota = pd.read_csv(DATA_PATH + 'confirmed_cases.csv')
df_confirmed_bogota['onset'] = pd.to_datetime(df_confirmed_bogota['onset'], errors='coerce')
df_waves = pd.read_csv(OUT_PATH + 'waves.csv')
df_waves['start_date'] = pd.to_datetime(df_waves['start_date'], errors='coerce')
df_waves['end_date'] = pd.to_datetime(df_waves['end_date'], errors='coerce')


def plot1(panel = None):
    df_counts = ut.counts('onset',df_confirmed_bogota)
    df_counts['cases_gs'] = ut.gaussian_smoothing(df_counts, 'cases', 'date', b)
    df_counts['cases_gs_diff'] = df_counts.cases_gs.diff()  
    df_counts = df_counts[df_counts.cases_gs_diff.notnull()]
    df_counts['cases_gs_diff_gs'] = ut.gaussian_smoothing(df_counts, 'cases_gs_diff', 'date', b)

    # Plot roots
    fig, ax1 = plt.subplots(figsize = (15,5))

    # Confirmed cases curve
    ln1 = ax1.plot(df_counts['date'], df_counts['cases'], color = colors[0], label = 'Confirmed cases')
    ln2 = ax1.plot(df_counts['date'], df_counts['cases_gs'], color= colors[2], ls = '--')
    ax1.set_ylabel('Confirmed cases')
    date_no_dup = []
    for d in range(len(df_waves)):
            ax1.axvline(x = df_waves['start_date'].iloc[d], color = 'black', alpha = 0.8, ls = '--')
            ax1.axvline(x = df_waves['end_date'].iloc[d], color = 'black', alpha = 0.8, ls = '--')
            date_no_dup.append(d)

    # First difference plot
    ax2 = ax1.twinx()
    ln3 = ax2.plot(df_counts['date'],df_counts['cases_gs_diff'],colors[1],label = 'Diff(Confirmed cases)')
    ln4 = ax2.plot(df_counts['date'],df_counts['cases_gs_diff_gs'],colors[2], ls = '--', label = 'Gaussian smoothing')
    ax2.axhline(y = 0, color = 'black', alpha = 0.8, ls = '--')
    ax2.set_ylabel('Diff(Confirmed cases)')
    ax2.spines.right.set_visible(True) #This was set as False by default in the .mpstyle file

    ##Legend
    lns = ln1 + ln2 + ln3 + ln4
    labs = [l.get_label() for l in lns]
    legend = ax1.legend(lns, labs, loc = 0)

    plt.savefig(FIG_PATH + 'waves.png')
    return fig, ax1

plot1()