# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 2022
@author: dsquevedo
@author: ntorresd
"""     

import yaml
import pandas as pd
import matplotlib.pyplot as plt

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

# Paths
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'rt')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'rt')

# Plot style
plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# Instantaneous reproduction number 
df_rt_all = pd.read_csv(OUT_PATH + 'rt_all_ages.csv')
df_rt_all['window_start'] = pd.to_datetime(df_rt_all['window_start'])
df_rt_all['window_end'] = pd.to_datetime(df_rt_all['window_end'])

df_rt_60p = pd.read_csv(OUT_PATH + 'rt_60_plus.csv')
df_rt_60p['window_start'] = pd.to_datetime(df_rt_60p['window_start'])
df_rt_60p['window_end'] = pd.to_datetime(df_rt_60p['window_end'])

# Instantaneous reproduction number R(t) plot
def plot_rt(ax):    
    ln1 = ax.plot(df_rt_all['window_end'], df_rt_all['Mean(R)'], color=colors[4], label = 'all')
    ln2 = ax.plot(df_rt_60p['window_end'], df_rt_60p['Mean(R)'], color=colors[1], label = '60+')
    ax.axhline(y=1, color='black', linestyle='--', linewidth=3)

    lns = ln1 + ln2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc = 0)

fig, ax = plt.subplots()
plot_rt(ax)
ax.set_xlabel('')
ax.set_ylabel('R(t)') 
ax.set_ylim(top=2.3)
ax.tick_params(axis='x', rotation=45)
fig.savefig(FIG_PATH + 'rt.png')
# fig.show()