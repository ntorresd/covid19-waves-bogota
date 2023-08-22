# -*- coding: utf-8 -*-
"""
Created on Mon 29 Nov 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""  

import warnings
warnings.filterwarnings('ignore')

import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt

config = yaml.load(open("config.yml", "r"))["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plots')
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
sys.path.append(SCRIPTS_PATH)

import results_severe_outcomes as results_severe_outcomes

fig, ax = plt.subplots(1,3, figsize=(12,5))
plt.rcParams["savefig.pad_inches"] = 0.4

results_severe_outcomes.plot_percentage(ax)
# results_severe_outcomes.plot_percentage_err(ax)

ax[0].set_ylabel('Porcentaje de ingreso a hospitalización')
ax[1].set_ylabel('Porcentaje de ingreso a UCI')
ax[2].set_ylabel('Porcentaje de fallecidos')
for axi in ax:
    axi.tick_params(axis='x', labelrotation=90)
    axi.set_xlabel('Grupo de edad')
ax[0].set_title('a.')
ax[1].set_title('b.')
ax[2].set_title('c.')
handles, labels = ax[2].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor = (0.79, -0.02), ncol = 6)

fig.savefig(FIG_PATH+'figure_4.png')
# fig.savefig(FIG_PATH+'figure_4_err.png')