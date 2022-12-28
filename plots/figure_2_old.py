# -*- coding: utf-8 -*-
"""
Created on Mon 29 2022
@author: dsquevedo
@author: ntorresd
"""  

import warnings
warnings.filterwarnings('ignore')

import sys
import yaml
import matplotlib.pyplot as plt

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plot_scripts')
sys.path.append(SCRIPTS_PATH)

import results_severe_outcomes as results_severe_outcomes

fig, ax = plt.subplots(2,3, figsize=(15,10))
fig.set_tight_layout(False)
plt.rcParams["savefig.pad_inches"] = 0.4

axi = ax[0]
results_severe_outcomes.plot_percentage(axi)

axi[0].set_ylabel('hospitalization percentage by age')
axi[1].set_ylabel('icu percentage by age')
axi[2].set_ylabel('death percentage by age')
for axii in axi:
    axii.tick_params(axis='x', labelrotation=90)
    axii.set_xlabel('age group')

handles, labels = axi[2].get_legend_handles_labels()
axi[2].legend(handles, labels, loc='upper left')

axi = ax[1]
vary = 'counts'
strat = 'wave'
results_severe_outcomes.plot_counts_histograms(axi)
for axii in axi:
    axii.set_xlabel(strat)
axi[0].set_ylabel('hospitalization counts by age')
axi[1].set_ylabel('icu counts by age')
axi[2].set_ylabel('death counts by age')
handles, labels = axi[2].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=len(labels))

fig.savefig(FIG_PATH+'figure_2.png')
# fig.show()