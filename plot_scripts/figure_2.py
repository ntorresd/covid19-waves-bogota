# -*- coding: utf-8 -*-
"""
Created on Mon 29 2022
@author: dsquevedo
@author: ntorresd
"""  
# %%
import sys
import yaml
import matplotlib.pyplot as plt

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plot_scripts')
sys.path.append(SCRIPTS_PATH)
# %%
import results_severe_outcomes as results_severe_outcomes
# %%
fig, ax = plt.subplots(2,3, figsize=(15,10))
fig.set_tight_layout(False)
plt.rcParams["savefig.pad_inches"] = 0.4

axi = ax[0]
results_severe_outcomes.plot_percentage(axi)

axi[0].set_title('hospitalization age percentage distribution')
axi[1].set_title('icu age percentage distribution')
axi[2].set_title('death age percentage distribution')

handles, labels = axi[1].get_legend_handles_labels()
axi[2].legend(handles, labels, loc='upper left')

axi = ax[1]
results_severe_outcomes.plot_counts_hist(axi)

axi[0].set_title('hospitalization counts by age')
axi[1].set_title('icu counts by age')
axi[2].set_title('death counts by age')

# fig.subplots_adjust(bottom=0.1, hspace=0.37)
handles, labels = axi[2].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=len(labels))
fig.savefig(FIG_PATH+'figure_2.png')

# results_severe_outcomes.plot_counts(ax[1])
# results_severe_outcomes.plot_counts_hist(ax[2])
# results_severe_outcomes.plot_proportions_hist(ax[3])
# %%
