# -*- coding: utf-8 -*-
"""
Created on Th December 1 2022
@author: dsquevedo
@author: ntorresd
"""  

import warnings
warnings.filterwarnings('ignore')

import sys
import yaml
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plot_scripts')
sys.path.append(SCRIPTS_PATH)

import results_severe_outcomes as results_severe_outcomes

fig = plt.figure(figsize=(15, 15))

gs = GridSpec(7, 3, figure=fig)

ax11 = fig.add_subplot(gs[0:3, 0])
ax12 = fig.add_subplot(gs[0:3, 1])
ax13 = fig.add_subplot(gs[0:3, 2])

ax_histograms = [ax11, ax12, ax13]
results_severe_outcomes.plot_proportions_hist(ax_histograms)

ax21 = fig.add_subplot(gs[3:5, 0])
ax22 = fig.add_subplot(gs[3:5, 1])
ax23 = fig.add_subplot(gs[3:5, 2])
ax31 = fig.add_subplot(gs[5:7, 0])
ax32 = fig.add_subplot(gs[5:7, 1])
axes_rates = [ax21, ax22, ax23, ax31, ax32]
ax33 = fig.add_subplot(gs[5:7, 2])

results_severe_outcomes.plot_rates(ax=ax21, var='CFR', var_name='CFR')
results_severe_outcomes.plot_rates(ax=ax22, var='HCR', var_name='HCR')
results_severe_outcomes.plot_rates(ax=ax23, var='HCR_I', var_name='HCR - ICU')

results_severe_outcomes.plot_rates(ax=ax31, var='HFR', var_name='HFR')
results_severe_outcomes.plot_rates(ax=ax32, var='HFR_I', var_name='HFR - ICU')

ax11.set_ylabel('hospitalization proportion')
ax12.set_ylabel('icu proportion')
ax13.set_ylabel('death proportion')
for axi in ax_histograms:
    axi.set_xlabel('wave')
handles, labels = ax13.get_legend_handles_labels()
ax13.legend(handles, labels, loc='upper right')

handles, labels = ax32.get_legend_handles_labels()
ax33.legend(handles, labels, loc='upper left')
ax33.axis('off')
for axi in axes_rates:
    axi.tick_params(axis='x', rotation=45)

fig.savefig(FIG_PATH+'figure_4.png')
# fig.show()