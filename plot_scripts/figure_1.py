# -*- coding: utf-8 -*-
"""
Created on Mon 28 2022
@author: dsquevedo
@author: ntorresd
"""  

import sys
import yaml
import matplotlib.pyplot as plt

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plot_scripts')
sys.path.append(SCRIPTS_PATH)

import results_rt as results_rt
import results_genomics as results_genomics
import population_pyramid as population_pyramid

fig, ax = plt.subplots(2,2, figsize=(15,10))

results_rt.plot_rt(ax=ax[0][0])
ax[0][0].tick_params(axis='x', rotation=45)

population_pyramid.plot_pyramid(ax=ax[0][1])
ax[0][1].set_xlim(left=-260000)

results_genomics.plot_prevalence(ax=ax[1][1])
ax[1][1].tick_params(axis='x', rotation=45)

fig.savefig(FIG_PATH+'figure_1.png')
fig.show()
