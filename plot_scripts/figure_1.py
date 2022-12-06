# -*- coding: utf-8 -*-
"""
Created on Mon 28 2022
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

import results_rt as results_rt
import results_genomics as results_genomics
# import results_waves as results_waves
import overview as overview

fig, ax = plt.subplots(2,2, figsize=(15,10))

axi = ax[0][0]
results_rt.plot_rt(ax=ax[0][0])
axi.set_xlabel('')
axi.set_ylabel('R(t)') 
axi.set_ylim(top=2.3)
axi.tick_params(axis='x', rotation=45)

axi = ax[0][1]
overview.plot_pyramid(ax=axi)
axi.set_xlim(left=-260000)
axi.set_xlabel('cases')
axi.set_ylabel('age group')
axi.legend()  

axi = ax[1][0]
overview.plot_cases_death_cum(ax=axi)

axi = ax[1][1]
results_genomics.plot_prevalence(ax=axi)
axi.tick_params(axis='x', rotation=45)
axi.set_xlabel('week')
axi.set_ylabel('prevalence')

fig.savefig(FIG_PATH+'figure_1.png')
# fig.show()
