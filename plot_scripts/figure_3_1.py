# -*- coding: utf-8 -*-
"""
Created on Wed 1 Dic 2022
@author: davidsantiagoquevedo
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

import results_genomics as rg

fig, ax = plt.subplots(1,2, figsize = (20, 10))
axi = ax[0]
rg.plot_multinomial(axi)
axi.set_xlabel('Week')
axi.set_ylabel('Prevalence')
axi.set_title('a.')

axi = ax[1]
rg.plot_heatmap(axi, n = 200)
axi.set_title('b.')
fig.savefig(FIG_PATH + 'figure_3_1.png')