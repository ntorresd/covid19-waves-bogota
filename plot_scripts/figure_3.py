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

import results_epidemiological_distributions as red

fig, ax = plt.subplots(2,5, figsize = (50, 10))
dist = 'hosp_stay'
red.plot_best_model_bar(dist, ax[0][0], 0, 'Hospital stay (days)', 'a.')
red.plot_violin(dist, 'Hospital stay (days)', 'f.', ax[1][0])
dist = 'icu_stay'
red.plot_best_model_bar(dist, ax[0][1], 1, 'ICU stay (days)', 'b.')
red.plot_violin(dist, 'ICU stay (days)', 'g.', ax[1][1])
dist = 'onset_hosp'
red.plot_best_model_bar(dist, ax[0][2], 2, 'Onset to hospitalization (days)', 'c.')
red.plot_violin(dist, 'Onset to hospitalization (days)', 'h.', ax[1][2])
dist = 'onset_icu'
red.plot_best_model_bar(dist, ax[0][3], 3, 'Onset to ICU entrance (days)', 'd.')
red.plot_violin(dist, 'Onset to ICU entrance (days)', 'i.', ax[1][3])
dist = 'onset_death'
red.plot_best_model_bar(dist, ax[0][4], 4, 'Onset to death (days)', 'e.')
red.plot_violin(dist, 'Onset to death (days)', 'j.', ax[1][4])
fig.savefig(FIG_PATH + 'figure_3.png')





fig, axes = plt.subplot_mosaic([['up', 'up'],['left1', 'right1'],
                                ['left2', 'right2'],['left3', 'right3']],
                              figsize = (10,10))
wt=0.1
ax = axes['up']
dist = 'onset_hosp'
w=-0.2
red.plot_best_model_bar_all(dist, ax, w, 3, 'a.')

dist = 'onset_icu'
w=-0.1
red.plot_best_model_bar_all(dist, ax, w, 2, 'a.')

dist = 'onset_death'
w= 0
red.plot_best_model_bar_all(dist, ax, w, 4, 'a.')

dist = 'hosp_stay'
w=0.1
red.plot_best_model_bar_all(dist, ax, w, 1, 'a.')
dist = 'icu_stay'
w=0.2
red.plot_best_model_bar_all(dist, ax, w, 0, 'a.')

ax = axes['left1']
dist = 'onset_hosp'
red.plot_violin(dist, 'Days', 'b.', ax)

ax = axes['right1']
dist = 'onset_icu'
red.plot_violin(dist, 'Days', 'c.', ax)

ax = axes['left2']
dist = 'onset_death'
red.plot_violin(dist, 'Days', 'd.', ax)

ax = axes['right2']
dist = 'hosp_stay'
red.plot_violin(dist, 'Days', 'e.', ax)

ax = axes['left3']
dist = 'icu_stay'
red.plot_violin(dist, 'Days', 'f.', ax)

ax = axes['right3']
ax.axis('off')

ax = axes['up']
ax.legend(['Onset to hospitalization', 
           'Onset to ICU entrance',
           'Onset to death',
           'Hospital stay',
           'ICU stay'], 
          bbox_to_anchor=(0.85, 0.25), bbox_transform=fig.transFigure, ncol=1)

fig.tight_layout()
fig.show()
fig.savefig(FIG_PATH + 'figure_3_v2.png')