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

fig, ax = plt.subplots(1,5, figsize = (50, 10))
dist = 'hosp_stay'
red.plot_best_model_bar(dist, ax[0], 0, 'Hospital stay (days)', 'a.')
dist = 'icu_stay'
red.plot_best_model_bar(dist, ax[1], 0, 'ICU stay (days)', 'b.')
dist = 'onset_hosp'
red.plot_best_model_bar(dist, ax[2], 0, 'Onset to hospitalization (days)', 'c.')
dist = 'onset_icu'
red.plot_best_model_bar(dist, ax[3], 0, 'Onset to ICU entrance (days)', 'd.')
dist = 'onset_death'
red.plot_best_model_bar(dist, ax[4], 0, 'Onset to death (days)', 'e.')
fig.savefig(FIG_PATH + 'Figure_3.png')


# Plot best distributions
fig, ax = plt.subplots(figsize = (10, 10))
dist = 'icu_stay'
red.plot_best_model_line(dist, ax, 0)
dist = 'hosp_stay'
red.plot_best_model_line(dist, ax, 1)
dist = 'onset_hosp'
red.plot_best_model_line(dist, ax, 2)
dist = 'onset_icu'
red.plot_best_model_line(dist, ax, 3)
dist = 'onset_death'
red.plot_best_model_line(dist, ax, 4)
handles, labels = ax.get_legend_handles_labels()
# remove the errorbars
handles = [h[0] for h in handles]
# use them in the legend
ax.legend(handles, labels, loc='upper right',numpoints=1)
fig.savefig(FIG_PATH + 'best_models.png')