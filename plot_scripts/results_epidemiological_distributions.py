# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 2022

@author: dsquevedo
@author: ntorres
"""
import yaml
import pandas as pd
import numpy as np
from scipy import stats
import scipy as scipy
import matplotlib.pyplot as plt
import seaborn as sns

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'epidemiological_distributions')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'epidemiological_distributions')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'epidemiological_distributions')

import sys
sys.path.append(UTILS_PATH)
import utilities as ut

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
"""
df_best_models = pd.read_csv(OUT_PATH + 'best_models.csv')
df_best_models = df_best_models.set_index(df_best_models.columns[0]).transpose()

##############################################################################
# 1. Prepare the data
all_dfs, columns = ut.prepare_confirmed_cases_data()

##############################################################################
# load the samples (models fits for every epidemiological distribution)
dist_posteriors = ut.load_samples()
"""
########################################################################
# Plot CDFs
def plot_cdf(df, epi_dist, max_val, ax, cdf_null_hyp, n, n_subset = None, subset = 'wave'):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    if n_subset:
        data = df.loc[(df[subset] == n_subset) & (df[epi_dist]<=max_val)][epi_dist]
        ext = f'[{n_subset}]'
    else:
        data = df.loc[df[epi_dist]<=max_val][epi_dist]   
        ext = ''
    data_sorted = np.sort(np.array(data.drop_duplicates()))
    #get cdf from obs data
    count, bins_count = np.histogram(data, bins=len(data_sorted))  
    # finding the PDF of the histogram using count values
    pdf_obs = count / sum(count)
    # using numpy np.cumsum to calculate the CDF
    # We can also find using the PDF values by looping and adding
    cdf_obs = np.cumsum(pdf_obs)    
    axt = ax.twinx()
    axt.plot(data_sorted, cdf_obs, linestyle='--', marker='',color='black')
    axt.plot(data_sorted, cdf_null_hyp,  linestyle='-',color=colors[n])
    """
    #get max different (K-S statistic)
    arr_dif_abs = np.abs(cdf_null_hyp-np.array(cdf_obs))
    dn_ks = max(arr_dif_abs)
    # Plot distances between cdf_null_hyp and observed data
    for x1,x2, y1, y2 in zip(data_sorted,data_sorted, cdf_obs, cdf_null_hyp):
        axt.plot([x1, x2], [y1, y2], color='green', alpha=0.4)
    """
    axt.spines.right.set_visible(True) #It was set as False by default in the .mpstyle file
    
########################################################################
# Plot function (PDFs and CDFs)
def plot_dist(df, epi_dist, max_val, ax, n_subset = None, subset = 'wave',
              dist_list = ['gamma', 'lognormal', 'weibull', 'exponential'],
              bin_unit = 1, title=None):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    best = df_best_models[df_best_models[epi_dist] == 0].index[0]
    if n_subset:
        data = df.loc[(df[subset ]== n_subset) & (df[epi_dist] <= max_val)][epi_dist]
        if max_val > data.max(): 
            max_val = data.max()
            data = df.loc[(df[subset] == n_subset) & (df[epi_dist] <= max_val)][epi_dist]
        ext = f'[{n_subset}]'
    else:
        data = df.loc[df[epi_dist]<=max_val][epi_dist]   
        if max_val > data.max(): 
            max_val = data.max()
            data = df.loc[(df[subset] == n_subset) & (df[epi_dist]<=max_val)][epi_dist]
        ext = ''
    data_sorted = np.sort(np.array(data.drop_duplicates()))
    params = dist_posteriors[epi_dist]      
    bins = int(max_val/bin_unit)
    ax.hist(data, density=True, bins = bins, alpha = 0.4, edgecolor = 'white')#, label=col)
    xx = np.linspace(min(data),max(data),1000)
    for n, dist in enumerate(dist_list):
        if dist == 'Gamma':
            a = params[dist]['alpha'+ext]
            b = params[dist]['beta'+ext]
            fit_stan = stats.gamma.pdf(xx, a = a, scale = 1/b, loc = 0)
            if dist == best:
                cdf_null_hyp = [stats.gamma.cdf(x, a = a, scale = 1/b, loc = 0) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        elif dist == 'Lognormal':
            s = params[dist]['sigma'+ext]
            mu = params[dist]['mu'+ext]
            fit_stan = stats.lognorm.pdf(xx, s = s, scale = np.exp(mu))
            if dist == best:
                cdf_null_hyp = [stats.lognorm.cdf(x, s = s, scale = np.exp(mu)) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        elif dist == 'Weibull':
            a = params[dist]['alpha'+ext]
            s = params[dist]['sigma'+ext]
            fit_stan = stats.weibull_min.pdf(xx, c = a, scale = s)
            if dist == best:
                cdf_null_hyp = [stats.weibull_min.cdf(x, c = a, scale = s) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n=n, n_subset = n_subset)
        elif dist == 'Exponential':
            b = params[dist]['beta'+ext]
            fit_stan = stats.expon.pdf(xx, loc = 0, scale = 1/b)
            if dist == best:
                cdf_null_hyp = [stats.expon.cdf(x, loc = 0 , scale = 1/b) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        elif dist == 'Gen Lognormal':
            mu = params[dist]['mu'+ext]
            s = params[dist]['sigma'+ext]
            g = params[dist]['g'+ext]
            fit_stan = ut.gln_pdf(xx, mu = mu, sigma = s, g = g)
            if dist == best:
                cdf_null_hyp = [ut.gln_cdf(x, mu = mu, sigma = s, g = g) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        #Fit
        ax.plot(xx, fit_stan, label = dist, color = colors[n])        
    ax.set_xlim([0,max_val])
    ax.set_xlabel(epi_dist.replace('_',' '))
    ax.set_ylabel('Probability')
    ax.set_title(title)
    return ax
########################################################################
# Plot epi distributions
"""
dist_list = ['Gamma', 'Lognormal', 'Weibull', 'Exponential', 'Gen Lognormal']
max_val_plot=60

########################################################################
########################################################################
# Hospital stay (general hospital bed and ICU bed)
fig, ax = plt.subplots(2, 4, figsize = (15, 7.5))
#ICU stay
var = 'icu_stay'
df = all_dfs[0] #df_icu_stay
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='a.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='b.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='c.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='d.')

#HOSP stay
var = 'hosp_stay'
df = all_dfs[1] #df_hosp_stay
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='e.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='f.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='g.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='h.')
handles, labels = ax[0][0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor = (0.8, -0.03), ncol = len(dist_list))
fig.savefig(FIG_PATH + 'icu_hosp_distributions.png')

########################################################################
########################################################################
# Onset to several outcomes
fig, ax = plt.subplots(3,4,figsize=(15, 10))
#ONSET HOSP
var = 'onset_hosp'
df = all_dfs[3] #df_onset_hosp
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='a.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='b.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='c.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='d.')
#ONSET ICU
var = 'onset_icu'
df = all_dfs[2] #df_onset_icu
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='e.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='f.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='g.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='h.')
#ONSET DEATH
var = 'onset_death'
df = all_dfs[4] #df_onset_death
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='i.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='j.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='k.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='l.')

handles, labels = ax[0][0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor = (0.8, -0.03), ncol = len(dist_list))
fig.savefig(FIG_PATH + 'onset_distributions.png')
"""

########################################################################
########################################################################
########################################################################
# Plot best model
def plot_best_model_line(dist, ax, n, title):
    df_best_models_sum = pd.read_csv(OUT_PATH + "best_fit_summary.csv")
    df_dist = df_best_models_sum[df_best_models_sum.dist == dist]
    df_dist = df_dist.set_index('stat')
    name = dist.replace('_', ' ')
    mean = np.array([df_dist['wave_1']['mean'], df_dist['wave_2']['mean'], 
                     df_dist['wave_3']['mean'], df_dist['wave_4']['mean']])
    err = abs(np.array([[df_dist['wave_1']['mean'] - df_dist['wave_1']['q025'], 
                     df_dist['wave_2']['mean'] - df_dist['wave_2']['q025'], 
                     df_dist['wave_3']['mean'] - df_dist['wave_3']['q025'], 
                     df_dist['wave_4']['mean'] - df_dist['wave_4']['q025']],
                    [df_dist['wave_1']['q975'] - df_dist['wave_1']['mean'], 
                     df_dist['wave_2']['q975'] - df_dist['wave_2']['mean'], 
                     df_dist['wave_3']['q975'] - df_dist['wave_3']['mean'], 
                     df_dist['wave_4']['q975'] - df_dist['wave_4']['mean']]]))
    ax.errorbar(['1', '2', '3', '4'], mean, err, 
                ls = '-', marker = 'o', 
                color = colors[n],
                label = name)
    ax.set_xlabel('Wave')
    ax.set_ylabel('Days')
    ax.set_title(title)
    return ax

def plot_best_model_bar(dist, ax, n, name_y, title):
    df_best_models_sum = pd.read_csv(OUT_PATH + "best_fit_summary.csv")
    df_dist = df_best_models_sum[df_best_models_sum.dist == dist]
    df_dist = df_dist.set_index('stat')
    name = dist.replace('_', ' ')
    mean = np.array([df_dist['wave_1']['mean'], df_dist['wave_2']['mean'], 
                     df_dist['wave_3']['mean'], df_dist['wave_4']['mean']])
    err = abs(np.array([[df_dist['wave_1']['mean'] - df_dist['wave_1']['q025'], 
                     df_dist['wave_2']['mean'] - df_dist['wave_2']['q025'], 
                     df_dist['wave_3']['mean'] - df_dist['wave_3']['q025'], 
                     df_dist['wave_4']['mean'] - df_dist['wave_4']['q025']],
                    [df_dist['wave_1']['q975'] - df_dist['wave_1']['mean'], 
                     df_dist['wave_2']['q975'] - df_dist['wave_2']['mean'], 
                     df_dist['wave_3']['q975'] - df_dist['wave_3']['mean'], 
                     df_dist['wave_4']['q975'] - df_dist['wave_4']['mean']]]))
    ax.bar(['1', '2', '3', '4'], mean, yerr = err,
    color = colors[n], label = name)
    ax.set_xlabel('Wave')
    ax.set_ylabel(name_y)
    ax.set_title(title)