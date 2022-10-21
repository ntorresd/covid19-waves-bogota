# -*- coding: utf-8 -*-
"""
Created on Thr Jul 19 2022
Adapted from: https://github.com/mrc-ide/Brazil_COVID19_distributions

@author: dsquevedo
@author: ntorres
"""
import pandas as pd
import numpy as np
import datetime as dt

import matplotlib.pyplot as plt
plt.style.use('../../plot_scripts/plot_style.mplstyle')

from scipy import stats
import scipy as scipy

##############################################################################
# Load, process, visualise for the separate models
DATA_PATH = '../../process_data/data/'
OUT_PATH = '../figures/'
SAMPLES_PATH = '../fitting_outputs/'
RESULTS_PATH = '../results/'
MAX_VAL = 133 
MIN_VAL = 1  # or 1

strat = 'wave'
stat = 'median'
##############################################################################
# load best models results
df_best_models = pd.read_csv(RESULTS_PATH + 'best_models.csv')
df_best_models = df_best_models.set_index(df_best_models.columns[0]).transpose()

##############################################################################
# load the observed data
drop_columns = ['start_date', 'end_date']

df_icu_stay = pd.read_csv(DATA_PATH + 'icu_stay_bog.csv')
df_icu_stay = df_icu_stay[(df_icu_stay['icu_stay'] > MIN_VAL)&(df_icu_stay['icu_stay'] <= MAX_VAL)]

df_hosp_stay = pd.read_csv(DATA_PATH + 'hosp_stay_bog.csv')
df_hosp_stay = df_hosp_stay[(df_hosp_stay['hosp_stay'] > MIN_VAL)&(df_hosp_stay['hosp_stay'].abs() <= MAX_VAL)]

df_onset_icu = pd.read_csv(DATA_PATH + 'onset_icu_bog.csv')
df_onset_icu = df_onset_icu[(df_onset_icu['onset_icu'] > MIN_VAL)&(df_onset_icu['onset_icu'] <= MAX_VAL)]

df_onset_hosp = pd.read_csv(DATA_PATH + 'onset_hosp_bog.csv')
df_onset_hosp = df_onset_hosp[(df_onset_hosp['onset_hosp'] > MIN_VAL)&(df_onset_hosp['onset_hosp'] <= MAX_VAL)]

df_onset_death = pd.read_csv(DATA_PATH + 'onset_death_bog.csv')
df_onset_death = df_onset_death[(df_onset_death['onset_death'] > MIN_VAL)&(df_onset_death['onset_death'] <= MAX_VAL)]

all_dfs = [df_icu_stay, df_hosp_stay, df_onset_icu, df_onset_hosp, df_onset_death]

# clean and preparation of data
for df in all_dfs:
    df.dropna(inplace=True)
    
strat_ages = df_onset_icu['age_group'].unique()
strat_sex = df_onset_icu['sex'].unique()
strat_wave= df_onset_icu['wave'].unique()

strat_sex.sort()
strat_ages.sort()
strat_wave.sort()

strat_sex_map = dict(zip(strat_sex, list(range(1, len(strat_sex)+1))))
strat_sex = list(range(1, len(strat_sex)+1))

strat_ages_map = dict(zip(strat_ages, list(range(1, len(strat_ages)+1))))
strat_ages = list(range(1, len(strat_ages)+1))

strat_wave_map = dict(zip(strat_wave, list(range(1, len(strat_wave)+1))))
strat_wave = list(range(1, len(strat_wave)+1))

if strat=='wave':
    strat_=strat_wave
elif strat=='age':
    strat_=strat_age
elif strat=='sex':
    strat_=strat_sex

columns = []
for df in all_dfs:
    df.dropna(inplace=True) # remove the rows with nan values
    col = str(df.columns[4])
    columns.append(col)
    df['age_group_id'] = df['age_group'].map(strat_ages_map)
    df['sex_id'] = df['sex'].map(strat_sex_map)
    df['wave_id'] = df['wave'].astype(int) 

##############################################################################
# print n samples and range of data
for df in all_dfs:
    col = str(df.columns[4])
    print(col, len(df[col].index), df[col].min(), '-', df[col].max())
##############################################################################
# load the samples (models fits for every epidemiological distribution)

dist_posteriors  = {'icu_stay':{},
                    'hosp_stay':{},
                    'onset_icu':{},
                    'onset_hosp':{},
                    'onset_death':{}
                   }

for col in columns:
    dist_posteriors[col].update({'Gamma': pd.read_csv(SAMPLES_PATH + col +'-samples-gamma_'+ strat + '.csv').agg(stat)})
    dist_posteriors[col].update({'Lognormal': pd.read_csv(SAMPLES_PATH + col +'-samples-logn_'+ strat + '.csv').agg(stat)})
    dist_posteriors[col].update({'Weibull': pd.read_csv(SAMPLES_PATH + col +'-samples-weibull_'+ strat + '.csv').agg(stat)})
    dist_posteriors[col].update({'Exponential': pd.read_csv(SAMPLES_PATH + col +'-samples-exponential_'+ strat + '.csv').agg(stat)})
    dist_posteriors[col].update({'Gen Lognormal': pd.read_csv(SAMPLES_PATH + col +'-samples-gln_'+ strat + '.csv').agg(stat)})
    dist_posteriors[col].update({'Gamma 3p': pd.read_csv(SAMPLES_PATH + col +'-samples-gamma3p_'+ strat + '.csv').agg(stat)})  
########################################################################
# Usefull statistical functions 
def gln_pdf(x, mu, sigma, g):
    """PDF of the generalised log-normal distribution"""
    k = g / (2**((g+1)/g) * sigma * scipy.special.gamma(1/g))
    return k/x * np.exp(-0.5 * np.abs((np.log(x)-mu)/sigma)**g)

def gln_cdf_help(x, mu, sigma, g):
    """CDF of the generalised log-normal distribution"""
    m = x
    result = scipy.integrate.quad(lambda x: gln_pdf(x,mu,sigma,g), 0, m)
    tmp = result[0]
    return tmp
    
def gln_cdf(x, mu, sigma, g):
    c = np.vectorize(gln_cdf_help)
    return c(x, mu, sigma, g)
########################################################################
# Plot CDFs
def plot_cdf(df, epi_dist, max_val, ax, cdf_null_hyp, n, n_strat = None):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    if n_strat:
        data = df.loc[(df[strat]==n_strat) & (df[epi_dist]<=max_val)][epi_dist]
        ext = f'[{n_strat}]'
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
def plot_dist(df, epi_dist, max_val, ax, n_strat = None,
              dist_list = ['gamma', 'lognormal', 'weibull', 'exponential'],
              bin_unit = 1, title=None):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    best = df_best_models[df_best_models[epi_dist]==0].index[0]
    if n_strat:
        data = df.loc[(df[strat]==n_strat) & (df[epi_dist]<=max_val)][epi_dist]
        if max_val > data.max(): 
            max_val = data.max()
            data = df.loc[(df[strat]==n_strat) & (df[epi_dist]<=max_val)][epi_dist]
        ext = f'[{n_strat}]'
    else:
        data = df.loc[df[epi_dist]<=max_val][epi_dist]   
        if max_val > data.max(): 
            max_val = data.max()
            data = df.loc[(df[strat]==n_strat) & (df[epi_dist]<=max_val)][epi_dist]
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
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_strat = n_strat)
        elif dist == 'Lognormal':
            s = params[dist]['sigma'+ext]
            mu = params[dist]['mu'+ext]
            fit_stan = stats.lognorm.pdf(xx, s = s, scale = np.exp(mu))
            if dist == best:
                cdf_null_hyp = [stats.lognorm.cdf(x, s = s, scale = np.exp(mu)) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_strat = n_strat)
        elif dist == 'Weibull':
            a = params[dist]['alpha'+ext]
            s = params[dist]['sigma'+ext]
            fit_stan = stats.weibull_min.pdf(xx, c = a, scale = s)
            if dist == best:
                cdf_null_hyp = [stats.weibull_min.cdf(x, c = a, scale = s) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n=n, n_strat = n_strat)
        elif dist == 'Exponential':
            b = params[dist]['beta'+ext]
            fit_stan = stats.expon.pdf(xx, loc = 0, scale = 1/b)
            if dist == best:
                cdf_null_hyp = [stats.expon.cdf(x, loc = 0 , scale = 1/b) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_strat = n_strat)
        elif dist == 'Gamma 3p':
            fit_stan=stats.gamma.pdf(xx, a=stan_stats[n_strat-1], scale=stan_stats[3+n_strat], loc=stan_stats[7+n_strat])
        elif dist == 'Gen Lognormal':
            mu = params[dist]['mu'+ext]
            s = params[dist]['sigma'+ext]
            g = params[dist]['g'+ext]
            fit_stan = gln_pdf(xx, mu = mu, sigma = s, g = g)
            if dist == best:
                cdf_null_hyp = [gln_cdf(x, mu = mu, sigma = s, g = g) for x in data_sorted]
                plot_cdf(df, epi_dist = var, max_val = max_val_plot,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_strat = n_strat)
        #Fit
        ax.plot(xx, fit_stan, label = dist, color = colors[n])        
    ax.set_xlim([0,max_val])
    ax.set_xlabel(epi_dist.replace('_',' '))
    ax.set_ylabel('Probability')
    ax.set_title(title)
    return ax

########################################################################
# Plots
dist_list = ['Gamma', 'Lognormal', 'Weibull', 'Exponential', 'Gen Lognormal']
max_val_plot=60

########################################################################
########################################################################
# Hospital stay (general hospital bed and ICU bed)
fig, ax = plt.subplots(2,4,figsize=(15, 7.5))
#ICU stay
var = 'icu_stay'
df = df_icu_stay
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][0], n_strat = 1, dist_list = dist_list,
          bin_unit = 1, title='a.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][1], n_strat = 2, dist_list = dist_list,
          bin_unit = 1, title='b.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][2], n_strat = 3, dist_list = dist_list,
          bin_unit = 1, title='c.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][3], n_strat = 4, dist_list = dist_list,
          bin_unit = 1, title='d.')
#HOSP stay
var = 'hosp_stay'
df = df_hosp_stay
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][0], n_strat = 1, dist_list = dist_list,
          bin_unit = 1, title='e.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][1], n_strat = 2, dist_list = dist_list,
          bin_unit = 1, title='f.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][2], n_strat = 3, dist_list = dist_list,
          bin_unit = 1, title='g.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][3], n_strat = 4, dist_list = dist_list,
          bin_unit = 1, title='h.')

handles, labels = ax[0][0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor = (0.8, -0.03), ncol = len(dist_list))
fig.savefig(OUT_PATH+'icu_hosp_distributions.png')

########################################################################
########################################################################
# Onset to several outcomes
fig, ax = plt.subplots(3,4,figsize=(15, 10))
#ONSET HOSP
var = 'onset_hosp'
df = df_onset_hosp
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][0], n_strat = 1, dist_list = dist_list,
          bin_unit = 1, title='a.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][1], n_strat = 2, dist_list = dist_list,
          bin_unit = 1, title='b.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][2], n_strat = 3, dist_list = dist_list,
          bin_unit = 1, title='c.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][3], n_strat = 4, dist_list = dist_list,
          bin_unit = 1, title='d.')
#ONSET ICU
var = 'onset_icu'
df = df_onset_icu
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][0], n_strat = 1, dist_list = dist_list,
          bin_unit = 1, title='e.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][1], n_strat = 2, dist_list = dist_list,
          bin_unit = 1, title='f.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][2], n_strat = 3, dist_list = dist_list,
          bin_unit = 1, title='g.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][3], n_strat = 4, dist_list = dist_list,
          bin_unit = 1, title='h.')
#ONSET DEATH
var = 'onset_death'
df = df_onset_death
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][0], n_strat = 1, dist_list = dist_list,
          bin_unit = 1, title='i.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][1], n_strat = 2, dist_list = dist_list,
          bin_unit = 1, title='j.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][2], n_strat = 3, dist_list = dist_list,
          bin_unit = 1, title='k.')
plot_dist(df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][3], n_strat = 4, dist_list = dist_list,
          bin_unit = 1, title='l.')

handles, labels = ax[0][0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor = (0.8, -0.03), ncol = len(dist_list))
fig.savefig(OUT_PATH+'onset_distributions.png')