# -*- coding: utf-8 -*-
"""
Created on Thr May 12 2022
Adapted from: https://github.com/mrc-ide/Brazil_COVID19_distributions

@author: dsquevedo
@author: ntorres
"""
import pandas as pd
import pystan
import numpy as np
import scipy
import matplotlib.pyplot as plt
import time
from scipy import stats

t0 = time.time()

DATA_PATH = '../../process_data/data/'
OUT_PATH = '../fitting_outputs/'
SEED = 4321
ITER = 2000
CHAINS = 4
MIN_VAL=0
MAX_VAL = 133 

strat='wave'

## Preparing data

drop_columns = ['Start_date', 'End_date']

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

# clean the data and prepare some the variables list 'columns'
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
    try:
        df.drop(columns = drop_columns, inplace = True)
    except:
        print('')
    col = str(df.columns[4])
    columns.append(col)
    df['age_group_id'] = df['age_group'].map(strat_ages_map)
    df['sex_id'] = df['sex'].map(strat_sex_map)
    df['wave_id'] = df['wave'].astype(int)

############################################################################## 
# Posteriors, district level
code_gln_district = """
functions{
 real custom_lpdf(real x, real mu, real sigma, real g)
    {
     real logK = log(g) - (g+1)/g*log(2)-log(sigma)-lgamma(1/g);
     real tmp = logK - log(x) - 0.5 * pow(fabs((log(x)-mu)/sigma),g);
     return tmp;
    }
}
data {
    int N;
    real y[N];
}
parameters {
    real<lower=0> mu;
    real<lower=0> sigma;
    real<lower=1> g;
}
model {
      for (i in 1:N){
        y[i] ~ custom(mu, sigma, g);
        }
      mu ~ normal(2,0.5);
      sigma ~ normal(0.5,0.5);
      g ~ normal(1.5,0.5);
}

"""
model_district = pystan.StanModel(model_code=code_gln_district)

def fit_district(values, list_of_params):
    stdata = values
    stan_data = {'N': len(stdata), 'y': stdata}
    fit = model_district.sampling(data=stan_data, iter=ITER, seed=SEED,
                                  chains=CHAINS, n_jobs=-1)
    print(fit)                            
    df = fit.to_dataframe()
    df = df[list_of_params]
    return df

def get_posteriors_gln(param_list):
    district_posteriors = {}
    for i in range(len(columns)):
        df = all_dfs[i]
        col = columns[i]
        print(col)
        vals = df[col].values
        # watch out here!!! we're shifting the data!!!!
        vals = vals + 0.5
        posterior = fit_district(vals, param_list)
        district_posteriors.update({col: posterior})       
        
    return district_posteriors

district_posteriors_gln = get_posteriors_gln(['mu', 'sigma', 'g'])

############################################################################
# Fitting distribution to the partially pooled data
code_pp_gln = """
functions{
  real custom_lpdf(real x, real mu, real sigma, real g)
    {
      real logK = log(g) - (g+1)/g*log(2)-log(sigma)-lgamma(1/g);
      real tmp = logK - log(x) - 0.5 * pow(fabs((log(x)-mu)/sigma),g);
      return tmp;
    }
}
data {
    int K; //  stratification (age, sex, waves)
    int N; // total number of observations
    real X[N]; // observations
    int INSERT_STRAT[N]; // index with the strat number for each observation
}
parameters {
    real<lower=0> mu[K];
    real<lower=0> sigma[K];
    real<lower=1> g[K];
    // hyperparameters
    real<lower=0> sigma_mu;
    real<lower=0> sigma_sigma;
    real<lower=0> sigma_g;
}

model {
    // likelihood
    for (i in 1:N){
            X[i] ~ custom(mu[INSERT_STRAT[i]], sigma[INSERT_STRAT[i]], g[INSERT_STRAT[i]]);
    }
    // priors
    mu ~ normal(INSERT_MU,sigma_mu);
    sigma ~ normal(INSERT_SIGMA,sigma_sigma);
    g ~ normal(INSERT_G,sigma_g);

    // hyperpriors
    sigma_mu ~ normal(2,0.5);
    sigma_sigma ~ normal(0.5,0.5);
    sigma_g ~ normal(1.5,0.5);

}
"""

def fit_partial_pooling(stan_code, df, col, mu, sigma, g, n_strats, strat_name):
    stan_code = stan_code.replace('INSERT_MU', str(mu))
    stan_code = stan_code.replace('INSERT_SIGMA', str(sigma))
    stan_code = stan_code.replace('INSERT_G', str(g))
    stan_code = stan_code.replace('INSERT_STRAT', str(strat_name))
    model = pystan.StanModel(model_code=stan_code)
    stan_pp_data = {'K': n_strats, 'N': df.shape[0], 
                    'X': df[col].values + 0.5,
                    strat_name : df[strat_name+'_id'].values}
    fit = model.sampling(data=stan_pp_data, iter=ITER, seed=SEED, 
                          chains=CHAINS, n_jobs=-1,
                          control={'adapt_delta': 0.8})
    print(fit)
    posterior_df = fit.to_dataframe()
    params_columns = posterior_df.columns.str.startswith('mu')\
                   + posterior_df.columns.str.startswith('sigma')\
                   + posterior_df.columns.str.startswith('g')\
                   + posterior_df.columns.str.startswith('sigma_')
    posterior_df = posterior_df.loc[:,params_columns]
    return posterior_df

strat_posteriors_gln = {}
for i in range(len(columns)):
    df = all_dfs[i]
    col = columns[i]
    print(col)
    stan_code = code_pp_gln
    mu = district_posteriors_gln[col]['mu'].values.mean()
    sigma = district_posteriors_gln[col]['sigma'].values.mean()
    g = district_posteriors_gln[col]['g'].values.mean()
    posterior = fit_partial_pooling(stan_code, df, col, mu, sigma, g, len(strat_), strat)
    # add national estimates
    posterior = pd.concat([posterior, district_posteriors_gln[col]], axis=1, sort=False)
    strat_posteriors_gln.update({col: posterior})
    # save the output
    posterior.to_csv(OUT_PATH + col +f'-samples-gln_{strat}.csv', index=False)
    strat_posteriors_gln.update({col: posterior})
    del posterior, df

###################################################################
print('GLN model fits done')
print('Time elapsed: ', round((time.time()-t0)/60,1), ' minutes')