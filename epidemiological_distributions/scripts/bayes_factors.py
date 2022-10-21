# -*- coding: utf-8 -*-
"""
Created on Thr May 12 2022
Adapted from: https://github.com/mrc-ide/Brazil_COVID19_distributions

@author: dsquevedo
@author: ntorres
"""
import pandas as pd
import numpy as np
from scipy import stats, special, integrate
import matplotlib.pyplot as plt
import time
import subprocess

t0 = time.time()

strat='wave'

##############################################################################
# Load, process, visualise for the separate models
DATA_PATH = '../../process_data/data/'
OUT_PATH = '../results/'
SAMPLES_PATH = '../fitting_outputs/'
MAX_VAL = 133 
MIN_VAL = 1  # or 1

##############################################################################

# load the data
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
# print n samples and range of data
for df in all_dfs:
    col = str(df.columns[4])
    print(col, len(df[col].index), df[col].min(), '-', df[col].max())
    

##############################################################################
columns_order = [0,1]
columns_ordered = columns#[columns[i] for i in columns_order]
number_samples = pd.DataFrame(columns = columns_ordered, index = strat_)

for df in all_dfs:
    count = df[strat+'_id'].value_counts().to_dict()
    col = df.columns[4]
    number_samples[col] = number_samples.index.map(count)
    
number_samples.fillna(0, inplace=True)
number_samples.to_csv(OUT_PATH + 'number_of_samples.csv')

number_samples_trans = number_samples.transpose()
number_samples_trans.to_csv(OUT_PATH + 'number_of_samples_wide.csv')

##############################################################################
# load the samples (models fits)

dist_posteriors_gamma = {}
dist_posteriors_wei = {}
dist_posteriors_lognorm = {}
dist_posteriors_exp = {}
dist_posteriors_gamma3p = {}
dist_posteriors_gln = {}

for i in range(len(columns)):
    col = columns[i]
    dist_posteriors_gamma.update({col: pd.read_csv(SAMPLES_PATH + col +'-samples-gamma_'+ strat + '.csv')})
    dist_posteriors_lognorm.update({col: pd.read_csv(SAMPLES_PATH + col +'-samples-logn_'+ strat + '.csv')})
    dist_posteriors_wei.update({col: pd.read_csv(SAMPLES_PATH + col +'-samples-weibull_'+ strat + '.csv')})
    dist_posteriors_exp.update({col: pd.read_csv(SAMPLES_PATH + col +'-samples-exponential_'+ strat + '.csv')})  
    dist_posteriors_gamma3p.update({col: pd.read_csv(SAMPLES_PATH + col +'-samples-gamma3p_'+ strat + '.csv')})
    try:
        dist_posteriors_gln.update({col: pd.read_csv(SAMPLES_PATH + col + '-samples-gln_' + strat + '.csv')})
    except:
        print(col, 'did not find samples for gln')
        dist_posteriors_gln.update({col: state_posteriors_gln[columns[0]] * np.nan})

##############################################################################
def gln_pdf(x, mu, sigma, g):
    """PDF of the generalised log-normal distribution"""
    k = g / (2**((g+1)/g) * sigma * special.gamma(1/g))
    return k/x * np.exp(-0.5 * np.abs((np.log(x)-mu)/sigma)**g)

def gln_lpdf(x, mu, sigma, g):
    """Log-PDF of the generalised log-normal distribution"""
    logk = np.log(g) - ( ((g+1)/g)*np.log(2) + np.log(sigma) + np.log(special.gamma(1/g)))
    return logk - np.log(x) - 0.5 * np.abs((np.log(x)-mu)/sigma)**g

def gln_cdf_help(x, mu, sigma, g):
    """CDF of the generalised log-normal distribution"""
    m = x
    result = integrate.quad(lambda x: gln_pdf(x,2,0.5,2.5), 0, m)
    tmp = result[0]
    return tmp
    
def gln_cdf(x, mu, sigma, g):
    c = np.vectorize(gln_cdf_help)
    return c(x, mu, sigma, g)
##############################################################################
# Laplace approx & Bayes Factors to find the best fit

def fit_posterior(distr, col):
    if distr == 'gamma':
        df = dist_posteriors_gamma[col].iloc[:,:-2]
        samples = df.values
        assert samples.shape[1] == len(strat_) * 2 + 2
    elif distr == 'lognormal':
        df = dist_posteriors_lognorm[col].iloc[:,:-2]
        samples = df.values
        assert samples.shape[1] == len(strat_) * 2 + 2
    elif distr == 'weibull':
        df = dist_posteriors_wei[col].iloc[:,:-2]
        samples = df.values
        assert samples.shape[1] == len(strat_) * 2 + 2
    elif distr == 'exponential':
        df = dist_posteriors_exp[col].iloc[:,:-1]
        samples = df.values
        assert samples.shape[1] == len(strat_) * 1 + 1
    elif distr == 'gamma3p':
        df = dist_posteriors_gamma3p[col].iloc[:,:-3]
        samples = df.values
        assert samples.shape[1] == len(strat_) * 3 + 3
    elif distr == 'gln':
        df = dist_posteriors_gln[col].iloc[:,:-3]
        samples = df.values
        assert samples.shape[1] == len(strat_) * 3 + 3
    else:
        print('Invalid distribution')
        return {}
    cov = np.cov(np.array(samples), rowvar = False)
    theta_mean = samples.mean(axis = 0)
    return {'mu': theta_mean, 'cov': cov, 'distr': distr}

columns_BF = columns

posterior_gamma = {}
posterior_wei = {}
posterior_lognorm = {}
posterior_exp = {}
posterior_gamma3p = {}
posterior_gln = {}

for col in columns_BF:
    posterior_gamma.update({col: fit_posterior('gamma', col)})
    posterior_lognorm.update({col: fit_posterior('lognormal', col)})
    posterior_wei.update({col: fit_posterior('weibull', col)})
    posterior_exp.update({col: fit_posterior('exponential', col)})
    posterior_gamma3p.update({col: fit_posterior('gamma3p', col)})
    posterior_gln.update({col: fit_posterior('gln', col)})
    

def Logf(posterior, col):
    col_id = columns.index(col)
    data = all_dfs[col_id]
    distr = posterior[col]['distr']
    theta = posterior[col]['mu']
    x = data[col].values + 0.5 
    strat_vals = data[strat+'_id'].values # int index for the state
    k = len(data[strat+'_id'].unique())
    
    log_prior = 0
    log_likelihood = 0
    
    if distr == 'gamma':
        sigma_alpha = theta[-2]
        sigma_beta = theta[-1]
        alpha_dist = dist_posteriors_gamma[col].alpha.values.mean()
        beta_dist =  dist_posteriors_gamma[col].beta.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_alpha)
        log_prior += stats.norm(0,1).logpdf(sigma_beta)
        for i in range(k):
            log_prior += stats.norm(alpha_dist,sigma_alpha).logpdf(theta[i]) # alphas
            log_prior += stats.norm(beta_dist,sigma_beta).logpdf(theta[i+k]) # betas
        for i in range(len(x)):
            alpha = theta[strat_vals[i]-1]
            beta = theta[strat_vals[i]-1 + k]
            log_likelihood += stats.gamma(a=alpha, scale = 1.0/beta).logpdf(x[i])

    elif distr == 'lognormal':
        sigma_mu = theta[-2]
        sigma_sigma = theta[-1]
        mu_dist = dist_posteriors_lognorm[col].mu.values.mean()
        sigma_dist = dist_posteriors_lognorm[col].sigma.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_mu)
        log_prior += stats.norm(0,1).logpdf(sigma_sigma)
        for i in range(k):
            log_prior += stats.norm(mu_dist,sigma_mu).logpdf(theta[i]) # mu
            log_prior += stats.norm(sigma_dist,sigma_sigma).logpdf(theta[i+k]) # sigmas
        for i in range(len(x)):
            mu = theta[strat_vals[i]-1]
            sigma = theta[strat_vals[i]-1 + k]
            log_likelihood += stats.lognorm(s=sigma, scale = np.exp(mu)).logpdf(x[i])  
            
    elif distr == 'weibull':
        sigma_alpha = theta[-2]
        sigma_sigma = theta[-1]
        alpha_dist = dist_posteriors_wei[col].alpha.values.mean()
        sigma_dist = dist_posteriors_wei[col].sigma.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_alpha)
        log_prior += stats.norm(0,1).logpdf(sigma_sigma)
        for i in range(k):
            log_prior += stats.norm(alpha_dist,sigma_alpha).logpdf(theta[i]) # alphas
            log_prior += stats.norm(sigma_dist,sigma_sigma).logpdf(theta[i+k]) # sigmas
        for i in range(len(x)):
            alpha = theta[strat_vals[i]-1]
            sigma = theta[strat_vals[i]-1 + k]
            log_likelihood += stats.weibull_min(c=alpha, scale = sigma).logpdf(x[i])
    
    elif distr == 'exponential':
        sigma_beta = theta[-2]
        beta_dist = dist_posteriors_exp[col].beta.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_beta)
        for i in range(k):
            log_prior += stats.norm(beta_dist,sigma_beta).logpdf(theta[i]) # betas
        for i in range(len(x)):
            beta = theta[strat_vals[i]-1]
            log_likelihood += stats.expon(scale = 1/beta).logpdf(x[i])
            
    elif distr == 'gamma3p':
        sigma_a = theta[-3]
        sigma_b = theta[-2]
        sigma_c = theta[-1]
        a_dist = dist_posteriors_gamma3p[col].a.values.mean()
        b_dist = dist_posteriors_gamma3p[col].b.values.mean()
        c_dist = dist_posteriors_gamma3p[col].c.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_a)
        log_prior += stats.norm(0,1).logpdf(sigma_b)
        log_prior += stats.norm(0,1).logpdf(sigma_c)
        for i in range(k):
            log_prior += stats.norm(a_dist,sigma_a).logpdf(theta[i]) # a
            log_prior += stats.norm(b_dist,sigma_b).logpdf(theta[i+k]) # b
            log_prior += stats.norm(c_dist,sigma_c).logpdf(theta[i+2*k]) # c
        for i in range(len(x)):
            a = theta[strat_vals[i]-1]
            b = theta[strat_vals[i]-1 + k]
            c = theta[strat_vals[i]-1 + 2*k]
            log_likelihood += stats.gamma(a=a, scale = b, loc=c).logpdf(x[i])
    elif distr == 'gln':
        sigma_mu = theta[-3]
        sigma_sigma = theta[-2]
        sigma_g = theta[-1]
        mu_dist = dist_posteriors_gln[col].mu.values.mean()
        sigma_dist = dist_posteriors_gln[col].sigma.values.mean()
        g_dist = dist_posteriors_gln[col].g.values.mean()
        log_prior += stats.norm(2,0.5).logpdf(sigma_mu)
        log_prior += stats.norm(0.5,0.5).logpdf(sigma_sigma)
        log_prior += stats.norm(1.5,0.5).logpdf(sigma_g)

        for i in range(k):
            log_prior += stats.norm(mu_dist,sigma_mu).logpdf(theta[i]) # mu
            log_prior += stats.norm(sigma_dist,sigma_sigma).logpdf(theta[i+k]) # sigmas
            log_prior += stats.norm(g_dist,sigma_g).logpdf(theta[i+2*k]) # g
        for i in range(len(x)):
            mu = theta[strat_vals[i]-1]
            sigma = theta[strat_vals[i]-1 + k]
            g = theta[strat_vals[i]-1 + 2*k]
            log_likelihood += gln_lpdf(x[i], mu, sigma, g)
    else:
        print('Invalid distribution')
        return np.nan

    return log_likelihood + log_prior    

def LogLaplaceCovariance(posterior, col):
    result = 1/2 * len(posterior[col]['mu']) * np.log(2*np.pi) 
    result += 1/2 * np.log(np.linalg.det(posterior[col]['cov']))
    result += posterior[col]['Logf']
    return result


# careful, this is slow
for col in columns_BF:
    print(col)
    posterior_gamma[col].update({'Logf': Logf(posterior_gamma, col)})
    posterior_gamma[col].update({'LogLaplace': LogLaplaceCovariance(posterior_gamma, col)})
    
    posterior_lognorm[col].update({'Logf': Logf(posterior_lognorm, col)})
    posterior_lognorm[col].update({'LogLaplace': LogLaplaceCovariance(posterior_lognorm, col)})
    
    
    posterior_wei[col].update({'Logf': Logf(posterior_wei, col)})
    posterior_wei[col].update({'LogLaplace': LogLaplaceCovariance(posterior_wei, col)})
    
    posterior_exp[col].update({'Logf': Logf(posterior_exp, col)})
    posterior_exp[col].update({'LogLaplace': LogLaplaceCovariance(posterior_exp, col)})
    
    posterior_gamma3p[col].update({'Logf': Logf(posterior_gamma3p, col)})
    posterior_gamma3p[col].update({'LogLaplace': LogLaplaceCovariance(posterior_gamma3p, col)})
    
    posterior_gln[col].update({'Logf': Logf(posterior_gln, col)})
    posterior_gln[col].update({'LogLaplace': LogLaplaceCovariance(posterior_gln, col)})
    
    
columns_models = ['Gamma', 'Lognormal', 'Weibull', 'Exponential', 'Gamma 3p', 'Gen Lognormal']
bayesfactorsBR = pd.DataFrame(columns = columns_models, index = columns_models)
BayesFactorsdist = {}
evidence_dict = {'Gamma': posterior_gamma, 'Lognormal': posterior_lognorm,
                 'Exponential': posterior_exp, 'Weibull': posterior_wei,
                 'Gamma 3p': posterior_gamma3p, 'Gen Lognormal': posterior_gln
                } 
                  
all_BF = {}
#row - col (negative: in favor of col; positive: in favor of row)
for col in columns_BF:
    for c1 in columns_models:
        for c2 in columns_models:
            bayesfactorsBR.loc[c1,c2] = 2*(evidence_dict[c1][col]['LogLaplace'] - evidence_dict[c2][col]['LogLaplace'])
    all_BF.update({col: bayesfactorsBR})
    print(col)
    print(bayesfactorsBR, '\n')
    df_bf=pd.DataFrame(bayesfactorsBR)
    df_bf.to_csv(f'{OUT_PATH}/bf_{col}.csv')
    

#Selecting best models based on BF
subprocess.call('python best_model.py', shell = True)

print('Time elapsed: ', round((time.time()-t0)/60,1), ' minutes')