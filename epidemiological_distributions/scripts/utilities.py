# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 2022

@author: dsquevedo
@author: ntorres
"""
import yaml
import pandas as pd
ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'epidemiological_distributions')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'epidemiological_distributions')
SCRIPTS_PATH = config['PATHS']['SCRIPTS_PATH'].format(dir = 'epidemiological_distributions')

SEED = config['MODELS']['SEED']
ITER = config['MODELS']['ITER']
CHAINS = config['MODELS']['CHAINS']
MIN_VAL = config['MODELS']['MIN_VAL']
MAX_VAL = config['MODELS']['MAX_VAL']

def prepare_confirmed_cases_data(strat = 'wave'):
    drop_columns = ['Start_date', 'End_date']

    df_icu_stay = pd.read_csv(DATA_PATH + 'icu_stay_bog.csv')
    df_icu_stay = df_icu_stay[(df_icu_stay['icu_stay'] > MIN_VAL) & (df_icu_stay['icu_stay'] <= MAX_VAL)]

    df_hosp_stay = pd.read_csv(DATA_PATH + 'hosp_stay_bog.csv')
    df_hosp_stay = df_hosp_stay[(df_hosp_stay['hosp_stay'] > MIN_VAL) & (df_hosp_stay['hosp_stay'].abs() <= MAX_VAL)]

    df_onset_icu = pd.read_csv(DATA_PATH + 'onset_icu_bog.csv')
    df_onset_icu = df_onset_icu[(df_onset_icu['onset_icu'] > MIN_VAL) & (df_onset_icu['onset_icu'] <= MAX_VAL)]

    df_onset_hosp = pd.read_csv(DATA_PATH + 'onset_hosp_bog.csv')
    df_onset_hosp = df_onset_hosp[(df_onset_hosp['onset_hosp'] > MIN_VAL) & (df_onset_hosp['onset_hosp'] <= MAX_VAL)]

    df_onset_death = pd.read_csv(DATA_PATH + 'onset_death_bog.csv')
    df_onset_death = df_onset_death[(df_onset_death['onset_death'] > MIN_VAL) & (df_onset_death['onset_death'] <= MAX_VAL)]

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

    if strat == 'wave':
        strat_= strat_wave
    elif strat == 'age':
        strat_= strat_ages
    elif strat == 'sex':
        strat_ = strat_sex

    columns = []
    for df in all_dfs:
        df.dropna(inplace=True) # remove the rows with nan values
        try:
            df.drop(columns = drop_columns, inplace = True)
        except:
            print('No columns to drop')
        col = str(df.columns[4])
        columns.append(col)
        df['age_group_id'] = df['age_group'].map(strat_ages_map)
        df['sex_id'] = df['sex'].map(strat_sex_map)
        df['wave_id'] = df['wave'].astype(int)
    
    return all_dfs, columns



def fit_district(values, list_of_params, model):
    stdata = values
    stan_data = {'N': len(stdata), 'y': stdata}
    fit = model.sampling(data = stan_data, iter = ITER, seed = SEED,
                        chains = CHAINS, n_jobs = -1)
    print(fit)                            
    df = fit.to_dataframe()
    df = df[list_of_params]
    return df

def get_posteriors_district(param_list, columns, all_dfs):
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
