# -*- coding: utf-8 -*-
"""
Created on Thr Jul 19 2022

@author: dsquevedo
@author: ntorres
"""
import warnings
warnings.filterwarnings('ignore')
import pandas as pd

DATA_PATH = '../results/'
OUT_PATH = '../results/'

ep_distributions  = {'icu_stay':{},
                     'hosp_stay':{},
                     'onset_icu':{},
                     'onset_hosp':{},
                     'onset_death':{}
                     }
best_models = pd.DataFrame({})

for dist in ep_distributions.keys():
    df = pd.read_csv(DATA_PATH+'bf_'+dist+'.csv')
    df = df.set_index(df.columns[0])
    ep_distributions[dist].update({'bf' : df,
                                   'best model' : df[df>=0].dropna()})
    best_models = pd.concat([best_models, ep_distributions[dist]['best model']])
cols = ['Epidemilogical distribution'] + best_models.columns.tolist()
best_models['Epidemilogical distribution'] = list(ep_distributions.keys()) 
best_models = best_models[cols]
best_models = best_models.set_index('Epidemilogical distribution')

best_models.to_csv(OUT_PATH+'best_models.csv')