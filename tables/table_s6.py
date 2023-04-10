# -*- coding: utf-8 -*-
'''
Created on Mon Nov 28 2022
@author: davidsantiagoquevedo
@author: ntorresd
'''  
import warnings
warnings.filterwarnings('ignore')

import sys
import yaml
import pandas as pd

config = yaml.load(open('config.yml', 'r'))['default']

RAT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'tables')

df_ratios = pd.read_csv(RAT_PATH + 'ratios.csv')
ratios = ['HCR', 'ICU-CR', 'CFR', 'HFR', 'ICU-FR']
df_final = df_ratios[['wave', 'age_group']]
df_final['wave'] = df_final['wave'].astype(int)
df_final.columns = ['Wave', 'Age group']

for col in ratios:
    df_final[col] = df_ratios[col].round(3).astype(str) \
                    + ' (' + df_ratios[col+'_lower'].round(3).astype(str) + ', ' \
                    + df_ratios[col+'_upper'].round(3).astype(str) + ')'

df_final.to_csv(OUT_PATH + 'table_s6.csv', index = False)