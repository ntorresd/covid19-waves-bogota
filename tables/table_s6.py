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

SEV_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'tables')

# Ratios
df_ratios = pd.read_csv(SEV_PATH + 'ratios.csv')
df_ratios['wave'] = df_ratios['wave'].astype(int)
ratios = ['HCR', 'ICU-CR', 'CFR', 'HFR', 'ICU-FR']

df_ratios_err = df_ratios[['wave', 'age_group']]
df_ratios_err.columns = ['Wave', 'Age group']

for col in ratios:
    df_ratios_err[col] = df_ratios[col].round(3).astype(str) \
                        + ' (' + df_ratios[col+'_lower'].round(3).astype(str) + ', ' \
                        + df_ratios[col+'_upper'].round(3).astype(str) + ')'

# Percentages
df_percentages = pd.read_csv(SEV_PATH + 'percentages.csv')
df_percentages['wave'] = df_percentages['wave'].astype(int)
percentages = ['HOSP-%', 'ICU-%', 'DEATH-%']

df_percentages_err = df_percentages[['wave', 'age_group']]
df_percentages_err.columns = ['Wave', 'Age group']

for col in percentages:
    df_percentages_err[col] = df_percentages[col].round(3).astype(str) \
                            + ' (' + df_percentages[col+'_lower'].round(3).astype(str) + ', ' \
                            + df_percentages[col+'_upper'].round(3).astype(str) + ')'

# Merge results
df_final = df_ratios_err.merge(df_percentages_err, how = 'left', on = ['Wave', 'Age group']).fillna('-')

df_final.to_csv(OUT_PATH + 'table_s6.csv', index = False)