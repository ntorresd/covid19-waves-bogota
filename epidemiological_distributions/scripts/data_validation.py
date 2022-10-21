# -*- coding: utf-8 -*-
"""
Created on Thr Jul 07 2022

@author: dsquevedo
@author: ntorres
"""
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from time import time
import datetime as dt
import matplotlib.pyplot as plt
plt.style.use('../../plot_scripts/plot_style.mplstyle')

DATA_PATH = '../../process_data/data/'
FIG_PATH = '../figures/'

# Reported occupancy in Bogotá
# Data sources: https://datosabiertos.bogota.gov.co, https://saludata.saludcapital.gov.co
df_hosp_oc_bogota = pd.read_csv(DATA_PATH+'hosp_occupancy_bog.csv',dtype=object)
df_icu_oc_bogota = pd.read_csv(DATA_PATH+'icu_occupancy_bog.csv',dtype=object)
    
## Data processing

df_hosp_oc_bogota['date']=pd.to_datetime(df_hosp_oc_bogota['date'],errors='coerce')
df_icu_oc_bogota['date']=pd.to_datetime(df_icu_oc_bogota['date'],errors='coerce')
df_hosp_oc_bogota['occupancy']=pd.to_numeric(df_hosp_oc_bogota['occupancy'])
df_icu_oc_bogota['occupancy']=pd.to_numeric(df_icu_oc_bogota['occupancy'])
df_hosp_oc_bogota=df_hosp_oc_bogota.sort_values(by='date')
df_icu_oc_bogota=df_icu_oc_bogota.sort_values(by='date')

# ICU and Hospital stay in Bogotá
df_hosp_stay_bogota = pd.read_csv(DATA_PATH+'hosp_stay_bog.csv')
df_icu_stay_bogota = pd.read_csv(DATA_PATH+'icu_stay_bog.csv')
df_hosp_stay_bogota['start_date']=pd.to_datetime(df_hosp_stay_bogota['start_date'],errors='coerce')
df_hosp_stay_bogota['end_date']=pd.to_datetime(df_hosp_stay_bogota['end_date'],errors='coerce')
df_icu_stay_bogota['start_date']=pd.to_datetime(df_icu_stay_bogota['start_date'],errors='coerce')
df_icu_stay_bogota['end_date']=pd.to_datetime(df_icu_stay_bogota['end_date'],errors='coerce')

## Occupancy calculated from confirmed cases
def calculate_occupancy(df,var_start,var_end,start_date,days):
    start_date=pd.to_datetime(start_date)
    n_days=pd.Series(np.linspace(0,days,days+1))
    days = n_days.apply(lambda x: start_date + dt.timedelta(days=x))
    df_occupancy={'date':[],'occupancy':[]}
    for date in days:
        occ=len(df[(df[var_start]<=date)&(df[var_end]>=date)])
        df_occupancy['date'].append(date)
        df_occupancy['occupancy'].append(occ)
    df_occupancy=pd.DataFrame(df_occupancy)
    df_occupancy=df_occupancy.sort_values(by='date')
    return df_occupancy
df_icu_occ_calc=calculate_occupancy(df_icu_stay_bogota,'start_date','end_date','2020-03-01',days=900)
df_hosp_occ_calc=calculate_occupancy(df_hosp_stay_bogota,'start_date','end_date','2020-03-01',days=900)

## Plots
fig,ax=plt.subplots(1,2,figsize=(14,7))
max_date=pd.to_datetime('2022-04-01')

axi=ax[0]
min_date=df_icu_oc_bogota.date.min()
df_hosp_oc_bogota.set_index('date')[['occupancy']].plot(ax=axi)
df_hosp_occ_calc.set_index('date')[['occupancy']].plot(ax=axi)
axi.set_title('a.')
axi.set_xlabel('')
axi.set_xlim([min_date,max_date])
axi.get_legend().remove()

axi=ax[1]
min_date=df_icu_oc_bogota.date.min()
min_date=df_icu_oc_bogota.date.min()
df_icu_oc_bogota.set_index('date')[['occupancy']].plot(ax=axi)
df_icu_occ_calc.set_index('date')[['occupancy']].plot(ax=axi)
axi.set_title('b.')
axi.set_xlabel('')
axi.set_xlim([min_date,max_date])
axi.get_legend().remove()

handles, labels = axi.get_legend_handles_labels()
fig.legend(handles, ['Reported - SaluData Bogotá','Confirmed cases'], bbox_to_anchor=(0.8,-0.02), ncol=2)
fig.savefig('../figures/icu_hosp_validation.png')