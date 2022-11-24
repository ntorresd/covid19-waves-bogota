# -*- coding: utf-8 -*-
"""
Created on Thr Jul 06 2022

@author: dsquevedo
@author: ntorresd
"""

# calculate counts and percentages by wave and age_group
def calculate_percentage(df, strat='wave', group='age_group', print_sum=False):
    df_percentage = df.groupby(by=[strat, group]).size().reset_index(name='counts')
    df_percentage['percentage'] = 0
    strat_list = df_percentage[strat].unique()
    for strat_ in strat_list:
        mask = (df_percentage[strat]==strat_)
        df_percentage.loc[mask, 'percentage'] = 100 * df_percentage.loc[mask, 'counts'] / df_percentage.loc[mask, 'counts'].sum()
        if print_sum:
            print('The sum over the percentages gives: ', df_percentage.loc[mask ,'percentage'].sum())
    return df_percentage
