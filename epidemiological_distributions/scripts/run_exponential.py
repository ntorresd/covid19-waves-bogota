# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 2022
Adapted from: https://github.com/mrc-ide/Brazil_COVID19_distributions

@author: dsquevedo
@author: ntorres
"""
import yaml
import pandas as pd
import pystan
import numpy as np
import scipy
import matplotlib.pyplot as plt
import time
from scipy import stats

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

