rm(list = ls())

library(rstudioapi) 
library(EpiEstim)
library(dplyr)
library(tidyverse)

library(config)
library(here)

# set the working directory 
setwd(here())

# paths
config <- config::get('PATHS')
DATA_PATH = config$DATA_PATH
OUT_PATH = config$OUT_PATH %>% str_replace('\\{dir\\}', 'rt')

# read data
df_confirmed_bogota <- read_csv(paste0(DATA_PATH, 'confirmed_cases.csv'))
df_confirmed_bogota <- df_confirmed_bogota[, c('onset', 'age', 'age_unit')]
names(df_confirmed_bogota)[names(df_confirmed_bogota) == 'onset'] <- 'dates'
df_confirmed_bogota$dates <- as.Date(df_confirmed_bogota$dates)

# define age groups
df_confirmed_bogota <- df_confirmed_bogota %>% 
  mutate(
  age_group = case_when(
    age_unit == 1 & age < 60 ~ '<60',
    age_unit == 2 ~ '<60',
    age_unit == 3 ~ '<60',
    age_unit == 1 & age >= 60 ~ '60+'
  ),
  age_group = factor(age_group, level=c('<60', '60+'))
)

# compute incidence
df_incidence <- df_confirmed_bogota %>% 
  count(dates, sort=TRUE, name='I')

mask = (df_confirmed_bogota$age_group == '60+')
df_incidence_60 <- df_confirmed_bogota[mask,] %>% 
  count(dates, sort=TRUE, name='I')
rm(mask)

####Completando las fechas faltantes####
complete_dates <- function(df){
  all_dates <- data.frame(dates = seq.Date(min(df$dates), 
                                           max(df$dates), 
                                           by = 'day')
                          )  
  
  df <- left_join(all_dates, df)
  df$I[is.na(df$I)] <- 0
  return(df)
} 

df_incidence <- df_incidence %>% complete_dates()
df_incidence_60 <- df_incidence_60 %>% complete_dates()

compute_rt <- function(df,
                       rt_days=14,
                       mean_si=6.48,
                       std_si=3.83){
  rt_days <- as.integer(rt_days)
  t_start <- seq_along(df$dates)[2:(length(df$dates) - rt_days)]
  t_end = t_start + rt_days
  rt_data  <- estimate_R(df, method = "parametric_si",
                         config = make_config(list(
                           mean_si = mean_si, std_si = std_si,
                           t_start= t_start, t_end = t_end )))
  
  df_rt <- rt_data$R
  df_rt$window_start <- min(df$dates) + df_rt$t_start
  df_rt$window_end <- min(df$dates) + df_rt$t_end
  df_rt <- df_rt %>% filter(window_end <= '2022-12-31')
  return(df_rt)
}

df_rt <- df_incidence %>% compute_rt()
df_rt_60 <- df_incidence_60 %>% compute_rt() 

# save data
write_csv(df_rt, paste0(OUT_PATH, 'rt_all_ages.csv'))
write_csv(df_rt_60, paste0(OUT_PATH, 'rt_60_plus.csv'))