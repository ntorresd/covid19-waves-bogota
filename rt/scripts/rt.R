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

# parameters 
incubation_period <- 5
rt_window <- 7

# read data
df_confirmed_bogota <- read_csv(paste0(DATA_PATH, 'confirmed_cases.csv'))
df_confirmed_bogota <- df_confirmed_bogota[, c('onset', 'age', 'age_unit')]
df_confirmed_bogota$onset <- as.Date(df_confirmed_bogota$onset)
df_confirmed_bogota$infection=df_confirmed_bogota$onset - incubation_period

# define age groups
df_confirmed_bogota <- df_confirmed_bogota %>% 
  mutate(
    age_group = case_when(
      age_unit == 1 & age < 60 ~ "<60",
      age_unit == 2 ~ "<60",
      age_unit == 3 ~ "<60",
      age_unit == 1 & age >= 60 ~ "60+"
    ),
    age_group = factor(age_group, level=c("<60", "60+"))
  )

# compute incidence

df_incidence <- df_confirmed_bogota %>% 
  group_by(infection) %>% summarise(I = n())

mask = (df_confirmed_bogota$age_group == "60+")
df_incidence_60p <- df_confirmed_bogota[mask, ] %>%
  group_by(infection) %>% summarise(I = n())
rm(mask)

# complete missing dates
complete_dates <- function(df){
  all_dates <- data.frame(infection = seq.Date(min(df$infection), 
                                               max(df$infection), 
                                               by = 'day')
  )  
  
  df <- left_join(all_dates, df)
  df$I[is.na(df$I)] <- 0
  return(df)
} 

df_incidence <- df_incidence %>% complete_dates()
df_incidence_60p <- df_incidence_60p %>% complete_dates()

# compute Rt
compute_rt <- function(df_incidence,
                       method = "parametric_si",
                       mean_si=6.48,
                       std_si=3.83){
  t_start <- seq(incubation_period, nrow(df_incidence) - rt_window)
  t_end <- t_start + rt_window
  rt_data <- estimate_R(df_incidence, 
                        method = method,
                        config = make_config(list(
                          mean_si = mean_si,
                          std_si = std_si,
                          t_start = t_start,
                          t_end = t_end)))
  
  df_rt <- rt_data$R
  df_rt$window_start <- min(df_incidence$infection) + df_rt$t_start 
  df_rt$window_end <- min(df_incidence$infection) + df_rt$t_end 

  return(df_rt)  
}

df_rt <- df_incidence %>% compute_rt() %>%
  filter(window_start >= as.Date("2020-03-01"))
df_rt_60p <- df_incidence_60p %>% compute_rt() %>%
  filter(window_start >= as.Date("2020-03-01"))

# remove atypical data first wave
remove_atypical_rt <- function(df_rt,
                               initial_date = as.Date("2020-03-01"),
                               final_date = as.Date("2020-09-25")) {
df_rt_wave <- df_rt %>%
  filter(window_start >= initial_date & window_start <= final_date)

typical_min <- quantile(df_rt_wave$`Mean(R)`)[[2]] - 1.5*IQR(df_rt_wave$`Mean(R)`)
typical_max <- quantile(df_rt_wave$`Mean(R)`)[[4]] + 1.5*IQR(df_rt_wave$`Mean(R)`)

df_rt_wave <- df_rt_wave %>%
  filter(`Mean(R)` >= typical_min & `Mean(R)` <= typical_max)

df_rt_typical <- df_rt %>%
  filter(window_start < initial_date | window_start > final_date) %>%
  full_join(df_rt_wave) %>%
  arrange(t_start)

return(df_rt_typical)
}

df_rt <- df_rt %>% remove_atypical_rt()
df_rt_60p <- df_rt_60p %>% remove_atypical_rt()

# save data
write_csv(df_rt, paste0(OUT_PATH, "rt_all_ages.csv"))
write_csv(df_rt_60p, paste0(OUT_PATH, "rt_60_plus.csv"))