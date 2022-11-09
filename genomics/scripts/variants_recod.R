rm(list= ls()) ## Borra los objetos anteriores

library(configr)
library(here)

#### PATHS ####
setwd(here())
DATA_PATH = "data/"
OUT_PATH = "genomics/outputs/"

#### LOAD FUNCTIONS ####
source("genomics/scripts/genomics_functions.R")
`%!in%` <- Negate(`%in%`)


#### DATES ####
initial_date = as.Date("2021-03-22") 

#### READ DATA ####

df_genomes_bog <- read_csv(paste0(DATA_PATH,'genomes_bog_raw.csv'))

df_genomes_bog$date <- df_genomes_bog$date %>% as.Date()
df_genomes_bog <- df_genomes_bog[df_genomes_bog$date >= initial_date, ]
colnames(df_genomes_bog) <- c('date', 'lineage')

load_filters()
dic_var_bog <- generate_dictionary(df_genomes_bog)
rm(list = ls(pat = "^filtro"))
df_genomes_bog <- lineage_codification(df_genomas = df_genomes_bog, dic_variants=dic_var_bog)

#### HISTOGRAM VARIANTS BOGOTÁ ####

df_grouped <- group_variant_epiweek(df_genomas_recod = df_genomes_bog, initial_date = initial_date)
df_grouped <- df_grouped %>% 
  mutate(lineage=factor(lineage, levels = c("Alpha", "Gamma", "Mu", "Other", "Delta", "Omicron"))) %>% 
  arrange(date, lineage)

write.csv(df_grouped, paste0(OUT_PATH, "/variants-ic-bog.csv"), row.names = FALSE)

unique(df_genomes_bog[(df_genomes_bog$lineage_recod=="Other") & (df_genomes_bog$date >="2021-12-16"),]$lineage)
unique(df_genomes_bog[(df_genomes_bog$lineage_recod=="Omicron"),]$lineage)
