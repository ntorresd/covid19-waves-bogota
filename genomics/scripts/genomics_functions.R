library(readxl)
library(openxlsx)
library(lubridate)
library(magrittr) #%>% function
library(MetBrewer)
library(tidyverse)
library(dplyr) #add_count
library(Hmisc)
library(reshape2)
library(ggplot2)
library(scales)
import::from(scales, rescale_none)

#### DICTIONARY FUNCTIONS #####
load_filters <- function(){
  ## Alpha ##
  filter_alpha <<- c("B.1.1.7", "^Q.")
  # ## Beta ## 
  # filter_beta <<- c("B.1.351", "^B.1.351.")
  ## Gamma ## 
  filter_gamma <<- c("P.1", "^P.1.")
  ## Delta ## 
  filter_delta <<- c("B.1.617.2", "^AY.")
  ## Lambda ## 
  # filter_lambda <<- c("C.37", "C.37.1")
  ## Mu ## 
  filter_mu <<- c("B.1.621", "B.1.621.1", "B.1.621.2", "BB.1", "BB.2")
  ## Omicron ##
  filter_omicron <<- c("B.1.1.529", "BA.1", "^BA.1.", "^BA.1.1","^BA.1.1", "BA.2","^BA.2.", 
                       "BA.3", "BA.4", "BA.5", "B.2.12.1")
}

add_variant_dictionary <- function(dic_variants,
                             list_lineages = list_lineages,
                             variant,
                             filter_variant){
  list_pos_variant <- grep(paste(filter_variant, collapse="|"), 
                             list_lineages)
  tabla_aux <- cbind(variant, list_lineages[list_pos_variant])
  dic_variants <- dic_variants %>% add_row(label = tabla_aux[,1],
                                             lineage = tabla_aux[,2])
  return (dic_variants)
}

generate_dictionary <- function(Base_genomas){
  list_lineages <- unique(Base_genomas$lineage) 
  dic_variants <- data.frame(label = character(),
                              lineage = character())
  
  dic_variants <- add_variant_dictionary (dic_variants = dic_variants,
                                     list_lineages = list_lineages,
                                     variant = "Alpha",
                                     filter_variant = filter_alpha)
  
  # dic_variants <- add_variant_dictionary (dic_variants = dic_variants,
  #                                        list_lineages = list_lineages,
  #                                        variant = "Beta",
  #                                        filter_variant = filter_beta)
  
  
  dic_variants <- add_variant_dictionary (dic_variants = dic_variants,
                                     list_lineages = list_lineages,
                                     variant = "Gamma",
                                     filter_variant = filter_gamma)
  
  
  dic_variants <- add_variant_dictionary (dic_variants = dic_variants,
                                     list_lineages = list_lineages,
                                     variant = "Delta",
                                     filter_variant = filter_delta)
  
  # dic_variants <- add_variant_dictionary (dic_variants = dic_variants,
  #                                        list_lineages = list_lineages,
  #                                        variant = "Lambda",
  #                                        filter_variant = filter_lambda)
  
  
  dic_variants <- add_variant_dictionary (dic_variants = dic_variants,
                                     list_lineages = list_lineages,
                                     variant = "Mu",
                                     filter_variant = filter_mu)
  
  
  dic_variants <- add_variant_dictionary (dic_variants = dic_variants,
                                     list_lineages = list_lineages,
                                     variant = "Omicron",
                                     filter_variant = filter_omicron)
  ## Others ##
  
  dic_variants <- left_join(data.frame(lineage = list_lineages), 
                             dic_variants, by = "lineage") %>% 
    # filter(lineage != "None") %>% filter(!is.na(lineage)) %>%
    mutate(label = replace_na(label, "Other")) %>% arrange(label)
  
  dic_variants <- dic_variants[, c(2,1)]
  return(dic_variants)
}

#### CODIFICATION FUNCTIONS AND GENOMES AGGRUPATIONS ####
filter_genomes <- function(df_genomas,
                            Variable_filtros = "Lugar",
                            filter_list){
  filter_positions <- grep(paste(filter_list, collapse="|"), df_genomas[, Variable_filtros][[1]])
  df_genomas <- df_genomas[filter_positions,]
  return (df_genomas)
}

lineage_codification <- function(df_genomas, 
                         dic_variants){
  df_genomas <- df_genomas %>% filter(lineage != "None") %>% filter(!is.na(lineage))
  df_genomas$lineage_recod = NA
  for (i in c(1:length(dic_variants$lineage))){
    df_genomas$lineage_recod[df_genomas$lineage == dic_variants$lineage[i]] = dic_variants$label[i]
  }
  return(df_genomas)
}

group_variant_epiweek <- function(df_genomas_recod, initial_date){
  df_genomas <- df_genomas_recod %>% add_column(epiweek = NA, .after = "date") %>% 
    add_column(year = NA, .after = "date")
  df_genomas$year <- format(df_genomas$date, format="%Y")
  df_genomas$epiweek <- epiweek(df_genomas$date)
  
  index <- (df_genomas$year == "2022" & df_genomas$date < as.Date("2022-01-31") & df_genomas$epiweek == 52)
  df_genomas$year[index] <- "2021"
  rm(index)
  
  df_aux <- aggregate(df_genomas$lineage_recod, by=list(df_genomas$year, 
                                                        df_genomas$epiweek,
                                                        df_genomas$lineage_recod), FUN = length)
  df_aux_tot <- aggregate(df_aux$x, by=list(df_aux$Group.1,
                                            df_aux$Group.2), FUN=sum) 
  
  df_aux <- df_aux[order(df_aux$Group.1, df_aux$Group.2), ]
  df_aux_tot <- df_aux_tot[order(df_aux_tot$Group.1, df_aux_tot$Group.2), ]
  
  date_list <- seq(initial_date, length.out = length(df_aux_tot$Group.1), by = 7)
  df_aux_tot <- data.frame(Group.0 = date_list, df_aux_tot)
  
  df_grouped <- left_join(df_aux, df_aux_tot, by=c('Group.1', 'Group.2'))
  df_grouped <- df_grouped[, c("Group.0", "Group.1", "Group.2", "Group.3", "x.x", "x.y")]
  colnames(df_grouped) <- c('date', 'year', 'week', 'lineage', 'n_seq_var', 'n_seq_tot')
  
  df_grouped <- cbind(df_grouped, (binconf(df_grouped$n_seq_var, df_grouped$n_seq_tot)))
  
  rm(df_aux, df_aux_tot)
  return(df_grouped)
}