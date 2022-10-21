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
  filtro_alpha <<- c("B.1.1.7", "^Q.")
  # ## Beta ## 
  # filtro_beta <<- c("B.1.351", "^B.1.351.")
  ## Gamma ## 
  filtro_gamma <<- c("P.1", "^P.1.")
  ## Delta ## 
  filtro_delta <<- c("B.1.617.2", "^AY.")
  ## Lambda ## 
  # filtro_lambda <<- c("C.37", "C.37.1")
  ## Mu ## 
  filtro_mu <<- c("B.1.621", "B.1.621.1", "B.1.621.2", "BB.1", "BB.2")
  ## Omicron ##
  filtro_omicron <<- c("B.1.1.529", "BA.1", "^BA.1.", "^BA.1.1","^BA.1.1", "BA.2","^BA.2.", 
                       "BA.3", "BA.4", "BA.5", "B.2.12.1")
}

add_variant_dictionary <- function(dic_variantes,
                             list_lineages = list_lineages,
                             variante,
                             filtro_variante){
  list_pos_variante <- grep(paste(filtro_variante, collapse="|"), 
                             list_lineages)
  tabla_aux <- cbind(variante, list_lineages[list_pos_variante])
  dic_variantes <- dic_variantes %>% add_row(label = tabla_aux[,1],
                                             lineage = tabla_aux[,2])
  return (dic_variantes)
}

generate_dictionary <- function(Base_genomas){
  list_lineages <- unique(Base_genomas$lineage) 
  dic_variantes <- data.frame(label = character(),
                              lineage = character())
  
  # dic_variantes <- add_variant_dictionary (dic_variantes = dic_variantes,
  #                                    list_lineages = list_lineages,
  #                                    variante = "Alpha",
  #                                    filtro_variante = filtro_alpha)
  
  # dic_variantes <- add_variant_dictionary (dic_variantes = dic_variantes,
  #                                        list_lineages = list_lineages,
  #                                        variante = "Beta",
  #                                        filtro_variante = filtro_beta)
  
  
  dic_variantes <- add_variant_dictionary (dic_variantes = dic_variantes,
                                     list_lineages = list_lineages,
                                     variante = "Gamma",
                                     filtro_variante = filtro_gamma)
  
  
  dic_variantes <- add_variant_dictionary (dic_variantes = dic_variantes,
                                     list_lineages = list_lineages,
                                     variante = "Delta",
                                     filtro_variante = filtro_delta)
  
  # dic_variantes <- add_variant_dictionary (dic_variantes = dic_variantes,
  #                                        list_lineages = list_lineages,
  #                                        Variante = "Lambda",
  #                                        filtro_variante = filtro_lambda)
  
  
  dic_variantes <- add_variant_dictionary (dic_variantes = dic_variantes,
                                     list_lineages = list_lineages,
                                     variante = "Mu",
                                     filtro_variante = filtro_mu)
  
  
  dic_variantes <- add_variant_dictionary (dic_variantes = dic_variantes,
                                     list_lineages = list_lineages,
                                     variante = "Omicron",
                                     filtro_variante = filtro_omicron)
  ## Others ##
  
  dic_variantes <- left_join(data.frame(lineage = list_lineages), 
                             dic_variantes, by = "lineage") %>% 
    # filter(lineage != "None") %>% filter(!is.na(lineage)) %>%
    mutate(label = replace_na(label, "Other")) %>% arrange(label)
  
  dic_variantes <- dic_variantes[, c(2,1)]
  return(dic_variantes)
}

#### CODIFICATION FUNCTIONS AND GENOMES AGGRUPATIONS ####
filter_genomes <- function(df_genomas,
                            Variable_filtros = "Lugar",
                            list_filtros){
  Filtros_pos <- grep(paste(list_filtros, collapse="|"), df_genomas[, Variable_filtros][[1]])
  df_genomas <- df_genomas[Filtros_pos,]
  return (df_genomas)
}

lineage_codification <- function(df_genomas, 
                         dic_variantes){
  df_genomas <- df_genomas %>% filter(lineage != "None") %>% filter(!is.na(lineage))
  df_genomas$lineage_recod = NA
  for (i in c(1:length(dic_variantes$lineage))){
    df_genomas$lineage_recod[df_genomas$lineage == dic_variantes$lineage[i]] = dic_variantes$label[i]
  }
  return(df_genomas)
}

group_variant_epiweek <- function(df_genomas_recod, fecha_inicial){
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
  
  date_list <- seq(fecha_inicial, length.out = length(df_aux_tot$Group.1), by = 7)
  df_aux_tot <- data.frame(Group.0 = date_list, df_aux_tot)
  
  df_agrupada <- left_join(df_aux, df_aux_tot, by=c('Group.1', 'Group.2'))
  df_agrupada <- df_agrupada[, c("Group.0", "Group.1", "Group.2", "Group.3", "x.x", "x.y")]
  colnames(df_agrupada) <- c('date', 'year', 'week', 'lineage', 'n_seq_var', 'n_seq_tot')
  
  df_agrupada <- cbind(df_agrupada, (binconf(df_agrupada$n_seq_var, df_agrupada$n_seq_tot)))
  
  rm(df_aux, df_aux_tot)
  return(df_agrupada)
}


#### PLOTTING FUNCTIONS ####
plot_prevalence_histogram <- function(df_agrupada,
                                      Fecha_corte = Fecha_GISAID,
                                      Titulo = "Prevalencia",
                                      Paleta = "Thomas"){
  Plot_prevalencia <- ggplot(df_agrupada, aes (x = as.factor(date), fill = lineage)) +
    ggtitle(Titulo) +
    geom_bar(aes (y = PointEst), stat = "identity", alpha=0.6, colour = "black")+
    geom_text(aes (y=1.035, label = n_seq_tot), hjust=0.5, alpha = 0.8, angle = 90) +
    geom_text(aes(y = PointEst, label = n_seq_var), position = position_stack(vjust = 0.5), angle = 90) +
    labs (caption=paste0("Corte al ", Fecha_corte,",\t epiweek = ", epiweek(Fecha_corte)),
          fill ="Variante", 
          x = "Fecha de recolección", 
          y= "Prevalencia por variante / Linaje viral") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
    theme(legend.position="bottom",
          axis.text.y = element_text(margin = margin(r = 0)),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_y_continuous(limits=c(0.045, 1.03), oob=rescale_none)+
    scale_fill_manual(values = met.brewer(Paleta, n = length(unique(df_agrupada$lineage))))
  # df_agrupada <- df_agrupada[(df_agrupada$n_seq_tot>=Secuencias_min),]
  # df_agrupada$date <- as.Date(df_agrupada$date)
  return (Plot_prevalencia)
}