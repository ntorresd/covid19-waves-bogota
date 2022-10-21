rm(list= ls()) ## Borra los objetos anteriores
# setwd("~/ASIS/Variantes/Prevalencia")

#### CARGAR FUNCIONES ####
source("genomics/functions/genomics_functions.R")
`%!in%` <- Negate(`%in%`)

#### PARAMETROS DE MODIFICACION ####
fecha_corte <-  '2022-08-17'
fecha_corte_ <- gsub("-","",fecha_corte)

fecha_inicial = as.Date("2021-08-22") 
fecha_inicial_hist <- as.Date("2021-03-26")

df_genomas_bog <- read_csv("data/genomes_bog_raw.csv")

df_genomas_bog$date <- df_genomas_bog$date %>% as.Date()
df_genomas_bog <- df_genomas_bog[df_genomas_bog$date >= fecha_inicial, ]
colnames(df_genomas_bog) <- c('date', 'lineage')

load_filters()
dic_var_bog <- generate_dictionary(df_genomas_bog)
rm(list = ls(pat = "^filtro"))
df_genomas_bog <- lineage_codification(df_genomas = df_genomas_bog, dic_variantes=dic_var_bog)

#### HISTOGRAMA VARIANTES BOGOTÁ ####
ruta_fig_hist_bog <- paste0('genomics/figures/prev_var_hist_bog_', fecha_corte_, '.png')

df_agrupada <- group_variant_epiweek(df_genomas_recod = df_genomas_bog, fecha_inicial = fecha_inicial_hist)
df_agrupada <- df_agrupada %>% 
  mutate(lineage=factor(lineage, levels = c("Alpha", "Gamma", "Mu", "Other", "Delta", "Omicron"))) %>% 
  arrange(date, lineage)

plot_prev_bog <- plot_prevalence_histogram(df_agrupada = df_agrupada,
                                           Fecha_corte = fecha_corte,
                                           Titulo = "Prevalencia por variantes - Bogotá",
                                           Paleta = "Thomas")

ggsave(plot = plot_prev_bog, ruta_fig_hist_bog, width = 35, height = 20, units = "cm")
write.csv(df_agrupada, paste0("genomics/outputs/variants-ic-bog_", fecha_corte_,".csv"), row.names = FALSE)

unique(df_genomas_bog[(df_genomas_bog$lineage_recod=="Other") & (df_genomas_bog$date >="2021-12-16"),]$lineage)
unique(df_genomas_bog[(df_genomas_bog$lineage_recod=="Omicron"),]$lineage)

#### HISTOGRAMA OMICRON BOGOTÁ BA.1/2/3 Y ^B.2. ####
ruta_fig_hist_omicron_bog <- paste0('genomics/figures/omicron_bog_', fecha_corte_, '.png')

df_omicron <- df_genomas_bog[df_genomas_bog$lineage_recod=="Omicron", c('date', 'lineage')]
# fecha_inicial_hist = min(df_genomas$date)
unique(df_omicron$lineage)
filtro_BA1 <- c("BA.1", "^BA.1.")
filtro_BA2 <- c("BA.2", "^BA.2.")
filtro_BA3 <- c("BA.3", "^BA.3.")
filtro_BA4 <- c("BA.4", "^BA.4")
filtro_BA5 <- c("BA.5", "^BA.5")
filtro_B2_12 <- c("^B.2.")

list_lineages <- unique(df_omicron$lineage) 
dic_omicron <- data.frame(label = character(),lineage = character())
dic_omicron <- add_variant_dictionary(dic_variantes = dic_omicron, 
                                list_lineages = list_lineages, 
                                variante="BA.1", 
                                filtro_variante=filtro_BA1)
dic_omicron <- add_variant_dictionary(dic_variantes = dic_omicron, 
                                list_lineages = list_lineages, 
                                variante="BA.2", 
                                filtro_variante=filtro_BA2)
# dic_omicron <- add_variant_dictionary(dic_variantes = dic_omicron, 
#                                 list_lineages = list_lineages, 
#                                 variante="BA.3", 
#                                 filtro_variante=filtro_BA3)
dic_omicron <- add_variant_dictionary(dic_variantes = dic_omicron, 
                                list_lineages = list_lineages, 
                                variante="BA.4", 
                                filtro_variante=filtro_BA4)
dic_omicron <- add_variant_dictionary(dic_variantes = dic_omicron, 
                                list_lineages = list_lineages, 
                                variante="BA.5", 
                                filtro_variante=filtro_BA5)

dic_omicron <- left_join(data.frame(lineage = list_lineages), 
                         dic_omicron, by = "lineage") %>% 
  mutate(label = replace_na(label, "Other")) %>% arrange(label)
dic_omicron <- dic_omicron[, c(2,1)]

df_omicron <- lineage_codification(df_genomas = df_omicron, dic_variantes=dic_omicron)
rm(list = ls(pat = "^filtro"))

df_agrupada <- group_variant_epiweek(df_genomas_recod = df_omicron, fecha_inicial = min(df_omicron$date))
plot_omicron <- plot_prevalence_histogram(df_agrupada = df_agrupada,
                                          Fecha_corte = fecha_corte,
                                          Titulo = "Subvariantes Omicron - Bogotá",
                                          Paleta = "Cross")
ggsave(plot = plot_omicron, ruta_fig_hist_omicron_bog,
       width = 35, height = 20, units = "cm")


#### HISTOGRAMA PREVALENCIA BOGOTÁ (CUSTOM)####
ruta_fig_hist_bog_custom <- paste0('genomics/figures/prev_var_hist_bog_custom_', fecha_corte_, '.png')

##FILTROS##
## Alpha ##
filtro_alpha <- c("B.1.1.7", "^Q.")
## Gamma ##
filtro_gamma <- c("P.1", "^P.1.")
## Delta ##
filtro_delta <- c("B.1.617.2", "^AY.")
## Mu ##
filtro_mu <- c("B.1.621", "B.1.621.1", "B.1.621.2", "BB.1", "BB.2")
## Omicron ##
filtro_ba1 <- c("BA.1", "^BA.1.")
filtro_ba2 <- c("BA.2", "^BA.2.")
# filtro_ba3 <- c("BA.3", "^BA.3.")
filtro_ba4 <- c("BA.4", "^BA.4.")
filtro_ba5 <- c("BA.5", "^BA.5.")
# filtro_b2_12 <- c("B.2.12", "^B.2.12.")

# list_filtros = list(filtro_alpha, filtro_gamma, filtro_delta, filtro_mu, filtro_ba1, filtro_ba2, filtro_ba4, filtro_ba5)
# list_variantes = list("Alpha", "Gamma", "Delta", "Mu", "BA.1", "BA.2", "BA.4", "BA.5")

list_filtros = list(filtro_gamma, filtro_delta, filtro_mu, filtro_ba1, filtro_ba2, filtro_ba4, filtro_ba5)
list_variantes = list("Gamma", "Delta", "Mu", "BA.1", "BA.2", "BA.4", "BA.5")

list_lineages <- unique(df_genomas_bog$lineage)
dic_custom <- data.frame(label = character(),lineage = character())
for (i in 1:length(list_filtros)){
  dic_custom <- add_variant_dictionary(dic_variantes = dic_custom,
                                  list_lineages = list_lineages,
                                  variante=list_variantes[[i]],
                                  filtro_variante=list_filtros[[i]])
}
dic_custom <- left_join(data.frame(lineage = list_lineages),
                         dic_custom, by = "lineage") %>%
  mutate(label = replace_na(label, "Other")) %>% arrange(label)
dic_custom <- dic_custom[, c(2,1)]

df_custom <- lineage_codification(df_genomas = df_genomas_bog, dic_variantes=dic_custom)

df_agrupada <- group_variant_epiweek(df_genomas_recod = df_custom, fecha_inicial = fecha_inicial_hist)
df_agrupada <- df_agrupada %>%
  mutate(lineage=factor(lineage, levels = c("Alpha", "Gamma", "Mu", "Other", "Delta", "BA.1", "BA.2", "BA.4", "BA.5"))) %>%
  arrange(date, lineage)
rm(list_filtros, list_variantes)
rm(list = ls(pat = "^filtro"))

plot_prev_custom <- plot_prevalence_histogram(df_agrupada = df_agrupada,
                                           Fecha_corte = fecha_corte,
                                           Titulo = "Prevalencia por variantes - Bogotá (custom)",
                                           Paleta = "Thomas")

ggsave(plot = plot_prev_custom, ruta_fig_hist_bog_custom,
       width = 35, height = 20, units = "cm")
# write.csv(df_agrupada, paste0("genomics/outputs/variants-ic-bog_custom_", fecha_corte_,".csv"), row.names = FALSE)

#### HISTOGRAMA OTHER BOGOTÁ ####
ruta_fig_hist_other_bog <- paste0('genomics/figures/other_bog_', fecha_corte_, '.png')

df_other <- df_genomas_bog[(df_genomas_bog$lineage_recod=="Other") & (df_genomas_bog$date >= "2022-03-27"), 
                           c('date', 'lineage')]
colnames(df_other) <- c("date", "lineage_recod")
unique(df_other$lineage_recod)
list_lineages <- unique(df_other$lineage_recod) 
df_agrupada <- group_variant_epiweek(df_genomas_recod = df_other, fecha_inicial = min(df_other$date))
plot_other <- plot_prevalence_histogram(df_agrupada = df_agrupada,
                                          Fecha_corte = fecha_corte,
                                          Titulo = "Otras Variantes - Bogotá",
                                          Paleta = "Cross")
ggsave(plot = plot_other, ruta_fig_hist_other_bog,
       width = 35, height = 20, units = "cm")
