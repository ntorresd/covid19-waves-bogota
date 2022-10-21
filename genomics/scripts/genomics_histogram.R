rm(list= ls()) ## Borra los objetos anteriores

#### CARGAR FUNCIONES ####
source("genomics/functions/genomics_functions.R")
`%!in%` <- Negate(`%in%`)

#### PARAMETROS DE MODIFICACION ####
download_date <-  '2022-08-17'
download_date_ <- gsub("-","",download_date)

initial_date = as.Date("2021-08-22") 
initial_date_hist <- as.Date("2021-03-26")

df_genomes_bog <- read_csv("data/genomes_bog_raw.csv")

df_genomes_bog$date <- df_genomes_bog$date %>% as.Date()
df_genomes_bog <- df_genomes_bog[df_genomes_bog$date >= initial_date, ]
colnames(df_genomes_bog) <- c('date', 'lineage')

load_filters()
dic_var_bog <- generate_dictionary(df_genomes_bog)
rm(list = ls(pat = "^filtro"))
df_genomes_bog <- lineage_codification(df_genomas = df_genomes_bog, dic_variants=dic_var_bog)

#### HISTOGRAMA VARIANTS BOGOTÁ ####
path_fig_hist_bog <- paste0('genomics/figures/prev_var_hist_bog_', download_date_, '.png')

df_grouped <- group_variant_epiweek(df_genomas_recod = df_genomes_bog, initial_date = initial_date_hist)
df_grouped <- df_grouped %>% 
  mutate(lineage=factor(lineage, levels = c("Alpha", "Gamma", "Mu", "Other", "Delta", "Omicron"))) %>% 
  arrange(date, lineage)

plot_prev_bog <- plot_prevalence_histogram(df_grouped = df_grouped,
                                           download_date = download_date,
                                           title = "Prevalence per variant - Bogotá",
                                           pallete = "Thomas")

ggsave(plot = plot_prev_bog, path_fig_hist_bog, width = 35, height = 20, units = "cm")
write.csv(df_grouped, paste0("genomics/outputs/variants-ic-bog_", download_date_,".csv"), row.names = FALSE)

unique(df_genomes_bog[(df_genomes_bog$lineage_recod=="Other") & (df_genomes_bog$date >="2021-12-16"),]$lineage)
unique(df_genomes_bog[(df_genomes_bog$lineage_recod=="Omicron"),]$lineage)

#### OMICRON HISTOGRAM BOGOTÁ BA.1/2/3 Y ^B.2. ####
path_fig_hist_omicron_bog <- paste0('genomics/figures/omicron_bog_', download_date_, '.png')

df_omicron <- df_genomes_bog[df_genomes_bog$lineage_recod=="Omicron", c('date', 'lineage')]
# initial_date_hist = min(df_genomas$date)
unique(df_omicron$lineage)
filter_BA1 <- c("BA.1", "^BA.1.")
filter_BA2 <- c("BA.2", "^BA.2.")
filter_BA3 <- c("BA.3", "^BA.3.")
filter_BA4 <- c("BA.4", "^BA.4")
filter_BA5 <- c("BA.5", "^BA.5")
filter_B2_12 <- c("^B.2.")

list_lineages <- unique(df_omicron$lineage) 
dic_omicron <- data.frame(label = character(),lineage = character())
dic_omicron <- add_variant_dictionary(dic_variants = dic_omicron, 
                                list_lineages = list_lineages, 
                                variant="BA.1", 
                                filter_variant=filter_BA1)
dic_omicron <- add_variant_dictionary(dic_variants = dic_omicron, 
                                list_lineages = list_lineages, 
                                variant="BA.2", 
                                filter_variant=filter_BA2)
# dic_omicron <- add_variant_dictionary(dic_variants = dic_omicron, 
#                                 list_lineages = list_lineages, 
#                                 variant="BA.3", 
#                                 filter_variant=filter_BA3)
dic_omicron <- add_variant_dictionary(dic_variants = dic_omicron, 
                                list_lineages = list_lineages, 
                                variant="BA.4", 
                                filter_variant=filter_BA4)
dic_omicron <- add_variant_dictionary(dic_variants = dic_omicron, 
                                list_lineages = list_lineages, 
                                variant="BA.5", 
                                filter_variant=filter_BA5)

dic_omicron <- left_join(data.frame(lineage = list_lineages), 
                         dic_omicron, by = "lineage") %>% 
  mutate(label = replace_na(label, "Other")) %>% arrange(label)
dic_omicron <- dic_omicron[, c(2,1)]

df_omicron <- lineage_codification(df_genomas = df_omicron, dic_variants=dic_omicron)
rm(list = ls(pat = "^filtro"))

df_grouped <- group_variant_epiweek(df_genomas_recod = df_omicron, initial_date = min(df_omicron$date))
plot_omicron <- plot_prevalence_histogram(df_grouped = df_grouped,
                                          download_date = download_date,
                                          title = "Subvariants Omicron - Bogotá",
                                          pallete = "Cross")
ggsave(plot = plot_omicron, path_fig_hist_omicron_bog,
       width = 35, height = 20, units = "cm")


#### HISTOGRAMA PREVALENCIA BOGOTÁ (CUSTOM)####
path_fig_hist_bog_custom <- paste0('genomics/figures/prev_var_hist_bog_custom_', download_date_, '.png')

##FILTROS##
## Alpha ##
filter_alpha <- c("B.1.1.7", "^Q.")
## Gamma ##
filter_gamma <- c("P.1", "^P.1.")
## Delta ##
filter_delta <- c("B.1.617.2", "^AY.")
## Mu ##
filter_mu <- c("B.1.621", "B.1.621.1", "B.1.621.2", "BB.1", "BB.2")
## Omicron ##
filter_ba1 <- c("BA.1", "^BA.1.")
filter_ba2 <- c("BA.2", "^BA.2.")
# filter_ba3 <- c("BA.3", "^BA.3.")
filter_ba4 <- c("BA.4", "^BA.4.")
filter_ba5 <- c("BA.5", "^BA.5.")
# filter_b2_12 <- c("B.2.12", "^B.2.12.")

# filter_list = list(filter_alpha, filter_gamma, filter_delta, filter_mu, filter_ba1, filter_ba2, filter_ba4, filter_ba5)
# list_variants = list("Alpha", "Gamma", "Delta", "Mu", "BA.1", "BA.2", "BA.4", "BA.5")

filter_list = list(filter_gamma, filter_delta, filter_mu, filter_ba1, filter_ba2, filter_ba4, filter_ba5)
list_variants = list("Gamma", "Delta", "Mu", "BA.1", "BA.2", "BA.4", "BA.5")

list_lineages <- unique(df_genomes_bog$lineage)
dic_custom <- data.frame(label = character(),lineage = character())
for (i in 1:length(filter_list)){
  dic_custom <- add_variant_dictionary(dic_variants = dic_custom,
                                  list_lineages = list_lineages,
                                  variant=list_variants[[i]],
                                  filter_variant=filter_list[[i]])
}
dic_custom <- left_join(data.frame(lineage = list_lineages),
                         dic_custom, by = "lineage") %>%
  mutate(label = replace_na(label, "Other")) %>% arrange(label)
dic_custom <- dic_custom[, c(2,1)]

df_custom <- lineage_codification(df_genomas = df_genomes_bog, dic_variants=dic_custom)

df_grouped <- group_variant_epiweek(df_genomas_recod = df_custom, initial_date = initial_date_hist)
df_grouped <- df_grouped %>%
  mutate(lineage=factor(lineage, levels = c("Alpha", "Gamma", "Mu", "Other", "Delta", "BA.1", "BA.2", "BA.4", "BA.5"))) %>%
  arrange(date, lineage)
rm(filter_list, list_variants)
rm(list = ls(pat = "^filtro"))

plot_prev_custom <- plot_prevalence_histogram(df_grouped = df_grouped,
                                           download_date = download_date,
                                           title = "Prevalence per variant - Bogotá (custom)",
                                           pallete = "Thomas")

ggsave(plot = plot_prev_custom, path_fig_hist_bog_custom,
       width = 35, height = 20, units = "cm")
# write.csv(df_grouped, paste0("genomics/outputs/variants-ic-bog_custom_", download_date_,".csv"), row.names = FALSE)

#### HISTOGRAMA OTHER BOGOTÁ ####
path_fig_hist_other_bog <- paste0('genomics/figures/other_bog_', download_date_, '.png')

df_other <- df_genomes_bog[(df_genomes_bog$lineage_recod=="Other") & (df_genomes_bog$date >= "2022-03-27"), 
                           c('date', 'lineage')]
colnames(df_other) <- c("date", "lineage_recod")
unique(df_other$lineage_recod)
list_lineages <- unique(df_other$lineage_recod) 
df_grouped <- group_variant_epiweek(df_genomas_recod = df_other, initial_date = min(df_other$date))
plot_other <- plot_prevalence_histogram(df_grouped = df_grouped,
                                          download_date = download_date,
                                          title = "Other variants - Bogotá",
                                          pallete = "Cross")
ggsave(plot = plot_other, path_fig_hist_other_bog,
       width = 35, height = 20, units = "cm")
