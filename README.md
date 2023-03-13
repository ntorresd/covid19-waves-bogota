# Four waves of SARS-CoV-2 in Bogotá: a detailed retrospective statistical comparison

In this repository we characterized and compared, using statistical tools, the first four consecutive epidemiological waves in the Bogotá, Colombia, that occurred between March 2020 and April 2022. We used the report of confirmed cases from the [District Health Secretary of Bogotá](https://datosabiertos.bogota.gov.co/dataset/numero-de-casos-confirmados-por-el-laboratorio-de-covid-19-bogota-d-c) , and the genomic surveillance data published by the [Global Initiative on Sharing All Influenza Data (GISAID)](https://gisaid.org/). We focused mainly on the estimation of:

1. The instantaneous reproduction number R(t).
2. The transmissibility advantage between variants.
3. The delay times for onset-to-hospitalization, onset-to-ICU, onset-to-death, hospital stay, and ICU stay. 
4. The characterization of severe outcomes by means of the proportion of cases using the general hospital and ICU services, as well, as deaths, and the severe ratios: Case Fatality Rate (CFR), Hospitalization Case Rate (HCR), and Hospitalization Fatality Rate (HFR) for each wave. 

## Data sources

1. Report of confirmed cases from the [District Health Secretary of Bogotá](https://datosabiertos.bogota.gov.co/dataset/numero-de-casos-confirmados-por-el-laboratorio-de-covid-19-bogota-d-c) (Private database - last update: 2022-08-02)

2. Genomic surveillance data published by the [Global Initiative on Sharing All Influenza Data (GISAID)](https://gisaid.org/) (Public database - last update: 2022-08-02)

## Methods

1. Reproduction number R(t): we estimated the time-varying instantaneous reproduction number R(t) using the epidemiological package for R: [EpiEstim](https://cran.r-project.org/web/packages/EpiEstim/index.html)

2. Transmissibility advantage: we evaluated the transmissibility advantage using a multinomial logistic regression with a single explanatory variable $t$ given by:

$$f(v,t)=\alpha + \beta_{v,0}t$$

In the previous expression $\alpha$ is the intercept of the model and $\beta_{v,0}$ is the variant-specific parameter for the time covariate, which was computed with respect to a reference (or pivot) variant. In order to compute the transmissibility advantage of $v$ with respect to $0$ from $\beta_{v,0}$ we used the following expression:

$$ T_{v,0} = \exp(\frac{\beta_{v,0}}{7} * g_0) $$

Where $g_0$ is the generation time of the pivot variant. Notice that we divided $\beta_{v,0}$ by 7 to convert time scale of the coefficients to daily.

With these coefficientes we computed the relative transmissibiliy advantage between two variants $w$ and $v$ as:

$$ T_{w,v} = \frac{T_{w,0}}{T_{v,0}}$$

The multinomial regressions were run in stan using the library [PyStan](https://pystan.readthedocs.io/en/latest/) for python.

3. Probability distributions of delay times: we used a bayesian hierarchical model adapted from [this repository](https://github.com/mrc-ide/Brazil_COVID19_distributions). We fitted initial parameters for the district level and then sample the parameters for each wave as follows:

$$q_{i,j} \sim N(q_{i, Bog},\sigma_i)$$

where $i=1,2,3, .., n$ runs over the $n$ parameters of the PDF, $j=1,2,3,4$ is the number of wave, $q_{i, Bog}$ is the value of the i-th parameter of the PDF estimated for Bogotá and $\sigma_i \sim N^+(0,1)$ is the standard deviation which is assumed to be distributed as a truncated normal distribution. We used a Hamiltonian Monte Carlo (HMC) algorithm implemented in Stan, setting four chains of 2000 iterations (1000 for warming up and 1000 for sampling). These models were also run using [PyStan](https://pystan.readthedocs.io/en/latest/).

4. Severe outcomes: we computed the proportion of severe outcomes (general hospitalization, ICU, and death) between the Covid-19 cases, with a confidence interval of 95% for binomial proportions; and the Case Fatality Rate (CFR), Hospitalization Case Rate (HCR), and Hospitalization Fatality Rate (HFR) for each wave, disaggregated by 10 years age-groups. Additionally, we calculated the contribution of waves to severe outcomes by age-groups.

## Repository description

All the folders, except [plots](/plots/) and [tables](/tables/), are organized following the same structure: a **scripts** subfolder that contains the necessary codes to run the models and the analysis, and an **outputs** subfolder that contains the results. 

1. The folder [epidemiological_distributions](/epidemiological_distributions/) contains the following scripts:
* Models implemented in stan for the distric level and the partial pooling for waves: [model_exponential_district.stan](/epidemiological_distributions/scripts/model_exponential_district.stan), [model_gamma_district.stan](/epidemiological_distributions/scripts/model_gamma_district.stan), [model_gln_district.stan](/epidemiological_distributions/scripts/model_gln_district.stan), [model_log_normal_district.stan](/epidemiological_distributions/scripts/model_log_normal_district.stan), [model_weibull_district.stan](/epidemiological_distributions/scripts/model_exponential_district.stan), [model_exponential_pool.stan](/epidemiological_distributions/scripts/model_exponential_pool.stan), [model_gamma_pool.stan](/epidemiological_distributions/scripts/model_gamma_pool.stan), [model_gln_pool.stan](/epidemiological_distributions/scripts/model_gln_pool.stan), [model_log_normal_pool.stan](/epidemiological_distributions/scripts/model_log_normal_pool.stan), [model_weibull_pool.stan](/epidemiological_distributions/scripts/model_exponential_pool.stan).

* Python scripts to run the models and extract the results: [run_exponential.py](/epidemiological_distributions/scripts/run_exponential.py), [run_gamma.py](/epidemiological_distributions/scripts/run_gamma.py), [run_gln.py](/epidemiological_distributions/scripts/run_gln.py), [run_log_normal.py](/epidemiological_distributions/scripts/run_log_normal.py), [run_weibull.py](/epidemiological_distributions/scripts/run_weibull.py)).

* [bayes_factor.py](/epidemiological_distributions/scripts/bayes_factors.py): python script to calculate the bayes factor.
* [summarize_results.py](/epidemiological_distributions/scripts/summarize_results.py): python script that sumarizes the main results of the bayesian inference.
* [utilities_epi_dist.py](/epidemiological_distributions/scripts/utilities_epi_dist.py): python script with functions used for preparing and cleaning the data, and tools for the statistical analysis .

2. The folder [genomics](/genomics/) contains the following scripts:
* [genomics_functions.R](/genomics/scripts/genomics_functions.R), [variants_record](/genomics/scripts/variants_recod.R): functions implemented in R to process the genomic data.
* [process_data.py](/genomics/scripts/process_data.py): python script to generate the inputs for the model.
* [multinomial_model.stan](/genomics/scripts/multinomial_model.stan): multinomial model implemented in stan. 
* [run_model.py](/genomics/scripts/run_model.py): python script to run the multinomial model.
* [process_results](/genomics/scripts/process_results.py): python script to process the results from the multinomial model.

3. The folder [rt](/rt/) contains the following script:
* [rt.R](/rt/scripts/rt.R): R script to estimate the Reproduction number. 

4. The folder [severe_outcomes](/severe_outcomes/) contains the following scripts:
* [percentages.py](/severe_outcomes/scripts/percentages.py): python script to calculate the percentages.
* [proportions.py](/severe_outcomes/scripts/proportions.py): python script to calculate the binomial proportions.
* [rates.py](/severe_outcomes/scripts/rates.py): python script to calculate the CFR, HCR, ICU-CR, HFR and ICU-FR.
* [utilities_severity.py](/severe_outcomes/scripts/utilities_severity.py): python script with tools used for the severity analysis.

5. The folder [waves](/waves/) contains the following scripts:
* [roots_confirmed_cases.py](/waves/scripts/roots_confirmed_cases.py): python script to find the roots of the epidemic curve using gaussian smoothing and interpolation.
* [process_waves.py](/waves/scripts/process_waves.py): python script to process the waves after visual inspection of the roots.
* [utilities_waves.py](/waves/scripts/utilities_waves.py): python script with tools used for determining the waves.

6. The folder [tables](/tables/) contains the following scripts:
* Python scripts to generate the tables of the supplementary materials: [table_s1.py](/tables/table_s1.py), [table_s2.py](/tables/table_s2.py), [table_s4.py](/tables/table_s4.py), [table_s5.py](/tables/table_s5.py). 

This scripts import and call functions from the utilities of each section and the results contained in the corresponding output subfolder.

7. The folder [plots](/tables/) contains the following scripts:
* [plot_style.mplstyle](/plots/plot_style.mplstyle): script with the styles used to generate all the figures.
* Python scripts to generate the figures included in the main document: [figure_1.py](/plots/figure_1.py), [figure_2.py](/plots/figure_2.py), [figure_3.py](/plots/figure_3.py), [figure_4.py](/plots/figure_4.py), [figure_5.py](/plots/figure_5.py)).

* [individual_plots.py](/plots/individual_plots.py): python script to generate the individual plots included in the supplementary materials and some visualizations of the analysis.
* Scripts with the visualization functions for every section of the analysis: [overview.py](/plots/overview.py), [results_epidemiological_distributions.py](/plots/results_epidemiological_distributions.py), [results_genomics.py](/plots/results_genomics.py), [results_rt.py](/plots/results_rt.py), [results_severe_outcomes.py](/plots/results_severe_outcomes.py), [results_waves.py](/plots/results_waves.py).

### Config file
This [YAML file](/config.yml) contains information about the models, the paths used in the scripts and the roots selected for the waves. It is called in the beginning of all the scripts.

## Authors and contributors
[davidsantiagoquevedo](https://github.com/davidsantiagoquevedo/), [ntorresd](https://github.com/ntorresd/), [cwhittaker1000](https://github.com/cwhittaker1000/), [zmcucunuba](https://github.com/zmcucunuba).