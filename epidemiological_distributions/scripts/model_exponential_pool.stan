// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
data {
    int K; // stratification (age, sex, waves)
    int N; // total number of observations
    real X[N]; // observations
    int wave[N]; // index with the strat number for each observation
    real beta_prior; // prior value of beta from district-level model
}
parameters {
    real<lower=0> beta[K];
    // hyperparameters
    real<lower=0> sigma_beta;
}
model {
    // likelihood
    for (i in 1:N){
            X[i] ~ exponential(beta[wave[i]]);
    }
    // priors
    beta ~ normal(beta_prior, sigma_beta);
    // hyperpriors
    sigma_beta ~ normal(0,1);
}

