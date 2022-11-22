// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
data {
    int K; // stratification (age, sex, waves)
    int N; // total number of observations
    real X[N]; // observations
    int wave[N]; // index with the strat number for each observation
    real mu_prior; // prior value of mu from district-level model
    real sigma_prior; // prior value of sigma from district-level model
}
parameters {
    real<lower=0> mu[K];
    real<lower=0> sigma[K];
    // hyperparameters
    real<lower=0> sigma_mu;
    real<lower=0> sigma_sigma;
}
model {
    // likelihood
    for (i in 1:N){
            X[i] ~ lognormal(mu[wave[i]], sigma[wave[i]]);
    }
    // priors
    mu ~ normal(mu_prior, sigma_mu);
    sigma ~ normal(sigma_prior, sigma_sigma);
    // hyperpriors
    sigma_mu ~ normal(0,1);
    sigma_sigma ~ normal(0,1);
}