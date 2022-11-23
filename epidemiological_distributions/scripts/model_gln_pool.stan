// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
functions{
  real custom_lpdf(real x, real mu, real sigma, real g)
    {
      real logK = log(g) - (g+1)/g*log(2)-log(sigma)-lgamma(1/g);
      real tmp = logK - log(x) - 0.5 * pow(fabs((log(x)-mu)/sigma),g);
      return tmp;
    }
}
data {
    int K; // stratification (age, sex, waves)
    int N; // total number of observations
    real X[N]; // observations
    int wave[N]; // index with the strat number for each observation
    real mu_prior; // prior value of mu from district-level model
    real sigma_prior; // prior value of sigma from district-level model
    real g_prior; // prior value of g from district-level model
}
parameters {
    real<lower=0> mu[K];
    real<lower=0> sigma[K];
    real<lower=1> g[K];
    // hyperparameters
    real<lower=0> sigma_mu;
    real<lower=0> sigma_sigma;
    real<lower=0> sigma_g;
}
model {
    // likelihood
    for (i in 1:N){
            X[i] ~ custom(mu[wave[i]], sigma[wave[i]], g[wave[i]]);
    }
    // priors
    mu ~ normal(mu_prior, sigma_mu);
    sigma ~ normal(sigma_prior, sigma_sigma);
    g ~ normal(g_prior, sigma_g);
    // hyperpriors
    sigma_mu ~ normal(2, 0.5);
    sigma_sigma ~ normal(0.5, 0.5);
    sigma_g ~ normal(1.5, 0.5);
}