// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
data {
    int N;
    vector[N] y;
}
parameters {
    real<lower=0> mu;
    real<lower=0> sigma;
}
model {
    mu ~ normal(0,1);
    sigma ~ normal(0,1);
    y ~ lognormal(mu, sigma);
}