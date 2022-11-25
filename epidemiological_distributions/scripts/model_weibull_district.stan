// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
data {
    int N;
    vector[N] y;
}
parameters {
    real<lower=0> alpha;
    real<lower=0> sigma;
}
model {
    alpha ~ normal(1,1);
    sigma ~ normal(1,0.5);
    y ~ weibull(alpha, sigma);
}