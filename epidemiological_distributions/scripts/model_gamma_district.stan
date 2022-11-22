// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
data {
    int N;
    vector[N] y;
}
parameters {
    real<lower=0> alpha;
    real<lower=0> beta;
}
model {
    alpha ~ normal(0,1);
    beta ~ normal(0,1);
    y ~ gamma(alpha, beta);
}