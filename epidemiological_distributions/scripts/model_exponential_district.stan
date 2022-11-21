// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
data {
    int N;
    vector[N] y;
}
parameters {
    real<lower=0> beta;
}
model {
    beta ~ normal(1,1);
    y ~ exponential(beta);
}

