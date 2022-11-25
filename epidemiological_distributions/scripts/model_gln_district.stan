// Created on Sun Nov 20 2022
// @author: dsquevedo
// @author: ntorres
functions{
 real custom_lpdf(real x, real mu, real sigma, real g)
    {
     real logK = log(g) - (g+1)/g*log(2) - log(sigma) - lgamma(1/g);
     real tmp = logK - log(x) - 0.5 * pow(fabs((log(x)- mu)/sigma),g);
     return tmp;
    }
}
data {
    int N;
    real y[N];
}
parameters {
    real<lower=0> mu;
    real<lower=0> sigma;
    real<lower=1> g;
}
model {
      for (i in 1:N){
        y[i] ~ custom(mu, sigma, g);
        }
      mu ~ normal(2, 0.5);
      sigma ~ normal(0.5, 0.5);
      g ~ normal(1.5, 0.5);
}