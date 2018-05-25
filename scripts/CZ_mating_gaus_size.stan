data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real level;
  real<lower=0, upper=10> lambda;
  real aver;
  real gamma_ta;
  real<lower=0, upper=10> stan_dev;
}

transformed parameters{
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = level + lambda * exp(-0.5*((ratio[i]-aver)/stan_dev)^2) + gamma_ta * ratio[i];
}

model {
  level ~ normal(0, 10);
  aver ~ normal(0, 10);
  gamma_ta ~ normal(0, 10);
  
  y ~ bernoulli_logit(y_hat);
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = bernoulli_logit_rng(level + lambda * exp(-0.5*((ratio[n]-aver)/stan_dev)^2) + gamma_ta * ratio[n]);
}
