data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real level;
  real aver;
  real gamma_ta;
  real<lower=0, upper=50> stan_dev;
}

model {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = level + (1 / sqrt(2 * pi() * stan_dev^2)) * exp(-0.5*((ratio[i]-aver)/stan_dev)^2) + gamma_ta * ratio[i];
  
  y ~ bernoulli_logit(y_hat);
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = bernoulli_logit_rng(level + (1 / sqrt(2 * pi() * stan_dev^2)) * exp(-0.5*((ratio[n]-aver)/stan_dev)^2) + gamma_ta * ratio[n]);
}
