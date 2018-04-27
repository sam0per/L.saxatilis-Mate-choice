data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real level;
  real aver;
  real gamma_ta;
  real<lower=0, upper=30> stan_dev;
  real<lower=0, upper=30> sigma_level;
  real<lower=0, upper=30> sigma_aver;
  real<lower=0, upper=30> sigma_gamma_ta;
}

model {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = level + (1 / sqrt(2 * pi() * stan_dev^2)) * exp(-0.5*((ratio[i]-aver)/stan_dev)^2) + gamma_ta * ratio[i];
  
  level ~ normal(0,sigma_level);
  aver ~ normal(0,sigma_aver);
  gamma_ta ~ normal(0,sigma_gamma_ta);

  sigma_level ~ cauchy(0, 10);
  sigma_aver ~ cauchy(0, 10);
  sigma_gamma_ta ~ cauchy(0, 10);
  
  y ~ bernoulli_logit(y_hat);
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = bernoulli_logit_rng(level + (1 / sqrt(2 * pi() * stan_dev^2)) * 
    exp(-0.5*((ratio[n]-aver)/stan_dev)^2) + gamma_ta * ratio[n]);
}
