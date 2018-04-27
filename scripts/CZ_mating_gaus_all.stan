data {
  int<lower=0> N;
  vector[N] ratio;
  vector[N] shape;
  int<lower=0> N_shore;
  int<lower=0> N_eco;
  int<lower=1, upper=N_shore> shore[N];
  int<lower=1, upper=N_eco> eco[N];
  int<lower=0, upper=1> y[N];
}

parameters {
  vector[N_shore] alpha;
  vector[N_eco] delta;
  real beta;
  vector[N_shore] alpha_sd;
  vector[N_eco] delta_sd;
  real beta_sd;
  vector[N_shore] alpha_av;
  vector[N_eco] delta_av;
  real beta_av;
  vector[N_shore] alpha_ta;
  vector[N_eco] delta_ta;
  real gamma_ta;
  real beta_ta;
  real<lower=0, upper=30> sigma_alpha;
  real<lower=0, upper=30> sigma_alpha_sd;
  real<lower=0, upper=30> sigma_alpha_av;
  real<lower=0, upper=30> sigma_alpha_ta;
  real<lower=0, upper=30> sigma_delta;
  real<lower=0, upper=30> sigma_delta_sd;
  real<lower=0, upper=30> sigma_delta_av;
  real<lower=0, upper=30> sigma_delta_ta;
  real<lower=0, upper=30> sigma_gamma_ta;
  real<lower=0, upper=30> sigma_beta;
  real<lower=0, upper=30> sigma_beta_sd;
  real<lower=0, upper=30> sigma_beta_av;
  real<lower=0, upper=30> sigma_beta_ta;
}

transformed parameters {
  vector[N] level;
  vector[N] aver;
  vector<lower=0, upper=30>[N] stan_dev;
  vector[N] rtail;
  for (i in 1:N) {
    level[i] = alpha[shore[i]] + delta[eco[i]] + beta * shape[i];
    stan_dev[i] = exp(alpha_sd[shore[i]] + delta_sd[eco[i]] + beta_sd * shape[i]);
    aver[i] = alpha_av[shore[i]] + delta_av[eco[i]] + beta_av * shape[i];
    rtail[i] = alpha_ta[shore[i]] + delta_ta[eco[i]] + gamma_ta * ratio[i] + beta_ta * shape[i];
  }
}

model {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = level[i] + (1 / sqrt(2 * pi() * stan_dev[i]^2)) * exp(-0.5*((ratio[i]-aver[i])/stan_dev[i])^2) + rtail[i];
  y ~ bernoulli_logit(y_hat);

  alpha ~ normal(0,sigma_alpha);
  alpha_sd ~ normal(0,sigma_alpha_sd);
  alpha_av ~ normal(0,sigma_alpha_av);
  alpha_ta ~ normal(0,sigma_alpha_ta);
  
  sigma_alpha ~ cauchy(0, 10);
  sigma_alpha_sd ~ cauchy(0, 10);
  sigma_alpha_av ~ cauchy(0, 10);
  sigma_alpha_ta ~ cauchy(0, 10);
  
  delta ~ normal(0,sigma_delta);
  delta_sd ~ normal(0,sigma_delta_sd);
  delta_av ~ normal(0,sigma_delta_av);
  delta_ta ~ normal(0,sigma_delta_ta);

  sigma_delta ~ cauchy(0, 10);
  sigma_delta_sd ~ cauchy(0, 10);
  sigma_delta_av ~ cauchy(0, 10);
  sigma_delta_ta ~ cauchy(0, 10);

  gamma_ta ~ normal(0,sigma_gamma_ta);
  beta ~ normal(0,sigma_beta);
  beta_sd ~ normal(0,sigma_beta_sd);
  beta_av ~ normal(0,sigma_beta_av);
  beta_ta ~ normal(0,sigma_beta_ta);

  sigma_gamma_ta ~ cauchy(0, 10);
  sigma_beta ~ cauchy(0, 10);
  sigma_beta_sd ~ cauchy(0, 10);
  sigma_beta_av ~ cauchy(0, 10);
  sigma_beta_ta ~ cauchy(0, 10);
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = bernoulli_logit_rng(level[n] + (1 / sqrt(2 * pi() * stan_dev[n]^2)) * 
    exp(-0.5*((ratio[n]-aver[n])/stan_dev[n])^2) + rtail[n]);
}
