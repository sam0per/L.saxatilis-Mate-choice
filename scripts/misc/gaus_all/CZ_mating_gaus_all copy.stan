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
  vector<lower=0, upper=20>[N] scale;
  vector[N] preference;
  vector<lower=0, upper=20>[N] choosiness;
  vector[N] asymmetry;
  vector[N] y_hat;
  for (i in 1:N) {
    level[i] = alpha1[shore[i]] + alpha2[eco[i]] + alpha3 * shape[i];
    scale[i] = exp(lambda1[shore[i]] + lambda2[eco[i]] + lambda3 * shape[i]);
    preference[i] = mu1[shore[i]] + mu2[eco[i]] + mu3 * shape[i];
    choosiness[i] = exp(sigma1[shore[i]] + sigma2[eco[i]] + sigma3 * shape[i]);
    asymmetry[i] = gamma1[shore[i]] + gamma2[eco[i]] + gamma3 * shape[i];
  }
  
  for (i in 1:N)
    y_hat[i] = level[i] + scale[i] * exp(-0.5*((ratio[i]-preference[i])/choosiness[i])^2) + asymmetry[i] * ratio[i];
}

model {
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
  
  y ~ bernoulli_logit(y_hat);
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = bernoulli_logit_rng(level[n] + (1 / sqrt(2 * pi() * stan_dev[n]^2)) * 
    exp(-0.5*((ratio[n]-aver[n])/stan_dev[n])^2) + rtail[n]);
}
