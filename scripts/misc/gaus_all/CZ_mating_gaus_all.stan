data {
  int<lower=0> N;
  vector[N] ratio;
  vector[N] shape;
  int<lower=1, upper=4> N_shore;
  int<lower=1, upper=2> N_eco;
  int<lower=1, upper=2> N_sex;
  int<lower=1, upper=N_shore> shore[N];
  int<lower=1, upper=N_eco> eco[N];
  int<lower=1, upper=N_sex> test[N];
  //int<lower=1, upper=N_sex> ref[N];
  int<lower=0, upper=1> y[N];
}

parameters {
  vector[N_shore] alpha1;
  //real alpha1_mu;
  real<lower=0> alpha1_sd;
  vector[N_eco] alpha2;
  //real alpha2_mu;
  real<lower=0> alpha2_sd;
  vector[N_sex] alpha3;
  //real alpha3_mu;
  real<lower=0> alpha3_sd;
  
  vector[N_shore] lambda1;
  //real lambda1_mu;
  real<lower=0> lambda1_sd;
  vector[N_eco] lambda2;
  //real lambda2_mu;
  real<lower=0> lambda2_sd;
  vector[N_sex] lambda3;
  //real lambda3_mu;
  real<lower=0> lambda3_sd;
  
  vector[N_shore] mu1;
  //real mu1_mu;
  real<lower=0> mu1_sd;
  vector[N_eco] mu2;
  //real mu2_mu;
  real<lower=0> mu2_sd;
  vector[N_sex] mu3;
  //real mu3_mu;
  real<lower=0> mu3_sd;
  
  vector[N_shore] sigma1;
  //real sigma1_mu;
  real<lower=0> sigma1_sd;
  vector[N_eco] sigma2;
  //real sigma2_mu;
  real<lower=0> sigma2_sd;
  vector[N_sex] sigma3;
  //real sigma3_mu;
  real<lower=0> sigma3_sd;
  
  vector[N_shore] gamma1;
  //real gamma1_mu;
  real<lower=0> gamma1_sd;
  vector[N_eco] gamma2;
  //real gamma2_mu;
  real<lower=0> gamma2_sd;
  vector[N_sex] gamma3;
  //real gamma3_mu;
  real<lower=0> gamma3_sd;
}

transformed parameters {
  vector[N] level;
  vector<lower=0>[N] scale;
  vector[N] preference;
  vector<lower=0>[N] choosiness;
  vector[N] asymmetry;
  vector[N] y_hat;
  for (i in 1:N) {
    level[i] = alpha1[shore[i]] + alpha2[eco[i]] + alpha3[test[i]] * shape[i];
    scale[i] = exp(lambda1[shore[i]] + lambda2[eco[i]] + lambda3[test[i]] * shape[i]);
    preference[i] = mu1[shore[i]] + mu2[eco[i]] + mu3[test[i]] * shape[i];
    choosiness[i] = exp(sigma1[shore[i]] + sigma2[eco[i]] + sigma3[test[i]] * shape[i]);
    asymmetry[i] = gamma1[shore[i]] + gamma2[eco[i]] + gamma3[test[i]] * shape[i];
  }
  for (i in 1:N)
    y_hat[i] = level[i] + scale[i] * exp(-0.5 * ((ratio[i] - preference[i]) / choosiness[i])^2) + asymmetry[i] * ratio[i];
}

model {
  alpha1 ~ normal(0, alpha1_sd);
  //alpha1_mu ~ normal(0, 1);
  alpha1_sd ~ cauchy(0, 25);
  alpha2 ~ normal(0, alpha2_sd);
  //alpha2_mu ~ normal(0, 1);
  alpha2_sd ~ cauchy(0, 25);
  alpha3 ~ normal(0, alpha3_sd);
  //alpha3_mu ~ normal(0, 1);
  alpha3_sd ~ cauchy(0, 25);
  
  lambda1 ~ normal(0, lambda1_sd);
  //lambda1_mu ~ normal(0, 1);
  lambda1_sd ~ cauchy(0, 25);
  lambda2 ~ normal(0, lambda2_sd);
  //lambda2_mu ~ normal(0, 1);
  lambda2_sd ~ cauchy(0, 25);
  lambda3 ~ normal(0, lambda3_sd);
  //lambda3_mu ~ normal(0, 1);
  lambda3_sd ~ cauchy(0, 25);
  
  mu1 ~ normal(0, mu1_sd);
  //mu1_mu ~ normal(0, 50);
  mu1_sd ~ cauchy(0, 25);
  mu2 ~ normal(0, mu2_sd);
  //mu2_mu ~ normal(0, 50);
  mu2_sd ~ cauchy(0, 25);
  mu3 ~ normal(0, mu3_sd);
  //mu3_mu ~ normal(0, 50);
  mu3_sd ~ cauchy(0, 25);
  
  sigma1 ~ normal(0, sigma1_sd);
  //sigma1_mu ~ normal(0, 50);
  sigma1_sd ~ cauchy(0, 25);
  sigma2 ~ normal(0, sigma2_sd);
  //sigma2_mu ~ normal(0, 50);
  sigma2_sd ~ cauchy(0, 25);
  sigma3 ~ normal(0, sigma3_sd);
  //sigma3_mu ~ normal(0, 50);
  sigma3_sd ~ cauchy(0, 25);
  
  gamma1 ~ normal(0, gamma1_sd);
  //gamma1_mu ~ normal(0, 50);
  gamma1_sd ~ cauchy(0, 25);
  gamma2 ~ normal(0, gamma2_sd);
  //gamma2_mu ~ normal(0, 50);
  gamma2_sd ~ cauchy(0, 25);
  gamma3 ~ normal(0, gamma3_sd);
  //gamma3_mu ~ normal(0, 50);
  gamma3_sd ~ cauchy(0, 25);
  
  y ~ bernoulli_logit(y_hat);
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n] = bernoulli_logit_rng(level[n] + scale[n] * exp(-0.5 * ((ratio[n] - preference[n]) / choosiness[n])^2) + asymmetry[n] * ratio[n]);
  }
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | level[n] + scale[n] * exp(-0.5 * ((ratio[n] - preference[n]) / choosiness[n])^2) + asymmetry[n] * ratio[n]);
  }
}
