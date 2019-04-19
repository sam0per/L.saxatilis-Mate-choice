data {
  int<lower=0> N;
  vector[N] ratio;
  vector[N] shape;
  int<lower=1> N_shore;
  int<lower=1> N_eco;
  int<lower=1, upper=N_shore> shore[N];
  int<lower=1, upper=N_eco> eco[N];
  int<lower=0, upper=1> y[N];
}

parameters {
  vector[N_shore] alpha1;
  //real alpha1_mu;
  real<lower=0> alpha1_sd;
  real alpha3;
  //real alpha3_mu;
  real<lower=0> alpha3_sd;
  
  vector[N_shore] lambda1;
  //real lambda1_mu;
  real<lower=0> lambda1_sd;
  real lambda3;
  //real lambda3_mu;
  real<lower=0> lambda3_sd;
  
  real preference;
  real<lower=0> choosiness;
  real asymmetry;
  
  //vector[N_shore] mu1;
  //real mu1_mu;
  //real<lower=0> mu1_sd;
  //real mu3;
  //real mu3_mu;
  //real<lower=0> mu3_sd;
  
  //vector[N_shore] sigma1;
  //real sigma1_mu;
  //real<lower=0> sigma1_sd;
  //real sigma3;
  //real sigma3_mu;
  //real<lower=0> sigma3_sd;
  
  //vector[N_shore] gamma1;
  //real gamma1_mu;
  //real<lower=0> gamma1_sd;
  //real gamma3;
  //real gamma3_mu;
  //real<lower=0> gamma3_sd;
}

transformed parameters {
  vector[N] level;
  vector<lower=0>[N] scale;
  //vector[N] preference;
  //vector<lower=0>[N] choosiness;
  //vector[N] asymmetry;
  real y_hat[N];
  level = alpha1[shore] + alpha3 * shape;
  scale = exp(lambda1[shore] + lambda3 * shape);
  //preference[i] = mu1[shore[i]] + mu3 * shape[i];
  //choosiness[i] = exp(sigma1[shore[i]] + sigma3 * shape[i]);
  //asymmetry[i] = gamma1[shore[i]] + gamma3 * shape[i];
    
  y_hat = level + scale * exp(-0.5 * ((ratio - preference) / choosiness)^2) + asymmetry * ratio;
}

model {
  alpha1 ~ normal(0, alpha1_sd);
  //alpha1_mu ~ normal(0, 1);
  alpha1_sd ~ cauchy(0, 20);
  alpha3 ~ normal(0, alpha3_sd);
  //alpha3_mu ~ normal(0, 1);
  alpha3_sd ~ cauchy(0, 20);
  
  lambda1 ~ normal(0, lambda1_sd);
  //lambda1_mu ~ normal(0, 1);
  lambda1_sd ~ cauchy(0, 20);
  lambda3 ~ normal(0, lambda3_sd);
  //lambda3_mu ~ normal(0, 1);
  lambda3_sd ~ cauchy(0, 20);
  
  //mu1 ~ normal(mu1_mu, mu1_sd);
  //mu1_mu ~ normal(0, 50);
  //mu1_sd ~ cauchy(0, 25);
  //mu3 ~ normal(mu3_mu, mu3_sd);
  //mu3_mu ~ normal(0, 50);
  //mu3_sd ~ cauchy(0, 25);
  
  //sigma1 ~ normal(sigma1_mu, sigma1_sd);
  //sigma1_mu ~ normal(0, 50);
  //sigma1_sd ~ cauchy(0, 25);
  //sigma3 ~ normal(sigma3_mu, sigma3_sd);
  //sigma3_mu ~ normal(0, 50);
  //sigma3_sd ~ cauchy(0, 25);
  
  //gamma1 ~ normal(gamma1_mu, gamma1_sd);
  //gamma1_mu ~ normal(0, 50);
  //gamma1_sd ~ cauchy(0, 25);
  //gamma3 ~ normal(gamma3_mu, gamma3_sd);
  //gamma3_mu ~ normal(0, 50);
  //gamma3_sd ~ cauchy(0, 25);
  
  //level ~ normal(pred_level, 10);
  //scale ~ lognormal(pred_scale, 0.5);
  preference ~ normal(0, 10);
  asymmetry ~ normal(0, 10);
  
  y ~ bernoulli_logit(y_hat);
}
