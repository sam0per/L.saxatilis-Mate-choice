data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
  matrix[N, K-1] M;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  //real<lower=0, upper=0.05> level;           //a
  //real<lower=0, upper=1> scale;           //b
  //real<lower=-0.5, upper=0.5> preference;   //c
  //real<lower=0, upper=3> choosiness;     //d
  //real<lower=0, upper=5> asymmetry;    //g

  vector<lower=-3, upper=3>[K] b_coeff;
  vector<lower=-3, upper=3>[K] c_coeff;
  real<lower=0.01, upper=3> d_intercept;
  vector<lower=-4, upper=4>[K-1] d_coeff;
  vector<lower=-3, upper=3>[K] g_coeff;
  real<lower=0> b_sigma;
  real<lower=0> c_sigma;
  real<lower=0> d_sigma;
  real<lower=0> g_sigma;
}

transformed parameters {
  vector[N] scale_preds;
  vector<lower=0, upper=1>[N] scale;
  vector[N] preference;
  vector<lower=0>[N] choosiness;
  vector<lower=0>[N] asymmetry;
  vector[N] y_hat;

  scale_preds = X * b_coeff;
  scale = inv_logit(scale_preds);
  preference = X * c_coeff;
  choosiness = exp(d_intercept + M * d_coeff);
  asymmetry = exp(X * g_coeff);

  for (i in 1:N) {
    y_hat[i] = scale[i] * exp(-0.5 * ((ratio[i] - preference[i]) / choosiness[i])^2) * (1 + erf(asymmetry[i] * (ratio[i] - preference[i]) / (1.414214 * choosiness[i])));
  }
}

model {
  //scale ~ beta(6, 10);
  //preference ~ normal(-0.1, 0.1);
  //choosiness ~ gamma(7, 10);
  //asymmetry ~ gamma(5, 3);

  b_sigma ~ cauchy(0, 1.5);
  c_sigma ~ cauchy(0, 1.5);
  d_sigma ~ cauchy(0, 1.5);
  g_sigma ~ cauchy(0, 1.5);
  b_coeff ~ normal(0, b_sigma);
  c_coeff ~ normal(0, c_sigma);
  d_coeff ~ normal(0, d_sigma);
  g_coeff ~ normal(0, g_sigma);

  y ~ bernoulli_logit(logit(y_hat));
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | scale[n] * exp(-0.5 * ((ratio[n] - preference[n]) / choosiness[n])^2) * (1 + erf(asymmetry[n] * (ratio[n] - preference[n]) / (1.414214 * choosiness[n]))));
  }
  for (s in 1:N) {
    y_rep[s] = bernoulli_logit_rng(scale[s] * exp(-0.5 * ((ratio[s] - preference[s]) / choosiness[s])^2) * (1 + erf(asymmetry[s] * (ratio[s] - preference[s]) / (1.414214 * choosiness[s]))));
  }
}
