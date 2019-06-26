data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  vector<lower=-3, upper=3>[K] b_coeff;
  vector<lower=-3, upper=3>[K] c_coeff;
  vector<lower=-3, upper=3>[K] d_coeff;
  vector<lower=-3, upper=3>[K] g_coeff;
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
  choosiness = exp(X * d_coeff);
  asymmetry = exp(X * g_coeff);

  for (i in 1:N) {
    y_hat[i] = scale[i] * exp(-0.5 * ((ratio[i] - preference[i]) / choosiness[i])^2) * (1 + erf(asymmetry[i] * (ratio[i] - preference[i]) / (1.414214 * choosiness[i])));
  }
}

model {
  b_coeff ~ normal(0, 2);
  c_coeff ~ normal(0, 2);
  d_coeff ~ normal(0, 2);
  g_coeff ~ normal(0, 2);

  y ~ bernoulli_logit(logit(y_hat));
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | logit(scale[n] * exp(-0.5 * ((ratio[n] - preference[n]) / choosiness[n])^2) * (1 + erf(asymmetry[n] * (ratio[n] - preference[n]) / (1.414214 * choosiness[n])))));
  }
  for (s in 1:N) {
    y_rep[s] = bernoulli_logit_rng(logit(scale[s] * exp(-0.5 * ((ratio[s] - preference[s]) / choosiness[s])^2) * (1 + erf(asymmetry[s] * (ratio[s] - preference[s]) / (1.414214 * choosiness[s])))));
  }
}
