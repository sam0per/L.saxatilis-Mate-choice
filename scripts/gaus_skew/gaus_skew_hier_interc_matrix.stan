data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
  //vector[N] X;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real<lower=-4, upper=4> alpha_intercept;
  real<lower=0, upper=1> b_intercept;
  real<lower=-2, upper=2> c_intercept;
  real<lower=0.0001, upper=3> d_intercept;
  vector[K] alpha_coeff;
  vector[K] b_coeff;
  vector[K] c_coeff;
  vector[K] d_coeff;
  //real d_coeff;
  //real<lower=0, upper=1> b_par;
  //real<lower=-10, upper=10> c_par;
  //real<lower=0.0001, upper=10> d_par;
  //real<lower=-10, upper=10> alpha_par;
}

transformed parameters {
  vector[N] alpha_hyp;
  vector[N] b_hyp;
  vector[N] c_hyp;
  vector[N] d_hyp;
  vector[N] y_hat;
  alpha_hyp = alpha_intercept + X * alpha_coeff;
  b_hyp = b_intercept + X * b_coeff;
  c_hyp = c_intercept + X * c_coeff;
  d_hyp = d_intercept + X * d_coeff;
  for (i in 1:N) {
    y_hat[i] = 0.01 + b_hyp[i] * exp(-0.5 * ((ratio[i] - c_hyp[i]) / d_hyp[i])^2) * (1 + erf(alpha_hyp[i] * (ratio[i] - c_hyp[i]) / (1.414214 * d_hyp[i])));
  }
}

model {
  alpha_coeff ~ normal(0, 1);
  b_coeff ~ normal(0, 1);
  c_coeff ~ normal(0, 1);
  d_coeff ~ normal(0, 1);
  alpha_hyp ~ gamma(12, 7);
  b_hyp ~ beta(2, 2);
  c_hyp ~ normal(-0.1, 0.5);
  d_hyp ~ gamma(7, 8);
  y ~ bernoulli_logit(logit(y_hat));
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | logit(0.01 + b_hyp[n] * exp(-0.5 * ((ratio[n] - c_hyp[n]) / d_hyp[n])^2) * (1 + erf(alpha_hyp[n] * (ratio[n] - c_hyp[n]) / (1.414214 * d_hyp[n])))));
  }
}
