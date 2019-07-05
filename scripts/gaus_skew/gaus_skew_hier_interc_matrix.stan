data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real<lower=0, upper=1> b_intercept;
  vector[K] b_coeff;
  real<lower=-10, upper=10> c_par;
  real<lower=0.0001, upper=10> d_par;
  real<lower=-10, upper=10> alpha_par;
}

transformed parameters {
  vector[N] b_hyp;
  vector[N] y_hat;
  b_hyp = b_intercept + X * b_coeff;
  for (i in 1:N) {
    y_hat[i] = 0.01 + b_hyp[i] * exp(-0.5 * ((ratio[i] - c_par) / d_par)^2) * (1 + erf(alpha_par * (ratio[i] - c_par) / (1.414214 * d_par)));
  }
}

model {
  c_par ~ normal(-0.1, 0.5);
  d_par ~ gamma(7, 8);
  alpha_par ~ gamma(12, 7);
  b_coeff ~ normal(0, 1);
  b_hyp ~ beta(2, 2);
  y ~ bernoulli_logit(logit(y_hat));
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | logit(0.01 + b_hyp[n] * exp(-0.5 * ((ratio[n] - c_par) / d_par)^2) * (1 + erf(alpha_par * (ratio[n] - c_par) / (1.414214 * d_par)))));
  }
}
