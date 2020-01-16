data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real<lower=0, upper=1> b0_par;
  real<lower=0, upper=1> b1_par;
  real<lower=-10, upper=10> c_par;
  real<lower=0, upper=10> d_par;
  real<lower=-10, upper=10> alpha_par;
}

transformed parameters {
  //vector[N] err_val;
  vector[N] y_hat;
  for (i in 1:N) {
    y_hat[i] = b0_par + b1_par * exp(-0.5 * ((ratio[i] - c_par) / d_par)^2) * (1 + erf(alpha_par * (ratio[i] - c_par) / (1.414214 * d_par)));
  }
  //y_prob = y_hat .* (1 + erf(err_val));
}

model {
  //level ~ normal(0, 10);
  //preference ~ normal(0, 10);
  //asymmetry ~ normal(0, 10);
  y ~ bernoulli_logit(logit(y_hat));
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | logit(b0_par + b1_par * exp(-0.5 * ((ratio[n] - c_par) / d_par)^2) * (1 + erf(alpha_par * (ratio[n] - c_par) / (1.414214 * d_par)))));
  }
  for (n in 1:N) {
    y_rep[n] = bernoulli_logit_rng(logit(b0_par + b1_par * exp(-0.5 * ((ratio[n] - c_par) / d_par)^2) * (1 + erf(alpha_par * (ratio[n] - c_par) / (1.414214 * d_par)))));
  }
}
