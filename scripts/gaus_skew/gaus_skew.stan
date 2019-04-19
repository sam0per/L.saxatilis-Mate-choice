data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real<lower=0, upper=1> level;
  real<lower=0, upper=1> scale;
  real<lower=-10, upper=10> preference;
  real<lower=0, upper=10> choosiness;
  real<lower=-10, upper=10> asymmetry;
}

transformed parameters{
  //vector[N] err_val;
  vector[N] y_hat;

  for (i in 1:N) {
    y_hat[i] = level + scale * exp(-0.5 * ((ratio[i] - preference) / choosiness)^2) * (1 + erf(asymmetry * (ratio[i] - preference) / (1.414214 * choosiness)));
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
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | level + scale * exp(-0.5 * ((ratio[n] - preference) / choosiness)^2) * (1 + erf(asymmetry * (ratio[n] - preference) / (1.414214 * choosiness))));
  }
}
