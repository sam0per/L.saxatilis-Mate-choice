data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real level;
  real<lower=0, upper=10> scale;
  real preference;
  real<lower=0, upper=10> choosiness;
  real asymmetry;
}

transformed parameters{
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = level + scale * exp(-0.5*((ratio[i]-preference)/choosiness)^2) + asymmetry * ratio[i];
}

model {
  level ~ normal(0, 5);
  preference ~ normal(0, 5);
  asymmetry ~ normal(0, 5);
  
  y ~ bernoulli_logit(y_hat);
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n] = bernoulli_logit_rng(level + scale * exp(-0.5*((ratio[n]-preference)/choosiness)^2) + asymmetry * ratio[n]);
  }
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | level + scale * exp(-0.5*((ratio[n]-preference)/choosiness)^2) + asymmetry * ratio[n]);
    }
}
