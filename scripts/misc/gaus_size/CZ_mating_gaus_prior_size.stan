data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  real level;
  real<lower=0, upper=20> scale;
  real preference;
  real<lower=0, upper=20> choosiness;
  real asymmetry;
}

transformed parameters{
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = level + scale * exp(-0.5*((ratio[i]-preference)/choosiness)^2) + asymmetry * ratio[i];
}

model {
  level ~ normal(0, 10);
  preference ~ normal(0, 10);
  asymmetry ~ normal(0, 10);
  
  if (!posterior_predictive) {
    y ~ bernoulli_logit(y_hat);
  }
}


