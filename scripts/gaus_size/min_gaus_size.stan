data {
  int<lower=0> N;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  //real level;
  real<lower=0, upper=5> scale;
  real preference;
  real<lower=0, upper=5> choosiness;
  real asymmetry;
}

transformed parameters{
  real y_hat[N];
  for (i in 1:N)
    y_hat[i] = scale * exp(-0.5*((ratio[i]-preference)/choosiness)^2) + asymmetry * ratio[i];
}

model {
  //level ~ normal(0, 10);
  preference ~ normal(0, 1);
  asymmetry ~ normal(0, 1);
  
  y ~ bernoulli_logit(y_hat);
}
