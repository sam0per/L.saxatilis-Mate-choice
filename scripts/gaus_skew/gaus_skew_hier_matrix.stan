data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
  vector[N] ratio;
  int<lower=0, upper=1> y[N];
}

parameters {
  //real<lower=0, upper=0.05> level;           //a
  real<lower=0, upper=1> scale;           //b
  //real<lower=-0.5, upper=0.5> preference;   //c
  //real<lower=0, upper=3> choosiness;     //d
  //real<lower=0, upper=5> asymmetry;    //g
  
  vector<lower=-3, upper=3>[K] c_coeff;
  vector<lower=-3, upper=3>[K] d_coeff;
  vector<lower=-3, upper=3>[K] g_coeff;
  real<lower=0> c_sigma;
  real<lower=0> d_sigma;
  real<lower=0> g_sigma;
}

transformed parameters {
  //vector[N] level;
  //vector<lower=0, upper=1>[N] scale_preds;
  //vector[N] scale_preds;
  //vector<lower=0, upper=1>[N] scale;
  //vector<lower=-10, upper=10>[N] preference_preds;
  //vector<lower=-20, upper=20>[N] preference;
  vector[N] preference;
  //vector<lower=0, upper=10>[N] choosiness_preds;
  vector<lower=0>[N] choosiness;
  vector<lower=0>[N] asymmetry;
  //vector<lower=-10, upper=10>[N] asymmetry_preds;
  vector[N] y_hat;
  
  //scale_preds = X * b_coeff;
  //scale = inv_logit(scale_preds);
  preference = X * c_coeff;
  choosiness = exp(X * d_coeff);
  asymmetry = X * g_coeff;
  
  for (i in 1:N) {
    y_hat[i] = scale * exp(-0.5 * ((ratio[i] - preference[i]) / choosiness[i])^2) * (1 + erf(asymmetry[i] * (ratio[i] - preference[i]) / (1.414214 * choosiness[i])));
  }
}

model {
  scale ~ beta(6, 10);
  //preference ~ normal(-0.1, 0.1);
  //choosiness ~ gamma(7, 10);
  //asymmetry ~ gamma(5, 3);
  
  c_sigma ~ cauchy(0, 1);
  d_sigma ~ cauchy(0, 1);
  g_sigma ~ cauchy(0, 1);
  c_coeff ~ normal(0, c_sigma);
  d_coeff ~ normal(0, d_sigma);
  g_coeff ~ normal(0, g_sigma);
  
  y ~ bernoulli_logit(logit(y_hat));
}
