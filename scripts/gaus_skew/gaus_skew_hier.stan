data {
  int<lower=0> N;
  vector[N] ratio;
  vector[N] shape;
  int<lower=1, upper=4> N_shore;
  int<lower=1, upper=4> N_eco;
  int<lower=1, upper=2> N_sex;
  int<lower=1, upper=N_shore> shore[N];
  int<lower=1, upper=N_eco> ref[N];
  int<lower=1, upper=N_sex> test[N];
  int<lower=0, upper=1> y[N];
}

parameters {
  real<lower=0, upper=1> level;
  //real<lower=0, upper=1> scale;
  //real<lower=-10, upper=10> preference;
  real<lower=0, upper=10> choosiness;
  real<lower=-10, upper=10> asymmetry;
  
  vector[N_shore] lambda1;
  //real lambda1_mu;
  real<lower=0> lambda1_sd;
  vector[N_eco] lambda2;
  //real lambda2_mu;
  real<lower=0> lambda2_sd;
  vector[N_sex] lambda3;
  //real lambda3_mu;
  real<lower=0> lambda3_sd;
  
  vector<lower=-5, upper=5>[N_shore] mu1;
  //real mu1_mu;
  real<lower=0> mu1_sd;
  vector<lower=-5, upper=5>[N_eco] mu2;
  //real mu2_mu;
  real<lower=0> mu2_sd;
  vector<lower=-5, upper=5>[N_sex] mu3;
  //real mu3_mu;
  real<lower=0> mu3_sd;
  
  //vector[N_shore] sigma1;
  //real sigma1_mu;
  //real<lower=0> sigma1_sd;
  //vector[N_eco] sigma2;
  //real sigma2_mu;
  //real<lower=0> sigma2_sd;
  //vector[N_sex] sigma3;
  //real sigma3_mu;
  //real<lower=0> sigma3_sd;
  
  //vector[N_shore] gamma1;
  //real gamma1_mu;
  //real<lower=0> gamma1_sd;
  //vector[N_eco] gamma2;
  //real gamma2_mu;
  //real<lower=0> gamma2_sd;
  //vector[N_sex] gamma3;
  //real gamma3_mu;
  //real<lower=0> gamma3_sd;
}

transformed parameters {
  //vector[N] level;
  //vector<lower=0, upper=1>[N] scale_preds;
  vector[N] scale_logit;
  vector<lower=0, upper=1>[N] scale;
  //vector<lower=-10, upper=10>[N] preference_preds;
  vector<lower=-20, upper=20>[N] preference;
  //vector<lower=0, upper=10>[N] choosiness_preds;
  //vector[N] choosiness;
  //vector<lower=-10, upper=10>[N] asymmetry_preds;
  vector[N] y_hat;
  for (i in 1:N) {
    //level[i] = alpha1[shore[i]] + alpha2[ref[i]] + alpha3[test[i]] * shape[i];
    scale_logit[i] = lambda1[shore[i]] + lambda2[ref[i]] + lambda3[test[i]] * shape[i];
    scale[i] = inv_logit(scale_logit[i]);
    //preference_preds[i] = mu1[shore[i]] + mu2[ref[i]] + mu3[test[i]] * shape[i];
    preference[i] = mu1[shore[i]] + mu2[ref[i]] + mu3[test[i]] * shape[i];
    //choosiness_preds[i] = exp(sigma1[shore[i]] + sigma2[ref[i]] + sigma3[test[i]] * shape[i]);
    //choosiness[i] = exp(sigma1[shore[i]] + sigma2[ref[i]] + sigma3[test[i]] * shape[i]);
    //asymmetry_preds[i] = gamma1[shore[i]] + gamma2[ref[i]] + gamma3[test[i]] * shape[i];

    y_hat[i] = level + scale[i] * exp(-0.5 * ((ratio[i] - preference[i]) / choosiness)^2) * (1 + erf(asymmetry * (ratio[i] - preference[i]) / (1.414214 * choosiness)));
  }
}

model {
  //alpha1 ~ normal(0, alpha1_sd);
  //alpha1_mu ~ normal(0, 1);
  //alpha1_sd ~ cauchy(0, 10);
  //alpha2 ~ normal(0, alpha2_sd);
  //alpha2_mu ~ normal(0, 1);
  //alpha2_sd ~ cauchy(0, 10);
  //alpha3 ~ normal(0, alpha3_sd);
  //alpha3_mu ~ normal(0, 1);
  //alpha3_sd ~ cauchy(0, 10);
  
  lambda1 ~ normal(0, lambda1_sd);
  //lambda1_mu ~ normal(0, 1);
  lambda1_sd ~ normal(0, 5);
  lambda2 ~ normal(0, lambda2_sd);
  //lambda2_mu ~ normal(0, 1);
  lambda2_sd ~ normal(0, 5);
  lambda3 ~ normal(0, lambda3_sd);
  //lambda3_mu ~ normal(0, 1);
  lambda3_sd ~ normal(0, 5);
  
  //mu1 ~ normal(0, mu1_sd);
  //mu1_mu ~ normal(0, 50);
  //mu1_sd ~ normal(0, 2);
  //mu2 ~ normal(0, mu2_sd);
  //mu2_mu ~ normal(0, 50);
  //mu2_sd ~ normal(0, 2);
  //mu3 ~ normal(0, mu3_sd);
  //mu3_mu ~ normal(0, 50);
  //mu3_sd ~ normal(0, 2);
  
  //sigma1 ~ normal(0, sigma1_sd);
  //sigma1 ~ normal(0, 5);
  //sigma1_mu ~ normal(0, 50);
  //sigma1_sd ~ normal(0, 5);
  //sigma2 ~ normal(0, sigma2_sd);
  //sigma2 ~ normal(0, 5);
  //sigma2_mu ~ normal(0, 50);
  //sigma2_sd ~ normal(0, 5);
  //sigma3 ~ normal(0, sigma3_sd);
  //sigma3 ~ normal(0, 5);
  //sigma3_mu ~ normal(0, 50);
  //sigma3_sd ~ normal(0, 5);
  
  //gamma1 ~ normal(0, gamma1_sd);
  //gamma1_mu ~ normal(0, 50);
  //gamma1_sd ~ normal(0, 1);
  //gamma2 ~ normal(0, gamma2_sd);
  //gamma2_mu ~ normal(0, 50);
  //gamma2_sd ~ normal(0, 1);
  //gamma3 ~ normal(0, gamma3_sd);
  //gamma3_mu ~ normal(0, 50);
  //gamma3_sd ~ normal(0, 1);
  
  //level ~ normal(0, 10);
  //level ~ normal(level_preds, 5);
  //scale_shape ~ lognormal(0, 1);
  //scale ~ beta(scale_preds, 1);
  //preference ~ normal(preference_preds, 10);
  //choosiness ~ gamma(choosiness_preds, 1);
  //asymmetry ~ normal(asymmetry_preds, 10);
  
  y ~ bernoulli_logit(logit(y_hat));
}
