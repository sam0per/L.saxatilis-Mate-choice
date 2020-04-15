rm(list = ls())

#############################
# Install required packages #
#############################

# List of packages for session
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "Rmisc", "boot", "pROC",
              "purrr", "reshape2", "gridExtra", "grid", "parallel", "projpred", "optparse")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

# library(pROC)
# library(rstan)
# library(dplyr)
# library(boot)

CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")

# gaus_size = readRDS("models/gaus_size/gaus_size_stan.rds")
skew_hier = readRDS("models/gaus_skew/gaus_skew_hier_BCDG_shore.rds")
# y_rep_size = summary(gaus_size, pars = c("y_rep"))$summary[,'mean']
y_rep_size = summary(skew_hier, pars = c("y_rep"))$summary[,'mean']
roc_obj_size = roc(CZ_data$mountYNcontact, y_rep_size, ci = TRUE)
auc(roc_obj_size)
ci.auc(roc_obj_size)

gaus_sex = readRDS("models/gaus_sex/gaus_sex.rds")
y_rep_sex = summary(gaus_sex, pars = c("y_rep"))$summary[,'mean']
roc_obj_sex = roc(CZ_data$mountYNcontact, y_rep_sex, ci = TRUE)
auc(roc_obj_sex)
ci.auc(roc_obj_sex)

roc.test(roc_obj_size, roc_obj_sex, method = "bootstrap")

# scale[i] = exp(lambda1[shore[i]] + lambda2[ref[i]] + lambda3[test[i]] * shape[i])
post_sex <- rstan::extract(gaus_sex)
names(post_sex)
dim(post_sex)
sex_scale = colMeans(post_sex$scale)
#sex_scale = rowMeans(post_sex$scale)
head(post_sex$lambda1)
colMeans(post_sex$lambda1)
dim(post_sex$lambda1)
as.numeric(CZ_data$shore)==1

int_lambda1 = sapply(1:4, function(x) ifelse(as.numeric(CZ_data$shore)==x, x, 0))
coef_lambda1 = sapply(1:4, function(x) int_lambda1[,x] * colMeans(post_sex$lambda1)[x]) %>%
  rowSums()
head(coef_lambda1)
tail(coef_lambda1)

int_lambda2 = sapply(1:4, function(x) ifelse(as.numeric(CZ_data$ref_eco_sex)==x, x, 0))
coef_lambda2 = sapply(1:4, function(x) int_lambda2[,x] * colMeans(post_sex$lambda2)[x]) %>%
  rowSums()
head(coef_lambda2)
tail(coef_lambda2)

int_lambda3 = sapply(1:2, function(x) ifelse(as.numeric(CZ_data$test_sex)==x, x, 0))
coef_lambda3 = sapply(1:2, function(x) int_lambda3[,x] * colMeans(post_sex$lambda3)[x] * CZ_data$shape) %>%
  rowSums()
head(coef_lambda3)
tail(coef_lambda3)

#res_lambda1 = sex_scale - exp(coef_lambda1 + coef_lambda2)
res_lambda1 = sex_scale - coef_lambda1
res_lambda1 = sex_scale - coef_lambda2
res_lambda1 = sex_scale - coef_lambda3
res_lambda1 = sex_scale - coef_lambda1 + coef_lambda2 + coef_lambda3


#Residual sum of squares:
RSS <- c(crossprod(res_lambda1))
#Mean squared error:
MSE <- RSS / length(res_lambda1)
#Root MSE:
RMSE <- sqrt(MSE)
#Pearson estimated residual variance (as returned by summary.lm):
sig2 <- RSS / res$df.residua

sqrt(mean(res_lambda1^2))

hist(post_sex$level)
length(post_sex$level)
hist(post_sex$scale)
dim(post_sex$scale)
length(rowMeans(post_sex$scale))
hist(rowMeans(post_sex$scale))
hist(post_sex$lambda1[,1])
hist(post_sex$lambda1[,2])

hist(sex_scale, breaks = 30, freq = FALSE)
hist(inv.logit(post_sex$level + sex_scale), breaks = 30)

points(x = density(sex_scale)$x, y = density(sex_scale)$y)
length(density(sex_scale)$x)
density(sex_scale)$x
hist(coef_lambda1 + coef_lambda2 + coef_lambda3, breaks = 30, freq = FALSE)
