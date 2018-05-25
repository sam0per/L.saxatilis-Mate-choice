rm(list = ls())

#############################
# Install required packages #
#############################
sessionInfo()
if (!require("rstan")) install.packages("rstan")
(.packages())
# List of packages for session
.packages = c("ggplot2", "plyr", "rms")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

load("Rpackages")
for (p in setdiff(packages, installed.packages()[,"Package"])){
  install.packages(p)
}
#############################
### read and adjust data ###
#############################

CZ_data = read.csv("data/CZ_mating_clean.csv",sep = ";")
CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
head(CZ_data)
summary(CZ_data)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#############################
# stan model with only size #
#############################
rm(CZ_mat_stan_size)
CZ_mat_stan_size = stan(file = "scripts/CZ_mating_gaus_size.stan", data = list(N = nrow(CZ_data),
                                                                               y = CZ_data$mountYN,
                                                                               ratio = CZ_data$size_ratio))

stan_dens(CZ_mat_stan_size, pars = c("level","lambda","aver","stan_dev","gamma_ta"))
CZ_size_params = round(summary(CZ_mat_stan_size, pars = c("level","lambda","aver","stan_dev","gamma_ta"),
                               probs=c(0.25, 0.975))$summary,2)
write.table(CZ_size_params, "tables/CZ_size_params.csv", row.names = TRUE, col.names = TRUE,sep = ";")

y_rep = data.frame(y_rep=summary(CZ_mat_stan_size, pars = c("y_rep"))$summary[,'mean'])
write.table(y_rep, "tables/CZ_size_mount_rep.csv", row.names = TRUE, col.names = TRUE,sep = ";")

CZ_data$y_rep = summary(CZ_mat_stan_size, pars = c("y_rep"))$summary[,'mean']
CZ_data$y_rep_se = summary(CZ_mat_stan_size, pars = c("y_rep"))$summary[,'se_mean']


library(boot)
y_hat = round(summary(CZ_mat_stan_size, pars = c("y_hat"))$summary, 2)
write.table(y_hat, "tables/CZ_size_mount_hat.csv", row.names = TRUE, col.names = TRUE,sep = ";")

CZ_logit = summary(CZ_mat_stan_size, pars = c("y_hat"))$summary[,'mean']
CZ_logit_se = summary(CZ_mat_stan_size, pars = c("y_hat"))$summary[,'se_mean']
CZ_data$preds = inv.logit(CZ_logit)
CZ_data$y_preds = rbinom(n = nrow(CZ_data),size = 1,prob = CZ_data$preds)
CZ_data$uci_preds = exp(CZ_logit+1.96*CZ_logit_se)/(1+exp(CZ_logit+1.96*CZ_logit_se))
CZ_data$lci_preds = exp(CZ_logit-1.96*CZ_logit_se)/(1+exp(CZ_logit-1.96*CZ_logit_se))

################################
################################
CZ_logit = CZ_size_params['level','mean'] + (1 / sqrt(2 * pi * CZ_size_params['stan_dev','mean']^2)) * 
  exp(-0.5*((CZ_data$size_ratio-CZ_size_params['aver','mean'])/CZ_size_params['stan_dev','mean'])^2) + 
  CZ_size_params['gamma_ta','mean'] * CZ_data$size_ratio

CZ_lci_logit = CZ_size_params['level','25%'] + (1 / sqrt(2 * pi * CZ_size_params['stan_dev','25%']^2)) * 
  exp(-0.5*((CZ_data$size_ratio-CZ_size_params['aver','25%'])/CZ_size_params['stan_dev','25%'])^2) + 
  CZ_size_params['gamma_ta','25%'] * CZ_data$size_ratio
CZ_data$uci_preds = inv.logit(CZ_lci_logit)
CZ_data$uci_preds = NULL
CZ_data$uci_preds = CZ_data$preds + CZ_data$y_rep_se
CZ_data$uci_preds = exp(CZ_logit+1.96*CZ_data$y_rep_se)/(1+exp(CZ_logit+1.96*CZ_data$y_rep_se))
CZ_data$uci_preds = CZ_data$preds+1.96*CZ_logit_se


CZ_uci_logit = CZ_size_params['level','97.5%'] + (1 / sqrt(2 * pi * CZ_size_params['stan_dev','97.5%']^2)) * 
  exp(-0.5*((CZ_data$size_ratio-CZ_size_params['aver','97.5%'])/CZ_size_params['stan_dev','97.5%'])^2) + 
  CZ_size_params['gamma_ta','97.5%'] * CZ_data$size_ratio
CZ_data$lci_preds = inv.logit(CZ_uci_logit)
CZ_data$lci_preds = NULL
CZ_data$lci_preds = CZ_data$preds - CZ_data$y_rep_se
CZ_data$lci_preds = exp(CZ_logit-1.96*CZ_data$y_rep_se)/(1+exp(CZ_logit-1.96*CZ_data$y_rep_se))
CZ_data$lci_preds = CZ_data$preds-1.96*CZ_logit_se

library(rstanarm)
launch_shinystan(CZ_mat_stan_size)

###################################
# plot observs and preds for size #
###################################                                                    
#install.packages("bayesplot")
library(bayesplot)
library(ggplot2)

y = CZ_data$mountYN
y_rep = extract(CZ_mat_stan_size, pars = 'y_rep', permuted = TRUE)$y_rep
dim(y_rep)
ppc_bars(y,y_rep)
breaks = c(-2,seq(-1.5,1.5,0.1),2)
bin = cut(CZ_data$size_ratio,breaks)
ppc_bars_grouped(y, y_rep, bin, prob = 0, freq = FALSE)

library(dplyr)
#install.packages("Rmisc")
library(Rmisc)
range(CZ_data$size_ratio)
breaks = c(-2,seq(-1.5,1.5,0.1),2)
CZ_data$bin = cut(CZ_data$size_ratio,breaks)
CZ_data_bin <- 
  CZ_data %>%
  group_by(bin) %>%
  dplyr::summarise(mount = mean(mountYN), 
                   uci_mount = CI(mountYN)['upper'], 
                   lci_mount = CI(mountYN)['lower'],
                   mean_ratio = mean(size_ratio),
                   y_rep = mean(y_rep),
                   preds_mount = mean(y_preds)) %>% 
  mutate(lci_mount = replace(lci_mount, which(lci_mount<0), 0))
summary(CZ_data_bin)

ggplot(data = CZ_data) +
  geom_vline(xintercept = 0) +
  geom_ribbon(aes(x = size_ratio,ymin = lci_preds, ymax = uci_preds), fill = "orange", alpha=0.3) +
  geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount, col="observations")) +
  geom_errorbar(data = CZ_data_bin, aes(x = mean_ratio, ymin = lci_mount, ymax = uci_mount),alpha = 0.2) +
  scale_colour_manual(values=c("blue","orange2")) +
  geom_line(aes(size_ratio,preds,col="predictions")) +
  labs(size="bin size",x="female - male size (log)",
       y="probability of mounting",col="") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"))

ggplot(CZ_data_bin,aes(mount,y_rep)) +
  geom_abline(slope = 1) +
  geom_point() +
  labs(x='observed mount proportion', y='predicted probability')

#######################
#######################
library(rstanarm)
library(arm)
yrep_prop <- sweep(y_rep, 2, trials, "/")
apply(y_rep,2,FUN = mean)
sweep(y_rep,2,bin,FUN = mean)
library(Matrix.utils)
aggregate.Matrix(y_rep,groupings = bin,fun = "mean")

ppc_error_binned(CZ_data_bin$mount,CZ_data_bin$preds_mount)

range(CZ_data$size_ratio)
breaks = c(-2,seq(-1.5,1.5,0.1),2)
bin = cut(CZ_data$size_ratio,breaks)
p_m_obs = aggregate(CZ_data$mountYN,by=list(bin),FUN=mean)[,2]
p_m_preds = aggregate(CZ_data$y_preds,by=list(bin),FUN=mean)[,2]
bin_mean = aggregate(CZ_data$size_ratio,by=list(bin),FUN=mean)[,2]
CZ_data$count = 1
bin_sum = aggregate(CZ_data$count,by=list(bin),FUN=sum)[,2]
CZ_data_bin = data.frame(cbind(mount=p_m_obs,mount_rep=p_m_preds,size_ratio=bin_mean,
                               count=bin_sum))
summary(CZ_data_bin)

ggplot(CZ_data_bin,aes(mount,preds_mount)) +
  geom_abline(slope = 1) +
  geom_point()

ggplot(CZ_data_bin,aes(preds_mount, y_rep)) +
  geom_abline(slope = 1) +
  geom_point()

CZ_data_bin %>%
  ggplot(aes(x = mean_ratio, y = mount, col="observations")) +
  geom_point() +
  geom_errorbar(aes(ymin = lci_mount, ymax = uci_mount),alpha = 0.2) +
  geom_vline(xintercept = 0) +
  scale_colour_manual(values=c("blue","orange")) +
  geom_line(data = CZ_data, aes(size_ratio,preds,col="predictions")) +
  geom_line(data = CZ_data, aes(size_ratio,lci_preds,col="predictions"),linetype='dashed',alpha=0.5) +
  geom_line(data = CZ_data, aes(size_ratio,uci_preds,col="predictions"),linetype='dashed',alpha=0.5) +
  labs(size="bin size",x="female - male size (log)",
       y="probability of mounting",col="") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"))


CZ_data_bin %>%
  ggplot(aes(x = mean_ratio, y = mount)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lci_mount, ymax = uci_mount), position = "dodge")


ggplot(data = CZ_data_bin,aes(size_ratio,mount,col="observations")) +
  geom_point(aes(size=CZ_data_bin$count)) +
  geom_vline(xintercept = 0) +
  scale_colour_manual(values=c("blue","orange")) +
  geom_point(aes(size_ratio,mount_rep,col="predictions",size=CZ_data_bin$count)) +
  labs(size="bin size",x="female - male size (log)",
       y="probability of mounting",col="") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"))
  
ggplot(data = CZ_data_bin,aes(size_ratio,mount,col="observations")) +
  geom_point(aes(size=CZ_data_bin$count)) +
  geom_vline(xintercept = 0) +
  scale_colour_manual(values=c("blue","orange")) +
  geom_point(data = CZ_data, aes(size_ratio,y_rep,col="predictions")) +
  labs(size="bin size",x="female - male size (log)",
       y="probability of mounting",col="") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"))

ggplot(data = CZ_data_bin,aes(size_ratio,mount,col="observations")) +
  geom_point(aes(size=CZ_data_bin$count)) +
  geom_vline(xintercept = 0) +
  scale_colour_manual(values=c("blue","orange")) +
  geom_line(data = CZ_data, aes(size_ratio,preds,col="predictions")) +
  labs(size="bin size",x="female - male size (log)",
       y="probability of mounting",col="") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"))


###################################
# apply size model to field distr #
###################################
devtools::install_github("rmcelreath/rethinking",force = TRUE)
library(rethinking)

CZA_data = read.csv("../2.mating/CZA/final_data/CZA_use_cleanup.csv", sep = ";")
with(CZA_data, plot(DistAlongPath, length_mm))

CZB_data = read.csv("../2.mating/CZB/final_data/CZA_use_cleanup.csv", sep = ";")

cline_2c3s <- function(phen,position,cl,cr,wl,wr,wave,crab,sc,sw,sh){
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  s_xl <- sqrt(sc^2 + 4*p_xl*(1-p_xl)*sh^2 + (p_xl^2)*(sw^2-sc^2))
  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing 
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  # combined cline
  z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  s_x[z_x < z_xl] <- s_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  return(minusll)
}

hist(rnorm(10000,mean = 2.4, sd = sqrt(0.2)))
crab <- data.frame(log_male=rnorm(10000,mean = 2.4,
                                  sd = sqrt(0.2)),
                   log_female=rep(rnorm(1000,mean = 2.5,
                                        sd = sqrt(0.3)),
                                  10),
                   ref_ecotype=rep("crab",10000))
wave <- data.frame(log_male=rnorm(10000,mean = 1.5,
                                  sd = sqrt(0.2)),
                   log_female=rep(rnorm(1000,mean = 1.7,
                                        sd = sqrt(0.3)),
                                  10),
                   ref_ecotype=rep("wave",10000))
hyb <- data.frame(log_male=rnorm(10000,mean = (2.4 + 1.5)/2,
                                 sd = sqrt(0.4)), # hyb male variance is closed to 0!
                  log_female=rep(rnorm(1000,mean = (2.5 + 1.7)/2,
                                       sd = sqrt(0.4)),
                                 10),
                  ref_ecotype=rep("hyb",10000))

# build the new dataframe from the wave and crab clines
CZ_cline <- merge(crab,wave,all=TRUE)
CZ_cline <- merge(CZ_cline,hyb,all=TRUE)
CZ_cline$newmale=NULL

crab$newmale=NULL
crab$size_ratio <- crab$log_female-crab$log_male
crab$logit = -1.23 + (1 / sqrt(2 * pi * 0.27^2)) * exp(-0.5*((crab$size_ratio-0.2)/0.27)^2) + 
  0.58 * crab$size_ratio
crab$preds = inv.logit(crab$logit)
crab$mount = rbinom(n = nrow(crab),size = 1,prob = crab$preds)

cztest <- CZ_cline %>% 
  nest(-ref_ecotype) %>%
  mutate(data = map(data, ~mutate(.x,
                                  ratio = log_female - log_male,
                                  logit = -1.23 + (1 / sqrt(2 * pi * 0.27^2)) * exp(-0.5*((ratio-0.2)/0.27)^2) + 
                                    0.58 * ratio,
                                  preds = inv.logit(logit),
                                  mount = rbinom(n = 10000,size = 1,prob = preds)))) %>%
  unnest() %>%
  group_by(ref_ecotype,log_female,mount) %>%
  dplyr::summarise(n=n(), mean_male = mean(log_male)) %>%
  filter(n<10)
  


mutate(YES = map(data, ~filter(.x, mount==1))) %>%
  mutate(NO = map(data, ~filter(.x, mount==0)))
  

cztest$N

head(crab)
CZidx = crab$log_female[crab$mount==0] %in% crab$log_female[crab$mount==1]
crab[CZidx,]

CZ_cline_YN[[i]]$`0`[CZ_cline_YN$`0`[, 2] %in% CZ_cline_YN$`1`[, 2], ]


with(crab, hist(log_female-log_male))
with(wave, hist(log_female-log_male))

eco_list = list(wave=wave,hyb=hyb,crab=crab)
#fem_list = list(wave$log_female,hyb$log_female,crab$log_female)
crabba = crab$log_female
for (i in nrow(crab)){
  male = sample(crab$log_male,size = 10,replace = TRUE)
  ratio = crab$log_female[i] - male
  logit = -1.23 + (1 / sqrt(2 * pi * 0.27^2)) * exp(-0.5*((ratio-0.2)/0.27)^2) + 
    0.58 * ratio
  pred = inv.logit(logit)
  mount = rbinom(n = 1,size = 1,prob = pred)
  print(mount)
}
  



print(sample(crab$log_male,size = 4,replace = TRUE))

CZ_idx <- eco_list %>%
  map_df(function(x) x$log_female - x$log_male) %>%
  map(function(x) -1.23 + (1 / sqrt(2 * pi * 0.27^2)) * exp(-0.5*((x-0.2)/0.27)^2) + 
        0.58 * x) %>%
  map(function(x) inv.logit(x)) %>%
  map_dfr(rbinom, n=10000, size=1)

test %>% nest(eco_list) 
CZ_idx <- eco_list %>%
  mutate(.data = map())
  map_df(function(x) x$log_female - x$log_male)
  


map(function(x) x == 1)

lapply(CZ_idx, select, 1)
lapply(CZ_idx, "[", 1:2)
CZ_Y = eco_list[CZ_idx]

eco_list_ratio = lapply(eco_list, function(x) x$log_female - x$log_male)
eco_list_logit = lapply(eco_list_ratio, function(x) -1.23 + (1 / sqrt(2 * pi * 0.27^2)) * exp(-0.5*((x-0.2)/0.27)^2) + 
                          0.58 * x)
eco_list_preds = lapply(eco_list_logit, function(x) inv.logit(x))
library(purrr)
library(magrittr)
eco_list_mount = map(eco_list_preds, rbinom, n=10000, size=1)
do.call('rbind',eco_list_mount) 
eco_list_mount %>%
  split(.)
split(unlist(eco_list_mount),as.factor(eco_list_mount))

lst <- list(a = data.frame(1:3, 4:6, 7), b = data.frame(1:3, 4:6, 8), 
            c = data.frame(1:3, 4:6, 7))
idx <- sapply(lst, function(x) x[1,3] == 7)
lst1 = lst[idx] 
lst2 = lst[!idx] 

crab$preds = eco_list_preds[[1]]
wave$preds = eco_list_preds[[2]]
hyb$preds = eco_list_preds[[3]]

eco_list = list(wave,hyb,crab)
head(eco_list)
lapply(1:length(eco_list), function(i) {
  y_preds <-  rbinom(n = nrow(eco_list[[i]]), size = 1, prob = eco_list[[i]]$preds)
})

map(mtcars, mean) %>% head
map_dbl(mtcars, mean) %>% head
mumu <- list(10, 100, -100)
sigmaa <- list(0.01, 1, 10)
map2(mumu, sigmaa, rnorm, n = 5)


for (i in seq_along(eco_list)){
  eco_list[[i]]$y_preds = rbinom(n = nrow(eco_list[[i]]), size = 1, prob = eco_list[[i]]$preds)
  eco_YN[[i]] <-  split(eco_list[[i]],as.factor(eco_list[[i]]$y_preds))
  
}


[CZ_cline_YN$`0`[, 2] %in% CZ_cline_YN$`1`[, 2], ] # NO mounts for pairs with mated females
CZ_cline_N[[i]] = CZ_cline_YN[[i]]$`0`[CZ_cline_YN$`0`[, 2] %in% CZ_cline_YN$`1`[, 2], ] # NO mounts for pairs with mated females
CZ_cline_Y = CZ_cline_YN$`1`
CZ_cline_clean = merge(CZ_cline_N,CZ_cline_Y,all = TRUE)
print(summary(CZ_cline_clean))

mapply(cbind, eco_list, function(x) rbinom(n = nrow(x),size = 1,prob = x$preds), SIMPLIFY=F)

crab$size_ratio <- CZ_cline$log_female-CZ_cline$log_male

CZ_cline_logit = -1.23 + (1 / sqrt(2 * pi * 0.27^2)) * exp(-0.5*((CZ_cline$size_ratio-0.2)/0.27)^2) + 
  0.58 * CZ_cline$size_ratio
CZ_cline$preds = inv.logit(CZ_cline_logit)
plot(CZ_cline$size_ratio,CZ_cline$preds)

CZ_cline$y_preds = rbinom(n = nrow(CZ_cline),size = 1,prob = CZ_cline$preds)

CZ_cline_YN <- split(CZ_cline,as.factor(CZ_cline$y_preds))
names(CZ_cline_YN$`0`)
names(CZ_cline_YN$`1`)
CZ_cline_N <- CZ_cline_YN$`0`[CZ_cline_YN$`0`[, 2] %in% CZ_cline_YN$`1`[, 2], ] # NO mounts for pairs with mated females
CZ_cline_Y <- CZ_cline_YN$`1`

CZ_cline_clean <- merge(CZ_cline_N,CZ_cline_Y,all = TRUE)

hist(CZ_cline_clean$log_female,breaks = 50)
hist(CZ_cline_clean$log_male,breaks = 50)

diff(by(CZ_cline_clean$log_male, CZ_cline_clean$y_preds, mean))
mean(CZ_cline_clean$log_male[CZ_cline_clean$y_preds==1])-
  mean(CZ_cline_clean$log_male[CZ_cline_clean$y_preds==0])
diff(by(CZ_cline_clean$log_male, sample(CZ_cline_clean$y_preds, length(CZ_cline_clean$y_preds), FALSE), mean))

male_rnd_diff <- replicate(1000, diff(by(CZ_cline_clean$log_male, sample(CZ_cline_clean$y_preds,
                                                                         length(CZ_cline_clean$y_preds), 
                                                                         FALSE), mean)))

hist(male_rnd_diff, xlim = c(-0.15, 0.15), col = "black", breaks = 100)
abline(v = diff(by(CZ_cline_clean$log_male, CZ_cline_clean$y_preds, mean)), col = "blue", lwd = 2)



CZ_cline_clean %>%
  group_by(as.factor(y_preds)) %>%
  dplyr::summarise(males = length(log_male))

CZ_cline_clean_male = data.frame(male_Y = sample(CZ_cline_clean$log_male[CZ_cline_clean$y_preds==1],size = 10000),
                                 male_N = sample(CZ_cline_clean$log_male[CZ_cline_clean$y_preds==0],size = 10000))
with(CZ_cline_clean_male, mean(male_Y - male_N))

ggplot(data = CZ_cline_clean_male, aes(male_Y, fill="mated")) +
  geom_histogram(binwidth = 0.1, stat = "identity") +
  geom_histogram(aes(male_N, fill="unmated"),binwidth = 0.1)

ggplot(data = CZ_cline_clean, aes(log_male, fill = as.factor(y_preds))) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = mean(CZ_cline_clean_male$male_Y), size=1.5, col='red') +
  geom_vline(xintercept = mean(CZ_cline_clean_male$male_N), size=1.5, col='blue') +
  scale_fill_manual(values = c('blue','red'), name='Males', labels=c('unmated','mated')) +
  theme(legend.title = element_text(size = 12,face = "bold"),
        axis.title = element_text(face = "bold")) +
  labs(x='male size (log)')

library(reshape2)
library(dplyr)

range(CZ_cline_clean$log_male)
breaks = c(-0.5,seq(0.5,3.5,0.15),5)
CZ_cline_clean$bin = cut(CZ_cline_clean$log_male,breaks)
CZ_cline_clean_bin <- 
  CZ_cline_clean %>%
  group_by(bin, ref_ecotype) %>%
  dplyr::summarise(mount = mean(y_preds),
                   uci_mount = CI(y_preds)['upper'], 
                   lci_mount = CI(y_preds)['lower'],
                   count = length(log_male),
                   mean_size = mean(log_male)) %>% 
  mutate(lci_mount = replace(lci_mount, which(lci_mount<0), 0)) %>%
  mutate(uci_mount = replace(uci_mount, which(uci_mount>1), 1))

factor(CZ_cline_clean_bin$ref_ecotype)
CZ_cline_clean_bin$ref_ecotype <- factor(CZ_cline_clean_bin$ref_ecotype, levels = c("wave","hyb","crab"))


ggplot(data = CZ_cline_clean_bin, aes(mean_size, mount)) +
  facet_wrap(~ref_ecotype) +
  geom_point() +
  geom_errorbar(aes(x = mean_size, ymin = lci_mount, ymax = uci_mount),alpha = 0.2) +
  labs(x='male size (log)', y='predicted probability of mounting') +
  theme(axis.title = element_text(face = "bold"))

hist(CZ_cline_clean$log_male, breaks = 50, col = as.factor(CZ_cline_clean$y_preds))
abline(v = mean(CZ_cline_clean_male$male_Y), col = 'red')
abline(v = mean(CZ_cline_clean_male$male_N), col = 'blue')

ggplot(subset(CZ_cline_clean, y_preds==1),aes(log_female,log_male, col=preds)) +
  geom_point() +
  geom_abline(slope = 1) + 
  labs(x='female size (log)', y='male size (log)', col='prob.') +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(size = 12,face = "bold"))

############################################
# stan model with shore, ecotype and shape #
############################################
stan_run_mat2 = stan(file = "stan/mat_gau2.stan",data = list(N = nrow(CZ_data),
                                                             y = CZ_data$mountYN,
                                                             ratio = CZ_data$size_ratio,
                                                             shape = CZ_data$Shape,
                                                             N_shore = length(unique(CZ_data$shore)),
                                                             N_eco = length(unique(CZ_data$ref_ecotype)),
                                                             shore = CZ_data$shore,
                                                             eco = CZ_data$ref_ecotype))

traceplot(stan_run_mat2, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"), inc_warmup = TRUE)
traceplot(stan_run_mat2, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"), inc_warmup = TRUE)
traceplot(stan_run_mat2, pars = c("delta","delta_sd","delta_av","delta_ta"), inc_warmup = TRUE)

CZ_gau2_alpha = round(summary(stan_run_mat2, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"))$summary,2)
CZ_gau2_delta = round(summary(stan_run_mat2, pars = c("delta","delta_sd","delta_av","delta_ta"))$summary,2)
CZ_gau2_beta = round(summary(stan_run_mat2, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"))$summary,2)

write.table(CZ_gau2_alpha, "CZ_ABCD/results/CZABCD_gau2_alpha.csv", row.names = TRUE, col.names = TRUE,sep = ";")
write.table(CZ_gau2_delta, "CZ_ABCD/results/CZABCD_gau2_delta.csv", row.names = TRUE, col.names = TRUE,sep = ";")
write.table(CZ_gau2_beta, "CZ_ABCD/results/CZABCD_gau2_beta.csv", row.names = TRUE, col.names = TRUE,sep = ";")

stan_dens(stan_run_mat2, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"))
stan_dens(stan_run_mat2, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"))
stan_dens(stan_run_mat2, pars = c("delta","delta_sd","delta_av","delta_ta"))

stan_plot(stan_run_mat2, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"))
stan_plot(stan_run_mat2, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"))
stan_plot(stan_run_mat2, pars = c("delta","delta_sd","delta_av","delta_ta"))




pp_check(stan_run_mat2,fun = "dens_overlay")
launch_shinystan(stan_run_mat2)
library(scales)
ggplot(CZ_bayes,aes(sim_le,y_rep,shape=ref_ecotype,col=Shape)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient2(low = "black",mid = "grey", high = "white") +
  labs(shape="ref. ecotype", col="shape", x="level",
       y="predicted probability", title = "Predicted curve intercept") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(log_female,log_male,shape=ref_ecotype,col=sim_le)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient2(low = "black", mid = "grey", high = "white", midpoint = -0.7) +
  labs(shape="ref. ecotype", col="level", x="female size (ln)",
       y="male size (ln)", title = "Predicted curve intercept") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))


ggplot(CZ_bayes,aes(sim_ta,y_rep,shape=ref_ecotype,col=Shape)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient2(low = "black",mid = "grey", high = "white") +
  labs(shape="ref. ecotype", col="shape", x="tail asymmetry",
       y="predicted probability", title = "Predicted tail asymmetry") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(log_female,log_male,shape=ref_ecotype,col=sim_ta)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient2(low = "black", mid = "grey", high = "white", midpoint = -0.5) +
  labs(shape="ref. ecotype", col="tail asymmetry", x="female size (ln)",
       y="male size (ln)", title = "Predicted tail asymmetry") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(log_female,log_male,shape=sex,col=sim_ta)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient2(low = "black", mid = "grey", high = "white", midpoint = -0.5) +
  labs(shape="sex", col="tail asymmetry", x="female size (ln)",
       y="male size (ln)", title = "Predicted tail asymmetry") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(CZ_bayes,aes(size_ratio,sim_sd,shape=ref_ecotype,col=Shape)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient(low = "black", high = "red") +
  labs(shape="ref. ecotype", col="shape", x="size ratio",
       y="sigma", title = "Predicted sexual selection") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(size_ratio,sim_sd,shape=sex,col=Shape)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient(low = "black", high = "red") +
  labs(shape="sex", col="shape", x="size ratio",
       y="sigma", title = "Predicted sexual selection") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(log_female,log_male,shape=ref_ecotype,col=sim_sd)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient2(low = "black",mid = "grey", high = "white",midpoint = 0.26) +
  labs(shape="ref. ecotype", col=expression(sigma~"(exp)"), x="female size",
       y="male size", title = "Predicted sexual selection") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(CZ_bayes,aes(log_female,log_male,shape=ref_ecotype,col=sim_av)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient2(low = "black",mid = "grey", high = "white",midpoint = 0.15) +
  labs(shape="ref. ecotype", col=expression(mu), x="female size",
       y="male size", title = "Predicted disassortative mating") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(CZ_bayes,aes(log_female,log_male,shape=ref_ecotype,col=y_rep)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  scale_color_gradient(low = "black", high = "red") +
  labs(shape="ref. ecotype", col="mount prob", x="female size",
       y="male size", title = "") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))


CZ_bayes$y_rep = summary(stan_run_mat2, pars = c("y_rep"))$summary[,'mean']

CZ_bayes$line = 0
ggplot(data = CZ_bayes,aes(size_ratio,y_rep,col=as.factor(ref_ecotype))) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref ecotype",size="bin size",x="shape test snails",
       y="predicted probability of mounting", title = "Predicted mounting probability") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(data = CZ_bayes,aes(Shape,y_rep,col=as.factor(ref_ecotype))) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref ecotype",size="bin size",x="ln(female size) - ln(male size)",
       y="predicted probability of mounting", title = "Predicted mounting probability") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))


ggplot(data = CZ_bayes,aes(size_ratio+Shape,y_rep,col=as.factor(ref_ecotype))) +
  #facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref ecotype",size="bin size",x="size_ratio + shape",
       y="predicted probability of mounting", title = "Predicted mounting probability") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

pars_sd_mean = summary(stan_run_mat2, pars = c('alpha_sd','delta_sd','beta_sd'))$summary[,'mean']
pars_sd = summary(stan_run_mat2, pars = c('alpha_sd','delta_sd','beta_sd'))$summary
write.table(pars_sd, "CZ_ABCD/results/CZABCD_gau_sd.csv", row.names = TRUE, col.names = TRUE,sep = ";")
#pars_sd_mean = pars_sd[,'mean']
library(dplyr)
group_by(CZ_bayes, shore) %>%
  summarise(
    count = n(),
    mean = mean(sim_sd, na.rm = TRUE),
    sd = sd(sim_sd, na.rm = TRUE)
  )
group_by(CZ_bayes, shore, ref_ecotype) %>%
  summarise(
    count = n(),
    mean = mean(sim_sd, na.rm = TRUE),
    sd = sd(sim_sd, na.rm = TRUE)
  )
library("ggpubr")
ggboxplot(CZ_bayes, x = "shore", y = "sim_sd", color = "ref_ecotype",
          palette = c("gold2", "blue")) +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(CZ_bayes,aes(Shape,sim_sd,col=ref_ecotype)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref. ecotype", x="Shape",
       y=expression("Trait variation in all mating pairs"~(sigma)), title = "Predicted sexual selection") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(size_ratio,Shape,shape=ref_ecotype,col=sim_sd)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.2) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_gradient(low = "black",high = "red") +
  labs(shape="ref. ecotype", x="female size - male size (ln)",
       col=expression(sigma), title = "Predicted sexual selection") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(CZ_bayes,aes(Shape,sim_le,col=ref_ecotype)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref. ecotype", x="Shape",
       y="Level or intercept", title = "Predicted curve intercept") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(size_ratio,Shape,shape=ref_ecotype,col=sim_le)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.2) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_gradient(low = "black",high = "red") +
  labs(shape="ref. ecotype", x="",
       col="Level or intercept", title = "Predicted curve intercept") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))


ggplot(CZ_bayes,aes(size_ratio,Shape,shape=ref_ecotype,col=sim_ta)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point() +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_gradient(low = "black", high = "red") +
  labs(shape="ref. ecotype", col="tail asymmetry", x="female size - male size (ln)",
       y="shape", title = "Predicted tail asymmetry") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))


pars_av_mean = summary(stan_run_mat2, pars = c('alpha_av','delta_av','beta_av'))$summary[,'mean']
pars_av = summary(stan_run_mat2, pars = c('alpha_av','delta_av','beta_av'))$summary
write.table(pars_av, "CZ_ABCD/results/CZABCD_gau_av.csv", row.names = TRUE, col.names = TRUE,sep = ";")
#pars_av_mean = pars_av[,'mean']

pars_ta_mean = summary(stan_run_mat2, pars = c('alpha_ta','delta_ta','gamma_ta','beta_ta'))$summary[,'mean']
sim_ta = pars_ta_mean[1:4][CZ_bayes$shore] + pars_ta_mean[5:6][CZ_bayes$ref_ecotype] + 
  pars_ta_mean[7] * CZ_bayes$size_ratio + pars_ta_mean[8] * CZ_bayes$Shape
CZ_bayes$sim_ta = sim_ta
ggboxplot(CZ_bayes, x = "shore", y = "sim_ta", color = "ref_ecotype",
          palette = c("gold2", "blue")) +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

pars_le_mean = summary(stan_run_mat2, pars = c('alpha','delta','beta'))$summary[,'mean']
sim_le = pars_le_mean[1:4][CZ_bayes$shore] + pars_le_mean[5:6][CZ_bayes$ref_ecotype] + 
  pars_le_mean[7] * CZ_bayes$Shape
CZ_bayes$sim_le = sim_le
ggboxplot(CZ_bayes, x = "shore", y = "sim_le", color = "ref_ecotype",
          palette = c("gold2", "blue")) +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

# exp in the stan code
sim_sd = exp(pars_sd_mean[1:4][CZ_bayes$shore] + pars_sd_mean[5:6][CZ_bayes$ref_ecotype] + pars_sd_mean[7] * CZ_bayes$Shape)
CZ_bayes$sim_sd = sim_sd
hist(CZ_bayes$sim_sd[CZ_bayes$mountYN==1])
ggplot(CZ_bayes,aes(sim_sd,fill=as.factor(ref_ecotype))) +
  geom_histogram(position = "dodge",bins = 10) +
  facet_wrap(~shore) +
  scale_fill_manual(values=c("gold2","blue")) +
  labs(fill="ref ecotype",x="Trait variation in all mating pairs",
       title = "Sexual selection") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

sim_av = pars_av_mean[1:4][CZ_bayes$shore] + pars_av_mean[5:6][CZ_bayes$ref_ecotype] + pars_av_mean[7] * CZ_bayes$Shape
CZ_bayes$sim_av = sim_av
hist(CZ_bayes$sim_av,col = CZ_bayes$shore,breaks = 25)
ggplot(subset(CZ_bayes,y_rep>=0.5),aes(sim_av,fill=as.factor(ref_ecotype))) +
  geom_histogram(position = "dodge",bins = 10) +
  facet_wrap(~shore) +
  scale_fill_manual(values=c("gold2","blue")) +
  labs(fill="ref ecotype",x="Size difference of successful pairs (ln)",
       title = "Disassortative mating") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(subset(CZ_bayes,mountYN==1),aes(size_ratio,fill=ref_ecotype)) +
  geom_histogram(position = "dodge",bins = 10) +
  facet_wrap(~shore) +
  scale_fill_manual(values=c("gold2","blue")) +
  labs(fill="ref ecotype",x="Size difference of successful pairs (ln)",
       title = "Disassortative mating") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(subset(CZ_bayes,y_rep>=0.5),aes(size_ratio,Shape,shape=ref_ecotype,col=sim_av)) +
  facet_wrap(~ as.factor(CZ_bayes$shore[CZ_bayes$y_rep>=0.5])) +
  geom_point(size=1.2) +
  geom_vline(xintercept = CZ_bayes$line[CZ_bayes$y_rep>=0.5]) +
  scale_color_gradient(low = "black",high = "red") +
  labs(shape="ref. ecotype", x="female size - male size (ln)",
       col=expression(mu), title = "Predicted disassortative mating") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))
ggplot(CZ_bayes,aes(Shape,sim_av,shape=ref_ecotype,col=y_rep)) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=1.5) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_gradient(low = "black",high = "red") +
  labs(shape="ref. ecotype",x="shape", y=expression(mu),
       col="mount prob", title = "Predicted disassortative mating") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

write.table(CZ_bayes, "CZ_ABCD/results/CZABCD_stan1.csv", row.names = FALSE, col.names = TRUE,sep = ";")
write.table(CZ_bayes, "CZ_ABCD/results/CZABCD_stan2.csv", row.names = FALSE, col.names = TRUE,sep = ";")

group_by(CZ_bayes, shore) %>%
  summarise(
    count = n(),
    mean = mean(sim_av, na.rm = TRUE),
    sd = sd(sim_av, na.rm = TRUE)
  )
group_by(subset(CZ_bayes,y_rep>=0.5), shore) %>%
  summarise(
    count = n(),
    mean = mean(sim_av, na.rm = TRUE),
    sd = sd(sim_av, na.rm = TRUE)
  )
group_by(subset(CZ_bayes,y_rep>=0.5),shore, ref_ecotype) %>%
  summarise(
    count = n(),
    mean = mean(sim_av, na.rm = TRUE),
    sd = sd(sim_av, na.rm = TRUE)
  )

group_by(subset(CZ_bayes,mountYN==1), shore) %>%
  summarise(
    count = n(),
    mean = mean(sim_av, na.rm = TRUE),
    sd = sd(sim_av, na.rm = TRUE)
  )

group_by(CZ_bayes, ref_ecotype) %>%
  summarise(
    count = n(),
    level = mean(sim_le, na.rm = TRUE),
    lev_sd = sd(sim_le, na.rm = TRUE),
    mu = mean(sim_av, na.rm = TRUE),
    mu_sd = sd(sim_av, na.rm = TRUE),
    sigma = mean(sim_sd, na.rm = TRUE),
    sig_sd = sd(sim_sd, na.rm = TRUE),
    tail = mean(sim_ta, na.rm = TRUE),
    tail_sd = sd(sim_ta, na.rm = TRUE)
  )

library(rstanarm)
CZ_av_aov = stan_aov(sim_av ~ shore + ref_ecotype + exp(Shape),data = CZ_bayes,prior = R2(0.5))
CZ_av_aov = stan_lm(sim_av ~ shore + ref_ecotype + exp(Shape),data = subset(CZ_bayes,y_rep>=0.5),prior = R2(0.5))
CZ_sd_aov = stan_aov(sim_sd ~ shore * ref_ecotype,data = CZ_bayes,prior = R2(0.5))

print(CZ_av_aov)
CZ_av_ci95 = posterior_interval(CZ_av_aov, prob = 0.95)
CZ_av_ci95 = round(CZ_av_ci95, 2)
write.table(CZ_av_ci95, "CZ_ABCD/results/CZABCD_av_ci95.csv", row.names = TRUE, col.names = TRUE,sep = ";")

library(car)
CZ_bayes_anova = aov(sim_av ~ shore * ref_ecotype, data = subset(CZ_bayes,y_rep>=0.5))
CZ_bayes_anova = aov(sim_av ~ shore + ref_ecotype, data = subset(CZ_bayes,y_rep>=0.5))
CZ_bayes_anova = aov(sim_av ~ shore, data = subset(CZ_bayes,y_rep>=0.5))
Anova(CZ_bayes_anova, type = "III")
group_by(subset(CZ_bayes,y_rep>=0.5), shore, ref_ecotype) %>%
  summarise(
    count = n(),
    mean = mean(sim_av, na.rm = TRUE),
    sd = sd(sim_av, na.rm = TRUE)
  )
TukeyHSD(CZ_bayes_anova,which = "shore")
TukeyHSD(CZ_bayes_anova)
# Homogeneity of variances
plot(CZ_bayes_anova,1)
leveneTest(sim_av ~ shore * ref_ecotype, data = CZ_bayes)
leveneTest(sim_av ~ shore, data = subset(CZ_bayes,y_rep>=0.5))
oneway.test(sim_av ~ shore, data = subset(CZ_bayes,y_rep>=0.5))

# Normality
plot(CZ_bayes_anova, 2)
# Extract the residuals
aov_residuals = residuals(object = CZ_bayes_anova)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals)
kruskal.test(sim_av ~ shore, data = subset(CZ_bayes,y_rep>=0.5))

stan_run_mat1 = stan(file = "stan/mat_gau1.stan",data = list(N = nrow(CZ_bayes),
                                                             y = CZ_bayes$mountYN,
                                                             ratio = CZ_bayes$size_ratio,
                                                             shape = CZ_bayes$Shape,
                                                             N_shore = length(unique(CZ_bayes$shore)),
                                                             N_eco = length(unique(CZ_bayes$ref_ecotype)),
                                                             shore = CZ_bayes$shore,
                                                             eco = CZ_bayes$ref_ecotype))

traceplot(stan_run_mat1, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"), inc_warmup = TRUE)
traceplot(stan_run_mat1, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"), inc_warmup = TRUE)
traceplot(stan_run_mat1, pars = c("delta","delta_sd","delta_av","delta_ta"), inc_warmup = TRUE)

print(stan_run_mat1, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"), digits_summary=3)
print(stan_run_mat1, pars = c("delta","delta_sd","delta_av","delta_ta"), digits_summary=3)
print(stan_run_mat1, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"), digits_summary=3)

stan_dens(stan_run_mat1, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"))
stan_dens(stan_run_mat1, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"))
stan_dens(stan_run_mat1, pars = c("delta","delta_sd","delta_av","delta_ta"))

stan_plot(stan_run_mat1, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"))
stan_plot(stan_run_mat1, pars = c("beta","beta_sd","beta_av","beta_ta","gamma_ta"))
stan_plot(stan_run_mat1, pars = c("delta","delta_sd","delta_av","delta_ta"))


launch_shinystan(stan_run_mat1)
print(stan_run_mat1, pars = c("y_rep"), digits_summary=3)
y_rep = extract(stan_run_mat1, pars = c('y_rep'))$y_rep
dim(y_rep)
CZ_bayes$y_rep = summary(stan_run_mat1, pars = c("y_rep"))$summary[,'mean']
write.table(CZ_bayes, "CZ_ABCD/results/CZABCD_stan1.csv", row.names = FALSE, col.names = TRUE,sep = ";")

CZ_bayes$line = 0
ggplot(data = CZ_bayes,aes(size_ratio,y_rep,col=as.factor(ref_ecotype))) +
  facet_wrap(~ as.factor(CZ_bayes$shore)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = CZ_bayes$line) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref ecotype",size="bin size",x="ln(female size) - ln(male size)",
       y="predicted probability of mounting", title = "Predicted mounting probability") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

pars_lev = summary(stan_run_mat1, pars = c('alpha','delta','beta'))$summary
pars_lev_mean = pars_lev[,'mean']

pars_sd_mean = summary(stan_run_mat1, pars = c('alpha_sd','delta_sd','beta_sd'))$summary[,'mean']
#pars_sd_mean = pars_sd[,'mean']

pars_av_mean = summary(stan_run_mat1, pars = c('alpha_av','delta_av','beta_av'))$summary[,'mean']
#pars_av_mean = pars_av[,'mean']

pars_ta = summary(stan_run_mat1, pars = c('alpha_ta','delta_ta','gamma_ta','beta_ta'))$summary
pars_ta_mean = pars_ta[,'mean']

sim_lev = pars_lev_mean[1:4][CZ_bayes$shore] + pars_lev_mean[5:6][CZ_bayes$ref_ecotype] + pars_lev_mean[7] * CZ_bayes$Shape

sim_sd = exp(pars_sd_mean[1:4][CZ_bayes$shore] + pars_sd_mean[5:6][CZ_bayes$ref_ecotype] + pars_sd_mean[7] * CZ_bayes$Shape) # exp in the stan code
CZ_bayes$sim_sd = sim_sd
hist(CZ_bayes$sim_sd[CZ_bayes$mountYN==1])
ggplot(CZ_bayes,aes(sim_sd,fill=as.factor(ref_ecotype))) +
  geom_histogram(position = "dodge",bins = 10) +
  facet_wrap(~shore) +
  scale_fill_manual(values=c("gold2","blue")) +
  labs(fill="ref ecotype",x="Trait variation in all mating pairs",
       title = "Sexual selection") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

sim_av = pars_av_mean[1:4][CZ_bayes$shore] + pars_av_mean[5:6][CZ_bayes$ref_ecotype] + pars_av_mean[7] * CZ_bayes$Shape
CZ_bayes$sim_av = sim_av
hist(CZ_bayes$sim_av,col = CZ_bayes$shore,breaks = 25)
ggplot(subset(CZ_bayes,y_rep>=0.5),aes(sim_av,fill=as.factor(ref_ecotype))) +
  geom_histogram(position = "dodge",bins = 10) +
  facet_wrap(~shore) +
  scale_fill_manual(values=c("gold2","blue")) +
  labs(fill="ref ecotype",x="Size difference of successful pairs (ln)",
       title = "Disassortative mating") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))

ggplot(subset(CZ_bayes,mountYN==1),aes(size_ratio,fill=ref_ecotype)) +
  geom_histogram(position = "dodge",bins = 10) +
  facet_wrap(~shore) +
  scale_fill_manual(values=c("gold2","blue")) +
  labs(fill="ref ecotype",x="Size difference of successful pairs (ln)",
       title = "Disassortative mating") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))


write.table(CZ_bayes, "CZ_ABCD/results/CZABCD_stan1.csv", row.names = FALSE, col.names = TRUE,sep = ";")

sim_ta = pars_ta_mean[1:4][CZ_bayes$shore] + pars_ta_mean[5:6][CZ_bayes$ref_ecotype] +
  pars_ta_mean[7] * CZ_bayes$size_ratio + pars_ta_mean[8] * CZ_bayes$Shape
library(boot)
CZ_bayes$y_sim_mean = inv.logit(sim_lev + (1 / sqrt(2 * pi * sim_sd^2)) * 
                                  exp(-0.5*((CZ_bayes$size_ratio-sim_av)/sim_sd)^2) + sim_ta)
y_sim = rbinom(nrow(CZ_bayes),size = 1, CZ_bayes$y_sim_mean)
y_sim = rbinom(nrow(CZ_bayes),size = 1, CZ_bayes$y_rep)
library(caret)
confusionMatrix(y_sim,CZ_bayes$mountYN)
library(loo)

stan_run_mat0 = stan(file = "stan/mat_gau.stan",data = list(N = nrow(CZ_bayes),
                                                            mount = CZ_bayes$mountYN,
                                                            ratio = CZ_bayes$size_ratio,
                                                            shape = CZ_bayes$Shape,
                                                            N_shore = length(unique(CZ_bayes$shore)),
                                                            shore = CZ_bayes$shore))

print(stan_run_mat0, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"), digits_summary=3)
print(stan_run_mat0, pars = c("beta","beta_sd","beta_av","beta_ta"), digits_summary=3)
traceplot(stan_run_mat0, pars = c("beta","beta_sd","beta_av","beta_ta"), inc_warmup = TRUE)
traceplot(stan_run_mat0, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"), inc_warmup = TRUE)
stan_dens(stan_run_mat0, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"))
stan_dens(stan_run_mat0, pars = c("beta","beta_sd","beta_av","beta_ta"))

pars_lev = summary(stan_run_mat0, pars = c('alpha','beta'))$summary
pars_lev_mean = pars_lev[,'mean']

pars_sd = summary(stan_run_mat0, pars = c('alpha_sd','beta_sd'))$summary
pars_sd_mean = pars_sd[,'mean']

pars_av = summary(stan_run_mat0, pars = c('alpha_av','beta_av'))$summary
pars_av_mean = pars_av[,'mean']

pars_ta = summary(stan_run_mat0, pars = c('alpha_ta','beta_ta'))$summary
pars_ta_mean = pars_ta[,'mean']

sim_lev = pars_lev_mean[1:4][CZ_bayes$shore] + pars_lev_mean[5] * CZ_bayes$Shape
sim_sd = exp(pars_sd_mean[1:4][CZ_bayes$shore] + pars_sd_mean[5] * CZ_bayes$Shape) # exp in the stan code
sim_av = pars_av_mean[1:4][CZ_bayes$shore] + pars_av_mean[5] * CZ_bayes$Shape
sim_ta = pars_ta_mean[1] * CZ_bayes$size_ratio + pars_ta_mean[2] * CZ_bayes$Shape

library(boot)
y_sim_mean = inv.logit(sim_lev + (1 / sqrt(2 * pi * sim_sd^2)) * 
                         exp(-0.5*((CZ_bayes$size_ratio-sim_av)/sim_sd)^2) + sim_ta)
CZ_bayes$y_sim = rbinom(nrow(CZ_bayes),size = 1, y_sim_mean)
summary(CZ_bayes)

library(caret)
confusionMatrix(CZ_bayes$y_sim,CZ_bayes$mountYN)

CZ_glm$line = 0
ggplot(data = CZ_glm,aes(size_ratio,y_probs,col=ref_ecotype)) +
  facet_wrap(~ as.factor(CZ_glm$shore)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = CZ_glm$line) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref ecotype",size="bin size",x="ln(female size) - ln(male size)",
       y="predicted probability of mounting", title = "Predicted mounting probability") +
  theme(legend.title = element_text(size = 11,face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 13,face = "bold",hjust = 0.5))


breaks <- c(-2,seq(-1.5,1.3,0.15),2)
bin <- cut(CZ_bayes$size_ratio,breaks)
p_m_obs <- aggregate(CZ_bayes$y_sim,by=list(bin,CZ_bayes$shore,CZ_bayes$ref_ecotype),FUN=mean)[,2:4]
bin_mean <- aggregate(CZ_bayes$size_ratio,by=list(bin,CZ_bayes$shore,CZ_bayes$ref_ecotype),FUN=mean)[,4]
CZ_bayes$count = 1
bin_sum <- aggregate(CZ_bayes$count,by=list(bin,CZ_bayes$shore,CZ_bayes$ref_ecotype),FUN=sum)[,4]
CZ_bayes_bin <- data.frame(cbind(shore=p_m_obs[,1],ref_ecotype=p_m_obs[,2],mount=p_m_obs[,3],
                                 size_ratio=bin_mean,count=bin_sum))
CZ_bayes_bin$shore = as.factor(CZ_bayes_bin$shore)
levels(CZ_bayes_bin$shore) = c("CZA","CZB","CZC","CZD")
CZ_bayes_bin$ref_ecotype = as.factor(CZ_bayes_bin$ref_ecotype)
levels(CZ_bayes_bin$ref_ecotype) = c("crab","wave")
library(ggplot2)
CZ_bayes_bin$line = 0
ggplot(data = CZ_bayes_bin,aes(size_ratio,mount,col=as.factor(ref_ecotype))) +
  geom_point(aes(size=CZ_bayes_bin$count)) +
  geom_vline(xintercept = CZ_bayes_bin$line) +
  facet_wrap(~ as.factor(CZ_bayes_bin$shore)) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref ecotype",size="bin size",x="ln(female size) - ln(male size)",
       y="proportion of mounts in each bin")

write.table(CZ_bayes, "CZ_ABCD/results/CZABCD_stan2.csv", row.names = FALSE, col.names = TRUE,
            sep = ";")




# this works
stan_code = '
data {
int<lower=0> N;
vector[N] ratio;
int<lower=0,upper=1> mount[N];
int<lower=0> N_shore;
int shore[N];
}
parameters {
vector[N_shore] beta0;
vector[N_shore] beta1;
vector<lower=0>[N_shore] sigma;
real<lower=0> sigma_beta0;
real<lower=0> sigma_beta1;

}
model {
for (i in 1:N)
mount[i] ~ bernoulli_logit(beta0[shore[i]] + (1 / sqrt(2 * pi() * sigma[shore[i]])^2) 
* exp(-0.5*((ratio[i]-beta1[shore[i]])/sigma[shore[i]])^2));

for (j in 1:N_shore) {
beta0[j] ~ normal(0, sigma_beta0);
beta1[j] ~ normal(0, sigma_beta1);
sigma[j] ~ cauchy(0, 10);
}
sigma_beta0 ~ cauchy(0, 10);
sigma_beta1 ~ cauchy(0, 10);

}
'
# Rhat
# distribution sigma and beta:
# beta 1 is the peak 
# sigma standard deviation
# beta0 shifts

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan_run = stan(data = list(N = nrow(CZ_bayes_test),
                            mount = CZ_bayes_test$mountYN,
                            ratio = CZ_bayes_test$size_ratio,
                            N_shore = length(unique(CZ_bayes_test$shore)),
                            shore = CZ_bayes_test$shore),
                model_code = stan_code)
plot(stan_run)
print(stan_run)
rstan::traceplot(stan_run)
library(rstanarm)
launch_shinystan(stan_run)
stan_dens(stan_run,pars = 'sigma')
stan_dens(stan_run,pars = 'beta1')
stan_dens(stan_run,pars = 'beta0')
pars = extract(stan_run)
pars_summ = summary(stan_run)$summary



# this works too
stan_code_mat = '
data {
  int<lower=0> N;
  vector[N] ratio;
  vector[N] shape;
  int<lower=0,upper=1> mount[N];
  int<lower=0> N_shore;
  int shore[N];
}
parameters {
  vector[N_shore] beta;
  vector[N_shore] beta_sd;
  vector[N_shore] beta_av;
  vector[N_shore] alpha;
  vector[N_shore] alpha_sd;
  vector[N_shore] alpha_av;
  vector[N_shore] alpha_ta;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta_sd;
  real<lower=0> sigma_beta_av;
  real<lower=0> sigma_alpha_sd;
  real<lower=0> sigma_alpha_av;
  real<lower=0> sigma_alpha_ta;
}

transformed parameters {
  vector[N] level;
  vector[N] aver;
  vector<lower=0>[N] stan_dev;
  vector[N] rtail;
  for (i in 1:N) {
    level[i] = alpha[shore[i]] + beta[shore[i]]*shape[i];
    stan_dev[i] = exp(alpha_sd[shore[i]] + beta_sd[shore[i]]*shape[i]);
    aver[i] = alpha_av[shore[i]] + beta_av[shore[i]]*shape[i];
    rtail[i] = alpha_ta[shore[i]];
  }
}

model {
  for (i in 1:N) {
    mount[i] ~ bernoulli_logit(level[i] + (1 / sqrt(2 * pi() * stan_dev[i]^2)) * 
      exp(-0.5*((ratio[i]-aver[i])/stan_dev[i])^2) + rtail[i]);
}
  for (j in 1:N_shore) {
    beta[j] ~ normal(0, sigma_beta);
    alpha[j] ~ normal(0, sigma_alpha);
    beta_sd[j] ~ normal(0, sigma_beta_sd);
    alpha_sd[j] ~ normal(0, sigma_alpha_sd);
    beta_av[j] ~ normal(0, sigma_beta_av);
    alpha_av[j] ~ normal(0, sigma_alpha_av);
    alpha_ta[j] ~ normal(0, sigma_alpha_ta);
  }
  sigma_beta ~ cauchy(0, 10);
  sigma_alpha ~ cauchy(0, 10);
  sigma_beta_sd ~ cauchy(0, 10);
  sigma_alpha_sd ~ cauchy(0, 10);
  sigma_beta_av ~ cauchy(0, 10);
  sigma_alpha_av ~ cauchy(0, 10);
  sigma_alpha_ta ~ cauchy(0, 10);
}
'
CZ_bayes = read.csv("CZ_ABCD/final_data/CZABCD_use_cleanup.csv",sep = ";")
CZ_bayes$ref_ecotype=as.integer(CZ_bayes$ref_ecotype) # 1 for crab and 2 for wave
CZ_bayes$shore=as.integer(CZ_bayes$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
hist(CZ_bayes$size_ratio)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan_run_mat = stan(data = list(N = nrow(CZ_bayes),
                                mount = CZ_bayes$mountYN,
                                ratio = CZ_bayes$size_ratio,
                                shape = CZ_bayes$Shape,
                                N_shore = length(unique(CZ_bayes$shore)),
                                shore = CZ_bayes$shore),
                    model_code = stan_code_mat)
summary(stan_run_mat)  # summarises parameter estimates
print(stan_run_mat)   # includes measures of fit (y_hat)
print(stan_run_mat, pars = c("alpha","alpha_ta"), digits_summary=3)
print(stan_run_mat, pars = c("beta","beta_sd","beta_av","beta_ta"), digits_summary=3)
stan_dens(stan_run_mat, pars = c("alpha","alpha_ta"))
#
#                             
# 
# 
# 
## they can be tight together but maybe sd should have a different beta because it cannot be negative
## explanatory variable can be added after the first x


## simple model
stan_code_mat = '
data {
  int<lower=0> N;
  vector[N] ratio;
  vector[N] shape;
  int<lower=0> N_shore;
  int shore[N];
  int<lower=0,upper=1> mount[N];
  
}
parameters {
  vector[N_shore] alpha;
  real beta;
  vector[N_shore] alpha_sd;
  real beta_sd;
  vector[N_shore] alpha_av;
  real beta_av;
  vector[N_shore] alpha_ta;
  real beta_ta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_alpha_sd;
  real<lower=0> sigma_beta_sd;
  real<lower=0> sigma_alpha_av;
  real<lower=0> sigma_beta_av;
  real<lower=0> sigma_alpha_ta;
  real<lower=0> sigma_beta_ta;
}

transformed parameters {
  vector[N] level;
  vector[N] aver;
  vector<lower=0>[N] stan_dev;
  vector[N] rtail;
  for (i in 1:N) {
    level[i] = alpha[shore[i]] + beta * shape[i];
    stan_dev[i] = exp(alpha_sd[shore[i]] + beta_sd * shape[i]);
    aver[i] = alpha_av[shore[i]] + beta_av * shape[i];
    rtail[i] = alpha_ta[shore[i]] + beta_ta * ratio[i];
  }
}

model {
  for (i in 1:N) {
    mount[i] ~ bernoulli_logit(level[i] + (1 / sqrt(2 * pi() * stan_dev[i]^2)) * 
      exp(-0.5*((ratio[i]-aver[i])/stan_dev[i])^2) + rtail[i]);
  }
  for (j in 1:N_shore) {
    alpha[j] ~ normal(0, sigma_alpha);
    alpha_sd[j] ~ normal(0, sigma_alpha_sd);
    alpha_av[j] ~ normal(0, sigma_alpha_av);
    alpha_ta[j] ~ normal(0, sigma_alpha_ta);
  }
  beta ~ normal(0, sigma_beta);
  beta_sd ~ normal(0, sigma_beta_sd);
  beta_av ~ normal(0, sigma_beta_av);
  beta_ta ~ normal(0, sigma_beta_ta);
  sigma_alpha ~ cauchy(0, 5);
  sigma_beta ~ cauchy(0, 5);
  sigma_alpha_sd ~ cauchy(0, 5);
  sigma_beta_sd ~ cauchy(0, 5);
  sigma_alpha_av ~ cauchy(0, 5);
  sigma_beta_av ~ cauchy(0, 5);
  sigma_alpha_ta ~ cauchy(0, 5);
  sigma_beta_ta ~ cauchy(0, 5);
}
'
CZ_glm = read.csv("CZ_ABCD/results/CZABCD_stan.csv",sep = ";")
summary(CZ_glm)

CZ_bayes = read.csv("CZ_ABCD/final_data/CZABCD_use_cleanup.csv",sep = ";")
head(CZ_bayes)
summary(CZ_bayes)
CZ_bayes$ref_ecotype=as.integer(CZ_bayes$ref_ecotype) # 1 for crab and 2 for wave
CZ_bayes$shore=as.integer(CZ_bayes$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan_run_mat = stan(data = list(N = nrow(CZ_bayes),
                                mount = CZ_bayes$mountYN,
                                ratio = CZ_bayes$size_ratio,
                                shape = CZ_bayes$Shape,
                                N_shore = length(unique(CZ_bayes$shore)),
                                shore = CZ_bayes$shore),
                    model_code = stan_code_mat)
class(stan_run_mat)
summary(stan_run_mat)  # summarises parameter estimates
print(stan_run_mat)   # includes measures of fit
print(stan_run_mat, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"), digits_summary=3)
print(stan_run_mat, pars = c("beta","beta_sd","beta_av","beta_ta"), digits_summary=3)
traceplot(stan_run_mat, pars = c("beta","beta_sd","beta_av","beta_ta"), inc_warmup = TRUE)
traceplot(stan_run_mat, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"), inc_warmup = TRUE)
stan_dens(stan_run_mat, pars = c("alpha","alpha_sd","alpha_av","alpha_ta"))
stan_dens(stan_run_mat, pars = c("beta","beta_sd","beta_av","beta_ta"))

lev_hat = head(summary(stan_run_mat, pars = 'level')$summary)

pars_lev = summary(stan_run_mat, pars = c('alpha','beta'))$summary
class(pars_lev)
head(pars_lev)
#pars_lev_mean = apply(pars_lev,2,'mean')
pars_lev_mean = pars_lev[,'mean']

pars_sd = summary(stan_run_mat, pars = c('alpha_sd','beta_sd'))$summary
pars_sd_mean = pars_sd[,'mean']

pars_av = summary(stan_run_mat, pars = c('alpha_av','beta_av'))$summary
pars_av_mean = pars_av[,'mean']

pars_ta = summary(stan_run_mat, pars = c('alpha_ta','beta_ta'))$summary
pars_ta_mean = pars_ta[,'mean']

sim_lev = pars_lev_mean[1:4][CZ_bayes$shore] + pars_lev_mean[5] * CZ_bayes$Shape
sim_sd = exp(pars_sd_mean[1:4][CZ_bayes$shore] + pars_sd_mean[5] * CZ_bayes$Shape) # exp in the stan code
sim_av = pars_av_mean[1:4][CZ_bayes$shore] + pars_av_mean[5] * CZ_bayes$Shape
sim_ta = pars_ta_mean[1:4][CZ_bayes$shore] + pars_ta_mean[5] * CZ_bayes$size_ratio

library(gtools)
y_sim_mean = inv.logit(sim_lev + (1 / sqrt(2 * pi * sim_sd^2)) * 
                         exp(-0.5*((CZ_bayes$size_ratio-sim_av)/sim_sd)^2) + sim_ta)
y_sim_mean = inv.logit(pars_lev_mean[1:4][CZ_bayes$shore] + pars_lev_mean[5] * CZ_bayes$Shape)
tail(y_sim_mean)
library(purrr)
rbernoulli()
install.packages('Rlab')
library(Rlab)
y_sim = rbern(nrow(CZ_bayes), y_sim_mean)
y_sim = rbinom(nrow(CZ_bayes),size = 1, y_sim_mean)
head(CZ_glm)
CZ_glm$y_sim = y_sim

breaks <- c(-2,seq(-1.5,1.3,0.15),2)
bin <- cut(CZ_glm$size_ratio,breaks)
p_m_obs <- aggregate(CZ_glm$y_sim,by=list(bin,CZ_glm$shore,CZ_glm$ref_ecotype),FUN=mean)[,2:4]
bin_mean <- aggregate(CZ_glm$size_ratio,by=list(bin,CZ_glm$shore,CZ_glm$ref_ecotype),FUN=mean)[,4]
CZ_glm$count = 1
bin_sum <- aggregate(CZ_glm$count,by=list(bin,CZ_glm$shore,CZ_glm$ref_ecotype),FUN=sum)[,4]
CZ_glm_bin <- data.frame(cbind(shore=p_m_obs[,1],ref_ecotype=p_m_obs[,2],mount=p_m_obs[,3],
                               size_ratio=bin_mean,count=bin_sum))
CZ_glm_bin$shore = as.factor(CZ_glm_bin$shore)
levels(CZ_glm_bin$shore) = c("CZA","CZB","CZC","CZD")
CZ_glm_bin$ref_ecotype = as.factor(CZ_glm_bin$ref_ecotype)
levels(CZ_glm_bin$ref_ecotype) = c("crab","wave")
library(ggplot2)
CZ_glm_bin$line = 0
ggplot(data = CZ_glm_bin,aes(size_ratio,mount,col=as.factor(ref_ecotype))) +
  geom_point(aes(size=CZ_glm_bin$count)) +
  geom_vline(xintercept = CZ_glm_bin$line) +
  facet_wrap(~ as.factor(CZ_glm_bin$shore)) +
  scale_color_manual(values=c("gold2","blue")) +
  labs(col="ref ecotype",size="bin size",x="ln(female size) - ln(male size)",
       y="proportion of mounts in each bin")

sim_mount = CZ_glm_bin$mount
plot(sim_mount,obs_mount)
abline(a=0, b=1)
confusionMatrix(CZ_glm$y_sim,CZ_glm$mountYN)[2:3]
head(CZ_bayes)
write.table(CZ_glm, "CZ_ABCD/results/CZABCD_stan.csv", row.names = FALSE, col.names = TRUE,sep = ";")
# pars = extract(stan_run, pars = 'beta_trt')$beta_trt
# beta_means = apply(pars,2,'mean')

# beta = apply(extract(stan_run, pars = 'beta')$beta, 2, 'mean')