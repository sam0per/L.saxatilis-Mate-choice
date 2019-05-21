rm(list = ls())

#############################
# Install required packages #
#############################
#if (!require("rstan")) install.packages("rstan")
#(.packages())

# List of packages for session
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc", "pander",
              "bbmle", "loo", "ggpubr", "cowplot", "purrr", "reshape2", "gridExtra", "grid", "arm", "parallel",
              "rstantools", "margins", "rstanarm", "projpred")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)

#load("Rpackages")
#for (p in setdiff(packages, installed.packages()[,"Package"])){
#  install.packages(p)
#}

#############################
### read and adjust data ###
#############################

#CZ_data = read.csv("data/CZ_mating_clean.csv",sep = ";")
CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
#CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
#CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
head(CZ_data)
summary(CZ_data)
range(CZ_data$shape)

######################################
# Marina's size variance test snails #
######################################
colnames(CZ_data)
CZ_data %>%
  #group_by(shore) %>%

if (CZ_data$shore == 'CZA') {
  filter(CZ_data$DistAlongPath <= 100 | CZ_data$DistAlongPath >= max(CZ_data$DistAlongPath)-30)
    #dplyr::summarise(MeanDistAlongPath = mean(DistAlongPath),
    #                 Count = n())
} else {
  filter(CZ_data$DistAlongPath <= 20 | CZ_data$DistAlongPath >= max(DistAlongPath)-20)
} %>%
  ggplot(., aes(DistAlongPath, log(length_mm))) +
  facet_wrap(~ shore, scale = 'free') +
  geom_point()

foo = CZ_data[CZ_data$shore == "CZA" & !(CZ_data$DistAlongPath > 120 & CZ_data$DistAlongPath < 300), ]
plot(foo$DistAlongPath, foo$length_mm)

shore = c("CZA", "CZB", "CZC", "CZD")
lapply(shore, function(x) {
  if (x == "CZA") {
    sample_n(CZ_data[CZ_data$shore == x & !(CZ_data$DistAlongPath > 120 | CZ_data$DistAlongPath < 300), ],
             size = 20, replace = TRUE)
    #nrow(CZ_data[CZ_data$shore == x & CZ_data$DistAlongPath < min(CZ_data$DistAlongPath) + 95, ])
  } else {
    sample_n(CZ_data[CZ_data$shore == x & !(CZ_data$DistAlongPath > 70 | CZ_data$DistAlongPath < 150), ],
             size = 20, replace = TRUE)
  }
}) %>%
  lapply(., "[", c("shore", "length_mm", "DistAlongPath", "test_sex")) %>%
  reduce(., full_join) %>%
  ggplot(., aes(DistAlongPath, log(length_mm))) +
  facet_wrap(~ shore, scale = 'free') +
  geom_point()



#############################
# stan model with only size #
#############################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)


dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio)
str(dat)
#dat$posterior_predictive = 1
#fit.prior.pred = stan(file = "scripts/CZ_mating_gaus_prior_size.stan", data = dat)

#writeLines(readLines("scripts/gaus_skew/gaus_skew_hier.stan"))
CZ_mat_stan_size = stan(file = "scripts/gaus_size/CZ_mating_gaus_size.stan", data = dat, iter = 8000, warmup = 2000, chains=2,
                        refresh=8000, control = list(adapt_delta = 0.90, max_treedepth = 15))
#gaus_all = stan(file="scripts/min_gaus_all.stan", data=dat, iter = 8000, warmup = 2000, chains=2, refresh=8000,
#                control = list(adapt_delta = 0.90, max_treedepth = 15))
#saveRDS(CZ_mat_stan_size, "tables/gaus_size/gaus_size_stan.rds")
saveRDS(CZ_mat_stan_size, "scripts/gaus_size/gaus_size_stan.rds")
CZ_mat_stan_size = readRDS("tables/gaus_size/gaus_size_stan.rds")

#############################
# plot post distr gaus size #
#############################
gaus_size_pars = c("level","scale","preference","choosiness","asymmetry")
gaus_size_parfig = c("preference","choosiness","asymmetry")

list_of_draws <- rstan::extract(CZ_mat_stan_size)
print(names(list_of_draws))
parfig = lapply(gaus_size_parfig, function(x) {
  ggplot() +
    geom_density(aes(list_of_draws[[x]]), fill='red', col='black') +
    labs(x="", title = x) +
    theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
})
names(parfig) = gaus_size_parfig
lev_sca_dens = ggplot() +
  geom_density(aes(inv.logit(list_of_draws$level + list_of_draws$scale)), fill='red', col='black') +
  labs(x="", title = "level + scale") +
  theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
parfig$lev_sca = lev_sca_dens

#level + scale * exp(-0.5 * ((x - preference)/choosiness)^2) + asymmetry * x
opt = function(asymmetry, scale, preference, x, choosiness){
  # derivative
  asymmetry - (scale * (x - preference) * exp(-(0.5 * (preference - x)^2)/choosiness^2))/choosiness^2
}
pars = list(asymmetry = list_of_draws$asymmetry, scale = list_of_draws$scale, preference = list_of_draws$preference,
            choosiness = list_of_draws$choosiness)
length(pars)
# find the root of the derivative
#str(xmin <- uniroot(opt, c(0, 1), tol = 0.0001, asymmetry = 2, scale = 7.2, preference = 0.1, choosiness = 1))

opt_draws = sapply(1:12000, function(z){
  uniroot(opt, c(0, 1), tol = 0.0001, asymmetry = pars[['asymmetry']][z], scale = pars[['scale']][z],
          preference = pars[['preference']][z], choosiness = pars[['choosiness']][z])$root
})

list_of_draws$optimum = opt_draws

opt_dens = ggplot() +
  geom_density(aes(list_of_draws$optimum), fill='red', col='black') +
  labs(x="", title = "optimum") +
  theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
mean(list_of_draws$optimum)
parfig$optimum = opt_dens


pdf("figures/gaus_size/gaus_size_pars_dens.pdf",width = 10, height = 7)
#do.call(ggarrange, parfig)
ggarrange(parfig$lev_sca, parfig$preference, parfig$choosiness, parfig$asymmetry, parfig$optimum)
dev.off()
#pdf("figures/gaus_size/gaus_size_opt_dens.pdf")
#ggarrange(parfig$optimum, widths = 3, heights = 2)
#stan_dens(CZ_mat_stan_size, pars = gaus_size_pars)
#dev.off()
pdf("figures/gaus_size/gaus_size_pars_plot.pdf")
stan_plot(CZ_mat_stan_size, pars = gaus_size_pars)
dev.off()


#########################
# pars values gaus size #
#########################
CZ_size_params = round(summary(CZ_mat_stan_size, pars = gaus_size_pars, probs=c(0.025, 0.975))$summary,2)
lev_sca_val = CZ_size_params['level',] + CZ_size_params['scale',]
lev_sca_val['mean'] = inv.logit(lev_sca_val['mean'])
list_of_draws$lev_sca = list_of_draws$level + list_of_draws$scale
lev_sca_val['se_mean'] = (sqrt(sd(list_of_draws$level)^2 + sd(list_of_draws$scale)^2))/
  sqrt(length(list_of_draws$lev_sca))
lev_sca_val['sd'] = sqrt(sd(list_of_draws$level)^2 + sd(list_of_draws$scale)^2)
xx <- seq(min(list_of_draws$lev_sca), max(list_of_draws$lev_sca), by = 0.001)
plot(density(list_of_draws$lev_sca))
dd = density(list_of_draws$lev_sca)
summary(dd$y)
fx <- splinefun(dd$x, dd$y) # interpolating function
pxx <- pmax(0, fx(xx)) # normalize so prob >0
# sample from the "empirical" distribution
samp <- sample(xx, 1e5, replace = TRUE, prob = pxx)
# and take sample quantiles
lev_sca_val['2.5%'] = inv.logit(quantile(samp, c(0.025, 0.975)))[1]
lev_sca_val['97.5%'] = inv.logit(quantile(samp, c(0.025, 0.975)))[2]

lev_sca_val['Rhat'] = 1
lev_sca_val = round(lev_sca_val, 2)
CZ_size_params = rbind(lev_sca_val, CZ_size_params)
rownames(CZ_size_params)[1] = "level + scale"
CZ_size_params

opt_val = data.frame(mean = mean(list_of_draws$optimum),
                     se_mean = sqrt(var(list_of_draws$optimum)/length(list_of_draws$optimum)),
                     sd = sd(list_of_draws$optimum),
                     `2.5%` = NA,
                     `97.5%` = NA,
                     n_eff = NA,
                     Rhat = NA)
xx <- seq(min(list_of_draws$optimum), max(list_of_draws$optimum), by = 0.001)
plot(density(list_of_draws$optimum))
dd = density(list_of_draws$optimum)
summary(dd$y)
fx <- splinefun(dd$x, dd$y) # interpolating function
pxx <- pmax(0, fx(xx)) # normalize so prob >0
# sample from the "empirical" distribution
samp <- sample(xx, 1e5, replace = TRUE, prob = pxx)
# and take sample quantiles
colnames(opt_val) = colnames(CZ_size_params)
opt_val['2.5%'] = quantile(samp, c(0.025, 0.975))[1]
opt_val['97.5%'] = quantile(samp, c(0.025, 0.975))[2]


rownames(opt_val) = "optimum"
opt_val[1:5] = round(opt_val[1:5], 2)
CZ_size_params = rbind(opt_val, CZ_size_params)
CZ_size_params

inv.logit((-6.94) + 7.17 * exp(-0.5 * (((-0.5) - 0.09)/0.97)^2) + 1.52 * (-0.5)) -
  inv.logit((-6.94) + 7.17 * exp(-0.5 * (((0.5) - 0.09)/0.97)^2) + 1.52 * (0.5))

inv.logit((-6.94) + 7.17 * exp((-0.5) * ((1 - 0.09)/0.97)^2) + 1.52 * 1)
inv.logit((-6.94) + 7.17 * exp((-0.5) * (((-0.31) - 0.09)/0.97)^2) + 1.52 * (-0.31))


CZ_size_params = rownames_to_column(as.data.frame(CZ_size_params), var="params")
write.table(CZ_size_params, "tables/gaus_size/gaus_size_params.csv", row.names = FALSE, col.names = TRUE,sep = ";")
#write.table(CZ_size_params[c(-2:-3),], "tables/gaus_size/gaus_size_tr_params.csv", row.names = FALSE, col.names = TRUE,sep = ";")
#sapply(CZ_size_params$params, function(x) rnorm(n = 1000, CZ_size_params$mean[x], sd = CZ_size_params$sd[x])) %>%
#  apply(., 2, hist)

#y_rep = data.frame(y_rep=summary(CZ_mat_stan_size, pars = c("y_rep"))$summary[,'mean'])
y_rep = data.frame(round(summary(CZ_mat_stan_size, pars = c("y_rep"))$summary, 2))
head(y_rep)
mean(y_rep[,'mean'])
y_rep = rownames_to_column(y_rep, var="rep")
write.table(y_rep, "tables/gaus_size/gaus_size_mount_yrep.csv", row.names = FALSE, col.names = TRUE,sep = ";")

CZ_data$y_rep = summary(CZ_mat_stan_size, pars = c("y_rep"))$summary[,'mean']
CZ_data$y_rep_se = summary(CZ_mat_stan_size, pars = c("y_rep"))$summary[,'se_mean']


y_hat = round(summary(CZ_mat_stan_size, pars = c("y_hat"))$summary, 2)
y_hat = rownames_to_column(as.data.frame(y_hat), var="hat")
write.table(y_hat, "tables/gaus_size/gaus_size_mount_yhat.csv", row.names = FALSE, col.names = TRUE,sep = ";")

CZ_logit = summary(CZ_mat_stan_size, pars = c("y_hat"))$summary[,'mean']
max(inv.logit(CZ_logit))
CZ_uci = summary(CZ_mat_stan_size, pars = c("y_hat"))$summary[,'97.5%']
CZ_lci = summary(CZ_mat_stan_size, pars = c("y_hat"))$summary[,'2.5%']
CZ_data$preds = inv.logit(CZ_logit) %>% round(3)
#CZ_data$preds = CZ_logit %>% round(3)
CZ_data$uci_preds = inv.logit(CZ_uci) %>% round(3)
CZ_data$lci_preds = inv.logit(CZ_lci) %>% round(3)
CZ_data$y_preds = rbinom(n = nrow(CZ_data),size = 1,prob = CZ_data$preds)
#CZ_data$y_preds = rbinom(n = nrow(CZ_data),size = 1,prob = inv.logit(CZ_data$preds))

write.table(CZ_data, "tables/gaus_size/gaus_size_mat.csv", row.names = FALSE, col.names = TRUE, sep = ";")



###################################
# plot observs and preds for size #
###################################
y = CZ_data$mountYNcontact
y_rep = rstan::extract(CZ_mat_stan_size, pars = 'y_rep', permuted = TRUE)$y_rep
dim(y_rep)
ppc_bars(y,y_rep)
range(CZ_data$size_ratio)
breaks = c(-2,seq(-1.5,1.5,0.1),2)
bin = cut(CZ_data$size_ratio,breaks)
pdf("figures/gaus_size/gaus_size_ppc_bars_grouped.pdf")
ppc_bars_grouped(y, y_rep, bin, prob = 0, freq = FALSE)
dev.off()

#CZ_data = read.csv("tables/CZs_gaus_size_mat.csv", sep = ";")
CZ_data$bin = cut(CZ_data$size_ratio,breaks)
CZ_data_bin =
  CZ_data %>%
  group_by(bin) %>%
  dplyr::summarise(mount = mean(mountYNcontact),
                   uci_mount = CI(mountYNcontact)['upper'],
                   lci_mount = CI(mountYNcontact)['lower'],
                   mean_ratio = mean(size_ratio),
                   y_rep = mean(y_rep),
                   preds_mount = mean(y_preds)) %>%
  mutate(lci_mount = replace(lci_mount, which(lci_mount<0), 0))
summary(CZ_data_bin)

pdf("figures/gaus_size/gaus_size_preds.pdf", width=8, height=7)
ggplot(data = CZ_data) +
  geom_vline(xintercept = 0) +
  geom_ribbon(aes(x = size_ratio,ymin = lci_preds, ymax = uci_preds), fill = "orange", alpha=0.3) +
  geom_errorbar(data = CZ_data_bin, aes(x = mean_ratio, ymin = lci_mount, ymax = uci_mount),alpha = 0.2) +
  scale_colour_manual(values=c("blue","orange2")) +
  geom_line(aes(size_ratio,preds,col="predictions")) +
  geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount, col="observations")) +
  labs(size="bin size",x="female - male size (ln)",
       y="probability of mounting",col="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.title = element_text(size = 11,face = "bold"), legend.position = c(0.05, 0.85),
        axis.title = element_text(face = "bold", size = 13)) +
  grids(linetype = "dashed")
dev.off()

pdf("figures/gaus_size/gaus_size_ppc.pdf")
ggplot(CZ_data_bin,aes(mount,y_rep)) +
  geom_abline(slope = 1, alpha=0.5) +
  geom_point() +
  labs(x='observed mount proportion', y='predicted probability') +
  theme(axis.title = element_text(face = "bold", size = 13))
dev.off()


###################################
##### cline analysis for size #####
###################################
CZA <- read.csv("../2.mating/CZA/final_data/CZA_use_cleanup.csv",sep = ";")
CZA$shore <- "CZA"
CZB <- read.csv("../2.mating/CZB/final_data/CZB_use_cleanup.csv",sep = ";")
CZB$shore <- "CZB"
CZC <- read.csv("../2.mating/CZC/final_data/CZC_use_cleanup.csv",sep = ";")
CZC$shore <- "CZC"
CZD <- read.csv("../2.mating/CZD/final_data/CZD_use_cleanup.csv",sep = ";")
CZD$shore <- "CZD"
CZ <- merge(CZB,CZA,all = TRUE)
CZ <- merge(CZ,CZC,all = TRUE)
CZ <- merge(CZ,CZD,all = TRUE)

CZall = CZ[!is.na(CZ$Shape) & CZ$log_female>0.5,]
write.table(CZall, "data/CZ_all_mating_clean.csv", row.names = FALSE, col.names = TRUE, sep = ";")

CZ_all = read.csv("data/CZ_all_mating_clean.csv", sep = ",")
# CZ_data = read.csv("tables/CZ_size_mating.csv", sep = ";")
# identical(sort(CZ_all$size_ratio), sort(CZ_data$size_ratio))
CZ_all %>%
  group_by(snail_ID, ref_ecotype) %>%
  dplyr::summarise(count=n()) %>%
  filter(count > 2)





cline_2c3s <- function(phen,position,sex,cl,cr,lwl,lwr,crab,wave,zs_c,zs_w,sc,sh,sw){
  wl = exp(lwl)
  wr = exp(lwr)
  # sc = exp(lsc)
  # sh = exp(lsh)
  # sw = exp(lsw)
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]
  s_xl <- sqrt(sc^2 + 4*p_xl*(1-p_xl)*sh^2 + (p_xl^2)*(sw^2-sc^2))

  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))

  # combined cline
  z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  s_x[z_x < z_xl] <- s_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  return(minusll)
}

cline_2c2s <- function(phen,position,sex,cl,cr,wl,wr,crab,wave,zs_c,zs_w,sc,sw){
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]
  s_xl <- sqrt(sc^2 + (p_xl)*(sw^2-sc^2))

  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  s_x <- sqrt(sc^2 + (p_x)*(sw^2-sc^2))

  # combined cline
  z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  s_x[z_x < z_xl] <- s_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=TRUE))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  return(minusll)
}

cline_2c1s <- function(phen, position, sex, cl, cr, lwl, lwr, crab, wave, zs_c, zs_w, ls_x){
  wl = exp(lwl)
  wr = exp(lwr)
  # s_x = exp(ls_x)
  s_x = ls_x
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]

  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]

  # combined cline
  z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  return(minusll)
}

rm(theta.init ,mle.cline.2c3s, p)
mle.cline.2c3s = list(CZA=NULL, CZB=NULL, CZC=NULL, CZD=NULL)

for (p in levels(CZ_all$shore)) {
  if (p=='CZA'){
    plot(CZ_all$LCmeanDist[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=130,cr=280,lwl=3,lwr=2.3,crab=-2.1,wave=-1.9,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.3,sw=0.2)
    mle.cline.2c3s$CZA = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$LCmeanDist[CZ_all$shore==p],
                                        sex=CZ_all$test_sex[CZ_all$shore==p]))
  }
  else if (p=='CZB'){
    plot(CZ_all$LCmeanDist[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=70,cr=150,lwl=1.6,lwr=3.9,crab=-2.5,wave=-1.5,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.3,sw=0.2)
    mle.cline.2c3s$CZB = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$LCmeanDist[CZ_all$shore==p],
                                        sex=CZ_all$test_sex[CZ_all$shore==p]))
  }
  else if (p=='CZC'){
    plot(CZ_all$LCmeanDist[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=50,cr=125,lwl=1.5,lwr=3,crab=-2.5,wave=-1.5,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.3,sw=0.2)
    mle.cline.2c3s$CZC = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$LCmeanDist[CZ_all$shore==p],
                                        sex=CZ_all$test_sex[CZ_all$shore==p]))
  }
  else {
    plot(CZ_all$LCmeanDist[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=80,cr=175,lwl=1.6,lwr=1.6,crab=-2.5,wave=-1.5,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.3,sw=0.2)
    mle.cline.2c3s$CZD = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$LCmeanDist[CZ_all$shore==p],
                                        sex=CZ_all$test_sex[CZ_all$shore==p]))
  }
}

lapply(mle.cline.2c3s, summary)

# ggarrange(tableGrob(round(summary(mle.cline.2c3s$CZA)@coef, 2)))
# ggarrange(tableGrob(round(summary(mle.cline.2c3s$CZB)@coef, 2)))
# ggarrange(tableGrob(round(summary(mle.cline.2c3s$CZC)@coef, 2)))
# ggarrange(tableGrob(round(summary(mle.cline.2c3s$CZD)@coef, 2)))

(CZ_cline_params = sapply(mle.cline.2c3s, function(x) round(coef(x), 2)))
(CZ_cline_params = abs(CZ_cline_params))

s_centre = function(sc, sh, sw) {
  # sc = exp(lsc)
  # sh = exp(lsh)
  # sw = exp(lsw)
  sqrt(sc^2 + sh^2 + 0.25*(sw^2-sc^2))
}
s_pars = sapply(c("CZA", "CZB", "CZC", "CZD"), function(x) {
  sh = s_centre(sc = CZ_cline_params["sc", x], sh = CZ_cline_params["sh", x], sw = CZ_cline_params["sw", x])
  # sc = exp(CZ_cline_params["lsc", x])
  # sw = exp(CZ_cline_params["lsw", x])
  sc = CZ_cline_params["sc", x]
  sw = CZ_cline_params["sw", x]
  return(round(rbind(sc, sh, sw), 2))
})

cline_pars = row.names(CZ_cline_params)
CZ_cline_params = rbind(CZ_cline_params[1:nrow(CZ_cline_params[c(-1:-3), ]), ], s_pars)
row.names(CZ_cline_params) = cline_pars
write.table(CZ_cline_params, "tables/clines/CZ_cline_params.csv", row.names = TRUE, col.names = TRUE, sep = ",")

# s_centre = function(sc, sh, sw) {sqrt(sc^2 + sh^2 + 0.25*(sw^2-sc^2))}
# s_centre(sc = exp(CZ_cline_params["lsc", "CZA"]), sh = exp(CZ_cline_params["lsh", "CZA"]),
#          sw = exp(CZ_cline_params["lsw", "CZA"]))

CZ_cline_params = rownames_to_column(as.data.frame(CZ_cline_params), var="params")
CZ_cline_params[(5:8),-1] = abs(CZ_cline_params[(5:8),-1])
CZ_cline_params[(-1:-2),-1] = round(exp(CZ_cline_params[(-1:-2),-1]), 2)
write.table(CZ_cline_params, "tables/clines/CZ_cline_params.csv", row.names = FALSE, col.names = TRUE, sep = ";")

CZ_cline_params[,1] = c("cl","cr","wl","wr","crab_mm","wave_mm","zs_c_mm","zs_w_mm","sc","sh","sw")
# s_centre <- sqrt(sc^2 + sh^2 + 0.25*(sw^2-sc^2))
round(sqrt(CZ_cline_params[9,-1]), 2)
round(sqrt(CZ_cline_params[11,-1]), 2)
CZ_cline_params[10,-1] = round(sapply(CZ_cline_params[-1], function(x) sqrt(x[9]^2 + 4*0.5*(1-0.5)*x[10]^2 +
                                                                          (0.5^2)*(x[11]^2-x[9]^2))), 2)
CZ_cline_params[11,-1] = round(sapply(CZ_cline_params[-1], function(x) sqrt(x[9]^2 + (x[11]^2-x[9]^2))), 2)

write.table(CZ_cline_params, "tables/clines/CZ_cline_params_mm.csv", row.names = FALSE, col.names = TRUE, sep = ";")



(CZ_cline_se = sapply(mle.cline.2c3s, function(x) round(sqrt(diag(vcov(x))), 2)))
CZ_cline_se = rownames_to_column(as.data.frame(CZ_cline_se), var="params")
write.table(CZ_cline_se, "tables/clines/CZ_cline_se.csv", row.names = FALSE, col.names = TRUE, sep = ";")

sapply(mle.cline.2c3s, function(x) summary(x))
sapply(mle.cline.2c3s, function(x) AIC(x))




CZ_all %>% group_by(shore, sex) %>% dplyr::summarise(mean_size=mean(log(length_mm)))
colnames(CZ_all)
CZ_all %>% group_by(shore) %>% dplyr::summarise(min=min(DistAlongPath), max=max(DistAlongPath))
breaks = c(0,seq(10,300,50),400)
CZ_all$bin = cut(CZ_all$DistAlongPath,breaks)
CZ_all_bin = CZ_all %>% group_by(shore, bin, sex) %>%
  dplyr::summarise(mean_size=mean(log(length_mm)))

CZ_all_bin$diff = c(0,round(diff(CZ_all_bin$mean_size, lag = 1),2))

pdf("figures/test_size_mfdiff_cline.pdf")
CZ_all_bin %>%
  dplyr::filter(row_number() %% 2 == 0) %>%
  ggplot() +
    geom_point(aes(bin,diff,col=shore)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(x="distance intervals", y="male - female size (ln)") +
    theme(axis.title = element_text(face = "bold"))
dev.off()

###################################
# apply size model to field distr #
###################################

# generate size distributions for each sex and ecotype using cline parameters
#s_xl <- sqrt(sc^2 + 4*0.5*(1-0.5)*sh^2 + (0.5^2)*(sw^2-sc^2))    # standard deviation at the zone centre
CZ_cline_params_mm = read.csv("tables/clines/CZ_cline_params_mm.csv", sep = ";")
ggarrange(tableGrob(CZ_cline_params_mm, rows = NULL))


CZ_cline_params = read.csv("tables/clines/CZ_cline_params.csv", sep = ";")
male_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[5]), sd = abs(x[9]))),2)
#apply(male_c, 2, hist)
female_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[5]+x[7]),
                                                               sd = abs(x[9]))),2)
#apply(female_c, 2, hist)


male_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs((x[5]+x[6])/2), sd = abs(x[10]))),2)
#apply(male_h, 2, hist)
female_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs((x[5]+x[6]+x[7]+x[8])/2),
                                                               sd = abs(x[10]))),2)
#apply(female_h, 2, hist)

male_w = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[6]), sd = abs(x[11]))),2)
#apply(male_w, 2, hist)
female_w = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[6]+x[8]),
                                                               sd = abs(x[11]))),2)
#apply(female_w, 2, hist)

# pair each female with every male within habitat and contact zone and compute mounting success YN
CZ_size_params = read.csv("tables/gaus_size/gaus_size_params.csv", sep = ";")


fem=list(crab=data.frame(female_c), hybrid=data.frame(female_h), wave=data.frame(female_w))
mal=list(crab=data.frame(male_c), hybrid=data.frame(male_h), wave=data.frame(male_w))

sim_mat = function(female, male) {
  bar = list()
  YN = data.frame()
  for (f in seq_along(female)) {
    success=FALSE
    i=1
    fem = female[f]
    while (!success) {
      m = sample(male, 1, replace = FALSE)
      p = inv.logit(CZ_size_params$mean[CZ_size_params$params=='level'] +
                      CZ_size_params$mean[CZ_size_params$params=='scale'] *
                      exp(-0.5 * (((fem - m) - CZ_size_params$mean[CZ_size_params$params=='preference'])/
                                    CZ_size_params$mean[CZ_size_params$params=='choosiness'])^2) +
                      CZ_size_params$mean[CZ_size_params$params=='asymmetry'] * (fem - m))
      s = rbinom(n = 1, size = 1, prob = p)
      YN[i,'male'] = m
      YN[i,'female'] = fem
      YN[i,'mountYN'] = s
      success = (s > 0)
      i = i + 1
    }
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}


res = lapply(names(fem), function(x) {
  sapply(names(fem$crab), function(y) sim_mat(female = fem[[x]][[y]], male = mal[[x]][[y]]))
})

shore = c("CZA", "CZB", "CZC", "CZD")
ecotype = c("crab", "hybrid", "wave")

eco_CZ = lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) do.call(rbind, res[[x]][,y]))
})

# save the simulated datasets of mate choice
lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) write.table(eco_CZ[[x]][[y]],
                                                   paste0("tables/gaus_size/sims/", shore[y], "_", ecotype[x], "_sim_YN.csv"),
                                                   row.names = FALSE, col.names = TRUE, sep = ";"))
})


####################
# sexual selection #
####################
rm(list = ls())
shore = c("CZA", "CZB", "CZC", "CZD")
ecotype = c("crab", "hybrid", "wave")


se = function(x) sqrt(var(x)/length(x))
lci = function(x) mean(x) - 1.96*se(x)
uci = function(x) mean(x) + 1.96*se(x)

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom", legend.justification = "center",
                                     legend.direction = "vertical"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

CZ_sim_dat = list()
CZ_sim_ss = list()
sim_fit_ci = list()


for (eco in ecotype) {
  if (eco=="wave") {
    sim_dat = list.files("tables/gaus_size/sims", pattern = paste0(eco, "_sim_YN"), full.names = TRUE)
    for (f in 1:length(sim_dat)) {
      CZ_sim_dat[[f]] = read.csv(sim_dat[f], sep = ";")
      CZ_sim_dat[[f]][,'male2'] = CZ_sim_dat[[f]][,'male']^2
      CZ_sim_ss[[f]] = glm(mountYN~male+male2, family = binomial(link = "logit"), data = CZ_sim_dat[[f]])
      sim_fit_ci[[f]] = tibble(fit=predict(CZ_sim_ss[[f]], se.fit=TRUE)$fit,
                               fit_se=predict(CZ_sim_ss[[f]], se.fit=TRUE)$se.fit,
                               fit_lci=exp(fit-1.96*fit_se)/(1+exp(fit-1.96*fit_se)),
                               fit_uci=exp(fit+1.96*fit_se)/(1+exp(fit+1.96*fit_se)),
                               male=CZ_sim_dat[[f]][,'male'],
                               female=CZ_sim_dat[[f]][,'female'],
                               mountYN=CZ_sim_dat[[f]][,'mountYN'],
                               ecotype=eco)
    }
    lapply(1:length(sim_fit_ci), function(x) write.table(sim_fit_ci[[x]], paste0("tables/gaus_size/sex_sel/", shore[x], "_",
                                                                                 eco, "_sim_ss_fit.csv"),
                                                       row.names = FALSE, col.names = TRUE, sep = ";"))
    sim_ss_coef = lapply(1:length(CZ_sim_ss), function(x) ggarrange(tableGrob(round(summary(CZ_sim_ss[[x]])$coef, 2))))
    lapply(1:length(sim_ss_coef), function(x) ggsave(filename = paste0("tables/gaus_size/sex_sel/", shore[x], "_",
                                                                       eco, "_sim_size_ss_coef.png"),
                                                     plot = sim_ss_coef[[x]]))


    CZ_sim_min = sapply(CZ_sim_dat, function(x) min(x[['male']]))
    CZ_sim_max = sapply(CZ_sim_dat, function(x) max(x[['male']]))
    CZ_sim_brk = lapply(1:length(sim_dat), function(x) c(CZ_sim_min[x]-0.1,
                                                         seq(CZ_sim_min[x]+0.1, CZ_sim_max[x]-0.1, 0.1),
                                                         CZ_sim_max[x]+0.1))
    CZ_sim_bin = lapply(1:length(sim_dat), function(x) cut(CZ_sim_dat[[x]][['male']], CZ_sim_brk[[x]]))
    CZ_sim_bin_dat = lapply(1:length(sim_dat), function(x) tibble(male=CZ_sim_dat[[x]][['male']],
                                                                  mountYN=CZ_sim_dat[[x]][['mountYN']],
                                                                  #preds=CZ_sim_dat[[x]][['preds']],
                                                                  bin=CZ_sim_bin[[x]])) %>%
      map(., function(x) group_by(x, bin)) %>%
      map(., function(x) summarise_all(x, .funs = c('mean', 'lci', 'uci'))) %>%
      map(., function(x) {
        mutate(x, mountYN_lci = replace(mountYN_lci, which(mountYN_lci<0), 0),
               mountYN_uci = replace(mountYN_uci, which(mountYN_uci>1), 1))
      })


    CZ_sim_ss_plot = lapply(1:length(sim_dat), function(x) {
      ggplot() +
        geom_ribbon(data = sim_fit_ci[[x]],
                    aes(x = male, ymin = fit_lci, ymax = fit_uci), fill = "orange", alpha=0.3) +
        geom_errorbar(data = CZ_sim_bin_dat[[x]], aes(x = male_mean, ymin = mountYN_lci, ymax = mountYN_uci),
                      alpha = 0.2) +
        geom_point(data = CZ_sim_bin_dat[[x]], aes(male_mean, mountYN_mean, col='binned data')) +
        #geom_point(data = CZ_sim_bin_dat[[x]], aes(male_mean, preds_mean), col='red') +
        geom_line(data = CZ_sim_ss[[x]], aes(male, fitted(CZ_sim_ss[[x]]), col='predictions')) +
        #geom_smooth(data = CZ_sim_ss[[x]], aes(male, fitted(CZ_sim_ss[[x]])), method = 'glm', formula=y~poly(x,2))
        labs(x="male size (ln)", y="predicted prob", col="", title=paste(shore[x], eco, sep = " ")) +
        scale_color_manual(values = c("blue", "orange"))
    })
    sim_ss_plot = do.call(grid_arrange_shared_legend, CZ_sim_ss_plot)
    ggsave(filename = paste0("figures/gaus_size/sex_sel/CZ_", eco, "_sim_size_ss.png"), plot = sim_ss_plot)
  }
}
dev.off()

rm(list = ls())
shore = c("CZA", "CZB", "CZC", "CZD")
ecotype = c("crab", "hybrid", "wave")

CZ_eco_chr = list(CZA=NULL, CZB=NULL, CZC=NULL, CZD=NULL)
for (s in seq_along(shore)) {
  for (e in seq_along(ecotype)) {
    CZ_eco_chr[[s]][e] = list.files("tables/gaus_size/sex_sel", pattern = paste(shore[s], ecotype[e], "sim", "ss", "fit",sep = "_"),
                                    full.names = TRUE)
  }
}

CZ_eco_df = lapply(names(CZ_eco_chr), function(x) {
  lapply(seq_along(ecotype), function(y) read.csv(CZ_eco_chr[[x]][y], sep = ";"))
}) %>%
  lapply(., function(x) reduce(x, full_join)) %>%
  map(., function(x) mutate(x, fit=round(inv.logit(fit),2)*100))


palette = colorRampPalette(colors=c("black", "blue", "red"))

n_fills = sapply(seq_along(CZ_eco_df), function(x) as.numeric(length(levels(factor(CZ_eco_df[[x]]$fit)))))
fills = sapply(seq_along(n_fills), function(x) palette(n_fills[x]))

CZ_eco_male_mx = map(CZ_eco_df, function(x) group_by(x, factor(ecotype), factor(fit))) %>%
  map(., function(x) dplyr::summarise(x, male_mean=mean(male))) %>%
  lapply(., setNames, nm = c("ecotype", "fit", "male_mean")) %>%
  map(., arrange, desc(as.integer(fit))) %>%
  map(., function(x) group_by(x, ecotype)) %>%
  lapply(., function(x) dplyr::summarise(x, male_mx = male_mean[which.max(x[[as.numeric(2)]])]))

#lapply(CZ_eco_male_mean, function(x) which.max(x[['fit']]))
#View(CZ_eco_male_mean[[1]])
#arrange(CZ_eco_male_mean[[1]], desc(fit))
#sample_n(CZ_eco_male_mean[[1]], size = 1)
#sapply(seq_along(CZ_eco_df), function(x) CZ_eco_df[[x]]$male[max(CZ_eco_df[[x]]$fit)])

CZ_eco_ss = lapply(seq_along(CZ_eco_df), function(x) {
  ggplot() +
    geom_dotplot(data = CZ_eco_df[[x]],
                 aes(x = factor(ecotype), y = male, fill=factor(fit), col=factor(fit)),
                 binwidth = 0.015, binaxis = "y", stackdir = "center",
                 stackratio = 0.5, dotsize = 1.5) +
    stat_summary(data = CZ_eco_df[[x]], aes(x = factor(ecotype), y = male),
                 fun.y = mean, geom = "point", size = 5, color="grey") +
    geom_point(data = CZ_eco_male_mx[[x]], aes(x=ecotype, y=male_mx), col="white", size=5, shape=1, stroke=2) +
    scale_fill_manual(values = fills[[x]]) +
    scale_color_manual(values = fills[[x]]) +
    labs(x="habitat", y="male size (ln)", title=shore[x]) +
    rremove("legend") +
    rremove("xlab")
})

lapply(seq_along(CZ_eco_ss), function(x) ggsave(filename = paste0("figures/gaus_size/sex_sel/", shore[x],
                                                                   "_sim_eco_ss.png"),
                                                 plot = CZ_eco_ss[[x]]))


CZ_male_mean = map(CZ_eco_df, function(x) group_by(x, factor(ecotype))) %>%
  map(., function(x) dplyr::summarise(x, male_mean=mean(male))) %>%
  lapply(., setNames, nm = c("ecotype", "male_mean"))

CZ_male_diff = map(seq_along(shore), function(x) full_join(CZ_male_mean[[x]], CZ_eco_male_mx[[x]])) %>%
  map(., function(x) mutate(x, male_diff=round(male_mean-male_mx, 2)))

lapply(seq_along(shore), function(x) write.table(CZ_male_diff[[x]], paste0("tables/gaus_size/sex_sel/", shore[x],
                                                                           "_sim_ss_mdiff.csv"),
                                                 row.names = FALSE, col.names = TRUE, sep = ";"))


######################
# assortative mating #
######################
CZ_eco_am = lapply(seq_along(CZ_eco_df), function(x) {
  ggplot() +
    geom_point(data = subset(CZ_eco_df[[x]], mountYN==1), aes(x = female, y = male, col=factor(ecotype))) +
    facet_wrap(~factor(ecotype)) +
    labs(x="female size (ln)", y="male size (ln)", title=shore[x]) +
    grids(linetype = "dashed") +
    #annotate("text", x=0, y=4, label=paste0('atop(bold("p_value is ',p_val,'"))'))
    rremove("legend")
})

lapply(seq_along(CZ_eco_am), function(x) ggsave(filename = paste0("figures/gaus_size/ass_mat/", shore[x],
                                                                  "_sim_eco_am.png"),
                                                plot = CZ_eco_am[[x]]))

CZ_eco_amr = lapply(CZ_eco_df, subset, mountYN==1) %>%
  map(., group_by, factor(ecotype)) %>%
  map(., function(x) dplyr::summarise(x, p_cor = round(cor(female, male, method = "pearson"), 2))) %>%
  lapply(., setNames, nm = c("ecotype", "p_cor"))

lapply(seq_along(shore), function(x) write.table(CZ_eco_amr[[x]], paste0("tables/gaus_size/ass_mat/", shore[x],
                                                                           "_sim_am_pcor.csv"),
                                                 row.names = FALSE, col.names = TRUE, sep = ";"))

lapply(CZ_eco_df, subset, mountYN==0) %>%
  map(., group_by, factor(ecotype)) %>%
  map(., function(x) dplyr::summarise(x, p_cor = round(cor(female, male, method = "pearson"), 2))) %>%
  lapply(., setNames, nm = c("ecotype", "p_cor"))





#######################################################
# stan model with shore, ecotype, sex and shape (all) #
#######################################################
rm(list = ls())
##########################
# terminal huluvu server #
##########################
library(rstan)
options(mc.cores = parallel::detectCores(logical = FALSE) - 15)
CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
CZ_data$test_sex=as.integer(CZ_data$test_sex) # 1 for female, 2 for male
#CZ_data$ref_eco_sex[CZ_data$test_sex=="female" & CZ_data$ref_ecotype=="crab"] = "crabm"
#CZ_data$ref_eco_sex[CZ_data$test_sex=="female" & CZ_data$ref_ecotype=="wave"] = "wavem"
#CZ_data$ref_eco_sex[CZ_data$test_sex=="male" & CZ_data$ref_ecotype=="crab"] = "crabf"
#CZ_data$ref_eco_sex[CZ_data$test_sex=="male" & CZ_data$ref_ecotype=="wave"] = "wavef"
#CZ_data$ref_eco_sex = as.factor(CZ_data$ref_eco_sex)
#write.table(CZ_data, "data/CZ_all_mating_clean.csv", row.names = FALSE, col.names = TRUE,sep = ";")
CZ_data$ref_eco_sex = as.integer(CZ_data$ref_eco_sex) # 1 for crab female, 2 for crab male, 3 for wave female, 4 for wave male
summary(CZ_data)
str(CZ_data)
dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio, shape = CZ_data$shape,
           N_shore = length(unique(CZ_data$shore)), N_eco = length(unique(CZ_data$ref_eco_sex)),
           N_sex = length(unique(CZ_data$test_sex)),
           shore = CZ_data$shore, ref = CZ_data$ref_eco_sex, test = CZ_data$test_sex)
rstan_options(auto_write = TRUE)


#gaus_all = stan(file="scripts/min_gaus_all.stan", data=dat, iter = 6000, warmup = 1500, chains=4,
#                refresh=6000, control = list(adapt_delta=0.90, max_treedepth = 15))
#gaus_all = stan(file="scripts/min_gaus_all.stan", data=dat, iter = 6000, warmup = 2000, chains=2, refresh=6000)
#gaus_all = stan(file="scripts/min_gaus_all.stan", data=dat, iter = 6000, warmup = 2000, chains=4, refresh=6000,
#                control = list(adapt_delta = 0.90, max_treedepth = 15))

gaus_all = stan(file="scripts/min_gaus_all.stan", data=dat, iter = 8000, warmup = 2000, chains=2, refresh=8000,
                control = list(adapt_delta = 0.90, max_treedepth = 15))
saveRDS(gaus_all, "scripts/gaus_sex/gaus_sex.rds")



#library(shinystan)
#launch_shinystan(gaus_all)
#rm(list = ls()[!ls() %in% c("CZ_data", "dat")])
#gaus_all = stan(file="scripts/min_gaus_all.stan", data=dat, iter = 6000, warmup = 1500, chains=2,refresh=6000)
#(CZ_all_pars = gaus_all@model_pars)
#library(parallel)
#n_pars = 1:30
#ls_pars = split(n_pars, ceiling(seq_along(n_pars)/(length(n_pars)/5)))
#all_pars_tbl = mclapply(ls_pars, function(x) round(summary(gaus_all, pars = CZ_all_pars[x],probs=c(0.025, 0.975))$summary,2),mc.cores = 5)

(CZ_all_pars = gaus_all@model_pars)
library(parallel)
n_pars = 5:22
ls_pars = split(n_pars, ceiling(seq_along(n_pars)/(length(n_pars)/3)))
ls_pars$`4` = 2:4
ls_pars$`5` = 1
ls_pars$`6` = seq(23,26)
all_pars_tbl = mclapply(ls_pars, function(x) round(summary(gaus_all, pars = CZ_all_pars[x],probs=c(0.025, 0.975))$summary,2),mc.cores = 6)
#all_pars_tbl = lapply(ls_pars, function(x) round(summary(gaus_all, pars = CZ_all_pars[x],probs=c(0.025, 0.975))$summary,2))
#rm(list = setdiff(ls(), c('CZ_data','dat')))

library(tibble)
all_pars_tbl = lapply(all_pars_tbl, function(x) rownames_to_column(as.data.frame(x), var="params"))
gaus_terms = c("preference", "choosiness", "asymmetry", "scale", "level", "hyperpars")
#lapply(seq_along(all_pars_tbl), function(x) write.table(all_pars_tbl[[x]],
#                                                        paste0("tables/gaus_all/gaus_all_", gaus_terms[x],
#                                                               "_coef.csv"),
#                                                        row.names = FALSE, col.names = TRUE, sep = ";"))
lapply(seq_along(all_pars_tbl), function(x) write.table(all_pars_tbl[[x]],
                                                        paste0("tables/gaus_sex/gaus_sex_", gaus_terms[x],
                                                               "_coef.csv"),
                                                        row.names = FALSE, col.names = TRUE, sep = ";"))

ls_pars$`6` = NULL
gaus_all_dens = mclapply(ls_pars, function(x) stan_dens(gaus_all, pars = CZ_all_pars[x]),mc.cores = 5)

gaus_all_plot = mclapply(ls_pars, function(x) stan_plot(gaus_all, pars = CZ_all_pars[x]),mc.cores = 4)

gaus_coef = c("mus", "sigmas", "gammas", "lambdas", "level")
#lapply(seq_along(gaus_all_dens), function(x) ggsave(filename = paste0("figures/gaus_all/gaus_all_",
#                                                                      gaus_coef[x], "_dens.png"),
#                                                    plot = gaus_all_dens[[x]]))
lapply(seq_along(gaus_all_dens), function(x) ggsave(filename = paste0("figures/gaus_sex/gaus_sex_",
                                                                      gaus_coef[x], "_dens.png"),
                                                    plot = gaus_all_dens[[x]]))
#lapply(seq_along(gaus_all_plot), function(x) ggsave(filename = paste0("figures/gaus_all/gaus_all_",
#                                                                      gaus_coef[x], "_plot.png"),
#                                                    plot = gaus_all_plot[[x]]))
lapply(seq_along(gaus_all_plot), function(x) ggsave(filename = paste0("figures/gaus_sex/gaus_sex_",
                                                                      gaus_coef[x], "_plot.png"),
                                                    plot = gaus_all_plot[[x]]))



#######################
# load gaus all model #
#######################
CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")

gaus_all = readRDS("tables/gaus_all/gaus_all_stanfit.rds")
gaus_all = readRDS("scripts/gaus_sex/gaus_sex.rds")
#round(summary(gaus_all, pars = "lambda1",probs=c(0.025, 0.975))$summary,2)

############################
# plot post distr gaus all #
############################
#gaus_size_pars = c("level","scale","preference","choosiness","asymmetry")
#gaus_size_parfig = c("preference","choosiness","asymmetry")
gaus_sex@model_pars

# plot and compute correlation matrix
mcmcpairs <- function( posterior , n=500 , cex=0.7 , pch=16 , adj=1 , ...) {
  panel.dens <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- density(x,adj=adj)
    y <- h$y
    y <- y/max(y)
    abline( v=0 , col="gray" , lwd=0.5 )
    lines( h$x , y )
  }
  panel.2d <- function( x , y , ... ) {
    i <- sample( 1:length(x) , size=n )
    abline( v=0 , col="gray" , lwd=0.5 )
    abline( h=0 , col="gray" , lwd=0.5 )
    dcols <- densCols( x[i] , y[i] )
    points( x[i] , y[i] , col=dcols , ... )
  }
  panel.cor <- function( x , y , ... ) {
    k <- cor( x , y )
    cx <- sum(range(x))/2
    cy <- sum(range(y))/2
    text( cx , cy , round(k,2) , cex=2*exp(abs(k))/exp(1) )
  }
  pairs( posterior , cex=cex , pch=pch , upper.panel=panel.2d , lower.panel=panel.cor , diag.panel=panel.dens , ... )
}
post <- rstan::extract(gaus_sex)
all_pars = do.call(cbind, lapply(1:24, function(x) post[[x]]))
dim(all_pars)
mcmcpairs(posterior = all_pars[,1:10])
mcmcpairs(posterior = all_pars[,11:20])
mcmcpairs(posterior = all_pars[,c(1:4, 11:14)])
mcmcpairs(posterior = all_pars[,c(5:8, 15:18)])
mcmcpairs(posterior = all_pars[,c(1:4, 15:18)])
mcmcpairs(posterior = all_pars[,21:33])
mcmcpairs(posterior = all_pars[,34:46])
mcmcpairs(posterior = all_pars[,47:59])
round(cor(all_pars), 2)
dim(cor(all_pars))

dim(post$level)
hist(apply(post$level, 2, mean))
hist(apply(post$scale, 2, mean))

cor(apply(post$level, 2, mean), apply(post$scale, 2, mean))
cor(apply(post$level, 2, mean), apply(post$preference, 2, mean))
cor(apply(post$scale, 2, mean), apply(post$preference, 2, mean))

list_of_draws <- rstan::extract(gaus_all)
print(names(list_of_draws))
parfig = lapply(gaus_size_parfig, function(x) {
  ggplot() +
    geom_density(aes(list_of_draws[[x]]), fill='red', col='black') +
    labs(x="", title = x) +
    theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
})
names(parfig) = gaus_size_parfig
lev_sca_dens = ggplot() +
  geom_density(aes(inv.logit(list_of_draws$level + list_of_draws$scale)), fill='red', col='black') +
  labs(x="", title = "level + scale") +
  theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
parfig$lev_sca = lev_sca_dens

#level + scale * exp(-0.5 * ((x - preference)/choosiness)^2) + asymmetry * x
opt = function(asymmetry, scale, preference, x, choosiness){
  # derivative
  asymmetry - (scale * (x - preference) * exp(-(0.5 * (preference - x)^2)/choosiness^2))/choosiness^2
}
pars = list(asymmetry = list_of_draws$asymmetry, scale = list_of_draws$scale, preference = list_of_draws$preference,
            choosiness = list_of_draws$choosiness)
length(pars)
# find the root of the derivative
#str(xmin <- uniroot(opt, c(0, 1), tol = 0.0001, asymmetry = 2, scale = 7.2, preference = 0.1, choosiness = 1))

opt_draws = sapply(1:16000, function(z){
  uniroot(opt, c(0, 1), tol = 0.0001, asymmetry = pars[['asymmetry']][z], scale = pars[['scale']][z],
          preference = pars[['preference']][z], choosiness = pars[['choosiness']][z])$root
})

list_of_draws$optimum = opt_draws

opt_dens = ggplot() +
  geom_density(aes(list_of_draws$optimum), fill='red', col='black') +
  labs(x="", title = "optimum") +
  theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
mean(list_of_draws$optimum)
parfig$optimum = opt_dens


pdf("figures/gaus_size/gaus_size_pars_dens.pdf",width = 10, height = 7)
#do.call(ggarrange, parfig)
ggarrange(parfig$lev_sca, parfig$preference, parfig$choosiness, parfig$asymmetry, parfig$optimum)
dev.off()
pdf("figures/gaus_size/gaus_size_opt_dens.pdf")
ggarrange(parfig$optimum, widths = 3, heights = 2)
#stan_dens(CZ_mat_stan_size, pars = gaus_size_pars)
dev.off()
pdf("figures/gaus_size/gaus_size_pars_plot.pdf")
stan_plot(CZ_mat_stan_size, pars = gaus_size_pars)
dev.off()


#########################
# pars values gaus size #
#########################


#level + scale * exp(-0.5 * ((x - preference)/choosiness)^2) + asymmetry * x
opt = function(asymmetry, scale, preference, x, choosiness){
  asymmetry - (scale * (x - preference) * exp(-(0.5 * (preference - x)^2)/choosiness^2))/choosiness^2
}
pars = list(asymmetry = c(1.5, 2, 1.7), scale = c(7.2, 7, 7.1), preference = c(0.1, 0.3, 0.2),
            choosiness = c(1, 1.5, 1.2))

str(xmin <- uniroot(opt, c(0, 1), tol = 0.0001, asymmetry = 2, scale = 7.2, preference = 0.1, choosiness = 1))

sapply(1:3, function(z){
  uniroot(opt, c(0, 2), tol = 0.0001, asymmetry = pars[['asymmetry']][z], scale = pars[['scale']][z],
          preference = pars[['preference']][z], choosiness = pars[['choosiness']][z])$root
})

#round(summary(gaus_all, pars = "level",probs=c(0.025, 0.975))$summary, 2)
#stan_plot(gaus_all, pars = "alpha1")
y_rep = data.frame(round(summary(gaus_all, pars = c("y_rep"))$summary,2))
y_rep = rownames_to_column(y_rep, var="rep")
head(y_rep)
write.table(y_rep, "tables/gaus_all/CZ_all_mount_yrep.csv", row.names = FALSE, col.names = TRUE,sep = ";")

CZ_data$y_rep = summary(gaus_all, pars = c("y_rep"))$summary[,'mean']
CZ_data$y_rep_se = summary(gaus_all, pars = c("y_rep"))$summary[,'se_mean']


y_hat = round(summary(gaus_all, pars = c("y_hat"))$summary, 2)
y_hat = rownames_to_column(as.data.frame(y_hat), var="hat")
write.table(y_hat, "tables/gaus_all/CZ_all_mount_yhat.csv", row.names = FALSE, col.names = TRUE,sep = ";")

CZ_logit = summary(gaus_all, pars = c("y_hat"))$summary[,'mean']
CZ_uci = summary(gaus_all, pars = c("y_hat"))$summary[,'97.5%']
CZ_lci = summary(gaus_all, pars = c("y_hat"))$summary[,'2.5%']
CZ_data$preds = inv.logit(CZ_logit) %>% round(3)
CZ_data$uci_preds = inv.logit(CZ_uci) %>% round(3)
CZ_data$lci_preds = inv.logit(CZ_lci) %>% round(3)
CZ_data$y_preds = rbinom(n = nrow(CZ_data),size = 1,prob = CZ_data$preds)

write.table(CZ_data, "tables/gaus_all/CZs_gaus_all_mat.csv", row.names = FALSE, col.names = TRUE, sep = ";")

###################################
# plot observs and preds for size #
###################################
CZ_data = read.csv("tables/gaus_all/CZs_gaus_all_mat.csv",sep = ";")
y = CZ_data$mountYN
y_rep = rstan::extract(gaus_all, pars = 'y_rep', permuted = TRUE)$y_rep
dim(y_rep)
ppc_bars(y,y_rep)
range(CZ_data$size_ratio)
breaks = c(-2,seq(-1.5,1.3,0.2),2)
bin = cut(CZ_data$size_ratio,breaks)
pdf("figures/gaus_all/gaus_all_ppc_bars_grouped.pdf")
ppc_bars_grouped(y, y_rep, bin, prob = 0, freq = FALSE)
dev.off()



CZ_data$bin = cut(CZ_data$size_ratio,breaks)
CZ_data_bin =
  CZ_data %>%
  group_by(bin, shore, ref_ecotype) %>%
  dplyr::summarise(mount = mean(mountYN),
                   uci_mount = CI(mountYN)['upper'],
                   lci_mount = CI(mountYN)['lower'],
                   mean_ratio = mean(size_ratio),
                   y_rep = mean(y_rep),
                   preds_mount = mean(y_preds),
                   mean_pred = mean(preds)) %>%
  mutate(lci_mount = replace(lci_mount, which(lci_mount<0), 0)) %>%
  mutate(uci_mount = replace(uci_mount, which(uci_mount>1), 1))
head(CZ_data_bin)
summary(CZ_data_bin)



pdf("figures/gaus_all/gaus_all_preds.pdf", width=8, height=7)
ggplot(data = CZ_data) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ shore) +
  geom_ribbon(aes(x = size_ratio,ymin = lci_preds, ymax = uci_preds), fill = "orange", alpha=0.3) +
  geom_errorbar(data = CZ_data_bin, aes(x = mean_ratio, ymin = lci_mount, ymax = uci_mount),alpha = 0.2) +
  #scale_colour_manual(values=c("blue", "black", "red")) +
  geom_line(data = CZ_data_bin, aes(mean_ratio, mean_pred,col=ref_ecotype)) +
  geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount, col=ref_ecotype)) +
  labs(size="bin size",x="female - male size (ln)",
       y="probability of mounting",col="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.title = element_text(size = 11,face = "bold"), legend.position = c(0.02, 0.95),
        axis.title = element_text(face = "bold", size = 13)) +
  grids(linetype = "dashed")
dev.off()

pdf("figures/gaus_all/gaus_all_ppc.pdf")
ggplot(CZ_data_bin,aes(mount,y_rep)) +
  geom_abline(slope = 1, alpha=0.5) +
  geom_point() +
  labs(x='observed mount proportion', y='predicted probability') +
  theme(axis.title = element_text(face = "bold", size = 13))
dev.off()


########################################################
###### approximate leave-one-out cross-validation ######
# http://mc-stan.org/loo/articles/loo2-with-rstan.html #
########################################################
CZ_mat_stan_size = readRDS("tables/gaus_size/gaus_size_stan.rds")
gaus_sex = readRDS("scripts/gaus_sex/gaus_sex.rds")
# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(CZ_mat_stan_size, merge_chains = FALSE)
#log_lik_2 <- extract_log_lik(gaus_all, merge_chains = FALSE)
log_lik_2 <- extract_log_lik(gaus_sex, merge_chains = FALSE)

# provide relative effective sample sizes
r_eff <- relative_eff(exp(log_lik_1))
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2, save_psis = TRUE)
print(loo_1)
plot(loo_1)
r_eff_2 <- relative_eff(exp(log_lik_2))
loo_2 <- loo(log_lik_2, r_eff = r_eff_2, cores = 2, save_psis = TRUE)
print(loo_2)
gaus_sex@model_pars
CZ_mat_stan_size@model_pars
plot(loo_2, label_points = TRUE)
comp <- compare(loo_1, loo_2)
print(comp)

# Marginal posterior predictive checks
y_rep1 = rstan::extract(CZ_mat_stan_size, pars = 'y_rep', permuted = TRUE)$y_rep
ppc_loo_pit_overlay(
  y = CZ_data$mountYN,
  yrep = y_rep1,
  lw = weights(loo_1$psis_object)
)

# The excessive number of values close to 1 indicates that the model is under-dispersed compared to the data,
# and we should consider a model that allows for greater dispersion.
ppc_loo_pit_overlay(
  y = CZ_data$mountYN,
  yrep = y_rep,
  lw = weights(loo_2$psis_object)
)

#loo_1 <- loo(log_lik_1, k_threshold=0.7, cores = 2)
#loo_2 <- loo(log_lik_2, k_threshold=0.7, cores = 2)
loo_list <- list(loo_1, loo_2)
loo_model_weights(loo_list)
loo_model_weights(loo_list, method = "pseudobma")
loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)


###################################
# apply full model to field distr #
###################################

# generate size distributions for each sex and ecotype using cline parameters
#s_xl <- sqrt(sc^2 + 4*0.5*(1-0.5)*sh^2 + (0.5^2)*(sw^2-sc^2))    # standard deviation at the zone centre


CZ_cline_params = read.csv("tables/clines/CZ_cline_params.csv", sep = ";")
male_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[5]), sd = abs(x[9]))),2)
#apply(male_c, 2, hist)
female_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[5]+x[7]),
                                                               sd = abs(x[9]))),2)
#apply(female_c, 2, hist)


male_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs((x[5]+x[6])/2), sd = abs(x[10]))),2)
#apply(male_h, 2, hist)
female_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs((x[5]+x[6]+x[7]+x[8])/2),
                                                               sd = abs(x[10]))),2)
#apply(female_h, 2, hist)

male_w = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[6]), sd = abs(x[11]))),2)
#apply(male_w, 2, hist)
female_w = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[6]+x[8]),
                                                               sd = abs(x[11]))),2)
#apply(female_w, 2, hist)

# pair each female with every male within habitat and contact zone and compute mounting success YN
CZ_size_params = read.csv("tables/gaus_size/gaus_size_params.csv", sep = ";")


fem=list(crab=data.frame(female_c), hybrid=data.frame(female_h), wave=data.frame(female_w))
mal=list(crab=data.frame(male_c), hybrid=data.frame(male_h), wave=data.frame(male_w))

sim_mat = function(female, male) {
  bar = list()
  YN = data.frame()
  for (f in seq_along(female)) {
    success=FALSE
    i=1
    fem = female[f]
    while (!success) {
      m = sample(male, 1, replace = FALSE)
      p = inv.logit(CZ_size_params$mean[CZ_size_params$params=='level'] +
                      CZ_size_params$mean[CZ_size_params$params=='scale'] *
                      exp(-0.5 * (((fem - m) - CZ_size_params$mean[CZ_size_params$params=='preference'])/
                                    CZ_size_params$mean[CZ_size_params$params=='choosiness'])^2) +
                      CZ_size_params$mean[CZ_size_params$params=='asymmetry'] * (fem - m))
      s = rbinom(n = 1, size = 1, prob = p)
      YN[i,'male'] = m
      YN[i,'female'] = fem
      YN[i,'mountYN'] = s
      success = (s > 0)
      i = i + 1
    }
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}


res = lapply(names(fem), function(x) {
  sapply(names(fem$crab), function(y) sim_mat(female = fem[[x]][[y]], male = mal[[x]][[y]]))
})

shore = c("CZA", "CZB", "CZC", "CZD")
ecotype = c("crab", "hybrid", "wave")

eco_CZ = lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) do.call(rbind, res[[x]][,y]))
})

# save the simulated datasets of mate choice
lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) write.table(eco_CZ[[x]][[y]],
                                                   paste0("tables/gaus_size/sims/", shore[y], "_", ecotype[x], "_sim_YN.csv"),
                                                   row.names = FALSE, col.names = TRUE, sep = ";"))
})

##################################
##################################
######  SKEW NORMAL MODEL  #######
##################################
library(shinystan)
#skew_hier = readRDS("models/gaus_skew/gaus_skew_hier.rds")
#skew_hier = readRDS("models/gaus_skew/gaus_skew_hier_d.rds")
skew_hier = readRDS("models/gaus_skew/gaus_skew_hier_BCDG.rds")
launch_shinystan(skew_hier)

# hier_draws$preference[1:10, 1:10]
# hist(hier_draws$preference[, 3707])
# hist(apply(hier_draws$preference, 1, mean))

skew_hier@model_pars
(hier_pars = skew_hier@model_pars[1:4])
# library(parallel)
(hier_pars_tbl = round(summary(skew_hier, pars = hier_pars, probs=c(0.025, 0.975))$summary,2))

# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
hier_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)
# hier_matrix = model.matrix(mountYNcontact ~ shore, data = CZ_data)
head(hier_matrix)
row.names(hier_pars_tbl)[grepl("c_coeff", row.names(hier_pars_tbl))] = paste0("c_", colnames(hier_matrix))
row.names(hier_pars_tbl)[grepl("d_coeff", row.names(hier_pars_tbl))] = paste0("d_", colnames(hier_matrix))
row.names(hier_pars_tbl)[grepl("g_coeff", row.names(hier_pars_tbl))] = paste0("g_", colnames(hier_matrix))
# row.names(hier_pars_tbl)[!grepl("c_|d_", row.names(hier_pars_tbl))] = c("b", "d", "g")
row.names(hier_pars_tbl)[!grepl("c_|d_|g_", row.names(hier_pars_tbl))] = c("b")

# library(tibble)
(hier_pars_df = rownames_to_column(as.data.frame(hier_pars_tbl), var="parameter"))

# write.table(hier_pars_df, paste0("tables/gaus_skew/gaus_skew_hier_coef.csv"),
#             row.names = FALSE, col.names = TRUE, sep = ";")
write.table(hier_pars_df, paste0("tables/gaus_skew/gaus_skew_hier_CDG_coef.csv"),
            row.names = FALSE, col.names = TRUE, sep = ";")


hier_draws = rstan::extract(skew_hier)
# names(hier_draws)[1:3] = row.names(hier_pars_tbl)[!grepl("c_", row.names(hier_pars_tbl))]
# hier_parfig = row.names(hier_pars_tbl)[!grepl("c_", row.names(hier_pars_tbl))]
# parfig_title = c("b", "d", "g")
# parfig = lapply(seq_along(hier_parfig), function(x) {
#   ggplot() +
#     geom_density(aes(hier_draws[[x]]), fill='red', col='black') +
#     labs(x="", title = parfig_title[x]) +
#     theme(axis.title = element_text(face = "bold", size = 14),
#           plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# })
parfig = ggplot() +
  geom_density(aes(hier_draws$scale), fill='red', col='black') +
  labs(x="", title = "b") +
  theme(axis.title = element_text(face = "bold", size = 14),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))


hier_hyp = hier_pars[-1]
hier_coeff = row.names(hier_pars_tbl)[grepl("c_|d_|g_", row.names(hier_pars_tbl))]
coeff_list = split(hier_coeff, ceiling(seq_along(hier_coeff)/length(colnames(hier_matrix))))
names(coeff_list) = hier_hyp
names(hier_draws)
coeff_list[["c_coeff"]]
seq_along(names(coeff_list))
head(hier_draws[["c_coeff"]])

coeff_draws = lapply(seq_along(hier_hyp) + 1, function(x) {
  hier_draws[[x]]
})


# colnames(coeff_draws[[3]]) = coeff_list[[3]]
rn_coeff_draws = lapply(seq_along(hier_hyp), function(x) {
  colnames(coeff_draws[[x]]) = coeff_list[[x]]
  coeff_draws[[x]]
})
head(rn_coeff_draws[[3]])
head(rn_coeff_draws[[1]][, 7])


# names(hier_draws$c_coeff) = hier_c
coeff_list[[1]][1]
# head(hier_draws$c_coeff[,"c_shoreCZB"])
coeff_parfig = lapply(seq_along(hier_hyp), function(x) {
  lapply(seq_along(colnames(hier_matrix)), function(y) {
    ggplot() +
      geom_density(aes(rn_coeff_draws[[x]][, y]), fill='red', col='black') +
      labs(x="", title = coeff_list[[x]][y]) +
      theme(axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  })
})
coeff_parfig[[2]][[3]]
lapply(seq_along(hier_hyp), function(x) {
  lapply(seq_along(colnames(hier_matrix)), function(y) {
    ggsave(filename = paste0("figures/gaus_skew/gaus_skew_hier_CDG_", coeff_list[[x]][y], "_dens.png"),
           plot = coeff_parfig[[x]][[y]])
  })
})
ggsave(filename = paste0("figures/gaus_skew/gaus_skew_hier_CDG_b_dens.png"),
       plot = parfig)

posterior = as.array(skew_hier)
dim(posterior)
dimnames(posterior)$parameters[1:25] = c("b", hier_coeff)
# stan_plot(skew_hier, pars = "scale", point_est = "mean") +
#   scale_y_discrete(labels = c(`scale` = 'b'))
mcmc_intervals(posterior, pars = hier_coeff[grepl("c_", x = hier_coeff)], point_est = "mean")
mcmc_intervals(posterior, pars = hier_coeff[grepl("d_", x = hier_coeff)], point_est = "mean")
mcmc_intervals(posterior, pars = hier_coeff[grepl("g_", x = hier_coeff)], point_est = "mean")
prx_hyp = c("c_", "d_", "g_")

parfig_plot = lapply(prx_hyp, function(hyp) {
  mcmc_intervals(posterior, pars = hier_coeff[grepl(hyp, x = hier_coeff)], point_est = "mean",
                 prob_outer = 0.95)
})
lapply(seq_along(prx_hyp), function(x) {
  ggsave(filename = paste0("figures/gaus_skew/gaus_skew_hier_CDG_", prx_hyp[x], "plot.png"),
         plot = parfig_plot[[x]])
})

color_scheme_set("mix-blue-red")
parfig_trace = lapply(prx_hyp, function(hyp) {
  mcmc_trace(posterior, pars = hier_coeff[grepl(hyp, x = hier_coeff)],
             facet_args = list(ncol = 1, strip.position = "left"))
})
lapply(seq_along(prx_hyp), function(x) {
  ggsave(filename = paste0("figures/gaus_skew/gaus_skew_hier_CDG_", prx_hyp[x], "trace.png"),
         plot = parfig_trace[[x]])
})
#all_parfig = append(parfig, c_parfig)

curve(dgamma(x, shape = 5, rate = 3), from = 0, to = 5)
# curve(dnorm(x, mean = -0.1, sd = 0.1), from = -0.5, to = 0.5)
# plot(dlnorm(1:100, meanlog = 0, sdlog = 0.5, log = FALSE))
#
# x = rlnorm(500,1,.6)
# grid = seq(0,8,.01)
#
# plot(grid,dlnorm(grid,0,4),type="l",xlab="x",ylab="f(x)")
#
# plot(grid,dbeta(grid,shape1 = 6,shape2 = 10),type="l",xlab="x",ylab="f(x)")
#
# plot(grid,dgamma(grid,shape = 5,rate = 2),type="l",xlab="x",ylab="f(x)")

grid = seq(0, 8, .01)
prior = data.frame(x=grid, y=dgamma(grid, shape = 5, rate = 2))
parfig[[3]] = parfig[[3]] +
  geom_line(data = prior, aes(x, y), linetype="dotdash")

round(summary(skew_hier, pars = "c_sigma", probs=c(0.025, 0.975))$summary,2)
c_grid = seq(-1, 1, .01)
c_prior = data.frame(x=c_grid, y=dnorm(c_grid, mean = 0, sd = 0.2))
prior_parfig = lapply(seq_along(c_parfig), function(x) {
  c_parfig[[x]] + geom_line(data = c_prior, aes(x, y), linetype="dotdash")
})
all_parfig = append(parfig, prior_parfig)

hier_all = c(hier_parfig, hier_c)
names(all_parfig) = hier_all

lapply(seq_along(hier_all), function(x) ggsave(filename = paste0("figures/gaus_skew/gaus_skew_hier_",
                                                                 hier_all[x], "_dens.png"),
                                                    plot = all_parfig[[x]]))

###################################
## projective variable selection ##
###################################
hier_draws = as.matrix(skew_hier)
head(hier_draws[, 1:6])
head(hier_draws[, 16:21])

sigma = hier_draws[, "c_sigma"]
c_coeff = hier_draws[, grepl("c_coeff", colnames(hier_draws))]
dim(c_coeff)
a = c_coeff[, 1]
b = c_coeff[, -1]
head(c_coeff)
c_hyp = hier_draws[, grepl("preference", colnames(hier_draws))]
dim(c_hyp)
y = apply(c_hyp, 1, mean)
hist(y)

y = c_hyp[1, ]
hist(y)
head(y)
hist(c_hyp_mean)
CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex:shape, data = CZ_data)
head(CZ_matrix)
x = CZ_matrix[, -1]
head(x)

# prediction with the reference model
predfun <- function(xt) t( b %*% t(xt) + a)
# initialize the reference model object. notice here that the first argument z
# denotes the features of the reference model, and x the features from which
# we want to select from
ref <- init_refmodel(x, y, gaussian(), predfun=predfun, dis=sigma)

cvs <- cv_varsel(ref)
varsel_plot(cvs, stats=c('elpd','rmse'))

sessionInfo()
