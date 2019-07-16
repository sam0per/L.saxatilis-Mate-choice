rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc", "pander",
              "purrr", "reshape2", "gridExtra", "grid", "arm", "parallel", "optparse", "pracma", "ggpubr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="input data", metavar="character"),
  make_option(c("-s", "--stanfile"), type="character", default=NULL,
              help="model written in Stan", metavar="character"),
  make_option(c("-o", "--output"), type = "character", default = "output",
              help = "prefix for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$stanfile)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input data and stan file).\n", call.=FALSE)
}

CZ_data = read.csv(opt$data, sep = ";")
pref_out = opt$output
# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
# pref_out = "gaus_skew/SKEW/gaus_skew"

############################
# stan model with skewness #
############################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE) - 2)


dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio)
# str(dat)
# dat$posterior_predictive = 1
# fit.prior.pred = stan(file = "scripts/CZ_mating_gaus_prior_size.stan", data = dat)

gaus_skew = rstan::stan(file = "../L.saxatilis-Mate-choice/scripts/gaus_skew/gaus_skew.stan",
                        data = dat, iter = 8000, warmup = 2000,
                        chains=4, refresh=8000,
                        control = list(adapt_delta = 0.90, max_treedepth = 15))
saveRDS(gaus_skew, "models/gaus_skew/SKEW/gaus_skew.rds")

#############################
# plot post distr gaus skew #
#############################
# gaus_skew = readRDS("models/gaus_skew/gaus_skew.rds")
gaus_skew@model_pars
gaus_size_pars = c("b","c","d","alpha")
# gaus_size_pars = c("b0","b1","c","d","alpha")
#gaus_size_parfig = c("preference","choosiness","asymmetry")

list_of_draws <- rstan::extract(gaus_skew)
names(list_of_draws) = c("b","c","d","alpha", "y_hat", "log_lik", "y_rep", "lp__")
# names(list_of_draws) = c("b0","b1","c","d","alpha", "y_hat", "log_lik", "y_rep", "lp__")
names(list_of_draws)
parfig = lapply(gaus_size_pars, function(x) {
  ggplot() +
    geom_density(aes(x = list_of_draws[[x]], y = ..scaled..), fill='red', col='black') +
    labs(x="", title = x) +
    theme(axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold",size = 15))
})
names(parfig) = gaus_size_pars


# find derivative of skew-normal in Wolphram|Alpha typing
# d[l + s * exp(-0.5 * ((x - p) / c)^2) * (1 + ERF(a * ((x - p) / (1.414214 * c)))) ,x]

opt = function(b, c, alpha, d, x){
  # derivative
  (b * exp(-(0.5 * (c - x)^2)/d^2) *
     (0.797884 * alpha * d * exp(-(0.5 * alpha^2 * (c - x)^2)/d^2) -
        (x - c) * (erf((0.707107 * alpha * (x - c))/d) + 1)))/d^2
}
pars = list(b = list_of_draws$b, c = list_of_draws$c, alpha = list_of_draws$alpha,
            d = list_of_draws$d)

# find the root of the derivative
#str(xmin <- uniroot(opt, c(0, 1), tol = 0.0001, scale = 0.44, preference = -0.07, asymmetry = 1.35, choosiness = 0.67))

opt_draws = sapply(1:24000, function(z){
  uniroot(opt, c(0, 1), tol = 0.0001, alpha = pars[['alpha']][z], b = pars[['b']][z],
          c = pars[['c']][z], d = pars[['d']][z])$root
})

list_of_draws$opt = opt_draws

opt_dens = ggplot() +
  geom_density(aes(x = list_of_draws$opt, y = ..scaled..), fill='red', col='black') +
  labs(x="", title = "optimum") +
  theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
mean(list_of_draws$opt)
parfig$opt = opt_dens


pdf("figures/gaus_skew/SKEW/gaus_skew_pars_dens.pdf",width = 10, height = 7)
#do.call(ggarrange, parfig)
ggarrange(parfig$b, parfig$c, parfig$d, parfig$alpha, parfig$opt)
dev.off()


#########################
# pars values gaus skew #
#########################
# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
# CZ_data = read.csv("data/CZ_all_mating_clean_copy.csv", sep = ";")
gaus_skew@model_pars
gaus_skew_pars = c("b_par","c_par","d_par","alpha_par")
# gaus_skew_pars = c("level","scale","preference","choosiness","asymmetry")

(skew_params = round(summary(gaus_skew, pars = gaus_skew_pars, probs=c(0.025, 0.975))$summary,2))
rownames(skew_params) = c("b","c","d","alpha")
# class(gaus_skew)

y_hat = summary(gaus_skew, pars = c("y_hat"))$summary[,'mean']
# CZ_erf = summary(gaus_skew, pars = c("err_val"))$summary[,'mean']
# y_prob = y_hat * (1 + erf(CZ_erf))
# plot(CZ_data$size_ratio, y_hat)

# err_val = skew_params['asymmetry', 'mean'] *
#   (CZ_data$size_ratio - skew_params['preference', 'mean']) / (1.414214 * skew_params['choosiness', 'mean'])
#
# y_hat = skew_params['level', 'mean'] + skew_params['scale', 'mean'] *
#   exp(-0.5 * ((CZ_data$size_ratio - skew_params['preference', 'mean']) /
#                 skew_params['choosiness', 'mean'])^2)
#
# y_prob = y_hat * (1 + erf(err_val))
#
# plot(CZ_data$size_ratio, y_prob)

(skew_params = rownames_to_column(as.data.frame(skew_params), var="parameter"))

opt_val = data.frame(parameter = 'optimum',
                     mean = mean(list_of_draws$opt),
                     se_mean = sqrt(var(list_of_draws$opt)/length(list_of_draws$opt)),
                     sd = sd(list_of_draws$opt),
                     `2.5%` = NA,
                     `97.5%` = NA,
                     n_eff = NA,
                     Rhat = NA)
xx <- seq(min(list_of_draws$opt), max(list_of_draws$opt), by = 0.001)
# plot(density(list_of_draws$opt))
dd = density(list_of_draws$opt)
summary(dd$y)
fx <- splinefun(dd$x, dd$y) # interpolating function
pxx <- pmax(0, fx(xx)) # normalize so prob >0
# sample from the "empirical" distribution
samp <- sample(xx, 1e5, replace = TRUE, prob = pxx)
# and take sample quantiles
colnames(opt_val) = colnames(skew_params)
opt_val['2.5%'] = quantile(samp, c(0.025, 0.975))[1]
opt_val['97.5%'] = quantile(samp, c(0.025, 0.975))[2]



opt_val[2:6] = round(opt_val[2:6], 2)
(skew_params = rbind(skew_params, opt_val))



write.table(skew_params, "tables/gaus_skew/SKEW/gaus_skew_params.csv", row.names = FALSE,
            col.names = TRUE,sep = ";")


CZ_data$stan_yhat = summary(gaus_skew, pars = c("y_hat"))$summary[,'mean']
CZ_data$stan_yhat_uci = summary(gaus_skew, pars = c("y_hat"))$summary[,'97.5%']
CZ_data$stan_yhat_lci = summary(gaus_skew, pars = c("y_hat"))$summary[,'2.5%']

CZ_data$stan_mountYN = rbinom(n = nrow(CZ_data),size = 1,prob = CZ_data$stan_yhat)

write.table(CZ_data, "tables/gaus_skew/SKEW/gaus_skew_mat.csv", row.names = FALSE, col.names = TRUE, sep = ";")



###################################
# plot observs and preds for skew #
###################################
CZ_data = read.csv("tables/gaus_skew/SKEW/gaus_skew_mat.csv", sep = ";")

# y = CZ_data$mountYNcontact
# y_rep = rstan::extract(CZ_mat_stan_size, pars = 'y_rep', permuted = TRUE)$y_rep
# dim(y_rep)
# ppc_bars(y,y_rep)
range(CZ_data$size_ratio)
breaks = c(-2,seq(-1.5,1.5,0.1),2)
# bin = cut(CZ_data$size_ratio,breaks)
# pdf("figures/gaus_size/gaus_size_ppc_bars_grouped.pdf")
# ppc_bars_grouped(y, y_rep, bin, prob = 0, freq = FALSE)
# dev.off()

CZ_data$bin = cut(CZ_data$size_ratio,breaks)
CZ_data_bin =
  CZ_data %>%
  group_by(bin) %>%
  dplyr::summarise(mount = mean(mountYNcontact),
                   uci_mount = CI(mountYNcontact)['upper'],
                   lci_mount = CI(mountYNcontact)['lower'],
                   mean_ratio = mean(size_ratio),
                   #y_rep = mean(y_rep),
                   preds_mount = mean(stan_mountYN)) %>%
  mutate(lci_mount = replace(lci_mount, which(lci_mount<0), 0))
# summary(CZ_data_bin)

pdf("figures/gaus_skew/SKEW/gaus_skew_preds.pdf", width=8, height=7)
ggplot(data = CZ_data) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(x = size_ratio,ymin = stan_yhat_lci, ymax = stan_yhat_uci), fill = "orange", alpha=0.3) +
  geom_errorbar(data = CZ_data_bin, aes(x = mean_ratio, ymin = lci_mount, ymax = uci_mount),alpha = 0.2) +
  scale_colour_manual(values=c("blue","orange2")) +
  geom_line(aes(size_ratio, stan_yhat, col="predictions")) +
  geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount, col="observations")) +
  labs(size="bin size",x="ln female size - ln male size",
       y="probability of mating",col="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.text = element_text(size = 15,face = "bold"), legend.position = 'none',
        axis.title = element_text(face = "bold", size = 15))
  # grids(linetype = "dashed")
dev.off()
