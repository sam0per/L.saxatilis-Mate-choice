rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc", "pander",
              "bbmle", "loo", "ggpubr", "cowplot", "purrr", "reshape2", "gridExtra", "grid", "arm", "parallel",
              "rstantools", "optparse", "VGAM")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="input data", metavar="character"),
  make_option(c("-s", "--stanfile"), type="character", default=NULL, 
              help="model written in Stan", metavar="character")) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$stanfile)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input data and stan file).\n", call.=FALSE)
}


CZ_data = read.csv(opt$data, sep = ";")
#CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
#CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
#head(CZ_data)
#summary(CZ_data)


############################
# stan model with skewness #
############################
rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores(logical = FALSE) - 15)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)

CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio)
#str(dat)
#dat$posterior_predictive = 1
#fit.prior.pred = stan(file = "scripts/CZ_mating_gaus_prior_size.stan", data = dat)

# gaus_skew = rstan::stan(file = opt$stanfile, data = dat, iter = 8000, warmup = 2000,
#                         chains=4, refresh=8000,
#                         control = list(adapt_delta = 0.90, max_treedepth = 15))
gaus_skew = rstan::stan(file = "scripts/gaus_skew/gaus_skew.stan", data = dat, iter = 8000, warmup = 2000,
                        chains=4, refresh=8000,
                        control = list(adapt_delta = 0.90, max_treedepth = 15))

saveRDS(gaus_skew, "models/gaus_skew/gaus_skew.rds")

#############################
# plot post distr gaus skew #
#############################
gaus_skew = readRDS("models/gaus_skew/gaus_skew.rds")
gaus_skew@model_pars
gaus_size_pars = c("level","scale","centre","choosiness","asymmetry")
#gaus_size_parfig = c("preference","choosiness","asymmetry")

list_of_draws <- rstan::extract(gaus_skew)
names(list_of_draws) = c("level","scale","centre","choosiness","asymmetry", "y_hat", "lp__")
print(names(list_of_draws))
parfig = lapply(gaus_size_pars, function(x) {
  ggplot() +
    geom_density(aes(list_of_draws[[x]]), fill='red', col='black') +
    labs(x="", title = x) +
    theme(axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold",size = 15))
})
names(parfig) = gaus_size_pars


# find derivative of skew-normal in Wolphram|Alpha typing
# d[l + s * exp(-0.5 * ((x - p) / c)^2) * (1 + ERF(a * ((x - p) / (1.414214 * c)))) ,x]

opt = function(scale, centre, asymmetry, choosiness, x){
  # derivative
  (scale * exp(-(0.5 * (centre - x)^2)/choosiness^2) *
     (0.797884 * asymmetry * choosiness * exp(-(0.5 * asymmetry^2 * (centre - x)^2)/choosiness^2) -
        (x - centre) * (erf((0.707107 * asymmetry * (x - centre))/choosiness) + 1)))/choosiness^2
}
pars = list(scale = list_of_draws$scale, centre = list_of_draws$centre, asymmetry = list_of_draws$asymmetry,
            choosiness = list_of_draws$choosiness)

# find the root of the derivative
#str(xmin <- uniroot(opt, c(0, 1), tol = 0.0001, scale = 0.44, preference = -0.07, asymmetry = 1.35, choosiness = 0.67))

opt_draws = sapply(1:24000, function(z){
  uniroot(opt, c(0, 1), tol = 0.0001, asymmetry = pars[['asymmetry']][z], scale = pars[['scale']][z],
          centre = pars[['centre']][z], choosiness = pars[['choosiness']][z])$root
})

list_of_draws$preference = opt_draws

opt_dens = ggplot() +
  geom_density(aes(list_of_draws$preference), fill='red', col='black') +
  labs(x="", title = "preference") +
  theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
mean(list_of_draws$preference)
parfig$preference = opt_dens


pdf("figures/gaus_skew/gaus_skew_pars_dens.pdf",width = 10, height = 7)
#do.call(ggarrange, parfig)
ggarrange(parfig$level, parfig$scale, parfig$centre, parfig$choosiness, parfig$asymmetry, parfig$preference)
dev.off()


#########################
# pars values gaus skew #
#########################
CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
gaus_size_pars = c("level","scale","preference","choosiness","asymmetry")
(skew_params = round(summary(gaus_skew, pars = gaus_size_pars, probs=c(0.025, 0.975))$summary,2))
class(gaus_skew)

y_hat = summary(gaus_skew, pars = c("y_hat"))$summary[,'mean']
#CZ_erf = summary(gaus_skew, pars = c("err_val"))$summary[,'mean']
#y_prob = y_hat * (1 + erf(CZ_erf))
plot(CZ_data$size_ratio, y_hat)

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
                     mean = mean(list_of_draws$optimum),
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
colnames(opt_val) = colnames(skew_params)
opt_val['2.5%'] = quantile(samp, c(0.025, 0.975))[1]
opt_val['97.5%'] = quantile(samp, c(0.025, 0.975))[2]



opt_val[2:6] = round(opt_val[2:6], 2)
(skew_params = rbind(opt_val, skew_params))



write.table(skew_params, "tables/gaus_skew/gaus_skew_params.csv", row.names = FALSE,
            col.names = TRUE,sep = ";")


CZ_data$stan_yhat = summary(gaus_skew, pars = c("y_hat"))$summary[,'mean']
CZ_data$stan_yhat_uci = summary(gaus_skew, pars = c("y_hat"))$summary[,'97.5%']
CZ_data$stan_yhat_lci = summary(gaus_skew, pars = c("y_hat"))$summary[,'2.5%']

CZ_data$stan_mountYN = rbinom(n = nrow(CZ_data),size = 1,prob = CZ_data$stan_yhat)

write.table(CZ_data, "tables/gaus_skew/gaus_skew_mat.csv", row.names = FALSE, col.names = TRUE, sep = ";")



###################################
# plot observs and preds for skew #
###################################                                                    
CZ_data = read.csv("tables/gaus_skew/gaus_skew_mat.csv", sep = ";")

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
summary(CZ_data_bin)

pdf("figures/gaus_skew/gaus_skew_preds.pdf", width=8, height=7)
ggplot(data = CZ_data) +
  geom_vline(xintercept = 0) +
  geom_ribbon(aes(x = size_ratio,ymin = stan_yhat_lci, ymax = stan_yhat_uci), fill = "orange", alpha=0.3) +
  geom_errorbar(data = CZ_data_bin, aes(x = mean_ratio, ymin = lci_mount, ymax = uci_mount),alpha = 0.2) +
  scale_colour_manual(values=c("blue","orange2")) +
  geom_line(aes(size_ratio, stan_yhat, col="predictions")) +
  geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount, col="observations")) +
  labs(size="bin size",x="ln(female) - ln(male size)",
       y="probability of mounting",col="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.title = element_text(size = 11,face = "bold"), legend.position = c(0.05, 0.85),
        axis.title = element_text(face = "bold", size = 13)) +
  grids(linetype = "dashed")
dev.off()


###################################
# apply skew model to field distr #
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
skew_params = read.csv("tables/gaus_skew/gaus_skew_params.csv", sep = ";")


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
      p = skew_params$mean[skew_params$parameter=='level'] + skew_params$mean[skew_params$parameter=='scale'] *
        exp(-0.5 * (((fem - m) - skew_params$mean[skew_params$parameter=='centre']) /
                      skew_params$mean[skew_params$parameter=='choosiness'])^2) *
        (1 + erf(skew_params$mean[skew_params$parameter=='asymmetry'] *
                   ((fem - m) - skew_params$mean[skew_params$parameter=='centre']) /
                   (1.414214 * skew_params$mean[skew_params$parameter=='choosiness'])))
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
                                                   paste0("tables/gaus_skew/sims/", shore[y], "_", ecotype[x], "_sim_YN.csv"),
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

# run the for loop with eco=="crab", eco=="hybrid", eco=="wave"
for (eco in ecotype) {
  if (eco=="wave") {
    sim_dat = list.files("tables/gaus_skew/BCDG/sims", pattern = paste0(eco, "_sim_YN"), full.names = TRUE)
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
    lapply(1:length(sim_fit_ci), function(x) write.table(sim_fit_ci[[x]], paste0("tables/gaus_skew/BCDG/sex_sel/", shore[x], "_",
                                                                                 eco, "_sim_ss_fit.csv"),
                                                         row.names = FALSE, col.names = TRUE, sep = ";"))
    sim_ss_coef = lapply(1:length(CZ_sim_ss), function(x) ggarrange(tableGrob(round(summary(CZ_sim_ss[[x]])$coef, 2))))
    lapply(1:length(sim_ss_coef), function(x) ggsave(filename = paste0("tables/gaus_skew/BCDG/sex_sel/", shore[x], "_",
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
        labs(x="ln(male size)", y="predicted prob", col="", title=paste(shore[x], eco, sep = " ")) +
        scale_color_manual(values = c("blue", "orange"))
    })
    sim_ss_plot = do.call(grid_arrange_shared_legend, CZ_sim_ss_plot)
    ggsave(filename = paste0("figures/gaus_skew/BCDG/sex_sel/CZ_", eco, "_sim_size_ss.png"), plot = sim_ss_plot)
  }
}
dev.off()

rm(list = ls())
shore = c("CZA", "CZB", "CZC", "CZD")
ecotype = c("crab", "hybrid", "wave")

CZ_eco_chr = list(CZA=NULL, CZB=NULL, CZC=NULL, CZD=NULL)
for (s in seq_along(shore)) {
  for (e in seq_along(ecotype)) {
    CZ_eco_chr[[s]][e] = list.files("tables/gaus_skew/BCDG/sex_sel", pattern = paste(shore[s], ecotype[e], "sim", "ss", "fit",sep = "_"),
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
    labs(x="habitat", y="ln(male size)", title=shore[x]) +
    rremove("legend") +
    rremove("xlab")
})
CZ_eco_ss[[1]] = CZ_eco_ss[[1]] + theme(axis.text.x=element_blank())
CZ_eco_ss[[2]] = CZ_eco_ss[[2]] + theme(axis.text.x=element_blank(),
                                        axis.title.y=element_blank())
CZ_eco_ss[[4]] = CZ_eco_ss[[4]] + theme(axis.title.y=element_blank())
png("figures/gaus_skew/BCDG/sex_sel/CZs_sim_eco_ss.png", width=580, height=500)
do.call(egg::ggarrange, CZ_eco_ss)
dev.off()

lapply(seq_along(CZ_eco_ss), function(x) ggsave(filename = paste0("figures/gaus_skew/BCDG/sex_sel/", shore[x],
                                                                  "_sim_eco_ss.png"),
                                                plot = CZ_eco_ss[[x]]))


CZ_male_mean = map(CZ_eco_df, function(x) group_by(x, factor(ecotype))) %>%
  map(., function(x) dplyr::summarise(x, male_mean=mean(male))) %>%
  lapply(., setNames, nm = c("ecotype", "male_mean"))

CZ_male_diff = map(seq_along(shore), function(x) full_join(CZ_male_mean[[x]], CZ_eco_male_mx[[x]])) %>%
  map(., function(x) mutate(x, male_diff=round(male_mean-male_mx, 2)))

lapply(seq_along(shore), function(x) write.table(CZ_male_diff[[x]], paste0("tables/gaus_skew/BCDG/sex_sel/", shore[x],
                                                                           "_sim_ss_mdiff.csv"),
                                                 row.names = FALSE, col.names = TRUE, sep = ";"))


######################
# assortative mating #
######################
CZ_eco_am = lapply(seq_along(CZ_eco_df), function(x) {
  ggplot() +
    geom_point(data = subset(CZ_eco_df[[x]], mountYN==1), aes(x = female, y = male, col=factor(ecotype))) +
    facet_wrap(~factor(ecotype)) +
    labs(x="ln(female size)", y="ln(male size)", title=shore[x]) +
    grids(linetype = "dashed") +
    #annotate("text", x=0, y=4, label=paste0('atop(bold("p_value is ',p_val,'"))'))
    rremove("legend")
})

lapply(seq_along(CZ_eco_am), function(x) ggsave(filename = paste0("figures/gaus_skew/BCDG/ass_mat/", shore[x],
                                                                  "_sim_eco_am.png"),
                                                plot = CZ_eco_am[[x]]))

CZ_eco_amr = lapply(CZ_eco_df, subset, mountYN==1) %>%
  map(., group_by, factor(ecotype)) %>%
  map(., function(x) dplyr::summarise(x, p_cor = round(cor(female, male, method = "pearson"), 2))) %>%
  lapply(., setNames, nm = c("ecotype", "p_cor"))

lapply(seq_along(shore), function(x) write.table(CZ_eco_amr[[x]], paste0("tables/gaus_skew/BCDG/ass_mat/", shore[x],
                                                                         "_sim_am_pcor.csv"),
                                                 row.names = FALSE, col.names = TRUE, sep = ";"))
names(CZ_eco_amr) = shore
(eco_pcor_table = data.frame(ecotype = ecotype, CZA = CZ_eco_amr$CZA$p_cor, CZB = CZ_eco_amr$CZB$p_cor,
                             CZC = CZ_eco_amr$CZC$p_cor, CZD = CZ_eco_amr$CZD$p_cor))
write.table(eco_pcor_table, file = "tables/gaus_skew/BCDG/ass_mat/CZs_sim_am_pcor.csv",
            row.names = FALSE, col.names = TRUE, sep = ";")

lapply(CZ_eco_df, subset, mountYN==0) %>%
  map(., group_by, factor(ecotype)) %>%
  map(., function(x) dplyr::summarise(x, p_cor = round(cor(female, male, method = "pearson"), 2))) %>%
  lapply(., setNames, nm = c("ecotype", "p_cor"))
