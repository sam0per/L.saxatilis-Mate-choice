rm(list = ls())

#############################
# Install required packages #
#############################
#if (!require("rstan")) install.packages("rstan")
#(.packages())

# List of packages for session
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc", "pander",
              "bbmle", "loo", "ggpubr", "cowplot", "purrr", "reshape2", "gridExtra", "grid")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

#load("Rpackages")
#for (p in setdiff(packages, installed.packages()[,"Package"])){
#  install.packages(p)
#}

###################################
# apply size model to field distr #
###################################
CZ_cline_params = read.csv("tables/CZ_cline_params.csv", sep = ";")
CZ_size_params = read.csv("tables/CZ_size_params.csv", sep = ";")

######
## crab habitat
######
# generate size distributions for each sex and ecotype using cline parameters
male_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[5]), sd = abs(x[9]))),2)
female_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs(x[5]+x[7]), sd = abs(x[9]))),2)

# pair each female with every male within habitat and contact zone and compute mounting success YN
CZ_cline_mountYN_c = matrix(nrow = nrow(female_c), ncol = nrow(male_c))
CZ_clines_YN_c = list()
for (i in 1:ncol(female_c)){
  for (j in 1:nrow(female_c)) {
    CZ_cline_mountYN_c[j,] = rbinom(n = nrow(female_c), size = 1,
                                    prob = inv.logit(CZ_size_params$mean[CZ_size_params$params=='level'] +
                                                       CZ_size_params$mean[CZ_size_params$params=='scale'] *
                                                       exp(-0.5 * (((female_c[j,i] - male_c[,i]) -
                                                                      CZ_size_params$mean[CZ_size_params$params=='preference'])/
                                                                     CZ_size_params$mean[CZ_size_params$params=='choosiness'])^2) +
                                                       CZ_size_params$mean[CZ_size_params$params=='asymmetry'] *
                                                       (female_c[j,i] - male_c[,i])))
    
  }
  CZ_clines_YN_c[[i]] = CZ_cline_mountYN_c
  names(CZ_clines_YN_c)[i] = colnames(female_c)[i]
}

# random sample of 1 successful male per female (variation in mounting success only in males)
CZ_clines_Yidx_c = lapply(CZ_clines_YN_c, function(x) which(x==1, arr.ind = TRUE)) %>%
  lapply(., function(x) as.data.frame(x[order(x[,'row']), ])) %>%
  map(., function(x) group_by(x, row)) %>% map(., function(x) sample_n(x, size = 1, replace = TRUE)) %>%
  lapply(., setNames, nm = c("row", "col")) %>%
  lapply(., function(x) data.frame(x[,'row'], x[,'col'], mount=1))

# random sample of 1 non-successful male per female (variation in mounting success only in males)
CZ_clines_Nidx_c = lapply(CZ_clines_YN_c, function(x) which(x==0, arr.ind = TRUE)) %>%
  lapply(., function(x) as.data.frame(x[order(x[,'row']), ])) %>%
  map(., function(x) group_by(x, row)) %>% map(., function(x) sample_n(x, size = 1, replace = TRUE)) %>%
  lapply(., function(x) data.frame(x[,'row'], x[,'col'], mount=0))

# merge the succ and non-succ mating pairs
YNidx_c = lapply(names(CZ_clines_YN_c), function(x) arrange(rbind(CZ_clines_Yidx_c[[x]], CZ_clines_Nidx_c[[x]]), row))
names(YNidx_c) = colnames(female_c)

# retrieve sizes (ln) of females and the sampled males from the generated male size distribution
CZ_dat_YN_c = lapply(names(YNidx_c), function(x) tibble(shore=x,
                                                        female=female_c[YNidx_c[[x]][['row']], x],
                                                        male=male_c[YNidx_c[[x]][['col']], x],
                                                        mountYN=YNidx_c[[x]][['mount']]))
names(CZ_dat_YN_c) = colnames(female_c)

# save the simulated datasets of mate choice
lapply(names(CZ_dat_YN_c), function(x) write.table(CZ_dat_YN_c[[x]], paste0("tables/", x, "_crab_sim_YN.csv"),
                                                   row.names = FALSE, col.names = TRUE, sep = ";"))


#################
## hybrid habitat
#################
# generate size distributions for each sex and ecotype using cline parameters
male_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs((x[5]+x[6])/2), sd = abs(x[10]))),2)
female_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = abs((x[5]+x[6]+x[7]+x[8])/2),
                                                               sd = sqrt(abs(x[10])))),2)
CZ_cline_mountYN_h = matrix(nrow = nrow(female_h), ncol = nrow(male_h))
CZ_clines_YN_h = list()
for (i in 1:ncol(female_h)){
  for (j in 1:nrow(female_h)) {
    CZ_cline_mountYN_h[j,] = rbinom(n = nrow(female_h), size = 1,
                                    prob = inv.logit(CZ_size_params$mean[CZ_size_params$params=='level'] +
                                                       CZ_size_params$mean[CZ_size_params$params=='scale'] *
                                                       exp(-0.5 * (((female_h[j,i] - male_h[,i]) -
                                                                      CZ_size_params$mean[CZ_size_params$params=='preference'])/
                                                                     CZ_size_params$mean[CZ_size_params$params=='choosiness'])^2) +
                                                       CZ_size_params$mean[CZ_size_params$params=='asymmetry'] *
                                                       (female_h[j,i] - male_h[,i])))
    
  }
  CZ_clines_YN_h[[i]] = CZ_cline_mountYN_h
  names(CZ_clines_YN_h)[i] = colnames(female_h)[i]
}

# random sample of 1 successful male per female (variation in mounting success only in males)
CZ_clines_Yidx_h = lapply(CZ_clines_YN_h, function(x) which(x==1, arr.ind = TRUE)) %>%
  lapply(., function(x) as.data.frame(x[order(x[,'row']), ])) %>%
  map(., function(x) group_by(x, row)) %>% map(., function(x) sample_n(x, size = 1, replace = TRUE)) %>%
  lapply(., setNames, nm = c("row", "col")) %>%
  lapply(., function(x) data.frame(x[,'row'], x[,'col'], mount=1))

# random sample of 1 non-successful male per female (variation in mounting success only in males)
CZ_clines_Nidx_h = lapply(CZ_clines_YN_h, function(x) which(x==0, arr.ind = TRUE)) %>%
  lapply(., function(x) as.data.frame(x[order(x[,'row']), ])) %>%
  map(., function(x) group_by(x, row)) %>% map(., function(x) sample_n(x, size = 1, replace = TRUE)) %>%
  lapply(., function(x) data.frame(x[,'row'], x[,'col'], mount=0))

# merge the succ and non-succ mating pairs
YNidx_h = lapply(names(CZ_clines_YN_h), function(x) arrange(rbind(CZ_clines_Yidx_h[[x]], CZ_clines_Nidx_h[[x]]), row))
names(YNidx_h) = colnames(female_h)

# retrieve sizes (ln) of females and the sampled males from the generated male size distribution
CZ_dat_YN_h = lapply(names(YNidx_h), function(x) tibble(shore=x,
                                                        female=female_h[YNidx_h[[x]][['row']], x],
                                                        male=male_h[YNidx_h[[x]][['col']], x],
                                                        mountYN=YNidx_h[[x]][['mount']]))
names(CZ_dat_YN_h) = colnames(female_h)

# save the simulated datasets of mate choice
lapply(names(CZ_dat_YN_h), function(x) write.table(CZ_dat_YN_h[[x]], paste0("tables/", x, "_hybrid_sim_YN.csv"),
                                                   row.names = FALSE, col.names = TRUE, sep = ";"))


#######################################
# stan model with only simulated size #
#######################################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))

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

shore = c("CZA", "CZB", "CZC", "CZD")
ecotype = c("crab", "hybrid", "wave")

CZ_sim_dat = list()
dat = list()
sim_mat_stan_size = list()
sim_size_params = list()


for (eco in ecotype) {
  if (eco=="crab") {
    sim_dat = list.files("tables", pattern = eco, full.names = TRUE)
    for (f in 1:length(sim_dat)) {
      CZ_sim_dat[[f]] = read.csv(sim_dat[f], sep = ";")
      dat[[f]] = list(N = nrow(CZ_sim_dat[[f]]), y=CZ_sim_dat[[f]][,'mountYN'],
                      ratio = CZ_sim_dat[[f]][,'female'] - CZ_sim_dat[[f]][,'male'])
      sim_mat_stan_size[[f]] = stan(file = "scripts/CZ_sim_gaus_size.stan", data = dat[[f]])
      sim_size_pars = sim_mat_stan_size[[f]]@model_pars[1:5]
      sim_size_params[[f]] = round(summary(sim_mat_stan_size[[f]], pars = sim_size_pars,probs=c(0.025, 0.975))$summary,2)
      sim_size_params[[f]] = rownames_to_column(as.data.frame(sim_size_params[[f]]), var="params")
      CZ_sim_dat[[f]][,'y_rep'] = summary(sim_mat_stan_size[[f]], pars = c("y_rep"))$summary[,'mean']
      CZ_sim_dat[[f]][,'y_rep_se'] = summary(sim_mat_stan_size[[f]], pars = c("y_rep"))$summary[,'se_mean']
      CZ_sim_dat[[f]][,'y_hat'] = summary(sim_mat_stan_size[[f]], pars = c("y_hat"))$summary[,'mean']
      CZ_sim_dat[[f]][,'preds'] = round(inv.logit(CZ_sim_dat[[f]][,'y_hat']), 2)
      CZ_sim_dat[[f]][,'uci'] = summary(sim_mat_stan_size[[f]], pars = c("y_hat"))$summary[,'97.5%']
      CZ_sim_dat[[f]][,'uci_preds'] = round(inv.logit(CZ_sim_dat[[f]][,'uci']), 2)
      CZ_sim_dat[[f]][,'lci'] = summary(sim_mat_stan_size[[f]], pars = c("y_hat"))$summary[,'2.5%']
      CZ_sim_dat[[f]][,'lci_preds'] = round(inv.logit(CZ_sim_dat[[f]][,'lci']), 2)
      write.table(CZ_sim_dat[[f]], paste0("tables/", shore[f], "_", eco,"_sim_YN.csv"),
                  row.names = FALSE, col.names = TRUE, sep = ";")
    }
    sim_size_dens = lapply(sim_mat_stan_size, function(x) stan_dens(x, pars = sim_size_pars))
    names(sim_size_dens) = shore
    lapply(names(sim_size_dens), function(x) ggsave(filename=paste0("figures/", x, "_", eco, "_sim_size_dens.png"),
                                                    plot=sim_size_dens[[x]]))
    sim_size_pars_tab = lapply(sim_size_params, tableGrob)
    names(sim_size_pars_tab) = shore
    lapply(names(sim_size_pars_tab), function(x) ggsave(filename=paste0("tables/", x, "_", eco, "_sim_size_pars.png"),
                                                        plot=sim_size_pars_tab[[x]]))
    
    sim_ss_mean = map(CZ_sim_dat, function(x) group_by(x, mountYN)) %>%
      map(., function(x) summarise_at(x, vars(male), mean))
    sim_ss = lapply(1:length(CZ_sim_dat), function(x) {
      ggplot() +
        geom_ribbon(data = CZ_sim_dat[[x]],
                    aes(x = exp(male), ymin = lci_preds, ymax = uci_preds, fill=as.factor(mountYN)),
                    alpha=0.3) +
        geom_point(data = CZ_sim_dat[[x]], aes(x=exp(male), y=preds, col=as.factor(mountYN))) +
        geom_vline(data = sim_ss_mean[[x]], aes(xintercept = exp(male),
                                                col=as.factor(mountYN))) +
        labs(title = paste(shore[x], eco, sep = " "), y="predicted probability", x="male size",
             fill="mount success", col="mount success") +
        scale_fill_manual(labels = c("no", "yes"), values = c("red", "blue")) +
        scale_color_manual(labels = c("no", "yes"), values = c("red", "blue"))
    })
    sim_ss_plot = do.call(grid_arrange_shared_legend, sim_ss)
    ggsave(filename = paste0("figures/CZ_", eco, "_sim_size_ss.png"), plot = sim_ss_plot)
  }
}




CZ_data = read.csv("data/CZ_mating_clean.csv",sep = ";")
CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD

qplot(CZ_data$shore)
qplot(CZ_data$shape)
qplot(rnorm(n = 100))
hist(rmultinom(n = 10, size = 4, prob = 0.2))


