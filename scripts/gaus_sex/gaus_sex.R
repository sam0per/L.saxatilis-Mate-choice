rm(list = ls())

writeLines(readLines("scripts/gaus_skew/gaus_skew_hier.stan"))

args = commandArgs(trailingOnly=TRUE)

library(rstan)
options(mc.cores = parallel::detectCores(logical = FALSE) - 15)
#CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")

CZ_data = read.csv(args[1], sep = ";")

CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
CZ_data$test_sex=as.integer(CZ_data$test_sex) # 1 for female, 2 for male
CZ_data$ref_eco_sex = as.integer(CZ_data$ref_eco_sex) # 1 for crab female, 2 for crab male, 3 for wave female, 4 for wave male
dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio, shape = CZ_data$shape,
           N_shore = length(unique(CZ_data$shore)), N_eco = length(unique(CZ_data$ref_eco_sex)),
           N_sex = length(unique(CZ_data$test_sex)),
           shore = CZ_data$shore, ref = CZ_data$ref_eco_sex, test = CZ_data$test_sex)
rstan_options(auto_write = TRUE)
gaus_all = stan(file="scripts/min_gaus_all.stan", data=dat, iter = 8000, warmup = 2000, chains=2, refresh=8000,
                control = list(adapt_delta = 0.90, max_treedepth = 15))
saveRDS(gaus_all, "scripts/gaus_sex/gaus_sex.rds")


#Rscript scripts/gaus_skew/gaus_skew_hier.R -d data/CZ_all_mating_clean.csv -s scripts/gaus_skew/gaus_skew_hier.stan -i 4000 -c 2
#gaus_all = readRDS("models/gaus_skew/gaus_skew_hier.rds")
(CZ_all_pars = gaus_all@model_pars)
#library(parallel)
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
gaus_terms = c("level&scale", "choosiness&asymmetry", "c_coeff&preference")
#lapply(seq_along(all_pars_tbl), function(x) write.table(all_pars_tbl[[x]],
#                                                        paste0("tables/gaus_all/gaus_all_", gaus_terms[x],
#                                                               "_coef.csv"),
#                                                        row.names = FALSE, col.names = TRUE, sep = ";"))
lapply(seq_along(all_pars_tbl), function(x) write.table(all_pars_tbl[[x]],
                                                        paste0("tables/gaus_sex/gaus_skew_hier_", gaus_terms[x],
                                                               "_coef.csv"),
                                                        row.names = FALSE, col.names = TRUE, sep = ";"))

ls_pars$`6` = NULL
gaus_all_dens = mclapply(ls_pars, function(x) stan_dens(gaus_all, pars = CZ_all_pars[x]),mc.cores = 5)

gaus_all_plot = mclapply(ls_pars, function(x) stan_plot(gaus_all, pars = CZ_all_pars[x]),mc.cores = 4)

gaus_coef = c("mus", "sigmas", "gammas", "lambdas", "level")
gaus_coef = c("level", "scale", "choosiness", "asymmetry", "c_coeff")
#gaus_coef = c("level", "scale", "choosiness", "asymmetry", "c_int", "c_CZB", "c_CZC", "c_CZD", "ref_ecotypewave", "test_sexmale", "shape", "test_sexmale:shape")

#lapply(seq_along(gaus_all_dens), function(x) ggsave(filename = paste0("figures/gaus_all/gaus_all_",
#                                                                      gaus_coef[x], "_dens.png"),
#                                                    plot = gaus_all_dens[[x]]))
lapply(seq_along(gaus_all_dens), function(x) ggsave(filename = paste0("figures/gaus_sex/gaus_skew_hier",
                                                                      gaus_coef[x], "_dens.png"),
                                                    plot = gaus_all_dens[[x]]))
#lapply(seq_along(gaus_all_plot), function(x) ggsave(filename = paste0("figures/gaus_all/gaus_all_",
#                                                                      gaus_coef[x], "_plot.png"),
#                                                    plot = gaus_all_plot[[x]]))
lapply(seq_along(gaus_all_plot), function(x) ggsave(filename = paste0("figures/gaus_sex/gaus_sex_",
                                                                      gaus_coef[x], "_plot.png"),
                                                    plot = gaus_all_plot[[x]]))

