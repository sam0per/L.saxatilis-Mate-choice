# Usage example:
# Rscript L.saxatilis-Mate-choice/scripts/mods_comp.R --modelone models/gaus_skew/SKEW/gaus_skew.rds \
#   --modeltwo models/gaus_skew/BCDG_shore/gaus_skew_hier_BCDG_shore.rds --modelthree models/gaus_skew/BCDG_eco/gaus_skew_hier_BCDG_eco.rds \
#   --modelfour models/gaus_skew/BCDG/gaus_skew_hier_BCDG.rds -d data/CZ_all_mating_clean.csv -o tables/mods_comp/comp_


rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc", "pander",
              "bbmle", "loo", "ggpubr", "cowplot", "purrr", "reshape2", "gridExtra", "grid", "arm", "parallel",
              "rstantools", "optparse", "pROC")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option("--modelone", type="character", default=NULL,
              help="first model for comparison", metavar="character"),
  make_option("--modeltwo", type="character", default=NULL,
              help="second model for comparison", metavar="character"),
  make_option("--modelthree", type="character", default=NULL,
              help="third model for comparison", metavar="character"),
  make_option("--modelfour", type="character", default=NULL,
              help="fourth model for comparison", metavar="character"),
  make_option("--nnhier", type="character", default=NULL,
              help="non-hierarchical model for comparison", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="mods_comp_out.csv",
              help="output file name [default = %default]", metavar="character"),
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="input data", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$modelone) | is.null(opt$modeltwo) | is.null(opt$data)) {
  print_help(opt_parser)
  stop("At least three arguments must be supplied (two models and one dataset).\n", call.=FALSE)
}

########################################################
###### approximate leave-one-out cross-validation ######
# http://mc-stan.org/loo/articles/loo2-with-rstan.html #
########################################################
CZ_data = read.csv(opt$data, sep = ";")
# mod1 = readRDS("models/gaus_skew/gaus_skew.rds")
mod1 = readRDS(opt$modelone)
# mod2 = readRDS("models/gaus_skew/gaus_skew_hier_BCDG_shore.rds")
mod2 = readRDS(opt$modeltwo)
mod3 = readRDS(opt$modelthree)
mod4 = readRDS(opt$modelfour)
nnhier = readRDS(opt$nnhier)

out_comp_str = lapply(c(opt$modelone,opt$modeltwo,opt$modelthree,opt$modelfour,opt$nnhier), function(x) {
  modsplit = strsplit(strsplit(basename(x), "[.]")[[1]][1], split = "_")[[1]]
  return(modsplit[length(modsplit)])
})

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(mod1, merge_chains = FALSE)
log_lik_2 <- extract_log_lik(mod2, merge_chains = FALSE)
log_lik_3 <- extract_log_lik(mod3, merge_chains = FALSE)
log_lik_4 <- extract_log_lik(mod4, merge_chains = FALSE)
log_lik_nnhier <- extract_log_lik(nnhier, merge_chains = FALSE)

# provide relative effective sample sizes
r_eff <- relative_eff(exp(log_lik_1))
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 4, save_psis = TRUE)
print(paste0("LOO results for model ", basename(opt$modelone)))
print(loo_1)

# split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
# out_comp = "tables/mods_comp/comp_skew_shore.csv"
# inp_mod1 = "models/gaus_skew/gaus_skew.rds"
out_mod1 = paste0(dirname(opt$out), "/loo_", strsplit(basename(opt$modelone), "[.]")[[1]][1], ".csv")
write.table(round(loo_1$estimates, 2), out_mod1, sep = ",", row.names = TRUE, col.names = TRUE)

# plot(loo_1)
r_eff_2 <- relative_eff(exp(log_lik_2))
loo_2 <- loo(log_lik_2, r_eff = r_eff_2, cores = 4, save_psis = TRUE)
print(paste0("LOO results for model ", basename(opt$modeltwo)))
print(loo_2)
out_mod2 = paste0(dirname(opt$out), "/loo_", strsplit(basename(opt$modeltwo), "[.]")[[1]][1], ".csv")
write.table(round(loo_2$estimates, 2), out_mod2, sep = ",", row.names = TRUE, col.names = TRUE)

r_eff_3 <- relative_eff(exp(log_lik_3))
loo_3 <- loo(log_lik_3, r_eff = r_eff_3, cores = 4, save_psis = TRUE)
print(paste0("LOO results for model ", basename(opt$modelthree)))
print(loo_3)
out_mod3 = paste0(dirname(opt$out), "/loo_", strsplit(basename(opt$modelthree), "[.]")[[1]][1], ".csv")
write.table(round(loo_3$estimates, 2), out_mod3, sep = ",", row.names = TRUE, col.names = TRUE)

r_eff_4 <- relative_eff(exp(log_lik_4))
loo_4 <- loo(log_lik_4, r_eff = r_eff_4, cores = 4, save_psis = TRUE)
print(paste0("LOO results for model ", basename(opt$modelfour)))
print(loo_4)
out_mod4 = paste0(dirname(opt$out), "/loo_", strsplit(basename(opt$modelfour), "[.]")[[1]][1], ".csv")
write.table(round(loo_4$estimates, 2), out_mod4, sep = ",", row.names = TRUE, col.names = TRUE)

r_eff_nnhier <- relative_eff(exp(log_lik_nnhier))
loo_nnhier <- loo(log_lik_nnhier, r_eff = r_eff_nnhier, cores = 4, save_psis = TRUE)
print(paste0("LOO results for model ", basename(opt$nnhier)))
print(loo_nnhier)
out_modnnhier = paste0(dirname(opt$out), "/loo_", strsplit(basename(opt$nnhier), "[.]")[[1]][1], ".csv")
write.table(round(loo_nnhier$estimates, 2), out_modnnhier, sep = ",", row.names = TRUE, col.names = TRUE)

print(paste0("Parameters in model ", basename(opt$modelone)))
mod1@model_pars
print(paste0("Parameters in model ", basename(opt$modeltwo)))
mod2@model_pars
print(paste0("Parameters in model ", basename(opt$modelthree)))
mod3@model_pars
print(paste0("Parameters in model ", basename(opt$modelfour)))
mod4@model_pars
print(paste0("Parameters in non-hierarchical model ", basename(opt$nnhier)))
nnhier@model_pars

# plot(loo_2, label_points = TRUE)
comp1_2 <- compare(loo_1, loo_2)
cat(paste0("estimated difference of expected leave-one-out prediction errors\nbetween ",
           basename(opt$modelone), " and ", basename(opt$modeltwo), " along with the standard error.\n",
           " Positive difference in elpd (and its scale relative to the standard error)\nindicates a preference for the second model.\n"))
print(comp1_2)
elpd_diff = comp1_2[1]
se_elpd_diff = comp1_2[2]
write.table(data.frame(elpd_diff, se_elpd_diff), file = paste0("tables/mods_comp/comp_", out_comp_str[[1]], "_", out_comp_str[[2]], ".csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")

comp1_3 <- compare(loo_1, loo_3)
cat(paste0("estimated difference of expected leave-one-out prediction errors\nbetween ",
           basename(opt$modelone), " and ", basename(opt$modelthree), " along with the standard error.\n",
           " Positive difference in elpd (and its scale relative to the standard error)\nindicates a preference for the second model.\n"))
print(comp1_3)
elpd_diff13 = comp1_3[1]
se_elpd_diff13 = comp1_3[2]
write.table(data.frame(elpd_diff13, se_elpd_diff13), file = paste0("tables/mods_comp/comp_", out_comp_str[[1]], "_", out_comp_str[[3]], ".csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")

comp1_4 <- compare(loo_1, loo_4)
cat(paste0("estimated difference of expected leave-one-out prediction errors\nbetween ",
           basename(opt$modelone), " and ", basename(opt$modelfour), " along with the standard error.\n",
           " Positive difference in elpd (and its scale relative to the standard error)\nindicates a preference for the second model.\n"))
print(comp1_4)
elpd_diff14 = comp1_4[1]
se_elpd_diff14 = comp1_4[2]
write.table(data.frame(elpd_diff14, se_elpd_diff14), file = paste0("tables/mods_comp/comp_", out_comp_str[[1]], "_", out_comp_str[[4]], ".csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")

comp2_3 <- compare(loo_2, loo_3)
cat(paste0("estimated difference of expected leave-one-out prediction errors\nbetween ",
           basename(opt$modeltwo), " and ", basename(opt$modelthree), " along with the standard error.\n",
           " Positive difference in elpd (and its scale relative to the standard error)\nindicates a preference for the second model.\n"))
print(comp2_3)
elpd_diff23 = comp2_3[1]
se_elpd_diff23 = comp2_3[2]
write.table(data.frame(elpd_diff23, se_elpd_diff23), file = paste0("tables/mods_comp/comp_", out_comp_str[[2]], "_", out_comp_str[[3]], ".csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")

comp2_4 <- compare(loo_2, loo_4)
cat(paste0("estimated difference of expected leave-one-out prediction errors\nbetween ",
           basename(opt$modeltwo), " and ", basename(opt$modelfour), " along with the standard error.\n",
           " Positive difference in elpd (and its scale relative to the standard error)\nindicates a preference for the second model.\n"))
print(comp2_4)
elpd_diff24 = comp2_4[1]
se_elpd_diff24 = comp2_4[2]
write.table(data.frame(elpd_diff24, se_elpd_diff24), file = paste0("tables/mods_comp/comp_", out_comp_str[[2]], "_", out_comp_str[[4]], ".csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")

comp3_4 <- compare(loo_3, loo_4)
cat(paste0("estimated difference of expected leave-one-out prediction errors\nbetween ",
           basename(opt$modelthree), " and ", basename(opt$modelfour), " along with the standard error.\n",
           " Positive difference in elpd (and its scale relative to the standard error)\nindicates a preference for the second model.\n"))
print(comp3_4)
elpd_diff34 = comp3_4[1]
se_elpd_diff34 = comp3_4[2]
write.table(data.frame(elpd_diff34, se_elpd_diff34), file = paste0("tables/mods_comp/comp_", out_comp_str[[3]], "_", out_comp_str[[4]], ".csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")

# y_rep_mod1 = summary(mod1, pars = c("y_rep"))$summary[,'mean']
# roc_obj_mod1 = roc(CZ_data$mountYNcontact, y_rep_mod1, ci = TRUE)
# cat(paste0("Area under the curve and 95% CI for ", basename(opt$modelone), ".\n"))
# auc(roc_obj_mod1)
# ci.auc(roc_obj_mod1)
#
# y_rep_mod2 = summary(mod2, pars = c("y_rep"))$summary[,'mean']
# roc_obj_mod2 = roc(CZ_data$mountYNcontact, y_rep_mod2, ci = TRUE)
# cat(paste0("Area under the curve and 95% CI for ", basename(opt$modeltwo), ".\n"))
# auc(roc_obj_mod2)
# ci.auc(roc_obj_mod2)
#
# cat(paste0("Bootstrap test for ROC comparison between ", basename(opt$modelone), " and ", basename(opt$modeltwo), ".\n"))
# roc.test(roc_obj_mod1, roc_obj_mod2, method = "bootstrap")

# Marginal posterior predictive checks
# y_rep1 = rstan::extract(mod1, pars = 'y_rep', permuted = TRUE)$y_rep
# ppc_loo_pit_overlay(
#  y = CZ_data$mountYN,
#  yrep = y_rep1,
#  lw = weights(loo_1$psis_object)
# )

# The excessive number of values close to 1 indicates that the model is under-dispersed compared to the data,
# and we should consider a model that allows for greater dispersion.
# y_rep2 = rstan::extract(mod2, pars = 'y_rep', permuted = TRUE)$y_rep
# ppc_loo_pit_overlay(
#  y = CZ_data$mountYN,
#  yrep = y_rep2,
#  lw = weights(loo_2$psis_object)
# )

# loo_1 <- loo(log_lik_1, k_threshold=0.7, cores = 2)
# loo_2 <- loo(log_lik_2, k_threshold=0.7, cores = 2)
# loo_list <- list(loo_1, loo_2)
# loo_model_weights(loo_list)
# loo_model_weights(loo_list, method = "pseudobma")
# loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
