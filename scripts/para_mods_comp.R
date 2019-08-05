# Usage example:
# Rscript L.saxatilis-Mate-choice/scripts/para_mods_comp.R --focus models/gaus_skew/B_all/gaus_skew_hier_B_all.rds --modelone models/gaus_skew/Bhyp_shore/gaus_skew_hier_Bhyp_shore.rds
#   -o tables/mods_comp/

rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot",
              "loo", "purrr", "reshape2", "arm", "parallel",
              "rstantools", "optparse")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(.packages, require, character.only=TRUE)

option_list = list(
  make_option("--focus", type="character", default=NULL,
              help="model that is compared against the others", metavar="character"),
  make_option("--modelone", type="character", default=NULL,
              help="a model for comparison", metavar="character"),
  make_option("--modeltwo", type="character", default=NULL,
              help="a second model for comparison", metavar="character"),
  make_option("--modelthree", type="character", default=NULL,
              help="a third model for comparison", metavar="character"),
  make_option("--modelfour", type="character", default=NULL,
              help="a fourth model for comparison", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./",
              help="output saved in current working directory [default = %default]", metavar="character"),
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="input data", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$modelone) | is.null(opt$focus)) {
  print_help(opt_parser)
  stop("At least two arguments must be supplied (the focus model and one model for comparison).\n", call.=FALSE)
}

########################################################
###### approximate leave-one-out cross-validation ######
# http://mc-stan.org/loo/articles/loo2-with-rstan.html #
########################################################
# CZ_data = read.csv(opt$data, sep = ";")
# mod1 = readRDS("models/gaus_skew/gaus_skew.rds")
cat("Reading model", opt$focus, "\n")
focus = readRDS(opt$focus)
cat("Reading model", opt$modelone, "\n")
mod1 = readRDS(opt$modelone)
# mod2 = readRDS("models/gaus_skew/gaus_skew_hier_BCDG_shore.rds")
# cat("Reading model", opt$modeltwo, "\n")
# mod2 = readRDS(opt$modeltwo)
# cat("Reading model", opt$modelthree, "\n")
# mod3 = readRDS(opt$modelthree)
# cat("Reading model", opt$modelfour, "\n")
# mod4 = readRDS(opt$modelfour)

# mod_str = c(opt$focus, opt$modelone, opt$modeltwo, opt$modelthree, opt$modelfour)
mod_str = c(opt$focus, opt$modelone)
out_comp_str = lapply(mod_str, function(x) {
  modsplit = strsplit(strsplit(basename(x), "[.]")[[1]][1], split = "_")[[1]]
  modstr = modsplit[(length(modsplit)-1):length(modsplit)]
  return(paste0(modstr[1], modstr[2]))
})

cat("Extracting pointwise log-likelihood and compute LOO", "\n")
# mod_ls = list(focus, mod1, mod2, mod3, mod4)
mod_ls = list(focus, mod1)
loglik_ls = parallel::mclapply(seq_along(mod_ls), function(m) {extract_log_lik(mod_ls[[m]], merge_chains = FALSE)}, mc.cores = length(mod_ls))

cat("Providing relative effective sample sizes", "\n")
reff_ls = parallel::mclapply(seq_along(mod_ls), function(m) {relative_eff(exp(loglik_ls[[m]]))}, mc.cores = length(mod_ls))
loo_ls = parallel::mclapply(seq_along(mod_ls), function(m) {loo(loglik_ls[[m]], r_eff = reff_ls[[m]], cores = 2, save_psis = TRUE)}, mc.cores = length(mod_ls))
# lapply(seq_along(mod_ls), function(m) {
#   cat("LOO results for model", basename(mod_str[m]), "\n")
#   print(loo_ls[[m]])
# })

# split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
# out_comp = "tables/mods_comp/comp_skew_shore.csv"
# inp_mod1 = "models/gaus_skew/gaus_skew.rds"
lapply(seq_along(mod_ls), function(m) {
  cat("LOO results for model", basename(mod_str[m]), "\n")
  print(loo_ls[[m]])
  out_loo_ls = paste0(opt$out, "loo_", out_comp_str[[m]], ".csv")
  cat("Writing output", out_loo_ls, "...\n")
  write.table(round(loo_ls[[m]]$estimates, 2), out_loo_ls, sep = ",", row.names = TRUE, col.names = TRUE)
})

# plot(loo_1)
# plot(loo_2, label_points = TRUE)

cat("Comparing focus model", opt$focus, "against the rest ...\n")
comp_ls = parallel::mclapply(2:length(mod_ls), function(m) {compare(loo_ls[[1]], loo_ls[[m]])}, mc.cores = length(mod_ls)-1)
# print(comp_ls)
lapply(seq_along(comp_ls), function(c) {
  cat("Comparing", out_comp_str[[1]], "vs", out_comp_str[[c+1]], "...\n")
  print(comp_ls[[c]])
  write.csv(tibble(x=c("elpd_diff","SE"), comp_ls[[c]]),
            paste0(opt$out, "comp_", out_comp_str[[1]], "_", out_comp_str[[c+1]], ".csv"), row.names=FALSE)
})

# cat("The following comparisons have been completed:\n",
#     out_comp_str[[1]], "vs", out_comp_str[[2]], "\n",
#     out_comp_str[[1]], "vs", out_comp_str[[3]], "\n",
#     out_comp_str[[1]], "vs", out_comp_str[[4]], "\n",
#     out_comp_str[[1]], "vs", out_comp_str[[5]], "\n")

cat("The following comparisons have been completed:\n",
    out_comp_str[[1]], "vs", out_comp_str[[2]], "\n")
