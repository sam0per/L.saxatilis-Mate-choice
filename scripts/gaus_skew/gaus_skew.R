rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc", "pander",
              "bbmle", "loo", "ggpubr", "cowplot", "purrr", "reshape2", "gridExtra", "grid", "arm", "parallel",
              "rstantools", "optparse")

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
# CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
# CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
# head(CZ_data)
# summary(CZ_data)


############################
# stan model with skewness #
############################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE) - 15)


dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio)
# str(dat)
# dat$posterior_predictive = 1
# fit.prior.pred = stan(file = "scripts/CZ_mating_gaus_prior_size.stan", data = dat)

gaus_skew = rstan::stan(file = opt$stanfile, data = dat, iter = 8000, warmup = 2000,
                        chains=4, refresh=8000,
                        control = list(adapt_delta = 0.90, max_treedepth = 15))
saveRDS(gaus_skew, "models/gaus_skew/gaus_skew.rds")
