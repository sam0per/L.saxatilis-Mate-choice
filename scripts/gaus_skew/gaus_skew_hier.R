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
              help="model written in Stan", metavar="character"),
  make_option(c("-i", "--iterations"), type = "integer", default = NULL,
              help = "number of MCMC iterations", metavar = "integer"),
  make_option(c("-c", "--chains"), type = "integer", default = 4,
              help = "number of MCMC chains [default: %default]", metavar = "integer"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$stanfile) | is.null(opt$iterations)) {
  print_help(opt_parser)
  stop("At least three arguments must be supplied (input data, stan file and MCMC iterations).\n", call.=FALSE)
}


CZ_data = read.csv(opt$data, sep = ";")

#########################################
# stan hierarchical model with skewness #
#########################################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE) - 15)

# names(CZ_data)
# head(CZ_data$ref_eco_sex)
# head(CZ_data$ref_ecotype)
# head(CZ_data$test_sex)
# CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_eco_sex + test_sex * shape, data = CZ_data)
CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex:shape, data = CZ_data)
# dim(CZ_matrix)

# CZ_data$ref_ecotype=as.integer(CZ_data$ref_ecotype) # 1 for crab and 2 for wave
# CZ_data$shore=as.integer(CZ_data$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
# CZ_data$test_sex=as.integer(CZ_data$test_sex) # 1 for female, 2 for male
# CZ_data$ref_eco_sex = as.integer(CZ_data$ref_eco_sex) # 1 for crab female, 2 for crab male, 3 for wave female, 4 for wave male
# dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio, shape = CZ_data$shape,
#            N_shore = length(unique(CZ_data$shore)), N_eco = length(unique(CZ_data$ref_eco_sex)),
#            N_sex = length(unique(CZ_data$test_sex)),
#            shore = CZ_data$shore, ref = CZ_data$ref_eco_sex, test = CZ_data$test_sex)

dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio,
           X = CZ_matrix, K = dim(CZ_matrix)[2])

gaus_skew = rstan::stan(file = opt$stanfile, data = dat, iter = opt$iterations, warmup = opt$iterations/4,
                        chains=opt$chains, refresh=opt$iterations,
                        control = list(adapt_delta = 0.90, max_treedepth = 15))

saveRDS(gaus_skew, "models/gaus_skew/gaus_skew_hier_CDG.rds")
