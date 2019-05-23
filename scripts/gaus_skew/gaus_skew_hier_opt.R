# Usage example: Rscript scripts/gaus_skew/gaus_skew_hier_opt.R -d data/CZ_all_mating_clean.csv
#                  -s models/gaus_skew/gaus_skew_hier_BCDG.rds -o gaus_skew/gaus_skew_hier_BCDG

rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc",
              "purrr", "reshape2", "gridExtra", "grid", "arm", "parallel","optparse", "pracma")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="input data", metavar="character"),
  make_option(c("-s", "--modelfile"), type="character", default=NULL,
              help="rds output for model written in Stan", metavar="character"),
  make_option(c("-o", "--output"), type = "character", default = "output",
              help = "prefix for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$modelfile)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input data and stan file).\n", call.=FALSE)
}


CZ_data = read.csv(opt$data, sep = ";")
pref_out = opt$output

#############################
# plot post distr gaus skew #
#############################
gaus_skew = readRDS(opt$modelfile)
# gaus_skew = readRDS("models/gaus_skew/gaus_skew_hier_BCDG.rds")
# gaus_size_pars = c("level","scale","centre","choosiness","asymmetry")


list_of_draws = rstan::extract(gaus_skew)
names(list_of_draws)[grepl("scale|preference|choosiness|asymmetry", names(list_of_draws))] = c("b_preds", "b", "c", "d", "alpha")
skew_hier_pars = c("b", "c", "d", "alpha")

# find derivative of skew-normal in Wolphram|Alpha typing
# d[l + s * exp(-0.5 * ((x - p) / c)^2) * (1 + ERF(a * ((x - p) / (1.414214 * c)))) ,x]

opt = function(b, c, d, alpha, x) {
  # derivative
  (b * exp(-(0.5 * (c - x)^2)/d^2) * (0.797884 * alpha * d * exp(-(0.5 * alpha^2 * (c - x)^2)/d^2) -
                                        (x - c) * (erf((0.707107 * alpha * (x - c))/d) + 1)))/d^2
}
opt_mx = matrix(nrow = nrow(list_of_draws$b), ncol = ncol(list_of_draws$b))
for (ind in 1:ncol(list_of_draws$b)) {
  pars = list(b = list_of_draws$b[, ind], c = list_of_draws$c[, ind],
              d = list_of_draws$d[, ind], alpha = list_of_draws$alpha[, ind])
  opt_draws = sapply(1:nrow(list_of_draws$b), function(z) {
    uniroot(opt, c(-2, 2), tol = 0.001, b = pars[['b']][z], c = pars[['c']][z], d = pars[['d']][z],
            alpha = pars[['alpha']][z])$root
  })
  opt_mx[, ind] = opt_draws
}

write.table(opt_mx, paste0("tables/", pref_out, "_opt.csv"), row.names = FALSE, col.names = TRUE, sep = ",")
