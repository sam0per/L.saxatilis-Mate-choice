# Usage example:
# Rscript L.saxatilis-Mate-choice/scripts/gaus_skew/gaus_skew_hier.R -d data/CZ_all_mating_clean.csv
#   -s L.saxatilis-Mate-choice/scripts/gaus_skew/gaus_skew_hier_matrix.stan -i 8000 -c 4 -p shore -m b_par
#   -o gaus_skew/BCDG_shore/gaus_skew_hier_BCDG_shore

rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "optparse", "tibble", "bayesplot")
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
              help = "number of MCMC chains [default: %default]", metavar = "integer"),
  make_option(c("-p", "--predictors"), type = "character", default = NULL,
              help = "select predictors of the model matrix [all, shore, ecotype, sex, islsex]", metavar = "character"),
  make_option(c("-m", "--modhyp"), type = "character", default = NULL,
              help = "choose which parameter to transform into hyper [b_par, c_par, d_par, alpha_par]", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "output",
              help = "prefix for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$stanfile) | is.null(opt$iterations)) {
  print_help(opt_parser)
  stop("At least three arguments must be supplied (input data, stan file and MCMC iterations).\n", call.=FALSE)
}

# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
CZ_data = read.csv(opt$data, sep = ";")
pref_out = opt$output
# pref_out = "gaus_skew/B_all/gaus_skew_hier_B_all"
pred_mx = opt$predictors

#########################################
# stan hierarchical model with skewness #
#########################################
rstan_options(auto_write = TRUE)
options(mc.cores = opt$chains)
# options(mc.cores = parallel::detectCores(logical = FALSE) - 15)
# options(mc.cores = parallel::detectCores(logical = FALSE) - 2)

if (pred_mx == "all") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
} else if (pred_mx == "shore") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore, data = CZ_data)[,-1]
} else if (pred_mx == "ecotype") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + shape, data = CZ_data)[,-1]
} else if (pred_mx == "sex") {
  CZ_matrix = matrix(model.matrix(mountYNcontact ~ test_sex, data = CZ_data)[,-1])
  colnames(CZ_matrix) = "test_sexmale"
} else if (pred_mx == "islsex") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore + test_sex, data = CZ_data)[,-1]
} else {
  print("Model matrix is missing")
}

dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio,
           X = CZ_matrix, K = dim(CZ_matrix)[2])
if (grepl(pattern = "alpha", x = opt$modhyp)) {
  hier_hyp_px = substr(opt$modhyp,1,6)
} else {
  hier_hyp_px = substr(opt$modhyp,1,2)
}
# hier_hyp_px = substr("b_par",1,2)
if (hier_hyp_px=="b_") {
  start_val = list(list(b_intercept=0.4, b_coeff=rep(0,ncol(CZ_matrix)), c_par=-0.17, d_par=0.85, alpha_par=2.32),
                   list(b_intercept=0.4, b_coeff=rep(0,ncol(CZ_matrix)), c_par=-0.17, d_par=0.85, alpha_par=2.32),
                   list(b_intercept=0.4, b_coeff=rep(0,ncol(CZ_matrix)), c_par=-0.17, d_par=0.85, alpha_par=2.32),
                   list(b_intercept=0.4, b_coeff=rep(0,ncol(CZ_matrix)), c_par=-0.17, d_par=0.85, alpha_par=2.32))
} else if (hier_hyp_px=="c_") {
  start_val = list(list(c_intercept=-0.17, c_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, d_par=0.85, alpha_par=2.32),
                   list(c_intercept=-0.17, c_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, d_par=0.85, alpha_par=2.32),
                   list(c_intercept=-0.17, c_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, d_par=0.85, alpha_par=2.32),
                   list(c_intercept=-0.17, c_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, d_par=0.85, alpha_par=2.32))
} else if (hier_hyp_px=="d_") {
  start_val = list(list(d_intercept=0.85, d_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, alpha_par=2.32),
                   list(d_intercept=0.85, d_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, alpha_par=2.32),
                   list(d_intercept=0.85, d_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, alpha_par=2.32),
                   list(d_intercept=0.85, d_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, alpha_par=2.32))
} else {
  start_val = list(list(alpha_intercept=2.32, alpha_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, d_par=0.85),
                   list(alpha_intercept=2.32, alpha_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, d_par=0.85),
                   list(alpha_intercept=2.32, alpha_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, d_par=0.85),
                   list(alpha_intercept=2.32, alpha_coeff=rep(0,ncol(CZ_matrix)), b_par=0.4, c_par=-0.17, d_par=0.85))
}


skew_hier = rstan::stan(file = opt$stanfile, data = dat, iter = opt$iterations, warmup = opt$iterations/4,
                        chains=opt$chains, refresh=opt$iterations, init = start_val,
                        control = list(adapt_delta = 0.95, max_treedepth = 15))

# skew_hier = rstan::stan(file = "../L.saxatilis-Mate-choice/scripts/gaus_skew/gaus_skew_hier_interc_matrix.stan",
#                         data = dat, iter = 8000, warmup = 2000,
#                         chains=4, refresh=8000, init = start_val,
#                         control = list(adapt_delta = 0.95, max_treedepth = 15))


cat("Saving", basename(pref_out), "Stan model file ...\n")

saveRDS(skew_hier, paste0("models/", pref_out, ".rds"))
# saveRDS(gaus_skew, "models/gaus_skew/gaus_skew_hier_BCDG_shore_eco.rds")

##########################################
# stan skewed hierarchical post analysis #
##########################################

if (sum(grepl(pattern = "intercept", skew_hier@model_pars)) == 1) {
  hier_pars = skew_hier@model_pars[1:5]
} else {
  hier_pars = skew_hier@model_pars[1:4]
}

hier_pars_tbl = round(summary(skew_hier, pars = hier_pars, probs=c(0.025, 0.975))$summary,2)

row.names(hier_pars_tbl)[grepl("b_coeff", row.names(hier_pars_tbl))] = paste0("b_", colnames(CZ_matrix))
row.names(hier_pars_tbl)[grepl("c_coeff", row.names(hier_pars_tbl))] = paste0("c_", colnames(CZ_matrix))
row.names(hier_pars_tbl)[grepl("d_coeff", row.names(hier_pars_tbl))] = paste0("d_", colnames(CZ_matrix))
row.names(hier_pars_tbl)[grepl("alpha_coeff", row.names(hier_pars_tbl))] = paste0("alpha_", colnames(CZ_matrix))

(hier_pars_df = rownames_to_column(as.data.frame(hier_pars_tbl), var="parameter"))
# hier_pars_df$parameter[grepl("ntercept", hier_pars_df$parameter)] = c("b_intercept", "c_intercept", "d_intercept", "g_intercept")

cat("Writing parameter descriptive stats in", paste0("tables/", pref_out, "_coef.csv"), "...\n")
write.table(hier_pars_df, paste0("tables/", pref_out, "_coef.csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")

hier_draws = rstan::extract(skew_hier)

hier_hyp = hier_pars[grepl(pattern = hier_hyp_px, hier_pars)]
hier_nnhyp = hier_pars[!grepl(pattern = paste0(hier_hyp_px, "coeff"), hier_pars)]
hier_coeff = hier_pars_df$parameter
hyp_draws = lapply(paste0(hier_hyp_px, "coeff"), function(x) {
  hier_draws[[x]]
})

hyp_draws_dens = lapply(seq_along(colnames(CZ_matrix)), function(y) {
  ggplot() +
    geom_density(aes(x = hyp_draws[[1]][, y], y = ..scaled..), fill='red', col='black') +
    labs(x="", title = paste0(hier_hyp_px, colnames(CZ_matrix)[y])) +
    theme(axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
})

# coeff_list = split(hier_coeff, ceiling(seq_along(hier_coeff)/length(colnames(CZ_matrix))))
# names(coeff_list) = hier_pars
# names(coeff_list) = hier_pars[grepl("coeff", hier_pars)]
# hier_draws$b_coeff = cbind(hier_draws$b_intercept, hier_draws$b_coeff)
# hier_draws$d_coeff = cbind(hier_draws$d_intercept, hier_draws$d_coeff)

# coeff_draws = lapply(hier_hyp, function(x) {
#   hier_draws[[x]]
# })

# rn_coeff_draws = lapply(seq_along(hier_hyp), function(x) {
#   colnames(coeff_draws[[x]]) = coeff_list[[x]]
#   coeff_draws[[x]]
# })

coeff_parfig = lapply(hier_nnhyp, function(x) {
  ggplot() +
    geom_density(aes(x = hier_draws[[x]], y = ..scaled..), fill='red', col='black') +
    labs(x="", title = x) +
    theme(axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
})


lapply(seq_along(colnames(CZ_matrix)), function(y) {
  cat("Saving density plot", paste0("figures/", pref_out, "_", hier_hyp_px, colnames(CZ_matrix)[y], "_dens.png"), "...\n")
  ggsave(filename = paste0("figures/", pref_out, "_", hier_hyp_px, colnames(CZ_matrix)[y], "_dens.png"),
         plot = hyp_draws_dens[[y]])
})
lapply(seq_along(hier_nnhyp), function(y) {
  cat("Saving density plot", paste0("figures/", pref_out, "_", hier_nnhyp[y], "_dens.png"), "...\n")
  ggsave(filename = paste0("figures/", pref_out, "_", hier_nnhyp[y], "_dens.png"),
         plot = coeff_parfig[[y]])
})


posterior = as.array(skew_hier)
# dim(posterior)
dimnames(posterior)$parameters[1:length(row.names(hier_pars_tbl))] = hier_coeff

# mcmc_intervals(posterior, pars = hier_coeff[grepl("b_", x = hier_coeff)], point_est = "mean")
# mcmc_intervals(posterior, pars = hier_coeff[grepl("c_", x = hier_coeff)], point_est = "mean")
# mcmc_intervals(posterior, pars = hier_coeff[grepl("d_", x = hier_coeff)], point_est = "mean")
# mcmc_intervals(posterior, pars = hier_coeff[grepl("g_", x = hier_coeff)], point_est = "mean")
# prx_hyp = c("b_", "c_", "d_", "alpha_")
parfig_plot = mcmc_intervals(posterior, pars = hier_coeff, point_est = "mean", prob_outer = 0.95)
cat("Saving dot plot", paste0("figures/", pref_out, "_plot.png"), "...\n")
ggsave(filename = paste0("figures/", pref_out, "_plot.png"), plot = parfig_plot)

# parfig_plot = lapply(prx_hyp, function(hyp) {
#   mcmc_intervals(posterior, pars = hier_coeff[grepl(hyp, x = hier_coeff)], point_est = "mean",
#                  prob_outer = 0.95)
# })
# lapply(seq_along(prx_hyp), function(x) {
#   cat("Saving dot plot", paste0("figures/", pref_out, "_", prx_hyp[x], "plot.png"), "...\n")
#   ggsave(filename = paste0("figures/", pref_out, "_", prx_hyp[x], "plot.png"),
#          plot = parfig_plot[[x]])
# })

# color_scheme_set("mix-blue-red")
# parfig_trace = lapply(prx_hyp, function(hyp) {
#   mcmc_trace(posterior, pars = hier_coeff[grepl(hyp, x = hier_coeff)],
#              facet_args = list(ncol = 1, strip.position = "left"))
# })
# lapply(seq_along(prx_hyp), function(x) {
#   cat("Saving trace plot", paste0("figures/", pref_out, "_", prx_hyp[x], "trace.png"), "...\n")
#   ggsave(filename = paste0("figures/", pref_out, "_", prx_hyp[x], "trace.png"),
#          plot = parfig_trace[[x]])
# })
