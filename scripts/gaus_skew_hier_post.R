rm(list = ls())

#############################
# Install required packages #
#############################

# List of packages for session
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "Rmisc",
              "purrr", "reshape2", "gridExtra", "grid", "parallel", "projpred", "optparse")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="input data", metavar="character"),
  make_option(c("-m", "--modelfile"), type="character", default=NULL, 
              help="rds model format written in Stan", metavar="character"),
  make_option(c("-o", "--output"), type = "character", default = "output",
              help = "prefix for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$modelfile)) {
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input data and input model.rds).\n", call.=FALSE)
}


CZ_data = read.csv(opt$data, sep = ";")
skew_hier = readRDS(opt$modelfile)
# "gaus_skew/gaus_skew_hier_BCDG"
pref_out = opt$output

#########################################
# stan hierarchical model with skewness #
#########################################
if ("d_intercept" %in% skew_hier@model_pars) {
  hier_pars = skew_hier@model_pars[1:5]
} else {
  hier_pars = skew_hier@model_pars[1:4]
}

(hier_pars_tbl = round(summary(skew_hier, pars = hier_pars, probs=c(0.025, 0.975))$summary,2))

if (length(row.names(hier_pars_tbl))==16) {
  hier_matrix = model.matrix(mountYNcontact ~ shore, data = CZ_data)
} else {
  hier_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)
}

row.names(hier_pars_tbl)[grepl("b_coeff", row.names(hier_pars_tbl))] = paste0("b_", colnames(hier_matrix))
row.names(hier_pars_tbl)[grepl("c_coeff", row.names(hier_pars_tbl))] = paste0("c_", colnames(hier_matrix))
row.names(hier_pars_tbl)[grepl("d_coeff", row.names(hier_pars_tbl))] = paste0("d_", colnames(hier_matrix)[-1])
row.names(hier_pars_tbl)[grepl("g_coeff", row.names(hier_pars_tbl))] = paste0("g_", colnames(hier_matrix))

(hier_pars_df = rownames_to_column(as.data.frame(hier_pars_tbl), var="parameter"))
hier_pars_df$parameter[grepl("ntercept", hier_pars_df$parameter)] = c("b_intercept", "c_intercept", "d_intercept", "g_intercept")

write.table(hier_pars_df, paste0("tables/", pref_out, "_coef.csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")


hier_draws = rstan::extract(skew_hier)

hier_hyp = hier_pars[grepl("coeff", hier_pars)]
hier_coeff = hier_pars_df$parameter
coeff_list = split(hier_coeff, ceiling(seq_along(hier_coeff)/length(colnames(hier_matrix))))
# names(coeff_list) = hier_pars
names(coeff_list) = hier_pars[grepl("coeff", hier_pars)]
hier_draws$d_coeff = cbind(hier_draws$d_coeff, hier_draws$d_intercept)

coeff_draws = lapply(hier_hyp, function(x) {
  hier_draws[[x]]
})

rn_coeff_draws = lapply(seq_along(hier_hyp), function(x) {
  colnames(coeff_draws[[x]]) = coeff_list[[x]]
  coeff_draws[[x]]
})

coeff_parfig = lapply(seq_along(hier_hyp), function(x) {
  lapply(seq_along(colnames(hier_matrix)), function(y) {
    ggplot() +
      geom_density(aes(rn_coeff_draws[[x]][, y]), fill='red', col='black') +
      labs(x="", title = coeff_list[[x]][y]) +
      theme(axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  })
})

lapply(seq_along(hier_hyp), function(x) {
  lapply(seq_along(colnames(hier_matrix)), function(y) {
    ggsave(filename = paste0("figures/", pref_out, "_", coeff_list[[x]][y], "_dens.png"),
           plot = coeff_parfig[[x]][[y]])
  })
}) 


posterior = as.array(skew_hier)
# dim(posterior)
dimnames(posterior)$parameters[1:length(row.names(hier_pars_tbl))] = hier_coeff

# mcmc_intervals(posterior, pars = hier_coeff[grepl("b_", x = hier_coeff)], point_est = "mean")
# mcmc_intervals(posterior, pars = hier_coeff[grepl("c_", x = hier_coeff)], point_est = "mean")
# mcmc_intervals(posterior, pars = hier_coeff[grepl("d_", x = hier_coeff)], point_est = "mean")
# mcmc_intervals(posterior, pars = hier_coeff[grepl("g_", x = hier_coeff)], point_est = "mean")
prx_hyp = c("b_", "c_", "d_", "g_")

parfig_plot = lapply(prx_hyp, function(hyp) {
  mcmc_intervals(posterior, pars = hier_coeff[grepl(hyp, x = hier_coeff)], point_est = "mean",
                 prob_outer = 0.95)
})
lapply(seq_along(prx_hyp), function(x) {
  ggsave(filename = paste0("figures/", pref_out, "_", prx_hyp[x], "plot.png"),
         plot = parfig_plot[[x]])
})

color_scheme_set("mix-blue-red")
parfig_trace = lapply(prx_hyp, function(hyp) {
  mcmc_trace(posterior, pars = hier_coeff[grepl(hyp, x = hier_coeff)],
             facet_args = list(ncol = 1, strip.position = "left"))
})
lapply(seq_along(prx_hyp), function(x) {
  ggsave(filename = paste0("figures/", pref_out, "_", prx_hyp[x], "trace.png"),
         plot = parfig_trace[[x]])
})
