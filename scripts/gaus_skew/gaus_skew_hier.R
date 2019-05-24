# Usage example:
# Rscript L.saxatilis-Mate-choice/scripts/gaus_skew/gaus_skew_hier.R -d data/CZ_all_mating_clean.csv
#   -s L.saxatilis-Mate-choice/scripts/gaus_skew/gaus_skew_hier_matrix.stan -i 8000 -c 4 -o gaus_skew/BCDG_shore/gaus_skew_hier_BCDG_shore

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
  make_option(c("-o", "--output"), type = "character", default = "output",
              help = "prefix for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$stanfile) | is.null(opt$iterations)) {
  print_help(opt_parser)
  stop("At least three arguments must be supplied (input data, stan file and MCMC iterations).\n", call.=FALSE)
}


CZ_data = read.csv(opt$data, sep = ";")
pref_out = opt$output

#########################################
# stan hierarchical model with skewness #
#########################################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE) - 15)

# names(CZ_data)
# head(CZ_data$ref_eco_sex)
# head(CZ_data$ref_ecotype)
# head(CZ_data$test_sex)

CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)
# CZ_matrix = model.matrix(mountYNcontact ~ shore, data = CZ_data)
# dM_matrix = CZ_matrix[, -1]
# CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype, data = CZ_data)


# CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex:shape, data = CZ_data)
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

skew_hier = rstan::stan(file = opt$stanfile, data = dat, iter = opt$iterations, warmup = opt$iterations/4,
                        chains=opt$chains, refresh=opt$iterations,
                        control = list(adapt_delta = 0.95, max_treedepth = 15))

cat("Saving", basename(pref_out), "Stan model file ...\n")

saveRDS(skew_hier, paste0("models/", pref_out, ".rds"))
# saveRDS(gaus_skew, "models/gaus_skew/gaus_skew_hier_BCDG_shore_eco.rds")

##########################################
# stan skewed hierarchical post analysis #
##########################################

if (sum(grepl(pattern = "b_intercept|d_intercept", skew_hier@model_pars)) == 2) {
  hier_pars = skew_hier@model_pars[1:6]
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
row.names(hier_pars_tbl)[grepl("d_coeff", row.names(hier_pars_tbl))] = paste0("d_", colnames(hier_matrix))
row.names(hier_pars_tbl)[grepl("g_coeff", row.names(hier_pars_tbl))] = paste0("g_", colnames(hier_matrix))

(hier_pars_df = rownames_to_column(as.data.frame(hier_pars_tbl), var="parameter"))
hier_pars_df$parameter[grepl("ntercept", hier_pars_df$parameter)] = c("b_intercept", "c_intercept", "d_intercept", "g_intercept")

cat("Writing parameter descriptive stats in", paste0("tables/", pref_out, "_coef.csv"), "...\n")
write.table(hier_pars_df, paste0("tables/", pref_out, "_coef.csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")


hier_draws = rstan::extract(skew_hier)

hier_hyp = hier_pars[grepl("coeff", hier_pars)]
hier_coeff = hier_pars_df$parameter
coeff_list = split(hier_coeff, ceiling(seq_along(hier_coeff)/length(colnames(hier_matrix))))
# names(coeff_list) = hier_pars
names(coeff_list) = hier_pars[grepl("coeff", hier_pars)]
# hier_draws$b_coeff = cbind(hier_draws$b_intercept, hier_draws$b_coeff)
# hier_draws$d_coeff = cbind(hier_draws$d_intercept, hier_draws$d_coeff)

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
    cat("Saving density plot", paste0("figures/", pref_out, "_", coeff_list[[x]][y], "_dens.png"), "...\n")
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
  cat("Saving dot plot", paste0("figures/", pref_out, "_", prx_hyp[x], "plot.png"), "...\n")
  ggsave(filename = paste0("figures/", pref_out, "_", prx_hyp[x], "plot.png"),
         plot = parfig_plot[[x]])
})

color_scheme_set("mix-blue-red")
parfig_trace = lapply(prx_hyp, function(hyp) {
  mcmc_trace(posterior, pars = hier_coeff[grepl(hyp, x = hier_coeff)],
             facet_args = list(ncol = 1, strip.position = "left"))
})
lapply(seq_along(prx_hyp), function(x) {
  cat("Saving trace plot", paste0("figures/", pref_out, "_", prx_hyp[x], "trace.png"), "...\n")
  ggsave(filename = paste0("figures/", pref_out, "_", prx_hyp[x], "trace.png"),
         plot = parfig_trace[[x]])
})
