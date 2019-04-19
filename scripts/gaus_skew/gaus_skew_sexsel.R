rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "tibble", "boot", "bayesplot", "Rmisc", "pander",
              "bbmle", "loo", "ggpubr", "cowplot", "purrr", "reshape2", "gridExtra", "grid",
              "arm", "parallel","optparse", "pracma")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="input data", metavar="character"),
  make_option(c("-s", "--modelfile"), type="character", default=NULL, 
              help="rds output for model written in Stan", metavar="character")) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$modelfile)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input data and stan file).\n", call.=FALSE)
}


CZ_data = read.csv(opt$data, sep = ";")


#############################
# plot post distr gaus skew #
#############################
gaus_skew = readRDS(opt$modelfile)
# gaus_skew = readRDS("models/gaus_skew/gaus_skew_hier_BCDG.rds")
# gaus_size_pars = c("level","scale","centre","choosiness","asymmetry")


list_of_draws = rstan::extract(gaus_skew)
# print(names(list_of_draws))
names(list_of_draws)[grepl("scale|preference|choosiness|asymmetry", names(list_of_draws))] = c("b_preds", "b", "c", "d", "g")
# names(list_of_draws) = c("level","scale","centre","choosiness","asymmetry", "y_hat", "lp__")
skew_hier_pars = c("b", "c", "d", "g")

skew_hier_distr = lapply(skew_hier_pars, function(x) {
  apply(list_of_draws[[x]], 2, mean)
})
names(skew_hier_distr) = skew_hier_pars
lapply(skew_hier_distr, mean)
lapply(skew_hier_pars, function(x) {
  xx = seq(min(skew_hier_distr[[x]]), max(skew_hier_distr[[x]]), by = 0.001)
  dd = density(skew_hier_distr[[x]])
  # interpolating function
  fx = splinefun(dd$x, dd$y)
  # normalize so prob >0
  pxx = pmax(0, fx(xx))
  # sample from the "empirical" distribution
  samp = sample(xx, 1e5, replace = TRUE, prob = pxx)
  # and take sample quantiles
  quantile(samp, c(0.025, 0.975))
})


parfig = lapply(skew_hier_pars, function(x) {
  ggplot() +
    geom_density(aes(skew_hier_distr[[x]]), fill='red', col='black') +
    labs(x="", title = x) +
    theme(axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold",size = 15))
})
names(parfig) = skew_hier_pars


# find derivative of skew-normal in Wolphram|Alpha typing
# d[l + s * exp(-0.5 * ((x - p) / c)^2) * (1 + ERF(a * ((x - p) / (1.414214 * c)))) ,x]

opt = function(b, c, d, g, x){
  # derivative
  (b * exp(-(0.5 * (c - x)^2)/d^2) * (0.797884 * g * d * exp(-(0.5 * g^2 * (c - x)^2)/d^2) -
                                        (x - c) * (erf((0.707107 * g * (x - c))/d) + 1)))/d^2
}
pars = list(b = skew_hier_distr$b, c = skew_hier_distr$c, d = skew_hier_distr$d, g = skew_hier_distr$g)

# find the root of the derivative
#str(xmin <- uniroot(opt, c(0, 1), tol = 0.0001, scale = 0.44, preference = -0.07, asymmetry = 1.35, choosiness = 0.67))


opt_draws = sapply(1:3709, function(z) {
  uniroot(opt, c(-1, 1), tol = 0.0001, b = pars[['b']][z], c = pars[['c']][z], d = pars[['d']][z],
          g = pars[['g']][z])$root
})

skew_hier_distr$o = opt_draws

opt_dens = ggplot() +
  geom_density(aes(skew_hier_distr$o), fill='red', col='black') +
  labs(x="", title = "o") +
  theme(axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 15))
# mean(list_of_draws$preference)
parfig$o = opt_dens


pdf("figures/gaus_skew/BCDG/gaus_skew_hier_BCDG_hyp_dens.pdf",width = 10, height = 7)
#do.call(ggarrange, parfig)
ggarrange(parfig$b, parfig$c, parfig$d, parfig$g, parfig$o)
dev.off()


############################################
# save Stan output for gaus_skew_hier_BCDG #
############################################
CZ_data$stan_yhat = summary(gaus_skew, pars = c("y_hat"))$summary[,'mean']
CZ_data$stan_yhat_uci = summary(gaus_skew, pars = c("y_hat"))$summary[,'97.5%']
CZ_data$stan_yhat_lci = summary(gaus_skew, pars = c("y_hat"))$summary[,'2.5%']
CZ_data$stan_mountYN = rbinom(n = nrow(CZ_data),size = 1,prob = CZ_data$stan_yhat)
write.table(CZ_data, "tables/gaus_skew/BCDG/gaus_skew_hier_BCDG_mat.csv", row.names = FALSE, col.names = TRUE, sep = ",")



###################################
# apply skew model to field distr #
###################################

# generate size distributions for each sex and ecotype using cline parameters
# s_xl <- sqrt(sc^2 + 4*0.5*(1-0.5)*sh^2 + (0.5^2)*(sw^2-sc^2))    # standard deviation at the zone centre
# CZ_cline_params = read.csv("tables/clines/CZ_cline_params.csv", sep = ";")
CZ_cline_params = read.csv("tables/clines/CZ_cline_params_mm_mf.csv")

# s_centre = function(sc, sh, sw) {sqrt(sc^2 + sh^2 + 0.25*(sw^2-sc^2))}

male_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = log(x[5]), sd = abs(x[13]))),2)
#apply(male_c, 2, hist)
female_c = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = log(x[10]),
                                                               sd = abs(x[13]))),2)
#apply(female_c, 2, hist)


male_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = log(x[6]), sd = abs(x[14]))),2)
#apply(male_h, 2, hist)
female_h = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = log(x[11]),
                                                               sd = abs(x[14]))),2)
#apply(female_h, 2, hist)

male_w = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = log(x[7]), sd = abs(x[15]))),2)
#apply(male_w, 2, hist)
female_w = round(sapply(CZ_cline_params[-1], function(x) rnorm(n = 1000, mean = log(x[12]),
                                                               sd = abs(x[15]))),2)
#apply(female_w, 2, hist)

# pair each female with every male within habitat and contact zone and compute mounting success YN
# skew_params = read.csv("tables/gaus_skew/gaus_skew_params.csv", sep = ";")


fem=list(crab=data.frame(female_c), hybrid=data.frame(female_h), wave=data.frame(female_w))
mal=list(crab=data.frame(male_c), hybrid=data.frame(male_h), wave=data.frame(male_w))

sim_mat = function(female, male) {
  bar = list()
  YN = data.frame()
  for (f in seq_along(female)) {
    success=FALSE
    i=1
    fem = female[f]
    while (!success) {
      m = sample(male, 1, replace = FALSE)
      p = mean(skew_hier_distr$b) * exp(-0.5 * (((fem - m) - mean(skew_hier_distr$c)) /
                                                  mean(skew_hier_distr$d))^2) *
        (1 + erf(mean(skew_hier_distr$g) * ((fem - m) - mean(skew_hier_distr$c)) /
                   (1.414214 * mean(skew_hier_distr$d))))
      s = rbinom(n = 1, size = 1, prob = p)
      YN[i,'male'] = m
      YN[i,'female'] = fem
      YN[i,'mountYN'] = s
      success = (s > 0)
      i = i + 1
    }
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}


res = lapply(names(fem), function(x) {
  sapply(names(fem$crab), function(y) sim_mat(female = fem[[x]][[y]], male = mal[[x]][[y]]))
})

shore = c("CZA", "CZB", "CZC", "CZD")
ecotype = c("crab", "hybrid", "wave")

eco_CZ = lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) do.call(rbind, res[[x]][,y]))
})

# save the simulated datasets of mate choice
lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) write.table(eco_CZ[[x]][[y]],
                                                   paste0("tables/gaus_skew/BCDG/sims/", shore[y], "_", ecotype[x], "_sim_YN.csv"),
                                                   row.names = FALSE, col.names = TRUE, sep = ";"))
})

