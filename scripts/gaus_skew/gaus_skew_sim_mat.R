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


# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ",")
# skew_pars = read.csv("tables/gaus_skew/SKEW/gaus_skew_params.csv", sep = ";")
# skew_pars = column_to_rownames(skew_pars, var = "parameter")
# colnames(CZ_data)

# CZ_cline_params = read.csv("tables/clines/CZ_cline_params.csv", row.names = 1)

cline_2c3s <- function(position, sex, cl, cr, wl, wr, crab, wave, zs_c, zs_w, sc, sh, sw) {
  # wl = exp(lwl)
  # wr = exp(lwr)
  # sc = exp(lsc)
  # sh = exp(lsh)
  # sw = exp(lsw)
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]
  s_xl <- sqrt(sc^2 + 4*p_xl*(1-p_xl)*sh^2 + (p_xl^2)*(sw^2-sc^2))
  
  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing 
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  # combined cline
  z_x[z_x > z_xl] <- z_xl[z_x > z_xl]
  s_x[z_x > z_xl] <- s_xl[z_x > z_xl]
  phen_cline = cbind(z_x, s_x, sex, position)
  return(phen_cline)
}
phen_cline = cline_2c3s(position = CZ_data$LCmeanDist[CZ_data$shore=="CZA"], sex = CZ_data$test_sex[CZ_data$shore=="CZA"],
                        cl = CZ_cline_params["cl", "CZA"], cr = CZ_cline_params["cr", "CZA"],
                        wl = CZ_cline_params["lwl", "CZA"], wr = CZ_cline_params["lwr", "CZA"],
                        crab = CZ_cline_params["crab", "CZA"], wave = CZ_cline_params["wave", "CZA"],
                        zs_c = CZ_cline_params["zs_c", "CZA"], zs_w = CZ_cline_params["zs_w", "CZA"],
                        sc = CZ_cline_params["sc", "CZA"], sh = CZ_cline_params["sh", "CZA"],
                        sw = CZ_cline_params["sw", "CZA"])
head(phen_cline)
tail(phen_cline)
phen_cline = as.data.frame(cbind(phen_cline, log_len=log(CZ_data[CZ_data$shore=="CZA", ]$length_mm)))

# plot(CZ_data$DistAlongPath[CZ_data$shore=="CZA"], z_xl)
ggplot(data = phen_cline, aes(position, log_len, col=as.factor(sex))) +
  geom_ribbon(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x), alpha=0.15) +
  geom_point() +
  geom_line(aes(position, z_x), size=2, alpha=0.6)
  # geom_errorbar(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x), width=0.25)

# with(data = CZ_data[CZ_data$shore=="CZA", ], plot(LCmeanDist, log(length_mm), pch=19, col=test_sex))
# points(CZ_data$LCmeanDist[CZ_data$shore=="CZA"], phen_cline[, 'z_x'])

islands = as.character(unique(CZ_data$shore))
CZs_phen_cline = lapply(islands, function(x) {
  phen_cline = cline_2c3s(position = CZ_data$LCmeanDist[CZ_data$shore==x], sex = CZ_data$test_sex[CZ_data$shore==x],
                          cl = CZ_cline_params["cl", x], cr = CZ_cline_params["cr", x],
                          wl = CZ_cline_params["lwl", x], wr = CZ_cline_params["lwr", x],
                          crab = CZ_cline_params["crab", x], wave = CZ_cline_params["wave", x],
                          zs_c = CZ_cline_params["zs_c", x], zs_w = CZ_cline_params["zs_w", x],
                          sc = CZ_cline_params["sc", x], sh = CZ_cline_params["sh", x],
                          sw = CZ_cline_params["sw", x])
  phen_cline = as.data.frame(cbind(phen_cline, log_len=log(CZ_data[CZ_data$shore==x, ]$length_mm)))
  return(phen_cline)
})

# CZ_cline_params["cl", "CZA"]-CZ_cline_params["lwl", "CZA"]/2

# CZ_sim[CZ_sim$position==round(cl-wl/2), ]

CZs_phen_cline = lapply(CZs_phen_cline, function(x) {
  x[, "sex"] = ifelse(x[, "sex"]==1, "female", "male")
  return(x)
})
lapply(seq_along(islands), function(x) {
  cline_pos = CZ_cline_params["cl", islands[x]]-CZ_cline_params["lwl", islands[x]]/2
  CZ_sim_sex = split(CZs_phen_cline[[x]], CZs_phen_cline[[x]][, "sex"])
  write.table(data.frame(male=as.numeric(), female=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
              file = paste0("tables/gaus_skew/SKEW/sims/", islands[x], "_", round(cline_pos), "_sim_YN.csv"), sep = ",",
              row.names = FALSE, col.names = TRUE)
  sim_mat(data = CZ_sim_sex, pos = cline_pos, isl = islands[x])
})


# phen_cline$sex = ifelse(phen_cline$sex==1, "female", "male")
CZ_sim_sex = split(phen_cline, phen_cline$sex)

# pos = (CZ_cline_params["cl", "CZA"]-CZ_cline_params["lwl", "CZA"]/2) - 1
# sum(round(CZ_sim_sex$female$position)==round(pos))
# sum(CZ_sim_sex$female$s_x > 0.5)
# sum(CZ_sim_sex$male$s_x > 0.5)
# rm(pos)

write.table(data.frame(male=as.numeric(), female=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
            file = "tables/gaus_skew/SKEW/sims/CZ_sim_YN.csv", sep = ",",
            row.names = FALSE, col.names = TRUE)
sim_mat = function(data, pos, isl) {
  bar = list()
  YN = data.frame()
  
  smp_pos = seq(round(pos)-2,round(pos)+2)
  fml_df = rbindlist(lapply(seq_along(smp_pos), function(z) {
    data$female[round(data$female$position)==smp_pos[z], ]
  }))
  ml_df = rbindlist(lapply(seq_along(smp_pos), function(z) {
    data$male[round(data$male$position)==smp_pos[z], ]
  }))
  
  # fml_m = mean(CZ_sim_sex$female[round(CZ_sim_sex$female$position)==seq(round(pos)-1,round(pos)+1), ]$z_x)
  # fml_sd = mean(CZ_sim_sex$female[round(CZ_sim_sex$female$position)==c(round(pos)-1,round(pos),round(pos)+1), ]$s_x)
  fml_m = mean(fml_df$z_x)
  fml_sd = mean(fml_df$s_x)
  ml_m = mean(ml_df$z_x)
  ml_sd = mean(ml_df$s_x)
  
  if (fml_sd > 0.4) {
    fml_dtr = rnorm(n = 10, mean = fml_m, sd = log(1.5))
  } else {
    fml_dtr = rnorm(n = 10, mean = fml_m, sd = fml_sd)
  }
  if (ml_sd > 0.4) {
    ml_dtr = rnorm(n = 10, mean = ml_m, sd = log(1.5))
  } else {
    ml_dtr = rnorm(n = 10, mean = ml_m, sd = ml_sd)
  }
  
  for (f in seq_along(fml_dtr)) {
    success=FALSE
    i=1
    fem = fml_dtr[f]
    while (!success) {
      m = sample(ml_dtr, 1, replace = FALSE)
      p = skew_pars["level","mean"] + skew_pars["scale","mean"] * exp(-0.5 * (((fem - m) - skew_pars["centre","mean"])
                                                                              / skew_pars["choosiness","mean"])^2) * 
        (1 + erf(skew_pars["asymmetry","mean"] * ((fem - m) - skew_pars["centre","mean"]) / (1.414214 * skew_pars["choosiness","mean"])))
      s = rbinom(n = 1, size = 1, prob = p)
      YN[i,'male'] = m
      YN[i,'female'] = fem
      YN[i,'mountYN'] = s
      success = (s > 0)
      i = i + 1
    }
    write.table(YN, file = paste0("tables/gaus_skew/SKEW/sims/", isl, "_", round(pos), "_sim_YN.csv"), append = TRUE,
                sep = ",", row.names = FALSE, col.names = FALSE)
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}

sim_mat(pos = (CZ_cline_params["cl", "CZA"]-CZ_cline_params["lwl", "CZA"]/2) - 1)

res = lapply(names(fem), function(x) {
  sapply(names(fem$crab), function(y) sim_mat(female = fem[[x]][[y]], male = mal[[x]][[y]]))
})