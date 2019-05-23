rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "optparse", "tibble", "bayesplot", "data.table", "purrr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="input data", metavar="character"),
  make_option(c("-m", "--modelpars"), type="character", default=NULL,
              help="mean estimates of the inferred parameters", metavar="character"),
  make_option(c("-c", "--clinepars"), type = "character", default = NULL,
              help = "cline parameter estimates", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "output_sim",
              help = "prefix for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$stanfile) | is.null(opt$iterations)) {
  print_help(opt_parser)
  stop("At least three arguments must be supplied (input data, stan file and MCMC iterations).\n", call.=FALSE)
}


CZ_data = read.csv(opt$data, sep = ",")
skew_pars = read.csv(opt$modelpars, sep = ";")
skew_pars = column_to_rownames(skew_pars, var = "parameter")
CZ_cline_params = read.csv(opt$clinepars, row.names = 1)
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

islands = as.character(unique(CZ_data$shore))
CZs_phen_cline = lapply(islands, function(x) {
  phen_cline = cline_2c3s(position = CZ_data$LCmeanDist[CZ_data$shore==x], sex = CZ_data$test_sex[CZ_data$shore==x],
                          cl = CZ_cline_params["cl", x], cr = CZ_cline_params["cr", x],
                          wl = exp(CZ_cline_params["lwl", x]), wr = exp(CZ_cline_params["lwr", x]),
                          crab = CZ_cline_params["crab", x], wave = CZ_cline_params["wave", x],
                          zs_c = CZ_cline_params["zs_c", x], zs_w = CZ_cline_params["zs_w", x],
                          sc = CZ_cline_params["sc", x], sh = CZ_cline_params["sh", x],
                          sw = CZ_cline_params["sw", x])
  phen_cline = as.data.frame(cbind(phen_cline, log_len=log(CZ_data[CZ_data$shore==x, ]$length_mm)))
  return(phen_cline)
})

CZs_phen_cline = lapply(CZs_phen_cline, function(x) {
  x[, "sex"] = ifelse(x[, "sex"]==1, "female", "male")
  return(x)
})

CZs_cline_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZs_phen_cline[[pl]], aes(position, log_len, col=sex)) +
    scale_color_manual(values = c("grey", "black")) +
    geom_ribbon(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x), fill="#e0ecf4", alpha=0.25) +
    geom_point() +
    geom_line(aes(position, z_x), size=1.3, alpha=0.7) +
    labs(title = islands[pl])
    # geom_errorbar(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x), width=0.25)
})


sim_mat = function(data, pos, isl, ce) {
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
    fml_dtr = rnorm(n = 200, mean = fml_m, sd = log(1.5))
  } else {
    fml_dtr = rnorm(n = 200, mean = fml_m, sd = fml_sd)
  }
  if (ml_sd > 0.4) {
    ml_dtr = rnorm(n = 200, mean = ml_m, sd = log(1.5))
  } else {
    ml_dtr = rnorm(n = 200, mean = ml_m, sd = ml_sd)
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
    write.table(YN, file = paste0("tables/gaus_skew/SKEW/sims/", isl, "_", ce, "_", round(pos), "_sim_YN.csv"), append = TRUE,
                sep = ",", row.names = FALSE, col.names = FALSE)
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}


CZs_mate_sim = function(s, centre, width) {
  map(seq_along(islands), function(x) {
    map(seq(from = 1, to = 39, by = 2), function(y) {
      # cline_pos = CZ_cline_params["cl", islands[x]] - y * CZ_cline_params["lwl", islands[x]]
      cline_pos = sum(CZ_cline_params[centre, islands[x]], s * y * CZ_cline_params[width, islands[x]])
      CZ_sim_sex = split(CZs_phen_cline[[x]], CZs_phen_cline[[x]][, "sex"])
      write.table(data.frame(male=as.numeric(), female=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
                  file = paste0("tables/gaus_skew/SKEW/sims/", islands[x], "_", centre, "_", round(cline_pos), "_sim_YN.csv"), sep = ",",
                  row.names = FALSE, col.names = TRUE)
      simYN = possibly(sim_mat, otherwise = "Missing snails")
      outYN = simYN(data = CZ_sim_sex, pos = cline_pos, isl = islands[x], ce = centre)
      return(outYN)
    })
  })
}
CZs_right_minus = CZs_mate_sim(s = -1, centre = "cr", width = "lwr")
CZs_right_plus = CZs_mate_sim(s = 1, centre = "cr", width = "lwr")
CZs_left_minus = CZs_mate_sim(s = -1, centre = "cl", width = "lwl")
CZs_left_plus = CZs_mate_sim(s = 1, centre = "cl", width = "lwl")

######################
# assortative mating #
######################
CZs_am = function(isls, mid_cline) {
  c_files = lapply(seq_along(isls), function(i) {
    lapply(seq_along(mid_cline), function(c) {
      list.files(path = "tables/gaus_skew/SKEW/sims", pattern = paste(isls[i], mid_cline[c], sep = "_"),
                 full.names = TRUE)
    })
  })
  c_df = lapply(seq_along(isls), function(i) {
    lapply(seq_along(mid_cline), function(c) {
      lapply(1:length(c_files[[i]][[c]]), function(f) {
        c_pos = strsplit(basename(c_files[[i]][[c]][f]), split = "_")[[1]][3]
        isl_df = read.csv(c_files[[i]][[c]][f])
        isl_df$male2 = isl_df$male^2
        p_cor = cor(isl_df[isl_df$mountYN==1, ]$female, isl_df[isl_df$mountYN==1, ]$male, method = "pearson")
        # ss_mod = glm(mountYN ~ male + male2, family = binomial(link = "logit"), data = isl_df)
        # p_ss_idx = which.max(fitted(ss_mod))
        pos_df = cbind(isl_df, position = as.integer(rep(c_pos, nrow(isl_df))), am_cor = rep(p_cor, nrow(isl_df)))
      })
    })
  })
  return(c_df)
}
CZs_am_df = CZs_am(isls = islands, mid_cline = c("cl", "cr"))



CZs_ss = function(cz_dat) {
  ss_mod = glm(mountYN ~ male + male2, family = binomial(link = "logit"), data = cz_dat)
  p_ss_idx = which.max(fitted(ss_mod))
  return(cbind(cz_dat, male_mx=rep(cz_dat$male[p_ss_idx], nrow(cz_dat)), male_av=rep(mean(cz_dat$male), nrow(cz_dat))))
}

CZs_am_ss = lapply(seq_along(islands), function(i) {
  lapply(seq_along(c("cl", "cr")), function(c) {
    lapply(seq_along(CZs_am_df[[i]][[c]]), function(d) {
      CZs_ss_poss = possibly(CZs_ss, otherwise = data.frame(male=as.numeric(), female=as.numeric(),
                                                            mountYN=as.integer(), male2=as.numeric(),
                                                            position=as.integer(), am_cor=as.numeric(),
                                                            male_mx=as.numeric(), male_av=as.numeric(),
                                                            stringsAsFactors=FALSE))
      CZs_ss_poss(cz_dat = CZs_am_df[[i]][[c]][[d]])
    })
  })
})


CZ_am_fun = function(data, isls, mid_cline) {
  mid_id = as.numeric(as.factor(mid_cline))
  CZ_clcr = rbind(rbindlist(data[[isls]][[mid_id[1]]]), rbindlist(data[[isls]][[mid_id[2]]]))
  return(CZ_clcr)
}

CZ_am_tot = lapply(seq_along(islands), function(i) {
  CZ_am_fun(data = CZs_am_ss, isls = i, mid_cline = c("cl", "cr"))
})

CZ_am_plot = function(data, x, y, yy, isls) {
  am_pl = ggplot(data = data, aes(x, y)) +
    geom_point(alpha=0.4) +
    geom_point(aes(x, yy), col="red", alpha=0.2) +
    geom_point(aes(x, am_cor), col="purple", size=2.5) +
    geom_point(aes(x, male_mx), col="green", size=2.5) +
    geom_point(aes(x, male_av), col="blue", size=2.5) +
    labs(x = "position", y = "ln size", title = isls)
  ggsave(filename = paste0("figures/gaus_skew/SKEW/ass_mat/", isls, "_sim_am.png"), plot = am_pl)
}

lapply(seq_along(islands), function(i) {
  CZ_am_plot(data = CZ_am_tot[[i]], x = CZ_am_tot[[i]]$position, y = CZ_am_tot[[i]]$male,
             yy = CZ_am_tot[[i]]$female, isls = islands[i])
})


####################
# sexual selection #
####################
