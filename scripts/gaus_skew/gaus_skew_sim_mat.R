rm(list = ls())

.packages = c("ggplot2", "dplyr", "rstan", "optparse", "tibble", "bayesplot", "data.table", "purrr",
              "pracma", "rgl", "parallel", "Rmisc")
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
              help = "directory for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$modelpars) | is.null(opt$clinepars)) {
  print_help(opt_parser)
  stop("At least three arguments must be supplied (input data, stan model and cline parameters).\n", call.=FALSE)
}

CZ_data = read.csv(opt$data, sep = ";")
skew_pars = read.csv(opt$modelpars, sep = ",")
skew_pars = column_to_rownames(skew_pars, var = "parameter")
CZ_cline_params = read.csv(opt$clinepars, row.names = 1)
pref_out = opt$output
# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
# skew_pars = read.csv("tables/gaus_skew/SKEW/gaus_skew_params.csv", sep = ";")
# skew_pars = column_to_rownames(skew_pars, var = "parameter")
# colnames(CZ_data)
# CZ_cline_params = read.csv("tables/clines/CZ_cline_params.csv", row.names = 1)
# pref_out = "gaus_skew/SKEW/sims/"
# pref_out = "gaus_skew/SKEW/sims2/"
# pref_out = "gaus_skew/SKEW/sims10000/"
dir.create(file.path("tables", pref_out))
dir.create(file.path("figures", pref_out))

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
# islands = "CZB"
# cline_2c3s(position = 3, sex = "female",
#            cl = CZ_cline_params["cl", islands], cr = CZ_cline_params["cr", islands],
#            wl = exp(CZ_cline_params["lwl", islands]), wr = exp(CZ_cline_params["lwr", islands]),
#            crab = CZ_cline_params["crab", islands], wave = CZ_cline_params["wave", islands],
#            zs_c = CZ_cline_params["zs_c", islands], zs_w = CZ_cline_params["zs_w", islands],
#            sc = CZ_cline_params["sc", islands], sh = CZ_cline_params["sh", islands],
#            sw = CZ_cline_params["sw", islands])
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
  x = mutate(x, figure="cline")
  return(x)
})

CZs_bin_cline = lapply(CZs_phen_cline, function(r) {
  brk = seq(from = as.integer(min(r[,"position"])-1), to = as.integer(max(r[,"position"])+1), by = 10)
  bin = cut(r[,"position"], brk)
  bin_dt = cbind(r, bin)
  log_len_m = aggregate(bin_dt[, c("log_len", "position")], list(bin_dt$bin, bin_dt$sex), CI)
  log_len_m = mutate(log_len_m, figure = "cline")
  # log_len_uCI = aggregate(bin_dt[, "log_len"], list(bin_dt$bin, bin_dt$sex), CI)
  # log_len_lCI = aggregate(bin_dt[, "log_len"], list(bin_dt$bin, bin_dt$sex), CI)['lower']
  return(log_len_m)
  })

CZs_cline_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZs_phen_cline[[pl]]) +
    geom_vline(xintercept = CZ_cline_params["cl", pl], linetype = "dashed") +
    geom_vline(xintercept = CZ_cline_params["cr", pl], linetype = "dashed") +
    scale_color_manual(values = c("lightcoral", "black")) +
    scale_fill_manual(values = c("lightcoral", "black")) +
    geom_ribbon(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x, fill=sex), alpha=0.15) +
    geom_point(data = CZs_bin_cline[[pl]], aes(x = CZs_bin_cline[[pl]]$position[, 'mean'],
                                               y = CZs_bin_cline[[pl]]$log_len[, 'mean'],
                                               col=CZs_bin_cline[[pl]]$Group.2)) +
    geom_line(aes(position, z_x, col=sex), size=1.3, alpha=0.7) +
    labs(x = paste0(islands[pl], " shore position"), y = 'ln shell size', fill='', col='') +
    theme(legend.position = 'top',
          legend.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(face = "bold", size = 15))
    # geom_errorbar(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x), width=0.25)
})
# lapply(seq_along(islands), function(s) {
#   ggsave(filename = paste0("figures/clines/", islands[s], "_size_sex.png"), plot = CZs_cline_plot[[s]])
# })

isl_pos = sapply(islands, function(x) {
  isl_c = round(c(CZ_cline_params["cl", x], CZ_cline_params["cr", x]))
  isl_rng = range(CZ_data[CZ_data$shore==x, ]$LCmeanDist)
  # str(isl_rng)
  wavel = sort(round(seq(from = isl_c[1], to = 0, by = -10)))
  crab = sort(round(seq(from = isl_c[2], to = isl_c[1]+5, by = -10)))
  waver = round(seq(from = isl_c[2]+10, to = isl_rng[2], by = 10))
  return(c(wavel, crab, waver))
})
# lapply(isl_pos, length)
# pos = 60
sim_mat = function(pos, isl) {
  bar = list()
  YN = data.frame()

  # smp_pos = seq(round(pos)-2,round(pos)+2)
  # smp_pos_m = smp_pos + dnorm(x = 1, mean = 0, sd = 1.5)
  # fml_df = rbindlist(lapply(seq_along(smp_pos), function(z) {
  #   data$female[round(data$female$position)==smp_pos[z], ]
  # }))
  # ml_df = rbindlist(lapply(seq_along(smp_pos), function(z) {
  #   data$male[round(data$male$position)==smp_pos[z], ]
  # }))

  # fml_m = mean(CZ_sim_sex$female[round(CZ_sim_sex$female$position)==seq(round(pos)-1,round(pos)+1), ]$z_x)
  # fml_sd = mean(CZ_sim_sex$female[round(CZ_sim_sex$female$position)==c(round(pos)-1,round(pos),round(pos)+1), ]$s_x)
  # fml_m = mean(fml_df$z_x)
  fml_m = as.numeric(cline_2c3s(position = pos, sex = "female",
                                cl = CZ_cline_params["cl", isl], cr = CZ_cline_params["cr", isl],
                                wl = exp(CZ_cline_params["lwl", isl]), wr = exp(CZ_cline_params["lwr", isl]),
                                crab = CZ_cline_params["crab", isl], wave = CZ_cline_params["wave", isl],
                                zs_c = CZ_cline_params["zs_c", isl], zs_w = CZ_cline_params["zs_w", isl],
                                sc = CZ_cline_params["sc", isl], sh = CZ_cline_params["sh", isl],
                                sw = CZ_cline_params["sw", isl])[,"z_x"])
  # fml_sd = mean(fml_df$s_x)
  fml_sd = as.numeric(cline_2c3s(position = pos, sex = "female",
                                 cl = CZ_cline_params["cl", isl], cr = CZ_cline_params["cr", isl],
                                 wl = exp(CZ_cline_params["lwl", isl]), wr = exp(CZ_cline_params["lwr", isl]),
                                 crab = CZ_cline_params["crab", isl], wave = CZ_cline_params["wave", isl],
                                 zs_c = CZ_cline_params["zs_c", isl], zs_w = CZ_cline_params["zs_w", isl],
                                 sc = CZ_cline_params["sc", isl], sh = CZ_cline_params["sh", isl],
                                 sw = CZ_cline_params["sw", isl])[,"s_x"])
  # ml_m = mean(ml_df$z_x)
  # ml_m = fml_m + rnorm(n=1,mean=0,sd=1.5)
  # ml_sd = mean(ml_df$s_x)

  # if (fml_sd > 0.4) {
  #   fml_dtr = rnorm(n = 1000, mean = fml_m, sd = log(1.5))
  # } else {
  fml_dtr = rnorm(n = 10000, mean = fml_m, sd = fml_sd)
  # }
  # if (ml_sd > 0.4) {
  #   ml_dtr = rnorm(n = 1000, mean = ml_m, sd = log(1.5))
  # } else {
  #   ml_dtr = rnorm(n = 1000, mean = ml_m, sd = ml_sd)
  # }
  
  for (f in seq_along(fml_dtr)) {
    success=FALSE
    i=1
    fem = fml_dtr[f]
    while (!success) {
      # m = sample(ml_dtr, 1, replace = FALSE)
      mpos = pos + rnorm(n=1, mean=0, sd=1.5)
      msize = as.numeric(cline_2c3s(position = mpos, sex = "male",
                                    cl = CZ_cline_params["cl", isl], cr = CZ_cline_params["cr", isl],
                                    wl = exp(CZ_cline_params["lwl", isl]), wr = exp(CZ_cline_params["lwr", isl]),
                                    crab = CZ_cline_params["crab", isl], wave = CZ_cline_params["wave", isl],
                                    zs_c = CZ_cline_params["zs_c", isl], zs_w = CZ_cline_params["zs_w", isl],
                                    sc = CZ_cline_params["sc", isl], sh = CZ_cline_params["sh", isl],
                                    sw = CZ_cline_params["sw", isl])[,"z_x"])
      msd = as.numeric(cline_2c3s(position = mpos, sex = "male",
                                  cl = CZ_cline_params["cl", isl], cr = CZ_cline_params["cr", isl],
                                  wl = exp(CZ_cline_params["lwl", isl]), wr = exp(CZ_cline_params["lwr", isl]),
                                  crab = CZ_cline_params["crab", isl], wave = CZ_cline_params["wave", isl],
                                  zs_c = CZ_cline_params["zs_c", isl], zs_w = CZ_cline_params["zs_w", isl],
                                  sc = CZ_cline_params["sc", isl], sh = CZ_cline_params["sh", isl],
                                  sw = CZ_cline_params["sw", isl])[,"s_x"])
      m = rnorm(n = 1, mean = msize, sd = msd)
      p = 0.01 + skew_pars["b","mean"] * exp(-0.5 * (((fem - m) - skew_pars["c","mean"])
                                                     / skew_pars["d","mean"])^2) *
        (1 + erf(skew_pars["alpha","mean"] * ((fem - m) - skew_pars["c","mean"]) / (1.414214 * skew_pars["d","mean"])))
      s = rbinom(n = 1, size = 1, prob = p)
      YN[i,'male'] = m
      YN[i,'female'] = fem
      YN[i,'sk_prob'] = p
      YN[i,'mountYN'] = s
      success = (s > 0)
      i = i + 1
      if (s > 0) {
        cat("male size", m, "mated female size", fem, ".\n")
      }
    }
    write.table(YN, file = paste0("tables/", pref_out, isl, "_", round(pos), "_sim_YN.csv"), append = TRUE,
                sep = ",", row.names = FALSE, col.names = FALSE)
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}

# CZs_mate_sim = function(centre, width) {
#   map(seq_along(islands), function(x) {
#     map(seq_along(c(5, 80)), function(y) {
#     # map(seq(from = 1, to = 20, by = 2), function(y) {
#       # cline_pos = CZ_cline_params["cl", islands[x]] - y * CZ_cline_params["lwl", islands[x]]
#       # cline_pos = sum(CZ_cline_params[centre, islands[x]], s * y * CZ_cline_params[width, islands[x]])
#       cline_pos = y
#       # CZ_sim_sex = split(CZs_phen_cline[[x]], CZs_phen_cline[[x]][, "sex"])
#       # write.table(data.frame(male=as.numeric(), female=as.numeric(), sk_prob=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
#       #             file = paste0("tables/", pref_out, islands[x], "_", centre, "_", round(cline_pos), "_sim_YN.csv"), sep = ",",
#       #             row.names = FALSE, col.names = TRUE)
#       write.table(data.frame(male=as.numeric(), female=as.numeric(), sk_prob=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
#                   file = paste0("tables/", pref_out, islands[x], "_", round(cline_pos), "_sim_YN.csv"), sep = ",",
#                   row.names = FALSE, col.names = TRUE)
#       simYN = possibly(sim_mat, otherwise = "Missing snails")
#       outYN = simYN(data = CZ_sim_sex, pos = cline_pos, isl = islands[x])
#       return(outYN)
#     })
#   })
# }
# CZs_right_minus = CZs_mate_sim(s = -1, centre = "cr", width = "lwr")
# CZs_right_plus = CZs_mate_sim(s = 1, centre = "cr", width = "lwr")
# CZs_left_minus = CZs_mate_sim(s = -1, centre = "cl", width = "lwl")
# CZs_left_plus = CZs_mate_sim(s = 1, centre = "cl", width = "lwl")

map(seq_along(islands), function(x) {
  map(isl_pos[[x]], function(y) {
    cline_pos = y
    write.table(data.frame(male=as.numeric(), female=as.numeric(), sk_prob=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
                file = paste0("tables/", pref_out, islands[x], "_", round(cline_pos), "_sim_YN.csv"), sep = ",",
                row.names = FALSE, col.names = TRUE)
    simYN = possibly(sim_mat, otherwise = "Missing snails")
    outYN = simYN(pos = cline_pos, isl = islands[x])
    # return(outYN)
    return(cat("island", islands[x], "at position", cline_pos, "has been simulated ...\n"))
  })
})

######################
# assortative mating #
######################
# rm(list = ls())
# option_list = list(
#   make_option(c("-d", "--data"), type="character", default=NULL,
#               help="input data", metavar="character"),
#   make_option(c("-o", "--output"), type = "character", default = "output_sim",
#               help = "directory for output files [default: %default]", metavar = "character"))
# opt_parser = OptionParser(option_list=option_list)
# opt = parse_args(opt_parser)
# CZ_data = read.csv(opt$data, sep = ";")
# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
# islands = as.character(unique(CZ_data$shore))
# pref_out = opt$output
# pref_out = "gaus_skew/SKEW/sims/"

CZs_am = function(isls) {
  c_files = lapply(seq_along(isls), function(i) {
    list.files(path = paste0("tables/", pref_out), pattern = isls[i],
               full.names = TRUE)
  })
  c_df = lapply(seq_along(isls), function(i) {
    lapply(1:length(c_files[[i]]), function(f) {
      c_pos = strsplit(basename(c_files[[i]][f]), split = "_")[[1]][2]
      isl_df = read.csv(c_files[[i]][f])
      isl_df$male2 = isl_df$male^2
      p_cor = cor(isl_df[isl_df$mountYN==1, ]$female, isl_df[isl_df$mountYN==1, ]$male, method = "pearson")
      # ss_mod = glm(mountYN ~ male + male2, family = binomial(link = "logit"), data = isl_df)
      # p_ss_idx = which.max(fitted(ss_mod))
      pos_df = cbind(isl_df, position = as.integer(rep(c_pos, nrow(isl_df))), am_cor = rep(p_cor, nrow(isl_df)))
    })
  })
  return(c_df)
}
CZs_am_df = CZs_am(isls = islands)
# CZs_am_df_one = read.csv("tables/gaus_skew/SKEW/sims/CZA_9_sim_YN.csv")
# CZs_am_df_two = read.csv("tables/gaus_skew/SKEW/sims/CZA_19_sim_YN.csv")
# head(CZs_am_df_one)
# with(data = CZs_am_df_one[CZs_am_df_one$mountYN==1, ], cor(female, male))
# with(data = CZs_am_df_two[CZs_am_df_two$mountYN==1, ], cor(female, male))
# hist(CZs_am_df_one$female)
# hist(CZs_am_df_two$female)
# coef(glm(mountYN ~ male + male2, family = binomial(link = "logit"), data = CZs_am_df[[1]][[3]]))["male2"]

CZs_ss = function(cz_dat) {
  ss_mod = glm(mountYN ~ male + male2, family = binomial(link = "logit"), data = cz_dat)
  dir_ss = coef(ss_mod)["male"]
  sta_ss = coef(ss_mod)["male2"]
  p_ss_idx = which.max(fitted(ss_mod))
  return(cbind(cz_dat, male_mx=rep(cz_dat$male[p_ss_idx], nrow(cz_dat)), male_av=rep(mean(cz_dat$male), nrow(cz_dat)),
               dss=rep(dir_ss, nrow(cz_dat)), sss=rep(sta_ss, nrow(cz_dat)), pss=fitted(ss_mod)))
}

CZs_am_ss = lapply(seq_along(islands), function(i) {
  lapply(seq_along(CZs_am_df[[i]]), function(d) {
    CZs_ss_poss = possibly(CZs_ss, otherwise = data.frame(male=as.numeric(), female=as.numeric(), sk_prob=as.numeric(),
                                                          mountYN=as.integer(), male2=as.numeric(),
                                                          position=as.integer(), am_cor=as.numeric(),
                                                          male_mx=as.numeric(), male_av=as.numeric(),
                                                          dss=as.numeric(), sss=as.numeric(), pss=as.numeric(),
                                                          stringsAsFactors=FALSE))
    CZs_ss_poss(cz_dat = CZs_am_df[[i]][[d]])
  })
})

CZ_am_fun = function(data, isls) {
  CZ_clcr = rbindlist(data[[isls]])
  return(CZ_clcr)
}

CZ_am_tot = lapply(seq_along(islands), function(i) {
  CZ_am_fun(data = CZs_am_ss, isls = i)
})

CZ_am_ss_fig = lapply(seq_along(islands), function(f) {
  am4fig = CZ_am_tot[[f]][, c("position", "am_cor")]
  am4fig = mutate(am4fig, figure="am")
  colnames(am4fig)[2] = "par_value"
  dss4fig = CZ_am_tot[[f]][, c("position", "dss")]
  dss4fig = mutate(dss4fig, figure="dss")
  colnames(dss4fig)[2] = "par_value"
  sss4fig = CZ_am_tot[[f]][, c("position", "sss")]
  sss4fig = mutate(sss4fig, figure="sss")
  colnames(sss4fig)[2] = "par_value"
  df4fig = rbind(am4fig, dss4fig, sss4fig)
  return(distinct(df4fig))
})
# sim_set = "sims2_amss"
# sim_set = "sims1_amss"
# lapply(seq_along(islands), function(t) {
#   write.csv(x = CZ_am_ss_fig[[t]],
#             file = paste0("tables/", dirname(pref_out), "/", sim_set, "/", islands[t], "_am_ss_pars.csv"),
#             row.names = FALSE)
# })
# lapply(seq_along(islands), function(pl) {
#   ggplot(data = CZ_am_ss_fig[[pl]]) +
#     facet_grid(rows = vars(figure), scales = "free") +
#     geom_point(aes(x = position, y = par_value))
# })
CZs_cline_am_ss_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZs_phen_cline[[pl]]) +
    facet_grid(rows = vars(figure), scales = "free") +
    geom_vline(xintercept = CZ_cline_params["cl", pl], linetype = "dashed") +
    geom_vline(xintercept = CZ_cline_params["cr", pl], linetype = "dashed") +
    scale_color_manual(values = c("red", "blue")) +
    scale_fill_manual(values = c("red", "blue")) +
    geom_ribbon(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x, fill=sex), alpha=0.15) +
    geom_point(data = CZs_bin_cline[[pl]], aes(x = CZs_bin_cline[[pl]]$position[, 'mean'],
                                               y = CZs_bin_cline[[pl]]$log_len[, 'mean'],
                                               col=CZs_bin_cline[[pl]]$Group.2)) +
    geom_line(aes(position, z_x, col=sex), size=1.3, alpha=0.7) +
    geom_point(data = CZ_am_ss_fig[[pl]], aes(x = position, y = par_value)) +
    labs(x = paste0(islands[pl], " shore position"), y = 'ln shell size', fill='', col='') +
    theme(legend.position = 'top',
          legend.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(face = "bold", size = 15))
  # geom_errorbar(aes(x=position, ymin=z_x-s_x, ymax=z_x+s_x), width=0.25)
})
lapply(seq_along(islands), function(s) {
  ggsave(filename = paste0("figures/gaus_skew/SKEW/sims/", islands[s], "_cline_am_ss.png"),
         plot = CZs_cline_am_ss_plot[[s]])
})

# CZ_am_ss_2_fig = read.csv("tables/gaus_skew/SKEW/sims2_amss/CZA_am_ss_pars.csv")
# CZ_am_ss_1_fig = CZ_am_ss_fig[[1]]
# pdf("misc/CZA_2sims_comp.pdf", width = 7, height = 7)
# ggplot(data = CZ_am_ss_1_fig) +
#   facet_grid(rows = vars(figure), scales = "free") +
#   geom_point(aes(x = position, y = par_value), col="black") +
#   geom_point(data = CZ_am_ss_2_fig, aes(x = position, y = par_value), col="orange") +
#   labs(x="CZA shore position")
# dev.off()
# 
# CZ_am_plot = function(data, x, y, yy, isls) {
#   am_pl = ggplot(data = data, aes(x, y)) +
#     geom_point(alpha=0.4) +
#     geom_point(aes(x, yy), col="red", alpha=0.2) +
#     geom_point(aes(x, am_cor), col="purple", size=2.5) +
#     # geom_point(aes(x, male_mx), col="green", size=2.5) +
#     # geom_point(aes(x, male_av), col="blue", size=2.5) +
#     labs(x = "position", y = "ln size", title = isls)
#   ggsave(filename = paste0("figures/gaus_skew/SKEW/ass_mat/", isls, "_sim_am.png"), plot = am_pl)
# }
# 
# lapply(seq_along(islands), function(i) {
#   CZ_am_plot(data = CZ_am_tot[[i]], x = CZ_am_tot[[i]]$position, y = CZ_am_tot[[i]]$male,
#              yy = CZ_am_tot[[i]]$female, isls = islands[i])
# })
# 
# 
# ############
# # 3D plots #
# ############
# CZ_sim_grid = mclapply(seq_along(islands), function(s) {
#   CZ_am_Y = CZ_am_tot[[s]]
#   # CZ_am_Y = split(CZ_am_tot[[s]], CZ_am_tot[[s]]$mountYN)$`1`
#   # CZ_am_Y$x = CZ_am_Y$position
#   CZ_am_Y$sratio = round(CZ_am_Y$female-CZ_am_Y$male, 1)
#   CZ_am_Y$male = round(CZ_am_Y$male, 1)
#   # plot(CZ_am_Y$position, CZ_am_Y$sratio)
#   # df_grid = expand.grid(position = seq(from = min(CZ_am_Y$position), to = max(CZ_am_Y$position), by = 10),
#   #                       sratio = unique(round(CZ_am_Y$sratio, 1)))
#   # plot(df_grid$position, df_grid$sratio)
#   # df_merge = merge(CZ_am_Y, df_grid, by = c("position", "sratio"), all = TRUE)
#   # plot(CZ_am_Y$position, CZ_am_Y$sratio)
#   cat("Adjusting", islands[s], "cline position onto a grid ...\n")
#   for (i in 1:nrow(CZ_am_Y)) {
#     i_rest = as.integer(CZ_am_Y$position[i]) %% 10
#     CZ_am_Y$xx[i] = CZ_am_Y$position[i] - i_rest
#   }
#   # plot(CZ_am_Y$xx, CZ_am_Y$male)
#   return(CZ_am_Y)
# }, mc.cores = 4)
# lapply(CZ_sim_grid, head)
# # plot(CZ_sim_grid[[1]]$position, CZ_sim_grid[[1]]$sratio)
# 
# # CZ_am_Y = split(CZ_am_tot[[1]], CZ_am_tot[[1]]$mountYN)$`1`
# # CZ_am_Y = CZ_am_tot[[1]]
# # CZ_am_Y = sample_n(tbl = CZ_am_Y, size = 1000, replace = FALSE)
# # CZ_am_Y$x = CZ_am_Y$position
# # CZ_am_Y$y = round(CZ_am_Y$female-CZ_am_Y$male, 1)
# # CZ_am_Y$y = round(CZ_am_Y$male, 1)
# # str(CZ_am_Y)
# # summary(CZ_am_Y)
# # plot(CZ_am_Y$x, CZ_am_Y$y)
# # for (i in 1:nrow(CZ_am_Y)) {
# #   i_rest = as.integer(CZ_am_Y$position[i]) %% 10
# #   CZ_am_Y$xx[i] = CZ_am_Y$position[i] - i_rest
# #   print(df_merge$xx[i])
# # }
# # plot(CZ_am_Y$xx, CZ_am_Y$y)
# #
# # df_grid  <- expand.grid(x = seq(from = min(CZ_am_Y$x), to = max(CZ_am_Y$x), by = 10), y = unique(round(CZ_am_Y$y, 1)))
# # head(df_grid)
# # plot(df_grid$x, df_grid$y)
# # plot(pp20$x, pp20$y)
# # head(CZ_am_Y)
# 
# # df_merge <- na.omit(merge(CZ_am_Y, df_grid, by = c("x", "y"), all = TRUE))
# # plot(df_merge$x, df_merge$y)
# # nrow(na.omit(df_merge))
# # for (i in 1:nrow(df_merge)) {
# #   i_rest = as.integer(df_merge$x[i]) %% 10
# #   df_merge$xx[i] = df_merge$x[i] - i_rest
# #   print(df_merge$xx[i])
# # }
# 
# # head(na.omit(df_merge))
# # plot(df_merge$xx, df_merge$y)
# # CZ_sim_grid = CZ_am_tot
# 
# # CZ_sim_grid = lapply(seq_along(islands), function(sp) {
# #   CZ_sim_grid[[sp]]$male = round(CZ_sim_grid[[sp]]$male, 1)
# #   return(CZ_sim_grid[[sp]])
# # })
# 
# df_am = lapply(seq_along(islands), function(sp) {
#   split(CZ_sim_grid[[sp]], CZ_sim_grid[[sp]]$mountYN)$`1`
# }) %>%
#   lapply(., function(g) {
#     group_by(g, xx, sratio)
#   }) %>%
#   lapply(., function(m) {
#     summarise_all(m, mean)
#   }) %>%
#   lapply(., function(p) {
#     ggplot(p, aes(xx, sratio)) +
#       geom_tile(aes(fill = am_cor)) +
#       scale_fill_viridis_c()
#   })
# df_am_c = lapply(seq_along(islands), function(cc) {
#   df_am[[cc]] + geom_vline(xintercept = as.numeric(CZ_cline_params["cl", islands[cc]]), linetype="dotted", color = "black", size=1.5) +
#   geom_vline(xintercept = as.numeric(CZ_cline_params["cr", islands[cc]]), linetype="dotted", color = "black", size=1.5)
# })
# 
# lapply(seq_along(islands), function(sv) {
#   ggsave(filename = paste0("figures/", pref_out, islands[sv], "_sim_am_grid.png"), plot = df_am_c[[sv]])
# })
# 
# 
# df_ss = lapply(CZ_sim_grid, function(g) {
#   group_by(g, xx, male)
# }) %>%
#   lapply(., function(m) {
#     summarise_all(m, mean)
#   })
# df_av = lapply(df_ss, function(g) {
#   group_by(g, xx)
# }) %>%
#   lapply(., function(m) {
#     summarise_all(m, mean)
#   })
# df_mx = lapply(df_ss, function(g) {
#   group_by(g, xx)
# }) %>%
#   lapply(., function(m) {
#     top_n(m, 1, sk_prob)
#   })
# df_ss_pl = lapply(seq_along(islands), function(p) {
#   ggplot(data = df_ss[[p]], aes(xx, male)) +
#     geom_tile(aes(fill = sk_prob)) +
#     scale_fill_viridis_c() +
#     geom_point(data = df_mx[[p]], aes(xx, male), col="red", size = 2.5) +
#     geom_point(data = df_av[[p]], aes(xx, male), col="white", size = 1.8) +
#     geom_vline(xintercept = as.numeric(CZ_cline_params["cl", islands[p]]), linetype="dotted", color = "black", size=1.5) +
#     geom_vline(xintercept = as.numeric(CZ_cline_params["cr", islands[p]]), linetype="dotted", color = "black", size=1.5)
# })
# lapply(seq_along(islands), function(sv) {
#   ggsave(filename = paste0("figures/", pref_out, islands[sv], "_sim_ss_grid.png"), plot = df_ss_pl[[sv]])
# })
# 
# 
# # Missing values
# # ggplot(df_merge, aes(x = xx, y = y)) +
# #   geom_tile(data = subset(df_merge, !is.na(am_cor)), aes(fill = am_cor)) +
# #   geom_tile(data = subset(df_merge,  is.na(am_cor)), aes(colour = "NA"),
# #             linetype = 0, fill = "black", alpha = 0.5)
# 
# 
# # ggplot(df_merge, aes(xx, y)) +
# #   geom_raster(aes(fill = am_cor), interpolate = TRUE) +
# #   scale_fill_viridis_c()
# #   scale_fill_continuous(low="thistle2", high="darkred",
# #                         guide="colorbar",na.value="thistle2")
