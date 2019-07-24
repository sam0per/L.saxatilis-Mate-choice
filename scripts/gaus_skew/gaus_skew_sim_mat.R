rm(list = ls())

# devtools::install_github("thomasp85/patchwork")
.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "dplyr", "rstan", "optparse", "tibble", "bayesplot", "data.table", "purrr",
              "pracma", "rgl", "parallel", "Rmisc", "bbmle")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="input data", metavar="character"),
  make_option(c("-m", "--modelpars"), type="character", default=NULL,
              help="mean estimates of the inferred parameters", metavar="character"),
  make_option(c("-n", "--numrun"), type = "integer", default = 3,
              help = "number of replicates for the full simulations [default: %default]", metavar = "integer"),
  make_option(c("-o", "--output"), type = "character", default = "output_sim",
              help = "directory for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$modelpars)) {
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input data and stan model).\n", call.=FALSE)
}

CZ_data = read.csv(opt$data, sep = ";")
skew_pars = read.csv(opt$modelpars, sep = ";")
skew_pars = column_to_rownames(skew_pars, var = "parameter")
# CZ_cline_params = read.csv(opt$clinepars, row.names = 1)
pref_out = opt$output
# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
# skew_pars = read.csv("tables/gaus_skew/SKEW/gaus_skew_params.csv", sep = ";")
# skew_pars = column_to_rownames(skew_pars, var = "parameter")
# colnames(CZ_data)
# CZ_cline_params = read.csv("tables/clines/CZ_cline_params.csv", row.names = 1)
# pref_out = "gaus_skew/SKEW/sims/"
# pref_out = "gaus_skew/SKEW/sims2/"
# pref_out = "gaus_skew/SKEW/sims10000/"
# pref_out = "gaus_skew/SKEW/sims_test/"
dir.create(file.path("tables", pref_out))
dir.create(file.path("figures", pref_out))
dir.create(file.path("tables", pref_out, "stats"))
dir.create(file.path("figures", pref_out, "stats"))
# islands = "CZB"
islands = as.character(unique(CZ_data$shore))

cline_2c4s <- function(phen,position,sex,cl,cr,lwl,lwr,crab,wave,zs_c,zs_w,sc,shl,sh,sw){
  wl = exp(lwl)
  wr = exp(lwr)
  # sc = exp(lsc)
  # sh = exp(lsh)
  # sw = exp(lsw)
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]
  s_xl <- sqrt(sc^2 + 4*p_xl*(1-p_xl)*shl^2 + (p_xl^2)*(sw^2-sc^2))
  
  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  # combined cline
  cond <- z_x < z_xl
  z_x[cond] <- z_xl[cond]
  s_x[cond] <- s_xl[cond]
  # z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  # s_x[z_x < z_xl] <- s_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  # phen_cline = data.frame(phen_cline = z_x, sd_cline = s_x, sex = sex, position = position)
  # return(phen_cline)
  return(minusll)
}

theta.init = list(CZA=list(cl=130, cr=280, lwl=3, lwr=2.3, crab=-2.1, wave=-1.9, zs_c=-0.1, zs_w=-0.1,
                           sc=0.2, shl=0.2, sh=0.2, sw=0.2),
                  CZB=list(cl=70, cr=150, lwl=1.6, lwr=3.9, crab=-2.5, wave=-1.5, zs_c=-0.1, zs_w=-0.1,
                           sc=0.2, shl=0.2, sh=0.2, sw=0.2),
                  CZC=list(cl=50, cr=125, lwl=1.5, lwr=3, crab=-2.5, wave=-1.5, zs_c=-0.1, zs_w=-0.1,
                           sc=0.2, shl=0.2, sh=0.2, sw=0.2),
                  CZD=list(cl=80, cr=175, lwl=1.6, lwr=1.6, crab=-2.5, wave=-1.5, zs_c=-0.1, zs_w=-0.1,
                           sc=0.2, shl=0.2, sh=0.2, sw=0.2))
# theta.init[[1]]
cline_pars = lapply(seq_along(islands), function(c) {
  cat("Fitting cline for island", islands[c], "...\n")
  mle.cline.2c4s = mle2(cline_2c4s, theta.init[[c]],
                        control=list(parscale=abs(unlist(theta.init[[c]]))),
                        data=list(phen=-log(CZ_data[CZ_data$shore==islands[c],]$length_mm),
                                  position=CZ_data[CZ_data$shore==islands[c],]$LCmeanDist,
                                  sex=CZ_data[CZ_data$shore==islands[c],]$test_sex))
  cline_est = round(coef(summary(mle.cline.2c4s)), 3)
  return(cline_est)
})

cline_sims <- function(phen,position,sex,cl,cr,lwl,lwr,crab,wave,zs_c,zs_w,sc,shl,sh,sw){
  wl = exp(lwl)
  wr = exp(lwr)
  # sc = exp(lsc)
  # sh = exp(lsh)
  # sw = exp(lsw)
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]
  s_xl <- sqrt(sc^2 + 4*p_xl*(1-p_xl)*shl^2 + (p_xl^2)*(sw^2-sc^2))
  
  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  # combined cline
  cond <- z_x < z_xl
  z_x[cond] <- z_xl[cond]
  s_x[cond] <- s_xl[cond]
  # z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  # s_x[z_x < z_xl] <- s_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  phen_cline = data.frame(phen_cline = z_x, sd_cline = s_x, sex = sex, position = position)
  return(phen_cline)
  # return(minusll)
}

CZs_phen_cline = lapply(seq_along(islands), function(x) {
  cat("Extracting fitted cline values for", islands[x], "...\n")
  cline_df = cline_sims(phen = -log(CZ_data[CZ_data$shore==islands[x],]$length_mm),
                        position = CZ_data[CZ_data$shore==islands[x],]$LCmeanDist,
                        sex = CZ_data[CZ_data$shore==islands[x],]$test_sex,
                        cl = cline_pars[[x]]['cl', 'Estimate'], cr = cline_pars[[x]]['cr', 'Estimate'],
                        lwl = cline_pars[[x]]['lwl', 'Estimate'], lwr = cline_pars[[x]]['lwr', 'Estimate'],
                        crab = cline_pars[[x]]['crab', 'Estimate'], wave = cline_pars[[x]]['wave', 'Estimate'],
                        zs_c = cline_pars[[x]]['zs_c', 'Estimate'], zs_w = cline_pars[[x]]['zs_w', 'Estimate'],
                        sc = cline_pars[[x]]['sc', 'Estimate'], shl = cline_pars[[x]]['shl', 'Estimate'],
                        sh = cline_pars[[x]]['sh', 'Estimate'], sw = cline_pars[[x]]['sw', 'Estimate'])
  clinefit_obs = as.data.frame(cbind(cline_df, log_len=log(CZ_data[CZ_data$shore==islands[x], ]$length_mm)))
  return(clinefit_obs)
})

CZs_phen_cline = lapply(CZs_phen_cline, function(x) {
  # x[, "sex"] = ifelse(x[, "sex"]==1, "female", "male")
  x = mutate(x, figure="cline")
  return(x)
})

CZs_bin_cline = lapply(CZs_phen_cline, function(r) {
  brk = seq(from = as.integer(min(r[,"position"])-1), to = as.integer(max(r[,"position"])+1), by = 10)
  bin = cut(r[,"position"], brk)
  bin_dt = cbind(r, bin)
  log_len_m = aggregate(bin_dt[, c("log_len", "position")], list(bin_dt$bin, bin_dt$sex), CI)
  log_len_m = mutate(log_len_m, figure = "cline")
  return(log_len_m)
  })

# CZs_cline_plot = lapply(seq_along(islands), function(pl) {
#   ggplot(data = CZs_phen_cline[[pl]]) +
#     geom_vline(xintercept = cline_pars[[pl]]['cl', 'Estimate'], linetype = "dashed") +
#     geom_vline(xintercept = cline_pars[[pl]]['cr', 'Estimate'], linetype = "dashed") +
#     scale_color_manual(values = c("red", "blue")) +
#     scale_fill_manual(values = c("red", "blue")) +
#     geom_ribbon(aes(x=position, ymin=abs(phen_cline)-sd_cline, ymax=abs(phen_cline)+sd_cline, fill=sex), alpha=0.15) +
#     geom_point(data = CZs_bin_cline[[pl]], aes(x = position[, 'mean'],
#                                                y = log_len[, 'mean'],
#                                                col = Group.2)) +
#     geom_line(aes(position, abs(phen_cline), col=sex), size=1.2, alpha=0.7) +
#     labs(x = paste0(islands[pl], " shore position"), y = 'ln shell size', fill='', col='') +
#     theme(legend.position = 'top',
#           legend.text = element_text(size = 13),
#           axis.title = element_text(face = "bold", size = 15))
# })
# lapply(seq_along(islands), function(s) {
#   ggsave(filename = paste0("figures/clines/", islands[s], "_size_sex.png"), plot = CZs_cline_plot[[s]])
# })

isl_pos = sapply(seq_along(islands), function(x) {
  cat("Chopping", islands[x], "transect into equally-distanced shore positions from\nthe left centre",
      cline_pars[[x]]['cl', 'Estimate'], "and the right centre", cline_pars[[x]]['cr', 'Estimate'], "...\n")
  isl_c = round(c(cline_pars[[x]]['cl', 'Estimate'], cline_pars[[x]]['cr', 'Estimate']))
  isl_rng = range(CZ_data[CZ_data$shore==islands[x], ]$LCmeanDist)
  # str(isl_rng)
  wavel = sort(round(seq(from = isl_c[1], to = 0, by = -10)))
  crab = sort(round(seq(from = isl_c[2], to = isl_c[1]+5, by = -10)))
  waver = round(seq(from = isl_c[2]+10, to = isl_rng[2], by = 10))
  return(c(wavel, crab, waver))
})
# lapply(isl_pos, length)
# pos = 9
# isl = "CZA"
# run = 2
sim_mat = function(pos, isl, run) {
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
  fml_cline = cline_sims(phen = -log(CZ_data[CZ_data$shore==islands[isl],]$length_mm),
                         position = pos, sex = 'female',
                         cl = cline_pars[[isl]]['cl', 'Estimate'], cr = cline_pars[[isl]]['cr', 'Estimate'],
                         lwl = cline_pars[[isl]]['lwl', 'Estimate'], lwr = cline_pars[[isl]]['lwr', 'Estimate'],
                         crab = cline_pars[[isl]]['crab', 'Estimate'], wave = cline_pars[[isl]]['wave', 'Estimate'],
                         zs_c = cline_pars[[isl]]['zs_c', 'Estimate'], zs_w = cline_pars[[isl]]['zs_w', 'Estimate'],
                         sc = cline_pars[[isl]]['sc', 'Estimate'], shl = cline_pars[[isl]]['shl', 'Estimate'],
                         sh = cline_pars[[isl]]['sh', 'Estimate'], sw = cline_pars[[isl]]['sw', 'Estimate'])
  # fml_m = as.numeric(cline_2c3s(position = pos, sex = "female",
  #                               cl = CZ_cline_params["cl", isl], cr = CZ_cline_params["cr", isl],
  #                               wl = exp(CZ_cline_params["lwl", isl]), wr = exp(CZ_cline_params["lwr", isl]),
  #                               crab = CZ_cline_params["crab", isl], wave = CZ_cline_params["wave", isl],
  #                               zs_c = CZ_cline_params["zs_c", isl], zs_w = CZ_cline_params["zs_w", isl],
  #                               sc = CZ_cline_params["sc", isl], sh = CZ_cline_params["sh", isl],
  #                               sw = CZ_cline_params["sw", isl])[,"z_x"])
  # fml_sd = as.numeric(cline_2c3s(position = pos, sex = "female",
  #                                cl = CZ_cline_params["cl", isl], cr = CZ_cline_params["cr", isl],
  #                                wl = exp(CZ_cline_params["lwl", isl]), wr = exp(CZ_cline_params["lwr", isl]),
  #                                crab = CZ_cline_params["crab", isl], wave = CZ_cline_params["wave", isl],
  #                                zs_c = CZ_cline_params["zs_c", isl], zs_w = CZ_cline_params["zs_w", isl],
  #                                sc = CZ_cline_params["sc", isl], sh = CZ_cline_params["sh", isl],
  #                                sw = CZ_cline_params["sw", isl])[,"s_x"])
  fml_dtr = rnorm(n = 1000, mean = abs(fml_cline[, 'phen_cline']), sd = fml_cline[, 'sd_cline'])
  for (f in seq_along(fml_dtr)) {
    success=FALSE
    i=1
    fem = fml_dtr[f]
    # fem = fml_dtr[2]
    while (!success) {
      # m = sample(ml_dtr, 1, replace = FALSE)
      mpos = pos + rnorm(n=1, mean=0, sd=1.5)
      mal_cline = cline_sims(phen = -log(CZ_data[CZ_data$shore==islands[isl],]$length_mm),
                             position = mpos, sex = 'male',
                             cl = cline_pars[[isl]]['cl', 'Estimate'], cr = cline_pars[[isl]]['cr', 'Estimate'],
                             lwl = cline_pars[[isl]]['lwl', 'Estimate'], lwr = cline_pars[[isl]]['lwr', 'Estimate'],
                             crab = cline_pars[[isl]]['crab', 'Estimate'], wave = cline_pars[[isl]]['wave', 'Estimate'],
                             zs_c = cline_pars[[isl]]['zs_c', 'Estimate'], zs_w = cline_pars[[isl]]['zs_w', 'Estimate'],
                             sc = cline_pars[[isl]]['sc', 'Estimate'], shl = cline_pars[[isl]]['shl', 'Estimate'],
                             sh = cline_pars[[isl]]['sh', 'Estimate'], sw = cline_pars[[isl]]['sw', 'Estimate'])
      m = rnorm(n = 1, mean = abs(mal_cline[, 'phen_cline']), sd = mal_cline[, 'sd_cline'])
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
        cat(islands[isl], "at position", pos, ": male size", m, "mated female size", fem, ".\n")
      }
    }
    write.table(YN, file = paste0("tables/", pref_out, islands[isl], "_", round(pos), "_sim_YN.csv"), append = TRUE,
                sep = ",", row.names = FALSE, col.names = FALSE)
    bar[[f]] = YN
    YN = data.frame()
  }
  bar = rbindlist(bar)
  Ycor = with(data = bar[bar$mountYN==1, ], cor(female, male, method = "pearson"))
  Ymean = with(data = bar[bar$mountYN==1, ], mean(male))
  Yvar = with(data = bar[bar$mountYN==1, ], var(male))
  YNmean = mean(bar$male)
  YNvar = var(bar$male)
  AM_SS = data.frame(am_r = Ycor, male_mated_mean = Ymean, male_mated_var = Yvar, male_all_mean = YNmean,
                     male_all_var = YNvar)
  write.table(AM_SS, file = paste0("tables/", pref_out, "stats/", islands[isl], "_", round(pos), "_AM_SS_stats_", run, ".csv"),
              append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
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
numrun = opt$numrun
# numrun = 3
map(1:numrun, function(n) {
  cat("Running run number", n, "...\n")
  map(seq_along(islands), function(x) {
    map(isl_pos[[x]], function(y) {
      cline_pos = y
      write.table(data.frame(male=as.numeric(), female=as.numeric(), sk_prob=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
                  file = paste0("tables/", pref_out, islands[x], "_", round(cline_pos), "_sim_YN.csv"), sep = ",",
                  row.names = FALSE, col.names = TRUE)
      write.table(data.frame(am_r=as.numeric(), male_mated_mean=as.numeric(), male_mated_var=as.numeric(),
                             male_all_mean=as.numeric(), male_all_var=as.numeric(), stringsAsFactors=FALSE),
                  file = paste0("tables/", pref_out, "stats/", islands[x], "_", round(cline_pos),
                                "_AM_SS_stats_", n, ".csv"), sep = ",", row.names = FALSE, col.names = TRUE)
      simYN = possibly(sim_mat, otherwise = "Missing snails")
      outYN = simYN(pos = cline_pos, isl = x, run = n)
      # return(outYN)
      return(cat("island", islands[x], "at position", cline_pos, "has been simulated ...\n"))
    })
  })
  return(cat("*** All encounters have been simulated successfully! ***\n"))
})
# CZ_simYN = read.csv(paste0("tables/", pref_out, "CZA", "_", "9", "_sim_YN.csv"))
# with(data = CZ_simYN[CZ_simYN$mountYN==1, ], cor(female, male, method = "pearson"))
# with(data = CZ_simYN[CZ_simYN$mountYN==1, ], mean(male))
# with(data = CZ_simYN[CZ_simYN$mountYN==1, ], var(male))
# mean(CZ_simYN$male)
# var(CZ_simYN$male)

# map(seq_along(islands), function(x) {
#   map(isl_pos[[x]], function(y) {
#     cline_pos = y
#     write.table(data.frame(male=as.numeric(), female=as.numeric(), sk_prob=as.numeric(), mountYN=as.integer(), stringsAsFactors=FALSE),
#                 file = paste0("tables/", pref_out, islands[x], "_", round(cline_pos), "_sim_YN.csv"), sep = ",",
#                 row.names = FALSE, col.names = TRUE)
#     simYN = possibly(sim_mat, otherwise = "Missing snails")
#     outYN = simYN(pos = cline_pos, isl = islands[x])
#     # return(outYN)
#     return(cat("island", islands[x], "at position", cline_pos, "has been simulated ...\n"))
#   })
# })

########################
# mating summary stats #
########################
CZs_join_runs = function(isls, pos) {
  run_fl = list.files(path = paste0("tables/", pref_out, "stats"), pattern = paste0(isls, "_", pos, "_"),
                      full.names = TRUE)
  lapply(1:numrun, function(r) {read.csv(run_fl[r])})
}
CZ_am_ss_fig = lapply(seq_along(islands), function(i) {
  cat("Calculating summary statistics for", islands[i], "...\n")
  rbindlist(lapply(isl_pos[[i]], function(p) {
    CZs_summ = apply(rbindlist(CZs_join_runs(isls = islands[i], pos = p)), 2, CI)
    CZs_stats = data.frame(AM = CZs_summ[, 'am_r'],
                           DSS = CZs_summ[, 'male_mated_mean']-CZs_summ[, 'male_all_mean'],
                           SSS = CZs_summ[, 'male_mated_var']-CZs_summ[, 'male_all_var'])
    CZs_dfplot = data.frame(position=rep(p, ncol(CZs_stats)),
                            low_val=t(CZs_stats)[, 'lower'],
                            mean_val=t(CZs_stats)[, 'mean'],
                            upp_val=t(CZs_stats)[, 'upper'],
                            figure=colnames(CZs_stats))
    return(CZs_dfplot)
  }))
})

# CZs_am = function(isls, run) {
#   lapply(seq_along(isls), function(i) {
#     lapply(seq_along(isl_pos[[i]]), function(p) {
#       c_files = list.files(path = paste0("tables/", pref_out, "stats"), pattern = paste(isls[i], isl_pos[[i]][p], sep = "_"),
#                            full.names = TRUE)
#       run_df = lapply(1:run, function(r) {
#         read.csv(c_files[r])
#       })
#     })
#   })
#   c_df = lapply(seq_along(isls), function(i) {
#     lapply(1:length(c_files[[i]]), function(f) {
#       c_pos = strsplit(basename(c_files[[i]][f]), split = "_")[[1]][2]
#       isl_df = read.csv(c_files[[i]][f])
#       isl_df$male2 = isl_df$male^2
#       p_cor = cor(isl_df[isl_df$mountYN==1, ]$female, isl_df[isl_df$mountYN==1, ]$male, method = "pearson")
#       # ss_mod = glm(mountYN ~ male + male2, family = binomial(link = "logit"), data = isl_df)
#       # p_ss_idx = which.max(fitted(ss_mod))
#       pos_df = cbind(isl_df, position = as.integer(rep(c_pos, nrow(isl_df))), am_cor = rep(p_cor, nrow(isl_df)))
#     })
#   })
#   return(c_df)
# }
# CZs_am_df = CZs_am(isls = islands, run = numrun)

# CZs_ss = function(cz_dat) {
#   ss_mod = glm(mountYN ~ male + male2, family = binomial(link = "logit"), data = cz_dat)
#   dir_ss = coef(ss_mod)["male"]
#   sta_ss = coef(ss_mod)["male2"]
#   p_ss_idx = which.max(fitted(ss_mod))
#   return(cbind(cz_dat, male_mx=rep(cz_dat$male[p_ss_idx], nrow(cz_dat)), male_av=rep(mean(cz_dat$male), nrow(cz_dat)),
#                dss=rep(dir_ss, nrow(cz_dat)), sss=rep(sta_ss, nrow(cz_dat)), pss=fitted(ss_mod)))
# }
#
# CZs_am_ss = lapply(seq_along(islands), function(i) {
#   lapply(seq_along(CZs_am_df[[i]]), function(d) {
#     CZs_ss_poss = possibly(CZs_ss, otherwise = data.frame(male=as.numeric(), female=as.numeric(), sk_prob=as.numeric(),
#                                                           mountYN=as.integer(), male2=as.numeric(),
#                                                           position=as.integer(), am_cor=as.numeric(),
#                                                           male_mx=as.numeric(), male_av=as.numeric(),
#                                                           dss=as.numeric(), sss=as.numeric(), pss=as.numeric(),
#                                                           stringsAsFactors=FALSE))
#     CZs_ss_poss(cz_dat = CZs_am_df[[i]][[d]])
#   })
# })
#
# CZ_am_fun = function(data, isls) {
#   CZ_clcr = rbindlist(data[[isls]])
#   return(CZ_clcr)
# }
#
# CZ_am_tot = lapply(seq_along(islands), function(i) {
#   CZ_am_fun(data = CZs_am_ss, isls = i)
# })
# lapply(CZ_am_tot, head)

# CZ_am_ss_fig = lapply(seq_along(islands), function(f) {
#   am4fig = CZ_am_tot[[f]][, c("position", "am_cor")]
#   am4fig = mutate(am4fig, figure="am")
#   colnames(am4fig)[2] = "par_value"
#   dss4fig = CZ_am_tot[[f]][, c("position", "dss")]
#   dss4fig = mutate(dss4fig, figure="dss")
#   colnames(dss4fig)[2] = "par_value"
#   sss4fig = CZ_am_tot[[f]][, c("position", "sss")]
#   sss4fig = mutate(sss4fig, figure="sss")
#   colnames(sss4fig)[2] = "par_value"
#   df4fig = rbind(am4fig, dss4fig, sss4fig)
#   return(distinct(df4fig))
# })

#########################
# figures summary stats #
#########################
cat("Preparing dataset for figure of the mating summary statistics ...\n")
CZ_dss_fig = lapply(CZ_am_ss_fig, function(x) {
  x[x$figure=='DSS', ]
})
CZ_sss_fig = lapply(CZ_am_ss_fig, function(x) {
  x[x$figure=='SSS', ]
})
CZ_am_fig = lapply(CZ_am_ss_fig, function(x) {
  x[x$figure=='AM', ]
})
CZs_cline_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZs_phen_cline[[pl]]) +
    facet_wrap(~figure, nrow = 1) +
    geom_vline(xintercept = cline_pars[[pl]]['cl', 'Estimate'], linetype = "dashed") +
    geom_vline(xintercept = cline_pars[[pl]]['cr', 'Estimate'], linetype = "dashed") +
    scale_color_manual(values = c("red", "blue")) +
    scale_fill_manual(values = c("red", "blue")) +
    geom_ribbon(aes(x=position, ymin=abs(phen_cline)-sd_cline, ymax=abs(phen_cline)+sd_cline, fill=sex), alpha=0.15) +
    geom_errorbar(data = CZs_bin_cline[[pl]], aes(x=position[, 'mean'],
                                                  ymin=log_len[, 'lower'],
                                                  ymax=log_len[, 'upper']), alpha=0.4, width=2) +
    geom_point(data = CZs_bin_cline[[pl]], aes(x = position[, 'mean'],
                                               y = log_len[, 'mean'],
                                               col=Group.2), size=0.7) +
    geom_line(aes(position, abs(phen_cline), col=sex), size=0.9, alpha=0.7) +
    labs(x = '', y = 'ln size', fill='', col='') +
    theme(legend.position = 'top',
          strip.text = element_text(face="bold", size=13),
          strip.background = element_rect(fill="lightblue", colour="black",size=1),
          legend.text = element_text(size = 13),
          axis.title.y = element_text(face = "bold", size = 9))
})
CZs_dss_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZ_dss_fig[[pl]]) +
    facet_wrap(~figure, nrow = 1) +
    geom_vline(xintercept = cline_pars[[pl]]['cl', 'Estimate'], linetype = "dashed") +
    geom_vline(xintercept = cline_pars[[pl]]['cr', 'Estimate'], linetype = "dashed") +
    geom_errorbar(aes(x=position, ymin=low_val, ymax=upp_val), alpha=0.4, width=2) +
    geom_point(aes(x = position, y = mean_val), size=0.8) +
    geom_line(aes(x = position, y = mean_val), size=0.5) +
    labs(x = '', y = paste0('mated males - all males\n(mean size)')) +
    ylim(c(-0.2, 0.2)) +
    theme(strip.text = element_text(face="bold", size=12),
          strip.background = element_rect(fill="lightblue", colour="black",size=1),
          axis.title.y = element_text(face = "bold", size = 9))
})
CZs_sss_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZ_sss_fig[[pl]]) +
    facet_wrap(~figure, nrow = 1) +
    geom_vline(xintercept = cline_pars[[pl]]['cl', 'Estimate'], linetype = "dashed") +
    geom_vline(xintercept = cline_pars[[pl]]['cr', 'Estimate'], linetype = "dashed") +
    geom_errorbar(aes(x=position, ymin=low_val, ymax=upp_val), alpha=0.4, width=2) +
    geom_point(aes(x = position, y = mean_val), size=0.8) +
    geom_line(aes(x = position, y = mean_val), size=0.5) +
    labs(x = '', y = paste0('mated males - all males\n(variance size)')) +
    ylim(c(-0.1, 0.1)) +
    theme(strip.text = element_text(face="bold", size=12),
          strip.background = element_rect(fill="lightblue", colour="black",size=1),
          axis.title.y = element_text(face = "bold", size = 9))
})
CZs_am_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZ_am_fig[[pl]]) +
    facet_wrap(~figure, nrow = 1) +
    geom_vline(xintercept = cline_pars[[pl]]['cl', 'Estimate'], linetype = "dashed") +
    geom_vline(xintercept = cline_pars[[pl]]['cr', 'Estimate'], linetype = "dashed") +
    geom_errorbar(aes(x=position, ymin=low_val, ymax=upp_val), alpha=0.4, width=2) +
    geom_point(aes(x = position, y = mean_val), size=0.8) +
    geom_line(aes(x = position, y = mean_val), size=0.5) +
    labs(x = paste0(islands[pl], " shore position"), y = 'r') +
    ylim(c(0,1)) +
    theme(strip.text = element_text(face="bold", size=12),
          strip.background = element_rect(fill="lightblue", colour="black",size=1),
          axis.title.y = element_text(face = "bold", size = 9),
          axis.title.x = element_text(face = "bold"))
})
CZs_cline_am_ss_plot = lapply(seq_along(islands), function(x) {
  CZs_cline_plot[[x]] + CZs_dss_plot[[x]] + CZs_sss_plot[[x]] + CZs_am_plot[[x]] + plot_layout(ncol = 1, heights = c(2,2,2,2))
})
lapply(seq_along(islands), function(s) {
  cat("Saving", paste0("figures/", pref_out, islands[s], "_cline_am_ss.png"), "...\n")
  ggsave(filename = paste0("figures/", pref_out, islands[s], "_cline_am_ss.png"),
         plot = CZs_cline_am_ss_plot[[s]])
})

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
