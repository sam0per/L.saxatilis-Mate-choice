rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "purrr", "reshape2", "pracma", "viridis", "data.table",
              "Cairo")
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
lapply(.packages, require, character.only=TRUE)

CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
#############################
######### Figure S1 #########
#############################
pref_out = lapply(c("B_all", "C_all", "D_all", "ALPHA_all"), function(m) {
  paste0("gaus_skew/", m, "/gaus_skew_hier_")
})

CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
gaus_skew = lapply(seq_along(pref_out), function(x) {
  gmod = readRDS(paste0("models/", pref_out[[x]], basename(dirname(pref_out[[x]])),".rds"))
  # print(gmod@model_pars[1:5])
  stan_val = summary(gmod, pars = gmod@model_pars[1:5])$summary
  if (gmod@model_pars[1] == "alpha_intercept") {
    hyp_px = substr(gmod@model_pars[1], start = 1, stop = 6)
    greek_px = expression("\u03b1")
  } else {
    hyp_px = substr(gmod@model_pars[1], start = 1, stop = 1)
    greek_px = hyp_px
  }
  row.names(stan_val)[grepl(paste0(hyp_px, "_", "coeff"), row.names(stan_val))] = paste0(hyp_px, "_", colnames(CZ_matrix))
  stan_val = rownames_to_column(as.data.frame(stan_val), var="parameter")
  stan_val$parameter = factor(stan_val$parameter,
                              levels = stan_val$parameter[order(sort(stan_val$parameter,
                                                                     decreasing = TRUE))])
  stan_val = mutate(.data = stan_val, fig = paste(greek_px, "hierarchical", sep = "-"))
  return(stan_val)
})
FigS1 = lapply(seq_along(pref_out), function(x) {
  ggplot(data = gaus_skew[[x]]) +
    facet_wrap(~fig) +
    geom_errorbarh(aes(xmax = `2.5%`, xmin = `97.5%`, y = parameter, height = 0, col = "95% CI"), size = 0.8) +
    geom_errorbarh(aes(xmax = `25%`, xmin = `75%`, y = parameter, height = 0, col = "50% CI"), size = 1.3) +
    geom_point(aes(x = mean, y = parameter, col = "odot"), size = 2.8) +
    geom_point(aes(x = mean, y = parameter, col = "idot"), size = 1.8) +
    scale_color_viridis_d(begin = 0, end = 1, option = "A", direction = -1) +
    theme(legend.position = "top", axis.title = element_blank(), legend.title = element_blank(),
          strip.text = element_text(face="bold", size=13),
          strip.background = element_rect(fill="lightblue", colour="black",size=1))
})
lapply(seq_along(pref_out), function(x) {
  ggsave(filename = paste0("figures/", pref_out[[x]], basename(dirname(pref_out[[x]])), "_coeff.pdf"),
         plot = FigS1[[x]], device = cairo_pdf, width = 7, height = 7)
})
#############################
######### Figure S2 #########
#############################
pref_out = "gaus_skew/BCDG/gaus_skew_hier_"
CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
gmod = readRDS(paste0("models/", pref_out, basename(dirname(pref_out)),".rds"))
gmod@model_pars
hyp_nm = c("b", "c", "d", "alpha")
gaus_skew = lapply(seq_along(hyp_nm), function(h) {
  coeff_nm = gmod@model_pars[grepl(pattern = paste0("^", hyp_nm[h]), x = gmod@model_pars)][-3]
  stan_val = summary(gmod, pars = coeff_nm)$summary
  if (coeff_nm[1] == "alpha_intercept") {
    # hyp_px = substr(gmod@model_pars[1], start = 1, stop = 6)
    greek_px = expression("\u03b1")
  } else {
    # hyp_px = substr(gmod@model_pars[1], start = 1, stop = 1)
    greek_px = hyp_nm[h]
  }
  row.names(stan_val)[grepl(paste0(hyp_nm[h], "_", "coeff"), row.names(stan_val))] = paste0(hyp_nm[h], "_", colnames(CZ_matrix))
  stan_val = rownames_to_column(as.data.frame(stan_val), var="parameter")
  stan_val$parameter = factor(stan_val$parameter,
                              levels = stan_val$parameter[order(sort(stan_val$parameter,
                                                                     decreasing = TRUE))])
  stan_val = mutate(.data = stan_val, fig = paste("Hyperparameter", greek_px, sep = " "))
  return(stan_val)
})
FigS2 = lapply(seq_along(hyp_nm), function(x) {
  ggplot(data = gaus_skew[[x]]) +
    facet_wrap(~fig) +
    geom_errorbarh(aes(xmax = `2.5%`, xmin = `97.5%`, y = parameter, height = 0, col = "95% CI"), size = 0.8) +
    geom_errorbarh(aes(xmax = `25%`, xmin = `75%`, y = parameter, height = 0, col = "50% CI"), size = 1.3) +
    geom_point(aes(x = mean, y = parameter, col = "odot"), size = 2.8) +
    geom_point(aes(x = mean, y = parameter, col = "idot"), size = 1.8) +
    scale_color_viridis_d(begin = 0, end = 1, option = "A", direction = -1) +
    theme(legend.position = "top", axis.title = element_blank(), legend.title = element_blank(),
          strip.text = element_text(face="bold", size=13),
          strip.background = element_rect(fill="lightblue", colour="black",size=1))
})
lapply(seq_along(hyp_nm), function(x) {
  ggsave(filename = paste0("figures/", pref_out, basename(dirname(pref_out)), "_", hyp_nm[x], "_coeff.pdf"),
         plot = FigS2[[x]], device = cairo_pdf, width = 7, height = 7)
})
##############################
##############################
######### island eff #########
##############################
pref_out = "gaus_skew/Bhyp_shore/gaus_skew_hier_"
pref_out = "gaus_skew/Chyp_shore/gaus_skew_hier_"
pref_out = "gaus_skew/Dhyp_shore/gaus_skew_hier_"
pref_out = "gaus_skew/ALPHAhyp_shore/gaus_skew_hier_"
gaus_skew = readRDS(paste0("models/", pref_out, basename(dirname(pref_out)),".rds"))
CZ_data$stan_yhat = summary(gaus_skew, pars = c("y_hat"))$summary[,'mean']
CZ_data$stan_yhat_uci = summary(gaus_skew, pars = c("y_hat"))$summary[,'97.5%']
CZ_data$stan_yhat_lci = summary(gaus_skew, pars = c("y_hat"))$summary[,'2.5%']
hyp_isl = tibble(size_ratio = CZ_data$size_ratio, shore = CZ_data$shore, ref_ecotype = CZ_data$ref_ecotype,
                 stan_yhat = CZ_data$stan_yhat, stan_yhat_lci = CZ_data$stan_yhat_lci,
                 stan_yhat_uci = CZ_data$stan_yhat_uci, hyp = "b-only-island")
hyp_isl = rbind(hyp_isl,
                tibble(size_ratio = CZ_data$size_ratio, shore = CZ_data$shore, ref_ecotype = CZ_data$ref_ecotype,
                       stan_yhat = CZ_data$stan_yhat, stan_yhat_lci = CZ_data$stan_yhat_lci,
                       stan_yhat_uci = CZ_data$stan_yhat_uci, hyp = "a-only-island"))
table(hyp_isl$hyp)
factor(hyp_isl$hyp)
# hyp_isl[hyp_isl$hyp=='hyperparameter b', ]$hyp = "b-only-island"
pdf(paste0("figures/gaus_skew/one_only_hyp_shore_fit.pdf"), width=8, height=6)
ggplot(data = hyp_isl, aes(x = size_ratio, y = stan_yhat)) +
  facet_wrap(~hyp) +
  geom_line(aes(col=shore), size=1.5) +
  # geom_line(aes(x = size_ratio, y = stan_yhat_lci), alpha=0.5, linetype = "dashed") +
  # geom_line(aes(x = size_ratio, y = stan_yhat_uci), alpha=0.5, linetype = "dashed")
  geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = shore), alpha=0.1) +
  labs(x="ln female size - ln male size", y="probability of mating", col="", fill ="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.text = element_text(size = 15,face = "bold"),
        legend.position = "top",
        axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))
dev.off()
#############################
#############################
######### Figure T1 #########
#############################
pref_out = "gaus_skew/Chyp_eco/gaus_skew_hier_"
pref_out = "gaus_skew/Bhyp_islsex/gaus_skew_hier_"
pdf(paste0("figures/", pref_out, basename(dirname(pref_out)),"_fit.pdf"), width=7, height=7)
ggplot(data = CZ_data, aes(x = size_ratio, y = stan_yhat)) +
  facet_wrap(~shore) +
  # geom_line(aes(col=test_sex), size=1.5) +
  geom_line(aes(col=ref_ecotype), size=1.5) +
  # geom_line(aes(x = size_ratio, y = stan_yhat_lci), alpha=0.5, linetype = "dashed") +
  # geom_line(aes(x = size_ratio, y = stan_yhat_uci), alpha=0.5, linetype = "dashed")
  # geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = test_sex), alpha=0.25) +
  geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = ref_ecotype), alpha=0.25) +
  # scale_color_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("darkgoldenrod3", "darkorchid4")) +
  labs(x="ln female size - ln male size", y="probability of mating", col="", fill ="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.position = 'top',
        strip.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill="lightblue", colour="black",size=1),
        legend.text = element_text(size = 13),
        axis.title = element_text(face = "bold", size = 15))
dev.off()
############################
############################
######### Figure 2 #########
############################
skew_pars = read.csv("tables/gaus_skew/SKEW/gaus_skew_params.csv", sep = ";")
(skew_pars = column_to_rownames(skew_pars, var = "parameter"))
pmat = function(b0, b1, c, d, alpha, dat) {
  b0 + b1 * exp(-0.5 * ((dat - c) / d)^2) * (1 + erf(alpha * (dat - c) / (1.414214 * d)))
}
range(pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
           alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio))
# max(pmat(b0 = 0, b1 = 0.36, c = 0.23, d = 0.74, alpha = 1.61, dat = CZ_data$size_ratio))
pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
     alpha = skew_pars["alpha","mean"], dat = -0.5)
CImat = read.csv("tables/gaus_skew/SKEW/gaus_skew_mat.csv", sep = ";")
max(CImat$stan_yhat_lci)
max(CImat$stan_yhat_uci)
max(CImat$stan_yhat)

pmat_df = tibble(p = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                          alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                 pp = pmat(b0 = 0.15, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                           alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                 ppp = pmat(b0 = 0.25, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                            alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                 size_ratio = CZ_data$size_ratio, par = "b0")
pmat_df = rbind(pmat_df,
                tibble(p = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                                alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       pp = pmat(b0 = 0.01, b1 = 0.25, c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                                 alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       ppp = pmat(b0 = 0.01, b1 = 0.5, c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                                  alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       size_ratio = CZ_data$size_ratio, par = "b1"))
pmat_df = rbind(pmat_df,
                tibble(p = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                                alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       pp = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = 0.4, d = skew_pars["d","mean"],
                                 alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       ppp = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = -0.6,d = skew_pars["d","mean"],
                                  alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       size_ratio = CZ_data$size_ratio, par = "c"))
pmat_df = rbind(pmat_df,
                tibble(p = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                                alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       pp = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"], d = 0.3,
                                 alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       ppp = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"], d = 1.4,
                                  alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       size_ratio = CZ_data$size_ratio, par = "d"))
pmat_df = rbind(pmat_df,
                tibble(p = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
                                alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
                       pp = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"], d = skew_pars["d","mean"],
                                 alpha = 0, dat = CZ_data$size_ratio),
                       ppp = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"], d = skew_pars["d","mean"],
                                  alpha = 4, dat = CZ_data$size_ratio),
                       size_ratio = CZ_data$size_ratio, par = "alpha"))

pref_out = "gaus_skew/SKEW/gaus_skew_"
pdf(paste0("figures/", pref_out, basename(dirname(pref_out)),"_pars.pdf"), width=8, height=6)
ggplot(data = pmat_df, aes(x = size_ratio, y = p)) +
  facet_wrap(~par, labeller = label_parsed) +
  geom_line(col='black') +
  geom_line(aes(x = size_ratio, y = pp), col='orange') +
  geom_line(aes(x = size_ratio, y = ppp), col='orange') +
  labs(x="ln female size - ln male size", y="probability of mating") +
  theme(axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))
dev.off()
############################
############################
######### Table S4 #########
############################
lapply(seq_along(islands), function(s) {
  cat("Saving", paste0("tables/clines/", islands[s], "_cline_params_2sh.csv"), "...\n")
  write.csv(x = cline_pars[[s]], file = paste0("tables/clines/", islands[s], "_cline_params_2sh.csv"), row.names = TRUE)
})
CZ_cline = read.csv("tables/clines/CZD_cline_params_2sh.csv", row.names = 1)
round(CZ_cline, 2)
# sqrt(0.118^2 + 0.188^2 + 0.25*(0.339^2-0.118^2))
s_centre = function(dt, sc, sh, sw) {
  sqrt(dt[as.character(sc),1]^2 + dt[as.character(sh),1]^2 + 0.25*(dt[as.character(sw),1]^2-dt[as.character(sc),1]^2))
}
CZ_cline["shl", 1] = round(s_centre(dt = CZ_cline, sc = "sc", sh = "shl", sw = "sw"), 2)
CZ_cline["sh", 1] = round(s_centre(dt = CZ_cline, sc = "sc", sh = "sh", sw = "sw"), 2)
CZ_cline["zs_c", 1] = round(CZ_cline["zs_c", 1]+CZ_cline["crab", 1], 2)
CZ_cline["zs_c", 2] = round(sqrt(CZ_cline["zs_c", 2]^2 + CZ_cline["crab", 2]^2), 2)
CZ_cline["zs_w", 1] = round(CZ_cline["zs_w", 1]+CZ_cline["wave", 1], 2)
CZ_cline["zs_w", 2] = round(sqrt(CZ_cline["zs_w", 2]^2 + CZ_cline["wave", 2]^2), 2)

CI_95 = function(dt, nm_par) {
  if (nm_par!="cl" & nm_par!="cr") {
    epar = exp(abs(dt[nm_par,1]))
    lci = exp(abs(dt[nm_par,1] - (2 * dt[nm_par,2])))
    uci = exp(abs(dt[nm_par,1] + (2 * dt[nm_par,2])))
  } else {
    lci = dt[nm_par,1] - (2 * dt[nm_par,2])
    uci = dt[nm_par,1] + (2 * dt[nm_par,2])
    epar = dt[nm_par,1]
  }
  paste0(round(epar, 2), " [", round(lci, 2), "- ", round(uci, 2), "]")
  # paste0(dt[as.character(epar),1], " (95% CI = ", lci, " to ", uci)
}
CI_95(dt = CZ_cline, nm_par = "wave")
lapply(row.names(CZ_cline), function(p) {
  CI_95(dt = CZ_cline, nm_par = p)
})
############################
############################
######### ESEB2019 #########
############################
gmod = readRDS(paste0("models/", pref_out, basename(dirname(pref_out)),".rds"))
CZ_data$stan_yhat = summary(gmod, pars = c("y_hat"))$summary[,'mean']
CZ_data$stan_yhat_uci = summary(gmod, pars = c("y_hat"))$summary[,'97.5%']
CZ_data$stan_yhat_lci = summary(gmod, pars = c("y_hat"))$summary[,'2.5%']

pdf(paste0("../../conf/ESEB/Turku2019/figures/Bhyp_shore_fit.pdf"), width=8, height=7)
ggplot(data = CZ_data, aes(x = size_ratio, y = stan_yhat)) +
  # facet_wrap(~hyp) +
  geom_line(aes(col=shore), size=2) +
  # geom_line(aes(x = size_ratio, y = stan_yhat_lci), alpha=0.5, linetype = "dashed") +
  # geom_line(aes(x = size_ratio, y = stan_yhat_uci), alpha=0.5, linetype = "dashed")
  geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = shore), alpha=0.15) +
  labs(x="ln female size - ln male size", y="probability of mating", col="", fill ="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,1)) +
  theme(legend.text = element_text(size = 15, face = "bold"),
        legend.position = "top",
        axis.title = element_text(face = "bold", size = 18))
dev.off()


CZs_dss_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZ_dss_fig[[pl]]) +
    facet_wrap(~figure, nrow = 1) +
    geom_vline(xintercept = cline_pars[[pl]]['cl', 'Estimate'], linetype = "dashed") +
    geom_vline(xintercept = cline_pars[[pl]]['cr', 'Estimate'], linetype = "dashed") +
    geom_errorbar(aes(x=position, ymin=low_val, ymax=upp_val), alpha=0.4, width=2) +
    geom_point(aes(x = position, y = mean_val), size=0.8) +
    geom_line(aes(x = position, y = mean_val), size=0.5) +
    labs(x = paste0(islands[pl], " shore position"), y = paste0('mated males - all males\n(mean size)')) +
    ylim(c(-0.2, 0.2)) +
    theme(strip.text = element_text(face="bold", size=12),
          strip.background = element_rect(fill="lightblue", colour="black",size=1),
          axis.title.y = element_text(face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 15))
})
CZs_sss_plot = lapply(seq_along(islands), function(pl) {
  ggplot(data = CZ_sss_fig[[pl]]) +
    facet_wrap(~figure, nrow = 1) +
    geom_vline(xintercept = cline_pars[[pl]]['cl', 'Estimate'], linetype = "dashed") +
    geom_vline(xintercept = cline_pars[[pl]]['cr', 'Estimate'], linetype = "dashed") +
    geom_errorbar(aes(x=position, ymin=low_val, ymax=upp_val), alpha=0.4, width=2) +
    geom_point(aes(x = position, y = mean_val), size=0.8) +
    geom_line(aes(x = position, y = mean_val), size=0.5) +
    labs(x = paste0(islands[pl], " shore position"), y = paste0('mated males - all males\n(variance size)')) +
    ylim(c(-0.1, 0.1)) +
    theme(strip.text = element_text(face="bold", size=12),
          strip.background = element_rect(fill="lightblue", colour="black",size=1),
          axis.title.y = element_text(face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 15))
})
lapply(seq_along(islands), function(s) {
  cat("Saving", paste0("figures/", pref_out, islands[s], "_sss.pdf"), "...\n")
  ggsave(filename = paste0("figures/", pref_out, islands[s], "_sss.pdf"),
         plot = CZs_sss_plot[[s]], device = "pdf", width = 8, height = 4)
})
############################
############################
############################