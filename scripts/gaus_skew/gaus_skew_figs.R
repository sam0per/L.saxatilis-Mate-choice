rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "purrr", "reshape2", "pracma", "viridis", "data.table")
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
lapply(.packages, require, character.only=TRUE)

CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")

pref_out = "gaus_skew/Bhyp_shore/gaus_skew_hier_"
pref_out = "gaus_skew/Chyp_shore/gaus_skew_hier_"
pref_out = "gaus_skew/Dhyp_shore/gaus_skew_hier_"
pref_out = "gaus_skew/ALPHAhyp_shore/gaus_skew_hier_"

pref_out = "gaus_skew/Bhyp_eco/gaus_skew_hier_"
pref_out = "gaus_skew/Chyp_eco/gaus_skew_hier_"

pref_out = "gaus_skew/Bhyp_islsex/gaus_skew_hier_"

pref_out = "gaus_skew/B_all/gaus_skew_hier_"

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
  } else {
    hyp_px = substr(gmod@model_pars[1], start = 1, stop = 2)
  }
  row.names(stan_val)[grepl(paste0(hyp_px, "coeff"), row.names(stan_val))] = paste0(hyp_px, colnames(CZ_matrix))
  stan_val = rownames_to_column(as.data.frame(stan_val), var="parameter")
  return(stan_val)
})

lapply(seq_along(pref_out), function(x) {
  gaus_skew[[x]]$parameter = factor(gaus_skew[[x]]$parameter,
                                    levels = gaus_skew[[x]]$parameter[order(sort(gaus_skew[[x]]$parameter, decreasing = TRUE))])
  return(gaus_skew[[x]])
})
# gaus_skew@model_pars[1:5]
# str(gaus_skew)
# stan_yhat = summary(gaus_skew, pars = gaus_skew@model_pars[1:5])$summary
# head(stan_yhat)
# CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
# colnames(CZ_matrix)
# row.names(stan_yhat)[grepl(paste0("b_", "coeff"), row.names(stan_yhat))] = paste0("b_", colnames(CZ_matrix))
# (stan_yhat = rownames_to_column(as.data.frame(stan_yhat), var="parameter"))
# colnames(stan_yhat)
# stan_yhat$parameter = factor(stan_yhat$parameter, levels = stan_yhat$parameter[order(sort(stan_yhat$parameter, decreasing = TRUE))])
ggplot(data = stan_yhat) +
  geom_errorbarh(aes(xmax = `2.5%`, xmin = `97.5%`, y = parameter, height = 0, col = "95% CI"), size = 0.8) +
  geom_errorbarh(aes(xmax = `25%`, xmin = `75%`, y = parameter, height = 0, col = "50% CI"), size = 1.3) +
  geom_point(aes(x = mean, y = parameter, col = "odot"), size = 2.8) +
  geom_point(aes(x = mean, y = parameter, col = "idot"), size = 1.8) +
  scale_color_viridis_d(begin = 0, end = 1, option = "A", direction = -1) +
  theme(legend.position = "top", axis.title = element_blank(), legend.title = element_blank())
  

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

# pref_out = lapply(c("Bhyp_shore", "Chyp_shore", "Dhyp_shore", "ALPHAhyp_shore"), function(m) {
#   paste0("gaus_skew/", m, "/gaus_skew_hier_")
# })
# colnames(CZ_data)
# gaus_skew = lapply(seq_along(pref_out), function(x) {
#   gmod = readRDS(paste0("models/", pref_out[[x]], basename(dirname(pref_out[[x]])),".rds"))
#   # print(gmod@model_pars[1:5])
#   ### TODO: sort shore names ###
#   isl = paste0("CZ", substr(basename(dirname(pref_out[[x]])), start = 1, stop = 1))
#   ##############################
#   # CZone = CZ_data[CZ_data$shore==isl, ]
#   CZone = data.frame(shore = isl, stan_yhat = summary(gmod, pars = c("y_hat"))$summary[,'mean'],
#                      stan_yhat_uci = summary(gmod, pars = c("y_hat"))$summary[,'97.5%'],
#                      stan_yhat_lci = summary(gmod, pars = c("y_hat"))$summary[,'2.5%'])
#   return(CZone)
# })
# lapply(gaus_skew, head)
# gaus_skew[[4]]$shore = "CZD"
# CZ_stan = as.data.frame(rbindlist(gaus_skew))
# head(CZ_stan)
# CZ_stan$size_ratio = CZ_data$size_ratio
# pdf(paste0("figures/gaus_skew/one_only_hyp_shore_fit.pdf"), width=8, height=6)
# ggplot(data = CZ_stan, aes(x = size_ratio, y = stan_yhat)) +
#   # facet_wrap(~shore) +
#   geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = shore), alpha=0.1) +
#   geom_line(aes(col=shore), size=1.5) +
#   # geom_line(aes(x = size_ratio, y = stan_yhat_lci), alpha=0.5, linetype = "dashed") +
#   # geom_line(aes(x = size_ratio, y = stan_yhat_uci), alpha=0.5, linetype = "dashed")
#   labs(x="ln female size - ln male size", y="probability of mating", col="", fill ="") +
#   scale_x_continuous(breaks = seq(-1.5,1.5,1)) +
#   theme(legend.text = element_text(size = 15),
#         legend.position = "top",
#         axis.title = element_text(face = "bold", size = 15))
# dev.off()

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
