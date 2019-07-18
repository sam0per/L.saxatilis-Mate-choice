rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "purrr", "reshape2", "pracma")
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


gaus_skew = readRDS(paste0("models/", pref_out, basename(dirname(pref_out)),".rds"))
# gaus_skew@model_pars
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

colnames(CZ_data)
pdf(paste0("figures/", pref_out, basename(dirname(pref_out)),"_fit.pdf"), width=8, height=7)
ggplot(data = CZ_data, aes(x = size_ratio, y = stan_yhat)) +
  facet_wrap(~shore) +
  geom_line(aes(col=test_sex), size=1.5) +
  # geom_line(aes(x = size_ratio, y = stan_yhat_lci), alpha=0.5, linetype = "dashed") +
  # geom_line(aes(x = size_ratio, y = stan_yhat_uci), alpha=0.5, linetype = "dashed")
  geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = test_sex), alpha=0.1) +
  labs(x="ln female size - ln male size", y="probability of mating", col="", fill ="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.text = element_text(size = 15,face = "bold"),
        axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))
dev.off()

skew_pars = read.csv("tables/gaus_skew/SKEW/gaus_skew_params.csv", sep = ";")
(skew_pars = column_to_rownames(skew_pars, var = "parameter"))
pmat = function(b0, b1, c, d, alpha, dat) {
  b0 + b1 * exp(-0.5 * ((dat - c) / d)^2) * (1 + erf(alpha * (dat - c) / (1.414214 * d)))
}
range(pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
           alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio))
# max(pmat(b0 = 0, b1 = 0.36, c = 0.23, d = 0.74, alpha = 1.61, dat = CZ_data$size_ratio))
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
