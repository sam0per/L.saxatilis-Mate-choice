rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "purrr", "reshape2", "pracma", "viridis", "data.table",
              "Cairo", "extrafont", "ggthemes")
.packagesdev = "thomasp85/patchwork"
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
############################
############################
######### Figure 1 #########
############################
# x <- seq(0, 1, length=1000)
# y <- dnorm(x, mean=0.5, sd=0.1)
# plot(x, dbeta(x, 15, 20), type="l")
# lines(x, y, lwd=1)
skew_pars = read.csv("tables/gaus_skew/SKEW/gaus_skew_params.csv", sep = ";")
(skew_pars = column_to_rownames(skew_pars, var = "parameter"))
pmat = function(b0, b1, c, d, alpha, dat) {
  b0 + b1 * exp(-0.5 * ((dat - c) / d)^2) * (1 + erf(alpha * (dat - c) / (1.414214 * d)))
}
range(pmat(b0 = 0.01, b1 = skew_pars["b1_par","mean"], c = skew_pars["c_par","mean"],d = skew_pars["d_par","mean"],
           alpha = skew_pars["alpha_par","mean"], dat = CZ_data$size_ratio))
# max(pmat(b0 = 0, b1 = 0.36, c = 0.23, d = 0.74, alpha = 1.61, dat = CZ_data$size_ratio))
pmat(b0 = 0.01, b1 = skew_pars["b1_par","mean"], c = skew_pars["c_par","mean"],d = skew_pars["d_par","mean"],
     alpha = skew_pars["alpha_par","mean"], dat = log(12.5)-log(5.2))
CImat = read.csv("tables/gaus_skew/SKEW/gaus_skew_mat.csv", sep = ";")
max(CImat$stan_yhat_lci)
max(CImat$stan_yhat_uci)
max(CImat$stan_yhat)
# exclude parameter b0
# pmat_df = tibble(p = pmat(b0 = 0.01, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
#                           alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
#                  pp = pmat(b0 = 0.15, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
#                            alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
#                  ppp = pmat(b0 = 0.25, b1 = skew_pars["b","mean"], c = skew_pars["c","mean"],d = skew_pars["d","mean"],
#                             alpha = skew_pars["alpha","mean"], dat = CZ_data$size_ratio),
#                  size_ratio = CZ_data$size_ratio, par = "b0")
pmat_pb = tibble(p = pmat(b0 = 0.01, b1 = 0.4, c = 0, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 pp = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 ppp = pmat(b0 = 0.01, b1 = 0.6, c = 0, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 size_ratio = seq(-2, to = 2, length.out = 1001), par = "b[1]")
b_1 = ggplot(data = pmat_pb, aes(x = size_ratio, y = pp)) +
  facet_wrap(~par, labeller = label_parsed) +
  geom_line(col='black') +
  geom_line(aes(x = size_ratio, y = p), col='orange') +
  geom_line(aes(x = size_ratio, y = ppp), col='orange') +
  labs(x="ln female size - ln male size", y="probability of mating") +
  theme(axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))
b_1

pmat_pc = tibble(p = pmat(b0 = 0.01, b1 = 0.5, c = -0.2, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 pp = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 ppp = pmat(b0 = 0.01, b1 = 0.5, c = 0.2, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 size_ratio = seq(-2, to = 2, length.out = 1001), par = "c")
c = ggplot(data = pmat_df, aes(x = size_ratio, y = pp)) +
  facet_wrap(~par, labeller = label_parsed) +
  geom_line(col='black') +
  geom_line(aes(x = size_ratio, y = p), col='orange') +
  geom_line(aes(x = size_ratio, y = ppp), col='orange') +
  labs(x="ln female size - ln male size", y="probability of mating") +
  theme(axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))
c

pmat_pd = tibble(p = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.4, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 pp = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 ppp = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.6, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 size_ratio = seq(-2, to = 2, length.out = 1001), par = "d")
d = ggplot(data = pmat_df, aes(x = size_ratio, y = pp)) +
  facet_wrap(~par, labeller = label_parsed) +
  geom_line(col='black') +
  geom_line(aes(x = size_ratio, y = p), col='orange') +
  geom_line(aes(x = size_ratio, y = ppp), col='orange') +
  labs(x="ln female size - ln male size", y="probability of mating") +
  theme(axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))

pmat_pa = tibble(p = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.5, alpha = -1, dat = seq(-2, to = 2, length.out = 1001)),
                 pp = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.5, alpha = 0, dat = seq(-2, to = 2, length.out = 1001)),
                 ppp = pmat(b0 = 0.01, b1 = 0.5, c = 0, d = 0.5, alpha = 1, dat = seq(-2, to = 2, length.out = 1001)),
                 size_ratio = seq(-2, to = 2, length.out = 1001), par = "alpha")
pmat_df = rbind(pmat_pb, pmat_pc, pmat_pd, pmat_pa)
# table(pmat_df$par)
# str(pmat_df)
pmat_df$par = factor(pmat_df$par, levels = c("b[1]","c","d","alpha"))
# pdf("manuscript/figures/FIG1_pars_eff.pdf", width = 25, height = 22)
# loadfonts(device="postscript")
# bitmap("manuscript/figures/FIG1_pars_eff.tiff",
#        height = 4, width = 5, units = 'in', res=600)
bitmap("manuscript/figures/FIG1_pars_eff_v2.tiff",
       height = 4, width = 5, units = 'in', res=600)
ggplot(data = pmat_df, aes(x = size_ratio, y = pp)) +
  facet_wrap(~par, labeller = label_parsed) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_line(col='black', size=1.1) +
  geom_line(aes(x = size_ratio, y = p), col='#1a9641', size=1.1) +
  geom_line(aes(x = size_ratio, y = ppp), col='orange', size=1.1) +
  labs(x="ln female size - ln male size", y="probability of mating") +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(face="italic", size=10),
        strip.background = element_rect(fill="lightblue", colour="black",size=1),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))
  
dev.off()


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
######### Figure 2 #########
############################
CZ_data = read.csv("tables/gaus_skew/SKEW/gaus_skew_mat.csv", sep = ";")

# y = CZ_data$mountYNcontact
# y_rep = rstan::extract(CZ_mat_stan_size, pars = 'y_rep', permuted = TRUE)$y_rep
# dim(y_rep)
# ppc_bars(y,y_rep)
range(CZ_data$size_ratio)
breaks = c(-2,seq(-1.5,1.5,0.1),2)
# bin = cut(CZ_data$size_ratio,breaks)
# pdf("figures/gaus_size/gaus_size_ppc_bars_grouped.pdf")
# ppc_bars_grouped(y, y_rep, bin, prob = 0, freq = FALSE)
# dev.off()

CZ_data$bin = cut(CZ_data$size_ratio,breaks)
CZ_data_bin =
  CZ_data %>%
  group_by(bin) %>%
  dplyr::summarise(mount = mean(mountYNcontact),
                   uci_mount = CI(mountYNcontact)['upper'],
                   lci_mount = CI(mountYNcontact)['lower'],
                   mean_ratio = mean(size_ratio),
                   #y_rep = mean(y_rep),
                   preds_mount = mean(stan_mountYN)) %>%
  mutate(lci_mount = replace(lci_mount, which(lci_mount<0), 0))
# summary(CZ_data_bin)

# pdf("figures/gaus_skew/SKEW/gaus_skew_preds.pdf", width=6, height=5)
# pdf("figures/gaus_skew/SKEW/gaus_skew_obs_bin.pdf", width=8, height=7)
bitmap("manuscript/figures/FIG2_SKEW.tiff",
       height = 4, width = 5, units = 'in', res=600)
ggplot(data = CZ_data) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(x = size_ratio,ymin = stan_yhat_lci, ymax = stan_yhat_uci), fill = "peachpuff") +
  geom_errorbar(data = CZ_data_bin, aes(x = mean_ratio, ymin = lci_mount, ymax = uci_mount)) +
  scale_colour_manual(values=c("blue","orange2")) +
  geom_line(aes(size_ratio, stan_yhat, col="predictions"), size=1.5) +
  geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount, col="observations"), size=2) +
  # geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount), col='blue', size=3.5) +
  labs(x="ln female size - ln male size",
       y="probability of mating", col="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.text = element_text(size = 15,face = "bold"), legend.position = 'none',
        axis.title = element_text(size = 12),
        # axis.ticks = element_line(size = 2),
        axis.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))
# grids(linetype = "dashed")
dev.off()
############################
############################
######### Figure 3 #########
############################
# dtmp = data.frame(x = rnorm(n = 10), y = rbinom(n = 10, size = 1, prob = 0.2))
# ggplot(data = dtmp, aes(x = x, y = y)) +
#   geom_rect(aes(xmin=0, xmax=0.5, ymin=0, ymax=1), fill="white") +
#   geom_point()
lapply(seq_along(islands), function(s) {
  cat("Saving", paste0("figures/", pref_out, islands[s], "_cline_am_ss.tiff"), "...\n")
  ggsave(filename = paste0("figures/", pref_out, islands[s], "_cline_am_ss.tiff"),
         plot = CZs_cline_am_ss_plot[[s]], device = "tiff", width = 4, height = 5, units = "in", dpi = 600)
})
################################
################################
######### Table S1, S2 #########
################################
install.packages("MuMIn")
library(MuMIn)
CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
top_glm = read.csv("tables/model_search/top_glm_2way_6predictors.csv")
CZ_data$size_ratio3 = CZ_data$size_ratio^3
call_top_glm = lapply(1:13, function(x) {
  glm(formula = as.character(top_glm$frml)[x], family = binomial(link = "logit"), data = CZ_data)
})
ave_top_glm = summary(model.avg(call_top_glm))
str(ave_top_glm)
write.csv(x = ave_top_glm[["coefmat.full"]], file = "tables/model_search/full_ave_top_glm_2way_6predictors.csv")
write.csv(x = ave_top_glm[["coefmat.subset"]], file = "tables/model_search/subset_ave_top_glm_2way_6predictors.csv")
##############################
##############################
######### S2 Fig S1a #########
##############################
pref_out = lapply(c("B_all", "C_all", "D_all", "ALPHA_all"), function(m) {
  paste0("gaus_skew/", m, "/gaus_skew_hier_")
})

CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
head(CZ_matrix)
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
isl = 1
hypp = c("Ball", "Call", "Dall", "ALPHAall")
gaus_skew[[isl]]
alpha_pars = gaus_skew[[4]]$parameter
# gaus_skew[[isl]]$fig = "b[1]-hierarchical\nIsland and test sex effect"
# gaus_skew[[isl]]$parameter = as.numeric(gaus_skew[[isl]]$parameter)
gaus_skew[[4]]$parameter = gaus_skew[[2]]$parameter
gaus_skew[[isl]]$parameter = gaus_skew[[3]]$parameter
gaus_skew[[isl]]$parameter = gsub(pattern = "d_", replacement = "b[1]_", x = gaus_skew[[isl]]$parameter)
gaus_skew[[isl]]$parameter = gsub(pattern = "b_", replacement = "c_", x = gaus_skew[[isl]]$parameter)
# gaus_skew[[1]]$parameter[1] = "b_intercept"
gaus_skew[[isl]]$parameter = factor(gaus_skew[[isl]]$parameter,
                                    levels = gaus_skew[[isl]]$parameter[order(sort(gaus_skew[[isl]]$parameter,
                                                                                   decreasing = TRUE))])
# gaus_skew[[1]]$parameter[1] = "b[1]-intercept"
str(gaus_skew[[isl]])
gaus_skew[[isl]]$fig = paste0(hypp[2], "-hierarchical\nIsland and transect sex effect")
tiff(filename = paste0("manuscript/figures/S2_FIG1a_",hypp[isl], ".tiff"), width = 2, height = 2, units = "in",
     res = 600)
# pdf(file = "manuscript/figures/S2_FIG1a_Call.pdf", width = 4, height = 4)
ggplot(data = gaus_skew[[isl]]) +
  facet_wrap(~fig) +
  # geom_segment(aes(x=-Inf, xend=Inf, y=10.5, yend=10.5), linetype = "dashed", size = 0.3, col = "grey80") +
  geom_segment(aes(x=0, xend=0, y=-Inf, yend=10.5), linetype = "dashed", size = 0.3, col = "grey80") +
  geom_errorbarh(aes(xmax = `2.5%`, xmin = `97.5%`, y = parameter, height = 0, col = "95% CI"), size = 0.4) +
  geom_errorbarh(aes(xmax = `25%`, xmin = `75%`, y = parameter, height = 0, col = "50% CI"), size = 0.5) +
  geom_point(aes(x = mean, y = parameter, col = "odot"), size = 0.2) +
  # geom_point(aes(x = mean, y = parameter, col = "idot"), size = 1.8) +
  scale_colour_manual(values = c("grey40", "grey60", "black")) +
  # scale_colour_brewer(type = "qual", palette = 8, direction = 1) +
  # scale_color_viridis_d(begin = 0, end = 1, option = "A", direction = -1) +
  # scale_y_continuous(breaks = seq(1,11,1)) +
  # scale_x_continuous(breaks = seq(-0.5,2,0.5)) +
  xlim(-1,3) +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_text(size = 3),
        axis.ticks = element_line(size = 0.3),
        strip.text = element_text(face="bold", size=5, margin = margin(t = 2, r = 2, b = 2, l = 2)),
        strip.background = element_rect(fill="lightblue", colour="black",size=0.6),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.3, linetype = "solid",
                                 colour = "black"))
# dev2bitmap(file = "manuscript/figures/S2_FIG1a_Call.pdf",
#            height = 4, width = 4, units = 'in', res=600, method = "pdf")
dev.off()

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
##############################
##############################
######### S2 Fig S1b #########
##############################
pref_out = "gaus_skew/BCDG/gaus_skew_hier_"
CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
head(CZ_matrix)
gmod = readRDS(paste0("models/", pref_out, basename(dirname(pref_out)),".rds"))
gmod@model_pars
hyp_nm = c("b", "c", "d", "alpha")
# h=1
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
isl = 4
hypp = c("BCDG_Bhyp", "BCDG_Chyp", "BCDG_Dhyp", "BCDG_ALPHAhyp")
gaus_skew[[isl]]
eff_var = list(b1 = c("b_intercept", "b_shoreCZB", "b_shoreCZC", "b_test_sexmale"),
               c = c("c_intercept", "c_shoreCZB", "c_shoreCZC", "c_ref_ecotypewave"),
               d = c("d_intercept"),
               alpha = c("alpha_intercept"))
# gaus_skew[[isl]]$fig = "c-hyperparameter\nIsland and ecotype effect"
# alpha_pars = gaus_skew[[4]]$parameter
# gaus_skew[[4]]$parameter = gaus_skew[[2]]$parameter
rownames(gaus_skew[[isl]]) = as.character(gaus_skew[[isl]][, 1])
gaus_skew[[isl]][, 1] = NULL
if (isl == 1) {
  gaus_skew[[isl]]$fig = "Full hierarchical - hyperparameter b[1]\nIsland and sex effects"
  gaus_skew[[isl]]$eff = 0
  gaus_skew[[isl]][grepl(pattern = "b_intercept|b_shoreCZB|b_shoreCZC|b_test_sexmale$",
                         x = rownames(gaus_skew[[isl]])), "eff"] = 1
} else if (isl == 2) {
  gaus_skew[[isl]]$fig = "Full hierarchical - hyperparameter c\nIsland and ecotype effects"
  gaus_skew[[isl]]$eff = 0
  gaus_skew[[isl]][grepl(pattern = "c_intercept|c_shoreCZB|c_shoreCZC|c_ref_ecotypewave",
                         x = rownames(gaus_skew[[isl]])), "eff"] = 1
} else {
  gaus_skew[[isl]]$fig = paste0("Full hierarchical - hyperparameter ", hyp_nm[isl], "\nNo additional effect")
  gaus_skew[[isl]]$eff = 0
  gaus_skew[[isl]][grepl(pattern = "d_intercept|alpha_intercept",
                         x = rownames(gaus_skew[[isl]])), "eff"] = 1
}
gaus_skew[[isl]]
gaus_skew[[isl]]$parameter = c(paste0(hyp_nm[isl],"0"), paste0(hyp_nm[isl],"1"), paste0(hyp_nm[isl],"2"),
                               paste0(hyp_nm[isl],"3"), paste0(hyp_nm[isl],"4"), paste0(hyp_nm[isl],"6"),
                               paste0(hyp_nm[isl],"5T"), paste0(hyp_nm[isl],"6:",hyp_nm[isl],"5T"))
gaus_skew[[isl]]$parameter = factor(gaus_skew[[isl]]$parameter)
gaus_skew[[isl]]$parameter = factor(gaus_skew[[isl]]$parameter, levels = rev(levels(gaus_skew[[isl]]$parameter)))
# gaus_skew[[isl]][8, "eff"] = 0
# gaus_skew[[isl]]$parameter = rownames(gaus_skew[[isl]])
gaus_skew[[isl]] = remove_rownames(gaus_skew[[isl]])
tiff(filename = paste0("manuscript/figures/S1_FIG1b_",hypp[isl], ".tiff"), width = 2, height = 2, units = "in",
     res = 600)
# pdf(file = "manuscript/figures/S2_FIG1a_Call.pdf", width = 4, height = 4)
ggplot(data = gaus_skew[[isl]]) +
  facet_wrap(~fig) +
  # geom_segment(aes(x=-Inf, xend=Inf, y=10.5, yend=10.5), linetype = "dashed", size = 0.3, col = "grey80") +
  geom_segment(aes(x=0, xend=0, y=-Inf, yend=7.5), linetype = "dashed", size = 0.3, col = "grey80") +
  geom_errorbarh(aes(xmax = `2.5%`, xmin = `97.5%`, y = parameter, height = 0), size = 0.4, col="grey60") +
  geom_errorbarh(aes(xmax = `25%`, xmin = `75%`, y = parameter, height = 0), size = 0.5, col="grey40") +
  geom_point(aes(x = mean, y = parameter, col = factor(eff)), size = 0.2) +
  # geom_point(aes(x = mean, y = parameter, col = "idot"), size = 1.8) +
  scale_colour_manual(values = c("black", "red")) +
  # scale_colour_brewer(type = "qual", palette = 8, direction = 1) +
  # scale_color_viridis_d(begin = 0, end = 1, option = "A", direction = -1) +
  # scale_y_continuous(breaks = seq(1,11,1)) +
  xlim(-1,2) +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_text(size = 3),
        axis.ticks = element_line(size = 0.3),
        strip.text = element_text(size=5, margin = margin(t = 2, r = 2, b = 2, l = 2)),
        strip.background = element_rect(fill="lightblue", colour="black",size=0.6),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.3, linetype = "solid",
                                 colour = "black"))
# dev2bitmap(file = "manuscript/figures/S2_FIG1a_Call.pdf",
#            height = 4, width = 4, units = 'in', res=600, method = "pdf")
dev.off()

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
##########################
##########################
##### S2 Fig S2-left #####
##########################
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
##########################
##########################
##### S2 Fig S2-left #####
##########################
pref_out = "gaus_skew/Chyp_isleco/gaus_skew_hier_"
pref_out = "gaus_skew/Bhyp_islsex/gaus_skew_hier_"
gaus_skew = readRDS(paste0("models/", pref_out, basename(dirname(pref_out)),".rds"))
CZ_data$stan_yhat = summary(gaus_skew, pars = c("y_hat"))$summary[,'mean']
CZ_data$stan_yhat_uci = summary(gaus_skew, pars = c("y_hat"))$summary[,'97.5%']
CZ_data$stan_yhat_lci = summary(gaus_skew, pars = c("y_hat"))$summary[,'2.5%']
# pdf(paste0("figures/", pref_out, basename(dirname(pref_out)),"_fit.pdf"), width=7, height=7)
tiff(filename = paste0("manuscript/figures/S2_FIGS2_left_", basename(dirname(pref_out)), ".tiff"),
     width = 2, height = 2, units = "in", res = 600)
# tiff(filename = paste0("manuscript/figures/S2_FIGS2_right_", basename(dirname(pref_out)), ".tiff"),
#      width = 2, height = 2, units = "in", res = 600)
ggplot(data = CZ_data, aes(x = size_ratio, y = stan_yhat)) +
  facet_wrap(~shore) +
  geom_line(aes(col=test_sex), size=0.5) +
  # geom_line(aes(col=ref_ecotype), size=0.5) +
  # geom_line(aes(x = size_ratio, y = stan_yhat_lci), alpha=0.5, linetype = "dashed") +
  # geom_line(aes(x = size_ratio, y = stan_yhat_uci), alpha=0.5, linetype = "dashed")
  geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = test_sex), alpha=0.25) +
  # geom_ribbon(aes(x = size_ratio, ymin = stan_yhat_lci, ymax = stan_yhat_uci, fill = ref_ecotype), alpha=0.25) +
  scale_color_manual(values = c("red", "blue")) +
  # scale_color_manual(values = c("darkgoldenrod3", "darkorchid4")) +
  labs(x="ln female size - ln male size", y="probability of mating", col="", fill ="") +
  scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.position = 'top', legend.key.size = unit(0.4,"line"),
        legend.justification="left", legend.margin=margin(6,-1,-4,8),
        legend.box.margin=margin(-10,-10,-10,-10), legend.spacing.x = unit(0.05, 'cm'),
        strip.text = element_text(size=4, margin = margin(t = 2, r = 2, b = 2, l = 2)),
        strip.background = element_rect(fill="lightblue", colour="black",size=0.6),
        legend.text = element_text(size = 3.5),
        axis.text = element_text(size = 3),
        axis.title = element_text(size = 4.5),
        axis.ticks = element_line(size = 0.3),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.3, linetype = "solid",
                                 colour = "black"))
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
#################################
#################################
######### Table dimorph #########
#################################
# CZ_am_ss_fig is on huluvu server
lapply(seq_along(islands), function(s) {
  cat("Saving", paste0("tables/", pref_out, "summary/", islands[s], "_cline_dimorph.csv"), "...\n")
  write.csv(x = CZ_am_ss_fig[[s]], file = paste0("tables/", pref_out, "summary/", islands[s], "_cline_dimorph.csv"),
            row.names = FALSE)
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

CZBdat = CZ_data[CZ_data$shore=='CZB', ]
range(CZBdat$LCmeanDist)
CZBdat50 = CZBdat[CZBdat$LCmeanDist < 53, ]
colnames(CZBdat50)
pdf("../../conf/ESEB/Turku2019/CZB_obs_fm_size_50.pdf", width=8, height=7)
ggplot(data = CZBdat50) +
  geom_histogram(aes(x = log(length_mm), y = (..count..)/sum(..count..), fill = test_sex),
                 col='black', position = 'dodge', binwidth = 0.2) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x='ln size', y='frequency', fill='') +
  theme(legend.position = 'top', legend.text = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 18),
        axis.ticks = element_line(size = 2),
        axis.text = element_text(size = 15))
dev.off()
log(mean(CZBdat50[CZBdat50$test_sex=='female', ]$length_mm))
log(mean(CZBdat50[CZBdat50$test_sex=='male', ]$length_mm))

# rbindlist_fread <- function(path, pattern = "*.csv") {
#   files = list.files(path, pattern, full.names = TRUE)
#   rbindlist(lapply(files, function(x) fread(x)))
# }

CZB6 = read.csv("tables/gaus_skew/SKEW/repsims/CZB_6_sim_YN.csv")
head(CZB6)
CZB6_1 = CZB6[CZB6$mountYN==1, ]
CZB6_CI = as.data.frame(rbind(CI(CZB6_1$male), CI(CZB6$male)))
CZB6_CI$fig = c('mated', 'all')
CZB6_1$fig = 'mated'
CZB6$fig = 'all'
CZB6 = rbind(CZB6, CZB6_1)
tail(CZB6)
pdf("../../conf/ESEB/Turku2019/CZB_sim_pos6_SS.pdf", width=8, height=7)
ggplot(data = CZB6) +
  geom_histogram(aes(x = male, y = (..count..)/sum(..count..), fill = fig),
                 col='black', binwidth = 0.15, position = 'dodge') +
  geom_vline(xintercept = CZB6_CI$mean[1], col="darkcyan", size=1.7) +
  geom_vline(xintercept = CZB6_CI$mean[2], col="blue", size=1.7) +
  # scale_fill_viridis_d() +
  # geom_histogram(data = CZB6[CZB6$mountYN==1, ], aes(x = male, y = (..count..)/sum(..count..), fill = 'mated'),
  #                col='black', binwidth = 0.2, position = 'dodge') +
  scale_fill_manual(values = c("blue", "darkcyan")) +
  labs(x='ln male size', y='frequency', fill='') +
  theme(legend.position = 'top', legend.text = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 18),
        axis.ticks = element_line(size = 2),
        axis.text = element_text(size = 15))
dev.off()

with(data = CZB6_1, cor(female, male))
pdf("../../conf/ESEB/Turku2019/figures/CZB_sim_pos6_AM.pdf", width=8, height=7)
ggplot(data = CZB6_1) +
  geom_abline(slope=1, linetype='dashed', size=2, alpha=0.6) +
  # geom_point(aes(x = male, y = female), size=3, col='black') +
  geom_point(aes(x = male, y = female), size=2, col='black') +
  labs(x='ln male size', y='ln female size') +
  theme(axis.title = element_text(face = "bold", size = 18),
        axis.ticks = element_line(size = 2),
        axis.text = element_text(size = 15)) +
  annotate("text", x = 1, y = 2.7, label = "paste(\"Pearson's \", ,italic(r), \" = 0.33\")",
           parse = TRUE, size = 7)
dev.off()
############################
############################
############################