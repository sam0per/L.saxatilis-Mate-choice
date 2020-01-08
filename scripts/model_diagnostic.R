rm(list = ls())

.packages = c("dplyr", "tibble", "boot", "Rmisc", "purrr", "reshape2", "parallel", "optparse", "pracma", "broom",
              "bayesplot", "ggplot2", "rstan", "hexbin", "data.table", "tools", "corrplot")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(.packages, require, character.only=TRUE)


# focus = readRDS("models/gaus_skew/Bhyp_shore/gaus_skew_hier_Bhyp_shore.rds")
focus = readRDS("models/gaus_skew/SKEW/gaus_skew.rds")
posterior_f <- as.array(focus)
par_n = dimnames(posterior_f)[['parameters']]
par_f = par_n[!grepl(pattern = "hyp|hat|lik|lp_", x = par_n)]
par_l = length(par_f)
par_comb = combn(1:par_l, m = 2)

mod_n = basename(dirname("models/gaus_skew/Bhyp_shore/gaus_skew_hier_Bhyp_shore.rds"))
pred_mx = strsplit(mod_n, split = "_")[[1]][2]
CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
if (pred_mx == "all") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
} else if (pred_mx == "shore") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore, data = CZ_data)[,-1]
} else if (pred_mx == "ecotype") {
  CZ_matrix = model.matrix(mountYNcontact ~ ref_ecotype + shape, data = CZ_data)[,-1]
} else if (pred_mx == "sex") {
  CZ_matrix = matrix(model.matrix(mountYNcontact ~ test_sex, data = CZ_data)[,-1])
  colnames(CZ_matrix) = "test_sexmale"
} else if (pred_mx == "isleco") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore + ref_ecotype + shape, data = CZ_data)[,-1]
} else if (pred_mx == "islsex") {
  CZ_matrix = model.matrix(mountYNcontact ~ shore + test_sex, data = CZ_data)[,-1]
} else if (pred_mx == "ecosex") {
  CZ_matrix = model.matrix(mountYNcontact ~ ref_ecotype + test_sex * shape, data = CZ_data)[,-1]
} else {
  print("Model name is missing")
}
colnames(CZ_matrix)
hyp_n = "b_"
par_f[grepl(pattern = paste0(hyp_n, "coeff"), x = par_f)] = paste0(hyp_n, colnames(CZ_matrix))
par_f
# dimnames(posterior_f)[['parameters']][1:10]
# posterior_f[1:10,2,1]
# plot(posterior_f[,4,1], posterior_f[,4,6])
# cor(posterior_f[,4,1], posterior_f[,4,6])
cor_dt = cor(posterior_f[,1,1:par_l])
rownames(cor_dt)[grepl(pattern = paste0(hyp_n, "coeff"), x = rownames(cor_dt))] = paste0(hyp_n, colnames(CZ_matrix))
colnames(cor_dt)[grepl(pattern = paste0(hyp_n, "coeff"), x = colnames(cor_dt))] = paste0(hyp_n, colnames(CZ_matrix))

dir.create("figures/mod_diagn/")
mod_bs = basename(file_path_sans_ext("models/gaus_skew/Bhyp_shore/gaus_skew_hier_Bhyp_shore.rds"))
svg(filename = paste0("figures/mod_diagn/", mod_bs, "_par_cor.svg"), width = 8, height = 8)
corrplot(cor_dt, method = "number")
dev.off()




# cor_dt = rbindlist(lapply(1:dim(par_comb)[2], function(x) {
#   p_combo = par_comb[, x]
#   par_x = p_combo[1]
#   par_y = p_combo[2]
#   ch1 = data.frame(cor_chain1 = cor(posterior_f[,1,par_x], posterior_f[,1,par_y]))
#   ch2 = data.frame(cor_chain2 = cor(posterior_f[,2,par_x], posterior_f[,2,par_y]))
#   ch3 = data.frame(cor_chain3 = cor(posterior_f[,3,par_x], posterior_f[,3,par_y]))
#   ch4 = data.frame(cor_chain4 = cor(posterior_f[,4,par_x], posterior_f[,4,par_y]))
#   ch0 = paste0(par_f[par_x], "_", par_f[par_y])
#   ch_dt = round(cbind(ch1, ch2, ch3, ch4), 3)
#   ch_dt = cbind(cor_pars = ch0, ch_dt)
#   return(ch_dt)
# }))
# 
# dir.create("tables/mod_diagn/")
# mod_bs = basename(file_path_sans_ext("models/gaus_skew/Bhyp_shore/gaus_skew_hier_Bhyp_shore.rds"))
# write.csv(x = cor_dt, file = paste0("tables/mod_diagn/", mod_bs, "_par_cor.csv"), row.names = FALSE)

