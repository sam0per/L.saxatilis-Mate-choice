start_val = list(list(b_intercept=0.4,
                      c_intercept=-0.17,
                      d_par=0.85,
                      alpha_par=2.32,
                      b_coeff=rep(0,ncol(CZX_matrix)),
                      c_coeff=rep(0,ncol(CZZ_matrix))),
                 list(b_intercept=0.4,
                      c_intercept=-0.17,
                      d_par=0.85,
                      alpha_par=2.32,
                      b_coeff=rep(0,ncol(CZX_matrix)),
                      c_coeff=rep(0,ncol(CZZ_matrix))),
                 list(b_intercept=0.4,
                      c_intercept=-0.17,
                      d_par=0.85,
                      alpha_par=2.32,
                      b_coeff=rep(0,ncol(CZX_matrix)),
                      c_coeff=rep(0,ncol(CZZ_matrix))),
                 list(b_intercept=0.4,
                      c_intercept=-0.17,
                      d_par=0.85,
                      alpha_par=2.32,
                      b_coeff=rep(0,ncol(CZX_matrix)),
                      c_coeff=rep(0,ncol(CZZ_matrix))))
dat = list(N = nrow(CZ_data), y = CZ_data$mountYNcontact, ratio = CZ_data$size_ratio,
           X = CZX_matrix, K = dim(CZX_matrix)[2], Z = CZZ_matrix, J = dim(CZZ_matrix)[2])
skew_hier = rstan::stan(file = "L.saxatilis-Mate-choice/scripts/gaus_skew/gaus_skew_hier_best_sub.stan",
                        data = dat, iter = 8000, warmup = 2000,
                        chains=4, refresh=8000, init = start_val,
                        control = list(adapt_delta = 0.95, max_treedepth = 15))
mod_str = c("models/gaus_skew/BCDG/gaus_skew_hier_BCDG.rds", "models/gaus_skew/BC/gaus_skew_hier_BC.rds")
lapply(seq_along(comp_ls), function(c) {
  cat("Comparing", out_comp_str[[1]], "vs", out_comp_str[[c+1]], "...\n")
  print(comp_ls[[c]], simplify = FALSE)
  write.csv(tibble(x=c("elpd_diff","SE"), comp_ls[[c]]),
            paste0("tables/mods_comp/comp_", out_comp_str[[1]], "_", out_comp_str[[c+1]], ".csv"), row.names=FALSE)
})