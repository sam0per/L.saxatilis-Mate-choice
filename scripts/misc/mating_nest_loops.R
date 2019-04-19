rm(list = ls())
df <- data.frame(fem=seq(1,4),
                 male=seq(5,8))

sample(seq(1,5), 4, replace = TRUE)

female_c = head(data.frame(female_c))
male_c = head(data.frame(male_c))
female_h = head(data.frame(female_h))
male_h = head(data.frame(male_h))
female_w = head(data.frame(female_w))
male_w = head(data.frame(male_w))



fem=list(crab=female_c, hybrid=female_h, wave=female_w)
mal=list(crab=male_c, hybrid=male_h, wave=male_w)

mat_f = function(female, male) {
  bar = list()
  YN = data.frame()
  for (f in seq_along(female)) {
    success=FALSE
    i=1
    fem = female[f]
    while (!success) {
      m = sample(male, 1, replace = FALSE)
      p = inv.logit(CZ_size_params$mean[CZ_size_params$params=='level'] +
                      CZ_size_params$mean[CZ_size_params$params=='scale'] *
                      exp(-0.5 * (((fem - m) - CZ_size_params$mean[CZ_size_params$params=='preference'])/
                                    CZ_size_params$mean[CZ_size_params$params=='choosiness'])^2) +
                      CZ_size_params$mean[CZ_size_params$params=='asymmetry'] * (fem - m))
      s = rbinom(n = 1, size = 1, prob = p)
      YN[i,'male'] = m
      YN[i,'female'] = fem
      YN[i,'mountYN'] = s
      success = (s > 0)
      i = i + 1
    }
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}


res = lapply(names(fem), function(x) {
  sapply(names(fem$crab), function(y) mat_f(female = fem[[x]][[y]], male = mal[[x]][[y]]))
})
do.call(rbind, res[[2]][,'CZD'])

eco_CZ = lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) do.call(rbind, res[[x]][,y]))
})
eco_CZ[[2]][[4]]

# save the simulated datasets of mate choice
lapply(seq_along(ecotype), function(x) {
  lapply(seq_along(shore), function(y) write.table(eco_CZ[[x]][[y]],
                                                   paste0("misc/", ecotype[x], "_", shore[y], "_sim_YN.csv"),
                                                   row.names = FALSE, col.names = TRUE, sep = ";"))
})





mat_f = function(female, male) {
  bar = list()
  YN = data.frame()
  for (f in seq_along(female)) {
    success=FALSE
    i=1
    femm = female[f]
    while (!success) {
      m = sample(male, 1)
      p = femm + m
      s = rbinom(n = 1, size = 1, prob = p/100)
      YN[i,'male'] = m
      YN[i,'female'] = femm
      YN[i,'mount'] = s
      success = (s > 0)
      i = i + 1
    }
    bar[[f]] = YN
    YN = data.frame()
  }
  return(bar)
}
mat_f(female = fem, male = mal)

for (f in seq_along(df$fem)) {
  success=FALSE
  i=1
  femm = df$fem[f]
  while (!success) {
    m = sample(df$male, 1)
    p = femm + m
    s = rbinom(n = 1, size = 1, prob = p/100)
    YN[i,'male'] = m
    YN[i,'female'] = femm
    YN[i,'mount'] = s
    success = (s > 0)
    i = i + 1
  }
  bar[[f]] = YN
  YN = data.frame()
}
bar
do.call(rbind, bar)

