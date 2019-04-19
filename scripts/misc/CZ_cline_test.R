library(bbmle)
install.packages('rdrop2')
library(rdrop2)
drop_auth()
drop_acc() %>% data.frame()
drop_dir("Samuel's studentship/CZ-mate-choice")
#drop_download("Samuel's studentship/CZ-mate-choice/data/CZ_all_mating_clean.csv")

CZ_all = drop_read_csv("Samuel's studentship/CZ-mate-choice/data/CZ_all_mating_clean.csv", sep=";")

cline_2c3s <- function(phen,position,sex,cl,cr,wl,wr,crab,wave,zs_c,zs_w,sc,sh,sw){
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # DECREASING
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]
  s_xl <- sqrt(sc^2 + 4*p_xl*(1-p_xl)*sh^2 + (p_xl^2)*(sw^2-sc^2))
  
  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # INCREASING 
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  # combined cline
  z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  s_x[z_x < z_xl] <- s_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  return(minusll)
}

# estimate cline parameters for the four CZs
mle.cline.2c3s = list(CZA=NULL, CZB=NULL, CZC=NULL, CZD=NULL)

for (p in levels(CZ_all$shore)) {
  if (p=='CZA'){
    plot(CZ_all$DistAlongPath[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=170,cr=280,wl=20,wr=10,crab=-2.1,wave=-1.9,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.2,sw=0.2)
    mle.cline.2c3s$CZA = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$DistAlongPath[CZ_all$shore==p],
                                        sex=CZ_all$sex[CZ_all$shore==p]))
  }
  else if (p=='CZB'){
    plot(CZ_all$DistAlongPath[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=70,cr=125,wl=5,wr=50,crab=-2.5,wave=-1.5,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.2,sw=0.2)
    mle.cline.2c3s$CZB = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$DistAlongPath[CZ_all$shore==p],
                                        sex=CZ_all$sex[CZ_all$shore==p]))
  }
  else if (p=='CZC'){
    plot(CZ_all$DistAlongPath[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=50,cr=125,wl=10,wr=20,crab=-2.5,wave=-1.5,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.2,sw=0.2)
    mle.cline.2c3s$CZC = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$DistAlongPath[CZ_all$shore==p],
                                        sex=CZ_all$sex[CZ_all$shore==p]))
  }
  else {
    plot(CZ_all$DistAlongPath[CZ_all$shore==p], log(CZ_all$length_mm[CZ_all$shore==p]))
    title(main = p)
    theta.init = list(cl=80,cr=165,wl=5,wr=10,crab=-2.5,wave=-1.5,zs_c=-0.1,zs_w=-0.1,sc=0.2,sh=0.2,sw=0.2)
    mle.cline.2c3s$CZD = mle2(cline_2c3s, theta.init,
                              control=list(parscale=abs(unlist(theta.init))),
                              data=list(phen=-log(CZ_all$length_mm[CZ_all$shore==p]),
                                        position=CZ_all$DistAlongPath[CZ_all$shore==p],
                                        sex=CZ_all$sex[CZ_all$shore==p]))
  }
}


# hybrid variances (sh) are odd
(CZ_cline_params = sapply(mle.cline.2c3s, function(x) round(coef(x), 2)))
(CZ_cline_se = sapply(mle.cline.2c3s, function(x) round(sqrt(diag(vcov(x))), 2)))
sapply(mle.cline.2c3s, function(x) summary(x))
sapply(mle.cline.2c3s, function(x) AIC(x))
