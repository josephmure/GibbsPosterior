#library(ggplot2)
#library(pscl)

library(Rcpp)
library(RcppGSL)

sourceCpp("../../RcppGibbsPosterior.cpp")

# NOMBRE_DIMENSIONS <- ncol(x_connus)
# REGULARITE <- scan("../regularite_vraie.txt")
# TYPE_NOYAU_MATERN_TYPE_PRIOR <- scan("../type_noyau_prior.txt", what="character")
# TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
# TYPE_PRIOR <- TYPE_NOYAU_MATERN_TYPE_PRIOR[2]
# NOMBRE_PAS_METROPOLIS <- 100
# ECART_TYPE_METROPOLIS <- 0.4
# 
# GERME_ALEATOIRE_C <- germe_aleatoire


