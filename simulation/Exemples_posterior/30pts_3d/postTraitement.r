library(pscl)
library(GenSA)

dossier_courant <- getwd()
setwd("../../../..")
source("MaxVraisemblanceIntegree.r")
source("densKernel.r")
setwd(dossier_courant)

NOMBRE_DIMENSIONS <- scan("../nombre_dimensions.txt")
TYPE_NOYAU_MATERN_TYPE_PRIOR <- c(as.matrix(read.table("../type_noyau_prior.txt")))
TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
#TYPE_PRIOR <- TYPE_NOYAU_MATERN_TYPE_PRIOR[2] ##Pas besoin ici.
REGULARITE <- read.table("../regularite_vraie.txt",as.is=TRUE)$V1

TEMPS_RECUIT <- 10

####Posterior <- ECHANTILLON_POSTERIOR


## Inutile de charger Posterior, car il est encore en memoire suite a scriptR.r. En outre, ce serait dangereux, car le MAP doit etre estime sur l'echantillon courant, pas la concatenation de tous les echantillons obtenus jusqu'ici !

#if(TYPE_NOYAU_MATERN == "geometrique")
#{
#	Posterior<-as.matrix(read.table("Posterior_anis_geom.txt", as.is=TRUE))
#}

#if (TYPE_NOYAU_MATERN == "tensorise")
#{
#	Posterior<-as.matrix(read.table("Posterior_tens.txt", as.is=TRUE))
#}

x_connus <- as.matrix(read.table("planXP.txt", as.is = TRUE))
y_connus <- scan("observations.txt")
tendance <- as.matrix(read.table("tendance.txt",as.is=TRUE))

opt <- optim(par=rep(0.5,NOMBRE_DIMENSIONS),fn = evalueOpposeLogVraisemblanceIntegreeAnisGeom, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 

if(opt$convergence !=0) print(paste("Attention ! Pas de convergence de optim pour estimer le MLE. Code :", opt$convergence))

argmax_vraisemblance_integree <- opt$par

write.matrix(argmax_vraisemblance_integree, "argmax_vraisemblance_integree.txt")

#recuit <- GenSA(par=rep(0.5,NOMBRE_DIMENSIONS), fn= densKernel, fenetre= 0.2, matrice_echantillon=Posterior, lower=rep(0,NOMBRE_DIMENSIONS), upper= rep(max(Posterior),NOMBRE_DIMENSIONS), control = list(max.time=TEMPS_RECUIT))

#mode_posterior <- recuit$par

#write.matrix(mode_posterior, "mode_posterior.txt")

opt_MAP <- optim(par=rep(0.5,NOMBRE_DIMENSIONS),fn = densKernel, fenetre= 0.2, matrice_echantillon=Posterior, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 

if(opt_MAP$convergence !=0) print(paste("Attention ! Pas de convergence de optim pour estimer le MAP. Code :", opt_MAP$convergence))

#mode_posterior_optim_classique <- opt_MAP$par
mode_posterior <- opt_MAP$par

#write.matrix(mode_posterior_optim_classique, "mode_posterior_optim_classique.txt")
write.matrix(mode_posterior, "mode_posterior.txt")
