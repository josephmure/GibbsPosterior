library(pscl)
#library(GenSA)


dossier_courant <- getwd()
setwd("../../../..")
source("densKernel.r")
setwd(dossier_courant)

# dossier_courant <- getwd()
# setwd("..")
# source("MaxVraisemblanceIntegree.r")
# source("densKernel.r")
# setwd(dossier_courant)
# 
# NOMBRE_DIMENSIONS <- ncol(x_connus)
# REGULARITE <- scan("regularite_vraie.txt")
# TYPE_NOYAU_MATERN_TYPE_PRIOR <- scan("type_noyau_prior.txt", what="character")
# TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
# REGULARITE <- read.table("regularite_vraie.txt",as.is=TRUE)$V1
# 
# TEMPS_RECUIT <- 10
# DEPART_opt_MAP <- rep(2, NOMBRE_DIMENSIONS)
# FENETRE_opt_MAP <- 1

#if(TYPE_NOYAU_MATERN == "geometrique")
#{
#	Posterior<-as.matrix(read.table("Posterior_anis_geom.txt", as.is=TRUE))
#}

#if (TYPE_NOYAU_MATERN == "tensorise")
#{
#	Posterior<-as.matrix(read.table("Posterior_tens.txt", as.is=TRUE))
#}



#recuit <- GenSA(par=rep(0.5,NOMBRE_DIMENSIONS), fn= densKernel, fenetre= 0.2, matrice_echantillon=Posterior, lower=rep(0,NOMBRE_DIMENSIONS), upper= rep(max(Posterior),NOMBRE_DIMENSIONS), control = list(max.time=TEMPS_RECUIT))

#mode_posterior <- recuit$par

#write.matrix(mode_posterior, "mode_posterior.txt")


trouveMAP <- function(FENETRE_opt_MAP, Posterior, NOMBRE_DIMENSIONS)
{
opt1 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = densKernel, fenetre= FENETRE_opt_MAP, matrice_echantillon=Posterior, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
opt2 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = densKernel, fenetre= FENETRE_opt_MAP, matrice_echantillon=Posterior, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
opt3 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = densKernel, fenetre= FENETRE_opt_MAP, matrice_echantillon=Posterior, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
opt4 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = densKernel, fenetre= FENETRE_opt_MAP, matrice_echantillon=Posterior, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
opt5 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = densKernel, fenetre= FENETRE_opt_MAP, matrice_echantillon=Posterior, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 


opt_value <- c(opt1$value,opt2$value,opt3$value,opt4$value,opt5$value)  
opt_par <- rbind(opt1$par,opt2$par,opt3$par,opt4$par,opt5$par) 

#if(opt_MAP$convergence !=0) print(paste("Attention ! Pas de convergence de optim pour estimer le MAP. Code :", opt_MAP$convergence))


#mode_posterior_optim_classique <- opt_MAP$par
mode_posterior <- opt_par[which.min(opt_value),]
}


#write.matrix(mode_posterior_optim_classique, "mode_posterior_optim_classique.txt")
#write.matrix(mode_posterior, "mode_posterior.txt")

#write.table(t(mode_posterior), append=TRUE, row.names=FALSE, col.names=FALSE, "mode_posterior.txt", sep= "\t") #mode_posterior est un vecteur. La transposition fait en sorte que ce soit une ligne qui soit ajoutee au fichier mode_posterior.txt
