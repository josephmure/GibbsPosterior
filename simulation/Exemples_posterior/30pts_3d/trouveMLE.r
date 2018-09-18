library(pscl)
#library(GenSA)


dossier_courant <- getwd()
setwd("../../../..")
source("MaxVraisemblanceIntegree.r")
setwd(dossier_courant)


#dossier_courant <- getwd()
#setwd("..")
#source("MaxVraisemblanceIntegree.r")
#source("densKernel.r")
#setwd(dossier_courant)

# REGULARITE <- scan("regularite_vraie.txt")
# TYPE_NOYAU_MATERN_TYPE_PRIOR <- scan("type_noyau_prior.txt", what="character")
# TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
# REGULARITE <- read.table("regularite_vraie.txt",as.is=TRUE)$V1
# 
# 
# x_connus <- as.matrix(read.table("planXP.txt", as.is = TRUE))
# y_connus <- scan("observations.txt")
# tendance <- as.matrix(read.table("tendance.txt",as.is=TRUE))
# 
# NOMBRE_DIMENSIONS <- ncol(x_connus)

trouveMLE <- function(NOMBRE_DIMENSIONS,TYPE_NOYAU_MATERN, REGULARITE, x_connus, y_connus, tendance)
{

if(TYPE_NOYAU_MATERN == "geometrique")
{
	opt1 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeAnisGeom, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS))
	opt2 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeAnisGeom, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS))
	opt3 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeAnisGeom, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS))
	opt4 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeAnisGeom, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS))
	opt5 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeAnisGeom, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS))
		}

if(TYPE_NOYAU_MATERN == "tensorise")
{
	print("Dans le fonction d'optimisation trouveMLE : cas tensorise")
	opt1 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeTens, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
	opt2 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeTens, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
	opt3 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeTens, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
	opt4 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeTens, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
	opt5 <- optim(par=runif(NOMBRE_DIMENSIONS,max=2),fn = evalueOpposeLogVraisemblanceIntegreeTens, x_connus= x_connus, y_connus= y_connus, regularite=REGULARITE, tendance=tendance, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 
	  }

#if(opt$convergence !=0) print(paste("Attention ! Pas de convergence de optim pour estimer le MLE. Code :", opt$convergence))

opt_value <- c(opt1$value,opt2$value,opt3$value,opt4$value,opt5$value)  
opt_par <- rbind(opt1$par,opt2$par,opt3$par,opt4$par,opt5$par)  


argmax_vraisemblance_integree <- opt_par[which.min(opt_value),]

argmax_vraisemblance_integree
}
#write.matrix(argmax_vraisemblance_integree, "argmax_vraisemblance_integree.txt")

