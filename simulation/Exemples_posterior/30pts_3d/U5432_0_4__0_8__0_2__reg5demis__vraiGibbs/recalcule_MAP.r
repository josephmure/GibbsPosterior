library(pscl)
library(GenSA)

dossier_courant <- getwd()
setwd("../../..")
source("MaxVraisemblanceIntegree.r")
source("densKernel.r")
setwd(dossier_courant)


NOMBRE_ECHANTILLONS <-500
NOMBRE_DIMENSIONS <- scan("nombre_dimensions.txt")

# Charge la matrice concatenant les echantillons a posteriori

Mat = as.matrix(read.table("matrice_echantillons_posterior.txt",as.is=TRUE))

TAILLE_ECHANTILLON <- nrow(Mat) / NOMBRE_ECHANTILLONS

liste_MAP <- matrix(NA,nrow=NOMBRE_ECHANTILLONS,ncol=NOMBRE_DIMENSIONS)

for(i in 1:NOMBRE_ECHANTILLONS)
{
	Posterior <- Mat[(TAILLE_ECHANTILLON * (i-1) + 1):(TAILLE_ECHANTILLON * i),]

opt_MAP <- optim(par=rep(0.5,NOMBRE_DIMENSIONS),fn = densKernel, fenetre= 0.2, matrice_echantillon=Posterior, method = "L-BFGS-B",lower = rep(0,NOMBRE_DIMENSIONS)) 

if(opt_MAP$convergence !=0) print(paste("Attention ! Pas de convergence de optim pour estimer le MAP. Code :", opt_MAP$convergence))


liste_MAP[i,] <- opt_MAP$par


write.matrix(liste_MAP[1:i,], "MAP.txt")

}
