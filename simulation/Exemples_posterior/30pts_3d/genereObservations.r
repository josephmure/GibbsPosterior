#library(distr)
library(pscl)
#library(ggplot2)
#library(rgl)

source("../../../../Matern_these.r")
source("../../../../Briques_these.r")
source("../../../../Krige_these.r")


REGULARITE <- read.table("../regularite_vraie.txt",as.is=TRUE)$V1
LONGUEUR_CORRELATION <- read.table("../longueur_correlation_vraie.txt",as.is=TRUE)$V1
NOMBRE_FONCTIONS_REELLES_TENDANCE <- length(FONCTIONS_REELLES)
tendance_reelle <- matrix(NA,nrow=NOMBRE_POINTS_PLANXP,ncol=NOMBRE_FONCTIONS_REELLES_TENDANCE)
MOYENNE <- 0
    if(NOMBRE_FONCTIONS_REELLES_TENDANCE>0)
    {        
        for(i in 1:NOMBRE_FONCTIONS_REELLES_TENDANCE)
        {
	    tendance_reelle[,i] <- apply(x_connus,1,FONCTIONS_REELLES[[i]])
        }
	MOYENNE <- c(tendance_reelle %*% BETA)
    }

TYPE_NOYAU_MATERN_TYPE_PRIOR <- c(as.matrix(read.table("../type_noyau_prior.txt")))
TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
TYPE_PRIOR <- TYPE_NOYAU_MATERN_TYPE_PRIOR[2]

GERME_ALEATOIRE_R <- germe_aleatoire


set.seed(GERME_ALEATOIRE_R) 

x_connus <- as.matrix(read.table("planXP.txt", as.is = TRUE))


if(TYPE_NOYAU_MATERN == "geometrique") 
{
	NOYAU <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE ) # anisotrope geometrique en fait !
} else if (TYPE_NOYAU_MATERN == "tensorise") 
{
	NOYAU <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE )
} else 
{
	print("Erreur dans genereObservations.r : le noyau doit-il etre anis. geometrique ou tensorise ?")
	NOYAU <- NULL
}

if(!is.null(NOYAU))
{
	write.table(c(TYPE_NOYAU_MATERN,TYPE_PRIOR),"type_noyau.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
	y_connus <- creeEchantillonNormal(x_connus = x_connus, noyau = NOYAU) + MOYENNE
	write.matrix(y_connus,"observations.txt",sep = "\t")
}


