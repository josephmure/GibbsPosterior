#library(distr)
library(pscl)
#library(ggplot2)
#library(rgl)

##source("../Matern_these.r")
##source("../Briques_these.r")
##source("../Krige_these.r")

LONGUEUR_CORRELATION <- c(0.4,0.8,0.2)



creeMoyenne <- function(x_connus,NOMBRE_POINTS_PLANXP, FONCTIONS_REELLES, BETA)
{
if(length(FONCTIONS_REELLES)==length(BETA)) {
tendance_reelle <- matrix(NA,nrow=NOMBRE_POINTS_PLANXP,ncol=length(BETA))
MOYENNE <- 0
    if(length(BETA)>0)
    {        
        for(i in 1:length(BETA))
        {
	    tendance_reelle[,i] <- apply(x_connus,1,FONCTIONS_REELLES[[i]])
        }
	MOYENNE <- c(tendance_reelle %*% BETA)
    }
} else {
print(paste0("Longueur de BETA = ",length(BETA)," alors que Longueur de FONCTIONS_REELLES = ",length(FONCTIONS_REELLES)))
}

MOYENNE
}


creeObservations <- function(x_connus,TYPE_NOYAU_MATERN, LONGUEUR_CORRELATION, REGULARITE, MOYENNE)
{
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
	#write.table(c(TYPE_NOYAU_MATERN,TYPE_PRIOR),"type_noyau.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
	y_connus <- creeEchantillonNormal(x_connus = x_connus, noyau = NOYAU) + MOYENNE
	#write.matrix(y_connus,"observations.txt",sep = "\t")
}
y_connus
}



