library(lhs)
library(pscl)

source("../../../../Briques_these.r")

NOMBRE_DIMENSIONS <- length( read.table("../longueur_correlation_vraie.txt",as.is=TRUE)$V1 ) ## Dimension = nombre de longueurs de correlation renseignees dans ../longueur_correlation_vraie.txt

#Cree FONCTIONS, la liste des fonctions servant a former les colonnes de tendance (H) ainsi que FONCTIONS_REELLES, la liste des fonctions dont la moyenne est REELLEMENT une combinaison lineaire.
source("../fonctions.r")

if(length(FONCTIONS_REELLES)==0) print("Attention ! FONCTIONS_REELLES ne devrait pas etre de longueur nulle : il faudrait au moins inclure la fonction constante_1 dans les fonctions dont la vraie fonction de moyenne est combinaison lineaire, quitte a l'annuler par BETA=0.")

if(length(FONCTIONS)<1)
{
	write.table("vide","tendance.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
} else
{
	tendance <- matrix(0,nrow=NOMBRE_POINTS_PLANXP, ncol=length(FONCTIONS))
	for(i in 1:length(FONCTIONS))
	{
		tendance[,i] <- apply(x_connus,1,FONCTIONS[[i]])
	}
	write.matrix(tendance,"tendance.txt",sep="\t")
}


BETA <- scan("../beta_vrai.txt") #le vecteur des coefficients lineaires reels des fonctions de FONCTIONS_REELLES