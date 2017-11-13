library(pscl)
NOMBRE_PROCESSUS <- 100

VECTEUR_GERMES_ALEATOIRES <- NULL

MLE <- NULL
MAP <- NULL
#MAP_optim_classique <- NULL

MATRICE_DISTANCES <- NULL
MATRICE_DISTANCES_FISHER <- NULL

MATRICE_ERREURS_PREDICTION <- NULL

## NOUVEAU : Suppression des fichiers au début
system("rm Posterior_anis_geom.txt")
system("rm Posterior_tens.txt")
## FIN NOUVEAU

for(germe_aleatoire in (GERME_ALEATOIRE_GLOBAL+1):(GERME_ALEATOIRE_GLOBAL+NOMBRE_PROCESSUS) )
{
	print(paste("Germe aleatoire", germe_aleatoire) )
	VECTEUR_GERMES_ALEATOIRES <- c(VECTEUR_GERMES_ALEATOIRES, germe_aleatoire)
	write.matrix(VECTEUR_GERMES_ALEATOIRES, "liste_germes_aleatoires.txt")
	source("../../nettoyage.r")
	source("../../generePlanXP.r")
	source("../../genereMatriceTendance.r")
	source("../../genereObservations.r")
	source("../../scriptR.r")
	source("../../postTraitement.r")


## INUTILE PUISQUE CES COMPARAISONS SONT FAITES APRES COUP
##	source("../../compareResultats.r")

##	MATRICE_DISTANCES <- rbind(MATRICE_DISTANCES, c( Distance_argmax_vraisemblance[1], Distance_mode_posterior[1]) )
##	write.matrix(MATRICE_DISTANCES,"matrice_distances.txt",sep = "\t")#

##	MATRICE_DISTANCES_FISHER <- rbind(MATRICE_DISTANCES_FISHER, c( Distance_argmax_vraisemblance[2], Distance_mode_posterior[2]) )
##	write.matrix(MATRICE_DISTANCES_FISHER,"matrice_distances_fisher.txt",sep = "\t")

##	MATRICE_ERREURS_PREDICTION <- rbind(MATRICE_ERREURS_PREDICTION, c(Ecart_prediction_argmax_vraisemblance, Ecart_prediction_mode_posterior) )
##	write.matrix(MATRICE_ERREURS_PREDICTION,"matrice_ecarts_predictions.txt",sep = "\t")
## FIN INUTILE


	MLE <- rbind(MLE, argmax_vraisemblance_integree)
	write.matrix(MLE,"liste_MLE.txt",sep = "\t")
	MAP <- rbind(MAP, mode_posterior)
	write.matrix(MAP,"liste_MAP.txt",sep = "\t")

##	MAP_optim_classique <- rbind(MAP_optim_classique, mode_posterior_optim_classique)
##	write.matrix(MAP_optim_classique,"liste_MAP_optim_classique.txt", sep = "\t")
}
