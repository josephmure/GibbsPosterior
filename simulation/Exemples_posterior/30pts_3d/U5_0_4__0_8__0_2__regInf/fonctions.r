##Définit la liste des fonctions de base intervenant en krigeage universel.

FONCTIONS <- NULL

FONCTIONS <- c(FONCTIONS,function(x) 1)


##Définit les fonctions formant une base de l'espace auquel appartient REELLEMENT la fonction moyenne

FONCTIONS_REELLES <- NULL

FONCTIONS_REELLES <- c(FONCTIONS_REELLES,function(x) 1)


##True value of the BETA parameter
BETA <- 5