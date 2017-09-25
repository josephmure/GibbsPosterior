##Définit la liste des fonctions de base intervenant en krigeage universel.

FONCTIONS <- NULL

FONCTIONS <- c(FONCTIONS, function(x) 1)

FONCTIONS <- c(FONCTIONS, function(x) x[1])

FONCTIONS <- c(FONCTIONS, function(x) x[2])

FONCTIONS <- c(FONCTIONS, function(x) x[3])

FONCTIONS <- c(FONCTIONS, function(x) x[4])

FONCTIONS <- c(FONCTIONS, function(x) x[5])

FONCTIONS <- c(FONCTIONS, function(x) x[6])

FONCTIONS <- c(FONCTIONS, function(x) x[7])


##Définit les fonctions formant une base de l'espace auquel appartient REELLEMENT la fonction moyenne

FONCTIONS_REELLES <- NULL


