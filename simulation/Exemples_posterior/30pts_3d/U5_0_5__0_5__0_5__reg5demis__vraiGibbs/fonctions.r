##Définit la liste des fonctions de base intervenant en krigeage universel.

FONCTIONS <- NULL

FONCTIONS <- c(FONCTIONS, function(x) 1)

##Définit la liste des fonctions de base dont la moyenne est reellement combinaison lineaire

FONCTIONS_REELLES <- NULL

FONCTIONS_REELLES <- c(FONCTIONS_REELLES, function(x) 1)
