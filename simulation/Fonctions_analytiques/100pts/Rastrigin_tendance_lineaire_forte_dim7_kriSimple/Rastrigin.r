## Fonction de Rastrigin issue du benchmark de Jean-Marc Martinez du GdR Mascot-Num
## http://gdr-mascotnum.math.cnrs.fr/data2/benchmarks/jm.pdf


Rastrigin_unPoint <- function(vecteur_point)
{
	sum(vecteur_point * vecteur_point - 10 * cos(2*pi*vecteur_point) + 10) + 1000 * sum(vecteur_point) ## ce dernier terme est artificiellement ajoute pour rendre la tendance lineaire pertinente
}

## La fonction d'Ackley est appliquee a chaque ligne de matrice_points
Rastrigin <- function(matrice_points)
{
    apply(matrice_points,1,Rastrigin_unPoint)
}


## Wrapper for Rastrigin, so that commandes_parallel.r does not have to know that Rastrigin is being used.
Fonction_emulee <- Rastrigin

## Nombres de dimensions
NOMBRE_DIMENSIONS <- 7
