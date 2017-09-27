## Fonction d'Ackley issue du benchmark de Jean-Marc Martinez du GdR Mascot-Num
## http://gdr-mascotnum.math.cnrs.fr/data2/benchmarks/jm.pdf

Ackley_unPoint <- function(vecteur_point)
{
    20 + exp(1) - 20*exp(-0.2*sqrt(mean(vecteur_point*vecteur_point))) - exp(mean(cos(2*pi*vecteur_point)))
}

## La fonction d'Ackley est appliquee a chaque ligne de matrice_points
Ackley <- function(matrice_points)
{
    apply(matrice_points,1,Ackley_unPoint)
}



## Wrapper for Ackley, so that commandes_parallel.r does not have to know that Ackley is being used.
Fonction_emulee <- Ackley

## Nombres de dimensions
NOMBRE_DIMENSIONS <- 7
