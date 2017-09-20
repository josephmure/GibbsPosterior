library(lhs)
library(pscl)


source("../../../../Briques_these.r")


generePlanXP <- function(GERME_ALEATOIRE_R, NOMBRE_POINTS_PLANXP, NOMBRE_DIMENSIONS, LHS)
{
  set.seed(GERME_ALEATOIRE_R)
  
  x_connus <- generePointsConnus(nombre_points = NOMBRE_POINTS_PLANXP, minima = rep(0,NOMBRE_DIMENSIONS), maxima = rep(1,NOMBRE_DIMENSIONS), lhs=LHS)
  
  x_connus
}
