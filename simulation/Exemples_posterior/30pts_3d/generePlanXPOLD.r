library(lhs)
library(pscl)

source("../../../../Briques_these.r")

GERME_ALEATOIRE_R <- germe_aleatoire
NOMBRE_DIMENSIONS <- length( read.table("../longueur_correlation_vraie.txt",as.is=TRUE)$V1 ) ## Dimension = nombre de longueurs de correlation renseignees dans ../longueur_correlation_vraie.txt

NOMBRE_POINTS_PLANXP <- 30
LHS <- FALSE

write.matrix(NOMBRE_DIMENSIONS,"../nombre_dimensions.txt")

set.seed(GERME_ALEATOIRE_R)

x_connus <- generePointsConnus(nombre_points = NOMBRE_POINTS_PLANXP, minima = rep(0,NOMBRE_DIMENSIONS), maxima = rep(1,NOMBRE_DIMENSIONS), lhs=LHS)
colnames(x_connus) <- NULL

#x_connus <- as.matrix(read.table("../../../x_connus20.txt", as.is = TRUE))

write.matrix(x_connus,"planXP.txt",sep="\t")

