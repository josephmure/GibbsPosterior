#library(distr)
library(pscl)
#library(ggplot2)
#library(rgl)

##source("../Matern_these.r")
##source("../Briques_these.r")
##source("../Krige_these.r")
source("../Ackley.r")

#x_connus <- as.matrix(read.table("planXP.txt", as.is = TRUE))
source("../../generePlanXP.r") ## genere le plan d'expériences

y_connus <- Fonction_emulee(x_connus)
#write.matrix(y_connus,"observations.txt",sep = "\t")

