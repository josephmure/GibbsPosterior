#library(distr)
library(pscl)
#library(ggplot2)
#library(rgl)

##source("../Matern_these.r")
##source("../Briques_these.r")
##source("../Krige_these.r")
source("../Rastrigin.r")

#x_connus <- as.matrix(read.table("planXP.txt", as.is = TRUE))
source("../../generePlanXP.r") ## genere le plan d'expériences

y_connus <- Rastrigin(x_connus)
#write.matrix(y_connus,"observations.txt",sep = "\t")

