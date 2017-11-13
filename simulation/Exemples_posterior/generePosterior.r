#library(distr)
library(lhs)
library(pscl)
library(ggplot2)
#library(rgl)

load("../../../../SamplingPosteriorSur.RData")
set.seed(1)

x_connus <- as.matrix(read.table("../../x_connus20.txt", as.is = TRUE))

REGULARITE <- 1.5
LONGUEUR_CORRELATION <- c(0.3,0.6)

NOYAU <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE )
y_connus <- creeEchantillonNormal(x_connus = x_connus, noyau = NOYAU)

write.matrix(y_connus,'y20.txt',sep = "\t")



PosteriorComplet <- Gibbs2DPosterior(x_connus = x_connus,y_connus = y_connus, longueur_correlation_init = c(0.5,0.5),regularite=REGULARITE , nb_iterations = 1000)
Posterior <- PosteriorComplet[,1:2]

write.matrix(PosteriorComplet, "PosteriorComplet.txt", sep= "\t")
write.matrix(Posterior, "Posterior.txt", sep= "\t")

plot(Posterior)
dev.copy(pdf,"TensReg1.5LC0.30.6nuage.pdf")
dev.off()

PosteriorDF  <- data.frame(Posterior )

ggplot(PosteriorDF ,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
dev.copy(pdf,"TensReg1.5LC0.30.6contour.pdf")
dev.off()

ggplot(PosteriorDF ,aes(x=X1,y=X2))+
    stat_density2d(aes(fill=..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="yellow")
dev.copy(pdf,"TensReg1.5LC0.30.6couleur.pdf")
dev.off()
