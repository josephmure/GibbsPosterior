
# Loi a posteriori

load("PosteriorNettoye.RData")

save.image("PosteriorNettoye.RData")

#load("PosteriorND.RData")

#save.image("PosteriorND.RData")


set.seed(8)

#library(distr)
library(lhs)
library(pscl)
library(ggplot2)
#library(rgl)

x_connus20<- generePointsConnus( nombre_points = 20, minima = c(0,0), maxima = c(1,1) ,lhs= TRUE)

plot(x_connus20,xlab="",ylab="")
dev.copy(pdf,"planXP.pdf")
dev.off()

plot(y_connus20)

NOYAU <- creeMaternTensorise(variance = 1, longueur = c(0.6,0.3,0.7), regularite = 1.5)

y_connus20 <- creeEchantillonNormal(x_connus = x_connus20, noyau = NOYAU)

# Required packages -- POURQUOI ?
library(mvtnorm)
library(ks)
evalueVraisemblanceIntegreeTens <- function( longueur_correlation, x_connus, y_connus, regularite)
    {
    Noyau <- creeMaternTensorise ( variance= 1, longueur= longueur_correlation, regularite= regularite)

    matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)

    Vraisemblance_integree <- det(matriceCorrelation)^(-1/2)* (sum(y_connus* solve(matriceCorrelation,y_connus)))^(-nrow(x_connus)/2)
    Vraisemblance_integree
}


evaluePosteriorCorrelCond <- function( x_connus, y_connus, regularite, longueur_correlation, indice_derive)
    {
    Noyau <- creeMaternTensorise ( variance= 1, longueur= longueur_correlation, regularite= regularite)

    matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)
##    inverseMatriceCorrelation <- solve(matriceCorrelation)
    
    matriceDerivee <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= creeMaternTensoriseDeriveLongueur( variance=1, longueur= Noyau$parametres$longueur, regularite= Noyau$parametres$regularite, indice_derive= indice_derive) )

    MatriceAuxiliaire <- try(solve(matriceCorrelation, matriceDerivee))
    if(class(MatriceAuxiliaire)=="try-error") 
    {
        res <- 0 
        print("INCALCULABLE")
    }
    else 
    {
        Prior_carre <- sum( diag(MatriceAuxiliaire %*% MatriceAuxiliaire) ) - 1/nrow(x_connus) * (sum(diag(MatriceAuxiliaire)))^2
        Prior <- sqrt(Prior_carre)
        Vraisemblance_integree <- det(matriceCorrelation)^(-1/2)* (sum(y_connus* solve(matriceCorrelation,y_connus)))^(-nrow(x_connus)/2)
        res <- Vraisemblance_integree * Prior
        if(is.nan(res)) res <- 0
    }
    res
}

densitePosterior <- function(t,indice_derive, longueur_correlation,x_connus,y_connus,regularite=1.5)
    {
    if(t<=0)res <- 0
    else
    {
        longueur_correlation[indice_derive] <- t

        res <- evaluePosteriorCorrelCond(x_connus = x_connus, y_connus = y_connus, regularite = regularite, longueur_correlation = longueur_correlation,indice_derive = indice_derive)
    }
    res
}

Metropolis <- function ( densite, sd_instrum, init, iterations,...)
{
    res <- NULL
    densinit <- densite(init,...)
    for(i in 1:iterations)
    {
            candidat <- rnorm(n = 1,mean = init, sd = sd_instrum)
            denscandidat <- densite(candidat,...)
            rapport <- denscandidat / densinit
            accepte <- rbinom(n = 1,size = 1, prob = min(1, rapport))
        if(accepte) 
        {
            init <- candidat
            densinit <- denscandidat
        }
  #      res <- c(res,init)
    }
   # res
    init
}

## genere une nouvelle coordonnée du posterior (krigeage simple et Matérn tensorisé)
generePosterior <- function(x_connus,y_connus, indice_derive, longueur_correlation,regularite)
{

#    dens_courante <- function (t)
#    {
#        densite2DPosterior(t, indice_derive = indice_derive, autre_longueur_correlation = autre_longueur_correlation,
#                           x_connus = x_connus, y_connus = y_connus, regularite = regularite)
#    }
    
    proposition <- Metropolis(densite = densitePosterior, sd_instrum = 0.2, init = 1, iterations = 20,
              indice_derive = indice_derive, longueur_correlation = longueur_correlation,
              x_connus = x_connus, y_connus = y_connus, regularite = regularite)
    
    list(proposition=proposition)
}



## Simule le posterior (krigeage simple) avec Matérn tensorisé
## probachoix : probabilité de choisir la 2e loi (cf genere2D)
GibbsPosterior <- function(x_connus, y_connus, longueur_correlation_init, regularite,nb_iterations)
    {
    lc <- longueur_correlation_init
    warnings <- 0
    PointsSimules <- matrix(NA,nrow=nb_iterations,ncol=ncol(x_connus))
    iter <- 0
    Plafond <- 2
#    for (iter in 1:nb_iterations)
    while(iter < nb_iterations)
    {
        iter <- iter + 1
        print(paste("iter =",iter))
        for(indice_derive in 1:length(longueur_correlation_init))
        {
            courant <- generePosterior(x_connus=x_connus,y_connus = y_connus,indice_derive=indice_derive, longueur_correlation=lc, regularite=regularite)
            lc[indice_derive] <- courant$proposition

            PointsSimules[iter,]<-lc 
        }
    }
    PointsSimules
}

evalueVraisemblanceIntegreeTens <- function( longueur_correlation, x_connus, y_connus, regularite)
    {
    Noyau <- creeMaternTensorise ( variance= 1, longueur= longueur_correlation, regularite= regularite)

    matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)

if(det(matriceCorrelation)!=0)    Vraisemblance_integree <- det(matriceCorrelation)^(-1/2)* (sum(y_connus* solve(matriceCorrelation,y_connus)))^(-nrow(x_connus)/2)
else Vraisemblance_integree <- 0
    Vraisemblance_integree
}

evalueVraisemblanceIntegreeAnisGeom <- function( longueur_correlation, x_connus, y_connus, regularite)
    {
    Noyau <- creeMaternIsotrope ( variance= 1, longueur= longueur_correlation, regularite= regularite)

    matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)

if(det(matriceCorrelation)!=0)    Vraisemblance_integree <- det(matriceCorrelation)^(-1/2)* (sum(y_connus* solve(matriceCorrelation,y_connus)))^(-nrow(x_connus)/2)
else Vraisemblance_integree <- 0
    Vraisemblance_integree

}

VraisemblanceIntegree2D<- function(t,autre_longueur_correlation,x_connus,y_connus,regularite=1.5)
    {
    taille <- length(t)
    res <- rep(NA,taille)

    for (i in 1:taille)
        {
        res[i] <- evalueVraisemblanceIntegree(x_connus = x_connus, y_connus = y_connus, regularite = regularite, longueur_correlation = c(t[i],autre_longueur_correlation))
        }
        

    res
}
