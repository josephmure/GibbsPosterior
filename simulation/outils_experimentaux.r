
# Loi a posteriori


#set.seed(8)

#library(distr)
library(lhs)
library(pscl)
library(ggplot2)
#library(rgl)



# Required packages -- POURQUOI ?
library(mvtnorm)
library(ks)


evaluePosteriorCorrelCondAnisGeom <- function( x_connus, y_connus, regularite, longueur_correlation, indice_derive)
    {
    Noyau <- creeMaternIsotrope ( variance= 1, longueur= longueur_correlation, regularite= regularite)

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

evalueInfoFisherAnisGeom <- function( x_connus, regularite, longueur_correlation, indice_derive1, indice_derive2, tendance=NULL)
{

  
  Noyau <- creeMaternIsotrope ( variance= 1, longueur= longueur_correlation, regularite= regularite)
  
  matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)
  ##    inverseMatriceCorrelation <- solve(matriceCorrelation)
  
  if(is.null(tendance)) injectionOrthogonal <- diag(rep(1,nrow(x_connus)))
  else   injectionOrthogonal <- qr.Q(qr(tendance),complete=TRUE)[,(1+ncol(tendance)):ncol(matriceCorrelation)]
  
  matriceDerivee1 <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= creeMaternIsotropeDeriveLongueur( variance=1, longueur= Noyau$parametres$longueur, regularite= Noyau$parametres$regularite, indice_derive= indice_derive1) )
  matriceDerivee2 <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= creeMaternIsotropeDeriveLongueur( variance=1, longueur= Noyau$parametres$longueur, regularite= Noyau$parametres$regularite, indice_derive= indice_derive2) )
  
  matriceCorrelation <- t(injectionOrthogonal) %*% matriceCorrelation %*% injectionOrthogonal
  matriceDerivee1 <- t(injectionOrthogonal) %*% matriceDerivee1 %*% injectionOrthogonal
  matriceDerivee2 <- t(injectionOrthogonal) %*% matriceDerivee2 %*% injectionOrthogonal
  
  
  MatriceAuxiliaire1 <- try(solve(matriceCorrelation, matriceDerivee1))
  MatriceAuxiliaire2 <- try(solve(matriceCorrelation, matriceDerivee2))
  

    res <- sum( diag(MatriceAuxiliaire1 %*% MatriceAuxiliaire2) ) - 1/nrow(x_connus) * (sum(diag(MatriceAuxiliaire1))) * (sum(diag(MatriceAuxiliaire2)))

  res
}

creeInfoFisherTriangulaire <- function(x_connus, regularite, longueur_correlation, tendance=NULL)
{
  res <- matrix(0,length(longueur_correlation),length(longueur_correlation))
  
  for(i in 1:length(longueur_correlation))
  {
    for(j in i:length(longueur_correlation))
    {
      res[i,j] <- evalueInfoFisherAnisGeom(x_connus= x_connus , regularite = regularite, longueur_correlation = longueur_correlation, indice_derive1 = i, indice_derive2 = j, tendance = tendance)
    }
  }
  
  res
}

creeInfoFisherDiagonale <- function(x_connus, regularite, longueur_correlation, tendance=NULL)
{
  res <- rep(0,length(longueur_correlation),length(longueur_correlation))
  
  for(i in 1:length(longueur_correlation))
  {
      res[i] <- evalueInfoFisherAnisGeom(x_connus= x_connus , regularite = regularite, longueur_correlation = longueur_correlation, indice_derive1 = i, indice_derive2 = i, tendance = tendance)

  }
  
  res
}



densitePosteriorAnisGeom <- function(t,indice_derive, longueur_correlation,x_connus,y_connus,regularite=1.5)
    {
    if(t<=0)res <- 0
    else
    {
        longueur_correlation[indice_derive] <- t

        res <- evaluePosteriorCorrelCondAnisGeom(x_connus = x_connus, y_connus = y_connus, regularite = regularite, longueur_correlation = longueur_correlation,indice_derive = indice_derive)
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
generePosteriorAnisGeom <- function(x_connus,y_connus, indice_derive, longueur_correlation,regularite, sd_instrum = 0.2, iterations = 20)
{

#    dens_courante <- function (t)
#    {
#        densite2DPosterior(t, indice_derive = indice_derive, autre_longueur_correlation = autre_longueur_correlation,
#                           x_connus = x_connus, y_connus = y_connus, regularite = regularite)
#    }
    
    proposition <- Metropolis(densite = densitePosteriorAnisGeom, sd_instrum = sd_instrum, init = 1, iterations = iterations,
              indice_derive = indice_derive, longueur_correlation = longueur_correlation,
              x_connus = x_connus, y_connus = y_connus, regularite = regularite)
    
    list(proposition=proposition)
}



## Simule le posterior (krigeage simple) avec Matérn tensorisé
## probachoix : probabilité de choisir la 2e loi (cf genere2D)
GibbsPosteriorAnisGeom <- function(x_connus, y_connus, longueur_correlation_init, regularite,nb_iterations)
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
            courant <- generePosteriorAnisGeom(x_connus=x_connus,y_connus = y_connus,indice_derive=indice_derive, longueur_correlation=lc, regularite=regularite)
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


source("../Krige_these.r")
source("../Matern_these.r")
source("../PosteriorNettoyeAnisGeom.r")

planXP <- matrix(scan("planXP.txt"),ncol=10,byrow=TRUE)
observations <- scan("observations.txt")
longueur_correlation <- rep(10,10)

sd_instrum=1
iterations=80

ddd <- NULL


for (i in 1:100) 
{
	ddd <- c(ddd,  generePosteriorAnisGeom(x_connus=planXP, y_connus=observations, indice_derive=1, longueur_correlation=longueur_correlation, regularite=2.5, sd_instrum=sd_instrum,iterations=iterations) )
	print(c(i,ddd[i]))
}




sd_instrum=1
iterations=40

ccc <- NULL


for (i in 1:100) 
{
	ccc <- c(ccc,  generePosteriorAnisGeom(x_connus=planXP, y_connus=observations, indice_derive=1, longueur_correlation=longueur_correlation, regularite=2.5, sd_instrum=sd_instrum,iterations=iterations) )
	print(c(i,ccc[i]))
}


XXX <- 0.5*0:100
YYY <- apply(matrix(XXX),1,densitePosteriorAnisGeom,indice_derive=1, longueur_correlation=c(20,rep(10,9)),x_connus=planXP,y_connus=observations,regularite=2.5)


