


#save.image("SamplingPosterior.RData")

##library(distr)
library(lhs)
library(pscl)
library(ggplot2)

x_connus_nuage<- generePointsConnus( nombre_points = 42, minima = c(0.1,0.1), maxima = c(0.9,0.9) ,lhs= TRUE)

plot(x_connus_nuage)

evaluePriorCorrelCond <- function( x_connus, regularite, longueur_correlation, indice_derive)
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
    else res <- sum( diag(MatriceAuxiliaire %*% MatriceAuxiliaire) ) - 1/nrow(x_connus) * (sum(diag(MatriceAuxiliaire)))^2
    sqrt(res)
}

evaluePriorCorrelCondIso <- function( x_connus, regularite, longueur_correlation, indice_derive)
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
    else res <- sum( diag(MatriceAuxiliaire %*% MatriceAuxiliaire) ) - 1/nrow(x_connus) * (sum(diag(MatriceAuxiliaire)))^2
    sqrt(res)
}

densite2D <- function(t,indice_derive,autre_longueur_correlation,x_connus,regularite=1.5)
    {
    taille <- length(t)
    res <- rep(NA,taille)
    print(autre_longueur_correlation)
    if(indice_derive==1)
        {
    for (i in 1:taille)
        {
        res[i] <- evaluePriorCorrelCond(x_connus = x_connus, regularite = regularite, longueur_correlation = c(t[i],autre_longueur_correlation),indice_derive = indice_derive)
    }
        }
    else
        {
        for (i in 1:taille)
        {
        res[i] <- evaluePriorCorrelCond(x_connus = x_connus, regularite = regularite, longueur_correlation = c(autre_longueur_correlation,t[i]),indice_derive = indice_derive)
        }
    }
    res
}

densite2DIso <- function(t,indice_derive,autre_longueur_correlation,x_connus,regularite=1.5)
    {
    taille <- length(t)
    res <- rep(NA,taille)
    if(indice_derive==1)
        {
    for (i in 1:taille)
        {
        res[i] <- evaluePriorCorrelCondIso(x_connus = x_connus, regularite = regularite, longueur_correlation = c(t[i],autre_longueur_correlation),indice_derive = indice_derive)
    }
        }
    else
        {
        for (i in 1:taille)
        {
        res[i] <- evaluePriorCorrelCondIso(x_connus = x_connus, regularite = regularite, longueur_correlation = c(autre_longueur_correlation,t[i]),indice_derive = indice_derive)
        }
    }
    res
}

## genere une nouvelle coordonnée du prior (krigeage simple et Matérn tensorisé)
genere2D <- function(x_connus,indice_derive, autre_longueur_correlation,regularite,probachoix)
    {

    
    accepte <- FALSE
    warnings <- 0
    
## Densité instrumentale = mélange de 2 lois
## 1e loi : épouse le comportement de la vraie densité au voisinage de 0   
    Ecarts <- abs(outer(x_connus[,indice_derive],x_connus[,indice_derive],"-"))
    Ecarts_sans_zeros <- Ecarts [Ecarts !=0]
    DistanceMinimale <- min(Ecarts_sans_zeros)
    alphaV0 <- 2
    betaV0 <- 2 * min(sqrt(regularite),1) * DistanceMinimale
    
## 2e loi : décroît assez lentement pour convenir à de nombreuses valeurs de autre_longueur_correlation  
    alphaVinf <- 0.1
    #betaVinf <- (1+alphaVinf)/(1+alphaV0) * betaV0
    betaVinf <- 0.01
    
#    ArgPlafond <- autre_longueur_correlation * 1.5
 #   if(ArgPlafond>1)
#        {
#    Vraie_dens_Argplafond <- densite2D(1, indice_derive=indice_derive, autre_longueur_correlation= autre_longueur_correlation, regularite=regularite, x_connus=x_connus) / ArgPlafond 
#    }
#    else
#        {
#            Vraie_dens_Argplafond <- densite2D(ArgPlafond, indice_derive=indice_derive, autre_longueur_correlation= autre_longueur_correlation, regularite=regularite, x_connus=x_connus) 
#}
#                Instr_dens_Argplafond <- (densigamma(ArgPlafond, alpha=alphaV0,beta=betaV0) + densigamma(ArgPlafond, alpha=alphaVinf,beta=betaVinf))/2
#
#    Plafond <- Vraie_dens_Argplafond / Instr_dens_Argplafond
#        print(paste("Plafond=",Plafond))
    Plafond <- 250
    Temps_plafond <- 0 ## pour débloquer l'algorithme quand la proba d'acceptation est trop faible
    
    while(!accepte)
    {
        if(Temps_plafond>100)
        {
            Plafond <- Plafond / 2
            Temps_plafond <- 0
            print("Plafond divisé par 2.")
        }
        choix <- rbinom(1,size=1,prob=probachoix)
#        print(paste("choix=",choix))
        if(choix==0) proposition <- rigamma(1,alpha=alphaV0,beta=betaV0)
        else proposition <- rigamma(1,alpha=alphaVinf,beta=betaVinf)

            
        densite <- densite2D(proposition, indice_derive=indice_derive, autre_longueur_correlation= autre_longueur_correlation, regularite=regularite, x_connus=x_connus)
        if(densite==0) warnings <- warnings + 0.01
            jugedepaix <- runif(1,max=Plafond)
        densiteInstrumentale <- (1-probachoix)*densigamma(proposition, alpha=alphaV0,beta=betaV0) + probachoix*densigamma(proposition, alpha=alphaVinf,beta=betaVinf)
        
#                    print("PROBACCEPTATION")
#            print(densite/densiteInstrumentale/Plafond)
#            print("")
        
        if(jugedepaix < densite/densiteInstrumentale) 
            {
                accepte <- TRUE
                print(proposition)
                if(Plafond < densite/densiteInstrumentale) warnings <- warnings + 1
            }
            
    }

        list(proposition=proposition, warnings=warnings, rapport= densite/densiteInstrumentale/Plafond)
}

## genere une nouvelle coordonnée du prior (krigeage simple et Matérn anisotrope géométrique)
genere2DIso <- function(x_connus,indice_derive, autre_longueur_correlation,regularite,probachoix)
    {

    
    accepte <- FALSE
    warnings <- 0
## Densité instrumentale = mélange de 2 lois
## 1e loi : épouse le comportement de la vraie densité au voisinage de 0
    Ecarts <- abs(outer(x_connus[,indice_derive],x_connus[,indice_derive],"-"))
    Ecarts_sans_zeros <- Ecarts [Ecarts !=0]
    DistanceMinimale <- min(Ecarts_sans_zeros)
    alphaV0 <- 2
    betaV0 <- 2 * min(sqrt(regularite),1) * DistanceMinimale
## 2e loi : décroît assez lentement pour convenir à de nombreuses valeurs de autre_longueur_correlation   
    alphaVinf <- 0.1
    #betaVinf <- (1+alphaVinf)/(1+alphaV0) * betaV0
    betaVinf <- 0.01
    
#    ArgPlafond <- autre_longueur_correlation * 1.5
 #   if(ArgPlafond>1)
#        {
#    Vraie_dens_Argplafond <- densite2DIso(1, indice_derive=indice_derive, autre_longueur_correlation= autre_longueur_correlation, regularite=regularite, x_connus=x_connus) / ArgPlafond 
#    }
#    else
#        {
#            Vraie_dens_Argplafond <- densite2DIso(ArgPlafond, indice_derive=indice_derive, autre_longueur_correlation= autre_longueur_correlation, regularite=regularite, x_connus=x_connus) 
#}
#                Instr_dens_Argplafond <- (densigamma(ArgPlafond, alpha=alphaV0,beta=betaV0) + densigamma(ArgPlafond, alpha=alphaVinf,beta=betaVinf))/2
#
#    Plafond <- Vraie_dens_Argplafond / Instr_dens_Argplafond
#        print(paste("Plafond=",Plafond))
    Plafond <- 600
    Temps_plafond <- 0 ## pour débloquer l'algorithme quand la proba d'acceptation est trop faible
    
    while(!accepte)
    {
        if(Temps_plafond>100)
        {
            Plafond <- Plafond / 2
            Temps_plafond <- 0
            print("Plafond divisé par 2.")
        }
        Temps_plafond <- Temps_plafond + 1
        choix <- rbinom(1,size=1,prob=probachoix)
#        print(paste("choix=",choix))
        if(choix==0) proposition <- rigamma(1,alpha=alphaV0,beta=betaV0)
        else proposition <- rigamma(1,alpha=alphaVinf,beta=betaVinf)

            
        densite <- densite2DIso(proposition, indice_derive=indice_derive, autre_longueur_correlation= autre_longueur_correlation, regularite=regularite, x_connus=x_connus)
        if(densite==0) warnings <- warnings + 0.01 ## La proposition est trop haute : vraie densité incalculable
            jugedepaix <- runif(1,max=Plafond)
        densiteInstrumentale <- (1-probachoix)*densigamma(proposition, alpha=alphaV0,beta=betaV0) + probachoix*densigamma(proposition, alpha=alphaVinf,beta=betaVinf)
        
#                    print("PROBACCEPTATION")
#            print(densite/densiteInstrumentale/Plafond)
#            print("")
        
        if(jugedepaix < densite/densiteInstrumentale) 
            {
                accepte <- TRUE
                print(proposition)
                if(Plafond < densite/densiteInstrumentale) warnings <- warnings + 1 #Rapport densités vraie/instrumentale > Plafond
            }
            
    }

        list(proposition=proposition, warnings=warnings, rapport= densite/densiteInstrumentale/Plafond)
}

## Simule le prior (krigeage simple) avec Matérn tensorisé
## probachoix : probabilité de choisir la 2e loi (cf genere2D)
Gibbs2D <- function(x_connus, longueur_correlation_init, regularite,nb_iterations,probachoix)
    {
    lc1 <- longueur_correlation_init[1]
    lc2 <- longueur_correlation_init[2]
    warnings <- 0
    PointsSimules <- matrix(NA,nrow=nb_iterations,ncol=5)
    for (iter in 1:nb_iterations)
        {
        courant <- genere2D(x_connus=x_connus,indice_derive=1, autre_longueur_correlation=lc2, regularite=regularite,probachoix=probachoix)
        lc1 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport1 <- courant$rapport
        courant <- genere2D(x_connus=x_connus,indice_derive=2, autre_longueur_correlation=lc1, regularite=regularite,probachoix=probachoix)
        lc2 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport2 <- courant$rapport
        PointsSimules[iter,]<-c(lc1,lc2,warnings,rapport1,rapport2)
    }
    print(c("warnings=",warnings))
    PointsSimules
}

## Simule le prior (krigeage simple) avec Matérn anisotrope géométrique
## probachoix : probabilité de choisir la 2e loi (cf genere2DIso)
Gibbs2DIso <- function(x_connus, longueur_correlation_init, regularite,nb_iterations,probachoix)
    {
    lc1 <- longueur_correlation_init[1]
    lc2 <- longueur_correlation_init[2]
    warnings <- 0
    PointsSimules <- matrix(NA,nrow=nb_iterations,ncol=5)
    for (iter in 1:nb_iterations)
        {
        courant <- genere2DIso(x_connus=x_connus,indice_derive=1, autre_longueur_correlation=lc2, regularite=regularite,probachoix=probachoix)
        lc1 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport1 <- courant$rapport
        courant <- genere2DIso(x_connus=x_connus,indice_derive=2, autre_longueur_correlation=lc1, regularite=regularite,probachoix=probachoix)
        lc2 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport2 <- courant$rapport
        PointsSimules[iter,]<-c(lc1,lc2,warnings,rapport1,rapport2)
    }
    print(c("warnings=",warnings))
    PointsSimules
}

## Comme Gibbs2D, mais simule en commençant par la 2e coordonnée
Gibbs2DAutreSens <- function(x_connus, longueur_correlation_init, nb_iterations)
    {
    lc1 <- longueur_correlation_init[1]
    lc2 <- longueur_correlation_init[2]
    warnings <- 0
    PointsSimules <- matrix(NA,nrow=nb_iterations,ncol=5)
    for (iter in 1:nb_iterations)
        {

        courant <- genere2D(x_connus=x_connus,indice_derive=2, autre_longueur_correlation=lc1)
        lc2 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport2 <- courant$rapport
        
        courant <- genere2D(x_connus=x_connus,indice_derive=1, autre_longueur_correlation=lc2)
        lc1 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport1 <- courant$rapport
        
        PointsSimules[iter,]<-c(lc1,lc2,warnings,rapport1,rapport2)
    }
    print(c("warnings=",warnings))
    PointsSimules
}

plot(iso3,xlim=c(0,2),ylim=c(0,2))
data.frame(iso3)

iso3 <- Gibbs2DIso(x_connus = x_connus_nuage,longueur_correlation_init = c(444.8566375,1603.6326257),regularite=1.5, probachoix=0.95, nb_iterations = 800)

essai <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(1,1),nb_iterations = 20)

essai2 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(0.59,5.6),nb_iterations = 20)

essai3 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(53,2.55),nb_iterations = 20)

essai4 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(197,0.07),nb_iterations = 20)

essai5 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(0.57,1.13),nb_iterations = 20)

essai6 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(4000,0.03),nb_iterations = 20)

new3 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(4761,100.99),nb_iterations = 100)

new4 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(1.06,0.16982359),nb_iterations = 10)

new5 <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(88,0.006157393),nb_iterations = 30)

big <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(14.66,42.11),nb_iterations = 700)

big[700,]

biggg <- Gibbs2D(x_connus = x_connus_nuage,longueur_correlation_init = c(0.0814573864098533,4.82550849474242),nb_iterations = 5)

bigAutreSens <- Gibbs2DAutreSens(x_connus = x_connus_nuage,longueur_correlation_init = c(1,1),nb_iterations = 1000)

big

ESSAI <-rbind(essai,essai2,essai3,essai4,essai5,essai6,new1[,1:2],new2[,1:2],new3[,1:2],new4[,1:2],new5[,1:2],big[,1:2])

ESSAI2<-matrix(ESSAI[rowSums(ESSAI)<100],ncol=2)
df2<-data.frame(ESSAI2)
nrow(ESSAI2)
ESSAI3<-matrix(ESSAI[apply(ESSAI,1,max)<10],ncol=2)
df3<-data.frame(ESSAI3)
nrow(ESSAI3)
ESSAI4<-matrix(ESSAI[apply(ESSAI,1,max)<1],ncol=2)
df4<-data.frame(ESSAI4)
nrow(ESSAI4)

df <- data.frame(x=rnorm(10000),y=rnorm(10000))
ggplot(df2,aes(x=X1,y=X2))+
    stat_density2d(aes(fill=..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="green")
dev.copy(pdf,"74_100.pdf")
dev.off()
ggplot(df2,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
dev.copy(pdf,"74_100_contour.pdf")
dev.off()
ggplot(df3,aes(x=X1,y=X2))+
    stat_density2d(aes(fill=..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="green")
dev.copy(pdf,"56_100.pdf")
dev.off()
ggplot(df3,aes(x=X1,y=X2))+geom_density2d()
dev.copy(pdf,"56_100_contour.pdf")
dev.off()
ggplot(df4,aes(x=X1,y=X2))+
    stat_density2d(aes(fill=..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="yellow")
dev.copy(pdf,"26_100.pdf")
dev.off()
ggplot(df4,aes(x=X1,y=X2))+geom_density2d()
dev.copy(pdf,"26_100_contour.pdf")
dev.off()

TRY <- rbind(try[,1:2],try2[,1:2],try3[,1:2],try4[,1:2],try5[,1:2])
TRY2<-matrix(TRY[rowSums(TRY)<100],ncol=2)
dftry2<-data.frame(TRY2)
nrow(TRY2)
TRY3<-matrix(TRY[apply(TRY,1,max)<10],ncol=2)
dftry3<-data.frame(TRY3)
nrow(TRY3)
TRY4<-matrix(TRY[apply(TRY,1,max)<1],ncol=2)
dftry4<-data.frame(TRY4)
nrow(TRY4)

ggplot(dftry2,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
ggplot(dftry3,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
ggplot(dftry4,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)

ISO<- rbind(iso3[167:637,1:2])
ISO2<-matrix(ISO[rowSums(ISO)<100],ncol=2)
dfiso2<-data.frame(ISO2)
nrow(ISO2)
ISO3<-matrix(ISO[apply(ISO,1,max)<10],ncol=2)
dfiso3<-data.frame(ISO3)
nrow(ISO3)
ISO4<-matrix(ISO[apply(ISO,1,max)<1],ncol=2)
dfiso4<-data.frame(ISO4)
nrow(ISO4)

ggplot(dfiso2,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
ggplot(dfiso3,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
plot(ISO3,xlim=c(0,3),ylim=c(0,3))
ggplot(dfiso4,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)

ISO<- rbind(iso1[,1:2],iso2[,1:2],iso3[1:166,1:2])
ISO2<-matrix(ISO[rowSums(ISO)<100],ncol=2)
dfiso2<-data.frame(ISO2)
nrow(ISO2)
ISO3<-matrix(ISO[apply(ISO,1,max)<10],ncol=2)
dfiso3<-data.frame(ISO3)
nrow(ISO3)
ISO4<-matrix(ISO[apply(ISO,1,max)<1],ncol=2)
dfiso4<-data.frame(ISO4)
nrow(ISO4)

ggplot(dfiso2,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
ggplot(dfiso3,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
ggplot(dfiso4,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)

bas<-bigAutreSens[,1:2]
dfa<-data.frame(bas)

bas2<-matrix(bas[rowSums(bas)<100],ncol=2)
dfa2<-data.frame(bas2)
nrow(bas2)

bas3<-matrix(bas[apply(bas,1,max)<10],ncol=2)
dfa3<-data.frame(bas3)
nrow(bas3)

bas4<-matrix(bas[apply(bas,1,max)<1],ncol=2)
dfa4<-data.frame(bas4)
nrow(bas4)


ggplot(dfa2,aes(x=X1,y=X2))+
    stat_density2d(aes(fill=..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="green")
dev.copy(pdf,"Simul2_80_100.pdf")
dev.off()
ggplot(dfa2,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
dev.copy(pdf,"Simul2_80_100_contour.pdf")
dev.off()
ggplot(dfa3,aes(x=X1,y=X2))+
    stat_density2d(aes(fill=..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="green")
dev.copy(pdf,"Simul2_60_100.pdf")
dev.off()
ggplot(dfa3,aes(x=X1,y=X2))+geom_density2d()
dev.copy(pdf,"Simul2_60_100_contour.pdf")
dev.off()
ggplot(dfa4,aes(x=X1,y=X2))+
    stat_density2d(aes(fill=..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="yellow")
dev.copy(pdf,"Simul2_27_100.pdf")
dev.off()
ggplot(dfa4,aes(x=X1,y=X2))+geom_density2d()
dev.copy(pdf,"Simul2_27_100_contour.pdf")
dev.off()

plot(rbind(essai,essai2,essai3,essai4,essai5,essai6,new1[,1:2],new2[,1:2],new3[,1:2],big[,1:2]),xlim=c(0,10),ylim=c(0,10))
