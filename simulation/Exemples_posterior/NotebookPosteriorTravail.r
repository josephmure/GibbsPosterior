
#save.image("SamplingPosterior.RData")

#load("SamplingPosterior.RData")



#library(distr)
library(lhs)
library(pscl)
library(ggplot2)
library(rgl)

#x_connus20<- generePointsConnus( nombre_points = 20, minima = c(0,0), maxima = c(1,1) ,lhs= TRUE)

#plot(x_connus20)

#plot(y_connus20)

#NOYAU <- creeMaternTensorise(variance = 1, longueur = c(0.3,0.6), regularite = 1.5)

#y_connus20 <- creeEchantillonNormal(x_connus = x_connus20, noyau = NOYAU)

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

evalueVraisemblanceIntegree <- function( x_connus, y_connus, regularite, longueur_correlation)
    {
    Noyau <- creeMaternTensorise ( variance= 1, longueur= longueur_correlation, regularite= regularite)

    matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)

    Vraisemblance_integree <- det(matriceCorrelation)^(-1/2)* (sum(y_connus* solve(matriceCorrelation,y_connus)))^(-nrow(x_connus)/2)
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

densite2DPosterior <- function(t,indice_derive,autre_longueur_correlation,x_connus,y_connus,regularite=1.5)
    {
    taille <- length(t)
    res <- rep(NA,taille)

    if(indice_derive==1)
        {
    for (i in 1:taille)
        {
        res[i] <- evaluePosteriorCorrelCond(x_connus = x_connus, y_connus = y_connus, regularite = regularite, longueur_correlation = c(t[i],autre_longueur_correlation),indice_derive = indice_derive)
    }
        }
    else
        {
        for (i in 1:taille)
        {
        res[i] <- evaluePosteriorCorrelCond(x_connus = x_connus, y_connus = y_connus, regularite = regularite, longueur_correlation = c(autre_longueur_correlation,t[i]),indice_derive = indice_derive)
        }
    }

    res
}

## genere une nouvelle coordonnée du posterior (krigeage simple et Matérn tensorisé)
genere2DPosterior <- function(x_connus,y_connus, indice_derive, autre_longueur_correlation,regularite)
    {

    
    accepte <- FALSE
    warnings <- 0
    puissance <- 0
    puissanceMax <- 4
    while(puissance < puissanceMax)
    {
        upper <- 10^puissance
        tol <- 0.1*upper
        opt <-optimize(densite2DPosterior,lower=0, upper=upper, tol=tol, maximum=TRUE, 
                       indice_derive=indice_derive,autre_longueur_correlation=autre_longueur_correlation,
                       x_connus=x_connus, y_connus = y_connus)

        if(opt$objective < (upper - tol) ) puissance<- puissanceMax
        puissance <- puissance + 1
    }

    Plafond <- 2

    alpha <- 0
    multiplicateur <- 1
        nb_passages <- 0
    while(alpha<=0)
        {
        nb_passages <- nb_passages + 1
        multiplicateur <- multiplicateur + 0.5
    dens<-densite2DPosterior(multiplicateur*opt$maximum, indice_derive=indice_derive,autre_longueur_correlation=autre_longueur_correlation,
                             x_connus=x_connus, y_connus = y_connus)
    rapport <- opt$objective/dens
    alpha <- log(rapport)/log(multiplicateur) - 1
    beta <- (alpha+1)*opt$maximum
    if(nb_passages>3) 
        {
            alpha <- 0.1
            beta <- 0.01
         }
    print(paste("alpha=",alpha,", beta=",beta))
        }

    
    while(!accepte)
    {
        proposition <- rigamma(1,alpha=alpha,beta=beta)            
        densite <- densite2DPosterior(proposition, indice_derive=indice_derive, autre_longueur_correlation= autre_longueur_correlation, regularite=regularite, x_connus=x_connus,y_connus = y_connus)
        

        if(densite==0) warnings <- warnings + 0.00001
        jugedepaix <- runif(1,max=Plafond)
        densiteInstrumentale <- densigamma(proposition, alpha=alpha,beta=beta) 

  
        
        if(densite!=0)
        {
        if(jugedepaix < densite/densiteInstrumentale/opt$objective*densigamma(opt$maximum,alpha=alpha,beta=beta)) 
            {
                accepte <- TRUE
                print(paste("proposition =",proposition))
                if(Plafond < densite/densiteInstrumentale/opt$objective*densigamma(opt$maximum,alpha=alpha,beta=beta)) warnings <- warnings + 1
            }
        }   
    }

        list(proposition=proposition, warnings=warnings, rapport= densite/densiteInstrumentale/opt$objective*densigamma(opt$maximum,alpha=alpha,beta=beta)/Plafond)
}

## Simule le posterior (krigeage simple) avec Matérn tensorisé
## probachoix : probabilité de choisir la 2e loi (cf genere2D)
Gibbs2DPosterior <- function(x_connus, y_connus, longueur_correlation_init, regularite,nb_iterations)
    {
    lc1 <- longueur_correlation_init[1]
    lc2 <- longueur_correlation_init[2]
    warnings <- 0
    PointsSimules <- matrix(NA,nrow=nb_iterations,ncol=5)
    for (iter in 1:nb_iterations)
        {
        print(paste("iter =",iter))
        courant <- genere2DPosterior(x_connus=x_connus,y_connus = y_connus,indice_derive=1, autre_longueur_correlation=lc2, regularite=regularite)
        lc1 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport1 <- courant$rapport
        courant <- genere2DPosterior(x_connus=x_connus,y_connus = y_connus,indice_derive=2, autre_longueur_correlation=lc1, regularite=regularite)
        lc2 <- courant$proposition
        warnings <- warnings + courant$warnings
        rapport2 <- courant$rapport
        PointsSimules[iter,]<-c(lc1,lc2,warnings,rapport1,rapport2)
    }
    print(c("warnings=",warnings))
    PointsSimules
}

trouveBeta <- function(alpha,theta1,theta2,indice_derive,autre_longueur_correlation,x_connus,y_connus,regularite)
{
	rapport <- densite2DPosterior(theta2,indice_derive,autre_longueur_correlation,x_connus,y_connus,regularite) / densite2DPosterior(theta1,indice_derive,autre_longueur_correlation,x_connus,y_connus,regularite)

	beta <- - ((alpha+1) * log(theta2/theta1) + log(rapport))/ (1/theta2 - 1/theta1)

	beta
}


# Z20 <- Gibbs2DPosterior(x_connus = x_connus20,y_connus = y_connus20, longueur_correlation_init = c(0.5,0.5),regularite=1.5, nb_iterations = 1000)
# 
# 
# Z20COMPLET <- rbind(Z20[,1:2])
# plot(Z20COMPLET,xlim=c(0,1),ylim=c(0,1))
# DFZ20COMPLET  <- data.frame(Z20COMPLET )
# ggplot(DFZ20COMPLET ,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
# ggplot(DFZ20COMPLET ,aes(x=X1,y=X2))+
#     stat_density2d(aes(fill=..level..), geom="polygon") +
#     scale_fill_gradient(low="blue", high="yellow")
# 
# #Zdisperse <- Gibbs2DPosterior(x_connus = x_connus_nuage,y_connus = y_connus_nuage_disperse, longueur_correlation_init = c(1,1),regularite=1.5, nb_iterations = 100)
# Zdisperse2 <- Gibbs2DPosterior(x_connus = x_connus_nuage,y_connus = y_connus_nuage_disperse, longueur_correlation_init = c(0.05304419,0.55978198),regularite=1.5, nb_iterations = 100)
# 
# ZDISPERSE <- rbind(Zdisperse[,1:2],Zdisperse2[,1:2])
# plot(ZDISPERSE,xlim=c(0,1),ylim=c(0,1))
# DFZDISPERSE  <- data.frame(ZDISPERSE )
# ggplot(DFZDISPERSE ,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
# ggplot(DFZDISPERSE ,aes(x=X1,y=X2))+
#     stat_density2d(aes(fill=..level..), geom="polygon") +
#     scale_fill_gradient(low="blue", high="yellow")
# 
# essai2 <- Gibbs2DPosterior(x_connus = x_connus_nuage,y_connus = y_connus_nuage, longueur_correlation_init = c(0.2719907,0.3117847),regularite=1.5, nb_iterations = 900)
# 
# ## avec 42 points LHS
# ESSAI <- rbind(essai[,1:2],essai2[,1:2])
# plot(ESSAI,xlim=c(0,1),ylim=c(0,1))
# DF <- data.frame(ESSAI)
# ggplot(DF,aes(x=X1,y=X2))+geom_density2d(contour=TRUE)
# ggplot(DF,aes(x=X1,y=X2))+
#     stat_density2d(aes(fill=..level..), geom="polygon") +
#     scale_fill_gradient(low="blue", high="yellow")
# 
# densNorm <- integrate(densite2DPosterior,lower=0, upper=Inf, rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)$value
# 
# system.time(integrate(densite2DPosterior,lower=0, upper=Inf, rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage))
# 
# integrate(densite2DPosterior,lower=0, upper=0.3, rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# 
# 
# system.time(integrate(densite2DPosterior,lower=0, upper=0.35, rel.tol=10, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage))
# 
# optimize(densite2DPosterior,lower=0, upper=1, tol=1, maximum=TRUE, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# 
# int <-integrate(densite2DPosterior,lower=0, upper=0.1, rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# int$value/densNorm
# 
# integrate(densite2DPosterior,lower=0, upper=0.285, rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# 
# 
# integrate(densite2DPosterior,lower=0, upper=0.2, rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# 
# evaluePosteriorCorrelCond
# 
# opt <-optimize(densite2DPosterior,lower=0, upper=1, tol=0.1, maximum=TRUE, indice_derive=1,autre_longueur_correlation=1,x_connus=x_connus_nuage, y_connus = y_connus_nuage_disperse)
# 
# dens<-densite2DPosterior(1.5*opt$maximum, indice_derive=1,autre_longueur_correlation=1,x_connus=x_connus_nuage, y_connus = y_connus_nuage_disperse)
# rapport= opt$objective/dens
# alpha <- log(rapport)/log(1.5) - 1
# alpha
# 
# abscisses <- 0.01*1:100
# invgamma <- densigamma(abscisses,alpha=alpha, beta=(alpha+1)*opt$maximum)
# vraiedens <- densite2DPosterior(abscisses,indice_derive=1,autre_longueur_correlation=1,x_connus=x_connus_nuage, y_connus = y_connus_nuage_disperse)
# vraiedens[10:100]/opt$objective/invgamma[10:100]*densigamma(opt$maximum,alpha=alpha,beta=(alpha+1)*opt$maximum)
# plot(abscisses,invgamma/densigamma(opt$maximum,alpha=alpha,beta=(alpha+1)*opt$maximum))
# points(abscisses,vraiedens/opt$objective,col=4)
# points(abscisses[10:100],vraiedens[10:100]/opt$objective/invgamma[10:100]*densigamma(opt$maximum,alpha=alpha,beta=(alpha+1)*opt$maximum),col=2)
# 
# dens<-densite2DPosterior(2*opt$maximum, indice_derive=1,autre_longueur_correlation=0.1,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# rapport= opt$objective/dens
# alpha <- log(rapport)/log(1.5) - 1
# alpha
# 
# abscisses <- 0.01*1:100
# invgamma <- densigamma(abscisses,alpha=alpha, beta=(alpha+1)*opt$maximum)
# vraiedens <- densite2DPosterior(abscisses,indice_derive=1,autre_longueur_correlation=0.1,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# plot(abscisses,invgamma/densigamma(opt$maximum,alpha=alpha,beta=(alpha+1)*opt$maximum))
# points(abscisses,vraiedens/opt$objective,col=4)
# points(abscisses[10:100],vraiedens[10:100]/opt$objective/invgamma[10:100]*densigamma(opt$maximum,alpha=alpha,beta=(alpha+1)*opt$maximum),col=2)
# 
# test <- function (t) integrate(densite2DPosterior,lower=0,upper=t,rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)$value/5.155683e-17
# densNorm<- function(t) densite2DPosterior(t, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)/5.155683e-17 
# 
# integrate(densite2DPosterior,lower=0, upper=Inf, rel.tol=1, indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# opt <-optimize(densNorm,lower=0, upper=1, tol=1, maximum=TRUE)
# opt
# nouveau <- test(opt$maximum)
# nouveau
# Point <- opt$maximum + (0.5-nouveau)/opt$objective
# Point
# nnouveau <- test(Point)
# nnouveau
# PPoint <- Point + (0.5-nnouveau)/densNorm(Point)
# PPoint
# nnnouveau <- test(PPoint)
# nnnouveau
# PPPoint <- PPoint + (0.5-nnnouveau)/densNorm(PPoint)
# PPPoint
# nnnnouveau <- test(PPPoint)
# nnnnouveau
# 
# testtt <- function(t) ( test(t) - 0.2 )^2
# 
# optimize(testtt,lower=0,upper=1,tol=0.01)
# 
# system.time(optimize(testtt,lower=0,upper=1,tol=0.01))
# 
# optimize(testtt,lower=0,upper=1,tol=0.2)
# 
# system.time(optimize(testtt,lower=0,upper=1,tol=0.1))
# 
# test(0.26940131435081)
# 
# test(0.315299338891701)
# 
# ?optimize
# 
# 
# abscisses <- c(0.01*1:100)
# ordonnees<-densite2DPosterior(abscisses,indice_derive=1,autre_longueur_correlation=0.2,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# ordonnees2<-densigamma(abscisses,alpha=0.1,beta=0.01)
# 
# plot(abscisses,ordonnees*10^10/ordonnees2,col=2)
# points(abscisses,ordonnees*10^10)
# lines(abscisses, ordonnees2)
# 
# 
# abscisses <- c(0.01*1:100)
# #ordonnees<-densite2DPosterior(abscisses,indice_derive=1,autre_longueur_correlation=0.1,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# ordonnees2<-densigamma(abscisses[10:100],alpha=2,beta=1.2)
# 
# plot(abscisses[10:100],ordonnees[10:100]/ordonnees2/4*10^22,col=2)
# points(abscisses[10:100],ordonnees[10:100]/4*10^22)
# lines(abscisses[10:100], ordonnees2)
# 
# integrate(f=densite2DPosterior, lower=0,upper=Inf, rel.tol=1,indice_derive=1,autre_longueur_correlation=0.1,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# 
# optimize(f=densite2DPosterior,maximum=TRUE, tol=1,lower=0,upper=1, indice_derive=1,autre_longueur_correlation=0.1,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# 
# ?integrate
# 
# 
# abscisses <- c(0.01*1:100)
# ordonnees<-VraisemblanceIntegree2D(abscisses,autre_longueur_correlation=0.1,x_connus=x_connus_nuage, y_connus = y_connus_nuage)
# plot(abscisses,ordonnees)
