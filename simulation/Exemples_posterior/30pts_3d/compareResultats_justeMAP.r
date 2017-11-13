#library(distr)
library(pscl)
#library(ggplot2)
#library(rgl)


source("../../../../Matern_these.r")
source("../../../../Briques_these.r")
source("../../../../Krige_these.r")
source("../../../../dichotomie.r")

REGULARITE <- read.table("../regularite_vraie.txt",as.is=TRUE)$V1
LONGUEUR_CORRELATION <- read.table("../longueur_correlation_vraie.txt",as.is=TRUE)$V1
LONGUEUR_CORRELATION_ARGMAX_VRAISEMBLANCE <- read.table("argmax_vraisemblance_integree.txt")$V1
LONGUEUR_CORRELATION_MODE_POSTERIOR <- read.table("mode_posterior.txt")$V1
NOMBRE_DIMENSIONS <- length( read.table("../longueur_correlation_vraie.txt",as.is=TRUE)$V1 ) ## Dimension = nombre de longueurs de correlation renseignees dans ../longueur_correlation_vraie.txt

TYPE_NOYAU_MATERN_TYPE_PRIOR <- c(as.matrix(read.table("../type_noyau_prior.txt")))
TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
#TYPE_PRIOR <- TYPE_NOYAU_MATERN_TYPE_PRIOR[2] ##Pas besoin ici.
CARDINAL_ENSEMBLE_TEST <- 100
POURCENTAGE_INTERVALLE_CONFIANCE <- 95
#POURCENTAGE_INTERVALLE_CONFIANCE_BAYESIEN <- 90

QNORM <- qnorm(1 - (1-POURCENTAGE_INTERVALLE_CONFIANCE/100)/2) ## QNORM est le quantile correspondant a la borne droite de l'intervalle de confiance/pari de la loi normale centree reduite entrant en jeu dans les calculs a suivre
x_connus <- as.matrix(read.table("planXP.txt", as.is = TRUE))
y_connus <- scan("observations.txt")
TAILLE_PLAN_XP <- nrow(x_connus)

tendance <- as.matrix(read.table("tendance.txt", as.is=TRUE))
#if(is.numeric(tendance)) 
#{
#	NOMBRE_FONCTIONS_TENDANCE <- ncol(tendance)
#}else 
#{
#	NOMBRE_FONCTIONS_TENDANCE <- 0
#}

NOMBRE_FONCTIONS_TENDANCE <- length(FONCTIONS)
NOMBRE_FONCTIONS_REELLES_TENDANCE <- length(FONCTIONS_REELLES)


if(TYPE_NOYAU_MATERN == "geometrique") 
{
	NOYAU <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE ) # anisotrope geometrique en fait !

	NOYAU_MODE_POSTERIOR <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION_MODE_POSTERIOR, regularite = REGULARITE ) # anisotrope geometrique en fait !
} else if (TYPE_NOYAU_MATERN == "tensorise") 
{
	NOYAU <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE )

	NOYAU_MODE_POSTERIOR <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION_MODE_POSTERIOR, regularite = REGULARITE ) # anisotrope geometrique en fait !
} else 
{
	print("Erreur dans compareResultats.r : le noyau doit-il etre anis. geometrique ou tensorise ?")
	NOYAU <- NULL

#	NOYAU_MODE_POSTERIOR <- NULL
}

if(!is.null(NOYAU)) ## Si on a bien precise la nature anisotrope geometrique ou tensorisee du noyau
{
    ## Matrice de correlation juste : correspondant aux vraies longueurs de correlation
    MATRICE_JUSTE <- creeMatriceCovariance(x1=x_connus, x2=x_connus, noyau= NOYAU)
    MATRICE_JUSTE_INVERSE <- solve(MATRICE_JUSTE) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.



    
    ## Matrice de correlation  correspondant aux longueurs de correlation estimees par MAP : max a posteriori 


    MATRICE_MODE_POSTERIOR <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= NOYAU_MODE_POSTERIOR)
    MATRICE_MODE_POSTERIOR_INVERSE <- solve(MATRICE_MODE_POSTERIOR) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.



    
    ## Recherche de la distance entre matrice de correl correspondant aux vraies longueurs de correlation et les matrices de correl obtenues avec MLE et MAP

    ## Transfo de Fisher appliquee aux matrices de correlation
	MATRICE_JUSTE_FISHER <- 1/2 * log( (1 +  MATRICE_JUSTE) / (1 - MATRICE_JUSTE) ) # transfo de Fisher
	diag(MATRICE_JUSTE_FISHER) <- 0
#print(MATRICE_JUSTE_FISHER)




	MATRICE_MODE_POSTERIOR_FISHER <- 1/2 * log( (1 +  MATRICE_MODE_POSTERIOR) / (1 - MATRICE_MODE_POSTERIOR) ) # transfo de Fisher
	diag(MATRICE_MODE_POSTERIOR_FISHER) <- 0
#print(MATRICE_MODE_POSTERIOR_FISHER)

    ## Distances "simple" et "de Fisher" entre matrices obtenues par MAP et la vraie matrice


	Distance_mode_posterior <- c( norm(MATRICE_MODE_POSTERIOR - MATRICE_JUSTE, "F"), norm(MATRICE_MODE_POSTERIOR_FISHER - MATRICE_JUSTE_FISHER, "F") ) # norme de Frobenius (ie euclidienne)
#print(Distance_mode_posterior)


    write.matrix(Distance_mode_posterior,"mode_posterior_distance.txt",sep = "\t")



    

# Ecart quadratique entre la prediction (moyenne conditionnelle sachant la valeur estimee des longueurs de correlation)
# et la moyenne conditionnelle sachant la vraie valeur des longueurs de correlation


## Generation de l'ensemble test
    Nouveaux_points <- generePointsConnus(nombre_points = CARDINAL_ENSEMBLE_TEST, minima = rep(0,NOMBRE_DIMENSIONS), maxima = rep(1,NOMBRE_DIMENSIONS), lhs=FALSE)
    Tous_les_points <- rbind(x_connus, Nouveaux_points)


    ## Matrice de variance marginale correspondant aux Nouveaux_points
    MATRICE_NOUVEAUX_JUSTE <- creeMatriceCovariance(x1= Nouveaux_points, x2= Nouveaux_points, noyau= NOYAU)
    ## INUTILE : tout ce qui compte, ce sont les éléments diagonaux, qui sont tous égaux à 1

    	#MATRICE_NOUVEAUX_MODE_POSTERIOR <- creeMatriceCovariance(x1= Nouveaux_points, x2= Nouveaux_points, noyau= NOYAU_MODE_POSTERIOR)

	MATRICE_NOUVEAUX_MODE_POSTERIOR <- 1
    
# Attention, non seulement la technique suivante est plus efficace que les commentaires, mais elle rend aussi les denominations plus justes. Par consequent,
# CORREL_CONNUS_NOUVEAUX_... actuel = t( CORREL_CONNUS_NOUVEAUX_...  ancien ) et Moyenne_conditionnelle_juste, ainsi que les predictions, sont des vecteurs lignes.
#MATRICE_JUSTE_TOTALE <- creeMatriceCovariance(x1= Tous_les_points, x2= Tous_les_points, noyau= NOYAU)
#CORREL_CONNUS_NOUVEAUX_JUSTE <- MATRICE_JUSTE_TOTALE[(nrow(x_connus)+1):(nrow(x_connus)+CARDINAL_ENSEMBLE_TEST),1:nrow(x_connus)]
CORREL_CONNUS_NOUVEAUX_JUSTE <- creeMatriceCovariance(x1= x_connus, x2= Nouveaux_points, noyau= NOYAU)


#MATRICE_MODE_POSTERIOR_TOTALE <- creeMatriceCovariance(x1= Tous_les_points, x2= Tous_les_points, noyau= NOYAU_MODE_POSTERIOR)
CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR <- creeMatriceCovariance(x1= x_connus, x2= Nouveaux_points, noyau= NOYAU_MODE_POSTERIOR)


### ATTENTION : JUSTE merite un traitement separe, car beta est connu.

    ## Calcul des moyennes marginales (qui dependent uniquemeent de beta, et non de sigma^2 ou theta)
    Moyenne_marginale_juste <- 0
    if(NOMBRE_FONCTIONS_REELLES_TENDANCE>0)
    {        
        for(i in 1:NOMBRE_FONCTIONS_REELLES_TENDANCE)
        {
            Moyenne_marginale_juste <- Moyenne_marginale_juste + BETA[i] * apply(Nouveaux_points,1,FONCTIONS_REELLES[[i]])
        }
        Moyenne_marginale_juste <- t(Moyenne_marginale_juste)

    ## Calcul des moyenne conditionnelles sachant les valeurs observees et les longeurs de correlation vraies ou estimees
    Moyenne_conditionnelle_juste <- t(y_connus - tendance_reelle %*% BETA) %*% MATRICE_JUSTE_INVERSE %*%  CORREL_CONNUS_NOUVEAUX_JUSTE + Moyenne_marginale_juste
    } else
	{
		Moyenne_conditionnelle_juste <- t(y_connus) %*% MATRICE_JUSTE_INVERSE %*% CORREL_CONNUS_NOUVEAUX_JUSTE
	}
    
    
    


    



        ## Calcul des vraies valeurs du processus gaussien aux points de l'ensemble test (Attention, ces valeurs sont presentees dans un vecteur LIGNE)
    Variance_conditionnelle_juste <-  MATRICE_NOUVEAUX_JUSTE - t(CORREL_CONNUS_NOUVEAUX_JUSTE) %*% MATRICE_JUSTE_INVERSE %*% CORREL_CONNUS_NOUVEAUX_JUSTE 
    ALEA <-   rnorm(CARDINAL_ENSEMBLE_TEST)
    Valeurs_nouveaux_points <- Moyenne_conditionnelle_juste + t( ALEA) %*% chol( Variance_conditionnelle_juste )
    #print((Valeurs_nouveaux_points- Moyenne_conditionnelle_juste) %*% solve(Variance_conditionnelle_juste, t(Valeurs_nouveaux_points - Moyenne_conditionnelle_juste) ) )
    #print(t(ALEA) %*% ALEA)
    #print(max(Variance_conditionnelle_juste) )# - t(chol(Variance_conditionnelle_juste)) %*% chol(Variance_conditionnelle_juste) )))



    ## Calcul de l'intervalle de pari correspondant aux vraies longueurs de correlation. Attention, chaque intervalle de pari est une COLONNE
    Intervalle_pari_juste <- rbind(Moyenne_conditionnelle_juste - QNORM * sqrt(diag(Variance_conditionnelle_juste)) , Moyenne_conditionnelle_juste + QNORM * sqrt(diag(Variance_conditionnelle_juste)))
    Longueur_intervalle_pari_juste <- Intervalle_pari_juste[2,] - Intervalle_pari_juste[1,]
    Longueur_moyenne_intervalle_pari_juste<- mean(Longueur_intervalle_pari_juste)



    ## Calcul de la frequence a laquelle la vraie valeur tombe dans l'intervalle de confiance correspondant aux vraies longueurs de correlation
    Appartient_intervalle_pari_juste <- (Intervalle_pari_juste[1,] <= Valeurs_nouveaux_points)*(Intervalle_pari_juste[2,] >= Valeurs_nouveaux_points)
    #print(Intervalle_pari_juste[1,] <= Valeurs_nouveaux_points)
    #print(Intervalle_pari_juste[2,] >= Valeurs_nouveaux_points)
    #print(Appartient_intervalle_pari_juste)
    Frequence_appartient_intervalle_pari_juste <- sum(Appartient_intervalle_pari_juste)/CARDINAL_ENSEMBLE_TEST    













    
### NOUVEAU : accomode la tendance
Moyenne_marginale <- 0
    
if(NOMBRE_FONCTIONS_TENDANCE>0)
{
	QR <- qr.Q(qr(tendance),complete=TRUE) #les NOMBRE_FONCTIONS_TENDANCE premieres colonnes de QR contiennent P, les dernieres W
	injectionOrthogonalTendance <- QR[,(NOMBRE_FONCTIONS_TENDANCE+1):ncol(QR)] #W
        injectionTendance <- QR[,1:NOMBRE_FONCTIONS_TENDANCE] #P


        ##Il faut redefinir tous les objets relatifs aux donnes connues



        TENDANCE_NOUVEAUX <-matrix(0,nrow=CARDINAL_ENSEMBLE_TEST, ncol=length(FONCTIONS))
	for(i in 1:NOMBRE_FONCTIONS_TENDANCE)
	{
		TENDANCE_NOUVEAUX[,i] <- apply(Nouveaux_points,1,FONCTIONS[[i]])
	}
        TENDANCE_NOUVEAUX_CONNUS <- TENDANCE_NOUVEAUX %*% solve(t(injectionTendance) %*% tendance, t(injectionTendance)) ## Attention, c'est H_0_0 %*% (P^T H)^(-1) %*% P^T = H_0_point %*% P^T
        





        ## Il faut redefinir la variance marginale des nouveaux points...
        #MATRICE_NOUVEAUX_MODE_POSTERIOR <- MATRICE_NOUVEAUX_MODE_POSTERIOR + TENDANCE_NOUVEAUX_CONNUS %*% MATRICE_MODE_POSTERIOR %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR - t(CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR) %*% t(TENDANCE_NOUVEAUX_CONNUS)

	## ... ou alors on peut être plus malin et se rappeler que en fin de compte seule la diagonale de cette matrice nous intéresse, ce qui permet de remplacer l'ancienne MATRICE_NOUVEAUX_MODE_POSTERIOR par la matrice remplie de 1.
        MATRICE_NOUVEAUX_MODE_POSTERIOR <- 1 + TENDANCE_NOUVEAUX_CONNUS %*% MATRICE_MODE_POSTERIOR %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR - t(CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR) %*% t(TENDANCE_NOUVEAUX_CONNUS)


        Moyenne_marginale <- t( TENDANCE_NOUVEAUX_CONNUS %*% y_connus )
        #print(Moyenne_marginale)
        
        ## Attention, on se met a redefinir des matrices. Leur ancienne signification est perdue a partir d'ici.


        CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR <-  t(injectionOrthogonalTendance) %*% (MATRICE_MODE_POSTERIOR %*% t(TENDANCE_NOUVEAUX_CONNUS) - CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR)

  
        y_connus <- - t(injectionOrthogonalTendance) %*% y_connus ## Attention, on se base sur MOINS W^T y
        



        MATRICE_MODE_POSTERIOR <- t(injectionOrthogonalTendance) %*% MATRICE_MODE_POSTERIOR %*% injectionOrthogonalTendance
        MATRICE_MODE_POSTERIOR_INVERSE <- solve(MATRICE_MODE_POSTERIOR)        
}



    

    
    #Aux_mode_posterior <-  solve(MATRICE_MODE_POSTERIOR, CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR)
    #Prediction_mode_posterior <- t(y_connus) %*% Aux_mode_posterior
    Prediction_mode_posterior <- t(y_connus) %*% MATRICE_MODE_POSTERIOR_INVERSE %*% CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR + Moyenne_marginale



## Ecarts au carre entre moyennes conditionnelles supposant pour les longueurs de correlation les valeurs du MAP et la vraie moyenne conditionnelle

    Ecart_prediction_mode_posterior_moyenne_conditionnelle_juste <- sum( (Prediction_mode_posterior - Moyenne_conditionnelle_juste)*(Prediction_mode_posterior - Moyenne_conditionnelle_juste) )




    ## Ecarts au carre entre moyennes conditionnelles supposant pour les longueurs de correlation les valeurs du MAP et la vraie valeur 

    Ecart_prediction_mode_posterior <- sum( (Prediction_mode_posterior - Valeurs_nouveaux_points)*(Prediction_mode_posterior - Valeurs_nouveaux_points) )

    



    ## Calcul de l'intervalle de pari correspondant aux longueurs de correlation donnees par le MAP. Attention, chaque intervalle de pari est une COLONNE
    Estimation_sigma_carre_mode_posterior<-  c( t(y_connus) %*% MATRICE_MODE_POSTERIOR_INVERSE %*% y_connus / length(y_connus) )
    print(paste("Estimation de sigma^2 (MAP) : ",Estimation_sigma_carre_mode_posterior))
    Variance_conditionnelle_mode_posterior<-  Estimation_sigma_carre_mode_posterior * ( MATRICE_NOUVEAUX_MODE_POSTERIOR - t(CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR) %*% MATRICE_MODE_POSTERIOR_INVERSE %*% CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR )
    Intervalle_pari_mode_posterior <- rbind(Prediction_mode_posterior - QNORM * sqrt(diag(Variance_conditionnelle_mode_posterior)) , Prediction_mode_posterior + QNORM * sqrt(diag(Variance_conditionnelle_mode_posterior)))
    Longueur_intervalle_pari_mode_posterior<- Intervalle_pari_mode_posterior[2,] - Intervalle_pari_mode_posterior[1,]
    Longueur_moyenne_intervalle_pari_mode_posterior<- mean(Longueur_intervalle_pari_mode_posterior)
    Appartient_intervalle_pari_mode_posterior<- (Intervalle_pari_mode_posterior[1,] <= Valeurs_nouveaux_points)*(Intervalle_pari_mode_posterior[2,] >= Valeurs_nouveaux_points)
    Frequence_appartient_intervalle_pari_mode_posterior<- sum(Appartient_intervalle_pari_mode_posterior)/CARDINAL_ENSEMBLE_TEST

    #print(Intervalle_pari_juste - rbind(Valeurs_nouveaux_points, Valeurs_nouveaux_points) )





    Longueurs_moyennes_intervalles_pari <- c(Longueur_moyenne_intervalle_pari_juste, 0, Longueur_moyenne_intervalle_pari_mode_posterior, 0, 0)
    
    Frequences_appartient_intervalle_pari <- c(Frequence_appartient_intervalle_pari_juste, 0, Frequence_appartient_intervalle_pari_mode_posterior, 0, 0)




Ecart_prediction_moyenne_conditionnelle_juste<- c(0, Ecart_prediction_mode_posterior_moyenne_conditionnelle_juste, 0)
write.matrix(Ecart_prediction_moyenne_conditionnelle_juste, "ecarts_predictions_moyenne_cond_juste.txt",sep = "\t")

    

Ecart_prediction <- c(0, Ecart_prediction_mode_posterior, 0)

    write.matrix(Ecart_prediction, "ecarts_predictions.txt",sep = "\t")


#        fonctionDensitePredictive <- function(t, numero_nouveau_point)    
#        {
#            mean( apply(Points_LC_parnorm[numero_nouveau_point,,],1,dnormParametresDevant, t=t) )
#        }

#    #xxx <- 0.05*-100:100
#    xxx <- Moyenne_conditionnelle_juste[1] + 0.01*-100:100
#    yyy <- apply(matrix(xxx),1,fonctionDensitePredictive,1)
#    #YYY <- apply(matrix(xxx),1,fonctionRepartitionPredictive,1)
#    zzz <- dnorm(xxx, mean= mean(Points_LC_parnorm[1,,1]), sd= sqrt( mean(Points_LC_parnorm[1,,2]) + var(Points_LC_parnorm[1,,1]) ) )
#    #ZZZ <- pnorm(xxx, mean= mean(Points_LC_parnorm[1,,1]), sd= sqrt( mean(Points_LC_parnorm[1,,2]) + var(Points_LC_parnorm[1,,1]) ) )
#    #ttt <- dnorm(xxx, mean= mean(Points_LC_parnorm[1,,1]), sd= sqrt( var(Points_LC_parnorm[1,,1]) ) )
##    plot(xxx, dnorm(xxx,mean=Moyenne_conditionnelle_juste[1], sd= sqrt(Variance_conditionnelle_juste[1,1])),type="l", xlab="Unobserved value", ylab="Predictive probability density", main="", lwd=2)
#    #lines(xxx, dnorm(xxx,mean=Prediction_argmax_vraisemblance[1], sd=sqrt(Variance_conditionnelle_argmax_vraisemblance[1,1])),col=3)
#    lines(xxx, dnorm(xxx,mean=Prediction_mode_posterior[1], sd=sqrt(Variance_conditionnelle_mode_posterior[1,1])),col=2, lty=2, lwd=2)
#    lines(xxx, yyy, col=4, lty=3, lwd=2)
#    lines(xxx, zzz, col=4, lty=3, lwd=1)
#    #abline(v=Intervalle_pari_full_bayesien_exact[1,1], col=4, lwd=2, lty=3)
#    #abline(v=Intervalle_pari_full_bayesien_exact[1,2], col=4, lwd=2, lty=3)
#    #abline(v=Intervalle_pari_full_bayesien[1,1], col=4, lty=3)
#    #abline(v=Intervalle_pari_full_bayesien[1,2], col=4, lty=3)

#    legende <- c("True","MAP","full-Bayes.", "approx. full-Bayes.")
#    col <- c(1,2,4,4)
#    lty <- c(1,2,3,3)
#    lwd <- c(2,2,2,1)

#    legend("topright",legende, col=col, lty=lty, lwd=lwd)

#    dev.copy(pdf,paste0("Graphes/dens_predictive_",numero_germe_aleatoire,".pdf"))
#    dev.off()
    
#    #lines(xxx, ttt, col=5)
#    #abline(v=mean(Points_LC_parnorm[1,,1], col=4))

#    #plot(xxx, ttt, col=5, type="l")
#    #lines(xxx, zzz, col=4)
    
#                                        #plot(xxx, (zzz-yyy)/yyy, type="l", col=4)
#    #plot(xxx, ZZZ-YYY, type="l", col=4)
#    #print(var(Points_LC_parnorm[1,,1]) / mean(Points_LC_parnorm[1,,2]) )

#    #print(Points_LC_parnorm[1,,2])
}

