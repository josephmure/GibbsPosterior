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
	NOYAU_ARGMAX_VRAISEMBLANCE <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION_ARGMAX_VRAISEMBLANCE, regularite = REGULARITE ) # anisotrope geometrique en fait !
	NOYAU_MODE_POSTERIOR <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION_MODE_POSTERIOR, regularite = REGULARITE ) # anisotrope geometrique en fait !
} else if (TYPE_NOYAU_MATERN == "tensorise") 
{
	NOYAU <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE )
	NOYAU_ARGMAX_VRAISEMBLANCE <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION_ARGMAX_VRAISEMBLANCE, regularite = REGULARITE ) # anisotrope geometrique en fait !
	NOYAU_MODE_POSTERIOR <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION_MODE_POSTERIOR, regularite = REGULARITE ) # anisotrope geometrique en fait !
} else 
{
	print("Erreur dans compareResultats.r : le noyau doit-il etre anis. geometrique ou tensorise ?")
	NOYAU <- NULL
#	NOYAU_ARGMAX_VRAISEMBLANCE <- NULL
#	NOYAU_MODE_POSTERIOR <- NULL
}

if(!is.null(NOYAU)) ## Si on a bien precise la nature anisotrope geometrique ou tensorisee du noyau
{
    ## Matrice de correlation juste : correspondant aux vraies longueurs de correlation
    MATRICE_JUSTE <- creeMatriceCovariance(x1=x_connus, x2=x_connus, noyau= NOYAU)
    MATRICE_JUSTE_INVERSE <- solve(MATRICE_JUSTE) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.



    
    ## Matrice de correlation  correspondant aux longueurs de correlation estimees par MLE : maximum de vraisemblance et MAP : max a posteriori respectivement
    MATRICE_ARGMAX_VRAISEMBLANCE <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= NOYAU_ARGMAX_VRAISEMBLANCE)
    MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE <- solve(MATRICE_ARGMAX_VRAISEMBLANCE) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.
    MATRICE_MODE_POSTERIOR <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= NOYAU_MODE_POSTERIOR)
    MATRICE_MODE_POSTERIOR_INVERSE <- solve(MATRICE_MODE_POSTERIOR) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.



    
    ## Recherche de la distance entre matrice de correl correspondant aux vraies longueurs de correlation et les matrices de correl obtenues avec MLE et MAP

    ## Transfo de Fisher appliquee aux matrices de correlation
	MATRICE_JUSTE_FISHER <- 1/2 * log( (1 +  MATRICE_JUSTE) / (1 - MATRICE_JUSTE) ) # transfo de Fisher
	diag(MATRICE_JUSTE_FISHER) <- 0
#print(MATRICE_JUSTE_FISHER)

	MATRICE_ARGMAX_VRAISEMBLANCE_FISHER <- 1/2 * log( (1 +  MATRICE_ARGMAX_VRAISEMBLANCE) / (1 - MATRICE_ARGMAX_VRAISEMBLANCE) ) # transfo de Fisher
	diag(MATRICE_ARGMAX_VRAISEMBLANCE_FISHER) <- 0
#print(MATRICE_ARGMAX_VRAISEMBLANCE_FISHER)


	MATRICE_MODE_POSTERIOR_FISHER <- 1/2 * log( (1 +  MATRICE_MODE_POSTERIOR) / (1 - MATRICE_MODE_POSTERIOR) ) # transfo de Fisher
	diag(MATRICE_MODE_POSTERIOR_FISHER) <- 0
#print(MATRICE_MODE_POSTERIOR_FISHER)

    ## Distances "simple" et "de Fisher" entre matrices obtenues par MLE et MAP et la vraie matrice
	Distance_argmax_vraisemblance <- c( norm(MATRICE_ARGMAX_VRAISEMBLANCE - MATRICE_JUSTE, "F"), norm(MATRICE_ARGMAX_VRAISEMBLANCE_FISHER - MATRICE_JUSTE_FISHER, "F") ) # norme de Frobenius (ie euclidienne)
#print(Distance_argmax_vraisemblance)
	Distance_mode_posterior <- c( norm(MATRICE_MODE_POSTERIOR - MATRICE_JUSTE, "F"), norm(MATRICE_MODE_POSTERIOR_FISHER - MATRICE_JUSTE_FISHER, "F") ) # norme de Frobenius (ie euclidienne)
#print(Distance_mode_posterior)

    write.matrix(Distance_argmax_vraisemblance,"argmax_vraisemblance_integree_distance.txt",sep = "\t")
    write.matrix(Distance_mode_posterior,"mode_posterior_distance.txt",sep = "\t")



    

# Ecart quadratique entre la prediction (moyenne conditionnelle sachant la valeur estimee des longueurs de correlation)
# et la moyenne conditionnelle sachant la vraie valeur des longueurs de correlation


## Generation de l'ensemble test
    Nouveaux_points <- generePointsConnus(nombre_points = CARDINAL_ENSEMBLE_TEST, minima = rep(0,NOMBRE_DIMENSIONS), maxima = rep(1,NOMBRE_DIMENSIONS), lhs=FALSE)
    Tous_les_points <- rbind(x_connus, Nouveaux_points)


    ## Matrice de variance marginale correspondant aux Nouveaux_points
    MATRICE_NOUVEAUX_JUSTE <- creeMatriceCovariance(x1= Nouveaux_points, x2= Nouveaux_points, noyau= NOYAU)
    ## INUTILE : tout ce qui compte, ce sont les éléments diagonaux, qui sont tous égaux à 1
	#MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- creeMatriceCovariance(x1= Nouveaux_points, x2= Nouveaux_points, noyau= NOYAU_ARGMAX_VRAISEMBLANCE)
    	#MATRICE_NOUVEAUX_MODE_POSTERIOR <- creeMatriceCovariance(x1= Nouveaux_points, x2= Nouveaux_points, noyau= NOYAU_MODE_POSTERIOR)
	MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- 1
	MATRICE_NOUVEAUX_MODE_POSTERIOR <- 1
    
# Attention, non seulement la technique suivante est plus efficace que les commentaires, mais elle rend aussi les denominations plus justes. Par consequent,
# CORREL_CONNUS_NOUVEAUX_... actuel = t( CORREL_CONNUS_NOUVEAUX_...  ancien ) et Moyenne_conditionnelle_juste, ainsi que les predictions, sont des vecteurs lignes.
#MATRICE_JUSTE_TOTALE <- creeMatriceCovariance(x1= Tous_les_points, x2= Tous_les_points, noyau= NOYAU)
#CORREL_CONNUS_NOUVEAUX_JUSTE <- MATRICE_JUSTE_TOTALE[(nrow(x_connus)+1):(nrow(x_connus)+CARDINAL_ENSEMBLE_TEST),1:nrow(x_connus)]
CORREL_CONNUS_NOUVEAUX_JUSTE <- creeMatriceCovariance(x1= x_connus, x2= Nouveaux_points, noyau= NOYAU)
#MATRICE_ARGMAX_VRAISEMBLANCE_TOTALE <- creeMatriceCovariance(x1= Tous_les_points, x2= Tous_les_points, noyau= NOYAU_ARGMAX_VRAISEMBLANCE)
CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- creeMatriceCovariance(x1= x_connus, x2= Nouveaux_points, noyau= NOYAU_ARGMAX_VRAISEMBLANCE)
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
        #MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE + TENDANCE_NOUVEAUX_CONNUS %*% MATRICE_ARGMAX_VRAISEMBLANCE %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE - t(CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE) %*% t(TENDANCE_NOUVEAUX_CONNUS)

	## ... ou alors on peut être plus malin et se rappeler que en fin de compte seule la diagonale de cette matrice nous intéresse, ce qui permet de remplacer l'ancienne MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE par la matrice remplie de 1.
	MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- 1 + TENDANCE_NOUVEAUX_CONNUS %*% MATRICE_ARGMAX_VRAISEMBLANCE %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE - t(CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE) %*% t(TENDANCE_NOUVEAUX_CONNUS)


        ## Il faut redefinir la variance marginale des nouveaux points...
        #MATRICE_NOUVEAUX_MODE_POSTERIOR <- MATRICE_NOUVEAUX_MODE_POSTERIOR + TENDANCE_NOUVEAUX_CONNUS %*% MATRICE_MODE_POSTERIOR %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR - t(CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR) %*% t(TENDANCE_NOUVEAUX_CONNUS)

	## ... ou alors on peut être plus malin et se rappeler que en fin de compte seule la diagonale de cette matrice nous intéresse, ce qui permet de remplacer l'ancienne MATRICE_NOUVEAUX_MODE_POSTERIOR par la matrice remplie de 1.
        MATRICE_NOUVEAUX_MODE_POSTERIOR <- 1 + TENDANCE_NOUVEAUX_CONNUS %*% MATRICE_MODE_POSTERIOR %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR - t(CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR) %*% t(TENDANCE_NOUVEAUX_CONNUS)


        Moyenne_marginale <- t( TENDANCE_NOUVEAUX_CONNUS %*% y_connus )
        #print(Moyenne_marginale)
        
        ## Attention, on se met a redefinir des matrices. Leur ancienne signification est perdue a partir d'ici.

        CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- t(injectionOrthogonalTendance) %*% (MATRICE_ARGMAX_VRAISEMBLANCE %*% t(TENDANCE_NOUVEAUX_CONNUS) - CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE)
        CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR <-  t(injectionOrthogonalTendance) %*% (MATRICE_MODE_POSTERIOR %*% t(TENDANCE_NOUVEAUX_CONNUS) - CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR)

  
        y_connus <- - t(injectionOrthogonalTendance) %*% y_connus ## Attention, on se base sur MOINS W^T y
        
        MATRICE_ARGMAX_VRAISEMBLANCE <- t(injectionOrthogonalTendance) %*% MATRICE_ARGMAX_VRAISEMBLANCE %*% injectionOrthogonalTendance
        MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE <- solve(MATRICE_ARGMAX_VRAISEMBLANCE)

        MATRICE_MODE_POSTERIOR <- t(injectionOrthogonalTendance) %*% MATRICE_MODE_POSTERIOR %*% injectionOrthogonalTendance
        MATRICE_MODE_POSTERIOR_INVERSE <- solve(MATRICE_MODE_POSTERIOR)        
}



    
    #Aux_argmax_vraisemblance <-  solve(MATRICE_ARGMAX_VRAISEMBLANCE, CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE)
    #Prediction_argmax_vraisemblance <- t(y_connus) %*% Aux_argmax_vraisemblance
    Prediction_argmax_vraisemblance <- t(y_connus) %*% MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE %*% CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE + Moyenne_marginale
    
    #Aux_mode_posterior <-  solve(MATRICE_MODE_POSTERIOR, CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR)
    #Prediction_mode_posterior <- t(y_connus) %*% Aux_mode_posterior
    Prediction_mode_posterior <- t(y_connus) %*% MATRICE_MODE_POSTERIOR_INVERSE %*% CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR + Moyenne_marginale



## Ecarts au carre entre moyennes conditionnelles supposant pour les longueurs de correlation les valeurs du MLE et du MAP respectivement et la vraie moyenne conditionnelle
Ecart_prediction_argmax_vraisemblance_moyenne_conditionnelle_juste <- sum( (Prediction_argmax_vraisemblance - Moyenne_conditionnelle_juste)*(Prediction_argmax_vraisemblance - Moyenne_conditionnelle_juste) )
    Ecart_prediction_mode_posterior_moyenne_conditionnelle_juste <- sum( (Prediction_mode_posterior - Moyenne_conditionnelle_juste)*(Prediction_mode_posterior - Moyenne_conditionnelle_juste) )




    ## Ecarts au carre entre moyennes conditionnelles supposant pour les longueurs de correlation les valeurs du MLE et du MAP respectivement et la vraie valeur 
    Ecart_prediction_argmax_vraisemblance <- sum( (Prediction_argmax_vraisemblance - Valeurs_nouveaux_points)*(Prediction_argmax_vraisemblance - Valeurs_nouveaux_points) )
    Ecart_prediction_mode_posterior <- sum( (Prediction_mode_posterior - Valeurs_nouveaux_points)*(Prediction_mode_posterior - Valeurs_nouveaux_points) )

    


                                        
    ## Calcul de l'intervalle de pari correspondant aux longueurs de correlation donnees par le MLE. Attention, chaque intervalle de pari est une COLONNE
    Estimation_sigma_carre_argmax_vraisemblance <-  c( t(y_connus) %*% MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE %*% y_connus / length(y_connus) )
    print(paste("Estimation de sigma^2 (MLE) : ",Estimation_sigma_carre_argmax_vraisemblance))
    Variance_conditionnelle_argmax_vraisemblance <-  Estimation_sigma_carre_argmax_vraisemblance * ( MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE - t(CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE) %*% MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE %*% CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE )

#print(diag( Estimation_sigma_carre_argmax_vraisemblance * ( MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE - t(CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE) %*% MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE %*% CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE ) ) -
#diag( Estimation_sigma_carre_argmax_vraisemblance * ( MATRICE_NOUVEAUX_ARGMAX_TRONQUE - t(CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE) %*% MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE %*% CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE ) ) )

    Intervalle_pari_argmax_vraisemblance <- rbind(Prediction_argmax_vraisemblance - QNORM * sqrt(diag(Variance_conditionnelle_argmax_vraisemblance)), Prediction_argmax_vraisemblance + QNORM * sqrt(diag(Variance_conditionnelle_argmax_vraisemblance)))
    Longueur_intervalle_pari_argmax_vraisemblance <- Intervalle_pari_argmax_vraisemblance[2,] - Intervalle_pari_argmax_vraisemblance[1,]
    Longueur_moyenne_intervalle_pari_argmax_vraisemblance<- mean(Longueur_intervalle_pari_argmax_vraisemblance)
    Appartient_intervalle_pari_argmax_vraisemblance <- (Intervalle_pari_argmax_vraisemblance[1,] <= Valeurs_nouveaux_points)*(Intervalle_pari_argmax_vraisemblance[2,] >= Valeurs_nouveaux_points)
    Frequence_appartient_intervalle_pari_argmax_vraisemblance <- sum(Appartient_intervalle_pari_argmax_vraisemblance)/CARDINAL_ENSEMBLE_TEST

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



    ## Intervalle de pari full-bayesien

    if(TYPE_NOYAU_MATERN == "geometrique") 
    {
        print("Debut inversions")

        if(NOMBRE_FONCTIONS_TENDANCE>0)
            {
        Tableau_matrices_correlation <- apply(ECHANTILLON_POSTERIOR,1,vecteurMatriceCorrelationAnisGeom, regularite=REGULARITE, x_connus=x_connus)
        Tableau_inverses_matrices_correlation <- apply(Tableau_matrices_correlation,2,vecteurInverseMatriceRestreinte,nb_lignes=TAILLE_PLAN_XP, injectionOrthogonalTendance=injectionOrthogonalTendance)
            } else
            {
                ## Tableau dont chaque COLONNE est la TOTALITE de l'inverse de la matrice de correlation correspondant a un jeu de longueurs de correlation
                Tableau_inverses_matrices_correlation <- apply(ECHANTILLON_POSTERIOR,1,vecteurInverseMatriceCorrelationAnisGeom, regularite=REGULARITE, x_connus=x_connus)
            }
        
        print("Fin inversions")
        ## Points_LC_parnrom est un Array en dimension 3
        ## dimension 1 : nouveaux points
        ## dimension 2 : jeux de longueurs de correlation obtenus par echantillonnage
        ## dimension 3 : moyenne au nouveau point sachant le jeu de longueurs de correlation ; variance etc.
        Points_LC_parnorm <- array(dim=c(CARDINAL_ENSEMBLE_TEST,nrow(ECHANTILLON_POSTERIOR), 2) )

        LC_correlConnusNouveaux <- array(dim=c(nrow(ECHANTILLON_POSTERIOR), TAILLE_PLAN_XP, CARDINAL_ENSEMBLE_TEST))

        print("Début corrélations")
        
        for(i in 1:nrow(ECHANTILLON_POSTERIOR) )
        {
            Noyau <- creeMaternIsotrope(variance=1, regularite=REGULARITE, longueur=ECHANTILLON_POSTERIOR[i,])
            LC_correlConnusNouveaux[i,,] <- creeMatriceCovariance(x1=x_connus, x2=Nouveaux_points, noyau=Noyau)
        }
        print("Fin corrélations")

    } else if(TYPE_NOYAU_MATERN == "tensorise")
    {
        if(NOMBRE_FONCTIONS_TENDANCE>0)
            {
        Tableau_matrices_correlation <- apply(ECHANTILLON_POSTERIOR,1,vecteurMatriceCorrelationTens, regularite=REGULARITE, x_connus=x_connus)
        Tableau_inverses_matrices_correlation <- apply(Tableau_matrices_correlation,2,vecteurInverseMatriceRestreinte,nb_lignes=TAILLE_PLAN_XP, injectionOrthogonalTendance=injectionOrthogonalTendance)
            } else
            {
                ## Tableau dont chaque COLONNE est la TOTALITE de l'inverse de la matrice de correlation correspondant a un jeu de longueurs de correlation
                Tableau_inverses_matrices_correlation <- apply(ECHANTILLON_POSTERIOR,1,vecteurInverseMatriceCorrelationTens, regularite=REGULARITE, x_connus=x_connus)
            }
        
        print("Fin inversions")
        ## Points_LC_parnrom est un Array en dimension 3
        ## dimension 1 : nouveaux points
        ## dimension 2 : jeux de longueurs de correlation obtenus par echantillonnage
        ## dimension 3 : moyenne au nouveau point sachant le jeu de longueurs de correlation ; variance etc.
        Points_LC_parnorm <- array(dim=c(CARDINAL_ENSEMBLE_TEST,nrow(ECHANTILLON_POSTERIOR), 2) )

        LC_correlConnusNouveaux <- array(dim=c(nrow(ECHANTILLON_POSTERIOR), TAILLE_PLAN_XP, CARDINAL_ENSEMBLE_TEST))
        ##LC_correlNouveauxNouveaux <- array(dim=c(nrow(ECHANTILLON_POSTERIOR), CARDINAL_ENSEMBLE_TEST, CARDINAL_ENSEMBLE_TEST)) ## NEW

        print("Début corrélations")
        
        for(i in 1:nrow(ECHANTILLON_POSTERIOR) )
        {
            Noyau <- creeMaternTensorise(variance=1, regularite=REGULARITE, longueur=ECHANTILLON_POSTERIOR[i,])
            LC_correlConnusNouveaux[i,,] <- creeMatriceCovariance(x1=x_connus, x2=Nouveaux_points, noyau=Noyau)
	    ##LC_correlNouveauxNouveaux[i,,] <- creeMatriceCovariance(x1=Nouveaux_points,x2=Nouveaux_points, noyau=Noyau) ## NEW
        }
        print("Fin corrélations")
    }
    
    

    if(NOMBRE_FONCTIONS_TENDANCE>0)
    {
        
#### A SIMPLIFIER        for(numero_nouveau_point in 1:CARDINAL_ENSEMBLE_TEST)
#### A SIMPLIFIER        {
            #print("Début corrélations")
            ## Tableau dont chaque COLONNE est le vecteur de correlation entre points du planXP et le nouveau point numer_nouveau_point correspondant a un jeu de longueurs de correlation
            #Tableau_correl_connus_nouveau <- apply(ECHANTILLON_POSTERIOR,1,vecteurCorrelConnusNouveauAnisGeom, regularite=REGULARITE, x_connus=x_connus, point_nouveau=Nouveaux_points[numero_nouveau_point,])
            #print("Fin corrélations")

            for(i in 1:nrow(ECHANTILLON_POSTERIOR))
            {
                matrice_correlation_inverse <-  matrix(Tableau_inverses_matrices_correlation[,i],nrow=TAILLE_PLAN_XP - NOMBRE_FONCTIONS_TENDANCE) ## rend a la matrice de correlation inverse sa structure. Attention, il s'agit en realite de la matrice restreinte W^T Sigma^(-1) W !!!

                matrice_correlation <- matrix(Tableau_matrices_correlation[,i],nrow=TAILLE_PLAN_XP)
                
                ##correl_connus_nouveau <- Tableau_correl_connus_nouveau[,i]            
#### A SIMPLIFIER                correl_connus_nouveau <- LC_correlConnusNouveaux[i,,numero_nouveau_point]
         	correl_connus_nouveau <- LC_correlConnusNouveaux[i,,]
		##correl_nouveau_nouveau <- LC_correlNouveauxNouveaux[i,,] ## NEW
       
### A SIMPLIFIER                variance_marginale <- 1 + t(TENDANCE_NOUVEAUX_CONNUS[numero_nouveau_point,]) %*% matrice_correlation %*% TENDANCE_NOUVEAUX_CONNUS[numero_nouveau_point,] - t(TENDANCE_NOUVEAUX_CONNUS[numero_nouveau_point,]) %*% correl_connus_nouveau - t(correl_connus_nouveau) %*% TENDANCE_NOUVEAUX_CONNUS[numero_nouveau_point,]

		## Attention, normalement, variance_marginale devrait etre egale à "matrice de variance (ou plutôt de corrélation) des nouveaux points + TENDANCE_NOUVEAUX_CONNUS %*% ...
		## Mais en pratique, on utilisera seulement les termes diagonaux de cette matrice, or une matrice de corrélation a tous ses coefficients diagonaux égaux à 1, donc on peut simplifier le calcul en remplaçant la matrice de corrélation par la matrice remplie de 1.
		variance_marginale <- 1 + TENDANCE_NOUVEAUX_CONNUS %*% matrice_correlation %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% correl_connus_nouveau - t(correl_connus_nouveau) %*% t(TENDANCE_NOUVEAUX_CONNUS)

                correl_connus_nouveau <- t(injectionOrthogonalTendance) %*% (matrice_correlation %*% t(TENDANCE_NOUVEAUX_CONNUS) - correl_connus_nouveau)
                
#### A SIMPLIFIER                moyenne_conditionnelle <-t(correl_connus_nouveau)%*% matrice_correlation_inverse %*% y_connus + Moyenne_marginale[numero_nouveau_point]
		moyenne_conditionnelle <-c( t(y_connus) %*% matrice_correlation_inverse %*% correl_connus_nouveau + Moyenne_marginale )

                variance_conditionnelle <- c( t(y_connus) %*% matrice_correlation_inverse %*% y_connus ) / TAILLE_PLAN_XP * diag(variance_marginale - t(correl_connus_nouveau) %*% matrice_correlation_inverse %*% correl_connus_nouveau)
                #print(paste("variance_conditionnelle=",variance_conditionnelle))
#### A SIMPLIFIER                Points_LC_parnorm[numero_nouveau_point,i,] <- c(moyenne_conditionnelle,variance_conditionnelle)
		Points_LC_parnorm[,i,] <- cbind(moyenne_conditionnelle,variance_conditionnelle)
#		print(max(abs(moyenne_conditionnelle - Moyenne_conditionnelle_juste)))
#print("Points_LC_parnorm[,i,1] = ")
#print(max(abs(Points_LC_parnorm[,i,1] - Moyenne_conditionnelle_juste)))
#print(moyenne_conditionnelle - Points_LC_parnorm[,i,1])

                #print(Points_LC_parnorm[numero_nouveau_point,i,])
            }
                #print(paste("dim_variance_conditionnelle=",dim(variance_conditionnelle)))

            #Points_LC_parnorm[,,2] <- sqrt(Points_LC_parnorm[,,2]) # pour que [,,2] contienne les ecarts-types conditionnels
            #print(Points_LC_parnorm[,,2])
#### A SIMPLIFIER        }
    }else
    {
        for(numero_nouveau_point in 1:CARDINAL_ENSEMBLE_TEST)
        {
            #print("Début corrélations")
            ## Tableau dont chaque COLONNE est le vecteur de correlation entre points du planXP et le nouveau point numer_nouveau_point correspondant a un jeu de longueurs de correlation
            #Tableau_correl_connus_nouveau <- apply(ECHANTILLON_POSTERIOR,1,vecteurCorrelConnusNouveauAnisGeom, regularite=REGULARITE, x_connus=x_connus, point_nouveau=Nouveaux_points[numero_nouveau_point,])
            #print("Fin corrélations")

            for(i in 1:nrow(ECHANTILLON_POSTERIOR))
            {
                matrice_correlation_inverse <-  matrix(Tableau_inverses_matrices_correlation[,i],nrow=TAILLE_PLAN_XP) ## rend a la matrice de correlation inverse sa structure
                                        #correl_connus_nouveau <- Tableau_correl_connus_nouveau[,i]
                correl_connus_nouveau <- LC_correlConnusNouveaux[i,,numero_nouveau_point]
                moyenne_conditionnelle <-t(correl_connus_nouveau)%*% matrice_correlation_inverse %*% y_connus
                variance_conditionnelle <- c( t(y_connus) %*% matrice_correlation_inverse %*% y_connus ) / TAILLE_PLAN_XP * (1 - t(correl_connus_nouveau) %*% matrice_correlation_inverse %*% correl_connus_nouveau)
                #print(paste("variance_conditionnelle=",variance_conditionnelle))
                Points_LC_parnorm[numero_nouveau_point,i,] <- c(moyenne_conditionnelle,variance_conditionnelle)
                #print(Points_LC_parnorm[numero_nouveau_point,i,])
            }
                #print(paste("dim_variance_conditionnelle=",dim(variance_conditionnelle)))

            #Points_LC_parnorm[,,2] <- sqrt(Points_LC_parnorm[,,2]) # pour que [,,2] contienne les ecarts-types conditionnels
            #print(Points_LC_parnorm[,,2])
        }
    }
    
   
        fonctionRepartitionPredictive <- function(t, numero_nouveau_point)    
        {
            mean( apply(Points_LC_parnorm[numero_nouveau_point,,],1,pnormParametresDevant, t=t) )
        }


    
    ## ATTENTION ! Les intervalles de pari sont ici les lignes. 
    Intervalle_pari_full_bayesien <- matrix(nrow=CARDINAL_ENSEMBLE_TEST,ncol=2)
    Intervalle_pari_full_bayesien_exact<- matrix(nrow=CARDINAL_ENSEMBLE_TEST,ncol=2)
#print("Moyenne des moyennes :")
    Moyennes_conditionnelles_bayesiennes <- apply(Points_LC_parnorm[,,1],1,mean)
#print(max(abs(Moyennes_conditionnelles_bayesiennes - Moyenne_conditionnelle_juste)))
    Moyennes_des_variances <- apply(Points_LC_parnorm[,,2],1,mean)
    Variances_des_moyennes <- apply(Points_LC_parnorm[,,1],1,var)
    Variances_conditionnelles_bayesiennes <- Moyennes_des_variances + Variances_des_moyennes

    Intervalle_pari_full_bayesien <- cbind(Moyennes_conditionnelles_bayesiennes - QNORM * sqrt(Variances_conditionnelles_bayesiennes) , Moyennes_conditionnelles_bayesiennes +  QNORM * sqrt(Variances_conditionnelles_bayesiennes) )
    
        for(numero_nouveau_point in 1:CARDINAL_ENSEMBLE_TEST)
        {

#            Intervalle_pari_full_bayesien[numero_nouveau_point,] <- c( mean(Points_LC_parnorm[numero_nouveau_point,,1]) - QNORM * sqrt( mean(Points_LC_parnorm[numero_nouveau_point,,2]) + var(Points_LC_parnorm[numero_nouveau_point,,1]) ),  mean(Points_LC_parnorm[numero_nouveau_point,,1]) + QNORM * sqrt( mean(Points_LC_parnorm[numero_nouveau_point,,2]) + var(Points_LC_parnorm[numero_nouveau_point,,1]) ) )
            
            
            Intervalle_pari_full_bayesien_exact[numero_nouveau_point,] <- dichotomieIntervallePari(fonctionRepartitionPredictive,numero_nouveau_point= numero_nouveau_point,  probaInf=(1-POURCENTAGE_INTERVALLE_CONFIANCE/100)/2, probaSup=1 - (1-POURCENTAGE_INTERVALLE_CONFIANCE/100)/2, pas = 1)
#            print(c(Intervalle_pari_full_bayesien[numero_nouveau_point,1],Valeurs_nouveaux_points[numero_nouveau_point],Intervalle_pari_full_bayesien[numero_nouveau_point,2]))
#            print(c(fonctionRepartitionPredictive(Intervalle_pari_full_bayesien[numero_nouveau_point,1],numero_nouveau_point=numero_nouveau_point),fonctionRepartitionPredictive(Intervalle_pari_full_bayesien[numero_nouveau_point,2],numero_nouveau_point=numero_nouveau_point)))
        }

        ## ATTENTION ! Les intervalles de pari sont ici les lignes. 
    Appartient_intervalle_pari_full_bayesien<- (Intervalle_pari_full_bayesien[,1] <= Valeurs_nouveaux_points)*(Intervalle_pari_full_bayesien[,2] >= Valeurs_nouveaux_points)
    Longueur_intervalle_pari_full_bayesien<- Intervalle_pari_full_bayesien[,2] - Intervalle_pari_full_bayesien[,1]
    Frequence_appartient_intervalle_pari_full_bayesien<- sum(Appartient_intervalle_pari_full_bayesien)/CARDINAL_ENSEMBLE_TEST       
    Longueur_moyenne_intervalle_pari_full_bayesien <- mean(Longueur_intervalle_pari_full_bayesien)


    Appartient_intervalle_pari_full_bayesien_exact<- (Intervalle_pari_full_bayesien_exact[,1] <= Valeurs_nouveaux_points)*(Intervalle_pari_full_bayesien_exact[,2] >= Valeurs_nouveaux_points)
    Longueur_intervalle_pari_full_bayesien_exact<- Intervalle_pari_full_bayesien_exact[,2] - Intervalle_pari_full_bayesien_exact[,1]
    Frequence_appartient_intervalle_pari_full_bayesien_exact<- sum(Appartient_intervalle_pari_full_bayesien_exact)/CARDINAL_ENSEMBLE_TEST       
    Longueur_moyenne_intervalle_pari_full_bayesien_exact<- mean(Longueur_intervalle_pari_full_bayesien_exact)

    Longueurs_moyennes_intervalles_pari <- c(Longueur_moyenne_intervalle_pari_juste, Longueur_moyenne_intervalle_pari_argmax_vraisemblance, Longueur_moyenne_intervalle_pari_mode_posterior, Longueur_moyenne_intervalle_pari_full_bayesien_exact, Longueur_moyenne_intervalle_pari_full_bayesien)
    
    Frequences_appartient_intervalle_pari <- c(Frequence_appartient_intervalle_pari_juste, Frequence_appartient_intervalle_pari_argmax_vraisemblance, Frequence_appartient_intervalle_pari_mode_posterior, Frequence_appartient_intervalle_pari_full_bayesien_exact, Frequence_appartient_intervalle_pari_full_bayesien)




## Ecart au carre entre moyenne conditionnelle bayesienne et la vraie moyenne conditionnelle
Ecart_prediction_bayesienne_moyenne_conditionnelle_juste<- sum( ( Moyennes_conditionnelles_bayesiennes - Moyenne_conditionnelle_juste)*( Moyennes_conditionnelles_bayesiennes - Moyenne_conditionnelle_juste) )

Ecart_prediction_moyenne_conditionnelle_juste<- c(Ecart_prediction_argmax_vraisemblance_moyenne_conditionnelle_juste, Ecart_prediction_mode_posterior_moyenne_conditionnelle_juste, Ecart_prediction_bayesienne_moyenne_conditionnelle_juste)
write.matrix(Ecart_prediction_moyenne_conditionnelle_juste, "ecarts_predictions_moyenne_cond_juste.txt",sep = "\t")

    
    ## Ecarts au carre entre moyenne conditionnelle bayesienne et la vraie valeur 
    Ecart_prediction_bayesienne <- sum( ( Moyennes_conditionnelles_bayesiennes - Valeurs_nouveaux_points)*( Moyennes_conditionnelles_bayesiennes - Valeurs_nouveaux_points) )

Ecart_prediction <- c(Ecart_prediction_argmax_vraisemblance, Ecart_prediction_mode_posterior, Ecart_prediction_bayesienne)

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

