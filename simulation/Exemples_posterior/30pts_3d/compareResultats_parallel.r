#library(distr)
library(pscl)
#library(ggplot2)
#library(rgl)
library(MASS)


#source("../../../../Matern_these.r")
#source("../../../../Briques_these.r")
#source("../../../../Krige_these.r")

#source("../../generePlanXP.r")
#source("../../genereMatriceTendance.r")

#source("../../../../dichotomie.r")
#source("../Rastrigin.r")

load("../planXPvariable/Simulations.RData")



#Pour que le germe aleatoire de chaque simulation soit connu :
if(is.null(Simulations[[1]]$germe_aleatoire)) {
for(i in 1:NOMBRE_PROCESSUS) Simulations[[i]]$germe_aleatoire <- GERME_ALEATOIRE_GLOBAL + i
}

CARDINAL_ENSEMBLE_TEST <- 1000#0
POURCENTAGE_INTERVALLE_CONFIANCE <- 95
#POURCENTAGE_INTERVALLE_CONFIANCE_BAYESIEN <- 90

NOMBRE_FONCTIONS_REELLES_TENDANCE <- length(FONCTIONS_REELLES)



QNORM <- qnorm(1 - (1-POURCENTAGE_INTERVALLE_CONFIANCE/100)/2) ## QNORM est le quantile correspondant a la borne droite de l'intervalle de confiance/pari de la loi normale centree reduite entrant en jeu dans les calculs a suivre

library(SparkR)


#Ce wrapper sert a etendre le timeout de SparkR (tres court par defaut) a 72h
connectBackend.orig <- getFromNamespace('connectBackend', pos='package:SparkR')
connectBackend.patched <- function(hostname, port, timeout = 3600*72) {
  connectBackend.orig(hostname, port, timeout)
}

assignInNamespace("connectBackend", value=connectBackend.patched, pos='package:SparkR')
#Fin wrapper

#sparkR.session("local[*]")
#nom_session_sparkR <- paste0("spark://",Sys.info()["nodename"],".athos.hbc.edf.fr:7077")
sparkR.session("spark://10.89.80.11:7077", appName="Evaluations")
#sparkR.session("spark://10.114.116.10:7077",appName="Evaluations")
#sparkR.session(appName="Evaluations")

evalueSimulations <- function(une_simulation)
{

 library(GibbsPosteriorRequired) 

  if(is.null(une_simulation$planXP)) { #Ici, planXP observations et tendance n'ont pas ete stockes dans une_simulation
    x_connus <- generePlanXP(une_simulation$germe_aleatoire,NOMBRE_POINTS_PLANXP = 100, NOMBRE_DIMENSIONS, LHS=FALSE)
    y_connus <- Fonction_emulee(x_connus)
    tendance <- genereMatriceTendance(x_connus,FONCTIONS)
 } else {
x_connus <- une_simulation$planXP ### Attention ! doit figurer dans Simulations
y_connus <- une_simulation$observations ### Attention ! doit figurer dans Simulations
tendance <- une_simulation$tendance ### Attention ! doit figurer dans Simulations
 }
  

  
NOMBRE_DIMENSIONS <- ncol(x_connus)
AMPLITUDE_OBSERVEE <- max(y_connus)-min(y_connus)
TAILLE_PLAN_XP <- nrow(x_connus)

NOMBRE_FONCTIONS_TENDANCE <- ncol(tendance)
if(nrow(tendance)*ncol(tendance)==1) NOMBRE_FONCTIONS_TENDANCE <- 0 # si la matrice est 1x1 (et alors son unique element est 0), c'est un code pour dire "krigeage simple, pas de fonction de tendance"


LONGUEUR_CORRELATION_ARGMAX_VRAISEMBLANCE <- une_simulation$MLE
LONGUEUR_CORRELATION_MODE_POSTERIOR <- une_simulation$MAP

if(nrow(une_simulation$Posterior)%%10 == 0) {
  ECHANTILLON_POSTERIOR <- une_simulation$Posterior[10 * (1:(nrow(une_simulation$Posterior)/10)),]
} else {
  print("L'echantillon de la loi a posteriori n'est pas divisible par 10.")
ECHANTILLON_POSTERIOR <- une_simulation$Posterior
}


if(TYPE_NOYAU_MATERN == "geometrique") 
{
  NOYAU <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE_VRAIE ) # anisotrope geometrique en fait !  
  NOYAU_ARGMAX_VRAISEMBLANCE <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION_ARGMAX_VRAISEMBLANCE, regularite = REGULARITE )
  NOYAU_MODE_POSTERIOR <- creeMaternIsotrope(variance = 1, longueur = LONGUEUR_CORRELATION_MODE_POSTERIOR, regularite = REGULARITE ) 
} else if (TYPE_NOYAU_MATERN == "tensorise") 
{
  NOYAU <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION, regularite = REGULARITE )
  NOYAU_ARGMAX_VRAISEMBLANCE <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION_ARGMAX_VRAISEMBLANCE, regularite = REGULARITE ) 
  NOYAU_MODE_POSTERIOR <- creeMaternTensorise(variance = 1, longueur = LONGUEUR_CORRELATION_MODE_POSTERIOR, regularite = REGULARITE )
} else 
{
  print("Erreur dans compareResultats.r : le noyau doit-il etre anis. geometrique ou tensorise ?")
}


## Matrice de correlation juste : correspondant aux vraies longueurs de correlation
MATRICE_JUSTE <- creeMatriceCovariance(x1=x_connus, x2=x_connus, noyau= NOYAU)
MATRICE_JUSTE_INVERSE <- solve(MATRICE_JUSTE) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.

## Matrice de correlation  correspondant aux longueurs de correlation estimees par MLE : maximum de vraisemblance et MAP : max a posteriori respectivement
MATRICE_ARGMAX_VRAISEMBLANCE <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= NOYAU_ARGMAX_VRAISEMBLANCE)
MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE <- solve(MATRICE_ARGMAX_VRAISEMBLANCE) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.
MATRICE_MODE_POSTERIOR <- creeMatriceCovariance(x1= x_connus, x2= x_connus, noyau= NOYAU_MODE_POSTERIOR)
MATRICE_MODE_POSTERIOR_INVERSE <- solve(MATRICE_MODE_POSTERIOR) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.






NOUVEAU_GERME_ALEATOIRE <- une_simulation$germe_aleatoire * 50  ##Important pour que les differentes machines en parallele ne choisissent pas toutes le meme germe aleatoire !! (*50 pour ne pas trop risquer de prendre des ensembles tests qui auraient ete d'apprentissage dans d'autres simulations)
set.seed(NOUVEAU_GERME_ALEATOIRE)


## Generation de l'ensemble test
Nouveaux_points <- generePointsConnus(nombre_points = CARDINAL_ENSEMBLE_TEST, minima = rep(0,NOMBRE_DIMENSIONS), maxima = rep(1,NOMBRE_DIMENSIONS), lhs=FALSE)
Tous_les_points <- rbind(x_connus, Nouveaux_points)







## Matrice de variance marginale correspondant aux Nouveaux_points
MATRICE_NOUVEAUX_JUSTE <- creeMatriceCovariance(x1= Nouveaux_points, x2= Nouveaux_points, noyau= NOYAU)









# Attention, non seulement la technique suivante est plus efficace que les commentaires, mais elle rend aussi les denominations plus justes. Par consequent,
# CORREL_CONNUS_NOUVEAUX_... actuel = t( CORREL_CONNUS_NOUVEAUX_...  ancien ) et Moyenne_conditionnelle_juste, ainsi que les predictions, sont des vecteurs lignes.

#CORREL_CONNUS_NOUVEAUX_JUSTE <- MATRICE_JUSTE_TOTALE[(nrow(x_connus)+1):(nrow(x_connus)+CARDINAL_ENSEMBLE_TEST),1:nrow(x_connus)]
CORREL_CONNUS_NOUVEAUX_JUSTE <- creeMatriceCovariance(x1= x_connus, x2= Nouveaux_points, noyau= NOYAU)
#MATRICE_ARGMAX_VRAISEMBLANCE_TOTALE <- creeMatriceCovariance(x1= Tous_les_points, x2= Tous_les_points, noyau= NOYAU_ARGMAX_VRAISEMBLANCE)
CORREL_CONNUS_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- creeMatriceCovariance(x1= x_connus, x2= Nouveaux_points, noyau= NOYAU_ARGMAX_VRAISEMBLANCE)
#MATRICE_MODE_POSTERIOR_TOTALE <- creeMatriceCovariance(x1= Tous_les_points, x2= Tous_les_points, noyau= NOYAU_MODE_POSTERIOR)
CORREL_CONNUS_NOUVEAUX_MODE_POSTERIOR <- creeMatriceCovariance(x1= x_connus, x2= Nouveaux_points, noyau= NOYAU_MODE_POSTERIOR)

### ATTENTION : JUSTE merite un traitement separe, car beta est connu.

    ## Calcul des moyennes marginales (qui dependent uniquemeent de beta, et non de sigma^2 ou theta)
    Moyenne_marginale_juste <- 0
    moyenne_reelle <- 0 #Dans le code parallelise, moyenne_reelle doit etre defini ici
    if(NOMBRE_FONCTIONS_REELLES_TENDANCE>0)
    {        
        for(i in 1:NOMBRE_FONCTIONS_REELLES_TENDANCE)
        {
            Moyenne_marginale_juste <- Moyenne_marginale_juste + BETA[i] * apply(Nouveaux_points,1,FONCTIONS_REELLES[[i]])
	    moyenne_reelle <- moyenne_reelle + BETA[i] * apply(x_connus,1,FONCTIONS_REELLES[[i]])
        }
        Moyenne_marginale_juste <- t(Moyenne_marginale_juste)

    ## Calcul des moyenne conditionnelles sachant les valeurs observees et les longeurs de correlation vraies ou estimees
    Moyenne_conditionnelle_juste <- t(y_connus - moyenne_reelle) %*% MATRICE_JUSTE_INVERSE %*%  CORREL_CONNUS_NOUVEAUX_JUSTE + Moyenne_marginale_juste
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

if(NOMBRE_FONCTIONS_TENDANCE>0){

  QR <- qr.Q(qr(tendance),complete=TRUE) #les NOMBRE_FONCTIONS_TENDANCE premieres colonnes de QR contiennent P, les dernieres W
  injectionOrthogonalTendance <- QR[,(NOMBRE_FONCTIONS_TENDANCE+1):ncol(QR)] #W
  injectionTendance <- QR[,1:NOMBRE_FONCTIONS_TENDANCE] #P


  ##Il faut redefinir tous les objets relatifs aux donnes connues



  TENDANCE_NOUVEAUX <-matrix(0,nrow=CARDINAL_ENSEMBLE_TEST, ncol=NOMBRE_FONCTIONS_TENDANCE)
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


  Estimation_beta_argmax_vraisemblance <- solve(t(tendance) %*% solve(MATRICE_ARGMAX_VRAISEMBLANCE) %*% tendance) %*% t(tendance) %*% solve(MATRICE_ARGMAX_VRAISEMBLANCE ) %*% y_connus

  Estimation_beta_mode_posterior <-  solve(t(tendance) %*% solve(MATRICE_MODE_POSTERIOR) %*% tendance) %*% t(tendance) %*% solve(MATRICE_MODE_POSTERIOR) %*% y_connus

  Estimations_beta <- cbind(Estimation_beta_argmax_vraisemblance,Estimation_beta_mode_posterior)

  y_connus <- - t(injectionOrthogonalTendance) %*% y_connus ## Attention, on se base sur MOINS W^T y

  MATRICE_ARGMAX_VRAISEMBLANCE_INJECTEE <- t(injectionOrthogonalTendance) %*% MATRICE_ARGMAX_VRAISEMBLANCE %*% injectionOrthogonalTendance
  MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE <- solve(MATRICE_ARGMAX_VRAISEMBLANCE_INJECTEE)

  MATRICE_MODE_POSTERIOR_INJECTEE <- t(injectionOrthogonalTendance) %*% MATRICE_MODE_POSTERIOR %*% injectionOrthogonalTendance
  MATRICE_MODE_POSTERIOR_INVERSE <- solve(MATRICE_MODE_POSTERIOR_INJECTEE)

}else { ##Cas du krigeage simple
  MATRICE_ARGMAX_VRAISEMBLANCE_INVERSE <- solve(MATRICE_ARGMAX_VRAISEMBLANCE)
  MATRICE_MODE_POSTERIOR_INVERSE <- solve(MATRICE_MODE_POSTERIOR)
  MATRICE_NOUVEAUX_ARGMAX_VRAISEMBLANCE <- 1 
  MATRICE_NOUVEAUX_MODE_POSTERIOR <- 1 
}

















## Calcul des moyenne conditionnelles sachant les valeurs observees et les longeurs de correlation vraies ou estimees


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
  ##mean( apply(Points_LC_parnorm[numero_nouveau_point,,],1,pnormParametresDevant, t=t) ) ##ancienne version
  mean( pnormParametresDevant(parametres =Points_LC_parnorm[numero_nouveau_point,,] , t=t) )
}



## ATTENTION ! Les intervalles de pari sont ici les lignes.
Intervalle_pari_full_bayesien <- matrix(nrow=CARDINAL_ENSEMBLE_TEST,ncol=2)
Intervalle_pari_full_bayesien_exact<- matrix(nrow=CARDINAL_ENSEMBLE_TEST,ncol=2)

Moyennes_conditionnelles_bayesiennes <- apply(Points_LC_parnorm[,,1],1,mean)
Moyennes_des_variances <- apply(Points_LC_parnorm[,,2],1,mean)
Variances_des_moyennes <- apply(Points_LC_parnorm[,,1],1,var)
Variances_conditionnelles_bayesiennes <- Moyennes_des_variances + Variances_des_moyennes

Intervalle_pari_full_bayesien <- cbind(Moyennes_conditionnelles_bayesiennes - QNORM * sqrt(Variances_conditionnelles_bayesiennes) , Moyennes_conditionnelles_bayesiennes +  QNORM * sqrt(Variances_conditionnelles_bayesiennes) )

for(numero_nouveau_point in 1:CARDINAL_ENSEMBLE_TEST)
{


  Intervalle_pari_full_bayesien_exact[numero_nouveau_point,] <- dichotomieIntervallePari(fonctionRepartitionPredictive,numero_nouveau_point= numero_nouveau_point,  probaInf=(1-POURCENTAGE_INTERVALLE_CONFIANCE/100)/2, probaSup=1 - (1-POURCENTAGE_INTERVALLE_CONFIANCE/100)/2, pas = 1)

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

Longueurs_moyennes_intervalles_pari_normalisees <- Longueurs_moyennes_intervalles_pari / AMPLITUDE_OBSERVEE

Frequences_appartient_intervalle_pari <- c(Frequence_appartient_intervalle_pari_juste, Frequence_appartient_intervalle_pari_argmax_vraisemblance, Frequence_appartient_intervalle_pari_mode_posterior, Frequence_appartient_intervalle_pari_full_bayesien_exact, Frequence_appartient_intervalle_pari_full_bayesien)

Estimations_sigma_carre <- c(Estimation_sigma_carre_argmax_vraisemblance,Estimation_sigma_carre_mode_posterior)




## Ecart au carre entre moyenne conditionnelle bayesienne et la vraie moyenne conditionnelle
Ecart_prediction_bayesienne_moyenne_conditionnelle_juste<- sum( ( Moyennes_conditionnelles_bayesiennes - Moyenne_conditionnelle_juste)*( Moyennes_conditionnelles_bayesiennes - Moyenne_conditionnelle_juste) )

Ecart_prediction_moyenne_conditionnelle_juste<- c(Ecart_prediction_argmax_vraisemblance_moyenne_conditionnelle_juste, Ecart_prediction_mode_posterior_moyenne_conditionnelle_juste, Ecart_prediction_bayesienne_moyenne_conditionnelle_juste)

Ecart_prediction_moyenne_conditionnelle_juste_normalise <- Ecart_prediction_moyenne_conditionnelle_juste / AMPLITUDE_OBSERVEE^2




## Ecarts au carre entre moyenne conditionnelle bayesienne et la vraie valeur
Ecart_prediction_bayesienne <- sum( ( Moyennes_conditionnelles_bayesiennes - Valeurs_nouveaux_points)*( Moyennes_conditionnelles_bayesiennes - Valeurs_nouveaux_points) )

Ecart_prediction <- c(Ecart_prediction_argmax_vraisemblance, Ecart_prediction_mode_posterior, Ecart_prediction_bayesienne)
Ecart_prediction_normalise <- Ecart_prediction / AMPLITUDE_OBSERVEE^2

list(ECART_PREDICTION = Ecart_prediction, ECART_PREDICTION_NORMALISE = Ecart_prediction_normalise, 
ECART_PREDICTION_MOYENNE_CONDITIONNELLE_JUSTE = Ecart_prediction_moyenne_conditionnelle_juste,
ECART_PREDICTION_MOYENNE_CONDITIONNELLE_JUSTE_NORMALISE = Ecart_prediction_moyenne_conditionnelle_juste_normalise,
ESTIMATIONS_SIGMA_CARRE = Estimations_sigma_carre, FREQUENCES_APPARTIENT_INTERVALLE_PARI = Frequences_appartient_intervalle_pari, LONGUEURS_MOYENNES_INTERVALLES_PARI = Longueurs_moyennes_intervalles_pari, LONGUEURS_MOYENNES_INTERVALLES_PARI_NORMALISEES = Longueurs_moyennes_intervalles_pari_normalisees)


}



Evaluations <- spark.lapply(list = Simulations,func = evalueSimulations)
#Evaluations <- lapply(X = Simulations,FUN = evalueSimulations)

sparkR.stop()

save.image(file = "Evaluations.RData",safe = TRUE)


MATRICE_FREQUENCES_APPARTIENT_INTERVALLE_PARI <- NULL
for(i in 1:length(Evaluations)) MATRICE_FREQUENCES_APPARTIENT_INTERVALLE_PARI <- rbind(MATRICE_FREQUENCES_APPARTIENT_INTERVALLE_PARI,Evaluations[[i]]$FREQUENCES_APPARTIENT_INTERVALLE_PARI)
write.matrix(MATRICE_FREQUENCES_APPARTIENT_INTERVALLE_PARI, "matrice_frequences_appartient_intervalle_pari.txt", sep = "\t")

MATRICE_LONGUEURS_MOYENNES_INTERVALLES_PARI <- NULL
for(i in 1:length(Evaluations)) MATRICE_LONGUEURS_MOYENNES_INTERVALLES_PARI <- rbind(MATRICE_LONGUEURS_MOYENNES_INTERVALLES_PARI,  Evaluations[[i]]$LONGUEURS_MOYENNES_INTERVALLES_PARI)
write.matrix(MATRICE_LONGUEURS_MOYENNES_INTERVALLES_PARI, "matrice_longueurs_moyennes_intervalles_pari.txt", sep = "\t")


# ## Ecriture des résultats
# 
# Prediction_argmax_vraisemblance <- c(Prediction_argmax_vraisemblance)
# Prediction_mode_posterior <- c(Prediction_mode_posterior)
# 
# Variances_conditionnelles_argmax_vraisemblance <- diag(Variance_conditionnelle_argmax_vraisemblance)
# Variances_conditionnelles_mode_posterior <- diag(Variance_conditionnelle_mode_posterior)
# observations <- scan("observations.txt")
# 
# write.csv(Valeurs_nouveaux_points,"Ensemble_test_valeurs.csv")
# write.csv(Nouveaux_points,"Ensemble_test.csv")
# 
# write.csv(Prediction_argmax_vraisemblance,"Prediction_argmax_vraisemblance.csv")
# write.csv(Variances_conditionnelles_argmax_vraisemblance, "Variances_conditionnelles_argmax_vraisemblance.csv")
# write.csv(t(Intervalle_pari_argmax_vraisemblance), "Intervalle_pari_argmax_vraisemblance.csv")
# write.csv(Appartient_intervalle_pari_argmax_vraisemblance, "Appartient_intervalle_pari_argmax_vraisemblance.csv")
# 
# write.csv(Prediction_mode_posterior, "Prediction_mode_posterior.csv")
# write.csv(Variances_conditionnelles_mode_posterior, "Variances_conditionnelles_mode_posterior.csv")
# write.csv(t(Intervalle_pari_mode_posterior), "Intervalle_pari_mode_posterior.csv")
# write.csv(Appartient_intervalle_pari_mode_posterior, "Appartient_intervalle_pari_mode_posterior.csv")
# 
# write.csv(Moyennes_conditionnelles_bayesiennes, "Moyennes_conditionnelles_bayesiennes.csv")
# write.csv(Variances_conditionnelles_bayesiennes, "Variances_conditionnelles_bayesiennes.csv")
# write.csv(Intervalle_pari_full_bayesien_exact, "Intervalle_pari_bayesien.csv")
# write.csv(Appartient_intervalle_pari_full_bayesien_exact, "Appartient_intervalle_pari_bayesien.csv")
# write.csv(Intervalle_pari_full_bayesien, "Intervalle_pari_bayesien_approx.csv")
# write.csv(Appartient_intervalle_pari_full_bayesien, "Appartient_intervalle_pari_bayesien_approx.csv")
# 
# 
# 
# ## Graphes
# 
# x11() #Necessaire avec Rscript, car bizarrement, sans ceci, Rscript ne cree pas correctement les graphes
# 
# plotMIN <- min(Valeurs_nouveaux_points)
# plotMAX <- max(Valeurs_nouveaux_points)
# 
# plot(Valeurs_nouveaux_points[Appartient_intervalle_pari_argmax_vraisemblance==1], Prediction_argmax_vraisemblance[Appartient_intervalle_pari_argmax_vraisemblance==1], pch=4, xlim=c(plotMIN,plotMAX), ylim=c(plotMIN,plotMAX), xlab="Values at test points", ylab= "Prediction (MLE)", col=4)
# 
# points(Valeurs_nouveaux_points[Appartient_intervalle_pari_argmax_vraisemblance==0], Prediction_argmax_vraisemblance[Appartient_intervalle_pari_argmax_vraisemblance==0], pch=4, col=2)
# 
# lines(c(plotMIN,plotMAX), c(plotMIN,plotMAX), lwd=3)
# points(observations, observations, col=3)
# 
# dev.copy(pdf,"perfMLE.pdf")
# dev.off()
# 
# 
# plot(Valeurs_nouveaux_points[Appartient_intervalle_pari_mode_posterior==1], Prediction_mode_posterior[Appartient_intervalle_pari_mode_posterior==1], pch=4, xlim=c(plotMIN,plotMAX), ylim=c(plotMIN,plotMAX), xlab="Values at test points", ylab= "Prediction (MAP)", col=4)
# 
# points(Valeurs_nouveaux_points[Appartient_intervalle_pari_mode_posterior==0], Prediction_mode_posterior[Appartient_intervalle_pari_mode_posterior==0], pch=4, col=2)
# 
# lines(c(plotMIN,plotMAX), c(plotMIN,plotMAX), lwd=3)
# points(observations, observations, col=3)
# 
# dev.copy(pdf,"perfMAP.pdf")
# dev.off()
# 
# 
# plot(Valeurs_nouveaux_points[Appartient_intervalle_pari_full_bayesien_exact==1], Moyennes_conditionnelles_bayesiennes[Appartient_intervalle_pari_full_bayesien_exact==1], pch=4, xlim=c(plotMIN,plotMAX), ylim=c(plotMIN,plotMAX), xlab="Values at test points", ylab= "Prediction (FPD)", col=4)
# 
# points(Valeurs_nouveaux_points[Appartient_intervalle_pari_full_bayesien_exact==0], Moyennes_conditionnelles_bayesiennes[Appartient_intervalle_pari_full_bayesien_exact==0], pch=4, col=2)
# 
# lines(c(plotMIN,plotMAX), c(plotMIN,plotMAX), lwd=3)
# points(observations, observations, col=3)
# 
# dev.copy(pdf,"perfFPD.pdf")
# dev.off()
# 
# 
# plot(Valeurs_nouveaux_points[Appartient_intervalle_pari_full_bayesien==1], Moyennes_conditionnelles_bayesiennes[Appartient_intervalle_pari_full_bayesien==1], pch=4, xlim=c(plotMIN,plotMAX), ylim=c(plotMIN,plotMAX), xlab="Values at test points", ylab= "Prediction (approximate FPD)", col=4)
# 
# points(Valeurs_nouveaux_points[Appartient_intervalle_pari_full_bayesien==0], Moyennes_conditionnelles_bayesiennes[Appartient_intervalle_pari_full_bayesien==0], pch=4, col=2)
# 
# lines(c(plotMIN,plotMAX), c(plotMIN,plotMAX), lwd=3)
# points(observations, observations, col=3)
# 
# dev.copy(pdf,"perfFPDapprox.pdf")
# dev.off()
# 
# 
# #    write.matrix(Ecart_prediction, "ecarts_predictions.txt",sep = "\t")
# 
# 
# #        fonctionDensitePredictive <- function(t, numero_nouveau_point)    
# #        {
# #            mean( apply(Points_LC_parnorm[numero_nouveau_point,,],1,dnormParametresDevant, t=t) )
# #        }
# 
# #    #xxx <- 0.05*-100:100
# #    xxx <- Moyenne_conditionnelle_juste[1] + 0.01*-100:100
# #    yyy <- apply(matrix(xxx),1,fonctionDensitePredictive,1)
# #    #YYY <- apply(matrix(xxx),1,fonctionRepartitionPredictive,1)
# #    zzz <- dnorm(xxx, mean= mean(Points_LC_parnorm[1,,1]), sd= sqrt( mean(Points_LC_parnorm[1,,2]) + var(Points_LC_parnorm[1,,1]) ) )
# #    #ZZZ <- pnorm(xxx, mean= mean(Points_LC_parnorm[1,,1]), sd= sqrt( mean(Points_LC_parnorm[1,,2]) + var(Points_LC_parnorm[1,,1]) ) )
# #    #ttt <- dnorm(xxx, mean= mean(Points_LC_parnorm[1,,1]), sd= sqrt( var(Points_LC_parnorm[1,,1]) ) )
# ##    plot(xxx, dnorm(xxx,mean=Moyenne_conditionnelle_juste[1], sd= sqrt(Variance_conditionnelle_juste[1,1])),type="l", xlab="Unobserved value", ylab="Predictive probability density", main="", lwd=2)
# #    #lines(xxx, dnorm(xxx,mean=Prediction_argmax_vraisemblance[1], sd=sqrt(Variance_conditionnelle_argmax_vraisemblance[1,1])),col=3)
# #    lines(xxx, dnorm(xxx,mean=Prediction_mode_posterior[1], sd=sqrt(Variance_conditionnelle_mode_posterior[1,1])),col=2, lty=2, lwd=2)
# #    lines(xxx, yyy, col=4, lty=3, lwd=2)
# #    lines(xxx, zzz, col=4, lty=3, lwd=1)
# #    #abline(v=Intervalle_pari_full_bayesien_exact[1,1], col=4, lwd=2, lty=3)
# #    #abline(v=Intervalle_pari_full_bayesien_exact[1,2], col=4, lwd=2, lty=3)
# #    #abline(v=Intervalle_pari_full_bayesien[1,1], col=4, lty=3)
# #    #abline(v=Intervalle_pari_full_bayesien[1,2], col=4, lty=3)
# 
# #    legende <- c("True","MAP","full-Bayes.", "approx. full-Bayes.")
# #    col <- c(1,2,4,4)
# #    lty <- c(1,2,3,3)
# #    lwd <- c(2,2,2,1)
# 
# #    legend("topright",legende, col=col, lty=lty, lwd=lwd)
# 
# #    dev.copy(pdf,paste0("Graphes/dens_predictive_",numero_germe_aleatoire,".pdf"))
# #    dev.off()
# 
# #    #lines(xxx, ttt, col=5)
# #    #abline(v=mean(Points_LC_parnorm[1,,1], col=4))
# 
# #    #plot(xxx, ttt, col=5, type="l")
# #    #lines(xxx, zzz, col=4)
# 
# #                                        #plot(xxx, (zzz-yyy)/yyy, type="l", col=4)
# #    #plot(xxx, ZZZ-YYY, type="l", col=4)
# #    #print(var(Points_LC_parnorm[1,,1]) / mean(Points_LC_parnorm[1,,2]) )
# 
# #    #print(Points_LC_parnorm[1,,2])
# 
# 
