## Calcul du prior quasi-Jeffreys


## Sortie : nablaG(H beta) quelle que soit la transformation utilisée
## Dépendances : aucune
## Attention : fonctionA et fonctionBprime sont des fonctions R --> R. En les exécutant sur un vecteur, on les applique en fait à chaque composante dudit vecteur.
calculeNablaG_moyenne <- function ( variance, moyenne, fonctionA, fonctionBprime, inverseMatriceCorrelation )
{
	# Moyenne <- Fonctions_points_connus %*% beta

##    print(fonctionA)
                                        #	print(fonctionBprime)
##    print(dim(fonctionA(moyenne)))
	- 1/variance * inverseMatriceCorrelation %*% fonctionA(moyenne) + fonctionBprime(moyenne)

	
}

## Sortie : nabla^2 G(H beta) quelle que soit la transformation utilisée
## Dépendances : aucune
## Attention : fonctionAprime et fonctionBseconde sont des fonctions R --> R. En les exécutant sur un vecteur, on les applique en fait à chaque composante dudit vecteur.
calculeNablaCarreG_moyenne <- function( variance, moyenne, fonctionAprime, fonctionBseconde, inverseMatriceCorrelation)
{
	terme1 <- -1/variance * inverseMatriceCorrelation %*% diag( fonctionAprime(c(moyenne)) ) ## pour que diag donne une matrice diagonale, il faut que diag prenne un vecteur en entrée

	terme1 + t(terme1) + diag( fonctionBseconde(moyenne) )
}

## Sortie : nabla^2 G(H beta) dans le cas particulier de la transformation de Box-Cox (utilise la programmation générique de calculeNablaG_moyenne)
## Dépendance : calculeNablaCarreG_moyenne
creeNablaCarreG_moyenne_boxCox <- function( variance, moyenne, lambda, inverseMatriceCorrelation)
{
	fonctionAprime <- function(t) log( abs( 1 + lambda * t) ) 
	fonctionBseconde <- function(t) - lambda / ( 1 + lambda * t )^2

	calculeNablaCarreG_moyenne( variance= variance, moyenne= moyenne, fonctionAprime= fonctionAprime, fonctionBseconde= fonctionBseconde, inverseMatriceCorrelation= inverseMatriceCorrelation )
}

## Sortie : nablaG(H beta) dans le cas particulier de la transformation de Box-Cox (utilise la programmation générique de calculeNablaG_moyenne)
## Dépendance : calculeNablaG_moyenne
creeNablaG_moyenne_boxCox <- function( variance, moyenne, lambda, inverseMatriceCorrelation)
{
	fonctionA <- function(t) 1/lambda * ( log( abs( 1 + lambda * t) ) * t + 1/lambda * log( abs( 1 + lambda * t) ) - t )
	fonctionBprime <- function(t) 1/ ( 1 + lambda * t )

	calculeNablaG_moyenne( variance= variance, moyenne= moyenne, fonctionA= fonctionA, fonctionBprime= fonctionBprime, inverseMatriceCorrelation= inverseMatriceCorrelation )
}

## Sortie : log(prior 2e ordre) quelle que soit la transformation utilisée
## Dépendances : aucune
calculePseudoPriorI_kriUniv <- function( variance, matriceCorrelation, Fonctions_points_connus, nablaG_moyenne)
{
#	variance <- noyau$parametres$variance
	nb_points <- nrow(matriceCorrelation)
	t_Fonctions_points_connus <- t(Fonctions_points_connus)
##        print("t_Fonctions_points_connus")
##        print(t_Fonctions_points_connus)
##        print(paste("variance = ",variance))
##        print(solve( matriceCorrelation, Fonctions_points_connus ))
##        print(   t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) )
	nablaG_moyenne <- matrix(nablaG_moyenne) ## pour le calcul matriciel, il faut que nablaG_moyenne soit une matrice à 1 colonne et non un vecteur
	auxiliaire0 <- 0 #nablaCarreG_moyenne %*% matriceCorrelation
	auxiliaire1 <- variance * c( t(nablaG_moyenne) %*% matriceCorrelation %*% nablaG_moyenne) ## c() car auxiliare doit être vu comme un scalaire et non comme une matrice 1x1
	auxiliaire2 <- 0 # variance^2 / 2 * sum( diag( auxiliaire0 %*% auxiliaire0 ) )
	auxiliaire3 <- 0 #- variance^2 / (2 * nb_points) * (sum ( diag ( auxiliaire0 ) ))^2
##        print(paste("auxiliaire1 = ", auxiliaire1))
##        print(paste("auxiliaire2 = ", auxiliaire2))
##        print(paste("auxiliaire3 = ", auxiliaire3))
##        print(paste("det = ", det( 1/variance *  t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) - 1/(auxiliaire1 + auxiliaire2 + auxiliaire3) * (t_Fonctions_points_connus %*% nablaG_moyenne) %*% (t(nablaG_moyenne) %*% Fonctions_points_connus)  ) ) )

        
	- log(variance) +
            1/2 * ( log( nb_points / 2 *  (auxiliaire1 + auxiliaire2 + auxiliaire3 ) ) +
                       log( det( 1/variance *  t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) - 1/(auxiliaire1 + auxiliaire2 + auxiliaire3) * (t_Fonctions_points_connus %*% nablaG_moyenne) %*% (t(nablaG_moyenne) %*% Fonctions_points_connus)  ) )
                   )
}

## Sortie : log(prior 2e ordre) quelle que soit la transformation utilisée
## Dépendances : aucune
calculePseudoPriorII_kriUniv <- function( variance, matriceCorrelation, Fonctions_points_connus, nablaG_moyenne, nablaCarreG_moyenne)
{
#	variance <- noyau$parametres$variance
	nb_points <- nrow(matriceCorrelation)
	t_Fonctions_points_connus <- t(Fonctions_points_connus)
##        print("t_Fonctions_points_connus")
##        print(t_Fonctions_points_connus)
##        print(paste("variance = ",variance))
##        print(solve( matriceCorrelation, Fonctions_points_connus ))
##        print(   t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) )
	nablaG_moyenne <- matrix(nablaG_moyenne) ## pour le calcul matriciel, il faut que nablaG_moyenne soit une matrice à 1 colonne et non un vecteur
	auxiliaire0 <- nablaCarreG_moyenne %*% matriceCorrelation
	auxiliaire1 <- variance * c( t(nablaG_moyenne) %*% matriceCorrelation %*% nablaG_moyenne) ## c() car auxiliare doit être vu comme un scalaire et non comme une matrice 1x1
	auxiliaire2 <- variance^2 / 2 * sum( diag( auxiliaire0 %*% auxiliaire0 ) )
	auxiliaire3 <- - variance^2 / (2 * nb_points) * (sum ( diag ( auxiliaire0 ) ))^2
##        print(paste("auxiliaire1 = ", auxiliaire1))
##        print(paste("auxiliaire2 = ", auxiliaire2))
##        print(paste("auxiliaire3 = ", auxiliaire3))
##        print(paste("det = ", det( 1/variance *  t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) - 1/(auxiliaire1 + auxiliaire2 + auxiliaire3) * (t_Fonctions_points_connus %*% nablaG_moyenne) %*% (t(nablaG_moyenne) %*% Fonctions_points_connus)  ) ) )

        
	- log(variance) +
            1/2 * ( log( nb_points / 2 *  (auxiliaire1 + auxiliaire2 + auxiliaire3 ) ) +
                       log( det( 1/variance *  t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) - 1/(auxiliaire1 + auxiliaire2 + auxiliaire3) * (t_Fonctions_points_connus %*% nablaG_moyenne) %*% (t(nablaG_moyenne) %*% Fonctions_points_connus)  ) )
                   )
}



## Sortie : log(prior 2e ordre) dans le cas de la transformation de Box-Cox
## Dépendances:	calculePseudoPrior_kriUniv
##		creeMatrices -> ... [Krige.r]
commande_calculePseudoPriorI_kriUniv_boxCox <- function( x_connus, noyau, liste_fonctions, beta, Matrices_utiles= NULL)
{
 if(is.null(Matrices_utiles)) Matrices_utiles <- creeMatrices( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions)

	moyenne <- Matrices_utiles$Fonctions_points_connus %*% beta

	nablaG_moyenne <- creeNablaG_moyenne_boxCox( variance= noyau$parametres$variance, moyenne= moyenne, lambda= noyau$boxCox, inverseMatriceCorrelation= Matrices_utiles$inverseMatriceCorrelation)

	calculePseudoPriorI_kriUniv( variance= noyau$parametres$variance, matriceCorrelation= Matrices_utiles$matriceCorrelation, Fonctions_points_connus= Matrices_utiles$Fonctions_points_connus, nablaG_moyenne= nablaG_moyenne )
}

## Sortie : log(prior 2e ordre) dans le cas de la transformation de Box-Cox
## Dépendances:	calculePseudoPrior_kriUniv
##		creeMatrices -> ... [Krige.r]
commande_calculePseudoPriorII_kriUniv_boxCox <- function( x_connus, noyau, liste_fonctions, beta, Matrices_utiles= NULL)
{
 if(is.null(Matrices_utiles)) Matrices_utiles <- creeMatrices( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions)

	moyenne <- Matrices_utiles$Fonctions_points_connus %*% beta

	nablaG_moyenne <- creeNablaG_moyenne_boxCox( variance= noyau$parametres$variance, moyenne= moyenne, lambda= noyau$boxCox, inverseMatriceCorrelation= Matrices_utiles$inverseMatriceCorrelation)
	nablaCarreG_moyenne <- creeNablaCarreG_moyenne_boxCox( variance= noyau$parametres$variance, moyenne= moyenne, lambda= noyau$boxCox, inverseMatriceCorrelation= Matrices_utiles$inverseMatriceCorrelation)

	calculePseudoPriorII_kriUniv( variance= noyau$parametres$variance, matriceCorrelation= Matrices_utiles$matriceCorrelation, Fonctions_points_connus= Matrices_utiles$Fonctions_points_connus, nablaG_moyenne= nablaG_moyenne, nablaCarreG_moyenne= nablaCarreG_moyenne )
}


## Sortie : log(prior 2e ordre) dans le cas de la transformation de Box-Cox
## Dépendances:	calculePseudoPrior_kriUniv
##		creeMatrices -> ... [Krige.r]
commande_calculePseudoPriorIII_kriUniv_boxCox <- function( noyau)
{

	- log( noyau$boxCox ) - 1/2 * log( noyau$parametres$variance)
}


## Sortie : liste comportant les pseudo-vraisemblances (en tenant compte du prior) atteintes en les minima trouvees, les warnings (0= tout va bien), les jeux de param argmax pseudo-vraisemblance, les jeux de param initiaux, MAIS PAS A^2 d'Anderson-Darling par valid directe et LOO
## NB : les numeros de lignes des matrices correspondent aux numeros d'entrees des vecteurs
## NB2: les initialisations se font par tirage aléatoire de loi uniforme déterminée par bornesInitialisation. Pour eviter les minima locaux, on reessaie nb_iterations fois.
## NB3: valable uniquement avec la transformation de Box-Cox
## Dependances:	lhs::geneticLHS
##		calculeFonctionAMinimiserVectoriel_BoxCox --> ... [Krige.r] [Matern.r] [Box-Cox.r]
## ATTENTION : ici, la regularite est supposee connue
trouveArgmaxLoiJointe_BoxCox <- function( x_connus, y_connus, regularite, liste_fonctions, bornesInitialisation= list( boxCox=list(minima=0,maxima=1), longueur= list(minima= c(0.001,0.001,0.001,0.001,0.001,0.001), maxima= c( 0.3, 1,2, 10, 6,6) ), beta= list(minima= c(0,0), maxima=c(1,1) ) ),  nb_iterations=100, avecPrior= TRUE)
    {
        nb_dim <-   length(bornesInitialisation$longueur$maxima )
        nb_regresseurs <- length(liste_fonctions)
        matrice_argmax_vraisemblance <-  matrix(rep(1000, nb_iterations * (1+nb_dim+nb_regresseurs)), ncol=1+nb_dim+nb_regresseurs)
        matrice_param_init <- matrix(rep(1000, nb_iterations * (1+nb_dim+nb_regresseurs)), ncol=1+nb_dim+nb_regresseurs)
        vect_pseudo_vraisemblance <- rep(1000, nb_iterations)
        vect_warning <- rep(0, nb_iterations)
##	matrice_validation <- matrix( rep(1000, nb_iterations*2)  , ncol=2)


        matrice_param_init <- geneticLHS(n=nb_iterations, k= 1 + nb_dim + nb_regresseurs)
        matrice_param_init[,1] <- bornesInitialisation$boxCox$minima + ( bornesInitialisation$boxCox$maxima - bornesInitialisation$boxCox$minima ) *matrice_param_init[,1]

        for( num_dim in 1:nb_dim ) matrice_param_init[,num_dim+1] <- bornesInitialisation$longueur$minima[num_dim]+ ( bornesInitialisation$longueur$maxima[num_dim]- bornesInitialisation$longueur$minima[num_dim]) *matrice_param_init[,num_dim+1]

        for( num_regresseur in 1:nb_regresseurs ) matrice_param_init[,(num_regresseur+1+nb_dim)] <- bornesInitialisation$beta$minima[num_regresseur]+ ( bornesInitialisation$beta$maxima[num_regresseur]- bornesInitialisation$beta$minima[num_regresseur]) *matrice_param_init[,(num_regresseur+1+nb_dim)]        
        
        for ( iter in 1:nb_iterations )
            {
                
                OPT <- try(optim( par= matrice_param_init[iter,] , fn= calculeFonctionAMinimiserVectoriel_BoxCox, x_connus= x_connus, y_connus= y_connus, regularite= regularite, liste_fonctions= liste_fonctions, beta= beta, avecPrior= avecPrior, method="BFGS" ) )

                print(OPT)
                if(class(OPT)=="try-error" )
                    {
                        vect_warning[iter] <- NA
                        matrice_argmax_vraisemblance[iter,] <- NA
                        vect_pseudo_vraisemblance[iter] <- NA
                    }
                else
                    {
                vect_warning[iter] <- OPT$convergence
                matrice_argmax_vraisemblance[iter,] <- OPT$par
                vect_pseudo_vraisemblance[iter] <- OPT$value
            }

##		matrice_validation[iter,] <- valide_parametres( x_connus= x_connus, y_connus= y_connus, parametres= list( regularite= Inf, longueur= matrice_argmax_vraisemblance[iter, -1] ), boxCox= matrice_argmax_vraisemblance[iter,1], esperanceFonctionsRegression= liste_fonctions )

		}

        list( vect_pseudo_vraisemblance= vect_pseudo_vraisemblance, vect_warning= vect_warning, matrice_argmax_vraisemblance= matrice_argmax_vraisemblance, matrice_param_init= matrice_param_init) ##, matrice_validation= matrice_validation)
    }



## Sortie : pseudo-vraisemblance des observations y_connus en fonction du noyau, de la base de fonctions d'esperance et du parametre de Box-Cox utilisés
## Dependances:	calculeFonctionAMinimiser_BoxCox --> ... [Krige.r] [Box-Cox.r]
##		creeMaternIsotrope [Matern.r]
## ATTENTION : ici, on impose a la regularite d'etre connue
calculeFonctionAMinimiserVectoriel_BoxCox<- function( vecteur_param_noyau, x_connus, y_connus, regularite, liste_fonctions, beta, avecPrior= TRUE)
    {
nb_param <- length(vecteur_param_noyau)
nb_regresseurs <- length(liste_fonctions)
        boxCox <- vecteur_param_noyau[1]
##        regularite <- Inf
longueur <- abs( vecteur_param_noyau[2:(nb_param-nb_regresseurs)] )
beta <- vecteur_param_noyau[(nb_param-nb_regresseurs+1):nb_param]
        noyau <- creeMaternIsotrope(variance=1, longueur= longueur, regularite= regularite)

        res= calculeFonctionAMinimiser_BoxCox( x_connus= x_connus, y_connus= y_connus, noyau= noyau, liste_fonctions = liste_fonctions,  boxCox= boxCox, beta= beta, avecPrior= avecPrior )
##        if(boxCox<0) res= 100000
        print(paste("boxCox = ",boxCox))
        print(paste("reg = ",regularite))
        print(paste("lon = ",longueur))
        print(res)
        print("-----------------------------------")
##        res
        res
    }



## Sortie : pseudo-vraisemblance des observations y_connus en fonction du noyau, de la base de fonctions d'esperance et du parametre de Box-Cox utilisés
## Dependances:	creeMatrices --> ... [Krige.r]
##		pseudoVraisemblance_BoxCox --> BoxCox [Box-Cox.r]
##		commande_calcule_PseudoPrior_kriUniv_BoxCox -> ...
calculeFonctionAMinimiser_BoxCox <-  function ( x_connus, y_connus, noyau, liste_fonctions,  boxCox, beta, avecPrior= TRUE )
    {
        Matrices <- creeMatrices ( x_connus= x_connus,  noyau= noyau, liste_fonctions= liste_fonctions )

        Num <- BoxCox(y_connus, lambda= boxCox)
        
        noyau$boxCox <- boxCox
#####   noyau$parametres$variance <- c(1/(nrow(x_connus) + length(liste_fonctions) + 1) * t(Num) %*% Matrices$matriceFormeQuadratique %*% Num) ## c() pour que Variance soit un scalaire, et non une matrice 1x1
        
        noyau$parametres$variance <- c(1/(nrow(x_connus) + length(liste_fonctions) + 1) * t(Num - Matrices$Fonctions_points_connus %*% beta) %*% Matrices$inverseMatriceCorrelation %*% (Num - Matrices$Fonctions_points_connus %*% beta)) ## c() pour que Variance soit un scalaire, et non une matrice 1x1




        
##        Beta_approche<- Matrices$matriceConstructionBeta %*% Num
##        print("Beta_approche_pourdefaux")
##        print(Beta_approche_pourdefaux)
##        Beta_approche <- c(1,2)
        
##        Moyenne <- Matrices$Fonctions_points_connus %*% Beta_approche
##        nablaG_Moyenne <- creeNablaG_moyenne_boxCox( moyenne= Moyenne, lambda = boxCox, matriceCorrelation= Matrices$matriceCorrelation)


##        Variance <- c(1/(nrow(x_connus) + length(liste_fonctions) + 1) * t(Num) %*% Matrices$matriceFormeQuadratique %*% Num) ## c() pour que Variance soit un scalaire, et non une matrice 1x1

        res <- pseudoVraisemblance_BoxCox( y_connus= y_connus, matriceCorrelation= Matrices$matriceCorrelation, variance= noyau$parametres$variance, boxCox= boxCox )
        if(avecPrior==1) res <- res - 2 / length(y_connus) * commande_calculePseudoPriorI_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= beta, Matrices_utiles= Matrices)
        if(avecPrior==2) res <- res - 2 / length(y_connus) * commande_calculePseudoPriorII_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= beta, Matrices_utiles= Matrices)
        if(avecPrior==3) res <- res - 2 / length(y_connus) * commande_calculePseudoPriorII_kriUniv_boxCox(  noyau= noyau)
        res
    }


## Sortie : pseudo-vraisemblance des observations y_connus en fonction de la matrice de corrélation  matriceCorrelation de y_connus et du paramètre lambda de Box-Cox (noté boxCox)
## Dépendances : BoxCox [Box-Cox.r]
pseudoVraisemblance_BoxCox<- function( y_connus, matriceCorrelation,variance,  boxCox ) ## INUTILE , Fonctions_points_connus, nablaG_moyenne)
    {
        nbPointsConnus <- nrow(matriceCorrelation) # = ncol(matriceCorrelation) car matriceCorrelation est carrée
        num <-  BoxCox(y_connus, lambda= boxCox)

	res <- 1/nbPointsConnus *  log(det(matriceCorrelation)) +log( variance ) - 2* (boxCox-1) / nbPointsConnus * sum(log(abs(y_connus)) )
       
 #####       res <- 1/nbPointsConnus *  log(det(matriceCorrelation)) +log( 1/ nbPointsConnus  * t( num ) %*% matriceFormeQuadratique %*% num  ) - 2* (boxCox-1) / nbPointsConnus * sum(log(abs(y_connus)) )

        ## INUTILE  + calculePseudoPrior_kriUniv( variance= variance, matriceCorrelation= matriceCorrelation, Fonctions_points_connus= Fonctions_points_connus, nablaG_moyenne= nablaG_moyenne)

        print(paste("pseudo-vraisemblance = ", res ))
        print(log(det(matriceCorrelation)))
        
        res

    }



## Sortie : liste de matrices nécessaires aux calculs
## Equivalences Bachoc : matriceFormeQuadratique <-> PI_theta, Fonctions_points_connus <-> H
## Dépendances : creeMatriceCovariance, creeMatriceValeursFonctions [Krige.r]
creeMatrices <- function(x_connus, noyau, liste_fonctions= NULL)
{
        noyau$parametres$variance <- 1 # Pour que la matrice soit de corrélation, non de covariance

        matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= noyau)

        inverseMatriceCorrelation <- solve(matriceCorrelation)

        if(is.null(liste_fonctions))
            
            {
                
##                matriceFormeQuadratique <- inverseMatriceCorrelation
		Fonctions_points_connus <- NULL
                matriceConstructionBeta <- NULL

            }

        else
            {
                Fonctions_points_connus <- creeMatriceValeursFonctions ( matrice_points = x_connus , liste_fonctions = liste_fonctions )

                aux <- solve( t( Fonctions_points_connus ) %*% inverseMatriceCorrelation %*% Fonctions_points_connus )

##                matriceConstructionBeta <- aux %*% t( Fonctions_points_connus ) %*% inverseMatriceCorrelation
                
##                matriceFormeQuadratique <- inverseMatriceCorrelation - inverseMatriceCorrelation %*% Fonctions_points_connus %*% matriceConstructionBeta

            }

        ##        list( matriceCorrelation= matriceCorrelation, matriceFormeQuadratique= matriceFormeQuadratique, Fonctions_points_connus= Fonctions_points_connus , matriceConstructionBeta= matriceConstructionBeta)
         list( matriceCorrelation= matriceCorrelation, inverseMatriceCorrelation= inverseMatriceCorrelation, Fonctions_points_connus= Fonctions_points_connus )
    }





## x_connus_nuage<- generePointsConnus( nombre_points = 25, minima = c(0.1,0.1), maxima = c(0.9,0.9) )
## x_connus_croix <- genereCroix( matrice_centres <- matrix(x_connus_nuage[1:2,],nrow=2), vecteur_nombre_points_par_dimension= c(4,4), vecteur_espacements= c(0.03, 0.03) )
## x_test_nuage <- generePointsConnus( nombre_points = c(25), minima = c(0,0), maxima = c(1,1) )
## x_test_croix <- genereCroix(  matrice_centres <- matrix(x_test_nuage[1,],nrow=1), vecteur_nombre_points_par_dimension=6, vecteur_espacements= c(0.03, 0.03) )


## x_complet <- rbind(x_connus_nuage, x_connus_croix)##, x_test_nuage, x_test_croix)

## y_complet <- apply( x_complet, 1, esperanceAffine(1,2) ) + creeEchantillonNormal(x_connus= x_complet, noyau= creeMaternIsotrope(variance=1, longueur= c(0.3,0.4), regularite= 1.5) )

## y_complet_boxCoxise <- BoxCox_reciproque(y_complet, lambda = 0.5)


## RESULTAT<- trouveArgmaxLoiJointe_BoxCox(x_connus = x_complet, y_connus = y_complet_boxCoxise, regularite = 1.5, liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = 2, bornesInitialisation = list( boxCox= list(minima=0, maxima=1), longueur= list(minima= c(0,0), maxima=c(1,1)) , beta= list(minima=c(0,0), maxima=c(5,5) ) ) )




####### OBSOLETE ########
##compare_Prior_sansPrior <-  function ( nb_iterations, x_connus, vrai_noyau, intercept_moyenne, pente_moyenne)
##    {
##        Matrice_res_sansPrior <- matrix(1000, nrow= nb_iterations, ncol= (1 + ncol(x_connus) + 2))
##        Matrice_res_avecPrior <- matrix(1000, nrow= nb_iterations, ncol= (1 + ncol(x_connus) + 2))
##        print(vrai_noyau)
##        for ( iter in 1:nb_iterations )
##            {
##                y_connus <- BoxCox_reciproque (  apply( x_connus, 1, esperanceAffine(a= intercept_moyenne, b= pente_moyenne) ) + creeEchantillonNormal(x_connus= x_connus, noyau= vrai_noyau )   , lambda= vrai_noyau$boxCox )
##
##                res_sansPrior <- trouveArgmaxLoiJointe_BoxCox(x_connus = x_connus, y_connus = y_connus, regularite = vrai_noyau$parametres$regularite, liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = 2, bornesInitialisation = list( boxCox= list(minima=0, maxima=1), longueur= list(minima= c(0,0), maxima=c(1,1)) , beta= list(minima=c(0,0), maxima=c(5,5) ) ), avecPrior= FALSE )

##                print(res_sansPrior)
##                Matrice_res_sansPrior[iter,] <- res_sansPrior$matrice_argmax_vraisemblance[1,]                
##                if( sd( res_sansPrior$matrice_argmax_vraisemblance[,1] ) > 0.1 )  Matrice_res_sansPrior[iter,] <- NA ## la valeur de lambda est notre juge de paix : comment savoir quelle est la bonne valeur si les optima obtenus sont trop disticts ?





##                 res_avecPrior <- trouveArgmaxLoiJointe_BoxCox(x_connus = x_connus, y_connus = y_connus, regularite = vrai_noyau$parametres$regularite, liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = 2, bornesInitialisation = list( boxCox= list(minima=0, maxima=1), longueur= list(minima= c(0,0), maxima=c(1,1)) , beta= list(minima=c(0,0), maxima=c(5,5) ) ), avecPrior= TRUE )

##                Matrice_res_avecPrior[iter,] <- res_avecPrior$matrice_argmax_vraisemblance[1,]
##                print( res_avecPrior$matrice_argmax_vraisemblance )
##                if( sd( res_sansPrior$matrice_argmax_vraisemblance[,1] ) > 0.1 )  Matrice_res_avecPrior[iter,] <- NA ## la valeur de lambda est notre juge de paix : comment savoir quelle est la bonne valeur si les optima obtenus sont trop disticts ?
##            }

##        list( Matrice_res_sansPrior= Matrice_res_sansPrior, Matrice_res_avecPrior= Matrice_res_avecPrior )
##    }

######## FIN OBSOLETE #####



## Sortie : liste de matrices permettant la comparaison entre la recherche de l'optimum de la vraisemblance et de la loi jointe (ie vraisemblance pondérée ou non par le prior au 2e ordre) (1er paramètres : lambda , 2e paramètre : norme 2 de l'erreur commise sur les longueurs de corrélation
## Hypothèses : moyenne affine
## nb_iterations désigne le nb de points de départ de l'optimisation
## Dépendances:	trouveArgmaxLoiJointe_BoxCox -->  lhs::geneticLHS ; calculeFonctionAMinimiserVectoriel_BoxCox --> ... [Krige.r] [Matern.r] [Box-Cox.r]		
represente_argoptima <- function( nb_iterations, x_connus, vrai_noyau, intercept_moyenne, pente_moyenne, pepite=0 )
    {
        Vraie_valeur <- c( vrai_noyau$boxCox, 0 )

         y_connus <- BoxCox_reciproque (  apply( x_connus, 1, esperanceAffine(a= intercept_moyenne, b= pente_moyenne) ) + creeEchantillonNormal(x_connus= x_connus, noyau= vrai_noyau, pepite=pepite )   , lambda= vrai_noyau$boxCox )

        Resultat_sansPrior <- trouveArgmaxLoiJointe_BoxCox(x_connus = x_connus, y_connus = y_connus, regularite = vrai_noyau$parametres$regularite, liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = nb_iterations , bornesInitialisation = list( boxCox= list(minima=0, maxima=1), longueur= list(minima= c(0,0), maxima=c(1,1)) , beta= list(minima=c(0,0), maxima=c(5,5) ) ), avecPrior= FALSE )

        Resultat_avecPriorI <- trouveArgmaxLoiJointe_BoxCox(x_connus = x_connus, y_connus = y_connus, regularite = vrai_noyau$parametres$regularite, liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = nb_iterations , bornesInitialisation = list( boxCox= list(minima=0, maxima=1), longueur= list(minima= c(0,0), maxima=c(1,1)) , beta= list(minima=c(0,0), maxima=c(5,5) ) ), avecPrior= TRUE )

        Resultat_avecPriorII <- trouveArgmaxLoiJointe_BoxCox(x_connus = x_connus, y_connus = y_connus, regularite = vrai_noyau$parametres$regularite, liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = nb_iterations , bornesInitialisation = list( boxCox= list(minima=0, maxima=1), longueur= list(minima= c(0,0), maxima=c(1,1)) , beta= list(minima=c(0,0), maxima=c(5,5) ) ), avecPrior= 2 )

        Resultat_avecPriorIII <- trouveArgmaxLoiJointe_BoxCox(x_connus = x_connus, y_connus = y_connus, regularite = vrai_noyau$parametres$regularite, liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = nb_iterations , bornesInitialisation = list( boxCox= list(minima=0, maxima=1), longueur= list(minima= c(0,0), maxima=c(1,1)) , beta= list(minima=c(0,0), maxima=c(5,5) ) ), avecPrior= 3 )

        if( sd( Resultat_sansPrior$matrice_argmax_vraisemblance[,1] ) > 0.1 ) print("Problème : je ne parviens pas à déterminer le maximum de vraisemblance.") ## la valeur de lambda est notre juge de paix : comment savoir quelle est la bonne valeur si les optima obtenus sont trop disticts ?


        Argoptima_locaux_sansPrior <- matrix(NA, nrow=nb_iterations, ncol= 2)
        for ( iter in 1:nb_iterations ) Argoptima_locaux_sansPrior[iter,]<- c( Resultat_sansPrior$matrice_argmax_vraisemblance[iter,1], sqrt( sum( (abs(Resultat_sansPrior$matrice_argmax_vraisemblance[iter,2:3]) - vrai_noyau$parametres$longueur)^2 )) )
        
        Argoptima_locaux_avecPriorI <- matrix(NA, nrow=nb_iterations, ncol= 2)
        for ( iter in 1:nb_iterations ) Argoptima_locaux_avecPriorI[iter,]<- c( Resultat_avecPriorI$matrice_argmax_vraisemblance[iter,1], sqrt( sum( (abs(Resultat_avecPriorI$matrice_argmax_vraisemblance[iter,2:3]) - vrai_noyau$parametres$longueur)^2 )) )
                                         
        Argoptima_locaux_avecPriorII <- matrix(NA, nrow=nb_iterations, ncol= 2)
        for ( iter in 1:nb_iterations ) Argoptima_locaux_avecPriorII[iter,]<- c( Resultat_avecPriorII$matrice_argmax_vraisemblance[iter,1], sqrt( sum( (abs(Resultat_avecPriorII$matrice_argmax_vraisemblance[iter,2:3]) - vrai_noyau$parametres$longueur)^2 )) )

        Argoptima_locaux_avecPriorIII <- matrix(NA, nrow=nb_iterations, ncol= 2)
        for ( iter in 1:nb_iterations ) Argoptima_locaux_avecPriorIII[iter,]<- c( Resultat_avecPriorIII$matrice_argmax_vraisemblance[iter,1], sqrt( sum( (abs(Resultat_avecPriorIII$matrice_argmax_vraisemblance[iter,2:3]) - vrai_noyau$parametres$longueur)^2 )) )

        
        list( Vraie_valeur= matrix(Vraie_valeur,nrow=1), Argoptima_locaux_sansPrior= Argoptima_locaux_sansPrior, Argoptima_locaux_avecPriorI= Argoptima_locaux_avecPriorI, Argoptima_locaux_avecPriorII= Argoptima_locaux_avecPriorII, Argoptima_locaux_avecPriorIII= Argoptima_locaux_avecPriorIII,  Numero_argoptimum_global_sansPrior= which.min(Resultat_sansPrior$vect_pseudo_vraisemblance), Numero_argoptimum_global_avecPriorI= which.min(Resultat_avecPriorI$vect_pseudo_vraisemblance), Numero_argoptimum_global_avecPriorII= which.min(Resultat_avecPriorII$vect_pseudo_vraisemblance) , Numero_argoptimum_global_avecPriorIII= which.min(Resultat_avecPriorIII$vect_pseudo_vraisemblance))

                                     }

x_connus_nuage<- generePointsConnus( nombre_points = 42, minima = c(0.1,0.1), maxima = c(0.9,0.9) ,lhs= TRUE)
x_connus_croix <- genereCroix( matrice_centres = matrix(x_connus_nuage[1,],nrow=1), vecteur_nombre_points_par_dimension = 4, vecteur_espacements = c(0.05,0.05) )
x_connus_complet <- rbind(x_connus_nuage, x_connus_croix)



res_argoptima <- represente_argoptima( nb_iterations= 5, x_connus= x_connus_complet, vrai_noyau= c( creeMaternIsotrope( variance= 1, longueur= c(0.3 , 0.4), regularite = 1.5) , boxCox= 0.5) , intercept_moyenne = 0, pente_moyenne = 0)

plot(res_argoptima$Argoptima_locaux_avecPriorI[- res_argoptima$Numero_argoptimum_global_avecPriorI,], main="50 points | var = 1 | correl = (0.3 , 0.4)", xlim=c(0,1), ylim=c(0,1), pch=3, lwd=2,xlab="lambda", ylab="norme erreur longueur de corrélation")                                
points(res_argoptima$Argoptima_locaux_avecPriorII[-res_argoptima$Numero_argoptimum_global_avecPriorII,],pch=3, lwd=2,col=3)   
points(res_argoptima$Argoptima_locaux_avecPriorIII[-res_argoptima$Numero_argoptimum_global_avecPriorII,],pch=3, lwd=2,col=5)                                      
points(res_argoptima$Argoptima_locaux_sansPrior[-res_argoptima$Numero_argoptimum_global_sansPrior,],pch=3, lwd=2,col=2)
points(res_argoptima$Vraie_valeur, pch=5,lwd=2,col=4)
points(matrix(res_argoptima$Argoptima_locaux_sansPrior[res_argoptima$Numero_argoptimum_global_sansPrior,] ,nrow=1) , pch=5,lwd=2, col=2)
points(matrix(res_argoptima$Argoptima_locaux_avecPriorI[res_argoptima$Numero_argoptimum_global_avecPriorI,] ,nrow=1) , pch=5,lwd=2, col=1)
points(matrix(res_argoptima$Argoptima_locaux_avecPriorII[res_argoptima$Numero_argoptimum_global_avecPriorII,] ,nrow=1) , pch=5,lwd=2, col=5)
legend("topright", legend= c("Paramètres réels","Avec prior ordre 1","Optimum", "Avec prior ordre 2","Optimum","Avec prior ordre 3","Optimum","Sans prior", "Optimum"), col=c(4,1,1,3,3,2,2,5,5),pch=c(5,3,5,3,5,3,5,3,5 ), pt.lwd=2 )


## dev.copy(pdf,"50Clambda0.5_var_1_correl0.3_0.4.pdf")
## dev.off()
## dev.copy(postscript,"50Clambda0.5_var_1_correl0.3_0.4.ps")
## dev.off()



















## Sortie : matrice des valeurs du (pseudo)prior (1er ordre) selon une coupe en lambda (lignes) et en l'un des paramètres régulant la moyenne (colonnes)
## Dépendances:	creeMaternIsotrope [Matern.r]
##		commande_calculePseudoPriorI/II_kriUniv_boxCox --> ... [Krige.r]
 coupePseudoPrior_kriUniv_boxCox <- function( x_connus, beta1= NULL, beta2= NULL, liste_fonctions, variance, longueur_correlation, regularite, lambda , avecPrior) 
 {
 	if(min(length(beta1), length(beta2))>1) print("erreur : a1 ou beta2 doit être un scalaire !")
 	 res<- matrix(0, nrow= length(lambda), ncol= max(length(beta1), length(beta2)) )

 	if(length(beta1)>1)
 	{print("beta1 est long")
 		for(l in 1:length(lambda))
 		{
 			for(b in 1:length(beta1))
 			{
				
 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation, regularite= regularite)

 				noyau$boxCox <- lambda[l]

if(avecPrior==1) res[l,b] <- commande_calculePseudoPriorI_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1[b],beta2) )


if(avecPrior==2) res[l,b] <- commande_calculePseudoPriorII_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1[b],beta2) )


 ##                                if(l==1 || l==length(lambda))
 ##                                    {
##                                         print(paste("beta1[b] =", beta1[b]))
                                         print(paste("lambda[l]=", lambda[l]))
                                         print(res[l,b])
  ##                                   }
                                        
 			}
 		}
 	}

 	if(length(beta2)>1)
 	{print("beta2 est long")
 		for(l in 1:length(lambda))
 		{
 			for(b in 1:length(beta2))
 			{
 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation, regularite= regularite)
 				noyau$boxCox <- lambda[l]
 				
 				if(avecPrior==1) res[l,b] <- commande_calculePseudoPriorI_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1,beta2[b]) )

 				if(avecPrior==2) res[l,b] <- commande_calculePseudoPriorII_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1,beta2[b]) )
 			}
 		}
 	}

 	res
 }




## Sortie : matrice des valeurs du (pseudo)prior (1er ordre) selon une coupe en lambda (lignes) et la longueur de correlation (colonnes) - le modele est monodimensionnel.
## Dépendances:	creeMaternIsotrope [Matern.r]
##		commande_calculePseudoPriorI/II_kriUniv_boxCox --> ... [Krige.r]
 coupePseudoPrior_lambda_correl_kriUniv_boxCox <- function( x_connus, beta1, beta2, liste_fonctions, variance, longueur_correlation, regularite, lambda, avecPrior ) 
 {
 	 res<- matrix(0, nrow= length(lambda), ncol= length(longueur_correlation) )

 		for(l in 1:length(lambda))
 		{
 			for(lc in 1:length(longueur_correlation))
 			{
				
 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation[lc], regularite= regularite)

 				noyau$boxCox <- lambda[l]

 				if(avecPrior==1) res[l,lc] <- commande_calculePseudoPriorI_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1,beta2) )


 				if(avecPrior==2) res[l,lc] <- commande_calculePseudoPriorII_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1,beta2) )

                print(paste("lambda[l]=", lambda[l]))

                                        
 			}
 		}

  	res
 }





## Sortie : matrice des valeurs de la (pseudo)loi jointe selon une coupe en lambda (lignes) et en l'un des paramètres régulant la moyenne (colonnes)
## Dépendances:	creeMaternIsotrope [Matern.r]
##		commande_calculePseudoPriorI/II_kriUniv_boxCox --> ... [Krige.r]
 coupePseudoLoiJointe_kriUniv_boxCox <- function( x_connus, y_connus, beta1= NULL, beta2= NULL, liste_fonctions, variance, longueur_correlation, regularite, lambda, avecPrior ) 
 {
 	if(min(length(beta1), length(beta2))>1) print("erreur : a1 ou beta2 doit être un scalaire !")
 	 res<- matrix(0, nrow= length(lambda), ncol= max(length(beta1), length(beta2)) )

 	if(length(beta1)>1)
 	{print("beta1 est long")
 		for(l in 1:length(lambda))
 		{
 			for(b in 1:length(beta1))
 			{
				
 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation, regularite= regularite)

 				noyau$boxCox <- lambda[l]

 				res[l,b] <- calculeFonctionAMinimiser_BoxCox( x_connus= x_connus, y_connus = y_connus, noyau= noyau, liste_fonctions= liste_fonctions, boxCox= noyau$boxCox, beta= c(beta1[b],beta2) , avecPrior= avecPrior)
 ##                                if(l==1 || l==length(lambda))
 ##                                    {
##                                         print(paste("beta1[b] =", beta1[b]))
                                         print(paste("lambda[l]=", lambda[l]))
                                         print(res[l,b])
  ##                                   }
                                        
 			}
 		}
 	}

 	if(length(beta2)>1)
 	{print("beta2 est long")
 		for(l in 1:length(lambda))
 		{
 			for(b in 1:length(beta2))
 			{
 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation, regularite= regularite)
 				noyau$boxCox <- lambda[l]
                                
 				res[l,b] <- calculeFonctionAMinimiser_BoxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, boxCox= noyau$boxCox, beta= c(beta1,beta2[b]) , avecPrior= avecPrior)
 			}
 		}
 	}

 	res
 }




## Sortie : matrice des valeurs de la (pseudo)loi jointe selon une coupe en lambda (lignes) et en la longueur de correlation (colonnes)
## Dépendances:	creeMaternIsotrope [Matern.r]
##		commande_calculePseudoPriorI/II_kriUniv_boxCox --> ... [Krige.r]
 coupePseudoLoiJointe_lambda_correl_kriUniv_boxCox <- function( x_connus, y_connus, beta1, beta2, liste_fonctions, variance, longueur_correlation, regularite, lambda, avecPrior ) 
 {
 	 res<- matrix(0, nrow= length(lambda), ncol= length(longueur_correlation) )

 		for(l in 1:length(lambda))
 		{
 			for(lc in 1:length(longueur_correlation))
 			{
				
 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation[lc], regularite= regularite)

 				noyau$boxCox <- lambda[l]
 				
 				print(paste("lambda[l]=", lambda[l]))
 				
 				res[l,lc] <- calculeFonctionAMinimiser_BoxCox( x_connus= x_connus, y_connus = y_connus, noyau= noyau, liste_fonctions= liste_fonctions, boxCox= noyau$boxCox, beta= c(beta1,beta2) , avecPrior= avecPrior)
                                        
 			}
 		}
 
 	res
 }


xxx <- generePointsConnus(5,lhs=TRUE)


cc0 <- 10:100*0.01
vvv0 <- -100:100*0.2  
lll0 <- 1:100 * 0.02

#lll0 <- lll0[-101] ## retrait de lambda=0 --> engendrerait des NA

 mat0 <- coupePseudoPrior_kriUniv_boxCox(x_connus = xxx, beta1= vvv0 , beta2=0, liste_fonctions = list(function (t) 1 , function(t) t), variance=1, longueur_correlation = 0.3, regularite = 1.5, lambda = lll0 , avecPrior=1)

image2D(x=lll0, y=vvv0, z=mat0, xlab= "lambda (paramètre de Box-Cox)", ylab= "intercept de la moyenne", contour=TRUE,  main= "Loi a priori 1 (log) | var = 1 | correl = 0.3 | pente = 0")


 mat00 <- coupePseudoPrior_lambda_correl_kriUniv_boxCox(x_connus = xxx, beta1= 10, beta2=0, liste_fonctions = list(function (t) 1 , function(t) t), variance=1, longueur_correlation = cc0, regularite = 1.5, lambda = lll0, avecPrior = 2 )

image2D(x=lll0, y=cc0, z=mat00, xlab= "lambda (paramètre de Box-Cox)", ylab= "longueur de correlation", contour=TRUE,  main= "Loi a priori 2 (log) | Box-Cox = 1 | var = 1 | correl = 0.3 | moyenne = 10")



dev.copy(pdf, "2eOrdre_var1_moyenne10.pdf")
dev.off()
dev.copy(postscript, "2eOrdre_var1_moyenne10.ps")
dev.off()


 yyy <- BoxCox_reciproque (  apply( xxx, 1, esperanceAffine(a= 10, b= 0) ) + creeEchantillonNormal(x_connus= xxx, noyau= creeMaternIsotrope( variance= 1, longueur= 0.3, regularite= 1.5) )   , lambda= 1 )


 mat1 <- coupePseudoLoiJointe_kriUniv_boxCox(x_connus = xxx, y_connus = yyy, beta1= vvv0 , beta2=0, liste_fonctions = list(function (t) 1 , function(t) t), variance=1, longueur_correlation = 0.3, regularite = 1.5, lambda = lll0, avecPrior = 1 )
 


image2D(x=lll0, y=vvv0, z=-mat1, xlab= "lambda (paramètre de Box-Cox)", ylab= "intercept de la moyenne", contour=TRUE,  main= "Loi jointe 1 (log) | Box-Cox = 1 | var = 1 | correl = 0.3 | moyenne = 10")


 mat10 <- coupePseudoLoiJointe_lambda_correl_kriUniv_boxCox(x_connus = xxx, y_connus = yyy, beta1= 10 , beta2=0, liste_fonctions = list(function (t) 1 , function(t) t), variance=1, longueur_correlation = cc0, regularite = 1.5, lambda = lll0, avecPrior = 2 )

image2D(x=lll0, y=cc0, z=-mat10, xlab= "lambda (paramètre de Box-Cox)", ylab= "longueur de correlation", contour=TRUE,  main= "Loi jointe 2 (log) | Box-Cox = 1 | var = 1 | correl = 0.3 | moyenne = 10")

dev.copy(pdf, "C_LJ2eOrdre_box1_var1_moyenne10.pdf")
dev.off()
dev.copy(postscript, "C_LJ2eOrdre_box1_var1_moyenne10.ps")
dev.off()











                                 
## mat0_trafique <- mat0
## mat0_trafique[mat0 < -log(5000)] <- -log(5000)
## image2D(x=lll0, y=vvv0, z=exp(-mat0_trafique), xlab= "lambda (paramètre de Box-Cox)", ylab= "intercept de la moyenne", contour=TRUE,  main= "Loi a priori - lambda = 0 retiré - sup artificiel = 5000")


## dev.copy(pdf,"lambda_negatif.pdf")
## dev.off()

## dev.copy(postscript,"lambda_negatif.ps")
## dev.off()



## image2D(x=lll0, y=vvv0, z=mat0_trafique, xlab= "lambda (paramètre de Box-Cox)", ylab= "intercept de la moyenne", contour=TRUE,  main= "Loi a priori ( - log) - lambda = 0 retiré - inf artificiel = -log(5000)")


## dev.copy(pdf,"lambda_negatif_log.pdf")
## dev.off()

## dev.copy(postscript,"lambda_negatif_log.ps")
## dev.off()


## vvv0 <- 1:100*0.1  
## lll0 <- 1:100 * 0.01

## lll0 <- lll0[-301] ## retrait de lambda=0 --> engendrerait des NA

## mat0 <- coupePseudoPrior_kriUniv_boxCox(x_connus = x_complet, beta1= vvv0 , beta2=0, liste_fonctions = list(function (t) 1 , function(t) t[1]), variance=1, longueur_correlation = c(0.3,0.4), regularite = 1.5, lambda = lll0 )


## image2D(x=lll0, y=vvv0, z=mat0, xlab= "lambda (paramètre de Box-Cox)", ylab= "intercept de la moyenne", contour=TRUE,  main= "Loi a priori ( - log) | lambda = 0 retiré | correl = 0.1 | pente = 0")


## dev.copy(pdf,"lambda_-3_3_log_correl_0.1_pente_0.pdf")
## dev.off()

## dev.copy(postscript,"lambda_-3_3_log_correl_0.1_pente_0.ps")
## dev.off()


## ##library(plot3D)
## ##image2D(x=vvv, y=vvv, z=mat1, xlab= "lambda (paramètre de Box-Cox)", ylab= "intercept de la moyenne", contour=TRUE, clim=c(-8,0), main= "Loi a priori ( - logarithme de)")









