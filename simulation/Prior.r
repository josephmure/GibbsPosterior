## Calcul du prior quasi-Jeffreys

## Sortie : nablaG(H beta) quelle que soit la transformation utilisée
## Dépendances : aucune
## Attention : fonctionA et fonctionBprime sont des fonctions R --> R. En les exécutant sur un vecteur, on les applique en fait à chaque composante dudit vecteur.
calculeNablaG_moyenne <- function ( moyenne, fonctionA, fonctionBprime, matriceCorrelation )
{
	# Moyenne <- Fonctions_points_connus %*% beta

##    print(fonctionA)
                                        #	print(fonctionBprime)
##    print(dim(fonctionA(moyenne)))
	2 * matriceCorrelation %*% fonctionA(moyenne) + fonctionBprime(moyenne)

	
}

## Sortie : nablaG(H beta) dans le cas particulier de la transformation de Box-Cox (utilise la programmation générique de calculeNablaG_moyenne)
## Dépendance : calculeNablaG_moyenne
creeNablaG_moyenne_boxCox <- function( moyenne, lambda, matriceCorrelation)
{
	fonctionA <- function(t) 1/lambda * ( log( abs( 1 + lambda * t) ) * t + 1/lambda * log( abs( 1 + lambda * t) ) - t )
	fonctionBprime <- function(t) 1/ ( 1 + lambda * t )

	calculeNablaG_moyenne( moyenne= moyenne, fonctionA= fonctionA, fonctionBprime= fonctionBprime, matriceCorrelation= matriceCorrelation )
}

## Sortie : log(prior 1er ordre) quelle que soit la transformation utilisée
## Dépendances : aucune
calculePseudoPrior_kriUniv <- function( variance, matriceCorrelation, Fonctions_points_connus, nablaG_moyenne)
{
#	variance <- noyau$parametres$variance
	t_Fonctions_points_connus <- t(Fonctions_points_connus)
	nablaG_moyenne <- matrix(nablaG_moyenne) ## pour le calcul matriciel, il faut que nablaG_moyenne soit une matrice à 1 colonne et non un vecteur
	auxiliaire <- c( t(nablaG_moyenne) %*% matriceCorrelation %*% nablaG_moyenne) ## c() car auxiliare doit être vu comme un scalaire et non comme une matrice 1x1
	#print("aux")
#print(auxiliaire)
	#print(dim(nablaG_moyenne))
                                        #print(dim(t_Fonctions_points_connus))

        print(paste("terme variance = ", 	- log(variance ) ) )
        print(paste("terme lambda = ", 1/2 *  log( variance * auxiliaire ) ) ) 
        print(paste("terme beta = ", 1/2 * log( det( 1/variance * ( t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) - 1/auxiliaire * (t_Fonctions_points_connus %*% nablaG_moyenne) %*% (t(nablaG_moyenne) %*% Fonctions_points_connus) ) ) ) ) )
	- log(variance) + 1/2 * ( log( variance * auxiliaire ) + log( det( 1/variance * ( t_Fonctions_points_connus %*% solve( matriceCorrelation, Fonctions_points_connus ) - 1/auxiliaire * (t_Fonctions_points_connus %*% nablaG_moyenne) %*% (t(nablaG_moyenne) %*% Fonctions_points_connus) ) ) ) )
}


## Sortie : log(prior 1er ordre) dans le cas de la transformation de Box-Cox
## Dépendances:	calculePseudoPrior_kriUniv
##		creeMatrices -> ... [Krige.r]
commande_calculePseudoPrior_kriUniv_boxCox <- function( x_connus, noyau, liste_fonctions, beta, Matrices_utiles= NULL)
{
 if(is.null(Matrices_utiles)) Matrices_utiles <- creeMatrices( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions)
 
	moyenne <- Matrices_utiles$Fonctions_points_connus %*% beta
	nablaG_moyenne <- creeNablaG_moyenne_boxCox( moyenne= moyenne, lambda= noyau$boxCox, matriceCorrelation= Matrices_utiles$matriceCorrelation)

	calculePseudoPrior_kriUniv( variance= noyau$parametres$variance, matriceCorrelation= Matrices_utiles$matriceCorrelation, Fonctions_points_connus= Matrices_utiles$Fonctions_points_connus, nablaG_moyenne= nablaG_moyenne )
}






## Sortie : liste comportant les pseudo-vraisemblances (en tenant compte du prior) atteintes en les minima trouvees, les warnings (0= tout va bien), les jeux de param argmax pseudo-vraisemblance, les jeux de param initiaux, MAIS PAS A^2 d'Anderson-Darling par valid directe et LOO
## NB : les numeros de lignes des matrices correspondent aux numeros d'entrees des vecteurs
## NB2: les initialisations se font par tirage aléatoire de loi uniforme déterminée par bornesInitialisation. Pour eviter les minima locaux, on reessaie nb_iterations fois.
## NB3: valable uniquement avec la transformation de Box-Cox
## Dependances:	lhs::geneticLHS
##		calculeFonctionAMinimiserVectoriel_BoxCox --> ... [Krige.r] [Matern.r] [Box-Cox.r]
## ATTENTION : ici, la regularite est supposee connue
trouveArgmaxLoiJointe_BoxCox <- function( x_connus, y_connus, regularite, beta, liste_fonctions, bornesInitialisation= list( boxCox=list(minima=0,maxima=1), longueur= list(minima= c(0.001,0.001,0.001,0.001,0.001,0.001), maxima= c( 0.3, 1,2, 10, 6,6) ) ),  nb_iterations=100)
    {
        nb_dim <-   length(bornesInitialisation$longueur$maxima )
        matrice_argmax_vraisemblance <- matrix(rep(1000, nb_iterations * (nb_dim+1)), ncol=nb_dim+1)
        matrice_param_init <- matrix(rep(1000, nb_iterations * (nb_dim+1)), ncol=nb_dim+1)
        vect_pseudo_vraisemblance <- rep(1000, nb_iterations)
        vect_warning <- rep(0, nb_iterations)
	matrice_validation <- matrix( rep(1000, nb_iterations*2)  , ncol=2)


        matrice_param_init <- geneticLHS(n=nb_iterations, k= 1 + nb_dim)
        matrice_param_init[,1] <- bornesInitialisation$boxCox$minima + ( bornesInitialisation$boxCox$maxima - bornesInitialisation$boxCox$minima ) *matrice_param_init[,1]

        for( num_dim in 1:nb_dim ) matrice_param_init[,num_dim+1] <- bornesInitialisation$longueur$minima[num_dim]+ ( bornesInitialisation$longueur$maxima[num_dim]- bornesInitialisation$longueur$minima[num_dim]) *matrice_param_init[,num_dim+1]
        
        for ( iter in 1:nb_iterations )
            {
                
                OPT <- try(optim( par= matrice_param_init[iter,] , fn= calculeFonctionAMinimiserVectoriel_BoxCox, x_connus= x_connus, y_connus= y_connus, regularite= regularite, liste_fonctions= liste_fonctions, beta= beta, method="BFGS" ) )

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

        list( vect_pseudo_vraisemblance= vect_pseudo_vraisemblance, vect_warning= vect_warning, matrice_argmax_vraisemblance= matrice_argmax_vraisemblance, matrice_param_init= matrice_param_init)##, matrice_validation= matrice_validation)
    }



## Sortie : pseudo-vraisemblance des observations y_connus en fonction du noyau, de la base de fonctions d'esperance et du parametre de Box-Cox utilisés
## Dependances:	calculeFonctionAMinimiser_BoxCox --> ... [Krige.r] [Box-Cox.r]
##		creeMaternIsotrope [Matern.r]
## ATTENTION : ici, on impose a la regularite d'etre connue
calculeFonctionAMinimiserVectoriel_BoxCox<- function( vecteur_param_noyau, x_connus, y_connus, regularite, liste_fonctions, beta)
    {
        boxCox <- vecteur_param_noyau[1]
##        regularite <- Inf
        longueur <- abs( vecteur_param_noyau[-1] )
        noyau <- creeMaternIsotrope(variance=1, longueur= longueur, regularite= regularite)

        res= calculeFonctionAMinimiser_BoxCox( x_connus= x_connus, y_connus= y_connus, noyau= noyau, liste_fonctions = liste_fonctions,  boxCox= boxCox, beta= beta )
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
calculeFonctionAMinimiser_BoxCox <-  function ( x_connus, y_connus, noyau, liste_fonctions,  boxCox, beta)
    {
        Matrices <- creeMatrices ( x_connus= x_connus,  noyau= noyau, liste_fonctions= liste_fonctions )

        Num <- BoxCox(y_connus, lambda= boxCox)
        
        noyau$boxCox <- boxCox
        noyau$parametres$variance <- c(1/(nrow(x_connus) + length(liste_fonctions) + 1) * t(Num) %*% Matrices$matriceFormeQuadratique %*% Num) ## c() pour que Variance soit un scalaire, et non une matrice 1x1
        
##print(noyau)


        
        Beta_approche<- Matrices$matriceConstructionBeta %*% Num
        print(paste("Beta_approche = ", Beta_approche))
##        Beta_approche <- c(1,2)
        
##        Moyenne <- Matrices$Fonctions_points_connus %*% Beta_approche
##        nablaG_Moyenne <- creeNablaG_moyenne_boxCox( moyenne= Moyenne, lambda = boxCox, matriceCorrelation= Matrices$matriceCorrelation)


##        Variance <- c(1/(nrow(x_connus) + length(liste_fonctions) + 1) * t(Num) %*% Matrices$matriceFormeQuadratique %*% Num) ## c() pour que Variance soit un scalaire, et non une matrice 1x1

        pseudoVraisemblance_BoxCox( y_connus= y_connus, matriceCorrelation= Matrices$matriceCorrelation, matriceFormeQuadratique= Matrices$matriceFormeQuadratique, variance= noyau$parametres$variance, boxCox= boxCox, Fonctions_points_connus=Matrices$Fonctions_points_connus, nablaG_moyenne= nablaG_Moyenne) - 2 / length(y_connus) * commande_calculePseudoPrior_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= Beta_approche, Matrices_utiles= Matrices)
    }


## Sortie : pseudo-vraisemblance des observations y_connus en fonction de la matrice de corrélation  matriceCorrelation de y_connus et du paramètre lambda de Box-Cox (noté boxCox)
## Dépendances : BoxCox [Box-Cox.r]
pseudoVraisemblance_BoxCox<- function( y_connus, matriceCorrelation, matriceFormeQuadratique,variance,  boxCox , Fonctions_points_connus, nablaG_moyenne)
    {
        nbPointsConnus <- nrow(matriceCorrelation) # = ncol(matriceCorrelation) car matriceCorrelation est carrée
        num <-  BoxCox(y_connus, lambda= boxCox)

       
        res <- 1/nbPointsConnus *  log(det(matriceCorrelation)) +log( 1/ nbPointsConnus  * t( num ) %*% matriceFormeQuadratique %*% num  ) - 2* (boxCox-1) / nbPointsConnus * sum(log(abs(y_connus)) )

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
                
                matriceFormeQuadratique <- inverseMatriceCorrelation
		Fonctions_points_connus <- NULL
                matriceConstructionBeta <- NULL

            }

        else
            {
                Fonctions_points_connus <- creeMatriceValeursFonctions ( matrice_points = x_connus , liste_fonctions = liste_fonctions )

                aux <- solve( t( Fonctions_points_connus ) %*% inverseMatriceCorrelation %*% Fonctions_points_connus )

                matriceConstructionBeta <- aux %*% t( Fonctions_points_connus ) %*% inverseMatriceCorrelation
                
                matriceFormeQuadratique <- inverseMatriceCorrelation - inverseMatriceCorrelation %*% Fonctions_points_connus %*% matriceConstructionBeta

            }

        list( matriceCorrelation= matriceCorrelation, matriceFormeQuadratique= matriceFormeQuadratique, Fonctions_points_connus= Fonctions_points_connus , matriceConstructionBeta= matriceConstructionBeta)
    }





## x_connus_nuage<- generePointsConnus( nombre_points = 25, minima = c(0.1,0.1), maxima = c(0.9,0.9) )
## x_connus_croix <- genereCroix( matrice_centres <- matrix(x_connus_nuage[1:2,],nrow=2), vecteur_nombre_points_par_dimension= c(4,4), vecteur_espacements= c(0.03, 0.03) )
## x_test_nuage <- generePointsConnus( nombre_points = c(25), minima = c(0,0), maxima = c(1,1) )
## x_test_croix <- genereCroix(  matrice_centres <- matrix(x_test_nuage[1,],nrow=1), vecteur_nombre_points_par_dimension=6, vecteur_espacements= c(0.03, 0.03) )


## x_complet <- rbind(x_connus_nuage, x_connus_croix)##, x_test_nuage, x_test_croix)

## y_complet <- apply( x_complet, 1, esperanceAffine(1,2) ) + creeEchantillonNormal(x_connus= x_complet, noyau= creeMaternIsotrope(variance=1, longueur= c(0.3,0.4), regularite= 1.5) )

## y_complet_boxCoxise <- BoxCox_reciproque(y_complet, lambda = 0.5)


RESULTAT_controle_sansPrior<- trouveArgmaxLoiJointe_BoxCox(x_connus = x_complet, y_connus = y_complet_boxCoxise, regularite = 1.5, beta= c(1,2), liste_fonctions = list(function(t) 1, function(t) t[1] ), nb_iterations = 2, bornesInitialisation = list( boxCox= list(minima=0, maxima=5), longueur= list(minima= c(0,0), maxima=c(1,1) ) ) )



































## ## Sortie : matrice des valeurs du prior (1er ordre) selon une coupe en lambda (lignes) et en l'un des paramètres régulant la moyenne (colonnes)
## ## Dépendances:	creeMaternIsotrope [Matern.r]
## ##			commande_calculePseudoPrior_kriUniv_boxCox --> ... [Krige.r]
## coupePseudoPrior_kriUniv_boxCox <- function( x_connus, beta1= NULL, beta2= NULL, liste_fonctions, variance, longueur_correlation, regularite, lambda )
## {
## 	if(min(length(beta1), length(beta2))>1) print("Erreur : beta1 ou beta2 doit être un scalaire !")
## 	res <- matrix(0, nrow= length(lambda), ncol= max(length(beta1), length(beta2)) )

## 	if(length(beta1)>1)
## 	{print("beta1 est long")
## 		for(l in 1:length(lambda))
## 		{
## 			for(b in 1:length(beta1))
## 			{
				
## 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation, regularite= regularite)

## 				noyau$boxCox <- lambda[l]

## 				res[l,b] <- commande_calculePseudoPrior_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1[b],beta2) )
## ##                                if(l==1 || l==length(lambda))
## ##                                    {
##                                         print(paste("beta1[b] =", beta1[b]))
##                                         print(paste("lambda[l]=", lambda[l]))
##                                         print(res[l,b])
##  ##                                   }
                                        
## 			}
## 		}
## 	}

## 	if(length(beta2)>1)
## 	{print("beta2 est long")
## 		for(l in 1:length(lambda))
## 		{
## 			for(b in 1:length(beta2))
## 			{
## 				noyau <- creeMaternIsotrope( variance= variance, longueur= longueur_correlation, regularite= regularite)
## 				noyau$boxCox <- lambda[l]
## 				res[l,b] <- commande_calculePseudoPrior_kriUniv_boxCox( x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions, beta= c(beta1,beta2[b]) )
## 			}
## 		}
## 	}

## 	res
## }








## xxx <- generePointsConnus(5)



## vvv0 <-0:100*0.1  
## lll0 <- -100:100 * 0.01

## lll0 <- lll0[-101] ## retrait de lambda=0 --> engendrerait des NA

## mat0 <- coupePseudoPrior_kriUniv_boxCox(x_connus = xxx, beta1= vvv0 , beta2=1, liste_fonctions = list(function (t) 1 , function(t) t), variance=1, longueur_correlation = 0.3, regularite = 1.5, lambda = lll0 )

## image2D(x=lll0, y=vvv0, z=exp(-mat0), xlab= "lambda (paramètre de Box-Cox)", ylab= "intercept de la moyenne", contour=TRUE,  main= "Loi a priori")

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









