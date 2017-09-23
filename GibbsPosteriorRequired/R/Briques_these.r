# Briques fondamentales pour les tests



# nombre_points est un entier, vecteurMinMax un vecteur à 2 composantes
posePointsUnif <- function ( nombre_points, vecteurMinMax)
{
	runif(nombre_points, min=vecteurMinMax[1], max= vecteurMinMax[2])
}


# nombre_points est le vecteur contenant le nombre de points dans chaque amas. Exemple: (80,10,10) --> 80 points répartis uniformément sur le domaine, 10 dans un 1er amas et 10 dans un 2e amas. Les centres des amas sont placés selon la loi uniforme sur le domaine (Mise à jour these 2015/10/16 : si point_central n'est pas NULL, le centre du dernier amas est point_central --> pas tiré uniformément sur le domaine.)
# minima et maxima st les vecteurs indiquant les bords du domaine. Exemple: minima=(0,1,-1) et maxima=(1,2,3). Alors le domaine est [0,1] x [1,2] x [-1,3]
# rayon est le rayon (en norme infinie) des amas. Exemple : si un amas a pour centre (1,2,0) et que le rayon est 0.1, alors l'amas recouvre la zone [0.9,1.1] x [1.9,2.1] x [-0.1,0.1].
# rayon n'est utilisé que si nombre_points est un vecteur de longueur >=2. Voilà pourquoi la valeur par défaut est NULL : si la longueur du vecteur nombre_points est >=2, cela siginifie qu'il y a au moins un amas, et que donc le rayon doit être spécifié.
# Dépendances : posePointsUnif
# Sortie : matrice dont chaque ligne représente les corrdonnées d'un point du plan d'expérience
# ATTENTION MAJ 2015/10/16 : si is.null(point_central)==FALSE, alors nombre_points DOIT être un vecteur d'au moins 2 entrées. La dernière entrée désigne le nombre de points amassés autour de point_central
generePointsConnus <- function ( nombre_points, minima=0, maxima=1, rayon=NULL, point_central= NULL, lhs= FALSE)
{
	if(lhs)
	{
		if(length(nombre_points) != 1) print("ATTENTION ! Avec geneticLHS, pas question d'amas de points !")		
		Points_connus <- geneticLHS( n= nombre_points, k= length(minima) , pop=1000, gen=10)		
        	for( num_dim in 1:length(minima) ) Points_connus[,num_dim] <- minima[num_dim]+ ( maxima[num_dim]- minima[num_dim]) * Points_connus[,num_dim] 
	}

	else
	{
			if ( length(minima) != length(maxima) ) print(" Erreur : A chaque minimum doit correspondre un maximum. ")


			Nbpts_glob_repartis <- nombre_points[1]
			Nb_centres_amas <- length(nombre_points) -1
			MinMax <- rbind ( minima, maxima )


			Points_connus <- apply ( MinMax , 2 , posePointsUnif , nombre_points=Nbpts_glob_repartis + Nb_centres_amas)

			if(is.null( dim(Points_connus) ) ) Points_connus <- matrix(Points_connus, nrow=1 ) # Points_connus doit être une matrice.

			if(!is.null(point_central)) 
			{
				Points_connus[nrow(Points_connus),] <- point_central
			}



			if(length(nombre_points) > 1)
			{
				Numeros_pts_centraux <- 1:Nb_centres_amas + Nbpts_glob_repartis
			  	for(numero_point_central in Numeros_pts_centraux )
					{	
						numero_amas <- numero_point_central - Nbpts_glob_repartis + 1
						MinMax <- rbind( Points_connus[numero_point_central,] - rayon , Points_connus[numero_point_central,] + rayon )
			      			Points_connus <- rbind(Points_connus , apply ( MinMax , 2 , posePointsUnif , nombre_points=nombre_points[numero_amas] - 1 ) )
					}
			}
	}

	Points_connus
}

# Sortie : matrices contenant les coordonnees (en ligne) de points formant des croix autour des points (lignes) de matrice_centres
# Dépendances : aucune
genereCroix <- function ( matrice_centres, vecteur_nombre_points_par_dimension, vecteur_espacements )
    {
        Nombre_dimensions <- ncol( matrice_centres )
        Nombre_centres <- length( vecteur_nombre_points_par_dimension )

##        if( length(vecteur_espacements) != Nombre_centres ) print("Combien de centres ?")
        
        Matrices_points <- matrix(0, nrow= sum(vecteur_nombre_points_par_dimension)*Nombre_dimensions , ncol= Nombre_dimensions )

        Numero_point_global <- 0
        
        for ( Numero_centre in 1:Nombre_centres)
            {
                print(Nombre_centres)
                Nombre_points_par_dimension <- vecteur_nombre_points_par_dimension[Numero_centre]
##                Espacement <- vecteur_espacements[Numero_centre]
                
                for( Numero_point in 0:(Nombre_points_par_dimension-1) )
                    {

                        Signe <- (-1)^Numero_point

                        for ( Numero_dimension in 1:Nombre_dimensions )
                            {
                                Ecart_centre <- (1 + (Numero_point %/% 2) ) * vecteur_espacements[Numero_dimension]
                                Numero_point_global <- Numero_point_global + 1
                                Matrices_points[Numero_point_global, ] <- matrice_centres[Numero_centre, ]
                             Matrices_points[Numero_point_global, Numero_dimension ] <-  Matrices_points[Numero_point_global, Numero_dimension ] + Signe * Ecart_centre
 ##                            print(Signe*Ecart_centre)
 ##                            print( Matrices_points )
                            }
##                        Ecart_centre <- (Numero_point %/% 2) * vecteur_espacements
##                        Matrices_points[Numero_point_global, ] <- matrice_centres[Numero_centre, ] + Signe * Ecart_centre
                    }
            }

        Matrices_points

    }


# Dépendances : creeMatriceCovarianceConditionnee [Krige.r]
# Sortie : valeurs du processus gaussien spécifié aux points (lignes) de x_connus
creeEchantillonNormal <- function ( x_connus , noyau , pepite=0 )
{

	# Préliminaires
	# x_connus doit être une matrice, y_connus un vecteur
	if ( is.null(dim(x_connus)) ) x_connus <- matrix(x_connus)
#print(x_connus)
	y_connus <- rnorm ( nrow ( x_connus ) )
	Connus_covariance_matrice <- creeMatriceCovarianceConditionnee ( x1=x_connus , x2=x_connus,  noyau=noyau, y=y_connus, pepite=pepite )
#print(dim(Connus_covariance_matrice))
#print(length(y_connus))
	t( chol ( Connus_covariance_matrice ) ) %*% y_connus
}


# Sortie : valeur A^2 du test d'Anderson-Darling effectué sur les entrées du vecteur valeurs.
Anderson_Darling <- function ( valeurs )
{
	valeurs <- sort(valeurs)
	coeffs <- 2 * seq(1, length(valeurs)) - 1
	- length(valeurs ) - mean ( coeffs * log ( pnorm(valeurs) * ( 1 - pnorm ( rev ( valeurs ) ) ) ) )
}


# Equivalent de Anderson_Darling pour le test de Cramér - von Mises.
Cramer_vonMises <- function (valeurs)
{
	valeurs <- sort ( valeurs )
	coeffs <- ( 2 * seq(1, length(valeurs)) - 1 ) / ( 2 * length(valeurs) )
	1/(12*length(valeurs)) * sum ( ( coeffs - pnorm(valeurs) )^2  )	
}
