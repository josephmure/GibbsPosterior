
# Validation directe

# Sortie : liste contenant (notations mémoire) 1) Matrice_auxiliaire = H ( t(H) K^{-1} H )^{-1} t(H) K^{-1} et 2) Covariance_valeurs= K-tilde = K - H ( t(H) K^{-1} H )^{-1} t(H) (en krigeage simple, 1) 0 et 2) K )
# Remarque : Il est ici tenu compte d'une éventuelle pépite.
# Dépendances : creeMatriceCovarianceConditionnee[Krige.r]
correlationValeurs <- function ( x_connus, x_nouveaux, y_connus, noyau, pepite=0, liste_fonctions = NULL ) 
{


	# Préliminaires
	# x_connus et x_nouveaux doivent être des matrices, y_connus un vecteur
	if ( is.null(dim(x_connus)) ) x_connus <- matrix(x_connus)
	dim(y_connus) <- NULL


	# Matrice des covariances entre points connus
	Connus_covariance_matrice <- creeMatriceCovariance ( x1=x_connus , x2=x_nouveaux,  noyau=noyau)


	
	# Créé Matrice_correlante ( matrice K-tilde-moins de Bachoc )
	if ( is.null(liste_fonctions ) )	# Krigeage simple
	{
		Matrice_auxiliaire <- 0 * Connus_covariance_matrice
		Covariance_valeurs <- Connus_covariance_matrice
	}

	else 					# Krigeage universel
	{

		Inv_connus_covariance_matrice <- solve(Connus_covariance_matrice )

		# Evaluation des fonctions aux points connus
		Fonctions_points_connus <- creeMatriceValeursFonctions ( matrice_points = x_connus , liste_fonctions = liste_fonctions )

		Matrice_correlante <-	Fonctions_points_connus %*% 
						solve( t( Fonctions_points_connus ) %*% Inv_connus_covariance_matrice %*% Fonctions_points_connus,
							t ( Fonctions_points_connus ) )
		
		Covariance_valeurs <-  Connus_covariance_matrice - Matrice_correlante
		Matrice_auxiliaire <- Matrice_correlante %*% Inv_connus_covariance_matrice
	}

	# Inversion de la matrice des coeffs diagonauxde Matrice_correlante

	# Détail technique : diag ( 1 2 ) = (1 4 )	diag( diag ( 1 2 ) ) = (1 0 )
	#			    3 4 			     3 4	0 4

#	Inv_diag <- solve ( diag ( diag ( Matrice_correlante ) ) ) 


	list ( Matrice_auxiliaire= Matrice_auxiliaire, Covariance_valeurs = Covariance_valeurs )

}






# Calcul des observations (recentrées / décorrélées)

# Sortie : liste contenant 1) y_centres= vecteur des observations recentrées 2) y_decorreles=  vecteur des observations décorrélées (en krigeage simple 1) est le vecteur des observations, puisqu'elles sont déjà centrées).
# Dépendances : correlationValeurs ->  creeMatriceCovarianceConditionnee[Krige.r]
valeurs <- function ( x_connus, y_connus , noyau, pepite=0, liste_fonctions = NULL ) 
{
	Infos_corr<-correlationValeurs(x_connus= x_connus, x_nouveaux= x_connus, y_connus = y_connus ,noyau= noyau,pepite= pepite,liste_fonctions= liste_fonctions)
	y_centres <- y_connus - Infos_corr$Matrice_auxiliaire %*% y_connus



	Nb_VPnonnulles <- length(y_connus) - length( liste_fonctions ) # Autant de vp nulles que d'éléments de liste_fonctions
        
	stockeSVD <-  svd(Infos_corr$Covariance_valeurs ,  nv= 0) # U == stockeSVD$u, D == stockeSVD$d


# Début de décorrélation : prémultiplication de y_centres par t(U)
	y_decorreles_zeros <- colSums ( stockeSVD$u * c(y_centres) )   # L'inverse de stockeSVD$u est t(stockesSVD$u). colSums(M*v) == t(M) %*% v, à condition que v soit du type VECTEUR, et non MATRICE A UNE COLONNE

# Fin de décorrélation : retrait des composantes constantes-nulles et réduction des autres
	VP_nonnulles <- stockeSVD$d[1:Nb_VPnonnulles] # vecteur des VP non nulles de K-tilde (K en krigeage simple)
	y_decorreles_VP <- y_decorreles_zeros[1:Nb_VPnonnulles] # Retrait des composantes contstantes-nulles
	y_decorreles <- solve( diag( sqrt ( VP_nonnulles ) ),  y_decorreles_VP ) # Réduction des composantes



        list ( y_centres = y_centres , y_decorreles = y_decorreles )
}



## Décorrélation par Cholesky des erreurs de krigeage sur un ensemble de test
## x_connus : matrice des points de l'ensemble d'apprentissage ( 1 ligne = 1 point )
## x_test : matrice des points de l'ensemble de test (1 ligne = 1 point)
## y_connus : vecteur des valeurs prises par le processus gaussien aux points de l'ensemble d'apprentissage
## y_test : vecteur des valeurs prises par le processus gaussien aux points de l'ensemble de test
## noyau : noyau du processus gaussien
## pepite : A PRIORI INUTILE
## liste_fonctions : liste contenant les fonctions de régression
## Sortie : erreurs décorrélées par Cholesky sur l'ensemble de test x_test sachant les valeurs y_connus prises sur l'ensemble d'apprentissage x_connus
erreursCholesky <- function ( x_connus, x_test, y_connus, y_test, noyau, pepite=0, liste_fonctions= NULL )
    {
        y_test_krige <- krige ( x_connus= x_connus, y_connus= y_connus, x_nouveaux= x_test, noyau= noyau, pepite= pepite, liste_fonctions= liste_fonctions )[1,]

        test_erreur_krige <- y_test_krige - y_test

        Variance_test <- correlationValeurs(x_connus= x_test, x_nouveaux= x_test, y_connus = 0 ,noyau= noyau,pepite= pepite,liste_fonctions= NULL)$Covariance_valeurs
	Variance_connus <- correlationValeurs(x_connus= x_connus, x_nouveaux= x_connus, y_connus = 0 ,noyau= noyau,pepite= pepite,liste_fonctions= NULL)$Covariance_valeurs
        Covariance_connus_test <- correlationValeurs(x_connus= x_connus, x_nouveaux= x_test, y_connus = 0 ,noyau= noyau,pepite= pepite,liste_fonctions= NULL)$Covariance_valeurs

        Variance_test_conditionnelle <- Variance_test - t(Covariance_connus_test) %*% solve( Variance_connus, Covariance_connus_test  )

        test_erreur_decorrelee <- solve( t( chol( Variance_test_conditionnelle ) ) , test_erreur_krige )

        test_erreur_decorrelee
    }


