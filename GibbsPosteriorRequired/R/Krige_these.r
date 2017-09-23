# Krigeage universel fréquentiste

####################################################################################################################

# Préparations : 

# En krigeage universel, on a besoin d'une base de fonctions qui engendrent l'espace
# vectoriel dans lequel vit la moyenne du processus gaussien. On a ensuite besoin de
# stocker sous forme matricielle les valeurs des fonctions de la base aux points qui
# nous intéressent.

# A chaque point (ie à chaque ligne de "matrice_points"), applique "fonction"
# Sortie : vecteur dont chaque entrée est le résultat de "fonction" appliqué à un point (ligne de "matrice_points")
appliqueFonction <- function ( matrice_points , fonction )
{
	apply ( matrice_points, 1, fonction)
}

# Applique chaque fonction de "liste_fonctions" à chacun des points (ligne de "matrice_points").
# Sortie : matrice dont la i-ème colonne donne les valeurs de la i-ième fonction de "liste_fonctions" aux différents points (lignes de "matrice_points")
# Dépendances : appliqueFonction
creeMatriceValeursFonctions <- function ( matrice_points , liste_fonctions )
{

	vapply ( liste_fonctions, appliqueFonction, rep(1, nrow( matrice_points ) ) , matrice_points = matrice_points )
#	équivalent:
#	sapply ( liste_fonctions, appliqueFonction , matrice_points = matrice_points )

#	vapply est plus économe en ressources que sapply, 
#	mais il faut lui préciser la structure du retour : 
#	le nb de lignes est celui de matrice_points (ie le nombre de points)
}


######################################################################################################################

# Krigeage :

# outer avec un seul vecteur (concaténation de deux ) en argument ,  + un entier indiquant la césure.
outerColonnesConcatenees <- function ( vecteurs_concatenes, taille_vecteur1, fonction )
{
	outer ( vecteurs_concatenes[1:taille_vecteur1] , vecteurs_concatenes [ (taille_vecteur1 +1 ): length( vecteurs_concatenes ) ] , fonction)
}


# Donne les écarts vectoriels entre points (lignes) de "matrice_coordonnees1" et de "matrice_coordonnees2".
# Sortie : matrice calculée comme suit : la différence entre la i-ème ligne de "matrice_coordonnees1" et la j-ième ligne de "matrice_coordonnees2" est la (i-1)*nrow(matrice_coordonnees2)+j -ième ligne de la matrice de sortie.
# REQUIS : "matrice_coordonnees1" et "matrice_coordonnes2" ont le même nb de colonnes (égal au nb de dimensions)
# Dépendances : outerColonnesConcatenees
creeEcart<- function (matrice_coordonnees1, matrice_coordonnees2)
{
	Nb_lignes <-  nrow (matrice_coordonnees1)
	Matrice_concatenee <- rbind ( matrice_coordonnees1 , matrice_coordonnees2 )
	apply ( Matrice_concatenee , 2 , outerColonnesConcatenees, taille_vecteur1 = Nb_lignes , fonction = "-" )
}

# Applique le noyau de covariance "noyau" aux écarts entre points recensés par les lignes de "matrice_ecarts"
# Sortie : vecteur des valeurs du noyau "noyau" pour chaque écart entre points.
appliqueNoyau <- function ( matrice_ecarts, noyau )
{
	apply ( matrice_ecarts, 1, noyau$nom, noyau$parametres )
}


# Sortie : matrice de covariance entre le processus gaussien gouverné par "noyau" pris aux points (lignes) de la matrice "x1" et le même processus gaussien pris aux points (lignes) de la matrice "x2".
# Dépendances : appliqueNoyau
creeMatriceCovariance <- function ( x1, x2, noyau )
{
	#Matrice des écarts entre points connus (nb colonnes = nb dimensions de l'espace, nb lignes = nb couples de points)
	Ecarts <- creeEcart(matrice_coordonnees1=x1, matrice_coordonnees2=x2)

	# Covariance entre points x1 et x2 (présentée sous forme de vecteur)
	Covariance_vecteur<- appliqueNoyau( matrice_ecarts=Ecarts , noyau=noyau )

	# Matrice des covariances entre points x1 et x2 (nb lignes : nb points x1, nb colonnes : nb points x2)
	Covariance_matrice <-  matrix ( Covariance_vecteur , nrow = nrow(x1) )

	Covariance_matrice
}


# Ajoute la pépite (préalablement mise à l'échelle de la variance du processus) à chaque élément diagonal de la matrice de variance-covariance.
# Sortie : "matrice" + "pepite_redimensionnee" x Identité
conditionneMatrice <- function ( matrice, pepite_redimensionnee )
{
	matrice + pepite_redimensionnee * diag ( nrow ( matrice ) )
}


# Sortie : Calcule la matrice de covariance du processus gaussien gouverné par "noyau" aux points (lignes) de la matrice "x1" et le même processus gaussien aux points de la matrice "x2", puis ajoute à la diagonale de cette matrice la pépite "pepite" mise à l'échelle.
# ATTENTION : "y" est censé être le vecteur des valeurs prises par le processus gaussien aux points (lignes) de la matrice "x1". Mais en réalité, il n'est utile que pour le calcul de var(y). Par ailleurs, il faut que "pepite"==0 si x2!=x1, car sinon, le résultat n'a aucun sens. Cette programmation est mauvaise.
# Dépendances : creeMatriceCovariance -> appliqueNoyau ; conditionneMatrice
creeMatriceCovarianceConditionnee <- function ( x1, x2, noyau, y, pepite )
{
	conditionneMatrice ( matrice= creeMatriceCovariance ( x1=x1, x2=x2, noyau=noyau ) , pepite_redimensionnee=var(y) * pepite )
}

#########################################################################################################################

# x_connus est une MATRICE (nb lignes = nb points, nb colonnes = nb dimensions de l'espace)
# y_connus est un VECTEUR de longueur nb points (OBSOLÈTE : Attention, il faut absolument que dim(y_connus) == NULL ! )
# x_nouveaux est une MATRICE (nb lignes = nb points, nb colonnes = nb dimensions de l'espace)
# De toute façon, x_connus, y_connus et x_nouveaux sont convertis à leurs bons types respectifs.
# noyau est une liste contenant toutes les informations
# pepite résoud les problèmes éventuels de conditionnement
# liste_fonctions est un liste contenant les fonctions de régression linéaire (krigeage universel). Si ==NULL, on est en krigeage simple.
# Sortie : prédiction de krigeage aux points (lignes) de la matrice "x_nouveaux", calculée à partir des valeurs "y_connus" du processus gouverné par "noyau" aux points (lignes) de "x_connus".
# Dépendances : creeMatriceCovarianceConditionnee -> conditionneMatrice ; creeMatriceCovariance -> appliqueNoyau
# Dépendances : creeMatriceValeursFonctions ->  appliqueFonction
krige <- function(x_connus, y_connus, x_nouveaux, noyau, pepite = 0, liste_fonctions = NULL )
{
	# Préliminaires
	# x_connus et x_nouveaux doivent être des matrices, y_connus un vecteur
	if ( is.null(dim(x_connus)) ) x_connus <- matrix(x_connus)
	if ( is.null(dim(x_nouveaux)) ) x_nouveaux <- matrix(x_nouveaux)
	dim(y_connus) <- NULL


	###################################################################################
	# Matrice des covariances entre points connus
	Connus_covariance_matrice <- creeMatriceCovarianceConditionnee ( x1=x_connus , x2=x_connus,  noyau=noyau, y=y_connus, pepite=pepite )
	# Reconditionnement éventuel
	#Connus_covariance_matrice <- Connus_covariance_matrice + var(y_connus) * pepite * diag(nrow(Connus_covariance_matrice))

	# Matrice des covariances entre points connus et points nouveaux (lignes : points connus, colonnes : points nouveaux)
	Nouveaux_covariance_matrice <-   creeMatriceCovariance ( x1=x_connus , x2=x_nouveaux , noyau=noyau )



	####################################################################################

	if ( is.null(liste_fonctions ) )	# Krigeage simple
	{
		# remplissage par 'placeholders'
		Fonctions_points_connus <- matrix(0)
		Fonctions_points_nouveaux <- matrix(0)
		Inverse_cov_estimation_coeffs <- matrix(1)
		Mat_auxiliaire <- matrix(0)
		Estimation_coeffs <- matrix(0)
	}

	else 					# Krigeage universel
	{
		# Evaluation des fonctions aux points connus et aux points nouveaux
		Fonctions_points_connus <- creeMatriceValeursFonctions ( matrice_points = x_connus , liste_fonctions = liste_fonctions )
		Fonctions_points_nouveaux <- creeMatriceValeursFonctions ( matrice_points = x_nouveaux , liste_fonctions = liste_fonctions )	

		##FAUX:	Mat_auxiliaire <- Fonctions_points_nouveaux - t( Fonctions_points_connus ) %*% solve ( ... )
		##car Fonctions_points_nouveaux est moralement t(h(x_new)) de F. Bachoc
		Mat_auxiliaire <- 	t (Fonctions_points_nouveaux) -
					t( Fonctions_points_connus ) %*% solve ( Connus_covariance_matrice , Nouveaux_covariance_matrice )

		# Inverse de la matrice de covariance de l'estimateur du vecteur des coordonnees de la fonction moyenne du processus gaussien
		Inverse_cov_estimation_coeffs <- t( Fonctions_points_connus ) %*% solve( Connus_covariance_matrice , Fonctions_points_connus )

		# Estimation des coefficients pondérant les fonctions
		Estimation_coeffs <- 	solve( Inverse_cov_estimation_coeffs  , 
					t( Fonctions_points_connus) %*% solve( Connus_covariance_matrice , y_connus ) )
		
	}

	##################################################################################
	# ESPÉRANCE conditionnelle au nouveau point sachant les points connus
	##FAUX:	Nouveaux_points_moyenne_cond <- colSums (Fonctions_points_nouveaux * Estimation_coeffs ) + colSums( ... )
	##car Fonctions_points_nouveaux est moralement t(h(x_new)) de F. Bachoc

	# Pour que colSums( A * B ) fonctionne, il est nécessaire que B soit de type vecteur
	Nouveaux_points_moyenne_cond <-	Fonctions_points_nouveaux %*% Estimation_coeffs + 
					colSums( Nouveaux_covariance_matrice * 
					c( solve(Connus_covariance_matrice,y_connus - Fonctions_points_connus %*% Estimation_coeffs ) ) )	

	# Du fait du produit matriciel, Nouveaux_points_moyenne_cond est une matrice. Or il doit être un vecteur.
	Nouveaux_points_moyenne_cond <- c( Nouveaux_points_moyenne_cond ) 


	###################################################################################
	# VARIANCE conditionnelle au nouveau point sachant les points connus
	Nouveaux_points_variance_cond<-	noyau$parametres$variance -
					colSums( Nouveaux_covariance_matrice * solve(Connus_covariance_matrice,Nouveaux_covariance_matrice) ) +
					colSums ( Mat_auxiliaire * solve ( Inverse_cov_estimation_coeffs , Mat_auxiliaire ) )
	

	###################################################################################
	# Présentation des résultats

	resultat <- rbind(Nouveaux_points_moyenne_cond, Nouveaux_points_variance_cond)

	rownames(resultat) <- c("Moyenne conditionnelle", "Variance conditionnelle")

	x_nouveaux_liste <- lapply(seq_len(nrow(x_nouveaux)), function(i) x_nouveaux[i,])
	colnames(resultat) <- x_nouveaux_liste

	resultat
}





