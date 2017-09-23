# ATTENTION !!! Le 19/11/2015, j'ai découvert une erreur dans la programmation de stage : dans la fonction "Matern", le code affirmait : " if( t==0 ) valeur = 1", alors qu'il aurait fallu écrire " if ( t==0 ) valeur = variance ".
# Cette erreur est corrigée dans cette nouvelle version de la programmation de stage.

# Un noyau est vu comme une liste contenant un nom et une liste de paramètres
# Même si MaternUnDemi, MaternTroisDemis et MaternCinqDemis fournissent des expressions analytiques intéressantes, il est probable que le gain de temps obtenu ainsi soit négligeable et ne justifie pas l'accroissement de la complexité de la programmation qu'elles induisent. Dorénavant, elles seront considérées comme obsolètes.



# Noyau trivial

#noyauTrivial <- function ( t, variance= 1 ) 
#{
#	variance
#}

#creeNoyauTrivial <- function( variance= 1 )
#{
#	list ( nom= noyauTrivial , parametres = list ( variance = variance ) )
#}

# Noyaux de Matérn

# Dimension 1 :

MaternUnDemi<- function(t, variance, longueur)
{
	variance *  exp(- sqrt(2) *abs(t) / longueur)
}

MaternTroisDemis<- function(t, variance, longueur)
{
	variance * ( 1 + sqrt(6) * abs(t) / longueur ) * exp(- sqrt(6) * abs(t) / longueur)
}

MaternCinqDemis<- function(t, variance, longueur)
{
	variance * ( 1 + sqrt(10) * abs(t) / longueur + 10/3 * t^2 / longueur^2) * exp(- sqrt(10) * abs(t) / longueur)
}

MaternInf<- function(t, variance, longueur)
{
	variance *  exp(- t^2 / longueur^2 )
}


Matern <- function ( t, variance, longueur, regularite)
{
#	if ( t == 0) valeur = variance
#	else
#	{
#		aux <- 2 * sqrt (regularite) * abs(t) / longueur
#		valeur <- variance / ( gamma(regularite) * 2^(regularite - 1) ) * aux^regularite * besselK(aux, regularite)
#	}
#	valeur

	t[t==0] <- NA
	aux <- 2 * sqrt (regularite) * abs(t) / longueur
	valeur <- variance / ( gamma(regularite) * 2^(regularite - 1) ) * aux^regularite * besselK(aux, regularite)
	valeur[is.na(valeur)] <- variance
	valeur
}

MaternDeriveLongueur <- function ( t, variance, longueur, regularite)
{
	if ( t == 0) valeur = 0
	else
	{
		aux <- 2 * sqrt (regularite) * abs(t) / longueur
		valeur <- variance / ( gamma(regularite) * 2^(regularite - 1) ) / longueur * aux^(regularite+1) * besselK(aux, regularite-1)
	}
	valeur
}



#####################################################################################################################

# Multidimensionnel isotrope :

MaternUnDemiIsotrope <- function ( t, parametres = list (variance = 1, longueur = 1) )
{
	MaternUnDemi ( sqrt ( sum ( t^2 / parametres$longueur^2 ) ), variance = parametres$variance, longueur = 1) 	
}


MaternTroisDemisIsotrope <- function ( t, parametres = list (variance = 1, longueur = 1))
{
	MaternTroisDemis ( sqrt ( sum ( t^2 / parametres$longueur^2 ) ), variance = parametres$variance, longueur = 1) 	
}


MaternCinqDemisIsotrope <- function ( t, parametres = list (variance = 1, longueur = 1))
{
	MaternCinqDemis ( sqrt ( sum ( t^2 / parametres$longueur^2 ) ), variance = parametres$variance, longueur = 1) 	
}

MaternInfIsotrope <- function ( t, parametres = list (variance = 1, longueur = 1))
{
	MaternInf ( sqrt ( sum ( t^2 / parametres$longueur^2 ) ), variance = parametres$variance, longueur = 1) 	
}

MaternIsotrope <- function ( t, parametres)
{
	Matern ( sqrt ( sum ( t^2 / parametres$longueur^2 ) ), variance = parametres$variance, longueur = 1, regularite = parametres$regularite ) 
}


# Astuce (à voir si c'est bien malin) : on remplace l'ancienne fonction creeMaternIsotrope par creeMaternIsotropePresoutenance
creeMaternIsotropePresoutenance <- function ( variance, longueur, regularite )
{
	if ( regularite == 0.5 )
	{
		list ( nom=MaternUnDemiIsotrope, parametres = list ( variance = variance , longueur = longueur ) )
	}

	else
	{

		if ( regularite == 1.5)
		{
			list ( nom=MaternTroisDemisIsotrope, parametres = list ( variance = variance , longueur = longueur ) )
		}

		else
		{
			if ( regularite == 2.5)
			{
				list ( nom=MaternCinqDemisIsotrope, parametres = list ( variance = variance , longueur = longueur ) )
			}
		

			else
			{
				if ( regularite == Inf)
				{
					list ( nom=MaternInfIsotrope, parametres = list ( variance = variance , longueur = longueur ) )
				}


				else
				{
					list ( nom= MaternIsotrope,parametres= list( variance= variance, longueur= longueur, regularite= regularite ) )
				}
			}
		 }
	 }
}




creeMaternIsotrope <- function ( variance, longueur, regularite )
{
	if ( regularite == Inf)
	{
		list ( nom=MaternInfIsotrope, parametres = list ( variance = variance , longueur = longueur ) )
	}


	else
        {
		list ( nom= MaternIsotrope,parametres= list( variance= variance, longueur= longueur, regularite= regularite ) )
	}

}









#############################################################################################################################
# Multidimensionnel tensorisé :

MaternUnDemiTensorise <- function ( t, parametres = list ( variance = 1, longueur = 1) )
{
	parametres$variance * prod ( MaternUnDemi ( t / parametres$longueur , variance = 1 , longueur = 1 ) )	
}

MaternTroisDemisTensorise <- function ( t, parametres = list ( variance = 1, longueur = 1) )
{
	parametres$variance * prod ( MaternTroisDemis ( t / parametres$longueur , variance = 1 , longueur = 1 ) )	
}

MaternCinqDemisTensorise <- function ( t, parametres = list ( variance = 1, longueur = 1) )
{
	parametres$variance * prod ( MaternCinqDemis ( t / parametres$longueur , variance = 1 , longueur = 1 ) )	
}

MaternInfTensorise <- function ( t, parametres = list ( variance = 1, longueur = 1) )
{
	parametres$variance * prod ( MaternInf ( t / parametres$longueur , variance = 1 , longueur = 1 ) )	
}

MaternTensorise <- function ( t, parametres )
{
	parametres$variance * prod ( Matern ( t , variance = 1 , longueur = parametres$longueur , regularite = parametres$regularite ) )
}



# Astuce (à voir si c'est bien malin) : on remplace l'ancienne fonction creeMaternTensorise par creeMaternTensorisePresoutenance
creeMaternTensorisePresoutenance <- function ( variance, longueur, regularite )
{
	if ( regularite == 0.5 )
	{
		list ( nom=MaternUnDemiTensorise, parametres = list ( variance = variance , longueur = longueur ) )
	}

	else
	{

		if ( regularite == 1.5)
		{
			list ( nom=MaternTroisDemisTensorise, parametres = list ( variance = variance , longueur = longueur ) )
		}

		else
		{
			if ( regularite == 2.5)
			{
				list ( nom=MaternCinqDemisTensorise, parametres = list ( variance = variance , longueur = longueur ) )
			}
		

			else
			{
				if ( regularite == Inf)
				{
					list ( nom=MaternInfTensorise, parametres = list ( variance = variance , longueur = longueur ) )
				}


				else
				{
					list ( nom= MaternTensorise,parametres= list( variance= variance, longueur= longueur, regularite= regularite ) )
				}
			}
		 }
	 }
}


creeMaternTensorise <- function ( variance, longueur, regularite )
{

	if ( regularite == Inf)
	{
		list ( nom=MaternInfTensorise, parametres = list ( variance = variance , longueur = longueur ) )
	}


	else
	{
		list ( nom= MaternTensorise, parametres= list(variance= variance, longueur= longueur, regularite= regularite) )
	}

}

## OBSOLÈTE ##
# Pour les applications automatisées, il faut que les fonctions de création de noyaux acceptent les hyperparamètres sous forme vectorielle.

creeMaternIsotropeVectoriel <- function ( vecteur_parametres )
{
	creeMaternIsotrope ( variance= vecteur_parametres[1], regularite= vecteur_parametres[2], longueur= vecteur_parametres[-c(1,2)] )
}

creeMaternTensoriseVectoriel <- function ( vecteur_parametres )
{
	creeMaternTensorise ( variance= vecteur_parametres[1], regularite= vecteur_parametres[2], longueur= vecteur_parametres[-c(1,2)] )
}
## FIN OBSOLÈTE ##





#########################################################################################################################################
# Certains algorithmes d'optimisation utilisent les dérivées partielles des noyaux de Matérn
# Attention, les fonctions ci-dessous ne permettent pas de gérer le cas où la régularité vaut + l'infini.

# t : vecteur contenant les coordonnées du points où on dérive
# Sortie : vecteur contenant les dérivées partielles du noyau de Matérn multidimensionnel isotrope par rapport aux différentes longueurs de corrélation
MaternIsotropeDeriveLongueur <- function ( t, parametres= list( variance=1, longueur, regularite , indice_derive) )
    {
        rayonModule <- sqrt( sum( t^2 / parametres$longueur^2 ) )

        aux <- 2 * sqrt (parametres$regularite) * rayonModule
        deriveeMatern1D <- -parametres$variance*2*sqrt(parametres$regularite) / gamma( parametres$regularite ) / 2^( parametres$regularite - 1 ) * aux^parametres$regularite * besselK( aux, parametres$regularite - 1 )

        if(rayonModule==0) res <- 0 else
	{
		if(length(parametres$longueur)==1) res <- - rayonModule / parametres$longueur * deriveeMatern1D else
        	res <- - t[parametres$indice_derive]^2 / parametres$longueur[parametres$indice_derive]^3 / rayonModule * deriveeMatern1D
	}
        
##        print(res)
        res
    }



creeMaternIsotropeDeriveLongueur <- function( variance, longueur, regularite, indice_derive)
    {
        list ( nom= MaternIsotropeDeriveLongueur, parametres= list (variance= variance, longueur= longueur, regularite= regularite, indice_derive= indice_derive ) )
    }


MaternTensoriseDeriveLongueur <- function( t, parametres )
    {
	parametres$variance * prod ( Matern ( t[ - parametres$indice_derive ] , variance = 1 , longueur = parametres$longueur[- parametres$indice_derive] , regularite = parametres$regularite ) ) * MaternDeriveLongueur( t[parametres$indice_derive], variance=1, longueur= parametres$longueur[parametres$indice_derive], regularite= parametres$regularite)
    }

creeMaternTensoriseDeriveLongueur <- function( variance, longueur, regularite, indice_derive)
    {
        list ( nom= MaternTensoriseDeriveLongueur, parametres= list (variance= variance, longueur= longueur, regularite= regularite, indice_derive= indice_derive ) )
    }



## Sortie : vecteur contenant les dérivées partielles (ou plutôt, les différences finies pour un certain pas ) du noyau de Matérn multidimensionnel isotrope par rapport à la régularité
MaternIsotropeDeriveRegularite <- function( t,  parametres= list( variance=1, longueur, regularite, pas)  )
    {
        parametresPlusPas <- parametres
        parametresPlusPas$regularite <- parametres$regularite + parametres$pas
        MaternRegPlusPas <- MaternIsotrope( t, parametres= parametresPlusPas )
        MaternReg <-  MaternIsotrope( t, parametres=parametres )

        ( MaternRegPlusPas - MaternReg ) / parametres$pas
    }

creeMaternIsotropeDeriveRegularite <- function( variance, longueur, regularite, pas )
    {
        list ( nom= MaternIsotropeDeriveRegularite, parametres= list (variance= variance, longueur= longueur, regularite= regularite, pas=pas ))
    }

