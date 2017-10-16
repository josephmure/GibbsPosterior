#source("Matern_these.r")
#source("Krige_these.r")



evalueVraisemblanceIntegreeAnisGeom <- function( longueur_correlation, x_connus, y_connus, regularite, tendance = matrix(0))
    {
    Noyau <- creeMaternIsotrope ( variance= 1, longueur= longueur_correlation, regularite= regularite)

    matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)

if(det(matriceCorrelation)!=0)
{
	tendanceInverseMatriceCorrelationTendance <- matrix(1)
	nbFonctionsTendance <- 0
	if(length(tendance)>1) ## en cas de krigeage simple : tendance est une matrice 1x1 contenant 0
	{
		nbFonctionsTendance <- ncol(tendance)
		inverseMatriceCorrelationTendance <- solve(matriceCorrelation,tendance)
		tendanceInverseMatriceCorrelationTendance <- t(tendance) %*% inverseMatriceCorrelationTendance
		ProjectionTendance <- tendance %*% solve(tendanceInverseMatriceCorrelationTendance, t(inverseMatriceCorrelationTendance))
		y_connus <- y_connus - ProjectionTendance %*% y_connus
	}
 	Vraisemblance_integree <- det(matriceCorrelation)^(-1/2)* det(tendanceInverseMatriceCorrelationTendance)^(-1/2) * (sum(y_connus* solve(matriceCorrelation,y_connus)))^(-(nrow(x_connus)-nbFonctionsTendance)/2)
}
else Vraisemblance_integree <- 0
    Vraisemblance_integree

}

evalueVraisemblanceIntegreeTens <- function( longueur_correlation, x_connus, y_connus, regularite, tendance = matrix(0))
    {
    Noyau <- creeMaternTensorise ( variance= 1, longueur= longueur_correlation, regularite= regularite)

    matriceCorrelation <- creeMatriceCovariance(x1= x_connus,x2=  x_connus, noyau= Noyau)

if(det(matriceCorrelation)!=0)
{
	tendanceInverseMatriceCorrelationTendance <- matrix(1)
	nbFonctionsTendance <- 0
	if(length(tendance)>1) ## en cas de krigeage simple : tendance est une matrice 1x1 contenant 0
	{
		nbFonctionsTendance <- ncol(tendance)
		inverseMatriceCorrelationTendance <- solve(matriceCorrelation,tendance)
		tendanceInverseMatriceCorrelationTendance <- t(tendance) %*% inverseMatriceCorrelationTendance
		ProjectionTendance <- tendance %*% solve(tendanceInverseMatriceCorrelationTendance, t(inverseMatriceCorrelationTendance))
		y_connus <- y_connus - ProjectionTendance %*% y_connus
	}
 	Vraisemblance_integree <- det(matriceCorrelation)^(-1/2)* det(tendanceInverseMatriceCorrelationTendance)^(-1/2) * (sum(y_connus* solve(matriceCorrelation,y_connus)))^(-(nrow(x_connus)-nbFonctionsTendance)/2)
}
else Vraisemblance_integree <- 0
    Vraisemblance_integree

}


## Attention, les fonctions d'optimisation preferent minimiser des logarithmes, plutot que de maximiser directement les fonctions

evalueOpposeLogVraisemblanceIntegreeAnisGeom <- function(longueur_correlation, x_connus, y_connus, regularite, tendance = "vide" )
{
    res <-  - log(evalueVraisemblanceIntegreeAnisGeom(longueur_correlation = longueur_correlation, x_connus=x_connus, y_connus= y_connus, regularite = regularite, tendance = tendance))
    if(res==Inf) res <- 1000
    res
}

evalueOpposeLogVraisemblanceIntegreeTens <- function(longueur_correlation, x_connus, y_connus, regularite, tendance = "vide" )
{
    res <-  - log(evalueVraisemblanceIntegreeTens(longueur_correlation = longueur_correlation, x_connus=x_connus, y_connus= y_connus, regularite = regularite, tendance = tendance))
    if(res==Inf) res <- 1000
    res
}
