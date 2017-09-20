## Transformation de Josselin


## Sortie : transformation de Josselin appliquée au vecteur de données
## Dépendances : aucune
## Attention, le vecteur_donnees doit nécessairement avoir des entrées strictement positives
transfoJosselin <- function(vecteur_donnees, parametre)
{
	if( parametre > 0) res <- 1/parametre * asinh( parametre * log(vecteur_donnees) )

	else if( parametre == 0 ) res <- log(vecteur_donnees)

	else res <- 1/parametre * sinh( parametre * log(vecteur_donnees) )

	res
}



deriveTransfoJosselin <- function(vecteur_donnees, parametre)
{
	if( parametre > 0) res <- 1/vecteur_donnees/ sqrt( ( parametre * log(vecteur_donnees) )^2 + 1)

	else if( parametre == 0 ) res <- 1/vecteur_donnees

	else res <- 1/vecteur_donnees * cosh( parametre * log(vecteur_donnees) )

	res
}
