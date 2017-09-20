calculeDistMinimale <- function(vecteur,matrice)
{
	vecteur <- matrix(rep(vecteur,nrow(matrice)),nrow=nrow(matrice),byrow=TRUE)
	matrice <- (matrice - vecteur)^2
	min(sqrt(apply(matrice,1,sum)))
}

calculeAutoDist <- function(matrice)
{
	res <- rep(NA,nrow(matrice))
	for (numeroLigne in 1:nrow(matrice))
	{
		res[numeroLigne] <- calculeDistMinimale(matrice[numeroLigne,],matrice[-numeroLigne,])
	}
	mean(res)
}
