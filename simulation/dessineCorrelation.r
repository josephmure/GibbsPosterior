source("Matern_these.r")

Nb_abscisses <- 100
Points <- matrix(NA,Nb_abscisses^2,2)
for(i in 1:Nb_abscisses)
{
	for(j in 1:Nb_abscisses)
	{
		Points[(i-1)*Nb_abscisses + j,] <- c(i,j)
	}
}
Points <- (Points / (Nb_abscisses+1) - 0.5)/2

dessineCorrelationAnisGeom <- function(Points, longueurs_correlation, variance=1, regularite=2.5, couleur="blue")
{
parametres <- list(variance= variance, regularite= regularite, longueurs= longueurs_correlation)

Correlation <- apply(Points,1,MaternIsotrope,parametres=parametres)

Data = data.frame(cbind(Points, Correlation))

gg <- ggplot() + stat_contour(data=Data,aes(x=V1,y=V2, z=Correlation, fill=..level..), geom="polygon") + scale_fill_continuous(name="Correlation", low="white", high=couleur) + theme_bw() + coord_fixed(xlim=c(min(Points[,1]),max(Points[,1])), ylim=c(min(Points[,2]),max(Points[,2])))

gg

}

dessineCorrelationTens <- function(Points, longueurs_correlation, variance=1, regularite=2.5, couleur="blue")
{
parametres <- list(variance= variance, regularite= regularite, longueurs= longueurs_correlation)

Correlation <- apply(Points,1,MaternTensorise,parametres=parametres)

Data = data.frame(cbind(Points, Correlation))

gg <- ggplot() + stat_contour(data=Data,aes(x=V1,y=V2, z=Correlation, fill=..level..), geom="polygon") + scale_fill_continuous(name="Correlation", low="white", high=couleur) + theme_bw() + coord_fixed(xlim=c(min(Points[,1]),max(Points[,1])), ylim=c(min(Points[,2]),max(Points[,2])))

gg
}

dessineContourTens <- function(Points, longueurs_correlation, variance=1, regularite=2.5)
{
parametres <- list(variance= variance, regularite= regularite, longueurs= longueurs_correlation)

Correlation <- apply(Points,1,MaternTensorise,parametres=parametres)

Data = data.frame(cbind(Points, Correlation))

gg <- ggplot() + geom_contour(data=Data,aes(x=V1,y=V2, z=Correlation)) + theme_bw() + coord_fixed(xlim=c(min(Points[,1]),max(Points[,1])), ylim=c(min(Points[,2]),max(Points[,2])))

gg
}

dessineContourAnisGeom <- function(Points, longueurs_correlation, variance=1, regularite=2.5)
{
parametres <- list(variance= variance, regularite= regularite, longueurs= longueurs_correlation)

Correlation <- apply(Points,1,MaternIsotrope,parametres=parametres)

Data = data.frame(cbind(Points, Correlation))

gg <- ggplot() + geom_contour(data=Data,aes(x=V1,y=V2, z=Correlation)) + theme_bw() + coord_fixed(xlim=c(min(Points[,1]),max(Points[,1])), ylim=c(min(Points[,2]),max(Points[,2])))

gg
}

## dataCoord jeu de points quelconque
## bb <-  bb + geom_point(data=dataCoord, aes(x=abscisses, y=ordonnees))
