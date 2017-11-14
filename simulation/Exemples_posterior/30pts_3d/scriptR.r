#library(ggplot2)
#library(pscl)

NOMBRE_DIMENSIONS <- scan("../nombre_dimensions.txt")
REGULARITE <- scan("../regularite_vraie.txt")
TYPE_NOYAU_MATERN_TYPE_PRIOR <- scan("../type_noyau_prior.txt", what="character")
TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
#TYPE_PRIOR <- TYPE_NOYAU_MATERN_TYPE_PRIOR[2] ##Pas besoin ici.
NOMBRE_POINTS_POSTERIOR_A_GENERER <- 1000
NOMBRE_PAS_METROPOLIS <- 20
ECART_TYPE_METROPOLIS <- 0.2

GERME_ALEATOIRE_C <- germe_aleatoire


#Creation de fichiers de configuration locaux pour le programme C++
vecteur_metaparametres <- c(NOMBRE_DIMENSIONS, REGULARITE, NOMBRE_POINTS_POSTERIOR_A_GENERER, NOMBRE_PAS_METROPOLIS, ECART_TYPE_METROPOLIS, length(FONCTIONS))
write.table(vecteur_metaparametres,"metaparametres.txt",row.names=FALSE,col.names=FALSE)



#Execution du programme C++
system(paste0("GSL_RNG_TYPE='taus' GSL_RNG_SEED=",GERME_ALEATOIRE_C, " ./GibbsPosterior") )

#Apres execution du programme C++, les fichiers de configuration locaux peuvent etre detruite
system("rm metaparametres.txt")
## system("rm planXP.txt") ## Ici, planXP a ete cree par le programme R, il est donc plus propre de le supprimer avec nettoyage.r

Posterior <- NULL

if(TYPE_NOYAU_MATERN == "geometrique")
{
	Posterior<-matrix(scan("point_posterior_anis_geom.txt"),ncol=NOMBRE_DIMENSIONS,byrow=TRUE)
	##system("rm point_posterior_anis_geom.txt")
	write.table(Posterior, append=TRUE, row.names=FALSE, col.names=FALSE, "Posterior_anis_geom.txt", sep= "\t")
}

if (TYPE_NOYAU_MATERN == "tensorise")
{
	Posterior<-matrix(scan("point_posterior_tens.txt"),ncol=NOMBRE_DIMENSIONS,byrow=TRUE)
	##system("rm point_posterior_tens.txt")
	write.table(Posterior, append=TRUE, row.names=FALSE, col.names=FALSE, "Posterior_tens.txt", sep= "\t")
}

###if(!is.null(Posterior)) write.matrix(Posterior, "Posterior.txt", sep= "\t")



