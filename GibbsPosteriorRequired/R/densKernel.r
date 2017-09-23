densKernel <- function(argument,fenetre,matrice_echantillon)
{
    if(length(argument)!=ncol(matrice_echantillon))
    {
        print("L'argument n'a pas la bonne taille.")
        print(paste("Taille de l'argument :",length(argument)))
        print(paste("Dimensions de la matrice_echantillon :", nrow(matrice_echantillon), ncol(matrice_echantillon)))
        res <- 0
    }
    
    else
    {
        matrice_echantillon <- matrix(rep(argument,nrow(matrice_echantillon)),ncol=length(argument),byrow=TRUE)- matrice_echantillon
        matrice_echantillon <- matrice_echantillon * matrice_echantillon / fenetre / fenetre
        
	vecteur_normes2 <- apply(matrice_echantillon,1,sum)
	vecteur_normes2 <- exp(-vecteur_normes2)

        res <- - sum(vecteur_normes2)
    }

    res
}

