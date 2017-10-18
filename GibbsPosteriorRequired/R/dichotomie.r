## ATTENTION : Matern_these.r est suppose charge

dichotomieIntervallePari <- function(fonctionRepartition, ..., probaInf, probaSup, pas = 1, tol = 0.001, borneInf_moins= -1, borneInf_plus= 0, borneSup_moins= 0, borneSup_plus= 1)
    {
        erreurSup = Inf
        erreurInf = Inf

        proba_moyenne <- -1
        
        while(erreurSup > tol)
        {
            if(fonctionRepartition(borneSup_plus,...)<probaSup)
            {
                borneSup_moins <- borneSup_plus
                borneSup_plus <- borneSup_plus + pas
            } else if(fonctionRepartition(borneSup_moins,...)>probaSup)
            {
                borneSup_plus <- borneSup_moins
                borneSup_moins <- borneSup_moins - pas
            } else
            {
                moyenne <- (borneSup_moins + borneSup_plus)/2
                proba_moyenne <- fonctionRepartition(moyenne,...)
                if(proba_moyenne < probaSup) borneSup_moins <- moyenne
                else borneSup_plus <- moyenne
            }
            erreurSup = abs(proba_moyenne-probaSup)
            #print(paste("erreurSup=",erreurSup))
        }

        borneSup <- moyenne

        proba_moyenne <- -1
        
        while(erreurInf > tol)
        {
            if(fonctionRepartition(borneInf_moins,...)>probaInf)
            {
                borneInf_plus <- borneInf_moins
                borneInf_moins <- borneInf_moins - pas
            } else if(fonctionRepartition(borneInf_plus,...)<probaInf)
            {
                borneInf_moins <- borneInf_plus
                borneInf_plus <- borneInf_plus + pas
            } else
            {
                moyenne <- (borneInf_moins + borneInf_plus)/2
                proba_moyenne <- fonctionRepartition(moyenne,...)
                if(proba_moyenne > probaInf) borneInf_plus <- moyenne
                else borneInf_moins <- moyenne
            }
            erreurInf <- abs(proba_moyenne-probaInf)
            #print(paste("erreurInf=",erreurInf))
        }    
        borneInf <- moyenne
        c(borneInf,borneSup)       
}


vecteurMatriceCorrelationAnisGeom <- function(longueurs_correlation, regularite, x_connus)
{
    Noyau = creeMaternIsotrope(variance = 1, longueur = longueurs_correlation, regularite = regularite )

    res <- c(creeMatriceCovariance(x1=x_connus, x2=x_connus, noyau= Noyau))
res
}

vecteurMatriceCorrelationTens <- function(longueurs_correlation, regularite, x_connus)
{
    Noyau = creeMaternTensorise(variance = 1, longueur = longueurs_correlation, regularite = regularite )

    res <- c(creeMatriceCovariance(x1=x_connus, x2=x_connus, noyau= Noyau))
res
}

vecteurInverseMatriceRestreinte <- function(vecteur_matrice,nb_lignes, injectionOrthogonalTendance)
{
	if(sum(vecteur_matrice-1)==0) res <- solve(t(injectionOrthogonalTendance) %*% injectionOrthogonalTendance) else # on verifie que la matrice de correlation n'est pas remplie de 1, auquel cas on decrete l'independance
	res <- solve(t(injectionOrthogonalTendance) %*% matrix(vecteur_matrice,nrow=nb_lignes) %*% injectionOrthogonalTendance)
res
}


vecteurInverseMatriceCorrelationAnisGeom <- function(longueurs_correlation, regularite, x_connus)
{
    Noyau = creeMaternIsotrope(variance = 1, longueur = longueurs_correlation, regularite = regularite )

    res <- c(solve(creeMatriceCovariance(x1=x_connus, x2=x_connus, noyau= Noyau)))
res
}

## Attention, on suppose que le point nouveau est unique
## SORTIE : vecteur de correlation entre les points connus (x_connus) et le point_nouveau
vecteurCorrelConnusNouveauAnisGeom <- function(longueurs_correlation, regularite, x_connus, point_nouveau)
{
    Noyau = creeMaternIsotrope(variance = 1, longueur = longueurs_correlation, regularite = regularite )
    
    res <- c(creeMatriceCovariance(x1=x_connus, x2=point_nouveau, noyau= Noyau))
res
}


vecteurInverseMatriceCorrelationTens <- function(longueurs_correlation, regularite, x_connus, injectionOrthogonalTendance)
{
    Noyau = creeMaternTensorise(variance = 1, longueur = longueurs_correlation, regularite = regularite )

    res <- c(solve(creeMatriceCovariance(x1=x_connus, x2=x_connus, noyau= Noyau)))
res
}

## Attention, on suppose que le point nouveau est unique
## SORTIE : vecteur de correlation entre les points connus (x_connus) et le point_nouveau
vecteurCorrelConnusNouveauTens <- function(longueurs_correlation, regularite, x_connus, point_nouveau)
{
    Noyau = creeMaternTensorise(variance = 1, longueur = longueurs_correlation, regularite = regularite )
    
    res <- c(creeMatriceCovariance(x1=x_connus, x2=point_nouveau, noyau= Noyau))
res
}

##pnormParametresDevant <- function(parametres,t) pnorm(t,mean=parametres[1],sd=sqrt(parametres[2])) ##ancienne version
pnormParametresDevant <- function(parametres,t) pnorm(t,mean=parametres[,1],sd=sqrt(parametres[,2]))

##dnormParametresDevant <- function(parametres,t) dnorm(t,mean=parametres[1],sd=sqrt(parametres[2])) ##ancienne version
dnormParametresDevant <- function(parametres,t) dnorm(t,mean=parametres[,1],sd=sqrt(parametres[,2]))
