## Optimisation des paramètres par Cholesky et Anderson-Darling

## Sortie : moyenne entre le résultat d'Anderson-Darling sur l'ensemble de test "nuage" x_test_nuage et l'ensemble de test "amas" (ou plutôt croix, en pratique) x_test_amas, pondérée par poids_nuage (compris entre 0 et 1)
## Dépendances : creeMaternIsotrope --> [Matern.r] ;
##               creeMatrices -->  [Estim_parametres_simple.r] ;
##               testeNormalite --> [Briques_these.r]
testeStructureCorrelation <- function( vect_param_forme, x_connus, x_test_nuage, x_test_amas, y_connus, y_test_nuage,y_test_amas, pepite= 0, liste_fonctions= NULL, poids_nuage, regularite= NULL)
    {
        ##  création du noyau
        if(is.null(regularite))
            {                
                reg <- abs( vect_param_forme[1] )
                long <- abs( vect_param_forme[-1] )
            }
        else
            {
                reg <- regularite
                long <- abs(vect_param_forme)
            }

        noyau <- creeMaternIsotrope( variance= 1, longueur= long, regularite= reg)
        var <- 1/ length(y_connus) * t(y_connus) %*% creeMatrices(x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions )$matriceFormeQuadratique %*% y_connus
        noyau$variance <- var
        print(paste("variance=",var))
        res_ad <- 0
        res_nuage <- testeNormalite(noyau= noyau, x_connus= x_connus, x_test= x_test_nuage, y_connus= y_connus, y_test= y_test_nuage, pepite= pepite, liste_fonctions= liste_fonctions )
        if(poids_nuage<1)
            {
                res_amas <-  testeNormalite(noyau= noyau, x_connus= x_connus, x_test= x_test_amas, y_connus= y_connus, y_test= y_test_amas, pepite= pepite, liste_fonctions= liste_fonctions )
                res_ad <- (1-poids_nuage)*res_amas
            }
        res_ad <- res_ad + poids_nuage*res_nuage ## + (1-poids_nuage)*res_amas
        print(res_ad)
        res_ad
    }

## Sortie : résultat d'Anderson-Darling sur l'ensemble de test soumis x_test
## Dépendances : Andeson-Darling --> [Briques_these.r]
testeNormalite <- function (noyau, variance, x_connus, x_test, y_connus, y_test, pepite= 0, liste_fonctions= NULL)
    {
        Anderson_Darling(  erreursCholesky( x_connus= x_connus, x_test= x_test, y_connus= y_connus, y_test= y_test, noyau= noyau, pepite= pepite, liste_fonctions= liste_fonctions ))
    }

        



## Sortie : liste comportant la matrice des paramètres initialisant l'algorithme d'optimisation, la matrice des paramètres trouvés, le vecteur des variances trouvées par maximum de vraisemblance correspondant aux paramètres trouvés, le vecteur comportant les A^2 d'Anderson-Darling pour les paramètres trouvés et le vecteur des warnings rencontrés (0 = tout va bien)
## NB : les numeros de lignes des matrices correspondent aux numeros d'entrees des vecteurs
## NB2: les initialisations se font par tirage aléatoire de loi uniforme déterminée par min_init et max_init. Pour éviter les minima locaux, on reessaie nb_points_initialisation fois.
## Dependances : tout ce qui precede
trouveStructureCorrelation <- function ( x_connus, x_test_nuage, x_test_amas, y_connus, y_test_nuage, y_test_amas, pepite= 0, liste_fonctions= NULL, nb_points_initialisation, min_init, max_init, poids_nuage, regularite=NULL)
    {
        nb_parametres <- length( min_init )

        if( length( max_init) != nb_parametres ) return ( "Combien de parametres ?" )
        
        Matrice_param <- matrix( nrow= nb_points_initialisation, ncol= nb_parametres )
        Matrice_param_init <- Matrice_param
        vect_anderson_darling <- vector(mode= "numeric", length= nb_points_initialisation )
        vect_warning <- vect_anderson_darling
        vect_variance <- vect_anderson_darling
        
        for (n in 1:nb_points_initialisation)
            {


                OPT <- try(solve(0))
                
                while( class(OPT)=="try-error" )

                    {
                        init <- NULL
                        for (i in 1:nb_parametres)
                            {
                                init[i] <- runif(1, min= min_init[i], max= max_init[i] )
                            }

                        Matrice_param_init[n,] <- init
                
                
                        OPT <- try( optim(par= init, fn= testeStructureCorrelation, x_connus= x_connus, x_test_nuage= x_test_nuage, x_test_amas= x_test_amas, y_connus= y_connus, y_test_nuage= y_test_nuage, y_test_amas= y_test_amas, pepite= pepite, liste_fonctions= liste_fonctions, poids_nuage= poids_nuage, regularite= regularite, method= "BFGS" ) )
                        ## , lower= min_init*0 )

                    }

                Matrice_param[n,] <- OPT$par
                vect_anderson_darling[n] <- OPT$value
                vect_warning[n] <- OPT$convergence
                
                if(is.null(regularite)) noyau <- creeMaternIsotrope( variance= 1, longueur=  abs(Matrice_param[n,-1]), regularite=  Matrice_param[n,1])
                else  noyau <- creeMaternIsotrope( variance= 1, longueur=  abs(Matrice_param[n,]), regularite= regularite )

                print(x_connus)
                print(noyau)
                print(liste_fonctions)

                print(creeMatrices(x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions ))
                var <- 1/ length(y_connus) * t(y_connus) %*% creeMatrices(x_connus= x_connus, noyau= noyau, liste_fonctions= liste_fonctions )$matriceFormeQuadratique %*% y_connus
                vect_variance[n] <- var
            }

        list(Matrice_param_init= Matrice_param_init, Matrice_param= Matrice_param, vect_variance= vect_variance, vect_anderson_darling= vect_anderson_darling, vect_warning= vect_warning )
    }

#### Création du processus gaussien

esperanceAffine <- function( a=0, b=0 )
    {
        function(t) a + b*t[1]
    }


x_connus_nuage<- generePointsConnus( nombre_points = 50, minima = c(0.1,0.1), maxima = c(0.9,0.9) )
x_connus_croix <- genereCroix( matrice_centres <- matrix(x_connus_nuage[1:2,],nrow=2), vecteur_nombre_points_par_dimension= c(6,6), vecteur_espacements= c(0.03, 0.03) )
x_test_nuage <- generePointsConnus( nombre_points = c(25), minima = c(0,0), maxima = c(1,1) )
x_test_croix <- genereCroix(  matrice_centres <- matrix(x_test_nuage[1,],nrow=1), vecteur_nombre_points_par_dimension=6, vecteur_espacements= c(0.03, 0.03) )

## Graphe du plan d'expérience

plot( x_connus_nuage, xlim= c(0,1), ylim= c(0,1), main= "Plan d'expériences", xlab= "", ylab="", pch=3)
points( x_connus_croix, pch= 3)
points( x_test_nuage, pch= 3, col= 4)
points( x_test_croix, pch= 3, col= 4)
legend("topright", legend= c("Ensemble d'apprentissage","Ensemble de test"), col=c(1,4),pch=3)

## dev.copy(pdf,"planXP.pdf")
## dev.off()

## NB : write.matrix nécessite library(MASS)

## write.matrix(x_connus_nuage,"planXP_apprentissage_nuage.dat", sep=";" )
## write.matrix(x_connus_croix,"planXP_apprentissage_croix.dat", sep=";" )
## write.matrix(x_test_nuage,"planXP_test_nuage.dat", sep=";" )
## write.matrix(x_test_croix,"planXP_test_croix.dat", sep=";" )

y_complet <- apply( rbind(x_connus_nuage, x_connus_croix, x_test_nuage, x_test_croix) , 1, esperanceAffine(1,2) ) + creeEchantillonNormal(x_connus=  rbind(x_connus_nuage, x_connus_croix, x_test_nuage, x_test_croix), noyau= creeMaternIsotrope(variance=1, longueur= c(0.3,0.4), regularite= 1.5) )             

y_connus_nuage <- y_complet[ 1:nrow(x_connus_nuage)]
y_connus_croix <- y_complet[ (nrow(x_connus_nuage)+1):(nrow(x_connus_nuage)+nrow(x_connus_croix)) ]
y_test_nuage <- y_complet[ (nrow(x_connus_nuage)+nrow(x_connus_croix)+1) : (nrow(x_connus_nuage)+nrow(x_connus_croix)+nrow(x_test_nuage)) ]
y_test_croix<- y_complet [ (nrow(x_connus_nuage)+nrow(x_connus_croix)+nrow(x_test_nuage)+1) : (nrow(x_connus_nuage)+nrow(x_connus_croix)+nrow(x_test_nuage)+nrow(x_test_croix)) ]

sc_separes<- trouveStructureCorrelation ( x_connus= rbind(x_connus_nuage, x_connus_croix), x_test_nuage, x_test_croix, y_connus= c(y_connus_nuage, y_connus_croix), y_test_nuage, y_test_croix, pepite= 0, liste_fonctions= list(function(t) 1, function(t) t[1] ), nb_points_initialisation=100, min_init=c(0,0), max_init=c(1,1),poids_nuage= 0.5, regularite= 1.5)


sc_confondus<-  trouveStructureCorrelation ( x_connus= rbind(x_connus_nuage, x_connus_croix), x_test_nuage=rbind(x_test_nuage, x_test_croix), x_test_amas=matrix(c(1,1), nrow=1) , y_connus= c(y_connus_nuage, y_connus_croix), y_test_nuage=c(y_test_nuage, y_test_croix), y_test_amas=0, pepite= 0, liste_fonctions= list(function(t) 1, function(t) t[1] ), nb_points_initialisation=100, min_init=c(0,0), max_init=c(1,1), poids_nuage=1, regularite=1.5)


# Graphe des paramètres trouvés

# Méthode avec ensembles de test séparés
plot(abs(sc_separes$Matrice_param_init), main="Paramètres estimés par nuage et croix d'apprentissage séparés", xlim=c(0,1), ylim=c(0,1), pch="+",xlab="Longueur de corrélation 1", ylab="Longueur de corrélation 2")
points(0.3,0.4,pch=5,col=4)
points(abs(sc_separes$Matrice_param), main="Paramètres estimés", xlim=c(0,1), ylim=c(0,1), pch="+",col=2)
points(abs( sc_separes$Matrice_param[ which.min(sc_separes$vect_anderson_darling) ,1 ] ), abs( sc_separes$Matrice_param [ which.min(sc_separes$vect_anderson_darling) ,2 ] ) , pch=5, col=2)

legend("topright", legend= c("Paramètres initiaux","Paramètres estimés","Optimum", "Paramètres réels"), col=c(1,2,2,4),pch=c(3,3,5,5 ))


## dev.copy(pdf,"nuage_croix_separes_0.3_0.4.pdf")
## dev.off()

## write.matrix(sc_separes$Matrice_param,"nuage_croix_separes_0.3_0.4_param_estimes.dat", sep=";" )
## write.matrix(sc_separes$Matrice_param_init,"nuage_croix_separes_0.3_0.4_param_init.dat", sep=";" )
## write.matrix(sc_separes$vect_anderson_darling,"nuage_croix_separes_0.3_0.4_anderson-darling.dat", sep=";" )
## write.matrix(sc_separes$vect_variance,"nuage_croix_separes_0.3_0.4_variance.dat", sep=";" )
## write.matrix(sc_separes$vect_warning,"nuage_croix_separes_0.3_0.4_warning.dat", sep=";") 



# Méthode avec ensembles de test confondus
plot(abs(sc_confondus$Matrice_param_init), main="Paramètres estimés par nuage et croix d'apprentissage confondus", xlim=c(0,1), ylim=c(0,1), pch="+",xlab="Longueur de corrélation 1", ylab="Longueur de corrélation 2")
points(0.3,0.4,pch=5,col=4)
points(abs(sc_confondus$Matrice_param), main="Paramètres estimés", xlim=c(0,1), ylim=c(0,1), pch="+",col=2)
points(abs( sc_confondus$Matrice_param[ which.min(sc_confondus$vect_anderson_darling) ,1 ] ), abs( sc_confondus$Matrice_param [ which.min(sc_confondus$vect_anderson_darling) ,2 ]), pch=5, col=2)

legend("topright", legend= c("Paramètres initiaux","Paramètres estimés","Optimum", "Paramètres réels"), col=c(1,2,2,4),pch=c(3,3,5,5 ))


## dev.copy(pdf,"nuage_croix_confondus_0.3_0.4.pdf")
## dev.off()


## write.matrix(sc_confondus$Matrice_param,"nuage_croix_confondus_0.3_0.4_param_estimes.dat", sep=";" )
## write.matrix(sc_confondus$Matrice_param_init,"nuage_croix_confondus_0.3_0.4_param_init.dat", sep=";" )
## write.matrix(sc_confondus$vect_anderson_darling,"nuage_croix_confondus_0.3_0.4_anderson-darling.dat", sep=";" )
## write.matrix(sc_confondus$vect_variance,"nuage_croix_confondus_0.3_0.4_variance.dat", sep=";" )
## write.matrix(sc_confondus$vect_warning,"nuage_croix_confondus_0.3_0.4_warning.dat", sep=";" )





