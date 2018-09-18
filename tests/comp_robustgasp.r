library(RobustGaSP)
#------------------------
# a 3 dimensional example
#------------------------
# dimensional of the inputs
dim_inputs <- 3
# number of the inputs
num_obs <- 30
# uniform samples of design
input <- matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs)
# Following codes use maximin Latin Hypercube Design, which is typically better than uniform
# library(lhs)
# input <- maximinLHS(n=num_obs, k=dim_inputs)  ##maximin lhd sample
####
# outputs from the 3 dim detpep10curv function
output = matrix(0,num_obs,1)
for(i in 1:num_obs){
output[i]<-dettepepel.3.data(input[i,])
}
# use constant mean basis, with no constraint on optimization
m1<- rgasp(design = input, response = output, lower_bound=FALSE, prior_choice = "ref_gamma", nugget.est = FALSE)
z1<- rgasp(design = input, response = output, lower_bound=FALSE, prior_choice = "ref_gamma", nugget.est = FALSE,zero.mean = "Yes")


library(GibbsPosterior)


REGULARITE = 2.5
NB_POINTS_A_GENERER = 500
ITERATIONS = 100
SD_INSTRUM = 0.2
NOMBRE_FONCTIONS_TENDANCE = 1
RANDOM_SEED = 42
WORDS = c("tensorise","REML")
POINT_DEPART <- c(1,1,1)

tendance = matrix(rep(1,length(output)))

g1 = GibbsPosteriorC(dim_inputs,REGULARITE,NB_POINTS_A_GENERER,ITERATIONS,SD_INSTRUM,NOMBRE_FONCTIONS_TENDANCE,RANDOM_SEED,WORDS,input,output,tendance,POINT_DEPART)

posterior = matrix(g1$posterior,ncol=dim_inputs,byrow=TRUE)

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


mode_GibbsRefPost <- optim(par=c(1,1,1), fn=densKernel, fenetre= max(apply(posterior,2,sd))/2, matrice_echantillon=posterior)

MLE <- optim(par = c(1,1,1), fn=evalueOpposeLogVraisemblanceIntegreeTens, x_connus= input, y_connus=output, regularite=REGULARITE, tendance=tendance)


test_points <- matrix(runif(num_obs*dim_inputs*5), num_obs*5,dim_inputs)
tendance_test_points <- matrix(rep(1,nrow(test_points)))

pred_m1 <- predict(m1,test_points,tendance_test_points)
pred_z1 <- predict(z1,test_points)
m1_corr_length = 1 / m1@beta_hat
z1_corr_length = 1 / z1@beta_hat

library(GibbsPosteriorRequired)

NOYAU_RobustGaSP_gamma <- creeMaternTensorise(variance = 1, longueur = m1_corr_length, regularite = 2.5 )  
NOYAU_RobustGaSP_gamma_KRISIMPLE <- creeMaternTensorise(variance = 1, longueur = m1_corr_length, regularite = 2.5 )  


MATRICE_RobustGaSP_gamma <- creeMatriceCovariance(x1= input, x2= input, noyau= NOYAU_RobustGaSP_gamma)
MATRICE_RobustGaSP_gamma_KRISIMPLE <- creeMatriceCovariance(x1= input, x2= input, noyau= NOYAU_RobustGaSP_gamma_)

MATRICE_RobustGaSP_gamma_INVERSE <- solve(MATRICE_RobustGaSP_gamma) ## Dommage qu'on ait a l'inverser, mais la matrice inverse intervient souvent, donc on n'a guere le choix.

CORREL_CONNUS_NOUVEAUX_RobustGaSP_gamma_KRISIMPLE <- creeMatriceCovariance(x1= input, x2= test_points, noyau= NOYAU_RobustGaSP_gamma)

  QR <- qr.Q(qr(tendance),complete=TRUE) #les NOMBRE_FONCTIONS_TENDANCE premieres colonnes de QR contiennent P, les dernieres W
  injectionOrthogonalTendance <- QR[,(NOMBRE_FONCTIONS_TENDANCE+1):ncol(QR)] #W
  injectionTendance <- QR[,1:NOMBRE_FONCTIONS_TENDANCE] #P
  
TENDANCE_NOUVEAUX_CONNUS <- tendance_test_points %*% solve(t(injectionTendance) %*% tendance, t(injectionTendance)) ## Attention, c'est H_0_0 %*% (P^T H)^(-1) %*% P^T = H_0_point %*% P^T

#MATRICE_NOUVEAUX_RobustGaSP_gamma <- 1 + TENDANCE_NOUVEAUX_CONNUS %*% MATRICE_RobustGaSP_gamma %*% t(TENDANCE_NOUVEAUX_CONNUS) - TENDANCE_NOUVEAUX_CONNUS %*% CORREL_CONNUS_NOUVEAUX_RobustGaSP_gamma - t(CORREL_CONNUS_NOUVEAUX_RobustGaSP_gamma) %*% t(TENDANCE_NOUVEAUX_CONNUS)

  CORREL_CONNUS_NOUVEAUX_RobustGaSP_gamma <- t(injectionOrthogonalTendance) %*% (MATRICE_RobustGaSP_gamma %*% t(TENDANCE_NOUVEAUX_CONNUS) - CORREL_CONNUS_NOUVEAUX_RobustGaSP_gamma_KRISIMPLE)
  Estimation_beta_RobustGaSP_gamma  <- solve(t(tendance) %*% solve(MATRICE_RobustGaSP_gamma) %*% tendance) %*% t(tendance) %*% solve(MATRICE_RobustGaSP_gamma) %*% output

  MATRICE_RobustGaSP_gamma_INJECTEE <- t(injectionOrthogonalTendance) %*% MATRICE_RobustGaSP_gamma %*% injectionOrthogonalTendance
  MATRICE_RobustGaSP_gamma_INVERSE <- solve(MATRICE_RobustGaSP_gamma_INJECTEE)
  MATRICE_RobustGaSP_gamma_INVERSE_KRISIMPLE <- solve(MATRICE_RobustGaSP_gamma)

    Moyenne_marginale <- t( TENDANCE_NOUVEAUX_CONNUS %*% output )
  y_connus <- - t(injectionOrthogonalTendance) %*% output ## Attention, on se base sur MOINS W^T y
  Prediction_RobustGaSP_gamma <- as.vector(t(y_connus) %*% MATRICE_RobustGaSP_gamma_INVERSE %*% CORREL_CONNUS_NOUVEAUX_RobustGaSP_gamma + Moyenne_marginale)
  Prediction_RobustGaSP_gamma_KRISIMPLE <- as.vector(t(output) %*% MATRICE_RobustGaSP_gamma_INVERSE_KRISIMPLE %*% CORREL_CONNUS_NOUVEAUX_RobustGaSP_gamma_KRISIMPLE)
 
  
#   MATRICE_NOUVEAUX_RobustGaSP_xi_KRISIMPLE <- 1  