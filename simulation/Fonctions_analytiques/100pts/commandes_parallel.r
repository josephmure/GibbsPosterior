library(pscl)
library(lhs)

source("../Rastrigin.r")

source("../../generePlanXP.r")
source("../../genereMatriceTendance.r")

library(Rcpp)
library(RcppGSL)
#sourceCpp("../../RcppGibbsPosterior.cpp")
library(GibbsPosteriorC)

source("../../trouveMLE.r")
source("../../postTraitement.r")

NOMBRE_PROCESSUS <-  500

#NOMBRE_DIMENSIONS <- 7
REGULARITE <- 2.5
NOMBRE_POINTS_POSTERIOR_A_GENERER <- 1000
NOMBRE_PAS_METROPOLIS <- 50 
ECART_TYPE_METROPOLIS <- 0.2



#TYPE_NOYAU_MATERN_TYPE_PRIOR <- scan("../type_noyau_prior.txt", what="character") ## interessant, mais en pratique c'est toujours la meme chose
TYPE_NOYAU_MATERN_TYPE_PRIOR <- c("geometrique","REML")
TYPE_NOYAU_MATERN <- TYPE_NOYAU_MATERN_TYPE_PRIOR[1]
TYPE_PRIOR <- TYPE_NOYAU_MATERN_TYPE_PRIOR[2]

#VECTEUR_GERMES_ALEATOIRES <- NULL
#GERME_ALEATOIRE_GLOBAL <- scan("germe_aleatoire.txt") # ???

MLE <- NULL
MAP <- NULL
MAP_optim_classique <- NULL


## Inutile : on n'ecrit plus de donnees intermediaires sur le disque dur
## NOUVEAU : Suppression des fichiers au début
#system("rm Posterior_anis_geom.txt")
#system("rm Posterior_tens.txt")
## FIN NOUVEAU

## "Pour du beurre", on genere des observations (et donc un plan d'experiences afin de connaitre la dimension)
#germe_aleatoire <- 0
#source("../genereObservations.r")
#germe_aleatoire <- NULL

#Simulations <- matrix(NA,nrow=NOMBRE_PROCESSUS,ncol=(2 + NOMBRE_POINTS_POSTERIOR_A_GENERER) * ncol(x_connus))

library(SparkR)


#Ce wrapper sert a etendre le timeout de SparkR (tres court par defaut) a 72h
connectBackend.orig <- getFromNamespace('connectBackend', pos='package:SparkR')
connectBackend.patched <- function(hostname, port, timeout = 3600*72) {
   connectBackend.orig(hostname, port, timeout)
}

assignInNamespace("connectBackend", value=connectBackend.patched, pos='package:SparkR')
#Fin wrapper

#sparkR.session("local[*]")
sparkR.session("spark://149.251.6.133:7077", appName="Simulations")


faitSimulations <- function(germe_aleatoire)
{
  library(GibbsPosterior)
  #  print(paste("Germe aleatoire", germe_aleatoire) )
  #	VECTEUR_GERMES_ALEATOIRES <- c(VECTEUR_GERMES_ALEATOIRES, germe_aleatoire)
  #	write.matrix(VECTEUR_GERMES_ALEATOIRES, "liste_germes_aleatoires.txt")
  #	source("../../nettoyage.r") #inutile : on n'ecrit plus de donnees intermediaires sur le disque dur
  #	source("../../generePlanXP.r") #inutile : deja contenu dans genereObservations.r
  x_connus <- generePlanXP(germe_aleatoire,NOMBRE_POINTS_PLANXP = 100, NOMBRE_DIMENSIONS, LHS=FALSE)
  y_connus <- Rastrigin(x_connus)
  #	source("../genereObservations.r")
  #	source("../../genereMatriceTendance.r")
  tendance <- genereMatriceTendance(x_connus,FONCTIONS)
  #	source("../../scriptR.r")
  
  #sourceCpp("../../RcppGibbsPosterior.cpp")
  #quarantedeux()
  
  point_posterior <- GibbsPosteriorC(ncol(x_connus), REGULARITE , NOMBRE_POINTS_POSTERIOR_A_GENERER, NOMBRE_PAS_METROPOLIS, ECART_TYPE_METROPOLIS, length(FONCTIONS),germe_aleatoire,TYPE_NOYAU_MATERN_TYPE_PRIOR, x_connus, y_connus, tendance)
  
  
  Posterior <- matrix(point_posterior,ncol=ncol(x_connus),byrow=TRUE)
  
  
  #source("../../trouveMLE.r")
  
  argmax_vraisemblance_integree <- trouveMLE(NOMBRE_DIMENSIONS= ncol(x_connus),TYPE_NOYAU_MATERN, REGULARITE,x_connus, y_connus,tendance)
  mode_posterior <- trouveMAP(FENETRE_opt_MAP = 1, Posterior, NOMBRE_DIMENSIONS = ncol(x_connus))  
  
  
  #	source("../../postTraitement.r")
  
  
  # 	MLE <- rbind(MLE, argmax_vraisemblance_integree)
  # 	write.matrix(MLE,"liste_MLE.txt",sep = "\t")
  # 	MAP <- rbind(MAP, mode_posterior)
  # 	write.matrix(MAP,"liste_MAP.txt",sep = "\t")
  
  #Simulations[germe_aleatoire,] <- c(argmax_vraisemblance_integree,mode_posterior,point_posterior)
  list(MLE = argmax_vraisemblance_integree,MAP = mode_posterior,Posterior = matrix(point_posterior,ncol=NOMBRE_DIMENSIONS,byrow=TRUE), planXP = x_connus, observations = y_connus, tendance = tendance) 


  }

Simulations <- spark.lapply((GERME_ALEATOIRE_GLOBAL+1):(GERME_ALEATOIRE_GLOBAL+NOMBRE_PROCESSUS) , faitSimulations)


sparkR.stop()

save.image(file = "Simulations.RData",safe = TRUE)
