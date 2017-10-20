load("Simulations_copie.RData")


source("../../postTraitement.r")

for(i in 1:length(Simulations))
{
  print(Simulations[[i]]$MAP)
  Simulations[[i]]$MAP <- trouveMAP(FENETRE_opt_MAP = 0.5, Posterior = Simulations[[i]]$Posterior, NOMBRE_DIMENSIONS = ncol(Simulations[[i]]$Posterior), NOMBRE_DEPARTS = 20)
  print(Simulations[[i]]$MAP)
  print("######################")
}