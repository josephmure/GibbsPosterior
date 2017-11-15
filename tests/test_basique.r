## Test elementaire pour voir si le paquet fonctionne
library(GibbsPosterior)

DIMENSION = 3
REGULARITE = 2.5
NB_POINTS_A_GENERER = 50
ITERATIONS = 100
SD_INSTRUM = 0.2
NB_FONCTIONS_TENDANCE = 1
RANDOM_SEED = 42
WORDS = c("geometrique","REML")
POINT_DEPART <- c(1,1,1)

planXP = matrix(runif(NB_POINTS_A_GENERER * DIMENSION),nrow=NB_POINTS_A_GENERER, ncol=DIMENSION)
observations = rnorm(NB_POINTS_A_GENERER)
tendance = matrix(rep(1,NB_POINTS_A_GENERER))

res = GibbsPosteriorC(DIMENSION,REGULARITE,NB_POINTS_A_GENERER,ITERATIONS,SD_INSTRUM,NB_FONCTIONS_TENDANCE,RANDOM_SEED,WORDS,planXP,observations,tendance,POINT_DEPART)

posterior = matrix(res$posterior,ncol=DIMENSION,byrow=TRUE)

posterior
