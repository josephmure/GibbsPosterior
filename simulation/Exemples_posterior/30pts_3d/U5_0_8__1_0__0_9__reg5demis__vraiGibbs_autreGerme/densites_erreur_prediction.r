LONGUEUR_CORRELATION <- scan("longueur_correlation_vraie.txt")

NOM_MLE <- "MLE"
NOM_MAP <- "MAP"
TITRE_GRAPHE <- paste("Gap between prediction and correct conditional mean :", "\n", LONGUEUR_CORRELATION[1], "-", LONGUEUR_CORRELATION[2], "-", LONGUEUR_CORRELATION[3])
LEGENDE_ABSCISSES <- "Error (quadratic norm)"
LEGENDE_ORDONNEES <- "Probabilty density"

Ecarts_predictions <- read.table("Ecarts_predictions.txt", as.is = TRUE)
    

Ecarts_predictions_MLE <- Ecarts_predictions[,1]
Ecarts_predictions_MAP <- Ecarts_predictions[,2]

Densite_MLE <- density(Ecarts_predictions_MLE)
Densite_MAP <- density(Ecarts_predictions_MAP)


##Creation du plot (densite de proba de l'erreur sur la prédiction (vis-à-vis de la juste moyenne conditionnelle)

 
plot(x=c(0,max(Ecarts_predictions)),y=c(0, max(c(Densite_MLE$y, Densite_MAP$y))), type="n", main=TITRE_GRAPHE, xlab=LEGENDE_ABSCISSES, ylab=LEGENDE_ORDONNEES)
lines(Densite_MLE)
lines(Densite_MAP, lwd=3)
legende <- c(NOM_MLE,NOM_MAP)
lwd <- c(1,3)
legend("topright", legende, lwd=lwd)
dev.copy(pdf,"Densite_ecarts_predictions.pdf")
dev.off()





