LONGUEUR_CORRELATION <- scan("longueur_correlation_vraie.txt")

NOM_MLE <- "MLE"
NOM_MAP <- "MAP"
TITRE_GRAPHE <- paste("Gap between the estimators and the correct correlation lengths :", "\n", LONGUEUR_CORRELATION[1], "-", LONGUEUR_CORRELATION[2], "-", LONGUEUR_CORRELATION[3])
LEGENDE_ABSCISSES_DISTANCE_SIMPLE <- "Error (simple distance)"
LEGENDE_ABSCISSES_DISTANCE_FISHER <- "Error (Fisher distance)"
LEGENDE_ORDONNEES <- "Probabilty density"

Distances_simples <- read.table("Distances_simples.txt", as.is = TRUE)
Distances_fisher <- read.table("Distances_fisher.txt", as.is=TRUE)


#dat_simples<- data.frame(Distance_Simple= c(Distances_simples$V1, Distances_simples$V2), Estimator=c( rep(NOM_MLE,length(Distances_simples$V1)), rep(NOM_MAP, length(Distances_simples$V2)) ))
    

#dat_fisher<- data.frame(Distance_Fisher= c(Distances_fisher$V1, Distances_fisher$V2), Estimator=c( rep(NOM_MLE,length(Distances_fisher$V1)), rep(NOM_MAP, length(Distances_fisher$V2)) ))
    

Distance_simple_MLE <- Distances_simples[,1]
Distance_simple_MAP <- Distances_simples[,2]

Distance_fisher_MLE <- Distances_fisher[,1]
Distance_fisher_MAP <- Distances_fisher[,2]

Densite_simple_MLE <- density(Distance_simple_MLE)
Densite_simple_MAP <- density(Distance_simple_MAP)

Densite_fisher_MLE <- density(Distance_fisher_MLE)
Densite_fisher_MAP <- density(Distance_fisher_MAP)

##Creation des plots (densite de proba de l'erreur sur l'estimation des longueurs de correlation

## Distance simple 
plot(x=c(0,max(Distances_simples)),y=c(0, max(c(Densite_simple_MLE$y, Densite_simple_MAP$y))), type="n", main=TITRE_GRAPHE, xlab=LEGENDE_ABSCISSES_DISTANCE_SIMPLE, ylab=LEGENDE_ORDONNEES)
lines(Densite_simple_MLE) #, sub=legende)
#abline(v=mean(Distance_simple_MLE))
lines(Densite_simple_MAP, lwd=3)
#abline(v=mean(Distance_simple_MAP), lwd=3)
legende <- c(NOM_MLE,NOM_MAP)
lwd <- c(1,3)
legend("topright", legende, lwd=lwd)
dev.copy(pdf,"Densite_dist_simple.pdf")
dev.off()

## Distance Fisher
plot(x=c(0,max(Distances_fisher)),y=c(0, max(c(Densite_fisher_MLE$y, Densite_fisher_MAP$y))), type="n", main=TITRE_GRAPHE, xlab=LEGENDE_ABSCISSES_DISTANCE_FISHER, ylab=LEGENDE_ORDONNEES)
lines(Densite_fisher_MLE) #, sub=legende)
lines(Densite_fisher_MAP, lwd=3)
legende <- c(NOM_MLE,NOM_MAP)
lwd <- c(1,3)
legend("topright", legende, lwd=lwd)
dev.copy(pdf,"Densite_dist_fisher.pdf")
dev.off()


#ggplot(dat_simples, aes(x=Distance_Simple, fill=Estimator)) +
#    geom_histogram(alpha=0.5, position="identity", binwidth = 1) +
#    ggtitle(TITRE_GRAPHE) +
#    scale_fill_manual(values = c("MLE" = "red","MAP" = "blue"))

       

#plot_fish <- ggplot(dat_fisher, aes(x=Distance_Fisher, fill=Estimator)) +
#    geom_histogram(alpha=0.5, position="identity", binwidth = 1.5) +
#    ggtitle(TITRE_GRAPHE) +
#    scale_fill_manual(values = c("MLE" = "red","MAP" = "blue"))





