library(ggplot2)

LONGUEUR_CORRELATION <- scan("longueur_correlation_vraie.txt")

NOM_MLE <- "MLE"
NOM_MAP <- "MAP"
TITRE_GRAPHE <- paste("Gap between the estimators and the correct correlation lengths :", LONGUEUR_CORRELATION[1], "-", LONGUEUR_CORRELATION[2], "-", LONGUEUR_CORRELATION[3])

Distances_simples <- read.table("Distances_simples.txt", as.is = TRUE)
Distances_fisher <- read.table("Distances_fisher.txt", as.is=TRUE)


dat_simples<- data.frame(Distance_Simple= c(Distances_simples$V1, Distances_simples$V2), Estimator=c( rep(NOM_MLE,length(Distances_simples$V1)), rep(NOM_MAP, length(Distances_simples$V2)) ))
    




dat_fisher<- data.frame(Distance_Fisher= c(Distances_fisher$V1, Distances_fisher$V2), Estimator=c( rep(NOM_MLE,length(Distances_fisher$V1)), rep(NOM_MAP, length(Distances_fisher$V2)) ))
    





plot_simple <- ggplot(dat_simples, aes(x=Distance_Simple, fill=Estimator)) +
    geom_histogram(alpha=0.5, position="identity") +
    ggtitle(TITRE_GRAPHE) +
    scale_fill_manual(values = c("MLE" = "red","MAP" = "blue"))

       

plot_fish <- ggplot(dat_fisher, aes(x=Distance_Fisher, fill=Estimator)) +
    geom_histogram(alpha=0.5, position="identity") +
    ggtitle(TITRE_GRAPHE) +
    scale_fill_manual(values = c("MLE" = "red","MAP" = "blue"))



plot_simple
dev.copy(pdf,"Histogramme_dist_simple.pdf")
dev.off()

plot_fish
dev.copy(pdf,"Histogramme_dist_fisher.pdf")
dev.off()

