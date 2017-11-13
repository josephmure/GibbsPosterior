library(pscl)


# Recupere les valeurs prises par le MLE

Mat1 = as.matrix(read.table("planXPvariable/Posterior_anis_geom.txt",as.is=TRUE))
Mat2 = as.matrix(read.table("planXPvariable_2/Posterior_anis_geom.txt",as.is=TRUE))
Mat3 = as.matrix(read.table("planXPvariable_3/Posterior_anis_geom.txt",as.is=TRUE))
Mat4 = as.matrix(read.table("planXPvariable_4/Posterior_anis_geom.txt",as.is=TRUE))
Mat5 = as.matrix(read.table("planXPvariable_5/Posterior_anis_geom.txt",as.is=TRUE))

Matrice_echantillons_posterior = rbind(Mat1,Mat2,Mat3,Mat4,Mat5)
colnames(Matrice_echantillons_posterior) = NULL

write.matrix(x=Matrice_echantillons_posterior,file= "matrice_echantillons_posterior.txt", sep= "\t")

