library(pscl)


# Recupere les valeurs prises par le MLE

Mat1 = as.matrix(read.table("planXPvariable/liste_MLE.txt",as.is=TRUE))
Mat2 = as.matrix(read.table("planXPvariable_2/liste_MLE.txt",as.is=TRUE))
Mat3 = as.matrix(read.table("planXPvariable_3/liste_MLE.txt",as.is=TRUE))
Mat4 = as.matrix(read.table("planXPvariable_4/liste_MLE.txt",as.is=TRUE))
Mat5 = as.matrix(read.table("planXPvariable_5/liste_MLE.txt",as.is=TRUE))

Matrice_liste_MLE = rbind(Mat1,Mat2,Mat3,Mat4,Mat5)
colnames(Matrice_liste_MLE) = NULL

write.matrix(x=Matrice_liste_MLE,file= "MLE.txt", sep= "\t")







# Recupere les valeurs prises par le MAP

Mat1 = as.matrix(read.table("planXPvariable/liste_MAP.txt",as.is=TRUE))
Mat2 = as.matrix(read.table("planXPvariable_2/liste_MAP.txt",as.is=TRUE))
Mat3 = as.matrix(read.table("planXPvariable_3/liste_MAP.txt",as.is=TRUE))
Mat4 = as.matrix(read.table("planXPvariable_4/liste_MAP.txt",as.is=TRUE))
Mat5 = as.matrix(read.table("planXPvariable_5/liste_MAP.txt",as.is=TRUE))

Matrice_liste_MAP = rbind(Mat1,Mat2,Mat3,Mat4,Mat5)
colnames(Matrice_liste_MAP) = NULL

write.matrix(x=Matrice_liste_MAP,file= "MAP.txt", sep= "\t")

