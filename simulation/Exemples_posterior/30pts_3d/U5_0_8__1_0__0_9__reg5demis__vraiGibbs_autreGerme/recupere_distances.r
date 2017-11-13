library(pscl)


# Pour les distances simples

Mat1 = as.matrix(read.table("planXPvariable/matrice_distances.txt",as.is=TRUE))
Mat2 = as.matrix(read.table("planXPvariable_2/matrice_distances.txt",as.is=TRUE))
Mat3 = as.matrix(read.table("planXPvariable_3/matrice_distances.txt",as.is=TRUE))
Mat4 = as.matrix(read.table("planXPvariable_4/matrice_distances.txt",as.is=TRUE))
Mat5 = as.matrix(read.table("planXPvariable_5/matrice_distances.txt",as.is=TRUE))

Matrice_distances_simples = rbind(Mat1,Mat2,Mat3,Mat4,Mat5)
colnames(Matrice_distances_simples) = NULL

write.matrix(x=Matrice_distances_simples,file= "Distances_simples.txt", sep= "\t")







# Pour les distances Fisher

Mat1 = as.matrix(read.table("planXPvariable/matrice_distances_fisher.txt",as.is=TRUE))
Mat2 = as.matrix(read.table("planXPvariable_2/matrice_distances_fisher.txt",as.is=TRUE))
Mat3 = as.matrix(read.table("planXPvariable_3/matrice_distances_fisher.txt",as.is=TRUE))
Mat4 = as.matrix(read.table("planXPvariable_4/matrice_distances_fisher.txt",as.is=TRUE))
Mat5 = as.matrix(read.table("planXPvariable_5/matrice_distances_fisher.txt",as.is=TRUE))

Matrice_distances_fisher = rbind(Mat1,Mat2,Mat3,Mat4,Mat5)
colnames(Matrice_distances_fisher) = NULL

write.matrix(x=Matrice_distances_fisher,file= "Distances_fisher.txt", sep= "\t")

