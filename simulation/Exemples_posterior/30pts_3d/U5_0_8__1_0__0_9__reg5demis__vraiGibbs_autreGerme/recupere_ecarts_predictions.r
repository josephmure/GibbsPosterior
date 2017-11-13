library(pscl)


# Ecarts entre la prediction des estimateurs MLE et MAP et la juste moyenne conditionnelle sachant les observations

Mat1 = as.matrix(read.table("planXPvariable/matrice_ecarts_predictions.txt",as.is=TRUE))
Mat2 = as.matrix(read.table("planXPvariable_2/matrice_ecarts_predictions.txt",as.is=TRUE))
Mat3 = as.matrix(read.table("planXPvariable_3/matrice_ecarts_predictions.txt",as.is=TRUE))
Mat4 = as.matrix(read.table("planXPvariable_4/matrice_ecarts_predictions.txt",as.is=TRUE))
Mat5 = as.matrix(read.table("planXPvariable_5/matrice_ecarts_predictions.txt",as.is=TRUE))

Matrice_ecarts_predictions = rbind(Mat1,Mat2,Mat3,Mat4,Mat5)
colnames(Matrice_ecarts_predictions) = NULL

write.matrix(x=Matrice_ecarts_predictions,file= "Ecarts_predictions.txt", sep= "\t")



