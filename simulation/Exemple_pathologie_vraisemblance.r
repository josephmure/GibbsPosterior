source("MaxVraisemblanceIntegree.r")
source("Briques_these.r")
source("Krige_these.r")
source("Matern_these.r")








Noyau = creeMaternIsotrope(variance=1,0.4,0.5)


x_connus = matrix(0.02 * (1:100)) 
y_connus = creeEchantillonNormal(x_connus,Noyau)

#x_connus = as.matrix(read.table("x_connus.txt")) 
#y_connus = as.matrix(read.table("y_connus.txt")) 

abscisses =  1:100 * 0.1 

likelihood = apply(matrix(abscisses),1,evalueVraisemblanceIntegreeAnisGeom, x_connus,y_connus,0.5,cbind(x_connus^0,x_connus)) 

likelihood = likelihood/min(likelihood)


plot(abscisses,likelihood,main="Exponential kernel (100 observations) and affine trend",xlab="theta",ylab="Likelihood",type="l",lwd=2) 

dev.copy(pdf,"Pathologie.pdf", title="Vraisemblance pathologique")
dev.off() 

library(MASS)
write.matrix(x_connus,"x_connus.txt")
write.matrix(y_connus,"y_connus.txt")

plot(x_connus,y_connus,type="l",xlab="",ylab="")

dev.copy(pdf,"Processus_pathologique.pdf",title="Processus gaussien pathologique")
dev.off()


