## Recherche les colonnes identiques dans 2 matrices différentes
## Correction de la fonction row.match du paquet "prodlim"
Row.match <- function (x, table, nomatch = NA) 
{
    if (class(table) == "matrix") 
        table <- as.data.frame(table)
    if (is.null(dim(x))) 
        x <- as.data.frame(matrix(x, nrow = 1))
    if (class(x) == "matrix") 
        x <- as.data.frame(x)
    cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
    ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
    match(cx, ct, nomatch = nomatch)
}

Mode <- function(x)
{
    ux <- unique(x) ## supprime les lignes inutiles
    tab <- tabulate(Row.match(x, ux)) ## pour chaque ligne de ux, indique combien de fois elle apparait dans x
    max <- max(tab)
    indices <- which(tab==max) ## vecteur des indices des lignes correspondant aux modes
    modes <- ux[indices,]
    print(paste("Nb points du mode :",max)) # combien de points coincident avec ce mode
    print(ux[which(tab>0.9*max),])
    modes
}


ModeArrondi<-function(x,nb_decimales=1)
{
    x <- round(x,nb_decimales)
    Mode(x)
}
