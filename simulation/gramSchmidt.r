gramSchmidt <- function (A, tol = .Machine$double.eps^0.5) 
{
    stopifnot(is.numeric(A), is.matrix(A))
    m <- nrow(A)
    n <- ncol(A)
    if (m < n) 
        stop("No. of rows of 'A' must be greater or equal no. of colums.")
    Q <- matrix(0, m, n)
    R <- matrix(0, n, n)
    for (k in 1:n) {
        Q[, k] <- A[, k]
        if (k > 1) {
            for (i in 1:(k - 1)) {
                R[i, k] <- t(Q[, i]) %*% Q[, k]
                Q[, k] <- Q[, k] - R[i, k] * Q[, i]
            }
        }
        R[k, k] <- sqrt( sum( Q[, k] ^2) )
        if (abs(R[k, k]) <= tol) Q[, k] <- 0 
        else Q[, k] <- Q[, k]/R[k, k]
    }
    return(list(Q = Q, R = R))
}
## Le resultat ($Q) est la transposee de la matrice Q_theta (avec des colonnes remplies de zeros).