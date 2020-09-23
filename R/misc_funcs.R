## misc useful functions

tr <- function(A) sum(diag(as.matrix(A)))

## numerically robust matrix inverse
rsolve <- function(mat, d.prop=1e-6, dmin=1e-9, dmax=1e9){
  o <- svd(mat)
  ## singular value threshold
  thresh <- min(max(sum(o$d) * d.prop, dmin), dmax)
  inv.mat <- o$v %*% diag(1/(pmax(o$d, thresh))) %*% t(o$u)
  return(inv.mat)
}

## generate Jordan blocks for a set of given eigenvalues
Jordan <- function(eigenvals){
    ## all Conjugate of eigenvals will be used 
    eigenvals2 <- unique(c(eigenvals, Conj(eigenvals)))
    ndim <- length(eigenvals2)
    eigenvals3 <- eigenvals2[Im(eigenvals2)>=0]
    J <- diag(ndim)
    dj <- 1
    for (lambda in eigenvals3){
        if (Im(lambda)==0){
           J[dj, dj] <- lambda
           dj <- dj+1
       } else {
           J[dj:(dj+1), dj:(dj+1)] <- matrix(c(Re(lambda), Im(lambda),
                                               -Im(lambda), Re(lambda)), nrow=2, byrow=T)
           dj <- dj+2
       }
    }
    return(Re(J))
}


