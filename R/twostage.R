## 02/08/2019: This is Xing's implementation of the two-stage method based on
## smoothing splines provided by the fda package
fsmooth <- function(Y, Ts, rough.pen=.001, norder=4, plot=FALSE, n.plot=5){
  ## represent the discrete data as functions
  mybasis <- create.bspline.basis(range(Ts), length(Ts)+norder-2, norder, Ts)
  mypar <- fdPar(mybasis, 2, lambda=rough.pen)  #under-smooth
  Xt <- smooth.basis(Ts, t(Y), mypar)[["fd"]]
  if (plot==TRUE) {
    n2 <- min(nrow(Y), n.plot)
    matplot(Ts, t(Y[1:n2,]), xlab="t", ylab="x(t)")
    plot(Xt[1:n2], add=TRUE)
  }
  return(Xt)
}

######################################################################
## 06/11/2020: Two relatively simple two-stage methods that are used
## in the identifiabililty paper. These two methods motivated the
## development of the wstar statistic (see IdentAnalysis1.r), which is
## an approximation of MSE(Ahat).
######################################################################
## Example ex:simple-2stage in the paper. 
twostage1 <- function(Y, tstep) {
  d <- nrow(Y); n <- ncol(Y)
  L <- cbind(-rbind(diag(n-1),0) +rbind(0, diag(n-1)),0)/tstep
  Ahat <- Y%*%L%*%t(Y) %*%rsolve(tcrossprod(Y))
  return(list(Ahat=Ahat, S=diag(n), L=L, x0=Y[1,]))
}

## Example ex:smooth-2stage in the paper.
twostage2 <- function(Y, Ts, nord=4, rough.pen=1e-3){
  ## represent the discrete data as functions
  mybasis <- create.bspline.basis(range(Ts), length(Ts)+nord-2, nord, Ts)
  ## Jb is the pairwise inprod of basis functions
  Jb <- inprod(mybasis, mybasis)
  ## matrix. rep. of D
  Dmat <- coef(deriv.fd(fd(diag(mybasis$nbasis), mybasis)))
  mypar <- fdPar(mybasis, 2, lambda=rough.pen)
  smod <- smooth.basis(Ts, t(Y), mypar); xt.hat <- smod[["fd"]]
  x0 <- drop(eval.fd(Ts[1], xt.hat))
  y2cMap <- smod[["y2cMap"]]
  S <- t(y2cMap)%*%Jb%*%y2cMap
  L <- t(y2cMap)%*%t(Dmat)%*%Jb%*%y2cMap
  Ahat <- Y%*%L%*%t(Y) %*%rsolve(Y%*%S%*%t(Y))
  return(list(Ahat=Ahat, S=S, L=L, xt.hat=xt.hat, x0=x0))
}

