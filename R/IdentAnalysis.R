############################################################
## Functions related to identifiability analysis
############################################################

## this is a convenient function that computes J.C.F. over the field
## of real numbers. Please note that this function is *not*
## numerically stable for all matrices (esp. large matrices). It
## assumes that there is no nilpotent cells in the J.C.F.  Returned
## values: K1 is the number of real eigenvalues; K2 is the number of
## pairs of complex eigenvalues. J is the pxp semi-diagonal matrix of
## Jordan blocks; Qmat is the matrix of generalized eigenvectors. Note
## that J is always organized in such way: the first K1 diagonal
## elements are real eigenvalues; the rest are 2x2 rotational matrices
## correspond with pairs of complex eigenvalues.
JordanReal <- function(A,tol=1e-6) {
  n <- nrow(A); m <- ncol(A)
  if (n !=m) stop("A must be a square matrix!")
  ## 1. J.C.F. on C
  ee <- eigen(A); lambdas <- ee$values
  rr <- Re(lambdas); ii <- Im(lambdas)
  xx <- Re(ee$vectors); yy <- Im(ee$vectors)
  ## 2. classify eigenvalues into real ones and *pairs* of complex
  ## ones
  real.idx <- which(abs(ii)<tol); comp.idx <- which(abs(ii)>=tol)
  K1 <- length(real.idx); K2 <- length(comp.idx)/2
  ## 3. the Jordan blocks and eigen.vecs; the real ones
  J <- Qmat <- matrix(0, n, n)
  if (K1>0) {
    for (k in 1:K1) {
      J[k,k] <- rr[real.idx[k]]
      Qmat[,k] <- xx[,real.idx[k]]
    }
  }
  ## 4. for complex pairs
  if (K2>0) {
    for (k in 1:K2) {
      k2 <- 2*k+K1-1
      J[k2,k2] <- J[k2+1,k2+1] <- rr[comp.idx[2*k]]
      J[k2,k2+1] <- ii[comp.idx[2*k-1]]
      J[k2+1,k2] <- ii[comp.idx[2*k]]
      Qmat[,k2] <- xx[,comp.idx[2*k-1]]
      Qmat[,k2+1] <- yy[,comp.idx[2*k-1]]
    }
  }
  return(list(J=J, Qmat=Qmat, Qinv=rsolve(Qmat), lambdas=lambdas, K1=K1, K2=K2))
}


## the minimum gap between eigenvalues
Lgapfun <- function(lambdas) {
  aa <- Re(lambdas); bb <- Im(lambdas)
  Lgap <- min(dist(cbind(aa,bb)))^2
  return(Lgap)
}

## this function computes the initial condition-based identifiability
## score (ICIS) and a few related quantitative measures of
## identifiability of an ODE system at initial conditon x0. n.digit is
## the number of precision digits to keep.
ICISAnalysis <- function(A, x0, n.digits=3) {
  ## sanity checks
  p <- nrow(A); p2 <- ncol(A); p3 <- length(x0)
  if (p !=p2) {
    stop("A must be a square matrix!")
  } else if (p != p3) {
    stop("Length of x0 must equal the dimension of A!")
  }
  ## 1. compute J.C.F. of A.
  jcf <- JordanReal(A); Qinv <- jcf$Qinv
  K1 <- jcf$K1; K2 <- jcf$K2; K <- K1+K2
  ## 2. Check if there are repeated real/complex eigenvalues. Note
  ## that eigenvalues returned by eigen() are *ordered*.
  lambdas <- round(eigen(A)$values, n.digits)
  ## calculate runs (repeats). To construct one equivalent network, We
  ## only need to use one of (possibly) many runs. I decide to use the
  ## first one.
  runs <- rle(lambdas)$lengths
  if (max(runs)==1) {
    Ident1 <- TRUE; RepEigenIdx <- NULL
  } else {
    Ident1 <- FALSE
    r0 <- min(which(runs>1)); r1 <- r0+runs[r0]-1
    RepEigenIdx <- r0:r1
  }
  ## also return Lgap, the minimum gap between eigenvalues, as an 
  ## identifiability measure, to the users
  Lgap <- Lgapfun(lambdas)
  ## 3. compute |w0,k|
  xtilde0 <- drop(Qinv %*% x0)
  w0k.norm <- rep(0, K)
  if (K1>0) w0k.norm[1:K1] <- round(abs(xtilde0[1:K1]),n.digits)
  if (K2>0) {
    for (k in 1:K2){
      k2 <- K1+2*k-1
      w0k.norm[K1+k] <- round(sqrt(sum(xtilde0[c(k2,k2+1)]^2)),n.digits)
    }
  }
  ICIS <- min(w0k.norm)
  ## 4. A is identifiable at x0 if: a) dlambda > tol, b) min(w0k.norm) > tol. 
  Ident2 <- ICIS>0
  ident <- Ident1 & Ident2
  ## 5. return useful results for constructing the identifiability
  ## class [A, x0]
  a0 <- ifelse(w0k.norm>0,0,1)
  if (K1>0) {
    aa.real <- a0[1:K1]
  } else {
    aa.real <- NULL
  }
  if (K2>0) {
    aa.comp <- rep(a0[(K1+1):K],each=2)
  } else {
    aa.comp <- NULL
  }
  aa <- c(aa.real, aa.comp)
  I0 <- diag(aa); Iplus <- diag(p)-I0
  return(list(Identifiable=ident, Ident1=Ident1, Ident2=Ident2, 
              ICIS=ICIS, w0k.norm=w0k.norm, I0=I0, Iplus=Iplus,
              Jordan=list(J=jcf$J, Qmat=jcf$Qmat, Qinv=jcf$Qinv, K1=K1, K2=K2,
              Lgap=Lgap, RepEigenIdx=RepEigenIdx)))
}

## EquivSystem computes [A, x0], the equivalence class of A given x0, and
## then creates an equivalent system A (for the same x0) based on the
## specified difference matrix Dmat (with the same dimensionality of A).
## Remarks: a) when Dmat is the zero matrix, EquivSystem(A, x0, 0)=A; b)
## when A is *identifiable* at x0, EquivSystem(A, x0, Dmat)=A for all input
## Dmat.  For unidentifiable cases, using a random Dmat a.s. will lead to a
## different system Atilde that is equivalent to A.
EquivSystem <- function(A, x0, Dmat, ...) {
  ## sanity checks
  if (!all(dim(A)==dim(Dmat))) stop("The dimension of A and Dmat must be equal!")
  p <- nrow(A); p2 <- ncol(A); p3 <- length(x0)
  if (p !=p2) {
    stop("A must be a square matrix!")
  } else if (p != p3) {
    stop("Length of x0 must equal the dimension of A!")
  }
  ## Identifiability analysis for (A, x0)
  o <- ICISAnalysis(A, x0, ...)
  Qmat <- o$Jordan$Qmat; Qinv <- o$Jordan$Qinv; I0 <- o$I0
  ## compute the equivalent network
  Atilde  <- A + Qmat%*%I0%*%Dmat%*%I0%*%Qinv
  return(Atilde)
}


## w* provides an approximation of MSE(A) based on the observed
## discrete data
wstarfun <- function(Y,S,L) {
  d <- nrow(Y)
  N <- rsolve(Y%*%S%*%t(Y)); A <- Y%*%L%*%t(Y)%*%N
  N2 <- N%*%N; tAA <- crossprod(A)
  Term1 <- sum(as.vector(t(N2))*as.vector(Y%*%L%*%L%*%t(Y) +d*crossprod(L%*%t(Y)) +Y%*%t(L)%*%t(L)%*%t(Y) -2*A%*%Y%*%S%*%L%*%t(Y) -2*tr(A)*Y%*%t(S)%*%L%*%t(Y) -2*Y%*%t(S)%*%t(L)%*%t(Y)%*%A +tAA%*%Y%*%S%*%S%*%t(Y) +tr(tAA)*crossprod(S%*%t(Y)) +Y%*%t(S)%*%t(S)%*%t(Y)%*%tAA))
  Term2 <- tr(crossprod(Y%*%L) -2*t(S)%*%t(Y)%*%t(A)%*%Y%*%L +crossprod(A%*%Y%*%S))*tr(N2)
  return(Term1+Term2)
}
## kappa is the metric proposed by Stanhope
kappafun <- function(Y) kappa(Y[, 1:nrow(Y)], exact=TRUE)
## tau is a generalization based on the smoothed 2stage
taufun <- function(Y,S) kappa(Y%*%S%*%t(Y), exact=TRUE)

## this is a wrapper that computes all three practical identifiability
## scores in one function
PISAnalysis <- function(Y,S,L){
  PIS <- wstarfun(Y,S,L)
  kappa <- kappafun(Y)
  SCN <- taufun(Y,S)
  return(c(SCN=SCN, PIS=PIS, kappa=kappa))
}
