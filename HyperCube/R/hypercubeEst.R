#' @export
hyercubeEst <-
function(x, y, V) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  sigma2 <- estSigma(x, y)
  cf <- canonicalForm(x, y, V)
  estrisk <- estRiskCanonical(cf, sigma2)
  list(etahat = cf$etahat, estsigma2 = sigma2, estrisk = estrisk)
}

canonicalForm <- 
function(x, y, V) {
  # Remember to check whether V is squre matrix, x is a vector
  
  p <- length(x)
  N <- t(x) %*% x
  Nsr <- expm::sqrtm(N)
  VNsr <- v %*% Nsr 
  
  U <- t(solve(Nsr, t(x)))
  z <- t(U) %*% y
  S <- t(VNsr) %*% solve( V %*% N %*% V + diag(p) - V %*% V , VNsr)
  etahat <- U %*% S %*% t(U) %*% y
  
  list(U = U, z = z, S = S, etahat = etahat, p = p)
}

estRiskCanonical <- 
function(canonicalform,  estsig) {
  (l2sq(canonicalform$z - canonicalform$S %*% canonicalform$z) + 
     (2 * tr(canonicalform$S) - canonicalform$p) * estsig ) /canonicalform$p
}

estSigma <- 
function(x, y) {
  qx <- qr(x)
  coef <- solve.qr(qx, y)
  df <- nrow(x)-ncol(x)
  sigma2 <- sum((y - x%*%coef)^2)/df
  sigma2
}



