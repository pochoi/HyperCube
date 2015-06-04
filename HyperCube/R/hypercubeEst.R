#' Hypercube Estimate
#' 
#' @export
hypercubeEst <-
function(X, y, V, ...) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  rX <- Matrix::rankMatrix(X)
  A <- hypercubeOp(X, V)
  rA <- Matrix::rankMatrix(A)
  etahat <- A %*% y
  coef <- MASS::ginv(X) %*% etahat
  sigma2 <- estSigma(X, y)
  estrisk <- estRisk(X, y, A, sigma2)
  residuals <- y - etahat
  
  list(coefficients = coef, fitted.values = etahat, residuals = residuals,
       estsigma2 = sigma2, estrisk = estrisk)
}

#' @export
hypercube <- function(x, ...) UseMethod("hypercube")

#' @export
#' @method hypercube default
hypercube.default <- function(X, y, V, ...)
{
  X <- as.matrix(X)
  y <- as.numeric(y)
  est <- hypercubeEst(X, y, V)
  est$call <- match.call()
  class(est) <- "hypercube"
  est
}

#' @export
#' @method hypercube formula
hypercube.formula <- function(formula, data, V, ...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  est <- hypercubeEst(X, y, V)
  est$call <- match.call()
  class(est) <- "hypercube"
  est
}




#' @export
#' @method print hypercube
print.hypercube <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

#' @export
hypercubeOp <- 
function(X, V) {
  # Check V is symmetric
  VXt <- V %*% t(X)
  t(VXt) %*% solve( VXt %*% t(VXt) + diag(dim(X)[2]) - V %*% V, VXt)
}

#' @export
estRisk <-
function(X, y, A, estsig) {
  p <- dim(X)[2]
  n <- length(y)
  (l2sq(y - A %*% y) + ( 2 * tr(A) - n) * estsig)/p
}

#' @export
estSigma <- 
function(x, y) {
  eta <- x %*% MASS::ginv(x) %*% y
  df <- nrow(x) - Matrix::rankMatrix(x)[1]
  # Check if df = 0, then use submodel
  sigma2 <- l2sq(y - eta)/df
  sigma2
}


#formula.hypercube
#model.frame.hypercube
#model.matrix.hypercube
#predict.hypercube
#print.hypercube
#residual.hypercube
#summary.hypercube
