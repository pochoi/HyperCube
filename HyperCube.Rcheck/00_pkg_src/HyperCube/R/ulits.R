# Trace of matrix
tr <- function(A) sum(diag(A))

# L2 norm of vector
l2sq <- function(v) sum(v^2)

#' Covert W matrix to V matrix
#' 
#' @description Convert W to V, as described in (1.6) in Beran (2014)
#' 
#' @param W a matrix, penalized least square
#' @return The function \code{plsW2V} returns a matrix V for the Hypercube Estimator.
#' 
#' @examples
#' # D_5 as in equation (3.10) in Beran (2014)
#' p <- 45
#' D <- diffMatrix(p, 5)
#' # The matrix W in equation (3.11) in Beran (2014)
#' W <- t(D) %*% D
#' # Convert W to V, as described in (1.6) in Beran (2014)
#' V <- plsW2V(W)
#' 
#' @export
plsW2V <- 
function(W) {
  # check W is square matrix
  solve(expm::sqrtm(diag(dim(W)[1])+W))
}

#' Difference Matrix
#' 
#' @description The difference matrix described in equation (3.10) in Beran (2014).
#' 
#' @param p number of coefficients
#' @param dth order of difference matrix
#' 
#' @examples
#' p <- 10
#' D <- diffMatrix(p, 5)
#' D
#' 
#' @export
diffMatrix <- 
function(p, dth) {
  D <- diag(p)
  dif <- function(p) cbind(diag(p-1),0) + cbind(0, -diag(p-1))
  for(k in seq(p, p-(dth-1),-1)){
    D <- dif(k) %*% D
  }
  D
}

