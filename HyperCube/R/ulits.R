tr <- function(A) sum(diag(A))
l2sq <- function(v) sum(v^2)

#' Generate V matrix
#' 
#' @description covert W matrix to V matrix
#' 
#' @param W a matrix, penalized least square
#' 
#' @export
plsW2V <- 
function(W) {
  # check W is square matrix
  solve(expm::sqrtm(diag(dim(W)[1])+W))
}

#' Difference Matrix
#' 
#' @param p number of coefficients
#' @param dth order of difference matrix
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

