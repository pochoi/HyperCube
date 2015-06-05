tr <- function(A) sum(diag(A))
l2sq <- function(v) sum(v^2)

#' @export
plsW2V <- 
function(W) {
  # check W is square matrix
  solve(expm::sqrtm(diag(dim(W)[1])+W))
}

#' @export
subX2V <- 
function(X, L){
  X0 <- X %*% L
  XplusX0 <- MASS::ginv(X) %*% X0
  XplusX0 %*% MASS::ginv(XplusX0)
}

#' @export
modelMatrix <-
function(formula, data) {
  mf <- model.frame(formula, data)
  model.matrix(attr(mf, "terms"), mf)
}

#' @export
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
  zetahat <- S %*% z
  etahat <- U %*% zetahat
  
  list(U = U, z = z, S = S, zetahat = zetahat, etahat = etahat, p = p)
}

#' @export
estRiskCanonical <- 
function(canonicalform,  estsig) {
  (l2sq(canonicalform$z - canonicalform$S %*% canonicalform$z) + 
     (2 * tr(canonicalform$S) - canonicalform$p) * estsig ) /canonicalform$p
}

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

#' @export
polyRegMatrix <- 
function(deg, x) {
  sapply(0:(deg-1), function(d) x^d)
}
