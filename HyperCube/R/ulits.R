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

#' @export
projectElement <- function(n){
  u <- rep(1,n)/sqrt(n)
  J <- outer(u,u)         
  H <- diag(n) - J
  list(J=J,H=H)
}

#' @export
projectMatrix <-
function(formula, data) {
  m <- model.frame(formula, data)
  fname <- names(m)
  fn <- length(fname)
  fnlevel <- sapply(fname, function(k) length(levels(m[,k])))
  fperm <- matrix(, nrow = 2^fn , ncol = fn)
  flist <- list()
  for(k in fname) flist[[k]] <- c("J", "H")
  fperm <- as.matrix(expand.grid(flist, stringsAsFactors = FALSE))
  proj <- list()
  for(k in 1:dim(fperm)[1]) {
    v <- fperm[k,]
    vname <- paste0(fname, 
                    sapply(fperm[k,], function(k) switch(k, J = 0, H = 1)), 
                    collapse=":")
    projtemp <- matrix(1,1,1)
    for(j in fname) {
      projtemp <- kronecker(projectElement(fnlevel[j])[[ v[j] ]], projtemp)
    }
    proj[[vname]] <- projtemp
  }
  attr(proj, "variables") <- fname
  attr(proj, "nlevels") <- fnlevel
  attr(proj, "dims") <- dim(projtemp)
  class(proj) <- c("hypercube.proj", class(proj))
  proj
}

#' @export
projectName <- 
function(proj, component) {
  fname <- attr(proj, "variables")
  nname <- length(fname)
  if(length(component) %% nname != 0) 
    stop("Incorrect component argument.")
  
  if(is.null(dim(component))) {
    component <- matrix(component, nrow = nname)
  } else if (dim(component)[1] != nname) {
    stop("Incorrect component argument.")
  }
  vname <- apply(component, 2, function(v) paste0(fname, v, collapse=":"))
  vname
}

#' @export
projectGet <-
function(proj, v) {
  fname <- attr(proj, "variables")
  if(typeof(v) == "character") {
    vtemp <- array(0, length(fname))
    names(vtemp) <- fname
    vtemp[v[v %in% fname]] <- 1
    v <- vtemp
  }
  vname <- paste0(fname, v, collapse=":")
  proj[[vname]]
}




#' @export
projectModel <- 
function(proj, component, weigths = NULL) {
  p <- length(attr(proj, "variables"))
  dims <- attr(proj, "dims")
  
  if(length(component) %% p != 0) {
    stop("incorrect model.")
  } else if(is.null(dim(component))) {
    component <- matrix(m, ncol = p , byrow = TRUE)
  }
  
  if(is.null(weights)) weights <- rep(1, dim(component)[1])
  
  ans <- matrix(0, nrow=dims[1], ncol=dims[2])
  for(k in 1:dim(component)[1]) {
    ans <- ans + weights[k] * projectGet(proj, component[k,])
  }
  ans
}

#' @export
projectFun <- 
function(proj, component = NULL) {
  
  
  if(is.null(component)) {
    project.matrix <- proj
  } else {
    project.matrix <- list()
    for(k in 1:dim(component)[1]) {
      project.matrix[[k]] <- projectGet(proj, component[k,])
      vnames
    }
  }
  project.name <- names(project.matrix)
  
  f <- function(weight) {
    Reduce(`+`, mapply(`*`, project.matrix, weight, SIMPLIFY = FALSE))
  }
  attr(f, "variable") <- names(project.matrix)
  f
}
