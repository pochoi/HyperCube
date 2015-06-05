
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
  class(proj) <- c("projection.hypercube", class(proj))
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
projectSub <-
function(proj, component) {  
  vname <- projectName(proj, component)
  subproj <- proj[vname]
  attr(subproj, "variables") <- attr(proj, "variables") 
  attr(subproj, "nlevels") <- attr(proj, "nlevels")
  attr(subproj, "dims") <- attr(proj, "dims")
  class(subproj) <- c("projection.hypercube", class(subproj))
  subproj
}

#' @export
projectWeight <- 
function(proj, component, weights = NULL) {
  p <- length(attr(proj, "variables"))
  dims <- attr(proj, "dims")  
  subproj <- projectSub(proj, component)
  if(is.null(weights)) weights <- rep(1, length(subproj))
  ans <- matrix(0, nrow=dims[1], ncol=dims[2])
  for(k in seq_along(subproj)) {
    ans <- ans + weights[k] * subproj[[k]]
  }
  attr(ans, "variable") <- attr(proj, "variables")
  attr(ans, "nlevels") <- attr(proj, "nlevels")
  attr(ans, "dims") <- attr(proj, "dims")
  attr(ans, "component") <- projectName(proj, component)
  attr(ans, "weights") <- weights
  ans
}

#' @export
projectFun <- 
function(proj, component = NULL) {
  if(is.null(component)) {
    subproj <- proj
  } else {
    subproj <- projectSub(proj, component)
  }
  
  f <- function(weight) {
    Reduce(`+`, mapply(`*`, subproj, weight, SIMPLIFY = FALSE))
  }
  attr(f, "variable") <- names(subproj)
  f
}
