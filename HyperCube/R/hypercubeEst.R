#' Hypercube Estimate
#' 
#' @param X design matrix
#' @param y observation
#' @param V sysmmetric matrix whose eigenvalues all lie in [0,1]
#' @param ... 
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @export
hypercubeEst <-
function(X, y, V, ...) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  #rX <- Matrix::rankMatrix(X)
  A <- hypercubeOp(X, V)
  #rA <- Matrix::rankMatrix(A)
  etahat <- A %*% y
  coef <- MASS::ginv(X) %*% etahat
#  sigma2 <- estSigma(mf)
#  estrisk <- estRisk(X, y, A, sigma2)
  residuals <- y - etahat
  
  list(coefficients = coef, fitted.values = etahat, 
       residuals = residuals
       #estsigma2 = sigma2, estrisk = estrisk
       )
}

#' Hypercube generic
#' 
#' @param ... Object to be hypercube
#'  @export
hypercube <- function(...) UseMethod("hypercube")

#' Hypercube Estimate
#' 
#' @param X design matrix
#' @param y observation
#' @param V sysmmetric matrix whose eigenvalues all lie in [0,1]
#' @param ... 
#' 
#' @describeIn hypercube 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
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

#' @describeIn hypercube 
#' 
#' @param formula formula to get estimate
#' @param data data you want to analysis
#' 
#' @export
#' @method hypercube formula
hypercube.formula <- function(formula, data, V, ...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  est <- hypercubeEst(X, y, V)
  est$mf
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

#' Hypercube Operator
#' 
#' @param X design matrix
#' @param V sysmmetric matrix whose eigenvalues all lie in [0,1]
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @export
hypercubeOp <- 
function(X, V) {
  # Check V is symmetric
  VXt <- V %*% t(X)
  t(VXt) %*% MASS::ginv( VXt %*% t(VXt) + diag(dim(X)[2]) - V %*% V) %*% VXt
}

#' Estimate Risk
#' 
#' @param X design matrix
#' @param y observation
#' @param A hypercuber operator
#' @param estsig estimated variance
#' 
#' @return The estimated risk
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @export
estRisk <-
function(X, y, A, estsig) {
  p <- dim(X)[2]
  n <- length(y)
  (l2sq(y - A %*% y) + ( 2 * tr(A) - n) * estsig)/p
}

#' Estimate Variance
#' 
#' @param mf model frame
#' 
#' @return The estimated variance
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @export
estSigma <- 
function(mf) {
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, data = mf)
  y <- model.response(mf)
  
  df <- nrow(X) - Matrix::rankMatrix(X)[1]

  if(df > 0) {
    eta <- X %*% MASS::ginv(X) %*% y
    sigma2 <- l2sq(y - eta)/df
    return(sigma2)
  } else {
    nvar <- length(all.vars(mt))-1
    if(nvar>1) {
      proj.formula <- reformulate(attr(mt,"term.labels"), 
                                  response = NULL, intercept = FALSE)
      proj <- projectMatrix(proj.formula, mf)
      component <- cbind(rep(0, nvar), diag(nvar))
      weights = rep(1, nvar+1) 
      P <- projectWeight(proj, component, weights) 
      XP <- X %*% P
      eta <- X %*% MASS::ginv(XP) %*% y
      rXP <- Matrix::rankMatrix(XP)[1]
      sigma2 <- l2sq(y - eta)/(nrow(X) - rXP)
      return(sigma2) 
    } else {
      stop("Please provide estimated sigma2.")
    }


  } 
  # Check if df = 0, then use submodel

  
  
}


#' Hypercube Optimization
#' 
#' @description give the component which minimizing the risk
#' 
#' @param formula formula
#' @param data data
#' @param sigma estimated variance
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @export
hypercubeOptimization <-
function(formula, data, sigma = NULL) {
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  
  mt <- attr(mf, "terms")
  proj.formula <- reformulate(attr(mt,"term.labels"), 
                              response = NULL, intercept = FALSE)
  
  proj <- projectMatrix(proj.formula, data)
  proj.fun <- projectFun(proj)
  nproj <- length(proj)

  if(is.null(sigma)) sigma <- estSigma(mf)

  projOptFun <- function(weights, estsigma) {
    proj <- proj.fun(weights)
    estRisk(X, y, hypercubeOp(X, proj), estsigma)
  }

  opt.ans <- optim(rep(0.5, nproj), 
                   projOptFun,
                   estsigma = sigma,
                   method = "L-BFGS-B", 
                   lower = rep(0, nproj), 
                   upper = rep(1, nproj)
                   )
  V <- proj.fun(opt.ans$par)
  est <- hypercube(formula, data, V)
  
  ans <- list( est = est, projcoef = opt.ans$par, risk = opt.ans$value)
  ans
}

#'@export
predict.hypercube <- 
function(object, newdata=NULL, ...) {
  if(is.null(newdata))
    y <- fitted(object)
  else{
    if(!is.null(object$formula)){
      ## model has been fitted using formula interface
      x <- model.matrix(object$formula, newdata)
    }
    else{
      x <- newdata
    }
    y <- as.vector(x %*% coef(object))
  }
  y
}


#print.hypercube
#residual.hypercube
#summary.hypercube



#formula.hypercube
#model.frame.hypercube
#model.matrix.hypercube