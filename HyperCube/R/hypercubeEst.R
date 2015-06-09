#' Estimation function for Hypercube Linear Model
#' 
#' This function calculates the hypercube estimates. Users should use \code{\link{hypercube}} which is a better user interface.
#' 
#' @param X design matrix
#' @param y observation
#' @param V sysmmetric matrix whose eigenvalues all lie in [0,1]. 
#' @param ... other optional arguments
#' @return A list containing the following:
#' \be{
#'   \item{coefficients}{Estimated coefficents \eqn{\hat{\beta}} under hypercube estimator}
#'   \item{fitted.values}{Fitted values (\eqn{\hat{\eta}})  under hypercube estimator}
#'   \item{residuals}{residuals, \eqn{y - \hat{\eta}}}
#'   \item{V}{Symmtric matrix with all eigenvalues lie in [0,1]. 
#'   The matrix V used the hypercube estimator.}
#' }
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @seealso \code{\link{hypercube}} is the user interfaces of hypercube estimators.
#'   
#' @examples
#' mf <- model.frame(weight ~ mother:infant -1, data = litter)
#' X <- model.matrix(attr(mf, "terms"), mf)
#' y <- model.response(mf)
#' litter.proj <- projectCreate( ~ mother:infant -1, data = litter)
#' V <- projectWeight(litter.proj, weights = c(1,1,1,0))
#' hcmod <- hypercubeEst(X, y, V)
#' hcmod$coefficients
#' 
#' @export
hypercubeEst <-
function(X, y, V, ...) {
  args <- list(...)
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  #rX <- Matrix::rankMatrix(X)
  A <- hypercubeOp(X, V)
  #rA <- Matrix::rankMatrix(A)
  etahat <- A %*% y
  coef <- MASS::ginv(X) %*% etahat
  rownames(coef) <- colnames(X)
  residuals <- y - etahat
  
  list(coefficients = coef, fitted.values = etahat, 
       residuals = residuals,
       V = V
       )
}

#' Hypercube Linear Model Fit
#' 
#' @note This package is still under development. For data with degree of freedom equal to zero, it may not be handled well. It is due the fact that it is not obvious to estimate the variance of the random error. Users need to provide estimated variance of random error in such cases.
#' 
#' @export
hypercube <- function(...) UseMethod("hypercube")

#' @describeIn hypercube Default method for \code{hypercube}
#' @param X design matrix
#' @param y observation
#' @param V sysmmetric matrix whose eigenvalues all lie in [0,1]
#' @param ...  other optional arguments
#' 
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @seealso \code{\link{hypercubeOptimization}}
#' @export
#' @method hypercube default
hypercube.default <- function(X, y, V, ...)
{
  X <- as.matrix(X)
  y <- as.numeric(y)
  est <- hypercubeEst(X, y, V)
  est$call <- match.call(expand.dots = TRUE)
  class(est) <- "hypercube"
  est
}

#' @describeIn hypercube Formula method for \code{hypercube}
#' 
#' @param formula formula to get estimate
#' @param data data you want to analysis
#' @return A object of the class "\code{hypercube}" containing the following:
#' \describe{
#'   \item{coefficients}{Estimated coefficents \eqn{\hat{\beta}} under hypercube estimator}
#'   \item{fitted.values}{Fitted values (\eqn{\hat{\eta}})  under hypercube estimator}
#'   \item{residuals}{residuals, \eqn{y - \hat{\eta}}}
#'   \item{estsigma2}{Estimated variance of \eqn{\epsilon}. 
#'    If the rank of the design matrix X is greater than the number of observation,
#'    the variance is estimated using least square regression. 
#'    If the rank is equal to the number of observation, 
#'    and if the model has more than one factor, 
#'    the variance is estimated using additive submodel fit.
#'    If neither of the above cases, 
#'    users need to provide their own estimated variance.}
#'   \item{estrisk}{Esitmated risk, which involved the estimated variance.}
#'   \item{V}{Symmtric matrix with all eigenvalues lie in [0,1]. 
#'   The matrix V used the hypercube estimator. }
#'   \item{modelframe}{The modelframe constructed from the argument formula and data.}
#' }
#' @export
#' @method hypercube formula
hypercube.formula <- function(formula, data, V, ...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  A <- hypercubeOp(X, V)
  est <- hypercubeEst(X, y, V)
  
  if("sigma2" %in% names(args)) {
    est$estsigma2 <- args[["sigma2"]]
  } else {
    est$estsigma2 <- estSigma(mf)$sigma
  }
  est$estrisk <- estRisk(X, y, A, est$estsigma2)
  est$modelframe <- mf
  est$formula <- formula
  est$call <- match.call(expand.dots = TRUE)
  class(est) <- "hypercube"
  est
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
    # If df > 0, then use least square fit
    eta <- X %*% MASS::ginv(X) %*% y
    sigma2 <- l2sq(y - eta)/df
    return(list(sigma2 = sigma2, eta = eta))
  } else {
    # If df = 0, and more then one factor, then use submodel fit
    nvar <- length(all.vars(mt))-1
    if(nvar>1) {
      proj.formula <- reformulate(attr(mt,"term.labels"), 
                                  response = NULL, intercept = FALSE)
      proj <- projectCreate(proj.formula, mf)
      component <- cbind(rep(0, nvar), diag(nvar))
      weights = rep(1, nvar+1) 
      P <- projectWeight(proj, component, weights) 
      XP <- X %*% P
      eta <- X %*% MASS::ginv(XP) %*% y
      rXP <- Matrix::rankMatrix(XP)[1]
      sigma2 <- l2sq(y - eta)/(nrow(X) - rXP)
      return(list(sigma2 = sigma2, eta = eta))
    } else {
      # At this moment, we cannot handle the case that df = 0 and only one factor.
      # Users need to provide estimated variance.
      stop("Please provide estimated sigma2.")
    }
  } 
}


#' Hypercube Optimization
#' 
#' @description Find the projection coefficients which minimizing the esitmated risk
#' 
#' @param formula formula
#' @param data data
#' @param sigma estimated variance
#' @return A list containing the following:
#' \describe{
#'   \item{est}{a object of the classs "\code{hypercube}".}
#'   \item{projcoef}{optimal coefficients for the projections.}
#'   \item{estrisk}{the minimum estimated risk.}
#' }
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @seealso \code{\link{hypercube}}
#' @examples
#' hcmodopt <- hypercubeOptimization( weight ~ mother:infant -1, data = litter)
#' hcmodopt$projcoef #projction coefficients
#' hcmodopt$estrisk #estimated risk
#' 
#' @export
hypercubeOptimization <-
function(formula, data, sigma = NULL) {
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  
  mt <- attr(mf, "terms")
  proj.formula <- reformulate(attr(mt,"term.labels"), 
                              response = NULL, intercept = FALSE)
  
  proj <- projectCreate(proj.formula, data)
  proj.fun <- projectFun(proj)
  nproj <- length(proj)

  if(is.null(sigma)) sigma <- estSigma(mf)$sigma2

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
  
  ans <- list( est = est, projcoef = opt.ans$par, estrisk = opt.ans$value)
  ans
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


#' Predicted values based on the HyperCube Linear Model Fits of the class "\code{hypecube}"
#' 
#' @param object Object of the class "\code{hypercube}"
#' @param newdata
#' 
#' @export
#' @method predict hypercube
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
#' Residuals of the class "\code{hypecube}"
#' 
#' @export
#' @method residuals hypercube
residuals.hypercube <- function(object, ...) {
  object$residuals
}

#' Fitted values of the class "\code{hypecube}".
#' @export
#' @method fitted hypercube
fitted.hypercube <- function(object, ...) {
  object$fitted.values
}

#' Fitted values of the class "\code{hypecube}"
#' 
#' @export
#' @method coefficients hypercube
coefficients.hypercube <- function(object, ...) {
  object$coefficients
}

#' Accessing Hypercube Estimator Fits
#' 
#' @export
#' @method summary hypercube
summary.hypercube <- function(object, ...) {
  X <- model.matrix(attr(object$modelframe, "terms"), object$modelframe)
  y <- model.response(object$modelframe)
  n <- nrow(X)
  r <- Matrix::rankMatrix(X)[1]
  df <- n - r
  
  if(df > 0) {
    estrisk <- object$estrisk
    fullestrisk <- estRisk(X, y, hypercubeOp(X, diag(ncol(X))), object$estsigma2)
  } else {
    estrisk <- object$estrisk
    fullestrisk <- NA
  }

  ans <- list(call = object$call,
              modelframe = object$modelframe,
              coefficients = object$coefficients,
              fitted.values = object$fitted.values,
              estrisk = estrisk,
              fullestrisk = fullestrisk
              )
  class(ans) <- "summary.hypercube"
  ans
}

#' @export
#' @method print summary.hypercube
print.summary.hypercube <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(paste("The estimated risk of hypercube estimation:", x$estrisk))
  cat("\n")
  if(!is.na(x$fullestrisk)) {
    cat(paste("The estimated risk of least square estimation: ", x$fullestrisk))
    cat("\n")
  }
}

