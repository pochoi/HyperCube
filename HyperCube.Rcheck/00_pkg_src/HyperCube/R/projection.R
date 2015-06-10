# As described in the example 2 in article Beran (2014),
# The ANOVA projections are constructed from 
# the matrix "J" and "H" (same notation in the Beran (2014)).
# This function makes "J" and "H".
# The parameter "n" is the number of level of the factor
# of the variable in data set.
projectElement <- function(n){
  u <- rep(1,n)/sqrt(n)
  J <- outer(u,u)         
  H <- diag(n) - J
  list(J=J,H=H)
}

# Generate all the possible permutations for the projection matrices,
# as described in the equation (5.4) in Beran (2014).
# the parameter "variables": the names of all covariate in the data set.
# the parameter "JH.flag": If TRUE, represent the permutation by "J" and "H",
#                          If FALSE, J is 0, H is 1.
projectPerm <- function(variables, JH.flag = TRUE) {
  flist <- list()
  for(k in variables) flist[[k]] <- if(JH.flag) c("J", "H") else 0:1
  fperm <- t(as.matrix(expand.grid(flist, stringsAsFactors = FALSE)))
  fperm
}

# Contrust the names of the projection matrix, given the components
# as described in the equation (5.4) in Beran (2014).
# For example, in (5.4) in Beran (2014), 
# the projection matrix "P_1" is corresponding to the component "c(0,0)";
# the projection matrix "P_2" is corresponding to the component "c(1,0)";
# the projection matrix "P_3" is corresponding to the component "c(0,1)";
# the projection matrix "P_4" is corresponding to the component "c(1,1)";
#
# We can get all the names of all components at once.
# projectName(variables, c(0,0,1,0,0,1,1,1))
# Or, 
# component <- cbind(c(0,0), c(1,0), c(0,1), c(1,1))
# projectName(variables, component)
projectName <- 
function(variables, component) {
  fname <- variables
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


#' Functions for object "\code{projection.hypercube}"
#' 
#' The standard ANOVA projections as described in equation (5.4) in Beran (2014).
#' 
#' @details
#' The standard ANOVA projections as described in equation (5.4) in Beran (2014). The set of projections are sysmmetric, idempotent, mutually orthogonal matrix. The linear combinations of the projection matrices are used as the matrix V (equation (1.3) in Beran (2014)) of the hypercube estimators.
#' Please provide the formula in the form, for example in the data "litter", "~ mother:infant -1". No responese and No intercept. Only the interaction of factors.
#'
#' @param formula input formula without response variable. For example, "~ mother:infant -1".
#' @param data input data
#' @return The function \code{projectCreate} creates a object of the "\code{projection.hyercube}", which is a list containing all the projection matrices. Also, the object "\code{projection.hyercube}" contains the following attributes:
#' \describe{
#'   \item{variables}{The names of the variables in the formula}
#'   \item{nlevels}{The numbers of levels of the factor in each variable}
#'   \item{dims}{The dimension of the projection matrices}
#'   \item{component}{The names of the projection matrices stored in the objects}
#' }
#' 
#' @examples 
#' proj <- projectCreate( ~ mother:infant -1, data = litter)
#' 
#' ## proj contains the projection matrices :
#' ## proj[["mother0:infant0"]], proj[["mother1:infant0"]], ...
#' ## To see all the names of the projection matrices:
#' attr(proj, "component")
#' 
#' ## If only the additive effects are needed, 
#' ## i.e. "mother0:infant0", "mother1:infant0" and "mother0:infant1",
#' component <- cbind(c(0,0), c(1,0), c(0,1))
#' proj.sub <- projectSub(proj, component)
#' attr(proj.sub, "component")
#' 
#' ## If we want 
#' ## P = 1 * proj[["mother0:infant0"]] + 0.5 * proj[["mother1:infant0"]]
#' component <- cbind(c(0,0), c(1,0), c(0,1))
#' weights <- c(1, 0.5, 0)
#' proj.weights <- projectWeight(proj, component, weights)
#' 
#' ## Create a function for more weighted projection matrices.
#' component <- cbind(c(0,0), c(1,0), c(0,1))
#' proj.fun <- projectFun(proj, component)
#' 
#' ## Same as proj.weights in above example
#' weights <- c(1, 0.5, 0)
#' proj.fun(weights)
#' 
#' ## Use the projection matrices for hypercube estimator
#' V <- proj.fun(weights)
#' hcmod <- hypercube(weight ~ mother:infant -1, data = litter, V = V)
#' summary(hcmod)
#' 
#' 
#' @references Beran, Rudolf. "Hypercube estimators: Penalized least squares, submodel selection, and numerical stability." Computational Statistics & Data Analysis 71 (2014): 654-666.
#' @export
projectCreate <-
function(formula, data) {
  m <- model.frame(formula, data)
  fname <- names(m)
  fn <- length(fname)
  fnlevel <- sapply(fname, function(k) length(levels(m[,k])))
  fperm <- projectPerm(fname, JH.flag = TRUE)
  
  component <- array(,dim(fperm)[2])
  proj <- list()
  for(k in 1:dim(fperm)[2]) {
    v <- fperm[,k]
    vname <- paste0(fname, 
                    sapply(v, function(k) switch(k, J = 0, H = 1)), 
                    collapse=":")
    projtemp <- matrix(1,1,1)
    for(j in fname) {
      projtemp <- kronecker(projectElement(fnlevel[j])[[ v[j] ]], projtemp)
    }
    component[k] <- vname
    proj[[vname]] <- projtemp
  }
  attr(proj, "variables") <- fname
  attr(proj, "nlevels") <- fnlevel
  attr(proj, "dims") <- dim(projtemp)
  attr(proj, "component") <- component
  class(proj) <- c("projection.hypercube", class(proj))
  proj
}

#' @describeIn projectCreate
#' Subset of object "\code{projection.hypercube}"
#' 
#' @param proj an object of class "\code{projection.hypercube}".
#' @param component a vector or a matrix to specify with components in the \code{proj} are used. See examples for how to use.
#' @return The function \code{projectSub} returns a object of "\code{projection.hypercube}", which contains only the projection matrices according to the argument \code{component}.
#' 
#' @export
projectSub <-
function(proj, component) {  
  vname <- projectName(attr(proj, "variables"), component)
  subproj <- proj[vname]
  attr(subproj, "variables") <- attr(proj, "variables") 
  attr(subproj, "nlevels") <- attr(proj, "nlevels")
  attr(subproj, "dims") <- attr(proj, "dims")
  attr(subproj, "component") <- vname
  class(subproj) <- c("projection.hypercube", class(subproj))
  subproj
}

#' @describeIn projectCreate
#' Weighted sum of matrix in object "\code{projection.hypercube}"
#' 
#' @param weights The weights for the weighted sum of projection matrix. The order of the weights should be the same as the order of the component.
#' @return The function \code{projectWeight} returns a projection matrix which is the weigthed sum of the projection matrices according to \code{component} and \code{weights}.
#' 
#' @export
projectWeight <- 
function(proj, component = NULL, weights = NULL) {
  p <- length(attr(proj, "variables"))
  dims <- attr(proj, "dims")
  
  # If components is NULL, we want all the components in argument proj
  if(is.null(component)) {
    component <- projectPerm(attr(proj, "variables"), JH.flag=FALSE)
  }
  
  # Pojection matrices according to the argument component
  subproj <- projectSub(proj, component)
  vname <- attr(subproj, "component")
  
  # If weights is NULL, set all weights as 1
  if(is.null(weights)) weights <- rep(1, length(subproj))
  
  # The weighted sum of projection matrices
  ans <- matrix(0, nrow=dims[1], ncol=dims[2])
  for(k in seq_along(subproj)) {
    ans <- ans + weights[k] * subproj[[k]]
  }
  # Attr stores some information about the projection matrices
  attr(ans, "variable") <- attr(proj, "variables")
  attr(ans, "nlevels") <- attr(proj, "nlevels")
  attr(ans, "dims") <- attr(proj, "dims")
  attr(ans, "component") <- vname
  attr(ans, "weights") <- weights
  ans
}

#' @describeIn projectCreate
#' Functions of object "\code{projection.hypercube}"
#' 
#' @return The function \code{projectFun} returns a function which computes and returns a weighted sum of projection matrices according to \code{component} in the argument of \code{projectFun}. The argument of the returned function is the \code{weight}. See example for how to use.
#' 
#' @export
projectFun <- 
function(proj, component = NULL) {
  # If NULL, use all the component in argument proj 
  if(is.null(component)) {
    subproj <- proj
  } else {
    subproj <- projectSub(proj, component)
  }
  
  # the function we want for the weighted projection matrices
  f <- function(weight) {
    Reduce(`+`, mapply(`*`, subproj, weight, SIMPLIFY = FALSE))
  }
  
  # Attr stores some information about the projection matrices
  attr(f, "variable") <- attr(proj, "variables")
  attr(f, "nlevels") <- attr(proj, "nlevels")
  attr(f, "dims") <- attr(proj, "dims")
  attr(f, "component") <- names(subproj)
  f
}