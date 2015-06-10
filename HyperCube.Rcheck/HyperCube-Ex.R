pkgname <- "HyperCube"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "HyperCube-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('HyperCube')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("HyperCube")
### * HyperCube

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hypercube
### Title: Hypercube Estimator Fits
### Aliases: hypercube hypercube.default hypercube.formula

### ** Examples

## Example 1 in Beran (2014)
## Fitting Canadian earning data with Hypercube Estimator

# The number of age, p, in Example 1 in Beran (2014).
p <- length(unique(canadian.earnings[,1]))
# D_5 as in equation (3.10) in Beran (2014)
D <- diffMatrix(p, 5)
# The parametor nu in equation (3.11) in Beran (2014)
nu <- 100
# The matrix W in equation (3.11) in Beran (2014)
W <- nu * t(D) %*% D
# Convert W to V, as described in (1.6) in Beran (2014)
V <- plsW2V(W)
# The variable age should be considered as a factor
canadian.earnings[,"age"] <- factor(canadian.earnings[,"age"])
# Hyperpercube Estimator Fit
hcmod <- hypercube( log.income ~ age -1, data=canadian.earnings, V)

# Plot of data
plot(as.numeric(as.character(canadian.earnings$age)),
     canadian.earnings$log.income,
     xlab = "age", ylab = "log(income)")
# Plot of fitted line
lines(levels(canadian.earnings$age), hcmod$coefficients)


## Example 2 in Beran (2014)
## Fitting rat litter data

# Projection matrices as decribed in equation (5.4) in Beran (2014)
litter.proj <- projectCreate( ~ mother:infant -1, data = litter)
# If only additive effect is consider,
# take V = P1 + P2 + P3 (notation in equation (5.4) in Beran (2014))
component <- cbind(c(0,0), c(1,0), c(0,1))
V <- projectWeight(litter.proj, component = component)
# Hypercube Estimator Fit
hcmod <- hypercube( weight ~ mother:infant -1, data = litter, V)

# Estimated Risk
summary(hcmod)

## Hypercube Estimator with optimal risk
##
hcmodopt <- hypercubeOptimization( weight ~ mother:infant -1,
                                   data = litter
                                   )
# The optimal projection coefficient which minimizes the risk.
hcmodopt$projcoef

# The minimum risk
hcmodopt$estrisk

# The Hypercube Estimator fit with the V of the optimal projection.
summary(hcmodopt$$est)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("HyperCube", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("diffMatrix")
### * diffMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: diffMatrix
### Title: Difference Matrix
### Aliases: diffMatrix

### ** Examples

p <- 10
D <- diffMatrix(p, 5)
D



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("diffMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hypercubeEst")
### * hypercubeEst

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hypercubeEst
### Title: Estimation on linear Mode for Hypercube
### Aliases: hypercubeEst

### ** Examples

## The use of the function \code{hypercubeEst}.
mf <- model.frame(weight ~ mother:infant -1, data = litter)
X <- model.matrix(attr(mf, "terms"), mf)
y <- model.response(mf)
litter.proj <- projectCreate( ~ mother:infant -1, data = litter)
V <- projectWeight(litter.proj, weights = c(1,1,1,0))
hcmod <- hypercubeEst(X, y, V)
hcmod$coefficients



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hypercubeEst", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hypercubeOptimization")
### * hypercubeOptimization

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hypercubeOptimization
### Title: Hypercube Optimization
### Aliases: hypercubeOptimization

### ** Examples

hcmodopt <- hypercubeOptimization( weight ~ mother:infant -1, data = litter)
hcmodopt$projcoef #projction coefficients
hcmodopt$estrisk #estimated risk



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hypercubeOptimization", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plsW2V")
### * plsW2V

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plsW2V
### Title: Covert W matrix to V matrix
### Aliases: plsW2V

### ** Examples

# D_5 as in equation (3.10) in Beran (2014)
p <- 45
D <- diffMatrix(p, 5)
# The matrix W in equation (3.11) in Beran (2014)
W <- t(D) %*% D
# Convert W to V, as described in (1.6) in Beran (2014)
V <- plsW2V(W)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plsW2V", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("projectCreate")
### * projectCreate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: projectCreate
### Title: Functions for object "'projection.hypercube'"
### Aliases: projectCreate projectFun projectSub projectWeight

### ** Examples

proj <- projectCreate( ~ mother:infant -1, data = litter)

## proj contains the projection matrices :
## proj[["mother0:infant0"]], proj[["mother1:infant0"]], ...
## To see all the names of the projection matrices:
attr(proj, "component")

## If only the additive effects are needed,
## i.e. "mother0:infant0", "mother1:infant0" and "mother0:infant1",
component <- cbind(c(0,0), c(1,0), c(0,1))
proj.sub <- projectSub(proj, component)
attr(proj.sub, "component")

## If we want
## P = 1 * proj[["mother0:infant0"]] + 0.5 * proj[["mother1:infant0"]]
component <- cbind(c(0,0), c(1,0), c(0,1))
weights <- c(1, 0.5, 0)
proj.weights <- projectWeight(proj, component, weights)

## Create a function for more weighted projection matrices.
component <- cbind(c(0,0), c(1,0), c(0,1))
proj.fun <- projectFun(proj, component)

## Same as proj.weights in above example
weights <- c(1, 0.5, 0)
proj.fun(weights)

## Use the projection matrices for hypercube estimator
V <- proj.fun(weights)
hcmod <- hypercube(weight ~ mother:infant -1, data = litter, V = V)
summary(hcmod)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("projectCreate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
