litter

proj <- projectMatrix( ~ mother:infant -1, data = litter)
component <- matrix(c(0,0,1,0,0,1,1,1),2,4)
V <- projectWeight(proj, component = component, weights = c(1,1,1,1))
estRisk()

hcmod <- hypercube(weight ~  mother:infant -1, data = litter, V)
class(hcmod)
lmmod <- lm(weight ~  mother:infant -1, data = litter)
cbind(lmmod$coefficients, hcmod$coefficients)


W <- as.matrix(expand.grid(0:1,0:1,0:1,0:1))
r <- apply(W, 1,
           function(w) {
             V <- projectWeight(proj, component = component, weights = w)
             mf <- model.frame(weight ~  mother:infant -1, data = litter)
             mt <- attr(mf, "terms")
             X <- model.matrix(mt, mf)
             y <- model.response(mf)
             sigma2 <- estSigma(mf)
             A <- hypercubeOp(X, V)
             estRisk(X, y, A, sigma2)
           }
           )
cbind(W, r)


mf <- model.frame(weight ~  mother:infant -1, data = litter)
mt <- attr(mf, "terms")
X <- model.matrix(mt, mf)
y <- model.response(mf)
sigma2 <- estSigma(mf)
A <- hypercubeOp(X, V)
estRisk(X, y, A, sigma2)


hcopt <- hypercubeOptimization(weight ~  mother:infant -1, data = litter)

