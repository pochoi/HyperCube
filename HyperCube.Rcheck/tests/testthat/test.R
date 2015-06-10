library(HyperCube)
context("Basic tests")

test_that("Check projectCreate", {
  proj <- projectCreate( ~ mother:infant -1, data = litter)
  expect_is(proj, "projection.hypercube")
})

test_that("Check hypercube, litter", {
  proj <- projectCreate( ~ mother:infant -1, data = litter)
  component <- matrix(c(0,0,1,0,0,1,1,1),2,4)
  V <- projectWeight(proj, component = component, weights = c(1,1,1,0))
  hcmod <- hypercube(weight ~  mother:infant -1, data = litter, V)
  expect_is(hcmod, "hypercube")
})

test_that("Check hypercube and lm, litter", {
  proj <- projectCreate( ~ mother:infant -1, data = litter)
  V <- projectWeight(proj, weights = c(1,1,1,1))
  hcmod <- hypercube(weight ~  mother:infant -1, data = litter, V)
  lmmod <- lm(weight ~  mother:infant -1, data = litter)
  expect_equal(object = as.numeric(hcmod$coefficients), 
               expected = as.numeric(lmmod$coefficients)
               )
})

test_that("Check hypercube, canadian.earnings", {
  p <- length(unique(canadian.earnings[,1]))
  D <- diffMatrix(p, 5)
  nu <- 100
  W <- nu * t(D) %*% D
  V <- plsW2V(W)
  canadian.earnings[,"age"] <- factor(canadian.earnings[,"age"])
  hcmod <- hypercube( log.income ~ age -1, data=canadian.earnings, V)
  expect_is(hcmod, "hypercube")
})

test_that("Check hypercubeOptimization", {
  hcmodopt <- hypercubeOptimization(weight ~  mother:infant -1, data = litter)
  expect_is(hcmodopt$est, "hypercube")
  expect_equal(object = hcmodopt$projcoef, 
               expected = c(0.997, 0.693, 0.000, 0.415), 
               tolerance = .001)
})

test_that("Check diffMatrix", {
  D <- diffMatrix(5,3)
  Dtest <- rbind(c(1,-3,3,-1,0), c(0,1,-3,3,-1))
  expect_equal(object = D, 
               expected = Dtest)
})