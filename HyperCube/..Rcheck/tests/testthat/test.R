library(HyperCube)
context("Basic tests")

test_that("Least Square estimation", {
  proj <- projectMatrix( ~ mother:infant -1, data = litter)
  component <- matrix(c(0,0,1,0,0,1,1,1),2,4)
  V <- projectWeight(proj, component = component, weights = c(1,1,1,1))

  hcmod <- hypercube(weight ~  mother:infant -1, data = litter, V)
  lmmod <- lm(weight ~  mother:infant -1, data = litter)
  
  expect_is(hcmod, "hypercube")  
  expect_equal(as.(lmmod$coefficients), as.array(hcmod$coefficients))
})