library(HyperCube)
context("Basic tests")

test_that("Check hypercube", {
  proj <- projectCreate( ~ mother:infant -1, data = litter)
  component <- matrix(c(0,0,1,0,0,1,1,1),2,4)
  V <- projectWeight(proj, component = component, weights = c(1,1,1,1))
  hcmod <- hypercube(weight ~  mother:infant -1, data = litter, V)
  expect_is(hcmod, "hypercube")
})

test_that("Check hypercubeOptimization", {
  hcmod <- hypercubeOptimization(weight ~  mother:infant -1, data = litter)
  expect_is(hcmod$est, "hypercube")
})

test_that("Check projectCreate", {
  proj <- projectCreate( ~ mother:infant -1, data = litter)
  expect_is(proj, "projection.hypercube")
})

