
test_that("HMC condition factor", {
  Q <- readRDS('fit.RDS')$mle$Q
  cf1 <- .print.mat.stats(Q=Q)
  cf2 <- .print.mat.stats(Qinv=solve(Q))
  cf2a <- .print.mat.stats(Qinv=solve(as.matrix(Q)))
  expect_equal(cf1,cf2,cf2a)
  expect_error(.print.mat.stats(Qinv=solve(Q), Q=Q), 'Only one of Q or')
})
