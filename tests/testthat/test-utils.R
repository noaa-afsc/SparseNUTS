
test_that("HMC condition factor", {
  Q <- readRDS('fit.RDS')$mle$Q
  cf1 <- .print.mat.stats(stats=list(Q=Q, max_cor=.1))
  cf2 <- .print.mat.stats(stats=list(Qinv=solve(Q), max_cor=.1))
  cf2a <- .print.mat.stats(stats=list(Qinv=solve(as.matrix(Q)), max_cor=.1))
  expect_equal(cf1,cf2,cf2a)
})

# rm(list=ls())
# devtools::load_all()
# TMB::runExample('simple')
# x <- .get_Q_stats(obj)
# Q <- x$Q
# y <- .get_Qinv_stats(obj)
#
# xx <- SparseNUTS:::.get_inputs(obj, FALSE, FALSE, 'dense', Q=NULL, Qinv=NULL)
# xx <- SparseNUTS:::.rotate_posterior()
#
# library(SparseNUTS)
# library(TMB)
# TMB::runExample('simple')
# mcmc1 <- sample_snuts(obj, metric='auto', chains=3, num_samples=500)
