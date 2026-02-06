
test_that("class methods working", {
  fit <- readRDS('fit.RDS')
  x <- as.data.frame(fit)
  expect_is(x, 'data.frame')
  summary(fit)
  print(fit)
})
