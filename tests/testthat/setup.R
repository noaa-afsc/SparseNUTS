
skip_consistency <- skip_reproducibility <- TRUE

# Only TMB is tested as of July 2025
skip_TMB <- FALSE
skip_RTMB <- FALSE

if(!skip_TMB){
  suppressWarnings(library(TMB))
  # setup simple object once, have to be careful b/c of the dynamic
  # links and I intentionally break the model downstream in tests
  # so that messes up obj. Hence the quick function to rebuild it
  # back to the MLE
  TMB::runExample('simple')
  obj0 <- obj
  par0 <- obj$env$parList()

  get_simple_obj <- function() {
    TMB::MakeADFun(data=obj0$env$data, parameters=par0, random=obj0$env$random,
                   dll=obj0$env$DLL, silent=TRUE)}

  obj <- get_simple_obj()


  # # dummy model for testing
  # pars <- list(x1=0,x2=0,x3=0)
  # f <- function(pars) {sum(-dnorm(c(pars$x1,pars$x2,pars$x3),0, 1,log=TRUE))}
  # f(pars)
  # obj <- RTMB::MakeADFun(func = f, parameters=pars, random=c('x2', 'x3'))
  # fit <- sample_snuts(obj, num_samples=300, num_warmup=100, chains=2, cores=1)
  # pairs(fit)
  # saveRDS(object=fit, file='tests/testthat/fit_snuts.RDS')
}

