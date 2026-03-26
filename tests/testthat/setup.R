
skip_consistency <- skip_reproducibility <- TRUE

# Only TMB is tested as of July 2025
skip_TMB <- FALSE
skip_RTMB <- FALSE

if(!skip_TMB){
  # suppressWarnings(library(TMB))
  # setup simple object once, have to be careful b/c of the dynamic
  # links and I intentionally break the model downstream in tests
  # so that messes up obj. Hence the quick function to rebuild it
  # back to the MLE
  TMB::runExample('simple')
  obj0 <- obj
  par0 <- obj$env$parList()

  get_tmb_obj <- function() {
    TMB::MakeADFun(data=obj0$env$data, parameters=par0, random=obj0$env$random,
                   DLL=obj0$env$DLL, silent=TRUE)}
}

if(!skip_RTMB){
  get_rtmb_obj <- function(){
    #library(RTMB)
    # dummy model for testing
    pars <- list(x1=0,x2=0,x3=0, xvec=c(1,1))
    f <- function(pars) {
      Sigma <- matrix(c(1,.8,.8,1), nrow=2)
      nll <- -RTMB::dmvnorm(pars$xvec, Sigma=Sigma, log=TRUE) +
        -sum(RTMB::dnorm(c(pars$x1,pars$x2,pars$x3),0, c(1,2,3),log=TRUE))
      return(nll)
    }
    obj <- RTMB::MakeADFun(func=f, parameters=pars, random=c('x2', 'x3'),
                           silent=TRUE)
    opt <- with(obj, nlminb(par,fn,gr))
    return(obj)
  }
#  obj <- get_rtmb_obj()
#  fit <- sample_snuts(obj, num_samples=400, num_warmup=100,
#                      chains=2, cores=1, init='random-t',
#                      seed=121414)
#  saveRDS(object=fit, file='tests/testthat/fit_snuts.RDS')

}
