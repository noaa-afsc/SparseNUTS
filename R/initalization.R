

#' Get a single initial value vector in untransformed model space
#' @param init The initial value strategy
#' @param obj2 The joint TMB model
#' @param inputs A list as returned by \code{.get_inputs}.
.get_inits <- function(init, obj2, inputs) {
  # only certain combinations of metrics and inputs can work
  metric <- inputs$metric
  if(metric=='stan' & init %in% c('random', 'random-t'))
    stop("'stan' metric not allowed with 'random' or 'random-t' init b/c no Qinv")
  if((is.null(inputs$Qinv | is.null(inputs$mle$est))) &
     init %in% c('random', 'random-t'))
    stop("Cannot use random inits b/c mode or Qinv does not exist, use a different init")
  # this will be the mode if skip_optimization=FALSE in .get_inputs
  lpb <- obj2$env$last.par.best
  if(inputs$laplace){
    # joint vector needs to collapse down to just the fixed
    # effects when doing ELA
    lpb <- lpb[!obj2$env$lrandom()]
  }
  npars <- length(lpb)
  for(ii in 1:10){
    inits <-
      switch(init,
             'last.par.best' = lpb,
             'random-t'      = inputs$mle$est + mymvnorm(inputs, df=2),
             'random'        = inputs$mle$est + mymvnorm(inputs, df=Inf),
             'unif'          = runif(n=npars, min=-2, max=2))
    inits <- as.numeric(inits)
    if(length(inits)!=npars)
      stop("Wrong vector length for inits:", length(inits), " when should be", npars)
    # if(init=='last.par.best')
    #   inits <- obj2$env$last.par.best
    # if(init=='random-t')
    #   inits <- as.numeric(inputs$mle$est + mvtnorm::rmvt(n=1, sigma=inputs$Qinv, df=2))
    # if(init=='random')
    #   inits <- as.numeric(inputs$mle$est + mvtnorm::rmvnorm(n=1, sigma=inputs$Qinv))
    # if(init=='unif')
    #   inits <- as.numeric(runif(n=length(rotation$x.cur), min=-2, max=2))
    success <- is.finite(obj2$fn(inits)) & is.finite(sum(obj2$gr(inits)))
    if(success) break
  }
  if(!success)
    stop(init, " inits resulted in NaN log-posterior after 10 tries, try another method or investigate model")
  return(inits)
}



#' Function to generate random initial values from a previous fit
#'
#' @param fit A fitted object from \code{\link{sample_snuts}}
#' @param chains The number of chains for the subsequent run, which
#'   determines the number to return.
#' @return A list of vectors which can be passed back into
#'   \code{\link{sample_snuts}}.
#' @export
sample_inits <- function(fit, chains){
  post <- extract_samples(fit)
  ind <- sample(1:nrow(post), size=chains)
  lapply(ind, function(i) as.numeric(post[i,]))
}



## Simulate a single draw from either a normal (df=Inf) or t distribution (df<Inf) using the sparse precision Q if available, otherwise the dense covariance.
## @param inputs A list returned by .get_inits
mymvnorm <- function(inputs, df=Inf){
  Q <- inputs$mle[['Q']]
  Qinv <- inputs$mle$Qinv
  if(!is.null(Q)){
    # use efficient precision sampling
    if(!is(Q, 'dsCMatrix'))
      stop("This function only works for dsCMatrix objects")
    L <- Matrix::Cholesky(Q, super=TRUE, LDL=FALSE)
    u <- matrix(rnorm(ncol(L)), ncol(L))
    ## NOTE: This code requires LDL=FALSE
    u <- Matrix::solve(L, u, system="Lt") ## Solve Lt^-1 %*% u
    u <- Matrix::solve(L, u, system="Pt") ## Multiply Pt %*% u
    u <- as.numeric(u) # mean-0 white noise with covar=Q^-1
    if(is.infinite(df)) return(u)
    # construct t from u via inverse gamma relationship
    g <- stats::rgamma(n=1, shape=df/2, rate=df/2)
    u/sqrt(g)
  } else if(!is.null(Qinv)) {
    if(is.infinite(df)){
      u <- mvtnorm::rmvnorm(n=1, sigma=Qinv)
    } else {
      u <- mvtnorm::rmvt(n=1, sigma=inputs$Qinv, df=df)
    }
    return(as.numeric(u))
  } else {
    stop("Neither Q nor Qinv available to simulate")
  }
}
