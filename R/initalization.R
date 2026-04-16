

#' Get a single initial value vector in untransformed model space
#' @param init The initial value strategy
#' @param obj2 The joint TMB model
#' @param inputs A list as returned by \code{.get_inputs}.
.get_inits <- function(init, obj2, inputs) {
  # only certain combinations of metrics and inputs can work
  metric <- inputs$metric
  if(init=='auto'){
    if(!is.null(inputs$mle[['Q']]) | !is.null(inputs$mle[['Qinv']])){
      init <- 'random-t'
    } else {
      init <- 'last.par.best'
    }
  }
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
  stopifnot(is.tmbfit(fit))
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
    u <- rmvnorm_Q(Q, nsim=1)
    u <- as.numeric(u)
    if(is.infinite(df)) return(u)
    # construct t from u via inverse gamma relationship
    g <- stats::rgamma(n=1, shape=df/2, rate=df/2)
    u/sqrt(g)
  } else if(!is.null(Qinv)) {
    if(is.infinite(df)){
      u <- mvtnorm::rmvnorm(n=1, sigma=Qinv)
    } else {
      u <- mvtnorm::rmvt(n=1, sigma=Qinv, df=df)
    }
    return(as.numeric(u))
  } else {
    stop("Neither Q nor Qinv available to simulate")
  }
}

#' Simulate draws from a MVN using the precision matrix Q
#' @param Q A sparse precision matrix
#' @param nsim How many replicates to simulate
rmvnorm_Q <- function(Q, nsim){
  if(!is(Q, 'dsCMatrix'))
    stop("This function only works for dsCMatrix objects")
  L <- Matrix::Cholesky(Q, super=TRUE, LDL=FALSE)
  u <- matrix(rnorm(nsim*ncol(L)), ncol(L))
  u <- Matrix::solve(L, u, system="Lt") ## Solve Lt^-1 %*% u
  u <- Matrix::solve(L, u, system="Pt") ## Multiply Pt %*% u
  u <- as.matrix(u) # mean-0 white noise with covar=Q^-1
  return(u)
}

#' Generate approximate posterior samples from the sparse
#' precision matrix Q, assuming multivariate normality.
#'
#' @param fit An object returned by \code{sample_snuts}
#' @param nsamples The number of samples to draw, defaulting to 4000.
#' @param seed An optional random seed
#' @return A data.frame of approximate samples
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#' x <- sample_Q(fit)
sample_Q <- function(fit, nsamples=4000, seed=NULL){
  stopifnot(is.tmbfit(fit))
  if(!is.null(seed)) set.seed(seed)
  Q <- fit$mle[['Q']]
  if(is.null(Q)){
    if(!is.null(fit$mle$Qinv)){
      out <- as.data.frame(mvtnorm::rmvnorm(n=nsamples, mean= fit$mle$est, sigma=fit$mle$Qinv))
    } else {
      stop("The fitted object does not contain Q or Qinv. Try rerunning with optimization turned on or specify Q directly.")
    }
  } else {
    if(is.null(fit$mle$est)) stop("No conditional mean found in fit$mle$est so cannot sample")
    out <- rmvnorm_Q(fit$mle$Q, nsim = nsamples)
    out <- as.data.frame(t(out + fit$mle$est))
    return(out)
  }
}
