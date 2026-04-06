



# #' Prepare inputs for sparse sampling
# #'
# #'  @param obj TMB object
# #'  @param skip_optimization Whether to skip or not
# #'  @param laplace Whether to due the LA or not
# #'  @param metric Which metric
# #'  @param Q Sparse precision
# #'  @param Qinv Inverse of Q
# #'  @param skip.cor Whether to skip calculating the correlation matrix
# #'  @return A list containing Q, Qinv, the mle list, and timings
.get_inputs <-  function(obj, skip_optimization, laplace, metric, Q, Qinv,
                         skip_cor) {
  time.opt <- time.Q <- time.Qinv <- 0
  if(metric=='stan'){
    parnames <- .make_unique_names(names(obj$env$last.par.best))
    # mle <- list(nopar=length(obj$env$last.par.best),
    #             parnames=parnames)
    out <- list(time.opt=time.opt,
                time.Qinv=time.Qinv, time.Q=time.Q,
                parnames=parnames,
                laplace=laplace, metric=metric)
    return(out)
  }

  hasRE <-  !is.null(obj$env$random)
  if(laplace & !hasRE)
    stop("No random effects found so laplace=TRUE fails, set to FALSE")
  if( (laplace | !hasRE) & metric=='sparse')
    stop("sparse metric only allowed with random effects
           and laplace=FALSE")

  if(!skip_optimization){
    nfe <- length(obj$par)
    nre <- length(obj$env$last.par.best)-nfe
    message("Optimizing marginal posterior with ", nfe, " fixed effects and ", nre, " random effects...")
    if(any(as.numeric(obj$gr(obj$par))==0))
      warning('Some gradients were identically 0 at initial optimization values. Typically this indicates a misspecified model. Investigate model structure and retry.')
    time.opt <-
      as.numeric(system.time(opt <- tryCatch(with(obj, nlminb(par, fn, gr)), error=function(e) 'error'))[3])
    if(is.character(opt))
      stop("Optimization failed. Try optimizing externally (potentially with TMBhelper::fit_tmb) and then setting 'skip_optimization=TRUE'.  Or specify 'metric='stan'' to bypass optimization altogether if no marginal mode is expected.")
    if(opt$convergence !=0) warning("Optimization convergence code indicated failed convergence")
    maxgrad <- max(abs(obj$gr(opt$par)))
    opt$maxgrad <- maxgrad
    if(maxgrad>.1)
      warning("Maximum absolute marginal gradient is large (", sprintf('%.3e', maxgrad), ") which may indicate Q is unreliable")
  } else {
    opt <- NULL
  }


  if(hasRE & !laplace){
    ## Make parameter names unique if vectors exist
    message("Getting Q and its stats...")
    time.Q <- as.numeric(system.time(
      stats <- .get_Q_stats(obj=obj, Q=Q))[3]
    )
  } else {
    ## fixed effect only model so get Qinv (covariance)
    message("Getting Qinv and stats for fixed effects...")
    time.Qinv <-
      as.numeric(system.time(
        stats <- .get_Qinv_stats(obj=obj, Qinv=Qinv))[3]
      )
    # otherwise will be joint
    if(laplace) stats$est <- opt$par
  }
  est <- stats$est
  ses <- stats$ses
  parnames <- .make_unique_names(names(est))
  if(!all(is.finite(ses))){
    if(metric %in% c('unit', 'auto')){
      warning("Some standard errors estimated to be NaN, so downstream plotting may be affected.")
    } else {
      stop("Some standard errors estimated to be NaN, use 'unit' or 'stan' metric for models without a mode or positive definite Hessian")
    }
  }
  mle <- list(nopar=length(est), est=est, se=ses,
              Q=stats[['Q']], Qinv=stats[['Qinv']])
  if(is.null(skip_cor)) skip_cor <- ifelse(mle$nopar <= 2000, FALSE, TRUE)
  if(mle$nopar==1) skip_cor <- TRUE
  stopifnot(is.logical(skip_cor))
  if(!skip_cor){
    if(!is.null(stats[['Qinv']])) {
      mle$cor <- cov2cor(stats[['Qinv']])
    } else if(!is.null(stats[['Q']])) {
      mle$cor <- cov2cor(as.matrix(Matrix::solve(stats[['Q']])))
    } else {
      stop("skip_cor is TRUE but no Q or Qinv available")
    }
    # known exactly now so use it
    stats$max_cor <- max(abs(mle$cor[lower.tri(mle$cor)]))
  }
  out <- list(mle=mle, time.opt=time.opt,
              time.Qinv=time.Qinv, time.Q=time.Q, parnames=parnames,
              laplace=laplace, metric=metric, max_cor=stats$max_cor,
              condition.factor=stats$condition.factor, opt=opt)
  return(out)
}
