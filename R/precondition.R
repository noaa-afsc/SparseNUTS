## NOTE: This file was moved over from the adnuts repo on
## 2025-12-04 and the git history can be found there.
## https://github.com/Cole-Monnahan-NOAA/adnuts/commit/33744b850b8ee0642b8c8095181ebb4b878e1495



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




#' Update algorithm for mass matrix.
#'
#' @param metric The metric to use
#' @param fn The current fn function.
#' @param gr The current gr function
#' @param inputs A list of inputs
#' @param y.cur The current parameter vector in unrotated (Y) space.
.rotate_posterior <- function(metric, fn, gr, inputs, y.cur){
  ## Rotation done using choleski decomposition
  Q <- inputs$mle[['Q']]
  Qinv <- inputs$mle[['Qinv']]
  if(metric=='dense'){
    if(is.null(Qinv)){
      if(is.null(Q)) stop("Neither Q nor Qinv available for dense preconditioner")
      Qinv <- solve(Q)
    }
    ## First case is a dense mass matrix
    M <- as.matrix(Qinv)
    # took this out b/c it was warning too often, better way to test?
    #if(!matrixcalc::is.symmetric.matrix(M) ||
    #  !matrixcalc::is.positive.definite(M))
    # warning("Estimated dense matrix was not positive definite so may be unreliable. Try different metric or turn on the laplace if there are random effects if it fails.")
    J <- NULL
    chd <- t(chol(M))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    ## Define rotated fn and gr functions
    fn2 <- function(x) fn(chd %*% x)
    gr2 <- function(x) {as.vector( gr(chd %*% x) %*% chd )}
    ## Now rotate back to "x" space using the new mass matrix M
    x.cur <- lapply(y.cur, \(y) as.numeric(chd.inv %*% y))
    finv <- function(x){
      t(chd %*% x)
    }
  } else if(metric=='diag'){
    ## diagonal but not unit
    if(is.null(inputs$mle$se)) stop("diag metric failed since SE unavailable")
    J <- NULL
    chd <- inputs$mle$se
    fn2 <- function(x) fn(chd * x)
    gr2 <- function(x) as.vector(gr(chd * x) ) * chd
    ## Now rotate back to "x" space using the new mass matrix M. M is a
    ## vector here. Note the big difference in efficiency without the
    ## matrix operations.
    x.cur <- lapply(y.cur, \(y) as.numeric((1/chd) * y))
    finv <- function(x) chd*x
  } else if(metric=='unit' | metric=='stan') {
    ## unit metric, change nothing
    fn2 <- function(x) fn(x)
    gr2 <- function(x) gr(x)
    x.cur <- y.cur
    finv <- function(x) x
    chd <- J <- NULL
  } else if(metric=='sparse-naive'){
    # This metric is carefully constructured to match the dense
    # metric up to numerical precision. But as it is slower it is
    # not typically used.
    # stopifnot(require(Matrix))
    if(!is(Q,"Matrix")) stop("Q is not a Matrix object, something went wrong")
    # M is actually Q, i.e., the inverse-mass
    # Antidiagonal matrix JJ = I
    J <- Matrix::sparseMatrix( i=1:nrow(Q), j=nrow(Q):1 )
    #chd <- Cholesky(M, super=FALSE, perm=FALSE)
    #chd <- Matrix::Cholesky(M, super=TRUE, perm=FALSE)
    chd <- Matrix::Cholesky(J%*%Q%*%J, super=TRUE, perm=FALSE) # perm
    Linv_times_x = function(chd,x){
      as.numeric(
        J%*% Matrix::solve(chd, Matrix::solve(chd, J%*%x, system="Lt"), system="Pt")
        )
    }
    x_times_Linv = function(chd,x){
      #x %*% chol()
      as.numeric(
        J%*%Matrix::solve(chd, Matrix::solve(chd, Matrix::t(x%*%J), system="L"),
                          system="Pt")
        )
    }
    fn2 <- function(x){
      Linv_x = Linv_times_x(chd, x)
      fn(Linv_x)
    }
    gr2 <- function(x){
      Linv_x = Linv_times_x(chd, x)
      grad = gr( Linv_x )
      as.vector(  x_times_Linv(chd, grad) )
    }
    ## Now rotate back to "x" space using the new mass matrix M
    #  solve(t(chol(solve(M)))) ~~ IS EQUAL TO ~~ J%*%chol(M)%*%J
    # J%*%chol(J%*%prec%*%J) %*% J%*%x
    x.cur <- lapply(y.cur, \(y) as.numeric(J%*%chol(J%*%Q%*%J) %*% J%*%y))
    finv <- function(x){
      t(as.numeric(
        J%*%Matrix::solve(chd, Matrix::solve(chd, J%*%x, system="Lt"), system="Pt")
      )
      )
    }
  } else if(metric=='sparse'){
    if(!is(Q,"Matrix")) stop("Q is not a Matrix object, something went wrong")
    # Do Cholesky on Q permuted directly
    J <- NULL
    chd <- Matrix::Cholesky(Q, super=TRUE, perm=TRUE)
    L <- as(chd, "sparseMatrix")
    perm <- chd@perm + 1L
    iperm <- Matrix::invPerm(perm)
    # Drop all numerical zeros and convert to triangular storage
    L <- Matrix::tril(Matrix::drop0(L)) ## class(L) == "dtCMatrix"
    Lt <- Matrix::t(L) ## class(Lt) == "dtCMatrix"
    x.cur <- lapply(y.cur, \(y) as.vector(Lt %*% y[perm]))
    fn2 <- function(x)  fn(Matrix::solve(Lt, x)[iperm])
    gr2 <- function(x){
      y <- Matrix::solve(Lt, x)[iperm]
      Matrix::solve(L, as.numeric(gr(y))[perm])
    }
    finv <- function(x)   as.numeric(Matrix::solve(Lt, x)[iperm])
  } else if(metric=='auto'){
    if(NROW(Qinv)==1 | NROW(Q)==1){
      message("diag metric selected b/c only 1 parameter")
      rdiag <-
        .rotate_posterior(metric='diag', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(rdiag)
    }

    if((is.null(Q) & is.null(Qinv)) |
       is.null(inputs$max_cor) |
       !is.finite(inputs$max_cor)){
      message("unit metric selected b/c Q and Qinv unavailable")
      runit <-
        .rotate_posterior(metric='unit', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(runit)
    }


    if(inputs$max_cor <= 0.2){
      message("diag metric selected b/c low max cor (", round(inputs$max_cor,4),")")
      rdiag <-
        .rotate_posterior(metric='diag', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(rdiag)
    }

    # high correlations exist
    # fixed effects only models or ELA
    if(is.null(Q) & !is.null(Qinv)){
      message("dense metric selected b/c no Q and high max cor (", round(inputs$max_cor,4),")")
      rdense <-
        .rotate_posterior(metric='dense', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(rdense)
    }

    # random effects models with sparse Q and high correlation
    if(!is.null(Q)){
      if(nrow(Q)>500){
        # high dimension, skip dense calcs since sparse almost always faster
        message("sparse metric selected b/c high dimesions and high max cor (", round(inputs$max_cor,4),")")
        rsparse <-
          .rotate_posterior(metric='sparse', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
        return(rsparse)
      } else {
        # low dimension sometimes dense is faster so check it
        if(!requireNamespace("microbenchmark", quietly=TRUE)){
          message("sparse metric selected b/c no timing available. Please install microbenchmark")
          rsparse <-
            .rotate_posterior(metric='sparse', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
          ## check for speed differences
          return(rsparse)
        } else {
          # when doing timing need to add random components to
          # input, otherwise TMB may skip calculations and throw
          # off benchmarking
          rsparse <-
            .rotate_posterior(metric='sparse', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
          if(is.null(Qinv)) inputs$Qinv <- as.matrix(Matrix::solve(Q))
          rdense <-
            .rotate_posterior(metric='dense', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
          npars <- length(rdense$x.cur[[1]])
          bench <- microbenchmark::microbenchmark(
            rdense$gr2(rdense$x.cur[[1]]+rnorm(npars, sd=1e-10)),
            rsparse$gr2(rsparse$x.cur[[1]]+rnorm(npars, sd=1e-10)),
            times = 500
          )
          tdense <- summary(bench)$median[1]
          tsparse <- summary(bench)$median[2]
          if(tdense < tsparse){
            message("dense metric selected b/c faster than sparse and high correlation (max=",
                    round(inputs$max_cor,4), ")")
            return(rdense)
          } else {
            message("sparse metric selected b/c faster than dense and high correlation (max=",
                    round(inputs$max_cor, 4), ")")
            return(rsparse)
          }
        }
      }
    }
  } else {
    stop("Invalid metric specified")
  }
  return(list(gr2=gr2, fn2=fn2, finv=finv, x.cur=x.cur, chd=chd, J=J, metric=metric))
}

