## NOTE: This file was copied from adnuts on 2025-12-04


#' Function to generate random initial values from a previous fit using
#' SparseNUTS
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


## Convert SparseNUTS fit (named list) into a \code{shinystan} object.
##
## @details The shinystan packages provides several conversion functions
##   for objects of different types, such as stanfit classes (Stan ouput)
##   and simple arrays. For the latter, option NUTS information, such as
##   \code{sampler_params} can be passed. This function essentially extends
##   the functionality of \code{as.shinystan} to work specifically with
##   fits from SparseNUTS (TMB or RTMB). The user can thus explore their model
##   with \code{launch_shinystan(.as.shinytmbfit(fit))} in the same way
##   that Stan models are examined.
## @param fit Output list from  \code{sample_snuts}.
## @return An S4 object of class shinystan. Depending on the algorithm
##   used, this list will have slight differences.
.as.shinytmb <- function(fit){
  if(fit$algorithm=="NUTS"){
    sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup, max_treedepth=max_treedepth,
             sampler_params=sampler_params, algorithm='NUTS', model_name=model))
  }  else {
    stop("Currently only works for NUTS chains")
  }
  # else if(fit$algorithm=="HMC"){
  #       sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup,
  #            sampler_params=sampler_params, algorithm='HMC', model_name=model))
  # } else {
  #   sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup,
  #            algorithm='RWM', model_name=model))
  # }
  return(invisible(sso))
}

#' Launch shinystan for a TMB fit.
#' @export
#' @param fit A named list returned by \code{sample_snuts}.
launch_shinytmb <- function(fit){
  if(!requireNamespace('shinystan', quietly = TRUE))
    stop("The shinystan package is required for this functionality")
  shinystan::launch_shinystan(.as.shinytmb(fit))
}


#' Extract posterior samples from a model fit.
#'
#' A helper function to extract posterior samples across multiple chains
#' into a single data.frame.
#'
#' @details This function is loosely based on the \pkg{rstan} function
#'   \code{extract}. Merging samples across chains should only be used for
#'   inference after appropriate diagnostic checks. Do not calculate
#'   diagnostics like Rhat or effective sample size after using this
#'   function, instead, use \code{\link[posterior]{summarize_draws}}. Likewise, warmup
#'   samples are not valid and should never be used for inference, but may
#'   be useful in some cases for diagnosing issues.
#'
#' @param fit A list returned by \code{sample_snuts}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @param inc_lp Whether to include a column for the log posterior density
#'   (last column). For diagnostics it can be useful.
#' @param as.list Whether to return the samples as a list (one element per
#'   chain). This could then be converted to a CODA mcmc object.
#' @param unbounded Boolean flag whether to return samples in
#'   unbounded (untransformed) space. This can be useful for model debugging.
#' @return If as.list is FALSE, an invisible data.frame containing samples
#'   (rows) of each parameter (columns). If multiple chains exist they will
#'   be rbinded together, maintaining order within each chain. If as.list
#'   is TRUE, samples are returned as a list of matrices.
#' @export
#' @examples
#' ## A previously run fitted TMB model
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#' post <- extract_samples(fit)
#' tail(apply(post, 2, median))
extract_samples <- function(fit, inc_warmup=FALSE, inc_lp=FALSE,
                            as.list=FALSE, unbounded=FALSE){
  if(!is.tmbfit(fit)) stop("fit is not a valid object")
  if(unbounded){
    x <- fit$samples_unbounded
    if(is.null(x))
      stop("No unbounded parameters in this fit")
  } else {
    x <- fit$samples
    if(is.null(x)) stop("No posterior samples found")
  }
  if(!is.array(x)) stop("fit$samples is not an array -- valid fit object?")
  ind <- if(inc_warmup) 1:dim(x)[1] else -(1:fit$warmup)
  ## Drop LP
  if(inc_lp){
    y <-  lapply(1:dim(x)[2], function(i) x[ind, i,])
  } else {
    y <-  lapply(1:dim(x)[2], function(i) x[ind, i, -dim(x)[3]])
    # case with only a single column (parameter) need to recover array
    if(is.null(dim(y[[1]]))){
      y <- lapply(y, function(ytmp)
        array(ytmp, dim=c(length(ytmp),1), dimnames=list(NULL, dimnames(x)[[3]][1])))
    }
  }
  if(as.list){
    return(invisible(y))
  } else {
    # if(length(y))==1) return(invisible(as.data.frame(y[[1]])))
    return(invisible(as.data.frame(do.call(rbind, y))))
  }
}


#' Extract sampler parameters from a fit.
#'
#' Extract information about NUTS trajectories, such as acceptance ratio
#' and treedepth, from a fitted object.
#'
#' @details Each trajectory (iteration) in NUTS has associated information
#'   about the trajectory: stepsize, acceptance ratio, treedepth, and number of
#'   leapfrog steps. This function extracts these into a data.frame, which
#'   may be useful for diagnosing issues in certain cases. In general, the
#'   user should not need to examine them, or preferably should via
#'   \code{\link{plot_sampler_params}} or  \code{\link{launch_shinytmb}}.
#'
#' @param fit A list returned by \code{sample_snuts}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @return An invisible data.frame containing samples (rows) of each
#'   parameter (columns). If multiple chains exist they will be rbinded
#'   together.
#' @seealso \code{\link{launch_shinytmb}}.
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#' sp <- extract_sampler_params(fit, inc_warmup=TRUE)
#' str(sp)
#'
extract_sampler_params <- function(fit, inc_warmup=FALSE){
  x <- fit$sampler_params
  if(!is.list(x)) stop("fit$sampler_parameters is not a list -- valid fit object?")
  if(inc_warmup){
    ind <- 1:dim(x[[1]])[1]
    its <- 1:length(ind)
  } else{
    ind <- -(1:fit$warmup)
    its <- (1:nrow(x[[1]]))[ind]
  }
  y <- do.call(rbind, lapply(1:length(x), function(i){
    # drop=FALSE to prevent error when there is a single sample
    if(length(its) != NROW(x[[i]][ind,,drop=FALSE]))
      warning("Length mismatch in extract_sampler_params: iterations=", length(its),
              " and draws=", NROW(x[[i]][ind,,drop=FALSE]))
    cbind(chain=i, iteration=its, x[[i]][ind,,drop=FALSE])
  }))
  return(invisible(as.data.frame(y)))
}


#' Function to take a character vector of parameter names and
#' force them to be unique by appending numbers in square
#' brackets as needed
#' @param x Character vector
.make_unique_names <- function(x){
  stopifnot(is.character(x))
  as.vector((unlist(sapply(unique(x),
                           function(y){
  temp <- x[x==y]
  if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
}))))}


#' Check NUTS diagnostics of a fitted model
#' @param fit A fitted SNUTS object
#' @param print Whether to print the results to console
#' @details This is a wrapper which calls the internal function \code{\link[StanEstimators]{check_hmc_diagnostics}}
#' @return A data.frame containing diagnostic results (invisibly)
#' @export
check_snuts_diagnostics <- function(fit, print=TRUE){
  stopifnot(SparseNUTS::is.tmbfit(fit))
  sp <- fit$sampler_params |> posterior::as_draws_df()
  sp <- sp[sp$.iteration > fit$warmup,]
  StanEstimators:::check_hmc_diagnostics(draws_df = sp,
        max_treedepth = as.numeric(fit$max_treedepth), print=print)
}



#' Extract posterior samples from a tmbfit object
#' @param x A fitted tmbfit object
#' @param invf The inverse function to decorrelate the parameters
#' @param parnames A vector of parameter names, excluding lp__
#' @param array Whether to return a data.frame (default) or array
#'   which is used in constructing other objects downstream
#' @export
get_post <- function(x, invf, parnames, array=FALSE) {
  p <- x@draws |> as.data.frame()
  #q <- subset(p, select=-c(lp__, .iteration, .draw, .chain))
  q <- p[, which(!names(p) %in% c('lp__', '.iteration', '.draw', '.chain')), drop=FALSE]
  names(q) <- parnames
  if(ncol(q)==1){
    q <- as.data.frame(apply(q, 1, invf)) |> cbind(p$lp__)
  } else {
    q <- as.data.frame(t(apply(q, 1, invf))) |> cbind(p$lp__)
  }
  colnames(q) <- c(parnames, 'lp__')
  ## build array
  if(array){
    samples <- array(NA, dim=c(max(p$.iter), max(p$.chain), 1 + length(parnames)),
                     dimnames = list(NULL, NULL, c(parnames, "lp__")))
    for(chain in 1:max(p$.chain))
      samples[,chain,] <- as.matrix(q[p$.chain==chain,])
    return(samples)
  }
  return(q)
}


#' Print matrix stats
#'
#' @param stats A list of matrix stats
#' @param getCF A flag whether to calculate the HMC condition factor.
#' @return Prints the HMC condition factor, sparsity, and max absolute correlation
#' @details The HMC condition factor used is from Langmore et al.
#'   (2019). A NULL value is returned if the matrix passed is a
#'   scalar or if getCF=FALSE.
.print.mat.stats <- function(stats, getCF){
  if(!is.null(stats$Qinv)){
    x <- stats$Qinv
    if(NROW(x)==1) return(NULL) # not a matrix!
    nm <- 'Qinv'
    if(getCF) ev <- eigen(x,TRUE)$value
  } else {
    x <- stats[['Q']]
    if(NROW(x)==1) return(NULL) # not a matrix!
    nm <- 'Q'
    # need eigen values for covariance not precision so do 1/ev
    if(getCF) ev <- 1/eigen(x,TRUE)$value
  }
  maxcor <- stats$max_cor
  pct.sparsity <- round(100*mean(x[lower.tri(x)] == 0),2)
  # see Langmore et al. https://arxiv.org/abs/1905.09813
  if(getCF){
  factor <- sum((max(ev)/ev)^4)^(1/4)
    message(nm, " is ", pct.sparsity,
            "% sparse | HMC condition factor=",round(factor,0),
            ' | max cor >=', round(maxcor,4))
  return(invisible(factor))
  } else {
    factor <- max(stats$se)/min(stats$se)
    message(nm, " is ", pct.sparsity,
            "% sparse | Ratio of marginal SDs=",round(factor,1),
            ' | Max abs cor >=', round(maxcor,4))
    return(NULL)
  }
}

#' Get the joint precision matrix Q from an optimized TMB or RTMB obj, along with the standard errors and lower bound of the absolute maximum correlation. Done using efficient sparse methods, specifically the Takahashi approach implemented in the Takahashi_Davis function in the \pkg{sparseinverse} package.
#'
#' @param obj An optimized TMB or RTMB object
#' @param Q An optional sparse precision matrix if previous calculated. Otherwise calculated internally.
#' @param getCF A flag to determine whether to calculate the HMC condition factor. Set to FASLE by default due to poor performance in high dimensions even for a sparse Q.
#'
#' @return A list containing a sparse matrix Q, the standard errors derived from the square root of the diagonal of the inverse of Q, and the lower bound on the maximum correlation among parameters.
#'
.get_Q_stats <- function(obj=NULL, Q=NULL, getCF=FALSE){
  isRTMB <- ifelse(obj$env$DLL=='RTMB', TRUE, FALSE)
  if(length(obj$env$random)==0){
    warning("Q not available for models without random effects")
    return(NULL)
  }
  if(is.null(Q)){
    if(isRTMB){
      Q <- RTMB::sdreport(obj, getJointPrecision=TRUE,                                              skip.delta.method=TRUE)$jointPrecision
    } else {
      Q <- TMB::sdreport(obj, getJointPrecision=TRUE,
                           skip.delta.method=TRUE)$jointPrecision
    }
  }
  d <- nrow(Q)
  est <- obj$env$last.par.best
  stopifnot(d==length(est))
  inverse_subset <- sparseinv::Takahashi_Davis(Q=Q)
  ses <- sqrt(Matrix::diag(inverse_subset))
  rows <- inverse_subset@i + 1
  cols <- rep(1:d, diff(inverse_subset@p))
  vals <- inverse_subset@x
  # Filter out diagonals
  ind <- rows != cols
  cors <-
    vals[ind]/(ses[rows[ind]]* ses[cols[ind]])
  max_cor <- max(abs(cors))
  parnames <- .make_unique_names(names(est))
  dimnames(Q) <- list(parnames, parnames)
  names(est) <- names(ses) <- parnames
  stopifnot(length(est)==nrow(Q))
  stats <- list(Q=Q, ses=ses, est=est, max_cor=max_cor)
  stats$condition.factor <- .print.mat.stats(stats=stats, getCF=getCF)
  return(stats)
}



#' Get the covariance matrix Qinv from an optimized TMB or RTMB obj, along with the standard errors and lower bound of the absolute maximum correlation.
#' @param obj An optimized TMB or RTMB object
#' @param Qinv An optional covariance matrix if previous calculated. Otherwise calculated internally.
#' @param getCF A flag whether to calculate the HMC condition factor. Set to FALSE by default because it is slow in high dimensions.
#'
#' @return A list containing a covariance matrix Qinv, the standard errors derived from the square root of the diagonal of Qinv, and the lower bound on the maximum correlation among parameters.
#'
.get_Qinv_stats <- function(obj=NULL, Qinv=NULL, getCF=FALSE){
  isRTMB <- ifelse(obj$env$DLL=='RTMB', TRUE, FALSE)
  if(is.null(Qinv)){
    if(isRTMB){
      Qinv <- RTMB::sdreport(obj, skip.delta.method=TRUE)$cov.fixed
    } else {
      Qinv <- TMB::sdreport(obj, skip.delta.method=TRUE)$cov.fixed
    }
  }
  npar <- nrow(Qinv)
  est <- obj$env$last.par.best
  if(npar>1){
    ses <- sqrt(diag(Qinv))
    cor <- stats::cov2cor(Qinv)
    max_cor <- max(abs(cor[lower.tri(cor, diag=FALSE)]))
  } else {
    ses <- sqrt(Qinv)
    max_cor <- NA
  }
  parnames <- .make_unique_names(names(obj$par))
  dimnames(Qinv) <- list(parnames, parnames)
  names(est) <- names(ses) <- parnames
  stats <- list(Qinv=Qinv, ses=ses, est=est, max_cor=max_cor)
  stats$condition.factor <- .print.mat.stats(stats=stats, getCF=getCF)
  return(stats)
  }

#' Calculate gradient timings on a model for different metrics
#'
#' @param obj A TMB object
#' @param metrics A character vector of different metrics to benchmark
#' @param times How many evaluations to do
#' @param model_name An optional character name for the model, if
#'   NULL will pull from the DLL name
#' @export
#' @return A data.frame containing the median gradient timing
#'   (time), the percent sparsity of \eqn{Q} and the dimension of
#'   the model (npar).
#'
benchmark_metrics <- function(obj, times=1000, metrics=NULL,
                              model_name=NULL){
  if(!requireNamespace(package='microbenchmark', quietly=TRUE))
    stop("The microbenchmark package is required for this function")
  hasRE <- length(obj$env$random)>0
  if(is.null(metrics)){
    metrics <- c('unit', 'diag', 'dense', 'sparse')
    # if no RE, sparse won't work
    if(!hasRE)  metrics <- c('unit', 'diag', 'dense')
  }
  if(!hasRE & any(c('sparse', 'sparse-naive') %in% metrics))
    stop('sparse metric not allowed for models without random effects')
  if(is.null(model_name)) model_name <- obj$env$DLL
  message("optimizing model ", model_name, "...")
  obj$env$beSilent()
  opt <- with(obj, nlminb(par, fn, gr))
  if(hasRE){
    Q <- .get_Q_stats(obj)[['Q']]
    M <- solve(as.matrix(Q))
  } else {
    M <- .get_Qinv_stats(obj)[['Qinv']]
    Q <- solve(M)
  }
  n <- length(obj$env$last.par.best)
  res <- lapply(metrics, function(metric) {
    out <- sample_snuts(obj, rotation_only = TRUE,
                        metric=metric, Q=Q, Qinv=M,
                        skip_optimization = TRUE)
    metric <- out$metric # in case auto
    x0 <- out$x.cur[[1]]
    # make sure to add tiny random component during benchmarking to
    # avoid TMB tricks of skipping calcs
    time <- summary(microbenchmark::microbenchmark(out$gr2(x0+rnorm(n, sd=1e-10)),
                                                   unit='ms', times=times))$median
    return(data.frame(model=model_name, metric=metric, time=time))
  })
  res <- do.call(rbind, res) |>
    cbind(pct.sparsity=round(100*mean(Q[lower.tri(Q)] == 0),2)) |>
    cbind(npar=length(obj$env$last.par.best))
  res
}

## Simulate a single draw from either a normal (df=Inf) or t distribution (df<Inf) using the sparse precision Q if available, otherwise the dense covariance.
## @param inputs A list returned by .get_inits
mymvnorm <- function(inputs, df=Inf){
  Q <- inputs[['Q']]
  Qinv <- inputs$Qinv
  if(!is.null(Q)){
    # use efficient precision sampling
    if(class(Q)!='dsCMatrix')
      stop("This function only works for dsCMatrix objects")
    L <- Matrix::Cholesky(Q, super=TRUE, LDL=FALSE)
    u <- matrix(rnorm(ncol(L)), ncol(L))
    ## NOTE: This code requires LDL=FALSE
    u <- Matrix::solve(L, u, system="Lt") ## Solve Lt^-1 %*% u
    u <- Matrix::solve(L, u, system="Pt") ## Multiply Pt %*% u
    u <- as.numeric(u) # mean-0 white noise with covar=Q^-1
    if(is.infinite(df)) return(u)
    # construct t from u via inverse gamma relationship
    g <- rgamma(n=1, shape=df/2, rate=df/2)
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


