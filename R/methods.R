


#' Constructor for the "tmbfit"  class
#' @param x Fitted object from \code{\link{sample_snuts}}
#' @return An object of class "tmbfit"
#' @export
tmbfit <- function(x){
  stopifnot(is.list(x))
  if(is.null(x$samples)) stop("Samples missing from fit")
  if(is.null(x$algorithm)) stop("Algorithm missing from fit")
  class(x) <- c('tmbfit', 'list')
  x
}

#' Check object of class tmbfit
#' @param x Returned list from \code{\link{sample_snuts}}
#' @export
is.tmbfit <- function(x) inherits(x, "tmbfit")


#' Construtor for tmbfit objects
#' @param x A fitted MCMC object
#' @param parnames A character vector of unique par names
#' @param mle A list of MLE parameters
#' @param invf The inverse function for the parameters
#' @param metric The metric used
#' @param skip_monitor Whether to skip monitor calculations
#' @param model A character giving the model name
#' @export
as.tmbfit <- function(x, parnames, mle, invf, metric,
                      skip_monitor,  model='anonymous'){
  post <- get_post(x, invf, parnames=parnames, array=TRUE)
  sp <- as.data.frame(x@diagnostics)
  spl <- list()
  for(chain in 1:max(sp$.chain)){
    spl[[chain]] <- as.matrix(sp[sp$.chain==chain,1:6])
  }
  timing <- sapply(x@timing, function(x) unlist(x))
  timing <- list(time.warmup=timing[1,],time.sampling=timing[2,],
                 time.total=timing[1,]+timing[2,])
  thin <- as.numeric(x@metadata$thin)
  warmup <- ceiling(as.numeric(x@metadata$num_warmup)/thin)
  iter <- ceiling(as.numeric(x@metadata$num_samples)/thin)
  if(dim(post)[1] != warmup + iter){
    stop("Error in output dimensions: iter=",iter, "; warmup=", warmup,
         "; nrow(post)=", nrow(post))
  }
  # make sure to call monitor *after* back-transforming the
  # parameters into the original, correlated space
  # mon <- cbind(variable=c(parnames, 'lp__'),
  #              rstan::monitor(post, warmup = warmup, print=FALSE))
  if(skip_monitor){
    mon <- NULL
  } else {
    mon <- posterior::summarise_draws(post[-(1:warmup),,])
  }
  x <- list(samples=post, sampler_params=spl, mle=mle,
            monitor=mon, model=model,
            metric=metric,
            par_names=parnames,
            max_treedepth=x@metadata$max_depth,
            warmup=warmup, iter=iter, thin=thin,
            timing=timing,
            algorithm='SNUTS')
  tmbfit(x)
}


#' Convert object of class tmbfit to data.frame. Calls
#' \code{\link{extract_samples}}
#'
#' @param x Fitted object from \code{\link{sample_snuts}}
#' @param row.names Ignored
#' @param optional Ignored
#' @param ... Ignored
#' @return A data frame with parameters as columns and samples as
#'   rows.
#' @details This calls the default settings of
#'   \code{\link{extract_samples}}, no warmup samples and no
#'   column for the log-posterior (lp__). Use this function
#'   directly for finer control.
#' @export
as.data.frame.tmbfit <-
  function(x, row.names=NULL, optional=FALSE, ...)
    extract_samples(x)

#' Plot object of class tmbfit
#' @param x Fitted object from \code{\link{sample_snuts}}
#' @param y Ignored
#' @param ... Ignored
#' @return Plot created
#' @method plot tmbfit
#' @export
plot.tmbfit <- function(x, y, ...) plot_marginals(x)

#' Print summary of object of class tmbfit
#' @param object Fitted object from \code{\link{sample_snuts}}
#' @param ... Ignored
#' @return Summary printed to screen
#' @method summary tmbfit
#' @export
summary.tmbfit <- function(object, ...) print(object)

#' Print summary of tmbfit object
#' @param x Fitted object from \code{\link{sample_snuts}}
#' @param ... Ignored
#' @return Summary printed to console
#' @method print tmbfit
#' @export
print.tmbfit <- function(x, ...){
  iter <- dim(x$samples)[1]
  chains <- dim(x$samples)[2]
  pars <- dim(x$samples)[3]-1
  samples <- (iter-x$warmup)*chains
  cat(paste0("Model '", x$model,"'", " has ", pars,
             " pars and was fit using ", x$algorithm,
             " with a '", x$metric, "' metric\n"))
  cat(paste0(chains," chain(s) of ", iter, " total iterations (", x$warmup, " warmup) were used\n"))
  rt <- c(sum(x$timing$time.total)/chains,  max(x$timing$time.total))
  ru <- 'seconds'
  if(rt[2]>60*60*24) {
    rt <- rt/(60*60*24); ru <- 'days'
  } else if(rt[2]>60*60) {
    rt <- rt/(60*60); ru <- 'hours'
  } else if(rt[2]>60){
    rt <- rt/60; ru <- 'minutes'
  }
  cat("Run time per chain: average=", round(rt[1],2), "and max=", round(rt[2],2), ru, '\n')
  if(!is.null(x$monitor)){
    minESS <- round(min(x$monitor$ess_bulk),1)
    maxRhat <- round(max(x$monitor$rhat),3)
    if(is.finite(minESS) & is.finite(maxRhat)){
      cat(paste0("Min bulk ESS=",
                 minESS,
                 " (",
                 round(100*minESS/samples,2),
                 "%) [", x$monitor$variable[which.min(x$monitor$ess_bulk)],
                 "] and maximum Rhat=", maxRhat, " [",
                 x$monitor$variable[which.max(x$monitor$rhat)],
                 ']\n'))
      if(minESS<200 | maxRhat > 1.1)
        cat('!! Warning: Signs of non-convergence found. Do not use for inference !!\n')
    } else {
      warning("ESS and Rhat calculations are not finite. Check model and rerun")
    }
  }
  if(x$algorithm=='SNUTS'){
    ndivs <- sum(extract_sampler_params(x)[,'divergent__'])
    cat(paste0("There were ", ndivs, " divergences after warmup\n"))
  }
}
