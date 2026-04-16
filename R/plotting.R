
#' Plot MLE vs MCMC marginal standard deviations for each
#' parameter
#'
#' @param fit A fitted object returned by
#'   \code{\link{sample_snuts}}
#' @param log Whether to plot the axes in log space (default TRUE).
#' @param plot Whether to plot it or not.
#' @details It can be helpful to compare uncertainty estimates
#'   between the two paradigms. This plots the marginal posterior
#'   standard deviation vs the frequentist standard error
#'   estimated from the .cor file. Large differences often
#'   indicate issues with one estimation method.
#' @return Invisibly returns data.frame with parameter name (row) and
#'   estimated uncertainties for each method (columns).
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#' x <- plot_uncertainties(fit, plot=FALSE)
#' head(x)
#' @export
plot_uncertainties <- function(fit, log=TRUE, plot=TRUE){
  stopifnot(is.tmbfit(fit))
  if(is.null(fit$mle$se))
    stop("MLE SEs not found so cannot plot them")
  sd.post <- apply(extract_samples(fit), 2, stats::sd)
  sd.mle <- fit$mle$se[1:length(sd.post)]
  pars <- fit$par_names[1:length(sd.post)]
  if(plot){
    plot(sd.post, sd.mle, log=ifelse(log, 'xy', ''),
         xlab='Posterior SD', ylab='MLE SE',
         main='Comparing Bayesian vs frequentist uncertainty estimates')
    abline(0,1)
  }
  df <- data.frame(par=pars, sd.post=sd.post, sd.mle=sd.mle)
  return(invisible(df))
}


#' Plot marginal distributions for a fitted model
#'
#' @param fit A fitted object returned by
#'   \code{\link{sample_snuts}}.
#' @param pars A numeric or character vector of parameters which
#'   to plot, for plotting a subset of the total (defaults to all)
#' @param order A string specifying the order to consider the
#'   parameters in. See \code{?pairs.tmbfit} for more
#'   information.
#' @param mfrow A custom grid size (vector of two) to be called
#'   as \code{par(mfrow)}, overriding the defaults.
#' @param add.mle Whether to add marginal normal distributions
#'   determined from the inverse Hessian file
#' @param add.monitor Whether to add ESS and Rhat information
#' @param breaks The number of breaks to use in \code{hist()},
#'   defaulting to 30
#' @export
#'
#' @details This function plots grid cells of all parameters
#'   in a model, comparing the marginal posterior histogram vs
#'   the asymptotic normal (red lines) from the inverse
#'   Hessian. Its intended use is to quickly gauge differences
#'   between frequentist and Bayesian inference on the same
#'   model.
#'
#' If \code{fit$monitor} exists the effective sample size
#' (ESS) and R-hat estimates are printed in the top right
#' corner. See
#' \url{https://mc-stan.org/rstan/reference/Rhat.html} for more
#' information. Generally Rhat>1.05 or ESS<100 (per chain)
#' suggest inference may be unreliable.
#'
#' This function is customized to work with multipage PDFs,
#' specifically:
#' \code{pdf('marginals.pdf', onefile=TRUE, width=7,height=5)}
#' produces a nice readable file.
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#' plot_marginals(fit, pars=1:2)
#' plot_marginals(fit)
#' plot_marginals(fit, pars=1:2)
#' plot_marginals(fit, pars=c(2,1))
#' plot_marginals(fit, pars=1:5, order='slow')
#' plot_marginals(fit, pars=1:2, order='fast')
#' plot_marginals(fit, pars=1:2, order='mismatch')
#'
plot_marginals <- function(fit, pars=NULL,  order='orig',
                           mfrow=NULL,
                           add.mle=TRUE, add.monitor=TRUE,
                           breaks=30){
  if(!is.tmbfit(fit)) stop("fit is not a valid object")
  if(!is.null(mfrow)) stopifnot(is.vector(mfrow) && length(mfrow)==2)
  stopifnot(add.mle %in% c(TRUE,FALSE))
  if(add.mle & is.null(fit$mle)) {
    add.mle <- FALSE
    warning("No MLE information found in fit$mle so cannot add")
  }
  if(!add.mle) fit$mle <- NULL
  if(!add.monitor) fit$monitor <- NULL
  par.old <- par()
  on.exit(par(mfrow=par.old$mfrow, mar=par.old$mar,
              mgp=par.old$mgp, oma=par.old$oma, tck=par.old$tck))
  posterior <- extract_samples(fit, inc_lp=TRUE)
  par.names <- names(posterior)
  if(is.null(pars)) {
    pars <- par.names
  }
  if(order!='orig'){
    pars <- pairs.tmbfit(fit, pars=pars, order=order, plot=FALSE)
  } else if(is.character(pars[1])){
    pars.ind <- match(x=pars, table=par.names)
    if(any(is.na(pars.ind))){
      warning("Some par names did not match -- dropped")
      print(pars[is.na(pars.ind)])
      pars.ind <- pars.ind[!is.na(pars.ind)]
    }
    pars <- pars.ind
  } else if(any(pars > NCOL(posterior))){
    warning("Some par numbers too big -- dropped")
    print(pars[pars > NCOL(posterior)])
    pars <- pars[ pars <=NCOL(posterior)]
  }
  n <- length(pars)
  stopifnot(is.numeric(pars[1]))
  stopifnot(ncol(posterior)>1)
  par(mar=c(1.5,0,.1,0), mgp=c(2,.4,0),
      oma=c(.25,.25,.25,.25), tck=-.02)
  if(!is.null(mfrow)){
    par(mfrow=mfrow)
  } else if(n>12){
    par(mfrow=c(4,4))
  } else if(n>9){
    par(mfrow=c(4,3))
  } else if(n>6){
    par(mfrow=c(3,3))
  } else if(n>4){
    par(mfrow=c(3,2))
  } else if(n>3){
    par(mfrow=c(2,2))
  } else {
    par(mfrow=c(1,n))
  }
  for(ii in pars){
    par <- par.names[ii]
    if(!is.null(fit$mle) & par!='lp__'){
      mle <- fit$mle$est[ii]
      se <-  fit$mle$se[ii]
      x1 <- seq(qnorm(.001, mle, se), qnorm(.999, mle, se), len=100)
      y1 <- dnorm(x1, mle, se)
    } else{
      x1 <- y1 <- NULL
    }
    tmp <- hist(posterior[,ii], plot=FALSE, breaks=breaks)
    x2 <- tmp$mids; y2 <- tmp$density
    plot(0,0, type='n', xlim=range(c(x1,x2)), yaxs='i',
         ylim=c(0, max(c(y1,y2))*1.3), axes=FALSE, ann=FALSE)
    hist(posterior[,ii], breaks=breaks, add=TRUE, yaxs='i', freq=FALSE, col=gray(.8))
    axis(1);  box(col=gray(.5));
    if(!is.null(fit$mle)) lines(x1,y1, col='red', lwd=2)
    if(!is.null(fit$monitor)){
      mon <- fit$monitor
      ## add ESS and Rhat to top right
      tmp <- par("usr"); xy <- c(.85,.88)
      text.x <- tmp[1]+xy[1]*diff(tmp[1:2])
      text.y <- tmp[3]+xy[2]*diff(tmp[3:4])
      label <- paste0('ESS=', round(mon[ii,'ess_bulk'],0), "\nRhat=", round(mon[ii,'rhat'],2))
      text(x=text.x, y=text.y, labels=label, cex=.8)
    }
    mtext(paste("",par), line=-1.6, adj=0, cex=.9)
  }
}



#' Plot adaptation metrics for a fitted model.
#'
#' @param fit A fitted object returned by
#' \code{\link{sample_snuts}}.
#' @param plot Whether to plot the results
#' @return Prints and invisibly returns a ggplot object
#'
#' @details This utility function quickly plots the adaptation output of NUTS
#' chains.
#' @importFrom rlang .data
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#' plot_sampler_params(fit)
plot_sampler_params <- function(fit, plot=TRUE){
  if(!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 package not found")
  sp <- SparseNUTS::extract_sampler_params(fit, inc_warmup=TRUE)
  sp.long <-
    data.frame(iteration=sp$iteration, chain=factor(sp$chain),
               value=c(sp$accept_stat__, log(sp$stepsize__),
                       sp$n_leapfrog__, sp$divergent__, sp$energy__),
               variable=rep(c('accept_stat', 'log_stepsize',
                              'n_leapfrog', 'divergent',
                              'energy'), each=nrow(sp)))
  g <- ggplot2::ggplot(sp.long, ggplot2::aes(.data$iteration, y=.data$value, color=.data$chain)) +
    ggplot2::geom_point(alpha=.5) +
    ggplot2::facet_wrap('variable', scales='free_y', ncol=1) + ggplot2::theme_bw()
  if(plot) print(g)
  return(invisible(g))
}




#' Make an image plot showing the correlation (lower triangle)
#' and sparsity (upper triangle).
#'
#' @param fit A fitted object
#' @param Q A sparse matrix. If NULL it will be extracted from
#'   \code{fit}.
#'
#' @details This function is used to visualize the sparsity and
#'   correlation patterns of the joint model. The upper triangle
#'   shows whether an element is 0 (white) or not (gray), while
#'   the lower triangle shows the correlation calculated from
#'   \code{cov2cor(solve(Q))}.
#' @return A plot created by \code{\link[Matrix]{image}}.
#' @export
plot_Q <- function(fit, Q=NULL){
  if(is.null(Q)){
    if(!is.tmbfit(fit)) stop("fit is not a valid fitted object")
    if(is.null(fit$mle$Q)) return(NULL)
    nn <- length(fit$par_names)
    if(is.null(fit$mle$Qinv)){
      corr <- stats::cov2cor(as.matrix(Matrix::solve(fit$mle$Q)))
    } else {
      corr <- stats::cov2cor(fit$mle$Qinv)
    }
    Q <- fit$mle$Q
  } else {
    corr <- stats::cov2cor(solve(Q))
  }
  Q[Q!=0] <- 1e-10
  Q[lower.tri(Q,TRUE)] <- corr[lower.tri(Q,TRUE)]
  Matrix::image(Q, useRaster=TRUE, at=seq(-1,1, len=50))
}
