## NOTE: This file was copied over from adnuts on 2025-12-04
#' Plot pairwise parameter posteriors and optionally the MLE points and
#' confidence ellipses.
#'
#' @param x A list as returned by \code{sample_nuts}.
#' @param pars A character vector of parameters or integers
#'   representing which parameters to subset. Useful if the model
#'   has a larger number of parameters and you just want to show
#'   a few key ones.
#' @param order The order to consider the parameters. Options are
#'   'orig' (default) to use the order declared in the model, or
#'   'slow' and 'fast' which are based on the effective sample
#'   sizes ordered by slowest or fastest mixing
#'   respectively. 'mismatch' sorts by parameters with
#'   large discrepancies between the MLE and posterior marginal
#'   variances, defined as the absolute relative difference of
#'   the MLE from the posterior i.e., abs((mle-post)/post).
#'   Finally, 'cor' orders by the largest maximum absolute
#'   pairwise posterior correlation (including lp__).
#'   See example for usage.
#' @param inc_warmup Whether to include the warmup samples or not
#'   (default).
#' @param diag What type of plot to include on the diagonal,
#'   options are 'acf' which plots the autocorrelation function
#'   \code{acf}, 'hist' shows marginal posterior histograms, and
#'   'trace' the trace plot.
#' @param acf.ylim If using the acf function on the diagonal,
#'   specify the y limit. The default is c(-1,1).
#' @param ymult A vector of length ncol(posterior) specifying how
#'   much room to give when using the hist option for the
#'   diagonal. For use if the label is blocking part of the
#'   plot. The default is 1.3 for all parameters.
#' @param axis.col Color of axes
#' @param label.cex Control size of outer and diagonal labels
#'   (default 1)
#' @param limits A list containing the ranges for each parameter
#'   to use in plotting.
#' @param add.mle Boolean whether to add 95% confidence ellipses
#' @param add.monitor Boolean whether to print effective sample
#' @param add.inits Boolean whether to add the initial values to the plot
#' @param point.col Color of posterior points. Default selects level of transparency depending on the number of samples.
#' @param point.pch Shape of posterior points. Defaults to 16 if unspecified.
#' @param unbounded Whether to use the bounded or unbounded
#'   version of the parameters.  size (ESS) and Rhat values on
#'   the diagonal.
#' @param ... Arguments to be passed to plot call in lower
#'   triangular panels (scatterplots).
#' @method pairs tmbfit
#' @return Produces a plot, and returns nothing.
#' @details This function is modified from the base \code{pairs}
#'   code to work specifically with fits from the \code{\link{sample_snuts}}
#'   function using the SNUTS algorithms. If an
#'   invertible Hessian was found (in \code{fit$mle}) then
#'   estimated covariances are available to compare and added
#'   automatically (red ellipses). Likewise, a "monitor" object
#'   from \code{\link[posterior]{summarize_draws}} is attached as \code{fit$monitor}
#'   and provides effective sample sizes (ESS) and Rhat
#'   values. The ESS are used to potentially order the parameters
#'   via argument \code{order}, but also printed on the diagonal.
#' @export
#' @author Cole Monnahan
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#' pairs(fit)
#' pairs(fit, pars=1:2)
#' pairs(fit, pars=c(2,1))
#' pairs(fit, pars=1:2, order='slow')
#' pairs(fit, pars=1:2, order='fast')
#' pairs(fit, pars=1:2, order='mismatch')
#'
pairs.tmbfit <- function(x,
                           pars=NULL,
                           order=c('orig', 'slow', 'fast', 'mismatch', 'cor'),
                           inc_warmup=FALSE,
                           diag=c("trace","acf","hist"),
                           acf.ylim=c(-1,1), ymult=NULL, axis.col=gray(.5),
                           label.cex=.8, limits=NULL,
                           add.mle=TRUE, add.monitor=TRUE, add.inits=FALSE,
                           point.col=NULL, point.pch=NULL,
                           unbounded=FALSE, ...){
  fit <- x
  if(unbounded | !add.mle){
    mle <- NULL
  } else {
    mle <- fit$mle
  }
  posterior <- extract_samples(fit, inc_lp=TRUE,
                               inc_warmup=inc_warmup,
                               unbounded=unbounded)
  if(!inc_warmup){
    chains <- rep(1:dim(fit$samples)[2],
                  each=dim(fit$samples)[1]-fit$warmup)
  } else {
    chains <- rep(1:dim(fit$samples)[2],
                  each=dim(fit$samples)[1])

  }
  divs <- if(fit$algorithm=="NUTS")
            extract_sampler_params(fit, inc_warmup=inc_warmup)$divergent__ else NULL
  ptcex <- .2
  divcex <- .75
  chaincols <- 1:length(unique(chains))
  ## reset to old par when exiting
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  diag <- match.arg(diag)
  order <- match.arg(order)
  par.names <- names(posterior)
  ess <- fit$monitor$ess_bulk
  Rhat <- fit$monitor$rhat
  if(add.inits){
    if(is.null(fit$inits))  stop("add.inits not possible when fit$inits slot is empty")
    if(!is.list(fit$inits)) stop("add.inits not possible when fit$inits is not a list")
    if(dim(fit$samples)[2] != length(fit$inits)) stop("Length of inits does not match number of chains")
  }
  if(!is.null(pars) & length(pars)<=1)
    stop("pars argument <=1, only makes sense for >=2")
  if(is.null(ess))
    warning("No monitor information found in fitted object so ESS and Rhat not available. See details of help.")

  # reorder the parameter name vector
  if(is.character(pars[1]) & order!='orig'){
    warning("Ignoring 'order' argument because parameter names supplied in 'pars'")
    order <- 'orig'
  }
  if(order=='orig'){
    ## do nothing
    ind <- seq_along(par.names)
  } else if(order %in% c('slow', 'fast')){
    if(is.null(ess))
      stop("No effective sample sizes found so cannot order by slow/fast.")
    ## Get slowest or fastest parameter names
      ind <- order(ess, decreasing=(order=='fast'))
  } else if(order=='mismatch'){
      if(is.null(fit$mle$se))
        stop("SEs unavailable so mismatch option fails")
      tmp <- plot_uncertainties(fit, log=FALSE, plot=FALSE)
      x <- abs((tmp$sd.mle-tmp$sd.post)/tmp$sd.post)
      ind <- order(x, decreasing=TRUE)
  } else if(order=='cor'){
      post.cor <- cor(posterior)
      diag(post.cor) <- 0 # zero out so can take max along rows
      max.cors <- sapply(1:ncol(post.cor),
                         function(i) post.cor[i,which.max(abs(post.cor[i,]))])
      ind <- order(abs(max.cors), decreasing=TRUE)
  }
  # everything below here matches by par name so this is the only
  # reordering needed
  par.names <- par.names[ind]

  ## if(!(NCOL(posterior) %in% c(mle$nopar, mle$nopar+1)))
  ##   stop("Number of parameters in posterior and mle not the same")
  ## pars will either be NULL, so use all parameters. OR a vector of
  ## indices OR a vector of characters. Want to force lp__ to be the very
  ## last one stylistically, and b/c there is no ellipse for it.
  if(is.null(pars)){
    ## Use all or first 5
    pars <- par.names[1:min(5,NCOL(posterior))]
  } else if(is.numeric(pars[1])){
    ## Index can be used instead of character names. Note this
    ## can be sorted from above
    pars <- par.names[pars]
  }
  ## Now pars is character and possibly reordered
  pars.bad <- match(x=pars, table=names(posterior))
  if(any(is.na(pars.bad))){
    warning("Some par names did not match -- dropped")
    print(pars.bad)
    pars <- pars[!is.na(pars.bad)]
  }
  ## Converts character to index which is used throughout to
  ## subset when looping
  pars.ind <- match(x=pars, table=names(posterior))
  n <- length(pars.ind)
  n.mle <- ifelse(is.null(mle$cor), 0, NROW(mle$cor))
  if(n==1) stop("This function is only meaningful for >1 parameter")
  if(is.null(ymult)) ymult <- rep(1.3, n)
  ## If no limits given, calculate the max range of the posterior samples and
  ## parameter confidence interval.
  if(is.null(limits)){
    limits <- list()
    for(i in 1:n){
      if(pars.ind[i]<=n.mle){
        limit.temp <- mle$est[pars.ind[i]] +
          c(-1,1)*1.96*mle$se[pars.ind[i]]
      } else {
        limit.temp <- c(NA,NA)
      }
      ## multiplier for the ranges, adjusts the whitespace around the
      ## plots
      draws.tmp <- posterior[,pars.ind[i]]
      if(add.inits) draws.tmp <- c(draws.tmp, sapply(fit$inits, \(x) x[i]))
      min.temp <- min(draws.tmp, limit.temp[1], na.rm=TRUE)
      max.temp <- max(draws.tmp, limit.temp[2], na.rm=TRUE)
      margin <- .15*(max.temp-min.temp)
      limits[[i]] <- c(min.temp-margin, max.temp+margin)
    }
  }
  ## Change posterior point look depending on how many samples. Makes
  ## easier to read.
  N <- NROW(posterior)
  if(!is.null(point.pch)){
    mypch <- point.pch
  } else {
    mypch <- 16
  }
  if(!is.null(point.col)){
    mycol <- point.col;
  } else {
    if(N>=200){
      mycol <- rgb(0,0,0,.6)
    } else if(N>=1000){
      mycol <- rgb(0,0,0,.4)
    } else if(N>=4000){
      mycol <- rgb(0,0,0,.25)
    } else if(N>=10000){
      mycol <- rgb(0,0,0,.05)
    }
  }
  if(is.null(divs)) divs <- rep(0, N)
  par(mfrow=c(n,n), mar=0*c(.1,.1,.1,.1), yaxs="i", xaxs="i", mgp=c(.25, .25,0),
      tck=-.02, cex.axis=.65, col.axis=axis.col, oma=c(2, 2, 2,2))
  temp.box <- function() box(col=axis.col, lwd=.5)
  ## Row and col here are not the posterior, but the matrix of pairwise
  ## combinations
  for(row in 1:n){
    for(col in 1:n){
      ii <- pars.ind[row]
      jj <- pars.ind[col]
      ## Diagonal, so add user choice
      if(row==col){
        if(diag=="hist"){
          h <- hist(posterior[,ii], plot=F)
          ## Annoyingling you can't pass NULL to xlim in hist. So
          ## have to split up for two cases depending on limits.
          if(is.null(limits)){
            hist(posterior[,ii], axes=F, freq=FALSE, ann=F,
                 ylim=c(0, ymult[row]*max(h$density)),
                 col=gray(.8), border=gray(.5))
          } else {
            ## Else use the user provided limits
            hist(posterior[,ii], axes=F, freq=FALSE, ann=F,
                 ylim=c(0, ymult[row]*max(h$density)),
                 col=gray(.8), border=gray(.5), xlim=limits[[row]])
          }
          temp.box()
        } else if(diag=="acf") {
          acf(posterior[,ii], axes=F, ann=F, ylim=acf.ylim)
          temp.box()
        } else if(diag=="trace") {
          ## Trace plots for each chain separately
          xlim <- c(1, length(chains[chains==1]))
          plot(x=0, y=0,  type="n", axes=FALSE,
               ann=FALSE, ylim=limits[[row]], xlim=xlim)
          for(ll in unique(chains)){
            lines(posterior[chains==ll,ii], col=chaincols[ll], lwd=.1)
          }
          temp.box()
        }
        ## Add ESS and Rhat info to diagonal
        if(!is.null(ess) & !is.null(Rhat) & add.monitor)
          mtext(paste0('ESS=', round(ess[ii], 0), " Rhat=", format(round(Rhat[ii],2),nsmall=2)),
                cex=.8*label.cex, line=-1)
      }
      ## If lower triangle and covariance known, add scatterplot
      if(row>col){
        par(xaxs="r", yaxs="r")
        plot(x=posterior[,jj], y=posterior[,ii], axes=FALSE, ann=FALSE,
             pch=mypch, cex=ptcex, col=mycol, xlim=limits[[col]],
             ylim=limits[[row]], ...)
        ## replot divegences on top so they are always visible
        points(x=posterior[which(divs==1),jj], y=posterior[which(divs==1),ii],
               pch=mypch, cex=divcex, col='green')
        ## can only add MLE stuff if not lp__ parameter which
        ## doesn'th ave one
        if(ii<=n.mle & jj <=n.mle){
          ## Add bivariate 95% normal levels from MLE
          points(x=mle$est[jj], y=mle$est[ii],
                 pch=16, cex=.5, col='red')
          if(add.inits)
            lapply(fit$inits, \(init) points(init[jj], init[ii], col='blue', cex=1, pch=16))
          ## Get points of a bivariate normal 95% confidence contour
          if(!requireNamespace("ellipse", quietly=TRUE)){
            warning("ellipse package needs to be installed to show ellipses")
          } else {
            ellipse.temp <- ellipse::ellipse(x=mle$cor[jj, ii],
                                    scale=mle$se[c(jj, ii)],
                                    centre= mle$est[c(jj, ii)], npoints=1000,
                                    level=.95)
            lines(ellipse.temp , lwd=.5, lty=1, col="red")
          }
        }
        par(xaxs="i", yaxs="i")
        temp.box()
      }
      if(row<col){
        ## If upper triangle add text showing the empirical correlation
        plot(0,0,type="n", xlim=c(0,1), ylim=c(0,1), axes=F,ann=F)
        temp.cor <- round(cor(posterior[,c(ii,jj)])[1,2],2)
        ## Set a minimum limit for this, so they're still
        ## visible, but still a function of correlation. This
        ## might need to be dynamic with n.
        legend("center", legend=NA, title=temp.cor,
               cex=(3*abs(temp.cor)+.25)*.5, bty='n')
        ## text(.5,.5, labels=temp.cor, cex=(3*abs(temp.cor)+.5)*.9,
        ##      col=1)
        temp.box()
      }
      ## Add special cases of axes on the ends
      if(row==n) {
        par( mgp=c(.05, ifelse(col %% 2 ==0, 0, .5),0) )
        axis(1, col=axis.col, lwd=.5)
      }
      if(col==1 & row >1) {
        par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      if(col==1 & row ==1){
        par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      if(row==1) mtext(pars[col], line=ifelse(col %% 2 ==1, .1, 1.1),
                       cex=label.cex)
      if(col==n)
        mtext(pars[row], side=4, line=ifelse(row %% 2 ==1, 0, 1), cex=label.cex)
    }
  }
}



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
#'
plot_marginals <- function(fit, pars=NULL, mfrow=NULL,
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
  posterior <- extract_samples(fit, inc_lp=FALSE)
  par.names <- names(posterior)
  if(is.null(pars)) pars <- par.names
  if(is.character(pars[1])){
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
    if(!is.null(fit$mle)){
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
      label <- paste0('ESS=', round(mon[ii,'ess_bulk'],2), "\nRhat=", round(mon[ii,'rhat'],3))
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
