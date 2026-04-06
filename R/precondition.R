## NOTE: This file was moved over from the adnuts repo on
## 2025-12-04 and the git history can be found there.
## https://github.com/Cole-Monnahan-NOAA/adnuts/commit/33744b850b8ee0642b8c8095181ebb4b878e1495

#' Update algorithm for mass matrix.
#'
#' @param metric The metric to use
#' @param fn The current fn function.
#' @param gr The current gr function
#' @param inputs A list of inputs
#' @param y.cur The current parameter vector in unrotated (Y) space.
.precondition <- function(metric, fn, gr, inputs, y.cur){
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
        .precondition(metric='diag', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(rdiag)
    }

    if((is.null(Q) & is.null(Qinv)) |
       is.null(inputs$max_cor) |
       !is.finite(inputs$max_cor)){
      message("unit metric selected b/c Q and Qinv unavailable")
      runit <-
        .precondition(metric='unit', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(runit)
    }


    if(inputs$max_cor <= 0.2){
      message("diag metric selected b/c low max cor (", round(inputs$max_cor,4),")")
      rdiag <-
        .precondition(metric='diag', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(rdiag)
    }

    # high correlations exist
    # fixed effects only models or ELA
    if(is.null(Q) & !is.null(Qinv)){
      message("dense metric selected b/c no Q and high max cor (", round(inputs$max_cor,4),")")
      rdense <-
        .precondition(metric='dense', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
      return(rdense)
    }

    # random effects models with sparse Q and high correlation
    if(!is.null(Q)){
      if(nrow(Q)>500){
        # high dimension, skip dense calcs since sparse almost always faster
        message("sparse metric selected b/c high dimesions and high max cor (", round(inputs$max_cor,4),")")
        rsparse <-
          .precondition(metric='sparse', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
        return(rsparse)
      } else {
        # low dimension sometimes dense is faster so check it
        if(!requireNamespace("microbenchmark", quietly=TRUE)){
          message("sparse metric selected b/c no timing available. Please install microbenchmark")
          rsparse <-
            .precondition(metric='sparse', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
          ## check for speed differences
          return(rsparse)
        } else {
          # when doing timing need to add random components to
          # input, otherwise TMB may skip calculations and throw
          # off benchmarking
          rsparse <-
            .precondition(metric='sparse', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
          if(is.null(Qinv)) inputs$Qinv <- as.matrix(Matrix::solve(Q))
          rdense <-
            .precondition(metric='dense', fn=fn, gr=gr, inputs=inputs, y.cur=y.cur)
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

