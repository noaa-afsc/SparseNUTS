# Package index

## All functions

- [`as.data.frame(`*`<tmbfit>`*`)`](https://noaa-afsc.github.io/SparseNUTS/reference/as.data.frame.tmbfit.md)
  :

  Convert object of class tmbfit to data.frame. Calls `extract_samples`

- [`as.tmbfit()`](https://noaa-afsc.github.io/SparseNUTS/reference/as.tmbfit.md)
  : Construtor for tmbfit objects

- [`benchmark_metrics()`](https://noaa-afsc.github.io/SparseNUTS/reference/benchmark_metrics.md)
  : Calculate gradient timings on a model for different metrics

- [`check_snuts_diagnostics()`](https://noaa-afsc.github.io/SparseNUTS/reference/check_snuts_diagnostics.md)
  : Check NUTS diagnostics of a fitted model

- [`.get_Q_stats()`](https://noaa-afsc.github.io/SparseNUTS/reference/dot-get_Q_stats.md)
  :

  Get the joint precision matrix Q from an optimized TMB or RTMB obj,
  along with the standard errors and lower bound of the absolute maximum
  correlation. Done using efficient sparse methods, specifically the
  Takahashi approach implemented in the Takahashi_Davis function in the
  sparseinverse package.

- [`.get_Qinv_stats()`](https://noaa-afsc.github.io/SparseNUTS/reference/dot-get_Qinv_stats.md)
  : Get the covariance matrix Qinv from an optimized TMB or RTMB obj,
  along with the standard errors and lower bound of the absolute maximum
  correlation.

- [`.get_inits()`](https://noaa-afsc.github.io/SparseNUTS/reference/dot-get_inits.md)
  : Get a single initial value vector in untransformed model space

- [`.make_unique_names()`](https://noaa-afsc.github.io/SparseNUTS/reference/dot-make_unique_names.md)
  : Function to take a character vector of parameter names and force
  them to be unique by appending numbers in square brackets as needed

- [`.precondition()`](https://noaa-afsc.github.io/SparseNUTS/reference/dot-precondition.md)
  : Update algorithm for mass matrix.

- [`.print.mat.stats()`](https://noaa-afsc.github.io/SparseNUTS/reference/dot-print.mat.stats.md)
  : Print matrix stats

- [`extract_sampler_params()`](https://noaa-afsc.github.io/SparseNUTS/reference/extract_sampler_params.md)
  : Extract sampler parameters from a fit.

- [`extract_samples()`](https://noaa-afsc.github.io/SparseNUTS/reference/extract_samples.md)
  : Extract posterior samples from a model fit.

- [`get_post()`](https://noaa-afsc.github.io/SparseNUTS/reference/get_post.md)
  : Extract posterior samples from a tmbfit object

- [`get_sparse_correlation()`](https://noaa-afsc.github.io/SparseNUTS/reference/get_sparse_correlation.md)
  : Extract a single correlation from a sparse precision matrix

- [`is.tmbfit()`](https://noaa-afsc.github.io/SparseNUTS/reference/is.tmbfit.md)
  : Check object of class tmbfit

- [`launch_shinytmb()`](https://noaa-afsc.github.io/SparseNUTS/reference/launch_shinytmb.md)
  : Launch shinystan for a TMB fit.

- [`pairs(`*`<tmbfit>`*`)`](https://noaa-afsc.github.io/SparseNUTS/reference/pairs.tmbfit.md)
  : Plot pairwise parameter posteriors and optionally the MLE points and
  confidence ellipses.

- [`plot(`*`<tmbfit>`*`)`](https://noaa-afsc.github.io/SparseNUTS/reference/plot.tmbfit.md)
  : Plot object of class tmbfit

- [`plot_Q()`](https://noaa-afsc.github.io/SparseNUTS/reference/plot_Q.md)
  : Make an image plot showing the correlation (lower triangle) and
  sparsity (upper triangle).

- [`plot_marginals()`](https://noaa-afsc.github.io/SparseNUTS/reference/plot_marginals.md)
  : Plot marginal distributions for a fitted model

- [`plot_sampler_params()`](https://noaa-afsc.github.io/SparseNUTS/reference/plot_sampler_params.md)
  : Plot adaptation metrics for a fitted model.

- [`plot_uncertainties()`](https://noaa-afsc.github.io/SparseNUTS/reference/plot_uncertainties.md)
  : Plot MLE vs MCMC marginal standard deviations for each parameter

- [`print(`*`<tmbfit>`*`)`](https://noaa-afsc.github.io/SparseNUTS/reference/print.tmbfit.md)
  : Print summary of tmbfit object

- [`rmvnorm_Q()`](https://noaa-afsc.github.io/SparseNUTS/reference/rmvnorm_Q.md)
  : Simulate draws from a MVN using the precision matrix Q

- [`sample_Q()`](https://noaa-afsc.github.io/SparseNUTS/reference/sample_Q.md)
  : Generate approximate posterior samples from the sparse precision
  matrix Q, assuming multivariate normality.

- [`sample_inits()`](https://noaa-afsc.github.io/SparseNUTS/reference/sample_inits.md)
  : Function to generate random initial values from a previous fit

- [`sample_snuts()`](https://noaa-afsc.github.io/SparseNUTS/reference/sample_snuts.md)
  : NUTS sampling for TMB models using a sparse metric (BETA).

- [`summary(`*`<tmbfit>`*`)`](https://noaa-afsc.github.io/SparseNUTS/reference/summary.tmbfit.md)
  : Print summary of object of class tmbfit

- [`tmbfit()`](https://noaa-afsc.github.io/SparseNUTS/reference/tmbfit.md)
  : Constructor for the "tmbfit" class
