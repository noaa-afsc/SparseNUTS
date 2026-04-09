# SparseNUTS

The goal of SparseNUTS is to provide an interface for performing
Bayesian analyses using the no-U-turn (NUTS) MCMC algorithm for TMB and
RTMB models.

SparseNUTS uses the NUTS algorithm in the Stan software as linked
through the
[StanEstimators](https://github.com/andrjohns/StanEstimators) R package.
SparseNUTS implements the sparse NUTS (SNUTS) algorithm of Monnahan et
al. (2026) by decorrelating and descaling the posterior distribution
prior to passing to Stan. For models with high correlations and sparse
precision matrices, SNUTS can substantially improve sampling efficiency
and thus reduce run times. It also allows easy comparisons between
posterior samples and asymptotic estimates (via the Laplace
approximation) of variances and correlations for parameters and derived
quantities.

Importantly, it works smiliarly for TMB and RTMB models, and can be run
on models from existing packages built on either of these platforms
(e.g., glmmTMB).

## Installation

SparseNUTS is not currently on CRAN. First, install the StanEstimators:

``` r
# we recommend running this is a fresh R session or restarting your current session
install.packages('StanEstimators', repos = c('https://andrjohns.r-universe.dev', 'https://cloud.r-project.org'))
```

Now install the latest version of SparseNUTS as

``` r
remotes::install_github("noaa-afsc/SparseNUTS")
```

## Usage

Let `obj` be a list as created by `MakeADFun` in TMB or RTMB, then basic
usage is:

``` r
fit <- sample_snuts(obj)    # sample using SNUTS
print(fit)                  # print basic info 
plot(fit)                   # plot marginal fits
pairs(fit, order='slow')    # pairs plot showing posterior vs asymptotic estimates
post <- as.data.frame(fit)  # get posterior samples into data.frame
```

See package
[vignettes](https://noaa-afsc.github.io/SparseNUTS/articles/SparseNUTS.html)
for more information.

## Known issues

In some rare cases sampling will stall for some chains prior to
sampling. I suspect this is due to how the initial step size is
determined in Stan. A simple workaround is to specify an initial
stepsize via `stepsize=1` for instance.

## NOAA Enterprise GitHub Disclaimer

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

## References

Monnahan, C. C., Kristensen, K., Thorson, J. T., & Carpenter, B. (2026).
Leveraging Sparsity to Improve No-U-Turn Sampling Efficiency for
Hierarchical Bayesian Models. arXiv preprint arXiv:2603.02437.
