# Sparse no-U-turn sampling for TMB and RTMB models

## Summary

`SparseNUTS` provides an interface for performing Bayesian analyses
using the no-U-turn (NUTS) MCMC algorithm (Hoffman and Gelman 2014) for
TMB (Kristensen et al. 2016) and RTMB models. It is closely related to
the `adnuts` package which provides functionality for ADMB (Fournier et
al. 2012) models.

For TMB models the NUTS algorithm is in the Stan software (Carpenter et
al. 2017) and linked through the
[StanEstimators](https://github.com/andrjohns/StanEstimators) R package.
`SparseNUTS` implements the sparse NUTS (SNUTS) algorithm of (Monnahan
et al. in prep) by decorrelating and descaling the posterior
distribution prior to passing to Stan. For models with high correlations
and sparse precision matrices, SNUTS can substantially improve sampling
efficiency. Importantly, it works smiliarly for TMB and RTMB models, and
can be run on models from existing packages built on either of these
platforms (e.g., `glmmTMB`).

This package was originally developed inside of the `adnuts` package but
was split off in late 2025 to have a separate package for TMB and ADMB
models. The `tmbstan` package also provides an interface to the Stan
software, but lacks the ability to decorrelate the target distribution
prior to sampling. `SparseNUTS` provides more flexible options related
to the mass matrix.

Basic usage for a TMB or RTMB object `obj` is shown as follows, with a
more in depth article on the website.

``` r
fit <- sample_snuts(obj)
print(fit)
plot(fit)
pairs(fit, order='slow')
post <- as.data.frame(fit)
launch_shinytmb(fit)
```

## References

Carpenter, B., A. Gelman, M. D. Hoffman, D. Lee, B. Goodrich, M.
Betancourt, A. Riddell, J. Q. Guo, P. Li, and A. Riddell. 2017. “Stan: A
Probabilistic Programming Language.” *Journal of Statistical Software*
76 (1): 1–29.

Fournier, D. A., H. J. Skaug, J. Ancheta, J. Ianelli, A. Magnusson, M.
N. Maunder, A. Nielsen, and J. Sibert. 2012. “AD Model Builder: Using
Automatic Differentiation for Statistical Inference of Highly
Parameterized Complex Nonlinear Models.” *Optimization Methods &
Software* 27 (2): 233–49.

Hoffman, M. D., and A. Gelman. 2014. “The No-u-Turn Sampler: Adaptively
Setting Path Lengths in Hamiltonian Monte Carlo.” *Journal of Machine
Learning Research* 15 (1): 1593–623.

Kristensen, Kasper, Anders Nielsen, Casper W. Berg, Hans Skaug, and
Bradley M. Bell. 2016. “TMB: Automatic Differentiation and Laplace
Approximation.” *Journal of Statistical Software* 70 (5): 21.

Monnahan, C. C., Thorson J. T., K. Kristensen, and B. Carpenter. in
prep. “Leveraging Sparsity to Improve No-u-Turn Sampling Efficiency for
Hierarchical Bayesian Models.” *arXiv Preprint*, in prep.
