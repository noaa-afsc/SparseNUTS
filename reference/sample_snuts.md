# NUTS sampling for TMB models using a sparse metric (BETA).

NUTS sampling for TMB models using a sparse metric (BETA).

## Usage

``` r
sample_snuts(
  obj,
  num_samples = 1000,
  num_warmup = NULL,
  chains = 4,
  cores = chains,
  thin = 1,
  adapt_stan_metric = NULL,
  control = NULL,
  seed = NULL,
  laplace = FALSE,
  init = c("last.par.best", "random", "random-t", "unif"),
  metric = c("auto", "unit", "diag", "dense", "sparse", "stan", "sparse-naive"),
  skip_optimization = FALSE,
  Q = NULL,
  Qinv = NULL,
  globals = NULL,
  model_name = NULL,
  refresh = NULL,
  print = TRUE,
  rotation_only = FALSE,
  lower = NULL,
  upper = NULL,
  iter = 2000,
  warmup = floor(iter/2),
  ...
)
```

## Arguments

- obj:

  The TMB object with random effects turned on and optimized.

- num_samples:

  The number of post-warmup iterations to run per chain.

- num_warmup:

  The number of warmup iterations to run per chain. The default of NULL
  indicates to automatically determine it based on other settings
  (recommended).

- chains:

  Number of chains

- cores:

  Number of parallel cores to use, defaults to `chains` so set this to 1
  to execute serial chains.

- thin:

  The thinning rate (defaults to 1).

- adapt_stan_metric:

  A boolean whether Stan should engage diagonal mass matrix adaptation.
  The default of NULL indicates to automatically select depending on
  other settings. See details.

- control:

  NUTS control list, currently available options are 'adapt_delta',
  'max_treedepth', and 'metric' which is the type of metric adaptation
  for Stan to do with options ('unit_e', 'diag_e', or 'dense_e'). For
  dense and sparse metrics this usually can be 'unit_e' to skip
  adaptation. NULL values (default) revert to `stan_sample` defaults.

- seed:

  Random number seed, used for generating initial values (if 'random")
  and for NUTS.

- laplace:

  Whether to leave the Laplace approximation on and only use NUTS to
  sample the fixed effects, or turn it off and sample from the joint
  parameter space (default).

- init:

  Either 'last.par.best' (default), 'random', 'random-t', or 'unif'. The
  former starts from the joint mode, while 'random' and 'random-t' draw
  from multivariate normal or multivariate t with 2 degrees of freedom
  distributions using the inverse joint precision matrix as a covariance
  matrix. 'random-t' is provided to allow for more dispersed initial
  values. 'unif' will draw U(-2,2) samples for all parameters, similar
  ot Stan's default behavior. If the joint NLL is undefined at the
  initial values then the model will exit and return the initial vector
  for further investigation by the user, if desired. Note that
  [`stan_sample`](https://andrjohns.github.io/StanEstimators/reference/stan_sample.html)
  only allows for the same init vector for all chains currently. If a
  seed is specified it will be set and thus the inits used will be
  reproducible. The inits are also returned in the 'inits' slot of the
  fitted object.

- metric:

  A character specifying which metric to use. Defaults to "auto" which
  uses an algorithm to select the best metric (see details), otherwise
  one of "sparse", "dense", "diag", "unit", "Stan", or "sparse-naive"
  can be specified.

- skip_optimization:

  Whether to skip optimization or not (default). If the model is already
  optimized and Q available these redundant steps can be skipped.

- Q:

  The sparse precision matrix. It is calculated internally if not
  specified (default).

- Qinv:

  The dense matrix \\Q^{-1}\\. It is calculated internally if not
  specified (default).

- globals:

  A named list of objects to pass to new R sessions when running in
  parallel and using RTMB. Typically this is the `data` object for now.

- model_name:

  An optional character giving the model name. If NULL it will use the
  DLL name which for RTMB models is just 'RTMB'. The name is used only
  for printing.

- refresh:

  How often to print updates to console (integer). 0 will turn off
  printing. The default is 100.

- print:

  Whether to print summary of run (default) or not

- rotation_only:

  Whether to return only the rotation object (for debugging purposes)

- lower, upper:

  A vector of lower and upper bounds. Only valid for 'stan' and 'unit'
  metrics. See details. Can be finite or infinite.

- iter:

  (Deprecated) Total iterations to run (warmup + sampling)

- warmup:

  (Deprecated) Total warmup iterations. Defaults to `iter`/2 based on
  Stan defaults, but when using dense, sparse, or diag metrics a much
  shorter warmup can be used (e.g., 150), especially if paired with a
  'unit_e' Stan metric. Use
  [`plot_sampler_params`](https://cole-monnahan-noaa.github.io/SparseNUTS/reference/plot_sampler_params.md)
  to investigate warmup performance and adjust as necessary for
  subsequent runs.

- ...:

  Additional arguments to pass to
  [`StanEstimators::stan_sample`](https://andrjohns.github.io/StanEstimators/reference/stan_sample.html).

## Value

A fitted MCMC object of class 'snutsfit'

## Details

The **TMB metric** is used to decorrelate and descale the posterior
distribution before sampling with Stan's algorithms. The chosen metric
and the reasoning for its selection are printed to the console before
sampling begins.

**Metric Options**

- `'auto'`: This is the default setting. It uses an internal algorithm
  to determine the optimal metric for the model. The choice depends on
  the availability of the precision matrix (\\Q\\) and/or the covariance
  matrix (\\M=Q^{-1}\\), the extent of parameter correlations, and the
  speed of gradient calculations.

- `'dense'` and `'sparse'`: Both of these options decorrelate and
  descale the posterior. However, the `'sparse'` metric is more
  computationally efficient for models with high-dimensional, sparse
  precision matrices. For models without random effects the `'sparse'`
  option is not available

- `'diag'`: This option only descales the posterior, using the marginal
  standard deviations derived from the covariance matrix \\M\\. It does
  not account for correlations.

- `'unit'`: This option uses an identity matrix, which is the default
  behavior in Stan. Unlike the `'Stan'` option below, the model is still
  optimized to find the mode, and the \\Q\\ matrix is calculated. This
  ensures that a full `mle` object (containing the mode, standard
  errors, and correlations) is returned.

- `'sparse-naive'`: This metric is constructed to be mathematically
  equivalent to `'dense'` but is often computationally faster. It is
  generally recommended only for testing and exploration.

- `'stan'`: This is a special flag that reverts the sampler to the
  standard Stan behavior. It skips the optimization and all \\Q\\ matrix
  calculations and ensures that Stan's mass matrix adaptation is engaged
  during warmup.

**Important Distinction**

Note that the `metric` parameter described here is specific to **TMB**
and is distinct from the Stan metric, which is controlled via the
`control` list argument in the sampling function. Stan by default adapts
a diagonal mass matrix (metric_e) using a series of expanding windows.
If \\Q\\ is a good estimate of the global precision then this is not
needed and disabling Stan's metric adaptation is recommended. This can
be done by setting `adapt_stan_metric=FALSE`. If left at NULL Stan's
adaptation will only be done for metrics 'stan' and 'unit' because those
two options do not descale the posterior. In this case, it is
recommended to use a longer warmup period to account for this adaptive
procedure.

Parameter bounds can be passed via the `lower` and `upper` vector
arguments. These can be finite or infinite, with length equal to the
number of parameters which will depend on whether `laplace` is TRUE or
not. However, bounds cannot be used with a metric that decorrelates or
descales the posterior because the bounds will be applied after this
operation which does not make sense. The 'stan' and 'unit' metrics do
not do this and so can be used. However, in general it is currently
recommended to do parameter transformations directly in the model and
add a Jacobian adjustment as necessary.
