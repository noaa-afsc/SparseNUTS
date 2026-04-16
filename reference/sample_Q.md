# Generate approximate posterior samples from the sparse precision matrix Q, assuming multivariate normality.

Generate approximate posterior samples from the sparse precision matrix
Q, assuming multivariate normality.

## Usage

``` r
sample_Q(fit, nsamples = 4000, seed = NULL)
```

## Arguments

- fit:

  An object returned by `sample_snuts`

- nsamples:

  The number of samples to draw, defaulting to 4000.

- seed:

  An optional random seed

## Value

A data.frame of approximate samples

## Examples

``` r
fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
x <- sample_Q(fit)
```
