# Plot MLE vs MCMC marginal standard deviations for each parameter

Plot MLE vs MCMC marginal standard deviations for each parameter

## Usage

``` r
plot_uncertainties(fit, log = TRUE, plot = TRUE)
```

## Arguments

- fit:

  A fitted object returned by `sample_admb`

- log:

  Whether to plot the axes in log space (default TRUE).

- plot:

  Whether to plot it or not.

## Value

Invisibly returns data.frame with parameter name (row) and estimated
uncertainties for each method (columns).

## Details

It can be helpful to compare uncertainty estimates between the two
paradigms. This plots the marginal posterior standard deviation vs the
frequentist standard error estimated from the .cor file. Large
differences often indicate issues with one estimation method.

## Examples

``` r
fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#> Warning: cannot open compressed file '', probable reason 'No such file or directory'
#> Error in gzfile(file, "rb"): cannot open the connection
x <- plot_uncertainties(fit, plot=FALSE)
#> Error: object 'fit' not found
head(x)
#> Error: object 'x' not found
```
