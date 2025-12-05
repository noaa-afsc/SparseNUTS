# Extract sampler parameters from a fit.

Extract information about NUTS trajectories, such as acceptance ratio
and treedepth, from a fitted object.

## Usage

``` r
extract_sampler_params(fit, inc_warmup = FALSE)
```

## Arguments

- fit:

  A list returned by `sample_admb`.

- inc_warmup:

  Whether to extract the warmup samples or not (default). Warmup samples
  should never be used for inference, but may be useful for diagnostics.

## Value

An invisible data.frame containing samples (rows) of each parameter
(columns). If multiple chains exist they will be rbinded together.

## Details

Each trajectory (iteration) in NUTS has associated information about the
trajectory: stepsize, acceptance ratio, treedepth, and number of
leapfrog steps. This function extracts these into a data.frame, which
may be useful for diagnosing issues in certain cases. In general, the
user should not need to examine them, or preferably should via
[`plot_sampler_params`](https://cole-monnahan-noaa.github.io/SparseNUTS/reference/plot_sampler_params.md)
or `launch_shinyadmb`.

## See also

`launch_shinyadmb`.

## Examples

``` r
fit <- readRDS(system.file('examples', 'fit.RDS', package='SparseNUTS'))
#> Warning: cannot open compressed file '', probable reason 'No such file or directory'
#> Error in gzfile(file, "rb"): cannot open the connection
sp <- extract_sampler_params(fit, inc_warmup=TRUE)
#> Error: object 'fit' not found
str(sp)
#> Error: object 'sp' not found
```
