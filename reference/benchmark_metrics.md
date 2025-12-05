# Calculate gradient timings on a model for different metrics

Calculate gradient timings on a model for different metrics

## Usage

``` r
benchmark_metrics(obj, times = 1000, metrics = NULL, model_name = NULL)
```

## Arguments

- obj:

  A TMB object

- times:

  How many evaluations to do

- metrics:

  A character vector of different metrics to benchmark

- model_name:

  An optional character name for the model, if NULL will pull from the
  DLL name

## Value

A data.frame containing the median gradient timing (time), the percent
sparsity of \\Q\\ and the dimension of the model (npar).
