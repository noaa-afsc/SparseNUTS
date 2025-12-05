# Convert object of class snutsfit to data.frame. Calls [`extract_samples`](https://cole-monnahan-noaa.github.io/SparseNUTS/reference/extract_samples.md)

Convert object of class snutsfit to data.frame. Calls
[`extract_samples`](https://cole-monnahan-noaa.github.io/SparseNUTS/reference/extract_samples.md)

## Usage

``` r
# S3 method for class 'snutsfit'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  Fitted object from `sample_rwm`

- row.names:

  Ignored

- optional:

  Ignored

- ...:

  Ignored

## Value

A data frame with parameters as columns and samples as rows.

## Details

This calls the default settings of
[`extract_samples`](https://cole-monnahan-noaa.github.io/SparseNUTS/reference/extract_samples.md),
no warmup samples and no column for the log-posterior (lp\_\_). Use this
function directly for finer control.
