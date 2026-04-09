# Function to generate random initial values from a previous fit

Function to generate random initial values from a previous fit

## Usage

``` r
sample_inits(fit, chains)
```

## Arguments

- fit:

  A fitted object from
  [`sample_snuts`](https://noaa-afsc.github.io/SparseNUTS/reference/sample_snuts.md)

- chains:

  The number of chains for the subsequent run, which determines the
  number to return.

## Value

A list of vectors which can be passed back into
[`sample_snuts`](https://noaa-afsc.github.io/SparseNUTS/reference/sample_snuts.md).
