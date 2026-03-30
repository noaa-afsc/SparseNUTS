# Extract a single correlation from a sparse precision matrix

Extract a single correlation from a sparse precision matrix

## Usage

``` r
get_sparse_correlation(Q, i, j)
```

## Arguments

- Q:

  The sparse precision

- i, j:

  The integer indices for the correlation

## Value

A single scalar correlation value

## Details

This function extracts correlations without having to invert Q using
sparse matrix routines. Used for plotting outputs only.
