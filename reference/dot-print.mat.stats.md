# Print matrix stats

Print matrix stats

## Usage

``` r
.print.mat.stats(stats, getCF)
```

## Arguments

- stats:

  A list of matrix stats

- getCF:

  A flag whether to calculate the HMC condition factor.

## Value

Prints the HMC condition factor, sparsity, and max absolute correlation

## Details

The HMC condition factor used is from Langmore et al. (2019). A NULL
value is returned if the matrix passed is a scalar or if getCF=FALSE.
