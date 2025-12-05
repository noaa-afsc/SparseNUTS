# Function to generate random initial values from a previous fit using SparseNUTS

Function to generate random initial values from a previous fit using
SparseNUTS

## Usage

``` r
sample_inits(fit, chains)
```

## Arguments

- fit:

  An outputted list from `sample_admb`

- chains:

  The number of chains for the subsequent run, which determines the
  number to return.

## Value

A list of lists which can be passed back into `sample_admb`.
