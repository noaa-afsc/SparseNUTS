# Extract posterior samples from a tmbfit object

Extract posterior samples from a tmbfit object

## Usage

``` r
get_post(x, invf, parnames, array = FALSE)
```

## Arguments

- x:

  A fitted tmbfit object

- invf:

  The inverse function to decorrelate the parameters

- parnames:

  A vector of parameter names, excluding lp\_\_

- array:

  Whether to return a data.frame (default) or array which is used in
  constructing other objects downstream
