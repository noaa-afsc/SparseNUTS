# Get the joint precision matrix Q from an optimized TMB or RTMB obj, along with the standard errors and lower bound of the absolute maximum correlation. Done using efficient sparse methods, specifically the Takahashi approach implemented in the Takahashi_Davis function in the sparseinverse package.

Get the joint precision matrix Q from an optimized TMB or RTMB obj,
along with the standard errors and lower bound of the absolute maximum
correlation. Done using efficient sparse methods, specifically the
Takahashi approach implemented in the Takahashi_Davis function in the
sparseinverse package.

## Usage

``` r
.get_Q_stats(obj = NULL, Q = NULL, getCF = FALSE)
```

## Arguments

- obj:

  An optimized TMB or RTMB object

- Q:

  An optional sparse precision matrix if previous calculated. Otherwise
  calculated internally.

- getCF:

  A flag to determine whether to calculate the HMC condition factor. Set
  to FASLE by default due to poor performance in high dimensions even
  for a sparse Q.

## Value

A list containing a sparse matrix Q, the standard errors derived from
the square root of the diagonal of the inverse of Q, and the lower bound
on the maximum correlation among parameters.
