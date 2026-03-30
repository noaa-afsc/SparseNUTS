# Get the covariance matrix Qinv from an optimized TMB or RTMB obj, along with the standard errors and lower bound of the absolute maximum correlation.

Get the covariance matrix Qinv from an optimized TMB or RTMB obj, along
with the standard errors and lower bound of the absolute maximum
correlation.

## Usage

``` r
.get_Qinv_stats(obj = NULL, Qinv = NULL, getCF = FALSE)
```

## Arguments

- obj:

  An optimized TMB or RTMB object

- Qinv:

  An optional covariance matrix if previous calculated. Otherwise
  calculated internally.

- getCF:

  A flag whether to calculate the HMC condition factor. Set to FALSE by
  default because it is slow in high dimensions.

## Value

A list containing a covariance matrix Qinv, the standard errors derived
from the square root of the diagonal of Qinv, and the lower bound on the
maximum correlation among parameters.
