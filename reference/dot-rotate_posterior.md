# Update algorithm for mass matrix.

Update algorithm for mass matrix.

## Usage

``` r
.rotate_posterior(metric, fn, gr, Q, Qinv, y.cur)
```

## Arguments

- metric:

  The metric to use

- fn:

  The current fn function.

- gr:

  The current gr function

- Q:

  The sparse precision matrix

- Qinv:

  The inverse of Q

- y.cur:

  The current parameter vector in unrotated (Y) space.
