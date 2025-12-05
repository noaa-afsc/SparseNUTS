# Make an image plot showing the correlation (lower triangle) and sparsity (upper triangle).

Make an image plot showing the correlation (lower triangle) and sparsity
(upper triangle).

## Usage

``` r
plot_Q(fit, Q = NULL)
```

## Arguments

- fit:

  A fitted object

- Q:

  A sparse matrix. If NULL it will be extracted from `fit`.

## Value

A plot created by
[`Matrix::image`](https://rdrr.io/r/graphics/image.html).

## Details

This function is used to visualize the sparsity and correlation patterns
of the joint model. The upper triangle shows whether an element is 0
(white) or not (gray), while the lower triangle shows the correlation
calculated from `cov2cor(solve(Q))`.
