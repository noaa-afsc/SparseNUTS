# Prepare inputs for sparse sampling

@param obj Object @param skip_optimization Whether to skip or not @param
laplace Whether to due the LA or not @param metric Which metric @param Q
Sparse precision @param Qinv Inverse of Q @return A list containing Q,
Qinv, the mle list, and timings

## Usage

``` r
.get_inputs(obj, skip_optimization, laplace, metric, Q, Qinv)
```
