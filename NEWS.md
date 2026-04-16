# SparseNUTS (development version)

# SparseNUTS 1.0.2 
- Changed default `init` option in `sample_snuts` to use 'auto' by default
- Added `sample_Q` function to efficiently sample approximate posterior draws from the joint precision matrix.
- Added `order` argument to `plot_marginals`
- Added `skip_monitor` argument to `sample_snuts`
- Improved console output
- Minor code and documentation improvements

# SparseNUTS 1.0.1
- Minimized dense matrix calculations throughout to reduce computational overhead of large models. This includes testing for the metric, calculating and returning correlations, sampling from Q or Qinv, etc.
- Added informative errors when the user forgets to pass global RTMB data or functions.
- Added checks and informative errors during marginal optimization.


# SparseNUTS 1.0.0

* Initial release
