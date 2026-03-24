# Plot Summary object for interval search models

Produces a plot for summary object of a fitted interval search model.
Plot 'visited' points against iteration steps. start.N points are
initial points selected before interval search starts.

## Usage

``` r
# S3 method for class 'sum.intsearch'
plot(x, type = "summary", startN = 21, ...)
```

## Arguments

- x:

  an object of class `sum.intsearch` as returned by the function
  [`summary.intsearch`](https://fbertran.github.io/c060/reference/summary.intsearch.md).

- type:

  type of plot to be drawn, `type="summary"` will plot the partial log
  likelihood deviance as a function of both tuning parameters \\\alpha\\
  and log \\\lambda\\. The final solution will be highlighted by solid
  red line. Alternativly, `type="points"` will draw the distribution of
  initial and visited points of the interval search plotted in
  chronological order.

- startN:

  number of initial points. Needed if `type="points"`

- ...:

  additional argument(s)

## References

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`EPSGO`](https://fbertran.github.io/c060/reference/EPSGO.md),
[`summary.intsearch`](https://fbertran.github.io/c060/reference/summary.intsearch.md)

## Author

Natalia Becker \\ <natalia.becker@dkfz.de>
