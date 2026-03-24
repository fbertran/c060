# Get coefficients for a model

Get coefficients for a model after applying interval search for tuning
parameters

## Usage

``` r
# S3 method for class 'sum.intsearch'
coef(object, ...)
```

## Arguments

- object:

  an object as returned by the function `summary.intsearch`.

- ...:

  additional argument(s)

## Value

named vector of non-zero coeficients for the optimal lambda

## References

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`EPSGO`](https://fbertran.github.io/c060/reference/EPSGO.md),
[`summary.intsearch`](https://fbertran.github.io/c060/reference/summary.intsearch.md),[`plot.sum.intsearch`](https://fbertran.github.io/c060/reference/plot.sum.intsearc.md)

## Author

Natalia Becker \\ <natalia.becker@dkfz.de>
