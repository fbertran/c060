# Summary method for interval search models

Produces a summary of a fitted interval search model

## Usage

``` r
# S3 method for class 'intsearch'
summary(
  object,
  digits = max(3, getOption("digits") - 3),
  verbose = TRUE,
  first.n = 5,
  ...
)
```

## Arguments

- object:

  an object of class `intsearch` as returned by the function `EPSGO`.

- digits:

  digits after the comma

- verbose:

  default set to TRUE.

- first.n:

  show first.n entries , default 5.

- ...:

  additional argument(s)

## Value

A list of following elements

- info:

  a data frame of four objects for optimal models  

  - alpha - a vector of alphas

  - lambda - a vector of penalization parameter lambda

  - deviances - a vector of deviances

  - n.features -a vector of number of features selected in each optimal
    model

- opt.alpha:

  an optimal value for alpha

- opt.lambda:

  an optimal value for lambda

- opt.error:

  an optimal value for error, hier minimal diviance

- opt.models:

  a list of optimal models with the same optimal error

## References

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`EPSGO`](https://fbertran.github.io/c060/reference/EPSGO.md),[`plot.sum.intsearch`](https://fbertran.github.io/c060/reference/plot.sum.intsearc.md)

## Author

Natalia Becker \\ <natalia.becker@dkfz.de>
