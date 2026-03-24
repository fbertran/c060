# Interface for determination of penalty lambda in penalized regression model via cross-validation

Determines the amount of shrinkage for a penalized regression model
fitted by glmnet via cross-validation, conforming to the calling
convention required by argument `complexity` in `peperr` call.

## Usage

``` r
complexity.glmnet(response, x, full.data, ...)
```

## Arguments

- response:

  a survival object (with `Surv(time, status)`, or a binary vector with
  entries 0 and 1).

- x:

  `n*p` matrix of covariates.

- full.data:

  data frame containing response and covariates of the full data set.

- ...:

  additional arguments passed to `cv.glmnet` call such as `family`.

## Value

Scalar value giving the optimal lambda.

## Details

Function is basically a wrapper for `cv.glmnet` of package `glmnet`. A
n-fold cross-validation (default n=10) is performed to determine the
optimal penalty lambda. For Cox PH regression models the deviance based
on penalized partial log-likelihood is used as loss function. For binary
endpoints other loss functions are available as well (see
`type.measure`). Deviance is default. Calling `peperr`, the default
arguments of `cv.glmnet` can be changed by passing a named list
containing these as argument `args.complexity`. Note that only penalized
Cox PH (`family="cox"`) and logistic regression models
(`family="binomial"`) are sensible for prediction error evaluation with
package `peperr`.

## References

Friedman, J., Hastie, T. and Tibshirani, R. (2008) *Regularization Paths
for Generalized Linear Models via Coordinate Descent*,
<https://web.stanford.edu/~hastie/Papers/glmnet.pdf>  
*Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010*  
<https://www.jstatsoft.org/v33/i01/>  
Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
*Regularization Paths for Cox's Proportional Hazards Model via
Coordinate Descent, Journal of Statistical Software, Vol. 39(5) 1-13*  
<https://www.jstatsoft.org/v39/i05/>  
Porzelius, C., Binder, H., and Schumacher, M. (2009) *Parallelized
prediction error estimation for evaluation of high-dimensional models,
Bioinformatics, Vol. 25(6), 827-829.*  
Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`peperr`](https://fbertran.github.io/peperr/reference/peperr.html),
[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html)

## Author

Thomas Hielscher <t.hielscher@dkfz.de>
