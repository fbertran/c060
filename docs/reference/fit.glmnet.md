# Interface function for fitting a penalized regression model with `glmnet`

Interface for fitting penalized regression models for binary of survival
endpoint using `glmnet`, conforming to the requirements for argument
`fit.fun` in `peperr` call.

## Usage

``` r
fit.glmnet(response, x, cplx, ...)
```

## Arguments

- response:

  a survival object (with `Surv(time, status)`, or a binary vector with
  entries 0 and 1).

- x:

  `n*p` matrix of covariates.

- cplx:

  lambda penalty value.

- ...:

  additional arguments passed to `glmnet` call such as `family`.

## Value

glmnet object

## Details

Function is basically a wrapper for `glmnet` of package glmnet. Note
that only penalized Cox PH (`family="cox"`) and logistic regression
models (`family="binomial"`) are sensible for prediction error
evaluation with package `peperr`.

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
[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)

## Author

Thomas Hielscher <t.hielscher@dkfz.de>
