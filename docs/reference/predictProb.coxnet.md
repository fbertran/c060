# Extract predicted survival probabilities from a glmnet fit

Extracts predicted survival probabilities from survival model fitted by
glmnet, providing an interface as required by `pmpec`.

## Usage

``` r
# S3 method for class 'coxnet'
predictProb(object, response, x, times, complexity, ...)
```

## Arguments

- object:

  a fitted model of class `glmnet`

- response:

  a two-column matrix with columns named 'time' and 'status'. The latter
  is a binary variable, with '1' indicating death, and '0' indicating
  right censored. The function
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) in package
  survival produces such a matrix

- x:

  `n*p` matrix of covariates.

- times:

  vector of evaluation time points.

- complexity:

  lambda penalty value.

- ...:

  additional arguments, currently not used.

## Value

Matrix with probabilities for each evaluation time point in `times`
(columns) and each new observation (rows).

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

`predictProb.glmnet`,[`peperr`](https://fbertran.github.io/peperr/reference/peperr.html),
[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)

## Author

Thomas Hielscher <t.hielscher@dkfz.de>
