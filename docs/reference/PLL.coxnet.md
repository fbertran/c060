# Predictive partial log-likelihood for glmnet Cox PH model fit

Extracts the predictive partial log-likelihood from a glmnet Cox PH
model fit.

## Usage

``` r
# S3 method for class 'coxnet'
PLL(object, newdata, newtime, newstatus, complexity, ...)
```

## Arguments

- object:

  fitted model of class `coxnet`.

- newdata:

  `n_new*p` matrix of covariates.

- newtime:

  `n_new`-vector of censored survival times.

- newstatus:

  `n_new`-vector of survival status, coded with 0 and .1

- complexity:

  lambda penalty value.

- ...:

  additional arguments, not used.

## Value

Vector of length `n_new`

## Details

Used by function `peperr`, if function `fit.glmnet` and `family="cox"`
is used for model fit, which gives a class `coxnet` object. This is
basically a wrapper based on the `coxnet.deviance` function from package
`glmnet`.

## References

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, *Journal of Statistical Software*, Volume 62(5), pages
1–22. https://doi.org/10.18637/jss.v062.i05.

## Author

Thomas Hielscher <t.hielscher@dkfz.de>
