# Function producing stratified/ balanced folds for cross validation

Get balanced folds for cross validation, which are used for tuning
penalization parameters

## Usage

``` r
balancedFolds(class.column.factor, cross.outer)
```

## Arguments

- class.column.factor:

  class labels of length n

- cross.outer:

  number of folds

## Value

- permutated.cut :

  vector of length n, indicating the fold belongs to

- model :

  model list

  - alpha - optimal alpha

  - lambda - optimal lambda

  - nfolds - cross-validation's folds

  - cvreg - `cv.glmnet` object for optimal alpha

  - fit - `glmnet` object for optimal alpha and optimal lambda

## References

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`EPSGO`](https://fbertran.github.io/c060/reference/EPSGO.md)

## Author

Natalia Becker natalia.becker at dkfz.de
