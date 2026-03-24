# Plot method for prediction error curves of a peperr object

Plots individual and aggregated prediction error estimates based on
bootstrap samples.

## Usage

``` r
Plot.peperr.curves(
  x,
  at.risk = TRUE,
  allErrors = FALSE,
  bootRuns = FALSE,
  bootQuants = TRUE,
  bootQuants.level = 0.95,
  leg.cex = 0.7,
  ...
)
```

## Arguments

- x:

  `peperr` object.

- at.risk:

  number at risk to be display. default is TRUE.

- allErrors:

  Display .632, no information and average out-of-bag error in addition.
  default is FALSE.

- bootRuns:

  Display individual out-of-bag bootstrap samples. default is FALSE.

- bootQuants:

  Display pointwise out-of-bag bootstrap quantiles as shaded area.
  default is TRUE.

- bootQuants.level:

  Quantile probabilities for pointwise out-of-bag bootstrap quantiles.
  default is 0.95, i.e. 2.5% and 97.5% quantiles.

- leg.cex:

  size of legend text

- ...:

  additional arguments, not used.

## Details

This function is literally taken from `plot.peperr` in the `peperr`
package. The display of prediction error curves is adapted to allow for
numbers at risk and pointwise bootstrap quantiles.

## References

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`peperr`](https://fbertran.github.io/peperr/reference/peperr.html)

## Author

Thomas Hielscher <t.hielscher@dkfz.de>

## Examples

``` r
if (FALSE) { # \dontrun{

# example from glmnet package
set.seed(10101)
library(glmnet)
library(survival)
library(peperr)

N=1000;p=30
nzc=p/3
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
fx=x[,seq(nzc)]%*%beta/3
hx=exp(fx)
ty=rexp(N,hx)
tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator
y=Surv(ty,1-tcens)

peperr.object <- peperr(response=y, x=x, 
                        fit.fun=fit.glmnet, args.fit=list(family="cox"), 
                        complexity=complexity.glmnet,  
                        args.complexity=list(family="cox",nfolds=10),
                        indices=resample.indices(n=N, method="sub632", sample.n=10))

# pointwise bootstrap quantiles and all error types
Plot.peperr.curves(peperr.object, allErrors=TRUE)

# individual bootstrap runs and selected error types
Plot.peperr.curves(peperr.object, allErrors=FALSE, bootRuns=TRUE)
} # }
```
