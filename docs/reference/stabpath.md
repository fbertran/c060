# Stability path for glmnet models

The function calculates the stability path for glmnet models, e.g. the
selection probabilities of the features along the range of
regularization parameters.

## Usage

``` r
stabpath(
  y,
  x,
  size = 0.632,
  steps = 100,
  weakness = 1,
  mc.cores = getOption("mc.cores", 2L),
  ...
)
```

## Arguments

- y:

  response variable. Like for the glment function: Quantitative for
  `family="gaussian"` or `family="poisson"` (non-negative counts). For
  `family="binomial"` should be either a factor with two levels, or a
  two-column matrix of counts or proportions. For
  `family="multinomial"`, can be a `nc>=2` level factor, or a matrix
  with `nc` columns of counts or proportions. For `family="cox"`, `y`
  should be a two-column matrix with columns named 'time' and 'status'.
  The latter is a binary variable, with '1' indicating death, and '0'
  indicating right censored. The function
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) in package
  `survival` produces such a matrix

- x:

  input matrix. Like for the glmnet function: of dimension nobs x nvars;
  each row is an observation vector. Can be in sparse matrix format
  (inherit from class `"sparseMatrix"` as in package `Matrix`; not yet
  available for `family="cox"`)

- size:

  proportion of samples drawn in every subsample used for the stability
  selection.

- steps:

  number of subsamples used for the stability selection.

- weakness:

  weakness parameter used for the randomised lasso as described in
  Meinshausen and B\\uhlmann (2010). For each subsample the features are
  reweighted by a random weight uniformly sampled in \[weakness,1\].
  This additional randomisation leads to a more consistent estimation of
  the stable set of features.

- mc.cores:

  number of cores used for the parallelization. For unix like system the
  parallelization is done by forking using the function `mclapply`. For
  windows systems socket cluster are used.

- ...:

  further arguments that are passed to the `glmnet` function.

## Value

an object of class "stabpath", which is a list of three objects

- fit:

  the fit object of class "glmnet" as returned from the glmnet function
  when applied to the complete data set.

- stabpath:

  a matrix which represents the stability path.

- qs:

  a vector holding the values of the average number of non-zero
  coefficients w.r.t to the lambdas in the regularization path.

## References

Meinshausen N. and B\\uhlmann P. (2010), *Stability Selection, Journal
of the Royal Statistical Society: Series B (Statistical Methodology)
Volume 72, Issue 4, pages 417–473.*  

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)`,`[`stabsel`](https://fbertran.github.io/c060/reference/stabsel.md)`,`[`plot.stabpath`](https://fbertran.github.io/c060/reference/plot.stabpath.md)

## Author

Martin Sill <m.sill@dkfz.de>

## Examples

``` r
if (FALSE) { # \dontrun{
#gaussian
set.seed(1234)
x <- matrix(rnorm(100*1000,0,1),100,1000)
y <- x[1:100,1:1000]%*% c(rep(2,5),rep(-2,5),rep(.1,990))
res <- stabpath(y,x,weakness=1,mc.cores=2)
plot(res)

#binomial
y=sample(1:2,100,replace=TRUE)
res <- stabpath(y,x,weakness=1,mc.cores=2,family="binomial")
plot(res)
    
#multinomial
y=sample(1:4,100,replace=TRUE)
res <- stabpath(y,x,weakness=1,mc.cores=2,family="multinomial")
plot(res)
    
#poisson
N=100; p=1000
nzc=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
f = x[,seq(nzc)]%*%beta
mu=exp(f)
y=rpois(N,mu)
res <- stabpath(y,x,weakness=1,mc.cores=2,family="poisson")
plot(res)

#Cox
library(survival)
set.seed(10101)
N=100;p=1000
nzc=p/3
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
fx=x[,seq(nzc)]%*%beta/3
hx=exp(fx)
ty=rexp(N,hx)
tcens=rbinom(n=N,prob=.3,size=1)
y=cbind(time=ty,status=1-tcens)
res <- stabpath(y,x,weakness=1,mc.cores=2,family="cox")
plot(res)
} # }
```
