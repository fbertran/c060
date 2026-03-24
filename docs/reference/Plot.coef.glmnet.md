# function to highlight the path of a pre-specified set of variables within the coefficient path

Creates several plots showing the coefficient path for the final model
of a cv.glmnet fit and highlights the path of a pre-specified set of
variables within the coefficient path.

## Usage

``` r
Plot.coef.glmnet(cvfit, betas)
```

## Arguments

- cvfit:

  an object of class "cv.glmnet" as returned by the function
  `cv.glmnet`.

- betas:

  a vector of names of variables; must be a subset of
  rownames(coef(cvfit)).

## Value

a list of four objects

- stable:

  a vector giving the positions of the estimated stable variables

- lambda:

  the penalization parameter used for the stability selection

- lpos:

  the position of the penalization parameter in the regularization path

- error:

  the desired type I error level w.r.t. to the chosen type I error rate

- type:

  the type I error rate

## References

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## Author

Manuela Zucknick \\ <m.zucknick@dkfz-heidelberg.de>

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1010)
n=1000;p=100
nzc=trunc(p/10)
x=matrix(rnorm(n*p),n,p)
beta=rnorm(nzc)
fx= x[,seq(nzc)] %*% beta
eps=rnorm(n)*5
y=drop(fx+eps)
px=exp(fx)
px=px/(1+px)
ly=rbinom(n=length(px),prob=px,size=1)
set.seed(1011)
cvob1=cv.glmnet(x,y)
Plot.coef.glmnet(cvob1, c("V1","V100"))
} # }
```
