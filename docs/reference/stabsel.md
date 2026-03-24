# function to estimate a stable set of variables

Given a desired type I error rate and a stability path calculated with
`stability.path` the function selects a stable set of variables.

## Usage

``` r
stabsel(x, error = 0.05, type = c("pfer", "pcer"), pi_thr = 0.6)
```

## Arguments

- x:

  an object of class "stabpath" as returned by the function `stabpath`.

- error:

  the desired type I error level w.r.t. to the chosen type I error rate.

- type:

  The type I error rate used for controlling the number falsely selected
  variables. If `type="pfer"` the per-family error rate is controlled
  and `error` corresponds to the expected number of type I errors.
  Selecting `type="pfer"` and `error` in the range of \$0 \> `error` \<
  1\$ will control the family-wise error rate, i.e. the probability that
  at least one variable in the estimated stable set has been falsely
  selected. If `type="pcer"` the per-comparison error rate is controlled
  and `error` corresponds to the expected number of type I errors
  divided by the number variables.

- pi_thr:

  the threshold used for the stability selection, should be in the range
  of \$0.5 \> pi_thr \< 1\$.

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

Meinshausen N. and B\\uhlmann P. (2010), *Stability Selection, Journal
of the Royal Statistical Society: Series B (Statistical Methodology)
Volume 72, Issue 4, pages 417–473.*  

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), *c060: Extended
Inference with Lasso and Elastic-Net Regularized Cox and Generalized
Linear Models, Journal of Statistical Software, Volume 62(5), pages
1–22.* https://doi.org/10.18637/jss.v062.i05.

## See also

[`plot.stabpath`](https://fbertran.github.io/c060/reference/plot.stabpath.md)`,`[`stabpath`](https://fbertran.github.io/c060/reference/stabpath.md)

## Author

Martin Sill \\ <m.sill@dkfz.de>

## Examples

``` r
if (FALSE) { # \dontrun{
#gaussian
set.seed(1234)
x=matrix(rnorm(100*1000,0,1),100,1000)
y <- x[1:100,1:1000]%*%c(rep(2,5),rep(-2,5),rep(.1,990))
res <- stabpath(y,x,weakness=1,mc.cores=2)
stabsel(res,error=0.05,type="pfer")
} # }
```
