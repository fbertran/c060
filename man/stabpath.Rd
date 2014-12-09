\name{stabpath}
\alias{stabpath}
\title{
Stability path for glmnet models
}
\description{
The function calculates the stability path for glmnet models, e.g. the selection probabilities of the features along the range of regularization parameters.
}
\usage{
stabpath(y,x,size=0.632,steps=100,weakness=1,mc.cores=getOption("mc.cores", 2L),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
response variable. Like for the glment function: Quantitative for \code{family="gaussian"} or
          \code{family="poisson"} (non-negative counts). For
          \code{family="binomial"} should be either a factor with two
          levels, or a two-column matrix of counts or proportions. For
          \code{family="multinomial"}, can be a \code{nc>=2} level factor, or a
          matrix with \code{nc} columns of counts or proportions. For
          \code{family="cox"}, \code{y} should be a two-column matrix with
          columns named 'time' and 'status'. The latter is a binary
          variable, with '1' indicating death, and '0' indicating right
          censored. The function \code{Surv()} in package \code{survival}
          produces such a matrix
}
  \item{x}{
input matrix. Like for the glmnet function:
of dimension nobs x nvars; each row is an
observation vector. Can be in sparse matrix format (inherit
from class \code{"sparseMatrix"} as in package \code{Matrix}; not yet
available for \code{family="cox"})
}
  \item{size}{
proportion of samples drawn in every subsample used for the stability selection.
}
  \item{steps}{
number of subsamples used for the stability selection.
}
  \item{weakness}{
weakness parameter used for the randomised lasso as described in Meinshausen and B\"uhlmann (2010). 
For each subsample the features are reweighted by a random weight uniformly sampled in [weakness,1].
This additional randomisation leads to a more consistent estimation of the stable set of features.
}
  \item{mc.cores}{
number of cores used for the parallelization. For unix like system the parallelization is done by forking using the function \code{mclapply}. For windows systems socket cluster are used. 
}
  \item{...}{
further arguments that are passed to the \code{glmnet} function. 
}
}
\value{
an object of class "stabpath", which is a list of three objects
  \item{fit}{
  the fit object of class "glmnet" as returned from the glmnet function when applied to the complete data set.
  }
  \item{stabpath}{
  a matrix which represents the stability path.
  }
  \item{qs}{
  a vector holding the values of the average number of non-zero coefficients w.r.t to the lambdas in the regularization path.
  }
}
\author{
Martin Sill 
\email{m.sill@dkfz.de}
}
\references{
Meinshausen N. and B\"uhlmann P. (2010), \emph{Stability Selection, Journal of the Royal Statistical Society: Series B (Statistical Methodology) Volume 72, Issue 4, pages 417--473.}\cr

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox and Generalized Linear Models, Journal of Statistical Software, Volume 62(5), pages 1--22.}
\url{http://www.jstatsoft.org/v62/i05/}
}
\seealso{ \code{\link{glmnet},\link{stabsel},\link{plot.stabpath}}
}
\examples{\dontrun{
#gaussian
set.seed(1234)
x <- matrix(rnorm(100*1000,0,1),100,1000)
y <- x[1:100,1:1000]\%*\% c(rep(2,5),rep(-2,5),rep(.1,990))
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
f = x[,seq(nzc)]\%*\%beta
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
fx=x[,seq(nzc)]\%*\%beta/3
hx=exp(fx)
ty=rexp(N,hx)
tcens=rbinom(n=N,prob=.3,size=1)
y=cbind(time=ty,status=1-tcens)
res <- stabpath(y,x,weakness=1,mc.cores=2,family="cox")
plot(res)
}
}
\keyword{stability selection}
