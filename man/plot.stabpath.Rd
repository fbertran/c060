\name{plot.stabpath}
\alias{plot.stabpath}
\title{
function to plot a stability path 
}
\description{
Given a desired family-wise error rate (FWER) and a stability path calculated with \code{stability.path} the function selects an stable set of features and plots the stability path and the corresponding regularization path.
}
\usage{
\method{plot}{stabpath}(x, fwer=0.5,pi_thr=0.6, xvar=c("lambda", "norm", "dev"), col.all="black", col.sel="red", ...)
}
\arguments{
  \item{x}{
 an object of class "stabpath" as returned by the function \code{stability.path}.
}
\item{fwer}{
 the desired family-wise error rate (fwer), e.g. the probability that at least one feature in the estimated set of stable features has been falsely selected. 
}
  \item{pi_thr}{
the threshold used for the stability selection, should be in the range of $0.5 > pi_thr < 1$.  
}
  \item{xvar}{
the variable used for the xaxis, e.g. for "lambda" the selection probabilities are plotted along the log of the penalization parameters,
for "norm" along the L1-norm and for "dev" along the fraction of explained deviance.
}
  \item{col.all}{
the color used for the variables that are not in the estimated stable set   
}
  \item{col.sel}{
the color used for the variables in the estimated stable set
}
 \item{...}{
further arguments that are passed to matplot
 }
}
\value{
a list of four objects
  \item{stable}{
  a vector giving the positions of the estimated stable variables 
  }
  \item{lambda}{
  the penalization parameter used for the stability selection 
  }
  \item{lpos}{
  the position of the penalization parameter in the regularization path
  }
  \item{fwer}{
  the controlled family-wise error rate (fwer), e.g. the probability that at least one feature in the estimated set of stable features     has been falsely selected.
  }
}
\author{
Martin Sill \
\email{m.sill@dkfz.de}
}
\references{
Meinshausen N. and B\"uhlmann P. (2010), Stability Selection, Journal of the Royal Statistical Society: Series B (Statistical Methodology) Volume 72, Issue 4, pages 417???473.
}
\seealso{ \code{\link{stability.selection},\link{stability.path}}
}
\examples{
\dontrun{
#gaussian
set.seed(1234)
x=matrix(rnorm(100*1000,0,1),100,1000)
y <- x[1:100,1:1000]\%*\%c(rep(2,5),rep(-2,5),rep(.1,990))
res <- stability.path(y,x,weakness=1,mc.cores=2)
plot(res,fwer=.5)
}
}
\keyword{stability selection}