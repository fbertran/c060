\name{stability.selection}
\alias{stability.selection}
\title{
function to estimate a stable set of variables  
}
\description{
Given a desired family-wise error rate (FWER) and a stability path calculated with \code{stability.path} the function selects a stable set of features.
}
\usage{
stability.selection(x,fwer,pi_thr=0.6)
}
%- maybe also 'usage' for other objects documented here.
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
\seealso{ \code{\link{plot.stabpath},\link{stability.path}}
}
\examples{
\dontrun{
#gaussian
set.seed(1234)
x=matrix(rnorm(100*1000,0,1),100,1000)
y <- x[1:100,1:1000]\%*\%c(rep(2,5),rep(-2,5),rep(.1,990))
res <- stability.path(y,x,weakness=1,mc.cores=2)
stability.selection(res,fwer=.5)
}}
\keyword{stability selection}
