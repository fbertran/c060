\name{stabsel}
\alias{stabsel}
\title{
function to estimate a stable set of variables  
}
\description{
Given a desired type I error rate and a stability path calculated with \code{stability.path} the function selects a stable set of variables.
}
\usage{
stabsel(x,error=0.05,type=c("pfer","pcer"),pi_thr=0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 an object of class "stabpath" as returned by the function \code{stabpath}.
}
  \item{error}{
the desired type I error level w.r.t. to the chosen type I error rate.
}
 \item{type}{
The type I error rate used for controlling the number falsely selected variables. If \code{type="pfer"} the per-family error rate is controlled and \code{error} corresponds to the expected number of type I errors.
Selecting \code{type="pfer"} and \code{error} in the range of $0 > \code{error} < 1$ will control the family-wise error rate, i.e. the probability that at least one variable in the estimated stable set has been falsely selected.
If \code{type="pcer"} the per-comparison error rate is controlled and \code{error} corresponds to the expected number of type I errors divided by the number variables. 
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
  \item{error}{the desired type I error level w.r.t. to the chosen type I error rate
  }
   \item{type}{the type I error rate 
  }
}
\author{
Martin Sill \
\email{m.sill@dkfz.de}
}
\references{
Meinshausen N. and B\"uhlmann P. (2010), \emph{Stability Selection, Journal of the Royal Statistical Society: Series B (Statistical Methodology) Volume 72, Issue 4, pages 417--473.}\cr

Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox and Generalized Linear Models, Journal of Statistical Software, Volume 62(5), pages 1--22.}
\url{http://www.jstatsoft.org/v62/i05/}
}
\seealso{ \code{\link{plot.stabpath},\link{stabpath}}
}
\examples{
\dontrun{
#gaussian
set.seed(1234)
x=matrix(rnorm(100*1000,0,1),100,1000)
y <- x[1:100,1:1000]\%*\%c(rep(2,5),rep(-2,5),rep(.1,990))
res <- stabpath(y,x,weakness=1,mc.cores=2)
stabsel(res,error=0.05,type="pfer")
}}
\keyword{stability selection}
