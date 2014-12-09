\name{fit.glmnet}
\alias{fit.glmnet}
\title{Interface function for fitting a penalized regression model with \code{glmnet}}
\description{
Interface for fitting penalized regression models for binary of survival endpoint using \code{glmnet}, conforming to the requirements for argument \code{fit.fun} in \code{peperr} call. 
}
\usage{
fit.glmnet(response, x, cplx, ...)
}
\arguments{
  \item{response}{a survival object (with \code{Surv(time, status)}, or a binary vector with entries 0 and 1).}
  \item{x}{\code{n*p} matrix of covariates.}
  \item{cplx}{lambda penalty value.}
  \item{\dots}{additional arguments passed to \code{glmnet} call such as \code{family}.}
}
\value{
glmnet object
}

\details{
Function is basically a wrapper for \code{glmnet} of package \pkg{glmnet}.
Note that only penalized Cox PH (\code{family="cox"}) and logistic regression models (\code{family="binomial"}) are sensible for prediction error
evaluation with package \code{peperr}.
}

\references{
   Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent},   \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
  Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
  \emph{Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent, Journal of Statistical Software, Vol. 39(5)
    1-13}\cr
  \url{http://www.jstatsoft.org/v39/i05/}\cr
   Porzelius, C., Binder, H., and Schumacher, M. (2009) 
   \emph{Parallelized prediction error estimation for evaluation of high-dimensional models,
   Bioinformatics, Vol. 25(6), 827-829.}\cr
Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox and Generalized Linear Models, Journal of Statistical Software, Volume 62(5), pages 1--22.}
\url{http://www.jstatsoft.org/v62/i05/}
}

\author{
Thomas Hielscher \
\email{t.hielscher@dkfz.de}
}
  
\seealso{ \code{\link[peperr]{peperr}}, \code{\link[glmnet]{glmnet}}}
\keyword{models} \keyword{penalized regression} 