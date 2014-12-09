\name{complexity.glmnet}
\alias{complexity.glmnet}
\title{Interface for determination of penalty lambda in penalized regression model via cross-validation}
\description{
Determines the amount of shrinkage for a penalized regression model fitted by glmnet via cross-validation, conforming to the calling convention required by argument \code{complexity} in \code{peperr} call. 
}
\usage{
complexity.glmnet(response, x, full.data, ...)
}
\arguments{
  \item{response}{a survival object (with \code{Surv(time, status)}, or a binary vector with entries 0 and 1).}
  \item{x}{\code{n*p} matrix of covariates.}
  \item{full.data}{data frame containing response and covariates of the full data set.}
  \item{\dots}{additional arguments passed to \code{cv.glmnet} call such as \code{family}.}
  }
\value{
Scalar value giving the optimal lambda.
}
\details{
Function is basically a wrapper for \code{cv.glmnet} of package \code{glmnet}. A n-fold cross-validation (default n=10) is performed to determine the optimal penalty lambda.
For Cox PH regression models the deviance based on penalized partial log-likelihood is used as loss function. For binary endpoints other loss functions are available as well (see \code{type.measure}). Deviance is default. Calling \code{peperr}, the default arguments of \code{cv.glmnet} can be changed by passing a named list containing these as argument \code{args.complexity}.
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

\seealso{\code{\link[peperr]{peperr}}, \code{\link[glmnet]{cv.glmnet}}}
\keyword{models} \keyword{penalized regression} 
