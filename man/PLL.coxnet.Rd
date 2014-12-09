\name{PLL.coxnet}
\alias{PLL.coxnet}
\title{Predictive partial log-likelihood for glmnet Cox PH model fit}
\description{
Extracts the predictive partial log-likelihood from a glmnet Cox PH model fit.
}
\usage{
\method{PLL}{coxnet}(object, newdata, newtime, newstatus, complexity, ...)
}
\arguments{
  \item{object}{fitted model of class \code{coxnet}.}
  \item{newdata}{\code{n_new*p} matrix of covariates.}
  \item{newtime}{\code{n_new}-vector of censored survival times.}
  \item{newstatus}{\code{n_new}-vector of survival status, coded with 0 and .1}
  \item{complexity}{lambda penalty value.}
  \item{\dots}{additional arguments, not used.}
}
\details{
Used by function \code{peperr}, if function \code{fit.glmnet} and \code{family="cox"} is used for model fit, which gives a class \code{coxnet} object.
This is basically a wrapper based on the \code{coxnet.deviance} function from package \code{glmnet}.
}
\value{
Vector of length \code{n_new}
}
\keyword{models} \keyword{penalized regression} \keyword{survival}

\author{
Thomas Hielscher \
\email{t.hielscher@dkfz.de}
}
\references{
Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox and Generalized Linear Models, Journal of Statistical Software, Volume 62(5), pages 1--22.}
\url{http://www.jstatsoft.org/v62/i05/}
}