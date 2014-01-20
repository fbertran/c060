\name{predictProb.glmnet}
\alias{predictProb.glmnet}
\title{Extract predicted survival probabilities from a glmnet fit}
\description{
Extracts predicted survival probabilities from survival model fitted by glmnet, providing an interface as required by \code{pmpec}. 
}
\usage{
\method{predictProb}{glmnet}(object, response, x, times, complexity, ...)
}
\arguments{
\item{object}{a fitted model of class \code{glmnet}.}
\item{response}{response variable. Quantitative for \code{family="gaussian"}, or \code{family="poisson"} (non-negative counts). For \code{family="binomial"} should be either a factor with two levels, or a two-column matrix of counts or proportions. For \code{family="multinomial"}, can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. }
\item{x}{\code{n*p} matrix of covariates.}
\item{times}{vector of evaluation time points.}
\item{complexity}{lambda penalty value.}
\item{\dots}{additional arguments, currently not used.}
}

\value{
Matrix with probabilities for each evaluation time point in \code{times} (columns) and each new observation (rows). 
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
   Bioinformatics, Vol. 25(6), 827-829.}
}

\author{
Thomas Hielscher \
\email{t.hielscher@dkfz.de}
}

\seealso{\code{\link[c060]{predictProb.coxnet}}, \code{\link[peperr]{peperr}}, \code{\link[glmnet]{glmnet}}}
\keyword{models} \keyword{penalized regression} \keyword{survival}