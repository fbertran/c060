\name{coef.sum.intsearch}
\alias{coef.sum.intsearch}

\title{
Get coefficients for a model 
}
\description{
Get coefficients for a model after applying interval search for tuning parameters 
}
\usage{
\method{coef}{sum.intsearch}(object,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{ an object as returned by the function \code{summary.intsearch}.}
\item{\dots}{additional argument(s)}
}
\value{
 named vector of non-zero coeficients for the optimal lambda}
\author{
Natalia Becker  \
\email{natalia.becker@dkfz.de}
}

\seealso{ \code{\link{EPSGO}}, \code{\link{summary.intsearch}},\code{\link{plot.sum.intsearch}}
}

\keyword{system}
\references{
Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox and Generalized Linear Models, Journal of Statistical Software, Volume 62(5), pages 1--22.}
\url{http://www.jstatsoft.org/v62/i05/}
}