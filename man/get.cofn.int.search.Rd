\name{get.cofn.int.search}
\alias{get.cofn.int.search}

\title{
Get coefficients for a model 
}
\description{
Get coefficients for a model after applying interval search for tuning parameters 
}
\usage{
get.cofn.int.search(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{model}{ an object as returned by the function \code{summary.int.search}.}
}


\value{
 named vector of non-zero coeficients for the optimal lambda}

\author{
Natalia Becker  \
\email{natalia.becker@dkfz.de}
}

\seealso{ \code{\link{EPSGO}}, \code{\link{summary.intsearch}}
}

\keyword{system}
