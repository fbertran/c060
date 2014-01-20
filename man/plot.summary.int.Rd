\name{plot.sum.intsearch}
\alias{plot.sum.intsearch}
\title{
Plot Summary object for interval search models
} 
\description{
Produces a plot for summary object of a fitted interval search model.
Plot 'visited' points against iteration steps. start.N points are initial points selected before interval search starts.
}
\usage{
\method{plot}{sum.intsearch}(x,type="summary",startN=21,... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 an object of class \code{sum.intsearch} as returned by the function \code{summary.intsearch}.
  }
  \item{type}{type of plot to be drawn, \code{type="summary"} will plot the partial log likelihood deviance as a function of   both tuning parameters \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} and log\ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}. The final solution will be highlighted by solid red line.
Alternativly, \code{type="points"} will draw the distribution of initial and visited points of the interval search plotted in chronological order.}
  \item{startN}{number of initial points. Needed if \code{type="points"}}
  \item{...}{additional argument(s)}
}


\author{
Natalia Becker  \
\email{natalia.becker@dkfz.de}
}

\seealso{ \code{\link{EPSGO}}, \code{\link{summary.intsearch}}
}

\keyword{plot}
