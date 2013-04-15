\name{aggregation.auc}
\alias{aggregation.auc}
\title{Determine the area under the ROC curve for a fitted model}
\description{
Evaluate the area under the ROC curve for a fitted model on new data. To be used as argument \code{aggregation.fun} in \code{peperr} call. 
}
\usage{
aggregation.auc(full.data=NULL, response, x, model, cplx=NULL,  type=c("apparent", "noinf"), 
   fullsample.attr = NULL, ...)
}
\arguments{
\item{full.data}{passed from \code{peperr}, but not used for calculation.}
\item{response}{vector of binary response.}
\item{x}{\code{n*p} matrix of covariates.}
\item{model}{model fitted as returned by a \code{fit.fun}, as used in a call to \code{peperr}.}
\item{cplx}{passed from \code{peperr}, but not necessary for calculation.}
\item{type}{character.}
\item{fullsample.attr}{passed from \code{peperr}, but not necessary for calculation.}
\item{\dots}{additional arguments, passed to \code{predict} function.}
}
\details{
Area under the ROC curve is calculated based on internal \code{glmnet:::auc} function from package \code{glmnet}.

}
\value{
Scalar, indicating the area under the ROC curve.
}

\author{
Thomas Hielscher \
\email{t.hielscher@dkfz.de}
}

\seealso{\code{\link[peperr]{peperr}}}

\keyword{models} \keyword{regression} \keyword{classification}