\name{balancedFolds }
\alias{balancedFolds }

\title{   Function producing stratified/ balanced folds for cross validation      }
\description{
  Get balanced folds for cross validation,  which are used for tuning penalization parameters
}

\usage{
balancedFolds(class.column.factor, cross.outer) 
}


\arguments{
  \item{class.column.factor}{class labels of length n  }
  \item{cross.outer}{  number of folds}
}



\value{
  \item{permutated.cut }{vector of length n, indicating the fold belongs to}
 \item{model }{ model list
    \itemize{
	  	\item alpha -  optimal alpha
		  \item lambda - optimal lambda
      \item nfolds - cross-validation's folds
		  \item cvreg -  \code{cv.glmnet} object for optimal alpha 
		  \item fit - \code{glmnet} object for optimal alpha and optimal lambda 
		} 
  }


} 
\author{Natalia Becker  natalia.becker at dkfz.de }

\seealso{\code{\link{EPSGO}}}

\keyword{models}
\keyword{multivariate}
\keyword{iteration}
\keyword{optimize}
\references{
Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox and Generalized Linear Models, Journal of Statistical Software, Volume 62(5), pages 1--22.}
\url{http://www.jstatsoft.org/v62/i05/}
}




