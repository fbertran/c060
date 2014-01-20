\name{tune.glmnet.interval}

\alias{tune.glmnet.interval}

\title{  Wrapper function for \code{glmnet} objects.    }
\description{
   Wrapper function for \code{glmnet} objects used by \code{epsgo} function.
   This function is  mainly used within the function \code{\link{epsgo}}  
}



\usage{
tune.glmnet.interval(parms, x, y,
                     weights, 
                     offset = NULL, 
                     lambda = NULL, 
                     type.measure = c("mse", "deviance", "class", "auc", "mae"),
                     seed=12345, 
                     nfolds = 10, 
                     foldid=NULL, 
                     grouped = TRUE, 
                     type.min=c("lambda.min", "lambda.1se"),
                     family,
                     verbose=FALSE,
                     \dots)
}
%tune.glmnet.interval: no visible binding for global variable
%  'parms.coding'


\arguments{
  \item{parms}{tuning parameter alpha for \code{glmnet} object}
  \item{x,y}{x is a matrix where each row refers to a sample a each
    column refers to a gene; y is a factor which includes the class for
    each sample}
   \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
   \item{offset}{A vector of length nobs that is included in the linear predictor (a nobs x nc matrix for the "multinomial" family). Useful for the "poisson" family (e.g. log of exposure time), or for refining a model by starting at a current fit. Default is NULL. If supplied, then values must also be supplied to the predict function.}
   \item{lambda}{A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. Supplying a value of lambda overrides this. WARNING: use with care. Do not supply a single value for lambda (for predictions after CV use predict() instead). Supply instead a decreasing sequence of lambda values. glmnet relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.}
   \item{type.measure}{loss to use for cross-validation. Currently five options, not all available for all models. The default is type.measure="deviance", which uses squared-error for gaussian models (a.k.a type.measure="mse" there), deviance for logistic and poisson regression, and partial-likelihood for the Cox model. type.measure="class" applies to binomial and multinomial logistic regression only, and gives misclassification error. type.measure="auc" is for two-class logistic regression only, and gives area under the ROC curve. type.measure="mse" or type.measure="mae" (mean absolute error) can be used by all models except the "cox"; they measure the deviation from the fitted mean to the response.}
  \item{seed}{seed}
  \item{nfolds}{number of cross-validation's folds, default 10.}
  \item{foldid}{an optional vector of values between 1 and nfold identifying what fold each observation is in. If supplied, nfold can be missing.}
  \item{grouped}{ This is an experimental argument, with default TRUE, and can be ignored by most users. For all models except the "cox", this refers to computing nfolds separate statistics, and then using their mean and estimated standard error to describe the CV curve. If grouped=FALSE, an error matrix is built up at the observation level from the predictions from the nfold fits, and then summarized (does not apply to type.measure="auc"). For the "cox" family, grouped=TRUE obtains the CV partial likelihood for the Kth fold by subtraction; by subtracting the log partial likelihood evaluated on the full dataset from that evaluated on the on the (K-1)/K dataset. This makes more efficient use of risk sets. With grouped=FALSE the log partial likelihood is computed only on the Kth fold}
  \item{type.min}{parameter for chosing optimal model: 'lambda.min'- value of lambda that gives minimum mean cross-validated error (cvm).
 'lambda.1se' - largest value of lambda such that error is within one standard error of the minimum.}
  \item{family}{family of the model, i.e. cox, glm,...} 
 \item{verbose}{verbose}   
  \item{\dots}{Further parameters}
}



\value{
  \item{q.val }{minimal value of Q.func on the interval defined by bounds. Here, q.val is  minimum mean cross-validate               d error (cvm)}
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

\seealso{\code{\link{epsgo}}}

\keyword{models}
\keyword{multivariate}
\keyword{iteration}
\keyword{optimize}





