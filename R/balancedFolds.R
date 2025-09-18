#' Function producing stratified/ balanced folds for cross validation
#' 
#' Get balanced folds for cross validation, which are used for tuning
#' penalization parameters
#' 
#' @export
#' @param class.column.factor class labels of length n
#' @param cross.outer number of folds
#' @return \item{permutated.cut }{vector of length n, indicating the fold
#' belongs to} \item{model }{ model list \itemize{ \item alpha - optimal alpha
#' \item lambda - optimal lambda \item nfolds - cross-validation's folds \item
#' cvreg - \code{cv.glmnet} object for optimal alpha \item fit - \code{glmnet}
#' object for optimal alpha and optimal lambda } }
#' @author Natalia Becker natalia.becker at dkfz.de
#' @seealso \code{\link[c060]{EPSGO}}
#' @references Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, Journal of Statistical Software, Volume
#' 62(5), pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords models multivariate iteration optimize
balancedFolds <- function(class.column.factor, cross.outer)
{
  #stolen from MCREstimate package
  # used for stratified(balanced) classification or regression
  # get balanced folds from pamr
  sampleOfFolds  <- get("balanced.folds",envir=asNamespace("pamr"))(class.column.factor, nfolds=cross.outer)
  permutated.cut <- rep(0,length(class.column.factor))
  for (sample in 1:cross.outer)
  {
    cat(sample,"\n")
    permutated.cut[sampleOfFolds[[sample]]] <- sample
  }
  return(permutated.cut)
}
