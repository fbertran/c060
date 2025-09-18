#' @keywords internal
#' @aliases c060-package c060 NULL
#' @author Thomas Hielscher \email{t.hielscher@@dkfz.de}, Manuela Zucknick, Natalia Becker, Frédéric Bertrand.
#' 
#' @seealso \code{\link[c060]{predictProb.coxnet}},
#' \code{\link[peperr]{peperr}}, \code{\link[glmnet]{glmnet}}
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#' \emph{Regularization Paths for Generalized Linear Models via Coordinate
#' Descent}, \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
#' \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}\cr
#' \url{https://www.jstatsoft.org/v33/i01/}\cr Simon, N., Friedman, J., Hastie,
#' T., Tibshirani, R. (2011) \emph{Regularization Paths for Cox's Proportional
#' Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol.
#' 39(5) 1-13}\cr \url{https://www.jstatsoft.org/v39/i05/}\cr Porzelius, C.,
#' Binder, H., and Schumacher, M. (2009) \emph{Parallelized prediction error
#' estimation for evaluation of high-dimensional models, Bioinformatics, Vol.
#' 25(6), 827-829.}\cr Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, Journal of Statistical Software, Volume
#' 62(5), pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords models penalized regression survival
#' 
"_PACKAGE"

#' @import glmnet survival pamr parallel peperr penalizedSVM lattice mlegp tgp
#' @importFrom penalized profL1
#' @importFrom penalized penalized
#' @importFrom methods is
#' @importFrom grDevices dev.off
#' @importFrom grDevices gray.colors
#' @importFrom grDevices pdf
#' @importFrom graphics abline
#' @importFrom graphics grid
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics matplot
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics polygon
#' @importFrom graphics text
#' @importFrom graphics axis
#' @importFrom stats approx
#' @importFrom stats coef
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats sd
#' 
NULL


