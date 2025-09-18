###########################################################################################################
#'
#' Get coefficients for a model
#' 
#' Get coefficients for a model after applying interval search for tuning
#' parameters
#' 
#' 
#' @param object an object as returned by the function \code{summary.intsearch}.
#' @param \dots additional argument(s)
#' @return named vector of non-zero coeficients for the optimal lambda
#' @author Natalia Becker \ \email{natalia.becker@@dkfz.de}
#' @seealso \code{\link[c060]{EPSGO}},
#' \code{\link[c060]{summary.intsearch}},\code{\link[c060]{plot.sum.intsearch}}
#' @references Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, Journal of Statistical Software, Volume
#' 62(5), pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords system
coef.sum.intsearch<-function(object,...){
  # get coef for a object from fit object after running interval search 
   
  f1 <- object$cvreg
  res<-f1$glmnet.fit
  cof <- as.vector(coef(res, s=object$lambda))
  names.cof <- rownames(res$beta)
  cofn <- cof[which(cof != 0)]
  names(cofn) <- names.cof[which(cof != 0)] 
  bet <- res$beta[match(names(cofn), rownames(res$beta)),]
  
  return(cofn) 
}  

