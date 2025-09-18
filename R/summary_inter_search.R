###########################################################################################################




#' Summary method for interval search models
#' 
#' Produces a summary of a fitted interval search model
#' 
#' 
#' @param object an object of class \code{intsearch} as returned by the
#' function \code{EPSGO}.
#' @param digits digits after the comma
#' @param verbose default set to TRUE.
#' @param first.n show first.n entries , default 5.
#' @param \dots additional argument(s)
#' @return A list of following elements \item{info}{ a data frame of four
#' objects for optimal models\cr \itemize{ \item alpha - a vector of alphas
#' \item lambda - a vector of penalization parameter lambda \item deviances - a
#' vector of deviances \item n.features -a vector of number of features
#' selected in each optimal model } } \item{opt.alpha}{ an optimal value for
#' alpha} \item{opt.lambda}{an optimal value for lambda} \item{opt.error}{ an
#' optimal value for error, hier minimal diviance} \item{opt.models}{ a list of
#' optimal models with the same optimal error}
#' @author Natalia Becker \ \email{natalia.becker@@dkfz.de}
#' @seealso \code{\link[c060]{EPSGO}},\code{\link[c060]{plot.sum.intsearch}}
#' @references Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, Journal of Statistical Software, Volume
#' 62(5), pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords summary
summary.intsearch<-function(object,digits = max(3, getOption("digits") - 3), verbose=TRUE, first.n=5, ...){
  fit <- object
  alphas <- fit$Xtrain[,1]
  lambdas <- unlist(sapply(sapply(fit$model, "[", "model"), "[", "lambda"))
  deviances <- fit$Ytrain
  # round problems!!!! take first from the fit 
  # number of selected features in the models; dfs
  tmp.models<-sapply(sapply(sapply(fit$model, "[", "model"), "[", "cvreg"), "[", "glmnet.fit")
  
  n.features<-mapply( function(List, lam) List$df[which(List$lambda %in% lam)], tmp.models, lambdas)
  
  # optimal models
  #print("chose the model with min num of FS ")      
  opt.models <- sapply(fit$model.list, "[", "model") [fit$Ytrain == fit$fmin ]
  
  opt.alpha <- opt.models[[1]]$alpha
  opt.lambda <- opt.models[[1]]$lambda
  opt.error <- fit$fmin 
  
  out <- list(info=data.frame(alpha=alphas,lambda=lambdas,deviance=deviances,n.features=n.features),
              opt.alpha=opt.alpha, opt.lambda=opt.lambda, opt.error=opt.error,
              opt.models=opt.models)
  class(out) <- "sum.intsearch"
  
  if(verbose){
    cat("Summary interval search \n\n")
    cat(paste("show the first", first.n,"out of",nrow(out$info),"entries\n"))
    print(out$info[1:first.n,])
    cat("\n..............................")
    
    cat("\n\n Optimal parameters found are: \n\n")
    cat(paste("alpha = ",round(out$opt.alpha,digits),
              "\t",
              "lambda = ",round(out$opt.lambda,digits),
              "deviance = ",round(out$opt.error,digits)))
  }
  invisible(out)
}
