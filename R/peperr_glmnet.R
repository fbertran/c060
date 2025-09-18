###############################################################
# baseline survival/ hazard Breslow estimator
# function essentially based on gbm::basehaz.gbm
###############################################################
basesurv <- function (response, lp, times.eval = NULL, centered = FALSE)
{
    if (is.null(times.eval)) times.eval <- sort(unique(response[,1]))
    
    t.unique <- sort(unique(response[,1][response[,2] == 1]))
    alpha    <- length(t.unique)

    for (i in 1:length(t.unique)) {
        alpha[i] <- sum(response[,1][response[,2] == 1] == t.unique[i])/sum(exp(lp[response[,1] >=  t.unique[i]]))
    }

    obj   <- approx(t.unique, cumsum(alpha), yleft=0, xout = times.eval, rule=2)

    if (centered) obj$y <- obj$y * exp(mean(lp))
    obj$z <- exp(-obj$y)

    names(obj) <- c("times","cumBaseHaz","BaseSurv")
    return(obj)
}

###############################################
### wrapper for glmnet
###############################################





#' Interface function for fitting a penalized regression model with
#' \code{glmnet}
#' 
#' Interface for fitting penalized regression models for binary of survival
#' endpoint using \code{glmnet}, conforming to the requirements for argument
#' \code{fit.fun} in \code{peperr} call.
#' 
#' Function is basically a wrapper for \code{glmnet} of package \pkg{glmnet}.
#' Note that only penalized Cox PH (\code{family="cox"}) and logistic
#' regression models (\code{family="binomial"}) are sensible for prediction
#' error evaluation with package \code{peperr}.
#' 
#' @param response a survival object (with \code{Surv(time, status)}, or a
#' binary vector with entries 0 and 1).
#' @param x \code{n*p} matrix of covariates.
#' @param cplx lambda penalty value.
#' @param \dots additional arguments passed to \code{glmnet} call such as
#' \code{family}.
#' @return glmnet object
#' @author Thomas Hielscher \email{t.hielscher@@dkfz.de}
#' @seealso \code{\link[peperr]{peperr}}, \code{\link[glmnet]{glmnet}}
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
#' @keywords models penalized regression
#' @export fit.glmnet
fit.glmnet <- function (response, x, cplx, ...) 
{
    #require(glmnet)
    res <- NULL
    tryerr <- try(res <- glmnet(y = response, x = data.matrix(x), lambda = cplx,  ...), silent=TRUE)

    if(!is(tryerr, 'try-error') && is(res,"coxnet")) {
          res$linear.predictor  <- as.numeric(predict(res, newx=data.matrix(x), type="link"))
          res$response          <- response
    }
    class(res) <- class(res)[1]
    res
}





#' Interface for determination of penalty lambda in penalized regression model
#' via cross-validation
#' 
#' Determines the amount of shrinkage for a penalized regression model fitted
#' by glmnet via cross-validation, conforming to the calling convention
#' required by argument \code{complexity} in \code{peperr} call.
#' 
#' Function is basically a wrapper for \code{cv.glmnet} of package
#' \code{glmnet}. A n-fold cross-validation (default n=10) is performed to
#' determine the optimal penalty lambda. For Cox PH regression models the
#' deviance based on penalized partial log-likelihood is used as loss function.
#' For binary endpoints other loss functions are available as well (see
#' \code{type.measure}). Deviance is default. Calling \code{peperr}, the
#' default arguments of \code{cv.glmnet} can be changed by passing a named list
#' containing these as argument \code{args.complexity}. Note that only
#' penalized Cox PH (\code{family="cox"}) and logistic regression models
#' (\code{family="binomial"}) are sensible for prediction error evaluation with
#' package \code{peperr}.
#' 
#' @param response a survival object (with \code{Surv(time, status)}, or a
#' binary vector with entries 0 and 1).
#' @param x \code{n*p} matrix of covariates.
#' @param full.data data frame containing response and covariates of the full
#' data set.
#' @param \dots additional arguments passed to \code{cv.glmnet} call such as
#' \code{family}.
#' @return Scalar value giving the optimal lambda.
#' @author Thomas Hielscher \email{t.hielscher@@dkfz.de}
#' @seealso \code{\link[peperr]{peperr}}, \code{\link[glmnet]{cv.glmnet}}
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
#' @keywords models penalized regression
#' @export complexity.glmnet
complexity.glmnet <- function (response, x, full.data, ...) 
{
    #require(glmnet)
    lambda <- NULL
    tryerr <- try(cv <- cv.glmnet(y = response, x = data.matrix(x),  ...), silent=TRUE)
    
    if(!is(tryerr, 'try-error')){
      lambda <-cv$lambda.min
    }    
    lambda
}





#' Extract predicted survival probabilities from a glmnet fit
#' 
#' Extracts predicted survival probabilities from survival model fitted by
#' glmnet, providing an interface as required by \code{pmpec}.
#' 
#' @aliases predictProb.glmnet
#' @param object a fitted model of class \code{glmnet}
#' @param response a two-column matrix with columns named 'time' and 'status'.
#' The latter is a binary variable, with '1' indicating death, and '0'
#' indicating right censored. The function \code{Surv()} in package survival
#' produces such a matrix
#' @param x \code{n*p} matrix of covariates.
#' @param times vector of evaluation time points.
#' @param complexity lambda penalty value.
#' @param \dots additional arguments, currently not used.
#' @return Matrix with probabilities for each evaluation time point in
#' \code{times} (columns) and each new observation (rows).
#' @author Thomas Hielscher \email{t.hielscher@@dkfz.de}
#' @seealso
#' \code{\link[c060]{predictProb.glmnet}},\code{\link[peperr]{peperr}},
#' \code{\link[glmnet]{glmnet}}
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
predictProb.coxnet <- predictProb.glmnet <- function (object, response, x, times, complexity,  ...) 
{
    #require(glmnet)    
    lp       <- as.numeric(predict(object, newx=data.matrix(x),s=complexity, type="link"))
    basesurv <- basesurv(object$response,object$linear.predictor, sort(unique(times)))
    p        <- exp(exp(lp) %*% -t(basesurv$cumBaseHaz))
    
    if (NROW(p) != NROW(x) || NCOL(p) != length(times)) 
        stop("Prediction failed")
    p
}





#' Predictive partial log-likelihood for glmnet Cox PH model fit
#' 
#' Extracts the predictive partial log-likelihood from a glmnet Cox PH model
#' fit.
#' 
#' Used by function \code{peperr}, if function \code{fit.glmnet} and
#' \code{family="cox"} is used for model fit, which gives a class \code{coxnet}
#' object. This is basically a wrapper based on the \code{coxnet.deviance}
#' function from package \code{glmnet}.
#' 
#' @export PLL.coxnet
#' @param object fitted model of class \code{coxnet}.
#' @param newdata \code{n_new*p} matrix of covariates.
#' @param newtime \code{n_new}-vector of censored survival times.
#' @param newstatus \code{n_new}-vector of survival status, coded with 0 and .1
#' @param complexity lambda penalty value.
#' @param \dots additional arguments, not used.
#' @return Vector of length \code{n_new}
#' @author Thomas Hielscher \email{t.hielscher@@dkfz.de}
#' @references Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, \emph{Journal of Statistical Software}, Volume
#' 62(5), pages 1--22. https://doi.org/10.18637/jss.v062.i05.
#' @keywords models penalized regression survival
PLL.coxnet <- function(object, newdata, newtime, newstatus, complexity, ...) 
{
   #require(glmnet)
   PLL <- glmnet::coxnet.deviance(pred = NULL, Surv(newtime,newstatus), x = data.matrix(newdata), offset = NULL, weights = NULL, beta = coef(object,s=complexity)) 
   PLL / -2
}


##################################################
###### classification aggregation functions   ####
##################################################

aggregation.misclass <- function (full.data = NULL, response, x, model, cplx = NULL, 
    type = c("apparent", "noinf"), fullsample.attr = NULL, ...) 
{
    data <- as.data.frame(x)
    data$response <- response

    if ("glmnet" %in% class(model)) {
        probs <- as.numeric(predict(model, newx = data.matrix(x), type="response", ...))
    }
    else if (is(model[1],"penfit")) {
      probs <- predict(model, data = data, penalized = x, ...)
    }
    else if (is(model[1],"glm")) {
        probs <- predict(model, newdata = data, type="response", ...)
    }
    else {
        probs <- predict(model, data = data, type = "response", 
            ...)
    }
    type <- match.arg(type)
    if (type == "apparent") {
        mr <- sum(abs(round(probs) - response))/length(response)
    }
    if (type == "noinf") {
        mr <- mean(abs((matrix(response, length(response), length(response), 
            byrow = TRUE) - round(probs))))
    }
    mr
}

aggregation.brier <- function (full.data = NULL, response, x, model, cplx = NULL, 
    type = c("apparent", "noinf"), fullsample.attr = NULL, ...) 
{
    data          <- as.data.frame(x)
    data$response <- response

    if ("glmnet" %in% class(model)) {
        probs <- as.numeric(predict(model, newx = data.matrix(x), type="response", ...))
    }
    else if (is(model[1],"penfit")) {
      probs <- predict(model, data = data, penalized = x, ...)
    }
    else if (is(model[1],"glm")) {
        probs <- predict(model, newdata = data, type="response", ...)
    }    
    else {
        probs <- predict(model, data = data, type = "response", 
            ...)
    }
    type <- match.arg(type)
    if (type == "apparent") {
        brier.score <- sum((probs - response)^2)/length(response)
    }
    if (type == "noinf") {
        brier.score <- mean((matrix(response, length(response), 
            length(response), byrow = TRUE) - probs)^2)
    }
    brier.score
}





#' Determine the area under the ROC curve for a fitted model
#' 
#' Evaluate the area under the ROC curve for a fitted model on new data. To be
#' used as argument \code{aggregation.fun} in \code{peperr} call.
#' 
#' Area under the ROC curve is calculated based on internal \code{glmnet:::auc}
#' function from package \code{glmnet}.
#' 
#' @param full.data passed from \code{peperr}, but not used for calculation.
#' @param response vector of binary response.
#' @param x \code{n*p} matrix of covariates.
#' @param model model fitted as returned by a \code{fit.fun}, as used in a call
#' to \code{peperr}.
#' @param cplx passed from \code{peperr}, but not necessary for calculation.
#' @param type character.
#' @param fullsample.attr passed from \code{peperr}, but not necessary for
#' calculation.
#' @param \dots additional arguments, passed to \code{predict} function.
#' @return Scalar, indicating the area under the ROC curve.
#' @author Thomas Hielscher \email{t.hielscher@@dkfz.de}
#' @seealso \code{\link[peperr]{peperr}}
#' @keywords models regression classification
#' @examples
#' 
#' \dontrun{
#' # binomial model - classification
#' 
#' library(c060)
#' library(gridExtra)
#' library(ggplot2)
#' 
#' set.seed(0815)
#' x <- matrix(rnorm(100*20),100,20)
#' y <- sample(0:1,100,replace=TRUE)
#' 
#' peperr_obj <- peperr(response=y, x=x, fit.fun=fit.glmnet, args.fit=list(family="binomial"),
#'            complexity=complexity.glmnet, args.complexity=list(nfolds=10, family="binomial"),
#'            trace=F, RNG="fixed",seed=0815,
#' #           aggregation.fun=c060:::aggregation.misclass,                  
#' #           aggregation.fun=c060:::aggregation.brier,                  
#'            aggregation.fun=c060:::aggregation.auc,                  
#'            indices=resample.indices(n=nrow(x), sample.n = 100, method = "sub632"))
#' 
#' tmp   <- data.frame(grp="",error=unlist(peperr_obj$sample.error)) 
#' errs  <- data.frame(error=c(perr(peperr_obj,"resample"),
#'          perr(peperr_obj,"632p"),perr(peperr_obj,"app"),
#'          perr(peperr_obj,"nullmodel")), col  = c("red","blue","green","brown"),
#'          row.names=c("mean\nout-of-bag",".632plus","apparent","null model"))
#'                  
#' p     <- ggplot(tmp, aes(grp,error))
#' pg    <- p + geom_boxplot(outlier.colour = rgb(0,0,0,0), outlier.size=0) +
#'          geom_jitter(position=position_jitter(width=.1)) + 
#'          theme_bw() + scale_y_continuous("AUC") +  scale_x_discrete("") +
#'          geom_hline(aes(yintercept=error, colour=col), data=errs, show_guide=T) + 
#'          scale_colour_identity("error type", guide = "legend", breaks=errs$col,
#'          labels=rownames(errs)) +
#'          ggtitle("AUC \n in bootstrap samples ")                       
#' 
#' p2     <- ggplot(data.frame(complx=peperr_obj$sample.complexity), aes(x=complx))
#' pg2    <- p2 + geom_histogram(binwidth = 0.02, fill = "white", colour="black") +
#'           theme_bw()+  xlab(expression(lambda)) +
#'           ylab("frequency") + 
#'           geom_vline(xintercept=peperr_obj$selected.complexity, colour="red") + 
#'           ggtitle("Selected complexity \n in bootstrap samples") +
#'           ggplot2::annotate("text", x = 0.12, y = -0.5,
#'           label = "full data", colour="red", size=4)
#' 
#' grid.arrange(pg2, pg, ncol=2)
#' 
#' }
#' @export aggregation.auc
aggregation.auc <- function (full.data = NULL, response, x, model, cplx = NULL, 
    type = c("apparent", "noinf"), fullsample.attr = NULL, ...) 
{
    data <- as.data.frame(x)
    data$response <- response
    if ("glmnet" %in% class(model)) {
        probs <- as.numeric(predict(model, newx = data.matrix(x), type="response", ...))
    }
    else if (is(model[1],"penfit")) {
      probs <- predict(model, data = data, penalized = x, ...)
    }    
    else if (is(model[1],"glm")) {
        probs <- predict(model, newdata = data, type="response", ...)
    }    
    else {
        probs <- predict(model, data = data, type = "response", 
            ...)
    }
    type <- match.arg(type)
    if (type == "apparent") {
        auc <- auc(response,probs)
    }
    if (type == "noinf") {
        resp.mat <- matrix(response, length(response),  length(response), byrow = TRUE)
        auc      <- mean(apply(resp.mat, 1, function(d) auc(d,probs)))
    }
    auc
}

########################
### plot pecs        ###
########################





#' Plot method for prediction error curves of a peperr object
#' 
#' Plots individual and aggregated prediction error estimates based on
#' bootstrap samples.
#' 
#' This function is literally taken from \code{plot.peperr} in the
#' \code{peperr} package. The display of prediction error curves is adapted to
#' allow for numbers at risk and pointwise bootstrap quantiles.
#' 
#' @param x \code{peperr} object.
#' @param at.risk number at risk to be display. default is TRUE.
#' @param allErrors Display .632, no information and average out-of-bag error
#' in addition. default is FALSE.
#' @param bootRuns Display individual out-of-bag bootstrap samples. default is
#' FALSE.
#' @param bootQuants Display pointwise out-of-bag bootstrap quantiles as shaded
#' area. default is TRUE.
#' @param bootQuants.level Quantile probabilities for pointwise out-of-bag
#' bootstrap quantiles. default is 0.95, i.e. 2.5\% and 97.5\% quantiles.
#' @param leg.cex size of legend text
#' @param \dots additional arguments, not used.
#' @author Thomas Hielscher \email{t.hielscher@@dkfz.de}
#' @seealso \code{\link[peperr]{peperr}}
#' @references Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, Journal of Statistical Software, Volume
#' 62(5), pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords models regression survival
#' @examples
#' 
#' \dontrun{
#' 
#' # example from glmnet package
#' set.seed(10101)
#' library(glmnet)
#' library(survival)
#' library(peperr)
#' 
#' N=1000;p=30
#' nzc=p/3
#' x=matrix(rnorm(N*p),N,p)
#' beta=rnorm(nzc)
#' fx=x[,seq(nzc)]%*%beta/3
#' hx=exp(fx)
#' ty=rexp(N,hx)
#' tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator
#' y=Surv(ty,1-tcens)
#' 
#' peperr.object <- peperr(response=y, x=x, 
#'                         fit.fun=fit.glmnet, args.fit=list(family="cox"), 
#'                         complexity=complexity.glmnet,  
#'                         args.complexity=list(family="cox",nfolds=10),
#'                         indices=resample.indices(n=N, method="sub632", sample.n=10))
#' 
#' # pointwise bootstrap quantiles and all error types
#' Plot.peperr.curves(peperr.object, allErrors=TRUE)
#' 
#' # individual bootstrap runs and selected error types
#' Plot.peperr.curves(peperr.object, allErrors=FALSE, bootRuns=TRUE)
#' }
#' 
#' @export Plot.peperr.curves
Plot.peperr.curves <- function(x,at.risk=TRUE,allErrors=FALSE,  bootRuns=FALSE, bootQuants=TRUE, bootQuants.level=0.95, leg.cex=0.7, ...) {
  
  #require(peperr)

  if (bootRuns) bootQuants <- FALSE
  
  plot(x$attribute, x$null.model, type = "n", las=1,
       col = "blue", xlab = "Evaluation time points", 
       ylab = "Prediction error", main = "Prediction error curves", 
       ylim = c(0, max(perr(x), x$full.apparent, x$null.model) + 0.1))
  
  if (length(x$sample.error) > 1 & bootRuns==TRUE) {
    for (i in 1:(length(x$sample.error))) {
      lines(x$attribute, x$sample.error[[i]], type = "l", col = "light grey", lty = 1)
    }
  }

  if (length(x$sample.error) > 1 & bootQuants==TRUE) {
    boots  <- do.call("rbind",x$sample.error)
    quants <- apply(boots, 2, function(d) quantile(d, probs=c((1-bootQuants.level)/2,1 - (1-bootQuants.level)/2)))
    polygon(c(x$attribute,rev(x$attribute)),c(quants[1,],rev(quants[2,])), col="light grey", border="light grey")
  }
  
  if (allErrors==FALSE) {
     lines(x$attribute, x$null.model, type = "l", col = "blue", lwd = 2, lty = 1)
     lines(x$attribute, perr(x, "632p"), type = "l", col= "black", lty = 1, lwd = 2)
     lines(x$attribute, x$full.apparent, type = "l", col = "red", lty = 1, lwd = 2)
     if (bootRuns==TRUE) {
       legend(x = "topleft", col = c("blue", "black", "red", "light grey"), lwd=c(2,2,2,1), cex=leg.cex,
            lty = 1, legend = c("Null model", ".632+ estimate", "Full apparent", "Bootstrap samples"))
     } else {
       legend(x = "topleft", col = c("blue", "black", "red"), lwd=c(2,2,2), cex=leg.cex,
              lty = 1, legend = c("Null model", ".632+ estimate", "Full apparent"))
     }   
  }

  if (allErrors==TRUE) {
    lines(x$attribute, x$null.model, type = "l", col = "blue", lwd = 2, lty = 1)
    lines(x$attribute, perr(x, "632p"), type = "l", col= "black", lty = 1, lwd = 2)
    lines(x$attribute, perr(x, "632"), type = "l", col= "brown", lty = 1, lwd = 2)
    lines(x$attribute, perr(x, "NoInf"), type = "l", col= "green", lty = 1, lwd = 2)
    lines(x$attribute, perr(x, "resample"), type = "l", col= "dark grey", lty = 1, lwd = 2)
    lines(x$attribute, x$full.apparent, type = "l", col = "red", lty = 1, lwd = 2)
    if (bootRuns==TRUE) {
      legend(x = "topleft", ncol=2, col = c("blue", "black","brown","green","dark grey","red", "light grey"), lwd=c(2,2,2,2,2,2,1), cex=leg.cex,
           lty = 1, legend = c("Null model", ".632+ estimate",".632 estimate", "No Information","Out-of-bag average","Full apparent", "Bootstrap samples"))
    } else {
      legend(x = "topleft", ncol=2, col = c("blue", "black","brown","green","dark grey","red"), lwd=c(2,2,2,2,2,2), cex=leg.cex,
             lty = 1, legend = c("Null model", ".632+ estimate",".632 estimate", "No Information","Out-of-bag average","Full apparent"))
    }
  }
  
  if (at.risk) {
     tmpxaxp   <- par("xaxp")
     tmpusr    <- par("usr")
     at.loc    <- seq(tmpxaxp[1],tmpxaxp[2],length=tmpxaxp[3]+1)
     n.at.risk <- summary(survfit(x$response ~ 1),times=at.loc)$n.risk     

     text(x=at.loc, y=tmpusr[3], labels=n.at.risk, cex=0.8, pos=3)
     text(x=tmpxaxp[2]+(tmpusr[2]-tmpxaxp[2])/2, y=tmpusr[3], labels="at\nrisk", cex=0.8, pos=3)
  }  
}

