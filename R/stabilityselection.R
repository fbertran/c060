#' Stability path for glmnet models
#' 
#' The function calculates the stability path for glmnet models, e.g. the
#' selection probabilities of the features along the range of regularization
#' parameters.
#' 
#' 
#' @param y response variable. Like for the glment function: Quantitative for
#' \code{family="gaussian"} or \code{family="poisson"} (non-negative counts).
#' For \code{family="binomial"} should be either a factor with two levels, or a
#' two-column matrix of counts or proportions. For \code{family="multinomial"},
#' can be a \code{nc>=2} level factor, or a matrix with \code{nc} columns of
#' counts or proportions. For \code{family="cox"}, \code{y} should be a
#' two-column matrix with columns named 'time' and 'status'. The latter is a
#' binary variable, with '1' indicating death, and '0' indicating right
#' censored. The function \code{Surv()} in package \code{survival} produces
#' such a matrix
#' @param x input matrix. Like for the glmnet function: of dimension nobs x
#' nvars; each row is an observation vector. Can be in sparse matrix format
#' (inherit from class \code{"sparseMatrix"} as in package \code{Matrix}; not
#' yet available for \code{family="cox"})
#' @param size proportion of samples drawn in every subsample used for the
#' stability selection.
#' @param steps number of subsamples used for the stability selection.
#' @param weakness weakness parameter used for the randomised lasso as
#' described in Meinshausen and B\"uhlmann (2010).  For each subsample the
#' features are reweighted by a random weight uniformly sampled in
#' [weakness,1]. This additional randomisation leads to a more consistent
#' estimation of the stable set of features.
#' @param mc.cores number of cores used for the parallelization. For unix like
#' system the parallelization is done by forking using the function
#' \code{mclapply}. For windows systems socket cluster are used.
#' @param ...  further arguments that are passed to the \code{glmnet} function.
#' @return an object of class "stabpath", which is a list of three objects
#' \item{fit}{ the fit object of class "glmnet" as returned from the glmnet
#' function when applied to the complete data set.  } \item{stabpath}{ a matrix
#' which represents the stability path.  } \item{qs}{ a vector holding the
#' values of the average number of non-zero coefficients w.r.t to the lambdas
#' in the regularization path.  }
#' @author Martin Sill \email{m.sill@@dkfz.de}
#' @seealso
#' \code{\link[glmnet]{glmnet},\link[c060]{stabsel},\link[c060]{plot.stabpath}}
#' @references Meinshausen N. and B\"uhlmann P. (2010), \emph{Stability
#' Selection, Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology) Volume 72, Issue 4, pages 417--473.}\cr
#' 
#' Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060:
#' Extended Inference with Lasso and Elastic-Net Regularized Cox and
#' Generalized Linear Models, Journal of Statistical Software, Volume 62(5),
#' pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords stability selection
#' @examples
#' \dontrun{
#' #gaussian
#' set.seed(1234)
#' x <- matrix(rnorm(100*1000,0,1),100,1000)
#' y <- x[1:100,1:1000]%*% c(rep(2,5),rep(-2,5),rep(.1,990))
#' res <- stabpath(y,x,weakness=1,mc.cores=2)
#' plot(res)
#' 
#' #binomial
#' y=sample(1:2,100,replace=TRUE)
#' res <- stabpath(y,x,weakness=1,mc.cores=2,family="binomial")
#' plot(res)
#'     
#' #multinomial
#' y=sample(1:4,100,replace=TRUE)
#' res <- stabpath(y,x,weakness=1,mc.cores=2,family="multinomial")
#' plot(res)
#'     
#' #poisson
#' N=100; p=1000
#' nzc=5
#' x=matrix(rnorm(N*p),N,p)
#' beta=rnorm(nzc)
#' f = x[,seq(nzc)]%*%beta
#' mu=exp(f)
#' y=rpois(N,mu)
#' res <- stabpath(y,x,weakness=1,mc.cores=2,family="poisson")
#' plot(res)
#' 
#' #Cox
#' library(survival)
#' set.seed(10101)
#' N=100;p=1000
#' nzc=p/3
#' x=matrix(rnorm(N*p),N,p)
#' beta=rnorm(nzc)
#' fx=x[,seq(nzc)]%*%beta/3
#' hx=exp(fx)
#' ty=rexp(N,hx)
#' tcens=rbinom(n=N,prob=.3,size=1)
#' y=cbind(time=ty,status=1-tcens)
#' res <- stabpath(y,x,weakness=1,mc.cores=2,family="cox")
#' plot(res)
#' }
#' 
#' @export stabpath
stabpath <- function(y,x,size=0.632,steps=100,weakness=1,mc.cores=getOption("mc.cores", 2L),...){
  fit <- glmnet(x,y,...)
  if(is(fit[1],"multnet")|is(fit[1],"lognet")) y <- as.factor(y)
  #if(is(fit[1],"lognet")) y <- as.logical(y) 
  p <- ncol(x)
  #draw subsets
  subsets <- sapply(1:steps,function(v){sample(1:nrow(x),nrow(x)*size)})
  
  # parallel computing depending on OS
  # UNIX/Mac
  if (.Platform$OS.type!="windows") {
    res <- mclapply(1:steps,mc.cores=mc.cores,glmnet.subset,subsets,x,y,lambda=fit$lambda,weakness,p,...)
  } else {
    # Windows  
    cl  <- makePSOCKcluster(mc.cores)
    clusterExport(cl,c("glmnet","drop0"))
    res <- parLapply(cl, 1:steps,glmnet.subset,subsets,x,y,lambda=fit$lambda,weakness,p,...)
    stopCluster(cl)
  }
  
  #merging
  res <- res[unlist(lapply(lapply(res,dim),function(x) x[2]==dim(res[[1]])[2]))]
  x <- as.matrix(res[[1]])
  qmat <- matrix(ncol=ncol(res[[1]]),nrow=length(res))
  qmat[1,] <- colSums(as.matrix(res[[1]]))
  for(i in 2:length(res)){
    qmat[i,] <- colSums(as.matrix(res[[i]]))
    x <- x + as.matrix(res[[i]])
  }
  x <- x/length(res)
  qs <- colMeans(qmat)
  out <- list(fit=fit,stabpath=x,qs=qs)	
  class(out) <- "stabpath" 
  return(out)
}

#internal function used by lapply 
glmnet.subset <- function(index,subsets,x,y,lambda,weakness,p,...){
  if(length(dim(y))==2|is(y,"Surv")){
    glmnet(x[subsets[,index],],y[subsets[,index],],lambda=lambda
           ,penalty.factor= 1/runif(p,weakness,1),...)$beta!=0
  }else{
    if(is.factor(y)&length(levels(y))>2){
      temp <- glmnet(x[subsets[,index],],y[subsets[,index]],lambda=lambda
                     ,penalty.factor= 1/runif(p,weakness,1),...)[[2]]
      temp <- lapply(temp,as.matrix)
      Reduce("+",lapply(temp,function(x) x!=0))
      
    }	
    else{
      glmnet(x[subsets[,index],],y[subsets[,index]],lambda=lambda
             ,penalty.factor= 1/runif(p,weakness,1),...)$beta!=0
    }
  }	
}


#performs error control and returns estimated set of stable variables and corresponding lambda




#' function to estimate a stable set of variables
#' 
#' Given a desired type I error rate and a stability path calculated with
#' \code{stability.path} the function selects a stable set of variables.
#' 
#' 
#' @param x an object of class "stabpath" as returned by the function
#' \code{stabpath}.
#' @param error the desired type I error level w.r.t. to the chosen type I
#' error rate.
#' @param type The type I error rate used for controlling the number falsely
#' selected variables. If \code{type="pfer"} the per-family error rate is
#' controlled and \code{error} corresponds to the expected number of type I
#' errors. Selecting \code{type="pfer"} and \code{error} in the range of $0 >
#' \code{error} < 1$ will control the family-wise error rate, i.e. the
#' probability that at least one variable in the estimated stable set has been
#' falsely selected. If \code{type="pcer"} the per-comparison error rate is
#' controlled and \code{error} corresponds to the expected number of type I
#' errors divided by the number variables.
#' @param pi_thr the threshold used for the stability selection, should be in
#' the range of $0.5 > pi_thr < 1$.
#' @return a list of four objects \item{stable}{ a vector giving the positions
#' of the estimated stable variables } \item{lambda}{ the penalization
#' parameter used for the stability selection } \item{lpos}{ the position of
#' the penalization parameter in the regularization path } \item{error}{the
#' desired type I error level w.r.t. to the chosen type I error rate }
#' \item{type}{the type I error rate }
#' @author Martin Sill \ \email{m.sill@@dkfz.de}
#' @seealso \code{\link{plot.stabpath},\link{stabpath}}
#' @references Meinshausen N. and B\"uhlmann P. (2010), \emph{Stability
#' Selection, Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology) Volume 72, Issue 4, pages 417--473.}\cr
#' 
#' Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060:
#' Extended Inference with Lasso and Elastic-Net Regularized Cox and
#' Generalized Linear Models, Journal of Statistical Software, Volume 62(5),
#' pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords stability selection
#' @examples
#' 
#' \dontrun{
#' #gaussian
#' set.seed(1234)
#' x=matrix(rnorm(100*1000,0,1),100,1000)
#' y <- x[1:100,1:1000]%*%c(rep(2,5),rep(-2,5),rep(.1,990))
#' res <- stabpath(y,x,weakness=1,mc.cores=2)
#' stabsel(res,error=0.05,type="pfer")
#' }
#' @export stabsel
stabsel <- function(x,error=0.05,type=c("pfer","pcer"),pi_thr=0.6){
  if(pi_thr <= 0.5 | pi_thr >= 1) stop("pi_thr needs to be > 0.5 and < 1!")
  if(is(x$fit[1],"multnet")){
    p <- dim(x$fit$beta[[1]])[1]
  }else{
    p <- dim(x$fit$beta)[1]
  }
  type <- match.arg(type)
  switch(type,
         "pcer"={
           if(error>=1 | error<=0)stop("pcer needs to be > 0 and < 1!")
           qv <- ceiling(sqrt(error* p * (2*pi_thr-1)*p)) },
         "pfer"={
           qv <- ceiling(sqrt(error * (2*pi_thr-1)*p)) }
  )
  if(x$qs[length(x$qs)]<=qv){ lpos <- length(x$qs)
  }else{
    lpos <- which(x$qs>qv)[1]
  }
  if(!is.na(lpos)){stable <- which(x$x[,lpos]>=pi_thr)}else{
    stable <- NA
  }
  out <- list(stable=stable,lambda=x$fit$lambda[lpos],lpos=lpos,error=error,type=type)
  return(out)
}

#' @exportS3Method
print.stabpath <- function(x,...){
  cat(" stabilitypath","\n",
      dim(x$x)[1],"variables","\n",
      dim(x$x)[2],"lambdas","\n")
}

#plot penalization and stability path 


#' function to plot a stability path
#' 
#' Given a desired family-wise error rate (FWER) and a stability path
#' calculated with \code{stability.path} the function selects an stable set of
#' features and plots the stability path and the corresponding regularization
#' path.
#' 
#' @exportS3Method
#' @param x an object of class "stabpath" as returned by the function
#' \code{stabpath}.
#' @param error the desired type I error level w.r.t. to the chosen type I
#' error rate.
#' @param type The type I error rate used for controlling the number falsely
#' selected variables. If \code{type="pfer"} the per-family error rate is
#' controlled and \code{error} corresponds to the expected number of type I
#' errors. Selecting \code{type="pfer"} and \code{error} in the range of 0 >
#' \code{error} < 1 will control the family-wise error rate, i.e. the
#' probability that at least one variable in the estimated stable set has been
#' falsely selected. If \code{type="pcer"} the per-comparison error rate is
#' controlled and \code{error} corresponds to the expected number of type I
#' errors divided by the number variables.
#' @param pi_thr the threshold used for the stability selection, should be in
#' the range of 0.5 > pi_thr < 1.
#' @param xvar the variable used for the xaxis, e.g. for "lambda" the selection
#' probabilities are plotted along the log of the penalization parameters, for
#' "norm" along the L1-norm and for "dev" along the fraction of explained
#' deviance.
#' @param col.all the color used for the variables that are not in the
#' estimated stable set
#' @param col.sel the color used for the variables in the estimated stable set
#' @param ...  further arguments that are passed to matplot
#' @return a list of four objects \item{stable}{ a vector giving the positions
#' of the estimated stable variables } \item{lambda}{ the penalization
#' parameter used for the stability selection } \item{lpos}{ the position of
#' the penalization parameter in the regularization path } \item{error}{the
#' desired type I error level w.r.t. to the chosen type I error rate }
#' \item{type}{the type I error rate }
#' @author Martin Sill \ \email{m.sill@@dkfz.de}
#' @seealso \code{\link{stabsel},\link{stabpath}}
#' @references Meinshausen N. and Buehlmann P. (2010), Stability Selection,
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology)
#' Volume 72, Issue 4, pages 417-473.\cr Sill M., Hielscher T., Becker N. and
#' Zucknick M. (2014), \emph{c060: Extended Inference with Lasso and
#' Elastic-Net Regularized Cox and Generalized Linear Models, Journal of
#' Statistical Software, Volume 62(5), pages 1--22.} 
#' https://doi.org/10.18637/jss.v062.i05.
#' @keywords stability selection
#' @examples
#' 
#' \dontrun{
#' #gaussian
#' set.seed(1234)
#' x=matrix(rnorm(100*1000,0,1),100,1000)
#' y <- x[1:100,1:1000]%*%c(rep(2,5),rep(-2,5),rep(.1,990))
#' res <- stabpath(y,x,weakness=1,mc.cores=2)
#' plot(res,error=.5,type='pfer')
#' }
#' 
plot.stabpath <- function(x,error=0.05,type=c("pfer","pcer"),pi_thr=0.6,xvar=c("lambda", "norm", "dev")
                          , col.all="black", col.sel="red",...){
  sel <- stabsel(x,error,type,pi_thr)
  if(is(x$fit[1],"multnet")){
    beta = as.matrix(Reduce("+",x$fit$beta))
  }else{
    beta = as.matrix(x$fit$beta)
  }  
  p <- dim(beta)[1]
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = apply(abs(beta), 2, sum)
    iname = "L1 Norm"
  }, lambda = {
    index = log(x$fit$lambda)
    iname = expression(paste("log ",lambda))
  }, dev = {
    index = x$fit$dev
    iname = "Fraction Deviance Explained"
  })
  #}
  #stability path
  cols <- rep(col.all,p)
  cols[sel$stable] <- col.sel
  lwds <- rep(1,p)
  lwds[sel$stable] <- 2
  if(!is(x$fit[1],"multnet")){
    par(mfrow=c(2,1))
    matplot(y=t(beta), x=index
            ,type="l",col=cols,lwd=lwds,lty=1,ylab=expression(paste(hat(beta)[i]))
            ,xlab=iname,main="Penalization Path",cex.lab=1,cex.axis=1,las=1,...)
  }
  matplot(y=as.matrix(t(x$x)), x=index
          ,type="l",col=cols,lwd=lwds,lty=1,ylab=expression(paste(hat(Pi)))
          ,xlab=iname,main="Stability Path",ylim=c(0,1),cex.lab=1,cex.axis=1,las=1,...)
  abline(h=pi_thr,col="darkred",lwd=1,lty=1)
  abline(v=index[sel$lpos],col="darkred",lwd=1,lty=1)
  #text(x=20,y=0.9,paste(expression(paste(lambda)),"=",paste(round(sel[[2]],digits=3)),sep=""),cex=0.75)
  return(sel)
}
NULL

