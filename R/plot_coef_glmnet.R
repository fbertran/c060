#' function to highlight the path of a pre-specified set of variables within
#' the coefficient path
#' 
#' Creates several plots showing the coefficient path for the final model of a
#' cv.glmnet fit and highlights the path of a pre-specified set of variables
#' within the coefficient path.
#' 
#' 
#' @param cvfit an object of class "cv.glmnet" as returned by the function
#' \code{cv.glmnet}.
#' @param betas a vector of names of variables; must be a subset of
#' rownames(coef(cvfit)).
#' @return a list of four objects \item{stable}{ a vector giving the positions
#' of the estimated stable variables } \item{lambda}{ the penalization
#' parameter used for the stability selection } \item{lpos}{ the position of
#' the penalization parameter in the regularization path } \item{error}{the
#' desired type I error level w.r.t. to the chosen type I error rate }
#' \item{type}{the type I error rate }
#' @author Manuela Zucknick \ \email{m.zucknick@@dkfz-heidelberg.de}
#' @references Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, Journal of Statistical Software, Volume
#' 62(5), pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords coefficient path
#' @examples
#' 
#' \dontrun{
#' set.seed(1010)
#' n=1000;p=100
#' nzc=trunc(p/10)
#' x=matrix(rnorm(n*p),n,p)
#' beta=rnorm(nzc)
#' fx= x[,seq(nzc)] %*% beta
#' eps=rnorm(n)*5
#' y=drop(fx+eps)
#' px=exp(fx)
#' px=px/(1+px)
#' ly=rbinom(n=length(px),prob=px,size=1)
#' set.seed(1011)
#' cvob1=cv.glmnet(x,y)
#' Plot.coef.glmnet(cvob1, c("V1","V100"))
#' }
#' @export Plot.coef.glmnet
Plot.coef.glmnet <- function(cvfit, betas){

op <- par(no.readonly = TRUE)
par(mar=c(4,4,2.5,1), mgp=c(2.5,1,0), mfrow=c(2,2))
fit <- cvfit$glmnet.fit
bet <- fit$beta[match(betas, rownames(fit$beta)),]

plot(fit, xvar="lambda", col="gray")
plotCoef(bet, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio, xvar = "lambda", add=TRUE, col="red")
abline(v=log(cvfit$lambda.min), lty=3)
abline(v=log(cvfit$lambda.1se), lty=3)

plotCoef(bet, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio, xvar = "lambda", add=FALSE, col="red")
abline(v=log(cvfit$lambda.min), lty=3)
abline(v=log(cvfit$lambda.1se), lty=3)

norm <- apply(abs(fit$beta), 2, sum)
plot(fit, xvar="norm", col="gray")
plotCoef(bet, xvar = "norm", add=TRUE, col="red",
                  norm = norm, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio)
abline(v=norm[match(cvfit$lambda.min, cvfit$lambda)], lty=3)
abline(v=norm[match(cvfit$lambda.1se, cvfit$lambda)], lty=3)

plot(fit, xvar="dev", col="gray")
plotCoef(bet, lambda = fit$lambda, df = fit$df, dev = fit$dev.ratio, xvar = "dev", add=TRUE, col="red")
abline(v=fit$dev.ratio[match(cvfit$lambda.min, cvfit$lambda)], lty=3)
abline(v=fit$dev.ratio[match(cvfit$lambda.1se, cvfit$lambda)], lty=3)

par(op)
}
