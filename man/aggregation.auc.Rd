\name{aggregation.auc}
\alias{aggregation.auc}
\title{Determine the area under the ROC curve for a fitted model}
\description{
Evaluate the area under the ROC curve for a fitted model on new data. To be used as argument \code{aggregation.fun} in \code{peperr} call. 
}
\usage{
aggregation.auc(full.data=NULL, response, x, model, cplx=NULL,  
type=c("apparent", "noinf"), fullsample.attr = NULL, \dots)
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

\examples{
\dontrun{
# binomial model - classification

library(c060)
library(gridExtra)
library(ggplot2)

set.seed(0815)
x <- matrix(rnorm(100*20),100,20)
y <- sample(0:1,100,replace=TRUE)

peperr_obj <- peperr(response=y, x=x, fit.fun=fit.glmnet, args.fit=list(family="binomial"),
           complexity=complexity.glmnet, args.complexity=list(nfolds=10, family="binomial"),
           trace=F, RNG="fixed",seed=0815,
#           aggregation.fun=c060:::aggregation.misclass,                  
#           aggregation.fun=c060:::aggregation.brier,                  
           aggregation.fun=c060:::aggregation.auc,                  
           indices=resample.indices(n=nrow(x), sample.n = 100, method = "sub632"))

tmp   <- data.frame(grp="",error=unlist(peperr_obj$sample.error)) 
errs  <- data.frame(error=c(perr(peperr_obj,"resample"),
         perr(peperr_obj,"632p"),perr(peperr_obj,"app"),
         perr(peperr_obj,"nullmodel")), col  = c("red","blue","green","brown"),
         row.names=c("mean\nout-of-bag",".632plus","apparent","null model"))
                 
p     <- ggplot(tmp, aes(grp,error))
pg    <- p + geom_boxplot(outlier.colour = rgb(0,0,0,0), outlier.size=0) +
         geom_jitter(position=position_jitter(width=.1)) + 
         theme_bw() + scale_y_continuous("AUC") +  scale_x_discrete("") +
         geom_hline(aes(yintercept=error, colour=col), data=errs, show_guide=T) + 
         scale_colour_identity("error type", guide = "legend", breaks=errs$col,
         labels=rownames(errs)) +
         ggtitle("AUC \n in bootstrap samples ")                       

p2     <- ggplot(data.frame(complx=peperr_obj$sample.complexity), aes(x=complx))
pg2    <- p2 + geom_histogram(binwidth = 0.02, fill = "white", colour="black") +
          theme_bw()+  xlab(expression(lambda)) +
          ylab("frequency") + 
          geom_vline(xintercept=peperr_obj$selected.complexity, colour="red") + 
          ggtitle("Selected complexity \n in bootstrap samples") +
          ggplot2::annotate("text", x = 0.12, y = -0.5,
          label = "full data", colour="red", size=4)

grid.arrange(pg2, pg, ncol=2)

}}

\keyword{models} \keyword{regression} \keyword{classification}