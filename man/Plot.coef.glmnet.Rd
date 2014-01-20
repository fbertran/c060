\name{Plot.coef.glmnet}
\alias{Plot.coef.glmnet}
\title{
function to highlight the path of a pre-specified set of variables within the coefficient path
}
\description{
Creates several plots showing the coefficient path for the final model of a cv.glmnet fit and highlights the path of a pre-specified set of variables within the coefficient path.
}
\usage{
Plot.coef.glmnet(cvfit, betas)
}
\arguments{
  \item{cvfit}{an object of class "cv.glmnet" as returned by the function \code{cv.glmnet}.}
  \item{betas}{a vector of names of variables; must be a subset of rownames(coef(cvfit)).}
}
\value{
a list of four objects
  \item{stable}{
  a vector giving the positions of the estimated stable variables 
  }
  \item{lambda}{
  the penalization parameter used for the stability selection 
  }
  \item{lpos}{
  the position of the penalization parameter in the regularization path
  }
  \item{error}{the desired type I error level w.r.t. to the chosen type I error rate
  }
   \item{type}{the type I error rate 
  }
}
\author{
Manuela Zucknick \
\email{m.zucknick@dkfz-heidelberg.de}
}
\examples{
\dontrun{
set.seed(1010)
n=1000;p=100
nzc=trunc(p/10)
x=matrix(rnorm(n*p),n,p)
beta=rnorm(nzc)
fx= x[,seq(nzc)] \%*\% beta
eps=rnorm(n)*5
y=drop(fx+eps)
px=exp(fx)
px=px/(1+px)
ly=rbinom(n=length(px),prob=px,size=1)
set.seed(1011)
cvob1=cv.glmnet(x,y)
plot.coef.glmnet(cvob1, c("V1","V100"))
}}
\keyword{coefficient path}