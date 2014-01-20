\name{epsgo}

\alias{epsgo}

\title{  Efficient Parameter Selection via Global Optimization  }
\description{
   Finds an optimal solution for the Q.func function.
}

\usage{
epsgo(Q.func, bounds,  round.n=5, parms.coding="none",
  fminlower=0, flag.find.one.min =FALSE,
  show=c("none", "final", "all"), N= NULL, maxevals = 500,
  pdf.name=NULL,  pdf.width=12,  pdf.height=12, my.mfrow=c(1,1),
  verbose=TRUE, seed=123,  \dots  )
}

\arguments{
  \item{Q.func}{ name of the function to be  minimized. }
  \item{bounds}{ bounds for parameters}
  \item{round.n}{ number of digits after comma, default: 5}           
  \item{parms.coding}{ parmeters coding: none  or log2, default: none.  }
  \item{fminlower}{ minimal value for the function Q.func, default is 0.     }
  \item{flag.find.one.min}{  do you want to find one min value and stop? Default: FALSE }
  \item{show}{  show plots of  DIRECT algorithm:    none, final iteration, all iterations. Default: none  }
  \item{N}{ define the number of start points, see details. }
  \item{maxevals}{ the maximum number of DIRECT function evaluations, default: 500.   }
  \item{pdf.name}{pdf name      }   
  \item{pdf.width}{ default: 12 }
  \item{pdf.height}{ default: 12 }
  \item{my.mfrow}{ default: c(1,1) }
  \item{verbose}{ verbose? default: TRUE. }
  \item{seed}{ seed }
  \item{\dots}{ additional argument(s) }
}





\value{
  \item{fmin }{minimal value of Q.func on the interval defined by bounds. }
  \item{xmin }{corresponding parameters for the minimum}
  \item{iter }{number of iterations}
  \item{neval }{  number of visited points }
  \item{maxevals }{  the maximum number of DIRECT function evaluations }
  \item{seed }{  seed}
  \item{bounds}{ bounds for parameters}
  \item{Q.func }{  name of the function to be  minimized. }
  \item{points.fmin }{  the set of points with the same fmin }
  \item{Xtrain }{  visited points }
  \item{Ytrain }{  the output of Q.func at visited points Xtrain }
  \item{gp.seed }{ seed for Gaussian Process }
  \item{model.list }{ detailed information of the search process }
} 


\details{
  if the number of start points (N)  is not defined by the user, it will be defined dependent on the dimensionality of the parameter space.
	N=10D+1, where  D is the number of parameters, but for high dimensional parameter space with more than 6 dimensions,  
  the initial set is restricted to 65. However for one-dimensional parameter space the N is set to 21 due to stability reasons.
 
   The idea of EPSGO (Efficient Parameter Selection via Global Optimization): Beginning
	from an intial Latin hypercube sampling containing N starting points we train
	an Online GP, look for the point with the maximal expected 	improvement, sample there and update the Gaussian Process(GP). Thereby
	it is not so important that GP really correctly 	models the error surface of the SVM in parameter space, but
	that it can give a us information about potentially interesting 	points in parameter space where we should sample next.
	We continue with sampling points until some convergence criterion is met.
  
  DIRECT is a sampling algorithm which requires no knowledge of the objective function gradient.
	Instead, the algorithm samples points in the domain, and uses the information it has obtained to decide where to
	search next. The DIRECT algorithm will globally converge to the maximal value of the objective function. The name
	DIRECT comes from the shortening of the phrase 'DIviding RECTangles', which describes the way the algorithm moves
	towards the optimum.  

	The code source was adopted from MATLAB originals, special thanks to Holger Froehlich.
}

	

\author{ Natalia Becker natalia.becker at dkfz.de }

\references{
Froehlich, H. and Zell, A. (2005) "Effcient parameter selection for support vector
machines in classification and regression via model-based global optimization"
\emph{In Proc. Int. Joint Conf. Neural Networks,  1431-1438 }.
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

 
nfolds = 10
set.seed(1234)
foldid <- balancedFolds(class.column.factor=y.classes, cross.outer=nfolds)

# y - binomial
y.classes<-ifelse(y>= median(y),1, 0)
bounds <- t(data.frame(alpha=c(0, 1)))
colnames(bounds)<-c("lower","upper")
 
fit <- epsgo(Q.func="tune.glmnet.interval", 
             bounds=bounds, 
             parms.coding="none", 
             seed = 1234, 
             show="none",
             fminlower = -100,
             x = x, y = y.classes, family = "binomial", 
             foldid = foldid,
             type.min = "lambda.1se",
             type.measure = "mse")
summary(fit)

# y - multinomial: low - low 25\%, middle - (25,75)-quantiles, high - larger 75\%.
y.classes<-ifelse(y <= quantile(y,0.25),1, ifelse(y >= quantile(y,0.75),3, 2))
bounds <- t(data.frame(alpha=c(0, 1)))
colnames(bounds)<-c("lower","upper")
 
fit <- epsgo(Q.func="tune.glmnet.interval", 
             bounds=bounds, 
             parms.coding="none", 
             seed = 1234, 
             show="none",
             fminlower = -100,
             x = x, y = y.classes, family = "multinomial", 
             foldid = foldid,
             type.min = "lambda.1se",
             type.measure = "mse")
summary(fit)

##poisson
N=500; p=20
nzc=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
f = x[,seq(nzc)]%*%beta
mu=exp(f)
y.classes=rpois(N,mu)

nfolds = 10
set.seed(1234)
foldid <- balancedFolds(class.column.factor=y.classes, cross.outer=nfolds)


fit <- epsgo(Q.func="tune.glmnet.interval", 
             bounds=bounds, 
             parms.coding="none", 
             seed = 1234, 
             show="none",
             fminlower = -100,
             x = x, y = y.classes, family = "poisson", 
             foldid = foldid,
             type.min = "lambda.1se",
             type.measure = "mse")
summary(fit)

#gaussian
set.seed(1234)
x=matrix(rnorm(100*1000,0,1),100,1000)
y <- x[1:100,1:1000]\%*\%c(rep(2,5),rep(-2,5),rep(.1,990))

foldid <- rep(1:10,each=10)

fit <- epsgo(Q.func="tune.glmnet.interval", 
             bounds=bounds, 
             parms.coding="none", 
             seed = 1234, 
             show="none",
             fminlower = -100,
             x = x, y = y, family = "gaussian", 
             foldid = foldid,
             type.min = "lambda.1se",
             type.measure = "mse")
summary(fit)  

# y - cox in vingette
}
}

\keyword{ models }
\keyword{ multivariate }
\keyword{ graphs }
\keyword{ iteration }
\keyword{ optimize }



