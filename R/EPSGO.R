###########################################################################################################
###########################################################################################################




#' Efficient Parameter Selection via Global Optimization
#' 
#' Finds an optimal solution for the Q.func function.
#' 
#' if the number of start points (N) is not defined by the user, it will be
#' defined dependent on the dimensionality of the parameter space.  N=10D+1,
#' where D is the number of parameters, but for high dimensional parameter
#' space with more than 6 dimensions, the initial set is restricted to 65.
#' However for one-dimensional parameter space the N is set to 21 due to
#' stability reasons.
#' 
#' The idea of EPSGO (Efficient Parameter Selection via Global Optimization):
#' Beginning from an intial Latin hypercube sampling containing N starting
#' points we train an Online GP, look for the point with the maximal expected
#' improvement, sample there and update the Gaussian Process(GP). Thereby it is
#' not so important that GP really correctly models the error surface of the
#' SVM in parameter space, but that it can give a us information about
#' potentially interesting points in parameter space where we should sample
#' next.  We continue with sampling points until some convergence criterion is
#' met.
#' 
#' DIRECT is a sampling algorithm which requires no knowledge of the objective
#' function gradient.  Instead, the algorithm samples points in the domain, and
#' uses the information it has obtained to decide where to search next. The
#' DIRECT algorithm will globally converge to the maximal value of the
#' objective function. The name DIRECT comes from the shortening of the phrase
#' 'DIviding RECTangles', which describes the way the algorithm moves towards
#' the optimum.
#' 
#' The code source was adopted from MATLAB originals, special thanks to Holger
#' Froehlich.
#' 
#' @param Q.func name of the function to be minimized.
#' @param bounds bounds for parameters
#' @param round.n number of digits after comma, default: 5
#' @param parms.coding parmeters coding: none or log2, default: none.
#' @param fminlower minimal value for the function Q.func, default is 0.
#' @param flag.find.one.min do you want to find one min value and stop?
#' Default: FALSE
#' @param show show plots of DIRECT algorithm: none, final iteration, all
#' iterations. Default: none
#' @param N define the number of start points, see details.
#' @param maxevals the maximum number of DIRECT function evaluations, default:
#' 500.
#' @param pdf.name pdf name
#' @param pdf.width default: 12
#' @param pdf.height default: 12
#' @param my.mfrow default: c(1,1)
#' @param verbose verbose? default: TRUE.
#' @param seed seed
#' @param \dots additional argument(s)
#' @return \item{fmin }{minimal value of Q.func on the interval defined by
#' bounds. } \item{xmin }{corresponding parameters for the minimum} \item{iter
#' }{number of iterations} \item{neval }{ number of visited points }
#' \item{maxevals }{ the maximum number of DIRECT function evaluations }
#' \item{seed }{ seed} \item{bounds}{ bounds for parameters} \item{Q.func }{
#' name of the function to be minimized. } \item{points.fmin }{ the set of
#' points with the same fmin } \item{Xtrain }{ visited points } \item{Ytrain }{
#' the output of Q.func at visited points Xtrain } \item{gp.seed }{ seed for
#' Gaussian Process } \item{model.list }{ detailed information of the search
#' process }
#' @author Natalia Becker natalia.becker at dkfz.de
#' @references Froehlich, H. and Zell, A. (2005) "Effcient parameter selection
#' for support vector machines in classification and regression via model-based
#' global optimization" \emph{In Proc. Int. Joint Conf. Neural Networks,
#' 1431-1438 }.\cr
#' 
#' Sill M., Hielscher T., Becker N. and Zucknick M. (2014), \emph{c060:
#' Extended Inference with Lasso and Elastic-Net Regularized Cox and
#' Generalized Linear Models, Journal of Statistical Software, Volume 62(5),
#' pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords models multivariate graphs iteration optimize
#' @examples
#' 
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
#' 
#' \donttest{
#' # y - binomial
#' y.classes<-ifelse(y>= median(y),1, 0)
#' set.seed(1234)
#' nfolds = 10
#' foldid <- balancedFolds(class.column.factor=y.classes, cross.outer=nfolds)
#' bounds <- t(data.frame(alpha=c(0, 1)))
#' colnames(bounds)<-c("lower","upper")
#'  
#' fit <- EPSGO(Q.func="tune.glmnet.interval", 
#'              bounds=bounds, 
#'              parms.coding="none", 
#'              seed = 1234, 
#'              show="none",
#'              fminlower = -100,
#'              x = x, y = y.classes, family = "binomial", 
#'              foldid = foldid,
#'              type.min = "lambda.1se",
#'              type.measure = "mse")
#' summary(fit)
#' }
#' 
#' \dontrun{
#' # y - multinomial: low - low 25%, middle - (25,75)-quantiles, high - larger 75%.
#' y.classes<-ifelse(y <= quantile(y,0.25),1, ifelse(y >= quantile(y,0.75),3, 2))
#' set.seed(1234)
#' nfolds = 10
#' foldid <- balancedFolds(class.column.factor=y.classes, cross.outer=nfolds)
#' bounds <- t(data.frame(alpha=c(0, 1)))
#' colnames(bounds)<-c("lower","upper")
#'  
#' fit <- EPSGO(Q.func="tune.glmnet.interval", 
#'              bounds=bounds, 
#'              parms.coding="none", 
#'              seed = 1234, 
#'              show="none",
#'              fminlower = -100,
#'              x = x, y = y.classes, family = "multinomial", 
#'              foldid = foldid,
#'              type.min = "lambda.1se",
#'              type.measure = "mse")
#' summary(fit)
#' }
#' 
#' \dontrun{
#' ##poisson
#' N=500; p=20
#' nzc=5
#' x=matrix(rnorm(N*p),N,p)
#' beta=rnorm(nzc)
#' f = x[,seq(nzc)]%*%beta
#' mu=exp(f)
#' y.classes=rpois(N,mu)
#' 
#' nfolds = 10
#' set.seed(1234)
#' foldid <- balancedFolds(class.column.factor=y.classes, cross.outer=nfolds)
#' 
#' fit <- EPSGO(Q.func="tune.glmnet.interval", 
#'              bounds=bounds, 
#'              parms.coding="none", 
#'              seed = 1234, 
#'              show="none",
#'              fminlower = -100,
#'              x = x, y = y.classes, family = "poisson", 
#'              foldid = foldid,
#'              type.min = "lambda.1se",
#'              type.measure = "mse")
#' summary(fit)
#' }
#' 
#' \dontrun{
#' #gaussian
#' set.seed(1234)
#' x=matrix(rnorm(100*1000,0,1),100,1000)
#' y <- x[1:100,1:1000]%*%c(rep(2,5),rep(-2,5),rep(.1,990))
#' 
#' foldid <- rep(1:10,each=10)
#' 
#' fit <- EPSGO(Q.func="tune.glmnet.interval", 
#'              bounds=bounds, 
#'              parms.coding="none", 
#'              seed = 1234, 
#'              show="none",
#'              fminlower = -100,
#'              x = x, y = y, family = "gaussian", 
#'              foldid = foldid,
#'              type.min = "lambda.1se",
#'              type.measure = "mse")
#' summary(fit)  
#' }
#' 
#' # y - cox in vignette
#' 
#' 
#' @export EPSGO
EPSGO<- function(
								# function to be minimized
								Q.func, 
								# bounds for parameters
								bounds, 
                # round.n -number of digits after comma
                round.n=5,
								# parms.coding="none", # or log2 
								parms.coding="none", # or log2 
								# min value for the function Q.func
								fminlower=0, 
								# do you want to find one min value and stop?
								flag.find.one.min =FALSE,
								# show plots ?   none, final iteration, all iterations 
								show=c("none", "final", "all"),
								# define the number of start points
								N= NULL, 
								#% maximum # of function evaluations 
								maxevals =   500,  
								# plot parameter
                 pdf.name=NULL, 
                 pdf.width=12, pdf.height=12,
                 my.mfrow=c(1,1),  
                # verbose?
                verbose=TRUE,
                seed=123, 
                ... ){
 

 
	# The EPSGO algorithm (from Holger's paper)
	#
	#EPSGO = Efficient Parameter Selection via Global Optimization

	#Input: function Q.func to measure gen. error
	#			  bounds: data frame with parameter bounds
	#								rownames: parameter names
	#								colnames: "lower" and "upper"
	#								Example: bounds=t(data.frame(lambda1=c(0, 10), lambda2=c(0,8)));colnames(bounds)<-c("lower", "upper")
	#								fminlower:  min Q value				
	#Output: Qmin, p*, number of visited points neval 
	
	## Scheme:
	#	1. number of tuning patameters
	#		 D = dim(P)
	#	2. create N = 10D sample points p1, ..., pN
	#		in [bounds] using Latin hypercube sampl.
	#	3. compute X= Q(p_i), i = 1, ...,N
	# 
	# 4. train Online GP
	#	5. number of visited points p_obs 
	#	 	 neval = N
	#	6.	REPEAT
	#					6.1  Find a new p, with max E[I(p)] ( the same as min -E[I(p)] )  
	#					?p = argmaxp E[I(p)] (computed by DIRECT)
	#					# Important! Direct.R calculate global min! --> change Problem function
	#					6.2 compute std. dev. and mean of E[I(p)]
	#					6.3 new p, new Q(p) 
	#						  Qnew = Q(?p)
	#					6.4 Add the new p ?
	#							if Qnew < Qmin
	#								Qmin = Qnew
	#								?p = ?p
	#							end
	#					6.5 update Online GP
	#							neval = neval + 1
	#			UNTIL convergence
	
	if (verbose){
	  print("parms.coding")
	  print(parms.coding)
	}
	
	
	###################################################################################################################
		## 1. ## number of tuning patameters D = dim(P)
	###################################################################################################################
	
	D<- nrow(bounds)
	
	###################################################################################################################
		## 2. ##  create N = 10D sample points p1, ..., pN
		#in [l, u] using Latin hypercube sampl.
	###################################################################################################################
	
	set.seed(seed)
	# in ego.m  X = start.points
	# N - number of start points in parameter space
	# wie in ego.m
	
	# define the number of start points 
	if (is.null(N)){
		ns<-c(21, 21, 33, 41, 51, 65, 65)  ##
		# N = 10D, D= number of parameters , but for high dim data with more than 6 dim, restrict the initial set of p to 65
		N<- ifelse ( D <= length(ns),ns[D], 65 )
		#N <- mult.factor * D
	}
	
	# start points X (= p_1,..., p_N)
	
	X<- lhs(N, bounds)
  
	# round.n -number of digits after comma
	X<-round(X,round.n)
	
	if (verbose) print(X)
	
	if(show !="none" & !is.null(pdf.name)) { 
		pdf(pdf.name, pdf.width, pdf.height)
		par(mfrow=my.mfrow)
	}

	if ( (show !="none") & (D<=2 )){
		# 1D plot
		if (D==1) {
			plot(X, xlab="Index", ylab=rownames(bounds)[1], col="orange", pch=20,
						 main=paste( "Latin hypercube sampling, n=",nrow(X)*D) )
			abline(h=seq(bounds[1,1],bounds[1,2],length=(nrow(X)+1)), lty=2, col=3)
		} else {
			# 2D plot
			plot(X, xlab=rownames(bounds)[1], ylab=rownames(bounds)[2], col="orange", pch=20, 
						 main=paste( "Latin hypercube sampling, n=",nrow(X)))
			abline(v=seq(bounds[1,1],bounds[1,2],length=(nrow(X)+1)), lty=2, col=3)
			abline(h=seq(bounds[2,1],bounds[2,2],length=(nrow(X)+1)), lty=2, col=3)
		}
	}


	###################################################################################################################
	## 3. ## compute Q(p_i), i = 1, ...,N
	###################################################################################################################
	
	# model.list<-apply(X, 1, eval(Q.func), x.svm=x.svm, y.svm=y.svm, maxIter=maxIter,parms.coding=parms.coding, inner.val.method=inner.val.method, cross.inner=cross.inner, seed=seed, verbose=verbose )
	model.list<-apply(X, 1, eval(Q.func), parms.coding=parms.coding,
									 maxevals=maxevals, seed=seed, verbose=verbose, ...)
	#debugging
  #save(model.list,file="model_list.RData") 
	
	
	# take the Q.values
	Q<- as.numeric(unlist( sapply(model.list, "[", "q.val")))
	Q.min <- min (Q, na.rm=T)
	min.p<- X[which.min(Q), , drop=FALSE]
	
	if (verbose) print(data.frame(X,Q))
	
	# delete the point(s) with no Q value(s)
	X <- X[!is.na(Q), ,drop=FALSE]
	Q <-Q[!is.na(Q)]
	
	# add start.q.values to the plot
	if (show !="none" & D ==1)  text( c(1: nrow(X)),X, labels=round(Q,6), pos=1, cex=0.5 )
	if (show !="none" & D ==2)  text( X, labels=signif(Q,6), pos=1, cex=0.5 )
	
	#		# 4. train Online GP
	
	# train Gaussian Process 
	# Input:  collection of random variables G(x)(here: x= points in tuning parameter space = start.points X ) 
	#
	###################################################################################################################
	## 4. ## train Gaussian Process
	###################################################################################################################
	
		gp.seed.new<- seed 
			
		# if we have tried 5 times and are still not able to fit --> break
			# the reason is in having a new point in Ytrain very close to one of the old ones.
			# --> matrix is singular! 
			
			if (exists("fit.gp")) rm(fit.gp)
			flag.fit.gp<- FALSE
      tmp.i<-1
      
     	while( ! flag.fit.gp   ) {
				
				if (tmp.i >5) {
					print(print( "At least one of the intial points is very close to the other points. It is not possible to fit the model via Gaussian Process. "))
					break  
				}
				try(fit.gp<-mlegp(X, Q,constantMean=0, seed=gp.seed.new,  verbose= 0  ))
				
				flag.fit.gp<- FALSE
      	# if fit.gp exists AND is not null
      	if (exists("fit.gp"))
      		if (!is.null(fit.gp))
      				flag.fit.gp<- TRUE
      
				# if fails to fit change seed
				if(!flag.fit.gp  ) {
					 print("fails to fit gp (fit.gp), change seed !")
					 gp.seed.new<- round(runif(1,min=1, max=10^3))
					 tmp.i<-tmp.i + 1
				}
			} # end of while 
			
			gp.seed<- gp.seed.new
		
	# define the exprected improvement function  --> function ExpImprovement.R
	# E(I(p)) = (Q_min - ^mu(p)) * PHI([Q_min - ^mu(p)] / ^sig(p) ) + 
	# 										^sig(p)* phi([Q_min - ^sig(p)] / ^sig(p) )^
	
	###################################################################################################################
	#	5. number of visited points p_obs 
	###################################################################################################################
	neval <- length(Q)
	
	###################################################################################################################
	##.6.	REPEAT - main block
	###################################################################################################################
	
	# Initialization for observed points in parameter space Xtrain in R^D and 
	# 							 for quality value Ytrain = Q		  
	Xtrain<- X
	Ytrain <- Q
	finished <- FALSE
	EImax <- Inf
	
	loop<-1
	fcalls<-length(Q)
	# Point with max E[I(p)]
	xmax=X[1, ,drop=FALSE]
	#fmin <- Inf
	# change to already calculated!
		
	# info for point with min Q.func
	fminold<-  Inf 
	fmin <- Q.min
	xmin<- min.p
	not_changed = 0
	
	set.seed(seed)
	
	while (!finished){ 
		if (verbose ) print(paste("loop", loop))
		# calculate Q for the new points 		
		if (loop >1) {
      EIold <- EImax	
		  fminold <- fmin
		
			
			#model.list.new<-apply(X, 1, eval(Q.func), x.svm=x.svm, y.svm=y.svm,
			#		 maxIter=maxIter,parms.coding=parms.coding, inner.val.method=inner.val.method,  cross.inner=cross.inner ,verbose=verbose )
      #model.list.new<-apply(X, 1, eval(Q.func), x=x, y=y, family=family, nfolds=nfolds,  maxit=maxit, seed=seed, verbose=verbose ,type.min=type.min)
			model.list.new<-apply(X, 1,  eval(Q.func), parms.coding=parms.coding,
									 								maxevals=maxevals, seed=seed, verbose=verbose, ...)
	   	# take the Q.values
      model.list<-c(model.list, model.list.new )
      Q<- as.numeric(unlist( sapply(model.list.new, "[", "q.val")))
			
			fcalls = fcalls + length(Q);
			Xtrain = rbind(Xtrain, X)
			Ytrain = c(Ytrain, Q )	
			
      if (verbose) print(data.frame(Xtrain,Ytrain))
			
			# fmin =  current min of Q.func at point xmin. 
			fmin<- min(Ytrain)
			xmin = Xtrain[which.min(Ytrain),]	
		}
		
		if (verbose) print(paste("fmin=",fmin))
		if (fmin < fminlower){ 
				break
			} 
		
		if (flag.find.one.min){
			### skip it (in case of having more than 1 points with global min )
			# break if reach the min value
			if (fmin <= fminlower){ 
				break
			} 
		}
		
		
		# break if no changes in the last 10 iterations
		if (fmin == fminold){
			not_changed = not_changed + 1
			if (not_changed >= 10){
				print("No changes in the last 10 iterations, break iterations")
				break
			}
		}else {
			not_changed = 0
		} 
		
		# train GP
		if (loop >1) {
			# sometimes error:  Error in solve.default(gp$invVarMatrix) :
	  	#					system is computationally singular: reciprocal condition number = 7.502e-17
			# Solution: change seed
						
			gp.seed.new<- c(seed + loop-1) 
			
			
			if (exists("fit.gp")) rm(fit.gp)
			
			# if we have tried 5 times and are still not able to fit --> break
			# the reason is in having a new point in Ytrain very close to one of the old ones.
			# --> matrix is singular! 
			
		  tmp.i<-1
		  flag.fit.gp<- FALSE
      
     	while(!flag.fit.gp   ) {
				
				if (tmp.i >5) {
					print(print( "The new point X is very close to the one of visited points."))
					finished <- TRUE
					break  
				}
				try(fit.gp<-mlegp(Xtrain, Ytrain,constantMean=0, seed=gp.seed.new,  verbose= 0  ))
				
				flag.fit.gp<- FALSE
      	# if fit.gp exists AND is not null
      	if (exists("fit.gp"))
      		if (!is.null(fit.gp))
      				flag.fit.gp<- TRUE
      
				# if fails to fit change seed
				if(!flag.fit.gp  ) {
					 print("fails to fit gp (fit.gp), change seed !")
					 gp.seed.new<- round(runif(1,min=1, max=10^3))
					 tmp.i<-tmp.i + 1
				}
			} # end of while  (!flag.fit.gp   ) 
			
			if (finished) break
			
			gp.seed<-c(gp.seed, gp.seed.new)
			#str(fit.gp)
		}
		
		if (show !="none" ) plot(fit.gp, main=paste("Gaussian Process", "iter ", loop))
			
		Problem<- list(f = "ExpImprovement" )
					
		#6.1  Find a new p, with max E[I(p)] ( the same as min -E[I(p)] ) 
		#[EImax xmax history] = Direct(Problem,bounds,options) 
		Dir.list<- Direct(Problem=Problem, 
								bounds=bounds, 
								# options 
								#% maximum of iterations 
								maxits =     50,   
								#% maximum # of function evaluations 
								maxevals = maxevals,       
								# the optimal value is unknown
								testflag= 0,
								# minimum value of function , min (ExpImprovement) = 0
								globalmin =  0, 
								#% print and plot iteration stats
								showits= show, 
								# verbose? 
								verbose=FALSE,
								tol= 0.01, 
								pdf.name=NULL, #"test.pdf" # 
								pdf.width=12, pdf.height=12,
								my.mfrow=c(1,1), #c(3,3),
								#
								# additional args for the problem function 
								fmin=fmin,
								fit.gp=fit.gp
								 )
		
		# E[I(p)], xmax - p with max EI
		EImax <- (-1 ) * Dir.list$minval
		xmax<-  Dir.list$final_point.xatmin
		History <-  Dir.list $ History
		
		if (verbose ){
			print(paste("EImax=", EImax ) )
			print( "xmax:")
			print( xmax)
			print( "history:")
			print( History)
		}
		
		# 6.2 compute std. dev. and mean of E[I(p)]
		### other random sample in the same size; Latin hypercube sampling over the 
		## whole parameter space. The idea is to make sure, that the expected improvement over the whole space is almost equally small. 
		Xsam <- lhs(N, bounds)
		EIall=rep(0,N)	
		for (i in 1:N ){		
			EIall[i] = - ExpImprovement(Xsam[i, , drop=FALSE], fmin, fit.gp, muX=NULL, muY=NULL)
		}
		
		EIstd =  sd( c(EIall, EImax) )
		EImean = mean(c(EIall, EImax))
		
		# if E[I(p)] is max by a random sample --> take it!
		if ( max(EIall) > EImax ){
			xmax = t(Xsam[which.max(EIall),, drop=FALSE])
			rownames(xmax)<- rownames(bounds)
			EImax <- max(EIall)
		}
		
		# 6.4 Add the new p ?
		# if the new p is one from the observed ones --> break
		# OR if the differences in functions is small,
		#  stop iterations
		
		# Is xmin already in the set of visited points?
		visited <-  Xtrain - matrix(xmax, nrow=nrow(Xtrain), ncol= ncol(Xtrain), byrow=TRUE) 
		visited<- any (	apply (visited == 0 , 1, all )  )
		
		if ( (visited )   | (abs(fmin - fminold) < 0.01) & ((EImax - EImean)^2 <= 0.1*EIstd)  ){
			if (visited) print( "The new point with min E[I(p)] is already in the set of visited points.")
			if ((abs(fmin - fminold) < 0.01) & ((EImax - EImean)^2 <= 0.1*EIstd) ) {
				print( "the differences in functions between 2 last iterations is small, stop iterations") 
			}
			finished = TRUE
			#		} else {
			#		# print("Take the second best point (candidate) ")
			#	
			#		# filter out E[I(p)]=fc==0, 
			#		cand<-Dir.list$fc[Dir.list$fc != 0 ]
			#		# mult by (-1) because of changing from max to min problem in Direct
			#		EImax2<-sort((-1) * cand, decr=T)[2]
			#	
			#		xmax <-  Dir.list$c[ , which((-1)*Dir.list$fc == EImax2)  ]
			#		X = t(xmax)
			#		EImax<- EImax2
			#		neval<- neval + 1
				
				
		}else{
			X = t(xmax)
      X<-round(X,round.n)
			rownames(X)<- NULL
			neval<- neval + 1
		}
		if (verbose ) print(paste("iteration :", loop, "  fmin = ", fmin ))
		print(paste("finished?", finished))
		print("X")
		print(X)
		loop = loop + 1
		
		print(data.frame(Xtrain,Ytrain))
	
	} # end of while (!finished) 
		
	if(show !="none" & !is.null(pdf.name)) dev.off()	
		
	# define the set of points with the same fmin
		tmp.set<-data.frame(Xtrain, f=Ytrain)
		points.fmin<- tmp.set[tmp.set$f == fmin, ,drop=FALSE  ]
	
  out <- list(fmin =fmin,
	     xmin = xmin, 
	     iter = loop,
	     neval =neval,
	     maxevals= maxevals,
	     seed =seed,
	     bounds = bounds,
	     Q.func =Q.func,
	     points.fmin =points.fmin,
	     Xtrain = Xtrain, 
	     Ytrain= Ytrain, 
	     gp.seed=gp.seed,
	     model.list = model.list 
	) 
  
  class(out) <- "intsearch"
	return(out)
	
}
