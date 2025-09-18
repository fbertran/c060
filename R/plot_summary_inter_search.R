#' Plot Summary object for interval search models
#' 
#' Produces a plot for summary object of a fitted interval search model. Plot
#' 'visited' points against iteration steps. start.N points are initial points
#' selected before interval search starts.
#' 
#' @name plot.sum.intsearc
#' @param x an object of class \code{sum.intsearch} as returned by the function
#' \code{\link[c060]{summary.intsearch}}.
#' @param type type of plot to be drawn, \code{type="summary"} will plot the
#' partial log likelihood deviance as a function of both tuning parameters
#' \eqn{\alpha} and log \eqn{\lambda}. The final solution will
#' be highlighted by solid red line. Alternativly, \code{type="points"} will
#' draw the distribution of initial and visited points of the interval search
#' plotted in chronological order.
#' @param startN number of initial points. Needed if \code{type="points"}
#' @param ... additional argument(s)
#' @author Natalia Becker \ \email{natalia.becker@@dkfz.de}
#' @export
#' @seealso \code{\link[c060]{EPSGO}}, \code{\link[c060]{summary.intsearch}}
#' @references Sill M., Hielscher T., Becker N. and Zucknick M. (2014),
#' \emph{c060: Extended Inference with Lasso and Elastic-Net Regularized Cox
#' and Generalized Linear Models, Journal of Statistical Software, Volume
#' 62(5), pages 1--22.} https://doi.org/10.18637/jss.v062.i05.
#' @keywords plot
plot.sum.intsearch <-function(x,type="summary",startN=21,...){
  if(type=="summary"){
  summary.int <- x
  summary.int$info$log_lambda <- log(summary.int$info$lambda)
  breaks                        <- do.breaks(range(summary.int$info$deviance), 20)
  n_init                        <- 21 # number of initial alpha values at iteration zero
  summary.int$info$cols       <- level.colors(summary.int$info$deviance,at=breaks, col.regions = gray.colors)
  n.features                  <- summary.int$info$n.features

  print(my.plot <- xyplot(log_lambda ~ alpha,
                          data = summary.int$info,
                          groups = summary.int$info$cols,
                          cex = 1, cex.axis=1.5,
                          col = "black",
                          jitter.y=T, amount=0.01,
                          ylab=list(expression(paste("log ",lambda)),cex=1.5),
                          xlab=list(expression(alpha),cex=1.5),
               #          scales=list(x=list(log=T, equispaced.log = FALSE)), # x axis on log-scale
                          panel = function(x, y, groups, ..., subscripts) {
                              fill <- groups[subscripts]
                              panel.grid(h = -1, v = -1)
                              panel.abline(h = log(summary.int$opt.lambda),
                                           v = summary.int$opt.alpha,
                                           col="red", lty = 1, lwd=2 )
                              panel.xyplot(x, y, pch = rep(c(22,21),c(n_init,nrow(summary.int$info)-n_init)),
                                       fill = fill, ...) ;
                              ltext(x=x, y=y, labels=n.features, pos=ifelse(y<0.1,3,4), offset=1.5, cex=1,col=1)
                          },
                          legend = list(top = list(fun = draw.colorkey,
                                                   args = list(key = list(space = "top",
                                                                          col = gray.colors,
                                                                          at = breaks),
                                                               draw = FALSE))),
                          main="Cross-validated partial log likelihood deviance",
                          scales=list(cex=1)
                          #sub="number of selected features are printed next to symbol \n rectangles show initial     alpha values"
  ))
  }
  if(type=="points"){
    # plot visited points vs. iteration steps
    summary.int <- x
    niter<- nrow(summary.int$info) - startN
    iter<-c(rep(0,startN), c(1:niter))
    plot(summary.int$info$alpha, iter, xlab=expression(alpha), ylab="Iteration", pch=20,cex=1.5,cex.axis=1.5)
    grid(NA, niter+1, lwd=2)
    abline(v=summary.int$opt.alpha, col="red")
  }
}

