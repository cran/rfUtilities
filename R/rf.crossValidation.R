#' @title Random Forest Classification Model Cross-validation 
#' @description Implements an n-Fold Cross-validation for Random Forests classification models
#'    
#' @param x random forest object
#' @param xdata x data used in model
#' @param p Percent data withhold
#' @param n Number of cross validations
#' @param seed Sets random seed in R global environment
#' @param plot plot cross-validation error statistic (TRUE/FALSE)
#' @param ... Additional arguments passed to Random Forests 
#'
#' @return A list class object with the following components:
#' @return   cv.Summary - Dataframe with summary statistics for error, pcc and ob
#' @return   Error.distribution - Vector of error values for each cv
#' @return   PCC.distribution - Vector of pcc values for each cv
#' @return   OOB.distribution - Vector of oob values for each cv
#'
#' @note
#' The crossvalidation statistics are based on the prediction error on the witheld data: 
#' percent correctly classified (PCC) and out-of-bag model error.  
#' The plot y axis represents the probability density function (see ?density)
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references
#' Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species Using Random Forest. Landscape Ecology 5:673-683.
#' @references
#' Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas connectivity in Yellowstone National Park with landscape genetics. Ecology 91:252-261
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#' 
#' @examples 
#' require(randomForest)
#'   data(iris)
#'     iris$Species <- as.factor(iris$Species)    	
#'       set.seed(1234)	
#' ( rf.mdl <- randomForest(iris[,1:4], iris[,"Species"], ntree=501) )
#'   ( rf.cv <- rf.crossValidation(rf.mdl, iris[,1:4], p=0.10, n=99, ntree=501) ) 	
rf.crossValidation <- function(x, xdata, p=0.10, n=99, seed=NULL, plot=TRUE, ...) {
    if (!inherits(x, "randomForest")) stop("x is not randomForest class object")
      if (!x$type == "classification") stop( "Random Forests Model is not classification") 
        if(!is.null(seed)) { set.seed(seed) }	
    sample.size = round( (length(x$y) * p) / length(x$classes), digits=0) 
      error <- vector()
        pcc <- vector()
		  oob <- vector()	
    for(i in 1:n) {	
	  samp.index <- vector()
        for(r in unique(x$y)) {			  
          cvalue <- which(x$y == r) 	  
            samp.index <- append(samp.index, sample(cvalue, sample.size)) 
        }
	    tx <- xdata[samp.index,]
	      ty <- x$y[samp.index]
	        mx <- xdata[-samp.index,]
	          my <- x$y[-samp.index] 
       rf.fit <- randomForest::randomForest(y=as.factor(my), x=mx, ytest=as.factor(ty), xtest=tx, ... )    
         error <- append(error, (median(rf.fit$test$err.rate[,1])) ) 
           pcc <- append(pcc, (sum(diag(rf.fit$test$confusion))/sum(rf.fit$test$confusion)) )
             oob <- append(oob, median(rf.fit$err.rate[,"OOB"]))	    	   
      }	 
    cv <- data.frame( c( summary(error), var(error) ), c( summary(pcc), var(pcc) ),
                      c( summary(oob), var(oob) ) )
      names(cv) <- c("ERROR","PCC", "OOB")
        rownames(cv)[7] <- "Variance" 	  	  
	if(plot == TRUE) {
	  par(mfrow=c(2,1))
        eden <- density(error)
         plot(eden, type="n", xlab="ERROR", ylab="", main="CROSS-VALIDATION ERROR DISTRIBUTION")
            polygon(eden, col = "blue2")
            lines(eden, lty=1)
              abline(v=min(error), lty=1)
               abline(v=max(error), lty=2)
              abline(v=mean(error), lty=3, col="red")
	    legend("topright", legend=c("min", "max", "mean"), 
	           lty=c(1,2,3), col=c("black", "black", "red"),
			   bg="white", bty="o", box.col="black", cex=0.75)			   
        pccden <- density(pcc)
        plot(pccden, type="n", xlab="PCC", ylab="", main="CROSS-VALIDATION PCC DISTRIBUTION")
            polygon(pccden, col = "blue2")
            lines(pccden, lty=1)
              abline(v=min(pcc), lty=1)
               abline(v=max(pcc), lty=2)
              abline(v=mean(pcc), lty=3, col="red")
	    legend("topleft", legend=c("min", "max", "mean"), 
	           lty=c(1,2,3), col=c("black", "black", "red"), 
			     bg="white", bty="o", box.col="black", cex=0.75)	
      }	  
    cat("\n MIN ERROR = ", min(error) )
      cat("\n MAX ERROR = ", max(error) )
        cat("\n MEAN ERROR = ", mean(error) )
          cat("\n MEDIAN ERROR = ", median(error) )
	        cat("\n ERROR VARIANCE = ", var(error) )	
    return( list( cv.Summary=cv, Error.distribution=error, PCC.distribution=pcc, 
	              OOB.distribution=oob) )		
}	
