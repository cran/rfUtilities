#' @title Random Forest Model Selection
#' @description Implements Murphy et al., (2010) Random Forests model selection approach. 
#' 
#' @param ydata Y Data for model
#' @param xdata X Data for model
#' @param imp.scale  Type of scaling for importance values (mir or se), default is mir
#' @param r Vector of importance percentiles to test i.e., c(0.1, 0.2, 0.5, 0.7, 0.9)
#' @param final Run final model with selected variables (TRUE/FALSE)
#' @param plot.imp Plot variable importance (TRUE/FALSE)
#' @param seed Sets random seed in the R global environment. This is hightly suggested.
#' @param parsimony  Threshold for competing model (0-1)
#' @param ... Arguments to pass to randomForest (e.g., ntree=1000, replace=TRUE, proximity=TRUE)
#'
#' @return A list class object with the following components:
#'  @return   rf.final - Final selected model, if final=TRUE(randomForest model object)
#'  @return   selvars - final selected variables (vector)
#'  @return   test - Validation parameters used on model selection (data.frame)
#'  @return   importance - Importance values for selected model (data.frame)
#'  @return   parameters - Variables used in each tested model (list)
#'
#' @note If you want to run classification, make sure that y is a factor, otherwise runs in regression mode
#' @note The mir scale option perfroms a row standardization and the se option performs normalization using The "standard errors" 
#' @note of the permutation-based importance measure. Both options result in a 0-1 range but "se" summs to 1.
#' @note The selection criteria are calculated as: mir = i/max(i) and se = (i / se) / ( sum(i) / se).
#' @note For regression the model selection criteria is; largest %variation explained, smallest MSE, and fewest parameters.
#' @note For classification; Smallest OOB error, smallest maximum within class error, and fewest parameters.
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @references
#' Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species Using Random Forest. Landscape Ecology 5:673-683.
#' @references
#' Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas connectivity in Yellowstone National Park with landscape genetics. Ecology 91:252-261
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer
#' @examples
#' # Classification on iris data
#' require(randomForest)
#' data(iris)
#'   iris$Species <- as.factor(iris$Species)
#' ( rf.class <- rf.modelSel(iris[,1:4], iris[,"Species"], seed=1234, imp.scale="mir") )
#' ( rf.class <- rf.modelSel(iris[,1:4], iris[,"Species"], seed=1234, imp.scale="mir", 
#'                           parsimony=0.03) )
#'     vars <- rf.class$PARAMETERS[[3]]
#'       ( rf.fit <- randomForest(x=iris[,vars], y=iris[,"Species"]) )
#'
#' # Regression on airquality data
#' data(airquality)
#'   airquality <- na.omit(airquality)
#' ( rf.regress <- rf.modelSel(airquality[,2:6], airquality[,1], imp.scale="se") )
#' ( rf.regress <- rf.modelSel(airquality[,2:6], airquality[,1], imp.scale="se", parsimony=0.03) )
#'     vars <- rf.regress$PARAMETERS[[3]]
#'       ( rf.fit <- randomForest(x=airquality[,vars], y=airquality[,1]) )
#'
rf.modelSel <- function(xdata, ydata, imp.scale="mir", r=c(0.25, 0.50, 0.75),  
                        final=FALSE, plot.imp=TRUE, seed=NULL, parsimony=NULL, ...) 
  {
 rf.ImpScale <- function (x, scale="mir") { 
  if (!inherits(x, "randomForest")) 
       stop(deparse(substitute(x)), " Must be a randomForest object")
	if(!is.null(seed)) { set.seed(seed) }
  if (x$type == "regression") {
   if (is.null(x$importanceSD) == TRUE | "%IncMSE" %in% 
       names(as.data.frame(x$importance)) == FALSE)
        stop("OBJECT DOES NOT CONTAIN PERMUTATED IMPORTANCE, PLEASE RUN 
              randomForest WITH importance=TRUE")   
	rf.imp <- x$importance[,"%IncMSE"]
    rf.impSD <- x$importanceSD
       rf.impSD[rf.impSD == 0] <- 0.000000001	
    if (scale == "mir") {
      i <- rf.imp / max(rf.imp) 
		}	  
    if (scale == "se") {
	  i <- ( rf.imp / rf.impSD ) / sum(rf.imp / rf.impSD, na.rm=TRUE)			 
        }
	 }
  if (x$type == "classification" | x$type == "unsupervised") {
   if (is.null(x$importanceSD) == TRUE | "MeanDecreaseAccuracy" %in% 
       names(as.data.frame(x$importance)) == FALSE)
        stop("OBJECT DOES NOT CONTAIN PERMUTATED IMPORTANCE, PLEASE RUN 
       randomForest WITH importance=TRUE") 
	rf.imp <- x$importance[,"MeanDecreaseAccuracy"]
    rf.impSD <- x$importanceSD[,"MeanDecreaseAccuracy"]
       rf.impSD[rf.impSD == 0] <- 0.000000001	
    if (scale == "mir") {
      i <- rf.imp / max(rf.imp) 
		}	  
    if (scale == "se") {
	  i <- ( rf.imp / rf.impSD ) / sum(rf.imp / rf.impSD, na.rm=TRUE)			 
        }
	 }
 	i <- as.data.frame(i)
	  names(i) <- c("importance") 
        row.names(i) <- names(rf.imp)	
   return( i )            
 }
 
RFtype <- is.factor(ydata) #TEST FOR FACTOR IN Y 
##CLASSIFICATION##
if (RFtype == "TRUE") {
    model.vars <- list()
    ln <- 0
    rf.all <- randomForest::randomForest(x=xdata, y=ydata, importance=TRUE, ...) 
      model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)  
       class.errors <- as.data.frame(rf.all$err.rate)
        class.errors <- na.omit(class.errors)  
         class.errors[class.errors == NaN] <- 0
          class.errors[class.errors == Inf] <- 1         
        i <- vector()
	      for ( l in 2:nlevels(as.factor(names(class.errors))) ) {              
            i <- append(i, median(class.errors[,l]))
          }        
        max.error = max(i) 
	imp <- rf.ImpScale(rf.all, scale=imp.scale) 
    results <- as.data.frame(array(0, dim=c( 0, 4 )))
      x <- c(0, (median(rf.all$err.rate[,"OOB"]) * 100), max.error * 100, dim(xdata)[2] )
    results <- rbind(results, x) 	 	 
     for (p in 1:length(r) ) {
		 t = quantile(imp[,1], probs=r[p])
         sel.imp <- subset(imp, importance >= t)
           sel.vars <- rownames(sel.imp)
     if (length( sel.vars ) > 1) {                             
         xdata.sub <- xdata[,sel.vars]       
      rf.model <- randomForest::randomForest(x=xdata.sub, y=ydata, importance=TRUE)          
           class.errors <- as.data.frame(rf.model$err.rate)
            class.errors <- na.omit(class.errors)  
             class.errors[class.errors == NaN] <- 0
              class.errors[class.errors == Inf] <- 1      
        i <- as.vector(array(0, dim=c((0),(1))))
       for ( l in 2:nlevels(as.factor(names(class.errors))) )
          {
          x.bar <- mean(class.errors[,l])              
            i <- as.vector(append(i, x.bar, after=length(i) ))
            }        
         max.error = max(i[2:length(i)] )     
         x <- c(t, median(rf.model$err.rate[,1]) * 100, max.error * 100, length(sel.vars) )
         results <- rbind(results, x)
		  model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)
         }
        }
  names(results) <- c("THRESHOLD", "OOBERROR", "CLASS.ERROR", "NPARAMETERS")
  results <- results[order(results$CLASS.ERROR, results$OOBERROR, results$NPARAMETERS),]
    if (is.null(parsimony) == FALSE) { 
	  if(parsimony < 0.00000001 | parsimony > 0.9) stop( "parsomony MUST RANGE 0-1")
        oob <- "TRUE"
        for(i in 2:nrow(results)) {
          if( abs((results[i,][2] - results[1,][2] ) / results[1,][2]) <= parsimony  &
              abs( (results[i,][3] - results[1,][3] ) / results[1,][3] ) <= parsimony ) {
            oob <- append(oob, "TRUE")
        	  } else {
        	oob <- append(oob, "FALSE")
            }
            final <- results[which( oob == "TRUE" ),]
        	  final <- final[final$NPARAMETERS == min(final$NPARAMETERS) ,]$THRESHOLD
        }
          } else {		
            final <- as.vector(results[,"THRESHOLD"])[1]
        }	
      sel.imp <- subset(imp, importance >= final)    
        sel.vars <- rownames(sel.imp)
          sel.post=which( results$NPARAMETERS == length(sel.vars) ) 
      results <- rbind(results[sel.post,],results[-sel.post,]) 	
  } # END OF CLASSIFICATION
  
##REGRESSION## 
if (RFtype == "FALSE") {
    model.vars <- list()
      ln <- 0      
    rf.all <- randomForest::randomForest(x=xdata, y=ydata, importance=TRUE, ...) 
	  model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)  
	   imp <- rf.ImpScale(rf.all, scale=imp.scale) 
      results <- as.data.frame(array(0, dim=c( 0, 4 )))
     x <- c(0, (median(rf.all$rsq)), mean(rf.all$mse), dim(xdata)[2] )
     results <- rbind(results, x)     
   for (p in 1:length(r) ) {
        t = quantile(imp[,1], probs=r[p])		 
        sel.vars <- rownames(subset(imp, importance >= t))  
     if (length( sel.vars ) > 1) {                             
      xdata.sub <- as.data.frame(xdata[,sel.vars]) 
      rf.model <- randomForest::randomForest(x=xdata.sub, y=ydata, importance=TRUE, ...)          
      x <- c(t, (median(rf.model$rsq)), mean(rf.model$mse), length(sel.vars) )
      results <- rbind(results, x)
	     model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)
      }
    }
   names(results) <- c("THRESHOLD", "VAREXP", "MSE", "NPARAMETERS")
     results <- results[order(-results$VAREXP, results$MSE, results$NPARAMETERS),]  
    if (!is.null(parsimony)) {
      if(parsimony < 0.00000001 | parsimony > 0.9) stop( "parsomony MUST RANGE 0-1")	
        oob <- "TRUE"
        for(i in 2:nrow(results)) {
          if( abs((results[i,][2] - results[1,][2] ) / results[1,][2]) <= parsimony  &
              abs( (results[i,][3] - results[1,][3] ) / results[1,][3] ) <= parsimony ) {
            oob <- append(oob, "TRUE")
        	  } else {
        	oob <- append(oob, "FALSE")
            }
            final <- results[which( oob == "TRUE" ),]
        	  final <- final[final$NPARAMETERS == min(final$NPARAMETERS) ,]$THRESHOLD
        }
          } else {		
            final <- as.vector(results[,"THRESHOLD"])[1]
        }	
      sel.imp <- subset(imp, importance >= final)    
        sel.vars <- rownames(sel.imp)
        sel.post=which( results$NPARAMETERS == length(sel.vars) ) 
      results <- rbind(results[sel.post,],results[-sel.post,]) 			
   } # END OF REGRESSION 
   
    if (plot.imp == TRUE) {
        if (imp.scale=="mir") {lable="Row Standardization Variable Importance"} 	
		  if (imp.scale=="se") {lable="Standardized Error Variable Importance"}
	  p <- as.matrix(subset(imp, importance >= final))    
       ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]])  
       dotchart(p[ord,1], main=lable, pch=19)
    }
	
    if (final == TRUE) {
       sub.xdata <- xdata[,sel.vars]  #SUBSET VARIABLES
        rf.final <- randomForest::randomForest(x=sub.xdata, y=ydata, importance=TRUE, ...)           
      ( list(MODEL=rf.final, SELVARS=sel.vars, TEST=results, IMPORTANCE=sel.imp, PARAMETERS=model.vars) )      
         } else {
      ( list(SELVARS=sel.vars, TEST=results, IMPORTANCE=sel.imp, PARAMETERS=model.vars) ) 
    }     
 }
