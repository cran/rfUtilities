% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf.significance.R
\name{rf.significance}
\alias{rf.significance}
\title{Random Forest model significance test}
\usage{
rf.significance(x, xdata, q = 0.99, p = 0.05, nperm = 999, ...)
}
\arguments{
\item{x}{randomForest class object}

\item{xdata}{Independent variables (x) used in model}

\item{q}{Quantile threshold to test classification models}

\item{p}{p-value to test for significance in regression models}

\item{nperm}{Number of permutations}

\item{...}{Additional Random Forests arguments}
}
\value{
A list class object with the following components:

For Regression problems:

RandRsquare Vector of random R-square values

Rsquare The R-square of the "true" model

Accept Is the model significant at specified p-value (TRUE/FALSE)

TestQuantile Quantile threshold used in significance plot

pValueThreshold Specified p-value

pValue p-values of randomizations

nPerm Number of permutations

For Classification problems:

RandOOB Vector of random out-of-bag (OOB) values

RandMaxError Maximum error of randomizations

test.OOB Error if the "true" model

Accept Is the model significant at specified p-value (TRUE/FALSE)

TestQuantile Quantile threshold used in significance plot

pValueThreshold Specified p-value

pValue p-values of randomizations

nPerm Number of permutations
}
\description{
Performs significance test for classification and regression Random Forests models.
}
\examples{
\dontrun{
# Regression
require(randomForest)
  set.seed(1234)	
    data(airquality)
    airquality <- na.omit(airquality)
 ( rf.mdl <- randomForest(x=airquality[,2:6], y=airquality[,1]) )
   ( rf.perm <- rf.significance(rf.mdl, airquality[,2:6], nperm=99, ntree=501) )
 
# Classification
require(randomForest)
  set.seed(1234)	
    data(iris)
      iris$Species <- as.factor(iris$Species) 
 ( rf.mdl <- randomForest(iris[,1:4], iris[,"Species"], ntree=501) )
   ( rf.perm <- rf.significance(rf.mdl, iris[,1:4], nperm=99, ntree=501) ) 
}

}
\references{
Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas connectivity in Yellowstone National Park with landscape genetics. Ecology 91:252-261

Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer
}
\author{
Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
}
