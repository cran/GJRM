\name{hazsurv.plot}
\alias{hazsurv.plot}

\title{Hazard and survival plots}

\usage{

hazsurv.plot(x, eq, newdata, type = "surv", t.range = NULL, intervals = TRUE, 
             n.sim = 100, prob.lev = 0.05, 
             shade = FALSE, ylim, ylab, xlab, ls = 100, baseline = FALSE, 
             pop.name = NULL, pop.min = NULL,  pop.max = NULL, pop.bin = NULL, 
             pop.build = FALSE, pop.grid = 200, 
             min.dn = 1e-200, min.pr = 1e-200, max.pr = 1, ...)



}

\arguments{

  \item{x}{A fitted \code{gamlss}/\code{gjrm} object.}

  \item{eq}{Equation number. This can be ignored for univariate models.}

  \item{newdata}{A data frame or list containing the values of the model covariates at which predictions are required. 
                 This must have one row and must be provided.}

  \item{type}{The plot to produce, either \code{"surv"}, \code{"hazard"} or \code{"cumhaz"}.}
  
  \item{t.range}{Time variable range to be considered for plotting. If \code{NULL} then it is determined automatically based on the observed data.}
  
  \item{intervals}{If \code{TRUE} then intervals are also produced.}

  \item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used
             for interval calculations.}
  
  \item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.}  
   
  \item{shade}{If \code{TRUE} then it produces shaded regions as confidence bands.}
   
  \item{ylim, ylab, xlab}{Usual plot arguments.}
  
  \item{ls}{Length of sequence to use for time variable.}
  
  \item{baseline}{If baseline is desired; this will set all covariate/smooth effects to zero.}
  
  \item{pop.name}{Variable name with respect to which you want to stratify.}
  
  \item{pop.min}{Min. value of interval for the stratification variable. If you would like to use an open interval then do not specify a 
                 value for \code{pop.max}; in this case you would have \code{pop.name} > \code{pop.min}.}
  
  \item{pop.max}{Max. value of interval for the stratification variable. For an open interval do not specify \code{pop.min} (\code{pop.name} < \code{pop.max}).}
  
  \item{pop.bin}{Either 0 or 1. This means \code{pop.name} is a binary variable and you want to stratify wrt either 0 or 1. Note that in 
                 this case \code{pop.min} and \code{pop.max} need not be used as the choice of 0 or 1 takes the place of the range of values in a binary setting.}
  
  \item{pop.build}{The default is \code{FALSE}, which is the setting in which the user provides the data. When \code{TRUE} we are 
                   assuming a setting in which the user provides a fixed set of features and would like the average net survival for
                   an individual whith pop.name in a certain range. Note that \code{pop.bin} can not be used in this setting as it would not 
                   make sense to do so (we would just have an individual net survival since it does not make sense to build a range
                   of values for a binary variable).}
  
  \item{pop.grid}{Default values is 200 and it is only used if \code{pop.build == TRUE}. It represents the number
                  of rows of the (internally built) dataset. In other terms it is the number of elements forming the sequence of
                  extremities \code{pop.min} and \code{pop.max}.}
  
  \item{min.dn, min.pr, max.pr}{Allowed minimum and maximum for estimated probabities and densities for survival, hazard and cumulative hazard calculations.}

  
  \item{...}{Other arguments to pass to plot.}

}



\value{It produces a plot or set of plots.}



\description{ 

This function produces either a survival or hazard plot. 

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\keyword{hplot} 
\keyword{smooth} 
\keyword{regression}
