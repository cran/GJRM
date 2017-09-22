\name{hazsurv.plot}
\alias{hazsurv.plot}

\title{Hazard and survival plots}

\usage{

hazsurv.plot(x, eq, newdata, type = "surv", intervals = TRUE, 
             n.sim = 100, prob.lev = 0.05, 
             shade = FALSE, ylim, ylab, xlab, ls = 100, ...)



}

\arguments{

  \item{x}{A fitted \code{gamlss}/\code{gjrm} object.}

  \item{eq}{Equation number. This can be ignored for univariate models.}

  \item{newdata}{A data frame or list containing the values of the model covariates at which predictions are required. 
                 This must have one row and must be provided.}

  \item{type}{The plot to produce, either \code{"surv"} or \code{"hazard"}.}

  \item{intervals}{If \code{TRUE} then intervals are also produced.}

  \item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used
             for interval calculations.}
  
  \item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.}  
   
  \item{shade}{If \code{TRUE} then it produces shaded regions as confidence bands.}
   
  \item{ylim, ylab, xlab}{Usual plot arguments.}
  
  \item{ls}{Length of sequence to use for time variable.}
  
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
