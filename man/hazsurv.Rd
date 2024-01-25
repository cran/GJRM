\name{hazsurv}
\alias{hazsurv}

\title{Post-estimation calculation of hazard, cumulative hazard and survival functions}



\description{This function produces estimated values, intervals and plots for the hazard, cumulative hazard and survival functions.}

\usage{

hazsurv(x, eq, newdata, type = "surv", t.range = NULL, t.vec = NULL, 
        intervals = TRUE, n.sim = 100, prob.lev = 0.05, shade = FALSE, 
        bars = FALSE, ylim, ylab, xlab, pch, ls = 100, baseline = FALSE,
        min.dn = 1e-200, min.pr = 1e-200, max.pr = 1, plot.out = TRUE, 
        print.progress = TRUE, ...)



}

\arguments{

  \item{x}{A fitted \code{gamlss}/\code{gjrm} object.}

  \item{eq}{Equation number. This can be ignored for univariate models.}

  \item{newdata}{A data frame or list containing the values of the model covariates at which predictions are required. This must always be provided.
  
  For the individual survival/hazard/cumulative hazard function, the data frame must have one row containing the values of the model covariates corresponding to the individual of interest. 
  
  For the (sub-)population survival/hazard/cumulative hazard function, the data frame must have as many rows as there are individuals in the (sub-)population of interest. Each row must contain the values of the model covariates of the corresponding individual.}

  \item{type}{Either \code{"surv"}, \code{"hazard"} or \code{"cumhaz"}. In the excess hazard setting these are, respectively, the net survival, the excess hazard and the cumulative excess hazard.}
  
  \item{t.range}{Time variable range. This must be a vector with only two elements: the minimum and maximum of the time range. If \code{NULL} then it is determined automatically based on the observed data.}
  
  \item{t.vec}{Vector of time values. This can also be a single time. Note you cannot provide both \code{t.range} and \code{t.vec} as they are two mutually exclusive ways of defining the time variable. If \code{NULL} then it is determined automatically based on the observed data. }
  
  \item{intervals}{If \code{TRUE} then intervals are also produced.}

  \item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used
             for interval calculations.}
  
  \item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.}  
   
  \item{shade}{If \code{TRUE} then it produces shaded regions as confidence bands.}
  
  \item{bars}{If \code{TRUE} then the confidence intervals are plotted as bars rather than continuous curves. If \code{t.vec} is used and only one time value is provided, this is the only possible plotting option for the confidence intervals. Note \code{shade} and \code{bars} are mutually exclusive.}
   
  \item{ylim, ylab, xlab, pch}{Usual plot arguments.}
  
  \item{ls}{Length of sequence to use for time variable.}
  
  \item{baseline}{If baseline is desired; this will set all covariate/smooth effects to zero.}
  
  \item{min.dn, min.pr, max.pr}{Allowed minimum and maximum for estimated probabities and densities for survival, hazard and cumulative hazard calculations.}

  \item{plot.out}{If \code{FALSE} then the function does not produce a plot. The default is \code{TRUE}.}

  \item{print.progress}{If \code{FALSE} then the function does not print progress made. The default is \code{TRUE}.}
  
  \item{...}{Other arguments to pass to plot.}

}



\value{It produces estimated values, intervals and plots for the hazard, cumulative hazard and survival functions.}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\keyword{hazsurv}
\keyword{hplot} 
\keyword{smooth} 
\keyword{regression}
