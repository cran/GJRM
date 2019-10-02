\name{resp.check}
\alias{resp.check}

\title{Plots for response variable}

\description{ 

It produces a histogram of the response along with the estimated density from the assumed distribution as well as a normal Q-Q plot for 
the (randomised) normalised quantile response. It also provides the log-likelihood for AIC calculation, for instance.

}

\usage{

resp.check(y, margin = "N", main = "Histogram and Density of Response",
           xlab = "Response", print.par = FALSE, plots = TRUE, 
           loglik = FALSE, os = FALSE,  
           intervals = FALSE, n.sim = 100, prob.lev = 0.05, 
           i.f = FALSE, ...)

}

\arguments{ 

\item{y}{Response.}
 
\item{margin}{The distributions allowed are: normal ("N"), log-normal ("LN"), generelised Pareto ("GP"), discrete generelised Pareto ("DGP"),
              Gumbel ("GU"), reverse Gumbel ("rGU"), logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"),
              Dagum ("DAGUM"), Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK"), Poisson ("PO"), zero truncated Poisson ("ZTP"), 
              negative binomial - type I ("NBI"), negative 
              binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG").} 

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{print.par}{If \code{TRUE} then the estimated parameters used to construct the plots are returned.}

\item{plots}{If \code{FALSE} then no plots are produced and only parameter estimates returned.}

\item{loglik}{If \code{TRUE} then it returns the logLik.}

\item{os}{If \code{TRUE} then the estimated parameters are returned on the original scale.}

\item{intervals}{If \code{TRUE} then intervals for the qqplot are produced.}   
  
\item{n.sim}{Number of replicate datasets used to simulate quantiles of the residual distribution.} 

\item{prob.lev}{Overall probability of the left and right tails of the probabilities' distribution used for interval calculations.} 

\item{i.f}{Internal fitting option. This is not for user purposes.}

\item{...}{Other graphics parameters to pass on to plotting commands.}

}


\details{

Prior to fitting a model with discrete and/or continuous margins, the distributions for the responses
may be chosen by looking at the histogram of the response along with the estimated density from the assumed distribution, and at the
normalised quantile responses. These will provide a rough guide to the adequacy of the chosen distribution.
The latter are defined as the quantile standard normal function of the cumulative distribution function of the response with scale and location
estimated by MLE. These should behave approximately as normally distributed variables (even though the original 
observations are not). Therefore, a normal Q-Q plot is appropriate here. 

If \code{loglik = TRUE} then this function also provides the log-likelihood for AIC calculation, for instance.

The shapiro test can also be performed. 

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{gjrm}}  
}




\keyword{generalised joint regression modelling} 
\keyword{regression} 
\keyword{Q-Q plot} 
\keyword{density plot}
\keyword{histogram}









