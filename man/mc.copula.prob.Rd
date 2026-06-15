\name{mc.copula.prob}
\alias{mc.copula.prob}

\title{Monte Carlo Copula Probabilities from a Fitted PCC Model}

\description{
\code{mc.copula.prob} computes probabilities of user-defined events
from a fitted PCC model using Monte Carlo simulation. 
Intervals are obtained via repeated posterior simulation of the model parameters.
}

\usage{
mc.copula.prob(obj, y.vec, newdata, event.fun, joint = TRUE,
               n.draw = 50000, n.sim = 100, prob.lev = 0.05)
}

\arguments{

\item{obj}{A fitted \code{gjrm.pcc} object as produced by the respective fitting function.}

\item{y.vec}{A vector containing the threshold values used inside \code{event.fun}.}

\item{newdata}{A data frame with one row, containing the covariate values at which the probability is to be evaluated.}

\item{event.fun}{A user-supplied function of the form \code{function(simY, y.vec)} returning a logical vector indicating whether 
                 the event of interest occurred for each simulated draw.}

\item{joint}{If \code{TRUE} then the calculation is done using the fitted joint model. If \code{FALSE} then 
             the calculation is done from univariate fits.}

\item{n.draw}{Number of Monte Carlo draws from the fitted model used to compute the point estimate.}

\item{n.sim}{Number of posterior simulations used to compute interval estimates.}

\item{prob.lev}{Overall probability mass in the two tails of the simulated distribution used to form the interval.}

}

\details{
The function first simulates \code{n.draw} draws from the fitted joint distribution
with fixed parameters to obtain the point estimate of the probability.

It then performs \code{n.sim} posterior simulations of the model parameters, and for each draw computes the corresponding probability.
The empirical distribution of these \code{n.sim} values is used to form interval estimates.

The event of interest is defined entirely by the user through \code{event.fun},
which receives the simulated responses and must return a logical vector.
This allows arbitrary probabilities of events to be evaluated.
}

\value{

A named numeric vector containing:

\item{pr}{The Monte Carlo estimate of the probability of the event.}

\item{2.5\%}{Lower interval bound based on the empirical bootstrap distribution.}

\item{97.5\%}{Upper interval bound based on the empirical bootstrap distribution.}


}

\author{
Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}
}

\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}, \code{\link{copula.prob}}
}

\keyword{generalised joint regression modelling} 
\keyword{copula}
\keyword{probability}
\keyword{bayesian posterior simulation}





