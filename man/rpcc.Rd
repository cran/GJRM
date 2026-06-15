\name{rpcc}
\alias{rpcc}

\title{Simulation from a Fitted Trivariate PCC Copula Model}

\description{
\code{rpcc} generates Monte Carlo draws from a fitted trivariate pair-copula construction (PCC) model estimated via \code{gjrm.pcc}. 
The function simulates the three response variables jointly, using the fitted marginal distributions and the estimated PCC copula structure. 
Posterior simulation of model parameters is available for uncertainty propagation.
}

\usage{
rpcc(obj, newdata, n.draw = 1000, sim.par = FALSE, independent = FALSE)
}

\arguments{

\item{obj}{A fitted \code{gjrm.pcc} object containing a trivariate PCC model.}

\item{newdata}{A data frame with exactly one row, containing the covariate values at which the simulation is to be performed.}

\item{n.draw}{Number of Monte Carlo draws to generate. Must be greater than 2.}

\item{sim.par}{If \code{TRUE}, model parameters are drawn from their posterior distribution before simulation, allowing uncertainty in the 
              fitted model to be propagated. If \code{FALSE}, parameters are fixed at their estimated values.}

\item{independent}{If \code{TRUE}, the three margins are simulated independently using uniform random 
                   variables. If \code{FALSE}, dependence is generated using the fitted PCC copula structure.}

}

\details{


The function returns a data frame with simulated values for the three response variables. 
It is primarily intended for use inside higher-level Monte Carlo procedures such as \code{mc.copula.prob}.

}

\value{
A data frame with three columns:

\item{Y1}{Simulated values for the first response.}
\item{Y2}{Simulated values for the second response.}
\item{Y3}{Simulated values for the third response.}

Each column contains \code{n.sim} simulated draws.
}

\author{
Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}
}

\seealso{
\code{\link{mc.copula.prob}}, 
\code{\link{copula.prob}}, 
\code{\link{gjrm}}
}

\keyword{simulation}
\keyword{bayesian posterior simulation}
\keyword{copula}
\keyword{pair-copula construction}
