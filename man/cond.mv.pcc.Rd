\name{cond.mv.pcc}
\alias{cond.mv.pcc}

\title{Conditional Mean/Variance from a PCC-based Trivariate Copula Model}

\description{
Function \code{cond.mv.pcc} calculates conditional means or variances from a trivariate copula model built via pair-copula constructions (PCC). Credible/confidence intervals are obtained using posterior simulation.
}

\usage{
cond.mv.pcc(x, eq, y1 = NULL, y2 = NULL, y3 = NULL, newdata, fun = "mean", 
            n.sim = 100, prob.lev = 0.05, method = "adaptive")
}

\arguments{

\item{x}{A fitted \code{cond.mv.pcc} object, as produced by a trivariate PCC copula model fitting function.}

\item{eq}{Equation of interest. The function infers the conditional structure based on this input.}

\item{y1}{Value of the first response variable. Required if the equation or conditioning involves it.}

\item{y2}{Value of the second response variable. Required if the equation or conditioning involves it.}

\item{y3}{Value of the third response variable. Required if the equation or conditioning involves it.}

\item{newdata}{A data frame with one row providing the covariate values for prediction. This must be supplied.}

\item{fun}{Specifies the functional form to compute: either \code{"mean"} (default) or \code{"variance"}.}

\item{n.sim}{Number of posterior draws from the distribution of the model parameters used for simulating the conditional mean or variance.}

\item{prob.lev}{The tail probability used for interval estimation. For example, \code{0.05} corresponds to a 95\% credible interval.}

\item{method}{The numerical integration method to be used. Options are:
\code{"adaptive"} (default), which performs adaptive quadrature with fallback to Gauss-Legendre integration if necessary, or 
\code{"quad"} (which is slower and tends to be less reliable), a transformation-based quadrature approach designed for integrals with infinite or semi-infinite limits. This method does not implement 
truncation or fallback. \code{"quad"} is often used as sensitivity analysis tool.}
}

\details{
\code{cond.mv.pcc} extends \code{cond.mv} to trivariate copula models using pair-copula constructions (PCC). It computes conditional means or variances by integrating over relevant conditional distributions. Uncertainty intervals are estimated using simulated parameter values from the posterior distribution.
}

\value{
\item{res}{A numeric vector containing: the lower bound of the interval, the estimated conditional mean or variance, and the upper bound.}

\item{prob.lev}{The probability level used for constructing the interval.}

\item{sim.mv}{A numeric vector of simulated conditional means or variances across posterior draws, used to derive the interval.}
}

\author{
Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}
}

\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}, \code{\link{cond.mv}}
}

\keyword{copulae}
\keyword{pair-copula construction}
\keyword{trivariate model}
\keyword{conditional mean}
\keyword{conditional variance}
\keyword{bayesian posterior simulation}
