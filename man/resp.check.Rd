\name{resp.check}
\alias{resp.check}

\title{Diagnostic plot for a variable}

\description{ 

It produces a normal Q-Q plot for 
the (randomised) normalised quantile response. It also provides the log-likelihood for AIC calculation, for instance. It is also used for internal purposes.

}

\usage{

resp.check(y, margin = "N", print.par = FALSE, plots = TRUE, 
           loglik = FALSE, os = FALSE, i.f = FALSE, 
           min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999,
           left.trunc = 0)

}

\arguments{ 

\item{y}{Response.}
 
\item{margin}{The distributions allowed are: normal ("N"), log-normal ("LN"), generelised Pareto ("GP"), discrete generelised Pareto ("DGP"),
              Gumbel ("GU"), reverse Gumbel ("rGU"), logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"),
              Dagum ("DAGUM"), Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK"), Poisson ("P"), zero truncated Poisson ("ZTP"), 
              negative binomial - type I ("NBI"), negative 
              binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG").} 

\item{print.par}{If \code{TRUE} then the estimated parameters used to construct the plot are returned.}

\item{plots}{If \code{FALSE} then no plot is produced and only parameter estimates returned.}

\item{loglik}{If \code{TRUE} then it returns the logLik.}

\item{os}{If \code{TRUE} then the estimated parameters are returned on the original scale.}

\item{i.f}{Internal fitting option. This is not for user purposes.}

\item{min.dn, min.pr, max.pr}{Allowed minimum and maximum for estimated probabities and densities for parameter estimation.}

\item{left.trunc}{Value of truncation at left. Currently done for count distributions only.}



}


\details{

Prior to fitting a model with discrete and/or continuous margins, the distributions for the outcome variables
may be chosen by checking the normalised quantile responses. These will provide a rough guide to the adequacy of the chosen distribution.
The latter are defined as the quantile standard normal function of the cumulative distribution function of the response with scale and location
estimated by MLE. These should behave approximately as normally distributed variables (even though the original 
observations are not). Therefore, a normal Q-Q plot is appropriate here. 

If \code{loglik = TRUE} then this function also provides the log-likelihood for AIC calculation, for instance.

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









