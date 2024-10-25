\name{k.tau}
\alias{k.tau}



\title{Kendall's tau a fitted joint model}

\description{ 
\code{k.tau} can be used to calculate the Kendall's tau from a fitted joint model with intervals obtained 
via posterior simulation.}



\usage{

k.tau(x, prob.lev = 0.05)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object as 
         produced by the respective fitting function.}

\item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.} 
  
}


\details{

This function calculates the Kendall's tau a fitted simultaneous model, with intervals 
obtained via posterior simulation. Note that this is derived under the assumption of continuous margins.

}



\value{

\item{res}{It returns the estimated tau with lower and upper interval limits.}
                        
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}  
}




\keyword{generalised joint regression modelling} 
\keyword{copula}
\keyword{Kendall's tau}
\keyword{marginal distribution}
\keyword{bayesian posterior simulation}





